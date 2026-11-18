"""Tier 1 master table and the 2-gate QC (Workflow step 11, `expression.py master`).

Tier 1 is pure raw: the SNV table (core_variants.tsv) joined with the AVI export and, per Atlas
track, a `{prefix}_raw` and `{prefix}_quantile` column. No aggregation happens here; the tissue
contrasts of Tier 2 (core/contrast.py) read this table.

Scores come from raw/scores_{scorer}.tsv and are matched to the track registry
(processed/track_registry.tsv) by (scorer, track_name, track_strand); a scorer that cohorts
dropped has no registry row, so its scores are left out. Gene-centric rows of another gene than
the locus's are left out too.

QC (two gates):
  Gate 1  coordinates and reference: expected row count and distinct positions, unique IDs,
          positions inside the window, valid ALT, REF equal to the window sequence and to the AVI
          REF, an AVI row for every SNV.
  Gate 2  attribution invariants: the AVI `top_modality` names the feature with the largest
          |fi_*| (`top_feature_key`), its value equals `top_feature_importance`, and the
          insertion/deletion attributions of an SNV are 0.

Outputs (processed/): master.parquet, track_scores_long.parquet, qc_report.json.

Library: MasterError, QcReport, build_master, run, load_master, load_registry, load_qc_report,
stage_metrics.
"""

from __future__ import annotations

import dataclasses
import json
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any

import polars as pl

import models
from core import locus as core_locus
from core import tracks as core_tracks
from core import variants as core_variants
from core.locus import LocusConfig
from models import AVI_FEATURES, AviRecord, TrackScoreRow
from study import load_decisions

WORKFLOW_STEP = 11
SCORES_STEP = 10
MASTER_FILE = "processed/master.parquet"
LONG_FILE = "processed/track_scores_long.parquet"
QC_FILE = "processed/qc_report.json"
REGISTRY_FILE = "processed/track_registry.tsv"
SCORES_GLOB = "scores_*.tsv"

# AVI display names of `top_modality` (alphagenome_atlas_avi.py AviFeature) -> feature key.
FEATURE_BY_DISPLAY = {
    "Splicing": "MERGED_SPLICING", "RNA-seq": "MAX_ABS_RNA_SEQ", "ATAC-seq": "MAX_ABS_ATAC",
    "DNASE-seq": "MAX_ABS_DNASE", "ChIP-TF": "MAX_ABS_CHIP_TF", "ChIP-Histone": "MAX_ABS_CHIP_HISTONE",
    "CAGE": "MAX_ABS_CAGE", "PRO-cap": "MAX_ABS_PROCAP", "Polyadenylation": "MAX_ABS_POLYADENYLATION",
    "3D Genome Contacts": "MAX_ABS_CONTACT_MAPS", "AlphaMissense": "ALPHAMISSENSE",
    "Cactus": "CACTUS_241_WAY", "Protein Termination": "PROTEIN_TERMINATION", "Start Lost": "START_LOST",
    "Stop Lost": "STOP_LOST", "PhastCons 470": "PHASTCONS_470_WAY", "Insertion": "IS_INSERTION",
    "Deletion": "IS_DELETION",
}
FI_COLUMNS = [f"fi_{k}" for k in AVI_FEATURES]
_LONG_COLUMNS = ["variant", "scorer", "track_index", "track_name", "track_strand", "ontology_curie",
                 "biosample_name", "biosample_type", "assay_title", "gtex_tissue", "data_source",
                 "gene_id", "gene_name", "gene_strand", "raw_score", "quantile_score", "column_prefix"]


class MasterError(Exception):
    """Invalid or incomplete Tier 1 input; the CLI prints it and exits with 1."""


@dataclasses.dataclass(frozen=True, kw_only=True)
class QcReport:
    """Result of the 2-gate QC.

    Attributes:
        gate1_pass: Gate 1 (coordinates and reference) passed.
        gate2_pass: Gate 2 (attribution invariants) passed.
        ref_mismatches: SNVs whose REF differs from the window sequence or from the AVI REF.
        missing_avi: SNVs without an AVI row.
        top_feature_mismatches: SNVs whose reported top feature is not the argmax of |fi_*|.
        indel_attributions: SNVs with a non-zero insertion or deletion attribution.
        total_variants: Rows checked.
        expected_variants: 3 x (window positions - excluded positions).
        distinct_positions: Distinct positions in the table.
        expected_positions: Window positions - excluded positions.
        out_of_bounds_positions: Rows outside the window.
        invalid_alt_alleles: Rows whose ALT is not a base different from REF.
        duplicated_variant_rows: Rows sharing a variant ID.
        quarantined_positions: Positions the researcher excluded.
        missing_track_scores: (SNV x registry track) cells without a score; informational.
        failed_variants: Variant IDs failing either gate.
    """

    gate1_pass: bool
    gate2_pass: bool
    ref_mismatches: int
    missing_avi: int
    top_feature_mismatches: int
    indel_attributions: int
    total_variants: int = 0
    expected_variants: int = 0
    distinct_positions: int = 0
    expected_positions: int = 0
    out_of_bounds_positions: int = 0
    invalid_alt_alleles: int = 0
    duplicated_variant_rows: int = 0
    quarantined_positions: tuple[int, ...] = ()
    missing_track_scores: int = 0
    failed_variants: tuple[str, ...] = ()

    @property
    def passed(self) -> bool:
        return self.gate1_pass and self.gate2_pass

    def to_dict(self, study_id: str) -> dict[str, Any]:
        return {
            "study_id": study_id,
            "overall_qc_passed": self.passed,
            "gates": [
                {"gate": "Gate 1: Coordinates & Reference", "passed": self.gate1_pass,
                 "total_variants": self.total_variants, "expected_variants": self.expected_variants,
                 "distinct_positions": self.distinct_positions, "expected_positions": self.expected_positions,
                 "out_of_bounds_positions": self.out_of_bounds_positions,
                 "ref_mismatches": self.ref_mismatches, "invalid_alt_alleles": self.invalid_alt_alleles,
                 "duplicated_variant_rows": self.duplicated_variant_rows, "missing_avi": self.missing_avi,
                 "quarantined_positions": list(self.quarantined_positions)},
                {"gate": "Gate 2: SNV Attribution Invariants", "passed": self.gate2_pass,
                 "top_feature_mismatches": self.top_feature_mismatches,
                 "indel_attributions": self.indel_attributions},
            ],
            "missing_track_scores": self.missing_track_scores,
            "failed_variants": list(self.failed_variants),
        }


# --- long table and wide pivot --------------------------------------------------------------


def long_scores(scores: Iterable[TrackScoreRow], registry: pl.DataFrame, variant_ids: Iterable[str],
                gene_id: str | None = None) -> pl.DataFrame:
    """Score rows matched to the registry, one row per (variant, registry track).

    Rows with an error, of a track outside the registry, of a variant outside `variant_ids` or of
    another gene than `gene_id` (versions ignored) are dropped.

    Raises:
        MasterError: two rows for the same (variant, track).
    """
    keep = set(variant_ids)
    rows = []
    for s in scores:
        if s.error or s.variant not in keep:
            continue
        if gene_id and s.gene_id and s.gene_id.split(".", 1)[0] != gene_id:
            continue
        rows.append({c: getattr(s, c, None) for c in _LONG_COLUMNS if c != "column_prefix"})
    schema = {c: (pl.Int64 if c == "track_index" else pl.Float64 if c.endswith("_score") else pl.String)
              for c in _LONG_COLUMNS if c != "column_prefix"}
    table = pl.DataFrame(rows, schema=schema)
    prefixes = registry.select("scorer", "track_name", "track_strand", "column_prefix")
    table = table.join(prefixes, on=["scorer", "track_name", "track_strand"], how="inner")
    dup = table.group_by("variant", "column_prefix").len().filter(pl.col("len") > 1)
    if dup.height:
        first = dup.row(0, named=True)
        raise MasterError(f"{dup.height} (variant, track) pair(s) have more than one score, e.g. "
                          f"{first['variant']} x {first['column_prefix']}; check raw/scores_*.tsv (Workflow step {SCORES_STEP}).")
    return table.select(_LONG_COLUMNS).sort("variant", "column_prefix")


def pivot_tracks(long_table: pl.DataFrame, registry: pl.DataFrame) -> pl.DataFrame:
    """One row per variant; `{prefix}_raw` then `{prefix}_quantile` per registry track, in registry order."""
    prefixes = registry["column_prefix"].to_list()
    raw = long_table.pivot(on="column_prefix", index="variant", values="raw_score")
    qnt = long_table.pivot(on="column_prefix", index="variant", values="quantile_score")
    raw = raw.rename({p: f"{p}_raw" for p in prefixes if p in raw.columns})
    qnt = qnt.rename({p: f"{p}_quantile" for p in prefixes if p in qnt.columns})
    wide = raw.join(qnt, on="variant", how="full", coalesce=True)
    ordered = ["variant"] + [f"{p}_{k}" for k in ("raw", "quantile") for p in prefixes if f"{p}_{k}" in wide.columns]
    return wide.select(ordered)


def avi_table(avi: Sequence[AviRecord]) -> pl.DataFrame:
    """AVI records as a table keyed by `variant_id` (the AVI variant string), with the top feature resolved."""
    rows = []
    for r in avi:
        key = FEATURE_BY_DISPLAY.get(r.top_modality)
        info = r.track.get(key) if key else None
        rows.append({
            "variant_id": r.variant, "avi_ref": r.ref, "avi_raw": r.avi_raw, "avi_quantile": r.avi_quantile,
            "avi_phred": r.avi_phred, "avi_top_percent": r.top_percentile,
            "top_feature_key": key, "top_feature": r.top_modality, "top_feature_value": r.top_feature_importance,
            "top_track_name": info.name if info else None, "top_track_biosample": info.biosample if info else None,
            "top_track_gene": info.gene if info else None,
            **{f"fi_{k}": float(v) for k, v in r.fi.items()},
        })
    schema = {"variant_id": pl.String, "avi_ref": pl.String, "top_feature_key": pl.String, "top_feature": pl.String,
              "top_track_name": pl.String, "top_track_biosample": pl.String, "top_track_gene": pl.String}
    return pl.DataFrame(rows, schema_overrides=schema)


# --- QC -------------------------------------------------------------------------------------


def run_qc(df: pl.DataFrame, locus: LocusConfig | None, quarantined: Iterable[int],
           missing_track_scores: int = 0) -> tuple[pl.DataFrame, QcReport]:
    """Appends `qc_pass_gate1` / `qc_pass_gate2` to the joined table and builds the report."""
    excluded = sorted(set(quarantined))
    if locus is not None:
        expected_positions = locus.window.length - len(excluded)
        in_window = pl.col("position").is_between(locus.window.start, locus.window.end)
    else:
        expected_positions = df["position"].n_unique()
        in_window = pl.lit(True)
    keys = pl.Series(list(AVI_FEATURES))
    fi_abs = pl.concat_list([pl.col(c).abs() for c in FI_COLUMNS])
    argmax = fi_abs.list.arg_max().map_elements(lambda i: keys[i], return_dtype=pl.String)
    is_snv = (pl.col("ref_allele").str.len_chars() == 1) & (pl.col("alt_allele").str.len_chars() == 1)
    checks = df.select(
        in_window.alias("in_window"),
        pl.col("avi_ref").is_not_null().alias("avi_ok"),
        (pl.col("ref_allele") == pl.col("avi_ref")).fill_null(True).alias("ref_ok"),
        (pl.col("alt_allele").is_in(list(core_variants.BASES)) & (pl.col("alt_allele") != pl.col("ref_allele"))).alias("alt_ok"),
        (~pl.col("variant_id").is_duplicated()).alias("unique"),
        (argmax == pl.col("top_feature_key")).fill_null(False).alias("top_ok"),
        ((pl.col("top_feature_value").abs() - pl.max_horizontal([pl.col(c).abs() for c in FI_COLUMNS])).abs() < 1e-9)
        .fill_null(False).alias("top_value_ok"),
        (~is_snv | ((pl.col("fi_IS_INSERTION") == 0) & (pl.col("fi_IS_DELETION") == 0))).fill_null(False).alias("indel_ok"),
    )
    g1 = checks.select(pl.all_horizontal("in_window", "avi_ok", "ref_ok", "alt_ok", "unique").alias("g1")).to_series()
    g2 = checks.select(pl.all_horizontal("top_ok", "top_value_ok", "indel_ok").alias("g2")).to_series()
    out = df.with_columns(g1.alias("qc_pass_gate1"), g2.alias("qc_pass_gate2"))
    total, distinct = df.height, df["position"].n_unique()
    gate1 = bool(g1.all()) and total == expected_positions * 3 and distinct == expected_positions
    report = QcReport(
        gate1_pass=gate1,
        gate2_pass=bool(g2.all()),
        ref_mismatches=int((~checks["ref_ok"]).sum()),
        missing_avi=int((~checks["avi_ok"]).sum()),
        top_feature_mismatches=int((~(checks["top_ok"] & checks["top_value_ok"])).sum()),
        indel_attributions=int((~checks["indel_ok"]).sum()),
        total_variants=total, expected_variants=expected_positions * 3, distinct_positions=distinct,
        expected_positions=expected_positions, out_of_bounds_positions=int((~checks["in_window"]).sum()),
        invalid_alt_alleles=int((~checks["alt_ok"]).sum()),
        duplicated_variant_rows=int((~checks["unique"]).sum()),
        quarantined_positions=tuple(excluded), missing_track_scores=missing_track_scores,
        failed_variants=tuple(out.filter(~(pl.col("qc_pass_gate1") & pl.col("qc_pass_gate2")))["variant_id"].to_list()),
    )
    return out, report


def build_master(
    variants: pl.DataFrame,
    avi: Sequence[AviRecord],
    scores: Iterable[TrackScoreRow],
    registry: pl.DataFrame,
    *,
    locus: LocusConfig | None = None,
    quarantined: Iterable[int] = (),
    gene_id: str | None = None,
) -> tuple[pl.DataFrame, QcReport]:
    """Tier 1 table (variants + AVI + wide track scores + QC flags) and the QC report.

    Args:
        variants: core_variants table (core.variants.enumerate_variants).
        avi: Records of raw/avi_region.parquet.
        scores: Rows of raw/scores_*.tsv.
        registry: processed/track_registry.tsv as a table.
        locus: Window for the Gate 1 bounds and counts; None skips the window checks.
        quarantined: Positions the researcher excluded (Gate 1 expects them absent).
        gene_id: Unversioned gene ID that gene-centric score rows must carry.

    Returns:
        (master, report); the long score table is `long_scores(...)` over the same inputs.
    """
    long_table = long_scores(scores, registry, variants["variant_id"].to_list(), gene_id)
    master, report, _ = _assemble(variants, avi, long_table, registry, locus, quarantined)
    return master, report


def _assemble(variants: pl.DataFrame, avi: Sequence[AviRecord], long_table: pl.DataFrame, registry: pl.DataFrame,
              locus: LocusConfig | None, quarantined: Iterable[int]) -> tuple[pl.DataFrame, QcReport, pl.DataFrame]:
    """Joins variants, AVI and the wide scores, then runs the QC.

    The AVI join is on position and ALT, not on the variant string, because that string carries the
    AVI REF: a REF disagreement then fails Gate 1 instead of looking like a missing AVI row.
    """
    avi_t = avi_table(avi).with_columns(
        pl.col("variant_id").str.extract(r":(\d+):", 1).cast(pl.Int64).alias("position"),
        pl.col("variant_id").str.extract(r">([ACGTN]+)$", 1).alias("alt_allele"),
    ).drop("variant_id")
    joined = (variants.join(avi_t, on=["position", "alt_allele"], how="left")
              .join(pivot_tracks(long_table, registry), left_on="variant_id", right_on="variant", how="left")
              .sort(["position", "alt_allele"]))
    missing = variants.height * registry.height - long_table.height
    flagged, report = run_qc(joined, locus, quarantined, missing_track_scores=missing)
    return flagged.drop("avi_ref"), report, long_table


# --- stage ----------------------------------------------------------------------------------


def load_registry(study_dir: Path) -> pl.DataFrame:
    """Reads processed/track_registry.tsv (written by the cohorts step)."""
    path = Path(study_dir) / REGISTRY_FILE
    if not path.is_file():
        raise MasterError(f"{REGISTRY_FILE} not found in {study_dir}; run expression.py cohorts (Workflow step 8).")
    df = pl.read_csv(path, separator="\t", infer_schema_length=None)
    missing = [c for c in core_tracks.REGISTRY_COLUMNS if c not in df.columns]
    if missing:
        raise MasterError(f"{REGISTRY_FILE} lacks the column(s) {missing} that the cohorts step writes; it was made by "
                          "another version of the skill or edited. Run expression.py cohorts again.")
    return df


def _load_scores(study_dir: Path, scorers: Sequence[str]) -> list[TrackScoreRow]:
    rows: list[TrackScoreRow] = []
    for scorer in scorers:
        path = study_dir / "raw" / f"scores_{scorer}.tsv"
        if not path.is_file():
            raise MasterError(f"raw/scores_{scorer}.tsv not found in {study_dir}; run Workflow step {SCORES_STEP} "
                              f"(atlas_query.py scores ... --scorers {scorer}).")
        rows.extend(models.load_raw(path))
    return rows


def run(study_dir: Path) -> QcReport:
    """Builds the Tier 1 table from core_variants.tsv, raw/avi_region.parquet and raw/scores_*.tsv.

    Writes master.parquet, track_scores_long.parquet and qc_report.json (also when QC fails, so the
    failing variants can be inspected).

    Raises:
        MasterError, VariantsError, core.locus.LocusError, models.DependencyOutputError, study.StudyError.
    """
    study_dir = Path(study_dir)
    decisions = load_decisions(study_dir)
    locus = core_locus.load_locus(study_dir)
    variants = core_variants.load_core_variants(study_dir)
    avi = core_variants._load_avi(study_dir)
    registry = load_registry(study_dir)
    scores = _load_scores(study_dir, sorted(set(registry["scorer"].to_list())))
    quarantined = core_variants.excluded_positions(decisions, include_rules=True)
    long_table = long_scores(scores, registry, variants["variant_id"].to_list(), locus.gene_id)
    master, report, long_table = _assemble(variants, avi, long_table, registry, locus, quarantined)
    (study_dir / "processed").mkdir(parents=True, exist_ok=True)
    master.write_parquet(study_dir / MASTER_FILE)
    long_table.write_parquet(study_dir / LONG_FILE)
    body = report.to_dict(locus.study_id) | {"rows": master.height, "columns": master.width,
                                             "n_long_scores": long_table.height}
    (study_dir / QC_FILE).write_text(json.dumps(body, indent=2) + "\n", encoding="utf-8")
    return report


def load_master(study_dir: Path) -> pl.DataFrame:
    """Reads processed/master.parquet."""
    path = Path(study_dir) / MASTER_FILE
    if not path.is_file():
        raise MasterError(f"{MASTER_FILE} not found in {study_dir}; run expression.py master (Workflow step {WORKFLOW_STEP}).")
    return pl.read_parquet(path)


def load_qc_report(study_dir: Path) -> dict[str, Any]:
    path = Path(study_dir) / QC_FILE
    if not path.is_file():
        raise MasterError(f"{QC_FILE} not found in {study_dir}; run expression.py master (Workflow step {WORKFLOW_STEP}).")
    return json.loads(path.read_text(encoding="utf-8"))


def stage_metrics(study_dir: Path) -> dict[str, int | float | str]:
    """Metrics of the master stage compared by `reproduce` (expected.json stages.master)."""
    study_dir = Path(study_dir)
    master = load_master(study_dir)
    qc = load_qc_report(study_dir)
    gate1, gate2 = qc["gates"]
    metrics: dict[str, int | float | str] = {
        "n_rows": master.height,
        "gate1_pass": "true" if gate1["passed"] else "false",
        "gate2_pass": "true" if gate2["passed"] else "false",
        "n_long_scores": qc["n_long_scores"],
    }
    for key, n in master["top_feature_key"].value_counts().iter_rows():
        metrics[f"top_feature.{key}"] = n
    return metrics
