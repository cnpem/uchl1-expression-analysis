"""Saturated SNV enumeration over the analysis window (Workflow step 11, `expression.py variants`).

Every window position gets the 3 SNVs of the other bases on the GRCh38 forward strand. REF comes
from the Ensembl window sequence (raw/sequence_*.txt); the Atlas AVI export (raw/avi_region.parquet)
is the universe the scores exist for, so the two are cross-checked and any disagreement becomes a
question instead of a silent fix. Codons are read in transcript orientation; a codon that leaves
exon 1 or the window keeps its CDS coordinates but no codon or amino-acid strings (those bases are
not in the window sequence).

Rule applied without asking (a technical check the researcher cannot judge position by position):
  variants.ref_mismatch:{pos}         the Atlas REF differs from the Ensembl sequence at pos
  variants.missing_in_atlas:{pos}     pos has SNVs without an AVI row
Each such position is excluded (`exclude_position`: all its SNVs leave every later table) and recorded
in decisions.yaml with decided_by: skill and the reason. When they exceed MAX_AUTO_EXCLUDED_FRACTION of
the window the stage stops instead: so many disagreements point to a capture problem (another window,
assembly or pull), not to isolated positions. An answer recorded earlier by the researcher
(`exclude_position` or `abort`) still applies. The site in a key is the 1-based position.

Output (processed/): core_variants.tsv (one row per SNV; `exon_number` is 1 inside exon 1 and empty
elsewhere because the GENCODE GTF carries exon numbers, not exon IDs).

Library: VariantsError, VariantsPending, enumerate_variants, excluded_positions, run,
load_core_variants, stage_metrics.
"""

from __future__ import annotations

import urllib.parse
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any

import polars as pl

import models
from core import locus as core_locus
from core.locus import LocusConfig
from models import AviRecord, EnsemblSequence
from study import DecisionEntry, Option, PendingDecision, load_decisions, load_study, sync_rules

WORKFLOW_STEP = 11
AVI_STEP = 9
AVI_FILE = "raw/avi_region.parquet"
CORE_VARIANTS_FILE = "processed/core_variants.tsv"
ATLAS_BASE_URL = "https://deepmind.google.com/science/alphagenome/atlas"
ENSEMBL_BASE_URL = "https://www.ensembl.org"
BASES = ("A", "C", "G", "T")
_COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}
_TRANSCRIPT_REGIONS = ("upstream_gene_variant", "5_prime_UTR_variant", "start_codon_variant",
                       "coding_sequence_variant")

CODON_TABLE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

CORE_SCHEMA: dict[str, Any] = {
    "chromosome": pl.String, "position": pl.Int64, "ref_allele": pl.String, "alt_allele": pl.String,
    "variant_id": pl.String, "ensembl_regulatory_id": pl.String, "ensembl_regulatory_biotype": pl.String,
    "exon_number": pl.Int64, "ensembl_transcript_region": pl.String, "genomic_domain": pl.String,
    "distance_to_tss": pl.Int64, "is_coding_cds": pl.Boolean, "is_start_codon": pl.Boolean,
    "atlas_url": pl.String, "ensembl_location_url": pl.String, "ensembl_regulatory_url": pl.String,
    "ensembl_transcript_url": pl.String, "cds_position": pl.Int64, "codon_number": pl.Int64,
    "codon_position": pl.Int64, "ref_codon": pl.String, "alt_codon": pl.String, "ref_aa": pl.String,
    "alt_aa": pl.String,
}
VARIANT_COLUMNS = tuple(CORE_SCHEMA)


class VariantsError(Exception):
    """Invalid or incomplete variant input; the CLI prints it and exits with 1."""


class VariantsPending(Exception):
    """Questions the researcher must answer before the variants can be enumerated.

    Attributes:
        items: The questions, for study.raise_pending.
    """

    def __init__(self, items: list[PendingDecision]) -> None:
        super().__init__(f"{len(items)} decision(s) needed")
        self.items = items


# --- annotation helpers ---------------------------------------------------------------------


def atlas_variant_url(variant_id: str) -> str:
    """`chr:pos:REF>ALT` -> Atlas variant page with the AVI section."""
    query = {"q": variant_id, "m": "variant", "lItems": "avi"}
    return f"{ATLAS_BASE_URL}?" + urllib.parse.urlencode(query, quote_via=urllib.parse.quote, safe=":,-")


def _species_path(species: str) -> str:
    return species[:1].upper() + species[1:]


def window_sequence(locus: LocusConfig, sequence: EnsemblSequence) -> str:
    """The window's bases (upper case) cut out of the raw sequence.

    Raises:
        VariantsError: the sequence is on another contig, does not cover the window, or holds a
            non-ACGT base inside it.
    """
    if sequence.chromosome != locus.chromosome:
        raise VariantsError(f"Sequence is on {sequence.chromosome}, the locus on {locus.chromosome}.")
    if sequence.start > locus.window.start or sequence.end < locus.window.end:
        raise VariantsError(
            f"Sequence {sequence.start}-{sequence.end} does not cover the window "
            f"{locus.window.start}-{locus.window.end}; run Workflow step 3 for a region that does.")
    offset = locus.window.start - sequence.start
    window = sequence.sequence[offset: offset + locus.window.length].upper()
    bad = sorted(set(window) - set(BASES))
    if bad:
        raise VariantsError(f"The window sequence holds non-ACGT bases {bad}; the SNV universe needs a resolved window.")
    return window


def _to_transcript(base: str, strand: int) -> str:
    return base if strand == 1 else _COMPLEMENT[base]


def _codon_fields(locus: LocusConfig, window: str, pos: int, alt: str) -> dict[str, Any]:
    s = locus.strand
    offset = (pos - locus.cds_start_genomic) * s  # 0-based CDS offset in transcript orientation
    codon_offset = offset - offset % 3
    genomic = [locus.cds_start_genomic + s * (codon_offset + i) for i in range(3)]
    fields: dict[str, Any] = {
        "cds_position": offset + 1, "codon_number": offset // 3 + 1, "codon_position": offset % 3 + 1,
        "ref_codon": None, "alt_codon": None, "ref_aa": None, "alt_aa": None,
    }
    if all(locus.first_exon.contains(g) and locus.window.contains(g) for g in genomic):
        ref_codon = "".join(_to_transcript(window[g - locus.window.start], s) for g in genomic)
        alt_codon = list(ref_codon)
        alt_codon[offset % 3] = _to_transcript(alt, s)
        alt_codon = "".join(alt_codon)
        fields |= {"ref_codon": ref_codon, "alt_codon": alt_codon,
                   "ref_aa": CODON_TABLE[ref_codon], "alt_aa": CODON_TABLE[alt_codon]}
    return fields


# --- decisions ------------------------------------------------------------------------------


MAX_AUTO_EXCLUDED_FRACTION = 0.01  # above 1% of the window's positions the stage asks first (capture problem?)
LIMIT_KEY = "variants.exclusion_limit"
RULE_PREFIXES = ("variants.ref_mismatch:", "variants.missing_in_atlas:")


def _limit_question(found: Sequence[dict[str, Any]], size: int, chromosome: str) -> PendingDecision:
    positions = sorted(f["pos"] for f in found)
    shown = ", ".join(f"{chromosome}:{p:,}" for p in positions[:10]) + (" …" if len(positions) > 10 else "")
    share = len(found) / size
    return PendingDecision(
        key=LIMIT_KEY,
        question=(f"At {len(found)} of the window's {size:,} positions ({share:.1%}) the AlphaGenome Atlas and the "
                  "reference sequence disagree, or the Atlas has no prediction. A few such positions are normal; this "
                  "many usually means the Atlas data and the sequence were captured for different windows. What should "
                  "happen?"),
        options=[Option(code="abort", label="Pause to capture the data again",
                        effect="The analysis pauses; the agent captures the Atlas data and the sequence again for the "
                               "study window and runs this step again."),
                 Option(code="exclude_all", label=f"Leave these {len(found)} positions out and continue",
                        effect=f"{3 * len(found):,} variants are not analysed; each position is listed in the report.")],
        recommended="abort", workflow_step=WORKFLOW_STEP,
        title="Many positions without matching predictions",
        why="Leaving out many positions removes variants from every result; if the cause is a capture problem, the "
            "rest of the data may be affected too.",
        identifiers={"positions": shown, "window positions": f"{size:,}"},
        agent_note="Before asking, check that raw/avi_region.parquet and the sequence file cover the window of "
                   "processed/locus.json (Workflow steps 3 and 9).")


def excluded_positions(decisions: Mapping[str, DecisionEntry], include_rules: bool = False) -> list[int]:
    """Positions excluded with `variants.ref_mismatch:{pos}` or `variants.missing_in_atlas:{pos}`.

    By default only the researcher's answers: the stage evaluates its rule again on every run, so a recapture
    that fixes a position brings it back. `include_rules` adds the positions the last run excluded by rule
    (decided_by: skill), i.e. every position missing from core_variants.tsv (what the master QC expects).
    """
    out: set[int] = set()
    for key, entry in decisions.items():
        if entry.decided_by == "skill" and not include_rules:
            continue
        kind, _, site = key.partition(":")
        if kind in ("variants.ref_mismatch", "variants.missing_in_atlas") and entry.value == "exclude_position":
            try:
                out.add(int(site))
            except ValueError:
                raise VariantsError(f"Decision key '{key}': the site must be a 1-based position.") from None
    return sorted(out)


def _answer(decisions: Mapping[str, DecisionEntry], key: str) -> str | None:
    entry = decisions.get(key)
    if entry is None or entry.decided_by == "skill":
        return None
    if entry.value == "exclude_position":
        return "exclude_position"
    if entry.value == "abort":
        raise VariantsError(f"Aborted by the researcher at '{key}'.")
    raise VariantsError(f"Decision '{key}' is {entry.value!r}; expected exclude_position or abort.")


# --- enumeration ----------------------------------------------------------------------------


def enumerate_variants(
    locus: LocusConfig,
    avi: Sequence[AviRecord],
    sequence: EnsemblSequence,
    *,
    decisions: Mapping[str, DecisionEntry] | None = None,
    species: str = "homo_sapiens",
    auto: list[dict[str, Any]] | None = None,
) -> pl.DataFrame:
    """All SNVs of the window (3 per position, REF from the sequence) minus the excluded positions.

    Args:
        locus: The window and annotation.
        avi: Records of raw/avi_region.parquet; they are cross-checked, not used as the universe.
        sequence: Ensembl sequence covering the window.
        decisions: Researcher answers; `variants.*` keys exclude positions.
        species: study.yaml species, for the Ensembl links.

        auto: When given, receives one {key, value, rule} per position excluded by the rule (the caller
            records them with study.record_rule).

    Returns:
        One row per SNV, columns in VARIANT_COLUMNS order, sorted by position then ALT.

    Raises:
        VariantsError: an aborting decision, more problem positions than MAX_AUTO_EXCLUDED_FRACTION of the
            window, or a sequence that cannot serve the window.
    """
    decisions = decisions or {}
    window = window_sequence(locus, sequence)
    have: dict[int, dict[str, str]] = {}
    for rec in avi:
        if locus.window.contains(rec.position):
            have.setdefault(rec.position, {})[rec.alt] = rec.ref
    excluded = set(excluded_positions(decisions))
    found: list[dict[str, Any]] = []
    for pos in range(locus.window.start, locus.window.end + 1):
        ref = window[pos - locus.window.start]
        seen = have.get(pos, {})
        wrong = sorted({r for r in seen.values() if r != ref})
        if wrong:
            if _answer(decisions, f"variants.ref_mismatch:{pos}") is None:
                found.append({"key": f"variants.ref_mismatch:{pos}", "pos": pos, "rule": (
                    f"At {locus.chromosome}:{pos:,} the reference sequence (Ensembl) has {ref} but the AlphaGenome Atlas "
                    f"has {', '.join(wrong)}; the position is left out, since its predictions refer to another base.")})
            continue
        absent = [a for a in BASES if a != ref and a not in seen]
        if absent and _answer(decisions, f"variants.missing_in_atlas:{pos}") is None:
            found.append({"key": f"variants.missing_in_atlas:{pos}", "pos": pos, "rule": (
                f"At {locus.chromosome}:{pos:,} ({ref}) the AlphaGenome Atlas has no prediction for "
                f"{', '.join(f'{ref}>{a}' for a in absent)}; the position is left out so every position is compared with "
                "all three possible base changes.")})
    size = locus.window.end - locus.window.start + 1
    if len(found) > MAX_AUTO_EXCLUDED_FRACTION * size:
        answer = decisions.get(LIMIT_KEY)
        accepted = answer is not None and answer.decided_by == "researcher" and str(answer.value) == "exclude_all"
        if answer is not None and str(answer.value) == "abort":
            raise VariantsError(f"Stopped by the researcher's answer to '{LIMIT_KEY}': capture the AlphaGenome Atlas "
                                "export (Workflow step 9) and the sequence (step 3) again for the study window.")
        if not accepted:
            raise VariantsPending([_limit_question(found, size, locus.chromosome)])
    for f in found:
        excluded.add(f["pos"])
        if auto is not None:
            auto.append({"key": f["key"], "value": "exclude_position", "rule": f["rule"]})

    transcript_url = f"{ENSEMBL_BASE_URL}/{_species_path(species)}/Transcript/Summary?t={locus.transcript_id}"
    base_url = f"{ENSEMBL_BASE_URL}/{_species_path(species)}"
    contig = locus.chromosome.removeprefix("chr")
    rows: list[dict[str, Any]] = []
    for pos in range(locus.window.start, locus.window.end + 1):
        if pos in excluded:
            continue
        ref = window[pos - locus.window.start]
        feature = core_locus.regulatory_feature_at(locus, pos)
        coding = core_locus.in_cds(locus, pos)
        for alt in (b for b in BASES if b != ref):
            vid = f"{locus.chromosome}:{pos}:{ref}>{alt}"
            row: dict[str, Any] = {
                "chromosome": locus.chromosome, "position": pos, "ref_allele": ref, "alt_allele": alt,
                "variant_id": vid,
                "ensembl_regulatory_id": feature.id if feature else None,
                "ensembl_regulatory_biotype": feature.biotype if feature else None,
                "exon_number": locus.first_exon.exon_number if locus.first_exon.contains(pos) else None,
                "ensembl_transcript_region": core_locus.transcript_region(locus, pos),
                "genomic_domain": core_locus.genomic_domain(locus, pos),
                "distance_to_tss": locus.distance_to_tss(pos),
                "is_coding_cds": coding,
                "is_start_codon": core_locus.in_start_codon(locus, pos),
                "atlas_url": atlas_variant_url(vid),
                "ensembl_location_url": f"{base_url}/Location/View?r={contig}:{max(1, pos - 25)}-{pos + 25}",
                "ensembl_regulatory_url": f"{base_url}/Regulation/Summary?rf={feature.id}" if feature else None,
                "ensembl_transcript_url": transcript_url,
                "cds_position": None, "codon_number": None, "codon_position": None,
                "ref_codon": None, "alt_codon": None, "ref_aa": None, "alt_aa": None,
            }
            if coding:
                row |= _codon_fields(locus, window, pos, alt)
            rows.append(row)
    return pl.DataFrame(rows, schema=CORE_SCHEMA).sort(["position", "alt_allele"])


# --- stage ----------------------------------------------------------------------------------


def _load_avi(study_dir: Path) -> list[AviRecord]:
    path = study_dir / AVI_FILE
    if not path.is_file():
        raise VariantsError(f"{AVI_FILE} not found in {study_dir}; run Workflow step {AVI_STEP} "
                            "(alphagenome_atlas_avi.py region ...).")
    return list(models.load_raw(path))


def run(study_dir: Path) -> pl.DataFrame:
    """Enumerates the variants from locus.json, raw/avi_region.parquet and the window sequence; writes core_variants.tsv.

    Raises:
        VariantsPending, VariantsError, core.locus.LocusError, models.DependencyOutputError, study.StudyError.
    """
    study_dir = Path(study_dir)
    study = load_study(study_dir)
    decisions = load_decisions(study_dir)
    locus = core_locus.load_locus(study_dir)
    if not locus.sequence_file:
        raise VariantsError("locus.json records no window sequence file; run expression.py locus (Workflow step 4).")
    sequence = models.load_raw(study_dir / locus.sequence_file)[0]
    avi = _load_avi(study_dir)
    auto: list[dict[str, Any]] = []
    try:
        df = enumerate_variants(locus, avi, sequence, decisions=decisions, species=study.data["species"], auto=auto)
    except (VariantsPending, VariantsError):
        sync_rules(study_dir, RULE_PREFIXES, [], WORKFLOW_STEP)  # the last run's exclusions no longer hold
        raise
    sync_rules(study_dir, RULE_PREFIXES, auto, WORKFLOW_STEP)
    (study_dir / "processed").mkdir(parents=True, exist_ok=True)
    df.write_csv(study_dir / CORE_VARIANTS_FILE, separator="\t")
    return df


def load_core_variants(study_dir: Path) -> pl.DataFrame:
    """Reads processed/core_variants.tsv with its declared column types."""
    path = Path(study_dir) / CORE_VARIANTS_FILE
    if not path.is_file():
        raise VariantsError(f"{CORE_VARIANTS_FILE} not found in {study_dir}; run expression.py variants "
                            f"(Workflow step {WORKFLOW_STEP}).")
    return pl.read_csv(path, separator="\t", schema_overrides=CORE_SCHEMA)


def stage_metrics(study_dir: Path) -> dict[str, int | float | str]:
    """Metrics of the variants stage compared by `reproduce` (expected.json stages.variants)."""
    study_dir = Path(study_dir)
    df = load_core_variants(study_dir)
    metrics: dict[str, int | float | str] = {
        "n_rows": df.height,
        "n_positions": df["position"].n_unique(),
        "distance_to_tss_min": int(df["distance_to_tss"].min()),
        "distance_to_tss_max": int(df["distance_to_tss"].max()),
    }
    counts = dict(df["ensembl_transcript_region"].value_counts().iter_rows())
    for region in sorted(set(_TRANSCRIPT_REGIONS) | set(counts)):
        metrics[f"region.{region}"] = int(counts.get(region, 0))
    ref_at = dict(df.select("position", "ref_allele").unique().iter_rows())
    metrics["ref_mismatches"] = sum(1 for r in _load_avi(study_dir) if r.position in ref_at and r.ref != ref_at[r.position])
    return metrics
