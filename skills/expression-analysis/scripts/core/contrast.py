"""Tier 2 target-vs-control contrasts (Workflow step 12, `expression.py contrast [--propose]`).

The contrast plan states which tracks enter each contrast and how they are aggregated. The agent
proposes it (`contrast --propose` -> processed/contrasts_proposed.yaml), the researcher approves it
and it is saved as contrasts.yaml; `contrast` then executes it on the Tier 1 table (read-only):
each tissue value is the mean (or median) of its selected tracks, and cohort statistics are over
tissues.

Researcher questions (answered with study.py decide, Workflow step 12):
  contrast.review:mixed_assays:{contrast}          selected tracks mix assay titles:
                                                   accept | abort
  contrast.review:mixed_sources:{contrast}         selected tracks mix data sources:
                                                   accept | abort
  contrast.review:single_tissue:{contrast}:{role}  a cohort of the contrast has one tissue:
                                                   accept | abort
  contrast.plan                                    approval of the plan itself: approved
`abort` stops the stage so the plan can be edited (e.g. `track_filter`). Without `contrast.plan`
the stage exits 1 (`study.require`).

Priority metrics (threshold-free), from tissue-level quantiles of the plan's `direction_contrast`
target tissues and the control tissues of `scope_contrast` (default: the same contrast):
  predicted_direction     Increase / Decrease when every target tissue has that sign, else Inconsistent
  specificity_scope       control tissues voting against that direction: Shared / Opposite / Mixed,
                          "No scope" when the direction is inconsistent or no control tissue has
                          a quantile (no control data)
  top variants            the first `report.top_n` rows ordered by consistent direction first, then
                          the largest |target - control| difference of the direction contrast

Outputs (processed/): contrasts_proposed.yaml (--propose); tissue_contrast.parquet,
tissue_contrast_meta.json, top_variants.tsv, top_variants.json (contrast).

Library: ContrastError, ContrastPending, TrackFilter, ContrastSpec, PriorityConfig, ContrastPlan,
propose_plan, review_points, apply_plan, top_variants, propose, run, stage_metrics.
"""

from __future__ import annotations

import dataclasses
import json
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import polars as pl
import yaml

from core import master as core_master
from study import PendingDecision, load_decisions, load_study, require

WORKFLOW_STEP = 12
GENE_EXPRESSION_SCORER = "RNA_SEQ"
PLAN_FILE = "contrasts.yaml"
PROPOSED_FILE = "processed/contrasts_proposed.yaml"
CONTRAST_FILE = "processed/tissue_contrast.parquet"
META_FILE = "processed/tissue_contrast_meta.json"
TOP_TSV_FILE = "processed/top_variants.tsv"
TOP_JSON_FILE = "processed/top_variants.json"
PLAN_APPROVED = ("approved", "accept")
STRANDS = ("sense", "antisense", "unstranded", "any")
AGGREGATIONS = ("mean", "median")
FILTER_FIELDS = ("assay_title", "data_source", "biosample_type")
TOP_FIELDS = ("biosample", "curie", "obo_url", "ols_url", "atlas_url")
_NAME_RE = re.compile(r"^[a-z][a-z0-9_]*$")

DIRECTION_UP = "Increase"
DIRECTION_DOWN = "Decrease"
DIRECTION_INCONSISTENT = "Inconsistent (within target)"
SCOPE_SHARED = "Shared (control same direction)"
SCOPE_OPPOSITE = "Opposite in control"
SCOPE_MIXED = "Mixed control"
SCOPE_NOT_CALLED = "No scope (inconsistent direction)"
SCOPE_NO_CONTROL = "No scope (no control data)"


class ContrastError(Exception):
    """Invalid plan or contrast input; the CLI prints it and exits with 1."""


class ContrastPending(Exception):
    """Review questions the researcher must answer before the plan can run.

    Attributes:
        items: The questions, for study.raise_pending.
    """

    def __init__(self, items: list[PendingDecision]) -> None:
        super().__init__(f"{len(items)} decision(s) needed")
        self.items = items


# --- plan -----------------------------------------------------------------------------------


def _unknown(where: str, d: Mapping[str, Any], allowed: Sequence[str]) -> None:
    extra = sorted(set(d) - set(allowed))
    if extra:
        raise ContrastError(f"{where}: unknown field '{extra[0]}' (allowed: {', '.join(allowed)}).")


@dataclasses.dataclass(frozen=True, kw_only=True)
class TrackFilter:
    """Optional allow-lists on registry attributes; None means no restriction.

    Attributes:
        assay_title: Allowed assay titles.
        data_source: Allowed data sources.
        biosample_type: Allowed biosample types.
    """

    assay_title: tuple[str, ...] | None = None
    data_source: tuple[str, ...] | None = None
    biosample_type: tuple[str, ...] | None = None

    def accepts(self, track: Mapping[str, Any]) -> bool:
        return all(allowed is None or track.get(f) in allowed for f in FILTER_FIELDS
                   for allowed in (getattr(self, f),))

    def to_dict(self) -> dict[str, list[str]]:
        return {f: list(getattr(self, f)) for f in FILTER_FIELDS if getattr(self, f) is not None}

    @classmethod
    def from_dict(cls, d: Mapping[str, Any] | None, *, where: str) -> TrackFilter:
        d = d or {}
        if not isinstance(d, Mapping):
            raise ContrastError(f"{where}: track_filter must be a mapping.")
        _unknown(f"{where}.track_filter", d, FILTER_FIELDS)
        return cls(**{k: (tuple(map(str, v)) if v is not None else None) for k, v in d.items()})


@dataclasses.dataclass(frozen=True, kw_only=True)
class ContrastSpec:
    """One target-vs-control contrast.

    Attributes:
        name: Column prefix of the contrast (`{name}_target_mean_raw`, ...).
        scorer: Atlas scorer.
        strand: Track strand class: sense, antisense, unstranded or any.
        track_filter: Extra registry allow-lists.
        target: Target tissue CURIEs, in plan order.
        control: Control tissue CURIEs, in plan order.
        tissue_aggregation: How a tissue's tracks are combined: mean or median.
    """

    name: str
    scorer: str
    strand: str = "any"
    track_filter: TrackFilter = dataclasses.field(default_factory=TrackFilter)
    target: tuple[str, ...]
    control: tuple[str, ...]
    tissue_aggregation: str = "mean"

    def __post_init__(self) -> None:
        if not _NAME_RE.match(self.name):
            raise ContrastError(f"contrast name must match {_NAME_RE.pattern}: {self.name!r}")
        if self.strand not in STRANDS:
            raise ContrastError(f"{self.name}: strand must be one of {STRANDS}, got {self.strand!r}.")
        if self.tissue_aggregation not in AGGREGATIONS:
            raise ContrastError(f"{self.name}: tissue_aggregation must be one of {AGGREGATIONS}.")
        if not self.target or not self.control:
            raise ContrastError(f"{self.name}: target and control need at least one CURIE each.")
        both = sorted(set(self.target) & set(self.control))
        if both:
            raise ContrastError(f"{self.name}: CURIEs in both target and control: {both}")

    def tracks(self, registry: Sequence[Mapping[str, Any]], role: str) -> dict[str, list[Mapping[str, Any]]]:
        """CURIE -> selected registry tracks, in plan CURIE order."""
        return {c: [t for t in registry
                    if t["scorer"] == self.scorer and t["ontology_curie"] == c
                    and (self.strand == "any" or t["strand_class"] == self.strand)
                    and self.track_filter.accepts(t)]
                for c in getattr(self, role)}

    def to_dict(self) -> dict[str, Any]:
        return {"name": self.name, "scorer": self.scorer, "strand": self.strand,
                "track_filter": self.track_filter.to_dict(), "target": list(self.target),
                "control": list(self.control), "tissue_aggregation": self.tissue_aggregation}


@dataclasses.dataclass(frozen=True, kw_only=True)
class PriorityConfig:
    """Which contrasts define the priority labels.

    Attributes:
        direction_contrast: Contrast whose target tissue quantiles define predicted_direction.
        scope_contrast: Contrast whose control tissues vote on specificity_scope (None = direction_contrast).
    """

    direction_contrast: str
    scope_contrast: str | None = None

    def to_dict(self) -> dict[str, Any]:
        out: dict[str, Any] = {"direction_contrast": self.direction_contrast}
        if self.scope_contrast is not None:
            out["scope_contrast"] = self.scope_contrast
        return out


@dataclasses.dataclass(frozen=True, kw_only=True)
class ContrastPlan:
    """The approved set of contrasts (contrasts.yaml).

    Attributes:
        study_id: Study the plan belongs to.
        contrasts: The contrasts, in column order.
        priority: Priority-label configuration.
    """

    study_id: str
    contrasts: tuple[ContrastSpec, ...]
    priority: PriorityConfig

    def __post_init__(self) -> None:
        names = [c.name for c in self.contrasts]
        if not names:
            raise ContrastError("The plan has no contrasts.")
        if len(set(names)) != len(names):
            raise ContrastError(f"duplicate contrast names: {names}")
        for ref in (self.priority.direction_contrast, self.priority.scope_contrast):
            if ref is not None and ref not in names:
                raise ContrastError(f"priority refers to unknown contrast {ref!r}")

    def contrast(self, name: str) -> ContrastSpec:
        return next(c for c in self.contrasts if c.name == name)

    def to_dict(self) -> dict[str, Any]:
        return {"study_id": self.study_id, "contrasts": [c.to_dict() for c in self.contrasts],
                "priority": self.priority.to_dict()}

    def validate_against(self, registry: Sequence[Mapping[str, Any]]) -> None:
        """Checks every scorer, strand, filter and CURIE against the track registry.

        Raises:
            ContrastError: a contrast selects no track for a CURIE or names an unknown scorer.
        """
        for spec in self.contrasts:
            if not any(t["scorer"] == spec.scorer for t in registry):
                raise ContrastError(f"{spec.name}: scorer {spec.scorer} has no track in the registry.")
            for role in ("target", "control"):
                empty = [c for c, ts in spec.tracks(registry, role).items() if not ts]
                if empty:
                    raise ContrastError(f"{spec.name}: no registry track for {role} CURIE(s) {empty} after strand/filters.")

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> ContrastPlan:
        if not isinstance(d, Mapping):
            raise ContrastError("The plan must be a mapping.")
        _unknown("plan", d, ("study_id", "contrasts", "priority"))
        for key in ("study_id", "contrasts", "priority"):
            if key not in d:
                raise ContrastError(f"plan: missing field '{key}'.")
        specs = []
        for i, c in enumerate(d["contrasts"] or []):
            where = f"contrasts[{i}]"
            if not isinstance(c, Mapping):
                raise ContrastError(f"{where}: must be a mapping.")
            _unknown(where, c, ("name", "scorer", "strand", "track_filter", "target", "control", "tissue_aggregation"))
            for key in ("name", "scorer", "target", "control"):
                if key not in c:
                    raise ContrastError(f"{where}: missing field '{key}'.")
            specs.append(ContrastSpec(
                name=str(c["name"]), scorer=str(c["scorer"]), strand=str(c.get("strand", "any")),
                track_filter=TrackFilter.from_dict(c.get("track_filter"), where=where),
                target=tuple(map(str, c["target"])), control=tuple(map(str, c["control"])),
                tissue_aggregation=str(c.get("tissue_aggregation", "mean"))))
        pr = d["priority"]
        if not isinstance(pr, Mapping) or "direction_contrast" not in pr:
            raise ContrastError("plan: priority.direction_contrast is required.")
        _unknown("priority", pr, ("direction_contrast", "scope_contrast"))
        return cls(study_id=str(d["study_id"]), contrasts=tuple(specs),
                   priority=PriorityConfig(direction_contrast=str(pr["direction_contrast"]),
                                           scope_contrast=pr.get("scope_contrast")))

    @classmethod
    def from_yaml(cls, path: str | Path) -> ContrastPlan:
        """Reads a plan; raises ContrastError (also for unreadable YAML)."""
        path = Path(path)
        try:
            data = yaml.safe_load(path.read_text(encoding="utf-8"))
        except (OSError, yaml.YAMLError) as exc:
            raise ContrastError(f"{path}: cannot read the plan ({exc}).") from None
        try:
            return cls.from_dict(data)
        except ContrastError as exc:
            raise ContrastError(f"{path}: {exc}") from None

    def to_yaml(self, header: str = "", labels: Mapping[str, str] | None = None) -> str:
        """The plan as YAML; `labels` (CURIE -> tissue name) become end-of-line comments."""
        body = yaml.safe_dump(self.to_dict(), sort_keys=False, allow_unicode=True)
        if labels:
            lines = []
            for line in body.splitlines():
                curie = line.strip().removeprefix("- ")
                lines.append(f"{line}  # {labels[curie]}" if line.strip().startswith("- ") and curie in labels else line)
            body = "\n".join(lines) + "\n"
        return header + body


# --- proposal and review --------------------------------------------------------------------


def _stranded(registry: Sequence[Mapping[str, Any]], scorer: str) -> bool:
    return any(t["strand_class"] != "unstranded" for t in registry if t["scorer"] == scorer)


def propose_plan(registry: pl.DataFrame, master: pl.DataFrame, *, resolution: Mapping[str, Any] | None = None,
                 study_id: str | None = None) -> tuple[ContrastPlan, dict[str, str]]:
    """Default plan and the scorers it could not include (scorer -> reason).

    One contrast per scorer; stranded scorers use sense-strand tracks; target and control CURIEs
    are the effective cohorts that have tracks for that scorer; every track is kept; tissues are
    averaged. The RNA_SEQ contrast (else the first) defines the priority labels. The order of
    scorers and CURIEs follows `resolution` (processed/ontology_resolution.json) when given, else
    the registry.

    Args:
        registry: processed/track_registry.tsv.
        master: Tier 1 table; every selected track must have its `{prefix}_raw` column.
        resolution: ontology_resolution.json content (effective scorers and CURIEs, study_id).
        study_id: Study ID when there is no resolution.

    Raises:
        ContrastError: no scorer has tracks in both cohorts, or a selected track is missing from the master.
    """
    rows = registry.to_dicts()
    if resolution is not None:
        scorers = list(resolution["effective_scorers"])
        target_all = list(resolution["effective_target_curies"])
        control_all = list(resolution["effective_control_curies"])
        study_id = resolution["study_id"]
    else:
        scorers = list(dict.fromkeys(r["scorer"] for r in rows))
        target_all = list(dict.fromkeys(r["ontology_curie"] for r in rows if r["cohort"] == "target"))
        control_all = list(dict.fromkeys(r["ontology_curie"] for r in rows if r["cohort"] == "control"))
    if study_id is None:
        raise ContrastError("propose_plan needs a resolution or a study_id.")
    specs: list[ContrastSpec] = []
    skipped: dict[str, str] = {}
    for scorer in scorers:
        strand = "sense" if _stranded(rows, scorer) else "unstranded"
        have = {t["ontology_curie"] for t in rows if t["scorer"] == scorer and t["strand_class"] == strand}
        target = [c for c in target_all if c in have]
        control = [c for c in control_all if c in have]
        if not target or not control:
            skipped[scorer] = f"no {strand} tracks in the {'target' if not target else 'control'} cohort"
            continue
        name = scorer.lower() + ("_sense" if strand == "sense" else "")
        specs.append(ContrastSpec(name=name, scorer=scorer, strand=strand, target=tuple(target), control=tuple(control)))
    if not specs:
        raise ContrastError("No contrast can be proposed: no scorer has tracks in both cohorts.")
    direction = next((s.name for s in specs if s.scorer == GENE_EXPRESSION_SCORER), specs[0].name)
    plan = ContrastPlan(study_id=study_id, contrasts=tuple(specs), priority=PriorityConfig(direction_contrast=direction))
    plan.validate_against(rows)
    missing = sorted({t["column_prefix"] for s in plan.contrasts for role in ("target", "control")
                      for ts in s.tracks(rows, role).values() for t in ts
                      if f"{t['column_prefix']}_raw" not in master.columns})
    if missing:
        raise ContrastError(f"Tier 1 has no column for track(s) {', '.join(missing)}; run expression.py master again.")
    return plan, skipped


def review_points(plan: ContrastPlan, registry: Sequence[Mapping[str, Any]]) -> list[PendingDecision]:
    """Questions about choices embedded in the plan, recomputed from the plan about to run."""
    labels = {t["ontology_curie"]: t.get("biosample_name") for t in registry}
    out: list[PendingDecision] = []
    for spec in plan.contrasts:
        selected = {role: spec.tracks(registry, role) for role in ("target", "control")}
        tracks = [t for role in selected.values() for ts in role.values() for t in ts]
        for field, kind, what in (("assay_title", "mixed_assays", "assay titles"),
                                  ("data_source", "mixed_sources", "data sources")):
            values = sorted({t.get(field) or "unknown" for t in tracks})
            if len(values) > 1:
                out.append(PendingDecision(
                    key=f"contrast.review:{kind}:{spec.name}",
                    question=(f"Contrast {spec.name} ({spec.scorer}) mixes {what} {values}. Accept, or abort and set "
                              f"`track_filter.{field}` in the plan?"),
                    options=["accept", "abort"], recommended="accept", workflow_step=WORKFLOW_STEP))
        for role, curies in selected.items():
            if len(curies) == 1:
                (curie,) = curies
                out.append(PendingDecision(
                    key=f"contrast.review:single_tissue:{spec.name}:{role}",
                    question=(f"Contrast {spec.name} ({spec.scorer}) has a single {role} tissue "
                              f"({labels.get(curie)}, {curie}); its {role} mean/max/min are that tissue's value. Accept?"),
                    options=["accept", "abort"], recommended="accept", workflow_step=WORKFLOW_STEP))
    return out


# --- aggregation ----------------------------------------------------------------------------


def tissue_value(prefixes: Sequence[str], kind: str = "raw", how: str = "mean") -> pl.Expr:
    cols = [pl.col(f"{p}_{kind}") for p in prefixes]
    return pl.mean_horizontal(cols) if how == "mean" else pl.concat_list(cols).list.median()


def tissue_columns(spec: ContrastSpec, registry: Sequence[Mapping[str, Any]], role: str) -> dict[str, list[str]]:
    """CURIE -> selected column prefixes for one cohort of a contrast."""
    return {c: [t["column_prefix"] for t in ts] for c, ts in spec.tracks(registry, role).items()}


def _add_top_tissue(df: pl.DataFrame, registry: Sequence[Mapping[str, Any]], spec: ContrastSpec) -> pl.DataFrame:
    """`top_{role}_tissue_*`: tissue with the largest |tissue raw value| in the direction contrast."""
    first: dict[str, Mapping[str, Any]] = {}
    for t in registry:
        first.setdefault(t["ontology_curie"], t)
    cols: list[pl.Expr] = []
    for role in ("target", "control"):
        groups = tissue_columns(spec, registry, role)
        curies = list(groups)
        values = pl.concat_list([tissue_value(p, "raw", spec.tissue_aggregation) for p in groups.values()])
        idx = values.list.eval(pl.element().abs()).list.arg_max()
        cols.append(values.list.get(idx).alias(f"top_{role}_tissue_raw_score"))
        meta = {c: {"biosample": first[c]["biosample_name"], "curie": c, "obo_url": first[c]["ontology_obo_url"],
                    "ols_url": first[c]["ontology_ols_url"], "atlas_url": first[c]["atlas_locus_track_url"]}
                for c in curies}
        for field in TOP_FIELDS:
            lookup = pl.Series([meta[c][field] for c in curies], dtype=pl.String)
            cols.append(idx.map_elements(lambda i, s=lookup: s[i], return_dtype=pl.String).alias(f"top_{role}_tissue_{field}"))
    return df.with_columns(cols)


def add_priority_metrics(df: pl.DataFrame, registry: Sequence[Mapping[str, Any]], plan: ContrastPlan) -> pl.DataFrame:
    """`predicted_direction`, `target_min_abs_quantile`, `specificity_scope` from tissue-level quantiles."""
    d_spec = plan.contrast(plan.priority.direction_contrast)
    s_spec = plan.contrast(plan.priority.scope_contrast or plan.priority.direction_contrast)
    q_target = [tissue_value(p, "quantile", d_spec.tissue_aggregation)
                for p in tissue_columns(d_spec, registry, "target").values()]
    q_control = [tissue_value(p, "quantile", s_spec.tissue_aggregation)
                 for p in tissue_columns(s_spec, registry, "control").values()]
    direction = (pl.when(pl.all_horizontal([q > 0 for q in q_target])).then(pl.lit(DIRECTION_UP))
                 .when(pl.all_horizontal([q < 0 for q in q_target])).then(pl.lit(DIRECTION_DOWN))
                 .otherwise(pl.lit(DIRECTION_INCONSISTENT)))
    df = df.with_columns(direction.alias("predicted_direction"),
                         pl.min_horizontal([q.abs() for q in q_target]).alias("target_min_abs_quantile"))
    sign = (pl.when(pl.col("predicted_direction") == DIRECTION_UP).then(1)
            .when(pl.col("predicted_direction") == DIRECTION_DOWN).then(-1).otherwise(None))
    same = pl.sum_horizontal([(q * sign > 0).cast(pl.Int32) for q in q_control])
    opposite = pl.sum_horizontal([(q * sign < 0).cast(pl.Int32) for q in q_control])
    # A null quantile votes neither way; with no control quantile at all there is nothing to call.
    any_control = pl.any_horizontal([q.is_not_null() for q in q_control]) if q_control else pl.lit(False)
    return df.with_columns(
        pl.when(sign.is_null()).then(pl.lit(SCOPE_NOT_CALLED))
        .when(~any_control).then(pl.lit(SCOPE_NO_CONTROL))
        .when(opposite == 0).then(pl.lit(SCOPE_SHARED))
        .when(opposite > same).then(pl.lit(SCOPE_OPPOSITE))
        .otherwise(pl.lit(SCOPE_MIXED)).alias("specificity_scope"))


def apply_plan(master: pl.DataFrame, plan: ContrastPlan, registry: pl.DataFrame) -> pl.DataFrame:
    """Tier 1 columns + per-contrast target/control mean, max, min and difference + priority labels.

    Adds `{contrast}_{target|control}_{mean|max|min}_raw`, `{contrast}_target_vs_control_diff`,
    `top_{target|control}_tissue_*` (direction contrast), `predicted_direction`,
    `target_min_abs_quantile` and `specificity_scope`.

    Args:
        master: Tier 1 table (read-only).
        plan: Approved plan.
        registry: processed/track_registry.tsv (the plan selects tracks from it).
    """
    rows = registry.to_dicts()
    plan.validate_against(rows)
    cols: list[pl.Expr] = []
    for spec in plan.contrasts:
        means = {}
        for role in ("target", "control"):
            values = [tissue_value(p, "raw", spec.tissue_aggregation) for p in tissue_columns(spec, rows, role).values()]
            means[role] = pl.mean_horizontal(values)
            cols += [means[role].alias(f"{spec.name}_{role}_mean_raw"),
                     pl.max_horizontal(values).alias(f"{spec.name}_{role}_max_raw"),
                     pl.min_horizontal(values).alias(f"{spec.name}_{role}_min_raw")]
        cols.append((means["target"] - means["control"]).alias(f"{spec.name}_target_vs_control_diff"))
    out = _add_top_tissue(master.with_columns(cols), rows, plan.contrast(plan.priority.direction_contrast))
    return add_priority_metrics(out, rows, plan)


def top_variants(contrast_df: pl.DataFrame, plan: ContrastPlan, n: int) -> pl.DataFrame:
    """The first `n` variants by the plan's priority metric.

    Order: consistent predicted direction first, then the largest |target - control| difference of
    the direction contrast, then variant ID. Columns: rank, `variant` (the key the Atlas links
    `table` command reads), avi_phred, avi_raw, top_modality, position, distance_to_tss,
    genomic_domain, ensembl_transcript_region, predicted_direction, specificity_scope and every
    contrast's difference.
    """
    diff = f"{plan.priority.direction_contrast}_target_vs_control_diff"
    diffs = [f"{c.name}_target_vs_control_diff" for c in plan.contrasts]
    ranked = (contrast_df
              .with_columns((pl.col("predicted_direction") != DIRECTION_INCONSISTENT).alias("_consistent"),
                            pl.col(diff).abs().alias("_priority"))
              .sort(["_consistent", "_priority", "variant_id"], descending=[True, True, False], nulls_last=True)
              .head(n))
    return ranked.select(
        pl.int_range(1, pl.len() + 1).alias("rank"), pl.col("variant_id").alias("variant"),
        "avi_phred", "avi_raw", pl.col("top_feature").alias("top_modality"), "position", "distance_to_tss",
        "genomic_domain", "ensembl_transcript_region", "predicted_direction", "specificity_scope", *diffs)


# --- stage ----------------------------------------------------------------------------------


def _need_inputs(study_dir: Path) -> tuple[pl.DataFrame, pl.DataFrame, dict[str, Any]]:
    from core import cohorts as core_cohorts

    qc = core_master.load_qc_report(study_dir)
    if not qc["overall_qc_passed"]:
        raise ContrastError("The Tier 1 QC did not pass (processed/qc_report.json); resolve it with expression.py master "
                            f"(Workflow step {core_master.WORKFLOW_STEP}) before contrasts.")
    return (core_master.load_master(study_dir), core_master.load_registry(study_dir),
            core_cohorts.load_resolution(study_dir))


def propose(study_dir: Path) -> tuple[Path, ContrastPlan, list[PendingDecision]]:
    """Writes processed/contrasts_proposed.yaml; returns it with the plan and the review points it would raise."""
    study_dir = Path(study_dir)
    master, registry, resolution = _need_inputs(study_dir)
    plan, skipped = propose_plan(registry, master, resolution=resolution)
    labels = {e["curie"]: (e["ols"]["label"] if e.get("ols") else e["curie"]) for e in resolution["curies"]}
    header = (f"# Proposed Contrast Plan for {plan.study_id} (generated). Review with the researcher, then save as "
              f"{PLAN_FILE} in the study directory.\n" + "".join(f"# Not included: {s}: {why}\n" for s, why in skipped.items()))
    path = study_dir / PROPOSED_FILE
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(plan.to_yaml(header, labels), encoding="utf-8")
    return path, plan, review_points(plan, registry.to_dicts())


_SHORT = {DIRECTION_INCONSISTENT: "Inconsistent", SCOPE_SHARED: "Shared", SCOPE_OPPOSITE: "Opposite",
          SCOPE_MIXED: "Mixed", SCOPE_NOT_CALLED: "No scope", SCOPE_NO_CONTROL: "No control"}


def _short(label: str) -> str:
    """Metric-key form of a priority label."""
    return _SHORT.get(label, label)


def run(study_dir: Path) -> pl.DataFrame:
    """Executes the approved plan on the Tier 1 table; writes the contrast table, its meta and the top variants.

    Exits with 1 (study.require) when `contrast.plan` has no answer.

    Raises:
        ContrastPending, ContrastError, MasterError, models.DependencyOutputError, study.StudyError.
    """
    study_dir = Path(study_dir)
    approval = require(study_dir, "contrast.plan", WORKFLOW_STEP)
    if str(approval) not in PLAN_APPROVED:
        raise ContrastError(f"Decision 'contrast.plan' is {approval!r}; expected 'approved'. Review the plan "
                            f"(Workflow step {WORKFLOW_STEP}) and record the approval with study.py decide.")
    study = load_study(study_dir)
    plan_path = study_dir / PLAN_FILE
    if not plan_path.is_file():
        raise ContrastError(f"{PLAN_FILE} not found in {study_dir}; run expression.py contrast --propose, review the "
                            f"proposal with the researcher and save it as {PLAN_FILE} (Workflow step {WORKFLOW_STEP}).")
    master, registry, resolution = _need_inputs(study_dir)
    plan = ContrastPlan.from_yaml(plan_path)
    if plan.study_id != study.study_id:
        raise ContrastError(f"{PLAN_FILE} is for study {plan.study_id!r}, not {study.study_id!r}.")
    rows = registry.to_dicts()
    plan.validate_against(rows)
    decisions = load_decisions(study_dir)
    pending, answered = [], []
    for q in review_points(plan, rows):
        entry = decisions.get(q.key)
        if entry is None:
            pending.append(q)
        elif entry.value == "abort":
            raise ContrastError(f"Aborted by the researcher at '{q.key}'; edit {PLAN_FILE} and run again.")
        elif entry.value == "accept":
            answered.append({"key": q.key, "value": entry.value})
        else:
            raise ContrastError(f"Decision '{q.key}' is {entry.value!r}; expected accept or abort.")
    if pending:
        raise ContrastPending(pending)

    tier2 = apply_plan(master, plan, registry)
    top = top_variants(tier2, plan, int(study.data["report"]["top_n"]))
    (study_dir / "processed").mkdir(parents=True, exist_ok=True)
    tier2.write_parquet(study_dir / CONTRAST_FILE)
    top.write_csv(study_dir / TOP_TSV_FILE, separator="\t")
    (study_dir / TOP_JSON_FILE).write_text(json.dumps(top.to_dicts(), indent=2) + "\n", encoding="utf-8")
    meta = {
        "rows": tier2.height, "plan_path": PLAN_FILE, "plan": plan.to_dict(),
        "review_decisions": answered + [{"key": "contrast.plan", "value": str(approval)}],
        "target_label": resolution["target_label"], "control_label": resolution["control_label"],
        "predicted_direction": dict(tier2["predicted_direction"].value_counts().iter_rows()),
        "specificity_scope": dict(tier2["specificity_scope"].value_counts().iter_rows()),
        "top_variants": top.height,
    }
    (study_dir / META_FILE).write_text(json.dumps(meta, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return tier2


def stage_metrics(study_dir: Path) -> dict[str, int | float | str]:
    """Metrics of the contrast stage compared by `reproduce` (expected.json stages.contrast)."""
    study_dir = Path(study_dir)
    path = study_dir / CONTRAST_FILE
    if not path.is_file() or not (study_dir / META_FILE).is_file():
        raise ContrastError(f"{CONTRAST_FILE} not found in {study_dir}; run expression.py contrast (Workflow step {WORKFLOW_STEP}).")
    df = pl.read_parquet(path)
    meta = json.loads((study_dir / META_FILE).read_text(encoding="utf-8"))
    metrics: dict[str, int | float | str] = {}
    for label, n in df["predicted_direction"].value_counts().iter_rows():
        metrics[f"direction.{_short(label)}"] = n
    for label, n in df["specificity_scope"].value_counts().iter_rows():
        metrics[f"scope.{_short(label)}"] = n
    for spec in meta["plan"]["contrasts"]:
        metrics[f"diff_sum.{spec['name']}"] = float(df[f"{spec['name']}_target_vs_control_diff"].sum())
    return metrics
