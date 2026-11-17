"""Live analysis of a study for the dashboard and `expression.py summarize`: tissues, effects, noise, filters.

Everything here is a pure function of the study's captured data (processed/master.parquet,
processed/track_registry.tsv, and raw/scores_RNA_SEQ_ACTIVE_baseline.tsv when the exploration captured a
baseline) and of the researcher's choices (Profile). Nothing is fixed in code: every threshold is a
Profile field whose value the researcher sets in the dashboard (or on the summarize command line); the
defaults are neutral (filters off) and only open the first view. Nothing is written to the study.

Scale of each scorer (alphagenome interpretation.md): RNA_SEQ is a natural-log fold change of the gene's
expression (0.001 pseudocount), CAGE/DNASE/ATAC a log2 ratio of the signal summed around the variant.
`pct` turns either into a percent change of what the scorer measures (SCALES names it). Percent
thresholds are converted on the side they apply to: a 30% decrease is log(0.7), a 30% increase log(1.3).

Score (Profile.score): "raw" is the scorer's own scale above; "quantile" is the Atlas quantile of the same
prediction (its rank among common variants, -1 to 1), comparable across scorers and tracks but not a
percent change, so the percent filters are skipped (the steps say so) and values are shown as they are.
Group (Profile.group_agg): the target and control values are the mean (default) or the median of their
tissues' effects.

Quantities (in the chosen score):
  tissue effect     mean score of the chosen tracks of that tissue (antisense tracks never enter)
  target_toward     smallest effect among target tissues in the direction of the question
                    (negative when some target tissue moves the other way)
  control_min/max   smallest and largest control tissue effect
  control_same      largest control effect in the same direction (0 when none)
  control_opposite  largest control effect in the opposite direction (0 when none)
  target_mean       mean (or median) of the target tissues' effects; control_mean likewise
  diff              target_mean minus control_mean
  diff_sigma        diff divided by the noise sigma at that magnitude (see noise_model)

Noise: for each ontology term with two or more chosen tracks of the scorer, the absolute difference
between the tracks (e.g. polyA vs total RNA-seq, GTEx vs ENCODE) over the window's SNVs is what two
measurements of the "same" tissue disagree by. It grows with the effect, so it is summarised per
magnitude bin (p95, p99 and a robust sigma); `noise_overall` gives the percentiles over all pairs. It is
pooled over every captured tissue with two tracks, which need not be the chosen ones. No pair -> no
noise model, and the sigma filter is skipped (with a note in the steps) instead of removing everything.

Library: Profile, SCORES, GROUP_AGGS, SCALES, pct, from_pct, uses_pct, unit, effect_text, load_study_data, tracks, tissue_table, reference_expression,
ref_by_tissue, noise_model, noise_overall, compute, expression_flags, apply, funnel.
"""

from __future__ import annotations

import dataclasses
import math
from pathlib import Path

import polars as pl

from core import master as core_master

BASELINE_NAMES = ("raw/scores_RNA_SEQ_ACTIVE_baseline.tsv", "atlas/scores_RNA_SEQ_ACTIVE_baseline.tsv")
PSEUDOCOUNT = 0.001  # RNA_SEQ LFC = ln((ALT + 0.001) / (REF + 0.001)) (alphagenome interpretation.md)
QUESTIONS = ("repression", "increase", "either")
SCORES = ("raw", "quantile")
GROUP_AGGS = ("mean", "median")
# scorer -> (log base, what a percent change refers to, with its article)
SCALES = {
    "RNA_SEQ": (math.e, "a expressão prevista do gene"),
    "CAGE": (2.0, "o sinal de início de transcrição em volta da variante"),
    "DNASE": (2.0, "a cromatina aberta (DNase) em volta da variante"),
    "ATAC": (2.0, "a cromatina aberta (ATAC) em volta da variante"),
}
UNIT = {"RNA_SEQ": "log-fold-change natural", "CAGE": "log2 do sinal local", "DNASE": "log2 do sinal local",
        "ATAC": "log2 do sinal local"}


def _base(scorer: str) -> float:
    return SCALES.get(scorer, (math.e, ""))[0]


def pct(value: float, scorer: str = "RNA_SEQ") -> float:
    """Percent change of what the scorer measures, from its raw (log) score."""
    return 100.0 * (_base(scorer) ** value - 1.0)


def from_pct(percent: float, scorer: str = "RNA_SEQ") -> float:
    """Raw (log) score of a signed percent change (+30 -> log 1.3, -30 -> log 0.7); 0 stays 0."""
    if not percent:
        return 0.0
    if percent <= -100:
        raise ValueError("a decrease of 100% or more has no log value")
    return math.log(1.0 + percent / 100.0, _base(scorer))


# --------------------------------------------------------------------------- choices


@dataclasses.dataclass
class Profile:
    """The researcher's analysis choices (dashboard controls or summarize options). Defaults are neutral."""

    scorer: str = "RNA_SEQ"
    target: list[str] = dataclasses.field(default_factory=list)      # ontology CURIEs
    control: list[str] = dataclasses.field(default_factory=list)
    shown: list[str] = dataclasses.field(default_factory=list)       # extra tissues shown as context
    target_label: str = "Alvo"         # group names shown in legends, figures and texts
    control_label: str = "Controle"
    color_by: str = "genomic_domain"   # colour of the window chart: genomic_domain | passa | sentido | ensembl_transcript_region | top_feature
    assays: list[str] = dataclasses.field(default_factory=list)      # assay titles kept; empty = all
    sources: list[str] = dataclasses.field(default_factory=list)     # data sources kept; empty = all
    question: str = "either"           # repression | increase | either
    min_change_pct: float = 0.0        # smallest % change in every target tissue, in the question's direction
    max_control_pct: float | None = None  # largest % change allowed in any control tissue, either way; None = off
    min_diff_sigma: float = 0.0        # diff / noise sigma in the question's direction (abs for "either"); 0 = off
    ref_floor: float | None = None     # reference expression floor (RNA_SEQ REF); None = off
    ref_floor_applies: str = "warn"    # warn | control | target | both
    domains: list[str] = dataclasses.field(default_factory=list)     # empty = all
    one_alt_per_position: bool = False
    y_metric: str = "target_mean"      # target_mean | target_toward | control_mean | diff | diff_sigma | avi_phred
    top_k: int | None = None           # display limit ('Quantas mostrar'), not a filter: see first_k
    score: str = "raw"                 # raw | quantile (Atlas quantile of the same prediction)
    group_agg: str = "mean"            # mean | median of the tissues of each group


def uses_pct(p: Profile) -> bool:
    """Percent changes (and the percent filters) exist only on the raw scale."""
    return p.score == "raw"


def unit(p: Profile) -> str:
    """Unit of the effect values shown for the profile."""
    return UNIT.get(p.scorer, "log") if uses_pct(p) else "quantil do Atlas, -1 a 1"


def effect_text(value: float | None, p: Profile, digits: int = 1) -> str:
    """An effect as the reader sees it: a signed percent on the raw scale, the quantile itself otherwise."""
    if value is None or value != value:
        return "—"
    if uses_pct(p):
        return f"{pct(value, p.scorer):+.{digits}f}%".replace(".", ",")
    return f"{value:+.3f}".replace(".", ",")


# --------------------------------------------------------------------------- data


@dataclasses.dataclass
class StudyData:
    study_dir: Path
    master: pl.DataFrame
    registry: pl.DataFrame
    baseline: pl.DataFrame | None   # RNA_SEQ_ACTIVE of one SNV of the gene, every track, when captured
    baseline_variant: str | None


def _find_baseline(study_dir: Path) -> Path | None:
    """The study's own copy first (init --from-exploration puts it in raw/), then the exploration next to it."""
    for base in (study_dir, study_dir.parent / "exploration", study_dir.parent):
        for name in BASELINE_NAMES:
            if (base / name).is_file():
                return base / name
    return None


def load_study_data(study_dir: Path) -> StudyData:
    study_dir = Path(study_dir)
    master = core_master.load_master(study_dir)
    registry = core_master.load_registry(study_dir)
    baseline, variant = None, None
    path = _find_baseline(study_dir)
    if path is not None:
        b = pl.read_csv(path, separator="\t", infer_schema_length=10000).filter(pl.col("scorer") == "RNA_SEQ_ACTIVE")
        if b.height:
            variant = b["variant"][0]
            baseline = b.filter(pl.col("variant") == variant)
    return StudyData(study_dir, master, registry, baseline, variant)


def tracks(data: StudyData, scorer: str, assays: list[str] | None = None,
           sources: list[str] | None = None, score: str = "raw") -> pl.DataFrame:
    """Tracks of a scorer in the registry, optionally restricted to assay titles and data sources.

    Stranded tracks on the opposite strand of the gene never enter (CAGE keeps only the gene's strand,
    as the contrast stage does), so a tissue's mean and its noise pairs compare like with like.
    """
    t = data.registry.filter(pl.col("scorer") == scorer)
    if "strand_class" in t.columns:
        t = t.filter(pl.col("strand_class") == "sense") if scorer == "CAGE" else t.filter(pl.col("strand_class") != "antisense")
    if assays:
        t = t.filter(pl.col("assay_title").is_in(assays))
    if sources:
        t = t.filter(pl.col("data_source").is_in(sources))
    return t.filter(pl.col("column_prefix").map_elements(lambda p: f"{p}_{score}" in data.master.columns,
                                                          return_dtype=pl.Boolean))


def reference_expression(data: StudyData) -> pl.DataFrame | None:
    """Predicted reference expression of the gene per RNA_SEQ track: REF = ACTIVE when LFC <= 0, else
    (ACTIVE + 0.001) * exp(-LFC) - 0.001, with ACTIVE from the baseline SNV and LFC from the same SNV in
    master. Exact by definition of the two scorers. None when no baseline was captured or when its SNV
    is outside the study's window (it is chosen in exon 1)."""
    if data.baseline is None or data.baseline_variant is None:
        return None
    row = data.master.filter(pl.col("variant_id") == data.baseline_variant)
    if row.height == 0:
        return None
    rna = data.registry.filter(pl.col("scorer") == "RNA_SEQ").select("track_name", "column_prefix", "ontology_curie")
    out = []
    act = {r["track_name"]: r["raw_score"] for r in data.baseline.iter_rows(named=True)}
    for r in rna.iter_rows(named=True):
        col = f"{r['column_prefix']}_raw"
        if r["track_name"] not in act or col not in row.columns:
            continue
        a, lfc = act[r["track_name"]], row[col][0]
        if a is None or lfc is None:
            continue
        ref = a if lfc <= 0 else (a + PSEUDOCOUNT) * math.exp(-lfc) - PSEUDOCOUNT
        out.append({"ontology_curie": r["ontology_curie"], "track_name": r["track_name"], "ref": max(ref, 0.0)})
    return pl.DataFrame(out) if out else None


def tissue_table(data: StudyData, scorer: str) -> pl.DataFrame:
    """One row per ontology term with tracks of the scorer: name, cohort in the study, track count,
    assays, sources, replicate available, and reference expression (min and max over its RNA_SEQ tracks)."""
    t = tracks(data, scorer)
    g = (t.group_by("ontology_curie")
         .agg(pl.col("biosample_name").first().alias("tissue"),
              pl.col("cohort").first().alias("study_cohort"),
              pl.len().alias("tracks"),
              pl.col("assay_title").unique().sort().str.join(", ").alias("assays"),
              pl.col("data_source").unique().sort().str.join(", ").alias("sources"))
         .with_columns((pl.col("tracks") >= 2).alias("replicate")))
    ref = reference_expression(data)
    if ref is not None:
        r = ref.group_by("ontology_curie").agg(pl.col("ref").min().alias("ref_min"), pl.col("ref").max().alias("ref_max"))
        g = g.join(r, on="ontology_curie", how="left")
    else:
        g = g.with_columns(pl.lit(None, dtype=pl.Float64).alias("ref_min"), pl.lit(None, dtype=pl.Float64).alias("ref_max"))
    return g.sort(["study_cohort", "tissue"])


# --------------------------------------------------------------------------- noise


def _pair_differences(data: StudyData, scorer: str, assays: list[str] | None, sources: list[str] | None,
                      score: str = "raw") -> pl.DataFrame | None:
    t = tracks(data, scorer, assays, sources, score)
    pairs = []
    for _, grp in t.group_by("ontology_curie"):
        cols = [f"{p}_{score}" for p in grp["column_prefix"].to_list()]
        for i in range(len(cols)):
            for j in range(i + 1, len(cols)):
                a, b = data.master[cols[i]], data.master[cols[j]]
                pairs.append(pl.DataFrame({"m": ((a.abs() + b.abs()) / 2), "d": (a - b).abs()}))
    return pl.concat(pairs).drop_nulls() if pairs else None


def noise_model(data: StudyData, scorer: str, assays: list[str] | None = None, sources: list[str] | None = None,
                bins: int = 5, score: str = "raw") -> pl.DataFrame | None:
    """Disagreement between tracks of the same ontology term, per magnitude bin.

    Columns: lo, hi (bin edges on the mean |effect| of the pair), n, p95, p99 of the absolute pair
    difference, and sigma = 1.4826 * its median (a robust scale for the difference between two
    measurements of the same tissue, the same kind of quantity as diff). None when no term has two tracks.
    """
    df = _pair_differences(data, scorer, assays, sources, score)
    if df is None or df.height == 0:
        return None
    edges = [df["m"].quantile(q) for q in [i / bins for i in range(bins + 1)]]
    edges[0], edges[-1] = 0.0, float("inf")
    rows = []
    for lo, hi in zip(edges, edges[1:]):
        s = df.filter((pl.col("m") >= lo) & (pl.col("m") < hi))["d"]
        if s.len() == 0:
            continue
        rows.append({"lo": lo, "hi": hi, "n": s.len(), "p95": s.quantile(0.95), "p99": s.quantile(0.99),
                     "sigma": 1.4826 * s.median()})
    return pl.DataFrame(rows) if rows else None


def noise_overall(data: StudyData, scorer: str, assays: list[str] | None = None,
                  sources: list[str] | None = None, score: str = "raw") -> tuple[float, float] | None:
    """(p95, p99) of the disagreement over all pairs of tracks of the same tissue (in the score chosen); None without pairs."""
    df = _pair_differences(data, scorer, assays, sources, score)
    if df is None or df.height == 0:
        return None
    return float(df["d"].quantile(0.95)), float(df["d"].quantile(0.99))


def _sigma_at(noise: pl.DataFrame | None, magnitude: pl.Expr) -> pl.Expr:
    if noise is None or noise.height == 0:
        return pl.lit(None, dtype=pl.Float64)
    expr = pl.lit(noise["sigma"][-1])
    for r in reversed(list(noise.iter_rows(named=True))):
        expr = pl.when(magnitude < r["hi"]).then(pl.lit(r["sigma"])).otherwise(expr)
    return expr


# --------------------------------------------------------------------------- compute and filter


def _question_sign(p: Profile, target_mean: pl.Expr) -> pl.Expr:
    if p.question == "repression":
        return pl.lit(-1.0)
    if p.question == "increase":
        return pl.lit(1.0)
    return pl.when(target_mean >= 0).then(1.0).otherwise(-1.0)


def compute(data: StudyData, p: Profile, noise: pl.DataFrame | None = None) -> pl.DataFrame:
    """Per-SNV quantities for the profile's tissues (see module docstring)."""
    t = tracks(data, p.scorer, p.assays or None, p.sources or None, p.score)
    by_tissue: dict[str, list[str]] = {}
    for r in t.iter_rows(named=True):
        by_tissue.setdefault(r["ontology_curie"], []).append(f"{r['column_prefix']}_{p.score}")
    keep = ["variant_id", "position", "ref_allele", "alt_allele", "genomic_domain", "distance_to_tss",
            "ensembl_transcript_region", "avi_phred", "top_feature", "atlas_url"]
    df = data.master.select([c for c in keep if c in data.master.columns])
    tissues = list(dict.fromkeys([*p.target, *p.control, *p.shown]))
    eff = {c: pl.mean_horizontal([data.master[col] for col in by_tissue[c]]) for c in tissues if c in by_tissue}
    df = df.with_columns([pl.Series(f"eff_{c}", data.master.select(eff[c]).to_series()) for c in eff])
    target = [c for c in p.target if c in eff]
    control = [c for c in p.control if c in eff]
    none_cols = ("control_mean", "control_min", "control_max", "control_same", "control_opposite")
    if not target:
        return df.with_columns(pl.lit(None, dtype=pl.Float64).alias(x) for x in
                               ("target_mean", "target_toward", *none_cols, "diff", "diff_sigma"))
    tm = _group([pl.col(f"eff_{c}") for c in target], p.group_agg)
    df = df.with_columns(tm.alias("target_mean"))
    df = df.with_columns(_question_sign(p, pl.col("target_mean")).alias("_s"))
    df = df.with_columns(pl.min_horizontal([pl.col(f"eff_{c}") * pl.col("_s") for c in target]).alias("target_toward"))
    if control:
        cols = [pl.col(f"eff_{c}") for c in control]
        df = df.with_columns(
            _group(cols, p.group_agg).alias("control_mean"),
            pl.min_horizontal(cols).alias("control_min"), pl.max_horizontal(cols).alias("control_max"),
            pl.max_horizontal([c * pl.col("_s") for c in cols]).clip(lower_bound=0).alias("control_same"),
            pl.max_horizontal([-c * pl.col("_s") for c in cols]).clip(lower_bound=0).alias("control_opposite"))
    else:
        df = df.with_columns(pl.lit(None, dtype=pl.Float64).alias(x) for x in none_cols)
    df = df.with_columns((pl.col("target_mean") - pl.col("control_mean")).alias("diff"))
    mag = (pl.col("target_mean").abs() + pl.col("control_mean").abs()) / 2
    df = df.with_columns((pl.col("diff") / _sigma_at(noise, mag)).alias("diff_sigma"))
    return df.drop("_s")


def _group(cols: list[pl.Expr], how: str) -> pl.Expr:
    """Mean or median of the group's tissue effects, per variant."""
    return pl.concat_list(cols).list.median() if how == "median" else pl.mean_horizontal(cols)


def _ref_by_tissue(data: StudyData) -> dict[str, float]:
    ref = reference_expression(data)
    if ref is None:
        return {}
    return {r["ontology_curie"]: r["ref"] for r in ref.group_by("ontology_curie").agg(pl.col("ref").mean()).iter_rows(named=True)}


def expression_flags(data: StudyData, p: Profile) -> dict[str, str]:
    """Tissue -> 'below floor', 'ok', 'unknown' (no REF) or 'floor off', for the profile's tissues."""
    ref = _ref_by_tissue(data)
    out = {}
    for c in dict.fromkeys([*p.target, *p.control]):
        if p.ref_floor is None:
            out[c] = "floor off"
        elif c not in ref:
            out[c] = "unknown"
        else:
            out[c] = "below floor" if ref[c] < p.ref_floor else "ok"
    return out


def apply(df: pl.DataFrame, data: StudyData, p: Profile, noise: pl.DataFrame | None = None
          ) -> tuple[pl.DataFrame, list[tuple[str, int]]]:
    """Applies the profile's filters in order; returns the kept rows and the count after each step.

    A filter that cannot be evaluated (no control tissue, no noise model) is skipped and the step says
    why, instead of removing every variant.
    """
    steps: list[tuple[str, int]] = [("todos os SNVs", df.height)]

    def step(label: str, cond: pl.Expr | None) -> None:
        nonlocal df
        if cond is not None:
            df = df.filter(cond)
        steps.append((label, df.height))

    has_control = "control_mean" in df.columns and df["control_mean"].null_count() < df.height
    if p.domains:
        step("domínios", pl.col("genomic_domain").is_in(p.domains))
    flags = expression_flags(data, p)
    if p.ref_floor is not None and p.ref_floor_applies in ("target", "both"):
        bad = [c for c in p.target if flags.get(c) == "below floor"]
        step("piso de expressão no alvo", pl.lit(not bad))
    pct_ok = uses_pct(p)
    if not pct_ok and (p.min_change_pct > 0 or p.max_control_pct is not None):
        steps.append(("mudança mínima no alvo / máxima no controle: ignorado (pontuação em quantil, sem %)", df.height))
    if not pct_ok and p.question != "either":
        # in quantiles the question still keeps only variants moving every target tissue the way asked
        step("todo tecido alvo no sentido da pergunta", pl.col("target_toward") > 0)
    elif p.min_change_pct > 0 or p.question != "either":
        # every target tissue must move in the question's direction by at least the % chosen, converted on
        # its own side (a 30% decrease is log 0.7, a 30% increase log 1.3); 0 and "either" skip the step
        up = from_pct(p.min_change_pct, p.scorer)
        down = -from_pct(-min(p.min_change_pct, 99.999), p.scorer)
        thr = pl.when(pl.col("target_mean") >= 0).then(up).otherwise(down) if p.question == "either" else \
            pl.lit(up if p.question == "increase" else down)
        step("mudança no alvo ≥ mínimo (todo tecido alvo, no sentido da pergunta)", pl.col("target_toward") >= thr)
    if p.max_control_pct is not None and pct_ok:
        if has_control:
            hi = from_pct(p.max_control_pct, p.scorer)
            lo = from_pct(-min(p.max_control_pct, 99.999), p.scorer)
            step("mudança no controle ≤ máximo (todo tecido controle, nos dois sentidos)",
                 (pl.col("control_min") >= lo) & (pl.col("control_max") <= hi))
        else:
            steps.append(("mudança no controle ≤ máximo: ignorado (sem tecido controle)", df.height))
    if p.min_diff_sigma > 0:
        if not has_control:
            steps.append(("alvo − controle acima do ruído: ignorado (sem tecido controle)", df.height))
        elif noise is None:
            steps.append(("alvo − controle acima do ruído: ignorado (sem medida de ruído)", df.height))
        else:
            # in the question's direction: an increase must be target above control, a repression target below
            # control; "either" accepts both. A variant that raises the target cannot pass a repression question
            s = p.min_diff_sigma
            if p.question == "increase":
                step(f"alvo − controle ≥ {s:g} σ do ruído (alvo acima)", pl.col("diff_sigma") >= s)
            elif p.question == "repression":
                step(f"alvo − controle ≥ {s:g} σ do ruído (alvo abaixo)", pl.col("diff_sigma") <= -s)
            else:
                step(f"alvo − controle ≥ {s:g} σ do ruído (nos dois sentidos)", pl.col("diff_sigma").abs() >= s)
    metric = p.y_metric if p.y_metric in df.columns else "target_mean"
    if p.one_alt_per_position:
        df = df.sort(pl.col(metric).abs(), descending=True, nulls_last=True).unique("position", keep="first", maintain_order=True)
        steps.append(("uma alternativa por posição", df.height))
    return df.sort(order_expr(metric), descending=True, nulls_last=True), steps


def order_expr(metric: str) -> pl.Expr:
    """Sort key of 'Ordenar por': the size of the change, except target_toward, which is already signed toward the
    question (larger = more in the direction asked)."""
    return pl.col(metric) if metric == "target_toward" else pl.col(metric).abs()


def first_k(df: pl.DataFrame, metric: str, k: int | None) -> pl.DataFrame:
    """The first `k` rows in the order of 'Ordenar por' (all when k is 0 or None). A display limit, not a filter."""
    metric = metric if metric in df.columns else "target_mean"
    df = df.sort(order_expr(metric), descending=True, nulls_last=True)
    return df.head(k) if k else df


def funnel(steps: list[tuple[str, int]]) -> pl.DataFrame:
    return pl.DataFrame({"etapa": [s for s, _ in steps], "SNVs": [n for _, n in steps]})


def ref_by_tissue(data: StudyData) -> dict[str, float]:
    """Ontology term -> mean predicted reference expression over its RNA_SEQ tracks (empty when unavailable)."""
    return _ref_by_tissue(data)
