"""Workbench data loading, filters, Figure 2 (lollipop canvas), Figure 3 (variant x tissue matrix),
Figure 4 (impact by tissue) and table helpers.

Reads a study directory (processed/tissue_contrast.parquet and its companions). Everything that
names a contrast, tissue or cohort comes from the approved contrast plan, the track registry and
study.yaml; nothing here is specific to one gene. `default_filter` is the variant set the workbench
shows when opened, and is what the `workbench` stage of `reproduce` counts.
"""

from __future__ import annotations

import dataclasses
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional, Sequence

import altair as alt
import polars as pl

from core import cohorts as core_cohorts
from core import contrast as core_contrast
from core import locus as core_locus
from core import master as core_master
from core.contrast import ContrastPlan, tissue_columns, tissue_value
from core.locus import LocusConfig
from study import StudyConfig, load_study
from viz.design_system import apply_print_style, assay_label, domain_color_scale
from viz.priority_metrics import DIRECTION_COLORS

TARGET_COLOR = "#2563eb"
CONTROL_COLOR = "#64748b"
SLIM_BASE_COLUMNS = ["variant_id", "predicted_direction", "target_min_abs_quantile", "specificity_scope", "position",
                     "distance_to_tss", "genomic_domain", "ensembl_transcript_region", "avi_phred", "top_feature_key"]
SLIM_LINK_COLUMNS = ["atlas_url", "ensembl_location_url"]


DEFAULT_MIN_QUANTILE = 0.99
"""Default of the workbench slider: |quantile| >= 0.99 in every target tissue of the direction contrast."""


@dataclass(frozen=True)
class WorkbenchData:
    """Everything the workbench reads from a study directory.

    Attributes:
        df: processed/tissue_contrast.parquet (Tier 2, with the priority labels).
        registry: processed/track_registry.tsv rows.
        plan: The contrast plan that the contrast step executed (from tissue_contrast_meta.json).
        locus: processed/locus.json.
        study: study.yaml.
        resolution: processed/ontology_resolution.json.
        qc_report: processed/qc_report.json.
    """

    df: pl.DataFrame
    registry: list[dict[str, Any]]
    plan: ContrastPlan
    locus: LocusConfig
    study: StudyConfig
    resolution: dict[str, Any]
    qc_report: dict[str, Any]

    @property
    def target_label(self) -> str:
        return self.resolution["target_label"]

    @property
    def control_label(self) -> str:
        return self.resolution["control_label"]

    @property
    def full_domain(self) -> tuple[int, int]:
        return (self.locus.window.start, self.locus.window.end)


def load_workbench_datasets(study_dir: Path | str) -> WorkbenchData:
    """Reads Tier 2, registry, the executed plan, locus, resolution and QC report from a study directory.

    Raises:
        core.contrast.ContrastError, core.master.MasterError, core.locus.LocusError, study.StudyError:
            a required processed file is missing (the message names the step that writes it).
    """
    study_dir = Path(study_dir)
    path = study_dir / core_contrast.CONTRAST_FILE
    meta_path = study_dir / core_contrast.META_FILE
    if not path.is_file() or not meta_path.is_file():
        raise core_contrast.ContrastError(
            f"{core_contrast.CONTRAST_FILE} not found in {study_dir}; run expression.py contrast "
            f"(Workflow step {core_contrast.WORKFLOW_STEP}) before the workbench.")
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    return WorkbenchData(
        df=pl.read_parquet(path),
        registry=core_master.load_registry(study_dir).to_dicts(),
        plan=ContrastPlan.from_dict(meta["plan"]),
        locus=core_locus.load_locus(study_dir),
        study=load_study(study_dir),
        resolution=core_cohorts.load_resolution(study_dir),
        qc_report=core_master.load_qc_report(study_dir),
    )


# --------------------------------------------------------------------------- filters


def apply_filters(df: pl.DataFrame, *, domains: Optional[Sequence[str]] = None, directions: Optional[Sequence[str]] = None,
                  scopes: Optional[Sequence[str]] = None, min_quantile: float = DEFAULT_MIN_QUANTILE,
                  metric_col: Optional[str] = None, min_abs_y: float = 0.0) -> pl.DataFrame:
    """The workbench filters: domain, direction and scope lists (None = all), |quantile| and |Y| thresholds."""
    cond = pl.col("target_min_abs_quantile") >= min_quantile
    if domains is not None:
        cond &= pl.col("genomic_domain").is_in(list(domains))
    if directions is not None:
        cond &= pl.col("predicted_direction").is_in(list(directions))
    if scopes is not None:
        cond &= pl.col("specificity_scope").is_in(list(scopes))
    if metric_col is not None:
        cond &= pl.col(metric_col).abs() >= min_abs_y
    return df.filter(cond)


def default_filter(master_contrast: pl.DataFrame) -> pl.DataFrame:
    """Variants shown by the workbench with its default filters (every domain, direction and scope;
    |quantile| >= 0.99 in every target tissue of the direction contrast)."""
    return apply_filters(master_contrast)


def domain_options(df: pl.DataFrame) -> list[str]:
    """Genomic-domain options of the workbench filter, in order of position."""
    return df.sort("position")["genomic_domain"].unique(maintain_order=True).to_list()


# --------------------------------------------------------------------------- metrics & table


def metric_options(data: WorkbenchData) -> dict[str, str]:
    """Display label -> column, one group per plan contrast (direction contrast first), then AVI Phred."""
    t, c = data.target_label, data.control_label
    names = [data.plan.priority.direction_contrast] + [x.name for x in data.plan.contrasts
                                                       if x.name != data.plan.priority.direction_contrast]
    options: dict[str, str] = {}
    for name in names:
        spec = data.plan.contrast(name)
        tag = f"{name} ({spec.scorer}{', ' + spec.strand if spec.strand not in ('any', 'unstranded') else ''})"
        options[f"{tag}: {t} − {c}"] = f"{name}_target_vs_control_diff"
        options[f"{tag}: {t} mean"] = f"{name}_target_mean_raw"
        options[f"{tag}: {c} mean"] = f"{name}_control_mean_raw"
    options["AVI Phred"] = "avi_phred"
    return options


def metric_axis_title(metric_col: str, data: WorkbenchData) -> list[str]:
    """Readable, self-contained axis title (two lines) for a Y metric column, for print and export."""
    if metric_col == "avi_phred":
        return ["AlphaGenome variant impact", "(AVI, Phred-scaled)"]
    for spec in data.plan.contrasts:
        assay = assay_label(spec.scorer)
        if spec.strand not in ("any", "unstranded"):
            assay = f"{assay} ({spec.strand} strand)"
        if metric_col == f"{spec.name}_target_vs_control_diff":
            return [f"Δ predicted {assay} effect", f"{data.target_label} − {data.control_label} (raw score)"]
        if metric_col == f"{spec.name}_target_mean_raw":
            return [f"Predicted {assay} effect", f"{data.target_label} mean (raw score)"]
        if metric_col == f"{spec.name}_control_mean_raw":
            return [f"Predicted {assay} effect", f"{data.control_label} mean (raw score)"]
    return [metric_col.replace("_", " ")]


def contrast_columns(plan: ContrastPlan) -> list[str]:
    return [col for spec in plan.contrasts for col in
            (f"{spec.name}_target_mean_raw", f"{spec.name}_control_mean_raw", f"{spec.name}_target_vs_control_diff")]


def format_slim_master_table(df: pl.DataFrame, plan: ContrastPlan, sort_by_col: str = "avi_phred",
                             sort_order: str = "descending") -> pl.DataFrame:
    cols = [c for c in [*SLIM_BASE_COLUMNS, *contrast_columns(plan), *SLIM_LINK_COLUMNS] if c in df.columns]
    out = df.select(cols).with_columns(pl.col(pl.Float64).round(4))
    return out.sort(sort_by_col, descending=sort_order == "descending", nulls_last=True) if sort_by_col in out.columns else out


# --------------------------------------------------------------------------- crop / selection state


def center_crop_on_positions(positions: Sequence[int], current_crop: Optional[tuple[int, int]],
                             full_domain: tuple[int, int], margin_bp: int = 10) -> tuple[int, int]:
    """Window that shows all `positions`: the current crop if they are already inside it; otherwise the crop
    shifted (same width) to contain them, widened only when they do not fit in that width."""
    if not positions:
        return current_crop or full_domain
    lo, hi = min(positions) - margin_bp, max(positions) + margin_bp
    if current_crop:
        c_lo, c_hi = current_crop
        if c_lo <= min(positions) and max(positions) <= c_hi:
            return (c_lo, c_hi)
        width = c_hi - c_lo
        if hi - lo <= width:
            center = (lo + hi) // 2
            lo, hi = center - width // 2, center - width // 2 + width
            if lo < full_domain[0]:
                lo, hi = full_domain[0], full_domain[0] + width
            if hi > full_domain[1]:
                lo, hi = full_domain[1] - width, full_domain[1]
    return (max(lo, full_domain[0]), min(hi, full_domain[1]))


# --------------------------------------------------------------------------- Figure 2

CROP_BRUSH_NAME = "crop_brush"
"""Name of the Figure 2 interval selection that `VariantCanvas` reports as the new window."""


def build_unified_lollipop_canvas(
    df: pl.DataFrame, metric_col: str, metric_label: str, color_col: str, *, max_lollipops: int = 800,
    selected_variant_ids: Optional[Sequence[str]] = None, crop_domain: Optional[tuple[int, int]] = None,
    full_domain: tuple[int, int], max_focus_labels: int = 5, target_label: str = "target", drag_mode: str = "pan",
    axis_title: Optional[Sequence[str]] = None, title: Optional[str] = None,
    subtitle: Optional[str | Sequence[str]] = None, color_title: Optional[str] = None,
    color_domain: Optional[Sequence[str]] = None, color_range: Optional[Sequence[str]] = None,
    faded: Optional[tuple[str, str]] = None,
    extra_tooltip: Optional[Sequence[tuple[str, str, Optional[str]]]] = None,
    x_title: str = "GRCh38 position (bp)", lang: str = "en", head_size: float = 90,
):
    """Figure 2: one lollipop per variant at its position, height = `metric_col`.

    The x domain is `crop_domain`. `drag_mode="pan"` lets a drag slide the view; `drag_mode="crop"`
    draws the `crop_brush` interval that `VariantCanvas` reports as the new window. Picked variants
    (`selected_variant_ids`) get a dashed red ring and a dashed vertical guide, and, when few, their
    ID at the top of the plot. Clicks are
    handled by `viz.variant_canvas.VariantCanvas`, because marimo does not report selections from
    layered charts.
    """
    data = (df.filter(pl.col(metric_col).is_not_null())
            .sort(pl.col(metric_col).abs(), descending=True).head(max_lollipops)
            .with_columns(pl.lit(0.0).alias("_zero"), pl.col(metric_col).alias("_y_val")))
    extra = list(extra_tooltip or [])
    keep = ["variant_id", "position", "distance_to_tss", "genomic_domain", "ensembl_transcript_region", "avi_phred",
            "predicted_direction", "specificity_scope", "_zero", "_y_val", color_col, *[c for c, _, _ in extra]]
    data = data.select(list(dict.fromkeys(c for c in keep if c in data.columns)))
    domain = list(crop_domain or full_domain)
    y_title = list(axis_title) if axis_title else metric_label
    x = alt.X("position:Q", title=x_title, scale=alt.Scale(domain=domain),
              axis=alt.Axis(format=",d", tickCount=6, labelOverlap="greedy", labelFlush=True))
    if color_domain is not None and color_range is not None:
        color = alt.Color(f"{color_col}:N", title=color_title or color_col,
                          scale=alt.Scale(domain=list(color_domain), range=list(color_range)))
    elif color_col == "predicted_direction":
        color = alt.Color(f"{color_col}:N", title="Predicted direction",
                          scale=alt.Scale(domain=list(DIRECTION_COLORS), range=list(DIRECTION_COLORS.values())))
    elif color_col == "genomic_domain":
        domains, colors = domain_color_scale(data["genomic_domain"].unique(maintain_order=True).to_list())
        color = alt.Color("genomic_domain:N", title="Genomic domain", scale=alt.Scale(domain=domains, range=colors))
    else:
        color = alt.Color(f"{color_col}:N", title=color_title or color_col.replace("_", " ").capitalize())
    if faded is not None and faded[0] in data.columns:
        # rows where faded[0] == faded[1] (e.g. variants that do not pass the filters) are drawn light grey
        is_faded = alt.datum[faded[0]] == faded[1]
        color = alt.condition(is_faded, alt.value("#c8c4c0"), color)
        fade = lambda solid: alt.condition(is_faded, alt.value(0.35), alt.value(solid))
    else:
        fade = lambda solid: alt.value(solid)
    tooltip = [alt.Tooltip("variant_id:N", title="Variant ID"), alt.Tooltip("position:Q", title="Position", format="d"),
               alt.Tooltip("distance_to_tss:Q", title="Distance to TSS (bp)"),
               alt.Tooltip("genomic_domain:N", title="Genomic Domain"),
               alt.Tooltip("ensembl_transcript_region:N", title="Consequence"),
               alt.Tooltip("_y_val:Q", title=metric_label, format=".4f"),
               alt.Tooltip("avi_phred:Q", title="AVI Phred", format=".2f"),
               alt.Tooltip("predicted_direction:N", title=f"Predicted {target_label} effect"),
               alt.Tooltip("specificity_scope:N", title="Specificity scope")]
    tooltip = [t for t in tooltip if t.to_dict()["field"] in data.columns]
    if lang == "pt":
        pt = {"Variant ID": "Variante", "Position": "Posição", "Distance to TSS (bp)": "Distância ao TSS (pb)",
              "Genomic Domain": "Domínio genômico", "Consequence": "Região do transcrito", "AVI Phred": "AVI Phred"}
        for tt in tooltip:
            tip_title = tt.to_dict().get("title")
            if tip_title in pt:
                tt.title = pt[tip_title]
    tooltip += [alt.Tooltip(f"{c}:{'Q' if fmt else 'N'}", title=t, **({"format": fmt} if fmt else {})) for c, t, fmt in extra]
    base = alt.Chart(data.to_pandas())
    spine = base.mark_tick(orient="horizontal", thickness=3, size=10, opacity=0.4, clip=True).encode(
        x=x, y=alt.Y("_zero:Q", title=y_title))
    stems = base.mark_rule(strokeWidth=1.8, clip=True).encode(x=x, y="_zero:Q", y2="_y_val:Q", color=color,
                                                             opacity=fade(0.7))
    if drag_mode == "crop":
        drag = alt.selection_interval(name=CROP_BRUSH_NAME, encodings=["x"], translate=False, zoom=False,
                                      mark=alt.BrushConfig(fill="#e11d48", fillOpacity=0.12, stroke="#e11d48", strokeWidth=1.5))
    else:
        drag = alt.selection_interval(bind="scales", encodings=["x", "y"], zoom=True)  # wheel zooms, drag pans, double-click resets
    heads = base.mark_circle(size=head_size, stroke="#ffffff", strokeWidth=0.8, clip=True, cursor="pointer").encode(
        x=x, y=alt.Y("_y_val:Q", title=y_title), color=color, opacity=fade(0.9), tooltip=tooltip)
    layers = [spine.add_params(drag), stems, heads]
    focus = data.filter(pl.col("variant_id").is_in(list(selected_variant_ids or [])))
    if not focus.is_empty():
        focus_base = alt.Chart(focus.to_pandas())
        layers.append(focus_base.mark_point(shape="circle", size=420, filled=False, strokeWidth=2.6, strokeDash=[4, 3],
                                            opacity=0.95, clip=True, stroke="#e11d48", cursor="pointer").encode(
            x=x, y="_y_val:Q", tooltip=[alt.Tooltip("variant_id:N", title="Picked variant (click to unpick)")]))
        layers.append(focus_base.mark_rule(strokeDash=[5, 4], strokeWidth=1.6, color="#e11d48", opacity=0.75, clip=True)
                      .encode(x=x))
        # One label per picked position (variants at the same position share it, e.g. "chr1:123 A>G, A>T"),
        # each label on its own row so labels at nearby positions do not overlap.
        groups = (focus.sort("position").with_columns(pl.col("variant_id").str.split(":").alias("_p"))
                  .group_by("position", maintain_order=True)
                  .agg(pl.col("_p").list.get(0).first().alias("_chr"),
                       pl.col("_p").list.get(-1).str.join(", ").alias("_alts")))
        groups = groups.with_columns((pl.col("_chr") + ":" + pl.col("position").cast(pl.Utf8) + " " + pl.col("_alts")).alias("_label"))
        if groups.height <= max_focus_labels:
            mid = (domain[0] + domain[1]) / 2
            for row, g in enumerate(groups.iter_rows(named=True)):
                side, dx = ("left", 6) if g["position"] < mid else ("right", -6)
                layers.append(alt.Chart(pl.DataFrame([g]).select("position", "_label").to_pandas()).mark_text(
                    y=4 + 18 * row, baseline="top", dx=dx, align=side, fontSize=13, fontWeight="bold", color="#e11d48",
                    clip=True).encode(x=x, text="_label:N"))
    chart = alt.layer(*layers).properties(width="container", height=440)
    if title:
        chart = chart.properties(title=alt.TitleParams(text=title, subtitle=(list(subtitle) if isinstance(subtitle, (list, tuple)) else subtitle) or alt.Undefined))
    return apply_print_style(chart)


# --------------------------------------------------------------------------- Figure 3


def tissue_effect_table(df: pl.DataFrame, variant_ids: Sequence[str], data: WorkbenchData,
                        contrasts: Optional[Sequence[core_contrast.ContrastSpec]] = None) -> pl.DataFrame:
    """Long table: one row per (variant, contrast, tissue) with the tissue raw and quantile values.

    `contrasts` restricts the table to some plan contrasts (default: all)."""
    sub = df.filter(pl.col("variant_id").is_in(list(variant_ids)))
    labels = {t["ontology_curie"]: t["biosample_name"] for t in data.registry}
    frames = []
    for spec in (contrasts if contrasts is not None else data.plan.contrasts):
        for role, cohort_label in (("target", data.target_label), ("control", data.control_label)):
            for curie, prefixes in tissue_columns(spec, data.registry, role).items():
                frames.append(sub.select(
                    "variant_id",
                    pl.lit(f"{labels.get(curie, curie)} ({curie})").alias("biosample"),
                    pl.lit(cohort_label).alias("cohort"),
                    pl.lit(f"{spec.name} ({spec.scorer})").alias("assay"),
                    tissue_value(prefixes, "raw", spec.tissue_aggregation).alias("raw_score"),
                    tissue_value(prefixes, "quantile", spec.tissue_aggregation).alias("quantile_score"),
                ))
    return pl.concat(frames) if frames else pl.DataFrame()


def contrast_options(data: WorkbenchData) -> dict[str, str]:
    """Display label -> contrast name, direction contrast first (same naming as the Y metric options)."""
    names = [data.plan.priority.direction_contrast] + [c.name for c in data.plan.contrasts
                                                       if c.name != data.plan.priority.direction_contrast]
    out = {}
    for name in names:
        spec = data.plan.contrast(name)
        strand = f", {spec.strand}" if spec.strand not in ("any", "unstranded") else ""
        out[f"{name} ({spec.scorer}{strand})"] = name
    return out


def contrast_of_metric(metric_col: str, data: WorkbenchData) -> Optional[str]:
    """Contrast whose columns the Y metric belongs to (None for non-contrast metrics such as AVI Phred)."""
    return next((c.name for c in data.plan.contrasts if metric_col.startswith(f"{c.name}_")), None)


def build_variant_tissue_matrix(
    df: pl.DataFrame, data: WorkbenchData, contrast: str, *, value: str = "quantile", max_variants: int = 60,
    order_by: str = "position", metric_col: Optional[str] = None,
    selected_variant_ids: Optional[Sequence[str]] = None,
) -> tuple[alt.Chart, dict[str, Any]]:
    """Figure 3: variants (rows) × the tissues of one plan contrast (columns).

    Columns are the target then control tissues of `contrast`, with the same tracks and tissue aggregation
    as Figure 4 (`tissue_effect_table`). `value` is `quantile` (fixed [-1, 1] scale) or `raw` (symmetric
    scale around 0). Rows are the first `max_variants` variants by genomic position, or by |`metric_col`|
    when `order_by == "metric"`. Target tissue labels are blue, control grey, picked variants' rows red.
    """
    if order_by == "metric" and metric_col:
        rows = df.sort(pl.col(metric_col).abs(), descending=True, nulls_last=True)
    else:
        rows = df.sort("position")
    rows = rows.head(max_variants)
    ids = rows["variant_id"].to_list()
    spec = data.plan.contrast(contrast)
    long = tissue_effect_table(rows, ids, data, [spec])
    col = "quantile_score" if value == "quantile" else "raw_score"
    tissue_order = list(dict.fromkeys(long["biosample"].to_list()))
    target_tissues = list(dict.fromkeys(long.filter(pl.col("cohort") == data.target_label)["biosample"].to_list()))
    color_scale = (alt.Scale(scheme="redblue", reverse=True, domain=[-1, 1]) if value == "quantile"
                   else alt.Scale(scheme="redblue", reverse=True, domainMid=0))
    picked = list(selected_variant_ids or [])
    label_color = f"indexof({picked!r}, datum.value) >= 0 ? '#e11d48' : '#334155'" if picked else "'#334155'"
    chart = alt.Chart(long.with_columns(pl.col(col).alias("value")).to_pandas()).mark_rect().encode(
        x=alt.X("biosample:N", sort=tissue_order, title=None,
                axis=alt.Axis(orient="top", labelAngle=-35, labelLimit=320, labelFontSize=13, labelColor=alt.ExprRef(
                    f"indexof({target_tissues!r}, datum.value) >= 0 ? '{TARGET_COLOR}' : '{CONTROL_COLOR}'"))),
        y=alt.Y("variant_id:N", sort=ids, title=None, axis=alt.Axis(labelColor=alt.ExprRef(label_color), labelFontSize=12)),
        color=alt.Color("value:Q", scale=color_scale,
                        title=f"{assay_label(spec.scorer)} {'quantile' if value == 'quantile' else 'raw score'}",
                        legend=alt.Legend(gradientLength=220, gradientThickness=14)),
        tooltip=["variant_id:N", "biosample:N", "cohort:N",
                 alt.Tooltip("raw_score:Q", format=".4f"), alt.Tooltip("quantile_score:Q", format=".4f")],
    ).properties(width=alt.Step(64), height=alt.Step(18))
    return apply_print_style(chart), {"n_variants": len(ids), "n_available": df.height, "contrast": contrast, "tissues": tissue_order,
                   "value": col}


def build_comparative_view(df: pl.DataFrame, selected_variant_ids: Sequence[str], data: WorkbenchData
                           ) -> tuple[alt.Chart, dict[str, Any]]:
    """Figure 4: per tissue and contrast, the raw effect of one variant or the mean ± SD of several."""
    long = tissue_effect_table(df, selected_variant_ids, data)
    n = len(set(selected_variant_ids))
    summary = long.group_by(["biosample", "cohort", "assay"], maintain_order=True).agg(
        pl.col("raw_score").mean().alias("raw_score"), pl.col("raw_score").std().alias("sd"),
        pl.col("raw_score").min().alias("min"), pl.col("raw_score").max().alias("max"))
    order = list(dict.fromkeys(summary["biosample"].to_list()))
    assays = list(dict.fromkeys(summary["assay"].to_list()))
    opacities = [round(0.95 - 0.35 * i / max(1, len(assays) - 1), 2) for i in range(len(assays))]
    tooltip = ["biosample:N", "cohort:N", "assay:N", alt.Tooltip("raw_score:Q", format=".4f")]
    if n > 1:
        tooltip += [alt.Tooltip("sd:Q", format=".4f"), alt.Tooltip("min:Q", format=".4f"), alt.Tooltip("max:Q", format=".4f")]
    chart = alt.Chart(summary.to_pandas()).mark_bar(cornerRadiusEnd=3).encode(
        y=alt.Y("biosample:N", sort=order, title=None),
        yOffset=alt.YOffset("assay:N", sort=assays),
        x=alt.X("raw_score:Q", title=("Predicted effect, mean of selected variants (raw score)" if n > 1
                                      else "Predicted effect (raw score)")),
        color=alt.Color("cohort:N", scale=alt.Scale(domain=[data.target_label, data.control_label],
                                                    range=[TARGET_COLOR, CONTROL_COLOR])),
        opacity=alt.Opacity("assay:N", scale=alt.Scale(domain=assays, range=opacities)),
        tooltip=tooltip,
    ).properties(width="container", height=max(320, 26 * len(order) * max(1, len(assays))))
    return apply_print_style(chart), {"n_variants": n, "tissues": order, "assays": assays}


def format_single_variant_tissue_table(df: pl.DataFrame, variant_id: str, data: WorkbenchData) -> pl.DataFrame:
    return (tissue_effect_table(df, [variant_id], data).drop("variant_id")
            .with_columns(pl.col("raw_score").round(5), pl.col("quantile_score").round(4)))


def format_tissue_legend_md(data: WorkbenchData) -> str:
    lines = [f"**Cohorts:** <span style='color:{TARGET_COLOR}'>■</span> {data.target_label} (target) · "
             f"<span style='color:{CONTROL_COLOR}'>■</span> {data.control_label} (control)", "",
             "**Contrasts (approved plan):**"]
    for spec in data.plan.contrasts:
        filt = {k: v for k, v in spec.track_filter.to_dict().items() if v}
        lines.append(f"* `{spec.name}`: {spec.scorer}, strand {spec.strand}, tissue {spec.tissue_aggregation}"
                     + (f", filter {filt}" if filt else ""))
    lines.append(f"* Priority labels from `{data.plan.priority.direction_contrast}`.")
    return "\n".join(lines)


# --------------------------------------------------------------------------- reproduce stage


def build_default_figures(data: WorkbenchData) -> dict[str, Any]:
    """Builds Figures 2-4 for the default state of the workbench (smoke test of the whole figure code).

    Returns:
        {"shown": rows after the default filters, "fig2": Vega-Lite dict, "fig3_meta": ..., "fig4_meta": ...}.

    Raises:
        ValueError: the default filters leave no variant.
    """
    shown = default_filter(data.df)
    if shown.is_empty():
        raise ValueError("The default filters leave no variant; lower the |quantile| threshold in the workbench.")
    metric_label, metric_col = next(iter(metric_options(data).items()))
    axis = metric_axis_title(metric_col, data)
    fig2 = build_unified_lollipop_canvas(shown, metric_col, metric_label, "genomic_domain", full_domain=data.full_domain,
                                         target_label=data.target_label, axis_title=axis)
    _, fig3_meta = build_variant_tissue_matrix(shown, data, data.plan.priority.direction_contrast)
    picked = shown["variant_id"].head(1).to_list()
    _, fig4_meta = build_comparative_view(data.df, picked, data)
    return {"shown": shown, "fig2": fig2.to_dict(), "fig3_meta": fig3_meta, "fig4_meta": fig4_meta}


def stage_metrics(study_dir: Path) -> dict[str, int | float | str]:
    """Metrics of the workbench stage compared by `reproduce` (expected.json stages.workbench)."""
    data = load_workbench_datasets(study_dir)
    return {
        "default_filter_count": default_filter(data.df).height,
        "domain_options": "|".join(domain_options(data.df)),
    }
