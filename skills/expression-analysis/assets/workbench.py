# /// script
# requires-python = ">=3.11"
# dependencies = ["marimo==0.25.2", "altair==6.3.0", "anywidget==0.11.0", "pandas==3.1.1", "polars==1.31.0", "pyarrow==26.0.1", "pyyaml==6.0.3", "matplotlib==3.11.2"]
# ///
"""Variant triage workbench (Workflow step 14), a marimo app over one study directory.

Open it with `expression.py dashboard --study-dir DIR [--port N] [--export-html FILE]`, or directly with
`marimo run assets/workbench.py -- --study-dir DIR`. It reads `processed/tissue_contrast.parquet`
and its companions from the study directory. Research use only.
"""

import marimo

__generated_with = "0.25.2"
app = marimo.App(width="full")


@app.cell
def _():
    import os
    import sys
    from pathlib import Path

    import marimo as mo
    import polars as pl

    # expression.py dashboard passes the skill's scripts/ folder; a direct `marimo run` finds it next to assets/.
    scripts_dir = Path(os.environ.get("EXPRESSION_ANALYSIS_SCRIPTS") or Path(__file__).resolve().parent.parent / "scripts")
    if str(scripts_dir) not in sys.path:
        sys.path.insert(0, str(scripts_dir))

    from viz.runtime import preload

    preload()  # one locked import of pandas / IPython / anywidget for all browser sessions

    from core.locus import footprint_blocks
    from viz.genomic_architecture import figure_svg
    from viz.priority_metrics import COLUMN_DOCS, column_docs
    from viz.variant_canvas import VariantCanvas, chart_spec
    from viz.workbench_helpers import (
        DEFAULT_MIN_QUANTILE,
        apply_filters,
        build_comparative_view,
        build_unified_lollipop_canvas,
        build_variant_tissue_matrix,
        center_crop_on_positions,
        contrast_of_metric,
        contrast_options,
        format_single_variant_tissue_table,
        format_slim_master_table,
        format_tissue_legend_md,
        load_workbench_datasets,
        metric_axis_title,
        metric_options,
    )

    _args = mo.cli_args()
    _study_dir = _args.get("study-dir") or _args.get("study_dir") or os.environ.get("STUDY_DIR")
    if not _study_dir:
        raise SystemExit("Pass the study directory: marimo run assets/workbench.py -- --study-dir DIR")
    data = load_workbench_datasets(Path(_study_dir))
    df = data.df
    docs = column_docs(data.target_label, data.control_label)
    return (
        COLUMN_DOCS,
        DEFAULT_MIN_QUANTILE,
        VariantCanvas,
        apply_filters,
        build_comparative_view,
        build_unified_lollipop_canvas,
        build_variant_tissue_matrix,
        center_crop_on_positions,
        chart_spec,
        contrast_of_metric,
        contrast_options,
        data,
        df,
        docs,
        figure_svg,
        footprint_blocks,
        format_single_variant_tissue_table,
        format_slim_master_table,
        format_tissue_legend_md,
        metric_axis_title,
        metric_options,
        mo,
        pl,
    )


@app.cell
def _(DEFAULT_MIN_QUANTILE, center_crop_on_positions, data, df, metric_options, mo, pl):
    # Shared state: picked variants (table rows and Figure 2 clicks) and the Figure 2 window
    get_focus_ids, set_focus_ids = mo.state([])
    # Self-loops let the buttons and the slider re-render their own cell after a change
    get_crop, set_crop = mo.state(tuple(data.full_domain), allow_self_loops=True)
    get_crop_mode, set_crop_mode = mo.state(False, allow_self_loops=True)

    def select_variants(ids):
        """Single entry point for table rows and canvas clicks: store the picks and re-centre the window on them."""
        ids = list(dict.fromkeys(ids))
        set_focus_ids(ids)
        if ids:
            positions = df.filter(pl.col("variant_id").is_in(ids))["position"].to_list()
            set_crop(center_crop_on_positions(positions, get_crop(), data.full_domain, margin_bp=10))

    metrics = metric_options(data)
    metric_dropdown = mo.ui.dropdown(options=metrics, value=next(iter(metrics)), label="Y metric")
    color_dropdown = mo.ui.dropdown(
        options=["genomic_domain", "predicted_direction", "specificity_scope", "ensembl_transcript_region", "top_feature_key"],
        value="genomic_domain", label="Color by")
    domains = df.sort("position")["genomic_domain"].unique(maintain_order=True).to_list()
    domain_filter = mo.ui.multiselect(options=domains, value=domains, label="Genomic domain")
    quantile_slider = mo.ui.slider(steps=sorted({0.0, 0.5, 0.9, 0.95, 0.99, 0.995, 0.999, DEFAULT_MIN_QUANTILE}), value=DEFAULT_MIN_QUANTILE,
                                   label=f"Min |quantile| in every {data.target_label} tissue", show_value=True)
    directions = df["predicted_direction"].unique().sort().to_list()
    direction_filter = mo.ui.multiselect(options=directions, value=directions, label="Predicted direction")
    scopes = df["specificity_scope"].unique().sort().to_list()
    scope_filter = mo.ui.multiselect(options=scopes, value=scopes, label="Specificity scope")
    n = df.height
    top_k_slider = mo.ui.slider(start=min(20, n), stop=n, step=20 if n > 20 else 1, value=n, label="Top-K by |Y|", show_value=True)
    table_sort_col = mo.ui.dropdown(options=["avi_phred", "target_min_abs_quantile", "position"], value="avi_phred", label="Sort table by")
    table_sort_order = mo.ui.dropdown(options=["descending", "ascending"], value="descending", label="Order")
    variant_multiselect = mo.ui.multiselect(options=df["variant_id"].to_list(), label="Compare variants", max_selections=50)
    return (
        color_dropdown,
        direction_filter,
        domain_filter,
        get_crop,
        get_crop_mode,
        get_focus_ids,
        metric_dropdown,
        quantile_slider,
        scope_filter,
        select_variants,
        set_crop,
        set_crop_mode,
        table_sort_col,
        table_sort_order,
        top_k_slider,
        variant_multiselect,
    )


@app.cell
def _(df, metric_dropdown, mo, pl):
    _max = df.select(pl.col(metric_dropdown.value).abs().max()).item() or 0.0
    y_threshold_slider = mo.ui.slider(start=0.0, stop=float(_max), step=float(_max) / 100 if _max else 0.01, value=0.0,
                                      label="Min |Y|", show_value=True)
    return (y_threshold_slider,)


@app.cell
def _(
    color_dropdown,
    data,
    direction_filter,
    domain_filter,
    format_tissue_legend_md,
    metric_dropdown,
    mo,
    quantile_slider,
    scope_filter,
    top_k_slider,
    y_threshold_slider,
):
    _qc = "passed" if data.qc_report["overall_qc_passed"] else "FAILED"
    mo.vstack([
        mo.md(f"# *{data.locus.gene_symbol}* variant triage workbench\n"
              f"`{data.locus.chromosome}:{data.locus.window.start:,}-{data.locus.window.end:,}` · "
              f"{data.df.height:,} SNVs · QC {_qc} · Research use only; not for clinical decisions."),
        mo.hstack([metric_dropdown, color_dropdown, y_threshold_slider], justify="start"),
        mo.hstack([quantile_slider, top_k_slider], justify="start"),
        mo.hstack([domain_filter, direction_filter, scope_filter], justify="start"),
        mo.accordion({"Cohorts and contrasts": mo.md(format_tissue_legend_md(data))}),
    ])
    return


@app.cell
def _(data, figure_svg, footprint_blocks, mo):
    _loc = data.locus
    mo.vstack([
        mo.md(f"## Figure 1 — Genomic architecture of the regulatory window\n"
              f"`{_loc.window_str}` ({_loc.window.length} bp, {3 * _loc.window.length:,} SNVs) around the TSS of "
              f"`{_loc.transcript_id}`: mutagenesis window, Ensembl Regulatory Build features and exon 1."),
        mo.Html(figure_svg(_loc, footprint_blocks(_loc))),
    ])
    return


@app.cell
def _(apply_filters, df, direction_filter, domain_filter, metric_dropdown, quantile_slider, scope_filter, y_threshold_slider):
    canvas_df = apply_filters(
        df, domains=domain_filter.value, directions=direction_filter.value, scopes=scope_filter.value,
        min_quantile=quantile_slider.value, metric_col=metric_dropdown.value, min_abs_y=y_threshold_slider.value)
    return (canvas_df,)


@app.cell
def _(canvas_df, metric_dropdown, pl, top_k_slider):
    active_df = canvas_df.sort(pl.col(metric_dropdown.value).abs(), descending=True, nulls_last=True).head(top_k_slider.value)
    return (active_df,)


@app.cell
def _(data, get_crop, get_crop_mode, mo, set_crop, set_crop_mode):
    _lo, _hi = data.full_domain
    crop_slider = mo.ui.range_slider(start=_lo, stop=_hi, step=1, value=list(get_crop()), label="Window (GRCh38)",
                                     show_value=False, full_width=True, on_change=lambda v: set_crop(tuple(v)))
    crop_button = mo.ui.button(label="Crop: drag over the plot" if get_crop_mode() else "Crop",
                               kind="danger" if get_crop_mode() else "neutral",
                               on_click=lambda _: set_crop_mode(not get_crop_mode()),
                               tooltip="Toggle crop mode: a drag selects the new window instead of panning")
    reset_button = mo.ui.button(label="Reset window", on_click=lambda _: set_crop(tuple(data.full_domain)),
                                disabled=tuple(get_crop()) == tuple(data.full_domain))
    return crop_button, crop_slider, reset_button


@app.cell
def _(
    VariantCanvas,
    active_df,
    build_unified_lollipop_canvas,
    chart_spec,
    color_dropdown,
    crop_button,
    crop_slider,
    data,
    get_crop,
    get_crop_mode,
    get_focus_ids,
    metric_axis_title,
    metric_dropdown,
    mo,
    quantile_slider,
    reset_button,
):
    _crop = tuple(get_crop())
    _mode = "crop" if get_crop_mode() else "pan"
    _axis = metric_axis_title(metric_dropdown.value, data)
    _chart = build_unified_lollipop_canvas(
        active_df, metric_dropdown.value, " ".join(_axis), color_dropdown.value,
        selected_variant_ids=get_focus_ids(), crop_domain=_crop, full_domain=data.full_domain,
        target_label=data.target_label, drag_mode=_mode, axis_title=_axis,
        title=f"{data.locus.gene_symbol}: {_axis[0]} across the regulatory window",
        subtitle=[f"{data.locus.chromosome}:{_crop[0]:,}-{_crop[1]:,} ({data.study.assembly}) · {active_df.height:,} SNVs shown · "
                  f"|quantile| ≥ {quantile_slider.value} in every {data.target_label} tissue",
                  "AlphaGenome predictions · research use only"])
    canvas = mo.ui.anywidget(VariantCanvas(spec=chart_spec(_chart), selected=list(get_focus_ids()), crop=list(_crop),
                                           file_name=f"{data.study.study_id}_figure2_{metric_dropdown.value}"))
    _tss = data.locus.tss
    _rel = sorted(((_crop[0] - _tss) * data.locus.strand, (_crop[1] - _tss) * data.locus.strand))
    _hint = "drag to select the new window" if _mode == "crop" else "drag to pan"
    mo.vstack([
        mo.md(f"## Figure 2 — {metric_dropdown.selected_key} along the window\n"
              f"{active_df.height:,} variants shown · window `{data.locus.chromosome}:{_crop[0]:,}-{_crop[1]:,}` "
              f"({_crop[1] - _crop[0] + 1} bp; {_rel[0]:+d} to {_rel[1]:+d} bp from the TSS) · {_hint}, hover for details, "
              "**click a variant to pick or unpick it** (dashed red ring; the window re-centres on picks). "
              "Save the figure from its **⋯** menu: SVG for papers, PNG (4× resolution) for posters."),
        mo.hstack([crop_button, reset_button], justify="start"),
        crop_slider,
        canvas,
    ])
    return (canvas,)


@app.cell
def _(canvas, data, get_crop, get_focus_ids, select_variants, set_crop, set_crop_mode):
    # Figure 2 -> shared state: clicked variants and a committed crop brush
    _picked = list(canvas.value.get("selected", []))
    if _picked != list(get_focus_ids()):
        select_variants(_picked)
    _new = tuple(canvas.value.get("crop", []) or ())
    if len(_new) == 2 and _new != tuple(get_crop()):
        set_crop((max(_new[0], data.full_domain[0]), min(_new[1], data.full_domain[1])))
        set_crop_mode(False)
    return


@app.cell
def _(
    active_df,
    docs,
    format_slim_master_table,
    data,
    get_focus_ids,
    mo,
    select_variants,
    table_sort_col,
    table_sort_order,
):
    slim_table_df = format_slim_master_table(active_df, data.plan, table_sort_col.value, table_sort_order.value)
    _ids = slim_table_df["variant_id"].to_list()
    kept_rows = [i for i, v in enumerate(_ids) if v in set(get_focus_ids())]
    variant_table = mo.ui.table(
        slim_table_df.to_dicts(),
        header_tooltip={k: v for k, v in docs.items() if k in slim_table_df.columns},
        selection="multi",
        initial_selection=kept_rows,
        page_size=12,
        on_change=lambda rows: select_variants([r["variant_id"] for r in rows or []]),
        label=f"{slim_table_df.height:,} variants after filters",
    )
    mo.vstack([mo.hstack([table_sort_col, table_sort_order], justify="start"), variant_table])
    return (variant_table,)


@app.cell
def _(active_df, contrast_of_metric, contrast_options, data, metric_dropdown, mo):
    _n = max(1, active_df.height)
    _contrasts = contrast_options(data)
    _default = contrast_of_metric(metric_dropdown.value, data) or data.plan.priority.direction_contrast
    matrix_contrast = mo.ui.dropdown(options=_contrasts, value=next(k for k, v in _contrasts.items() if v == _default),
                                     label="Track (contrast)")
    matrix_value = mo.ui.dropdown(options={"Quantile": "quantile", "Raw score": "raw"}, value="Quantile", label="Value")
    matrix_order = mo.ui.dropdown(options={"Genomic position": "position", "|Y metric|": "metric"}, value="Genomic position",
                                  label="Rows ordered by")
    matrix_rows = mo.ui.slider(start=1, stop=min(_n, 200), step=1, value=min(_n, 60), label="Max rows", show_value=True)
    return matrix_contrast, matrix_order, matrix_rows, matrix_value


@app.cell
def _(
    active_df,
    build_variant_tissue_matrix,
    data,
    get_crop,
    get_focus_ids,
    matrix_contrast,
    matrix_order,
    matrix_rows,
    matrix_value,
    metric_dropdown,
    mo,
    pl,
):
    _lo, _hi = get_crop()
    _in_window = active_df.filter(pl.col("position").is_between(_lo, _hi))
    if _in_window.is_empty():
        _out = mo.md("## Figure 3 — Variant × tissue matrix\nNo variant in the current filters and window.")
    else:
        _chart, _meta = build_variant_tissue_matrix(
            _in_window, data, matrix_contrast.value, value=matrix_value.value, max_variants=matrix_rows.value,
            order_by=matrix_order.value,
            metric_col=metric_dropdown.value, selected_variant_ids=get_focus_ids())
        _out = mo.vstack([
            mo.md(f"## Figure 3 — Variant × tissue matrix: {matrix_contrast.selected_key}\n"
                  f"{_meta['n_variants']} of {_meta['n_available']} variants in the current filters and window (rows) × "
                  f"the tissues of the selected contrast (columns), {matrix_value.selected_key.lower()} values. "
                  f"Tissue labels: blue = {data.target_label}, grey = {data.control_label}. Picked variants have red labels; pick them in Figure 2 or the table above."),
            mo.hstack([matrix_contrast, matrix_value, matrix_order, matrix_rows], justify="start"),
            _chart,
        ])
    _out
    return


@app.cell
def _(active_df, get_focus_ids, variant_multiselect):
    chosen_vids = list(get_focus_ids()) or list(variant_multiselect.value) or active_df["variant_id"].head(1).to_list()
    return (chosen_vids,)


@app.cell
def _(chosen_vids, mo, variant_multiselect):
    detail_variant_dropdown = mo.ui.dropdown(options=chosen_vids, value=chosen_vids[0] if chosen_vids else None,
                                             label="Per-tissue table for")
    mo.hstack([variant_multiselect, detail_variant_dropdown], justify="start")
    return (detail_variant_dropdown,)


@app.cell
def _(
    build_comparative_view,
    chosen_vids,
    data,
    detail_variant_dropdown,
    df,
    format_single_variant_tissue_table,
    format_tissue_legend_md,
    mo,
    pl,
):
    if not chosen_vids:
        _out = mo.md("No variant matches the filters.")
    else:
        bio_chart, _meta = build_comparative_view(df, chosen_vids, data)
        _vid = detail_variant_dropdown.value or chosen_vids[0]
        _row = df.filter(pl.col("variant_id") == _vid).row(0, named=True)
        _mode = "single variant" if _meta["n_variants"] == 1 else f"mean ± SD of {_meta['n_variants']} variants (SD in tooltip)"
        tissue_table = mo.ui.table(format_single_variant_tissue_table(df, _vid, data).to_dicts(), selection=None,
                                   label=f"Per-tissue values for {_vid}")
        _out = mo.vstack([
            mo.md(f"## Figure 4 — Variant impact by cell/tissue ({_mode})\n"
                  f"**{_vid}** · {_row['genomic_domain']} · {_row['predicted_direction']} · {_row['specificity_scope']} · "
                  f"[AlphaGenome Atlas]({_row['atlas_url']}) · [Ensembl]({_row['ensembl_location_url']})"),
            bio_chart,
            tissue_table,
            mo.md(format_tissue_legend_md(data)),
        ])
    _out
    return


if __name__ == "__main__":
    app.run()
