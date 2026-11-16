# /// script
# requires-python = ">=3.11"
# dependencies = ["marimo==0.25.2", "altair==6.3.0", "anywidget==0.11.0", "pandas==3.1.1", "polars==1.31.0", "pyarrow==26.0.1", "pyyaml==6.0.3", "matplotlib==3.11.2", "vl-convert-python==1.9.0.post1"]
# ///
"""Analysis dashboard of an expression study (marimo). Research use only.

See first, refine later: the page opens on the result for the study's own groups with no filter, a
sentence saying what it shows, and only three choices in view (question, target tissues, control
tissues). Filters and advanced options (group names, tracks, expression floor, presentation) are in
sections that open on demand. Nothing here calls the network or writes to the study; the agent reproduces
any view with `expression.py summarize` (same computation, same options) to explain it.

Open: `expression.py dashboard --study-dir DIR` (app) or `... --edit` (marimo editor; edits this file).
"""

import marimo

__generated_with = "0.25.2"
app = marimo.App(width="full")


@app.cell
def _():
    import dataclasses
    import os
    import sys
    from pathlib import Path

    import altair as alt
    import marimo as mo
    import polars as pl

    scripts_dir = Path(os.environ.get("EXPRESSION_ANALYSIS_SCRIPTS") or Path(__file__).resolve().parent.parent / "scripts")
    if str(scripts_dir) not in sys.path:
        sys.path.insert(0, str(scripts_dir))
    import study as study_mod
    from core import analysis as A
    from core import locus as core_locus
    from viz import dashboard_helpers as H
    from viz.genomic_architecture import figure_svg
    from viz.variant_canvas import VariantCanvas, chart_spec
    from viz.workbench_helpers import build_unified_lollipop_canvas, center_crop_on_positions

    alt.data_transformers.disable_max_rows()
    _args = mo.cli_args()
    study_dir = Path(_args.get("study-dir") or _args.get("study_dir") or os.environ.get("STUDY_DIR") or ".").resolve()
    data = A.load_study_data(study_dir)
    locus = core_locus.load_locus(study_dir)
    study_cfg = study_mod.load_study(study_dir)
    full_domain = (locus.window.start, locus.window.end)
    ref_tissue = A.ref_by_tissue(data)

    def save_buttons(chart, name, width=1100):
        """PNG (4x, for posters) and SVG (vector, for papers) of a chart as drawn now (filters, window and
        highlights; the mouse zoom of the browser view is not included)."""
        import vl_convert as vlc

        def _spec():
            try:
                c = chart.properties(width=width) if chart.to_dict().get("width") == "container" else chart
            except Exception:
                c = chart
            return c.to_dict()

        return mo.hstack([
            mo.download(data=lambda: vlc.vegalite_to_png(_spec(), scale=4), filename=f"{name}.png", mimetype="image/png",
                        label="Salvar PNG"),
            mo.download(data=lambda: vlc.vegalite_to_svg(_spec()).encode("utf-8"), filename=f"{name}.svg",
                        mimetype="image/svg+xml", label="Salvar SVG"),
        ], justify="start")
    return (A, H, VariantCanvas, alt, build_unified_lollipop_canvas, center_crop_on_positions, chart_spec, core_locus,
            data, dataclasses, figure_svg, full_domain, locus, mo, pl, ref_tissue, save_buttons, study_cfg, study_dir)


@app.cell
def _(A, data, mo):
    # Advanced: scorer (shown in the "Avançado" section; read here because it sets the tissues and tracks)
    scorers = sorted(set(data.registry["scorer"].to_list()))
    scorer_pick = mo.ui.dropdown(options=scorers, value="RNA_SEQ" if "RNA_SEQ" in scorers else scorers[0],
                                 label="Tipo de predição")
    return (scorer_pick,)


@app.cell
def _(A, data, study_cfg):
    # Starting point: the study's own groups, no filter
    _reg = data.registry
    _coh = study_cfg.data.get("cohorts", {})
    base = A.Profile(target=sorted(set(_reg.filter(_reg["cohort"] == "target")["ontology_curie"].to_list())),
                     control=sorted(set(_reg.filter(_reg["cohort"] == "control")["ontology_curie"].to_list())),
                     target_label=_coh.get("target", {}).get("label") or "Alvo",
                     control_label=_coh.get("control", {}).get("label") or "Controle")
    return (base,)


@app.cell
def _(A, H, data, scorer_pick):
    tissues = A.tissue_table(data, scorer_pick.value)
    descriptions = H.describe_tissues(data, scorer_pick.value)
    _t = A.tracks(data, scorer_pick.value)
    assays = sorted(set(_t["assay_title"].drop_nulls().to_list()))
    sources = sorted(set(_t["data_source"].drop_nulls().to_list()))
    labels = {r["ontology_curie"]: f"{r['tissue']} ({r['tracks']} {'medições' if r['tracks'] > 1 else 'medição'})"
              for r in tissues.iter_rows(named=True)}
    return assays, descriptions, labels, sources, tissues


@app.cell
def _(H, assays, base, data, labels, mo, sources):
    _opts = {v: k for k, v in labels.items()}
    _pick = lambda xs: [labels[x] for x in xs if x in labels]
    # --- the three main choices
    question = mo.ui.radio(options={"Diminui no alvo": "repression", "Aumenta no alvo": "increase",
                                    "Qualquer mudança": "either"},
                           value={"repression": "Diminui no alvo", "increase": "Aumenta no alvo"}.get(
                               base.question, "Qualquer mudança"), label="O que você procura?", inline=True)
    target_pick = mo.ui.multiselect(options=_opts, value=_pick(base.target), label="Onde (tecidos do alvo)")
    control_pick = mo.ui.multiselect(options=_opts, value=_pick(base.control), label="Comparado com (tecidos do controle)")
    # --- filters
    min_target_pct = mo.ui.number(start=0.0, stop=500.0, step=1.0, value=base.min_change_pct,
                                  label="Mudança mínima em cada tecido do alvo (%)")
    use_max_control = mo.ui.switch(value=base.max_control_pct is not None, label="Exigir pouca mudança no controle")
    max_control_pct = mo.ui.number(start=0.0, stop=500.0, step=1.0,
                                   value=base.max_control_pct if base.max_control_pct is not None else 2.0,
                                   label="Mudança máxima no controle (%)")
    min_sigma = mo.ui.number(start=0.0, stop=20.0, step=0.5, value=base.min_diff_sigma,
                             label="Alvo − controle acima do ruído (× σ; 0 = desligado)")
    _domains = data.master.sort("position")["genomic_domain"].unique(maintain_order=True).to_list()
    domain_pick = mo.ui.multiselect(options=_domains, value=[d for d in base.domains if d in _domains] or _domains,
                                    label="Regiões da janela (vazio = todas)")
    # --- advanced
    target_name = mo.ui.text(value=base.target_label, label="Nome do grupo alvo")
    control_name = mo.ui.text(value=base.control_label, label="Nome do grupo controle")
    shown_pick = mo.ui.multiselect(options=_opts, value=_pick(base.shown), label="Tecidos de contexto (só exibidos)")
    assay_pick = mo.ui.multiselect(options=assays, value=[a for a in base.assays if a in assays] or assays, label="Ensaios que entram (vazio = todos)")
    source_pick = mo.ui.multiselect(options=sources, value=[s for s in base.sources if s in sources] or sources, label="Fontes que entram (vazio = todas)")
    use_floor = mo.ui.switch(value=base.ref_floor is not None, label="Usar piso de expressão de referência")
    ref_floor = mo.ui.number(start=0.0, stop=100.0, step=0.005, value=base.ref_floor if base.ref_floor is not None else 0.01,
                             label="Piso (REF)")
    floor_applies = mo.ui.dropdown(options={"Só avisar": "warn", "Tirar controles abaixo do piso": "control",
                                            "Não analisar se o alvo estiver abaixo": "target", "Os dois": "both"},
                                   value={"control": "Tirar controles abaixo do piso", "target": "Não analisar se o alvo estiver abaixo",
                                          "both": "Os dois"}.get(base.ref_floor_applies, "Só avisar"), label="O que fazer com o piso")
    one_alt = mo.ui.switch(value=base.one_alt_per_position, label="Uma alternativa por posição")
    _ymetric = H.Y_METRICS
    y_metric = mo.ui.dropdown(options=_ymetric, value=next((k for k, v in _ymetric.items() if v == base.y_metric), "Mudança no alvo"),
                              label="Ordenar por")
    _colors = H.color_options()
    color_by = mo.ui.dropdown(options=_colors, value=next((k for k, v in _colors.items() if v == base.color_by), "Domínio genômico"),
                              label="Cor no gráfico da janela")
    show_mode = mo.ui.dropdown(options={"Só as que passam nos filtros": "kept", "Todas (as que não passam em cinza)": "all"},
                               value="Só as que passam nos filtros", label="Variantes nos gráficos")
    top_k = mo.ui.number(start=0, stop=data.master.height, step=1, value=base.top_k or 0,
                         label="Quantas mostrar (0 = todas)")
    return (assay_pick, color_by, control_name, control_pick, domain_pick, floor_applies, max_control_pct, min_sigma,
            min_target_pct, one_alt, question, ref_floor, show_mode, shown_pick, source_pick, target_name, target_pick, top_k,
            use_floor, use_max_control, y_metric)


@app.cell
def _(A, assay_pick, color_by, control_name, control_pick, domain_pick, floor_applies, max_control_pct,
      min_sigma, min_target_pct, one_alt, question, ref_floor, scorer_pick, shown_pick, source_pick, target_name,
      target_pick, top_k, use_floor, use_max_control, y_metric):
    profile = A.Profile(
        scorer=scorer_pick.value,
        target=list(target_pick.value), control=list(control_pick.value), shown=list(shown_pick.value),
        target_label=target_name.value.strip() or "Alvo", control_label=control_name.value.strip() or "Controle",
        color_by=color_by.value, assays=list(assay_pick.value), sources=list(source_pick.value), question=question.value,
        min_change_pct=float(min_target_pct.value or 0.0),
        max_control_pct=float(max_control_pct.value or 0.0) if use_max_control.value else None,
        min_diff_sigma=float(min_sigma.value), ref_floor=float(ref_floor.value) if use_floor.value else None,
        ref_floor_applies=floor_applies.value, domains=list(domain_pick.value), one_alt_per_position=one_alt.value,
        y_metric=y_metric.value, top_k=int(top_k.value) or None)
    return (profile,)


@app.cell
def _(A, H, data, dataclasses, pl, profile, show_mode):
    noise = A.noise_model(data, profile.scorer, profile.assays or None, profile.sources or None)
    _p = profile
    if _p.ref_floor is not None and _p.ref_floor_applies in ("control", "both"):
        _flags = A.expression_flags(data, _p)
        _p = dataclasses.replace(_p, control=[c for c in _p.control if _flags.get(c) != "below floor"])
    effects = A.compute(data, _p, noise)
    kept, steps = A.apply(effects, data, _p, noise)
    flags = A.expression_flags(data, profile)
    metric = profile.y_metric if profile.y_metric in effects.columns else "target_mean"
    # What the charts and the table show: the variants that pass (or all, the rest in grey), the first K in the
    # order chosen. K limits what is shown; the headline and the filter steps still count every variant that passes.
    shown_df = H.display_frame(effects, kept, profile, show_mode.value)
    return effects, flags, kept, metric, noise, shown_df, steps


# ============================================================================ page


@app.cell
def _(A, data, locus, mo, scorer_pick):
    _what = A.SCALES.get(scorer_pick.value, (None, "o sinal previsto"))[1]
    mo.vstack([
        mo.md(f"# *{locus.gene_symbol}*: quais variantes mudam {_what}\n"
              f"Janela `{locus.chromosome}:{locus.window.start:,}-{locus.window.end:,}` · {data.master.height:,} variantes testadas "
              "(todas as trocas de uma base) · predições do AlphaGenome, uso somente em pesquisa."),
        mo.accordion({"Como usar": mo.md(
            "A página já abre com o resultado para os grupos do estudo, sem filtros. Mude as três escolhas logo abaixo e veja o "
            "resultado se atualizar. Filtros e opções avançadas ficam nos blocos mais abaixo, que se abrem quando você quiser.\n\n"
            f"As mudanças são em **% d{_what}** (por exemplo −20% = a variante reduz em um quinto). "
            "O tipo de predição muda em *Avançado*: RNA-seq mede a expressão do gene inteiro; CAGE, DNase e ATAC medem o "
            "sinal local em volta da variante.")}),
    ])
    return


@app.cell
def _(control_pick, descriptions, mo, question, show_mode, target_pick, top_k, y_metric):
    _desc = lambda xs: "  \n".join(f"- {descriptions.get(c, '')}" for c in xs if descriptions.get(c))
    mo.vstack([
        question,
        mo.hstack([target_pick, control_pick], justify="start", gap=2),
        mo.hstack([y_metric, top_k, show_mode], justify="start", gap=2),
        mo.accordion({"De onde vêm as medições destes tecidos": mo.md(
            "**Alvo**  \n" + (_desc(target_pick.value) or "nenhum tecido escolhido")
            + "\n\n**Controle**  \n" + (_desc(control_pick.value) or "nenhum tecido escolhido")
            + "\n\nGTEx e ENCODE são os projetos que mediram a expressão usada para treinar o modelo: GTEx mede tecidos de adultos "
            "(mistura de tipos celulares); ENCODE inclui tipos celulares isolados, linhagens e alguns tecidos. polyA mede sobretudo "
            "mRNA maduro; total inclui RNA ainda não processado."),
                      "Como mostrar": mo.md(
            "**Ordenar por** define a ordem da tabela, quais entram em *Quantas mostrar* e a altura de cada haste no "
            "gráfico da janela. *Mudança no alvo* é a média dos tecidos do alvo; *no sentido procurado, no tecido que "
            "menos muda* favorece variantes que mudam todos os tecidos do alvo; *diferença alvo − controle* favorece "
            "efeitos próprios do alvo. **Quantas mostrar** desenha e lista só as N primeiras nessa ordem (0 = todas); o resumo "
            "e as etapas dos filtros continuam contando todas as que passam. **Variantes nos gráficos**: com *Todas*, as que "
            "não passam nos filtros aparecem em cinza, e *Quantas mostrar* vale sobre todas.")}),
    ])
    return


@app.cell
def _(H, effects, kept, labels, mo, profile, ref_tissue):
    if not profile.target:
        _out = mo.callout(mo.md("Escolha pelo menos um tecido do alvo."), kind="warn")
    else:
        _out = mo.callout(mo.md(H.headline(kept, effects, profile, labels, ref_tissue)), kind="info")
    _out
    return


@app.cell
def _(center_crop_on_positions, data, full_domain, mo, pl):
    get_focus_ids, set_focus_ids = mo.state([])
    get_crop, set_crop = mo.state(tuple(full_domain), allow_self_loops=True)
    get_crop_mode, set_crop_mode = mo.state(False, allow_self_loops=True)

    def select_variants(ids):
        """Única entrada para destaques (gráficos, tabela, identificador): guarda e centraliza a janela."""
        ids = list(dict.fromkeys(ids))
        set_focus_ids(ids)
        if ids:
            positions = data.master.filter(pl.col("variant_id").is_in(ids))["position"].to_list()
            set_crop(center_crop_on_positions(positions, get_crop(), full_domain, margin_bp=10))
    return get_crop, get_crop_mode, get_focus_ids, select_variants, set_crop, set_crop_mode, set_focus_ids


@app.cell
def _(full_domain, get_crop, get_crop_mode, get_focus_ids, mo, set_crop, set_crop_mode, set_focus_ids):
    crop_slider = mo.ui.range_slider(start=full_domain[0], stop=full_domain[1], step=1, value=list(get_crop()),
                                     label="Trecho da janela", show_value=False, full_width=True, on_change=lambda v: set_crop(tuple(v)))
    crop_button = mo.ui.button(label="Recortar: arraste sobre o gráfico" if get_crop_mode() else "Recortar",
                               kind="danger" if get_crop_mode() else "neutral",
                               on_click=lambda _: set_crop_mode(not get_crop_mode()),
                               tooltip="Arrastar no gráfico escolhe o trecho a mostrar")
    reset_button = mo.ui.button(label="Janela inteira", on_click=lambda _: set_crop(tuple(full_domain)),
                                disabled=tuple(get_crop()) == tuple(full_domain))
    clear_button = mo.ui.button(label="Limpar destaques", on_click=lambda _: set_focus_ids([]), disabled=not get_focus_ids())
    return clear_button, crop_button, crop_slider, reset_button


@app.cell
def _(mo):
    point_size = mo.ui.slider(start=5, stop=120, step=5, value=40, label="Tamanho dos pontos", show_value=True)
    return (point_size,)


@app.cell
def _(A, H, VariantCanvas, build_unified_lollipop_canvas, chart_spec, clear_button, color_by, crop_button, crop_slider,
      figure_svg, full_domain, get_crop, get_crop_mode, get_focus_ids, locus, metric, mo, one_alt, pl, point_size, profile,
      reset_button, save_buttons, shown_df, study_cfg, y_metric):
    from core.locus import footprint_blocks
    _crop = tuple(get_crop())
    _mode = "crop" if get_crop_mode() else "pan"
    _df = shown_df
    if one_alt.value:
        _df = _df.sort(pl.col(metric).abs(), descending=True, nulls_last=True).unique("position", keep="first")
    _dom, _rng = H.fixed_colors(color_by.value)
    _unit = "" if metric in ("avi_phred", "diff_sigma") else f" ({A.UNIT.get(profile.scorer, 'log')})"
    _chart = build_unified_lollipop_canvas(
        _df, metric, y_metric.selected_key, color_by.value, max_lollipops=max(1, _df.height),
        selected_variant_ids=get_focus_ids(), crop_domain=_crop, full_domain=full_domain, target_label=profile.target_label,
        drag_mode=_mode, axis_title=[y_metric.selected_key + _unit], color_title=color_by.selected_key,
        color_domain=_dom, color_range=_rng, x_title="Posição no cromossomo (GRCh38)", lang="pt",
        head_size=point_size.value * 2, faded=("passa", "não passa"),
        extra_tooltip=[("target_mean", f"Mudança no {profile.target_label} ({A.UNIT.get(profile.scorer, 'log')})", ".4f"),
                       ("control_mean", f"Mudança no {profile.control_label} ({A.UNIT.get(profile.scorer, 'log')})", ".4f"),
                       ("top_feature", "Modalidade AVI", None), ("passa", "Filtros", None)])
    canvas = mo.ui.anywidget(VariantCanvas(spec=chart_spec(_chart), selected=list(get_focus_ids()), crop=list(_crop),
                                           file_name=f"{study_cfg.data['study_id']}_janela_{metric}"))
    _rel = sorted(((_crop[0] - locus.tss) * locus.strand, (_crop[1] - locus.tss) * locus.strand))
    mo.vstack([
        mo.md(f"## Ao longo da janela\n{_df.height:,} variantes" + (f" (as {profile.top_k:,} primeiras em *Ordenar por*)"
              if profile.top_k and _df.height >= profile.top_k else "") + f" · trecho de {_crop[1] - _crop[0] + 1} pb "
              f"({_rel[0]:+d} a {_rel[1]:+d} pb do início do gene). Clique numa variante para destacá-la."),
        mo.accordion({"Mapa da janela (gene, promotor, éxon 1)": mo.Html(figure_svg(locus, footprint_blocks(locus))),
                      "Como usar este gráfico": mo.md(
                          "Cada haste é uma variante; para cima aumenta, para baixo diminui. **Clique** para destacar (anel vermelho; "
                          "clique de novo para tirar). **Roda do mouse** aproxima, **arrastar** desloca e **clique duplo** volta. Pontos sobrepostos: diminua o "
                          "*Tamanho dos pontos*, ligue *Uma alternativa por posição* (fica só a troca de maior efeito em cada "
                          "posição) ou aproxime com a roda. Para ver só um trecho, clique em "
                          "*Recortar* e arraste sobre ele, ou use a barra; *Janela inteira* volta. O trecho escolhido vale também "
                          "para o gráfico de comparação e para a tabela. Para salvar a figura, use os botões *Salvar PNG* e *Salvar SVG* abaixo dele.")}),
        mo.hstack([crop_button, reset_button, clear_button, color_by], justify="start"),
        mo.hstack([point_size, one_alt], justify="start"),
        crop_slider,
        canvas,
        save_buttons(_chart, f"{study_cfg.data['study_id']}_janela"),
    ])
    return (canvas,)


@app.cell
def _(canvas, full_domain, get_crop, get_crop_mode, get_focus_ids, select_variants, set_crop, set_crop_mode):
    _picked = list(canvas.value.get("selected", []))
    if _picked != list(get_focus_ids()):
        select_variants(_picked)
    _new = tuple(canvas.value.get("crop", []) or ())
    if get_crop_mode() and len(_new) == 2 and _new != tuple(get_crop()):
        set_crop((max(_new[0], full_domain[0]), min(_new[1], full_domain[1])))
        set_crop_mode(False)
    return


@app.cell
def _(A, H, VariantCanvas, alt, chart_spec, color_by, get_crop, get_focus_ids, mo, pl, point_size, profile, save_buttons,
      shown_df, study_cfg):
    _lo, _hi = get_crop()
    _schart = H.target_control_scatter(shown_df.filter(pl.col("position").is_between(_lo, _hi)), profile, color_by.value,
                                       color_by.selected_key, point_size=point_size.value, focus_ids=list(get_focus_ids()))
    if _schart is None:
        _schart = alt.Chart(alt.Data(values=[{"t": "Nenhuma variante com valores no alvo e no controle."}])).mark_text().encode(
            text="t:N").properties(width=560, height=60)
    scatter = mo.ui.anywidget(VariantCanvas(spec=chart_spec(_schart),
                                            selected=list(get_focus_ids()), crop=[],
                                            file_name=f"{study_cfg.data['study_id']}_alvo_x_controle"))
    mo.vstack([
        mo.md(f"## O efeito é próprio do alvo?\nCada ponto é uma variante: horizontal = mudança no **{profile.target_label}**, "
              f"vertical = no **{profile.control_label}**. **Roda do mouse** aproxima, **arrastar** desloca, "
              "**clique duplo** volta ao enquadramento original; clique num ponto para destacá-lo."),
        mo.accordion({"Como ler": mo.md(
            "- **Sobre a diagonal tracejada**: muda igual nos dois grupos (efeito geral, não próprio do alvo).\n"
            "- **Perto da linha horizontal do zero**, longe do centro: muda no alvo e quase nada no controle (próprio do alvo).\n"
            "- **Em quadrantes opostos**: os grupos mudam em sentidos contrários; só tem sentido se o gene for expresso nos dois.\n\n"
            "As cores seguem o seletor \"Cor no gráfico da janela\" (ao lado de Recortar) e são as mesmas nos dois gráficos. "
            "Com \"Variantes nos gráficos: todas\", as que não passam nos filtros ficam claras.")}),
        scatter,
        save_buttons(_schart, f"{study_cfg.data['study_id']}_alvo_x_controle", width=560),
    ])
    return (scatter,)


@app.cell
def _(get_focus_ids, scatter, select_variants):
    _picked = list(scatter.value.get("selected", []))
    if _picked != list(get_focus_ids()):
        select_variants(_picked)
    return


@app.cell
def _(A, get_crop, get_focus_ids, labels, mo, pl, profile, select_variants, show_mode, shown_df):
    _base = A.SCALES.get(profile.scorer, (2.718281828459045, ""))[0]
    _lo, _hi = get_crop()
    _rows = shown_df.filter(pl.col("position").is_between(_lo, _hi))
    _eff = [c for c in _rows.columns if c.startswith("eff_")]
    _pct = lambda c: ((pl.lit(_base).pow(pl.col(c)) - 1) * 100).round(1)
    _show = _rows.select([pl.col("variant_id").alias("variante"), pl.col("genomic_domain").alias("região"),
                          pl.col("distance_to_tss").alias("dist. início do gene"),
                          _pct("target_mean").alias(f"% {profile.target_label}"), _pct("control_mean").alias(f"% {profile.control_label}"),
                          *[_pct(c).alias(f"% {labels.get(c[4:], c[4:]).split(' (')[0]}") for c in _eff],
                          pl.col("top_feature").alias("modalidade AVI"), pl.col("atlas_url").alias("Atlas"),
                          *([pl.col("passa").alias("filtros")] if show_mode.value == "all" else [])])
    _ids = _show["variante"].to_list()
    variant_table = mo.ui.table(_show.to_dicts(), selection="multi",
                                initial_selection=[i for i, v in enumerate(_ids) if v in set(get_focus_ids())], page_size=12,
                                on_change=lambda rows: select_variants([r["variante"] for r in rows or []]),
                                label=f"{_show.height:,} variantes no trecho atual")
    _splice = _rows.height and "Splicing" in _rows.head(10)["top_feature"].to_list()
    mo.vstack([
        mo.md(f"## Variantes\nMudança prevista em % d{A.SCALES.get(profile.scorer, (None, 'o sinal previsto'))[1]}, "
              "por grupo e por tecido. Marque linhas para destacar."
              + ("  \n**Atenção**: há variantes de modalidade AVI *Splicing* entre as primeiras: o efeito delas é provavelmente "
                 "de splicing, não de transcrição." if _splice else "")),
        variant_table,
    ])
    return


@app.cell
def _(get_focus_ids, kept, mo, select_variants):
    _opts = list(dict.fromkeys([*get_focus_ids(), *kept["variant_id"].to_list()]))
    pick_by_id = mo.ui.multiselect(options=_opts, value=list(get_focus_ids()), label="Destacar pelo identificador",
                                   max_selections=50, on_change=lambda v: select_variants(list(v)))
    return (pick_by_id,)


@app.cell
def _(H, data, get_focus_ids, labels, mo, pick_by_id, pl, profile, save_buttons, study_cfg):
    _ids = list(get_focus_ids())
    if not _ids:
        _body = mo.md("Destaque uma ou mais variantes (clique nos gráficos, marque na tabela ou escolha pelo identificador).")
    else:
        _long = H.track_effects(data, profile, _ids, labels)
        _links = data.master.filter(pl.col("variant_id").is_in(_ids)).select("variant_id", "atlas_url", "genomic_domain", "top_feature")
        _md = "  \n".join(f"**{r['variant_id']}** · {r['genomic_domain']} · modalidade AVI {r['top_feature']} · "
                          f"[ver no AlphaGenome Atlas]({r['atlas_url']})" for r in _links.iter_rows(named=True))
        _body = mo.vstack([
            mo.md(_md),
            mo.md("Uma barra por medição. Barras do mesmo tecido que discordam muito indicam medição instável naquele tecido."),
            *([_bars := H.variant_tissue_bars(_long, profile),
               save_buttons(_bars, f"{study_cfg.data['study_id']}_variantes_por_tecido", width=900)]
              if _long.height else [mo.md("Sem medições para os tecidos escolhidos.")]),
        ])
    mo.vstack([mo.md(f"## Variantes destacadas ({len(_ids)})"), pick_by_id, _body])
    return


@app.cell
def _(A, mo, steps):
    mo.md("---")
    return


@app.cell
def _(A, H, data, domain_pick, max_control_pct, min_sigma, min_target_pct, mo, profile, steps, use_max_control):
    _ov = A.noise_overall(data, profile.scorer, profile.assays or None, profile.sources or None)
    _p95 = f"{H._n(A.pct(_ov[0], profile.scorer), 1)}%" if _ov is not None else "não disponível (nenhum tecido com duas medições)"
    mo.accordion({"Filtrar resultados": mo.vstack([
        mo.md("Todos começam desligados. Como referência, duas medições do **mesmo** tecido costumam discordar em até "
              f"**{_p95}** neste estudo (95% dos casos): mudanças menores que isso podem ser só diferença entre medições."),
        min_target_pct,
        mo.md("Cada tecido do alvo precisa mudar pelo menos isso, no sentido escolhido. Com vários tecidos, vale o que muda menos."),
        mo.hstack([use_max_control, max_control_pct], justify="start"),
        mo.md("Para efeitos próprios do alvo: nenhum tecido do controle pode mudar mais que isso."),
        min_sigma,
        mo.md("A diferença entre alvo e controle precisa superar este múltiplo do ruído do estudo (2 = além do comum; 3 = bem além). "
              "O ruído é estimado juntando todos os tecidos com duas ou mais medições; um tecido com uma só medição usa esse ruído "
              "geral. Sem controle ou sem ruído, o filtro fica de fora (a tabela de etapas diz quando)."),
        domain_pick,
        mo.md("**Quantas variantes passam em cada etapa:**"),
        mo.ui.table(A.funnel(steps).to_pandas(), selection=None, pagination=False),
    ])})
    return


@app.cell
def _(A, H, assay_pick, control_name, descriptions, flags, floor_applies, labels, mo, noise, pl,
      profile, ref_floor, scorer_pick, shown_pick, source_pick, target_name, tissues, use_floor):
    _t = tissues.select(pl.col("tissue").alias("tecido"), pl.col("ontology_curie").alias("termo"),
                        pl.col("study_cohort").replace({"target": "alvo", "control": "controle"}).alias("na captura"),
                        pl.col("tracks").alias("medições"), pl.col("ref_min").round(4).alias("REF mín"),
                        pl.col("ref_max").round(4).alias("REF máx"),
                        pl.col("ontology_curie").replace_strict(descriptions, default="").alias("de onde vêm as medições"))
    if noise is None:
        _noise_md = "Nenhum tecido tem duas medições com os ensaios e fontes escolhidos: não há como medir o ruído."
    else:
        _rows = [f"| até {r['hi']:.3f} | {r['p95']:.3f} | {r['p99']:.3f} | {r['sigma']:.3f} |" if r["hi"] != float("inf")
                 else f"| acima | {r['p95']:.3f} | {r['p99']:.3f} | {r['sigma']:.3f} |" for r in noise.iter_rows(named=True)]
        _noise_md = (f"| tamanho da mudança ({A.UNIT.get(profile.scorer, 'log')}, em módulo) | p95 | p99 | σ |\n"
                     "|---|---|---|---|\n" + "\n".join(_rows))
    _warn = [f"- **{labels.get(c, c)}**: abaixo do piso de expressão" for c, f in flags.items() if f == "below floor"]
    mo.accordion({"Avançado": mo.vstack([
        mo.md("### Grupos"),
        mo.hstack([target_name, control_name], justify="start"),
        shown_pick,
        mo.md(H.group_legend_md(profile, labels)),
        mo.md("### Medições que entram"),
        mo.md("Misturar polyA com total, ou ENCODE com GTEx, junta medições que às vezes discordam muito."),
        mo.hstack([assay_pick, source_pick], justify="start"),
        mo.md("### Expressão de referência"),
        mo.md("REF = expressão do gene prevista no tecido, sem a variante. Perto de 0,001 (o pseudocount do modelo) o gene está "
              "praticamente desligado e as mudanças previstas ficam instáveis. Na ativação ectópica, use \"Só avisar\"."),
        mo.hstack([use_floor, ref_floor, floor_applies], justify="start"),
        mo.md("\n".join(_warn)),
        mo.md("### Tecidos capturados"),
        mo.ui.table(_t.to_pandas(), selection=None, pagination=False),
        mo.md("### Ruído do estudo"),
        mo.md("Quanto duas medições do mesmo tecido discordam, por tamanho da mudança (na escala bruta do modelo). σ é a escala típica "
              "usada no filtro \"acima do ruído\"."),
        mo.md(_noise_md),
        mo.md("### Tipo de predição"),
        scorer_pick,
        mo.md("RNA_SEQ = expressão do gene; CAGE = início de transcrição; DNASE/ATAC = cromatina aberta."),
    ])})
    return


@app.cell
def _(kept, mo):
    _n = max(1, kept.height)
    matrix_order = mo.ui.dropdown(options={"Posição na janela": "position", "Métrica do gráfico da janela": "metric"},
                                  value="Posição na janela", label="Ordem das linhas")
    matrix_rows = mo.ui.slider(start=1, stop=min(_n, 200), step=1, value=min(_n, 60), label="Linhas", show_value=True)
    return matrix_order, matrix_rows


@app.cell
def _(H, get_crop, get_focus_ids, kept, labels, matrix_order, matrix_rows, metric, mo, pl, profile, save_buttons, study_cfg):
    _lo, _hi = get_crop()
    _in = kept.filter(pl.col("position").is_between(_lo, _hi))
    if _in.is_empty():
        _body = mo.md("Nenhuma variante passa nos filtros no trecho atual.")
    else:
        _chart, _meta = H.tissue_matrix(_in, profile, labels, max_rows=matrix_rows.value, order_by=matrix_order.value,
                                        metric_col=metric, picked=get_focus_ids())
        _body = mo.vstack([
            mo.md(f"{_meta['n_variants']} de {_meta['n_available']} variantes (linhas) × tecidos (colunas). Vermelho = aumenta, "
                  "azul = diminui. Tecido em azul = alvo, cinza = controle, bege = contexto; variantes destacadas em vermelho."),
            mo.hstack([matrix_order, matrix_rows], justify="start"),
            _chart,
            save_buttons(_chart, f"{study_cfg.data['study_id']}_matriz")])
    mo.accordion({"Comparar todos os tecidos (matriz variante × tecido)": _body})
    return


if __name__ == "__main__":
    app.run()
