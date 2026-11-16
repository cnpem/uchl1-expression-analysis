"""Figures and texts of the analysis dashboard (assets/dashboard.py) built from an analysis profile.

The legacy workbench helpers read a fixed contrast plan; these read the profile's groups instead (target,
control and context tissues chosen in the dashboard), so every figure follows what the researcher picked.

Library: GROUP_COLORS, group_of, color_options, add_display_columns, track_effects, tissue_matrix,
variant_tissue_bars, group_legend_md.
"""

from __future__ import annotations

from typing import Any, Optional, Sequence

import altair as alt
import polars as pl

from core import analysis as A
from viz.design_system import apply_print_style
from viz.workbench_helpers import CONTROL_COLOR, TARGET_COLOR

CONTEXT_COLOR = "#a8a29e"
GROUP_COLORS = {"alvo": TARGET_COLOR, "controle": CONTROL_COLOR, "contexto": CONTEXT_COLOR}


def group_of(profile: A.Profile, curie: str) -> str:
    if curie in profile.target:
        return "alvo"
    if curie in profile.control:
        return "controle"
    return "contexto"


def color_options() -> dict[str, str]:
    """Display label -> column used to colour Figure 2."""
    return {"Domínio genômico": "genomic_domain", "Passa nos filtros": "passa", "Sentido do efeito no alvo": "sentido",
            "Região do transcrito": "ensembl_transcript_region", "Modalidade principal do AVI": "top_feature"}


def add_display_columns(effects: pl.DataFrame, kept_ids: Sequence[str]) -> pl.DataFrame:
    """Columns used for colour and tooltips: whether the variant passes the filters and the sign in the target."""
    return effects.with_columns(
        pl.when(pl.col("variant_id").is_in(list(kept_ids))).then(pl.lit("passa")).otherwise(pl.lit("não passa")).alias("passa"),
        pl.when(pl.col("target_mean") > 0).then(pl.lit("aumenta no alvo"))
        .when(pl.col("target_mean") < 0).then(pl.lit("diminui no alvo")).otherwise(pl.lit("sem efeito")).alias("sentido"))


def fixed_colors(color_col: str) -> tuple[Optional[list[str]], Optional[list[str]]]:
    if color_col == "passa":
        return ["passa", "não passa"], ["#2a78d6", "#d6d3d1"]
    if color_col == "sentido":
        return ["aumenta no alvo", "diminui no alvo", "sem efeito"], ["#c2410c", "#1d4ed8", "#a8a29e"]
    return None, None


def track_effects(data: A.StudyData, profile: A.Profile, variant_ids: Sequence[str],
                  labels: dict[str, str]) -> pl.DataFrame:
    """Long table: one row per (variant, track) of the profile's tissues, with group, assay, source and value."""
    t = A.tracks(data, profile.scorer, profile.assays or None, profile.sources or None)
    tissues = list(dict.fromkeys([*profile.target, *profile.control, *profile.shown]))
    t = t.filter(pl.col("ontology_curie").is_in(tissues))
    sub = data.master.filter(pl.col("variant_id").is_in(list(variant_ids)))
    rows = []
    for r in t.iter_rows(named=True):
        col = f"{r['column_prefix']}_raw"
        if col not in sub.columns:
            continue
        for vid, val in zip(sub["variant_id"].to_list(), sub[col].to_list()):
            rows.append({"variant_id": vid, "curie": r["ontology_curie"],
                         "tecido": labels.get(r["ontology_curie"], r["ontology_curie"]).split(" (")[0],
                         "grupo": group_of(profile, r["ontology_curie"]),
                         "trilha": f"{r['assay_title']} · {r['data_source']}", "efeito": val})
    return pl.DataFrame(rows) if rows else pl.DataFrame()


def tissue_matrix(effects: pl.DataFrame, profile: A.Profile, labels: dict[str, str], *, max_rows: int = 60,
                  order_by: str = "position", metric_col: Optional[str] = None,
                  picked: Optional[Sequence[str]] = None) -> tuple[Any, dict[str, Any]]:
    """Figure 3: variants (rows) x the profile's tissues (columns: target, control, context), LFC per tissue."""
    tissues = [c for c in dict.fromkeys([*profile.target, *profile.control, *profile.shown]) if f"eff_{c}" in effects.columns]
    rows = (effects.sort(pl.col(metric_col).abs(), descending=True, nulls_last=True) if order_by == "metric" and metric_col
            else effects.sort("position")).head(max_rows)
    ids = rows["variant_id"].to_list()
    names = {c: labels.get(c, c).split(" (")[0] for c in tissues}
    long = rows.select(["variant_id", *[f"eff_{c}" for c in tissues]]).unpivot(index="variant_id", variable_name="col",
                                                                                value_name="efeito")
    long = long.with_columns(pl.col("col").str.slice(4).alias("curie")).with_columns(
        pl.col("curie").replace_strict(names, default=None).alias("tecido"),
        pl.col("curie").replace_strict({c: group_of(profile, c) for c in tissues}, default="contexto").alias("grupo"))
    order = [names[c] for c in tissues]
    by_group = {g: [names[c] for c in tissues if group_of(profile, c) == g] for g in GROUP_COLORS}
    label_expr = (f"indexof({by_group['alvo']!r}, datum.value) >= 0 ? '{TARGET_COLOR}' : "
                  f"(indexof({by_group['controle']!r}, datum.value) >= 0 ? '{CONTROL_COLOR}' : '{CONTEXT_COLOR}')")
    pk = list(picked or [])
    row_color = f"indexof({pk!r}, datum.value) >= 0 ? '#e11d48' : '#334155'" if pk else "'#334155'"
    lim = float(long["efeito"].abs().max() or 0.01)
    chart = alt.Chart(long.to_pandas()).mark_rect().encode(
        x=alt.X("tecido:N", sort=order, title=None,
                axis=alt.Axis(orient="top", labelAngle=-35, labelLimit=320, labelFontSize=13, labelColor=alt.ExprRef(label_expr))),
        y=alt.Y("variant_id:N", sort=ids, title=None, axis=alt.Axis(labelColor=alt.ExprRef(row_color), labelFontSize=12)),
        color=alt.Color("efeito:Q", scale=alt.Scale(scheme="redblue", reverse=True, domain=[-lim, lim]),
                        title="Efeito previsto (LFC)", legend=alt.Legend(gradientLength=220)),
        tooltip=["variant_id:N", "tecido:N", "grupo:N", alt.Tooltip("efeito:Q", format=".4f")],
    ).properties(width=alt.Step(70), height=alt.Step(18))
    return apply_print_style(chart), {"n_variants": len(ids), "n_available": effects.height}


def variant_tissue_bars(tracks_long: pl.DataFrame, profile: A.Profile) -> Any:
    """Figure 4: for the picked variants, the effect in every track of every tissue, coloured by group.

    One bar per track (polyA, total, GTEx...), so the agreement between tracks of the same tissue, which
    is what the noise model measures, is visible; several variants show their mean with min-max in the tooltip.
    """
    s = (tracks_long.group_by(["tecido", "grupo", "trilha"], maintain_order=True)
         .agg(pl.col("efeito").mean().alias("efeito"), pl.col("efeito").min().alias("mínimo"),
              pl.col("efeito").max().alias("máximo"), pl.len().alias("variantes")))
    order = list(dict.fromkeys(s.sort(pl.col("grupo").replace_strict({"alvo": 0, "controle": 1, "contexto": 2}))["tecido"].to_list()))
    chart = alt.Chart(s.to_pandas()).mark_bar(cornerRadiusEnd=3).encode(
        y=alt.Y("tecido:N", sort=order, title=None),
        yOffset=alt.YOffset("trilha:N"),
        x=alt.X("efeito:Q", title="Efeito previsto (LFC)"),
        color=alt.Color("grupo:N", scale=alt.Scale(domain=list(GROUP_COLORS), range=list(GROUP_COLORS.values())),
                        legend=alt.Legend(title="Grupo")),
        opacity=alt.Opacity("trilha:N", legend=alt.Legend(title="Trilha")),
        tooltip=["tecido:N", "grupo:N", "trilha:N", alt.Tooltip("efeito:Q", format=".4f"), alt.Tooltip("mínimo:Q", format=".4f"),
                 alt.Tooltip("máximo:Q", format=".4f"), "variantes:Q"],
    ).properties(width="container", height=max(260, 24 * max(1, s.height)))
    return apply_print_style(chart)


def group_legend_md(profile: A.Profile, labels: dict[str, str]) -> str:
    name = lambda c: labels.get(c, c)
    lines = [f"<span style='color:{TARGET_COLOR}'>■</span> **{profile.target_label}** (alvo): "
             + (", ".join(name(c) for c in profile.target) or "nenhum tecido escolhido"),
             f"<span style='color:{CONTROL_COLOR}'>■</span> **{profile.control_label}** (controle): "
             + (", ".join(name(c) for c in profile.control) or "nenhum tecido escolhido")]
    if profile.shown:
        lines.append(f"<span style='color:{CONTEXT_COLOR}'>■</span> **Contexto** (só exibidos): "
                     + ", ".join(name(c) for c in profile.shown))
    return "  \n".join(lines)


# --------------------------------------------------------------------------- plain-language texts

_SOURCE = {"gtex": "GTEx (tecido de adulto post-mortem, média de muitos doadores; mistura de tipos celulares)",
           "encode": "ENCODE", "fantom": "FANTOM5"}
_ASSAY = {"polyA plus RNA-seq": "polyA (RNA com cauda poli-A, sobretudo mRNA maduro)",
          "total RNA-seq": "total (inclui RNA ainda não processado e não codificante)"}
_KIND = {"tissue": "tecido", "primary_cell": "célula primária", "cell_line": "linhagem celular",
         "in_vitro_differentiated_cells": "célula diferenciada em laboratório"}


def describe_tissues(data: A.StudyData, scorer: str) -> dict[str, str]:
    """Ontology term -> one plain sentence on its measurements and how much they agree in this window."""
    t = A.tracks(data, scorer)
    out: dict[str, str] = {}
    for curie, grp in t.group_by("ontology_curie"):
        curie = curie[0] if isinstance(curie, tuple) else curie
        parts = []
        for r in grp.iter_rows(named=True):
            src = _SOURCE.get(r["data_source"], r["data_source"] or "fonte não informada")
            assay = _ASSAY.get(r["assay_title"], r["assay_title"] or "")
            kind = _KIND.get(r["biosample_type"], r["biosample_type"] or "")
            parts.append(f"{src}, {kind}{', ' + assay if assay else ''}")
        text = f"{grp.height} {'medições' if grp.height > 1 else 'medição'}: " + "; ".join(parts) + "."
        cols = [f"{p}_raw" for p in grp["column_prefix"].to_list() if f"{p}_raw" in data.master.columns]
        if len(cols) >= 2:
            # the least concordant pair: one unstable measurement is enough to make the tissue unstable
            rs = [data.master.select(pl.corr(a, b)).item() for i, a in enumerate(cols) for b in cols[i + 1:]]
            rs = [r for r in rs if r is not None and r == r]
            if rs:
                r = min(rs)
                how = "concordam bem" if r >= 0.9 else ("concordam em parte" if r >= 0.6 else "concordam pouco")
                which = "" if len(cols) == 2 else ", o par que menos concorda"
                text += (f" Neste gene, as medições {how} (correlação {_n(r, 2)}{which})"
                         + ("; compare com cuidado." if r < 0.6 else "."))
        out[curie] = text
    return out


def _n(x: float, digits: int = 1, sign: bool = False) -> str:
    """Number in Brazilian format (comma decimal, dot thousands)."""
    s = f"{x:+,.{digits}f}" if sign else f"{x:,.{digits}f}"
    return s.replace(",", "§").replace(".", ",").replace("§", ".")


def headline(kept: pl.DataFrame, effects: pl.DataFrame, profile: A.Profile, labels: dict[str, str],
             ref_by_tissue: dict[str, float], top: int = 3) -> str:
    """One paragraph answering the current question with the current choices (Portuguese)."""
    q = {"repression": "diminuem", "increase": "aumentam", "either": "mudam"}[profile.question]
    total = effects.height
    lines = [f"**{_n(kept.height, 0)} de {_n(total, 0)} variantes** passam nas escolhas atuais "
             f"(alvo: {profile.target_label}; controle: {profile.control_label})."]
    if kept.height and "target_mean" in kept.columns:
        col = "target_toward" if profile.question != "either" else "target_mean"
        best = kept.sort(pl.col(col).abs() if profile.question == "either" else pl.col(col), descending=True, nulls_last=True).head(top)
        items = []
        for r in best.iter_rows(named=True):
            t = A.pct(r["target_mean"], profile.scorer) if r["target_mean"] is not None else float("nan")
            c = A.pct(r["control_mean"], profile.scorer) if r.get("control_mean") is not None else None
            items.append(f"`{r['variant_id']}` ({_n(t, 1, True)}% no alvo" + (f", {_n(c, 1, True)}% no controle)" if c is not None else ")"))
        what = A.SCALES.get(profile.scorer, (None, "o sinal previsto"))[1]
        lines.append(f"Variantes que mais {q} {what} no alvo: " + "; ".join(items) + ".")
    # REF is the gene's predicted expression: the warning concerns RNA_SEQ effects on a gene that is off.
    # Reference level: the researcher's floor when set, else 0.01 (ten times the model's pseudocount, 0.001).
    floor = profile.ref_floor if profile.ref_floor is not None else 0.01
    low = ([labels.get(c, c).split(" (")[0] for c in profile.target if c in ref_by_tissue and ref_by_tissue[c] < floor]
           if profile.scorer == "RNA_SEQ" else [])
    if low:
        lines.append(f"**Atenção**: o gene tem expressão prevista quase nula em {', '.join(low)} (REF < {_n(floor, 3)}; "
                     "o pseudocount do modelo é 0,001): efeitos previstos ali são instáveis; só interprete como possível "
                     "ativação ectópica.")
    return "  \n".join(lines)
