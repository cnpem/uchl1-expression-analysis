"""Dashboard figures for the report (`expression.py figures`, Workflow step 14).

The dashboard saves nothing. When the researcher settles on a view, the agent passes the same choices to
`expression.py figures`, which runs the dashboard's computation (core.analysis) and draws the same charts
(dashboard_helpers, workbench_helpers) as static PNG files. `report.md` then shows them with the choices
that produced them, so every figure in the report can be redrawn and checked in the dashboard.

Outputs (figures/dashboard/): janela.png (variants along the window), alvo_x_controle.png (target vs
control), variantes.png (per-measurement bars of the highlighted variants), matriz.png (variants x
tissues) and choices.json (the choices, the headline, the filter steps and the files).

Library: FIG_DIR, CHOICES_FILE, write_figures, load_choices, describe_choices.
"""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import polars as pl

from core import analysis as A
from viz import dashboard_helpers as H
from viz.workbench_helpers import build_unified_lollipop_canvas

FIG_DIR = "figures/dashboard"
CHOICES_FILE = f"{FIG_DIR}/choices.json"
QUESTION_PT = {"repression": "diminui no alvo", "increase": "aumenta no alvo", "either": "qualquer mudança"}
SHOW_PT = {"kept": "só as que passam nos filtros", "all": "todas (as que não passam em cinza)"}


def _png(chart: Any, path: Path, width: int = 1100) -> None:
    import vl_convert as vlc

    spec = chart.to_dict()
    if spec.get("width") == "container":
        spec = chart.properties(width=width).to_dict()
    path.write_bytes(vlc.vegalite_to_png(spec, scale=2))


def write_figures(study_dir: Path, data: A.StudyData, profile: A.Profile, effects: pl.DataFrame, kept: pl.DataFrame,
                  steps: Sequence[tuple[str, int]], *, show_mode: str, color_col: str, highlight: Sequence[str],
                  labels: Mapping[str, str], window: tuple[int, int], command: str,
                  opening: bool = False) -> dict[str, Any]:
    """Draws the dashboard's charts for these choices into figures/dashboard/ and writes choices.json.

    `opening` marks the study's opening view (its groups and question, no filter), drawn by `report` when the
    researcher has not chosen a view yet; a later `figures` run replaces it."""
    out = Path(study_dir) / FIG_DIR
    out.mkdir(parents=True, exist_ok=True)
    metric = profile.y_metric if profile.y_metric in effects.columns else "target_mean"
    metric_label = next((k for k, v in H.Y_METRICS.items() if v == metric), metric)
    color_label = next((k for k, v in H.color_options().items() if v == color_col), color_col)
    shown = H.display_frame(effects, kept, profile, show_mode)
    picked = [v for v in highlight if v in set(effects["variant_id"].to_list())]
    files: dict[str, str] = {}

    dom, rng = H.fixed_colors(color_col)
    unit = "" if metric in ("avi_phred", "diff_sigma") else f" ({A.unit(profile)})"
    if shown.height:
        lolli = build_unified_lollipop_canvas(
            shown, metric, metric_label, color_col, max_lollipops=shown.height, selected_variant_ids=picked,
            crop_domain=window, full_domain=window, target_label=profile.target_label, drag_mode="pan",
            axis_title=[metric_label + unit], color_title=color_label, color_domain=dom, color_range=rng,
            x_title="Posição no cromossomo (GRCh38)", lang="pt", head_size=70, faded=("passa", "não passa"),
            extra_tooltip=[("top_feature", "Modalidade AVI", None), ("passa", "Filtros", None)])
        _png(lolli, out / "janela.png")
        files["janela"] = f"{FIG_DIR}/janela.png"
        scatter = H.target_control_scatter(shown, profile, color_col, color_label, focus_ids=picked, zoom=False)
        if scatter is not None:
            _png(scatter, out / "alvo_x_controle.png", width=560)
            files["alvo_x_controle"] = f"{FIG_DIR}/alvo_x_controle.png"
    if picked:
        long = H.track_effects(data, profile, picked, dict(labels))
        if long.height:
            _png(H.variant_tissue_bars(long, profile), out / "variantes.png", width=900)
            files["variantes"] = f"{FIG_DIR}/variantes.png"
    if kept.height:
        matrix, _meta = H.tissue_matrix(kept, profile, dict(labels), max_rows=40, order_by="metric", metric_col=metric,
                                        picked=picked)
        _png(matrix, out / "matriz.png")
        files["matriz"] = f"{FIG_DIR}/matriz.png"
    for name in ("janela", "alvo_x_controle", "variantes", "matriz"):  # a previous run's figure would mislead
        if name not in files and (out / f"{name}.png").is_file():
            (out / f"{name}.png").unlink()

    record = {
        "command": command,
        "opening": opening,
        "profile": dataclasses.asdict(profile),
        "display": {"show_mode": show_mode, "color_by": color_col, "highlight": picked, "metric_label": metric_label},
        "labels": {c: labels.get(c, c) for c in [*profile.target, *profile.control, *profile.shown]},
        "headline": H.headline(kept, effects, profile, dict(labels), A.ref_by_tissue(data),
                               missing=A.without_tracks(data, profile)),
        "steps": [list(s) for s in steps],
        "n_shown": shown.height,
        "figures": files,
    }
    (Path(study_dir) / CHOICES_FILE).write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return record


def load_choices(study_dir: Path) -> dict[str, Any] | None:
    path = Path(study_dir) / CHOICES_FILE
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def describe_choices(record: Mapping[str, Any]) -> list[str]:
    """The choices behind the figures, in the words of the dashboard (Portuguese, for report.md)."""
    p, d, lab = record["profile"], record["display"], record["labels"]
    name = lambda c: f"{lab.get(c, c).split(' (')[0]} (`{c}`)"
    pct = lambda v: "desligado" if not v else f"{v:g}%"
    rows = [
        ("O que você procura?", QUESTION_PT.get(p["question"], p["question"])),
        (f"Onde ({p['target_label']})", ", ".join(name(c) for c in p["target"]) or "—"),
        (f"Comparado com ({p['control_label']})", ", ".join(name(c) for c in p["control"]) or "—"),
        ("Tipo de predição", p["scorer"]),
        ("Pontuação", "quantil do Atlas" if p.get("score") == "quantile" else "valor bruto do modelo"),
        ("Como juntar os tecidos de cada grupo", "mediana" if p.get("group_agg") == "median" else "média"),
        ("Mudança mínima em cada tecido do alvo", pct(p["min_change_pct"])),
        ("Mudança máxima no controle", "desligado" if p["max_control_pct"] is None else f"{p['max_control_pct']:g}%"),
        ("Alvo − controle acima do ruído", "desligado" if not p["min_diff_sigma"] else f"{p['min_diff_sigma']:g} σ"),
        ("Regiões da janela", ", ".join(p["domains"]) or "todas"),
        ("Ensaios / fontes", f"{', '.join(p['assays']) or 'todos'} / {', '.join(p['sources']) or 'todas'}"),
        ("Piso de expressão (REF)", "desligado" if p["ref_floor"] is None else f"{p['ref_floor']:g} ({p['ref_floor_applies']})"),
        ("Uma alternativa por posição", "sim" if p["one_alt_per_position"] else "não"),
        ("Ordenar por", d["metric_label"]),
        ("Quantas mostrar", f"{p['top_k']}" if p["top_k"] else "todas"),
        ("Variantes nos gráficos", SHOW_PT.get(d["show_mode"], d["show_mode"])),
        ("Variantes destacadas", ", ".join(f"`{v}`" for v in d["highlight"]) or "nenhuma"),
    ]
    return ["| Escolha | Valor |", "|---|---|", *[f"| {k} | {v} |" for k, v in rows]]
