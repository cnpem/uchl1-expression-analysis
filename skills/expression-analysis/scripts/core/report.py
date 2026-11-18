"""Research report of a study (Workflow step 13, `expression.py report`).

Writes <study-dir>/report.md in Portuguese with seven sections: Estudo, Decisões, Análise do
pesquisador (dashboard), Resultados da triagem fixa, Skills utilizadas, Limitações and Citações.
Everything comes from files of the study directory: study.yaml, decisions.yaml, manifest.json,
processed/*, figures/dashboard/ (the researcher's view drawn by `expression.py figures`) and, when
present, raw/atlas_links.md and raw/interpret_*.txt (optional: included only if produced).
The skills listed are the ones recorded in manifest.json, so the report names what was actually
used. Nothing here is specific to one gene. Research use only: the report describes predicted
molecular impact and never a diagnosis.

Library: ReportError, write_report, render_report.
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import polars as pl

from core import cohorts as core_cohorts
from core import contrast as core_contrast
from core import locus as core_locus
from core import master as core_master
from study import load_decisions, load_study, option_codes, option_label

REPORT_FILE = "report.md"
MANIFEST_FILE = "manifest.json"
ATLAS_LINKS_FILE = "raw/atlas_links.md"
INTERPRET_GLOB = "interpret_*.txt"
MAX_INTERPRET_LINES = 60
PRIMARY_CITATION = ("Avsec, Ž. et al. (2026). Advancing regulatory variant effect prediction with AlphaGenome. "
                    "Nature. (AlphaGenome model and Atlas predictions used throughout this study.)")

# Citation of each dependency skill (key = `skill` field of manifest.json steps).
SKILL_CITATIONS: dict[str, str] = {
    "alphagenome-variant-impact-score": (
        "DeepMind science skills, `alphagenome-variant-impact-score` (GENCODE gene structure, AlphaGenome Atlas "
        "variant impact scores and track catalogue). https://github.com/google-deepmind/science-skills. "
        "Underlying model: Avsec et al. (2026)."),
    "alphagenome-single-variant-analysis": (
        "DeepMind science skills, `alphagenome-single-variant-analysis` (per-variant AlphaGenome predictions). "
        "https://github.com/google-deepmind/science-skills. Underlying model: Avsec et al. (2026)."),
    "alphagenome-atlas-website-links": (
        "DeepMind science skills, `alphagenome-atlas-website-links` (Atlas deep links). "
        "https://github.com/google-deepmind/science-skills."),
    "embl-ebi-ols": (
        "DeepMind science skills, `embl-ebi-ols`: EMBL-EBI Ontology Lookup Service (OLS4) for UBERON/CL terms. "
        "https://github.com/google-deepmind/science-skills; https://www.ebi.ac.uk/ols4."),
    "ensembl-database": (
        "DeepMind science skills, `ensembl-database`: Ensembl REST API (Regulatory Build, reference sequence). "
        "https://github.com/google-deepmind/science-skills; https://rest.ensembl.org."),
    "alphagenome": (
        "K-Dense scientific agent skills, `alphagenome` (AlphaGenome API: per-track scores). "
        "https://github.com/K-Dense-AI/scientific-agent-skills. Underlying model: Avsec et al. (2026)."),
}
LOCAL_SKILL = "expression-analysis"


class ReportError(Exception):
    """A study file the report needs is missing or unreadable; the CLI prints it and exits with 1."""


def _read_json(path: Path, step: str) -> Any:
    if not path.is_file():
        raise ReportError(f"{path.name} not found in {path.parent}; {step}.")
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ReportError(f"{path}: invalid JSON ({exc}).") from None


def _cell(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value).replace("|", "\\|").replace("\n", " ")


def _table(headers: Sequence[str], rows: Sequence[Sequence[Any]]) -> list[str]:
    out = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    out += ["| " + " | ".join(_cell(c) for c in row) + " |" for row in rows]
    return out


def _skills_used(manifest: Mapping[str, Any]) -> list[tuple[str, list[str]]]:
    """Skills of the manifest steps, in first-use order, with the step names that used them."""
    used: dict[str, list[str]] = {}
    for step in manifest.get("steps", []):
        used.setdefault(step["skill"], []).append(step["step"])
    return list(used.items())


# --- sections -------------------------------------------------------------------------------


def _section_study(study: Any, locus: core_locus.LocusConfig, resolution: Mapping[str, Any],
                   manifest: Mapping[str, Any]) -> list[str]:
    cfg = study.data
    window = cfg["window"]
    strand = "+" if locus.strand == 1 else "-"
    gene_version = f".{locus.gene_version}" if locus.gene_version is not None else ""
    tx_version = f".{locus.transcript_version}" if locus.transcript_version is not None else ""
    lines = [
        "## Estudo", "",
        f"- Gene: **{locus.gene_symbol}** (`{locus.gene_id}{gene_version}`), {study.assembly}, fita {strand}.",
        f"- Transcrito: `{locus.transcript_id}{tx_version}` ({locus.transcript_type}"
        f"{', MANE Select' if locus.is_mane_select else ''}); TSS {locus.chromosome}:{locus.tss:,}; "
        f"início do CDS {locus.cds_start_genomic:,}.",
        f"- Janela de mutagênese saturada: `{locus.window_str}` ({locus.window.length} pb, "
        f"{3 * locus.window.length:,} SNVs; política `{locus.window_policy}`).",
        "- Features do Ensembl Regulatory Build na janela: "
        + (", ".join(f"`{f.id}` ({f.biotype}, {f.start:,}-{f.end:,})" for f in locus.features) or "nenhuma") + ".",
        f"- Coorte alvo: **{resolution['target_label']}**, CURIEs efetivos: "
        f"{', '.join(f'`{c}`' for c in resolution['effective_target_curies'])}.",
        f"- Coorte controle: **{resolution['control_label']}**, CURIEs efetivos: "
        f"{', '.join(f'`{c}`' for c in resolution['effective_control_curies'])}.",
        f"- Scorers do Atlas: {', '.join(resolution['effective_scorers'])}"
        + (f" (descartados: {', '.join(resolution['dropped_scorers'])})" if resolution["dropped_scorers"] else "") + ".",
    ]
    dropped = [c for c in resolution["curies"] if c["status"] != "kept"]
    for c in dropped:
        label = (c.get("ols") or {}).get("label") or c["curie"]
        extra = f" (substituído por `{c['replaced_by']}`)" if c.get("replaced_by") else ""
        lines.append(f"- CURIE `{c['curie']}` ({label}), coorte {c['role']}: {c['status']}{extra}.")
    origin = manifest.get("origin", {})
    lines += ["", f"Origem dos dados brutos: `{origin.get('kind', 'desconhecida')}`"
              + (f". {origin['note']}" if origin.get("note") else "."), ""]
    return lines


def _section_decisions(study_dir: Path, meta: Mapping[str, Any]) -> list[str]:
    decisions = load_decisions(study_dir)
    lines = ["## Decisões", ""]
    if not decisions:
        return lines + ["Nenhuma decisão registrada em `decisions.yaml`.", ""]
    ordered = sorted(decisions.items(), key=lambda kv: (kv[1].workflow_step, kv[0]))
    asked = [(k, e) for k, e in ordered if e.decided_by == "researcher"]
    rules = [(k, e) for k, e in ordered if e.decided_by == "skill"]

    def shown(options: Any, code: Any) -> str:
        if isinstance(code, list):
            return ", ".join(map(str, code))
        label = option_label(options, code)
        return f"{label} (`{code}`)" if label else f"`{code}`"

    def origin(e: Any) -> str:
        if not e.said:
            return "—"
        words = " ".join(str(e.said).split())
        words = words if len(words) <= 160 else words[:157] + "…"
        return f"Delegada ao agente: “{words}”" if e.delegated else f"“{words}”"

    if asked:
        lines += ["### Escolhas do pesquisador", "",
                  "Respostas registradas em `decisions.yaml`, com as opções oferecidas e a recomendada; o código entre "
                  "crases é o valor gravado. \"—\" quando a pergunta não tinha opções fixas. *Nas palavras do "
                  "pesquisador* cita de onde veio cada resposta; \"Delegada ao agente\" quando o pesquisador pediu ao "
                  "agente que escolhesse:", ""]
        rows = [(e.workflow_step, f"`{key}`",
                 "; ".join(shown(e.options, c) for c in option_codes(e.options)) if e.options else "—",
                 shown(e.options, e.recommended) if e.recommended is not None else "—",
                 shown(e.options, e.value), origin(e), e.note or "—", e.answered_at[:10]) for key, e in asked]
        lines += _table(["Passo", "Decisão", "Opções", "Recomendada", "Escolhida", "Nas palavras do pesquisador", "Nota",
                         "Data"], rows)
        delegated = [key for key, e in asked if e.delegated]
        if delegated:
            lines += ["", f"**{len(delegated)} resposta(s) delegada(s) ao agente:** "
                      + ", ".join(f"`{k}`" for k in delegated) + ". O pesquisador pediu que o agente escolhesse; as "
                      "palavras dele estão na tabela."]
    redefined = [(k, e) for k, e in ordered if e.previous]
    if redefined:
        lines += ["", "### Redefinições", "",
                  "Decisões que o pesquisador mudou depois de respondidas (por exemplo, os grupos trocados depois dos "
                  "resultados). Os valores anteriores ficam em `decisions.yaml` (`previous`); a tabela acima mostra o "
                  "atual:", ""]
        rows = []
        for key, e in redefined:
            for old in e.previous:
                rows.append((f"`{key}`", shown(old.get("options"), old.get("value")), old.get("note") or "—",
                             str(old.get("answered_at", ""))[:10] or "—"))
            rows.append((f"`{key}`", f"**{shown(e.options, e.value)}** (atual)", e.note or "—", e.answered_at[:10]))
        lines += _table(["Decisão", "Valor", "Nota", "Data"], rows)
    if rules:
        lines += ["", "### Regras aplicadas pela skill", "",
                  "Escolhas feitas por regra fixa, sem pergunta (verificações técnicas e o plano da triagem), "
                  "registradas em `decisions.yaml` com `decided_by: skill`:", ""]
        rows = [(e.workflow_step, f"`{key}`", shown(None, e.value), e.rule or "—", e.answered_at[:10]) for key, e in rules]
        lines += _table(["Passo", "Item", "Aplicado", "Regra", "Data"], rows)
    plan = meta["plan"]
    lines += ["", "Plano de contrastes da triagem (`contrasts.yaml`):", ""]
    rows = [(c["name"], c["scorer"], c["strand"], c["tissue_aggregation"], len(c["target"]), len(c["control"]))
            for c in plan["contrasts"]]
    lines += _table(["Contraste", "Scorer", "Fita", "Agregação por tecido", "Tecidos alvo", "Tecidos controle"], rows)
    lines += ["", f"Rótulos de prioridade a partir do contraste `{plan['priority']['direction_contrast']}`.", ""]
    return lines


FIGURE_TEXT = {
    "janela": ("Variantes ao longo da janela",
               "Uma haste por variante, na sua posição; a altura é a medida escolhida em *Ordenar por*. As que não passam "
               "nas escolhas ficam em cinza (com *Variantes nos gráficos: todas*); as destacadas têm um anel vermelho."),
    "alvo_x_controle": ("O efeito é próprio do alvo?",
                        "Cada ponto é uma variante: mudança no alvo (horizontal) contra o controle (vertical). Sobre a "
                        "diagonal tracejada o efeito é igual nos dois grupos; perto da linha horizontal do zero e longe do "
                        "centro, é próprio do alvo; em quadrantes opostos, os grupos mudam em sentidos contrários."),
    "variantes": ("Variantes destacadas, medição por medição",
                  "Uma barra por medição (ensaio e fonte) de cada tecido. Barras do mesmo tecido que discordam muito "
                  "indicam um tecido instável para a comparação."),
    "matriz": ("Variantes × tecidos",
               "As primeiras variantes que passam (linhas, na ordem escolhida) em cada tecido (colunas): vermelho aumenta, "
               "azul diminui."),
}


def _section_dashboard(study_dir: Path) -> list[str]:
    """The researcher's view of the dashboard, drawn by `expression.py figures` with the choices behind it."""
    from viz import report_figures

    rec = report_figures.load_choices(study_dir)
    opening = bool(rec and rec.get("opening"))
    lines = ["## Visão inicial do estudo (dashboard)" if opening else "## Análise do pesquisador (dashboard)", ""]
    if opening:
        lines += ["Figuras da visão com que o dashboard abre: os grupos e a pergunta do estudo, sem filtros. Quando o "
                  "pesquisador escolher outra visão no dashboard, `expression.py figures` com as escolhas dele substitui "
                  "esta seção.", ""]
    if rec is None:
        return lines + ["Nenhuma visão do dashboard foi registrada para o relatório. Com as escolhas do pesquisador, "
                        "`expression.py figures --study-dir D …` desenha os gráficos do dashboard em `figures/dashboard/`; "
                        "depois, `expression.py report` os inclui aqui.", ""]
    master = study_dir / core_master.MASTER_FILE
    if master.is_file() and master.stat().st_mtime > (study_dir / report_figures.CHOICES_FILE).stat().st_mtime:
        lines += ["> **Atenção**: as figuras abaixo são anteriores à última tabela de variantes (`master.parquet`); "
                  "gere-as de novo com o comando abaixo antes de usar este relatório.", ""]
    lines += ["As figuras desta seção foram desenhadas com as mesmas contas e gráficos do dashboard, para as escolhas "
              "abaixo. Para conferir qualquer uma, abra o dashboard (`expression.py dashboard`) e faça as mesmas escolhas.",
              "", *report_figures.describe_choices(rec), "", "**Resultado com essas escolhas**", "", rec["headline"], "",
              "Variantes em cada etapa dos filtros:", ""]
    lines += _table(["Etapa", "Variantes"], [(s, n) for s, n in rec["steps"]])
    n = 0
    for name, (title, caption) in FIGURE_TEXT.items():
        path = rec["figures"].get(name)
        if not path or not (study_dir / path).is_file():
            continue
        n += 1
        lines += ["", f"**Figura D{n}. {title}{"" if title.endswith("?") else "."}** {caption}", "", f"![Figura D{n}: {title}]({path})"]
    lines += ["", f"Comando que gerou as figuras: `{rec['command']}`", ""]
    return lines


def _section_results(study_dir: Path, study: Any, meta: Mapping[str, Any], qc: Mapping[str, Any],
                     df: pl.DataFrame, top: list[dict[str, Any]]) -> list[str]:
    from viz import workbench_helpers

    gates = qc["gates"]
    lines = ["## Resultados da triagem fixa", "",
             f"Tier 1 (`master.parquet`): {qc['rows']:,} SNVs, {qc['n_long_scores']:,} pontuações por trilha; QC "
             f"{'aprovado' if qc['overall_qc_passed'] else 'REPROVADO'} (porta 1: {'ok' if gates[0]['passed'] else 'falhou'}, "
             f"porta 2: {'ok' if gates[1]['passed'] else 'falhou'}).",
             "",
             f"Tier 2 (`tissue_contrast.parquet`): {meta['target_label']} contra {meta['control_label']}. "
             "Direção prevista (todos os tecidos alvo com o mesmo sinal no contraste de direção) e escopo "
             "(votos dos tecidos controle):", ""]
    lines += _table(["Rótulo", "Valor", "SNVs"],
                    [("Direção", k, v) for k, v in sorted(meta["predicted_direction"].items(), key=lambda kv: -kv[1])]
                    + [("Escopo", k, v) for k, v in sorted(meta["specificity_scope"].items(), key=lambda kv: -kv[1])])
    shown = workbench_helpers.default_filter(df)
    lines += ["", f"Com o filtro fixo da triagem (|quantil| >= {workbench_helpers.DEFAULT_MIN_QUANTILE} em todos os "
              f"tecidos alvo) restam {shown.height:,} de {df.height:,} SNVs. Estes rótulos e este filtro são uma triagem "
              "heurística; os grupos, a pergunta e os filtros do estudo são escolhidos pelo pesquisador no dashboard.",
              "", "SNVs por domínio genômico:", ""]
    counts = dict(df.group_by("genomic_domain").len().rows())
    lines += _table(["Domínio", "SNVs"], [(d, counts[d]) for d in workbench_helpers.domain_options(df)])
    diffs = [c["name"] for c in meta["plan"]["contrasts"]]
    lines += ["", f"Principais variantes (top {len(top)}, ordenadas por direção consistente e maior "
              f"|{meta['target_label']} - {meta['control_label']}| em `{meta['plan']['priority']['direction_contrast']}`; "
              "`top_variants.tsv`):", ""]
    headers = ["#", "Variante", "AVI Phred", "Modalidade AVI", "Dist. TSS", "Domínio", "Direção", "Escopo",
               *[f"Δ {n}" for n in diffs]]
    rows = [(t["rank"], f"`{t['variant']}`", t["avi_phred"], t["top_modality"], t["distance_to_tss"], t["genomic_domain"],
             t["predicted_direction"], t["specificity_scope"], *[t[f"{n}_target_vs_control_diff"] for n in diffs])
            for t in top]
    lines += _table(headers, rows)
    fig1 = _figure1_svg(study_dir)
    lines += ["", "Δ = média alvo menos média controle do escore bruto previsto (comparável apenas dentro de um scorer). "
              "A Figura 1 mostra a janela (versão interativa em `figures/fig1_genomic_architecture.html`); os gráficos "
              "interativos estão no dashboard (`expression.py dashboard`).", ""]
    if fig1:
        lines += ["**Figura 1. Arquitetura genômica da janela.** Gene, transcritos, elementos do Ensembl Regulatory Build "
                  "e a janela analisada.", "", f"![Figura 1: Arquitetura genômica da janela]({fig1})", ""]
    links = study_dir / ATLAS_LINKS_FILE
    if links.is_file():
        lines += ["### Links do Atlas", "", _tight_tables(links.read_text(encoding="utf-8").strip()), ""]
    interprets = sorted((study_dir / "raw").glob(INTERPRET_GLOB))
    if interprets:
        lines += ["### Interpretação por variante (AlphaGenome)", ""]
        for path in interprets:
            text = path.read_text(encoding="utf-8", errors="replace").splitlines()
            clipped = len(text) > MAX_INTERPRET_LINES
            lines += [f"`{path.name}`" + (f" (primeiras {MAX_INTERPRET_LINES} de {len(text)} linhas)" if clipped else ""),
                      "", "```", *text[:MAX_INTERPRET_LINES], "```", ""]
    return lines


FIG1_HTML = "figures/fig1_genomic_architecture.html"
FIG1_SVG = "figures/fig1_genomic_architecture.svg"


def _figure1_svg(study_dir: Path) -> str | None:
    """Figure 1 as an image for report.md: the SVG drawn inside the locus stage's HTML, written next to it (and
    written again whenever the HTML is newer). None when there is no Figure 1."""
    html, svg = study_dir / FIG1_HTML, study_dir / FIG1_SVG
    if not html.is_file():
        return None
    if not svg.is_file() or svg.stat().st_mtime < html.stat().st_mtime:
        text = html.read_text(encoding="utf-8")
        a, b = text.find("<svg"), text.find("</svg>")
        if a < 0 or b < 0:
            return None
        body = text[a:b + 6]
        if "xmlns=" not in body[:300]:
            body = body.replace("<svg", '<svg xmlns="http://www.w3.org/2000/svg"', 1)
        svg.write_text(body + "\n", encoding="utf-8")
    return FIG1_SVG


def _tight_tables(text: str) -> str:
    """Drops blank lines between rows of a markdown table (they split it, and it no longer renders as a table)."""
    lines = text.splitlines()
    out = []
    for i, line in enumerate(lines):
        if not line.strip() and out and out[-1].lstrip().startswith("|") and \
                any(nxt.strip() for nxt in lines[i + 1:i + 2]) and lines[i + 1].lstrip().startswith("|"):
            continue
        out.append(line)
    return "\n".join(out)


def _section_skills(manifest: Mapping[str, Any]) -> list[str]:
    lines = ["## Skills utilizadas", "",
             "Lidas de `manifest.json` (cada passo registra a skill e o comando que gerou os arquivos de `raw/`):", ""]
    used = _skills_used(manifest)
    lines += _table(["Skill", "Passos do workflow"],
                    [(f"`{skill}`", ", ".join(f"`{s}`" + (f" (x{steps.count(s)})" if steps.count(s) > 1 else "")
                                              for s in dict.fromkeys(steps))) for skill, steps in used])
    lines += ["", f"Versão da skill `{manifest.get('skill', LOCAL_SKILL)}`: {manifest.get('skill_version', 'desconhecida')}.", ""]
    return lines


def _section_limits(resolution: Mapping[str, Any], meta: Mapping[str, Any]) -> list[str]:
    lines = [
        "## Limitações", "",
        "- **Uso em pesquisa apenas.** São previsões de impacto molecular (AlphaGenome Atlas) para priorizar experimentos; "
        "não são medições nem achados clínicos, e não devem ser usadas como diagnóstico.",
        "- **Escores brutos** são comparáveis dentro de um scorer, não entre scorers; use quantis para comparar ensaios.",
        "- **Escopo da janela:** só a janela analisada foi mutada; enhancers distais e íntrons ficam de fora.",
        "- **Transcrito:** as anotações seguem o transcrito selecionado; outros transcritos que iniciam na janela estão em "
        "`processed/transcripts_in_window.tsv`.",
        "- **Colunas `top_track_*`** do master ficam nulas quando a exportação de AVI não traz colunas de trilha; "
        "a triagem não depende delas.",
    ]
    if resolution["dropped_scorers"]:
        lines.append(f"- Scorers sem trilhas em uma das coortes foram descartados: {', '.join(resolution['dropped_scorers'])}.")
    unserved = [c["curie"] for c in resolution["curies"] if c["status"] == "dropped"]
    if unserved:
        lines.append("- CURIEs sem trilhas próprias no catálogo foram descartados das coortes: "
                     + ", ".join(f"`{c}`" for c in unserved) + ".")
    notes = [n["note"] for n in meta.get("plan_notes", [])]
    if notes:
        lines.append("- Observações do plano de contrastes da triagem (ensaios ou fontes misturados, grupo com um só "
                     "tecido): " + " ".join(notes))
    return lines + [""]


def _section_citations(manifest: Mapping[str, Any]) -> list[str]:
    lines = ["## Citações", "", f"- {PRIMARY_CITATION}"]
    for skill, _ in _skills_used(manifest):
        if skill in SKILL_CITATIONS:
            lines.append(f"- {SKILL_CITATIONS[skill]}")
        elif skill.split("/")[0] == LOCAL_SKILL:
            lines.append(f"- Script local da skill `{skill}` (expression-analysis, versão {manifest.get('skill_version', '?')}).")
        else:
            lines.append(f"- Skill `{skill}` (ver o `SKILL.md` da skill para a citação).")
    return lines + [""]


# --- entry point ----------------------------------------------------------------------------


def render_report(study_dir: Path) -> str:
    """Report text of a study that completed the contrast step.

    Raises:
        ReportError: manifest.json or a processed file is missing.
        study.StudyError, core.* errors: a required study file is invalid.
    """
    study_dir = Path(study_dir)
    study = load_study(study_dir)
    locus = core_locus.load_locus(study_dir)
    resolution = core_cohorts.load_resolution(study_dir)
    qc = core_master.load_qc_report(study_dir)
    manifest = _read_json(study_dir / MANIFEST_FILE, "record the workflow steps with manifest.py add")
    meta = _read_json(study_dir / core_contrast.META_FILE, f"run expression.py contrast (Workflow step {core_contrast.WORKFLOW_STEP})")
    top = _read_json(study_dir / core_contrast.TOP_JSON_FILE, f"run expression.py contrast (Workflow step {core_contrast.WORKFLOW_STEP})")
    contrast_path = study_dir / core_contrast.CONTRAST_FILE
    if not contrast_path.is_file():
        raise ReportError(f"{core_contrast.CONTRAST_FILE} not found in {study_dir}; run expression.py contrast "
                          f"(Workflow step {core_contrast.WORKFLOW_STEP}).")
    df = pl.read_parquet(contrast_path)
    lines = [f"# Relatório: {locus.gene_symbol} ({study.study_id})", "",
             "> Uso em pesquisa apenas. Predições de impacto molecular, não diagnóstico.", ""]
    lines += _section_study(study, locus, resolution, manifest)
    lines += _section_decisions(study_dir, meta)
    lines += _section_dashboard(study_dir)
    lines += _section_results(study_dir, study, meta, qc, df, top)
    lines += _section_skills(manifest)
    lines += _section_limits(resolution, meta)
    lines += _section_citations(manifest)
    return "\n".join(lines)


def write_report(study_dir: Path) -> Path:
    """Writes <study_dir>/report.md and returns its path."""
    study_dir = Path(study_dir)
    path = study_dir / REPORT_FILE
    path.write_text(render_report(study_dir), encoding="utf-8")
    return path
