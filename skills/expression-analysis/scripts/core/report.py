"""Research report of a study (Workflow step 13, `expression.py report`).

Writes <study-dir>/report.md in Portuguese with six sections: Estudo, Decisões, Resultados,
Skills utilizadas, Limitações and Citações. Everything comes from files of the study directory:
study.yaml, decisions.yaml, manifest.json, processed/* and, when present, raw/atlas_links.md and
raw/interpret_*.txt (both optional: the report includes them only if the researcher produced them).
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
from study import load_decisions, load_study

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
    lines += ["Respostas do pesquisador registradas em `decisions.yaml` (passo do workflow, chave, opções oferecidas, "
              "recomendada, valor escolhido); \"—\" quando a pergunta foi registrada sem opções:", ""]
    rows = [(e.workflow_step, f"`{key}`", "; ".join(f"`{o}`" for o in e.options) if e.options else "—",
             f"`{e.recommended}`" if e.recommended is not None else "—",
             f"`{e.value}`" if not isinstance(e.value, list) else ", ".join(map(str, e.value)),
             e.answered_at[:10]) for key, e in sorted(decisions.items(), key=lambda kv: (kv[1].workflow_step, kv[0]))]
    lines += _table(["Passo", "Decisão", "Opções", "Recomendada", "Valor", "Data"], rows)
    plan = meta["plan"]
    lines += ["", "Plano de contrastes aprovado (`contrasts.yaml`):", ""]
    rows = [(c["name"], c["scorer"], c["strand"], c["tissue_aggregation"], len(c["target"]), len(c["control"]))
            for c in plan["contrasts"]]
    lines += _table(["Contraste", "Scorer", "Fita", "Agregação por tecido", "Tecidos alvo", "Tecidos controle"], rows)
    lines += ["", f"Rótulos de prioridade a partir do contraste `{plan['priority']['direction_contrast']}`.", ""]
    return lines


def _section_results(study_dir: Path, study: Any, meta: Mapping[str, Any], qc: Mapping[str, Any],
                     df: pl.DataFrame, top: list[dict[str, Any]]) -> list[str]:
    from viz import workbench_helpers

    gates = qc["gates"]
    lines = ["## Resultados", "",
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
    lines += ["", f"Com os filtros padrão do workbench (|quantil| >= {workbench_helpers.DEFAULT_MIN_QUANTILE} em todos os "
              f"tecidos alvo) restam {shown.height:,} de {df.height:,} SNVs.", "", "SNVs por domínio genômico:", ""]
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
    lines += ["", "Δ = média alvo menos média controle do escore bruto previsto (comparável apenas dentro de um scorer). "
              "A Figura 1 (`figures/fig1_genomic_architecture.html`) mostra a janela; as Figuras 2 a 4 estão no workbench "
              "(`expression.py dashboard`).", ""]
    links = study_dir / ATLAS_LINKS_FILE
    if links.is_file():
        lines += ["### Links do Atlas", "", links.read_text(encoding="utf-8").strip(), ""]
    interprets = sorted((study_dir / "raw").glob(INTERPRET_GLOB))
    if interprets:
        lines += ["### Interpretação por variante (AlphaGenome)", ""]
        for path in interprets:
            text = path.read_text(encoding="utf-8", errors="replace").splitlines()
            clipped = len(text) > MAX_INTERPRET_LINES
            lines += [f"`{path.name}`" + (f" (primeiras {MAX_INTERPRET_LINES} de {len(text)} linhas)" if clipped else ""),
                      "", "```", *text[:MAX_INTERPRET_LINES], "```", ""]
    return lines


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
    reviews = [d["key"] for d in meta.get("review_decisions", []) if d["key"].startswith("contrast.review:")]
    if reviews:
        lines.append("- Pontos de revisão do plano de contrastes aceitos pelo pesquisador (fontes ou ensaios misturados, "
                     "tecido único): " + ", ".join(f"`{k.removeprefix('contrast.review:')}`" for k in reviews) + ".")
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
