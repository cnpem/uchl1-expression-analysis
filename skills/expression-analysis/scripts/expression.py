#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["pyyaml", "matplotlib", "polars", "pyarrow", "pandas", "altair"]
# ///
"""CLI for the expression-analysis skill.

Subcommands:
  check      Verify prerequisites (uv, dependency skills, licence files, API keys).
             Prints what is missing and the install command; never installs anything.
  locus      Workflow step 4: transcript, TSS, exon 1 and regulatory window from raw/gtf.json,
             raw/regulatory.json and raw/sequence_*.txt. Writes processed/locus.json,
             processed/transcripts_in_window.tsv and figures/fig1_genomic_architecture.html;
             open questions go to processed/pending_decisions.json (exit 1).
  cohorts    Workflow step 8: UBERON/CL cohorts matched with the Atlas track catalogue. Writes
             processed/track_registry.tsv, processed/ontology_resolution.json and
             processed/tissue_cohorts.md; open questions (cohort.unserved:*,
             cohort.scorer_without_cohort:*) go to processed/pending_decisions.json (exit 1).
  variants   Workflow step 11: saturated SNV table of the window (REF from the Ensembl sequence,
             cross-checked with the AVI export). Writes processed/core_variants.tsv; REF
             disagreements and SNVs without an AVI row (variants.*) go to
             processed/pending_decisions.json (exit 1).
  master     Workflow step 11: Tier 1 table (variants + AVI + per-track raw/quantile) with the
             2-gate QC. Writes processed/master.parquet, track_scores_long.parquet, qc_report.json;
             exit 1 when a gate fails.
  contrast   Workflow step 12: Tier 2 target-vs-control contrasts. --propose writes
             processed/contrasts_proposed.yaml and lists the review points. Without --propose it
             needs decision contrast.plan and contrasts.yaml; writes processed/tissue_contrast.parquet,
             tissue_contrast_meta.json, top_variants.tsv and top_variants.json; open review points
             (contrast.review:*) go to processed/pending_decisions.json (exit 1).
  report     Workflow step 13: writes <study-dir>/report.md (Estudo, Decisoes, Resultados, Skills
             utilizadas, Limitacoes, Citacoes) from processed/*, decisions.yaml and manifest.json.
  dashboard  Workflow step 14: opens the marimo triage workbench (assets/workbench.py) on a study
             directory; --export-html writes a static HTML with Figures 1-4 instead.
  reproduce  Re-run the registered stages on a bundled reference study and compare with
             its expected.json.

Workflow helper (no analysis; it spares the agent hand arithmetic):
  regions    Regions the dependency commands need: the regulatory search region of step 3 (gene
             span of raw/gtf.json over all transcripts, starts normalised to 1-based, +/-
             window.regulatory_search_upstream_bp on both sides, so either strand is covered) and,
             once processed/locus.json exists, the analysis window of steps 9-10 and max_window.
             JSON on stdout; --only NAME prints one value (chr:start-end, or the integer max_window).
"""

import argparse
import json
import os
import re
import shutil
import sys
from pathlib import Path
from typing import Callable

import models
import reference
import study
from core import cohorts as core_cohorts
from core import contrast as core_contrast
from core import locus as core_locus
from core import master as core_master
from core import report as core_report
from core import variants as core_variants
from viz import genomic_architecture
from viz import workbench_helpers

SKILL_DIR = Path(__file__).resolve().parent.parent


def _stage_error(exc: Exception) -> None:
    print(f"Error: {exc}", file=sys.stderr)
    sys.exit(1)


def run_locus(study_dir: Path) -> None:
    """Workflow step 4: builds the locus and window, writes locus.json, the transcript table and Figure 1."""
    try:
        locus = core_locus.run(study_dir)
        figure = genomic_architecture.write_figure(study_dir, locus)
    except core_locus.LocusPending as exc:
        study.raise_pending(study_dir, exc.items)
    except (core_locus.LocusError, models.DependencyOutputError, study.StudyError) as exc:
        _stage_error(exc)
    codon = ("start codon reads ATG" if locus.start_codon_verified
             else "start codon not covered by a raw sequence, not checked")
    print(f"Locus: {locus.transcript_id} ({'+' if locus.strand == 1 else '-'} strand), TSS {locus.chromosome}:{locus.tss}, "
          f"window {locus.window_str} ({locus.window.length} bp, {len(locus.features)} regulatory feature(s)); "
          f"{codon}; window sequence from {locus.sequence_file}. Wrote {core_locus.LOCUS_FILE}, {core_locus.TRANSCRIPTS_FILE}, "
          f"{figure.relative_to(study_dir)}")


def locus_cmd(args: argparse.Namespace) -> None:
    run_locus(Path(args.study_dir))


def run_cohorts(study_dir: Path) -> None:
    """Workflow step 8: resolves the cohorts against the track catalogue and writes the track registry."""
    try:
        resolution = core_cohorts.run(study_dir)
    except core_cohorts.CohortPending as exc:
        study.raise_pending(study_dir, exc.items)
    except (core_cohorts.CohortError, core_locus.LocusError, models.DependencyOutputError, study.StudyError) as exc:
        _stage_error(exc)
    metrics = core_cohorts.stage_metrics(study_dir)
    counts = ", ".join(f"{s} {metrics[f'n_tracks.{s}']}" for s in resolution["effective_scorers"])
    dropped = resolution["dropped_scorers"]
    print(f"Cohorts: target {len(resolution['effective_target_curies'])} CURIE(s), control "
          f"{len(resolution['effective_control_curies'])} CURIE(s); tracks {counts}"
          + (f"; dropped scorer(s) {', '.join(dropped)}" if dropped else "")
          + f". Wrote {core_cohorts.REGISTRY_FILE}, {core_cohorts.RESOLUTION_FILE}, {core_cohorts.TISSUE_DOC_FILE}")


def cohorts_cmd(args: argparse.Namespace) -> None:
    run_cohorts(Path(args.study_dir))


def run_variants(study_dir: Path) -> None:
    """Workflow step 11: enumerates the SNVs of the window and writes core_variants.tsv."""
    try:
        df = core_variants.run(study_dir)
    except core_variants.VariantsPending as exc:
        study.raise_pending(study_dir, exc.items)
    except (core_variants.VariantsError, core_locus.LocusError, models.DependencyOutputError, study.StudyError) as exc:
        _stage_error(exc)
    print(f"Variants: {df.height} SNV(s) at {df['position'].n_unique()} position(s). Wrote {core_variants.CORE_VARIANTS_FILE}")


def variants_cmd(args: argparse.Namespace) -> None:
    run_variants(Path(args.study_dir))


def run_master(study_dir: Path) -> None:
    """Workflow step 11: Tier 1 master table and the 2-gate QC; exits 1 when a gate fails."""
    try:
        report = core_master.run(study_dir)
    except (core_master.MasterError, core_variants.VariantsError, core_locus.LocusError,
            models.DependencyOutputError, study.StudyError) as exc:
        _stage_error(exc)
    status = "passed" if report.passed else "FAILED"
    print(f"Master: {report.total_variants} row(s), QC {status} (gate 1 {report.gate1_pass}, gate 2 {report.gate2_pass}; "
          f"ref mismatches {report.ref_mismatches}, missing AVI {report.missing_avi}, top-feature mismatches "
          f"{report.top_feature_mismatches}, indel attributions {report.indel_attributions}). "
          f"Wrote {core_master.MASTER_FILE}, {core_master.LONG_FILE}, {core_master.QC_FILE}")
    if not report.passed:
        print(f"Error: QC failed for {len(report.failed_variants)} variant(s); see {core_master.QC_FILE} (failed_variants).",
              file=sys.stderr)
        sys.exit(1)


def master_cmd(args: argparse.Namespace) -> None:
    run_master(Path(args.study_dir))


def run_contrast(study_dir: Path) -> None:
    """Workflow step 12: executes the approved contrast plan (needs decision contrast.plan)."""
    try:
        df = core_contrast.run(study_dir)
    except core_contrast.ContrastPending as exc:
        study.raise_pending(study_dir, exc.items)
    except (core_contrast.ContrastError, core_master.MasterError, core_cohorts.CohortError, models.DependencyOutputError,
            study.StudyError) as exc:
        _stage_error(exc)
    print(f"Contrast: {df.height} row(s). Wrote {core_contrast.CONTRAST_FILE}, {core_contrast.META_FILE}, "
          f"{core_contrast.TOP_TSV_FILE}, {core_contrast.TOP_JSON_FILE}")


def contrast_cmd(args: argparse.Namespace) -> None:
    study_dir = Path(args.study_dir)
    if not args.propose:
        run_contrast(study_dir)
        return
    try:
        path, plan, points = core_contrast.propose(study_dir)
    except (core_contrast.ContrastError, core_master.MasterError, core_cohorts.CohortError, models.DependencyOutputError,
            study.StudyError) as exc:
        _stage_error(exc)
    print(f"Proposed plan: {', '.join(c.name for c in plan.contrasts)}; priority from {plan.priority.direction_contrast}. "
          f"Wrote {path.relative_to(study_dir)}")
    for q in points:
        print(f"  review: {q.key}: {q.question}")


def run_workbench_metrics(study_dir: Path) -> None:
    """Workflow step 14 (headless): loads the study like the workbench and builds Figures 2-4 for its default state."""
    try:
        data = workbench_helpers.load_workbench_datasets(study_dir)
        built = workbench_helpers.build_default_figures(data)
    except (core_contrast.ContrastError, core_master.MasterError, core_cohorts.CohortError, core_locus.LocusError,
            models.DependencyOutputError, study.StudyError, ValueError) as exc:
        _stage_error(exc)
    print(f"Workbench: {built['shown'].height} of {data.df.height} variant(s) under the default filters; "
          f"Figure 3 {built['fig3_meta']['n_variants']} row(s) x {len(built['fig3_meta']['tissues'])} tissue(s).")


def report_cmd(args: argparse.Namespace) -> None:
    study_dir = Path(args.study_dir)
    try:
        path = core_report.write_report(study_dir)
    except (core_report.ReportError, core_contrast.ContrastError, core_master.MasterError, core_cohorts.CohortError,
            core_locus.LocusError, models.DependencyOutputError, study.StudyError) as exc:
        _stage_error(exc)
    print(f"Report: wrote {path}")


WORKBENCH_APP = SKILL_DIR / "assets" / "workbench.py"
WORKBENCH_DEPS = ["marimo>=0.25", "altair>=5.4", "anywidget", "pandas", "polars", "pyarrow", "pyyaml", "matplotlib"]


def dashboard_cmd(args: argparse.Namespace) -> None:
    """Opens the marimo workbench, or exports it to static HTML (marimo runs in an ephemeral uv environment)."""
    import subprocess

    study_dir = Path(args.study_dir).resolve()
    try:
        workbench_helpers.load_workbench_datasets(study_dir)
    except (core_contrast.ContrastError, core_master.MasterError, core_cohorts.CohortError, core_locus.LocusError,
            study.StudyError) as exc:
        _stage_error(exc)
    if shutil.which("uv") is None:
        _stage_error(RuntimeError("'uv' is not on PATH; install it with the 'uv' skill."))
    # marimo writes __marimo__/session/ next to the notebook file. Launch it through a symlink inside
    # the study directory so that folder lands in the study, never inside the skill's assets/.
    # (The app resolves the symlink to find scripts/.)
    app_dir = study_dir / ".workbench"
    app_dir.mkdir(exist_ok=True)
    app = app_dir / WORKBENCH_APP.name
    if app.is_symlink() or app.exists():
        app.unlink()
    app.symlink_to(WORKBENCH_APP)
    base = ["uv", "run", "--no-project", "--python", "3.13", *[a for d in WORKBENCH_DEPS for a in ("--with", d)], "marimo"]
    if args.export_html:
        out = Path(args.export_html).resolve()
        out.parent.mkdir(parents=True, exist_ok=True)
        cmd = [*base, "export", "html", str(app), "-o", str(out), "--force", "--", "--study-dir", str(study_dir)]
    else:
        cmd = [*base, "run", str(app), "--port", str(args.port), "--", "--study-dir", str(study_dir)]
    sys.exit(subprocess.call(cmd, cwd=study_dir))


# --------------------------------------------------------------------------- workflow helper

REGION_NAMES = ("regulatory_search", "window", "max_window")


def _region_entry(chromosome: str, start: int, end: int, rule: str, used_by: str) -> dict:
    return {"region": models.interval_str(chromosome, start, end), "chromosome": chromosome, "start": start,
            "end": end, "length": end - start + 1, "rule": rule, "used_by": used_by}


def study_regions(study_dir: Path) -> dict:
    """Regions of the Workflow's dependency commands, from study.yaml, raw/gtf.json and processed/locus.json.

    regulatory_search: gene span over every transcript of raw/gtf.json (after models.normalize_gtf, so
    1-based) +/- window.regulatory_search_upstream_bp on both sides; symmetric, so it holds the
    upstream flank on either strand. window and max_window: from processed/locus.json, or None
    before the locus step.

    Raises:
        core_locus.LocusError, models.DependencyOutputError, study.StudyError, ValueError.
    """
    study_dir = Path(study_dir)
    config = study.load_study(study_dir)
    gtf = study_dir / core_locus.GTF_FILE
    if not gtf.is_file():
        raise core_locus.LocusError(f"{core_locus.GTF_FILE} not found in {study_dir}; run Workflow step 2 "
                                    "(alphagenome-variant-impact-score gtf).")
    transcripts = [models.normalize_gtf(t) for t in models.load_raw(gtf)]
    chromosomes = sorted({t.chromosome for t in transcripts})
    if len(chromosomes) != 1:
        raise core_locus.LocusError(f"{core_locus.GTF_FILE}: expected the transcripts of one gene on one contig, "
                                    f"found {len(transcripts)} transcript(s) on {chromosomes or 'none'}")
    upstream = int((config.data.get("window") or {}).get("regulatory_search_upstream_bp") or 0)
    start = max(1, min(t.start for t in transcripts) - upstream)
    end = max(t.end for t in transcripts) + upstream
    out = {
        "study_id": config.study_id,
        "gene_symbol": config.data.get("gene_symbol"),
        "regulatory_search": _region_entry(
            chromosomes[0], start, end,
            f"gene span of {core_locus.GTF_FILE} ({len(transcripts)} transcripts, 1-based) +/- "
            f"window.regulatory_search_upstream_bp ({upstream}) on both sides",
            "step 3 regulatory_features.py overlap --region and ensembl-database get-sequence"),
        "window": None,
        "max_window": None,
    }
    if (study_dir / core_locus.LOCUS_FILE).is_file():
        locus = core_locus.load_locus(study_dir)
        out["window"] = _region_entry(locus.chromosome, locus.window.start, locus.window.end,
                                      f"{core_locus.LOCUS_FILE} window (policy {locus.window_policy})",
                                      "step 9 region --region and step 10 atlas_query.py scores --interval")
        out["max_window"] = study.max_window(study_dir)
    return out


def regions_cmd(args: argparse.Namespace) -> None:
    try:
        regions = study_regions(Path(args.study_dir))
    except (core_locus.LocusError, models.DependencyOutputError, study.StudyError, ValueError) as exc:
        _stage_error(exc)
    if args.only and regions[args.only] is None:
        _stage_error(core_locus.LocusError(f"{args.only} needs {core_locus.LOCUS_FILE}; run expression.py locus "
                                           f"(Workflow step {core_locus.WORKFLOW_STEP}) first."))
    if args.only == "max_window":
        print(regions["max_window"])
    elif args.only:
        print(regions[args.only]["region"])
    else:
        print(json.dumps(regions, indent=2))


# Stage registry for `reproduce`. Each task that implements a stage appends
# (name, run_fn, metrics_fn): run_fn(study_dir) runs it, metrics_fn(study_dir) returns its
# {metric: value} dict for comparison with expected.json. Order = execution order.
STAGES: list[tuple[str, Callable[[Path], None], Callable[[Path], dict]]] = [
    ("locus", run_locus, core_locus.stage_metrics),
    ("cohorts", run_cohorts, core_cohorts.stage_metrics),
    ("variants", run_variants, core_variants.stage_metrics),
    ("master", run_master, core_master.stage_metrics),
    ("contrast", run_contrast, core_contrast.stage_metrics),
    ("workbench", run_workbench_metrics, workbench_helpers.stage_metrics),
]

# name -> (origin, install command). Install commands are only printed, never run.
_DEEPMIND = "google-deepmind/science-skills"
_KDENSE = "K-Dense-AI/scientific-agent-skills"
DEPENDENCIES = {
    "uv": ("DeepMind", f"npx skills add {_DEEPMIND} --skill uv"),
    "credentials": ("DeepMind", f"npx skills add {_DEEPMIND} --skill credentials"),
    "alphagenome-variant-impact-score": ("DeepMind", f"npx skills add {_DEEPMIND} --skill alphagenome-variant-impact-score"),
    "alphagenome-single-variant-analysis": ("DeepMind", f"npx skills add {_DEEPMIND} --skill alphagenome-single-variant-analysis"),
    "alphagenome-atlas-website-links": ("DeepMind", f"npx skills add {_DEEPMIND} --skill alphagenome-atlas-website-links"),
    "embl-ebi-ols": ("DeepMind", f"npx skills add {_DEEPMIND} --skill embl-ebi-ols"),
    "ensembl-database": ("DeepMind", f"npx skills add {_DEEPMIND} --skill ensembl-database"),
    "alphagenome": ("K-Dense", f"npx skills add {_KDENSE} --skill alphagenome"),
}
REQUIRED_KEYS = ["ALPHAGENOME_API_KEY"]
OPTIONAL_KEYS = []
_LICENSE_RE = re.compile(r"\.licenses/([A-Za-z0-9_.-]+_LICENSE\.txt)")


def _search_dirs() -> list[Path]:
    """Directories searched for installed skills, in priority order."""
    cwd, home = Path.cwd(), Path.home()
    return [
        SKILL_DIR.parent,
        cwd / ".claude" / "skills",
        cwd / ".agents" / "skills",
        home / ".claude" / "skills",
        home / ".agents" / "skills",
    ]


def _find_skill(name: str) -> Path | None:
    for d in _search_dirs():
        candidate = d / name / "SKILL.md"
        if candidate.is_file():
            return candidate
    return None


def _key_present(name: str) -> bool:
    """True if the variable is set; the value is never stored or printed."""
    if os.environ.get(name):
        return True
    env_file = Path.home() / ".env"
    try:
        for line in env_file.read_text().splitlines():
            line = line.strip()
            if line.startswith("export "):
                line = line[len("export "):]
            key, sep, value = line.partition("=")
            if sep and key.strip() == name and value.strip().strip("\"'"):
                return True
    except OSError:
        pass
    return False


def check(args: argparse.Namespace) -> None:
    errors: list[str] = []
    warnings: list[str] = []
    installed: dict[str, Path] = {}

    if shutil.which("uv") is None:
        errors.append(
            "Missing prerequisite: 'uv' is not on PATH. "
            "Install it with the 'uv' skill (see below) or from https://docs.astral.sh/uv/."
        )
    for name, (origin, command) in DEPENDENCIES.items():
        found = _find_skill(name)
        if found is None:
            errors.append(
                f"Missing dependency skill '{name}' (origin: {origin}).\n"
                f"  Install command (not executed): {command}"
            )
        else:
            installed[name] = found

    licenses: dict[str, str] = {}
    for name, path in installed.items():
        for match in _LICENSE_RE.findall(path.read_text(errors="replace")):
            licenses.setdefault(match, name)
    for filename, owner in sorted(licenses.items()):
        if not (Path.cwd() / ".licenses" / filename).is_file():
            errors.append(
                f"Missing licence file '.licenses/{filename}' required by skill '{owner}'. "
                f"Follow the Prerequisites section of '{owner}' with the researcher."
            )

    for key in REQUIRED_KEYS:
        if not _key_present(key):
            errors.append(
                f"Missing API key {key}. Set it with the 'credentials' skill "
                "(stored in ~/.env) or export it in the environment."
            )
    for key in OPTIONAL_KEYS:
        if not _key_present(key):
            warnings.append(f"Warning: {key} is not set; requests will run at the lower rate limit.")

    for w in warnings:
        print(w, file=sys.stderr)
    if errors:
        print(f"{len(errors)} prerequisite(s) missing:", file=sys.stderr)
        for e in errors:
            print(f"- {e}", file=sys.stderr)
        sys.exit(1)
    print("All prerequisites satisfied")


def reproduce(args: argparse.Namespace) -> None:
    reference_dir = SKILL_DIR / "assets" / "reference" / args.reference
    if not reference_dir.is_dir():
        available = sorted(p.name for p in (SKILL_DIR / "assets" / "reference").glob("*") if p.is_dir())
        print(f"Unknown reference '{args.reference}'. Available: {', '.join(available) or 'none'}.", file=sys.stderr)
        sys.exit(1)
    keep = Path(args.keep) if args.keep else None
    # Capture the per-scorer track counts while the cohorts stage runs (the temporary study is
    # removed afterwards) so the summary line can state the total number of tracks.
    track_counts: dict[str, int] = {}

    def _cohorts_metrics(study_dir: Path) -> dict:
        metrics = core_cohorts.stage_metrics(study_dir)
        track_counts.update({k.split(".", 1)[1]: v for k, v in metrics.items() if k.startswith("n_tracks.")})
        return metrics

    stages = [(n, r, _cohorts_metrics if n == "cohorts" else m) for n, r, m in STAGES]
    code = reference.run_reference(reference_dir, stages, args.until, keep)
    if code == 0 and track_counts:
        parts = ", ".join(f"{s} {n}" for s, n in track_counts.items())
        print(f"Summary: {sum(track_counts.values())} tracks in total ({parts}). "
              "Report these numbers, and the headline results above, in the answer.")
    sys.exit(code)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("check", help="Verify prerequisites; print install commands for what is missing")
    p.set_defaults(func=check)

    p = sub.add_parser("locus", help="Workflow step 4: locus, TSS and regulatory window")
    p.add_argument("--study-dir", required=True, help="Study directory")
    p.set_defaults(func=locus_cmd)

    p = sub.add_parser("cohorts", help="Workflow step 8: UBERON/CL cohorts and Atlas track registry")
    p.add_argument("--study-dir", required=True, help="Study directory")
    p.set_defaults(func=cohorts_cmd)

    p = sub.add_parser("variants", help="Workflow step 11: saturated SNV table")
    p.add_argument("--study-dir", required=True, help="Study directory")
    p.set_defaults(func=variants_cmd)

    p = sub.add_parser("master", help="Workflow step 11: Tier 1 master table and QC gates")
    p.add_argument("--study-dir", required=True, help="Study directory")
    p.set_defaults(func=master_cmd)

    p = sub.add_parser("contrast", help="Workflow step 12: Tier 2 contrasts (--propose writes the proposal)")
    p.add_argument("--study-dir", required=True, help="Study directory")
    p.add_argument("--propose", action="store_true", help="Write processed/contrasts_proposed.yaml and stop")
    p.set_defaults(func=contrast_cmd)

    p = sub.add_parser("report", help="Workflow step 13: research report (report.md)")
    p.add_argument("--study-dir", required=True, help="Study directory")
    p.set_defaults(func=report_cmd)

    p = sub.add_parser("dashboard", help="Workflow step 14: marimo triage workbench")
    p.add_argument("--study-dir", required=True, help="Study directory")
    p.add_argument("--port", type=int, default=2718, help="Port of the marimo server (default 2718)")
    p.add_argument("--export-html", help="Write a static HTML export to this path instead of serving")
    p.set_defaults(func=dashboard_cmd)

    p = sub.add_parser("reproduce", help="Re-run stages on a bundled reference study")
    p.add_argument("--reference", required=True, help="Reference study name")
    p.add_argument("--until", help="Last stage to run")
    p.add_argument("--keep", help="Keep the working copy in this directory")
    p.set_defaults(func=reproduce)

    p = sub.add_parser("regions", help="Helper: regulatory search region, window and max_window (JSON)",
                       description="Prints the regions the Workflow's dependency commands need, so no coordinate "
                                   "is worked out by hand: regulatory_search (step 3: gene span of raw/gtf.json "
                                   "over all transcripts, starts normalised to 1-based, +/- "
                                   "window.regulatory_search_upstream_bp on both sides) and, after the locus "
                                   "step, window (steps 9-10, from processed/locus.json) and max_window "
                                   "(= study.py window --print-max-window).")
    p.add_argument("--study-dir", required=True, help="Study directory")
    p.add_argument("--only", choices=REGION_NAMES,
                   help="Print only this value: chr:start-end for a region, the integer for max_window")
    p.set_defaults(func=regions_cmd)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
