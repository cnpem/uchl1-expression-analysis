#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["pyyaml==6.0.3", "matplotlib==3.11.2", "polars==1.31.0", "pyarrow==26.0.1", "pandas==3.1.1", "altair==6.3.0", "vl-convert-python==1.9.0.post1"]
# ///
"""CLI for the expression-analysis skill.

Subcommands:
  check      Verify prerequisites (uv, dependency skills, licence files, API keys).
             Prints what is missing and the install command; never installs anything. --json
             lists them as questions for the setup.
  locus      Workflow step 4: transcript, TSS, exon 1 and regulatory window from raw/gtf.json,
             raw/regulatory.json and raw/sequence_*.txt. Writes processed/locus.json,
             processed/transcripts_in_window.tsv and figures/fig1_genomic_architecture.html;
             open questions (locus.transcript, locus.feature_*, locus.no_feature_chain) go to
             processed/pending_decisions.json (exit 1); the rest of the window follows its policy
             (rules recorded in decisions.yaml).
  cohorts    Workflow step 8: UBERON/CL cohorts matched with the Atlas track catalogue. Writes
             processed/track_registry.tsv, processed/ontology_resolution.json and
             processed/tissue_cohorts.md; open questions (cohort.unserved:*, cohort.no_rna_seq:*)
             go to processed/pending_decisions.json (exit 1); a prediction type without tracks in
             one group is left out by rule.
  variants   Workflow step 11: saturated SNV table of the window (REF from the Ensembl sequence,
             cross-checked with the AVI export). Writes processed/core_variants.tsv; positions
             with REF disagreements or SNVs without an AVI row are left out by rule, unless they
             exceed 1% of the window (question variants.exclusion_limit, exit 1).
  master     Workflow step 11: Tier 1 table (variants + AVI + per-track raw/quantile) with the
             2-gate QC. Writes processed/master.parquet, track_scores_long.parquet, qc_report.json;
             exit 1 when a gate fails.
  contrast   Workflow step 12: Tier 2 target-vs-control contrasts (the report's fixed triage).
             Generates contrasts.yaml from the cohorts when absent (rule contrast.plan) and writes
             processed/tissue_contrast.parquet, tissue_contrast_meta.json, top_variants.tsv and
             top_variants.json; no question. --propose only writes processed/contrasts_proposed.yaml
             for a maintainer.
  report     Workflow step 13: writes <study-dir>/report.md (Estudo, Decisoes, Analise do pesquisador,
             Resultados da triagem, Skills utilizadas, Limitacoes, Citacoes) from processed/*,
             figures/dashboard/, decisions.yaml and manifest.json.
  summarize  Text summary of the dashboard view for given choices (same computation; for the agent to
             explain what the researcher sees). Writes nothing.
  figures    Workflow step 14: the dashboard's charts for given choices as PNG in figures/dashboard/,
             with choices.json; report.md then shows them with the choices behind them.
  dashboard  Workflow step 14: opens the analysis dashboard (assets/dashboard.py) on a study directory;
             --export-html writes a static HTML of its opening view instead. --edit (maintainers) opens
             the marimo editor on the skill's file; --legacy opens the previous triage workbench.
  reproduce  Re-run the registered stages on a bundled reference study and compare with
             its expected.json.

Captures and exploration (each runs a dependency skill with pinned packages and records itself in
manifest.json; --dry-run prints the commands):
  explore    Workflow step 0: read-only exploration of a gene before the study exists -> OUT/exploration/
             (summary.md for the researcher, exploration.json, raw/, atlas/, manifest.json).
  fetch      fetch gtf | regulatory | sequence | ols | ols-term | tracks | avi | scores | links | interpret
             --study-dir D: one dependency capture of Workflow steps 2-13, written under the name the next
             stage reads.

Every stage subcommand above appends its own command line, outputs and the decision keys it used to
manifest.json when it succeeds (report and dashboard included).

Workflow helper (no analysis; it spares the agent hand arithmetic):
  regions    Regions the dependency commands need: the regulatory search region of step 3 (gene
             span of raw/gtf.json over all transcripts, starts normalised to 1-based, +/-
             window.regulatory_search_upstream_bp on both sides, so either strand is covered) and,
             once processed/locus.json exists, the analysis window of steps 9-10 and max_window.
             JSON on stdout; --only NAME prints one value (chr:start-end, or the integer max_window).
"""

import argparse
import datetime
import json
import os
import shutil
import sys
from pathlib import Path
from typing import Any, Callable

import deps
import explore as core_explore
import fetch
import manifest
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


def _record(study_dir: Path, step: str, outputs: list[str], areas: tuple[str, ...] = (),
            versions: dict[str, str] | None = None) -> None:
    """Appends the stage that just ran (this exact command line) to the study manifest.

    The decision keys are the recorded answers of the stage's areas (e.g. `locus.`), so the manifest says
    which researcher decisions the outputs rest on. Called by the CLI only: `reproduce` runs the same
    stages on a copy of the reference and records nothing.
    """
    # The researcher's answers under the stage's areas, plus the skill's rules other than per-position exclusions
    # (those can be hundreds; decisions.yaml and the report list each one)
    decisions = study.load_decisions(study_dir) if areas else {}
    keys = [k for k, e in decisions.items() if k.startswith(areas)
            and not (e.decided_by == "skill" and k.startswith(core_variants.RULE_PREFIXES))]
    present = [o for o in outputs if (study_dir / o).is_file()]
    try:
        manifest.add_step(study_dir, manifest.ManifestStep(
            step=step, skill="expression-analysis", command=fetch.invocation(), outputs=present,
            source_versions={"skill_version": manifest._skill_identity()[1], **(versions or {})},
            recorded_at=manifest._now(), decision_keys=keys))
    except manifest.ManifestError as exc:
        print(f"Warning: not recorded in manifest.json ({exc})", file=sys.stderr)
        return
    print(f"Recorded step {step!r} in {study_dir / manifest.MANIFEST_NAME}")


def run_locus(study_dir: Path) -> None:
    """Workflow step 4: builds the locus and window, writes locus.json, the transcript table and Figure 1."""
    try:
        locus = core_locus.run(study_dir)
        figure = genomic_architecture.write_figure(study_dir, locus)
        study.clear_pending(study_dir, ("locus.",))
    except core_locus.LocusPending as exc:
        study.raise_pending(study_dir, exc.items, ("locus.",))
    except (core_locus.LocusError, models.DependencyOutputError, study.StudyError) as exc:
        _stage_error(exc)
    codon = ("start codon reads ATG" if locus.start_codon_verified
             else "start codon not covered by a raw sequence, not checked")
    print(f"Locus: {locus.transcript_id} ({'+' if locus.strand == 1 else '-'} strand), TSS {locus.chromosome}:{locus.tss}, "
          f"window {locus.window_str} ({locus.window.length} bp, {len(locus.features)} regulatory feature(s)); "
          f"{codon}; window sequence from {locus.sequence_file}. Wrote {core_locus.LOCUS_FILE}, {core_locus.TRANSCRIPTS_FILE}, "
          f"{figure.relative_to(study_dir)}")


def locus_cmd(args: argparse.Namespace) -> None:
    study_dir = Path(args.study_dir)
    run_locus(study_dir)
    _record(study_dir, "4-locus", [core_locus.LOCUS_FILE, core_locus.TRANSCRIPTS_FILE,
                                   "figures/fig1_genomic_architecture.html"], ("locus.",))


def run_cohorts(study_dir: Path) -> None:
    """Workflow step 8: resolves the cohorts against the track catalogue and writes the track registry."""
    try:
        resolution = core_cohorts.run(study_dir)
        study.clear_pending(study_dir, ("cohort.",))
    except core_cohorts.CohortPending as exc:
        study.raise_pending(study_dir, exc.items, ("cohort.",))
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
    study_dir = Path(args.study_dir)
    run_cohorts(study_dir)
    _record(study_dir, "8-cohorts", [core_cohorts.REGISTRY_FILE, core_cohorts.RESOLUTION_FILE,
                                     core_cohorts.TISSUE_DOC_FILE], ("cohort.",))


def run_variants(study_dir: Path) -> None:
    """Workflow step 11: enumerates the SNVs of the window and writes core_variants.tsv."""
    try:
        df = core_variants.run(study_dir)
        study.clear_pending(study_dir, ("variants.",))
    except core_variants.VariantsPending as exc:
        study.raise_pending(study_dir, exc.items, ("variants.",))
    except (core_variants.VariantsError, core_locus.LocusError, models.DependencyOutputError, study.StudyError) as exc:
        _stage_error(exc)
    print(f"Variants: {df.height} SNV(s) at {df['position'].n_unique()} position(s). Wrote {core_variants.CORE_VARIANTS_FILE}")


def variants_cmd(args: argparse.Namespace) -> None:
    study_dir = Path(args.study_dir)
    run_variants(study_dir)
    _record(study_dir, "11-variants", [core_variants.CORE_VARIANTS_FILE], ("variants.",))


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
    study_dir = Path(args.study_dir)
    run_master(study_dir)
    _record(study_dir, "11-master", [core_master.MASTER_FILE, core_master.LONG_FILE, core_master.QC_FILE])


def run_contrast(study_dir: Path) -> None:
    """Workflow step 12: runs the contrast plan of the triage (generated from the cohorts when the study has none)."""
    try:
        df = core_contrast.run(study_dir)
        study.clear_pending(study_dir, ("contrast.",))
    except core_contrast.ContrastPending as exc:
        study.raise_pending(study_dir, exc.items, ("contrast.",))
    except (core_contrast.ContrastError, core_master.MasterError, core_cohorts.CohortError, models.DependencyOutputError,
            study.StudyError) as exc:
        _stage_error(exc)
    print(f"Contrast: {df.height} row(s). Wrote {core_contrast.CONTRAST_FILE}, {core_contrast.META_FILE}, "
          f"{core_contrast.TOP_TSV_FILE}, {core_contrast.TOP_JSON_FILE}")


def contrast_cmd(args: argparse.Namespace) -> None:
    study_dir = Path(args.study_dir)
    if not args.propose:
        run_contrast(study_dir)
        _record(study_dir, "12-contrast", ["contrasts.yaml", core_contrast.CONTRAST_FILE, core_contrast.META_FILE,
                                           core_contrast.TOP_TSV_FILE, core_contrast.TOP_JSON_FILE], ("contrast.",))
        return
    try:
        path, plan, points = core_contrast.propose(study_dir)
    except (core_contrast.ContrastError, core_master.MasterError, core_cohorts.CohortError, models.DependencyOutputError,
            study.StudyError) as exc:
        _stage_error(exc)
    print(f"Proposed plan: {', '.join(c.name for c in plan.contrasts)}; priority from {plan.priority.direction_contrast}. "
          f"Wrote {path.relative_to(study_dir)}")
    for q in points:
        print(f"  note: {q.question}")
    _record(study_dir, "12-contrast-propose", [path.relative_to(study_dir).as_posix()])


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


def _opening_figures(study_dir: Path) -> list[str]:
    """Draws the study's opening view (its groups and question, no filter) for report.md when the researcher has not
    drawn a view with `figures`, or when that opening view is older than the variant table. Returns the files."""
    from viz import report_figures

    rec = report_figures.load_choices(study_dir)
    choices = study_dir / report_figures.CHOICES_FILE
    master = study_dir / core_master.MASTER_FILE
    stale = rec is not None and rec.get("opening") and master.is_file() and \
        master.stat().st_mtime > choices.stat().st_mtime
    if rec is not None and not stale:
        return []
    args = argparse.Namespace(
        study_dir=str(study_dir), scorer="RNA_SEQ", target=None, control=None, target_label=None, control_label=None,
        question=None, score="raw", group_agg="mean", min_change_pct=0.0, max_control_pct=None, min_sigma=0.0,
        assay=None, source=None, domain=None, one_alt=False, ref_floor=None, floor_applies="warn",
        y_metric="target_mean", top_k=60, color_by="genomic_domain")
    v = _dashboard_view(args)
    locus = core_locus.load_locus(study_dir)
    rec = report_figures.write_figures(
        study_dir, v["data"], v["profile"], v["effects"], v["kept"], v["steps"], show_mode="all",
        color_col="genomic_domain", highlight=v["kept"].head(5)["variant_id"].to_list(), labels=v["labels"],
        window=(locus.window.start, locus.window.end), command=fetch.invocation() + " (visão inicial)", opening=True)
    return [*rec["figures"].values(), report_figures.CHOICES_FILE]


def report_cmd(args: argparse.Namespace) -> None:
    study_dir = Path(args.study_dir)
    try:
        drawn = _opening_figures(study_dir)
    except (core_master.MasterError, core_locus.LocusError, study.StudyError, OSError, ValueError) as exc:
        drawn = []
        print(f"Warning: the opening view's figures could not be drawn ({exc}); report.md will have no dashboard figures.",
              file=sys.stderr)
    if drawn:
        print(f"Figures of the opening view (groups and question of the study, no filter): {len(drawn) - 1} PNG in "
              "figures/dashboard/.")
    try:
        path = core_report.write_report(study_dir)
    except (core_report.ReportError, core_contrast.ContrastError, core_master.MasterError, core_cohorts.CohortError,
            core_locus.LocusError, models.DependencyOutputError, study.StudyError) as exc:
        _stage_error(exc)
    print(f"Report: wrote {path}")
    fig1 = study_dir / core_report.FIG1_SVG
    _record(study_dir, "13-report", ["report.md", *drawn, *([core_report.FIG1_SVG] if fig1.is_file() else [])])


WORKBENCH_APP = SKILL_DIR / "assets" / "workbench.py"      # legacy triage workbench (reproduce stage metrics)
DASHBOARD_APP = SKILL_DIR / "assets" / "dashboard.py"      # analysis dashboard (groups, question, filters chosen on the page)
WORKBENCH_DEPS = list(deps.WORKBENCH_DEPS)
DEFAULT_PORT = 2718


def _wait_until_up(url: str, proc: Any, timeout_s: float = 600.0) -> bool:
    """True once marimo answers on url/health; False if the process exits first or the wait times out."""
    import time
    import urllib.request

    deadline = time.monotonic() + timeout_s
    while time.monotonic() < deadline:
        if proc.poll() is not None:
            return False
        try:
            with urllib.request.urlopen(f"{url}/health", timeout=2) as resp:
                if resp.status == 200:
                    return True
        except OSError:
            pass
        time.sleep(1)
    return False


def _free_port(preferred: int) -> int:
    """`preferred` when nothing listens on it, else a free port chosen by the system."""
    import socket

    for port in (preferred, 0):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            try:
                s.bind(("127.0.0.1", port))
            except OSError:
                continue
            return s.getsockname()[1]
    raise RuntimeError("no free TCP port on 127.0.0.1")


def _place_app(app_dir: Path, source: Path) -> Path:
    """The notebook inside the study (marimo writes __marimo__/ next to it, never inside the skill).

    A symlink where the system allows one; otherwise a copy (Windows without developer mode refuses
    symlinks). Either way the app finds scripts/ through EXPRESSION_ANALYSIS_SCRIPTS.
    """
    app_dir.mkdir(exist_ok=True)
    app = app_dir / source.name
    if app.is_symlink() or app.exists():
        app.unlink()
    try:
        app.symlink_to(source)
    except OSError:
        shutil.copyfile(source, app)
    return app


def dashboard_cmd(args: argparse.Namespace) -> None:
    """Opens the analysis dashboard (default), the legacy workbench (--legacy), the dashboard in the marimo
    editor (--edit), or exports a static HTML (marimo runs in an ephemeral uv environment)."""
    import subprocess

    study_dir = Path(args.study_dir).resolve()
    source = WORKBENCH_APP if args.legacy else DASHBOARD_APP
    try:
        if args.legacy:
            workbench_helpers.load_workbench_datasets(study_dir)
        else:
            core_master.load_master(study_dir)
            core_master.load_registry(study_dir)
    except (core_contrast.ContrastError, core_master.MasterError, core_cohorts.CohortError, core_locus.LocusError,
            study.StudyError) as exc:
        _stage_error(exc)
    if shutil.which("uv") is None:
        _stage_error(RuntimeError("'uv' is not on PATH; install it with the 'uv' skill."))
    env = {**os.environ, "EXPRESSION_ANALYSIS_SCRIPTS": str(SKILL_DIR / "scripts"), "PYTHONIOENCODING": "utf-8"}
    base = ["uv", "run", "--no-project", "--python", deps.DEP_PYTHON, *deps.uv_with(deps.WORKBENCH_DEPS), "marimo"]
    packages = ", ".join(deps.WORKBENCH_DEPS)
    step = "14-workbench" if args.legacy else "14-dashboard"
    if args.export_html:
        app = _place_app(study_dir / ".workbench", source)
        out = Path(args.export_html).resolve()
        out.parent.mkdir(parents=True, exist_ok=True)
        cmd = [*base, "export", "html", "--no-sandbox", str(app), "-o", str(out), "--force", "--", "--study-dir", str(study_dir)]
        code = subprocess.call(cmd, cwd=study_dir, env=env)
        if code == 0:
            inside = out.is_relative_to(study_dir)
            versions = {"marimo_command": deps.display_command(cmd), "packages": packages}
            if not inside:
                versions["export_outside_study"] = str(out)
            _record(study_dir, step, [out.relative_to(study_dir).as_posix()] if inside else [], versions=versions)
        sys.exit(code)
    env["EXPRESSION_DASHBOARD_LINK"] = "1"  # live page: follows .dashboard/agent_view.json, writes researcher_view.json
    port = _free_port(args.port)
    if port != args.port:
        print(f"Port {args.port} is busy; using {port}.")
    url = f"http://localhost:{port}"
    if args.edit:
        # The editor opens the skill's own file, so changes to cells, charts and defaults are kept in the skill
        # (and in its git history); the study is passed as an argument.
        cmd = [*base, "edit", "--no-sandbox", "--no-token", str(source), "--headless", "--port", str(port), "--",
               "--study-dir", str(study_dir)]
        print(f"Dashboard editor: open {url} in a browser. Edits are saved to {source} (Ctrl+C here stops it).", flush=True)
        try:
            sys.exit(subprocess.call(cmd, cwd=study_dir, env=env))
        except KeyboardInterrupt:
            sys.exit(130)
    app = _place_app(study_dir / ".workbench", source)
    cmd = [*base, "run", "--no-sandbox", str(app), "--headless", "--port", str(port), "--", "--study-dir", str(study_dir)]
    print(f"Dashboard: starting on {url} (the first start can take a minute while the packages install)…", flush=True)
    code = 1
    try:
        proc = subprocess.Popen(cmd, cwd=study_dir, env=env)
        # Only say it is up once marimo answers: an agent that starts it in the background reads this line.
        if _wait_until_up(url, proc):
            print(f"Dashboard READY: open {url} in a browser (Ctrl+C here stops it).", flush=True)
        else:
            print(f"Dashboard FAILED: marimo exited with code {proc.poll()} before answering on {url}; nothing is "
                  "being served. Do not give the researcher the link; read the error above.", file=sys.stderr, flush=True)
        code = proc.wait()
    except KeyboardInterrupt:
        code = 130
    finally:
        # Recorded after the session ends, with its exit code, so a launch that failed is not reported as served.
        _record(study_dir, step, [], versions={"url": url, "marimo_command": deps.display_command(cmd),
                                               "packages": packages, "exit_code": str(code)})
    sys.exit(code)


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


# --------------------------------------------------------------------------- captures (fetch) and exploration

FETCH_KINDS = ("gtf", "regulatory", "sequence", "ols", "ols-term", "tracks", "scores", "avi", "links", "interpret")


def _study_curies(study_dir: Path, config: study.StudyConfig) -> list[str]:
    """CURIEs whose tracks the study needs: study.yaml cohorts, the cohort decisions and their replace:/add: answers."""
    curies = [c for side in ("target", "control") for c in config.data["cohorts"][side]["curies"]]
    try:
        curies += sorted(core_cohorts.chosen_curies(study.load_decisions(study_dir)))
    except core_cohorts.CohortError:
        pass  # cohort.target / cohort.control not decided yet
    return list(dict.fromkeys(curies))


def _queried(study_dir: Path) -> set[tuple[str, str]]:
    return core_cohorts._queried_pairs(study_dir)


def fetch_cmd(args: argparse.Namespace) -> None:
    """Runs one dependency capture of the Workflow into the study and records it in manifest.json."""
    study_dir = Path(args.study_dir)
    ws = fetch.Workspace(root=study_dir, invoked_by=fetch.invocation(), dry_run=args.dry_run)
    try:
        config = study.load_study(study_dir)
        gene = config.gene_symbol
        kind = args.kind
        if kind == "gtf":
            fetch.gtf(ws, gene)
        elif kind in ("regulatory", "sequence"):
            region = args.region or study_regions(study_dir)["regulatory_search"]["region"]
            (fetch.regulatory if kind == "regulatory" else fetch.sequence)(ws, region)
        elif kind == "ols":
            if not args.query:
                raise fetch.FetchError("fetch ols needs --query TEXT (repeat it for several tissues or cell types)")
            for text in args.query:
                path = fetch.ols_search(ws, text, args.ontology, args.rows)
                if not args.dry_run:
                    _print_hits(path)
        elif kind == "ols-term":
            if not args.curie:
                raise fetch.FetchError("fetch ols-term needs --curie CURIE (repeatable)")
            for curie in args.curie:
                fetch.ols_term(ws, curie)
        elif kind == "tracks":
            scorers = args.scorer or list(config.data["scorers"])
            if args.catalogue:
                for s in scorers:
                    if args.force or not (study_dir / "raw" / f"tracks_{s}_catalogue.parquet").is_file():
                        fetch.catalogue(ws, s)
            curies = args.curie or _study_curies(study_dir, config)
            done = set() if args.force else _queried(study_dir)
            for c in curies:
                for s in scorers:
                    exists = (study_dir / "raw" / f"tracks_{s}_{fetch.curie_file(c)}.parquet").is_file()
                    if not args.force and (exists or (s, c) in done):
                        print(f"Skip {s} x {c}: already captured (use --force to capture again).")
                        continue
                    fetch.tracks(ws, c, s)
        elif kind in ("avi", "scores"):
            regions = study_regions(study_dir)
            if regions["window"] is None:
                raise core_locus.LocusError("the window needs processed/locus.json; run expression.py locus first.")
            window, max_window = regions["window"]["region"], int(regions["max_window"])
            if kind == "avi":
                fetch.avi_region(ws, window, max_window)
            else:
                resolution = core_cohorts.load_resolution(study_dir)
                curies = resolution["effective_target_curies"] + resolution["effective_control_curies"]
                for s in args.scorer or resolution["effective_scorers"]:
                    fetch.scores(ws, window, max_window, s, curies, gene)
        elif kind == "links":
            fetch.atlas_links(ws)
        elif kind == "interpret":
            if not args.ontology_id:
                raise fetch.FetchError("fetch interpret needs --ontology-id CURIE (a target CURIE the researcher chose)")
            top = json.loads((study_dir / core_contrast.TOP_JSON_FILE).read_text(encoding="utf-8"))
            variants = args.variant or [r["variant"] for r in top][: args.top]
            if not args.dry_run:
                (study_dir / "processed" / "interpret_variants.txt").write_text("\n".join(variants) + "\n",
                                                                               encoding="utf-8")
            label = args.tissue or config.data["cohorts"]["target"]["label"]
            for v in variants:
                fetch.interpret(ws, v, args.ontology_id, gene, label)
    except (fetch.FetchError, deps.DependencyError, core_locus.LocusError, core_cohorts.CohortError,
            core_contrast.ContrastError, models.DependencyOutputError, study.StudyError, manifest.ManifestError,
            OSError, ValueError) as exc:
        _stage_error(exc)


def _print_hits(path: Path) -> None:
    """The OLS hits as label, CURIE and the start of the definition, for showing to the researcher."""
    hits = models.load_raw(path)[0].terms
    seen = set()
    print(f"OLS hits in {path.name}:")
    for h in hits:
        if h.obo_id in seen:
            continue
        seen.add(h.obo_id)
        definition = (h.description[0] if h.description else "")[:110]
        print(f"  {h.obo_id:<16} {h.label} -- {definition}")


def explore_cmd(args: argparse.Namespace) -> None:
    """Workflow step 0: read-only exploration into <out>/exploration (no study yet)."""
    if args.activity_only:
        try:
            path = core_explore.write_activity(Path(args.out) / "exploration", args.floor)
        except (core_explore.ExploreError, OSError) as exc:
            _stage_error(exc)
        print(f"Wrote {path} (no capture). Show it to the researcher to choose or change the target and control groups.")
        return
    if not args.gene:
        _stage_error(ValueError("--gene is required (except with --activity-only)"))
    try:
        summary = core_explore.run(Path(args.out), args.gene.strip(), args.tissue or [],
                                   tuple(args.scorers), args.flank_bp, args.floor, fetch.invocation(), args.dry_run)
    except (core_explore.ExploreError, fetch.FetchError, deps.DependencyError, core_locus.LocusError,
            models.DependencyOutputError, manifest.ManifestError, OSError, ValueError) as exc:
        _stage_error(exc)
    if not args.dry_run:
        print(f"Exploration: wrote {summary} and {summary.with_name(core_explore.EXPLORATION_FILE)}. Show summary.md "
              "to the researcher and ask for the window, target, control and scorers before study.py init "
              "--from-exploration.")


def _choice_args(p: argparse.ArgumentParser) -> None:
    """The dashboard's analysis choices as options (summarize and figures)."""
    p.add_argument("--scorer", default="RNA_SEQ", help="RNA_SEQ (gene expression; default), CAGE, DNASE or ATAC")
    p.add_argument("--target", action="append", help="Target CURIE (repeatable; default: the study's target tissues)")
    p.add_argument("--control", action="append", help="Control CURIE (repeatable; default: the study's control tissues)")
    p.add_argument("--shown", action="append",
                   help="Context tissue CURIE (repeatable): shown in the tables and figures, outside the group values")
    p.add_argument("--target-label", help="Name of the target group (default: study.yaml)")
    p.add_argument("--control-label", help="Name of the control group (default: study.yaml)")
    p.add_argument("--question", choices=["repression", "increase", "either"],
                   help="Direction in the target (default: the study's decision analysis.question, else either)")
    p.add_argument("--score", choices=["raw", "quantile"], default="raw",
                   help="Value compared: the model's raw score (default; percent changes) or the Atlas quantile")
    p.add_argument("--group-agg", choices=["mean", "median"], default="mean",
                   help="How each group's tissues are combined: mean (default) or median")
    p.add_argument("--min-change-pct", type=float, default=0.0, help="Smallest change in every target tissue, in %% (0 = off)")
    p.add_argument("--max-control-pct", type=float, help="Largest change allowed in the control, in %% (default: no limit)")
    p.add_argument("--min-sigma", type=float, default=0.0, help="Target - control beyond this many noise sigmas, in the question's direction (both for either; 0 = off)")
    p.add_argument("--assay", action="append", help="Assay titles kept (repeatable; default all)")
    p.add_argument("--source", action="append", help="Data sources kept, e.g. encode, gtex (repeatable; default all)")
    p.add_argument("--domain", action="append", help="Genomic domains kept (repeatable; default all)")
    p.add_argument("--one-alt", action="store_true", help="Keep one alternative per position (the largest effect)")
    p.add_argument("--ref-floor", type=float, help="Reference expression floor (RNA_SEQ REF; default: off)")
    p.add_argument("--floor-applies", choices=["warn", "control", "target", "both"], default="warn",
                   help="What the floor does: warn (default), drop controls below it, skip a target below it, or both")
    p.add_argument("--y-metric", default="target_mean",
                   choices=["target_mean", "target_toward", "control_mean", "diff", "diff_sigma", "avi_phred"],
                   help="Order of the list (the dashboard's 'Ordenar por')")


AGENT_VIEW = ".dashboard/agent_view.json"            # written by `view`; the open dashboard follows it
RESEARCHER_VIEW = ".dashboard/researcher_view.json"  # written by the open dashboard: what is on screen now
# option (argparse dest) -> key in the view files (the fields of core.analysis.Profile, plus how to show)
VIEW_KEYS = {"scorer": "scorer", "target": "target", "control": "control", "shown": "shown", "target_label": "target_label",
             "control_label": "control_label", "question": "question", "score": "score", "group_agg": "group_agg",
             "min_change_pct": "min_change_pct", "max_control_pct": "max_control_pct", "min_sigma": "min_diff_sigma",
             "assay": "assays", "source": "sources", "domain": "domains", "one_alt": "one_alt_per_position",
             "ref_floor": "ref_floor", "floor_applies": "ref_floor_applies", "y_metric": "y_metric", "top_k": "top_k",
             "show": "show_mode", "color_by": "color_by", "highlight": "highlight"}


def _researcher_view(study_dir: Path) -> dict:
    """What the researcher has on screen in the open dashboard (written by the page on every change)."""
    path = Path(study_dir) / RESEARCHER_VIEW
    if not path.is_file():
        _stage_error(RuntimeError(
            f"{path} not found: the dashboard has not been opened on this study since it started writing its view. "
            "Open it with expression.py dashboard (and the page in a browser), then run this again."))
    try:
        view = json.loads(path.read_text(encoding="utf-8"))
    except ValueError as exc:
        _stage_error(RuntimeError(f"{path} could not be read ({exc}); run this again"))
    if not isinstance(view, dict) or not isinstance(view.get("choices"), dict):
        _stage_error(RuntimeError(f"{path} has no choices; reload the dashboard page"))
    return view


def _use_dashboard_choices(args: argparse.Namespace) -> None:
    """--from-dashboard: the options take the values the researcher has on screen (other choice options are ignored)."""
    choices = _researcher_view(Path(args.study_dir))["choices"]
    for dest, key in VIEW_KEYS.items():
        if key in choices:
            value = choices[key]
            if dest == "top_k":
                value = value or 0
            setattr(args, dest, value)
    print(f"Choices taken from the open dashboard ({RESEARCHER_VIEW}).", file=sys.stderr)


def view_cmd(args: argparse.Namespace) -> None:
    """Without choices: prints what the researcher has on screen. With choices: sends them to the open dashboard,
    on top of what is on screen (--reset: on top of the study's opening view); the page follows without reloading."""
    study_dir = Path(args.study_dir).resolve()
    if not (study_dir / study.STUDY_FILE).is_file():
        _stage_error(RuntimeError(f"{study_dir} is not a study directory"))
    given = {VIEW_KEYS[d]: v for d, v in vars(args).items() if d in VIEW_KEYS}
    if args.no_max_control:
        given["max_control_pct"] = None
    if args.no_floor:
        given["ref_floor"] = None
    if args.clear_highlight:
        given["highlight"] = []
    for flag, key in (("no_control", "control"), ("no_shown", "shown"), ("all_domains", "domains"),
                      ("all_assays", "assays"), ("all_sources", "sources")):
        if getattr(args, flag):
            given[key] = []  # an empty selection: no control or context tissue; every region, assay or source
    if args.full_window:
        given["window"] = "full"
    elif args.window:
        given["window"] = sorted(args.window)
    if args.reset:  # the opening view also clears the highlights and shows the whole window
        given.setdefault("highlight", [])
        given.setdefault("window", "full")
        sent_reset = not (set(given) - {"highlight", "window"})
    else:
        sent_reset = False
    if not given and not args.reset:
        view = _researcher_view(study_dir)
        print(json.dumps(view, indent=1, ensure_ascii=False))
        print(f"\n{view.get('n_pass')} of {view.get('n_total')} variant(s) pass these choices (as of {view.get('updated_at')}). "
              "To explain or draw exactly this view: expression.py summarize / figures --from-dashboard.", file=sys.stderr)
        return
    choices: dict = {}
    seen = (study_dir / RESEARCHER_VIEW).is_file()
    if not args.reset and seen:
        choices = dict(_researcher_view(study_dir)["choices"])
        choices.pop("highlight", None)  # the researcher's highlights stay as they are unless the agent sends some
    choices.update(given)
    out = study_dir / AGENT_VIEW
    out.parent.mkdir(exist_ok=True)
    body = {"sent_at": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "message": args.message or "", "choices": choices}
    out.write_text(json.dumps(body, indent=1), encoding="utf-8")  # ASCII JSON, written in place: the page watches this file
    changed = "the study's opening view (highlights cleared, whole window)" if sent_reset else \
        ", ".join(f"{k}={json.dumps(v, ensure_ascii=False)}" for k, v in given.items())
    print(f"Sent to the dashboard ({AGENT_VIEW}): {changed}"
          + (" on top of the study's opening view." if args.reset or not seen else " on top of what is on screen."))
    print("The open page follows by itself (no reload). Read it back with expression.py view."
          if seen else "No open dashboard was seen on this study yet; the view applies when it is opened.")


def _decided_question(study_dir: Path) -> str:
    """The researcher's question (decision analysis.question), or "either" for a study without it."""
    entry = study.load_decisions(Path(study_dir)).get("analysis.question")
    return entry.value if entry is not None and entry.value in ("repression", "increase", "either") else "either"


def _dashboard_view(args: argparse.Namespace) -> dict:
    """The dashboard's computation for the choices in `args` (shared by summarize and figures)."""
    import dataclasses

    from core import analysis as A

    study_dir = Path(args.study_dir)
    if getattr(args, "from_dashboard", False):
        _use_dashboard_choices(args)
    try:
        data = A.load_study_data(study_dir)
        config = study.load_study(study_dir)
    except (core_master.MasterError, study.StudyError) as exc:
        _stage_error(exc)
    reg = data.registry
    cohorts = config.data.get("cohorts", {})
    profile = A.Profile(
        scorer=args.scorer,
        target=args.target or sorted(set(reg.filter(reg["cohort"] == "target")["ontology_curie"].to_list())),
        control=args.control or sorted(set(reg.filter(reg["cohort"] == "control")["ontology_curie"].to_list())),
        target_label=args.target_label or cohorts.get("target", {}).get("label") or "Alvo",
        control_label=args.control_label or cohorts.get("control", {}).get("label") or "Controle",
        shown=getattr(args, "shown", None) or [],
        assays=args.assay or [], sources=args.source or [], question=args.question or _decided_question(study_dir),
        score=args.score, group_agg=args.group_agg,
        min_change_pct=args.min_change_pct, max_control_pct=args.max_control_pct,
        min_diff_sigma=args.min_sigma, ref_floor=args.ref_floor, ref_floor_applies=args.floor_applies,
        domains=args.domain or [], one_alt_per_position=args.one_alt, y_metric=args.y_metric,
        top_k=getattr(args, "top_k", None) or None, color_by=getattr(args, "color_by", "genomic_domain"))
    noise = A.noise_model(data, profile.scorer, profile.assays or None, profile.sources or None, score=profile.score)
    # Same order of operations as the dashboard: controls below the expression floor leave before computing
    used = profile
    if profile.ref_floor is not None and profile.ref_floor_applies in ("control", "both"):
        flags = A.expression_flags(data, profile)
        used = dataclasses.replace(profile, control=[c for c in profile.control if flags.get(c) != "below floor"])
    effects = A.compute(data, used, noise)
    kept, steps = A.apply(effects, data, used, noise)
    tissues = A.tissue_table(data, profile.scorer)
    labels = {r["ontology_curie"]: r["tissue"] for r in tissues.iter_rows(named=True)}
    return {"study_dir": study_dir, "data": data, "config": config, "profile": profile, "used": used, "noise": noise,
            "effects": effects, "kept": kept, "steps": steps, "labels": labels}


def figures_cmd(args: argparse.Namespace) -> None:
    """Draws the dashboard's charts for the given choices into figures/dashboard/ for report.md."""
    from viz import report_figures

    v = _dashboard_view(args)
    study_dir, kept = v["study_dir"], v["kept"]
    highlight = list(args.highlight or [])
    if not highlight and args.auto_highlight:
        highlight = kept.head(args.auto_highlight)["variant_id"].to_list()
    locus = core_locus.load_locus(study_dir)
    try:
        rec = report_figures.write_figures(
            study_dir, v["data"], v["profile"], v["effects"], kept, v["steps"], show_mode=args.show,
            color_col=args.color_by, highlight=highlight, labels=v["labels"],
            window=(locus.window.start, locus.window.end), command=fetch.invocation())
    except (OSError, ValueError) as exc:
        _stage_error(exc)
    print(f"Figures: {', '.join(rec['figures']) or 'none'} for {kept.height} variant(s) that pass "
          f"({rec['n_shown']} shown). Wrote {report_figures.CHOICES_FILE}; run expression.py report to include them.")
    _record(study_dir, "14-figures", [*rec["figures"].values(), report_figures.CHOICES_FILE])


def summarize_cmd(args: argparse.Namespace) -> None:
    """Text summary of the dashboard view for the same choices (the agent reads it to explain the dashboard).

    Same computation as assets/dashboard.py (core.analysis); defaults are the dashboard's opening view:
    the study's target and control tissues, any direction, no filter. Writes nothing.
    """
    import polars as pl

    from core import analysis as A
    from viz import dashboard_helpers as H

    v = _dashboard_view(args)
    data, config, profile, noise = v["data"], v["config"], v["profile"], v["noise"]
    effects, kept, steps, labels = v["effects"], v["kept"], v["steps"], v["labels"]
    used = v["used"]
    scorer = profile.scorer
    desc = H.describe_tissues(data, scorer)
    what = A.SCALES.get(scorer, (None, "o sinal previsto"))[1]
    out = [f"# {config.gene_symbol}: resumo da visão do dashboard ({scorer}; % d{what})", "",
           H.headline(kept, effects, profile, labels, A.ref_by_tissue(data), missing=A.without_tracks(data, profile)),
           "", "## Grupos e medições", ""]
    for role, curies, name in (("Alvo", profile.target, profile.target_label), ("Controle", profile.control, profile.control_label)):
        out.append(f"**{role} ({name})**")
        out += [f"- {labels.get(c, c)} ({c}): {desc.get(c, 'sem medições deste tipo')}" for c in curies] or ["- nenhum"]
        out.append("")
    dropped = [c for c in profile.control if c not in used.control]
    if dropped:
        out += [f"Controles abaixo do piso de expressão ({H._n(profile.ref_floor, 3)}), fora do cálculo: "
                + ", ".join(labels.get(c, c) for c in dropped) + ".", ""]
    ref = A.ref_by_tissue(data)
    out += ["## Expressão de referência (REF) e piso", ""]
    if not ref:
        out += ["Sem REF neste estudo (não há a variante de referência da exploração em raw/), então o piso não se aplica.", ""]
    else:
        top = max(ref.values())
        floor = profile.ref_floor
        out += ["REF é a expressão do gene prevista no tecido, sem variante (escala do AlphaGenome; o pseudocount do modelo é "
                "0,001). Onde o gene está quase desligado, mudanças pequenas viram % grandes e instáveis. Compare os tecidos "
                "pela coluna '% do maior REF': a exploração chama de 'não expresso' o que fica abaixo de 5% do maior REF "
                f"(aqui, {H._n(0.05 * top, 3)}); 0,01 (10× o pseudocount) é só o limite absoluto mínimo. O piso (`--ref-floor N`) usa "
                "o valor de REF: `--floor-applies warn` só avisa, `control` "
                "tira do cálculo os controles abaixo dele, `target` não deixa passar nada se um tecido do alvo estiver "
                "abaixo, `both` faz os dois.", "",
                "| tecido | grupo | REF | % do maior REF | com o piso atual |", "|---|---|---|---|---|"]
        for role, curies in (("alvo", profile.target), ("controle", profile.control)):
            for c in curies:
                r = ref.get(c)
                if r is None:
                    out.append(f"| {labels.get(c, c)} | {role} | — | — | sem REF |")
                    continue
                state = ("piso desligado" if floor is None else
                         ("abaixo do piso" + (", fora do cálculo" if role == "controle" and profile.ref_floor_applies in ("control", "both") else
                                              ", nenhuma variante passa" if role == "alvo" and profile.ref_floor_applies in ("target", "both") else
                                              ", só aviso")) if r < floor else "acima do piso")
                out.append(f"| {labels.get(c, c)} | {role} | {H._n(r, 4)} | {H._n(100 * r / top, 1)}% | {state} |")
        out.append("")
    out += ["## Ruído entre medições do mesmo tecido", ""]
    overall = A.noise_overall(data, scorer, profile.assays or None, profile.sources or None, score=profile.score)
    if overall is None:
        out.append("Não disponível (nenhum tecido com duas medições nas escolhas atuais); o filtro de ruído fica de fora.")
    else:
        size = (lambda v: f"{H._n(A.pct(v, scorer), 1)}%") if A.uses_pct(profile) else (lambda v: f"{H._n(v, 3)} (quantil)")
        out.append(f"Duas medições do mesmo tecido discordam em até {size(overall[0])} (p95) e "
                   f"{size(overall[1])} (p99), juntando todos os tecidos com duas ou mais medições.")
    out += ["", "## Etapas dos filtros", "", "| etapa | variantes |", "|---|---|"]
    out += [f"| {s} | {n} |" for s, n in steps]
    # Order: the one of 'Ordenar por', as in the dashboard's table and headline (kept comes sorted by it)
    top = kept.head(args.top)
    eff = [c for c in top.columns if c.startswith("eff_")]
    out += ["", f"## Primeiras {top.height} variantes", "",
            (lambda m: f"| variante | região | {m} alvo | {m} controle | "
             + " | ".join(f"{m} {labels.get(c[4:], c[4:])}" for c in eff))("%" if A.uses_pct(profile) else "quantil")
            + " | modalidade AVI |",
            "|" + "---|" * (5 + len(eff))]
    fmt = (lambda v: "—" if v is None else H._n(A.pct(v, scorer), 1, True)) if A.uses_pct(profile) else \
        (lambda v: "—" if v is None else H._n(v, 3, True))
    for r in top.iter_rows(named=True):
        out.append(f"| {r['variant_id']} | {r['genomic_domain']} | {fmt(r['target_mean'])} | {fmt(r['control_mean'])} | "
                   + " | ".join(fmt(r[c]) for c in eff) + f" | {r['top_feature']} |")
    out += ["", "Predições do AlphaGenome; uso somente em pesquisa. Variantes de modalidade AVI *Splicing* "
            "provavelmente agem por splicing, não por transcrição."]
    print("\n".join(out))


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

DEPENDENCIES = deps.DEPENDENCIES
REQUIRED_KEYS = ["ALPHAGENOME_API_KEY"]
OPTIONAL_KEYS = []


def _key_present(name: str) -> bool:
    """True if the project's .env sets the variable; the value is never stored or printed."""
    return deps.key_present(name)


SETUP_INSTALL = ("npx skills add romulorcruz/alphagenome-skills-setup --skill alphagenome-skills-setup "
                 "--agent <agent> --copy -y")


SKILL_PURPOSE = {"uv": "to run its scripts", "credentials": "to store the AlphaGenome key without showing it",
                 "alphagenome-single-variant-analysis": "to interpret the top variants",
                 "alphagenome-atlas-website-links": "to link the top variants to the AlphaGenome Atlas"}
SETUP_SKILL = "alphagenome-skills-setup"  # prepares the environment step by step, with the researcher's approval


def _setup_item(key: str, title: str, question: str, why: str, identifiers: dict, options: list[dict],
                recommended: str | None, handled_by: str) -> dict:
    """One missing prerequisite in the shape of a study question (study.py pending --json), so the agent and the
    setup skill present it like any other question: in the researcher's language, labels and effects, no codes."""
    return {"key": key, "title": title, "question": question, "why": why, "identifiers": identifiers,
            "options": options, "recommended": recommended, "handled_by": handled_by}


def prerequisite_items() -> list[dict]:
    """What is missing for this skill to run, one item per prerequisite (empty when all is in place)."""
    items: list[dict] = []
    installed: dict[str, Path] = {}
    later = {"code": "later", "label": "Not now", "effect": "Nothing is installed; the analysis cannot run until this is in place."}
    if shutil.which("uv") is None:
        items.append(_setup_item(
            "setup.uv", "Install uv (runs the analysis scripts)",
            "uv is not installed. It runs every script of the analysis in its own environment with fixed package "
            "versions, without changing the computer's Python. Install it?",
            "Without uv no script of the analysis can run.", {"program": "uv", "source": "https://docs.astral.sh/uv/"},
            [{"code": "install", "label": "Install uv", "effect": "The 'uv' skill installs it for this user account."},
             later], "install", "uv"))
    for name, (origin, command) in DEPENDENCIES.items():
        found = deps.find_skill(name)
        if found is None:
            items.append(_setup_item(
                f"setup.skill:{name}", f"Install the skill {name}",
                f"The analysis uses the skill {name} ({origin}) {SKILL_PURPOSE.get(name, 'to fetch part of its data')}. "
                "Install it in this project?",
                "Installing a skill brings third-party code and, for some, terms of use; each one is approved on its own.",
                {"skill": name, "origin": origin, "install command": command},
                [{"code": "install", "label": f"Install {name}", "effect": f"Runs: {command}"}, later],
                "install", SETUP_SKILL))
        else:
            installed[name] = found
    searched = ", ".join(str(d) for d in deps.licence_dirs(installed))
    for filename, owner in deps.missing_licences(installed):
        items.append(_setup_item(
            f"setup.licence:{owner}", f"Terms of use of {owner}",
            f"The skill {owner} has terms of use that must be accepted before it is used. The agent shows them to you "
            "now. Do you accept them?",
            "The notice file records that these terms were seen and accepted; one acceptance never covers another skill.",
            {"skill": owner, "notice file": f".licenses/{filename}"},
            [{"code": "accept", "label": "I have read and accept the terms",
              "effect": f"The agent creates .licenses/{filename} recording the acceptance."},
             {"code": "decline", "label": "I do not accept",
              "effect": f"Nothing is recorded; the steps that use {owner} cannot run."}],
            None, SETUP_SKILL))
        items[-1]["agent_note"] = (f"Show the terms from the Prerequisites section of {owner}'s SKILL.md before asking "
                                   f"(looked for the notice in: {searched}).")
    for key in REQUIRED_KEYS:
        if not _key_present(key):
            items.append(_setup_item(
                f"setup.key:{key}", "AlphaGenome access key",
                f"The AlphaGenome key ({key}) is not set. It is needed to fetch the predictions. Set it now?",
                "The key identifies the researcher to the AlphaGenome service; it is personal, kept in the project's "
                ".env (ignored by git), and never shown.",
                {"variable": key, "stored in": ".env of the project",
                 "get one": "https://alphagenome.google/api"},
                [{"code": "set", "label": "Set the key",
                  "effect": "The setup skill opens the key page and the .env file for the researcher to paste it."}, later],
                "set", SETUP_SKILL))
    return items


def check(args: argparse.Namespace) -> None:
    items = prerequisite_items()
    for f in deps.env_unprotected():
        print(f"Warning: {f} holds keys but git does not ignore it; add '.env' to .gitignore before any commit "
              "(an AlphaGenome key is personal and must not be shared).", file=sys.stderr)
    for key in OPTIONAL_KEYS:
        if not _key_present(key):
            print(f"Warning: {key} is not set; requests will run at the lower rate limit.", file=sys.stderr)
    if args.json:
        print(json.dumps({"ready": not items, "setup_skill": SETUP_SKILL,
                          "setup_skill_installed": deps.find_skill(SETUP_SKILL) is not None, "items": items},
                         indent=2, ensure_ascii=False))
        sys.exit(1 if items else 0)
    if not items:
        print("All prerequisites satisfied")
        return
    print(f"{len(items)} prerequisite(s) missing:", file=sys.stderr)
    for it in items:
        print(f"- {it['title']}", file=sys.stderr)
    if deps.find_skill(SETUP_SKILL) is not None:
        print(f"Next: hand the setup to the skill '{SETUP_SKILL}' (it asks the researcher item by item and installs "
              "only what they approve). Do not install anything from this skill.", file=sys.stderr)
    else:
        print(f"Next: the skill '{SETUP_SKILL}' installs these with the researcher's approval, and it is not installed. "
              f"Ask the researcher whether to install it (repository private for now: they need GitHub access), then: "
              f"{SETUP_INSTALL} (with <agent> the id of their agent app, e.g. claude-code, opencode, antigravity). Do not "
              "install the dependency skills from this skill.", file=sys.stderr)
    sys.exit(1)


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
    deps.use_utf8()
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("explore", help="Workflow step 0: read-only exploration before the study exists",
                       description="Gene structure, Regulatory Build, candidate windows, OLS terms per tissue text, Atlas "
                                   "tracks per term x scorer and the gene's predicted activity per tissue "
                                   "(RNA_SEQ_ACTIVE), written to OUT/exploration/ with summary.md (Portuguese).")
    p.add_argument("--gene", help="HGNC symbol (required, except with --activity-only)")
    p.add_argument("--out", required=True, help="Folder for the exploration (OUT/exploration/); later the study's parent")
    p.add_argument("--activity-only", action="store_true",
                   help="Only (re)write OUT/exploration/activity.md from the captured baseline: every catalogue term with "
                        "the gene's predicted activity, for choosing or changing the groups. No network")
    p.add_argument("--tissue", action="append", help="A tissue or cell type in the researcher's words (repeatable)")
    p.add_argument("--scorers", nargs="+", default=list(core_explore.DEFAULT_SCORERS), help="Atlas scorers to survey")
    p.add_argument("--flank-bp", type=int, default=core_explore.DEFAULT_FLANK_BP,
                   help="Search Regulatory Build features this far around the gene (default 2000, the study default)")
    p.add_argument("--floor", type=float, default=core_explore.DEFAULT_FLOOR,
                   help="Relative activity below which the gene reads as not expressed (default 0.05 of its maximum)")
    p.add_argument("--dry-run", action="store_true", help="Print the first command and stop")
    p.set_defaults(func=explore_cmd)

    p = sub.add_parser("fetch", help="Run one dependency capture into the study and record it in manifest.json",
                       description="Runs a dependency skill script with pinned packages, writes the file the next stage "
                                   "reads and records the exact command. Kinds: gtf (step 2), regulatory and sequence "
                                   "(step 3), ols (5), ols-term (6), tracks [--catalogue] (7), avi (9), scores (10), "
                                   "links and interpret (13).")
    p.add_argument("kind", choices=FETCH_KINDS)
    p.add_argument("--study-dir", required=True, help="Study directory")
    p.add_argument("--region", help="regulatory/sequence: override the region (default: regions regulatory_search)")
    p.add_argument("--query", action="append", help="ols: tissue or cell text (repeatable); file raw/ols_search_<slug>.json")
    p.add_argument("--ontology", default=fetch.DEFAULT_OLS_ONTOLOGIES, help="ols: ontologies (default uberon,cl,efo,clo)")
    p.add_argument("--rows", type=int, default=10, help="ols: hits to keep (default 10)")
    p.add_argument("--curie", action="append", help="ols-term / tracks: CURIE (repeatable; tracks default: the study's)")
    p.add_argument("--scorer", action="append", help="tracks / scores: scorer (repeatable; default: the study's)")
    p.add_argument("--catalogue", action="store_true", help="tracks: also the whole catalogue of each scorer")
    p.add_argument("--force", action="store_true", help="tracks: capture again what was already captured")
    p.add_argument("--ontology-id", help="interpret: target CURIE for the model-API interpretation")
    p.add_argument("--tissue", help="interpret: tissue label for the figures (default: the target cohort label)")
    p.add_argument("--top", type=int, default=5, help="interpret: number of top variants (default 5)")
    p.add_argument("--variant", action="append", help="interpret: a chr:pos:REF>ALT instead of the top list (repeatable)")
    p.add_argument("--dry-run", action="store_true", help="Print the commands; run and record nothing")
    p.set_defaults(func=fetch_cmd)

    p = sub.add_parser("summarize", help="Text summary of the dashboard view for given choices (for the agent)",
                       description="Same computation as the dashboard. Defaults: the study's target and control tissues, any "
                                   "direction, no filter. Writes nothing; prints Markdown.")
    p.add_argument("--study-dir", required=True, help="Study directory")
    _choice_args(p)
    p.add_argument("--from-dashboard", action="store_true",
                   help="Use the choices the researcher has on screen in the open dashboard (other choice options are ignored)")
    p.add_argument("--top", type=int, default=10, help="Variants listed (default 10)")
    p.set_defaults(func=summarize_cmd)

    p = sub.add_parser("view", help="Link with the open dashboard: read what the researcher has on screen, or send a view to it",
                       description="Without choice options: prints what the researcher has on screen now. With them: the open "
                                   "dashboard moves to these choices, on top of what is on screen (--reset: on top of the "
                                   "study's opening view), without reloading the page.")
    p.add_argument("--study-dir", required=True, help="Study directory")
    _choice_args(p)
    p.add_argument("--top-k", type=int, help="'Quantas mostrar' (0 = all)")
    p.add_argument("--show", choices=["kept", "all"], help="'Variantes nos gráficos'")
    p.add_argument("--color-by", choices=["genomic_domain", "passa", "sentido", "ensembl_transcript_region", "top_feature"],
                   help="'Cor no gráfico da janela'")
    p.add_argument("--highlight", action="append", help="Variant ID to highlight (repeatable)")
    for action in p._actions:  # only the options given are sent: no default takes part
        if action.dest not in ("help", "study_dir"):
            action.default = argparse.SUPPRESS
    p.add_argument("--no-max-control", action="store_true", help="Switch off the limit of change in the control")
    p.add_argument("--no-floor", action="store_true", help="Switch off the reference expression floor")
    p.add_argument("--clear-highlight", action="store_true", help="Remove the highlights")
    p.add_argument("--no-control", action="store_true", help="No control tissue (the control filters are then listed as ignored)")
    p.add_argument("--no-shown", action="store_true", help="No context tissue")
    p.add_argument("--all-domains", action="store_true", help="Every region of the window")
    p.add_argument("--all-assays", action="store_true", help="Every assay")
    p.add_argument("--all-sources", action="store_true", help="Every data source")
    p.add_argument("--window", nargs=2, type=int, metavar=("START", "END"),
                   help="Stretch of the window to show (GRCh38 positions), as the page's 'Recortar'")
    p.add_argument("--full-window", action="store_true", help="Show the whole window again ('Janela inteira')")
    p.add_argument("--reset", action="store_true",
                   help="Start from the study's opening view instead of what is on screen (also clears the highlights "
                        "and shows the whole window, unless --highlight or --window is given)")
    p.add_argument("--message", help="One line shown to the researcher on the page, saying what this view is")
    p.set_defaults(func=view_cmd)

    p = sub.add_parser("figures", help="Workflow step 14: the dashboard's charts for given choices, as PNG for report.md",
                       description="Same choices as summarize (and the dashboard), plus how to show them. Writes "
                                   "figures/dashboard/*.png and choices.json; expression.py report then includes them.")
    p.add_argument("--study-dir", required=True, help="Study directory")
    _choice_args(p)
    p.add_argument("--top-k", type=int, default=0, help="'Quantas mostrar': draw only the first K in that order (0 = all)")
    p.add_argument("--show", choices=["kept", "all"], default="kept",
                   help="'Variantes nos gráficos': only those that pass (default) or all, the rest in grey")
    p.add_argument("--color-by", default="genomic_domain",
                   choices=["genomic_domain", "passa", "sentido", "ensembl_transcript_region", "top_feature"],
                   help="'Cor no gráfico da janela' (default genomic_domain)")
    p.add_argument("--highlight", action="append", help="Variant ID to highlight (repeatable), as in the dashboard")
    p.add_argument("--auto-highlight", type=int, default=5,
                   help="Without --highlight, highlight the first N that pass, in that order (default 5; 0 = none)")
    p.add_argument("--from-dashboard", action="store_true",
                   help="Use the choices, the way of showing and the highlights the researcher has on screen in the open dashboard")
    p.set_defaults(func=figures_cmd)

    p = sub.add_parser("check", help="Verify prerequisites; print install commands for what is missing")
    p.add_argument("--json", action="store_true",
                   help="Print the missing prerequisites as questions (same shape as study.py pending --json)")
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

    p = sub.add_parser("dashboard", help="Workflow step 14: interactive analysis dashboard (marimo)")
    p.add_argument("--study-dir", required=True, help="Study directory")
    p.add_argument("--port", type=int, default=DEFAULT_PORT,
                   help="Preferred port of the marimo server (default 2718; a free one is used when it is busy)")
    p.add_argument("--export-html", help="Write a static HTML export to this path instead of serving")
    p.add_argument("--edit", action="store_true",
                   help="Maintainers: open the dashboard in the marimo editor (edits are saved to the skill's assets/dashboard.py)")
    p.add_argument("--legacy", action="store_true", help="Open the previous triage workbench (assets/workbench.py)")
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
