#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["pyyaml==6.0.3", "matplotlib==3.11.2", "polars==1.31.0", "pyarrow==26.0.1", "pandas==3.1.1", "altair==6.3.0"]
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
  summarize  Text summary of the dashboard view for given choices (same computation; for the agent to
             explain what the researcher sees). Writes nothing.
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
import json
import os
import shutil
import sys
from pathlib import Path
from typing import Callable

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
    keys = [k for k in study.load_decisions(study_dir) if k.startswith(areas)] if areas else []
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
    study_dir = Path(args.study_dir)
    run_locus(study_dir)
    _record(study_dir, "4-locus", [core_locus.LOCUS_FILE, core_locus.TRANSCRIPTS_FILE,
                                   "figures/fig1_genomic_architecture.html"], ("locus.",))


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
    study_dir = Path(args.study_dir)
    run_cohorts(study_dir)
    _record(study_dir, "8-cohorts", [core_cohorts.REGISTRY_FILE, core_cohorts.RESOLUTION_FILE,
                                     core_cohorts.TISSUE_DOC_FILE], ("cohort.",))


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
        print(f"  review: {q.key}: {q.question}")
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


def report_cmd(args: argparse.Namespace) -> None:
    study_dir = Path(args.study_dir)
    try:
        path = core_report.write_report(study_dir)
    except (core_report.ReportError, core_contrast.ContrastError, core_master.MasterError, core_cohorts.CohortError,
            core_locus.LocusError, models.DependencyOutputError, study.StudyError) as exc:
        _stage_error(exc)
    print(f"Report: wrote {path}")
    _record(study_dir, "13-report", ["report.md"])


WORKBENCH_APP = SKILL_DIR / "assets" / "workbench.py"      # legacy triage workbench (reproduce stage metrics)
DASHBOARD_APP = SKILL_DIR / "assets" / "dashboard.py"      # analysis dashboard (groups, question, filters chosen on the page)
WORKBENCH_DEPS = list(deps.WORKBENCH_DEPS)
DEFAULT_PORT = 2718


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
    print(f"Dashboard: open {url} in a browser (Ctrl+C here stops it).", flush=True)
    code = 1
    try:
        code = subprocess.call(cmd, cwd=study_dir, env=env)
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


def summarize_cmd(args: argparse.Namespace) -> None:
    """Text summary of the dashboard view for the same choices (the agent reads it to explain the dashboard).

    Same computation as assets/dashboard.py (core.analysis); defaults are the dashboard's opening view:
    the study's target and control tissues, any direction, no filter. Writes nothing.
    """
    import dataclasses

    import polars as pl

    from core import analysis as A
    from viz import dashboard_helpers as H

    study_dir = Path(args.study_dir)
    try:
        data = A.load_study_data(study_dir)
        config = study.load_study(study_dir)
    except (core_master.MasterError, study.StudyError) as exc:
        _stage_error(exc)
    reg = data.registry
    cohorts = config.data.get("cohorts", {})
    scorer = args.scorer
    toward = args.question != "either"
    profile = A.Profile(
        scorer=scorer,
        target=args.target or sorted(set(reg.filter(reg["cohort"] == "target")["ontology_curie"].to_list())),
        control=args.control or sorted(set(reg.filter(reg["cohort"] == "control")["ontology_curie"].to_list())),
        target_label=args.target_label or cohorts.get("target", {}).get("label") or "Alvo",
        control_label=args.control_label or cohorts.get("control", {}).get("label") or "Controle",
        assays=args.assay or [], sources=args.source or [], question=args.question,
        min_change_pct=args.min_change_pct, max_control_pct=args.max_control_pct,
        min_diff_sigma=args.min_sigma, ref_floor=args.ref_floor, ref_floor_applies=args.floor_applies,
        domains=args.domain or [], one_alt_per_position=args.one_alt, y_metric=args.y_metric, top_k=None)
    tissues = A.tissue_table(data, scorer)
    labels = {r["ontology_curie"]: r["tissue"] for r in tissues.iter_rows(named=True)}
    desc = H.describe_tissues(data, scorer)
    noise = A.noise_model(data, scorer, profile.assays or None, profile.sources or None)
    # Same order of operations as the dashboard: controls below the expression floor leave before computing
    used = profile
    if profile.ref_floor is not None and profile.ref_floor_applies in ("control", "both"):
        flags = A.expression_flags(data, profile)
        used = dataclasses.replace(profile, control=[c for c in profile.control if flags.get(c) != "below floor"])
    effects = A.compute(data, used, noise)
    kept, steps = A.apply(effects, data, used, noise)
    what = A.SCALES.get(scorer, (None, "o sinal previsto"))[1]
    out = [f"# {config.gene_symbol}: resumo da visão do dashboard ({scorer}; % d{what})", "",
           H.headline(kept, effects, profile, labels, A.ref_by_tissue(data)), "", "## Grupos e medições", ""]
    for role, curies, name in (("Alvo", profile.target, profile.target_label), ("Controle", profile.control, profile.control_label)):
        out.append(f"**{role} ({name})**")
        out += [f"- {labels.get(c, c)} ({c}): {desc.get(c, 'sem medições deste tipo')}" for c in curies] or ["- nenhum"]
        out.append("")
    dropped = [c for c in profile.control if c not in used.control]
    if dropped:
        out += [f"Controles abaixo do piso de expressão ({H._n(profile.ref_floor, 3)}), fora do cálculo: "
                + ", ".join(labels.get(c, c) for c in dropped) + ".", ""]
    out += ["## Ruído entre medições do mesmo tecido", ""]
    overall = A.noise_overall(data, scorer, profile.assays or None, profile.sources or None)
    if overall is None:
        out.append("Não disponível (nenhum tecido com duas medições nas escolhas atuais); o filtro de ruído fica de fora.")
    else:
        out.append(f"Duas medições do mesmo tecido discordam em até {H._n(A.pct(overall[0], scorer), 1)}% (p95) e "
                   f"{H._n(A.pct(overall[1], scorer), 1)}% (p99), juntando todos os tecidos com duas ou mais medições.")
    out += ["", "## Etapas dos filtros", "", "| etapa | variantes |", "|---|---|"]
    out += [f"| {s} | {n} |" for s, n in steps]
    # Order: the chosen metric (as in the dashboard's table); with the default metric and a direction, the
    # headline's order (the change in the direction asked, in the target tissue that changes least)
    col = "target_toward" if (toward and args.y_metric == "target_mean") else args.y_metric
    signed = col == "target_toward"
    top = kept.sort(pl.col(col) if signed else pl.col(col).abs(), descending=True, nulls_last=True).head(args.top)
    eff = [c for c in top.columns if c.startswith("eff_")]
    out += ["", f"## Primeiras {top.height} variantes", "",
            "| variante | região | % alvo | % controle | " + " | ".join(f"% {labels.get(c[4:], c[4:])}" for c in eff) + " | modalidade AVI |",
            "|" + "---|" * (5 + len(eff))]
    fmt = lambda v: "—" if v is None else H._n(A.pct(v, scorer), 1, True)
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
    """True if the variable is set; the value is never stored or printed."""
    if os.environ.get(name):
        return True
    env_file = Path.home() / ".env"
    try:
        for line in env_file.read_text(encoding="utf-8", errors="replace").splitlines():
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
        found = deps.find_skill(name)
        if found is None:
            errors.append(
                f"Missing dependency skill '{name}' (origin: {origin}).\n"
                f"  Install command (not executed): {command}"
            )
        else:
            installed[name] = found

    searched = ", ".join(str(d) for d in deps.licence_dirs(installed))
    for filename, owner in deps.missing_licences(installed):
        errors.append(
            f"Missing licence file '.licenses/{filename}' required by skill '{owner}' (looked in {searched}). "
            f"Follow the Prerequisites section of '{owner}' with the researcher, one notice at a time."
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
    deps.use_utf8()
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("explore", help="Workflow step 0: read-only exploration before the study exists",
                       description="Gene structure, Regulatory Build, candidate windows, OLS terms per tissue text, Atlas "
                                   "tracks per term x scorer and the gene's predicted activity per tissue "
                                   "(RNA_SEQ_ACTIVE), written to OUT/exploration/ with summary.md (Portuguese).")
    p.add_argument("--gene", required=True, help="HGNC symbol")
    p.add_argument("--out", required=True, help="Folder for the exploration (OUT/exploration/); later the study's parent")
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
    p.add_argument("--scorer", default="RNA_SEQ", help="RNA_SEQ (gene expression; default), CAGE, DNASE or ATAC")
    p.add_argument("--target", action="append", help="Target CURIE (repeatable; default: the study's target tissues)")
    p.add_argument("--control", action="append", help="Control CURIE (repeatable; default: the study's control tissues)")
    p.add_argument("--target-label", help="Name of the target group (default: study.yaml)")
    p.add_argument("--control-label", help="Name of the control group (default: study.yaml)")
    p.add_argument("--question", choices=["repression", "increase", "either"], default="either", help="Direction in the target")
    p.add_argument("--min-change-pct", type=float, default=0.0, help="Smallest change in every target tissue, in %% (0 = off)")
    p.add_argument("--max-control-pct", type=float, help="Largest change allowed in the control, in %% (default: no limit)")
    p.add_argument("--min-sigma", type=float, default=0.0, help="Target - control above this many noise sigmas (0 = off)")
    p.add_argument("--assay", action="append", help="Assay titles kept (repeatable; default all)")
    p.add_argument("--source", action="append", help="Data sources kept, e.g. encode, gtex (repeatable; default all)")
    p.add_argument("--domain", action="append", help="Genomic domains kept (repeatable; default all)")
    p.add_argument("--one-alt", action="store_true", help="Keep one alternative per position (the largest effect)")
    p.add_argument("--ref-floor", type=float, help="Reference expression floor (RNA_SEQ REF; default: off)")
    p.add_argument("--floor-applies", choices=["warn", "control", "target", "both"], default="warn",
                   help="What the floor does: warn (default), drop controls below it, skip a target below it, or both")
    p.add_argument("--y-metric", default="target_mean",
                   choices=["target_mean", "target_toward", "control_mean", "diff", "diff_sigma", "avi_phred"],
                   help="Order of the list (the dashboard's 'Altura no gráfico da janela e ordem da tabela')")
    p.add_argument("--top", type=int, default=10, help="Variants listed (default 10)")
    p.set_defaults(func=summarize_cmd)

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
