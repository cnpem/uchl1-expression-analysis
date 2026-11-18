"""Dependency captures of the Workflow, run and recorded by the skill (`expression.py fetch ...`, `explore`).

Each function runs one dependency skill script with its pinned packages, writes the output under the
exact name the next stage reads (raw/gtf.json, raw/ols_search_<slug>.json, raw/ols_term_<PREFIX>_<digits>.json,
raw/tracks_<S>_<PREFIX>_<digits>.parquet, ...), and appends the command it ran, the package versions and the
`expression.py` call that ran it to the folder's manifest.json. The agent chooses the parameters; nothing is
typed into a shell twice, so the same call works in PowerShell, bash and zsh and the manifest records what
actually ran.

A Workspace is a study directory, or the exploration folder (`<out>/exploration/`, no study.yaml; its
manifest carries its own id). With `dry_run` the commands are printed and nothing runs or is recorded.

Library: Workspace, FetchError, slug, curie_file, gtf, regulatory, sequence, ols_search, ols_term, tracks,
catalogue, avi_region, scores, atlas_tracks, baseline, atlas_links, interpret.
"""

from __future__ import annotations

import dataclasses
import datetime
import json
import re
import sys
from pathlib import Path

import deps
import manifest

AVI = ("alphagenome-variant-impact-score", "scripts/alphagenome_atlas_avi.py")
ATLAS_QUERY = ("alphagenome", "scripts/atlas_query.py")
OLS_SEARCH = ("embl-ebi-ols", "scripts/search_ols.py")
OLS_TERM = ("embl-ebi-ols", "scripts/get_term.py")
ENSEMBL = ("ensembl-database", "scripts/ensembl_api.py")
LINKS = ("alphagenome-atlas-website-links", "scripts/alphagenome_atlas_links.py")
SVA_INTERPRET = ("alphagenome-single-variant-analysis", "scripts/interpret_splicing.py")
SVA_VISUALIZE = ("alphagenome-single-variant-analysis", "scripts/visualize_variant_effects.py")
REGULATORY_SCRIPT = deps.SKILL_DIR / "scripts" / "regulatory_features.py"
GENE_CENTRIC = ("RNA_SEQ", "RNA_SEQ_ACTIVE", "POLYADENYLATION", "SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS")
DEFAULT_OLS_ONTOLOGIES = "uberon,cl,efo,clo"


class FetchError(Exception):
    """A capture failed or cannot start; the CLI prints it and exits with 1."""


def slug(text: str) -> str:
    """File-name slug of a search text: lower case, runs of other characters as one '_'."""
    return re.sub(r"[^a-z0-9]+", "_", text.strip().lower()).strip("_") or "query"


def curie_file(curie: str) -> str:
    return curie.replace(":", "_", 1)


def _today() -> str:
    return datetime.date.today().isoformat()


@dataclasses.dataclass
class Workspace:
    """Where captures go and are recorded.

    Attributes:
        root: Study directory, or the exploration folder.
        study_id: Manifest id for a folder without study.yaml (exploration); None reads study.yaml.
        invoked_by: The `expression.py` command line that asked for the capture (recorded with it).
        dry_run: Print the commands only.
    """

    root: Path
    study_id: str | None = None
    invoked_by: str = ""
    dry_run: bool = False

    def path(self, relative: str) -> Path:
        return self.root / relative

    def record(self, step: str, skill: str, command: str, outputs: list[str], versions: dict[str, str],
               decision_keys: list[str] | None = None) -> None:
        if self.dry_run:
            return
        versions = {**versions, "uv": deps.uv_version()}
        if self.invoked_by:
            versions["invoked_by"] = self.invoked_by
        manifest.add_step(self.root, manifest.ManifestStep(
            step=step, skill=skill, command=command, outputs=outputs, source_versions=versions,
            recorded_at=manifest._now(), decision_keys=list(decision_keys or [])), study_id=self.study_id)
        print(f"Recorded step {step!r} in {self.root / manifest.MANIFEST_NAME} ({', '.join(outputs) or 'no output'})")


def _rel(ws: Workspace, path: Path) -> str:
    return path.relative_to(ws.root).as_posix()


def _avi_argv(*args: str) -> list[str]:
    return ["uv", "run", "--no-project", "--python", deps.DEP_PYTHON, *deps.uv_with(deps.AVI_DEPS),
            str(deps.skill_script(*AVI)), *args]


def _atlas_query_argv(*args: str) -> list[str]:
    env = deps.env_file()
    pre = ["--env-file", env.as_posix()] if env else []  # uv drops backslashes from this path on Windows
    if not deps.key_present("ALPHAGENOME_API_KEY"):
        raise FetchError("ALPHAGENOME_API_KEY is not set in the project's .env; the researcher sets it with the "
                         "alphagenome-skills-setup skill (expression.py check names what is missing).")
    return ["uv", "run", *pre, *deps.uv_with(deps.ATLAS_QUERY_PINS), str(deps.skill_script(*ATLAS_QUERY)), *args]


def _env_args() -> list[str]:
    """`--env-file <the project's .env>` for uv: the key reaches the scripts from the project's .env only (K-Dense
    reads the environment; DeepMind would otherwise look in ~/.env); the value never appears in the command."""
    env = deps.env_file()
    return ["--env-file", env.as_posix()] if env else []  # uv drops backslashes from this path on Windows


def _polite_argv(script: tuple[str, str], *args: str) -> list[str]:
    return ["uv", "run", *deps.uv_with(deps.POLITE_HTTP_PINS), str(deps.skill_script(*script)), *args]


def _exec(ws: Workspace, argv: list[str], *, polite: bool = False) -> deps.DependencyRun | None:
    if ws.dry_run:
        print(deps.display_command(argv))
        return None
    run = deps.run(argv, env_extra=deps.polite_http_env() if polite else None)
    if run.returncode != 0:
        tail = (run.stderr.strip() or run.stdout.strip())[-1500:]
        raise FetchError(f"command failed with exit {run.returncode}: {deps.display_command(argv)}\n{tail}")
    return run


def _pins(packages: tuple[str, ...]) -> str:
    return ", ".join(packages)


def _out(ws: Workspace, relative: str) -> Path:
    path = ws.path(relative)
    if not ws.dry_run:
        path.parent.mkdir(parents=True, exist_ok=True)
    return path


# --------------------------------------------------------------------------- step 2


def gtf(ws: Workspace, gene: str) -> Path:
    """GENCODE v46 structure of `gene` -> raw/gtf.json (the `Loading GTF from ...` line before the JSON is dropped)."""
    out = _out(ws, "raw/gtf.json")
    argv = _avi_argv("gtf", "--gene", gene, "--all_transcripts", "--exons", "--cds", "--utr", "--format", "json")
    run = _exec(ws, argv)
    if run is None:
        return out
    lines = run.stdout.splitlines()
    start = next((i for i, line in enumerate(lines) if line.strip() in ("[", "[]")), None)
    if start is None:
        raise FetchError(f"gtf printed no JSON list for {gene!r}; output began with: {run.stdout[:300]!r}")
    text = "\n".join(lines[start:]) + "\n"
    records = json.loads(text)
    if not records:
        raise FetchError(f"gtf found no transcript for {gene!r} in GENCODE v46; check the HGNC symbol.")
    out.write_text(text, encoding="utf-8")
    ws.record("2-gtf", AVI[0], deps.display_command(argv) + f" > {_rel(ws, out)}  # stdout; lines before the JSON "
              f"list ({start}) dropped by expression.py fetch", [_rel(ws, out)],
              {"gencode": "v46", "python": deps.DEP_PYTHON, "packages": _pins(deps.AVI_DEPS)})
    return out


# --------------------------------------------------------------------------- step 3


def regulatory(ws: Workspace, region: str) -> Path:
    out = _out(ws, "raw/regulatory.json")
    argv = ["uv", "run", str(REGULATORY_SCRIPT), "overlap", "--region", region, "--output", str(out)]
    run = _exec(ws, argv)
    if run is None:
        return out
    m = re.search(r"Ensembl release (\S+)\)", run.stdout)
    print(run.stdout.strip())
    ws.record("3-regulatory", "expression-analysis/regulatory_features.py", deps.display_command(argv), [_rel(ws, out)],
              {"ensembl": m.group(1) if m else "unknown", "assembly": "GRCh38",
               "region_rule": "gene span of raw/gtf.json (normalised) +/- regulatory_search_upstream_bp"})
    return out


def sequence(ws: Workspace, region: str) -> Path:
    out = _out(ws, f"raw/sequence_{region.replace(':', '_').replace('-', '_')}.txt")
    argv = _polite_argv(ENSEMBL, "get-sequence", region, "--species", "human", "--output", str(out))
    if _exec(ws, argv, polite=True) is None:
        return out
    if not out.is_file():
        raise FetchError(f"get-sequence wrote no {out}")
    ws.record("3-sequence", ENSEMBL[0], deps.display_command(argv), [_rel(ws, out)],
              {"assembly": "GRCh38", "packages": _pins(deps.POLITE_HTTP_PINS)})
    return out


# --------------------------------------------------------------------------- steps 5-6


def ols_search(ws: Workspace, query: str, ontology: str = DEFAULT_OLS_ONTOLOGIES, rows: int = 10) -> Path:
    """OLS search for one tissue or cell text -> raw/ols_search_<slug>.json (the name cohorts reads)."""
    out = _out(ws, f"raw/ols_search_{slug(query)}.json")
    argv = _polite_argv(OLS_SEARCH, "--query", query, "--ontology", ontology, "--defining", "--rows", str(rows),
                        "--output", str(out))
    if _exec(ws, argv, polite=True) is None:
        return out
    data = json.loads(out.read_text(encoding="utf-8"))
    if data.get("status") != "success":
        raise FetchError(f"OLS search {query!r} failed: {data.get('message')}")
    ws.record("5-ols-search", OLS_SEARCH[0], deps.display_command(argv), [_rel(ws, out)],
              {"ols": "OLS4", "packages": _pins(deps.POLITE_HTTP_PINS)})
    return out


def ols_term(ws: Workspace, curie: str) -> Path:
    """OLS term with children and all descendants -> raw/ols_term_<PREFIX>_<digits>.json."""
    out = _out(ws, f"raw/ols_term_{curie_file(curie)}.json")
    argv = _polite_argv(OLS_TERM, "--obo_id", curie, "--relations", "hierarchicalChildren,hierarchicalDescendants",
                        "--output", str(out))
    if _exec(ws, argv, polite=True) is None:
        return out
    ws.record("6-ols-term", OLS_TERM[0], deps.display_command(argv), [_rel(ws, out)],
              {"ols": "OLS4", "packages": _pins(deps.POLITE_HTTP_PINS)})
    return out


# --------------------------------------------------------------------------- step 7


def tracks(ws: Workspace, curie: str, scorer: str) -> Path | None:
    """Atlas tracks of one CURIE x scorer -> raw/tracks_<S>_<PREFIX>_<digits>.parquet; a zero-track pair is
    recorded without output (that is how cohorts tells it from a pair never queried)."""
    out = _out(ws, f"raw/tracks_{scorer}_{curie_file(curie)}.parquet")
    argv = _avi_argv("metadata", "--tracks", "--ontology", curie, "--scorer", scorer, "--output", str(out))
    run = _exec(ws, argv)
    if run is None:
        return out
    versions = {"atlas_pull": _today(), "python": deps.DEP_PYTHON, "packages": _pins(deps.AVI_DEPS)}
    if out.is_file():
        ws.record("7-tracks", AVI[0], deps.display_command(argv), [_rel(ws, out)], versions)
        return out
    if "No tracks matched" not in run.stdout:
        raise FetchError(f"metadata --tracks wrote no {out.name} and did not report 'No tracks matched': {run.stdout[-500:]}")
    print(f"No {scorer} track for {curie}; recorded as a zero-track pair.")
    ws.record("7-tracks", AVI[0], deps.display_command(argv), [], {**versions, "result": "no_tracks"})
    return None


def catalogue(ws: Workspace, scorer: str) -> Path:
    out = _out(ws, f"raw/tracks_{scorer}_catalogue.parquet")
    argv = _avi_argv("metadata", "--tracks", "--scorer", scorer, "--output", str(out))
    if _exec(ws, argv) is None:
        return out
    ws.record("7-tracks-catalogue", AVI[0], deps.display_command(argv), [_rel(ws, out)],
              {"atlas_pull": _today(), "python": deps.DEP_PYTHON, "packages": _pins(deps.AVI_DEPS)})
    return out


# --------------------------------------------------------------------------- steps 9-10


def avi_region(ws: Workspace, region: str, max_window: int) -> Path:
    out = _out(ws, "raw/avi_region.parquet")
    argv = _avi_argv("region", "--region", region, "--top_k", "20", "--min_phred", "0", "--max_window_size",
                     str(max_window), "--include_features", "--include_track_info", "--output", str(out))
    if _exec(ws, argv) is None:
        return out
    ws.record("9-avi-region", AVI[0], deps.display_command(argv), [_rel(ws, out)],
              {"atlas_pull": _today(), "python": deps.DEP_PYTHON, "packages": _pins(deps.AVI_DEPS)})
    return out


def scores(ws: Workspace, region: str, max_window: int, scorer: str, curies: list[str], gene: str | None) -> Path:
    """Per-track scores of every SNV of the window for one scorer -> raw/scores_<S>.tsv (`--gene` only for
    gene-centric scorers, which are the only ones that carry a gene)."""
    out = _out(ws, f"raw/scores_{scorer}.tsv")
    gene_args = ["--gene", gene] if gene and scorer in GENE_CENTRIC else []
    argv = _atlas_query_argv("scores", "--interval", region, "--max-window", str(max_window), "--scorers", scorer,
                             "--ontology", *curies, *gene_args, "-o", str(out))
    if _exec(ws, argv) is None:
        return out
    ws.record("10-scores", ATLAS_QUERY[0], deps.display_command(argv), [_rel(ws, out)],
              {"atlas_pull": _today(), "packages": _pins(deps.ATLAS_QUERY_PINS)})
    return out


# --------------------------------------------------------------------------- exploration only


def atlas_tracks(ws: Workspace) -> Path:
    """The whole Atlas track catalogue with data source and assay (K-Dense `tracks`) -> atlas/atlas_tracks.tsv."""
    out = _out(ws, "atlas/atlas_tracks.tsv")
    argv = _atlas_query_argv("tracks", "-o", str(out))
    if _exec(ws, argv) is None:
        return out
    ws.record("0-atlas-catalogue", ATLAS_QUERY[0], deps.display_command(argv), [_rel(ws, out)],
              {"atlas_pull": _today(), "packages": _pins(deps.ATLAS_QUERY_PINS)})
    return out


def baseline(ws: Workspace, variant: str, gene: str) -> Path:
    """RNA_SEQ_ACTIVE of one SNV inside the gene, every track: the predicted activity of the gene per tissue."""
    out = _out(ws, "atlas/scores_RNA_SEQ_ACTIVE_baseline.tsv")
    argv = _atlas_query_argv("scores", "--variant", variant, "--scorers", "RNA_SEQ_ACTIVE", "--gene", gene, "-o", str(out))
    if _exec(ws, argv) is None:
        return out
    ws.record("0-baseline", ATLAS_QUERY[0], deps.display_command(argv), [_rel(ws, out)],
              {"atlas_pull": _today(), "packages": _pins(deps.ATLAS_QUERY_PINS)})
    return out


# --------------------------------------------------------------------------- step 13


def atlas_links(ws: Workspace) -> Path:
    src = ws.path("processed/top_variants.json")
    if not src.is_file() and not ws.dry_run:
        raise FetchError("processed/top_variants.json not found; run expression.py contrast first (Workflow step 12).")
    out = _out(ws, "raw/atlas_links.md")
    argv = ["uv", "run", "--no-project", "--python", deps.DEP_PYTHON, *deps.uv_with(deps.LINKS_DEPS),
            str(deps.skill_script(*LINKS)), "table", "--input", str(src)]
    run = _exec(ws, argv)
    if run is None:
        return out
    out.write_text(run.stdout, encoding="utf-8")
    ws.record("13-atlas-links", LINKS[0], deps.display_command(argv) + f" > {_rel(ws, out)}  # stdout",
              [_rel(ws, out)], {"python": deps.DEP_PYTHON, "packages": _pins(deps.LINKS_DEPS)})
    return out


_VARIANT_RE = re.compile(r"^(chr[0-9A-Za-z]+):(\d+):([ACGT]+)>([ACGT]+)$")


def interpret(ws: Workspace, variant: str, ontology: str, gene: str, tissue: str) -> list[Path]:
    """Model-API interpretation (text) and REF/ALT figures of one variant (alphagenome-single-variant-analysis)."""
    m = _VARIANT_RE.match(variant)
    if m is None:
        raise FetchError(f"variant {variant!r}: expected chr:pos:REF>ALT")
    chrom, pos, ref, alt = m.groups()
    tag = f"{chrom}_{pos}_{ref}_{alt}"
    out = _out(ws, f"raw/interpret_{tag}.txt")
    common = ["--chrom", chrom, "--pos", pos, "--ref", ref, "--alt", alt]
    argv = ["uv", "run", *_env_args(), "--no-project", "--python", deps.DEP_PYTHON, *deps.uv_with(deps.SVA_DEPS),
            str(deps.skill_script(*SVA_INTERPRET)), *common, "--ontology_id", ontology]
    written: list[Path] = []
    run = _exec(ws, argv)
    if run is not None:
        out.write_text(run.stdout, encoding="utf-8")
        ws.record("13-interpret", SVA_INTERPRET[0], deps.display_command(argv) + f" > {_rel(ws, out)}  # stdout",
                  [_rel(ws, out)], {"python": deps.DEP_PYTHON, "packages": _pins(deps.SVA_DEPS)})
        written.append(out)
    vdir = _out(ws, f"figures/variants/{tag}")
    argv = ["uv", "run", *_env_args(), "--no-project", "--python", deps.DEP_PYTHON, *deps.uv_with(deps.SVA_DEPS),
            str(deps.skill_script(*SVA_VISUALIZE)), *common, "--gene", gene, "--tissue", tissue, "--ontology",
            ontology, "--tracks", "regulatory", "--output_dir", str(vdir)]
    if _exec(ws, argv) is not None:
        files = sorted(p for p in vdir.rglob("*") if p.is_file())
        ws.record("13-visualize", SVA_VISUALIZE[0], deps.display_command(argv), [_rel(ws, p) for p in files],
                  {"python": deps.DEP_PYTHON, "packages": _pins(deps.SVA_DEPS)})
        written.extend(files)
    return written


def invocation() -> str:
    """This process's command line, as `uv run <script> <args>` (recorded with every capture)."""
    return deps.display_command(["uv", "run", sys.argv[0], *sys.argv[1:]])
