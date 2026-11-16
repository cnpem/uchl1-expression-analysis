"""Dependency skills, pinned runtimes and the environment around them (library, no CLI).

Everything a fetch or a check needs to find and run a dependency skill the same way on Windows
(PowerShell), macOS and Linux (bash, zsh): where the skills are installed, which licence notices
they ask for and where those notices may live, the exact package versions each dependency script
runs with, where `~/.env` is, and how its command is written into manifest.json.

Pinned versions: resolved with `uv pip compile` on 2026-11-15 and tested on Windows ARM64 (x86-64
Python under emulation) and with the bundled reference. polars stays at 1.31.0 because newer
releases read the CPU features through `platform.machine()`, which reports ARM64 under emulation,
and then refuse to start (`unknown feature flag: 'sse3'`). pysam has no Windows wheel; the
DeepMind pyproject marks it `sys_platform != 'win32'` and the AVI script imports it only for VCF
output, which this skill never asks for.

Library: DEPENDENCIES, search_dirs, find_skill, skill_script, licence_dirs, missing_licences,
env_file, use_utf8, polite_http_env, DependencyRun, run, display_command, uv_version.
"""

from __future__ import annotations

import dataclasses
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent.parent
SKILL_NAME = SKILL_DIR.name

# Python used for the pyproject-based DeepMind scripts (their pyproject asks for >= 3.14, where
# h5py has no wheel), and the exact packages each dependency script runs with.
DEP_PYTHON = "3.13"
AVI_DEPS = ("alphagenome==2.0.0", "absl-py==2.5.1", "anndata==0.13.4", "numpy==2.5.4", "pandas==2.3.3",
            "polars==1.31.0", "pyarrow==26.0.1", "pysam==0.24.1; sys_platform != 'win32'",
            "python-dotenv==1.2.4", "tabulate==0.10.0")
LINKS_DEPS = ("alphagenome==2.0.0", "numpy==2.5.4", "pandas==2.3.3", "pyarrow==26.0.1", "python-dotenv==1.2.4")
# alphagenome 2.0.0 removed dna_client, which the model-API scripts import.
SVA_DEPS = ("alphagenome==0.11.0", "absl-py==2.5.1", "python-dotenv==1.2.4", "matplotlib==3.11.2",
            "numpy==2.5.4", "pandas==3.1.1", "pyarrow==26.0.1", "scipy==1.18.1", "seaborn==0.13.2")
# Scripts with their own PEP 723 header get the same packages pinned with --with.
ATLAS_QUERY_PINS = ("alphagenome==2.0.0",)
POLITE_HTTP_PINS = ("polite-http==0.1.4",)
WORKBENCH_DEPS = ("marimo==0.25.2", "altair==6.3.0", "anywidget==0.11.0", "pandas==3.1.1", "polars==1.31.0",
                  "pyarrow==26.0.1", "pyyaml==6.0.3", "matplotlib==3.11.2",
                  "vl-convert-python==1.9.0.post1")

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
_LICENSE_RE = re.compile(r"\.licenses/([A-Za-z0-9_.-]+_LICENSE\.txt)")


class DependencyError(Exception):
    """A dependency skill, script or input is missing; the CLI prints it and exits with 1."""


def use_utf8() -> None:
    """Writes stdout/stderr as UTF-8, so labels with accents or symbols never stop a run on a cp1252 console."""
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, ValueError):
            pass


# --------------------------------------------------------------------------- locations


def search_dirs() -> list[Path]:
    """Directories searched for installed skills, in priority order."""
    cwd, home = Path.cwd(), Path.home()
    return [
        SKILL_DIR.parent,
        cwd / ".claude" / "skills",
        cwd / ".agents" / "skills",
        home / ".claude" / "skills",
        home / ".agents" / "skills",
    ]


def find_skill(name: str) -> Path | None:
    """The folder of an installed skill (the first search dir holding <name>/SKILL.md), or None."""
    for d in search_dirs():
        if (d / name / "SKILL.md").is_file():
            return d / name
    return None


def skill_script(name: str, relative: str) -> Path:
    """Path of a dependency script; raises DependencyError when the skill or script is missing."""
    folder = find_skill(name)
    if folder is None:
        raise DependencyError(f"Dependency skill '{name}' not found in {', '.join(str(d) for d in search_dirs())}. "
                              "Run expression.py check and install it with the researcher's approval.")
    path = folder / relative
    if not path.is_file():
        raise DependencyError(f"{path} not found: the installed '{name}' skill has no {relative}.")
    return path


def _workspace_of(skill_folder: Path) -> Path | None:
    """The workspace that holds `<ws>/.claude/skills/<name>` or `<ws>/.agents/skills/<name>`."""
    skills = skill_folder.parent
    if skills.name == "skills" and skills.parent.name in (".claude", ".agents"):
        return skills.parent.parent
    if skills.name == "skills":  # a repository layout: <repo>/skills/<name>
        return skills.parent
    return None


def licence_dirs(installed: dict[str, Path] | None = None) -> list[Path]:
    """Folders that may hold `.licenses/`: the working directory first, then the workspaces of this skill and of
    each dependency skill found, so `check` gives the same answer from any folder of the workspace."""
    roots = [Path.cwd()]
    for folder in [SKILL_DIR, *(p.parent if p.name == "SKILL.md" else p for p in (installed or {}).values())]:
        ws = _workspace_of(folder)
        if ws is not None:
            roots.append(ws)
    out: list[Path] = []
    for root in roots:
        d = (root / ".licenses").resolve()
        if d not in out:
            out.append(d)
    return out


def missing_licences(installed: dict[str, Path]) -> list[tuple[str, str]]:
    """(notice file, owner skill) of every licence notice a dependency asks for that no licence dir holds."""
    wanted: dict[str, str] = {}
    for name, path in installed.items():
        skill_md = path if path.name == "SKILL.md" else path / "SKILL.md"
        for match in _LICENSE_RE.findall(skill_md.read_text(encoding="utf-8", errors="replace")):
            wanted.setdefault(match, name)
    dirs = licence_dirs(installed)
    return [(f, o) for f, o in sorted(wanted.items()) if not any((d / f).is_file() for d in dirs)]


def env_file() -> Path | None:
    """`~/.env`, expanded here so no shell has to expand `~`; None when it does not exist."""
    path = Path(os.path.expanduser("~")) / ".env"
    return path if path.is_file() else None


def polite_http_env() -> dict[str, str]:
    """Environment for the scripts built on polite-http (OLS, Ensembl).

    polite-http paces requests through a lock file per host holding the `time.monotonic()` of the next
    free slot. A file written by a process with another monotonic clock (another Python build, or before a
    reboot) can hold a slot hours ahead, and every later call then sleeps silently. The skill keeps its own
    lock folder and resets a reservation that lies more than ten minutes ahead of this clock.
    """
    lock_dir = Path.home() / ".cache" / "alphagenome-gene-skills" / "polite-http"
    lock_dir.mkdir(parents=True, exist_ok=True)
    now = time.monotonic()
    for lock in lock_dir.glob("polite-http-*.lock"):
        try:
            slot = float(lock.read_text(encoding="ascii").strip() or 0)
        except (OSError, ValueError):
            continue
        if slot > now + 600:
            lock.write_text("0", encoding="ascii")
    return {"POLITE_HTTP_LOCK_DIR": str(lock_dir)}


# --------------------------------------------------------------------------- running


def uv_version() -> str:
    try:
        out = subprocess.run(["uv", "--version"], capture_output=True, text=True, check=False).stdout.split()
        return out[1] if len(out) > 1 else "unknown"
    except OSError:
        return "unknown"


_SAFE = re.compile(r"^[\w@%+=:,./-]+$")  # a backslash (Windows path) is quoted: bash would eat it bare


def _quote(arg: str) -> str:
    if _SAFE.match(arg):
        return arg
    if not re.search(r'["$`]|\\$', arg):
        return f'"{arg}"'  # double quotes mean the same in bash, zsh and PowerShell for these characters
    return shlex.quote(arg)


def display_command(argv: list[str]) -> str:
    """The command as one line that runs unchanged in bash, zsh and PowerShell (quoted only where needed)."""
    return " ".join(_quote(str(a)) for a in argv)


def uv_with(packages: tuple[str, ...]) -> list[str]:
    return [a for p in packages for a in ("--with", p)]


@dataclasses.dataclass(frozen=True)
class DependencyRun:
    """A finished dependency command.

    Attributes:
        argv: The command as run.
        returncode: Its exit code.
        stdout, stderr: Its output, decoded as UTF-8.
        seconds: Wall time.
    """

    argv: list[str]
    returncode: int
    stdout: str
    stderr: str
    seconds: float


def run(argv: list[str], *, env_extra: dict[str, str] | None = None, cwd: Path | None = None,
        echo: bool = True) -> DependencyRun:
    """Runs a dependency command, echoing its stderr, and returns its output (stdout kept for the caller)."""
    if shutil.which(argv[0]) is None:
        raise DependencyError(f"'{argv[0]}' is not on PATH; install it with the 'uv' skill.")
    env = {**os.environ, "PYTHONIOENCODING": "utf-8", "PYTHONUTF8": "1", **(env_extra or {})}
    if echo:
        print(f"$ {display_command(argv)}", file=sys.stderr, flush=True)
    t0 = time.monotonic()
    proc = subprocess.run(argv, capture_output=True, env=env, cwd=cwd, check=False)
    out = proc.stdout.decode("utf-8", errors="replace")
    err = proc.stderr.decode("utf-8", errors="replace")
    if echo and err.strip():
        print(err.rstrip(), file=sys.stderr)
    return DependencyRun(list(map(str, argv)), proc.returncode, out, err, time.monotonic() - t0)
