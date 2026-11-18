#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["pyyaml==6.0.3"]
# ///
"""Provenance manifest of a study directory (manifest.json).

Subcommands:
  add   Append one step (dependency command, its outputs, source versions, decision keys).
        Safe for concurrent callers (subagents): the read-modify-write runs under a file lock
        (<study-dir>/.manifest.lock) and the file is replaced atomically.
  show  Print the manifest as a table.

Library: ManifestStep and add_step(study_dir, step).

manifest.json format:
  {"study_id", "skill", "skill_version", "origin": {"kind": "live" | "snapshot", "note"},
   "steps": [{"step", "skill", "command", "outputs", "source_versions", "recorded_at",
              "decision_keys"}]}

`skill` and `skill_version` come from the frontmatter (name, metadata.version) of the SKILL.md
of the skill this script lives in; `study_id` comes from <study-dir>/study.yaml.
Output paths are relative to the study directory and must exist when the step is recorded.
A command that legitimately writes nothing (for example a catalogue query with no match) is
recorded with "outputs": [].
"""

import argparse
import contextlib
import dataclasses
import datetime
import json
import os
import re
import sys
import tempfile
from pathlib import Path
from typing import Any

import yaml

try:
    import fcntl
except ImportError:  # Windows: no advisory lock; writes stay atomic
    fcntl = None

SKILL_DIR = Path(__file__).resolve().parent.parent
MANIFEST_NAME = "manifest.json"
LOCK_NAME = ".manifest.lock"
ORIGIN_KINDS = ("live", "snapshot")
# Version of the study's file formats (manifest.json, study.yaml, decisions.yaml, pending_decisions.json), written
# when the study is created. A study with a higher number was made by a newer version of the skill: it is refused
# instead of being read wrongly. Studies made before the field existed are format 1.
FORMAT_VERSION = 1


class ManifestError(Exception):
    """Invalid manifest input; the CLI prints it and exits with 1."""


def _skill_identity() -> tuple[str, str]:
    """Returns (name, metadata.version) from this skill's SKILL.md frontmatter."""
    path = SKILL_DIR / "SKILL.md"
    match = re.match(r"^---\n(.*?)\n---\n", path.read_text(encoding="utf-8"), re.S)
    if not match:
        raise ManifestError(f"{path}: missing YAML frontmatter")
    front = yaml.safe_load(match.group(1)) or {}
    name = front.get("name")
    version = (front.get("metadata") or {}).get("version")
    if not name or not version:
        raise ManifestError(f"{path}: frontmatter needs 'name' and 'metadata.version'")
    return str(name), str(version)


def _study_id(study_dir: Path) -> str:
    path = study_dir / "study.yaml"
    if not path.is_file():
        raise ManifestError(f"{path} not found; create the study first (study.py init).")
    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    study_id = data.get("study_id") if isinstance(data, dict) else None
    if not study_id:
        raise ManifestError(f"{path}: missing 'study_id'")
    return str(study_id)


def _parse_pairs(pairs: list[str]) -> dict[str, str]:
    out: dict[str, str] = {}
    for pair in pairs:
        key, sep, value = pair.partition("=")
        if not sep or not key.strip():
            raise ManifestError(f"--source-version expects key=value, got {pair!r}")
        out[key.strip()] = value.strip()
    return out


def _check_outputs(study_dir: Path, outputs: list[str]) -> list[str]:
    checked = []
    for output in outputs:
        rel = Path(output)
        if rel.is_absolute() or ".." in rel.parts:
            raise ManifestError(f"output {output!r} must be a path relative to the study directory")
        if not (study_dir / rel).is_file():
            raise ManifestError(f"output {output!r} does not exist in {study_dir}")
        checked.append(rel.as_posix())
    return checked


@dataclasses.dataclass(frozen=True, kw_only=True)
class ManifestStep:
    """One recorded dependency command.

    Attributes:
        step: Step label, e.g. "2-gtf".
        skill: Dependency skill that produced the outputs.
        command: Exact command that was run.
        outputs: Paths relative to the study directory.
        source_versions: Versions of the data sources used.
        recorded_at: UTC ISO-8601 time the step was recorded.
        decision_keys: Decision keys (decisions.yaml) the step depended on.
    """

    step: str
    skill: str
    command: str
    outputs: list[str]
    source_versions: dict[str, str]
    recorded_at: str
    decision_keys: list[str]


def _now() -> str:
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


@contextlib.contextmanager
def _locked(study_dir: Path):
    """Exclusive lock on the study's manifest, held across read-modify-write."""
    with open(study_dir / LOCK_NAME, "a+") as handle:
        if fcntl is not None:
            fcntl.flock(handle, fcntl.LOCK_EX)
        try:
            yield
        finally:
            if fcntl is not None:
                fcntl.flock(handle, fcntl.LOCK_UN)


def _atomic_write(path: Path, text: str) -> None:
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=path.name + ".", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(text)
        os.replace(tmp, path)
    except BaseException:
        with contextlib.suppress(OSError):
            os.unlink(tmp)
        raise


def check_format(data: dict, path: Path) -> None:
    """Refuses a file written in a newer format than this skill reads (absent field = format 1)."""
    found = data.get("format_version", 1) if isinstance(data, dict) else 1
    if not isinstance(found, int) or found > FORMAT_VERSION:
        raise ManifestError(f"{path}: format_version {found!r} is newer than this skill reads ({FORMAT_VERSION}); the "
                            "study was made by a newer version of the skill. Update the skill before using this study.")


def _load_or_create(study_dir: Path, origin_kind: str | None, origin_note: str | None,
                    study_id: str | None = None) -> dict:
    path = study_dir / MANIFEST_NAME
    if path.is_file():
        manifest = json.loads(path.read_text(encoding="utf-8"))
        check_format(manifest, path)
        origin = manifest.get("origin", {})
        if origin_kind is not None and origin_kind != origin.get("kind"):
            raise ManifestError(f"{path}: origin is already {origin.get('kind')!r}; it is set only on creation")
        if origin_note is not None and origin_note != origin.get("note"):
            raise ManifestError(f"{path}: origin note is already set; it is set only on creation")
        return manifest
    name, version = _skill_identity()
    return {
        "study_id": study_id or _study_id(study_dir),
        "skill": name,
        "skill_version": version,
        "format_version": FORMAT_VERSION,
        "origin": {"kind": origin_kind or "live", "note": origin_note or ""},
        "steps": [],
    }


def add_step(
    study_dir: Path,
    step: ManifestStep,
    *,
    origin_kind: str | None = None,
    origin_note: str | None = None,
    study_id: str | None = None,
) -> None:
    """Appends a step to <study_dir>/manifest.json, creating the manifest if needed.

    Concurrent callers are serialised by a file lock and the file is replaced atomically, so
    no step is lost and readers never see a partial file. Output paths must already exist.

    Args:
        study_dir: Study directory (holds study.yaml).
        step: The step to record; `outputs` are validated as relative paths that exist.
        origin_kind: "live" or "snapshot"; only honoured when the manifest is created.
        origin_note: Origin note; only honoured when the manifest is created.
        study_id: Identifier for a manifest created outside a study (the exploration folder has no
            study.yaml); by default it is read from study.yaml.

    Raises:
        ManifestError: On a missing study directory, bad output path or origin conflict.
    """
    study_dir = Path(study_dir)
    if not study_dir.is_dir():
        raise ManifestError(f"study directory {study_dir} not found")
    outputs = _check_outputs(study_dir, step.outputs)
    with _locked(study_dir):
        manifest = _load_or_create(study_dir, origin_kind, origin_note, study_id)
        manifest["steps"].append({**dataclasses.asdict(step), "outputs": outputs})
        _atomic_write(study_dir / MANIFEST_NAME, json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")


def add(args: argparse.Namespace) -> None:
    step = ManifestStep(
        step=args.step,
        skill=args.skill,
        command=args.command,
        outputs=list(args.output),
        source_versions=_parse_pairs(args.source_version),
        recorded_at=_now(),
        decision_keys=list(args.decision_key),
    )
    add_step(Path(args.study_dir), step, origin_kind=args.origin, origin_note=args.origin_note)
    print(f"Recorded step {args.step!r} ({len(step.outputs)} output(s)) in {Path(args.study_dir) / MANIFEST_NAME}")


def _clip(text: str, width: int) -> str:
    text = " ".join(str(text).split())
    return text if len(text) <= width else text[: width - 3] + "..."


def show(args: argparse.Namespace) -> None:
    study_dir = Path(args.study_dir)
    path = study_dir / MANIFEST_NAME
    if not path.is_file():
        raise ManifestError(f"{path} not found; no step has been recorded yet (manifest.py add).")
    manifest: dict[str, Any] = json.loads(path.read_text(encoding="utf-8"))
    check_format(manifest, path)
    origin = manifest.get("origin", {})
    print(f"study {manifest.get('study_id')} | {manifest.get('skill')} {manifest.get('skill_version')}"
          f" | origin: {origin.get('kind')} ({_clip(origin.get('note', ''), 70)})")
    rows = [("step", "skill", "outputs", "decisions", "recorded_at", "command")]
    for s in manifest.get("steps", []):
        rows.append((str(s.get("step")), str(s.get("skill")), str(len(s.get("outputs", []))),
                     ",".join(s.get("decision_keys", [])) or "-", str(s.get("recorded_at")),
                     _clip(s.get("command", ""), 60)))
    widths = [max(len(r[i]) for r in rows) for i in range(len(rows[0]))]
    for r in rows:
        print("  ".join(c.ljust(w) for c, w in zip(r, widths)).rstrip())
    print(f"{len(rows) - 1} step(s)")


def main() -> None:
    for stream in (sys.stdout, sys.stderr):  # UTF-8 on any console (cp1252 on Windows)
        stream.reconfigure(encoding="utf-8", errors="replace")
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command_name", required=True)

    p = sub.add_parser("add", help="Append one provenance step to manifest.json")
    p.add_argument("--study-dir", required=True, help="Study directory (holds study.yaml)")
    p.add_argument("--step", required=True, help="Step label, e.g. 2-gtf")
    p.add_argument("--skill", required=True, help="Dependency skill that produced the outputs")
    p.add_argument("--command", required=True, help="Exact command that was run")
    p.add_argument("--output", action="append", default=[], help="Output path relative to the study dir (repeatable)")
    p.add_argument("--source-version", action="append", default=[], metavar="K=V", help="Source version (repeatable)")
    p.add_argument("--decision-key", action="append", default=[], metavar="KEY", help="Decision key used (repeatable)")
    p.add_argument("--origin", choices=ORIGIN_KINDS, help="Origin kind; only when the manifest is created")
    p.add_argument("--origin-note", help="Origin note; only when the manifest is created")
    p.set_defaults(func=add)

    p = sub.add_parser("show", help="Print manifest.json as a table")
    p.add_argument("--study-dir", required=True, help="Study directory")
    p.set_defaults(func=show)

    args = parser.parse_args()
    try:
        args.func(args)
    except (ManifestError, OSError, json.JSONDecodeError, yaml.YAMLError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
