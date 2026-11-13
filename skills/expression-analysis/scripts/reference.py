#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
"""Reference reproduction harness: re-run the stages on a bundled reference study and compare.

Used by `expression.py reproduce` and `splicing.py reproduce`. It is a library: the CLIs pass
their registered stages to run_reference().

A stage is a tuple (name, run_fn, metrics_fn):
  run_fn(study_dir: Path) -> None          runs the stage inside the (temporary) study directory,
                                           reading only raw/, decisions.yaml and contrasts.yaml
  metrics_fn(study_dir: Path) -> dict      returns {metric: int | float | str} for that stage

expected.json format:
  {"reference": "name", "sources": {...},
   "stages": {"<stage>": {"<metric>": {"value": V, "tol": T, "headline": true}}}}
`tol` absent means exact equality (integers and text); with `tol`, |observed - value| <= tol.
`headline: true` marks a main result of the reference: after a passing run, the observed value of
each such metric of the stages that ran is printed in a "Headline results" block
(`stage.metric = value`), so the agent can report the numbers and not only "it matches".
"""

import dataclasses
import json
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Any, Callable

Stage = tuple[str, Callable[[Path], None], Callable[[Path], dict]]


@dataclasses.dataclass(frozen=True, kw_only=True)
class Mismatch:
    """One metric whose observed value differs from expected.json.

    Attributes:
        stage: Stage name.
        metric: Metric name ("<stage>" itself when the whole stage has no expected entry).
        expected: Value in expected.json (None when missing there).
        observed: Value the stage produced (None when the stage did not report it).
        tol: Absolute tolerance, or None for exact equality.
    """

    stage: str
    metric: str
    expected: Any
    observed: Any
    tol: float | None


def compare(expected: dict, observed: dict[str, dict]) -> list[Mismatch]:
    """Compares observed stage metrics with a parsed expected.json.

    Only stages present in `observed` are checked (the others were not run). A stage that ran
    but has no entry in expected.json is a mismatch, as is an expected metric the stage did not
    report. Metrics reported but not listed in expected.json are ignored.

    Args:
        expected: Parsed expected.json (uses its "stages" mapping).
        observed: {stage: {metric: value}} from the metrics functions.

    Returns:
        The mismatches, in stage then metric order; empty when everything matches.
    """
    stages = expected.get("stages", {})
    out: list[Mismatch] = []
    for stage, metrics in observed.items():
        if stage not in stages:
            out.append(Mismatch(stage=stage, metric=stage, expected=None, observed="(stage ran)", tol=None))
            continue
        for metric, spec in stages[stage].items():
            want, tol = spec["value"], spec.get("tol")
            got = metrics.get(metric)
            if got is None:
                out.append(Mismatch(stage=stage, metric=metric, expected=want, observed=None, tol=tol))
            elif tol is None or isinstance(want, str) or isinstance(got, str):
                if got != want:
                    out.append(Mismatch(stage=stage, metric=metric, expected=want, observed=got, tol=tol))
            elif not abs(got - want) <= tol:
                out.append(Mismatch(stage=stage, metric=metric, expected=want, observed=got, tol=tol))
    return out


def _format_value(value: Any) -> str:
    if isinstance(value, float) and not isinstance(value, bool):
        return f"{value:.6f}" if value == 0 or 1e-4 <= abs(value) < 1e6 else f"{value:.6g}"
    return str(value)


def headline_lines(expected: dict, observed: dict[str, dict]) -> list[str]:
    """`stage.metric = value` (observed) for every metric expected.json marks `headline: true`.

    Only stages present in `observed` (the ones that ran) are listed, in expected.json order.
    """
    lines = []
    for stage, metrics in expected.get("stages", {}).items():
        if stage not in observed:
            continue
        for metric, spec in metrics.items():
            if isinstance(spec, dict) and spec.get("headline") is True:
                lines.append(f"  {stage}.{metric} = {_format_value(observed[stage].get(metric))}")
    return lines


def _print_mismatches(mismatches: list[Mismatch]) -> None:
    rows = [("stage", "metric", "expected", "observed", "tolerance")]
    rows += [
        (m.stage, m.metric, repr(m.expected), repr(m.observed), "exact" if m.tol is None else str(m.tol))
        for m in mismatches
    ]
    widths = [max(len(r[i]) for r in rows) for i in range(5)]
    for r in rows:
        print("  ".join(c.ljust(w) for c, w in zip(r, widths)).rstrip(), file=sys.stderr)


def run_reference(
    reference_dir: Path,
    stages: list[Stage],
    until: str | None,
    keep: Path | None,
) -> int:
    """Copies a reference study, runs the stages in order up to `until`, and compares.

    Args:
        reference_dir: assets/reference/<name> of the skill (holds expected.json).
        stages: Registered stages, in execution order.
        until: Last stage to run; None runs all of them.
        keep: Directory for the working copy; None uses a temporary directory that is removed.

    Returns:
        Process exit code: 0 when every compared metric matches, 1 otherwise (divergence,
        unknown or unimplemented stage, no stages, missing reference).
    """
    names = [s[0] for s in stages]
    if until is not None and until not in names:
        print(f"Stage '{until}' is not implemented yet.", file=sys.stderr)
        return 1
    if not stages:
        print("No stages are implemented yet; nothing to reproduce.", file=sys.stderr)
        return 1
    expected_path = reference_dir / "expected.json"
    if not expected_path.is_file():
        print(f"Reference {reference_dir} has no expected.json.", file=sys.stderr)
        return 1
    expected = json.loads(expected_path.read_text(encoding="utf-8"))
    to_run = stages if until is None else stages[: names.index(until) + 1]

    if keep is not None:
        if keep.exists() and any(keep.iterdir()):
            print(f"--keep directory {keep} is not empty.", file=sys.stderr)
            return 1
        work = keep
        shutil.copytree(reference_dir, work, dirs_exist_ok=True)
        cleanup = None
    else:
        cleanup = tempfile.mkdtemp(prefix="reproduce-")
        work = Path(cleanup) / reference_dir.name
        shutil.copytree(reference_dir, work)
    try:
        observed: dict[str, dict] = {}
        for name, run_fn, metrics_fn in to_run:
            try:
                run_fn(work)
            except SystemExit as exc:
                print(f"Stage '{name}' failed (exit {exc.code}).", file=sys.stderr)
                return 1
            observed[name] = metrics_fn(work)
        mismatches = compare(expected, observed)
        if mismatches:
            print(f"Reference '{reference_dir.name}' diverges from expected.json:", file=sys.stderr)
            _print_mismatches(mismatches)
            return 1
        print(f"Reference '{reference_dir.name}': {len(to_run)} stage(s) match expected.json ({', '.join(observed)}).")
        headlines = headline_lines(expected, observed)
        if headlines:
            print("Headline results (observed; report these numbers in the answer):")
            for line in headlines:
                print(line)
        if keep is not None:
            print(f"Working copy kept in {work}")
        return 0
    finally:
        if cleanup is not None:
            shutil.rmtree(cleanup, ignore_errors=True)
