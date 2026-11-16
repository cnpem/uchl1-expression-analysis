"""Saved views of the dashboard (`analyses/<slug>.yaml` in the study).

A view is the researcher's choices on the dashboard, saved only when they ask, with a name and a note:
the analysis choices (`core.analysis.Profile`: groups, question, filters, order, how many) and how they
are shown (variants that do not pass, colour, point size, highlighted variants). The dashboard saves and
opens views; `expression.py summarize --view`, `figures --view` and `report` read the same file, so what
the researcher saw, what the agent explains and what the report shows come from one record.

Views hold choices among the parameters the skill defines; they never define new parameters. A file
with an unknown field, or a value outside the allowed ones, is refused when it is read.

Library: VIEWS_DIR, View, slug, save_view, load_view, list_views, ViewError.
"""

from __future__ import annotations

import dataclasses
import datetime as dt
import re
import unicodedata
from pathlib import Path
from typing import Any

import yaml

from core import analysis as A

VIEWS_DIR = "analyses"
DISPLAY_DEFAULTS: dict[str, Any] = {"show_mode": "kept", "color_by": "genomic_domain", "point_size": 40, "highlight": []}
SHOW_MODES = ("kept", "all")
COLOR_COLUMNS = ("genomic_domain", "passa", "sentido", "ensembl_transcript_region", "top_feature")
Y_METRICS = ("target_mean", "target_toward", "control_mean", "diff", "diff_sigma", "avi_phred")
FLOOR_APPLIES = ("warn", "control", "target", "both")
_PROFILE_FIELDS = {f.name: f for f in dataclasses.fields(A.Profile)}


class ViewError(Exception):
    """A view that cannot be saved or read (bad name, unknown field, value out of range)."""


@dataclasses.dataclass
class View:
    name: str
    slug: str
    note: str
    saved_at: str
    created_at: str
    skill_version: str
    created_by: str
    profile: A.Profile
    display: dict[str, Any]
    path: Path


def _now() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def slug(name: str) -> str:
    """File name of a view: lower-case ASCII letters, digits and hyphens ("Repressão no timo" -> "repressao-no-timo")."""
    text = unicodedata.normalize("NFKD", name).encode("ascii", "ignore").decode().lower()
    out = re.sub(r"[^a-z0-9]+", "-", text).strip("-")[:60].strip("-")
    if not out:
        raise ViewError(f"the name {name!r} has no letters or digits; give the view a name")
    return out


def _skill_version() -> str:
    try:
        import manifest

        return manifest._skill_identity()[1]
    except Exception:  # the dashboard runs without the manifest module on some paths
        return "unknown"


def _check_profile(raw: Any, where: str) -> A.Profile:
    if not isinstance(raw, dict):
        raise ViewError(f"{where}: 'profile' must be a mapping")
    unknown = sorted(set(raw) - set(_PROFILE_FIELDS))
    if unknown:
        raise ViewError(f"{where}: unknown analysis choice(s) {unknown}; a view can only hold the skill's own choices")
    p = A.Profile(**raw)
    checks = [
        (p.question in A.QUESTIONS, f"question must be one of {list(A.QUESTIONS)}"),
        (p.y_metric in Y_METRICS, f"y_metric must be one of {list(Y_METRICS)}"),
        (p.color_by in COLOR_COLUMNS, f"color_by must be one of {list(COLOR_COLUMNS)}"),
        (p.ref_floor_applies in FLOOR_APPLIES, f"ref_floor_applies must be one of {list(FLOOR_APPLIES)}"),
        (isinstance(p.scorer, str) and p.scorer, "scorer must be a prediction type such as RNA_SEQ"),
        (all(isinstance(x, list) and all(isinstance(v, str) for v in x)
             for x in (p.target, p.control, p.shown, p.assays, p.sources, p.domains)),
         "target, control, shown, assays, sources and domains must be lists of text"),
        (isinstance(p.min_change_pct, (int, float)) and p.min_change_pct >= 0, "min_change_pct must be >= 0"),
        (p.max_control_pct is None or (isinstance(p.max_control_pct, (int, float)) and p.max_control_pct >= 0),
         "max_control_pct must be empty or >= 0"),
        (isinstance(p.min_diff_sigma, (int, float)) and p.min_diff_sigma >= 0, "min_diff_sigma must be >= 0"),
        (p.ref_floor is None or (isinstance(p.ref_floor, (int, float)) and p.ref_floor >= 0), "ref_floor must be empty or >= 0"),
        (p.top_k is None or (isinstance(p.top_k, int) and p.top_k >= 0), "top_k must be empty or a whole number >= 0"),
        (isinstance(p.one_alt_per_position, bool), "one_alt_per_position must be true or false"),
        (not set(p.target) & set(p.control), "a tissue cannot be in both the target and the control group"),
    ]
    for ok, msg in checks:
        if not ok:
            raise ViewError(f"{where}: {msg}")
    return p


def _check_display(raw: Any, where: str) -> dict[str, Any]:
    raw = raw or {}
    if not isinstance(raw, dict):
        raise ViewError(f"{where}: 'display' must be a mapping")
    unknown = sorted(set(raw) - set(DISPLAY_DEFAULTS))
    if unknown:
        raise ViewError(f"{where}: unknown display setting(s) {unknown}")
    d = {**DISPLAY_DEFAULTS, **raw}
    if d["show_mode"] not in SHOW_MODES:
        raise ViewError(f"{where}: show_mode must be one of {list(SHOW_MODES)}")
    if d["color_by"] not in COLOR_COLUMNS:
        raise ViewError(f"{where}: color_by must be one of {list(COLOR_COLUMNS)}")
    if not isinstance(d["point_size"], int) or not 1 <= d["point_size"] <= 500:
        raise ViewError(f"{where}: point_size must be a whole number between 1 and 500")
    if not isinstance(d["highlight"], list) or not all(isinstance(v, str) for v in d["highlight"]):
        raise ViewError(f"{where}: highlight must be a list of variant IDs")
    return d


CREATED_BY = ("researcher", "agent")


def save_view(study_dir: Path, name: str, profile: A.Profile, display: dict[str, Any] | None = None,
              note: str = "", created_by: str = "researcher") -> View:
    """Writes analyses/<slug>.yaml (replacing a view of the same name, keeping its creation date).

    `created_by` is "researcher" (saved on the dashboard) or "agent" (saved by the agent with
    `expression.py views --save`, after the researcher agreed); the dashboard and the report show it."""
    if created_by not in CREATED_BY:
        raise ViewError(f"created_by must be one of {list(CREATED_BY)}")
    name = name.strip()
    s = slug(name)
    path = Path(study_dir) / VIEWS_DIR / f"{s}.yaml"
    raw_profile = dataclasses.asdict(profile)
    _check_profile(raw_profile, path.name)
    disp = _check_display({k: v for k, v in (display or {}).items() if k in DISPLAY_DEFAULTS}, path.name)
    created = _now()
    if path.is_file():
        try:
            created = (yaml.safe_load(path.read_text(encoding="utf-8")) or {}).get("created_at") or created
        except yaml.YAMLError:
            pass
    now = _now()
    body = {"name": name, "note": note.strip(), "created_at": created, "saved_at": now,
            "skill_version": _skill_version(), "created_by": created_by, "profile": raw_profile, "display": disp}
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("# Saved view of the expression-analysis dashboard (choices only; see references/formats.md)\n"
                    + yaml.safe_dump(body, sort_keys=False, allow_unicode=True, width=1000), encoding="utf-8")
    return load_view(study_dir, s)


def load_view(study_dir: Path, name: str) -> View:
    """Reads a view by its name or file name; refuses unknown fields and values out of range."""
    path = Path(study_dir) / VIEWS_DIR / f"{slug(name)}.yaml"
    if not path.is_file():
        known = [v.name for v in list_views(study_dir)]
        raise ViewError(f"no saved view {name!r} in {path.parent}" + (f"; saved views: {known}" if known else ""))
    try:
        raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except yaml.YAMLError as exc:
        raise ViewError(f"{path.name}: invalid YAML ({exc})") from None
    allowed = {"name", "note", "created_at", "saved_at", "skill_version", "created_by", "profile", "display"}
    unknown = sorted(set(raw) - allowed)
    if unknown:
        raise ViewError(f"{path.name}: unknown field(s) {unknown}")
    return View(name=str(raw.get("name") or path.stem), slug=path.stem, note=str(raw.get("note") or ""),
                saved_at=str(raw.get("saved_at") or ""), created_at=str(raw.get("created_at") or ""),
                skill_version=str(raw.get("skill_version") or ""),
                created_by=str(raw.get("created_by") or "researcher"),
                profile=_check_profile(raw.get("profile"), path.name),
                display=_check_display(raw.get("display"), path.name), path=path)


def list_views(study_dir: Path) -> list[View]:
    """The study's saved views, by name; a file that cannot be read is skipped (load_view names the problem)."""
    out = []
    for path in sorted((Path(study_dir) / VIEWS_DIR).glob("*.yaml")):
        try:
            out.append(load_view(study_dir, path.stem))
        except ViewError:
            continue
    return sorted(out, key=lambda v: v.name.lower())
