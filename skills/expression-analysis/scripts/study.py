#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["pyyaml==6.0.3"]
# ///
"""Study configuration, researcher decisions and pending questions of a study directory.

The kind of study ("expression" or "splicing") is read from the name in the SKILL.md of the
skill this script lives in, so the same file serves both skills.

Subcommands:
  init      Create the study directory and study.yaml from assets/study_template.yaml, or, with
            --from-exploration DIR, from the options the researcher chose in an `expression.py explore`
            run (window, target and control CURIEs, labels, scorers), copying its captures into raw/
            with their manifest steps (no network call is repeated).
  validate  Check study.yaml against the fields of the skill's kind.
  decide    Record a researcher's answer in decisions.yaml and clear it from the pending list.
  require   Exit 1 with the missing-decision message unless the key has an answer.
  pending   List the open questions (processed/pending_decisions.json).
  window    --print-max-window: max(1000, length of the longest window to score).

Files: study.yaml, decisions.yaml ({key: {value, workflow_step, asked_at,
answered_at, note[, options, recommended]}}), processed/pending_decisions.json (list of {key,
question, options, recommended, workflow_step}). Decision keys look like
<area>.<topic>[:<subject>] (KEY_PATTERN); `decide` refuses any other key unless a stage raised it
verbatim (a pending key) or it was decided before, and warns about a key that is neither pending
nor decided before (unless --option says the agent asked it).
`--key` is repeatable: every key given gets the same value. `decide` copies `options` and
`recommended` of the matching pending entry into the decision, so decisions.yaml shows what was
offered next to what was chosen. For a question the agent asks itself (no pending entry, e.g.
cohort.target), `--option` (repeatable) and `--recommended` record what was offered; they take
precedence over a pending entry. A key decided with neither has no offered fields (unless an
earlier decision of the same key had them). An answer `replace:{curie}` to `cohort.unserved:{old}`
also writes {curie} in place of {old} in study.yaml, with a comment keeping the original.

Library: StudyConfig.from_dict, load_study, saturation_flank_bp, DecisionEntry, load_decisions,
require, PendingDecision, raise_pending, max_window.
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
from typing import Any, Literal, NoReturn

import yaml

try:
    import fcntl
except ImportError:  # Windows: no advisory lock; writes stay atomic
    fcntl = None

SKILL_DIR = Path(__file__).resolve().parent.parent
STUDY_FILE = "study.yaml"
DECISIONS_FILE = "decisions.yaml"
PENDING_FILE = "processed/pending_decisions.json"
LOCUS_FILE = "processed/locus.json"
STUDY_SUBDIRS = ("raw", "processed", "figures")
MIN_WINDOW = 1000
DECISION_FIELDS = ("value", "workflow_step", "asked_at", "answered_at", "note")
OFFERED_FIELDS = ("options", "recommended")  # optional in decisions.yaml, always together
TRACE_FIELDS = ("value_label", "decided_by", "rule", "previous", "said", "delegated")  # optional (absent in older studies)
PREVIOUS_FIELDS = ("value", "value_label", "answered_at", "note", "options", "recommended", "decided_by", "said", "delegated")
STEP1_DELEGABLE = {"window": "study.window", "scorers": "study.scorers", "question": "analysis.question"}
SAID_HELP = ("the researcher's own words for this answer, copied from the conversation (not paraphrased; several "
             "messages may be joined with ' / '). An answer the researcher gave before being asked counts: quote it. "
             "If they told you to choose, quote that and pass --delegated")
DECIDED_BY = ("researcher", "skill")
# Plain names of the prediction types (Atlas scorers) for texts the researcher reads
PREDICTION_NAMES = {"RNA_SEQ": "gene expression (RNA-seq)", "RNA_SEQ_ACTIVE": "gene activity (RNA-seq)",
                    "CAGE": "transcription start (CAGE)", "DNASE": "open chromatin (DNase)",
                    "ATAC": "open chromatin (ATAC)", "SPLICE_SITES": "splice sites", "SPLICE_SITE_USAGE": "splice-site use",
                    "SPLICE_JUNCTIONS": "splice junctions", "POLYADENYLATION": "polyadenylation"}


def prediction_names(scorers: Any) -> str:
    """'gene expression (RNA-seq), transcription start (CAGE) and open chromatin (DNase)'."""
    names = [PREDICTION_NAMES.get(s, s) for s in scorers]
    return names[0] if len(names) == 1 else ", ".join(names[:-1]) + " and " + names[-1] if names else "none"
DEFAULT_SATURATION_FLANK_BP = 100  # splicing: bp on each side of every exon-intron boundary scanned
KEY_PATTERN = re.compile(r"^[a-z_]+\.[a-z0-9_]+(?::\S+)?\Z")  # <area>.<topic>[:<subject>], no whitespace
_REGION_RE = re.compile(r"^\s*([^:\s]+):(\d+)-(\d+)\s*$")

Kind = Literal["expression", "splicing"]
KINDS = {"expression-analysis": "expression", "splicing-analysis": "splicing"}


class StudyError(Exception):
    """Invalid study input; the CLI prints it and exits with 1."""


def _fail(message: str) -> NoReturn:
    print(message, file=sys.stderr)
    sys.exit(1)


def _now() -> str:
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# --- file helpers ---------------------------------------------------------------------------


@contextlib.contextmanager
def _locked(study_dir: Path):
    """Exclusive lock for read-modify-write of decisions and pending files."""
    with open(study_dir / ".study.lock", "a+") as handle:
        if fcntl is not None:
            fcntl.flock(handle, fcntl.LOCK_EX)
        try:
            yield
        finally:
            if fcntl is not None:
                fcntl.flock(handle, fcntl.LOCK_UN)


def _atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=path.name + ".", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(text)
        os.replace(tmp, path)
    except BaseException:
        with contextlib.suppress(OSError):
            os.unlink(tmp)
        raise


def _read_yaml(path: Path) -> Any:
    if not path.is_file():
        raise StudyError(f"{path} not found")
    try:
        return yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise StudyError(f"{path}: invalid YAML ({exc})") from exc


def skill_kind() -> Kind:
    """Kind of the skill this script lives in, from the name in its SKILL.md frontmatter."""
    path = SKILL_DIR / "SKILL.md"
    match = re.match(r"^---\n(.*?)\n---\n", path.read_text(encoding="utf-8"), re.S)
    name = (yaml.safe_load(match.group(1)) or {}).get("name") if match else None
    if name not in KINDS:
        raise StudyError(f"{path}: frontmatter name {name!r} is not one of {sorted(KINDS)}")
    return KINDS[name]


# --- study.yaml -----------------------------------------------------------------------------

# Schema: a type (or tuple of types) for a leaf; {"fields": {...}, "required": {...}} for a
# mapping; ("list", item_schema) for a list.
_STR, _INT, _NUM = str, int, (int, float)


def _obj(fields: dict, required: tuple = ()) -> dict:
    return {"fields": fields, "required": set(required)}


_REPORT = _obj({"top_n": _INT})
_SCHEMAS: dict[str, dict] = {
    "expression": _obj(
        {
            "study_id": _STR,
            "species": _STR,
            "assembly": _STR,
            "gene_symbol": _STR,
            "window": _obj(
                {
                    "policy": _STR,
                    "explicit_interval": (_STR, type(None)),
                    "regulatory_search_upstream_bp": _INT,
                    "max_feature_gap_bp": _INT,
                    "feature_biotypes": ("list", _STR),
                },
                required=("policy",),
            ),
            "cohorts": _obj(
                {
                    side: _obj({"label": _STR, "curies": ("list", _STR)}, required=("label", "curies"))
                    for side in ("target", "control")
                },
                required=("target", "control"),
            ),
            "scorers": ("list", _STR),
            "report": _REPORT,
        },
        required=("study_id", "species", "assembly", "gene_symbol", "window", "cohorts", "scorers"),
    ),
    "splicing": _obj(
        {
            "study_id": _STR,
            "assembly": _STR,
            "gene_symbol": _STR,
            "transcript_id": (_STR, type(None)),
            "introns_of_interest": ("list", _INT),
            "isoforms": _obj({"reference": (_STR, type(None)), "alternative": (_STR, type(None))}),
            "saturation_region": (_STR, type(None)),
            "saturation": _obj({"flank_bp": _INT}),
            "distance_bins": ("list", _INT),
            "manual_elements": (
                "list",
                _obj(
                    {"label": _STR, "start": _INT, "end": _INT, "citation": _STR},
                    required=("label", "start", "end", "citation"),
                ),
            ),
            "cohort": _obj(
                {
                    "sources": ("list", _STR),
                    "class_rules": _STR,
                    "include_indels": bool,
                    "gnomad_filter": (_STR, type(None)),
                }
            ),
            "composite": _obj({"junction_weight": _NUM}),
            "threshold": _obj({"fixed": ("list", _obj({"label": _STR, "value": _NUM}, required=("label", "value")))}),
            "view": _obj({"start": _INT, "end": _INT}),
            "report": _REPORT,
        },
        required=("study_id", "assembly", "gene_symbol", "introns_of_interest"),
    ),
}


def _check(value: Any, schema: Any, path: str) -> None:
    if isinstance(schema, dict):
        if not isinstance(value, dict):
            raise StudyError(f"{path}: expected a mapping, got {type(value).__name__}")
        for key in value:
            if key not in schema["fields"]:
                raise StudyError(f"{path}.{key}: unknown field" if path else f"{key}: unknown field")
        for key in sorted(schema["required"]):
            if key not in value:
                raise StudyError(f"{path}.{key}: required field is missing" if path else f"{key}: required field is missing")
        for key, item in value.items():
            _check(item, schema["fields"][key], f"{path}.{key}" if path else key)
    elif isinstance(schema, tuple) and schema and schema[0] == "list":
        if not isinstance(value, list):
            raise StudyError(f"{path}: expected a list, got {type(value).__name__}")
        for i, item in enumerate(value):
            _check(item, schema[1], f"{path}[{i}]")
    else:
        types = schema if isinstance(schema, tuple) else (schema,)
        if isinstance(value, bool) and bool not in types:
            raise StudyError(f"{path}: expected {_names(types)}, got bool")
        if not isinstance(value, types):
            raise StudyError(f"{path}: expected {_names(types)}, got {type(value).__name__}")


def _names(types: tuple) -> str:
    return " or ".join("null" if t is type(None) else t.__name__ for t in types)


@dataclasses.dataclass(frozen=True, kw_only=True)
class StudyConfig:
    """Validated contents of study.yaml.

    Attributes:
        kind: "expression" or "splicing".
        study_id: Study identifier (lowercase gene symbol by default).
        assembly: Genome assembly.
        gene_symbol: HGNC gene symbol.
        data: The full validated mapping, as written in study.yaml (nested fields by key).
    """

    kind: Kind
    study_id: str
    assembly: str
    gene_symbol: str
    data: dict[str, Any]

    @classmethod
    def from_dict(cls, d: Any, *, kind: Kind) -> "StudyConfig":
        """Validates the schema fields of `kind` and rejects unknown ones.

        Args:
            d: Parsed study.yaml.
            kind: "expression" or "splicing".

        Raises:
            StudyError: Naming the offending field (dotted path) for an unknown, missing or
                wrongly typed field.
        """
        if kind not in _SCHEMAS:
            raise StudyError(f"unknown study kind {kind!r}")
        _check(d, _SCHEMAS[kind], "")
        for key in ("study_id", "assembly", "gene_symbol"):
            if not d[key].strip():
                raise StudyError(f"{key}: must not be empty")
        flank = (d.get("saturation") or {}).get("flank_bp")
        if flank is not None and flank < 1:
            raise StudyError(f"saturation.flank_bp: must be at least 1, got {flank}")
        return cls(kind=kind, study_id=d["study_id"], assembly=d["assembly"],
                   gene_symbol=d["gene_symbol"], data=d)


def load_study(study_dir: Path) -> StudyConfig:
    """Reads and validates <study_dir>/study.yaml for the kind of this skill."""
    return StudyConfig.from_dict(_read_yaml(Path(study_dir) / STUDY_FILE), kind=skill_kind())


def saturation_flank_bp(config: StudyConfig) -> int:
    """Splicing: study.yaml saturation.flank_bp, else DEFAULT_SATURATION_FLANK_BP."""
    flank = (config.data.get("saturation") or {}).get("flank_bp")
    return DEFAULT_SATURATION_FLANK_BP if flank is None else int(flank)


# --- decisions ------------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True, kw_only=True)
class DecisionEntry:
    """One researcher decision in decisions.yaml.

    Attributes:
        value: The answer (any YAML value).
        workflow_step: Workflow step number that asked the question.
        asked_at: ISO-8601 time the question was raised.
        answered_at: ISO-8601 time the answer was recorded.
        note: Free text with the question, options and rationale.
        options: The answers the pending question offered, or None when the key was decided
            without a pending question (then decisions.yaml has neither this nor `recommended`).
            Each is {code, label} (older studies: the code alone).
        recommended: The code the pending question recommended (None when it had none, or
            when `options` is None).
        value_label: The label of the chosen option as it was offered, for the report (None when
            the answer is not one of the offered options, e.g. a list of CURIEs).
        decided_by: "researcher" (an answer to a question) or "skill" (a fixed rule the skill applied
            without asking; `rule` says which).
        rule: For decided_by "skill", the rule in plain words.
        previous: Earlier answers to the same key, oldest first, when the decision was redefined (e.g. the
            groups changed after the results). Each keeps its value, label, time, note and the options it was
            chosen from, so the report shows what changed and why. None when never redefined.
        said: The researcher's own words the answer comes from, copied from the conversation (None in older
            studies and for the skill's rules). The report shows it next to the answer, so a choice the agent made
            without the researcher has nothing to quote and stands out.
        delegated: True when the researcher told the agent to choose (e.g. "use the recommended window"); `said`
            then quotes that instruction and the report lists the answer as delegated.
    """

    value: Any
    workflow_step: int
    asked_at: str
    answered_at: str
    note: str
    options: list[Any] | None = None
    recommended: Any = None
    value_label: str | None = None
    decided_by: str = "researcher"
    rule: str | None = None
    previous: list[dict[str, Any]] | None = None
    said: str | None = None
    delegated: bool = False


def load_decisions(study_dir: Path) -> dict[str, DecisionEntry]:
    """Reads decisions.yaml; an absent file means no decisions yet."""
    path = Path(study_dir) / DECISIONS_FILE
    if not path.is_file():
        return {}
    raw = _read_yaml(path) or {}
    if not isinstance(raw, dict):
        raise StudyError(f"{path}: expected a mapping of decision key to entry")
    out: dict[str, DecisionEntry] = {}
    for key, entry in raw.items():
        extra = set(entry) - set(DECISION_FIELDS) if isinstance(entry, dict) else set()
        offered = extra & set(OFFERED_FIELDS)
        if (not isinstance(entry, dict) or not set(DECISION_FIELDS) <= set(entry)
                or not extra <= set(OFFERED_FIELDS) | set(TRACE_FIELDS) or offered not in (set(), set(OFFERED_FIELDS))):
            raise StudyError(f"{path}: decision {key!r} must have the fields {list(DECISION_FIELDS)}, optionally both "
                             f"of {list(OFFERED_FIELDS)} and any of {list(TRACE_FIELDS)}")
        if offered and not isinstance(entry["options"], list):
            raise StudyError(f"{path}: decision {key!r}: options must be a list")
        if "previous" in entry and not (isinstance(entry["previous"], list)
                                        and all(isinstance(x, dict) and "value" in x for x in entry["previous"])):
            raise StudyError(f"{path}: decision {key!r}: previous must be a list of earlier answers (each with a value)")
        if entry.get("decided_by", "researcher") not in DECIDED_BY:
            raise StudyError(f"{path}: decision {key!r}: decided_by must be one of {list(DECIDED_BY)}")
        body = {f: _text(entry[f]) if f.endswith("_at") else entry[f] for f in DECISION_FIELDS}
        out[str(key)] = DecisionEntry(**body, **{f: entry[f] for f in extra})
    return out


def option_codes(options: Any) -> list[str]:
    """Codes of offered options, in the new ({code, label}) or the old (code only) form."""
    return [str(o["code"]) if isinstance(o, dict) else str(o) for o in options or []]


def option_label(options: Any, code: Any) -> str | None:
    """Label of the option whose code is `code` (None when absent or when options carry no label)."""
    for o in options or []:
        if isinstance(o, dict) and str(o.get("code")) == str(code):
            return o.get("label")
    return None


def _text(value: Any) -> str:
    """YAML turns unquoted timestamps into datetimes; keep them as ISO text."""
    return value.isoformat() if isinstance(value, (datetime.date, datetime.datetime)) else str(value)


def _write_decisions(study_dir: Path, decisions: dict[str, DecisionEntry]) -> None:
    body = {}
    for key, entry in decisions.items():
        body[key] = dataclasses.asdict(entry)
        if entry.options is None:  # decided without a pending question: no offered fields
            for f in OFFERED_FIELDS:
                del body[key][f]
        for f in ("value_label", "rule", "previous", "said"):  # optional trace fields stay out when empty (older files unchanged)
            if body[key][f] is None:
                del body[key][f]
        if not entry.delegated:
            del body[key]["delegated"]
        if entry.decided_by == "researcher":
            del body[key]["decided_by"]
    _atomic_write(study_dir / DECISIONS_FILE, yaml.safe_dump(body, sort_keys=False, allow_unicode=True, width=1000))


def require(study_dir: Path, key: str, step: int) -> Any:
    """Returns the recorded value for `key`, or exits 1 with the missing-decision message.

    Args:
        study_dir: Study directory.
        key: Decision key, e.g. "locus.transcript".
        step: Workflow step that obtains the answer, named in the message.
    """
    entry = load_decisions(study_dir).get(key)
    if entry is None:
        _fail(f"Missing decision '{key}'. Run Workflow step {step} and record the answer with study.py decide.")
    return entry.value


def sync_rules(study_dir: Path, prefixes: tuple[str, ...], rules: list[dict[str, Any]], step: int) -> None:
    """Makes the skill's rules under `prefixes` match `rules` ({key, value, rule[, note]}): adds or updates them
    and removes the ones that no longer apply (e.g. after a recapture), in one read and at most one write.
    Researcher answers are never touched; an unchanged rule keeps its first timestamps."""
    study_dir = Path(study_dir)
    wanted = {r["key"]: r for r in rules}
    now = _now()
    with _locked(study_dir):
        decisions = load_decisions(study_dir)
        changed = False
        for k in [k for k, e in decisions.items()
                  if e.decided_by == "skill" and k.startswith(prefixes) and k not in wanted]:
            del decisions[k]
            changed = True
        for key, r in wanted.items():
            old = decisions.get(key)
            if old is not None and old.decided_by == "researcher":
                continue  # a researcher's answer is never overwritten by a rule
            note = r.get("note", "")
            if old is not None and old.value == r["value"] and old.rule == r["rule"] and old.note == note:
                continue
            decisions[key] = DecisionEntry(value=r["value"], workflow_step=step, asked_at=now, answered_at=now,
                                           note=note, decided_by="skill", rule=r["rule"])
            changed = True
        if changed:
            _write_decisions(study_dir, decisions)


def record_rule(study_dir: Path, key: str, value: Any, step: int, rule: str, note: str = "") -> None:
    """Records one choice the skill made by a fixed rule, without asking (decided_by: skill).

    Kept in decisions.yaml next to the researcher's answers so the report and the manifest show every
    choice that shaped the result. Re-running a stage with the same value keeps the first timestamps.
    """
    sync_rules(study_dir, (), [{"key": key, "value": value, "rule": rule, "note": note}], step)


# --- pending decisions ----------------------------------------------------------------------


@dataclasses.dataclass(frozen=True, kw_only=True)
class Option:
    """One answer to a question.

    Attributes:
        code: What is recorded in decisions.yaml and read by the stages (stable; never shown to the
            researcher).
        label: The answer in plain words.
        effect: What happens if it is chosen (counts, what is kept or dropped, what to do next).
    """

    code: str
    label: str
    effect: str = ""


def _as_option(o: Any) -> Option:
    if isinstance(o, Option):
        return o
    if isinstance(o, dict):
        return Option(code=str(o["code"]), label=str(o.get("label") or o["code"]), effect=str(o.get("effect") or ""))
    return Option(code=str(o), label=str(o))


@dataclasses.dataclass(frozen=True, kw_only=True)
class PendingDecision:
    """A question a stage found and the researcher must answer.

    The texts are plain English with no jargon; the agent presents them in the researcher's language
    (the language of the conversation) and records the answer with the option's code.

    Attributes:
        key: Decision key (<area>.<topic>[:<subject>]).
        question: The question, in plain words.
        options: The accepted answers (Option, or a bare code in older callers).
        recommended: The code of the suggested answer, or None.
        workflow_step: Workflow step to run to obtain the answer.
        title: A short name for the question.
        why: Why the answer matters for the conclusions, in one sentence.
        identifiers: What the question is about, for traceability (CURIEs, transcript and feature IDs,
            coordinates, counts); shown next to the names and kept with the answer.
        agent_note: An instruction for the agent (what to run first, how to present a long list); never
            part of what the researcher is asked.
    """

    key: str
    question: str
    options: list[Any]
    recommended: str | None
    workflow_step: int
    title: str = ""
    why: str = ""
    identifiers: dict[str, Any] = dataclasses.field(default_factory=dict)
    agent_note: str = ""

    def __post_init__(self) -> None:
        object.__setattr__(self, "options", [_as_option(o) for o in self.options])

    @property
    def codes(self) -> list[str]:
        return [o.code for o in self.options]


def _load_pending(study_dir: Path) -> list[dict]:
    path = Path(study_dir) / PENDING_FILE
    if not path.is_file():
        return []
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise StudyError(f"{path}: invalid JSON ({exc})") from exc
    if not isinstance(data, list):
        raise StudyError(f"{path}: expected a list")
    for item in data:  # older files list bare codes
        item["options"] = [dataclasses.asdict(_as_option(o)) for o in item.get("options") or []]
    return data


def _write_pending(study_dir: Path, items: list[dict]) -> None:
    _atomic_write(study_dir / PENDING_FILE, json.dumps(items, indent=2, ensure_ascii=False) + "\n")


def clear_pending(study_dir: Path, prefixes: tuple[str, ...], keep: tuple[str, ...] = ()) -> None:
    """Drops pending questions under `prefixes` (a stage's areas) except the keys in `keep`: questions the stage
    no longer raises (answered elsewhere, made moot by another answer, or now a rule) must not be asked."""
    study_dir = Path(study_dir)
    with _locked(study_dir):
        items = _load_pending(study_dir)
        left = [p for p in items if not p["key"].startswith(prefixes) or p["key"] in keep]
        if len(left) != len(items):
            _write_pending(study_dir, left)


def raise_pending(study_dir: Path, items: list[PendingDecision], prefixes: tuple[str, ...] = ()) -> NoReturn:
    """Merges `items` into processed/pending_decisions.json and exits 1.

    Items with a key already pending replace the older entry; with `prefixes` (the stage's areas), older
    entries under them that the stage no longer raises are dropped. Exit message:
    "{n} decision(s) needed: run study.py pending." with n the total pending count.
    """
    study_dir = Path(study_dir)
    clear_pending(study_dir, prefixes, tuple(i.key for i in items)) if prefixes else None
    with _locked(study_dir):
        merged = {p["key"]: p for p in _load_pending(study_dir)}
        for item in items:
            previous = merged.get(item.key, {})
            merged[item.key] = {**dataclasses.asdict(item), "asked_at": previous.get("asked_at", _now())}
        _write_pending(study_dir, list(merged.values()))
    _fail(f"{len(merged)} decision(s) needed: run study.py pending.")


# --- max window -----------------------------------------------------------------------------


def _locus(study_dir: Path) -> dict | None:
    path = study_dir / LOCUS_FILE
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise StudyError(f"{path}: invalid JSON ({exc})") from exc


def _saturation_windows(study_dir: Path, config: StudyConfig) -> list[dict] | None:
    """Splicing: the saturation windows processed/locus.json recorded, checked against study.yaml.

    Returns None when locus.json does not exist yet.

    Raises:
        StudyError: locus.json predates the saturation windows, or study.yaml saturation_region /
            saturation.flank_bp changed since the locus step wrote them.
    """
    locus = _locus(study_dir)
    if locus is None:
        return None
    plan = locus.get("saturation")
    if not isinstance(plan, dict) or not plan.get("windows"):
        raise StudyError(f"{study_dir / LOCUS_FILE} has no saturation windows (written by an older version); "
                         "run splicing.py locus again (Workflow step 4).")
    flank = saturation_flank_bp(config)
    if plan.get("rule") != "flank" or plan.get("flank_bp") != flank:
        raise StudyError(f"{study_dir / LOCUS_FILE} holds saturation windows for {plan.get('rule')!r} "
                         f"(flank_bp {plan.get('flank_bp')}), but study.yaml now asks for +/- {flank} bp around "
                         "each splice site; run splicing.py locus again (Workflow step 4).")
    return plan["windows"]


def max_window(study_dir: Path) -> int:
    """Largest window to score: max(1000, longest window length), the AlphaGenome --max-window value.

    Expression: `window.length` of processed/locus.json. Splicing: the explicit
    `saturation_region` of study.yaml (one interval), else the longest of the saturation windows
    in processed/locus.json (+/- saturation.flank_bp around each splice site of the introns of
    interest, merged where they overlap; written by the locus step).

    Raises:
        StudyError: When neither source exists yet (the locus step has not run), or when
            locus.json is stale for the study.yaml saturation settings.
    """
    study_dir = Path(study_dir)
    config = load_study(study_dir)
    length: int | None = None
    if config.kind == "expression":
        locus = _locus(study_dir)
        if locus is not None:
            length = int(locus["window"]["length"])
    else:
        region = config.data.get("saturation_region")
        if region and region.strip():
            match = _REGION_RE.match(region)
            if not match:
                raise StudyError(f"saturation_region: expected CHR:START-END, got {region!r}")
            start, end = int(match.group(2)), int(match.group(3))
            if end < start:
                raise StudyError(f"saturation_region: end {end} is before start {start}")
            length = end - start + 1
        else:
            windows = _saturation_windows(study_dir, config)
            if windows:
                length = max(int(w["end"]) - int(w["start"]) + 1 for w in windows)
    if length is None:
        raise StudyError(
            f"No window available for {study_dir}: {LOCUS_FILE} does not exist and study.yaml "
            "has no saturation_region. Run the locus step first (Workflow step 4)."
        )
    return max(MIN_WINDOW, length)


# --- exploration and CURIE replacement --------------------------------------------------------


_REPLACE_ANSWER = re.compile(r"^replace:(\S+)")
_EXPLORATION_RAW = (re.compile(r"^gtf\.json$"), re.compile(r"^regulatory\.json$"), re.compile(r"^sequence_.+\.txt$"),
                    re.compile(r"^ols_search_.+\.json$"), re.compile(r"^ols_term_.+\.json$"))
_CATALOGUE = re.compile(r"^tracks_([A-Z0-9_]+)_catalogue\.parquet$")


def apply_replacement(study_dir: Path, key: str, value: Any) -> str | None:
    """After `cohort.unserved:{curie}` = replace:{new}, writes {new} in place of {curie} in study.yaml.

    study.yaml stays the request the study actually ran on, and a comment at its end keeps the original CURIE
    and the decision that changed it. Returns a one-line report, or None when nothing changed.
    """
    if not key.startswith("cohort.unserved:"):
        return None
    m = _REPLACE_ANSWER.match(str(value).strip())
    if m is None:
        return None
    old, new = key.split(":", 1)[1], m.group(1)
    path = Path(study_dir) / STUDY_FILE
    text = path.read_text(encoding="utf-8")
    token = re.compile(rf"(?<![\w:]){re.escape(old)}(?![\w:])")
    lines = text.splitlines(keepends=True)
    changed = False
    for i, line in enumerate(lines):
        body = line.split("#", 1)[0]
        if ("curies:" in body or body.lstrip().startswith("- ")) and token.search(body):
            lines[i] = token.sub(new, body, count=1) + line[len(body):]
            changed = True
    if not changed:
        return None
    note = f"# {key}: {old} replaced by {new} on {_now()} (decisions.yaml keeps the original cohort answer)\n"
    new_text = "".join(lines)
    new_text = new_text if new_text.endswith("\n") else new_text + "\n"
    StudyConfig.from_dict(yaml.safe_load(new_text), kind=skill_kind())
    _atomic_write(path, new_text + note)
    return f"study.yaml: {old} replaced by {new} in the cohort CURIEs (original kept in a comment)."


def apply_cohort_answer(study_dir: Path, key: str, value: Any, label: str | None = None) -> str | None:
    """After `cohort.target` / `cohort.control` is decided (again), writes the CURIEs (and the group name, when
    given) into that group of study.yaml, so the request matches the decision without editing the file by hand.

    Returns a one-line report, or None when the key is not a group or nothing changed.
    """
    if key not in ("cohort.target", "cohort.control") or not isinstance(value, list):
        return None
    side = key.split(".", 1)[1]
    path = Path(study_dir) / STUDY_FILE
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    start = next((i for i, ln in enumerate(lines) if re.match(rf"^  {side}:\s*(#.*)?$", ln)), None)
    if start is None:
        return None
    changed = []
    for i in range(start + 1, len(lines)):
        ln = lines[i]
        if re.match(r"^ {0,2}\S", ln):
            break
        body, comment = (ln.split("#", 1) + [""])[:2]
        comment = ("  #" + comment.rstrip("\n")) if comment else ""
        m = re.match(r"^(\s*)(curies|label):", body)
        if m and m.group(2) == "curies":
            lines[i] = f"{m.group(1)}curies: {_flow([str(v) for v in value])}{comment}\n"
            changed.append("CURIEs")
        elif m and m.group(2) == "label" and label:
            lines[i] = f"{m.group(1)}label: {json.dumps(label, ensure_ascii=False)}{comment}\n"
            changed.append("name")
    if not changed:
        return None
    new_text = "".join(lines)
    StudyConfig.from_dict(yaml.safe_load(new_text), kind=skill_kind())
    _atomic_write(path, new_text)
    return f"study.yaml: {side} group updated ({' and '.join(changed)}) from decision '{key}'."


def _exploration_root(path: Path) -> Path:
    path = Path(path)
    for candidate in (path, path / "exploration"):
        if (candidate / "exploration.json").is_file():
            return candidate
    raise StudyError(f"no exploration.json in {path} or {path / 'exploration'}; run expression.py explore first.")


def _fill(text: str, pattern: str, value: str) -> str:
    new, n = re.subn(pattern, lambda m: m.group(1) + value, text, count=1, flags=re.M)
    if n != 1:
        raise StudyError(f"study template: line {pattern!r} not found")
    return new


def _flow(items: list[str]) -> str:
    return "[" + ", ".join(items) + "]"


# The answers a study cannot start without (Workflow step 1). The researcher gives each one; the agent never
# fills one from the recommendation. init --from-exploration refuses to run without all of them and records them
# in decisions.yaml with the options offered, next to cohort.target and cohort.control.
REQUIRED_ANSWERS = {
    "window": ("study.window", "--window", "which window (region) to analyse, by its name in summary.md"),
    "target": ("cohort.target", "--target", "which tissues or cells form the target group"),
    "control": ("cohort.control", "--control", "which tissues or cells form the control group"),
    "scorers": ("study.scorers", "--scorers", "which kinds of prediction to compare (gene expression is always included)"),
    "question": ("analysis.question", "--question", "what is looked for in the target: a decrease, an increase or any change"),
}
QUESTION_OPTIONS = [{"code": "repression", "label": "Diminui a expressão no alvo"},
                    {"code": "increase", "label": "Aumenta a expressão no alvo"},
                    {"code": "either", "label": "Qualquer mudança no alvo"}]
SCORER_LABELS = {"RNA_SEQ": "Expressão do gene (RNA-seq)", "CAGE": "Início de transcrição (CAGE)",
                 "DNASE": "Cromatina aberta (DNase)", "ATAC": "Cromatina aberta (ATAC)"}


def missing_answers(args: argparse.Namespace) -> list[str]:
    """The step-1 answers not given on the command line, as '--flag: what to ask'."""
    return [f"{flag}: {what}" for attr, (_, flag, what) in REQUIRED_ANSWERS.items() if not getattr(args, attr, None)]


def _record_answers(study_dir: Path, args: argparse.Namespace, data: dict, scorers: list[str]) -> None:
    """Writes the researcher's step-1 answers (window, kinds of prediction, question) with what was offered."""
    now = _now()
    note = args.answer_note or "Resposta do pesquisador no passo 1 (a partir da exploração)"
    delegated = {STEP1_DELEGABLE[d] for d in (args.delegated or [])}
    windows = [{"code": w["id"], "label": w["label"]} for w in data["windows"]]
    explored = [{"code": s, "label": SCORER_LABELS.get(s, s)} for s in data["scorers"]]
    entries = {
        "study.window": (args.window, windows, data.get("recommended_window")),
        "study.scorers": (scorers, explored, None),
        "analysis.question": (args.question, QUESTION_OPTIONS, None),
    }
    with _locked(study_dir):
        decisions = load_decisions(study_dir)
        for key, (value, options, rec) in entries.items():
            decisions[key] = DecisionEntry(value=value, workflow_step=1, asked_at=now, answered_at=now, note=note,
                                           options=options, recommended=rec,
                                           value_label=option_label(options, value) if isinstance(value, str) else None,
                                           previous=_history(decisions.get(key), value),
                                           said=args.said.strip(), delegated=key in delegated)
        _write_decisions(study_dir, decisions)


def init_from_exploration(args: argparse.Namespace) -> None:
    """study.yaml from the options the researcher chose in the exploration, and the exploration captures in raw/."""
    import shutil

    import manifest  # same scripts/ folder; only this path needs it

    if skill_kind() != "expression":
        raise StudyError("--from-exploration is for expression studies")
    missing = missing_answers(args)
    if not (args.said or "").strip():
        missing.append(f"--said: {SAID_HELP}")
    if missing:
        raise StudyError("the study needs the researcher's answers to step 1, and these were not given:\n  "
                         + "\n  ".join(missing) + "\nAsk the researcher (all of them in one message, with the options of "
                         "summary.md); never fill an answer from the recommendation. Nothing was created.")
    root = _exploration_root(Path(args.from_exploration))
    data = json.loads((root / "exploration.json").read_text(encoding="utf-8"))
    windows = {w["id"]: w for w in data["windows"]}
    if args.window not in windows:
        raise StudyError(f"--window {args.window!r} is not one of the exploration windows {list(windows)}")
    both = sorted(set(args.target) & set(args.control))
    if both:
        raise StudyError(f"{both} in both --target and --control; a tissue cannot be its own control")
    scorers = list(args.scorers)
    if "RNA_SEQ" not in scorers:
        raise StudyError("RNA_SEQ must be among the scorers (the expression contrast)")
    missing = [s for s in scorers if not (root / "raw" / f"tracks_{s}_catalogue.parquet").is_file()]
    if missing:
        raise StudyError(f"the exploration has no catalogue for {missing}; run explore with --scorers including them")
    gene = data["gene"]
    study_dir = Path(args.study_dir)
    target = study_dir / STUDY_FILE
    if target.exists():
        raise StudyError(f"{target} already exists; refusing to overwrite an existing study.")
    window = windows[args.window]
    policy = window["policy"]
    text = (SKILL_DIR / "assets" / "study_template.yaml").read_text(encoding="utf-8")
    text = text.replace("{{STUDY_ID}}", (args.study_id or gene.lower())).replace("{{GENE_SYMBOL}}", gene)
    text = _fill(text, r"^(\s*policy:\s*)\S+", policy["policy"])
    text = _fill(text, r"^(\s*explicit_interval:\s*)\S+", policy.get("explicit_interval") or "null")
    text = _fill(text, r"^(\s*regulatory_search_upstream_bp:\s*)\d+", str(data["flank_bp"]))
    text = _fill(text, r"^(\s*max_feature_gap_bp:\s*)\d+", str(policy.get("max_feature_gap_bp", 0)))
    text = _fill(text, r"^(\s*label:\s*)Target", json.dumps(args.target_label or "Target", ensure_ascii=False))
    text = _fill(text, r"^(\s*curies:\s*)\[\]", _flow(args.target))
    text = _fill(text, r"^(\s*label:\s*)Control", json.dumps(args.control_label or "Control", ensure_ascii=False))
    text = _fill(text, r"^(\s*curies:\s*)\[\]", _flow(args.control))
    text = _fill(text, r"^(scorers:\s*)\[[^\]]*\]", _flow(scorers))
    text = text.replace("# Ask the researcher for the gene, target and control tissues, scorers and window policy\n"
                        "# (Workflow step 1), then run `study.py validate`.",
                        f"# Written by study.py init --from-exploration from {root} (window option {args.window}).")
    StudyConfig.from_dict(yaml.safe_load(text), kind="expression")

    study_dir.mkdir(parents=True, exist_ok=True)
    for sub in STUDY_SUBDIRS:
        (study_dir / sub).mkdir(exist_ok=True)
    target.write_text(text, encoding="utf-8")
    invocation = " ".join(__import__("shlex").quote(a) for a in ["uv", "run", sys.argv[0], *sys.argv[1:]])
    manifest.add_step(study_dir, manifest.ManifestStep(
        step="1-study", skill="expression-analysis", command=invocation, outputs=[STUDY_FILE],
        source_versions={"exploration": str(root.resolve()), "window_option": args.window},
        recorded_at=manifest._now(), decision_keys=["study.window", "study.scorers", "analysis.question"]),
        origin_kind="live", origin_note=args.origin_note or f"Created from the exploration in {root}")
    _record_answers(study_dir, args, data, scorers)

    explored = json.loads((root / "manifest.json").read_text(encoding="utf-8"))["steps"]
    copied = 0
    untraced: list[str] = []
    for src in sorted((root / "raw").iterdir()):
        cat = _CATALOGUE.match(src.name)
        if not (any(r.match(src.name) for r in _EXPLORATION_RAW) or (cat and cat.group(1) in scorers)):
            continue
        rel = f"raw/{src.name}"
        step = next((s for s in explored if rel in s.get("outputs", [])), None)
        if step is None:  # e.g. the error reply of a search that failed: it cannot be traced, so it stays out
            untraced.append(src.name)
            continue
        shutil.copy2(src, study_dir / rel)
        manifest.add_step(study_dir, manifest.ManifestStep(
            step=step["step"], skill=step["skill"], command=step["command"], outputs=[rel],
            source_versions={**step.get("source_versions", {}), "copied_from": str((root / rel).resolve()),
                             "captured_at": step.get("recorded_at", "")},
            recorded_at=manifest._now(), decision_keys=[]))
        copied += 1
    # The RNA_SEQ_ACTIVE scores of the exploration's probe SNV: the dashboard reconstructs the reference
    # expression (REF) from them, to tell where the gene's predicted level is low. Optional (older explorations have none).
    rel_src = "atlas/scores_RNA_SEQ_ACTIVE_baseline.tsv"
    step = next((s for s in explored if rel_src in s.get("outputs", [])), None)
    if (root / rel_src).is_file() and step is not None:
        rel = "raw/scores_RNA_SEQ_ACTIVE_baseline.tsv"
        shutil.copy2(root / rel_src, study_dir / rel)
        manifest.add_step(study_dir, manifest.ManifestStep(
            step=step["step"], skill=step["skill"], command=step["command"], outputs=[rel],
            source_versions={**step.get("source_versions", {}), "copied_from": str((root / rel_src).resolve()),
                             "captured_at": step.get("recorded_at", "")},
            recorded_at=manifest._now(), decision_keys=[]))
        copied += 1
    print(f"Created {target} (study_id {args.study_id or gene.lower()}, window {window['region']} [{args.window}], "
          f"scorers {scorers}, question {args.question}) and copied {copied} exploration file(s) into raw/ with their "
          "manifest steps. Recorded study.window, study.scorers and analysis.question in decisions.yaml. Record the "
          "groups with study.py decide (cohort.target, cohort.control with the options shown), then run study.py validate."
          + (f" Left out (no manifest step in the exploration, e.g. a failed search): {', '.join(untraced)}." if untraced else ""))


# --- CLI ------------------------------------------------------------------------------------


def init(args: argparse.Namespace) -> None:
    if args.from_exploration:
        init_from_exploration(args)
        return
    if not args.gene:
        raise StudyError("--gene is required (or --from-exploration)")
    study_dir = Path(args.study_dir)
    target = study_dir / STUDY_FILE
    if target.exists():
        raise StudyError(f"{target} already exists; refusing to overwrite an existing study.")
    template = SKILL_DIR / "assets" / "study_template.yaml"
    if not template.is_file():
        raise StudyError(f"{template} not found")
    gene = args.gene.strip()
    if not gene:
        raise StudyError("--gene must not be empty")
    text = template.read_text(encoding="utf-8")
    text = text.replace("{{STUDY_ID}}", gene.lower()).replace("{{GENE_SYMBOL}}", gene)
    study_dir.mkdir(parents=True, exist_ok=True)
    for sub in STUDY_SUBDIRS:
        (study_dir / sub).mkdir(exist_ok=True)
    target.write_text(text, encoding="utf-8")
    print(f"Created {target} (study_id {gene.lower()}); edit it with the researcher, then run study.py validate.")


def validate(args: argparse.Namespace) -> None:
    config = load_study(Path(args.study_dir))
    print(f"{Path(args.study_dir) / STUDY_FILE}: valid {config.kind} study '{config.study_id}'")


def _offered(args: argparse.Namespace, pending_entry: dict | None, earlier: DecisionEntry | None) -> dict:
    """The options/recommended to store with a decision: --option/--recommended, else the pending entry's,
    else an earlier decision's of the same key, else none."""
    if args.option:
        if (earlier is not None and earlier.options and all(isinstance(o, dict) for o in earlier.options)
                and set(args.option) == set(option_codes(earlier.options))):
            # the same question asked again with its codes only: keep the labels it was offered with
            return {"options": earlier.options, "recommended": args.recommended}
        return {"options": list(args.option), "recommended": args.recommended}
    if pending_entry is not None:  # what the question offered (code and label), next to the answer
        return {"options": [{"code": o["code"], "label": o["label"]} for o in pending_entry.get("options") or []],
                "recommended": pending_entry.get("recommended")}
    if earlier is not None and earlier.options is not None:  # a revised answer to the same question
        return {"options": earlier.options, "recommended": earlier.recommended}
    return {}


_ENSEMBL_ID = re.compile(r"^(ENS[A-Z]*\d+)(?:\.\d+)?$")


def _match_code(raw: str, parsed: Any, codes: list[str]) -> tuple[bool, Any]:
    """Whether the answer `raw` is one of `codes`, and the value to record.

    The text is compared before YAML parsing (so `yes`, `no` or `01` stay codes). A code with a `<…>`
    placeholder (e.g. `fixed:<value>`) takes any text after its prefix. Ensembl IDs match with or without
    their version and are recorded as offered. Lists and mappings (e.g. CURIEs) are not checked.
    """
    if raw in codes:
        return True, raw
    for c in codes:
        if "<" in c and c.endswith(">"):
            prefix = c[:c.index("<")]
            if raw.startswith(prefix) and len(raw) > len(prefix):
                return True, raw
    m = _ENSEMBL_ID.match(raw)
    if m:
        same = [c for c in codes if (mc := _ENSEMBL_ID.match(c)) and mc.group(1) == m.group(1)]
        if len(same) == 1:
            return True, same[0]
    if isinstance(parsed, (list, dict)):
        return True, parsed
    return False, parsed


def _history(old: DecisionEntry | None, answer: Any) -> list[dict[str, Any]] | None:
    """The earlier answers to keep when a key is decided again: the old entry is appended when the value
    changes (a redefinition), so nothing the researcher decided before is lost."""
    if old is None:
        return None
    hist = list(old.previous or [])
    if old.value != answer:
        snap = {f: getattr(old, f) for f in PREVIOUS_FIELDS}
        hist.append({f: v for f, v in snap.items()
                     if v is not None and v != "" and not (f == "decided_by" and v == "researcher")
                     and not (f == "delegated" and v is False)})
    return hist or None


def decide(args: argparse.Namespace) -> None:
    study_dir = Path(args.study_dir)
    if not study_dir.is_dir():
        raise StudyError(f"study directory {study_dir} not found")
    keys = list(dict.fromkeys(args.key))
    if not (args.said or "").strip():
        raise StudyError(f"--said is required: {SAID_HELP}. If the researcher has not answered, ask them; never record "
                         "an answer of your own. Nothing was recorded.")
    if args.recommended is not None and not args.option:
        raise StudyError("--recommended needs the options it was chosen from: pass each with --option")
    groups = [k for k in keys if k in ("cohort.target", "cohort.control")]
    if groups and args.recommended is not None:
        raise StudyError(f"--recommended does not apply to {', '.join(groups)}: the exploration recommends one term per "
                         "text the researcher gave, not a group. Pass the terms offered with --option and leave "
                         "--recommended out. Nothing was recorded.")
    if args.option and args.recommended is not None and args.recommended not in args.option:
        raise StudyError(f"--recommended {args.recommended!r} is not one of the --option values {args.option}")
    try:
        value = yaml.safe_load(args.value)
    except yaml.YAMLError:
        value = args.value
    now = _now()
    with _locked(study_dir):
        pending = {p["key"]: p for p in _load_pending(study_dir)}
        decisions = load_decisions(study_dir)
        bad = [k for k in keys if k not in pending and k not in decisions and not KEY_PATTERN.match(k)]
        if bad:
            raise StudyError(
                f"invalid decision key(s) {', '.join(repr(k) for k in bad)}: a key is <area>.<topic>[:<subject>] "
                "(lower-case area and topic, no whitespace), as printed by study.py pending; pass one --key per key. "
                "Nothing was recorded.")
        # --option says the agent asked the question itself, so no stage had it pending
        unknown = [] if args.option else [k for k in keys if k not in pending and k not in decisions]
        for key in keys:
            previous = pending.get(key)
            codes = option_codes(previous.get("options")) if previous else []
            answer = value
            if codes:
                ok, answer = _match_code(args.value.strip(), value, codes)
                if not ok:
                    listed = "; ".join(f"{o['code']} = {o['label']}" for o in previous["options"])
                    raise StudyError(f"'{args.value}' is not an answer to '{key}'. Record one of the codes ({listed}). "
                                     "Nothing was recorded.")
            offered = _offered(args, previous, decisions.get(key))
            decisions[key] = DecisionEntry(
                previous=_history(decisions.get(key), answer),
                value=answer,
                workflow_step=args.step,
                asked_at=previous["asked_at"] if previous and "asked_at" in previous else now,
                answered_at=now,
                note=args.note or "",
                value_label=option_label(offered.get("options"), answer),
                said=args.said.strip(),
                delegated=bool(args.delegated),
                **offered,
            )
        _write_decisions(study_dir, decisions)
        if any(k in pending for k in keys):
            _write_pending(study_dir, [p for k, p in pending.items() if k not in keys])
    for key in unknown:
        print(f"Warning: '{key}' is neither pending (processed/pending_decisions.json) nor decided before; recorded "
              "anyway. Check the spelling against study.py pending unless the question was asked by the agent "
              "(e.g. cohort.target).", file=sys.stderr)
    for key in keys:
        print(f"Recorded decision '{key}' = {json.dumps(value, ensure_ascii=False)} in {study_dir / DECISIONS_FILE}")
        if (study_dir / STUDY_FILE).is_file():
            for report in (apply_replacement(study_dir, key, value),
                           apply_cohort_answer(study_dir, key, value, getattr(args, "group_label", None))):
                if report:
                    print(report)


def require_cmd(args: argparse.Namespace) -> None:
    value = require(Path(args.study_dir), args.key, args.step)
    print(json.dumps(value, ensure_ascii=False))


def pending(args: argparse.Namespace) -> None:
    items = _load_pending(Path(args.study_dir))
    if args.json:
        print(json.dumps(items, indent=2, ensure_ascii=False))
        return
    if not items:
        print("No pending decisions.")
        return
    print(f"{len(items)} question(s) for the researcher. Ask one at a time (questions of the same kind may be asked "
          "together), in the researcher's language: the title, the question, why it matters, the names with their "
          "identifiers, each option's label and what happens, and which one is recommended and why. With more than about "
          "six options, show the first ones and offer the rest on request. Never show codes or keys; record the answer "
          "with the code of the chosen option. Follow the 'For the agent' lines before asking.")
    for n, item in enumerate(items, 1):
        print(f"\n[{n}] {item.get('title') or item['key']}")
        print(f"    Question: {item.get('question')}")
        if item.get("why"):
            print(f"    Why it matters: {item['why']}")
        if item.get("identifiers"):
            print("    About: " + "; ".join(f"{k}: {v}" for k, v in item["identifiers"].items()))
        if item.get("agent_note"):
            print(f"    For the agent: {item['agent_note']}")
        print("    Options:")
        for o in item["options"]:
            rec = " (recommended)" if o["code"] == item.get("recommended") else ""
            effect = f": {o['effect']}" if o.get("effect") else ""
            print(f"      - {o['label']}{rec}{effect}  [code: {o['code']}]")
        print(f"    Record (agent): study.py decide --study-dir D --key \"{item['key']}\" --value CODE "
              f"--step {item.get('workflow_step')} --said \"<the researcher's words>\" --note \"<their reason>\"")
    print("\nSeveral keys with the same answer go in one decide call, one --key each.")


def window(args: argparse.Namespace) -> None:
    if not args.print_max_window:
        raise StudyError("nothing to do: pass --print-max-window")
    print(max_window(Path(args.study_dir)))


def main() -> None:
    for stream in (sys.stdout, sys.stderr):  # UTF-8 on any console (cp1252 on Windows)
        stream.reconfigure(encoding="utf-8", errors="replace")
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command_name", required=True)

    def study_dir_arg(p: argparse.ArgumentParser) -> None:
        p.add_argument("--study-dir", required=True, help="Study directory")

    p = sub.add_parser("init", help="Create a study directory and study.yaml from the template")
    study_dir_arg(p)
    p.add_argument("--gene", help="Gene symbol, e.g. the HGNC symbol (taken from the exploration with --from-exploration)")
    p.add_argument("--from-exploration", metavar="DIR",
                   help="Build study.yaml from an `expression.py explore` folder and copy its captures into raw/")
    p.add_argument("--window", help="--from-exploration: window option id from summary.md (e.g. regulatory_unit)")
    p.add_argument("--target", action="append", help="--from-exploration: target CURIE (repeatable)")
    p.add_argument("--control", action="append", help="--from-exploration: control CURIE (repeatable)")
    p.add_argument("--target-label", help="--from-exploration: short name of the target group")
    p.add_argument("--control-label", help="--from-exploration: short name of the control group")
    p.add_argument("--scorers", nargs="+", help="--from-exploration: kinds of prediction the researcher chose (required)")
    p.add_argument("--question", choices=["repression", "increase", "either"],
                   help="--from-exploration: what the researcher looks for in the target (required)")
    p.add_argument("--answer-note", help="--from-exploration: note stored with the step-1 answers")
    p.add_argument("--said", help=f"--from-exploration (required): {SAID_HELP}")
    p.add_argument("--delegated", nargs="+", choices=sorted(STEP1_DELEGABLE),
                   help="--from-exploration: answers the researcher explicitly told you to choose (quoted in --said)")
    p.add_argument("--study-id", help="--from-exploration: study id (default: lower-case gene symbol)")
    p.add_argument("--origin-note", help="--from-exploration: origin note of the manifest (who asked what, when)")
    p.set_defaults(func=init)

    p = sub.add_parser("validate", help="Validate study.yaml")
    study_dir_arg(p)
    p.set_defaults(func=validate)

    p = sub.add_parser("decide", help="Record a researcher decision")
    study_dir_arg(p)
    p.add_argument("--key", required=True, action="append",
                   help="Decision key, <area>.<topic>[:<subject>]; repeat it to give several keys the same value")
    p.add_argument("--value", required=True, help="Answer; parsed as YAML (e.g. '[3]'), else kept as text")
    p.add_argument("--step", required=True, type=int, help="Workflow step that asked")
    p.add_argument("--note", help="Question, options and rationale")
    p.add_argument("--option", action="append",
                   help="An option offered for a question the agent asked itself (no pending entry); repeat it "
                        "for each option. Stored as `options` next to the answer")
    p.add_argument("--recommended", help="The option recommended to the researcher (one of --option)")
    p.add_argument("--group-label", help="With cohort.target / cohort.control: the group's short name, written to study.yaml")
    p.add_argument("--said", help=f"Required: {SAID_HELP}")
    p.add_argument("--delegated", action="store_true",
                   help="The researcher explicitly told you to choose this one (their words in --said)")
    p.set_defaults(func=decide)

    p = sub.add_parser("require", help="Exit 1 with the missing-decision message if the key has no answer")
    study_dir_arg(p)
    p.add_argument("--key", required=True)
    p.add_argument("--step", required=True, type=int, help="Workflow step that obtains the answer")
    p.set_defaults(func=require_cmd)

    p = sub.add_parser("pending", help="List the open questions for the researcher")
    study_dir_arg(p)
    p.add_argument("--json", action="store_true", help="Print the questions as JSON (title, why, identifiers, options)")
    p.set_defaults(func=pending)

    p = sub.add_parser("window", help="Window sizing")
    study_dir_arg(p)
    p.add_argument("--print-max-window", action="store_true", help="Print max(1000, longest window length)")
    p.set_defaults(func=window)

    args = parser.parse_args()
    try:
        args.func(args)
    except (StudyError, OSError, KeyError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
