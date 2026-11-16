"""UBERON/CL cohorts and the Atlas track registry (Workflow step 8, `expression.py cohorts`).

The target and control CURIEs the researcher chose (decisions `cohort.target`, `cohort.control`)
are matched with the Atlas track catalogue (raw/tracks_{scorer}_{curie}.parquet), the OLS records
(raw/ols_search_*.json, raw/ols_term_*.json) give each CURIE its label and definition, and every
CURIE or scorer that would leave a cohort empty becomes a question instead of an empty cohort.

Researcher questions (answered with study.py decide, Workflow step 8):
  cohort.unserved:{curie}             the CURIE has no track for any requested scorer:
                                      drop_curie | replace:{descendant CURIE}
  cohort.scorer_without_cohort:{s}    scorer s has no track in the target or in the control cohort:
                                      drop_scorer
  cohort.no_rna_seq:{role}            the target or control cohort keeps no RNA_SEQ track:
                                      add:{descendant CURIE} | abort
Any question also accepts `abort`, which stops the stage with exit 1. RNA_SEQ is the scorer of the
expression contrast, so a cohort without an RNA_SEQ track cannot be fixed by dropping the scorer.
`add:{curie}` adds an OLS descendant of one of that cohort's CURIEs to the cohort (the CURIEs it
already has stay, for the other scorers); the recommended one is the descendant with the most
RNA_SEQ tracks in the catalogue.

`replace:{curie}` and `add:{curie}` take OLS descendants (from raw/ols_term_{curie}.json, relation
hierarchicalDescendants, else hierarchicalChildren). When the whole-catalogue files
raw/tracks_{scorer}_catalogue.parquet were captured (Workflow step 7), the questions list only the
descendants they serve (any requested scorer for replace, RNA_SEQ for add, with track counts);
without them every descendant is listed and the question says so. After choosing one, pull that
CURIE's tracks (Workflow step 7) and run this stage again. A replacement that turns out to have no
tracks raises its own `cohort.unserved:{descendant}` question. `study.py decide` also writes a
`replace:` answer into study.yaml (this stage reads the CURIEs from decisions.yaml either way).

Outputs (processed/): track_registry.tsv, ontology_resolution.json, tissue_cohorts.md.

Library: CohortError, CohortPending, resolve_cohorts, chosen_curies, run, load_resolution, stage_metrics.
"""

from __future__ import annotations

import csv
import json
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import models
from core import locus as core_locus
from core import tracks as core_tracks
from models import AtlasTrack
from study import DecisionEntry, PendingDecision, StudyConfig, load_decisions, load_study

WORKFLOW_STEP = 8
OLS_STEP = 5
TRACKS_STEP = 7
GENE_EXPRESSION_SCORER = "RNA_SEQ"
REGISTRY_FILE = "processed/track_registry.tsv"
RESOLUTION_FILE = "processed/ontology_resolution.json"
TISSUE_DOC_FILE = "processed/tissue_cohorts.md"
MANIFEST_FILE = "manifest.json"
_QUERY_RE = re.compile(r"--ontology\s+(\S+)\s+--scorer\s+(\S+)")
_REPLACE_RE = re.compile(r"^replace:(\S+)(?:\s.*)?$")
_ADD_RE = re.compile(r"^add:(\S+)(?:\s.*)?$")
NO_RNA_SEQ_KEY = "cohort.no_rna_seq"


class CohortError(Exception):
    """Invalid or incomplete cohort input; the CLI prints it and exits with 1."""


class CohortPending(Exception):
    """Questions the researcher must answer before the cohorts can be built.

    Attributes:
        items: The questions, for study.raise_pending.
    """

    def __init__(self, items: list[PendingDecision]) -> None:
        super().__init__(f"{len(items)} decision(s) needed")
        self.items = items


# --- raw inputs -----------------------------------------------------------------------------


def _curie_file(curie: str) -> str:
    return curie.replace(":", "_", 1)


def _queried_pairs(study_dir: Path) -> set[tuple[str, str]]:
    """(scorer, CURIE) pairs the manifest records a catalogue command for (zero-track pairs have no file)."""
    path = study_dir / MANIFEST_FILE
    if not path.is_file():
        return set()
    try:
        steps = json.loads(path.read_text(encoding="utf-8")).get("steps", [])
    except (json.JSONDecodeError, AttributeError) as exc:
        raise CohortError(f"{path}: invalid manifest ({exc})") from None
    pairs: set[tuple[str, str]] = set()
    for step in steps:
        match = _QUERY_RE.search(step.get("command", ""))
        if match and "metadata --tracks" in step["command"]:
            pairs.add((match.group(2), match.group(1)))
    return pairs


def _load_tracks(study_dir: Path, curies: Sequence[str], scorers: Sequence[str]) -> dict[tuple[str, str], list[AtlasTrack]]:
    """Catalogue tracks per (scorer, CURIE). A pair with neither a file nor a manifest step was never queried."""
    queried = _queried_pairs(study_dir)
    out: dict[tuple[str, str], list[AtlasTrack]] = {}
    missing: list[str] = []
    for curie in curies:
        for scorer in scorers:
            path = study_dir / "raw" / f"tracks_{scorer}_{_curie_file(curie)}.parquet"
            if path.is_file():
                out[(scorer, curie)] = list(models.load_raw(path))
            elif (scorer, curie) in queried:
                out[(scorer, curie)] = []
            else:
                missing.append(f"raw/{path.name}")
    if missing:
        shown = ", ".join(missing[:6]) + (f" and {len(missing) - 6} more" if len(missing) > 6 else "")
        raise CohortError(
            f"Track catalogue missing for {len(missing)} CURIE x scorer pair(s): {shown}. Run Workflow step "
            f"{TRACKS_STEP} (expression.py fetch tracks --study-dir D, which queries every missing pair), then run "
            "cohorts again."
        )
    return out


def _catalogue_counts(study_dir: Path, scorer: str) -> dict[str, int] | None:
    """Tracks per CURIE in raw/tracks_{scorer}_catalogue.parquet (the whole catalogue); None when not captured."""
    path = study_dir / "raw" / f"tracks_{scorer}_catalogue.parquet"
    if not path.is_file():
        return None
    counts: dict[str, int] = {}
    for track in models.load_raw(path):
        counts[track.ontology_curie] = counts.get(track.ontology_curie, 0) + 1
    return counts


def _served_by(study_dir: Path, scorers: Sequence[str]) -> dict[str, list[str]] | None:
    """Scorers that serve each CURIE, from raw/tracks_{scorer}_catalogue.parquet; None unless every scorer has one."""
    served: dict[str, set[str]] = {}
    for scorer in scorers:
        path = study_dir / "raw" / f"tracks_{scorer}_catalogue.parquet"
        if not path.is_file():
            return None
        for track in models.load_raw(path):
            served.setdefault(track.ontology_curie, set()).add(scorer)
    return {c: [s for s in scorers if s in ss] for c, ss in served.items()}


def _ols_info(study_dir: Path) -> tuple[dict[str, dict[str, Any]], dict[str, list[tuple[str, str]]]]:
    """OLS records by CURIE, and descendants (CURIE, label) by term CURIE.

    Priority for a CURIE: its own get_term record, then a search hit, then a mention in a relation list.
    """
    info: dict[str, dict[str, Any]] = {}
    descendants: dict[str, list[tuple[str, str]]] = {}

    def put(curie: str, record: dict[str, Any]) -> None:
        info.setdefault(curie, record)

    raw = study_dir / "raw"
    terms = [models.load_raw(p)[0].term for p in sorted(raw.glob("ols_term_*.json"))]
    for term in terms:
        put(term.obo_id, _record(term.obo_id, term.label, term.iri, term.ontology_name, term.description, term.synonyms))
        related = term.hierarchicalDescendants if term.hierarchicalDescendants is not None else term.hierarchicalChildren
        descendants[term.obo_id] = [(r.obo_id, r.label) for r in related or ()]
    for path in sorted(raw.glob("ols_search_*.json")):
        for hit in models.load_raw(path)[0].terms:
            put(hit.obo_id, _record(hit.obo_id, hit.label, hit.iri, hit.ontology_name, hit.description, hit.exact_synonyms))
    for term in terms:
        for rel in (term.hierarchicalChildren or ()) + (term.hierarchicalDescendants or ()):
            put(rel.obo_id, _record(rel.obo_id, rel.label, rel.iri, rel.obo_id.split(":")[0].lower(), (), ()))
    return info, descendants


def _record(curie: str, label: str, iri: str, ontology: str, description: Sequence[str], synonyms: Sequence[str]) -> dict[str, Any]:
    return {"obo_id": curie, "label": label, "iri": iri, "ontology_name": ontology,
            "description": description[0] if description else None, "synonyms": list(synonyms)[:6]}


# --- decisions ------------------------------------------------------------------------------


def _cohort_curies(decisions: Mapping[str, DecisionEntry]) -> dict[str, list[str]]:
    """Target and control CURIEs from the decisions (Workflow step 5), which settle the request of study.yaml."""
    out: dict[str, list[str]] = {}
    for role in ("target", "control"):
        entry = decisions.get(f"cohort.{role}")
        if entry is None:
            raise CohortError(
                f"Missing decision 'cohort.{role}'. Run Workflow step {OLS_STEP} and record the answer with study.py decide.")
        value = entry.value
        if not isinstance(value, list) or not value or not all(isinstance(v, str) for v in value):
            raise CohortError(f"decision 'cohort.{role}' must be a non-empty list of CURIEs, got {value!r}")
        if len(set(value)) != len(value):
            raise CohortError(f"decision 'cohort.{role}' lists a CURIE twice: {value}")
        out[role] = list(value)
    both = sorted(set(out["target"]) & set(out["control"]))
    if both:
        raise CohortError(
            f"CURIE(s) {', '.join(both)} are in both cohort.target and cohort.control; a tissue cannot be its own "
            f"control. Choose different CURIEs for one cohort (Workflow step {OLS_STEP}) and record them with study.py decide.")
    return out


def _choice(decisions: Mapping[str, DecisionEntry], key: str) -> str | None:
    entry = decisions.get(key)
    return None if entry is None else str(entry.value).strip()


# --- resolution -----------------------------------------------------------------------------


def resolve_cohorts(
    study: StudyConfig,
    tracks: Mapping[tuple[str, str], Sequence[AtlasTrack]],
    ols: Mapping[str, Mapping[str, Any]],
    descendants: Mapping[str, Sequence[tuple[str, str]]],
    decisions: Mapping[str, DecisionEntry],
    served: Mapping[str, Sequence[str]] | None = None,
    rna_seq_counts: Mapping[str, int] | None = None,
) -> dict[str, Any]:
    """Matches the chosen cohorts with the catalogue and applies the researcher's answers.

    Args:
        study: StudyConfig (labels of the cohorts, requested scorers).
        tracks: Catalogue tracks per (scorer, CURIE) for every CURIE the cohorts name; replacements named
            by a decision are looked up in the same mapping.
        ols: OLS record per CURIE (obo_id, label, iri, ontology_name, description, synonyms).
        descendants: (CURIE, label) pairs of each term's OLS descendants.
        decisions: decisions.yaml entries.
        served: Scorers serving each CURIE in the whole catalogue (raw/tracks_{scorer}_catalogue.parquet); the
            replacement options of an unserved CURIE are limited to its descendants found here. None (no
            catalogue captured) offers every OLS descendant and says so in the question.
        rna_seq_counts: RNA_SEQ tracks per CURIE in the whole catalogue (raw/tracks_RNA_SEQ_catalogue.parquet), for
            the options of cohort.no_rna_seq:{role}; None (not captured) offers every OLS descendant.

    Returns:
        The ontology resolution (written to processed/ontology_resolution.json): requested and effective
        scorers, per-CURIE status (kept, dropped, replaced), track counts and the answers that shaped it.

    Raises:
        CohortPending: Questions to answer (cohort.unserved:*, cohort.scorer_without_cohort:*, cohort.no_rna_seq:*).
        CohortError: An `abort` answer, an invalid answer, a CURIE without an OLS record, or a cohort still
            without RNA_SEQ after its cohort.no_rna_seq answer.
    """
    requested = list(study.data["scorers"])
    chosen = _cohort_curies(decisions)
    pending: list[PendingDecision] = []
    used: list[str] = []
    entries: list[dict[str, Any]] = []
    seen: set[str] = {c for cs in chosen.values() for c in cs}

    def entry(curie: str, role: str) -> dict[str, Any]:
        if curie not in ols:
            raise CohortError(
                f"No OLS record for {curie} in raw/. Run Workflow step {OLS_STEP} (search_ols.py) or step 6 (get_term.py) for it.")
        rec = ols[curie]
        per_scorer = {s: len(tracks.get((s, curie), ())) for s in requested}
        listed = [t for s in requested for t in tracks.get((s, curie), ())]
        return {
            "curie": curie, "role": role, "status": "kept", "replaced_by": None,
            "ols": dict(rec),
            "obo_url": core_tracks.obo_purl(curie),
            "ols_url": core_tracks.ols_term_url(rec["ontology_name"], rec["iri"]),
            "tracks_per_scorer": per_scorer,
            "tracks": [{"scorer": t.scorer, "track_index": t.track_index, "name": t.track_name, "strand": t.strand,
                        "biosample_name": t.biosample_name, "biosample_type": t.biosample_type} for t in listed],
        }

    def expand(curie: str, role: str) -> None:
        item = entry(curie, role)
        entries.append(item)
        if sum(item["tracks_per_scorer"].values()):
            return
        key = f"cohort.unserved:{curie}"
        answer = _choice(decisions, key)
        if answer is None:
            pending.append(_unserved_question(item, requested, descendants.get(curie), served))
            return
        used.append(key)
        if answer == "drop_curie":
            item["status"] = "dropped"
        elif answer == "abort":
            raise CohortError(f"Aborted by the researcher at '{key}'.")
        elif (m := _REPLACE_RE.match(answer)):
            target = m.group(1)
            if target not in {c for c, _ in descendants.get(curie, ())}:
                raise CohortError(
                    f"Decision '{key}': {target} is not an OLS descendant of {curie} in raw/ols_term_{_curie_file(curie)}.json.")
            if target in seen:
                raise CohortError(f"Decision '{key}': {target} is already in a cohort; choose another replacement or drop_curie.")
            seen.add(target)
            item["status"], item["replaced_by"] = "replaced", target
            expand(target, role)
        else:
            raise CohortError(f"Decision '{key}' is {answer!r}; expected drop_curie, replace:CURIE or abort.")

    # Replacement CURIEs are not known until the decisions are read, so tracks are looked up lazily.
    for role in ("target", "control"):
        for curie in chosen[role]:
            expand(curie, role)

    def has_rna_seq(role: str) -> bool:
        return any(e["tracks_per_scorer"].get(GENE_EXPRESSION_SCORER) for e in entries
                   if e["status"] == "kept" and e["role"] == role)

    if GENE_EXPRESSION_SCORER in requested:
        waiting = {q.key for q in pending}
        for role in ("target", "control"):
            if has_rna_seq(role):
                continue
            role_curies = [e["curie"] for e in entries if e["role"] == role]
            if any(f"cohort.unserved:{c}" in waiting for c in role_curies):
                continue  # answer the unserved question first: a replacement may bring RNA_SEQ
            key = f"{NO_RNA_SEQ_KEY}:{role}"
            answer = _choice(decisions, key)
            if answer is None:
                pending.append(_no_rna_seq_question(role, entries, descendants, rna_seq_counts, seen))
                continue
            used.append(key)
            if answer == "abort":
                raise CohortError(f"Aborted by the researcher at '{key}'.")
            m = _ADD_RE.match(answer)
            if m is None:
                raise CohortError(f"Decision '{key}' is {answer!r}; expected add:CURIE or abort.")
            target = m.group(1)
            offered = {c for curie in role_curies for c, _ in descendants.get(curie, ())}
            if target not in offered:
                raise CohortError(f"Decision '{key}': {target} is not an OLS descendant of the {role} cohort's CURIEs "
                                  f"({', '.join(role_curies)}) in raw/ols_term_*.json.")
            if target in seen:
                raise CohortError(f"Decision '{key}': {target} is already in a cohort; choose another descendant.")
            seen.add(target)
            expand(target, role)
            for e in entries:
                if e["curie"] == target:
                    e["added_by"] = key

    kept = [e for e in entries if e["status"] == "kept"]
    dropped_scorers: list[str] = []
    adding = {q.key.split(":", 1)[1] for q in pending if q.key.startswith(f"{NO_RNA_SEQ_KEY}:")}
    for scorer in requested:
        empty = [role for role in ("target", "control")
                 if not any(e["tracks_per_scorer"].get(scorer) for e in kept if e["role"] == role)]
        if not empty or scorer == GENE_EXPRESSION_SCORER or set(empty) <= adding:
            continue  # (a descendant added for RNA_SEQ may bring this scorer's tracks too: ask that first)
        key = f"cohort.scorer_without_cohort:{scorer}"
        answer = _choice(decisions, key)
        if answer is None:
            pending.append(_scorer_question(scorer, empty, kept))
            continue
        used.append(key)
        if answer == "drop_scorer":
            dropped_scorers.append(scorer)
        elif answer == "abort":
            raise CohortError(f"Aborted by the researcher at '{key}'.")
        else:
            raise CohortError(f"Decision '{key}' is {answer!r}; expected drop_scorer or abort.")
    if pending:
        raise CohortPending(pending)

    if GENE_EXPRESSION_SCORER in requested:
        for role in ("target", "control"):
            if not has_rna_seq(role):
                key = f"{NO_RNA_SEQ_KEY}:{role}"
                raise CohortError(
                    f"The {role} cohort still keeps no {GENE_EXPRESSION_SCORER} track after '{key}' = "
                    f"{_choice(decisions, key)!r}: its tracks (Workflow step 7) have no {GENE_EXPRESSION_SCORER} track. "
                    f"Record another answer to '{key}' with study.py decide, or choose other CURIEs for cohort.{role} "
                    f"(Workflow step {OLS_STEP}).")
    effective = [s for s in requested if s not in dropped_scorers]
    if not effective:
        raise CohortError("Every requested scorer was dropped; nothing is left to contrast.")
    cohorts = study.data["cohorts"]
    return {
        "study_id": study.study_id,
        "requested_scorers": requested,
        "target_label": cohorts["target"]["label"],
        "control_label": cohorts["control"]["label"],
        "curies": entries,
        "effective_scorers": effective,
        "dropped_scorers": dropped_scorers,
        "effective_target_curies": [e["curie"] for e in kept if e["role"] == "target"],
        "effective_control_curies": [e["curie"] for e in kept if e["role"] == "control"],
        "decisions": [{"key": k, "value": decisions[k].value} for k in sorted(used)],
    }


def _unserved_question(item: Mapping[str, Any], requested: Sequence[str],
                       descendants: Sequence[tuple[str, str]] | None,
                       served: Mapping[str, Sequence[str]] | None) -> PendingDecision:
    curie, role, label = item["curie"], item["role"], item["ols"]["label"]
    options = ["drop_curie"]
    unique = dict(descendants or ())
    if not unique:
        hint = (f" No descendants are recorded: run Workflow step 6 (get_term.py --relations "
                f"hierarchicalChildren,hierarchicalDescendants) for {curie} to see replacements.")
    elif served is None:
        ordered = sorted(unique.items(), key=lambda kv: (not kv[0].startswith("UBERON:"), kv[1].lower(), kv[0]))
        options += [f"replace:{c} ({l})" for c, l in ordered]
        hint = (f" The Atlas catalogue was not captured (raw/tracks_SCORER_catalogue.parquet, Workflow step 7), so all "
                f"{len(ordered)} OLS descendant(s) are listed, served or not; after choosing one pull its catalogue and "
                "run cohorts again.")
    else:
        ordered = sorted(((c, l) for c, l in unique.items() if served.get(c)), key=lambda kv: (kv[1].lower(), kv[0]))
        options += [f"replace:{c} ({l}; served by {', '.join(served[c])})" for c, l in ordered]
        hint = (f" {len(ordered)} of its {len(unique)} OLS descendant(s) are served by at least one requested scorer and "
                "are listed as replace:CURIE options; after choosing one pull its catalogue (Workflow step 7) and run "
                "cohorts again." if ordered else
                f" None of its {len(unique)} OLS descendant(s) is served by a requested scorer.")
    return PendingDecision(
        key=f"cohort.unserved:{curie}",
        question=(f"{role.capitalize()} cohort term {curie} ({label}) has no AlphaGenome track for any requested "
                  f"scorer {list(requested)}. Drop it from the cohort?{hint}"),
        options=options,
        recommended="drop_curie",
        workflow_step=WORKFLOW_STEP,
    )


def _no_rna_seq_question(role: str, entries: Sequence[Mapping[str, Any]],
                         descendants: Mapping[str, Sequence[tuple[str, str]]],
                         counts: Mapping[str, int] | None, seen: set[str]) -> PendingDecision:
    """cohort.no_rna_seq:{role}: add an OLS descendant served by RNA_SEQ, or abort."""
    mine = [e for e in entries if e["role"] == role]
    curies = [e["curie"] for e in mine]
    unique: dict[str, str] = {}
    for curie in curies:
        for c, label in descendants.get(curie, ()):
            if c not in seen:
                unique.setdefault(c, label)
    tracks = ", ".join(f"{e['curie']} ({e['ols']['label']}): "
                       + ", ".join(f"{s} {n}" for s, n in e["tracks_per_scorer"].items()) for e in mine)
    scorer = GENE_EXPRESSION_SCORER
    options: list[str] = []
    recommended = "abort"
    if not unique:
        missing = [c for c in curies if c not in descendants]
        hint = (f" No OLS descendants are recorded for {', '.join(missing or curies)}: run Workflow step 6 (get_term.py "
                "--relations hierarchicalChildren,hierarchicalDescendants) for them, then run cohorts again to see "
                f"descendants served by {scorer}.")
    elif counts is None:
        ordered = sorted(unique.items(), key=lambda kv: (kv[1].lower(), kv[0]))
        options = [f"add:{c} ({l})" for c, l in ordered]
        recommended = None
        hint = (f" The {scorer} catalogue was not captured (raw/tracks_{scorer}_catalogue.parquet, Workflow step 7), so "
                f"all {len(ordered)} OLS descendant(s) are listed, served or not, and none is recommended; capture it "
                "and run cohorts again to see which are served.")
    else:
        ordered = sorted(((c, l) for c, l in unique.items() if counts.get(c)),
                         key=lambda kv: (-counts[kv[0]], kv[1].lower(), kv[0]))
        options = [f"add:{c} ({l}; {counts[c]} {scorer} tracks)" for c, l in ordered]
        if options:
            recommended = options[0]
        hint = (f" {len(ordered)} of their {len(unique)} OLS descendant(s) have {scorer} tracks in the catalogue and are "
                f"listed as add:CURIE options, most tracks first (recommended: the first). After choosing one, pull its "
                "tracks (Workflow step 7) and run cohorts again." if ordered else
                f" None of their {len(unique)} OLS descendant(s) has a {scorer} track in the catalogue; choose other "
                f"CURIEs for cohort.{role} (Workflow step {OLS_STEP}) or abort.")
    return PendingDecision(
        key=f"{NO_RNA_SEQ_KEY}:{role}",
        question=(f"The {role} cohort keeps no {scorer} track ({tracks}), and {scorer} is the expression contrast. Add an "
                  f"OLS descendant that {scorer} serves to the {role} cohort (its CURIEs stay for the other scorers), or "
                  f"abort?{hint}"),
        options=[*options, "abort"],
        recommended=recommended,
        workflow_step=WORKFLOW_STEP,
    )


def _scorer_question(scorer: str, empty: Sequence[str], kept: Sequence[Mapping[str, Any]]) -> PendingDecision:
    by_curie = ", ".join(f"{e['curie']} {e['tracks_per_scorer'].get(scorer, 0)}" for e in kept)
    return PendingDecision(
        key=f"cohort.scorer_without_cohort:{scorer}",
        question=(f"Scorer {scorer} has no track in the {' and '.join(empty)} cohort, so no target-vs-control contrast "
                  f"is possible for it. Drop {scorer} from this study? (tracks by CURIE: {by_curie})"),
        options=["drop_scorer"],
        recommended="drop_scorer",
        workflow_step=WORKFLOW_STEP,
    )


# --- stage ----------------------------------------------------------------------------------


def chosen_curies(decisions: Mapping[str, DecisionEntry]) -> set[str]:
    """Every CURIE the cohorts may use: cohort.target, cohort.control, and the replace:/add: answers."""
    chosen = _cohort_curies(decisions)
    needed = set(chosen["target"]) | set(chosen["control"])
    for key, entry in decisions.items():
        value = str(entry.value).strip()
        if key.startswith("cohort.unserved:"):
            m = _REPLACE_RE.match(value)
        elif key.startswith(f"{NO_RNA_SEQ_KEY}:"):
            m = _ADD_RE.match(value)
        else:
            m = None
        if m:
            needed.add(m.group(1))
    return needed


def run(study_dir: Path) -> dict[str, Any]:
    """Builds the cohorts and registry from raw/ and decisions.yaml; writes the three processed files.

    Raises:
        CohortPending, CohortError, core.locus.LocusError, models.DependencyOutputError, study.StudyError.
    """
    study_dir = Path(study_dir)
    study = load_study(study_dir)
    decisions = load_decisions(study_dir)
    locus = core_locus.load_locus(study_dir)
    ols, descendants = _ols_info(study_dir)
    needed = chosen_curies(decisions)
    tracks = _load_tracks(study_dir, sorted(needed), study.data["scorers"])
    resolution = resolve_cohorts(study, tracks, ols, descendants, decisions, _served_by(study_dir, study.data["scorers"]),
                                 _catalogue_counts(study_dir, GENE_EXPRESSION_SCORER))
    catalogue = [t for ts in tracks.values() for t in ts]
    registry = core_tracks.build_track_registry(
        catalogue, study, resolution=resolution, chromosome=locus.chromosome, gene_strand=locus.strand,
        window=(locus.window.start, locus.window.end))
    out = study_dir / "processed"
    out.mkdir(parents=True, exist_ok=True)
    registry.write_csv(study_dir / REGISTRY_FILE, separator="\t")
    (study_dir / RESOLUTION_FILE).write_text(json.dumps(resolution, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (study_dir / TISSUE_DOC_FILE).write_text(core_tracks.render_tissue_cohorts_markdown(resolution, registry), encoding="utf-8")
    return resolution


def load_resolution(study_dir: Path) -> dict[str, Any]:
    """Reads processed/ontology_resolution.json (written by the cohorts step)."""
    path = Path(study_dir) / RESOLUTION_FILE
    if not path.is_file():
        raise CohortError(f"{RESOLUTION_FILE} not found in {study_dir}; run expression.py cohorts (Workflow step {WORKFLOW_STEP}).")
    return json.loads(path.read_text(encoding="utf-8"))


def stage_metrics(study_dir: Path) -> dict[str, int | float | str]:
    """Metrics of the cohorts stage compared by `reproduce` (expected.json stages.cohorts)."""
    study_dir = Path(study_dir)
    resolution = load_resolution(study_dir)
    path = study_dir / REGISTRY_FILE
    if not path.is_file():
        raise CohortError(f"{REGISTRY_FILE} not found in {study_dir}; run expression.py cohorts (Workflow step {WORKFLOW_STEP}).")
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    metrics: dict[str, int | float | str] = {
        f"n_tracks.{s}": sum(1 for r in rows if r["scorer"] == s) for s in resolution["effective_scorers"]}
    metrics["target_curies"] = ",".join(sorted(resolution["effective_target_curies"]))
    metrics["control_curies"] = ",".join(sorted(resolution["effective_control_curies"]))
    metrics["dropped_scorers"] = ",".join(sorted(resolution["dropped_scorers"]))
    return metrics
