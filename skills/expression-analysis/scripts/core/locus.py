"""Locus and regulatory window of the study gene (Workflow step 4, `expression.py locus`).

Builds the `regulatory_unit` window, evaluates locus decisions, and computes strand-aware
position annotations from the dependency outputs in raw/:

  raw/gtf.json          GENCODE v46 transcripts (alphagenome-variant-impact-score `gtf`). The file
                        writes 0-based starts; every function here that takes GtfTranscript records
                        takes them AS LOADED from raw/ and applies models.normalize_gtf once itself.
  raw/regulatory.json   Ensembl Regulatory Build features (regulatory_features.py overlap).
  raw/sequence_*.txt    Reference sequence (ensembl-database get-sequence). The interval is read
                        from each file name, never recomputed; the file that covers the window is
                        recorded and the start codon must read ATG (reverse complement on '-').

Coordinates are 1-based closed. "Oriented" coordinates are pos * strand, so that upstream is a
smaller oriented value on either strand. The TSS is the 5' end of exon 1: its start on '+', its
end on '-'.

Window policies (study.yaml window.policy):
  regulatory_unit  The `feature_biotypes` feature that contains or abuts the TSS, extended
                   upstream through features that touch, overlap or lie within
                   max_feature_gap_bp, through the 3' end of exon 1.
  explicit         window.explicit_interval (chr:start-end) as given; it must contain the TSS.

Researcher questions (answered with study.py decide, Workflow step 4):
  locus.transcript                  0 or several MANE Select transcripts
  locus.multiple_tss_features       several features contain or abut the TSS: union | abort
  locus.feature_overlap:{a}+{b}     chained features overlap: union | abort
  locus.feature_gap:{a}+{b}         chained features are separated by a gap: bridge_gap | abort
  locus.exon1_overlap:{feature}     a chained feature extends into exon 1: keep_full_unit | abort
  locus.no_feature_chain            no feature contains or abuts the TSS: set_explicit_interval | abort

Library: LocusConfig, build_locus, check_sequence, transcripts_in_window, footprint_blocks,
snv_counts, classify_block, transcript_region, genomic_domain, run, load_locus, stage_metrics.
"""

from __future__ import annotations

import csv
import dataclasses
import json
import re
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

import models
from models import EnsemblSequence, GtfExon, GtfTranscript, RegulatoryFeature
from study import DecisionEntry, Option, PendingDecision, StudyConfig, load_decisions, load_study, sync_rules

WORKFLOW_STEP = 4
GTF_FILE = "raw/gtf.json"
REGULATORY_FILE = "raw/regulatory.json"
SEQUENCE_GLOB = "sequence_*.txt"
LOCUS_FILE = "processed/locus.json"
TRANSCRIPTS_FILE = "processed/transcripts_in_window.tsv"
POLICIES = ("regulatory_unit", "explicit")
DEFAULT_BIOTYPES = ("promoter", "enhancer")
TRANSCRIPT_KEY = "locus.transcript"
NO_CHAIN_KEY = "locus.no_feature_chain"
MULTIPLE_TSS_KEY = "locus.multiple_tss_features"
TRANSCRIPT_COLUMNS = (
    "transcript_id", "transcript_version", "transcript_type", "is_selected", "is_mane_select", "strand",
    "tss", "tss_offset_bp", "first_exon_start", "first_exon_end",
)
# Footprint kind -> name used in the snvs.* metrics.
METRIC_KIND = {"5_prime_utr": "utr5"}
STANDARD_METRIC_KINDS = ("enhancer", "promoter", "utr5", "cds")

_INTERVAL_RE = re.compile(r"^\s*(chr[0-9A-Za-z]+):(\d+)-(\d+)\s*$")
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


class LocusError(Exception):
    """An input prevents building the locus; the CLI prints the message and exits with 1."""


class LocusPending(Exception):
    """Researcher decisions are missing; `items` go to study.raise_pending."""

    def __init__(self, items: list[PendingDecision]):
        self.items = items
        super().__init__(", ".join(item.key for item in items))


# --------------------------------------------------------------------------- data classes


@dataclasses.dataclass(frozen=True, kw_only=True)
class Interval:
    """1-based closed genomic interval.

    Attributes:
        start: First base.
        end: Last base (>= start).
    """

    start: int
    end: int

    def __post_init__(self) -> None:
        if self.start < 1 or self.end < self.start:
            raise ValueError(f"invalid 1-based closed interval {self.start}..{self.end}")

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    def contains(self, pos: int) -> bool:
        return self.start <= pos <= self.end

    def to_dict(self) -> dict[str, Any]:
        extra = {k: v for k, v in dataclasses.asdict(self).items() if k not in ("start", "end")}
        return {**extra, "start": self.start, "end": self.end, "length": self.length}

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]):
        body = dict(d)
        length = body.pop("length", None)
        item = cls(**body)
        if length is not None and length != item.length:
            raise ValueError(f"length {length} inconsistent with {item.start}..{item.end}")
        return item


@dataclasses.dataclass(frozen=True, kw_only=True)
class ExonRef(Interval):
    """Exon 1 of the selected transcript (normalised GENCODE coordinates).

    Attributes:
        start: First base, 1-based.
        end: Last base, 1-based.
        exon_number: Exon number in 5'->3' transcript order (1).
    """

    exon_number: int


@dataclasses.dataclass(frozen=True, kw_only=True)
class FeatureRef(Interval):
    """Ensembl Regulatory Build feature overlapping the window.

    Attributes:
        start: First base, 1-based.
        end: Last base, 1-based.
        id: Stable ID (ENSR...).
        biotype: Regulatory Build `description` (promoter, enhancer, CTCF_binding_site, ...).
    """

    id: str
    biotype: str


@dataclasses.dataclass(frozen=True, kw_only=True)
class LocusConfig:
    """Machine-derived locus contract, written to processed/locus.json.

    Attributes:
        study_id: study.yaml study_id.
        gene_id: Unversioned Ensembl gene ID of the selected transcript.
        gene_version: GENCODE gene version (None when the ID has none).
        gene_symbol: HGNC symbol (study.yaml gene_symbol).
        transcript_id: Unversioned Ensembl transcript ID of the selected transcript.
        transcript_version: GENCODE transcript version (None when the ID has none).
        transcript_type: GENCODE transcript biotype.
        is_mane_select: True when the selected transcript is MANE Select.
        assembly: Genome assembly (study.yaml).
        chromosome: Contig with chr prefix.
        strand: 1 or -1.
        tss: Transcription start site, the 5' end of exon 1 (1-based).
        cds_start_genomic: First CDS base in transcript orientation (first base of the start codon).
        start_codon: 1-based closed genomic span (lo, hi) of the start codon.
        first_exon: Exon 1 of the selected transcript.
        window: Saturation-mutagenesis window.
        window_policy: study.yaml window.policy that produced the window.
        features: Regulatory Build features (any biotype) overlapping the window, by start.
        decisions: {decision key: value} of the researcher answers the window depends on.
        sequence_file: raw/ file whose sequence covers the window (None before check_sequence).
        start_codon_verified: True when a raw sequence covered the start codon and read ATG.
    """

    study_id: str
    gene_id: str
    gene_version: int | None
    gene_symbol: str
    transcript_id: str
    transcript_version: int | None
    transcript_type: str
    is_mane_select: bool
    assembly: str
    chromosome: str
    strand: int
    tss: int
    cds_start_genomic: int
    start_codon: tuple[int, int]
    first_exon: ExonRef
    window: Interval
    window_policy: str
    features: tuple[FeatureRef, ...]
    decisions: Mapping[str, Any]
    sequence_file: str | None = None
    start_codon_verified: bool = False

    def __post_init__(self) -> None:
        if self.strand not in (1, -1):
            raise ValueError(f"strand must be 1 or -1, got {self.strand!r}")
        if not self.window.contains(self.tss):
            raise ValueError(f"TSS {self.tss} lies outside the window {self.window.start}..{self.window.end}")
        lo, hi = self.start_codon
        expected = ((self.cds_start_genomic, self.cds_start_genomic + 2) if self.strand == 1
                    else (self.cds_start_genomic - 2, self.cds_start_genomic))
        if (lo, hi) != expected:
            raise ValueError(f"start_codon {self.start_codon} inconsistent with CDS start {self.cds_start_genomic}")
        for f in self.features:
            if f.end < self.window.start or f.start > self.window.end:
                raise ValueError(f"regulatory feature {f.id} does not overlap the window")

    def distance_to_tss(self, pos: int) -> int:
        """Strand-aware offset of `pos` from the TSS (negative = upstream)."""
        return (pos - self.tss) * self.strand

    def upstream_of_tss(self, pos: int) -> bool:
        return self.distance_to_tss(pos) < 0

    @property
    def window_str(self) -> str:
        return models.interval_str(self.chromosome, self.window.start, self.window.end)

    def to_dict(self) -> dict[str, Any]:
        out: dict[str, Any] = {}
        for field in dataclasses.fields(self):
            value = getattr(self, field.name)
            if isinstance(value, Interval):
                value = value.to_dict()
            elif field.name == "features":
                value = [f.to_dict() for f in value]
            elif field.name == "start_codon":
                value = list(value)
            elif field.name == "decisions":
                value = dict(value)
            out[field.name] = value
        return out

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> LocusConfig:
        """Rebuilds a LocusConfig from processed/locus.json; raises ValueError/KeyError/TypeError."""
        names = {f.name for f in dataclasses.fields(cls)}
        unknown = sorted(set(d) - names)
        if unknown:
            raise ValueError(f"unknown field {unknown[0]!r}")
        body = dict(d)
        body["first_exon"] = ExonRef.from_dict(d["first_exon"])
        body["window"] = Interval.from_dict(d["window"])
        body["features"] = tuple(FeatureRef.from_dict(f) for f in d["features"])
        body["start_codon"] = tuple(d["start_codon"])
        body["decisions"] = dict(d["decisions"])
        return cls(**body)


@dataclasses.dataclass(frozen=True, kw_only=True)
class FootprintBlock:
    """Contiguous run of window positions with the same annotation.

    Attributes:
        kind: Regulatory biotype, "5_prime_utr", "cds", "unannotated" (upstream, no feature) or
            "other" (downstream of exon 1, explicit windows only).
        feature_id: Regulatory feature ID for regulatory kinds, else None.
        start: First base, 1-based.
        end: Last base, 1-based.
    """

    kind: str
    feature_id: str | None
    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    @property
    def snvs(self) -> int:
        return 3 * self.length


# --------------------------------------------------------------------------- helpers


def strip_version(stable_id: str) -> str:
    return stable_id.split(".", 1)[0]


def _version(stable_id: str) -> int | None:
    m = re.match(r"^[^.]+\.(\d+)", stable_id)
    return int(m.group(1)) if m else None


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _strand(t: GtfTranscript) -> int:
    if t.strand == "+":
        return 1
    if t.strand == "-":
        return -1
    raise LocusError(f"{GTF_FILE}: transcript {t.transcript_id} has strand {t.strand!r}, expected '+' or '-'")


def _contig_without_prefix(chromosome: str) -> str:
    name = chromosome[3:] if chromosome.startswith("chr") else chromosome
    return "MT" if name == "M" else name


def _oriented(start: int, end: int, strand: int) -> tuple[int, int]:
    a, b = start * strand, end * strand
    return (a, b) if a <= b else (b, a)


def _first_exon(t: GtfTranscript) -> GtfExon:
    """Exon 1 of a normalised transcript, checked against coordinate order and the transcript span."""
    if not t.exons:
        raise LocusError(f"{GTF_FILE}: transcript {t.transcript_id} has no exons (capture gtf with --exons)")
    strand = _strand(t)
    ordered = sorted(t.exons, key=lambda e: e.start, reverse=strand == -1)
    numbered = [e for e in t.exons if e.exon_number == 1]
    if len(numbered) != 1 or numbered[0] != ordered[0]:
        raise LocusError(f"{GTF_FILE}: exon 1 of {t.transcript_id} is not its 5'-most exon on strand {t.strand}")
    first = ordered[0]
    if (strand == 1 and first.start != t.start) or (strand == -1 and first.end != t.end):
        raise LocusError(f"{GTF_FILE}: the 5' end of exon 1 of {t.transcript_id} does not match the transcript span "
                         f"{t.start}-{t.end}")
    return first


def _tss(t: GtfTranscript) -> int:
    """TSS of a normalised transcript: start on '+', end on '-'."""
    return t.start if _strand(t) == 1 else t.end


def _cds_length(t: GtfTranscript) -> int:
    return sum(s.end - s.start + 1 for s in t.cds)


def _question(key: str, question: str, options: Iterable[Option], recommended: str | None, *, title: str,
              why: str, identifiers: dict[str, Any]) -> PendingDecision:
    return PendingDecision(key=key, question=question, options=list(options), recommended=recommended,
                           workflow_step=WORKFLOW_STEP, title=title, why=why, identifiers=identifiers)


WINDOW_WHY = ("The window is the set of variants every conclusion rests on: what it leaves out is never scored, "
              "and every base it adds costs three more variants to capture.")
STOP_WINDOW = Option(code="abort", label="Pause to choose another window",
                     effect="The analysis pauses here. You then choose another window from the exploration page (for "
                            "example the start of the gene +/- 500 bp); the agent sets it and the analysis continues "
                            "from this step.")


def _span(a: int, b: int) -> str:
    return f"{b - a + 1:,} bp, {3 * (b - a + 1):,} variants"


def _resolve(triggers: list[PendingDecision], decisions: Mapping[str, DecisionEntry]) -> dict[str, Any]:
    """Recorded answers for `triggers`; raises LocusPending for missing ones, LocusError on abort."""
    missing = [t for t in triggers if t.key not in decisions]
    if missing:
        raise LocusPending(missing)
    chosen: dict[str, Any] = {}
    for t in triggers:
        value = decisions[t.key].value
        if str(value) not in t.codes:
            raise LocusError(f"Decision '{t.key}' = {value!r} is not one of {t.codes}; record it again with "
                             f"study.py decide (Workflow step {WORKFLOW_STEP}).")
        chosen[t.key] = str(value)
    aborted = sorted(k for k, v in chosen.items() if v == "abort")
    if aborted:
        raise LocusError(f"Stopped: the researcher answered 'abort' to {aborted}. Change window.policy in "
                         "study.yaml (e.g. explicit with window.explicit_interval) or record another answer.")
    return chosen


# --------------------------------------------------------------------------- transcript


def select_transcript(transcripts: list[GtfTranscript], gene_symbol: str,
                      decisions: Mapping[str, DecisionEntry]) -> GtfTranscript:
    """The single MANE Select transcript, or the one named by the `locus.transcript` decision.

    Args:
        transcripts: Normalised records of raw/gtf.json.
        gene_symbol: study.yaml gene_symbol; every record must carry it.
        decisions: Recorded decisions; `locus.transcript` is read only when there is not exactly
            one MANE Select transcript (0 or several, e.g. several gene IDs share the symbol).

    Raises:
        LocusPending: The decision is needed and missing.
        LocusError: Empty or foreign gtf.json, abort, unknown transcript or no CDS.
    """
    if not transcripts:
        raise LocusError(f"{GTF_FILE} has no transcripts; run Workflow step 2 (alphagenome-variant-impact-score gtf).")
    names = sorted({t.gene_name for t in transcripts})
    if names != [gene_symbol]:
        raise LocusError(f"{GTF_FILE} holds transcripts of {names}, but study.yaml asks for {gene_symbol!r}; "
                         "capture gtf again for the study gene (Workflow step 2).")
    mane = [t for t in transcripts if t.is_mane_select]
    if len(mane) == 1:
        chosen = mane[0]
    else:
        coding = sorted((t for t in transcripts if t.cds), key=lambda t: (-_cds_length(t), t.transcript_id))
        pool = mane or coding
        candidates = [t.transcript_id for t in pool]
        entry = decisions.get(TRANSCRIPT_KEY)
        if entry is None:
            gene_ids = sorted({t.gene_id for t in transcripts})
            genes = f" The symbol {gene_symbol} matches {len(gene_ids)} genes ({', '.join(gene_ids)})." if len(gene_ids) > 1 else ""
            basis = (f"{len(mane)} transcripts are marked MANE Select (the reference transcript agreed by NCBI and "
                     "Ensembl)" if mane else
                     "no transcript is marked MANE Select (the reference agreed by NCBI and Ensembl does not cover "
                     "every gene); the protein-coding transcripts are listed, longest coding sequence first")
            options = [Option(
                code=t.transcript_id,
                label=f"{t.transcript_id} ({'MANE Select, ' if t.is_mane_select else ''}{t.transcript_type.replace('_', '-')}, "
                      f"coding sequence {_cds_length(t):,} bp)",
                effect=f"The gene starts at {t.chromosome}:{_tss(t):,} (strand {'+' if _strand(t) == 1 else '-'}); "
                       "the window, exon 1 and the start codon follow this transcript.") for t in pool]
            options.append(Option(code="abort", label="Stop to check the gene",
                                  effect="Nothing is computed; check the gene symbol or capture the gene structure again."))
            raise LocusPending([_question(
                TRANSCRIPT_KEY,
                f"Which transcript of {gene_symbol} defines where the gene starts, its first exon and its start codon? "
                f"{basis.capitalize()}.{genes}",
                options, candidates[0] if candidates else "abort",
                title=f"Reference transcript of {gene_symbol}",
                why="The start of the transcript places the window and labels each variant (promoter, 5' UTR, coding); "
                    "another transcript can move the window and change those labels.",
                identifiers={"gene": gene_symbol, "gene IDs": ", ".join(gene_ids), "source": GTF_FILE})])
        value = str(entry.value)
        if value == "abort":
            raise LocusError(f"Stopped: the researcher answered 'abort' to {TRANSCRIPT_KEY}.")
        matches = [t for t in transcripts if t.transcript_id == value
                   or strip_version(t.transcript_id) == strip_version(value)]
        if len(matches) != 1:
            raise LocusError(f"Decision '{TRANSCRIPT_KEY}' = {value!r} names no single transcript of {GTF_FILE}.")
        chosen = matches[0]
    if not chosen.cds:
        raise LocusError(f"Transcript {chosen.transcript_id} ({chosen.transcript_type}) has no CDS in {GTF_FILE}; the "
                         "window and the variant annotation need a start codon (capture gtf with --cds --utr).")
    return chosen


# --------------------------------------------------------------------------- window policy


def derive_regulatory_unit_window(
    tx: GtfTranscript,
    features: Iterable[RegulatoryFeature],
    biotypes: Iterable[str],
    max_gap_bp: int,
    chromosome: str,
    decisions: Mapping[str, DecisionEntry] | None = None,
    rules: list[dict[str, Any]] | None = None,
) -> tuple[Interval | None, list[RegulatoryFeature], list[PendingDecision]]:
    """Applies the `regulatory_unit` policy to a normalised transcript.

    Args:
        decisions: Recorded answers. The chain stops at the first feature_overlap/feature_gap question
            without an answer (the next run asks the following one) and at an `exclude_feature` answer. None
            (the exploration) chains every feature it can, collecting every question.
        rules: When given, receives the choices the policy makes by definition ({key, value, rule}):
            several regions at the start of the gene are all kept (`locus.multiple_tss_features`), and a
            region that reaches into exon 1 keeps the whole of exon 1 (`locus.exon1_overlap:{feature}`).
            They are not questions: the researcher chose this policy, whose window is the regulatory
            regions at the start of the gene through the end of exon 1.

    Returns:
        (window, chained features, questions). With no feature at the TSS the window is None and
        the only question is `locus.no_feature_chain`.
    """
    exploring = decisions is None
    decisions = decisions or {}
    feat = lambda f: f"{f.id} ({f.description}, {chromosome}:{f.start:,}-{f.end:,})"
    s = _strand(tx)
    ex1 = _first_exon(tx)
    tss_o = _tss(tx) * s
    ex1_o = _oriented(ex1.start, ex1.end, s)
    wanted = set(biotypes)
    candidates = [f for f in features if f.description in wanted]
    oriented = {f.id: _oriented(f.start, f.end, s) for f in candidates}

    anchors = [f for f in candidates if oriented[f.id][0] <= tss_o <= oriented[f.id][1] + 1 + max_gap_bp]
    if not anchors:
        return None, [], [_question(
            NO_CHAIN_KEY,
            f"No {' or '.join(sorted(wanted))} region of the Ensembl Regulatory Build contains or touches the start of "
            f"{tx.transcript_id} ({chromosome}:{_tss(tx):,}), so the window cannot be built from regulatory regions. "
            "Which window should the study use?",
            [Option(code="set_explicit_interval", label="Use a window from the exploration page",
                    effect="You choose one of the exploration's windows that does not depend on regulatory regions "
                           "(for example the start of the gene +/- 500 bp); the agent sets it and the analysis continues "
                           "from this step."),
             STOP_WINDOW], "set_explicit_interval",
            title="No regulatory region at the start of the gene", why=WINDOW_WHY,
            identifiers={"transcript": tx.transcript_id, "gene start": f"{chromosome}:{_tss(tx):,}",
                         "regions searched": ", ".join(sorted(wanted))})]

    triggers: list[PendingDecision] = []
    found_rules: list[dict[str, Any]] = []
    chain = [min(anchors, key=lambda f: oriented[f.id][0])]
    if len(anchors) > 1:
        lo = min(f.start for f in anchors)
        hi = max(f.end for f in anchors)
        listed = "; ".join(feat(f) for f in sorted(anchors, key=lambda g: g.start))
        touching = all(oriented[f.id][0] <= tss_o <= oriented[f.id][1] + 1 for f in anchors)
        if touching:  # by definition of the window: every region at the start of the gene is part of it
            found_rules.append({"key": MULTIPLE_TSS_KEY, "value": "union", "rule": (
                f"{len(anchors)} regulatory regions contain or touch the start of the gene ({chromosome}:{_tss(tx):,}): "
                f"{listed}. All of them are kept, so the window starts from {chromosome}:{lo:,}-{hi:,}.")})
        else:  # some reach the start of the gene only through the allowed gap: a real choice
            triggers.append(_question(
                MULTIPLE_TSS_KEY,
                f"{len(anchors)} regulatory regions lie at or near the start of the gene ({chromosome}:{_tss(tx):,}), "
                f"some only within the {max_gap_bp:,} bp the study allows between regions: {listed}. Should the window "
                "start from all of them together?",
                [Option(code="union", label="Use all of them together",
                        effect=f"The core of the window spans {chromosome}:{lo:,}-{hi:,} ({_span(lo, hi)})."),
                 STOP_WINDOW], "union",
                title="Several regulatory regions near the start of the gene", why=WINDOW_WHY,
                identifiers={"regions": ", ".join(sorted(f.id for f in anchors)),
                             "gene start": f"{chromosome}:{_tss(tx):,}"}))
        chain = sorted(anchors, key=lambda f: oriented[f.id][0])

    used = {f.id for f in chain}
    while True:
        head_o = oriented[chain[0].id][0]
        nxt = [f for f in candidates
               if f.id not in used and oriented[f.id][0] < head_o and oriented[f.id][1] >= head_o - 1 - max_gap_bp]
        if not nxt:
            break
        f = min(nxt, key=lambda g: oriented[g.id][0])
        gap = head_o - oriented[f.id][1] - 1
        added = oriented[chain[0].id][0] - oriented[f.id][0]  # bases the window gains upstream with f
        if gap != 0:
            key = f"locus.feature_{'overlap' if gap < 0 else 'gap'}:{f.id}+{chain[0].id}"
            how = (f"overlaps {chain[0].id} by {-gap:,} bp" if gap < 0 else
                   f"is {gap:,} bp upstream of {chain[0].id} (regions up to {max_gap_bp:,} bp apart are joined)")
            include = Option(code="union" if gap < 0 else "bridge_gap", label=f"Include {f.id} in the window",
                             effect=f"The window grows by {added:,} bp upstream ({3 * added:,} more variants)"
                                    + ("." if gap < 0 else ", including the gap between the two regions."))
            exclude = Option(code="exclude_feature", label=f"Leave {f.id} out",
                             effect=f"The window starts at {chain[0].id}; {f.id} and any region beyond it are not scored.")
            triggers.append(_question(
                key, f"The regulatory region {feat(f)} {how}. Should it be part of the window?",
                [include, exclude, STOP_WINDOW], include.code,
                title=f"Extend the window to {f.id}?", why=WINDOW_WHY,
                identifiers={"region": feat(f), "next to": feat(chain[0])}))
            entry = decisions.get(key)
            if (entry is None and not exploring) or (entry is not None and str(entry.value) == "exclude_feature"):
                break  # the regions beyond depend on this answer: ask them only once it is known
        chain.insert(0, f)
        used.add(f.id)

    unit_lo = min(oriented[f.id][0] for f in chain)
    unit_hi = ex1_o[1]
    a, b = sorted((unit_lo * s, unit_hi * s))
    for f in chain:
        lo, hi = oriented[f.id]
        if hi >= ex1_o[0]:
            overlap = min(hi, ex1_o[1]) - ex1_o[0] + 1
            found_rules.append({"key": f"locus.exon1_overlap:{f.id}", "value": "keep_full_unit", "rule": (
                f"Region {feat(f)} overlaps the first {overlap:,} bp of exon 1 of {tx.transcript_id} "
                f"({chromosome}:{ex1.start:,}-{ex1.end:,}); the window keeps all of exon 1, and those bases count as both "
                f"{f.description} and exon 1. Window {chromosome}:{a:,}-{b:,} ({_span(a, b)}).")})
    if rules is not None:
        rules.extend(found_rules)
    return Interval(start=a, end=b), chain, triggers


def parse_explicit_interval(interval: str | None, chromosome: str) -> Interval:
    m = _INTERVAL_RE.match(interval or "")
    if m is None:
        raise LocusError(f"study.yaml window.explicit_interval {interval!r}: expected chr:start-end (1-based closed).")
    if m.group(1) != chromosome:
        raise LocusError(f"study.yaml window.explicit_interval {interval!r} is not on {chromosome}.")
    try:
        return Interval(start=int(m.group(2)), end=int(m.group(3)))
    except ValueError as exc:
        raise LocusError(f"study.yaml window.explicit_interval {interval!r}: {exc}") from None


# --------------------------------------------------------------------------- assembly


def build_locus(
    transcripts: list[GtfTranscript],
    regulatory: list[RegulatoryFeature],
    study: StudyConfig,
    decisions: dict[str, DecisionEntry],
    rules: list[dict[str, Any]] | None = None,
) -> LocusConfig:
    """Derives the locus and window of the study gene.

    `rules`, when given, receives the window choices made by the policy's definition (recorded by `run`).

    Args:
        transcripts: Records of raw/gtf.json AS LOADED (0-based starts); normalised here once.
        regulatory: Records of raw/regulatory.json.
        study: Validated study.yaml.
        decisions: decisions.yaml.

    Raises:
        LocusPending: A researcher question has no recorded answer (all open window questions are
            raised together; the transcript question first, since the window depends on it).
        LocusError: Inconsistent input, an `abort` answer, or a policy the inputs cannot satisfy.
    """
    normalised = [models.normalize_gtf(t) for t in transcripts]
    tx = select_transcript(normalised, study.gene_symbol, decisions)
    strand = _strand(tx)
    ex1 = _first_exon(tx)
    tss = _tss(tx)
    chromosome = tx.chromosome
    contig = _contig_without_prefix(chromosome)
    on_contig = [f for f in regulatory if f.seq_region_name == contig]

    policy_cfg = study.data.get("window") or {}
    policy = policy_cfg.get("policy")
    chosen: dict[str, Any] = {}
    if policy == "explicit":
        window = parse_explicit_interval(policy_cfg.get("explicit_interval"), chromosome)
    elif policy == "regulatory_unit":
        found: list[dict[str, Any]] = []
        window, _chain, triggers = derive_regulatory_unit_window(
            tx, on_contig, policy_cfg.get("feature_biotypes") or DEFAULT_BIOTYPES,
            int(policy_cfg.get("max_feature_gap_bp") or 0), chromosome, decisions, found)
        chosen = _resolve(triggers, decisions)
        for r in found:  # an older study may hold the researcher's answer to what is now a rule
            entry = decisions.get(r["key"])
            if entry is not None and entry.decided_by == "researcher":
                if str(entry.value) == "abort":
                    raise LocusError(f"Stopped: the researcher answered 'abort' to {r['key']}. Change the window in "
                                     "study.yaml or record another answer.")
                chosen[r["key"]] = str(entry.value)
            else:
                chosen[r["key"]] = r["value"]
        if rules is not None:
            rules.extend(found)
        if window is None:
            raise LocusError(f"Decision '{NO_CHAIN_KEY}' = set_explicit_interval: set window.policy: explicit and "
                             "window.explicit_interval in study.yaml, then run expression.py locus again.")
    else:
        raise LocusError(f"study.yaml window.policy {policy!r} is not one of {list(POLICIES)}.")
    if sum(t.is_mane_select for t in normalised) != 1:  # select_transcript used the recorded answer
        chosen = {TRANSCRIPT_KEY: str(decisions[TRANSCRIPT_KEY].value), **chosen}
    if not window.contains(tss):
        raise LocusError(f"The window {chromosome}:{window.start}-{window.end} does not contain the TSS "
                         f"{chromosome}:{tss} of {tx.transcript_id}; widen window.explicit_interval.")

    cds_start = min(seg.start for seg in tx.cds) if strand == 1 else max(seg.end for seg in tx.cds)
    start_codon = (cds_start, cds_start + 2) if strand == 1 else (cds_start - 2, cds_start)
    in_window = sorted(
        (FeatureRef(id=f.id, biotype=f.description, start=f.start, end=f.end)
         for f in on_contig if f.end >= window.start and f.start <= window.end),
        key=lambda f: (f.start, f.end, f.id))
    return LocusConfig(
        study_id=study.study_id,
        gene_id=strip_version(tx.gene_id),
        gene_version=_version(tx.gene_id),
        gene_symbol=study.gene_symbol,
        transcript_id=strip_version(tx.transcript_id),
        transcript_version=_version(tx.transcript_id),
        transcript_type=tx.transcript_type,
        is_mane_select=tx.is_mane_select,
        assembly=study.assembly,
        chromosome=chromosome,
        strand=strand,
        tss=tss,
        cds_start_genomic=cds_start,
        start_codon=start_codon,
        first_exon=ExonRef(exon_number=ex1.exon_number, start=ex1.start, end=ex1.end),
        window=window,
        window_policy=policy,
        features=tuple(in_window),
        decisions=chosen,
    )


def check_sequence(locus: LocusConfig, sequences: list[tuple[str, EnsemblSequence]]) -> LocusConfig:
    """Records the raw sequence that covers the window and checks that the start codon reads ATG.

    Args:
        locus: Built locus.
        sequences: (study-relative file name, record) of every raw/sequence_*.txt; the interval of
            each record is the one in its file name.

    Returns:
        The locus with `sequence_file` and `start_codon_verified` set.

    Raises:
        LocusError: No file covers the window, or the start codon does not read ATG (which is how
            a 1-bp coordinate shift or a 5'-incomplete CDS shows up).
    """
    same = [(name, s) for name, s in sequences
            if (s.chromosome if s.chromosome.startswith("chr") else "chr" + s.chromosome) == locus.chromosome]
    covering = [name for name, s in same if s.start <= locus.window.start and locus.window.end <= s.end]
    if not covering:
        have = [name for name, _ in sequences] or "none"
        raise LocusError(f"No raw/{SEQUENCE_GLOB} covers the window {locus.window_str} (found: {have}); run Workflow "
                         "step 3 (ensembl-database get-sequence) for a region that contains it.")
    lo, hi = locus.start_codon
    verified = False
    for name, s in same:
        if s.start <= lo and hi <= s.end:
            codon = s.sequence[lo - s.start:hi - s.start + 1].upper()
            if locus.strand == -1:
                codon = reverse_complement(codon)
            if codon != "ATG":
                raise LocusError(f"The start codon of {locus.transcript_id} at {locus.chromosome}:{lo}-{hi} reads "
                                 f"{codon!r} in {name}, expected 'ATG' (coordinate shift, or a 5'-incomplete CDS: "
                                 f"choose another transcript with {TRANSCRIPT_KEY}).")
            verified = True
            break
    return dataclasses.replace(locus, sequence_file=covering[0], start_codon_verified=verified)


def transcripts_in_window(transcripts: list[GtfTranscript], locus: LocusConfig) -> list[dict[str, Any]]:
    """Transcripts of the selected gene whose TSS lies in the window, plus the selected one.

    Args:
        transcripts: Records of raw/gtf.json AS LOADED (normalised here).
        locus: Built locus.

    Returns:
        Rows with TRANSCRIPT_COLUMNS; `tss_offset_bp` is the strand-aware distance to the selected
        TSS (negative = upstream). Sorted by offset, the selected transcript first on ties.
    """
    rows = []
    for raw in transcripts:
        t = models.normalize_gtf(raw)
        if strip_version(t.gene_id) != locus.gene_id or not t.exons:
            continue
        tss = _tss(t)
        tid = strip_version(t.transcript_id)
        selected = tid == locus.transcript_id
        if not (selected or locus.window.contains(tss)):
            continue
        first = next((e for e in t.exons if e.exon_number == 1), None)
        rows.append({
            "transcript_id": tid,
            "transcript_version": _version(t.transcript_id),
            "transcript_type": t.transcript_type,
            "is_selected": selected,
            "is_mane_select": t.is_mane_select,
            "strand": _strand(t),
            "tss": tss,
            "tss_offset_bp": locus.distance_to_tss(tss),
            "first_exon_start": first.start if first else None,
            "first_exon_end": first.end if first else None,
        })
    return sorted(rows, key=lambda r: (r["tss_offset_bp"], not r["is_selected"], r["transcript_id"]))


# --------------------------------------------------------------------------- position annotation


def in_cds(locus: LocusConfig, pos: int) -> bool:
    """Position lies in exon 1 at or after the CDS start (transcript orientation)."""
    return (locus.first_exon.contains(pos)
            and locus.distance_to_tss(pos) >= locus.distance_to_tss(locus.cds_start_genomic))


def in_start_codon(locus: LocusConfig, pos: int) -> bool:
    lo, hi = locus.start_codon
    return lo <= pos <= hi


def regulatory_feature_at(locus: LocusConfig, pos: int) -> FeatureRef | None:
    """Regulatory Build feature containing `pos`; the shortest one wins when features overlap."""
    containing = [f for f in locus.features if f.contains(pos)]
    return min(containing, key=lambda f: (f.length, f.id)) if containing else None


def classify_block(locus: LocusConfig, pos: int) -> tuple[str, str | None]:
    """Footprint kind and feature ID: regulatory biotype, `5_prime_utr`, `cds`, `unannotated` or `other`."""
    if locus.first_exon.contains(pos):
        return ("cds" if in_cds(locus, pos) else "5_prime_utr", None)
    if locus.upstream_of_tss(pos):
        feature = regulatory_feature_at(locus, pos)
        return (feature.biotype, feature.id) if feature else ("unannotated", None)
    return ("other", None)


def transcript_region(locus: LocusConfig, pos: int) -> str:
    """Sequence Ontology-style consequence label for an SNV at `pos`."""
    if locus.first_exon.contains(pos):
        if in_start_codon(locus, pos):
            return "start_codon_variant"
        return "coding_sequence_variant" if in_cds(locus, pos) else "5_prime_UTR_variant"
    if locus.upstream_of_tss(pos):
        return "upstream_gene_variant"
    return "other"


def _display(biotype: str) -> str:
    text = biotype.replace("_", " ")
    return text[:1].upper() + text[1:]


def genomic_domain(locus: LocusConfig, pos: int) -> str:
    """Single display label shared by every table and figure."""
    if locus.first_exon.contains(pos):
        if in_start_codon(locus, pos):
            return "ATG Start Codon"
        return "Exon 1 (CDS)" if in_cds(locus, pos) else "Exon 1 (5' UTR)"
    if locus.upstream_of_tss(pos):
        feature = regulatory_feature_at(locus, pos)
        return f"{feature.id} ({_display(feature.biotype)})" if feature else "Unannotated"
    return "Other"


def footprint_blocks(locus: LocusConfig) -> list[FootprintBlock]:
    """Partitions the window into contiguous blocks in genomic order; lengths sum to the window length."""
    blocks: list[FootprintBlock] = []
    run_start = locus.window.start
    current = classify_block(locus, run_start)
    for pos in range(locus.window.start + 1, locus.window.end + 2):
        nxt = classify_block(locus, pos) if pos <= locus.window.end else None
        if nxt != current:
            blocks.append(FootprintBlock(kind=current[0], feature_id=current[1], start=run_start, end=pos - 1))
            run_start, current = pos, nxt
    return blocks


def snv_counts(locus: LocusConfig) -> dict[str, int]:
    """SNVs (3 per base) per footprint kind, named as in METRIC_KIND, plus `total`."""
    counts = {k: 0 for k in STANDARD_METRIC_KINDS}
    for block in footprint_blocks(locus):
        name = METRIC_KIND.get(block.kind, block.kind)
        counts[name] = counts.get(name, 0) + block.snvs
    counts["total"] = 3 * locus.window.length
    return counts


# --------------------------------------------------------------------------- stage I/O


def _load(study_dir: Path, relative: str, step: int, producer: str) -> list[Any]:
    path = study_dir / relative
    if not path.is_file():
        raise LocusError(f"{relative} not found in {study_dir}; run Workflow step {step} ({producer}).")
    return models.load_raw(path)


def _load_sequences(study_dir: Path) -> list[tuple[str, EnsemblSequence]]:
    paths = sorted((study_dir / "raw").glob(SEQUENCE_GLOB))
    if not paths:
        raise LocusError(f"No raw/{SEQUENCE_GLOB} in {study_dir}; run Workflow step 3 (ensembl-database "
                         "get-sequence for the regulatory search region).")
    return [(str(p.relative_to(study_dir)), models.load_raw(p)[0]) for p in paths]


def _write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    def cell(v: Any) -> str:
        return "" if v is None else ("true" if v is True else "false" if v is False else str(v))

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(TRANSCRIPT_COLUMNS)
        for row in rows:
            writer.writerow([cell(row[c]) for c in TRANSCRIPT_COLUMNS])


def run(study_dir: Path) -> LocusConfig:
    """Builds the locus from raw/ and decisions.yaml; writes locus.json and transcripts_in_window.tsv.

    Raises:
        LocusPending, LocusError, models.DependencyOutputError, study.StudyError.
    """
    study_dir = Path(study_dir)
    study = load_study(study_dir)
    decisions = load_decisions(study_dir)
    transcripts = _load(study_dir, GTF_FILE, 2, "alphagenome-variant-impact-score gtf")
    regulatory = _load(study_dir, REGULATORY_FILE, 3, "regulatory_features.py overlap")
    sequences = _load_sequences(study_dir)
    try:
        rules: list[dict[str, Any]] = []
        locus = build_locus(transcripts, regulatory, study, decisions, rules)
    except ValueError as exc:  # normalize_gtf or a LocusConfig invariant
        raise LocusError(f"Cannot build the locus: {exc}") from None
    sync_rules(study_dir, (MULTIPLE_TSS_KEY, "locus.exon1_overlap:"), rules, WORKFLOW_STEP)
    locus = check_sequence(locus, sequences)
    (study_dir / "processed").mkdir(parents=True, exist_ok=True)
    (study_dir / LOCUS_FILE).write_text(json.dumps(locus.to_dict(), indent=2) + "\n", encoding="utf-8")
    _write_tsv(study_dir / TRANSCRIPTS_FILE, transcripts_in_window(transcripts, locus))
    return locus


def load_locus(study_dir: Path) -> LocusConfig:
    """Reads processed/locus.json (written by the locus step)."""
    path = Path(study_dir) / LOCUS_FILE
    if not path.is_file():
        raise LocusError(f"{LOCUS_FILE} not found in {study_dir}; run expression.py locus (Workflow step {WORKFLOW_STEP}).")
    try:
        return LocusConfig.from_dict(json.loads(path.read_text(encoding="utf-8")))
    except (json.JSONDecodeError, KeyError, TypeError, ValueError) as exc:
        raise LocusError(f"{path}: invalid locus ({exc}); run expression.py locus again.") from None


def load_transcripts_in_window(study_dir: Path) -> list[dict[str, str]]:
    """Rows of processed/transcripts_in_window.tsv as text, or [] when the file is absent."""
    path = Path(study_dir) / TRANSCRIPTS_FILE
    if not path.is_file():
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def stage_metrics(study_dir: Path) -> dict[str, int | float | str]:
    """Metrics of the locus stage compared by `reproduce` (expected.json stages.locus)."""
    locus = load_locus(study_dir)
    metrics: dict[str, int | float | str] = {
        "transcript_id": locus.transcript_id,
        "strand": locus.strand,
        "tss": locus.tss,
        "cds_start": locus.cds_start_genomic,
        "window_start": locus.window.start,
        "window_end": locus.window.end,
        "window_length": locus.window.length,
        "n_regulatory_in_window": len(locus.features),
    }
    metrics.update({f"snvs.{kind}": n for kind, n in snv_counts(locus).items()})
    return metrics
