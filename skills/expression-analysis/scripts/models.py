#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["pyarrow==26.0.1"]
# ///
"""Typed views of the dependency outputs kept in a study's raw/ folder.

Every raw file is parsed by the dataclass of its file-name pattern (RAW_PATTERNS). Field names
are the keys (or columns) the dependency writes, so a renamed, added or dropped key raises
DependencyOutputError(file, field, skill) instead of silently changing an analysis. Keys a
dependency legitimately omits are listed in each class's OPTIONAL set and documented there.
Coordinates are reported exactly as the dependency wrote them; the helpers at the end convert
gnomAD `chrom-pos-ref-alt` ids and ClinVar `allele_info` to 1-based `chr:pos:REF>ALT`.
See references/dependency_outputs.md for the producing command of each file.

Subcommands:
  validate --dir D   Parse every file in D; exit 1 and print file, field and skill on the
                     first problem of each file.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import json
import math
import re
import sys
import types
from collections.abc import Callable, Mapping
from pathlib import Path
from typing import Any, ClassVar

SKILL_AVI = "alphagenome-variant-impact-score"
SKILL_KDENSE = "alphagenome"
SKILL_OLS = "embl-ebi-ols"
SKILL_ENSEMBL = "ensembl-database"
SKILL_REGULATORY = "expression-analysis (regulatory_features.py, Ensembl REST)"
SKILL_CLINVAR = "clinvar-database"
SKILL_GNOMAD = "gnomad-database"
SKILL_SVA = "alphagenome-single-variant-analysis"
SKILL_LINKS = "alphagenome-atlas-website-links"
SKILL_SNAPSHOT = "splicing-analysis snapshot (2026-10)"
SKILL_DBSNP = "dbsnp-database"

# The 18 AVI feature attributions (fi_<NAME>, track_*_<NAME>) of alphagenome_atlas_avi.py.
AVI_FEATURES = (
    "MERGED_SPLICING", "MAX_ABS_ATAC", "MAX_ABS_CONTACT_MAPS", "MAX_ABS_DNASE",
    "MAX_ABS_CHIP_TF", "MAX_ABS_CHIP_HISTONE", "MAX_ABS_CAGE", "MAX_ABS_PROCAP",
    "MAX_ABS_RNA_SEQ", "MAX_ABS_POLYADENYLATION", "ALPHAMISSENSE", "CACTUS_241_WAY",
    "PROTEIN_TERMINATION", "START_LOST", "STOP_LOST", "PHASTCONS_470_WAY",
    "IS_INSERTION", "IS_DELETION",
)


class DependencyOutputError(Exception):
    """A dependency output does not have the expected shape.

    Attributes:
      file: Path of the offending file.
      field: Key, column or line that failed (nested keys joined with '.', records as [i]).
      skill: Dependency skill that produced the file.
      reason: Short explanation.
    """

    def __init__(self, file: str, field: str, skill: str, reason: str = ""):
        self.file, self.field, self.skill, self.reason = str(file), str(field), skill, reason
        super().__init__(f"{self.file}: field '{self.field}' ({self.skill}): {reason or 'invalid'}")


# --------------------------------------------------------------------------------------------
# Value converters. Each takes (value, source) and raises ValueError/TypeError on a bad type.
# Text inputs (TSV/CSV) arrive as strings; an empty cell is a missing value.
# --------------------------------------------------------------------------------------------

def _is_missing(v: Any) -> bool:
    return v is None or (isinstance(v, str) and v == "") or (isinstance(v, float) and math.isnan(v))


def _str(v: Any, source: str = "") -> str:
    if not isinstance(v, str):
        raise TypeError(f"expected text, got {type(v).__name__}")
    return v


def _ostr(v: Any, source: str = "") -> str | None:
    return None if _is_missing(v) else _str(v)


def _text(v: Any, source: str = "") -> str:
    """Required, non-empty text (an empty TSV/CSV cell or null is a missing value)."""
    if _is_missing(v):
        raise ValueError("missing value")
    return _str(v)


def _int(v: Any, source: str = "") -> int:
    if isinstance(v, bool):
        raise TypeError("expected an integer, got a boolean")
    if isinstance(v, int):
        return v
    if isinstance(v, float) and v.is_integer():
        return int(v)
    if isinstance(v, str) and re.fullmatch(r"-?\d+(\.0+)?", v.strip()):
        return int(float(v))
    raise TypeError(f"expected an integer, got {v!r}")


def _oint(v: Any, source: str = "") -> int | None:
    return None if _is_missing(v) else _int(v)


def _float(v: Any, source: str = "") -> float:
    """Required number; null, an empty cell and NaN are missing values."""
    if _is_missing(v) or (isinstance(v, str) and v.strip().lower() == "nan"):
        raise ValueError("missing value")
    if isinstance(v, bool):
        raise TypeError("expected a number, got a boolean")
    if isinstance(v, (int, float)):
        return float(v)
    if isinstance(v, str):
        try:
            return float(v)
        except ValueError as exc:
            raise TypeError(f"expected a number, got {v!r}") from exc
    raise TypeError(f"expected a number, got {type(v).__name__}")


def _ofloat(v: Any, source: str = "") -> float | None:
    """Optional number; null, an empty cell and NaN become None."""
    if _is_missing(v) or (isinstance(v, str) and v.strip().lower() == "nan"):
        return None
    return _float(v)


def _bool(v: Any, source: str = "") -> bool:
    if isinstance(v, bool):
        return v
    if isinstance(v, str) and v in ("True", "False", "true", "false"):
        return v.lower() == "true"
    raise TypeError(f"expected a boolean, got {v!r}")


def _str_tuple(v: Any, source: str = "") -> tuple[str, ...]:
    if not isinstance(v, list):
        raise TypeError(f"expected a list, got {type(v).__name__}")
    return tuple(_str(x) for x in v)


def _ostr_tuple(v: Any, source: str = "") -> tuple[str, ...] | None:
    return None if v is None else _str_tuple(v)


def _mapping(v: Any, source: str = "") -> Mapping[str, Any]:
    if not isinstance(v, dict):
        raise TypeError(f"expected an object, got {type(v).__name__}")
    return types.MappingProxyType(dict(v))


def _omapping(v: Any, source: str = "") -> Mapping[str, Any] | None:
    return None if v is None else _mapping(v)


def _mappings(v: Any, source: str = "") -> tuple[Mapping[str, Any], ...]:
    """A list of objects kept as read-only mappings (nested content not modelled field by field)."""
    if not isinstance(v, list):
        raise TypeError(f"expected a list, got {type(v).__name__}")
    return tuple(_mapping(x) for x in v)


def _many(cls: type) -> Callable[[Any, str], tuple]:
    def convert(v: Any, source: str) -> tuple:
        if not isinstance(v, list):
            raise TypeError(f"expected a list, got {type(v).__name__}")
        out = []
        for i, item in enumerate(v):
            try:
                out.append(cls.from_dict(item, source=source))
            except DependencyOutputError as exc:
                raise DependencyOutputError(source, f"[{i}].{exc.field}", exc.skill, exc.reason) from None
        return tuple(out)
    return convert


def _one(cls: type, *, nullable: bool = False) -> Callable[[Any, str], Any]:
    def convert(v: Any, source: str) -> Any:
        if v is None and nullable:
            return None
        return cls.from_dict(v, source=source)
    return convert


class _Record:
    """Strict mapping-to-dataclass constructor shared by the raw-file dataclasses.

    Subclasses set SKILL (the producing dependency), CONVERT (field -> converter),
    OPTIONAL (file keys the dependency may omit) and, rarely, ALIASES (file key -> field
    name, for keys that are not Python identifiers).
    """

    SKILL: ClassVar[str] = ""
    CONVERT: ClassVar[Mapping[str, Callable[[Any, str], Any]]] = {}
    OPTIONAL: ClassVar[frozenset[str]] = frozenset()
    ALIASES: ClassVar[Mapping[str, str]] = {}

    @classmethod
    def _required(cls, d: Mapping[str, Any]) -> set[str]:
        return set(cls._file_keys()) - set(cls.OPTIONAL)

    @classmethod
    def _file_keys(cls) -> dict[str, str]:
        inverse = {field: key for key, field in cls.ALIASES.items()}
        return {inverse.get(name, name): name for name in cls.CONVERT}

    @classmethod
    def from_dict(cls, d: Mapping[str, Any], *, source: str):
        """Builds the record from one mapping; raises DependencyOutputError on any mismatch."""
        if not isinstance(d, Mapping):
            raise DependencyOutputError(source, cls.__name__, cls.SKILL,
                                        f"expected an object, got {type(d).__name__}")
        keys = cls._file_keys()
        unknown = sorted(set(d) - set(keys))
        if unknown:
            raise DependencyOutputError(source, unknown[0], cls.SKILL, "unknown field")
        missing = sorted(cls._required(d) - set(d))
        if missing:
            raise DependencyOutputError(source, missing[0], cls.SKILL, "missing field")
        kwargs = {}
        for key, name in keys.items():
            if key not in d:
                continue
            try:
                kwargs[name] = cls.CONVERT[name](d[key], source)
            except DependencyOutputError as exc:
                field = exc.field if exc.field.startswith("[") else f".{exc.field}"
                raise DependencyOutputError(source, f"{key}{field}", exc.skill, exc.reason) from None
            except (TypeError, ValueError) as exc:
                raise DependencyOutputError(source, key, cls.SKILL, str(exc)) from None
        return cls(**kwargs)


def _record(cls):
    return dataclasses.dataclass(frozen=True, kw_only=True)(cls)


# --------------------------------------------------------------------------------------------
# alphagenome-variant-impact-score: gtf (gtf.json, gtf_junctions.json)
# --------------------------------------------------------------------------------------------

@_record
class GtfSegment(_Record):
    """One CDS or UTR interval of a `gtf --format json` transcript.

    Attributes:
      start: GENCODE v46 Feather `Start` as written, 0-based (normalize_gtf gives 1-based).
      end: Last base, 1-based.
      width: end - start + 1 as written (one more than the length; see normalize_gtf).
    """

    start: int
    end: int
    width: int
    SKILL: ClassVar[str] = SKILL_AVI
    CONVERT: ClassVar = {"start": _int, "end": _int, "width": _int}


@_record
class GtfExon(_Record):
    """One exon of a `gtf --exons --format json` transcript.

    Attributes:
      exon_number: Exon index in 5'->3' transcript order.
      chromosome: Contig with chr prefix.
      start: Feather `Start` as written, 0-based (normalize_gtf gives 1-based).
      end: Last base, 1-based.
      width: end - start + 1 as written (one more than the length).
      strand: '+' or '-'.
      acceptor: 5' exon boundary: start on '+' (0-based as written), end on '-'.
      donor: 3' exon boundary: end on '+', start on '-' (0-based as written).
      overlaps_variant: True when a --variant query falls inside the exon.
    """

    exon_number: int
    chromosome: str
    start: int
    end: int
    width: int
    strand: str
    acceptor: int
    donor: int
    overlaps_variant: bool
    SKILL: ClassVar[str] = SKILL_AVI
    CONVERT: ClassVar = {"exon_number": _int, "chromosome": _str, "start": _int, "end": _int,
                         "width": _int, "strand": _str, "acceptor": _int, "donor": _int,
                         "overlaps_variant": _bool}


@_record
class GtfJunction(_Record):
    """One canonical intron or exon-skipping junction of a `gtf` transcript.

    Attributes:
      type: 'Canonical Intron {i}' or 'Exon {k} Skipping'.
      upstream_exon: Flanking exon number, 5' side.
      downstream_exon: Flanking exon number, 3' side.
      junction_start: Last base of the upstream exon, 1-based (= 0-based intron start; Atlas encoding).
      junction_end: Last base of the intron, 1-based (= 0-based half-open intron end).
      junction_str: 'junction_start:junction_end'.
      intron_length: junction_end - junction_start - 1 as written (one less than the intron length).
      score_id: Atlas ScoreId (only for canonical introns of a --variant query, else null).
      atlas_url: Atlas track-predictions link (same condition as score_id).

    Optional keys: score_id and atlas_url are absent on 'Exon {k} Skipping' entries.
    """

    type: str
    upstream_exon: int
    downstream_exon: int
    junction_start: int
    junction_end: int
    junction_str: str
    intron_length: int
    score_id: str | None = None
    atlas_url: str | None = None
    SKILL: ClassVar[str] = SKILL_AVI
    CONVERT: ClassVar = {"type": _str, "upstream_exon": _int, "downstream_exon": _int,
                         "junction_start": _int, "junction_end": _int, "junction_str": _str,
                         "intron_length": _int, "score_id": _ostr, "atlas_url": _ostr}
    OPTIONAL: ClassVar = frozenset({"score_id", "atlas_url"})


@_record
class GtfTranscript(_Record):
    """One transcript of `gtf --format json` (top-level records of gtf.json and gtf_junctions.json).

    Attributes:
      gene_name: HGNC symbol.
      gene_id: Versioned Ensembl gene ID (GENCODE v46).
      transcript_id: Versioned Ensembl transcript ID.
      transcript_type: GENCODE transcript biotype.
      chromosome: Contig with chr prefix.
      start: Feather transcript `Start` as written, 0-based (normalize_gtf gives 1-based).
      end: Last base, 1-based.
      strand: '+' or '-'.
      is_mane_select: True for the MANE Select transcript.
      num_exons: Number of exons.
      exons: Exons in 5'->3' order.
      junctions: Canonical introns, then exon-skipping junctions (always written).
      cds: CDS intervals (empty unless the transcript is coding).
      utr: UTR intervals, 5' and 3' together.
    """

    gene_name: str
    gene_id: str
    transcript_id: str
    transcript_type: str
    chromosome: str
    start: int
    end: int
    strand: str
    is_mane_select: bool
    num_exons: int
    exons: tuple[GtfExon, ...]
    junctions: tuple[GtfJunction, ...]
    cds: tuple[GtfSegment, ...]
    utr: tuple[GtfSegment, ...]
    SKILL: ClassVar[str] = SKILL_AVI
    CONVERT: ClassVar = {"gene_name": _str, "gene_id": _str, "transcript_id": _str,
                         "transcript_type": _str, "chromosome": _str, "start": _int, "end": _int,
                         "strand": _str, "is_mane_select": _bool, "num_exons": _int,
                         "exons": _many(GtfExon), "junctions": _many(GtfJunction),
                         "cds": _many(GtfSegment), "utr": _many(GtfSegment)}


# --------------------------------------------------------------------------------------------
# alphagenome-variant-impact-score: region / query (avi_region.parquet, avi_query.json)
# --------------------------------------------------------------------------------------------

@_record
class AviTrackInfo(_Record):
    """Track behind one AVI feature (`track_{idx,name,biosample,gene}_<FEATURE>` columns).

    Attributes:
      idx: Atlas track index (track_idx_<FEATURE>).
      name: Track name (track_name_<FEATURE>).
      biosample: Biosample name (track_biosample_<FEATURE>).
      gene: Gene symbol (track_gene_<FEATURE>); absent for features without a gene.

    Optional keys: gene (written only when the attributed track has a gene).
    """

    idx: int
    name: str
    biosample: str
    gene: str | None = None
    SKILL: ClassVar[str] = SKILL_AVI
    CONVERT: ClassVar = {"idx": _int, "name": _str, "biosample": _str, "gene": _ostr}
    OPTIONAL: ClassVar = frozenset({"gene"})


_AVI_CORE = {"rank": _oint, "variant": _text, "chromosome": _text, "position": _int, "ref": _text,
             "alt": _text, "avi_phred": _float, "avi_raw": _float, "avi_quantile": _float,
             "top_percentile": _float, "top_modality": _text, "top_feature_importance": _float}
_AVI_TRACK_RE = re.compile(r"^track_(idx|name|biosample|gene)_(.+)$")


@_record
class AviRecord(_Record):
    """One scored SNV of `region` (avi_region.parquet) or `query --format json` (avi_query.json).

    Attributes:
      rank: 1-based rank by Phred (region exports only).
      variant: chr:pos:REF>ALT, 1-based.
      chromosome: Contig with chr prefix.
      position: 1-based position.
      ref: Reference base.
      alt: Alternate base.
      avi_phred: Calibrated Phred score.
      avi_raw: Raw AVI model output.
      avi_quantile: Tail quantile (1 - CDF).
      top_percentile: Genome-wide top percentile.
      top_modality: Display name of the top feature.
      top_feature_importance: Attribution of the top feature.
      fi: The 18 `fi_<FEATURE>` attributions, keyed by FEATURE (all 18 required).
      track: `track_*_<FEATURE>` provenance keyed by FEATURE (query with --include_track_info).

    Optional keys: rank (absent in query output) and every track_* column (region never writes
    them; query writes them only with --include_track_info).

    Required values: avi_phred, avi_raw, avi_quantile, top_percentile, top_modality,
    top_feature_importance and every fi_* must hold a value (null, empty or NaN fails). The script
    can leave them null when the Atlas returns no quantiles or no attributions, but none of the
    outputs recorded on 2026-11-12 does (the expression reference avi_region.parquet and the
    splicing live-capture avi_query.json), so a null there is reported instead of passed on.
    """

    rank: int | None = None
    variant: str
    chromosome: str
    position: int
    ref: str
    alt: str
    avi_phred: float
    avi_raw: float
    avi_quantile: float
    top_percentile: float
    top_modality: str
    top_feature_importance: float
    fi: Mapping[str, float]
    track: Mapping[str, AviTrackInfo]
    SKILL: ClassVar[str] = SKILL_AVI

    @classmethod
    def from_dict(cls, d: Mapping[str, Any], *, source: str) -> AviRecord:
        if not isinstance(d, Mapping):
            raise DependencyOutputError(source, cls.__name__, cls.SKILL, "expected an object")
        kwargs: dict[str, Any] = {}
        fi: dict[str, float] = {}
        tracks: dict[str, dict[str, Any]] = {}
        for key, value in d.items():
            try:
                if key in _AVI_CORE:
                    kwargs[key] = _AVI_CORE[key](value)
                elif key.startswith("fi_") and key[3:] in AVI_FEATURES:
                    fi[key[3:]] = _float(value)
                elif (m := _AVI_TRACK_RE.match(key)) and m.group(2) in AVI_FEATURES:
                    tracks.setdefault(m.group(2), {})[m.group(1)] = value
                else:
                    raise DependencyOutputError(source, key, cls.SKILL, "unknown field")
            except (TypeError, ValueError) as exc:
                raise DependencyOutputError(source, key, cls.SKILL, str(exc)) from None
        missing = sorted({k for k in _AVI_CORE if k != "rank"} - set(kwargs))
        missing += [f"fi_{f}" for f in AVI_FEATURES if f not in fi]
        if missing:
            raise DependencyOutputError(source, missing[0], cls.SKILL, "missing field")
        track = {}
        for feature, parts in tracks.items():
            parts = {k: v for k, v in parts.items() if not (k == "gene" and _is_missing(v))}
            if all(_is_missing(v) for v in parts.values()):
                continue  # null cells of a feature without track info in a table export
            try:
                track[feature] = AviTrackInfo.from_dict(parts, source=source)
            except DependencyOutputError as exc:
                raise DependencyOutputError(source, f"track_{exc.field}_{feature}", cls.SKILL, exc.reason) from None
        return cls(**kwargs, fi=types.MappingProxyType(fi), track=types.MappingProxyType(track))


# --------------------------------------------------------------------------------------------
# alphagenome-variant-impact-score: metadata --tracks (tracks_{SCORER}_{CURIE}.parquet)
# --------------------------------------------------------------------------------------------

@_record
class AtlasTrack(_Record):
    """One Atlas catalogue track of `metadata --tracks --ontology C --scorer S --output`.

    Attributes:
      scorer: Atlas scorer name.
      track_index: Track index within the scorer.
      track_name: Track name.
      biosample_name: Biosample name.
      biosample_type: tissue, cell_line, primary_cell, ...
      ontology_curie: UBERON/CL/EFO CURIE.
      strand: '+', '-' or '.'.
      transcription_factor: ChIP-TF target (CHIP_TF tracks only).
      histone_mark: Histone mark (CHIP_HISTONE tracks only).
      assay: Assay label (only when the catalogue row has one).

    Optional keys: transcription_factor, histone_mark and assay, written only when the scorer's
    catalogue has that column.
    """

    scorer: str
    track_index: int
    track_name: str
    biosample_name: str
    biosample_type: str
    ontology_curie: str
    strand: str
    transcription_factor: str | None = None
    histone_mark: str | None = None
    assay: str | None = None
    SKILL: ClassVar[str] = SKILL_AVI
    CONVERT: ClassVar = {"scorer": _str, "track_index": _int, "track_name": _str,
                         "biosample_name": _str, "biosample_type": _str, "ontology_curie": _str,
                         "strand": _str, "transcription_factor": _ostr, "histone_mark": _ostr,
                         "assay": _ostr}
    OPTIONAL: ClassVar = frozenset({"transcription_factor", "histone_mark", "assay"})


# --------------------------------------------------------------------------------------------
# alphagenome (K-Dense): atlas_query.py scores (scores_*.tsv) and score_variants.py (indels)
# --------------------------------------------------------------------------------------------

@_record
class TrackScoreRow(_Record):
    """One row of `atlas_query.py scores` (variant x track, x gene or junction where applicable).

    Attributes:
      variant: chr:pos:REF>ALT, 1-based.
      scorer: Atlas scorer name (empty on error rows).
      gene_id: Versioned gene ID (gene-centric scorers).
      gene_name: Gene symbol (gene-centric scorers).
      gene_strand: Gene strand (gene-centric scorers).
      junction_start: Junction start as served (SPLICE_JUNCTIONS).
      junction_end: Junction end as served (SPLICE_JUNCTIONS).
      track_index: Track index within the scorer.
      track_name: Track name.
      track_strand: Track strand.
      ontology_curie: Track ontology CURIE.
      biosample_name: Biosample name.
      biosample_type: Biosample type.
      assay_title: Assay title.
      gtex_tissue: GTEx tissue (GTEx tracks).
      transcription_factor: ChIP-TF target (CHIP_TF).
      histone_mark: Histone mark (CHIP_HISTONE).
      data_source: encode, gtex, fantom, ...
      raw_score: Raw score on the scorer's scale.
      quantile_score: Calibrated quantile when served.
      error: Per-variant lookup failure (--variant/--input runs only).

    Optional keys: every column except variant and scorer; the script writes the union of the
    track-metadata, gene and junction columns the server returned. track_index and raw_score
    are required unless the file has an error column.

    Required values: variant always; on every row whose error cell is empty, scorer, track_index
    and raw_score must hold a value (an empty cell fails). Only error rows may leave them empty.
    """

    variant: str
    scorer: str
    gene_id: str | None = None
    gene_name: str | None = None
    gene_strand: str | None = None
    junction_start: int | None = None
    junction_end: int | None = None
    track_index: int | None = None
    track_name: str | None = None
    track_strand: str | None = None
    ontology_curie: str | None = None
    biosample_name: str | None = None
    biosample_type: str | None = None
    assay_title: str | None = None
    gtex_tissue: str | None = None
    transcription_factor: str | None = None
    histone_mark: str | None = None
    data_source: str | None = None
    raw_score: float | None = None
    quantile_score: float | None = None
    error: str | None = None
    SKILL: ClassVar[str] = SKILL_KDENSE
    CONVERT: ClassVar = {"variant": _text, "scorer": _str, "gene_id": _ostr, "gene_name": _ostr,
                         "gene_strand": _ostr, "junction_start": _oint, "junction_end": _oint,
                         "track_index": _oint, "track_name": _ostr, "track_strand": _ostr,
                         "ontology_curie": _ostr, "biosample_name": _ostr,
                         "biosample_type": _ostr, "assay_title": _ostr, "gtex_tissue": _ostr,
                         "transcription_factor": _ostr, "histone_mark": _ostr,
                         "data_source": _ostr, "raw_score": _ofloat, "quantile_score": _ofloat,
                         "error": _ostr}
    OPTIONAL: ClassVar = frozenset(set(CONVERT) - {"variant", "scorer"})

    @classmethod
    def _required(cls, d: Mapping[str, Any]) -> set[str]:
        return {"variant", "scorer"} | (set() if "error" in d else {"track_index", "raw_score"})

    @classmethod
    def from_dict(cls, d: Mapping[str, Any], *, source: str) -> TrackScoreRow:
        record = super().from_dict(d, source=source)
        if record.error is None:
            for name in ("scorer", "track_index", "raw_score"):
                if getattr(record, name) in (None, ""):
                    raise DependencyOutputError(source, name, cls.SKILL, "missing value on a row without error")
        return record


@_record
class ModelScoreRow(_Record):
    """One row of `score_variants.py` (the official `variant_scorers.tidy_scores` table).

    Attributes:
      variant_id: chr:pos:REF>ALT, 1-based (indels allowed).
      scored_interval: Model input interval 'chr:start-end:strand'.
      gene_id: Gene ID (gene-centric scorers).
      gene_name: Gene symbol (gene-centric scorers).
      gene_type: Gene biotype (gene-centric scorers).
      gene_strand: Gene strand (gene-centric scorers).
      junction_Start: Junction start (junction scorer).
      junction_End: Junction end (junction scorer).
      output_type: Model output type, e.g. SPLICE_JUNCTIONS.
      variant_scorer: Scorer configuration string.
      track_name: Track name.
      track_strand: Track strand.
      assay_title: Column 'Assay title' (renamed: not a Python identifier).
      ontology_curie: Track ontology CURIE.
      biosample_name: Biosample name.
      biosample_type: Biosample type.
      biosample_life_stage: Life stage.
      gtex_tissue: GTEx tissue.
      data_source: Data source.
      raw_score: Raw score.
      quantile_score: Calibrated quantile (human recommended scorers).

    Optional keys: everything except variant_id, scored_interval, output_type, variant_scorer and
    raw_score; gene, junction and track-metadata columns depend on the scorers requested.
    """

    variant_id: str
    scored_interval: str
    gene_id: str | None = None
    gene_name: str | None = None
    gene_type: str | None = None
    gene_strand: str | None = None
    junction_Start: int | None = None
    junction_End: int | None = None
    output_type: str
    variant_scorer: str
    track_name: str | None = None
    track_strand: str | None = None
    assay_title: str | None = None
    ontology_curie: str | None = None
    biosample_name: str | None = None
    biosample_type: str | None = None
    biosample_life_stage: str | None = None
    gtex_tissue: str | None = None
    data_source: str | None = None
    raw_score: float
    quantile_score: float | None = None
    SKILL: ClassVar[str] = SKILL_KDENSE
    CONVERT: ClassVar = {"variant_id": _str, "scored_interval": _str, "gene_id": _ostr,
                         "gene_name": _ostr, "gene_type": _ostr, "gene_strand": _ostr,
                         "junction_Start": _oint, "junction_End": _oint, "output_type": _str,
                         "variant_scorer": _str, "track_name": _ostr, "track_strand": _ostr,
                         "assay_title": _ostr, "ontology_curie": _ostr, "biosample_name": _ostr,
                         "biosample_type": _ostr, "biosample_life_stage": _ostr,
                         "gtex_tissue": _ostr, "data_source": _ostr, "raw_score": _float,
                         "quantile_score": _ofloat}
    ALIASES: ClassVar = {"Assay title": "assay_title"}
    OPTIONAL: ClassVar = frozenset(
        {"gene_id", "gene_name", "gene_type", "gene_strand", "junction_Start", "junction_End",
         "track_name", "track_strand", "Assay title", "ontology_curie", "biosample_name",
         "biosample_type", "biosample_life_stage", "gtex_tissue", "data_source", "quantile_score"})


# --------------------------------------------------------------------------------------------
# embl-ebi-ols: search_ols.py (ols_search_*.json), get_term.py (ols_term_*.json)
# --------------------------------------------------------------------------------------------

def _check_status(d: Any, source: str, skill: str) -> None:
    if isinstance(d, Mapping) and d.get("status") not in (None, "success"):
        raise DependencyOutputError(source, "status", skill, f"{d.get('status')}: {d.get('message', '')}")


@_record
class OlsPagination(_Record):
    """Pagination block of search_ols.py.

    Attributes:
      start: Offset requested.
      rows: Rows requested.
      has_more: True when more hits exist after this page.
    """

    start: int
    rows: int
    has_more: bool
    SKILL: ClassVar[str] = SKILL_OLS
    CONVERT: ClassVar = {"start": _int, "rows": _int, "has_more": _bool}


@_record
class OlsSearchHit(_Record):
    """One term of search_ols.py output.

    Attributes:
      iri: Term IRI.
      label: Term label.
      description: Definitions.
      ontology_name: Ontology id (e.g. uberon).
      ontology_prefix: Ontology prefix (e.g. UBERON).
      obo_id: CURIE.
      short_form: CURIE with '_'.
      type: Entity type (class, property, ...).
      is_defining_ontology: As served (often false even with --defining).
      exact_synonyms: Exact synonyms.
    """

    iri: str
    label: str
    description: tuple[str, ...]
    ontology_name: str
    ontology_prefix: str
    obo_id: str
    short_form: str
    type: str
    is_defining_ontology: bool
    exact_synonyms: tuple[str, ...]
    SKILL: ClassVar[str] = SKILL_OLS
    CONVERT: ClassVar = {"iri": _str, "label": _str, "description": _str_tuple,
                         "ontology_name": _str, "ontology_prefix": _str, "obo_id": _str,
                         "short_form": _str, "type": _str, "is_defining_ontology": _bool,
                         "exact_synonyms": _str_tuple}


@_record
class OlsSearchResult(_Record):
    """The whole search_ols.py output (one per ols_search_*.json).

    Attributes:
      status: 'success' (an error file raises DependencyOutputError on status).
      total_found: Hits in OLS.
      results_count: Hits in this page.
      pagination: Page block.
      terms: Hits of this page.
    """

    status: str
    total_found: int
    results_count: int
    pagination: OlsPagination
    terms: tuple[OlsSearchHit, ...]
    SKILL: ClassVar[str] = SKILL_OLS
    CONVERT: ClassVar = {"status": _str, "total_found": _int, "results_count": _int,
                         "pagination": _one(OlsPagination), "terms": _many(OlsSearchHit)}

    @classmethod
    def from_dict(cls, d: Mapping[str, Any], *, source: str) -> OlsSearchResult:
        _check_status(d, source, cls.SKILL)
        return super().from_dict(d, source=source)


@_record
class OlsRelatedTerm(_Record):
    """One related term listed under a relation of get_term.py.

    Attributes:
      iri: Term IRI.
      label: Term label.
      obo_id: CURIE.
    """

    iri: str
    label: str
    obo_id: str
    SKILL: ClassVar[str] = SKILL_OLS
    CONVERT: ClassVar = {"iri": _str, "label": _str, "obo_id": _str}


_OLS_RELATIONS = ("parents", "children", "ancestors", "descendants", "hierarchicalParents",
                  "hierarchicalChildren", "hierarchicalAncestors", "hierarchicalDescendants")


@_record
class OlsTerm(_Record):
    """The `term` object of get_term.py.

    Attributes:
      iri: Term IRI.
      label: Term label.
      description: Definitions.
      obo_id: CURIE.
      ontology_name: Ontology id.
      ontology_prefix: Ontology prefix.
      is_defining_ontology: As served.
      is_obsolete: Obsolete flag.
      has_children: Has is-a children.
      is_root: Root flag.
      short_form: CURIE with '_'.
      synonyms: Synonyms.
      annotation: Free-form OLS annotation object.
      in_subset: Subsets.
      parents: Relation list (only when requested with --relations).
      children: Relation list (idem).
      ancestors: Relation list (idem).
      descendants: Relation list (idem).
      hierarchicalParents: Relation list (idem).
      hierarchicalChildren: Relation list (idem).
      hierarchicalAncestors: Relation list (idem).
      hierarchicalDescendants: Relation list (idem).
      graph: Graph JSON (only with --relations graph).

    Optional keys: the relation lists and graph, written only for the relations requested.
    """

    iri: str
    label: str
    description: tuple[str, ...]
    obo_id: str
    ontology_name: str
    ontology_prefix: str
    is_defining_ontology: bool
    is_obsolete: bool
    has_children: bool
    is_root: bool
    short_form: str
    synonyms: tuple[str, ...]
    annotation: Mapping[str, Any]
    in_subset: tuple[str, ...]
    parents: tuple[OlsRelatedTerm, ...] | None = None
    children: tuple[OlsRelatedTerm, ...] | None = None
    ancestors: tuple[OlsRelatedTerm, ...] | None = None
    descendants: tuple[OlsRelatedTerm, ...] | None = None
    hierarchicalParents: tuple[OlsRelatedTerm, ...] | None = None
    hierarchicalChildren: tuple[OlsRelatedTerm, ...] | None = None
    hierarchicalAncestors: tuple[OlsRelatedTerm, ...] | None = None
    hierarchicalDescendants: tuple[OlsRelatedTerm, ...] | None = None
    graph: Mapping[str, Any] | None = None
    SKILL: ClassVar[str] = SKILL_OLS
    CONVERT: ClassVar = {"iri": _str, "label": _str, "description": _str_tuple, "obo_id": _str,
                         "ontology_name": _str, "ontology_prefix": _str,
                         "is_defining_ontology": _bool, "is_obsolete": _bool,
                         "has_children": _bool, "is_root": _bool, "short_form": _str,
                         "synonyms": _str_tuple, "annotation": _mapping,
                         "in_subset": lambda v, source="": () if v is None else _str_tuple(v),  # OLS sends null for none
                         **{rel: _many(OlsRelatedTerm) for rel in _OLS_RELATIONS},
                         "graph": _omapping}
    OPTIONAL: ClassVar = frozenset({*_OLS_RELATIONS, "graph"})


@_record
class OlsTermResult(_Record):
    """The whole get_term.py output (one per ols_term_*.json).

    Attributes:
      status: 'success'.
      term: The term.
    """

    status: str
    term: OlsTerm
    SKILL: ClassVar[str] = SKILL_OLS
    CONVERT: ClassVar = {"status": _str, "term": _one(OlsTerm)}

    @classmethod
    def from_dict(cls, d: Mapping[str, Any], *, source: str) -> OlsTermResult:
        _check_status(d, source, cls.SKILL)
        return super().from_dict(d, source=source)


# --------------------------------------------------------------------------------------------
# ensembl-database get-sequence (sequence_{chr}_{start}_{end}.txt) and regulatory_features.py
# --------------------------------------------------------------------------------------------

@_record
class EnsemblSequence(_Record):
    """Genomic sequence of `ensembl_api.py get-sequence REGION --output` (plain text, no header).

    The dependency writes only the bases, so the region comes from the file name
    sequence_{chromosome}_{start}_{end}.txt (1-based closed, as requested).

    Attributes:
      chromosome: Contig with chr prefix (from the file name).
      start: First base, 1-based (from the file name).
      end: Last base, 1-based (from the file name).
      sequence: Bases A/C/G/T/N; length must be end - start + 1.
    """

    chromosome: str
    start: int
    end: int
    sequence: str
    SKILL: ClassVar[str] = SKILL_ENSEMBL
    CONVERT: ClassVar = {"chromosome": _str, "start": _int, "end": _int, "sequence": _str}

    @classmethod
    def from_dict(cls, d: Mapping[str, Any], *, source: str) -> EnsemblSequence:
        record = super().from_dict(d, source=source)
        if not re.fullmatch(r"[ACGTNacgtn]*", record.sequence):
            raise DependencyOutputError(source, "sequence", cls.SKILL, "non-ACGTN characters")
        if len(record.sequence) != record.end - record.start + 1:
            raise DependencyOutputError(source, "sequence", cls.SKILL,
                                        f"length {len(record.sequence)} != {record.end - record.start + 1}")
        return record


@_record
class RegulatoryFeature(_Record):
    """One Ensembl Regulatory Build feature of `regulatory_features.py overlap` (REST, as served).

    Attributes:
      id: Stable ID (ENSR...).
      feature_type: Always 'regulatory'.
      start: First base, 1-based.
      end: Last base, 1-based.
      strand: 0, 1 or -1.
      description: Biotype (promoter, enhancer, CTCF_binding_site, open_chromatin_region, ...).
      extended_start: Extended start, 1-based.
      extended_end: Extended end, 1-based.
      seq_region_name: Chromosome without chr prefix.
      source: 'Ensembl'.
    """

    id: str
    feature_type: str
    start: int
    end: int
    strand: int
    description: str
    extended_start: int
    extended_end: int
    seq_region_name: str
    source: str
    SKILL: ClassVar[str] = SKILL_REGULATORY
    CONVERT: ClassVar = {"id": _str, "feature_type": _str, "start": _int, "end": _int,
                         "strand": _int, "description": _str, "extended_start": _int,
                         "extended_end": _int, "seq_region_name": _str, "source": _str}


# --------------------------------------------------------------------------------------------
# clinvar-database: clinvar_api.py search / summary / evidence
# --------------------------------------------------------------------------------------------

@_record
class ClinVarSearchResult(_Record):
    """Output of `clinvar_api.py search` (clinvar_search.json).

    Attributes:
      total_count: Matches in ClinVar.
      fetched_count: IDs retrieved (must equal total_count for a complete pull).
      variant_ids: ClinVar Variation IDs.
    """

    total_count: int
    fetched_count: int
    variant_ids: tuple[str, ...]
    SKILL: ClassVar[str] = SKILL_CLINVAR
    CONVERT: ClassVar = {"total_count": _int, "fetched_count": _int, "variant_ids": _str_tuple}


@_record
class ClinVarGene(_Record):
    """One gene of a ClinVar summary.

    Attributes:
      gene_id: NCBI Gene ID.
      symbol: Gene symbol.
      strand: '+', '-' or ''.
    """

    gene_id: str
    symbol: str
    strand: str
    SKILL: ClassVar[str] = SKILL_CLINVAR
    CONVERT: ClassVar = {"gene_id": _str, "symbol": _str, "strand": _str}


@_record
class ClinVarSummary(_Record):
    """One record of `clinvar_api.py summary` (clinvar_summary[_NNN].json is a list of these).

    Attributes:
      variant_id: ClinVar Variation ID.
      title: Preferred title (HGVS for small variants).
      clinical_significance: Classification ('Unknown' when absent).
      review_status: Review status ('Unknown' when absent).
      last_evaluated: 'YYYY/MM/DD HH:MM' or 'Unknown'.
      phenotypes: Trait names.
      genes: Genes the record spans (CNVs list hundreds).
      variation_type: e.g. 'single nucleotide variant', 'Deletion', 'copy number gain'.
      molecular_consequences: e.g. 'intron variant'.
    """

    variant_id: str
    title: str
    clinical_significance: str
    review_status: str
    last_evaluated: str
    phenotypes: tuple[str, ...]
    genes: tuple[ClinVarGene, ...]
    variation_type: str
    molecular_consequences: tuple[str, ...]
    SKILL: ClassVar[str] = SKILL_CLINVAR
    CONVERT: ClassVar = {"variant_id": _str, "title": _str, "clinical_significance": _str,
                         "review_status": _str, "last_evaluated": _str, "phenotypes": _str_tuple,
                         "genes": _many(ClinVarGene), "variation_type": _str,
                         "molecular_consequences": _str_tuple}


@_record
class ClinVarAllele(_Record):
    """`allele_info` of `clinvar_api.py evidence` (first GRCh38 SequenceLocation of the record).

    Attributes:
      chromosome: Chromosome without chr prefix.
      position_start: Start (written as text; converted to int).
      position_stop: Stop (idem).
      reference_allele: Reference allele or null.
      alternate_allele: Alternate allele or null.
      cytogenetic_band: Cytogenetic band.
      dbsnp_rsid: 'rs…' or null.
    """

    chromosome: str | None
    position_start: int | None
    position_stop: int | None
    reference_allele: str | None
    alternate_allele: str | None
    cytogenetic_band: str | None
    dbsnp_rsid: str | None
    SKILL: ClassVar[str] = SKILL_CLINVAR
    CONVERT: ClassVar = {"chromosome": _ostr, "position_start": _oint, "position_stop": _oint,
                         "reference_allele": _ostr, "alternate_allele": _ostr,
                         "cytogenetic_band": _ostr, "dbsnp_rsid": _ostr}


@_record
class ClinVarCondition(_Record):
    """One condition of a ClinVar evidence record.

    Attributes:
      name: Preferred name ('Unknown' when absent).
      medgen_cui: MedGen CUI.
      omim_id: OMIM ID.
      orphanet_id: Orphanet ID.
      hpo_terms: HPO IDs.
    """

    name: str
    medgen_cui: str | None
    omim_id: str | None
    orphanet_id: str | None
    hpo_terms: tuple[str, ...]
    SKILL: ClassVar[str] = SKILL_CLINVAR
    CONVERT: ClassVar = {"name": _str, "medgen_cui": _ostr, "omim_id": _ostr,
                         "orphanet_id": _ostr, "hpo_terms": _str_tuple}


@_record
class ClinVarConsequence(_Record):
    """One functional consequence (Sequence Ontology) of a ClinVar evidence record.

    Attributes:
      value: Consequence label.
      sequence_ontology_id: SO ID or null.
    """

    value: str
    sequence_ontology_id: str | None
    SKILL: ClassVar[str] = SKILL_CLINVAR
    CONVERT: ClassVar = {"value": _str, "sequence_ontology_id": _ostr}


@_record
class ClinVarStructuralVariant(_Record):
    """`structural_variant_details` of a ClinVar evidence record (CNVs only, else null).

    Attributes:
      outer_start: Outer start or null.
      inner_start: Inner start or null.
      inner_stop: Inner stop or null.
      outer_stop: Outer stop or null.
      copy_number: Absolute copy number or null.
    """

    outer_start: int | None
    inner_start: int | None
    inner_stop: int | None
    outer_stop: int | None
    copy_number: str | None
    SKILL: ClassVar[str] = SKILL_CLINVAR
    CONVERT: ClassVar = {"outer_start": _oint, "inner_start": _oint, "inner_stop": _oint,
                         "outer_stop": _oint, "copy_number": _ostr}


@_record
class ClinVarSubmission(_Record):
    """One submitter record of a ClinVar evidence record.

    Attributes:
      submitter_name: Submitter.
      classification: Submitted classification.
      curator_notes: Comments joined with '; ' or null.
      assertion_criteria: Assertion method or review status, or null.
      date_last_evaluated: 'YYYY-MM-DD' or null.
    """

    submitter_name: str
    classification: str
    curator_notes: str | None
    assertion_criteria: str | None
    date_last_evaluated: str | None
    SKILL: ClassVar[str] = SKILL_CLINVAR
    CONVERT: ClassVar = {"submitter_name": _str, "classification": _str, "curator_notes": _ostr,
                         "assertion_criteria": _ostr, "date_last_evaluated": _ostr}


@_record
class ClinVarEvidence(_Record):
    """Output of `clinvar_api.py evidence --variant_id ID` (clinvar_evidence_{ID}.json).

    Attributes:
      variant_id: ClinVar Variation ID (must match the file name).
      allele_info: Location and alleles.
      conditions: Conditions.
      functional_consequences: SO consequences.
      structural_variant_details: CNV details or null.
      citation_references: PubMed IDs.
      submissions: Submitter records.
    """

    variant_id: str
    allele_info: ClinVarAllele
    conditions: tuple[ClinVarCondition, ...]
    functional_consequences: tuple[ClinVarConsequence, ...]
    structural_variant_details: ClinVarStructuralVariant | None
    citation_references: tuple[str, ...]
    submissions: tuple[ClinVarSubmission, ...]
    SKILL: ClassVar[str] = SKILL_CLINVAR
    CONVERT: ClassVar = {"variant_id": _str, "allele_info": _one(ClinVarAllele),
                         "conditions": _many(ClinVarCondition),
                         "functional_consequences": _many(ClinVarConsequence),
                         "structural_variant_details": _one(ClinVarStructuralVariant, nullable=True),
                         "citation_references": _str_tuple,
                         "submissions": _many(ClinVarSubmission)}


# --------------------------------------------------------------------------------------------
# gnomad-database: search_variants.py --chrom --start --end (gnomad_region.json, or
# gnomad_region_{i}.json for the tiles of a region split into queries of at most 15 kb)
# --------------------------------------------------------------------------------------------

@_record
class GnomadPopulation(_Record):
    """One genetic ancestry group of a gnomAD frequency block.

    Attributes:
      id: Group id (afr, nfe, XX, ...).
      ac: Allele count.
      an: Allele number.
      homozygote_count: Homozygotes.
      hemizygote_count: Hemizygotes (chrX/Y; null elsewhere).
    """

    id: str
    ac: int
    an: int
    homozygote_count: int | None
    hemizygote_count: int | None
    SKILL: ClassVar[str] = SKILL_GNOMAD
    CONVERT: ClassVar = {"id": _str, "ac": _int, "an": _int, "homozygote_count": _oint,
                         "hemizygote_count": _oint}


@_record
class GnomadFrequency(_Record):
    """The exome or genome block of a gnomAD variant.

    Attributes:
      ac: Allele count.
      an: Allele number.
      af: Allele frequency or null.
      homozygote_count: Homozygotes.
      hemizygote_count: Hemizygotes.
      populations: Per-group counts.
    """

    ac: int
    an: int
    af: float | None
    homozygote_count: int | None
    hemizygote_count: int | None
    populations: tuple[GnomadPopulation, ...]
    SKILL: ClassVar[str] = SKILL_GNOMAD
    CONVERT: ClassVar = {"ac": _int, "an": _int, "af": _ofloat, "homozygote_count": _oint,
                         "hemizygote_count": _oint, "populations": _many(GnomadPopulation)}


@_record
class GnomadRegionVariant(_Record):
    """One variant of data.region.variants in gnomad_region.json (the GraphQL answer as served).

    The query selects no filters (PASS/AC0/...) field, so the file cannot tell PASS variants apart.

    Attributes:
      variant_id: 'chrom-pos-ref-alt' (VCF-style, no chr prefix; see gnomad_to_variant).
      rsids: dbSNP IDs or null.
      consequence: Most severe consequence or null.
      exome: Exome frequencies or null.
      genome: Genome frequencies or null.
    """

    variant_id: str
    rsids: tuple[str, ...] | None
    consequence: str | None
    exome: GnomadFrequency | None
    genome: GnomadFrequency | None
    SKILL: ClassVar[str] = SKILL_GNOMAD
    CONVERT: ClassVar = {"variant_id": _str, "rsids": _ostr_tuple, "consequence": _ostr,
                         "exome": _one(GnomadFrequency, nullable=True),
                         "genome": _one(GnomadFrequency, nullable=True)}


# --------------------------------------------------------------------------------------------
# dbsnp-database: dbsnp_cli.py resolve-rsid (dbsnp_rsid_{N}.json)
# --------------------------------------------------------------------------------------------

@_record
class DbsnpSpdi(_Record):
    """SPDI of one dbSNP allele (NCBI Variation Services, as served).

    Attributes:
      seq_id: RefSeq sequence accession (NC_0000NN.v; see refseq_chromosome).
      position: 0-based inter-residue position of the first deleted base.
      deleted_sequence: Reference bases replaced; NCBI expands indels to the whole ambiguous
        (repeat) context, so a 1-bp deletion in a run can read e.g. 'AAA' -> 'AA'.
      inserted_sequence: Bases inserted in their place (equal to deleted_sequence for the
        reference allele).
    """

    seq_id: str
    position: int
    deleted_sequence: str
    inserted_sequence: str
    SKILL: ClassVar[str] = SKILL_DBSNP
    CONVERT: ClassVar = {"seq_id": _str, "position": _int, "deleted_sequence": _str,
                         "inserted_sequence": _str}


@_record
class DbsnpAlleleSpec(_Record):
    """The `allele` object of a dbSNP placement allele.

    Attributes:
      spdi: The allele as SPDI.
    """

    spdi: DbsnpSpdi
    SKILL: ClassVar[str] = SKILL_DBSNP
    CONVERT: ClassVar = {"spdi": _one(DbsnpSpdi)}


@_record
class DbsnpAllele(_Record):
    """One allele of a dbSNP placement.

    Attributes:
      allele: SPDI wrapper.
      hgvs: Genomic HGVS ('NC_…:g.N=' for the reference allele, 3'-shifted for indels).
    """

    allele: DbsnpAlleleSpec
    hgvs: str
    SKILL: ClassVar[str] = SKILL_DBSNP
    CONVERT: ClassVar = {"allele": _one(DbsnpAlleleSpec), "hgvs": _str}

    @property
    def is_reference(self) -> bool:
        """True for the reference allele (HGVS ending in '=')."""
        return self.hgvs.endswith("=")


@_record
class DbsnpPlacement(_Record):
    """One top-level placement of an rsID on the requested assembly.

    Attributes:
      seq_id: RefSeq sequence accession.
      alleles: Reference and alternate alleles.
    """

    seq_id: str
    alleles: tuple[DbsnpAllele, ...]
    SKILL: ClassVar[str] = SKILL_DBSNP
    CONVERT: ClassVar = {"seq_id": _str, "alleles": _many(DbsnpAllele)}


@_record
class DbsnpRsidPlacements(_Record):
    """Output of `dbsnp_cli.py resolve-rsid rsN --output dbsnp_rsid_{N}.json`.

    The script writes no file when the rsID has no current snapshot (merged or withdrawn); it
    prints a JSON error and exits 1 instead.

    Attributes:
      rsid: Numeric rsID without 'rs' (must match the file name).
      assembly: RefSeq assembly accession (GCF_000001405.40 = GRCh38).
      placements: Top-level placements on that assembly (usually one).
    """

    rsid: str
    assembly: str
    placements: tuple[DbsnpPlacement, ...]
    SKILL: ClassVar[str] = SKILL_DBSNP
    CONVERT: ClassVar = {"rsid": _str, "assembly": _str, "placements": _many(DbsnpPlacement)}


# --------------------------------------------------------------------------------------------
# ensembl-database: ensembl_api.py vep HGVS (vep_clinvar_{id}.json)
# --------------------------------------------------------------------------------------------

@_record
class VepResult(_Record):
    """One result of `ensembl_api.py vep HGVS --output vep_clinvar_{id}.json` (Ensembl VEP REST).

    The file is the JSON list VEP returns; the script writes {} when VEP returns nothing, read
    as no result. Only the location and alleles are modelled field by field; the consequence
    blocks (with the plugin fields the script requests) are kept as read-only mappings.

    Attributes:
      input: The HGVS string sent.
      id: Identifier VEP echoes (the input for HGVS queries).
      assembly_name: Assembly (GRCh38).
      seq_region_name: Chromosome without chr prefix.
      start: First reference base, 1-based (for an insertion, start = end + 1).
      end: Last reference base, 1-based.
      strand: Strand of allele_string (1 or -1).
      allele_string: 'REF/ALT' with '-' for an empty allele.
      most_severe_consequence: Most severe SO term over all transcripts.
      transcript_consequences: Per-transcript blocks (mappings), when served.
      colocated_variants: Known variants at the position (mappings), when served.
      regulatory_feature_consequences: Regulatory blocks (mappings), when served.
      motif_feature_consequences: Motif blocks (mappings), when served.
      intergenic_consequences: Intergenic blocks (mappings), when served.
    """

    input: str
    assembly_name: str
    seq_region_name: str
    start: int
    end: int
    strand: int
    allele_string: str
    most_severe_consequence: str
    id: str | None = None
    transcript_consequences: tuple[Mapping[str, Any], ...] | None = None
    colocated_variants: tuple[Mapping[str, Any], ...] | None = None
    regulatory_feature_consequences: tuple[Mapping[str, Any], ...] | None = None
    motif_feature_consequences: tuple[Mapping[str, Any], ...] | None = None
    intergenic_consequences: tuple[Mapping[str, Any], ...] | None = None
    SKILL: ClassVar[str] = SKILL_ENSEMBL
    OPTIONAL: ClassVar = frozenset({"id", "transcript_consequences", "colocated_variants",
                                    "regulatory_feature_consequences", "motif_feature_consequences",
                                    "intergenic_consequences"})
    CONVERT: ClassVar = {"input": _str, "id": _ostr, "assembly_name": _str, "seq_region_name": _str,
                         "start": _int, "end": _int, "strand": _int, "allele_string": _str,
                         "most_severe_consequence": _str, "transcript_consequences": _mappings,
                         "colocated_variants": _mappings, "regulatory_feature_consequences": _mappings,
                         "motif_feature_consequences": _mappings, "intergenic_consequences": _mappings}


# --------------------------------------------------------------------------------------------
# alphagenome-single-variant-analysis: interpret_splicing.py (stdout, interpret_*.txt)
# --------------------------------------------------------------------------------------------

@_record
class SpliceSiteChange(_Record):
    """A 'Position P (Offset k): REF=x, ALT=y, Delta=z' line (|Delta| > 0.05).

    Attributes:
      position: Genomic position, 1-based.
      offset: Offset from the variant (-2..2).
      REF: Splice-site usage, REF allele.
      ALT: Splice-site usage, ALT allele.
      Delta: ALT - REF.
      flag: Text after '->' on the next line ('Loss of strong splice site at P!' or 'Gain …').
    """

    position: int
    offset: int
    REF: float
    ALT: float
    Delta: float
    flag: str | None = None
    SKILL: ClassVar[str] = SKILL_SVA
    CONVERT: ClassVar = {"position": _int, "offset": _int, "REF": _float, "ALT": _float,
                         "Delta": _float, "flag": _ostr}
    OPTIONAL: ClassVar = frozenset({"flag"})


@_record
class JunctionDelta(_Record):
    """A junction overlapping the variant: '  S-E (Strand s): REF=x, ALT=y, Delta=z'.

    Attributes:
      start: Junction start as served.
      end: Junction end as served.
      strand: '+', '-' or '.'.
      REF: Junction score, REF allele.
      ALT: Junction score, ALT allele.
      Delta: ALT - REF.
      flag: Text after '->' on the next line, e.g. 'Skipping/Loss of canonical intron (Exon skipping?)'.
    """

    start: int
    end: int
    strand: str
    REF: float
    ALT: float
    Delta: float
    flag: str | None = None
    SKILL: ClassVar[str] = SKILL_SVA
    CONVERT: ClassVar = {"start": _int, "end": _int, "strand": _str, "REF": _float, "ALT": _float,
                         "Delta": _float, "flag": _ostr}
    OPTIONAL: ClassVar = frozenset({"flag"})


@_record
class CrypticJunction(_Record):
    """A 'Potential Cryptic Junctions (Gain in ALT)' line: '  S-E: REF=x -> ALT=y'.

    Attributes:
      start: Junction start as served.
      end: Junction end as served.
      REF: Junction score, REF allele.
      ALT: Junction score, ALT allele.
    """

    start: int
    end: int
    REF: float
    ALT: float
    SKILL: ClassVar[str] = SKILL_SVA
    CONVERT: ClassVar = {"start": _int, "end": _int, "REF": _float, "ALT": _float}


@_record
class TopJunctionChange(_Record):
    """A row of the 'Top 5 Most Changed Junctions' pandas table.

    Attributes:
      index: Row label printed by pandas.
      start: Junction start as served.
      end: Junction end as served.
      score_REF: Junction score, REF allele.
      score_ALT: Junction score, ALT allele.
      delta: score_ALT - score_REF.
    """

    index: int
    start: int
    end: int
    score_REF: float
    score_ALT: float
    delta: float
    SKILL: ClassVar[str] = SKILL_SVA
    CONVERT: ClassVar = {"index": _int, "start": _int, "end": _int, "score_REF": _float,
                         "score_ALT": _float, "delta": _float}


@_record
class SpliceInterpretation(_Record):
    """Parsed stdout of interpret_splicing.py (the script has no --output; stdout is saved).

    Attributes:
      variant: chr:pos:REF>ALT from the 'Quantifying splicing changes for … in CURIE...' line.
      ontology_id: CURIE from the same line.
      splice_site_usage: '--- Splice Site Usage Analysis ---' lines.
      junctions_overlapping_variant: 'Junctions overlapping variant (P):' lines.
      potential_cryptic_junctions: 'Potential Cryptic Junctions (Gain in ALT):' lines.
      top_changed_junctions: 'Top 5 Most Changed Junctions:' table rows.
      mean_rna_coverage_ref: 'Mean RNA Coverage (+/- 50bp): REF=…'.
      mean_rna_coverage_alt: '… ALT=…'.
      log2_fold_change: 'Log2 Fold Change at variant: …' (absent when REF coverage is 0).
      messages: Warning/Info/'missing' lines printed by the script.
    """

    variant: str
    ontology_id: str
    splice_site_usage: tuple[SpliceSiteChange, ...]
    junctions_overlapping_variant: tuple[JunctionDelta, ...]
    potential_cryptic_junctions: tuple[CrypticJunction, ...]
    top_changed_junctions: tuple[TopJunctionChange, ...]
    mean_rna_coverage_ref: float | None
    mean_rna_coverage_alt: float | None
    log2_fold_change: float | None
    messages: tuple[str, ...]
    SKILL: ClassVar[str] = SKILL_SVA
    CONVERT: ClassVar = {"variant": _str, "ontology_id": _str,
                         "splice_site_usage": _many(SpliceSiteChange),
                         "junctions_overlapping_variant": _many(JunctionDelta),
                         "potential_cryptic_junctions": _many(CrypticJunction),
                         "top_changed_junctions": _many(TopJunctionChange),
                         "mean_rna_coverage_ref": _ofloat, "mean_rna_coverage_alt": _ofloat,
                         "log2_fold_change": _ofloat, "messages": _str_tuple}


_NUM = r"-?(?:\d+\.?\d*(?:[eE][-+]?\d+)?|nan|inf)"
_INTERPRET_LINES = [
    ("header", re.compile(r"^Quantifying splicing changes for (\S+) in (\S+?)\.\.\.$")),
    ("section", re.compile(r"^--- (Splice Site Usage Analysis|Junction Analysis|RNA-seq Analysis) ---$")),
    ("site", re.compile(rf"^Position (\d+) \(Offset (-?\d+)\): REF=({_NUM}), ALT=({_NUM}), Delta=({_NUM})$")),
    ("overlap_head", re.compile(r"^Junctions overlapping variant \((\d+)\):$")),
    ("none", re.compile(r"^  None found\.$")),
    ("junction", re.compile(rf"^  (\d+)-(\d+) \(Strand ([-+.])\): REF=({_NUM}), ALT=({_NUM}), Delta=({_NUM})$")),
    ("flag", re.compile(r"^  -> (.+)$")),
    ("cryptic_head", re.compile(r"^Potential Cryptic Junctions \(Gain in ALT\):$")),
    ("cryptic", re.compile(rf"^  (\d+)-(\d+): REF=({_NUM}) -> ALT=({_NUM})$")),
    ("top_head", re.compile(r"^Top 5 Most Changed Junctions:$")),
    ("top_columns", re.compile(r"^\s+start\s+end\s+score_REF\s+score_ALT\s+delta$")),
    ("top_row", re.compile(rf"^(\d+)\s+({_NUM})\s+({_NUM})\s+({_NUM})\s+({_NUM})\s+({_NUM})$")),
    ("rna", re.compile(rf"^Mean RNA Coverage \(\+/- 50bp\): REF=({_NUM}), ALT=({_NUM})$")),
    ("log2fc", re.compile(rf"^Log2 Fold Change at variant: ({_NUM})$")),
    ("message", re.compile(r"^(Warning: .*|Info: .*|Splice site usage track missing\.|"
                           r"No junctions found or failed to parse\.|RNA-seq track missing\.)$")),
]


def parse_interpret_text(text: str, source: str) -> dict[str, Any]:
    """Turns interpret_splicing.py stdout into the mapping SpliceInterpretation.from_dict takes."""
    out: dict[str, Any] = {"splice_site_usage": [], "junctions_overlapping_variant": [],
                           "potential_cryptic_junctions": [], "top_changed_junctions": [],
                           "mean_rna_coverage_ref": None, "mean_rna_coverage_alt": None,
                           "log2_fold_change": None, "messages": []}
    last: dict[str, Any] | None = None
    for n, line in enumerate(text.splitlines(), start=1):
        if not line.strip():
            continue
        for kind, pattern in _INTERPRET_LINES:
            m = pattern.match(line)
            if m:
                break
        else:
            raise DependencyOutputError(source, f"line {n}", SKILL_SVA, f"unrecognised line: {line[:80]!r}")
        g = m.groups()
        if kind == "header":
            out["variant"], out["ontology_id"] = g
        elif kind == "site":
            last = dict(zip(("position", "offset", "REF", "ALT", "Delta"), g))
            out["splice_site_usage"].append(last)
        elif kind == "junction":
            last = dict(zip(("start", "end", "strand", "REF", "ALT", "Delta"), g))
            out["junctions_overlapping_variant"].append(last)
        elif kind == "flag":
            if last is None or "flag" in last:
                raise DependencyOutputError(source, f"line {n}", SKILL_SVA, "'->' line without a record")
            last["flag"] = g[0]
        elif kind == "cryptic":
            out["potential_cryptic_junctions"].append(dict(zip(("start", "end", "REF", "ALT"), g)))
        elif kind == "top_row":
            out["top_changed_junctions"].append(
                dict(zip(("index", "start", "end", "score_REF", "score_ALT", "delta"), g)))
        elif kind == "rna":
            out["mean_rna_coverage_ref"], out["mean_rna_coverage_alt"] = g
        elif kind == "log2fc":
            out["log2_fold_change"] = g[0]
        elif kind == "message":
            out["messages"].append(g[0])
        if kind not in ("site", "junction", "flag"):
            last = None
    if "variant" not in out:
        raise DependencyOutputError(source, "line 1", SKILL_SVA, "missing 'Quantifying splicing changes' header")
    return out


# --------------------------------------------------------------------------------------------
# Curated cohort and saturation snapshot schemas
# --------------------------------------------------------------------------------------------

@_record
class CuratedCohortRow(_Record):
    """One row of snapshot_cohort.csv (curated variant cohort).

    Attributes:
      record_id: Cohort record id.
      variant_key: 'chrom-pos-ref-alt'.
      chrom: Chromosome without chr prefix.
      pos: 1-based position.
      ref: Reference allele.
      alt: Alternate allele.
      variant_type: SNV, deletion, ...
      region_annotation: e.g. 'Intron N'.
      sources: Sources joined with ';'.
      clinvar_variation_id: ClinVar Variation ID (written as a float) or null.
      clinvar_title: ClinVar title or null.
      clinvar_significance: ClinVar classification or null.
      clinvar_review_status: Review status or null.
      clinvar_phenotypes: Phenotypes or null.
      molecular_consequences: ClinVar consequences or null.
      category: Automatic class.
      manual_review_required: Review flag.
      final_category: Class after review.
      manual_evidence_type: Manual curation field (empty in the snapshot).
      manual_phenotype: Idem.
      manual_citation: Idem.
      manual_curator: Idem.
      manual_curation_date: Idem.
      manual_notes: Idem.
      analysis_split: candidate, ...
      curation_timestamp: ISO timestamp.
      gnomad_dataset: gnomAD dataset or null.
      gnomad_rsids: rsIDs or null.
      gnomad_consequence: gnomAD consequence or null.
      gnomad_ac: Allele count (float) or null.
      gnomad_an: Allele number (float) or null.
      gnomad_af: Allele frequency or null.
      gnomad_hom: Homozygotes (float) or null.
      gnomad_hemi: Hemizygotes (float) or null.
      gnomad_filters: gnomAD filters (PASS, ...) or null.
      background_assumption: Control rationale or null.
      qc_status: QC status.
      qc_reasons: QC reasons.
      SPLICE_JUNCTIONS: max |delta| SPLICE_JUNCTIONS.
      SPLICE_SITES: max |delta| SPLICE_SITES.
      SPLICE_SITE_USAGE: max |delta| SPLICE_SITE_USAGE.
      alphagenome_splicing: Composite splicing score.
      alphagenome_status: Scoring status.
      dist_nearest_ss: Distance to the nearest canonical splice site (bp).
      nearest_splice_site: Nearest canonical site label.
      distance_zone: Distance bin label.
    """

    record_id: str
    variant_key: str
    chrom: str
    pos: int
    ref: str
    alt: str
    variant_type: str
    region_annotation: str
    sources: str
    clinvar_variation_id: int | None
    clinvar_title: str | None
    clinvar_significance: str | None
    clinvar_review_status: str | None
    clinvar_phenotypes: str | None
    molecular_consequences: str | None
    category: str
    manual_review_required: bool
    final_category: str
    manual_evidence_type: str | None
    manual_phenotype: str | None
    manual_citation: str | None
    manual_curator: str | None
    manual_curation_date: str | None
    manual_notes: str | None
    analysis_split: str
    curation_timestamp: str
    gnomad_dataset: str | None
    gnomad_rsids: str | None
    gnomad_consequence: str | None
    gnomad_ac: float | None
    gnomad_an: float | None
    gnomad_af: float | None
    gnomad_hom: float | None
    gnomad_hemi: float | None
    gnomad_filters: str | None
    background_assumption: str | None
    qc_status: str
    qc_reasons: str
    SPLICE_JUNCTIONS: float
    SPLICE_SITES: float
    SPLICE_SITE_USAGE: float
    alphagenome_splicing: float
    alphagenome_status: str
    dist_nearest_ss: int
    nearest_splice_site: str
    distance_zone: str
    SKILL: ClassVar[str] = SKILL_SNAPSHOT
    CONVERT: ClassVar = {
        "record_id": _str, "variant_key": _str, "chrom": _str, "pos": _int, "ref": _str,
        "alt": _str, "variant_type": _str, "region_annotation": _str, "sources": _str,
        "clinvar_variation_id": _oint, "clinvar_title": _ostr, "clinvar_significance": _ostr,
        "clinvar_review_status": _ostr, "clinvar_phenotypes": _ostr,
        "molecular_consequences": _ostr, "category": _str, "manual_review_required": _bool,
        "final_category": _str, "manual_evidence_type": _ostr, "manual_phenotype": _ostr,
        "manual_citation": _ostr, "manual_curator": _ostr, "manual_curation_date": _ostr,
        "manual_notes": _ostr, "analysis_split": _str, "curation_timestamp": _str,
        "gnomad_dataset": _ostr, "gnomad_rsids": _ostr, "gnomad_consequence": _ostr,
        "gnomad_ac": _ofloat, "gnomad_an": _ofloat, "gnomad_af": _ofloat, "gnomad_hom": _ofloat,
        "gnomad_hemi": _ofloat, "gnomad_filters": _ostr, "background_assumption": _ostr,
        "qc_status": _str, "qc_reasons": _str, "SPLICE_JUNCTIONS": _float, "SPLICE_SITES": _float,
        "SPLICE_SITE_USAGE": _float, "alphagenome_splicing": _float, "alphagenome_status": _str,
        "dist_nearest_ss": _int, "nearest_splice_site": _str, "distance_zone": _str}


@_record
class SaturationRow(_Record):
    """One SNV of snapshot_saturation.parquet (saturation scores over the reference window).

    Attributes:
      variant: chr:pos:REF>ALT, 1-based.
      variant_key: 'chrom-pos-ref-alt'.
      chrom: Chromosome with chr prefix.
      pos: 1-based position.
      ref: Reference base.
      alt: Alternate base.
      splice_site_score: Splice-site score used by the original analysis.
      SPLICE_SITES: max |delta| SPLICE_SITES.
      SPLICE_SITE_USAGE: max |delta| SPLICE_SITE_USAGE.
      SPLICE_JUNCTIONS: max |delta| SPLICE_JUNCTIONS.
      alphagenome_splicing: Composite splicing score.
      avi_score: AVI raw score.
      avi_raw: AVI raw score.
      avi_phred: AVI Phred.
      avi_quantile: AVI tail quantile.
      is_damaging_tau_1_00: alphagenome_splicing >= 1.0.
      is_pathogenic_tau_2_65: alphagenome_splicing >= 2.65.
    """

    variant: str
    variant_key: str
    chrom: str
    pos: int
    ref: str
    alt: str
    splice_site_score: float
    SPLICE_SITES: float
    SPLICE_SITE_USAGE: float
    SPLICE_JUNCTIONS: float
    alphagenome_splicing: float
    avi_score: float
    avi_raw: float
    avi_phred: float
    avi_quantile: float
    is_damaging_tau_1_00: bool
    is_pathogenic_tau_2_65: bool
    SKILL: ClassVar[str] = SKILL_SNAPSHOT
    CONVERT: ClassVar = {
        "variant": _str, "variant_key": _str, "chrom": _str, "pos": _int, "ref": _str, "alt": _str,
        "splice_site_score": _float, "SPLICE_SITES": _float, "SPLICE_SITE_USAGE": _float,
        "SPLICE_JUNCTIONS": _float, "alphagenome_splicing": _float, "avi_score": _float,
        "avi_raw": _float, "avi_phred": _float, "avi_quantile": _float,
        "is_damaging_tau_1_00": _bool, "is_pathogenic_tau_2_65": _bool}


# --------------------------------------------------------------------------------------------
# File readers and the name-pattern table
# --------------------------------------------------------------------------------------------

def _read_json(path: Path, skill: str) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise DependencyOutputError(str(path), "(file)", skill, f"not readable JSON: {exc}") from None


def _read_table(path: Path, skill: str) -> list[dict[str, Any]]:
    try:
        if path.suffix == ".parquet":
            import pyarrow.parquet as pq  # noqa: PLC0415 - only needed for parquet files
            return pq.read_table(path).to_pylist()
        delimiter = "\t" if path.suffix == ".tsv" else ","
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            rows = list(reader)
            if reader.fieldnames is None:
                raise ValueError("empty file")
            if any(None in row for row in rows):
                raise ValueError("a row has more cells than the header")
            return rows
    except (OSError, UnicodeDecodeError, ValueError) as exc:
        raise DependencyOutputError(str(path), "(file)", skill, f"not a readable table: {exc}") from None
    except Exception as exc:  # pyarrow raises its own error types
        raise DependencyOutputError(str(path), "(file)", skill, f"not a readable table: {exc}") from None


def _records(cls: type, items: Any, source: str) -> list[Any]:
    if not isinstance(items, list):
        raise DependencyOutputError(source, "(file)", cls.SKILL, f"expected a list, got {type(items).__name__}")
    out = []
    for i, item in enumerate(items):
        try:
            out.append(cls.from_dict(item, source=source))
        except DependencyOutputError as exc:
            raise DependencyOutputError(source, f"[{i}].{exc.field}", exc.skill, exc.reason) from None
    return out


def _load_json_list(cls: type) -> Callable[[Path, re.Match], list[Any]]:
    return lambda path, m: _records(cls, _read_json(path, cls.SKILL), str(path))


def _load_json_one(cls: type) -> Callable[[Path, re.Match], list[Any]]:
    return lambda path, m: [cls.from_dict(_read_json(path, cls.SKILL), source=str(path))]


def _load_table(cls: type) -> Callable[[Path, re.Match], list[Any]]:
    return lambda path, m: _records(cls, _read_table(path, cls.SKILL), str(path))


def _load_tracks(path: Path, m: re.Match) -> list[Any]:
    rows = _records(AtlasTrack, _read_table(path, SKILL_AVI), str(path))
    scorer, curie = m.group("scorer"), m.group("curie").replace("_", ":", 1)
    for i, row in enumerate(rows):
        if row.scorer != scorer:
            raise DependencyOutputError(str(path), f"[{i}].scorer", SKILL_AVI, f"{row.scorer!r} != file name {scorer!r}")
        if row.ontology_curie != curie:
            raise DependencyOutputError(str(path), f"[{i}].ontology_curie", SKILL_AVI,
                                        f"{row.ontology_curie!r} != file name {curie!r}")
    return rows


def _load_sequence(path: Path, m: re.Match) -> list[Any]:
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        raise DependencyOutputError(str(path), "(file)", SKILL_ENSEMBL, f"not readable text: {exc}") from None
    d = {"chromosome": m.group("chrom"), "start": m.group("start"), "end": m.group("end"),
         "sequence": text.strip()}
    return [EnsemblSequence.from_dict(d, source=str(path))]


def _load_evidence(path: Path, m: re.Match) -> list[Any]:
    record = ClinVarEvidence.from_dict(_read_json(path, SKILL_CLINVAR), source=str(path))
    if record.variant_id != m.group("id"):
        raise DependencyOutputError(str(path), "variant_id", SKILL_CLINVAR,
                                    f"{record.variant_id!r} != file name {m.group('id')!r}")
    return [record]


def _load_gnomad(path: Path, m: re.Match) -> list[Any]:
    source = str(path)
    d = _read_json(path, SKILL_GNOMAD)
    if not isinstance(d, dict):
        raise DependencyOutputError(source, "(file)", SKILL_GNOMAD, "expected an object")
    if "errors" in d:
        raise DependencyOutputError(source, "errors", SKILL_GNOMAD, f"GraphQL errors: {str(d['errors'])[:200]}")
    for key, expected in (("", {"data"}), ("data", {"region"}), ("data.region", {"variants"})):
        node = d
        for part in filter(None, key.split(".")):
            node = node[part]
        if not isinstance(node, dict):
            raise DependencyOutputError(source, key or "(file)", SKILL_GNOMAD, "expected an object")
        extra, missing = sorted(set(node) - expected), sorted(expected - set(node))
        if extra or missing:
            field = f"{key}.{(extra or missing)[0]}".lstrip(".")
            raise DependencyOutputError(source, field, SKILL_GNOMAD, "unknown field" if extra else "missing field")
    rows = d["data"]["region"]["variants"]
    try:
        return _records(GnomadRegionVariant, rows, source)
    except DependencyOutputError as exc:
        raise DependencyOutputError(source, f"data.region.variants{exc.field}", exc.skill, exc.reason) from None


def _load_dbsnp(path: Path, m: re.Match) -> list[Any]:
    record = DbsnpRsidPlacements.from_dict(_read_json(path, SKILL_DBSNP), source=str(path))
    if record.rsid != m.group("rsid"):
        raise DependencyOutputError(str(path), "rsid", SKILL_DBSNP, f"{record.rsid!r} != file name {m.group('rsid')!r}")
    return [record]


def _load_vep(path: Path, m: re.Match) -> list[Any]:
    """vep_clinvar_{id}.json: the VEP list; {} (what the script writes for no result) is no record."""
    d = _read_json(path, SKILL_ENSEMBL)
    if d == {}:
        return []
    return _records(VepResult, d, str(path))


def _load_interpret(path: Path, m: re.Match) -> list[Any]:
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        raise DependencyOutputError(str(path), "(file)", SKILL_SVA, f"not readable text: {exc}") from None
    return [SpliceInterpretation.from_dict(parse_interpret_text(text, str(path)), source=str(path))]


def _load_text(path: Path, m: re.Match) -> list[Any]:
    """atlas_links.md: no dataclass, the report embeds the text as is."""
    try:
        return [path.read_text(encoding="utf-8")]
    except (OSError, UnicodeDecodeError) as exc:
        raise DependencyOutputError(str(path), "(file)", SKILL_LINKS, f"not readable text: {exc}") from None


# (file-name pattern, dataclass name for messages, skill, loader). First match wins.
RAW_PATTERNS: list[tuple[re.Pattern, str, str, Callable[[Path, re.Match], list[Any]]]] = [
    (re.compile(r"^gtf\.json$"), "GtfTranscript", SKILL_AVI, _load_json_list(GtfTranscript)),
    (re.compile(r"^gtf_junctions\.json$"), "GtfTranscript", SKILL_AVI, _load_json_list(GtfTranscript)),
    (re.compile(r"^avi_region\.parquet$"), "AviRecord", SKILL_AVI, _load_table(AviRecord)),
    (re.compile(r"^avi_query\.json$"), "AviRecord", SKILL_AVI, _load_json_list(AviRecord)),
    (re.compile(r"^tracks_(?P<scorer>[A-Z0-9_]+?)_(?P<curie>[A-Za-z]+_\d+)\.parquet$"), "AtlasTrack",
     SKILL_AVI, _load_tracks),
    (re.compile(r"^tracks_[A-Z0-9_]+_catalogue\.parquet$"), "AtlasTrack", SKILL_AVI, _load_table(AtlasTrack)),
    (re.compile(r"^model_scores_indels\.tsv$"), "ModelScoreRow", SKILL_KDENSE, _load_table(ModelScoreRow)),
    (re.compile(r"^scores_[A-Za-z0-9_]+\.tsv$"), "TrackScoreRow", SKILL_KDENSE, _load_table(TrackScoreRow)),
    (re.compile(r"^ols_search_[A-Za-z0-9_.-]+\.json$"), "OlsSearchResult", SKILL_OLS,
     _load_json_one(OlsSearchResult)),
    (re.compile(r"^ols_term_[A-Za-z]+_\d+\.json$"), "OlsTermResult", SKILL_OLS, _load_json_one(OlsTermResult)),
    (re.compile(r"^sequence_(?P<chrom>chr[0-9A-Za-z]+)_(?P<start>\d+)_(?P<end>\d+)\.txt$"), "EnsemblSequence",
     SKILL_ENSEMBL, _load_sequence),
    (re.compile(r"^regulatory\.json$"), "RegulatoryFeature", SKILL_REGULATORY,
     _load_json_list(RegulatoryFeature)),
    (re.compile(r"^clinvar_search\.json$"), "ClinVarSearchResult", SKILL_CLINVAR,
     _load_json_one(ClinVarSearchResult)),
    (re.compile(r"^clinvar_summary(_\d+)?\.json$"), "ClinVarSummary", SKILL_CLINVAR,
     _load_json_list(ClinVarSummary)),
    (re.compile(r"^clinvar_evidence_(?P<id>\d+)\.json$"), "ClinVarEvidence", SKILL_CLINVAR, _load_evidence),
    (re.compile(r"^gnomad_region\.json$"), "GnomadRegionVariant", SKILL_GNOMAD, _load_gnomad),
    (re.compile(r"^gnomad_region_(?P<i>\d+)\.json$"), "GnomadRegionVariant", SKILL_GNOMAD, _load_gnomad),
    (re.compile(r"^dbsnp_rsid_(?P<rsid>\d+)\.json$"), "DbsnpRsidPlacements", SKILL_DBSNP, _load_dbsnp),
    (re.compile(r"^vep_clinvar_(?P<id>\d+)\.json$"), "VepResult", SKILL_ENSEMBL, _load_vep),
    (re.compile(r"^interpret_[A-Za-z0-9_.>:-]+\.txt$"), "SpliceInterpretation", SKILL_SVA, _load_interpret),
    (re.compile(r"^atlas_links\.md$"), "text (no dataclass)", SKILL_LINKS, _load_text),
    (re.compile(r"^snapshot_cohort\.csv$"), "CuratedCohortRow", SKILL_SNAPSHOT, _load_table(CuratedCohortRow)),
    (re.compile(r"^snapshot_saturation\.parquet$"), "SaturationRow", SKILL_SNAPSHOT,
     _load_table(SaturationRow)),
]


def _match(path: Path) -> tuple[re.Match, str, str, Callable[[Path, re.Match], list[Any]]]:
    for pattern, kind, skill, loader in RAW_PATTERNS:
        m = pattern.match(path.name)
        if m:
            return m, kind, skill, loader
    raise DependencyOutputError(str(path), "(file name)", "unknown",
                                "no dataclass for this file-name pattern (see RAW_PATTERNS)")


def load_raw(path: Path) -> list[Any]:
    """Parses one raw file with the dataclass of its name pattern; returns its top-level records.

    List files (gtf, regulatory, avi_query, clinvar_summary, tables, gnomAD variants) return one
    record per entry or row; single-object files (OLS, clinvar_search, clinvar_evidence,
    sequence, interpret) return a one-element list; atlas_links.md returns [text].
    """
    path = Path(path)
    m, _kind, _skill, loader = _match(path)
    return loader(path, m)


# --------------------------------------------------------------------------------------------
# Coordinate helpers (1-based chr:pos:REF>ALT)
# --------------------------------------------------------------------------------------------

_ALLELE_RE = re.compile(r"^[ACGTN]+$")


def variant_str(chromosome: str, position: int, ref: str, alt: str) -> str:
    """Formats a 1-based variant as chr:pos:REF>ALT (adds the chr prefix, MT -> chrM)."""
    chrom = str(chromosome)
    chrom = chrom if chrom.startswith("chr") else "chr" + ("M" if chrom == "MT" else chrom)
    ref, alt = ref.upper(), alt.upper()
    if int(position) < 1 or not _ALLELE_RE.match(ref) or not _ALLELE_RE.match(alt):
        raise ValueError(f"invalid variant {chromosome}:{position}:{ref}>{alt}")
    return f"{chrom}:{int(position)}:{ref}>{alt}"


def parse_variant(text: str) -> tuple[str, int, str, str]:
    """Splits 'chr:pos:REF>ALT' into (chromosome, position, ref, alt)."""
    m = re.fullmatch(r"(chr[0-9A-Za-z]+):(\d+):([ACGTNacgtn]+)>([ACGTNacgtn]+)", text.strip())
    if not m:
        raise ValueError(f"expected chr:pos:REF>ALT, got {text!r}")
    return m.group(1), int(m.group(2)), m.group(3).upper(), m.group(4).upper()


def interval_str(chromosome: str, start: int, end: int) -> str:
    """Formats a 1-based closed interval as chr:start-end."""
    chrom = str(chromosome)
    chrom = chrom if chrom.startswith("chr") else "chr" + ("M" if chrom == "MT" else chrom)
    if int(start) < 1 or int(end) < int(start):
        raise ValueError(f"invalid interval {chromosome}:{start}-{end}")
    return f"{chrom}:{int(start)}-{int(end)}"


def gnomad_to_variant(variant_id: str) -> str:
    """gnomAD '1-1000-AG-A' (VCF-style, 1-based, anchored indels) -> 'chr1:1000:AG>A'."""
    parts = variant_id.split("-")
    if len(parts) != 4 or not parts[1].isdigit():
        raise ValueError(f"expected chrom-pos-ref-alt, got {variant_id!r}")
    return variant_str(parts[0], int(parts[1]), parts[2], parts[3])


def clinvar_allele_to_variant(allele: ClinVarAllele) -> str | None:
    """ClinVar allele_info -> chr:pos:REF>ALT, or None when it does not name one genomic allele.

    Returns None when the location has no alleles (the location is then the gene range or a
    CNV span), when an allele is '-' or empty (an unanchored indel needs the preceding base,
    e.g. from dbSNP SPDI), or when the span does not match the REF length.
    """
    a = allele
    if not a.chromosome or a.position_start is None or a.position_stop is None:
        return None
    ref, alt = (a.reference_allele or "").upper(), (a.alternate_allele or "").upper()
    if not _ALLELE_RE.match(ref) or not _ALLELE_RE.match(alt):
        return None
    if a.position_stop - a.position_start + 1 != len(ref):
        return None
    return variant_str(a.chromosome, a.position_start, ref, alt)


def refseq_chromosome(accession: str) -> str | None:
    """RefSeq chromosome accession -> chr name ('NC_000023.11' -> 'chrX'); None for other sequences."""
    m = re.fullmatch(r"NC_0+(\d+)\.\d+", accession.strip())
    if not m:
        return None
    number = int(m.group(1))
    named = {23: "X", 24: "Y", 12920: "M"}
    if number in named:
        return "chr" + named[number]
    return f"chr{number}" if 1 <= number <= 22 else None


def _segment_1based(s: GtfSegment) -> GtfSegment:
    start = s.start + 1
    return dataclasses.replace(s, start=start, width=s.end - start + 1)


def _exon_1based(e: GtfExon) -> GtfExon:
    if e.strand == "+" and (e.acceptor, e.donor) == (e.start, e.end):
        start = e.start + 1
        return dataclasses.replace(e, start=start, width=e.end - start + 1, acceptor=start, donor=e.end)
    if e.strand == "-" and (e.acceptor, e.donor) == (e.end, e.start):
        start = e.start + 1
        return dataclasses.replace(e, start=start, width=e.end - start + 1, acceptor=e.end, donor=start)
    raise ValueError(f"exon {e.exon_number}: acceptor/donor {e.acceptor}/{e.donor} do not follow the "
                     f"recorded gtf convention for strand {e.strand!r} (start {e.start}, end {e.end})")


def normalize_gtf(t: GtfTranscript) -> GtfTranscript:
    """Returns a copy of a raw `gtf` transcript with 1-based closed starts (apply once, to raw/).

    Convention of the raw file: `alphagenome_atlas_avi.py gtf` writes the GENCODE v46 Feather
    `Start` column as `start`, and that column is 0-based (the alphagenome package's
    `gene_annotation.extract_tss` documents "Start == End, 0-based"); `end` is the 1-based last
    base. Evidence recorded on 2026-11-12 (see references/dependency_outputs.md): every exon start
    of both reference genes is one less than the Ensembl REST 1-based start, and `width` is one
    more than the Ensembl length, while the ends agree.

    What changes:
      - transcript, exon, CDS and UTR `start` become start + 1; `width` = end - start + 1;
      - exon `acceptor`/`donor` (the 5'/3' exon boundaries): the script writes acceptor = start,
        donor = end on '+' and acceptor = end, donor = start on '-'; the value equal to the raw
        start moves with it, the one equal to `end` is kept. Any other pattern raises ValueError;
      - junction coordinates are NOT shifted: junction_start is the last base of the upstream
        exon and junction_end the last base of the intron (1-based), i.e. the 0-based half-open
        intron that Atlas SPLICE_JUNCTIONS uses, so they keep matching the Atlas scores. The
        1-based intron is junction_start + 1 .. junction_end; `intron_length` is recomputed as
        junction_end - junction_start (the raw value is one less).
    The TSS of the normalised transcript is `start` on '+' and `end` on '-'.
    """
    junctions = tuple(dataclasses.replace(j, intron_length=j.junction_end - j.junction_start)
                      for j in t.junctions)
    return dataclasses.replace(
        t, start=t.start + 1, exons=tuple(_exon_1based(e) for e in t.exons), junctions=junctions,
        cds=tuple(_segment_1based(s) for s in t.cds), utr=tuple(_segment_1based(s) for s in t.utr))


# --------------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------------

def validate(args: argparse.Namespace) -> None:
    directory = Path(args.dir)
    if not directory.is_dir():
        print(f"Error: {directory} is not a directory", file=sys.stderr)
        sys.exit(1)
    files = sorted(p for p in directory.iterdir() if p.is_file() and not p.name.startswith("."))
    errors = 0
    for path in files:
        try:
            _m, kind, _skill, _loader = _match(path)
            records = load_raw(path)
        except DependencyOutputError as exc:
            errors += 1
            print(f"FAIL file={exc.file} field={exc.field} skill={exc.skill}: {exc.reason}", file=sys.stderr)
            continue
        print(f"ok   {path.name}: {len(records)} {kind}")
    if errors:
        print(f"{errors} of {len(files)} file(s) failed validation in {directory}", file=sys.stderr)
        sys.exit(1)
    print(f"All {len(files)} file(s) in {directory} match their dataclasses")


def main() -> None:
    for stream in (sys.stdout, sys.stderr):  # UTF-8 on any console (cp1252 on Windows)
        stream.reconfigure(encoding="utf-8", errors="replace")
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("validate", help="Validate every raw file of a directory")
    p.add_argument("--dir", required=True, help="Directory of raw dependency outputs (study raw/)")
    p.set_defaults(func=validate)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
