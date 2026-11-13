"""Atlas track registry and tissue documentation (Workflow step 8, `expression.py cohorts`).

Builds the track registry from the pre-fetched Atlas track catalogue
(`raw/tracks_{scorer}_{curie}.parquet`) so later score tables can be matched to it by track name.

The catalogue has no assay, data-source or GTEx columns, so those three are read from the track
name, which the Atlas builds as `[<assay prefix>] <CURIE> [gtex <Tissue>] <assay>`:
  * `gtex <Tissue> <assay>`  -> data_source gtex, gtex_tissue <Tissue>
  * an `hCAGE` track          -> data_source fantom
  * any other track           -> data_source encode
The score tables carry the same three columns from the server; `track_name` is the join key, so
a divergence shows up when the scores are merged (Workflow step 11).

Library: parse_track_name, strand_class, atlas_locus_url, build_track_registry,
render_tissue_cohorts_markdown, REGISTRY_COLUMNS.
"""

from __future__ import annotations

import re
import urllib.parse
from collections.abc import Mapping, Sequence
from typing import Any

import polars as pl

from models import AtlasTrack

REGISTRY_COLUMNS = (
    "scorer", "track_name", "track_strand", "strand_class", "ontology_curie", "ontology_label",
    "biosample_name", "biosample_type", "assay_title", "data_source", "gtex_tissue", "cohort",
    "tissue_key", "ontology_obo_url", "ontology_ols_url", "atlas_locus_track_url", "column_prefix",
)
STRAND_ORDER = {"sense": 0, "antisense": 1, "unstranded": 2}
ATLAS_BASE_URL = "https://deepmind.google.com/science/alphagenome/atlas"
ATLAS_MODALITY_LABELS = {
    "RNA_SEQ": "RNA-seq", "DNASE": "DNase", "ATAC": "ATAC-seq", "CHIP_TF": "ChIP-TF",
    "CHIP_HISTONE": "ChIP-Histone", "CAGE": "CAGE", "PROCAP": "PRO-cap",
    "POLYADENYLATION": "Polyadenylation", "SPLICE_JUNCTIONS": "Splice junctions",
    "SPLICE_SITE_USAGE": "Splice site usage",
}
_CURIE_RE = re.compile(r"^[A-Za-z][A-Za-z0-9]*:\d+$")


def obo_purl(curie: str) -> str:
    """OBO PURL of a CURIE (PREFIX:digits -> http://purl.obolibrary.org/obo/PREFIX_digits)."""
    prefix, local = curie.split(":", 1)
    return f"http://purl.obolibrary.org/obo/{prefix}_{local}"


def ols_term_url(ontology_name: str, iri: str) -> str:
    """OLS4 page of a term."""
    return f"https://www.ebi.ac.uk/ols4/ontologies/{ontology_name}/classes?iri={urllib.parse.quote(iri, safe='')}"


def tissue_key(curie: str) -> str:
    """Lowercase CURIE with non-alphanumerics as underscores (PREFIX:digits -> prefix_digits)."""
    return re.sub(r"[^a-z0-9]+", "_", curie.lower()).strip("_")


def strand_class(track_strand: str | None, gene_strand: int) -> str:
    """sense / antisense relative to the gene, or unstranded ('.' or missing)."""
    if track_strand in (None, "", "."):
        return "unstranded"
    return "sense" if (track_strand == "+") == (gene_strand == 1) else "antisense"


def parse_track_name(track_name: str) -> dict[str, str | None]:
    """Assay title, data source and GTEx tissue from an Atlas track name (see module docstring)."""
    tokens = [t for t in track_name.split() if not _CURIE_RE.match(t)]
    gtex_tissue = None
    source = "encode"
    if len(tokens) >= 3 and tokens[0] == "gtex":
        source, gtex_tissue, tokens = "gtex", tokens[1], tokens[2:]
    elif tokens and tokens[0] == "hCAGE":
        source = "fantom"
    return {"assay_title": " ".join(tokens) or None, "data_source": source, "gtex_tissue": gtex_tissue}


def atlas_locus_url(chromosome: str, start: int, end: int, *, scorer: str, biosample: str | None) -> str:
    """Atlas locus page for a 1-based closed interval, filtered to one scorer and biosample."""
    label = ATLAS_MODALITY_LABELS.get(scorer.upper(), scorer)
    parts = ([f"BIOSAMPLE_NAME:{biosample}"] if biosample else []) + [f"SCORER_MODALITY:{label}"]
    query = {"q": f"{chromosome}:{start}-{end}", "m": "locus", "lItems": f"avi,section:{scorer.upper()}",
             "f": ",".join(parts)}
    return f"{ATLAS_BASE_URL}?" + urllib.parse.urlencode(query, quote_via=urllib.parse.quote, safe=":,-")


def build_track_registry(
    tracks: Sequence[AtlasTrack],
    study: Any,
    *,
    resolution: Mapping[str, Any],
    chromosome: str,
    gene_strand: int,
    window: tuple[int, int],
) -> pl.DataFrame:
    """One row per Atlas track of the effective cohorts and scorers.

    Args:
        tracks: Catalogue tracks (any scorer or CURIE; the others are ignored).
        study: StudyConfig; its `scorers` fix the requested scorers.
        resolution: Output of cohorts.resolve_cohorts (effective scorers, CURIE roles and labels).
        chromosome: Locus chromosome with chr prefix.
        gene_strand: 1 or -1.
        window: Locus window (start, end), 1-based closed, for the Atlas links.

    Returns:
        A frame with REGISTRY_COLUMNS, ordered by scorer, strand class, CURIE and track name;
        `column_prefix` numbers the tracks of each scorer in that order.

    Raises:
        ValueError: A track whose scorer was not requested in study.yaml.
    """
    requested = set(study.data["scorers"])
    stray = sorted({t.scorer for t in tracks} - requested)
    if stray:
        raise ValueError(f"tracks of scorer(s) {stray} are not in study.yaml scorers {sorted(requested)}")
    roles = {c["curie"]: c for c in resolution["curies"] if c["status"] == "kept"}
    keep = set(resolution["effective_scorers"])
    records: list[dict[str, Any]] = []
    for t in sorted(tracks, key=lambda t: t.track_index):
        if t.scorer not in keep or t.ontology_curie not in roles:
            continue
        info = roles[t.ontology_curie]
        records.append({
            "scorer": t.scorer,
            "track_name": t.track_name,
            "track_strand": t.strand,
            "strand_class": strand_class(t.strand, gene_strand),
            "ontology_curie": t.ontology_curie,
            "ontology_label": info["ols"]["label"],
            "biosample_name": t.biosample_name,
            "biosample_type": t.biosample_type,
            **parse_track_name(t.track_name),
            "cohort": info["role"],
            "tissue_key": tissue_key(t.ontology_curie),
            "ontology_obo_url": info["obo_url"],
            "ontology_ols_url": info["ols_url"],
            "atlas_locus_track_url": atlas_locus_url(chromosome, window[0], window[1], scorer=t.scorer,
                                                     biosample=t.biosample_name),
        })
    records.sort(key=lambda r: (r["scorer"], STRAND_ORDER[r["strand_class"]], r["ontology_curie"], r["track_name"]))
    ordinal: dict[str, int] = {}
    for r in records:
        n = ordinal.get(r["scorer"], 0)
        ordinal[r["scorer"]] = n + 1
        r["column_prefix"] = f"{r['scorer'].lower()}_t{n}_{r['tissue_key']}"
    schema = {c: pl.Utf8 for c in REGISTRY_COLUMNS}
    return pl.DataFrame(records, schema=schema).select(REGISTRY_COLUMNS)


def render_tissue_cohorts_markdown(resolution: Mapping[str, Any], registry: pl.DataFrame) -> str:
    """Per-study tissue documentation: what each cohort CURIE is (OLS) and which tracks represent it.

    Every statement comes from a recorded dependency output (OLS records, the track catalogue).
    """
    lines = [
        f"# Tissue cohorts - {resolution['study_id']}",
        "",
        "Generated from the EMBL-EBI OLS4 term records in raw/, the Atlas track catalogue and the "
        "cohort decisions of this study.",
        f"Target cohort label: **{resolution['target_label']}**; control cohort label: "
        f"**{resolution['control_label']}**.",
        "",
        "| CURIE | Tissue (OLS label) | Cohort | Status | Atlas tracks |",
        "| :--- | :--- | :--- | :--- | :--- |",
    ]
    for c in resolution["curies"]:
        n = sum(c["tracks_per_scorer"].values())
        lines.append(f"| [`{c['curie']}`]({c['ols_url']}) | {c['ols']['label']} | {c['role']} | "
                     f"{_status_text(c)} | {n} |")
    for c in resolution["curies"]:
        ols = c["ols"]
        lines += ["", f"## {ols['label']} (`{c['curie']}`) - {c['role']}", ""]
        lines.append(f"* **Definition (OLS4):** {ols['description'] or 'none provided'}")
        if ols["synonyms"]:
            lines.append(f"* **Synonyms:** {', '.join(ols['synonyms'])}")
        lines.append(f"* **Ontology:** {ols['ontology_name']}; IRI <{ols['iri']}>; [OLS]({c['ols_url']})")
        lines.append("* **Atlas tracks per scorer:** " + ", ".join(f"{s} {n}" for s, n in c["tracks_per_scorer"].items()))
        if c["status"] != "kept":
            lines.append(f"* **Status:** {_status_text(c)} (decision `cohort.unserved:{c['curie']}`).")
        rows = registry.filter(pl.col("ontology_curie") == c["curie"]).to_dicts() if c["status"] == "kept" else []
        if rows:
            lines += ["", "| Column prefix | Scorer | Strand class | Assay | Source | Atlas biosample | GTEx tissue |",
                      "| :--- | :--- | :--- | :--- | :--- | :--- | :--- |"]
            lines += [f"| `{t['column_prefix']}` | {t['scorer']} | {t['strand_class']} | {t['assay_title'] or ''} | "
                      f"{t['data_source'] or ''} | {t['biosample_name'] or ''} | {t['gtex_tissue'] or ''} |"
                      for t in rows]
    return "\n".join(lines) + "\n"


def _status_text(c: Mapping[str, Any]) -> str:
    if c["status"] == "kept":
        return f"used (added by decision `{c['added_by']}`)" if c.get("added_by") else "used"
    if c["status"] == "dropped":
        return "dropped (no tracks)"
    return f"replaced by {c['replaced_by']}"
