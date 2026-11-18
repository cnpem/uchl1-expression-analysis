"""Exploration before the study exists (Workflow step 0, `expression.py explore`).

Cheap, read-only queries that show the researcher what the data allow before anything is fixed in
study.yaml: the gene structure, the Regulatory Build around it, 3-5 candidate windows derived from those
features, the ontology terms that match each tissue or cell named, which of them have Atlas tracks (per
scorer, with data source and assay), and how strongly the gene itself is predicted to be active in each
tissue (`RNA_SEQ_ACTIVE` of one probe SNV inside the gene, the only variant scored; no variant of a window is
scored). Everything goes to `<out>/exploration/` with its own
manifest; `study.py init --from-exploration` reuses the captures without repeating a call.

Outputs (in `<out>/exploration/`):
  raw/                 files with the names the study stages read (gtf.json, regulatory.json, sequence_*.txt,
                       ols_search_<slug>.json, ols_term_*.json, tracks_<S>_catalogue.parquet)
  atlas/               atlas_tracks.tsv (whole catalogue with data_source and assay_title) and
                       scores_RNA_SEQ_ACTIVE_baseline.tsv (one SNV inside the gene, every RNA-seq track)
  exploration.json     the options in machine-readable form (read by study.py init --from-exploration)
  summary.md           the options for the researcher, in Portuguese, with a recommendation and the reason
  manifest.json        every command, its versions and outputs

Network calls: 1 GENCODE gtf, 2 Ensembl (Regulatory Build, sequence), 1 OLS search per tissue, at most 1 OLS
term per tissue whose best term has no RNA-seq track, 1 Atlas catalogue per scorer, 1 Atlas catalogue with
sources, 1 Atlas score query of one SNV.

Library: candidate_windows, term_table, baseline_activity, run.
"""

from __future__ import annotations

import json
import re
from collections import Counter
from pathlib import Path
from typing import Any

import polars as pl

import fetch
import models
import nomenclature as N
import terms as T
from core import locus as core_locus

STEP = 0
DEFAULT_SCORERS = ("RNA_SEQ", "CAGE", "DNASE", "ATAC")
DEFAULT_FLANK_BP = 2000
DEFAULT_FLOOR = 0.05  # terms whose expression level is below 5% of the gene's highest term are flagged (relative cut)
TSS_HALF_WIDTH = 500
EXPLORATION_FILE = "exploration.json"
EXPLORATION_FORMAT = 1  # format_version of exploration.json; study.py init refuses a newer one
SUMMARY_FILE = "summary.md"
ACTIVITY_FILE = "activity.md"
BASELINE_FILE = "atlas/scores_RNA_SEQ_ACTIVE_baseline.tsv"
_KIND_PT = {"tissue": "tecido", "primary_cell": "célula primária", "in_vitro_differentiated_cells": "célula diferenciada in vitro",
            "cell_line": "linhagem celular", "organoid": "organoide"}
STUDY_RAW = (re.compile(r"^gtf\.json$"), re.compile(r"^regulatory\.json$"), re.compile(r"^sequence_.+\.txt$"),
             re.compile(r"^ols_search_.+\.json$"), re.compile(r"^ols_term_.+\.json$"),
             re.compile(r"^tracks_[A-Z0-9_]+_catalogue\.parquet$"))
ATLAS_LIMIT_BP = 1000  # the Atlas scripts refuse wider windows unless --max-window is raised (the skill does it)


# What each optional capture gives, and what the exploration shows without it (Portuguese, for summary.md).
OPTIONAL_CAPTURES = {
    "regulatory": {"service": "Ensembl REST (Regulatory Build)",
                   "effect": "sem os elementos regulatórios, só a janela TSS ± 500 pb é oferecida"},
    "sequence": {"service": "Ensembl REST (sequência de referência)",
                 "effect": "sem a sequência, a variante de sondagem não pode ser montada e o nível de expressão "
                           "previsto por tecido não é mostrado"},
}


def _unavailable(capture: str, exc: Exception) -> dict[str, str]:
    """A capture that failed, as recorded in exploration.json and shown in summary.md."""
    first = next((line for line in str(exc).splitlines()[1:] if line.strip()), str(exc).splitlines()[0])
    print(f"Note: {OPTIONAL_CAPTURES[capture]['service']} could not be reached; the exploration goes on without it.")
    return {"capture": capture, **OPTIONAL_CAPTURES[capture], "error": first.strip()[:300]}


class ExploreError(Exception):
    """An exploration input is missing or inconsistent; the CLI prints it and exits with 1."""


# --------------------------------------------------------------------------- locus and windows


def _gene_locus(gtf_path: Path, gene: str) -> tuple[Any, list[Any], str | None]:
    """(selected transcript, all normalised transcripts, note when no single MANE Select transcript exists)."""
    transcripts = [models.normalize_gtf(t) for t in models.load_raw(gtf_path)]
    try:
        return core_locus.select_transcript(transcripts, gene, {}), transcripts, None
    except core_locus.LocusPending as exc:
        q = exc.items[0]
        chosen = next(t for t in transcripts if t.transcript_id == q.recommended)
        return chosen, transcripts, (f"no single MANE Select transcript; {chosen.transcript_id} (longest CDS) is used "
                                     f"here and the study asks locus.transcript ({', '.join(q.codes[:-1])})")


def search_region(transcripts: list[Any], flank: int) -> tuple[str, int, int]:
    """Gene span over every transcript +/- flank on both sides (the study's regulatory_search rule)."""
    chroms = sorted({t.chromosome for t in transcripts})
    if len(chroms) != 1:
        raise ExploreError(f"the gene's transcripts are on {chroms}; expected one contig")
    return chroms[0], max(1, min(t.start for t in transcripts) - flank), max(t.end for t in transcripts) + flank


def _features_report(features: list[Any], window: tuple[int, int], tss: int, strand: int) -> tuple[list[str], list[str]]:
    inside, outside = [], []
    for f in sorted(features, key=lambda f: f.start):
        mid = (f.start + f.end) // 2
        dist = (mid - tss) * strand
        text = f"{f.description} {f.id} ({f.start}-{f.end}, {'+' if dist >= 0 else ''}{dist} pb do TSS)"
        if f.start >= window[0] and f.end <= window[1]:
            inside.append(text)
        elif f.end >= window[0] and f.start <= window[1]:
            inside.append(text + " [parcial]")
        else:
            outside.append(text)
    return inside, outside


def candidate_windows(tx: Any, features: list[Any], region: tuple[str, int, int], flank: int) -> list[dict[str, Any]]:
    """3-5 windows derived from the features around the TSS, each with how the study would express it.

    Generic rules, in order (duplicates dropped): the regulatory unit the locus stage builds by default
    (feature chain at the TSS through exon 1); that chain with gaps up to the search flank bridged; the
    TSS feature(s) alone; TSS +/- 500 bp; every promoter/enhancer of the search region with exon 1.
    """
    chrom = tx.chromosome
    strand = core_locus._strand(tx)
    tss = core_locus._tss(tx)
    ex1 = core_locus._first_exon(tx)
    biotypes = core_locus.DEFAULT_BIOTYPES
    flank_gap = flank
    out: list[dict[str, Any]] = []

    def add(wid: str, label: str, lo: int, hi: int, policy: dict[str, Any], why: str, questions: list[str]) -> None:
        lo, hi = min(lo, tss), max(hi, tss)
        if any(w["start"] == lo and w["end"] == hi for w in out):
            return
        bp = hi - lo + 1
        inside, outside = _features_report(features, (lo, hi), tss, strand)
        out.append({"id": wid, "label": label, "region": f"{chrom}:{lo}-{hi}", "start": lo, "end": hi, "bp": bp,
                    "snvs": 3 * bp, "policy": policy, "why": why, "study_questions": questions,
                    "features_included": inside, "features_excluded": outside,
                    "needs_max_window": bp > ATLAS_LIMIT_BP})

    unit, chain, triggers = core_locus.derive_regulatory_unit_window(tx, features, biotypes, 0, chrom)
    if unit is not None:
        add("regulatory_unit", "Unidade regulatória do TSS + éxon 1", unit.start, unit.end,
            {"policy": "regulatory_unit", "max_feature_gap_bp": 0},
            "promotor/enhancer que contém o TSS, estendido por elementos contíguos a montante, até o fim do éxon 1 "
            "(padrão da skill)", [q.key for q in triggers])
    bridged, chain2, trig2 = core_locus.derive_regulatory_unit_window(tx, features, biotypes, flank_gap, chrom)
    if bridged is not None:
        add("unit_bridged", "Unidade regulatória com lacunas a montante incluídas", bridged.start, bridged.end,
            {"policy": "regulatory_unit", "max_feature_gap_bp": flank_gap},
            "inclui promotores/enhancers a montante separados por lacunas (até a borda da região buscada)",
            [q.key for q in trig2])
    if chain:
        anchor = [f for f in chain if f.start <= tss + 1 and f.end >= tss - 1] or chain[-1:]
        add("tss_feature", "Só o elemento que contém o TSS", min(f.start for f in anchor), max(f.end for f in anchor),
            {"policy": "explicit"}, "o núcleo do promotor, menor custo; deixa de fora o resto do éxon 1", [])
    add("tss_500", f"TSS ± {TSS_HALF_WIDTH} pb", tss - TSS_HALF_WIDTH, tss + TSS_HALF_WIDTH, {"policy": "explicit"},
        "janela fixa, independente da anotação regulatória (útil quando não há elemento no TSS)", [])
    pe = [f for f in features if f.description in biotypes]
    if pe:
        add("all_features", "Todos os promotores/enhancers da região buscada + éxon 1",
            min([f.start for f in pe] + [ex1.start]), max([f.end for f in pe] + [ex1.end]), {"policy": "explicit"},
            "cobre todos os elementos regulatórios anotados em volta do gene; custo maior", [])
    for w in out:
        if w["policy"]["policy"] == "explicit":
            w["policy"]["explicit_interval"] = w["region"]
    return out


# --------------------------------------------------------------------------- terms and tracks


def _norm(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", str(text).lower()).strip()


def _same_name(a: str, b: str) -> bool:
    """Same name written another way: case, hyphens and spaces do not count ("MCF-7" = "MCF7", "T-cell" = "T cell")."""
    return _norm(a).replace(" ", "") == _norm(b).replace(" ", "")


def _catalogue(atlas_tracks: Path, scorers: tuple[str, ...]) -> pl.DataFrame:
    df = pl.read_csv(atlas_tracks, separator="\t", infer_schema_length=0)
    for col in ("assay_title", "data_source", "biosample_name", "biosample_type"):
        if col not in df.columns:
            df = df.with_columns(pl.lit(None, dtype=pl.Utf8).alias(col))
    return df.filter(pl.col("scorer").is_in(list(scorers)) & pl.col("ontology_curie").is_not_null()
                     & (pl.col("ontology_curie") != ""))


def _served(cat: pl.DataFrame, curie: str, scorers: tuple[str, ...]) -> dict[str, Any]:
    rows = cat.filter(pl.col("ontology_curie") == curie)
    per = {}
    for s in scorers:
        r = rows.filter(pl.col("scorer") == s)
        if r.height:
            per[s] = {"tracks": r.height,
                      "sources": dict(Counter(x or "?" for x in r["data_source"].to_list())),
                      "assays": dict(Counter(x or "?" for x in r["assay_title"].to_list()))}
    names = sorted({n for n in rows["biosample_name"].to_list() if n})
    return {"per_scorer": per, "biosample_names": names}


def term_table(ols_file: Path, text: str, cat: pl.DataFrame, scorers: tuple[str, ...], top: int = 6) -> dict[str, Any]:
    """OLS candidates for one tissue text with their Atlas tracks, the recommended term and catalogue terms whose
    biosample name contains the text (catches cell lines and subtypes the search ranks low)."""
    hits = models.load_raw(ols_file)[0].terms
    seen: dict[str, Any] = {}
    for h in hits:
        if h.obo_id and h.obo_id not in seen and re.match(r"^[A-Za-z]+:\d+$", h.obo_id):
            seen[h.obo_id] = h
    candidates = []
    for curie, h in list(seen.items())[:top]:
        candidates.append({"curie": curie, "label": h.label, "ontology": curie.split(":")[0],
                           "definition": (h.description[0] if h.description else ""), **_served(cat, curie, scorers)})
    want = _norm(text)
    # The Atlas catalogue names its samples; a sample with the researcher's own name (written any way) is the term
    # they mean, even when the ontology search ranks it low, misses it, or knows it no more (cell lines).
    in_cat = {}
    for curie, name in cat.select("ontology_curie", "biosample_name").unique().iter_rows():
        if name and _same_name(name, text) and curie not in in_cat:
            in_cat[curie] = name
    for curie, name in in_cat.items():
        if curie not in seen:
            candidates.insert(0, {"curie": curie, "label": name, "ontology": curie.split(":")[0],
                                  "definition": "", "source": "atlas_catalogue", **_served(cat, curie, scorers)})
    by_cat = [c for c in candidates if c["curie"] in in_cat and c["per_scorer"].get("RNA_SEQ")]
    exact = [c for c in candidates if _same_name(c["label"], text) and c["per_scorer"].get("RNA_SEQ")] or \
        [c for c in candidates if _norm(c["label"]) == want]
    recommended = (by_cat or exact or [c for c in candidates if c["per_scorer"].get("RNA_SEQ")] or candidates or [None])[0]
    words = want.split()
    related = []
    if words:
        names = cat.select("ontology_curie", "biosample_name").unique()
        for curie, name in names.iter_rows():
            if curie in seen or curie in in_cat or not name:
                continue
            if all(w in _norm(name).split() for w in words):
                related.append({"curie": curie, "label": name, **_served(cat, curie, scorers)})
    related.sort(key=lambda r: (-len(r["per_scorer"]), -sum(v["tracks"] for v in r["per_scorer"].values()), r["curie"]))
    return {"text": text, "ols_file": ols_file.name, "candidates": candidates,
            "recommended": recommended["curie"] if recommended else None, "related_in_catalogue": related[:8],
            "descendants": []}


def served_descendants(term_file: Path, cat: pl.DataFrame, scorers: tuple[str, ...]) -> list[dict[str, Any]]:
    term = models.load_raw(term_file)[0].term
    rel = term.hierarchicalDescendants if term.hierarchicalDescendants is not None else term.hierarchicalChildren
    out = []
    for r in rel or ():
        info = _served(cat, r.obo_id, scorers)
        if info["per_scorer"]:
            out.append({"curie": r.obo_id, "label": r.label, **info})
    out.sort(key=lambda d: (-len(d["per_scorer"]), -d["per_scorer"].get("RNA_SEQ", {}).get("tracks", 0), d["label"]))
    return out


# --------------------------------------------------------------------------- baseline expression


def baseline_variant(tx: Any, sequence_file: Path) -> str:
    """One SNV in the middle of exon 1 (inside the gene), REF from the Ensembl sequence."""
    seq = models.load_raw(sequence_file)[0]
    ex1 = core_locus._first_exon(tx)
    pos = (ex1.start + ex1.end) // 2
    ref = seq.sequence[pos - seq.start].upper()
    if ref not in "ACGT":
        raise ExploreError(f"reference base at {tx.chromosome}:{pos} is {ref!r}; cannot build a baseline SNV")
    alt = next(b for b in "ACGT" if b != ref)
    return f"{tx.chromosome}:{pos}:{ref}>{alt}"


def baseline_activity(path: Path, floor: float) -> dict[str, Any]:
    """Mean RNA_SEQ_ACTIVE per ontology term, relative to the gene's highest term."""
    df = pl.read_csv(path, separator="\t", infer_schema_length=10000)
    df = df.filter(pl.col("raw_score").is_not_null())
    if "ontology_curie" not in df.columns or df.height == 0:
        raise ExploreError(f"{path.name}: no RNA_SEQ_ACTIVE score returned")
    names = df.group_by("ontology_curie").agg(pl.col("biosample_name").first().alias("name"),
                                              pl.col("raw_score").mean().alias("activity"), pl.len().alias("tracks"))
    top = names["activity"].max() or 0.0
    names = names.with_columns((pl.col("activity") / top if top > 0 else pl.lit(0.0)).alias("relative"))
    by = {r["ontology_curie"]: {"name": r["name"], "activity": r["activity"], "relative": r["relative"],
                                "tracks": r["tracks"], "expressed": r["relative"] >= floor}
          for r in names.iter_rows(named=True)}
    ranked = sorted(by.items(), key=lambda kv: -kv[1]["activity"])
    return {"max": top, "floor": floor, "by_curie": by, "top": [c for c, _ in ranked[:10]],
            "median_relative": float(names["relative"].median() or 0.0)}


def write_activity(root: Path, floor: float = 0.05) -> Path:
    """exploration/activity.md: every term of the Atlas catalogue with the gene's predicted expression level, highest
    first, with its kind of sample, assays, sources and number of measurements. The list to choose (or change)
    the target and control groups from, without new captures or ad-hoc code."""
    root = Path(root)
    df = pl.read_csv(root / BASELINE_FILE, separator="\t", infer_schema_length=10000).filter(pl.col("raw_score").is_not_null())
    if df.height == 0:
        raise ExploreError(f"{BASELINE_FILE}: no RNA_SEQ_ACTIVE score returned")
    t = (df.group_by("ontology_curie")
         .agg(pl.col("biosample_name").first().alias("name"), pl.col("biosample_type").first().alias("kind"),
              pl.col("raw_score").mean().alias("activity"), pl.len().alias("tracks"),
              pl.col("assay_title").drop_nulls().unique().sort().str.join(" + ").alias("assays"),
              pl.col("data_source").drop_nulls().unique().sort().str.join(", ").alias("sources"))
         .sort("activity", descending=True))
    top = float(t["activity"].max() or 0.0)
    variant = df["variant"][0] if "variant" in df.columns else "?"
    L = [f"# {T.cap(T.label('expression_level'))} do gene em todos os tecidos e células do catálogo ({t.height})", "",
         f"{T.cap(T.label('expression_level'))} do gene: pontuação `RNA_SEQ_ACTIVE` de uma {T.label('probe_snv')} ({variant}), isto é, o maior "
         "entre os alelos REF e ALT da cobertura média de RNA-seq prevista nos éxons do gene. A tabela traz a média das "
         f"trilhas de cada termo e o valor relativo ao termo de maior média ({N.dec(top)}). Abaixo de {_pct(floor)} desse "
         f"máximo, a coluna *≥ {_pct(floor)} do máximo* diz **não**: é um corte desta skill, relativo ao máximo do próprio gene e sem "
         "calibração; não quer dizer que o transcrito esteja ausente. Onde o valor é muito baixo, a mudança em % de uma "
         "variante pode ser grande mesmo com uma mudança absoluta mínima. Use esta lista para escolher, ou trocar depois "
         "dos resultados, os tecidos do alvo e do controle; um termo novo precisa de captura (passos 7 a 10).", "",
         "O valor é igual ao da sequência de referência quando a variante de sondagem não aumenta a expressão, e maior "
         "caso contrário; a exploração não mede essa diferença. Para os tecidos que já estão num estudo, `summarize` "
         "mostra o nível com o alelo de referência (REF), reconstruído das duas pontuações. As trilhas de um termo podem "
         "ser de ensaios (polyA, total) e fontes (ENCODE, GTEx) diferentes, cujos valores não são diretamente "
         "comparáveis: as colunas *Ensaios* e *Fontes* dizem o que entrou em cada média.", "",
         f"| # | Termo | Nome | Tipo | {T.cap(T.short('expression_level'))} | {T.cap(T.short('share_of_maximum'))} | "
         f"≥ {_pct(floor)} do máximo | Ensaios | Fontes | {T.cap(T.label('track'))}s |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for i, r in enumerate(t.iter_rows(named=True), 1):
        rel = r["activity"] / top if top > 0 else 0.0
        L.append(f"| {i} | {r['ontology_curie']} | {r['name']} | {_KIND_PT.get(r['kind'], r['kind'] or '—')} | "
                 f"{N.dec(r['activity'])} | {_pct(rel)} | {'sim' if rel >= floor else 'não'} | {r['assays'] or '—'} | "
                 f"{r['sources'] or '—'} | {r['tracks']} |")
    out = root / ACTIVITY_FILE
    out.write_text("\n".join(L) + "\n", encoding="utf-8")
    return out


# --------------------------------------------------------------------------- run


def run(out_dir: Path, gene: str, tissues: list[str], scorers: tuple[str, ...] = DEFAULT_SCORERS,
        flank: int = DEFAULT_FLANK_BP, floor: float = DEFAULT_FLOOR, invoked_by: str = "",
        dry_run: bool = False) -> Path:
    """Runs the exploration into <out_dir>/exploration and returns the summary path."""
    if "RNA_SEQ" not in scorers:
        raise ExploreError("RNA_SEQ must be among the scorers: it is the expression contrast of the study.")
    root = Path(out_dir) / "exploration"
    root.mkdir(parents=True, exist_ok=True)
    ws = fetch.Workspace(root=root, study_id=f"{gene.lower()}-exploration", invoked_by=invoked_by, dry_run=dry_run)
    raw = root / "raw"

    def have(name: str) -> bool:  # captures already made are reused (re-running explore costs nothing)
        return (root / name).is_file() and not dry_run

    if not have("raw/gtf.json"):
        fetch.gtf(ws, gene)
    if dry_run:
        print("(dry run: the remaining commands need the gene structure and are not shown)")
        return root / SUMMARY_FILE
    tx, transcripts, tx_note = _gene_locus(raw / "gtf.json", gene)
    chrom, s0, e0 = search_region(transcripts, flank)
    region = models.interval_str(chrom, s0, e0)
    # Ensembl (Regulatory Build and reference sequence) is optional here: when it cannot be reached, the
    # exploration still shows the terms, the tracks and the TSS window, says what is missing, and a later
    # explore fills it in (captures already made are reused). A study needs both, so init refuses until then.
    unavailable: list[dict[str, str]] = []
    if not have("raw/regulatory.json"):
        try:
            fetch.regulatory(ws, region)
        except fetch.FetchError as exc:
            (raw / "regulatory.json").unlink(missing_ok=True)  # a partial file would be reused as a capture
            unavailable.append(_unavailable("regulatory", exc))
    seq_name = f"raw/sequence_{region.replace(':', '_').replace('-', '_')}.txt"
    if not have(seq_name):
        try:
            fetch.sequence(ws, region)
        except fetch.FetchError as exc:
            (root / seq_name).unlink(missing_ok=True)
            unavailable.append(_unavailable("sequence", exc))
    contig = core_locus._contig_without_prefix(chrom)
    features = ([f for f in models.load_raw(raw / "regulatory.json") if f.seq_region_name == contig]
                if (raw / "regulatory.json").is_file() else [])
    windows = candidate_windows(tx, features, (chrom, s0, e0), flank)

    for s in scorers:
        if not have(f"raw/tracks_{s}_catalogue.parquet"):
            fetch.catalogue(ws, s)
    if not have("atlas/atlas_tracks.tsv"):
        fetch.atlas_tracks(ws)
    cat = _catalogue(root / "atlas" / "atlas_tracks.tsv", tuple(scorers))

    terms = []
    for text in tissues:
        f = f"raw/ols_search_{fetch.slug(text)}.json"
        if not have(f):
            fetch.ols_search(ws, text)
        table = term_table(root / f, text, cat, tuple(scorers))
        rec = next((c for c in table["candidates"] if c["curie"] == table["recommended"]), None)
        if rec is not None and not rec["per_scorer"].get("RNA_SEQ"):
            tf = f"raw/ols_term_{fetch.curie_file(rec['curie'])}.json"
            try:  # a term the ontology no longer knows must not stop the other tissues
                if not have(tf):
                    fetch.ols_term(ws, rec["curie"])
                table["descendants"] = served_descendants(root / tf, cat, tuple(scorers))
            except (fetch.FetchError, models.DependencyOutputError, OSError, ValueError) as exc:
                table["descendants_error"] = f"{rec['curie']}: {str(exc).splitlines()[0][:200]}"
        terms.append(table)

    variant = baseline_variant(tx, root / seq_name) if (root / seq_name).is_file() else None
    activity = None  # the probe SNV needs the reference base, which comes from the sequence
    if variant is not None:
        if not have("atlas/scores_RNA_SEQ_ACTIVE_baseline.tsv"):
            fetch.baseline(ws, variant, gene)
        activity = baseline_activity(root / "atlas" / "scores_RNA_SEQ_ACTIVE_baseline.tsv", floor)

    rec_curies = [t["recommended"] for t in terms if t["recommended"]]
    rec_tracks = {s: sum(cat.filter((pl.col("ontology_curie") == c) & (pl.col("scorer") == s)).height for c in rec_curies)
                  for s in scorers}
    for w in windows:
        w["rows_per_scorer"] = {s: w["snvs"] * n for s, n in rec_tracks.items()}
    data = {
        "format_version": EXPLORATION_FORMAT,
        "gene": gene, "transcript": tx.transcript_id, "transcript_note": tx_note, "chromosome": chrom,
        "strand": tx.strand, "tss": core_locus._tss(tx), "gene_span": [min(t.start for t in transcripts),
                                                                        max(t.end for t in transcripts)],
        "n_transcripts": len(transcripts), "search_region": region, "flank_bp": flank, "scorers": list(scorers),
        "features": [{"id": f.id, "biotype": f.description, "start": f.start, "end": f.end,
                      "distance_to_tss": ((f.start + f.end) // 2 - core_locus._tss(tx)) * core_locus._strand(tx)}
                     for f in sorted(features, key=lambda f: f.start)],
        "windows": windows, "recommended_window": windows[0]["id"] if windows else None,
        "terms": terms, "baseline_variant": variant, "activity": activity,
        "recommended_tracks_per_scorer": rec_tracks, "unavailable": unavailable,
    }
    (root / EXPLORATION_FILE).write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    summary = root / SUMMARY_FILE
    summary.write_text(render_summary(data), encoding="utf-8")
    if activity is not None:
        write_activity(root, floor)
    return summary


# --------------------------------------------------------------------------- summary (Portuguese)


def _pct(x: float | None) -> str:
    return "—" if x is None else f"{100 * x:.1f}%".replace(".", ",")


def _n(x: int) -> str:
    """Integer with the Portuguese thousands separator."""
    return f"{x:,}".replace(",", ".")


def _scorer_cell(info: dict[str, Any], s: str) -> str:
    v = info["per_scorer"].get(s)
    if not v:
        return "0"
    src = "/".join(sorted(v["sources"]))
    return f"{v['tracks']} ({src})"


def _assays(info: dict[str, Any]) -> str:
    a = info["per_scorer"].get("RNA_SEQ", {}).get("assays", {})
    return ", ".join(f"{k} ×{n}" for k, n in sorted(a.items())) or "—"


def _act(activity: dict[str, Any] | None, curie: str) -> str:
    if activity is None:
        return "indisponível (sem a sequência de referência)"
    a = activity["by_curie"].get(curie)
    if a is None:
        return "sem trilha de RNA-seq"
    flag = "" if a["expressed"] else f" ⚠ abaixo de {_pct(activity['floor'])} do máximo do gene"
    return f"{N.dec(a['activity'])} ({_pct(a['relative'])} do máximo){flag}"


def render_summary(d: dict[str, Any]) -> str:
    act = d["activity"]
    scorers = d["scorers"]
    L: list[str] = []
    L.append(f"# Exploração do gene {N.gene(d['gene'])}: opções antes de definir o estudo\n")
    L.append("Esta página resume consultas **só de leitura** (estrutura do gene, elementos regulatórios, termos de "
             "ontologia, catálogo de trilhas do AlphaGenome Atlas e o nível de expressão previsto do gene em cada tecido). "
             "Só uma variante de sondagem foi pontuada, para ler o nível do gene em cada tecido; nenhuma variante da "
             "janela foi pontuada ainda. Com base nela, o pesquisador escolhe a janela, os tecidos alvo e "
             "controle e os tipos de predição (scorers); só então o estudo é criado. Tudo aqui é predição de modelo, para "
             "priorizar experimentos, não medida nem diagnóstico.\n")
    L.append("## 1. O gene\n")
    L.append(f"- Transcrito usado: **{d['transcript']}** (fita {d['strand']}), TSS em {d['chromosome']}:{d['tss']}; "
             f"o gene tem {d['n_transcripts']} transcrito(s) em GENCODE v46 (a anotação das predições do Atlas).")
    if d["transcript_note"]:
        L.append(f"- Atenção: {d['transcript_note']}.")
    L.append(f"- Região buscada para elementos regulatórios: {d['search_region']} (o gene ± {d['flank_bp']} pb).\n")
    missing = {u["capture"]: u for u in d.get("unavailable", [])}
    if missing:
        L.append("> **Capturas que faltam.** " + " ".join(
            f"{u['service']} não respondeu ({u['error'] or 'sem detalhe'}): {u['effect']}." for u in missing.values())
            + " O estudo só pode ser criado quando elas existirem: rode a exploração de novo mais tarde, com o mesmo "
              "comando (o que já foi capturado é reaproveitado).\n")
    L.append("## 2. Elementos regulatórios em volta do gene (Ensembl Regulatory Build)\n")
    if "regulatory" in missing:
        L.append("Indisponível nesta exploração (ver o aviso acima).")
    else:
        L.append("| Elemento | Tipo | Posição | Distância ao TSS (− = a montante) |\n|---|---|---|---|")
        for f in d["features"]:
            L.append(f"| {f['id']} | {f['biotype']} | {f['start']}-{f['end']} | {f['distance_to_tss']:+d} pb |")
    L.append("")
    L.append("## 3. Janelas possíveis\n")
    L.append("Cada posição da janela gera 3 variantes (SNVs). O custo cresce com o número de SNVs × trilhas: a coluna "
             "\"linhas\" estima o tamanho da consulta ao Atlas para os termos recomendados abaixo. Janelas acima de "
             f"{ATLAS_LIMIT_BP} pb são aceitas (a skill ajusta o limite), mas custam proporcionalmente mais.\n")
    L.append("| Opção | Região | Tamanho | SNVs | Linhas por scorer (termos recomendados) | Inclui | Deixa de fora |\n"
             "|---|---|---|---|---|---|---|")
    for w in d["windows"]:
        rows = ", ".join(f"{s} {_n(n)}" for s, n in w["rows_per_scorer"].items() if n)
        inc = "; ".join(w["features_included"]) or "nenhum elemento anotado"
        exc = "; ".join(w["features_excluded"]) or "—"
        star = " **(recomendada)**" if w["id"] == d["recommended_window"] else ""
        L.append(f"| `{w['id']}`{star}: {w['label']} | {w['region']} | {_n(w['bp'])} pb | {_n(w['snvs'])} | {rows or '—'} "
                 f"| {inc} | {exc} |")
    if d["windows"]:
        w0 = d["windows"][0]
        why = (" É a menor janela que cobre o elemento onde a transcrição começa e o éxon 1 inteiro, que é onde variantes "
               "reguladoras do promotor costumam agir." if w0["id"] == "regulatory_unit" else "")
        if "regulatory" in missing:
            why = (" Recomendação provisória: sem o Regulatory Build, é a única janela que pode ser montada; as janelas "
                   "baseadas nos elementos regulatórios aparecem quando a exploração for refeita.")
        L.append(f"\n**Recomendação:** `{w0['id']}` — {w0['why']}.{why}")
        if w0["study_questions"]:
            n = len(w0["study_questions"])
            L.append(f"Ao criar o estudo com ela, {'uma pergunta' if n == 1 else f'{n} perguntas'} sobre a extensão da "
                     "janela (incluir ou não elementos reguladores vizinhos) será feita na etapa da janela.")
    L.append("")
    L.append("## 4. Tecidos e células: termos de ontologia e trilhas disponíveis\n")
    L.append("Para cada texto citado, os termos candidatos (UBERON = tecido, CL = tipo celular, EFO/CLO = linhagens "
             "celulares), com o número de trilhas do Atlas por scorer (fonte entre parênteses), o protocolo de RNA-seq e "
             "o nível de expressão previsto do gene naquele termo. Um termo sem trilha não pode entrar na comparação; nesses "
             "casos aparecem os termos mais específicos (descendentes) que têm trilhas.\n")
    for t in d["terms"]:
        L.append(f"### \"{t['text']}\"\n")
        L.append("| Termo | Nome | " + " | ".join(scorers) + " | RNA-seq (protocolo) | Nível de expressão previsto |\n"
                 "|---|---|" + "---|" * len(scorers) + "---|---|")
        for c in t["candidates"]:
            star = " **(recomendado)**" if c["curie"] == t["recommended"] else ""
            L.append(f"| {c['curie']}{star} | {c['label']} | " + " | ".join(_scorer_cell(c, s) for s in scorers)
                     + f" | {_assays(c)} | {_act(act, c['curie'])} |")
        rec = next((c for c in t["candidates"] if c["curie"] == t["recommended"]), None)
        if rec is not None:
            n_rna = rec["per_scorer"].get("RNA_SEQ", {}).get("tracks", 0)
            reason = ("é a amostra do catálogo do Atlas com o nome citado" if rec.get("source") == "atlas_catalogue" else
                      "o nome coincide com o texto citado" if _same_name(rec["label"], t["text"]) else
                      "é o primeiro resultado com trilhas de RNA-seq")
            L.append(f"\nRecomendado: **{rec['curie']} ({rec['label']})**, porque {reason}"
                     + (f"; tem {n_rna} trilha(s) de RNA-seq" if n_rna else "; **não tem trilha de RNA-seq**") + ".")
            if n_rna == 1:
                L.append("Só uma trilha de RNA-seq: não há outra trilha do mesmo termo para comparar; considere juntar "
                         "termos irmãos ou descendentes da tabela.")
        if t.get("descendants_error"):
            L.append(f"\nA ontologia não devolveu os termos mais específicos do recomendado ({t['descendants_error']}); "
                     "os outros tecidos seguiram normalmente.")
        if t["descendants"]:
            L.append("\nDescendentes com trilhas (substitutos possíveis do termo recomendado):\n")
            L.append("| Termo | Nome | " + " | ".join(scorers) + " | Nível de expressão previsto |\n|---|---|"
                     + "---|" * len(scorers) + "---|")
            for x in t["descendants"][:10]:
                L.append(f"| {x['curie']} | {x['label']} | " + " | ".join(_scorer_cell(x, s) for s in scorers)
                         + f" | {_act(act, x['curie'])} |")
        if t["related_in_catalogue"]:
            L.append("\nOutros termos do catálogo cujo nome contém o texto (subtipos, linhagens):\n")
            L.append("| Termo | Nome no catálogo | " + " | ".join(scorers) + " | Nível de expressão previsto |\n|---|---|"
                     + "---|" * len(scorers) + "---|")
            for x in t["related_in_catalogue"]:
                L.append(f"| {x['curie']} | {x['label']} | " + " | ".join(_scorer_cell(x, s) for s in scorers)
                         + f" | {_act(act, x['curie'])} |")
        L.append("")
    L.append(f"## 5. {T.cap(T.label('expression_level'))} do gene nesses tecidos (`RNA_SEQ_ACTIVE`)\n")
    if act is None:
        L.append("Indisponível nesta exploração: a variante de sondagem precisa da base de referência, que vem da "
                 "sequência do Ensembl (ver o aviso no início). A coluna de nível da seção 4 também fica vazia.\n")
    else:
        L.append(f"Uma consulta ao Atlas com uma {T.label('probe_snv')} dentro do gene ({d['baseline_variant']}) devolve, para "
                 f"cada trilha de RNA-seq, o {T.label('expression_level')} do gene: o maior entre os dois alelos dessa variante. "
                 "Os valores abaixo são relativos ao termo de maior nível médio do gene no catálogo "
                 f"({N.dec(act['max'])}). Abaixo de {_pct(act['floor'])} desse máximo o termo é marcado: é um corte desta skill, "
                 "relativo ao máximo do próprio gene e sem calibração; não quer dizer que o transcrito esteja ausente. Onde o "
                 "nível é muito baixo, a mudança em % de uma variante pode ser grande mesmo com uma mudança absoluta mínima, e "
                 "não deve ser lida como efeito específico do tecido. As trilhas de um termo podem ser de ensaios e fontes "
                 "diferentes, cujos valores não são diretamente comparáveis.\n")
        L.append(f"A lista completa, com os {len(act['by_curie'])} tecidos e células do catálogo, tipo de amostra, ensaios e "
                 f"fontes, está em `{ACTIVITY_FILE}`.\n")
        L.append(f"| Termos de maior {T.label('expression_level')} | {T.cap(T.short('expression_level'))} | "
                 f"{T.cap(T.short('share_of_maximum'))} |\n|---|---|---|")
        for c in act["top"]:
            a = act["by_curie"][c]
            L.append(f"| {c} ({a['name']}) | {N.dec(a['activity'])} | {_pct(a['relative'])} |")
        low = [(t["text"], t["recommended"]) for t in d["terms"]
               if t["recommended"] in act["by_curie"] and not act["by_curie"][t["recommended"]]["expressed"]]
        L.append("")
        if low:
            L.append(f"**Aviso:** abaixo de {_pct(act['floor'])} do máximo do gene: "
                     + "; ".join(f"\"{txt}\" ({c})" for txt, c in low)
                     + ". Se algum deles for o tecido alvo, uma mudança grande em % ali pode ser uma mudança absoluta mínima "
                       "(e o modelo erra mais na especificidade de tecido para genes de nível baixo); considere um "
                       "termo da lista acima com nível mais alto. Como controle, um termo de nível muito baixo pode dominar o "
                       "valor do grupo com mudanças grandes em % que são mudanças absolutas mínimas; o painel tem um piso do "
                       "nível REF para tirá-lo do cálculo.")
        else:
            L.append(f"Todos os termos recomendados com RNA-seq estão acima de {_pct(act['floor'])} do máximo do gene.")
        L.append("")
    L.append(f"## 6. {T.cap(T.label('scorer'))} (scorers)\n")
    L.append("| Scorer | O que mede | Termos recomendados com trilha |\n|---|---|---|")
    meaning = {"RNA_SEQ": "expressão (RNA)", "CAGE": "início de transcrição", "DNASE": "cromatina acessível (DNase)",
               "ATAC": "cromatina acessível (ATAC)"}
    for s in scorers:
        served = [t["recommended"] for t in d["terms"] if t["recommended"] and any(
            c["curie"] == t["recommended"] and c["per_scorer"].get(s) for c in t["candidates"])]
        L.append(f"| {s} | {meaning.get(s, s)} | {len(served)} de {len(d['terms'])} ({', '.join(served) or '—'}) |")
    L.append("\n**Recomendação:** RNA_SEQ é obrigatório (é o contraste de expressão). Os demais valem a pena quando "
             "servem ao menos um termo do alvo e um do controle; um tipo de previsão sem trilhas em um dos grupos é "
             "deixado de fora automaticamente e registrado no relatório.\n")
    L.append("## 7. Limites\n")
    L.append("- O Atlas só tem **substituições de uma base (SNVs)** do genoma de referência GRCh38. Inserções e deleções "
             "(indels) conhecidas, ou uma lista fornecida pelo pesquisador, podem ser pontuadas pelo modelo "
             "(`score_variants.py` da skill `alphagenome`, API do modelo, ~10-20 s por variante); varredura saturada "
             "de indels fica fora por custo.\n"
             "- Só a janela escolhida é mutada: elementos distantes (enhancers a dezenas de kb, íntrons) ficam fora.\n"
             "- Termos com uma única trilha não têm outra para comparar. Duas trilhas do mesmo termo diferem em ensaio "
             "(polyA, total) ou fonte (ENCODE, GTEx): não são réplicas. Essas diferenças entre alvo e controle são "
             "declaradas no relatório.\n"
             "- Tudo é predição de modelo para priorizar experimentos; não é medida, nem classificação clínica.\n")
    L.append("## 8. O que o pesquisador decide agora\n")
    L.append("1. A janela (seção 3).\n2. Quais termos formam o grupo **alvo** e quais o **controle** (seção 4), e um "
             "nome curto para cada grupo.\n3. Os scorers (seção 6).\n\nDepois disso o estudo é criado com "
             "`study.py init --from-exploration`, reaproveitando todos os arquivos já baixados.\n")
    L.append("Arquivos desta exploração: `raw/` (respostas das skills de dependência), `atlas/` (catálogo e nível de expressão "
             "previsto), `exploration.json`, `manifest.json` (comandos exatos e versões).\n")
    return "\n".join(L)


def study_raw_files(root: Path, scorers: list[str]) -> list[Path]:
    """Exploration captures a study reads, limited to the catalogues of the chosen scorers."""
    out = []
    for p in sorted((root / "raw").iterdir()):
        m = re.match(r"^tracks_([A-Z0-9_]+)_catalogue\.parquet$", p.name)
        if m and m.group(1) not in scorers:
            continue
        if any(r.match(p.name) for r in STUDY_RAW):
            out.append(p)
    return out
