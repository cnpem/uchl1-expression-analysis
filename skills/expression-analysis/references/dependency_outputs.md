# Dependency outputs: expression-analysis

What each file in a study's `raw/` folder is, which dependency skill writes it, and the fields
`scripts/models.py` expects. Field names are the ones the dependency writes; `models.py validate
--dir <study>/raw` fails with file, field and skill when a dependency adds, drops or renames one.
Formats were recorded on 2026-11-12 (reference study `assets/reference/uchl1/`, manifest there).

Placeholders: `G` gene symbol, `R` regulatory search region (gene span ± `window.regulatory_search_upstream_bp`),
`W` analysis window `chr:start-end` (from `processed/locus.json`), `M` = `study.py window --print-max-window`,
`C` CURIE, `S` Atlas scorer. Each file is written by `expression.py fetch <kind>` (or `explore`), which
runs the command shown with the pinned packages of `scripts/deps.py` and records it in `manifest.json`;
the commands below show what runs, not something to type.

## Dependency versions tested

Versions are the `computedHash` in `skills-lock.json` (the skills have no release numbers) plus the
version fields their own files declare.

| Skill | Source | Lock hash | Version declared | Runtime used on 2026-11-12 |
|---|---|---|---|---|
| `alphagenome-variant-impact-score` | google-deepmind/science-skills | `e9980588bc1e…` | SKILL.md `metadata.publisher: google` (no version); `pyproject.toml` 0.1.0, `requires-python >=3.14` | Python 3.13 (installed), `alphagenome` 2.0.0 |
| `alphagenome` | K-Dense-AI/scientific-agent-skills | `ff10c6aa5025…` | SKILL.md `metadata.version: "1.2"`, upstream `alphagenome 0.9.0` | PEP 723 `alphagenome>=0.9.0` → 2.0.0 |
| `embl-ebi-ols` | google-deepmind/science-skills | `3c05b0a57043…` | publisher google | PEP 723 |
| `ensembl-database` | google-deepmind/science-skills | `da1c6d0f87ef…` | publisher google | PEP 723; Ensembl REST release 116 |
| `alphagenome-single-variant-analysis` | google-deepmind/science-skills | `c61d75fd3f60…` | publisher google; `pyproject.toml` 0.1.0 | `alphagenome` 0.11.0 (pinned `<2`, see below) |
| `alphagenome-atlas-website-links` | google-deepmind/science-skills | `5a50f0be962b…` | publisher google; `pyproject.toml` 0.1.0, `requires-python >=3.14` | Python 3.13, `alphagenome` 2.0.0 |

**Running the pyproject-based DeepMind scripts.** `alphagenome_atlas_avi.py` and
`alphagenome_atlas_links.py` have no PEP 723 header; their skill's `pyproject.toml` holds the
dependencies. `uv run` from the skill folder would write `.venv/` and `uv.lock` into the installed
skill, and `uv run --with-requirements pyproject.toml` is refused by uv, so the fetches pass the same
list explicitly, with exact versions (`deps.AVI_DEPS`, `deps.LINKS_DEPS`), under
`uv run --no-project --python 3.13`. Versions resolved on 2026-11-15: `alphagenome` 2.0.0, `polars`
1.31.0 (newer releases refuse to start on Windows ARM64 under x86-64 emulation: `unknown feature flag:
'sse3'`), `pandas` 2.3.3. `pysam` carries the marker `sys_platform != 'win32'`, as in the DeepMind
pyproject: it has no Windows wheel and the AVI script imports it only for VCF output. Scripts with
their own PEP 723 header (`atlas_query.py`, `search_ols.py`, `get_term.py`, `ensembl_api.py`) get
their packages pinned with `--with` (`alphagenome==2.0.0`, `polite-http==0.1.4`).

**`alphagenome` 2.0.0.** An unpinned `alphagenome` now resolves to 2.0.0, which no longer exposes
`alphagenome.models.dna_client` (`ImportError: cannot import name 'dna_client'`). The Atlas
client (`alphagenome.atlas`) and `alphagenome.data.genome` still work, so `alphagenome_atlas_avi.py`,
`atlas_query.py` and `alphagenome_atlas_links.py` run; `interpret_splicing.py` (and every model-API
script) needs `uv run --with "alphagenome<2" …` (resolved 0.11.0).

**Keys.** The key lives in the workspace's `.env` or in `~/.env` (`deps.env_file()` picks the first that sets
it). The DeepMind scripts load `~/.env` themselves; the fetch also passes the chosen file with `--env-file`. The
K-Dense scripts read only the
environment, so the fetch runs them as `uv run --env-file <home>/.env …` with the path expanded in
Python and written with `/` (uv drops the backslashes of a Windows path given to `--env-file`). uv
loads the file into the process environment; no shell sources it, and the key never appears on a
command line or in the manifest.

**polite-http.** The OLS and Ensembl scripts pace requests through a lock file per host holding the
monotonic time of the next free slot. A file left by a process with another clock can hold a slot
hours ahead, and the script then sleeps without a message. The fetches set `POLITE_HTTP_LOCK_DIR` to
`~/.cache/alphagenome-gene-skills/polite-http` and reset a slot more than ten minutes ahead.

## `gtf.json` — gene structure (Workflow step 2)

```bash
alphagenome_atlas_avi.py gtf --gene G --all_transcripts --exons --cds --utr --format json   # fetch gtf
```

`--format json` prints the JSON to stdout and ignores `--output`; the script also prints one
`Loading GTF from …` line to stdout first, which `fetch gtf` drops in Python (nothing else is changed).
The first run downloads the 318 MB GENCODE v46 Feather to `~/.cache/alphagenome/`.

List of `GtfTranscript`:

| Field | Type | Note |
|---|---|---|
| `gene_name`, `gene_id`, `transcript_id`, `transcript_type`, `chromosome`, `strand` | str | IDs are versioned (GENCODE v46) |
| `start`, `end` | int | see "Coordinates" below |
| `is_mane_select` | bool | |
| `num_exons` | int | |
| `exons` | list of `GtfExon` | `exon_number, chromosome, start, end, width, strand, acceptor, donor` (int/str), `overlaps_variant` (bool) |
| `junctions` | list of `GtfJunction` | always written, even without `--junctions` |
| `cds`, `utr` | list of `GtfSegment` | `start, end, width` (int); 5' and 3' UTR in one list; GTF `CDS` excludes the stop codon |

**Coordinates and `normalize_gtf`.** Convention of the raw file (kept as written in `raw/` and in
the dataclasses):

- `start` of transcripts, exons, CDS and UTR segments is the GENCODE v46 Feather `Start` column as
  stored, which is **0-based**; `end` is the **1-based** last base. So `start` is one below the
  1-based first base and `width = end - start + 1` is one more than the feature length. The skill's
  SKILL.md calls the output "1-based closed"; for starts it is not.
- `acceptor`/`donor` are the 5'/3' exon boundaries: the script writes acceptor = `start`,
  donor = `end` on `+`, and acceptor = `end`, donor = `start` on `−`, so the one equal to `start`
  carries the same 0-based shift.
- Junctions are **not** shifted: `junction_start` is the last base of the upstream exon and
  `junction_end` the last base of the intron (both 1-based), i.e. the 0-based half-open intron,
  which is how Atlas `SPLICE_JUNCTIONS` reports `junction_start`/`junction_end`. The 1-based intron
  is `junction_start + 1 .. junction_end`; `intron_length = junction_end - junction_start - 1` as
  written is one less than its length.

Evidence (2026-11-12): the AlphaGenome package documents the Feather `Start` as 0-based
(`alphagenome.data.gene_annotation.extract_tss`: "Start == End, 0-based"). In both bundled
references every MANE exon `start` in `raw/gtf.json` is exactly one below the Ensembl REST
(release 116) first base of the same exon, every `end` equals Ensembl's, and the expression
reference's TSS and exon-1 length match its reference case only after the shift (+1 / −1 bp).
The splicing live capture shows the reference transcript's canonical intron-3 junction with the
same `junction_start`/`junction_end` pair in `gtf_junctions.json` and in Atlas
`SPLICE_JUNCTIONS` rows (`scores_cohort.tsv`). Concrete coordinates stay in
`assets/reference/<gene>/` (no gene values in this document).

Use: stages call `models.normalize_gtf(transcript)` once on each record read from `raw/`. It
returns a copy with 1-based closed starts for the transcript, exons, CDS and UTR segments,
`width = end - start + 1`, the shifted `acceptor`/`donor`, junction coordinates unchanged and
`intron_length = junction_end - junction_start`; it raises `ValueError` if an exon's
acceptor/donor do not follow the convention above. After it, the TSS is `start` on `+` and `end`
on `−`. Only `+`-strand genes were recorded so far; the `−` branch follows the script's code.

## Exploration-only files (`OUT/exploration/atlas/`, step 0)

- `atlas_tracks.tsv` — `atlas_query.py tracks -o …` (K-Dense, no `--scorer`): the whole Atlas track
  catalogue, one row per scorer × track: `scorer`, `track_index`, `name`, `strand`, `ontology_curie`,
  `biosample_name`, `biosample_type`, `assay_title`, `gtex_tissue`, `data_source` (columns present only
  when the server sends them). It adds the source and assay that the AVI catalogue lacks.
- `scores_RNA_SEQ_ACTIVE_baseline.tsv` — `atlas_query.py scores --variant V --scorers RNA_SEQ_ACTIVE
  --gene G -o …` for one SNV in the middle of exon 1 (REF from the sequence): one row per RNA-seq track
  with the gene's predicted activity (`raw_score`); same columns as `scores_{S}.tsv`.

Neither is copied into a study; the files of `OUT/exploration/raw/` are, with their manifest steps.

## `regulatory.json` — Ensembl Regulatory Build (step 3)

```bash
uv run skills/expression-analysis/scripts/regulatory_features.py overlap --region R --output raw/regulatory.json
```

Own client (`EnsemblRegulatoryClient`, 15 requests/s with a file lock at
`~/.cache/alphagenome-gene-skills/ensembl.lock`); `GET /overlap/region/human/{R}?feature=regulatory`,
written as served. The success line reports the Ensembl release (`/info/data`) for the manifest.

List of `RegulatoryFeature`: `id` (str), `feature_type` (str, `regulatory`), `start`, `end` (int,
1-based), `strand` (int: 0, 1, −1), `description` (str: promoter, enhancer, CTCF_binding_site, …),
`extended_start`, `extended_end` (int), `seq_region_name` (str, no `chr`), `source` (str).

## `sequence_{chr}_{start}_{end}.txt` — reference sequence (step 3)

```bash
uv run .claude/skills/ensembl-database/scripts/ensembl_api.py get-sequence R --species human --output raw/sequence_{chr}_{start}_{end}.txt
```

`get-sequence` writes plain text (only the bases, no FASTA header), not JSON, so the region goes
into the file name. `EnsemblSequence`: `chromosome`, `start`, `end` (from the name, 1-based
closed), `sequence` (str); `models.py` checks the length equals `end - start + 1`.

## `ols_search_{slug}.json` — UBERON/CL candidates (step 5)

```bash
uv run .claude/skills/embl-ebi-ols/scripts/search_ols.py --query "TISSUE" --ontology uberon --defining --rows 5 --output raw/ols_search_{slug}.json
```

`OlsSearchResult`: `status` (str, `success`), `total_found`, `results_count` (int),
`pagination` {`start`, `rows` (int), `has_more` (bool)}, `terms`: list of `OlsSearchHit`
{`iri`, `label`, `ontology_name`, `ontology_prefix`, `obo_id`, `short_form`, `type` (str),
`description`, `exact_synonyms` (list of str), `is_defining_ontology` (bool; served `false` even
with `--defining`)}. A plain tissue word may not return the intended term in the first rows
(e.g. a muscle query returned only sub-structures); search again with the full label and keep
both files. On failure the script writes `{"status": "error", "message": …}`, which `models.py`
reports on `status`.

## `ols_term_{PREFIX}_{digits}.json` — children of an unserved CURIE (step 6, optional)

```bash
uv run .claude/skills/embl-ebi-ols/scripts/get_term.py --obo_id C --relations hierarchicalChildren,hierarchicalDescendants --output raw/ols_term_{PREFIX}_{digits}.json
```

`OlsTermResult`: `status` (str), `term`: `OlsTerm` {`iri`, `label`, `obo_id`, `ontology_name`,
`ontology_prefix`, `short_form` (str), `description`, `synonyms`, `in_subset` (list of str),
`annotation` (object, free form), `is_defining_ontology`, `is_obsolete`, `has_children`,
`is_root` (bool)}, plus one list of {`iri`, `label`, `obo_id`} per requested relation
(`hierarchicalChildren`, …; optional keys, written only when requested) and `graph` (object,
only with `--relations graph`). `hierarchicalChildren` is one level; a tissue such as a cortical lobe is several levels below its
parent, so step 6 requests `hierarchicalDescendants` too (`cohorts` offers them as replacements).

## `tracks_{S}_{PREFIX}_{digits}.parquet` — Atlas track catalogue (step 7)

```bash
uv run --no-project --python 3.13 <AVI_DEPS pinned> .claude/skills/alphagenome-variant-impact-score/scripts/alphagenome_atlas_avi.py \
  metadata --tracks --ontology C --scorer S --output raw/tracks_{S}_{PREFIX}_{digits}.parquet
```

One command per CURIE × scorer. When nothing matches, the script prints
`No tracks matched the specified metadata filter criteria.` and writes no file; the manifest
records that command with `"outputs": []`, which is how a zero-track pair is told apart from a
pair never queried. `AtlasTrack` rows: `scorer`, `track_name`, `biosample_name`,
`biosample_type`, `ontology_curie`, `strand` (str), `track_index` (int); optional
`transcription_factor`, `histone_mark`, `assay` (only for scorers whose catalogue has them).
`models.py` checks that `scorer` and `ontology_curie` match the file name.

### `tracks_{S}_catalogue.parquet` — whole catalogue per scorer (step 7, once per study scorer)

```bash
uv run --no-project --python 3.13 <AVI_DEPS pinned> .claude/skills/alphagenome-variant-impact-score/scripts/alphagenome_atlas_avi.py \
  metadata --tracks --scorer S --output raw/tracks_{S}_catalogue.parquet
```

Same `AtlasTrack` rows as above for every track of the scorer (`--output` writes all rows; `--top_n`
only limits stdout). `models.py` checks the columns but not a CURIE, since the file spans all of them.
`cohorts` reads it to offer, for an unserved CURIE, only the OLS descendants that at least one study
scorer serves; without it the question lists every descendant and says so. A manifest step records it.

## `avi_region.parquet` — saturation scan with AVI (step 9)

```bash
uv run --no-project --python 3.13 <AVI_DEPS pinned> .claude/skills/alphagenome-variant-impact-score/scripts/alphagenome_atlas_avi.py \
  region --region W --top_k 20 --min_phred 0 --max_window_size M --include_features --include_track_info --output raw/avi_region.parquet
```

All 3 × N SNVs, ranked. `AviRecord`: `rank`, `position` (int), `variant`, `chromosome`, `ref`,
`alt`, `top_modality` (str), `avi_phred`, `avi_raw`, `avi_quantile`, `top_percentile`,
`top_feature_importance` (float), and the 18 `fi_<FEATURE>` columns (float; all required). All of
these must hold values (null, empty or NaN fails); none of the recorded outputs leaves one null.
`region` never resolves track provenance, so `--include_track_info` adds no `track_*` column here
(they appear only in `query` output, `avi_query.json`).

## `scores_{S}.tsv` — per-tissue track scores (step 10)

```bash
uv run --env-file <home>/.env --with alphagenome==2.0.0 .claude/skills/alphagenome/scripts/atlas_query.py scores \
  --interval W --max-window M --scorers S --ontology C1 C2 … [--gene G] -o raw/scores_{S}.tsv
```

`--gene G` only for gene-centric scorers (`RNA_SEQ`, `RNA_SEQ_ACTIVE`, `POLYADENYLATION`,
`SPLICE_*`), as the original pipeline did. One row per variant × track (× gene, × junction).
`TrackScoreRow`: `variant`, `scorer` (str, required); `track_index` (int) and `raw_score` (float)
required unless the file has an `error` column, and on every row with an empty `error` cell
`scorer`, `track_index` and `raw_score` must hold values; optional, written only when the server returns
them: `gene_id`, `gene_name`, `gene_strand` (gene-centric), `junction_start`, `junction_end` (int,
`SPLICE_JUNCTIONS`), `track_name`, `track_strand`, `ontology_curie`, `biosample_name`,
`biosample_type`, `assay_title`, `gtex_tissue`, `transcription_factor`, `histone_mark`,
`data_source` (str), `quantile_score` (float), `error` (str, `--variant/--input` runs only).
Empty cells are missing values. Columns seen on 2026-11-12: `RNA_SEQ` has the gene columns and
`gtex_tissue`; `CAGE`, `DNASE`, `ATAC` have neither.

## `interpret_{variant}.txt` — splicing interpretation of the top variants (step 13)

```bash
uv run --no-project --python 3.13 <SVA_DEPS pinned, alphagenome==0.11.0> .claude/skills/alphagenome-single-variant-analysis/scripts/interpret_splicing.py \
  --chrom CHR --pos POS --ref REF --alt ALT --ontology_id C > raw/interpret_{CHR}_{POS}_{REF}_{ALT}.txt
```

**Known issue:** an unpinned `alphagenome` resolves to 2.0.0, which removed
`alphagenome.models.dna_client`; without `--with "alphagenome<2"` the script stops with
`ImportError: cannot import name 'dna_client' from 'alphagenome.models'` (0.11.0 was used).

The script has no `--output` and prints text (not JSON), so stdout is saved. `models.py` parses
it into `SpliceInterpretation` (`variant`, `ontology_id`, `splice_site_usage`
[`position`, `offset`, `REF`, `ALT`, `Delta`, `flag`], `junctions_overlapping_variant`
[`start`, `end`, `strand`, `REF`, `ALT`, `Delta`, `flag`], `potential_cryptic_junctions`
[`start`, `end`, `REF`, `ALT`], `top_changed_junctions` [`index`, `start`, `end`, `score_REF`,
`score_ALT`, `delta`], `mean_rna_coverage_ref`, `mean_rna_coverage_alt`, `log2_fold_change`,
`messages`); an unrecognised line fails validation with its line number. Format recorded from the
splicing capture; the overlapping-junction list can run to thousands of lines.

## `atlas_links.md` — Atlas deep links (step 13)

```bash
uv run --no-project --python 3.13 <LINKS_DEPS pinned> .claude/skills/alphagenome-atlas-website-links/scripts/alphagenome_atlas_links.py \
  table --input FILE.json > raw/atlas_links.md
```

`table --input` reads a JSON list of records with `variant` and, when present, `avi_phred`,
`avi_raw`, `top_modality`; it does not read TSV. A Markdown table on stdout (no `--output`).
No dataclass: the report embeds the text as is.
