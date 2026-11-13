# Workflow reference: expression-analysis

The complete commands of the 14 Workflow steps of `SKILL.md`. Formats and fields of every raw
file are in `dependency_outputs.md`; how to read the results is in `interpretation.md`.

## Contents

- [Conventions](#conventions): variables, recording a step, researcher decisions
- [Step 1. Prerequisites and study](#step-1-prerequisites-and-study)
- [Step 2. Gene structure](#step-2-gene-structure)
- [Step 3. Regulatory Build and sequence](#step-3-regulatory-build-and-sequence)
- [Step 4. Locus and window](#step-4-locus-and-window)
- [Step 5. UBERON/CL terms](#step-5-uberoncl-terms)
- [Step 6. Descendants (optional)](#step-6-descendants-optional)
- [Step 7. Track catalogue](#step-7-track-catalogue)
- [Step 8. Cohorts](#step-8-cohorts)
- [Step 9. Saturation and AVI](#step-9-saturation-and-avi)
- [Step 10. Scores per tissue](#step-10-scores-per-tissue)
- [Step 11. Variants and Tier 1](#step-11-variants-and-tier-1)
- [Step 12. Contrast (Tier 2)](#step-12-contrast-tier-2)
- [Step 13. Interpretation and report](#step-13-interpretation-and-report)
- [Step 14. Workbench](#step-14-workbench)
- [When a command fails](#when-a-command-fails)

## Conventions

Run the commands from the directory that holds the skills (the one with `.claude/skills/`). The
loops below work in bash and zsh. Set these once per shell:

```bash
SKILL=.claude/skills/expression-analysis   # path to this skill directory
DEPS=.claude/skills                        # where `expression.py check` found the dependency skills
D=studies/mygene-expression                # study directory chosen by the researcher, outside the skill
G=MYGENE                                   # HGNC symbol, as in study.yaml gene_symbol
AVI_DEPS="alphagenome>=0.9.0,absl-py,anndata,numpy>=2.0,pandas<3.0,polars,pyarrow,pysam,python-dotenv,tabulate"
LINKS_DEPS="alphagenome>=0.9.0,numpy>=2.0,pandas<3.0,pyarrow,python-dotenv"
SVA_DEPS="alphagenome<2,absl-py,python-dotenv,matplotlib,numpy,pandas,pyarrow,scipy,seaborn"
```

Values in `<…>` and CURIEs such as `UBERON:AAAAAAA` are placeholders for the study's own values;
a `--decision-key` placeholder stands for one `--decision-key KEY` per question answered in that
step (none when the stage asked nothing).

Why the dependency lists: `alphagenome_atlas_avi.py` and `alphagenome_atlas_links.py` have no
PEP 723 header, and their skills' `pyproject.toml` asks for Python ≥ 3.14, where `h5py` has no
wheel. Running them from the skill folder would also write `.venv/` into the installed skill. So
they run with `--no-project --python 3.13` and the list passed explicitly. `SVA_DEPS` is the
`alphagenome-single-variant-analysis` list with `alphagenome<2`, because `alphagenome` 2.0.0
removed `dna_client`, which its model-API scripts import.

**Recording a step.** Every step ends with `manifest.py add`: the report lists the skills used from
`manifest.json`, the dependencies ask to be named in any output, and the manifest is what lets
someone else re-run the study. Record the command exactly as it ran: wrap `--command` in double
quotes so `$D`, `$G` and the paths expand, and escape what must stay literal (`\"` for inner double
quotes, `\$` for a `$` the command itself contains). Outputs are paths relative to `$D` and must
exist; a command that legitimately wrote nothing is recorded without `--output`. The first `add`
creates the manifest and is the only one that takes `--origin`/`--origin-note`.

**Researcher decisions.** A stage that meets an open choice writes it to
`processed/pending_decisions.json` and exits 1 with `N decision(s) needed: run study.py pending.`

```bash
uv run $SKILL/scripts/study.py pending --study-dir $D
# show each question, its options and the recommendation to the researcher; then, with THEIR answer:
uv run $SKILL/scripts/study.py decide --study-dir $D --key 'locus.exon1_overlap:FEATURE' --value keep_full_unit \
  --step 4 --note "Question as shown; options; the researcher's reason"
```

`--value` is parsed as YAML, so a list is written `'[UBERON:0000955, UBERON:0002107]'`. Quote keys
(they contain `:`). Then run the stage again; it reads `decisions.yaml` and goes on. A missing
answer that a stage needs stops it with `Missing decision '{key}'. Run Workflow step {N} …`.
`decide` copies the question's `options` and `recommended` from `pending_decisions.json` into the
`decisions.yaml` entry, so the record shows what was offered next to what was chosen. For a
question you ask yourself (no pending entry, e.g. `cohort.target` in step 5), pass what you offered
with `--option` (once per option) and `--recommended`. If the researcher delegated the choice
("take the recommended options"), record the recommendation and write that in the note.

A key must look like `<area>.<topic>[:<subject>]` (as `study.py pending` prints it): `decide`
refuses anything else with exit 1 and records nothing, and warns about a key that is neither pending
nor decided before. Several keys with the same answer go in one call, one `--key` each. Write the
keys out as separate words; never loop over an unquoted variable that holds several keys (`for K in
$KEYS` gets one newline-joined word in zsh). To take the keys from a file, read it line by line:

```bash
uv run $SKILL/scripts/study.py decide --study-dir $D --key 'variants.ref_mismatch:<POS_A>' \
  --key 'variants.ref_mismatch:<POS_B>' --value exclude_position --step 11 --note "Both positions excluded by the researcher"
while IFS= read -r K; do                                   # keys.txt: one key per line (bash and zsh)
  uv run $SKILL/scripts/study.py decide --study-dir $D --key "$K" --value exclude_position --step 11 --note "…"
done < keys.txt
```

## Step 1. Prerequisites and study

```bash
uv run $SKILL/scripts/expression.py check
```

- `All prerequisites satisfied` (exit 0): go on.
- Exit 1: the list names each missing item. For a **skill**, show the researcher the skill, its
  origin and the printed install command, and ask about one skill at a time. Run the install only
  for the skill they approved (or let them run it), then run `check` again. For a **licence
  notice**, open that dependency's `SKILL.md` Prerequisites and go through the notice with the
  researcher; create the `.licenses/` file only after their explicit OK. For a **key**, follow the
  `credentials` skill. For **uv**, the `uv` skill.

Create the study and fill it in with the researcher:

```bash
uv run $SKILL/scripts/study.py init --study-dir $D --gene $G
```

Ask, and write the answers into `$D/study.yaml`: the gene (HGNC symbol); the target and control
tissues or cell types in words and a short label for each group (`cohorts.*.label`; the CURIEs
come in step 5); the scorers (default `RNA_SEQ, CAGE, DNASE, ATAC`; `RNA_SEQ` is required because
it is the expression contrast); the window policy (`regulatory_unit`, built from the Regulatory
Build around the TSS, or `explicit` with `explicit_interval: chrN:START-END`);
`regulatory_search_upstream_bp` (default 2000) and `report.top_n` (default 20).

```bash
uv run $SKILL/scripts/study.py validate --study-dir $D
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 1-study --skill expression-analysis \
  --command "uv run $SKILL/scripts/study.py init --study-dir $D --gene $G" --output study.yaml \
  --origin live --origin-note "Requested by <researcher> on <date>: <question in one line>"
```

## Step 2. Gene structure

```bash
uv run --no-project --python 3.13 --with "$AVI_DEPS" $DEPS/alphagenome-variant-impact-score/scripts/alphagenome_atlas_avi.py \
  gtf --gene $G --all_transcripts --exons --cds --utr --format json | sed -n '/^\[$/,$p' > $D/raw/gtf.json
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 2-gtf --skill alphagenome-variant-impact-score \
  --command "uv run --no-project --python 3.13 --with \"$AVI_DEPS\" $DEPS/alphagenome-variant-impact-score/scripts/alphagenome_atlas_avi.py gtf --gene $G --all_transcripts --exons --cds --utr --format json | sed -n '/^\[\$/,\$p' > $D/raw/gtf.json" \
  --output raw/gtf.json --source-version gencode=v46 --source-version python=3.13
```

- GENCODE v46 is the annotation of the Atlas scores, which is why exons come from here and not
  from Ensembl REST. The first run downloads the 318 MB GENCODE Feather to `~/.cache/alphagenome/`.
- `--format json` prints to stdout and ignores `--output`; the `sed` filter drops the
  `Loading GTF from …` line printed before the JSON.
- Starts are 0-based as written (ends 1-based). Keep the file unedited: every stage applies
  `models.normalize_gtf`.

## Step 3. Regulatory Build and sequence

```bash
uv run $SKILL/scripts/expression.py regions --study-dir $D          # JSON: regulatory_search (and window after step 4)
R=$(uv run $SKILL/scripts/expression.py regions --study-dir $D --only regulatory_search)
uv run $SKILL/scripts/regulatory_features.py overlap --region $R --output $D/raw/regulatory.json
uv run $DEPS/ensembl-database/scripts/ensembl_api.py get-sequence $R --species human \
  --output $D/raw/sequence_$(echo $R | tr ':-' '__').txt
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 3-regulatory --skill expression-analysis/regulatory_features.py \
  --command "uv run $SKILL/scripts/regulatory_features.py overlap --region $R --output $D/raw/regulatory.json" \
  --output raw/regulatory.json --source-version ensembl=<release printed by overlap> \
  --source-version region_rule="gene span of raw/gtf.json (normalised) +/- regulatory_search_upstream_bp"
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 3-sequence --skill ensembl-database \
  --command "uv run $DEPS/ensembl-database/scripts/ensembl_api.py get-sequence $R --species human --output $D/raw/sequence_$(echo $R | tr ':-' '__').txt" \
  --output raw/sequence_$(echo $R | tr ':-' '__').txt --source-version ensembl=<release> --source-version assembly=GRCh38
```

- `regions` computes the gene span over every transcript of `raw/gtf.json` (starts normalised to
  1-based) ± `window.regulatory_search_upstream_bp` on both sides, so the upstream flank is covered
  on either strand and nobody adds coordinates by hand.
- The features decide the window (step 4); the sequence gives REF for every SNV and checks the
  start codon. The sequence file name carries the region (`get-sequence` writes bare bases).

## Step 4. Locus and window

```bash
uv run $SKILL/scripts/expression.py locus --study-dir $D
```

Questions (answer with the researcher, `--step 4`, then run `locus` again):

| Key | Asks | Options |
|---|---|---|
| `locus.transcript` | 0 or several MANE Select transcripts | a transcript ID |
| `locus.exon1_overlap:{feature}` | a chained feature runs into exon 1 | `keep_full_unit`, `abort` |
| `locus.multiple_tss_features` | several features contain or abut the TSS | `union`, `abort` |
| `locus.feature_overlap:{a}+{b}`, `locus.feature_gap:{a}+{b}` | how to chain neighbouring features | `union` / `bridge_gap`, `abort` |
| `locus.no_feature_chain` | no feature contains or abuts the TSS | `set_explicit_interval` (then `window.policy: explicit` in `study.yaml`), `abort` |

The success line gives the transcript, strand, TSS, window and whether the start codon reads ATG
(reverse complement on the − strand). Open `figures/fig1_genomic_architecture.html` with the
researcher before going on: the window is the set of variants the whole study is about.

```bash
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 4-locus --skill expression-analysis \
  --command "uv run $SKILL/scripts/expression.py locus --study-dir $D" \
  --output processed/locus.json --output processed/transcripts_in_window.tsv \
  --output figures/fig1_genomic_architecture.html --decision-key '<each locus.* key answered>'
uv run $SKILL/scripts/expression.py regions --study-dir $D     # now also prints window and max_window
```

## Step 5. UBERON/CL terms

One search per tissue or cell type the researcher named (`--ontology cl` for cell types); the
slug is the query in lower case with `_` for spaces:

```bash
uv run $DEPS/embl-ebi-ols/scripts/search_ols.py --query "cerebral cortex" --ontology uberon --defining --rows 5 \
  --output $D/raw/ols_search_cerebral_cortex.json
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 5-ols-search --skill embl-ebi-ols \
  --command "uv run $DEPS/embl-ebi-ols/scripts/search_ols.py --query \"cerebral cortex\" --ontology uberon --defining --rows 5 --output $D/raw/ols_search_cerebral_cortex.json" \
  --output raw/ols_search_cerebral_cortex.json --source-version ols=OLS4
```

Show the researcher each hit's `label`, `obo_id` and `description`, and say which hit you
recommend for each tissue they named (the one whose label and definition match their words). A
plain word may not return the intended term in the first rows (a muscle query returned only
sub-structures): search again with the full label and keep both files. The researcher picks the
CURIEs; record them with the hits you offered (`--option` once per hit, `--recommended` the one you
suggested) and copy them into `study.yaml` (`cohorts.target.curies`, `cohorts.control.curies`),
then validate:

```bash
uv run $SKILL/scripts/study.py decide --study-dir $D --key cohort.target --value '[UBERON:AAAAAAA, UBERON:BBBBBBB]' \
  --option 'UBERON:AAAAAAA (<label>)' --option 'UBERON:BBBBBBB (<label>)' --option 'UBERON:DDDDDDD (<label>)' \
  --recommended 'UBERON:AAAAAAA (<label>)' --step 5 --note "Target cohort chosen by the researcher from the OLS hits; why"
uv run $SKILL/scripts/study.py decide --study-dir $D --key cohort.control --value '[UBERON:CCCCCCC]' \
  --option 'UBERON:CCCCCCC (<label>)' --option 'UBERON:EEEEEEE (<label>)' --recommended 'UBERON:CCCCCCC (<label>)' \
  --step 5 --note "Control cohort chosen by the researcher; why"
uv run $SKILL/scripts/study.py validate --study-dir $D
```

Substituting a term yourself (a cell type for a tissue the Atlas does not serve, say) is the
researcher's choice, not yours: step 8 asks it as `cohort.no_rna_seq:{role}` or
`cohort.unserved:{curie}` with the served options.

A CURIE cannot be in both cohorts (a tissue is not its own control); `cohorts` stops if it is.

## Step 6. Descendants (optional)

For each chosen CURIE that step 7 finds without tracks for every scorer (all its pairs printed
`No tracks matched …`) or without any `RNA_SEQ` track (its `RNA_SEQ` pair printed `No tracks
matched …`):

```bash
C=UBERON:AAAAAAA
uv run $DEPS/embl-ebi-ols/scripts/get_term.py --obo_id $C --relations hierarchicalChildren,hierarchicalDescendants \
  --output $D/raw/ols_term_${C/:/_}.json
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 6-ols-children --skill embl-ebi-ols \
  --command "uv run $DEPS/embl-ebi-ols/scripts/get_term.py --obo_id $C --relations hierarchicalChildren,hierarchicalDescendants --output $D/raw/ols_term_${C/:/_}.json" \
  --output raw/ols_term_${C/:/_}.json --source-version ols=OLS4
```

Why both relations: `hierarchicalChildren` is one level, and a served tissue can sit several
levels below the one asked for. Step 8 offers only the descendants that the step-7 catalogue
serves (for `cohort.no_rna_seq`, those with `RNA_SEQ` tracks, most first). Steps 6 and 7 can run
in either order; step 6 is only needed for CURIEs without tracks or without `RNA_SEQ` tracks.

## Step 7. Track catalogue

Every CURIE of `cohort.target` and `cohort.control` (and any replacement chosen in step 8) ×
every scorer of `study.yaml`. A pair without tracks prints `No tracks matched the specified
metadata filter criteria.` and writes no file; it is still recorded (without `--output`), which is
how a zero-track pair is told apart from a pair never queried.

```bash
for C in UBERON:AAAAAAA UBERON:BBBBBBB UBERON:CCCCCCC; do      # every CURIE of both cohorts
  for S in RNA_SEQ CAGE DNASE ATAC; do                       # the scorers of study.yaml
    f=tracks_${S}_${C/:/_}.parquet
    cmd="uv run --no-project --python 3.13 --with \"$AVI_DEPS\" $DEPS/alphagenome-variant-impact-score/scripts/alphagenome_atlas_avi.py metadata --tracks --ontology $C --scorer $S --output $D/raw/$f"
    uv run --no-project --python 3.13 --with "$AVI_DEPS" $DEPS/alphagenome-variant-impact-score/scripts/alphagenome_atlas_avi.py \
      metadata --tracks --ontology $C --scorer $S --output $D/raw/$f
    if [ -f "$D/raw/$f" ]; then
      uv run $SKILL/scripts/manifest.py add --study-dir $D --step 7-tracks --skill alphagenome-variant-impact-score \
        --command "$cmd" --output raw/$f --source-version atlas_pull=$(date +%F)
    else
      uv run $SKILL/scripts/manifest.py add --study-dir $D --step 7-tracks --skill alphagenome-variant-impact-score \
        --command "$cmd" --source-version atlas_pull=$(date +%F) --source-version result=no_tracks
    fi
  done
done
```

Then the whole catalogue of each scorer, once per study. `cohorts` uses it to offer, for a CURIE
without tracks, only descendants that some study scorer serves:

```bash
for S in RNA_SEQ CAGE DNASE ATAC; do                         # the scorers of study.yaml
  uv run --no-project --python 3.13 --with "$AVI_DEPS" $DEPS/alphagenome-variant-impact-score/scripts/alphagenome_atlas_avi.py \
    metadata --tracks --scorer $S --output $D/raw/tracks_${S}_catalogue.parquet
  uv run $SKILL/scripts/manifest.py add --study-dir $D --step 7-tracks-catalogue --skill alphagenome-variant-impact-score \
    --command "uv run --no-project --python 3.13 --with \"$AVI_DEPS\" $DEPS/alphagenome-variant-impact-score/scripts/alphagenome_atlas_avi.py metadata --tracks --scorer $S --output $D/raw/tracks_${S}_catalogue.parquet" \
    --output raw/tracks_${S}_catalogue.parquet --source-version atlas_pull=$(date +%F)
done
```

## Step 8. Cohorts

```bash
uv run $SKILL/scripts/expression.py cohorts --study-dir $D
```

| Key | Asks | Options |
|---|---|---|
| `cohort.unserved:{curie}` | the CURIE has no track for any study scorer | `drop_curie`, `replace:{descendant CURIE}` (served descendants listed), `abort` |
| `cohort.scorer_without_cohort:{scorer}` | a scorer has no track in the target or the control cohort | `drop_scorer`, `abort` |
| `cohort.no_rna_seq:{role}` | the target or control cohort keeps no `RNA_SEQ` track (the expression contrast) | `add:{descendant CURIE}` (OLS descendants with `RNA_SEQ` tracks in the catalogue, most tracks first and recommended), `abort` |

After `replace:{curie}` or `add:{curie}`, pull that CURIE's tracks (step 7, all scorers) and run
`cohorts` again. `add:` keeps the cohort's CURIEs (they still serve the other scorers) and adds the
descendant. A cohort without any `RNA_SEQ` track cannot be fixed by dropping a scorer; when no
descendant is listed, step 6 has not run for its CURIEs, or none is served: choose other CURIEs
(step 5) or `abort`. Show `processed/tissue_cohorts.md` (tissues, labels, track counts) to the
researcher.

```bash
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 8-cohorts --skill expression-analysis \
  --command "uv run $SKILL/scripts/expression.py cohorts --study-dir $D" \
  --output processed/track_registry.tsv --output processed/ontology_resolution.json \
  --output processed/tissue_cohorts.md --decision-key '<each cohort.* key answered>'
```

## Step 9. Saturation and AVI

```bash
W=$(uv run $SKILL/scripts/expression.py regions --study-dir $D --only window)
M=$(uv run $SKILL/scripts/study.py window --study-dir $D --print-max-window)
uv run --no-project --python 3.13 --with "$AVI_DEPS" $DEPS/alphagenome-variant-impact-score/scripts/alphagenome_atlas_avi.py \
  region --region $W --top_k 20 --min_phred 0 --max_window_size $M --include_features --include_track_info \
  --output $D/raw/avi_region.parquet
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 9-avi-region --skill alphagenome-variant-impact-score \
  --command "uv run --no-project --python 3.13 --with \"$AVI_DEPS\" $DEPS/alphagenome-variant-impact-score/scripts/alphagenome_atlas_avi.py region --region $W --top_k 20 --min_phred 0 --max_window_size $M --include_features --include_track_info --output $D/raw/avi_region.parquet" \
  --output raw/avi_region.parquet --source-version atlas_pull=$(date +%F)
```

`--min_phred 0` keeps all 3 × N SNVs (the file is the universe of variants, not a top list;
`--top_k` only limits what is printed). `--max_window_size` is `max(1000, window length)`: the
script refuses windows above 1000 bp by default.

## Step 10. Scores per tissue

One call per kept scorer, with every CURIE of both cohorts (including replacements; a CURIE
without tracks simply returns no rows) and `--gene $G` only for the gene-centric `RNA_SEQ`:

```bash
CURIES="UBERON:AAAAAAA UBERON:BBBBBBB UBERON:CCCCCCC"
uv run --env-file ~/.env $DEPS/alphagenome/scripts/atlas_query.py scores --interval $W --max-window $M --scorers RNA_SEQ \
  --ontology $(echo $CURIES) --gene $G -o $D/raw/scores_RNA_SEQ.tsv
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 10-scores --skill alphagenome \
  --command "uv run --env-file ~/.env $DEPS/alphagenome/scripts/atlas_query.py scores --interval $W --max-window $M --scorers RNA_SEQ --ontology $CURIES --gene $G -o $D/raw/scores_RNA_SEQ.tsv" \
  --output raw/scores_RNA_SEQ.tsv --source-version atlas_pull=$(date +%F)
for S in CAGE DNASE ATAC; do                                 # the other scorers kept in step 8
  uv run --env-file ~/.env $DEPS/alphagenome/scripts/atlas_query.py scores --interval $W --max-window $M --scorers $S \
    --ontology $(echo $CURIES) -o $D/raw/scores_$S.tsv
  uv run $SKILL/scripts/manifest.py add --study-dir $D --step 10-scores --skill alphagenome \
    --command "uv run --env-file ~/.env $DEPS/alphagenome/scripts/atlas_query.py scores --interval $W --max-window $M --scorers $S --ontology $CURIES -o $D/raw/scores_$S.tsv" \
    --output raw/scores_$S.tsv --source-version atlas_pull=$(date +%F)
done
```

`$(echo $CURIES)` splits the list into separate arguments in bash and zsh alike. The K-Dense
script reads the key only from the environment, hence `--env-file ~/.env`. Capture steps 9 and 10
together: scores of another window or pull do not line up with the AVI table.

## Step 11. Variants and Tier 1

```bash
uv run $SKILL/scripts/expression.py variants --study-dir $D
uv run $SKILL/scripts/expression.py master --study-dir $D
```

`variants` asks `variants.ref_mismatch:{pos}` (the Atlas REF differs from the Ensembl sequence) and
`variants.missing_in_atlas:{pos}` (SNVs without an AVI row): `exclude_position` or `abort`, with
the researcher. `master` exits 1 when a QC gate fails and lists the variants in
`processed/qc_report.json` (`failed_variants`); the usual cause is raw files from different windows
or pulls. Find the cause rather than continuing, because Tier 2 would inherit it.

```bash
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 11-variants --skill expression-analysis \
  --command "uv run $SKILL/scripts/expression.py variants --study-dir $D" --output processed/core_variants.tsv
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 11-master --skill expression-analysis \
  --command "uv run $SKILL/scripts/expression.py master --study-dir $D" --output processed/master.parquet \
  --output processed/track_scores_long.parquet --output processed/qc_report.json
```

## Step 12. Contrast (Tier 2)

```bash
uv run $SKILL/scripts/expression.py contrast --study-dir $D --propose
```

It writes `processed/contrasts_proposed.yaml` and prints the review points. Go through the plan
with the researcher: one contrast per scorer (strand, track filter, target and control CURIEs,
tissue aggregation) and `priority.direction_contrast`. Save the version they approve, with their
edits, as `contrasts.yaml`, and record the approval:

```bash
cp $D/processed/contrasts_proposed.yaml $D/contrasts.yaml        # then apply the researcher's edits, if any
uv run $SKILL/scripts/study.py decide --study-dir $D --key contrast.plan --value approved --step 12 \
  --note "Plan reviewed with the researcher: <what was changed, or 'accepted as proposed'>"
uv run $SKILL/scripts/expression.py contrast --study-dir $D
```

Review questions `contrast.review:mixed_assays:{contrast}`, `contrast.review:mixed_sources:{contrast}`
and `contrast.review:single_tissue:{contrast}:{role}` take `accept` or `abort` (abort, edit
`contrasts.yaml`, e.g. its `track_filter`, and run again).

```bash
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 12-contrast --skill expression-analysis \
  --command "uv run $SKILL/scripts/expression.py contrast --study-dir $D" \
  --output contrasts.yaml --output processed/tissue_contrast.parquet --output processed/tissue_contrast_meta.json \
  --output processed/top_variants.tsv --output processed/top_variants.json \
  --decision-key contrast.plan --decision-key '<each contrast.review:* key answered>'
```

## Step 13. Interpretation and report

Atlas deep links of the top variants (`table` reads JSON only and has no `--output`):

```bash
uv run --no-project --python 3.13 --with "$LINKS_DEPS" $DEPS/alphagenome-atlas-website-links/scripts/alphagenome_atlas_links.py \
  table --input $D/processed/top_variants.json > $D/raw/atlas_links.md
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 13-atlas-links --skill alphagenome-atlas-website-links \
  --command "uv run --no-project --python 3.13 --with \"$LINKS_DEPS\" $DEPS/alphagenome-atlas-website-links/scripts/alphagenome_atlas_links.py table --input $D/processed/top_variants.json > $D/raw/atlas_links.md" \
  --output raw/atlas_links.md --source-version python=3.13
```

Optional, per variant, with `alphagenome-single-variant-analysis` (model API, about 20 s per call;
ask the researcher how many of the top variants and which target CURIE). The text of
`interpret_splicing.py` goes into the report; `visualize_variant_effects.py` draws REF/ALT tracks:

```bash
CT=UBERON:AAAAAAA          # a target CURIE chosen by the researcher
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="variant") c=i; next} {print $c}' $D/processed/top_variants.tsv | head -n 5 > $D/processed/interpret_variants.txt
for v in $(cat $D/processed/interpret_variants.txt); do
  chrom=${v%%:*}; rest=${v#*:}; pos=${rest%%:*}; alleles=${rest#*:}; ref=${alleles%>*}; alt=${alleles#*>}
  out=interpret_${chrom}_${pos}_${ref}_${alt}.txt
  uv run --with "alphagenome<2" $DEPS/alphagenome-single-variant-analysis/scripts/interpret_splicing.py \
    --chrom $chrom --pos $pos --ref $ref --alt $alt --ontology_id $CT > $D/raw/$out
  uv run $SKILL/scripts/manifest.py add --study-dir $D --step 13-interpret --skill alphagenome-single-variant-analysis \
    --command "uv run --with 'alphagenome<2' $DEPS/alphagenome-single-variant-analysis/scripts/interpret_splicing.py --chrom $chrom --pos $pos --ref $ref --alt $alt --ontology_id $CT > $D/raw/$out" \
    --output raw/$out --source-version alphagenome='<2'
  vdir=figures/variants/${chrom}_${pos}_${ref}_${alt}
  uv run --no-project --python 3.13 --with "$SVA_DEPS" $DEPS/alphagenome-single-variant-analysis/scripts/visualize_variant_effects.py \
    --chrom $chrom --pos $pos --ref $ref --alt $alt --gene $G --tissue "<target label>" --ontology $CT \
    --tracks regulatory --output_dir $D/$vdir
  uv run $SKILL/scripts/manifest.py add --study-dir $D --step 13-visualize --skill alphagenome-single-variant-analysis \
    --command "uv run --no-project --python 3.13 --with \"$SVA_DEPS\" $DEPS/alphagenome-single-variant-analysis/scripts/visualize_variant_effects.py --chrom $chrom --pos $pos --ref $ref --alt $alt --gene $G --tissue \"<target label>\" --ontology $CT --tracks regulatory --output_dir $D/$vdir" \
    $(cd $D && for f in $vdir/*; do [ -f "$f" ] && echo --output $f; done) --source-version alphagenome='<2'
done
```

`$(cd $D && for f in …; do echo --output $f; done)` turns the files the run wrote into one
`--output` each (`manifest.py` accepts files, not folders); the same form serves any step that
writes many files. `processed/interpret_variants.txt` keeps the list of variants sent to this
step. Then the report, which reads the manifest, so record the steps above first:

```bash
uv run $SKILL/scripts/expression.py report --study-dir $D
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 13-report --skill expression-analysis \
  --command "uv run $SKILL/scripts/expression.py report --study-dir $D" --output report.md
```

Read `report.md` back before handing it over: every claim should be a predicted effect with its
tissue and scorer, the skills section should match the manifest, and nothing should read as a
diagnosis (`interpretation.md`).

## Step 14. Workbench

```bash
uv run $SKILL/scripts/expression.py dashboard --study-dir $D                  # serves on http://localhost:2718
uv run $SKILL/scripts/expression.py dashboard --study-dir $D --export-html $D/figures/workbench.html
uv run $SKILL/scripts/manifest.py add --study-dir $D --step 14-workbench --skill expression-analysis \
  --command "uv run $SKILL/scripts/expression.py dashboard --study-dir $D --export-html $D/figures/workbench.html" \
  --output figures/workbench.html
```

marimo runs in an ephemeral `uv` environment, so the first run needs network access. The defaults
and figures are explained in `interpretation.md`. Only the export writes a file to record.
marimo's session folder (`__marimo__/`) is written under `$D/.workbench/`, never inside the skill.

Reproducing the reference (`expression.py reproduce --reference uchl1`): after a pass it prints a
"Headline results" block and a `Summary:` line with the total number of tracks (22). Report both in the
answer to the researcher; a reproduction that is not reported cannot be checked.

## When a command fails

| Message | Meaning | Do |
|---|---|---|
| `N decision(s) needed: run study.py pending.` | an open question | Conventions: show, ask, `decide`, run again |
| `Missing decision '{key}'. Run Workflow step {N} …` | a stage needs an earlier answer | go back to step N |
| `invalid decision key(s) …` (`decide`) | a key with spaces, newlines or another shape; nothing was recorded | copy the key from `study.py pending`; one `--key` per key (Conventions) |
| `… not found …; run Workflow step N` | a raw or processed input is missing | run that step |
| `FAIL file=… field=… skill=…` / `DependencyOutputError` | a dependency changed its output format | stop; compare with `dependency_outputs.md`; tell the researcher |
| `interval is … above --max-window` / refused window | the window is wider than 1000 bp | pass `study.py window --print-max-window` |
| `ImportError: cannot import name 'dna_client'` | `alphagenome` 2.0.0 | add `--with "alphagenome<2"` |
| `h5py` build failure in the links or AVI script | Python 3.14 from the skill's pyproject | `--no-project --python 3.13 --with "$…_DEPS"` |
| `QC failed for N variant(s)` (`master`) | raw files disagree | check that steps 3, 9 and 10 used the same window and pull |
| `Reference 'uchl1': … mismatch` (`reproduce`) | a code or data change moved a reference number | read the stage/metric table it prints before changing anything |
