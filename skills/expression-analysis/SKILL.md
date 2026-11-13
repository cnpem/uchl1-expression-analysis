---
name: expression-analysis
description: "Compares AlphaGenome-predicted regulatory and expression effects of every SNV in a human gene's regulatory window between target and control tissues or cell types (UBERON/CL), and opens an interactive marimo workbench. Routes to the DeepMind science skills for GENCODE v46 gene structure, OLS ontology terms, Atlas AVI and track catalogue, and to AlphaGenome Atlas track scores. Use this skill whenever someone asks which variants in a gene change expression in one tissue versus another, wants a saturation-mutagenesis or promoter/enhancer variant triage, or mentions tissue-specific regulatory effects with AlphaGenome — even if they do not say 'expression analysis'. Research use only."
metadata:
  version: 0.1.0
---

# expression-analysis

## Overview

For one human gene (GRCh38), this skill takes every possible SNV of the gene's regulatory window
(by default the promoter/enhancer unit that holds the TSS, through the end of exon 1), collects the
AlphaGenome Atlas predictions for each of them, and asks: which variants change the predicted
signal (RNA-seq, CAGE, DNase, ATAC) in the target tissues, and do they do so differently in the
control tissues? The answer comes as a Tier 1 table (one row per SNV), a Tier 2 tissue contrast,
a ranked top list with Atlas links, a Portuguese `report.md` and an interactive marimo workbench.

How the work is split, and why:

- **Dependency skills fetch, this skill analyses.** Gene structure, ontology terms, track
  catalogue and scores come from the DeepMind science skills and the K-Dense `alphagenome` skill.
  Their outputs are saved unedited in the study's `raw/` folder, and every command is recorded in
  `manifest.json`. The stages of `scripts/expression.py` read only the study folder, so a study
  can be re-run offline and every number traces back to a recorded command.
- **The researcher decides, the agent asks.** When the data leave a choice open (which
  transcript, a regulatory feature that runs into exon 1, a tissue without Atlas tracks, a
  contrast that mixes assays), the stage stops with exit 1 and `N decision(s) needed`, and writes
  the question with its options and a recommendation. The agent shows it to the researcher and
  records the researcher's answer; it does not pick one itself, because these choices change
  which variants and tissues the conclusions rest on.
- **Research use only.** Everything here is a predicted molecular effect used to prioritise
  experiments. It is not a measurement, not a clinical classification and never a diagnosis.

Read `references/interpretation.md` before reporting results, and `references/workflow.md` for
the complete command of every step.

## Prerequisites

1. **`uv`** on PATH (the `uv` skill). Every script runs with `uv run`, which installs its PEP 723
   dependencies into a cached environment, so nothing is installed into the system Python.
2. **The dependency skills** of the next section. `expression.py check` looks for them next to
   this skill, in `.claude/skills/`, `.agents/skills/`, `~/.claude/skills/` and
   `~/.agents/skills/`.
3. **Licence notices.** `embl-ebi-ols`, `ensembl-database` and
   `alphagenome-single-variant-analysis` ask for a notice file in `.licenses/` of the workspace
   after the researcher has been shown the terms of use. `check` lists the missing ones; follow
   the Prerequisites section of that dependency together with the researcher, and create the
   file only after the researcher says so, because the file records that they saw the terms.
4. **`ALPHAGENOME_API_KEY`** in `~/.env`, set with the `credentials` skill. The DeepMind scripts
   read `~/.env` themselves; the K-Dense scripts read only the environment, so they run as
   `uv run --env-file ~/.env …`. The key stays out of command lines and of the conversation:
   do not source, print or copy `~/.env`; `check` only tests whether the key exists.
5. **A study directory** chosen by the researcher, outside this skill's folder (the skill folder
   is read-only, and a study holds a researcher's data and decisions).

## Dependencies

Cited by the `name` of each skill. Versions tested on 2026-11-12 are in
`references/dependency_outputs.md` (with each raw file's command and fields).

| Skill | Origin | Steps | Used for |
|---|---|---|---|
| `uv` | DeepMind | all | running every script with `uv run` |
| `credentials` | DeepMind | 1 | storing `ALPHAGENOME_API_KEY` in `~/.env` |
| `alphagenome-variant-impact-score` | DeepMind | 2, 7, 9 | `gtf` (GENCODE v46 structure), `metadata --tracks` (track catalogue), `region` (AVI of every SNV) |
| `embl-ebi-ols` | DeepMind | 5, 6 | `search_ols.py` (UBERON/CL terms), `get_term.py` (descendants) |
| `ensembl-database` | DeepMind | 3 | `get-sequence` (reference bases: REF check, start codon) |
| `alphagenome` | K-Dense | 10 | `atlas_query.py scores` (per-track scores per tissue) |
| `alphagenome-single-variant-analysis` | DeepMind | 13 | `interpret_splicing.py`, `visualize_variant_effects.py` on the top variants |
| `alphagenome-atlas-website-links` | DeepMind | 13 | `table --input` (Atlas deep links for the top variants) |

**Installation is the researcher's decision.** Installing a skill brings third-party code onto
their machine and, for several of these, terms of use. So no script of this skill installs
anything: `expression.py check` prints each missing skill with its origin and its install command
(`npx skills add google-deepmind/science-skills --skill <name>` or
`npx skills add K-Dense-AI/scientific-agent-skills --skill alphagenome`) and stops. Show that
list to the researcher and ask about one skill at a time; an approval covers only the skill it
names. Licence notices and keys follow each dependency's own Prerequisites, with the researcher.

## Quick Start

Variables used in every command (set them once per shell; `references/workflow.md` explains them):

```bash
SKILL=.claude/skills/expression-analysis   # this skill's folder
DEPS=.claude/skills                        # folder holding the dependency skills
D=studies/mygene-expression                # study directory chosen by the researcher
G=MYGENE                                   # HGNC symbol of the gene
```

Reproduce the bundled reference study (offline, a few seconds; six stages compared with
`assets/reference/uchl1/expected.json`), then open its workbench:

```bash
uv run $SKILL/scripts/expression.py reproduce --reference uchl1 --keep /tmp/uchl1-demo
uv run $SKILL/scripts/expression.py dashboard --study-dir /tmp/uchl1-demo
```

In the answer to the researcher, report the "Headline results" block that `reproduce` prints after a
pass, and the total number of tracks from its `Summary:` line (22 for UCHL1). A reproduction that is
not reported cannot be checked.

A new study, in outline only. These five lines are not a full run: the dependency captures
(steps 2-3, 5-7 and 9-10) and the stages between them (`cohorts`, `variants`, `master`,
`contrast`, `report`) are left out here and are in `references/workflow.md`, step by step:

```bash
uv run $SKILL/scripts/expression.py check                       # step 1: prerequisites (prints, never installs)
uv run $SKILL/scripts/study.py init --study-dir $D --gene $G    # step 1: study.yaml from the template
uv run $SKILL/scripts/expression.py regions --study-dir $D      # step 3: regions for the captures (after step 2)
uv run $SKILL/scripts/expression.py locus --study-dir $D        # step 4; then cohorts, variants, master, contrast, report
uv run $SKILL/scripts/expression.py dashboard --study-dir $D    # step 14: workbench
```

## Utility Scripts

All in `scripts/`, run as `uv run $SKILL/scripts/<script> <subcommand> --help` for the options.
Stage subcommands take `--study-dir D` and read and write only inside `D`.

| Command | Does |
|---|---|
| `expression.py check` | Checks `uv`, dependency skills, licence notices, `ALPHAGENOME_API_KEY`; prints install commands, runs none |
| `expression.py regions [--only NAME]` | Prints the regulatory search region (step 3), the window and `max_window` (steps 9-10), so no coordinate is computed by hand |
| `expression.py locus` | Step 4: transcript, TSS, exon 1 and window → `processed/locus.json`, `transcripts_in_window.tsv`, Figure 1 |
| `expression.py cohorts` | Step 8: CURIEs × Atlas track catalogue → `processed/track_registry.tsv`, `ontology_resolution.json` |
| `expression.py variants` / `master` | Step 11: SNV table (REF from the sequence) and Tier 1 table with the two QC gates |
| `expression.py contrast [--propose]` | Step 12: proposes the contrast plan, then runs the approved one → tissue contrast, top variants |
| `expression.py report` | Step 13: `report.md` (study, decisions, results, skills used, limits, citations) |
| `expression.py dashboard [--port N] [--export-html P]` | Step 14: marimo workbench, or a static HTML export |
| `expression.py reproduce --reference uchl1 [--until S] [--keep DIR]` | Re-runs the stages on the bundled reference and compares with `expected.json` |
| `study.py init / validate / pending / decide / require / window` | Study file, researcher decisions, open questions, `--print-max-window` |
| `manifest.py add / show` | Records a step (command, outputs, versions, decision keys); prints the manifest |
| `models.py validate --dir D/raw` | Checks every raw file against the fields the dependency wrote on 2026-11-12 |
| `regulatory_features.py overlap --region R --output P` | Ensembl Regulatory Build features of a region (this skill's only own HTTP client) |

## Workflow

Each step ends with `manifest.py add` (the report reads the skills used from the manifest, and the
dependencies ask to be named in any output). Complete commands: `references/workflow.md`.

**Decisions (every stage).** When a stage exits 1 with `N decision(s) needed: run study.py
pending.`, run `study.py pending --study-dir $D`, show each question, its options and the
recommendation to the researcher, record the researcher's answer with
`study.py decide --study-dir $D --key K --value V --step N --note "…"` (one `--key` per key when
several share an answer), and run the stage again. The recommendation is the stage's suggestion,
not an answer. `decide` copies the options and the recommendation into `decisions.yaml` next to the
answer, so the record shows what was offered; for a question you ask yourself (the CURIEs of step
5), pass them with `--option` (repeated) and `--recommended`. When the researcher has delegated
("take the recommended options"), record the recommendation and say so in the note. The final
answer lists every decision raised (key, options, recommended, chosen, reason) and, in a separate
list, the choices you made that no stage raised (e.g. a search word or a scorer list), so the
researcher can tell their decisions from yours.

1. **Prerequisites and study.** `expression.py check`; if anything is missing, show the list and
   ask (one skill at a time, see Dependencies). `study.py init`, then ask the researcher for the
   gene, the target and control tissues or cell types (in words; CURIEs come in step 5), the
   scorers and the window policy; edit `study.yaml` and run `study.py validate`.
2. **Gene structure.** `alphagenome-variant-impact-score` `gtf --all_transcripts --exons --cds
   --utr --format json` → `raw/gtf.json`. Use this GTF because it is GENCODE v46, the annotation
   the Atlas scores were computed on; its starts are 0-based as written and the stages normalise
   them, so the file is kept as is.
3. **Regulatory Build and sequence.** `regions --only regulatory_search` gives the gene span
   ± `regulatory_search_upstream_bp`; `regulatory_features.py overlap` and `ensembl-database
   get-sequence` of that region → `raw/regulatory.json`, `raw/sequence_*.txt`. The features build
   the window; the sequence checks REF alleles and the start codon.
4. **Locus and window.** `expression.py locus` → `processed/locus.json`, Figure 1. Questions
   `locus.*` (transcript, exon-1 overlap, gaps between features). Show Figure 1 to the researcher.
5. **UBERON/CL terms.** `embl-ebi-ols` `search_ols.py` for each tissue named; show labels and
   definitions with the hit you recommend; the researcher picks the CURIEs → decisions
   `cohort.target`, `cohort.control` recorded with the hits offered (`--option`, `--recommended`),
   and the same CURIEs in `study.yaml`. Never substitute another term on your own.
6. **Descendants (optional).** `get_term.py --relations hierarchicalChildren,hierarchicalDescendants`
   for each CURIE that step 7 finds without tracks or without `RNA_SEQ` tracks, so that step 8 can
   offer a served descendant (a cortical lobe sits several levels below its parent, so children
   alone miss it).
7. **Track catalogue.** `metadata --tracks --ontology C --scorer S --output
   raw/tracks_{S}_{C}.parquet` for every CURIE × scorer (no match writes no file; record it
   anyway), plus `metadata --tracks --scorer S --output raw/tracks_{S}_catalogue.parquet` for the
   whole catalogue of each scorer, which step 8 needs to offer only descendants that have tracks.
8. **Cohorts.** `expression.py cohorts` → track registry. Questions `cohort.unserved:{curie}`
   (drop, or replace by a served descendant, then step 7 for it),
   `cohort.scorer_without_cohort:{scorer}` and `cohort.no_rna_seq:{role}` (add a descendant with
   `RNA_SEQ` tracks, the one with most recommended, then step 7 for it; or abort).
9. **Saturation and AVI.** `region --region W --max_window_size M --min_phred 0
   --include_features --include_track_info` with `W` from `regions --only window` and `M` from
   `study.py window --print-max-window` → `raw/avi_region.parquet`, the universe of SNVs.
10. **Scores per tissue.** K-Dense `atlas_query.py scores --interval W --max-window M --scorers S
    --ontology C…` (plus `--gene G` for `RNA_SEQ`) once per scorer → `raw/scores_{S}.tsv`.
11. **Variants and Tier 1.** `expression.py variants` (questions `variants.ref_mismatch:*`,
    `variants.missing_in_atlas:*`), then `expression.py master` with two QC gates (coordinates
    and REF; AVI attribution invariants). A failed gate stops the study until its cause is found.
12. **Contrast (Tier 2).** `expression.py contrast --propose`; review the plan and its review
    points with the researcher; save the approved plan as `contrasts.yaml`; decision
    `contrast.plan`; `expression.py contrast` (questions `contrast.review:*`).
13. **Interpretation and report.** Atlas links `table --input processed/top_variants.json >
    raw/atlas_links.md`; optionally `alphagenome-single-variant-analysis` on the top variants with
    a target CURIE; `expression.py report` → `report.md`.
14. **Workbench.** `expression.py dashboard` (or `--export-html`) for the researcher's triage.

Copy this checklist into the conversation and tick it as you go:

```
expression-analysis progress:
- [ ] Step 1: check passes (missing skills asked one at a time); study.yaml validated; manifest started
- [ ] Step 2: raw/gtf.json (GENCODE v46) recorded
- [ ] Step 3: regions printed; raw/regulatory.json and raw/sequence_*.txt recorded
- [ ] Step 4: locus run; locus.* answered by the researcher; Figure 1 shown
- [ ] Step 5: OLS searches recorded; cohort.target and cohort.control decided by the researcher, with options
- [ ] Step 6: descendants pulled for CURIEs without tracks or without RNA_SEQ tracks (or not needed)
- [ ] Step 7: tracks per CURIE x scorer and per-scorer catalogues recorded (no-match pairs too)
- [ ] Step 8: cohorts run; cohort.unserved / scorer_without_cohort / no_rna_seq answered
- [ ] Step 9: raw/avi_region.parquet with --max_window_size from study.py window
- [ ] Step 10: raw/scores_{scorer}.tsv for every kept scorer, same window and pull
- [ ] Step 11: variants and master run; QC gates passed
- [ ] Step 12: plan reviewed with the researcher; contrasts.yaml saved; contrast.plan approved; contrast run
- [ ] Step 13: atlas_links.md (and interpretations) recorded; report.md written and read back
- [ ] Step 14: workbench opened or exported for the researcher
```

## Rate Limiting

Only `regulatory_features.py` makes its own requests: 15 requests/s (the Ensembl REST limit the
`ensembl-database` skill states), with a file lock at `~/.cache/alphagenome-gene-skills/ensembl.lock`
so concurrent subagents share the budget, and retries with exponential backoff on 429/5xx.
Every other request comes from a dependency script, which applies its own limits (Ensembl
15 requests/s; the Atlas client's `--workers`, default 8). Each loop iteration in the Workflow is a
separate process, so a limiter cannot see the others: run the loops one at a time rather than
spreading them over parallel subagents. The model API calls of step 13 take about 20 s each;
run them sequentially.

## Common Mistakes

- **Annotating exons from Ensembl REST instead of `raw/gtf.json`.** The Atlas scores were computed
  on GENCODE v46; Ensembl release 116 can differ in transcripts and boundaries. Ensembl is used
  only for the Regulatory Build, the sequence and IDs.
- **"Fixing" the GTF by hand.** The `gtf` output writes starts 0-based (ends 1-based); every
  stage applies `models.normalize_gtf` once. Editing `raw/` shifts every coordinate twice and
  breaks `models.py validate`.
- **Omitting `--max_window_size` / `--max-window` for a window above 1000 bp.** Both scripts refuse
  wider windows by default; pass the value of `study.py window --print-max-window`.
- **Mixing scores of different Atlas pulls or windows.** `avi_region.parquet` and every
  `scores_{S}.tsv` describe one window captured together; a re-capture of one file means
  re-capturing all, or the QC gates of `master` fail (and if they did not, tissues would be
  compared across different data). The manifest records each pull date.
- **Presenting a filter as a clinical threshold.** |quantile| ≥ 0.99, AVI Phred 20 or a large Δ
  rank predicted effects; none is a pathogenicity cut-off. The same holds for any τ quoted from
  another study.
- **`get_term.py` with `hierarchicalChildren` only.** Deep descendants (a frontal cortex under the
  cerebral cortex) are then never offered; request `hierarchicalChildren,hierarchicalDescendants`.
- **Skipping the per-scorer catalogue (`tracks_{S}_catalogue.parquet`).** Without it the
  `cohort.unserved` question lists every descendant, most of them without tracks.
- **`--gene` on CAGE, DNASE or ATAC.** Only gene-centric scorers (`RNA_SEQ`, `RNA_SEQ_ACTIVE`,
  `POLYADENYLATION`, `SPLICE_*`) carry a gene; use `--gene G` for those only.
- **`gtf --output` or `--format json` without the `sed` filter.** `--format json` prints to stdout
  and ignores `--output`, after one `Loading GTF from …` line that the filter drops.
- **An unpinned `alphagenome` for model-API scripts.** `alphagenome` 2.0.0 removed `dna_client`;
  run `interpret_splicing.py` and `visualize_variant_effects.py` with `alphagenome<2`.
- **Running the links skill through its `pyproject.toml`.** It asks for Python ≥ 3.14, where
  `h5py` has no wheel; use `uv run --no-project --python 3.13 --with "$LINKS_DEPS"`, and give
  `table` a JSON file (it does not read TSV).
- **Choosing for the researcher, or editing `decisions.yaml` by hand.** Record answers with
  `study.py decide`, which also clears the pending question and keeps the timestamps.
- **Installing a dependency, or writing a `.licenses/` file, without the researcher's explicit
  approval for that skill.** Ask one skill at a time.
- **Sourcing or printing `~/.env`.** Use `uv run --env-file ~/.env`; keys never reach the screen.

## Outputs

```
<study_dir>/                         # chosen by the researcher, never inside the skill
  study.yaml                         # the request (gene, cohorts, scorers, window policy, top_n)
  decisions.yaml                     # the researcher's answers, with timestamps and notes
  contrasts.yaml                     # the approved contrast plan (step 12)
  manifest.json                      # provenance: every command, output, version and decision key
  raw/                               # dependency outputs, unedited (gtf.json, regulatory.json,
                                     #   sequence_*.txt, ols_*.json, tracks_*.parquet,
                                     #   avi_region.parquet, scores_{scorer}.tsv, atlas_links.md,
                                     #   interpret_*.txt)
  processed/                         # this skill's outputs: locus.json, transcripts_in_window.tsv,
                                     #   track_registry.tsv, ontology_resolution.json, tissue_cohorts.md,
                                     #   core_variants.tsv, master.parquet, track_scores_long.parquet,
                                     #   qc_report.json, contrasts_proposed.yaml, tissue_contrast.parquet,
                                     #   tissue_contrast_meta.json, top_variants.tsv/.json,
                                     #   interpret_variants.txt (variants sent to step 13),
                                     #   pending_decisions.json
  figures/                           # fig1_genomic_architecture.html, workbench.html (export),
                                     #   variants/<chr>_<pos>_<ref>_<alt>/ (REF/ALT plots, step 13)
  report.md                          # Portuguese report: study, decisions, results, skills used, limits, citations
```

## Citation

- Avsec, Ž. et al. (2026). Advancing regulatory variant effect prediction with AlphaGenome.
  *Nature* 649:1206. (AlphaGenome model and Atlas predictions.)
- Cheng et al. (2026). AlphaGenome Atlas report. medRxiv, doi:10.64898/2026.09.16.26363192.
- DeepMind science skills (`alphagenome-variant-impact-score`, `alphagenome-single-variant-analysis`,
  `alphagenome-atlas-website-links`, `embl-ebi-ols`, `ensembl-database`, `credentials`, `uv`):
  https://github.com/google-deepmind/science-skills.
- Kassis, T. et al. (2026). Scientific Agent Skills: A Library of Procedural Knowledge for Research
  Agents. arXiv:2609.00065 (the K-Dense `alphagenome` skill).
- Data behind the dependency outputs: GENCODE v46, Ensembl (Regulatory Build, sequence), OLS
  (UBERON, CL); cite them as their skills ask. `report.md` lists the skills actually used, read
  from `manifest.json`.
