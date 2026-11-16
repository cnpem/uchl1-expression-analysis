# Workflow reference: expression-analysis

The complete commands of the Workflow steps of `SKILL.md`. Formats and fields of every raw file are
in `dependency_outputs.md`; how to read the results is in `interpretation.md`.

## Contents

- [Conventions](#conventions): paths, what each command records, researcher decisions
- [Step 0. Exploration](#step-0-exploration)
- [Step 1. Study from the exploration](#step-1-study-from-the-exploration)
- [Steps 2-3. Gene structure, Regulatory Build and sequence](#steps-2-3-gene-structure-regulatory-build-and-sequence)
- [Step 4. Locus and window](#step-4-locus-and-window)
- [Steps 5-6. Terms and descendants](#steps-5-6-terms-and-descendants)
- [Step 7. Track catalogue](#step-7-track-catalogue)
- [Step 8. Cohorts](#step-8-cohorts)
- [Steps 9-10. Saturation, AVI and scores](#steps-9-10-saturation-avi-and-scores)
- [Step 11. Variants and Tier 1](#step-11-variants-and-tier-1)
- [Step 12. Contrast (Tier 2)](#step-12-contrast-tier-2)
- [Step 13. Interpretation and report](#step-13-interpretation-and-report)
- [Step 14. Dashboard](#step-14-dashboard)
- [When a command fails](#when-a-command-fails)

## Conventions

Run the commands from the researcher's workspace (the folder whose `.claude/skills/` or
`.agents/skills/` holds the dependency skills). `SKILL` is this skill's folder, `OUT` the folder
the researcher chose, `D` the study directory (normally `OUT/<gene>-study`); write them out in full.
Every command below is one line that runs unchanged in PowerShell, bash and zsh: no variables,
`$(...)`, pipes or `~` to expand, because a shell difference would otherwise change what runs and
what the manifest says ran.

**What gets recorded.** Each `expression.py fetch`, each stage (`locus` … `contrast`, `report`,
`dashboard`) and `study.py init --from-exploration` appends one step to `manifest.json` when it
succeeds: the exact command (for a fetch, the dependency command it ran with its pinned packages,
plus the `expression.py` call as `invoked_by`), the outputs, the versions and the decision keys the
stage used. The report lists the skills used from the manifest, and the dependencies ask to be named
in any output. `fetch … --dry-run` prints the commands and records nothing. `manifest.py add` is
only for a step run outside the skill; give it the command exactly as it ran.

**Researcher decisions.** A stage that meets an open choice writes it to
`processed/pending_decisions.json` and exits 1 with `N decision(s) needed: run study.py pending.`
Each question has a title, the question, why it matters, what it is about (names with their
identifiers: CURIEs, transcript and region IDs, coordinates, counts) and options with a label, the
effect of choosing it and a code. Present it in the researcher's language (the language of the
conversation), with the names and identifiers as given, and never show the codes: they are what
`decide` records and what the stages read. `study.py pending --json` gives the same in JSON.
`decide` refuses a code the question did not offer (a code such as `fixed:<value>` takes any text after
its prefix; Ensembl IDs match with or without their version) and records the option's label next to it
(`value_label`), so `decisions.yaml` and the report show what the researcher saw.

**Rules the skill applies without asking.** Checks the researcher cannot judge, and choices that do not
change the analysis they steer, are applied by fixed rules and recorded in `decisions.yaml` with
`decided_by: skill` and the rule in words: the `regulatory_unit` window keeps every regulatory region
at the start of the gene and all of exon 1, which is how that window is defined (step 4); a prediction
type without tracks in one group is left out (step 8), positions where the Atlas and the reference sequence disagree are left out (step 11; above
1% of the window the researcher is asked first), and the triage's contrast plan is generated from the
groups (step 12). The rules are evaluated again on every run (a recapture that fixes a position brings
it back) and never overwrite a researcher's answer. Report them to the researcher with the results.

```bash
uv run SKILL/scripts/study.py pending --study-dir D
uv run SKILL/scripts/study.py decide --study-dir D --key "locus.exon1_overlap:FEATURE" --value keep_full_unit --step 4 --note "Question as shown; the researcher's reason"
```

`--value` is parsed as YAML, so a list is written `"[UBERON:0000955, UBERON:0002107]"`. Quote keys
(they contain `:`). Several keys with the same answer go in one call, one `--key` each. Then run the
stage again. `decide` copies the question's `options` and `recommended` into `decisions.yaml`, so the
record shows what was offered next to what was chosen; for a question you ask yourself, pass what you
offered with `--option` (once per option) and `--recommended`. If the researcher delegated the choice,
record the recommendation and write that in the note. A key must look like
`<area>.<topic>[:<subject>]` (as `pending` prints it); `decide` refuses anything else.

## Step 0. Exploration

```bash
uv run SKILL/scripts/expression.py check
uv run SKILL/scripts/expression.py explore --gene MYGENE --tissue "<words of the researcher>" --tissue "<another>" --out OUT
```

- `check` exit 1 lists what is missing and hands it to the `alphagenome-skills-setup` skill, which asks
  the researcher item by item (skills, terms of use, the key in the project's `.env`) and installs only
  what they approve. If that skill is not installed, `check` prints how to install it; ask first. Install
  nothing from this skill. Run `check` again until it passes.
- `explore` asks nothing and writes `OUT/exploration/` (see `SKILL.md` Outputs). Options:
  `--scorers` (default `RNA_SEQ CAGE DNASE ATAC`), `--flank-bp` (Regulatory Build search around the
  gene, default 2000, the study default), `--floor` (relative activity below which the gene reads as
  not expressed, default 0.05 of its maximum). Re-running it reuses what was captured, so adding a
  tissue costs one OLS search.
- Network: one call each for GENCODE, the Regulatory Build and the sequence; one OLS search per
  tissue, one OLS descendants query per tissue whose best term has no RNA-seq track; one Atlas
  catalogue per scorer plus one with sources and assays; one Atlas score query of a single SNV. No
  window is scored.

Show `summary.md` to the researcher (it is written for someone who does not program). What it holds,
and why each part matters for the choice:

| Section | Holds | Why |
|---|---|---|
| 3. Windows | 3-5 windows derived from the features (regulatory unit at the TSS through exon 1, with gaps bridged, the TSS feature alone, TSS ± 500 bp, all promoters/enhancers of the region), with size, SNVs = 3 × bp, Atlas rows per scorer and features included or left out | the window is the set of variants every conclusion rests on, and its size is the cost |
| 4. Terms | per tissue text: OLS candidates (UBERON, CL, EFO, CLO) with tracks per scorer, data source and RNA-seq protocol; catalogue terms whose name contains the text (subtypes, cell lines); served descendants of a term without RNA-seq tracks | a term without tracks cannot be contrasted; one track means no replicate |
| 5. Baseline | `RNA_SEQ_ACTIVE` of the gene per tissue, relative to its maximum, with a warning below the floor | a gene that is silent in the target makes predicted fold changes large and unreliable there |
| 6. Scorers | which recommended terms each scorer serves | a scorer without tracks in one group is dropped later |
| 7. Limits | SNVs only (indels through the model API, `score_variants.py` of the `alphagenome` skill), window scope | what the study cannot answer |

### Changing the groups after the results

When the researcher wants other target or control tissues (for example, controls where the gene is also
expressed), that is a redefinition of the study, and it follows the same path as the definition. Choose from
the exploration, never from ad-hoc code:

1. Show `OUT/exploration/activity.md`. It lists every catalogue term with the gene's predicted activity
   (relative to the most active), kind of sample, assays, sources and number of measurements. Without it (an
   older exploration), run `expression.py explore --out OUT --activity-only`; it captures nothing. For
   tissues the researcher names in their own words, run `explore` again with those `--tissue` texts: it
   reuses the captures and only searches the new terms.
2. Discuss the choice with the researcher. Prefer terms where the gene is expressed. Keep the kinds of
   sample and the assays balanced with the target. Say what the new comparison asks, which may differ from
   the original question.
3. Record it: `study.py decide --key cohort.control` (or `cohort.target`) with the new list, the terms
   offered (`--option`) and the recommended one, and a note saying why the study was redefined. Then
   `fetch tracks`, `cohorts`, `fetch scores`, `variants`, `master`, `contrast`, `report`, as in steps 7-13.
   Tissues taken out of a group can stay visible as context tissues on the dashboard.

## Step 1. Study from the exploration

Ask the researcher, with the words of the summary page: which window (section 3, by its name; you
pass its id to `--window`), which terms form the target and the control groups (section 4) with a short
name for each group, and which kinds of prediction to compare (section 6; gene expression, `RNA_SEQ`,
is always included). Then:

```bash
uv run SKILL/scripts/study.py init --study-dir D --from-exploration OUT --window regulatory_unit --target CL:AAAAAAA --target UBERON:BBBBBBB --control UBERON:CCCCCCC --target-label "<label>" --control-label "<label>" --scorers RNA_SEQ CAGE --origin-note "Requested by <researcher> on <date>: <question in one line>"
uv run SKILL/scripts/study.py decide --study-dir D --key cohort.target --value "[CL:AAAAAAA, UBERON:BBBBBBB]" --option "CL:AAAAAAA (<label>)" --option "UBERON:BBBBBBB (<label>)" --option "UBERON:DDDDDDD (<label>)" --recommended "CL:AAAAAAA (<label>)" --step 1 --note "Chosen by the researcher from the exploration; why"
uv run SKILL/scripts/study.py decide --study-dir D --key cohort.control --value "[UBERON:CCCCCCC]" --option "UBERON:CCCCCCC (<label>)" --recommended "UBERON:CCCCCCC (<label>)" --step 1 --note "Chosen by the researcher; why"
uv run SKILL/scripts/study.py validate --study-dir D
```

`init --from-exploration` writes `study.yaml` (window policy from the option: `regulatory_unit`
with its `max_feature_gap_bp`, or `explicit` with the interval) and copies `gtf.json`,
`regulatory.json`, the sequence, the OLS files and the catalogues of the chosen scorers into `raw/`,
each with the manifest step that captured it (`copied_from`), so steps 2, 3, 5 and part of 7 need no
new call. Values in `<…>` and CURIEs such as `UBERON:AAAAAAA` are placeholders. Substituting a term
yourself is the researcher's choice, not yours; a CURIE cannot be in both cohorts.

A study without an exploration still works: `study.py init --study-dir D --gene G`, edit `study.yaml`
with the researcher, then the fetches of steps 2-7.

## Steps 2-3. Gene structure, Regulatory Build and sequence

Only for a study not created from an exploration (or to capture again):

```bash
uv run SKILL/scripts/expression.py fetch gtf --study-dir D
uv run SKILL/scripts/expression.py fetch regulatory --study-dir D
uv run SKILL/scripts/expression.py fetch sequence --study-dir D
```

- `gtf` is GENCODE v46, the annotation of the Atlas scores, which is why exons come from here and not
  from Ensembl REST. The first run downloads the 318 MB GENCODE Feather to `~/.cache/alphagenome/`.
  The script prints a `Loading GTF from …` line before the JSON; the fetch drops it and keeps the JSON
  as written (starts 0-based; every stage applies `models.normalize_gtf`).
- `regulatory` and `sequence` use the region `expression.py regions --only regulatory_search` prints:
  the gene span over every transcript ± `window.regulatory_search_upstream_bp` on both sides, so the
  upstream flank is covered on either strand. The sequence file name carries the region.

## Step 4. Locus and window

```bash
uv run SKILL/scripts/expression.py locus --study-dir D
```

Questions (answer with the researcher, `--step 4`, then run `locus` again). Labels are shortened here;
`pending` prints them in full with the effect of each option (sizes and variant counts):

| Key | When | Options (label → code) |
|---|---|---|
| `locus.transcript` | 0 or several MANE Select transcripts | each transcript (type, coding length, where the gene starts) → its ID; "Stop to check the gene" → `abort` |
| `locus.feature_overlap:{a}+{b}`, `locus.feature_gap:{a}+{b}` | a further region overlaps or lies close upstream | "Include {a} in the window" → `union` / `bridge_gap`; "Leave {a} out" → `exclude_feature` (the window starts at {b}; nothing beyond is asked); stop → `abort` |
| `locus.no_feature_chain` | no region contains or touches the start of the gene | "Use a window from the exploration page" → `set_explicit_interval` (then set `window.policy: explicit` and `window.explicit_interval` in `study.yaml` to the window they pick); "Pause to choose another window" → `abort` |

"Pause" answers end the stage; to go on, set the window the researcher picks in `study.yaml` and record
a new answer.

Rules (not asked; `decided_by: skill`): `locus.multiple_tss_features` (several regions at the start of
the gene are all kept, `union`) and `locus.exon1_overlap:{feature}` (a region that reaches into exon 1:
the window keeps all of exon 1, `keep_full_unit`). The researcher chose this when choosing the
`regulatory_unit` window; to change it, choose another window. An answer recorded in an older study
still applies.

The exploration's summary already names the questions its window option will raise. The success line
gives the transcript, strand, TSS, window and whether the start codon reads ATG. Open
`figures/fig1_genomic_architecture.html` with the researcher before going on.

## Steps 5-6. Terms and descendants

The exploration captured the OLS searches. For a tissue the researcher adds later, or descendants of
a CURIE that step 7 finds without tracks or without RNA_SEQ tracks:

```bash
uv run SKILL/scripts/expression.py fetch ols --study-dir D --query "cerebral cortex"
uv run SKILL/scripts/expression.py fetch ols-term --study-dir D --curie UBERON:AAAAAAA
```

`fetch ols` writes `raw/ols_search_<slug>.json` (slug: lower case, `_` for other characters), searches
UBERON, CL, EFO and CLO by default (`--ontology` to narrow) and prints each hit's CURIE, label and
definition to show the researcher. A plain word may not return the intended term in the first rows:
search again with the full label (both files are kept). `fetch ols-term` writes
`raw/ols_term_<PREFIX>_<digits>.json` with `hierarchicalChildren` and `hierarchicalDescendants`,
because a served tissue can sit several levels below the one asked for.

## Step 7. Track catalogue

```bash
uv run SKILL/scripts/expression.py fetch tracks --study-dir D
```

It queries every CURIE of the study (`study.yaml` cohorts, the cohort decisions and their `replace:` /
`add:` answers) × every scorer of `study.yaml`, skipping pairs already captured (`--force` captures
again; `--curie`, `--scorer` narrow it; `--catalogue` also fetches `tracks_{S}_catalogue.parquet`,
which an exploration-based study already has). A pair without tracks writes no file and is recorded
with `result=no_tracks`, which is how `cohorts` tells it from a pair never queried.

## Step 8. Cohorts

```bash
uv run SKILL/scripts/expression.py cohorts --study-dir D
```

| Key | When | Options (label → code) |
|---|---|---|
| `cohort.unserved:{curie}` | a term has no AlphaGenome prediction | "Remove {name} from the group" → `drop_curie` (recommended while the group keeps a term with RNA-seq, otherwise the best replacement is); "Use {more specific term} instead" → `replace:{CURIE}` (terms with predictions, most RNA-seq data first, marked tissue or cell type); "Pause to choose other tissues" → `abort` |
| `cohort.no_rna_seq:{role}` | no term of a group has an RNA-seq prediction | "Add {more specific term} to the group for gene expression" → `add:{CURIE}` (most RNA-seq data first, recommended); "Pause" → `abort` |

A term chosen as a replacement in one question is not offered again in the others on the next run.
Questions a stage no longer raises (answered elsewhere, made moot, or now rules) leave `pending` when the
stage runs again. For the window, the questions about neighbouring regions come one at a time, from the
start of the gene outwards, since each depends on the answer before it. When
a question's "For the agent" line asks for `fetch ols-term`, run it and `cohorts` again before asking.

Rule (not asked): `cohort.scorer_without_cohort:{scorer}`, a prediction type without measurements in
one group is left out (`drop_scorer`, `decided_by: skill`); an `abort` recorded earlier still stops.

Recording `replace:{curie}` also writes the new CURIE into `study.yaml` (a comment at its end keeps the
original and the decision), so `study.yaml` stays the request the study ran on. After `replace:` or
`add:`, run `fetch tracks` (it picks up the new CURIE) and `cohorts` again. A cohort without any
`RNA_SEQ` track cannot be fixed by dropping a scorer; when no descendant is listed, run `fetch
ols-term` for its CURIEs, or choose other CURIEs. Show `processed/tissue_cohorts.md` to the researcher.

## Steps 9-10. Saturation, AVI and scores

```bash
uv run SKILL/scripts/expression.py fetch avi --study-dir D
uv run SKILL/scripts/expression.py fetch scores --study-dir D
```

- `avi` writes `raw/avi_region.parquet` for the window of `processed/locus.json`, every SNV kept
  (`--min_phred 0`; the file is the universe of variants) and `--max_window_size` set to
  `max(1000, window length)` (the script refuses wider windows otherwise).
- `scores` writes `raw/scores_{S}.tsv` for each scorer kept by `cohorts`, with every kept CURIE of
  both cohorts, `--gene` only for gene-centric scorers (`RNA_SEQ`, `RNA_SEQ_ACTIVE`,
  `POLYADENYLATION`, `SPLICE_*`), and the key from the project's `.env` (the K-Dense script reads only the
  environment; the fetch passes the file to `uv`, so the key never reaches a command line).
- Capture both together: scores of another window or pull do not line up with the AVI table.

## Step 11. Variants and Tier 1

```bash
uv run SKILL/scripts/expression.py variants --study-dir D
uv run SKILL/scripts/expression.py master --study-dir D
```

`variants` leaves out, by rule, positions where the Atlas REF differs from the Ensembl sequence
(`variants.ref_mismatch:{pos}`) or where SNVs have no AVI row (`variants.missing_in_atlas:{pos}`), and
records each with `decided_by: skill`. Above 1% of the window it asks first
(`variants.exclusion_limit`: "Pause to capture the data again", recommended, → `abort`; "Leave these
positions out and continue" → `exclude_all`): so many disagreements usually mean the AVI export and
the sequence describe different windows, so check that first (steps 3 and 9). `master` exits 1 when a QC gate fails and lists the variants in
`processed/qc_report.json` (`failed_variants`); the usual cause is raw files from different windows
or pulls. Find the cause rather than continuing, because Tier 2 would inherit it.

## Step 12. Contrast (Tier 2)

```bash
uv run SKILL/scripts/expression.py contrast --study-dir D
```

The contrast plan feeds the fixed triage of `report.md` and `reproduce`; the analysis the researcher
steers is the dashboard (step 14), where groups, assays, sources and filters are chosen. So nothing is
asked here: without `D/contrasts.yaml` the stage generates it from the groups (one contrast per
prediction type with tracks in both) and records `contrast.plan` with `decided_by: skill`; a
`contrasts.yaml` already in the study is used as it is. Mixed assays or sources and a group with one
tissue are written as notes to `processed/tissue_contrast_meta.json` (`plan_notes`) and the report.
`contrast --propose` writes `processed/contrasts_proposed.yaml` for a maintainer who wants to edit the
plan (`track_filter`, aggregation) before running it.

## Step 13. Interpretation and report

```bash
uv run SKILL/scripts/expression.py fetch links --study-dir D
uv run SKILL/scripts/expression.py fetch interpret --study-dir D --ontology-id UBERON:AAAAAAA --top 5
uv run SKILL/scripts/expression.py report --study-dir D
```

- `links` writes `raw/atlas_links.md` (Atlas deep links of `processed/top_variants.json`).
- `interpret` is optional and uses the model API (about 20 s per variant): ask the researcher how many
  of the top variants and which target CURIE. For each variant it writes `raw/interpret_<variant>.txt`
  (`interpret_splicing.py`) and `figures/variants/<variant>/` (`visualize_variant_effects.py`, tissue
  label `--tissue`, default the target label), and lists the variants in
  `processed/interpret_variants.txt`. `--variant chr:pos:REF>ALT` (repeatable) replaces the top list.
- `report` reads the manifest, so run it last. Read `report.md` back before handing it over: every
  claim should be a predicted effect with its tissue and scorer, the skills section should match the
  manifest, and nothing should read as a diagnosis (`interpretation.md`).

## Step 14. Dashboard

```bash
uv run SKILL/scripts/expression.py dashboard --study-dir D
uv run SKILL/scripts/expression.py summarize --study-dir D --question repression --top 10
uv run SKILL/scripts/expression.py dashboard --study-dir D --export-html D/figures/dashboard.html
uv run SKILL/scripts/expression.py views --study-dir D
uv run SKILL/scripts/expression.py summarize --study-dir D --view "<saved view name>"
uv run SKILL/scripts/expression.py figures --study-dir D --view "<saved view name>"
uv run SKILL/scripts/expression.py report --study-dir D
```

`dashboard` serves on port 2718, or on a free port when 2718 is taken, and prints the URL to give the
researcher. marimo runs in an ephemeral `uv` environment (the first run needs network access) and
writes its session folder under `D/.workbench/`, never inside the skill. Open the link once, after the URL
is printed: two pages opened while it starts can break the first session (a matplotlib import error in the
server output; reload the page). If you restart it, tell the
researcher to reload the page: an old page is refused by the new server ("Unauthorized").

The dashboard writes to the study only when the researcher saves a view (*Visões salvas*), as
`analyses/<name>.yaml`; `views` lists them. To explain what the researcher sees, run `summarize --view <name>`
on their saved view, or `summarize` with the same choices (`--target`/`--control` CURIEs, `--question`, `--min-change-pct`, `--max-control-pct`,
`--min-sigma`, `--ref-floor`, `--floor-applies`, `--assay`, `--source`, `--domain`, `--one-alt`,
`--y-metric`, `--scorer`): it runs the dashboard's
computation and prints the headline, where each tissue's measurements come from, the noise between
measurements of the same tissue, the filter steps and the top variants. `references/dashboard.md`
describes every part of the page in the words the page uses, what to say when presenting it, and each
field of *Avançado*.

For each saved view wanted in the report, `figures --view <name>` draws it as PNG in
`figures/dashboard/<name>/` with `choices.json`, and `report` run again adds one subsection per view to
*Análise do pesquisador (dashboard)*: the name and note, the choices, the result, the filter steps and the
figures. Views hold only choices among the skill's parameters; a file with anything else is refused.

The export is the opening view (the study's groups, no filter), not what the researcher chose on the page.
`--edit` is for maintainers: it opens the marimo editor on the skill's own `assets/dashboard.py`, so
changes land in the skill, not in the study; never use it for a researcher's session. `--legacy` opens
the previous triage workbench.

Reproducing the reference (`expression.py reproduce --reference uchl1`): after a pass it prints a
"Headline results" block and a `Summary:` line with the total number of tracks (22). Report both.

## When a command fails

| Message | Meaning | Do |
|---|---|---|
| `N decision(s) needed: run study.py pending.` | an open question | Conventions: show, ask, `decide`, run again |
| `Missing decision '{key}'. Run Workflow step {N} …` | a stage needs an earlier answer | go back to step N |
| `invalid decision key(s) …` (`decide`) | a key with spaces or another shape; nothing was recorded | copy the key from `study.py pending` |
| `… not found …; run Workflow step N` | a raw or processed input is missing | run that step's fetch |
| `Dependency skill '…' not found` (`fetch`) | the skill is not installed where `check` looks | `expression.py check`; install only with the researcher's approval |
| `command failed with exit N: …` (`fetch`) | the dependency script failed; its last output follows | read it; network, key or an argument (e.g. an unknown gene symbol) |
| `FAIL file=… field=… skill=…` / `DependencyOutputError` | a dependency changed its output format | stop; compare with `dependency_outputs.md`; tell the researcher |
| `ALPHAGENOME_API_KEY is not set in the project's .env` | no key in the project | the `alphagenome-skills-setup` skill |
| `Warning: … .env holds keys but git does not ignore it` | a workspace `.env` could be committed | add `.env` to `.gitignore` now |
| `QC failed for N variant(s)` (`master`) | raw files disagree | check that steps 3, 9 and 10 used the same window and pull |
| An OLS or Ensembl fetch that never ends | a stale `polite-http` reservation | the fetch resets the skill's own lock folder; a call made outside the skill uses the system temp folder |
| `'X' is not an answer to 'K'` (`decide`) | the value is not one of the codes the question offered | record one of the codes it lists (`pending` shows label and code) |
| `Reference 'uchl1': … mismatch` (`reproduce`) | a code or data change moved a reference number | read the stage/metric table it prints before changing anything |
