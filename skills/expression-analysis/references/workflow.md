# Workflow reference: expression-analysis

The complete commands of the Workflow steps of `SKILL.md`. Formats and fields of every raw file are
in `dependency_outputs.md`; how to read the results is in `interpretation.md`.

## Contents

- [Progress checklist](#progress-checklist)
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
- [Common mistakes](#common-mistakes)

## Progress checklist

Copy this checklist into the conversation and tick it as you go:

```
expression-analysis progress:
- [ ] Step 0: check passes (anything missing set up by alphagenome-skills-setup); explore run; summary.md shown
- [ ] Step 1: window, target, control, labels, scorers and question given by the researcher (early answers matched and shown; the rest asked together; none filled from the recommendation unless explicitly delegated); init --from-exploration with --said; cohort decisions recorded
- [ ] Requests made ahead of their step listed, and each done at its step
- [ ] Step 4: locus run; locus.* answered by the researcher; Figure 1 shown
- [ ] Steps 6-8: fetch tracks; cohorts run; cohort.* answered (study.yaml updated on replace)
- [ ] Steps 9-10: fetch avi and fetch scores, same window and pull
- [ ] Step 11: variants and master run; positions left out by rule noted; QC gates passed
- [ ] Step 12: contrast run
- [ ] Step 13: atlas links (and interpretations) fetched; report.md written and read back
- [ ] Step 14: dashboard opened and presented (choices, how to show, filters, Avançado); questions answered from `summarize`; the researcher's view drawn with `figures` and `report` run again
- [ ] Final answer: researcher decisions (delegated ones named), the skill's rules and agent choices in separate lists; requests still to do
- [ ] Every number explained with the glossary's name, its unit and status, from a command; levels given next to large percentages; open points said as open
```

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
`processed/pending_decisions.json` and exits 3 with `N decision(s) needed: run study.py pending.`
(exit 1 is an error, 0 done, 130 interrupted).
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
uv run SKILL/scripts/study.py decide --study-dir D --key "locus.exon1_overlap:FEATURE" --value keep_full_unit --step 4 --said "<the researcher's words>" --note "Question as shown; the researcher's reason"
```

`--value` is parsed as YAML, so a list is written `"[UBERON:0000955, UBERON:0002107]"`. Quote keys
(they contain `:`). Several keys with the same answer go in one call, one `--key` each. Then run the
stage again. `decide` copies the question's `options` and `recommended` into `decisions.yaml`, so the
record shows what was offered next to what was chosen; for a question you ask yourself, pass what you
offered with `--option` (once per option) and `--recommended`. `--said` (required) quotes the
researcher's words the answer comes from, as they wrote them; the report shows them next to the answer. If
the researcher explicitly told you to choose this question, record the recommendation, quote that
instruction in `--said` and add `--delegated` (a general "go ahead" is not delegation: see `conversation.md`). A key must look like
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
- `explore` asks nothing and writes `OUT/exploration/` (see `formats.md`, *Folder layout*). Options:
  `--scorers` (default `RNA_SEQ CAGE DNASE ATAC`), `--flank-bp` (Regulatory Build search around the
  gene, default 2000, the study default), `--floor` (share of the gene's maximum below which a term is
  flagged, default 0.05). Re-running it reuses what was captured, so adding a
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
| 4. Terms | per tissue text: OLS candidates (UBERON, CL, EFO, CLO) with tracks per scorer, data source and RNA-seq protocol; catalogue terms whose name contains the text (subtypes, cell lines); served descendants of a term without RNA-seq tracks | a term without tracks cannot be contrasted; one track leaves nothing to compare it with (two tracks of a term differ in assay or source: they are not replicates) |
| 5. Baseline | `RNA_SEQ_ACTIVE` of the gene per tissue, relative to its maximum, with a warning below the floor | where the gene's predicted level is very low, a large % can be a minimal absolute change |
| 6. Scorers | which recommended terms each scorer serves | a scorer without tracks in one group is dropped later |
| 7. Limits | SNVs only (indels through the model API, `score_variants.py` of the `alphagenome` skill), window scope | what the study cannot answer |

**What to ask, and when.** The gene is its HGNC symbol (`explore --gene`). When the researcher names a
disease, a protein or an alias ("PKU", "phenylalanine hydroxylase"), say which gene you take it to be (*PAH*)
and confirm before exploring; a disease with several genes is a question, not a choice for you. Write the
symbol in italics (`glossary.md`, *Writing conventions*).

Ask the researcher only for the gene and the tissues or cell types in their own words, for both sides: where
they want to know whether the variants change expression (target) and what to compare it with (control), or
that they want the exploration to show candidates first. Ask even when the request names some tissues: a
request that names only the target, or names tissues only as context, would be explored again once the groups
come. Then `expression.py explore`. If the first request already brings other answers (a window, the groups,
the question), keep them for step 1 and say so; do not act on them yet. Tissues that reach you after the
exploration (new groups at step 1, a redefinition) go through `explore` again with every text (the other
captures are reused), and the researcher sees the terms and the gene's predicted expression level found for
them before `init` or `decide`.

Show `summary.md`: windows with size, SNVs and cost; candidate terms with tracks per scorer, source and assay;
served descendants of terms without tracks; the gene's predicted expression level per tissue, with a warning
where it is below the floor. Explain that level with *The exploration's expression level* in
`interpretation.md` (a prediction of the probe SNV's larger allele, no official unit, not TPM, comparable only
within one kind of track). Do not ask for the study's choices before this page is shown, because they depend
on it. Asking before the data are in hand makes the researcher choose blind, and every later surprise (a
term without tracks, a tissue where the gene's predicted level is very low) becomes another interruption.

**When Ensembl does not answer**, `explore` still writes the page, without the regulatory elements and the
expression level, and opens it with *Capturas que faltam*. Show it as it is, say what is missing and that the
study waits for it (`init` refuses), and run the same `explore` again later; nothing already captured is
fetched twice. Do not fill the gap with another source or your own code.

### Changing the groups after the results

When the researcher wants other target or control tissues (for example, controls where the gene is also
expressed), that is a redefinition of the study, and it follows the same path as the definition. Choose from
the exploration, never from ad-hoc code; never rank tissues with your own code:

1. Show `OUT/exploration/activity.md`. It lists every catalogue term with the gene's predicted expression level
   (as a share of the gene's maximum), kind of sample, assays, sources and number of tracks. Without it (an
   older exploration), run `expression.py explore --out OUT --activity-only`; it captures nothing. For
   tissues the researcher names in their own words, run `explore` again with those `--tissue` texts: it
   reuses the captures and only searches the new terms.
2. Discuss the choice with the researcher. Prefer terms where the gene is expressed. Keep the kinds of
   sample and the assays balanced with the target. Say what the new comparison asks, which may differ from
   the original question.
3. Record it: `study.py decide --key cohort.control` (or `cohort.target`) with the new list, the terms
   offered (`--option`) and the recommended one, the group's new name (`--group-label`), the researcher's words (`--said`) and a note saying
   why the study was redefined. `decide` writes the new terms and name into `study.yaml`; do not edit it.
   Then `fetch tracks`, `cohorts`, `fetch scores`, `variants`, `master`, `contrast`, `report`, as in steps
   7-13. The contrast plan the skill generated follows the new groups by itself; do not delete it.
   Tissues taken out of a group can stay visible as context tissues on the dashboard. `decide` keeps the
   earlier answer in `previous` (value, note, date, options), and `report.md` lists it under
   *Redefinições*; never edit `decisions.yaml` to replace an answer.
4. Close it in the report. When the researcher settles on the redefined analysis, run `figures` with its
   choices and `report` again, and read back *Redefinições* and the dashboard section: what changed, why,
   and what the new comparison shows.

## Step 1. Study from the exploration

Ask the researcher, with the words of the summary page, for the five answers the study needs: which
window (section 3, by its name; you pass its id to `--window`), which terms form the target and the control
groups (section 4) with a short name for each group, which kinds of prediction to compare (section 6; gene
expression, `RNA_SEQ`, is always included) and what they look for in the target (`--question`: `repression`
for a decrease, `increase`, or `either`). Ask for the missing ones together; an answer the researcher did not
give is asked, never taken from the recommendation. `init` stops and lists any that is missing. If the
researcher answers only some (for example, only the groups), ask the rest before going on; never fill an
answer from the recommendation, even when the researcher seems in a hurry. Answers given before this step
count once they are matched to the options and the researcher has seen the match (`conversation.md`, *The
researcher talks freely*). Assays, sources, the score (raw or quantile) and how a group's tissues are combined
are not asked: everything is captured and the researcher changes them on the dashboard (*Avançado*).
`study.py init --from-exploration` refuses to run without `--window`, `--target`, `--control`, `--scorers`,
`--question` and `--said` (their words, joined with ` / ` when they came in several messages; `--delegated
window` etc. for an answer they explicitly left to you). Then:

```bash
uv run SKILL/scripts/study.py init --study-dir D --from-exploration OUT --window regulatory_unit --target CL:AAAAAAA --target UBERON:BBBBBBB --control UBERON:CCCCCCC --target-label "<label>" --control-label "<label>" --scorers RNA_SEQ CAGE --question repression --said "<the researcher's words / joined across messages>" --origin-note "Requested by <researcher> on <date>: <question in one line>"
uv run SKILL/scripts/study.py decide --study-dir D --key cohort.target --value "[CL:AAAAAAA, UBERON:BBBBBBB]" --option "CL:AAAAAAA (<label>)" --option "UBERON:BBBBBBB (<label>)" --option "UBERON:DDDDDDD (<label>)" --step 1 --said "<the researcher's words>" --note "Chosen by the researcher from the exploration; why"
uv run SKILL/scripts/study.py decide --study-dir D --key cohort.control --value "[UBERON:CCCCCCC]" --option "UBERON:CCCCCCC (<label>)" --step 1 --said "<the researcher's words>" --note "Chosen by the researcher; why"
uv run SKILL/scripts/study.py validate --study-dir D
```

The groups take no `--recommended` (`decide` refuses it): the page recommends one term per text the researcher
gave, not a group. Groups or tissues given only now, after the exploration, go through `explore` again first
(every text; the other captures are reused), so the researcher sees the terms and the gene's predicted expression level found.

`init --from-exploration` records `study.window`, `study.scorers` and `analysis.question` in
`decisions.yaml` with the options offered (the report lists them with the cohort answers; the dashboard opens
on that question), writes `study.yaml` (window policy from the option: `regulatory_unit`
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

If the researcher's text finds nothing, another spelling of the same name ("MCF7" for "MCF-7") is fine: say so
and put it in the `--note` of the cohort decision. A different term than the one they named is their choice:
ask.

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

Rule (not asked): `cohort.scorer_without_cohort:{scorer}`, a prediction type without tracks in
one group is left out (`drop_scorer`, `decided_by: skill`); an `abort` recorded earlier still stops.
Rule (not asked): `cohort.term_without_ols:{curie}`, a study term the Atlas catalogue serves but the ontology
no longer knows (cell lines renamed or retired, e.g. some EFO cell lines) enters with the catalogue's sample
name (`atlas_catalogue`).

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

Read `dashboard.md` first. With the researcher: run `expression.py dashboard` and give them the link (tell
them to reload the page if you restart it), with the short presentation of that file: the three choices, the
row of how to show, where filters and *Avançado* are. When they open *Avançado* or ask about a field, explain
it with the field table there, including what the page does without saying (empty = all, ignored filters,
display limits). When they ask about what they see, read it with `expression.py view` and run
`expression.py summarize --from-dashboard`, then explain from its numbers; do not describe results you have
not computed. To show them a view, say what you will change and why, then send it with
`expression.py view <choices> --message "…"`: the open page follows, and they can change it further. When they
settle on a view, run `expression.py figures --from-dashboard`, then `expression.py report` again, so
`report.md` shows that view's figures with the choices behind them; the dashboard itself saves no analysis.

```bash
uv run SKILL/scripts/expression.py dashboard --study-dir D
uv run SKILL/scripts/expression.py view --study-dir D
uv run SKILL/scripts/expression.py summarize --study-dir D --from-dashboard --top 10
uv run SKILL/scripts/expression.py view --study-dir D --min-sigma 1 --message "<what this view is>"
uv run SKILL/scripts/expression.py figures --study-dir D --from-dashboard
uv run SKILL/scripts/expression.py summarize --study-dir D --question repression --top 10
uv run SKILL/scripts/expression.py dashboard --study-dir D --export-html D/figures/dashboard.html
uv run SKILL/scripts/expression.py figures --study-dir D --question repression --target CL:AAAAAAA --control UBERON:BBBBBBB --top-k 100 --show all
uv run SKILL/scripts/expression.py report --study-dir D
```

`dashboard` serves on port 2718, or on a free port when 2718 is taken, and prints the URL to give the
researcher. marimo runs in an ephemeral `uv` environment (the first run needs network access) and
writes its session folder under `D/.workbench/`, never inside the skill. Open the link once, after the URL
is printed: two pages opened while it starts can break the first session (a matplotlib import error in the
server output; reload the page). If you restart it, tell the
researcher to reload the page: an old page is refused by the new server ("Unauthorized").

While it is open, the page and you share two files in `D/.dashboard/`: the page writes what the researcher
has on screen (`researcher_view.json`) and follows the view you send (`agent_view.json`). `view` without
choices prints the screen; with choices, the page moves to that view without reloading (only what you pass
changes; `--reset` starts from the opening view). `summarize` and `figures` take `--from-dashboard` to use
the choices on screen. The page must be open in a browser for the link to work; `references/dashboard.md`
(*Link with the page*) has the rules. Nothing else is written to the study.

To explain what the researcher sees, run `summarize --from-dashboard`, or, with the page closed, with
the same choices (`--target`/`--control` CURIEs, `--question`, `--min-change-pct`, `--max-control-pct`,
`--min-sigma`, `--ref-floor`, `--floor-applies`, `--assay`, `--source`, `--domain`, `--one-alt`,
`--y-metric`, `--scorer`): it runs the dashboard's
computation and prints the headline, where each tissue's tracks come from, the disagreement between
tracks of the same tissue, the filter steps and the top variants. `references/dashboard.md`
describes every part of the page in the words the page uses, what to say when presenting it, and each
field of *Avançado*.

When the researcher settles on a view, `figures` (same choices, plus `--top-k`, `--show`, `--color-by`,
`--highlight`) draws it as PNG in `figures/dashboard/` with `choices.json`, and `report` run again adds
the section *Análise do pesquisador (dashboard)* to `report.md`: the choices, the result, the filter steps
and the figures. A new `figures` run replaces the previous view.

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
| `--said is required` (`decide`, `init`) | no words of the researcher were passed for the answer | quote what they wrote; if they have not answered, ask (never record an answer of your own) |
| `Reference 'uchl1': … mismatch` (`reproduce`) | a code or data change moved a reference number | read the stage/metric table it prints before changing anything |

## Common mistakes

- **Asking for tissues, scorers and window before the exploration.** The researcher then chooses
  without knowing which terms have tracks or where the gene is expressed; run `explore` first.
- **Typing dependency commands by hand.** Use `expression.py fetch`: it pins the packages, writes the
  file name the stages read and records the command; a hand-typed command is easy to misname
  (`ols_t_cell.json` instead of `ols_search_t_cell.json` breaks `cohorts`) or to paraphrase.
- **Annotating exons from Ensembl REST instead of `raw/gtf.json`.** The Atlas scores were computed
  on GENCODE v46; Ensembl is used only for the Regulatory Build, the sequence and IDs.
- **Editing `raw/`.** The `gtf` starts are 0-based as written; every stage normalises them once.
- **Mixing scores of different Atlas pulls or windows.** `avi_region.parquet` and every
  `scores_{S}.tsv` describe one window captured together; re-capture all of them or none.
- **Reading a predicted effect in a tissue where the gene's predicted level is very low as tissue-specific.**
  A large % there can be a minimal absolute change; check the exploration's activity table and give the level.
- **Defining a score from memory, or naming a quantity your own way.** What each number is, its unit and its
  status are in `glossary.md`; use the Portuguese name written there, and say "not established" for
  an entry marked undecided. The full definition of a score is in the `alphagenome` dependency skill
  (`interpretation.md`, *Where each score is defined*).
- **Explaining the dashboard from memory.** Its numbers depend on the researcher's current choices; read the
  screen with `view`, run `summarize --from-dashboard` and quote it, and use the wording of `dashboard.md`.
- **Moving the researcher's page without saying so, or assuming a view arrived.** Say what you will change and
  why before `view`, give it a `--message`, and read the screen back; with the page closed the view only waits.
- **Presenting a filter as a clinical threshold.** A minimum % change, a multiple of the between-track disagreement, AVI
  Phred 20 or a large Δ rank predicted effects; none is a pathogenicity cut-off.
- **Writing your own scripts in the study folder.** Ad-hoc code in `D/` is not in the manifest and is
  mistaken for the skill's output; use the commands here, and keep out-of-scope work outside `D/`
  (`conversation.md`, *Out of scope*).
- **Choosing for the researcher, or editing `decisions.yaml` by hand.** Record answers with
  `study.py decide`, which also clears the pending question and keeps the timestamps.
- **Showing codes instead of the question.** `keep_full_unit`, `drop_curie` or `abort` mean nothing to
  the researcher; ask with the labels and effects, in their language, and keep the identifiers.
- **Installing a dependency, writing a `.licenses/` file or handling a key from this skill.** That is
  `alphagenome-skills-setup`'s job; doing it here repeats its questions.
- **Printing or committing a `.env`, or putting the key elsewhere.** The key lives only in the project's
  `.env`, listed in `.gitignore`; the fetch commands load it and it never reaches the screen, the manifest or git.
