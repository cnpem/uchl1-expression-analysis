---
name: expression-analysis
description: "Explores a human gene with AlphaGenome before any study is fixed (gene structure, Regulatory Build, candidate windows with their cost, UBERON/CL/EFO terms with Atlas tracks per tissue, the gene's predicted baseline expression per tissue), then compares the predicted effects of every SNV in the chosen window between target and control tissues or cell types, with a Portuguese report and an interactive dashboard where the researcher picks groups, question and filters. Routes to the DeepMind science skills (GENCODE, OLS, Atlas AVI and track catalogue, Ensembl) and to the K-Dense AlphaGenome Atlas scores. Use it whenever someone asks which variants in a gene change expression in one tissue versus another, whether a gene is expressed in a tissue according to AlphaGenome, which tissues or windows can be studied for a gene, wants a saturation-mutagenesis or promoter/enhancer variant triage, or mentions tissue-specific regulatory effects with AlphaGenome, even without saying 'expression analysis'. Research use only."
metadata:
  version: 0.3.0
---

# expression-analysis

## Overview

For one human gene (GRCh38), this skill takes every possible SNV of the gene's regulatory window
(by default the promoter/enhancer unit that holds the TSS, through the end of exon 1), collects the
AlphaGenome Atlas predictions for each of them, and asks: which variants change the predicted
signal (RNA-seq, CAGE, DNase, ATAC) in the target tissues, and do they do so differently in the
control tissues? The researcher answers it in an interactive marimo dashboard, choosing the groups,
the question and every filter on the data already captured. Behind it: a Tier 1 table (one row per
SNV), a fixed triage (Tier 2 contrast and ranked top list with Atlas links) and a Portuguese
`report.md`.

It starts with an **exploration** (`expression.py explore`): cheap, read-only queries that show the
researcher, in a Portuguese `summary.md`, which windows are possible and what each costs, which
ontology terms match the tissues they named and whether the Atlas has tracks for them, and whether
the gene is expressed there at all. The researcher chooses from that page; only then is the study
created. Asking before the data are in hand makes the researcher choose blind, and every later
surprise (a term without tracks, a tissue where the gene is silent) becomes another interruption.

How the work is split, and why:

- **Dependency skills fetch, this skill analyses.** Gene structure, ontology terms, track
  catalogue and scores come from the DeepMind science skills and the K-Dense `alphagenome` skill.
  `expression.py fetch <kind>` runs each of them with pinned package versions, writes the output
  unedited under the name the next stage reads, and records the exact command in `manifest.json`.
  The stages read only the study folder, so a study can be re-run offline and every number traces
  back to a recorded command. You choose the parameters; the script runs and records them, so no
  shell syntax has to be adapted between PowerShell, bash and zsh.
- **The researcher decides, the agent asks.** When the data leave a choice open that changes which
  variants and tissues the conclusions rest on (which transcript, how far the window goes, a tissue
  without Atlas measurements), the stage stops with exit 1 and `N decision(s) needed`, and writes the
  question with its options, their effects and a recommendation. Ask it in the researcher's language
  and record their answer; do not pick one yourself. Technical checks and choices that do not change
  the steered analysis are applied by fixed rules and recorded as the skill's (`decided_by: skill`).
  The identifiers (CURIEs, transcript and region IDs, coordinates) always stay next to the names, so
  every choice can be traced.
- **Research use only.** Everything here is a predicted molecular effect used to prioritise
  experiments. It is not a measurement, not a clinical classification and never a diagnosis.

Read `references/workflow.md` for the complete commands and the questions of each step,
`references/interpretation.md` before reporting results, `references/dashboard.md` before opening the
dashboard or answering a question about it, `references/formats.md` for the fields of the study files
(exploration, questions, decisions, contrast plan), and `references/dependency_outputs.md` when a
dependency output fails validation.

## Prerequisites

The **`alphagenome-skills-setup`** skill installs and configures everything this skill needs, with the
researcher's approval for each item; this skill installs nothing and asks nothing about installation.
Run `expression.py check`. When it reports something missing, hand the setup to that skill. If it is not
installed, `check` prints how to install it; ask the researcher first (the repository is private for now).
Do not install dependency skills, accept terms of use or handle keys from this skill: doing it twice,
from two skills, asks the researcher the same questions twice.

What `check` tests:

1. **`uv`** on PATH. Every script runs with `uv run`, which installs its pinned PEP 723 dependencies into a
   cached environment, so nothing is installed into the system Python.
2. **The dependency skills** of the next section, looked for next to this skill, in the workspace's
   `.claude/skills/` and `.agents/skills/`, and in `~/.claude/skills/` and `~/.agents/skills/`.
3. **Terms of use** of `embl-ebi-ols`, `ensembl-database` and `alphagenome-single-variant-analysis`,
   recorded in the workspace's `.licenses/` once the researcher accepted them.
4. **`ALPHAGENOME_API_KEY` in the project's `.env`**, the only place the key is read from (not `~/.env`,
   not environment variables), so each study runs with its own project's key. The key is personal (the
   AlphaGenome terms forbid sharing it, even inside the institution), so `.env` must be in `.gitignore`;
   `check` warns when it is not. The fetch commands hand that file to the scripts that need it. Do not
   source, print or copy a `.env`.
5. **An output folder** chosen by the researcher, outside this skill's folder (the skill folder is
   read-only; the exploration and the study hold the researcher's data and decisions).

## Dependencies

Cited by the `name` of each skill. Versions tested are in `references/dependency_outputs.md`;
the pinned package lists are in `scripts/deps.py`.

| Skill | Origin | Steps | Used for |
|---|---|---|---|
| `uv` | DeepMind | all | running every script with `uv run` |
| `credentials` | DeepMind | prerequisites | required by the DeepMind skills; the key itself is set by `alphagenome-skills-setup` in the project's `.env` |
| `alphagenome-variant-impact-score` | DeepMind | 0, 2, 7, 9 | `gtf` (GENCODE v46 structure), `metadata --tracks` (track catalogue), `region` (AVI of every SNV) |
| `embl-ebi-ols` | DeepMind | 0, 5, 6 | `search_ols.py` (UBERON/CL/EFO/CLO terms), `get_term.py` (descendants) |
| `ensembl-database` | DeepMind | 0, 3 | `get-sequence` (reference bases: REF check, start codon) |
| `alphagenome` | K-Dense | 0, 10 | `atlas_query.py tracks` (sources, assays), `scores` (per-track scores; `RNA_SEQ_ACTIVE` baseline) |
| `alphagenome-single-variant-analysis` | DeepMind | 13 | `interpret_splicing.py`, `visualize_variant_effects.py` on the top variants |
| `alphagenome-atlas-website-links` | DeepMind | 13 | `table --input` (Atlas deep links for the top variants) |

**Installation belongs to `alphagenome-skills-setup`.** Installing a skill brings third-party code and,
for several of these, terms of use; that skill asks the researcher about each one. No script of this
skill installs anything.

## Quick Start

`SKILL` stands for this skill's folder, `OUT` for the folder the researcher chose.

Reproduce the bundled reference study (offline, a few seconds; six stages compared with
`assets/reference/uchl1/expected.json`), then open its dashboard:

```bash
uv run SKILL/scripts/expression.py reproduce --reference uchl1 --keep OUT/uchl1-demo
uv run SKILL/scripts/expression.py dashboard --study-dir OUT/uchl1-demo
```

Report the "Headline results" block that `reproduce` prints after a pass, and the total number of
tracks from its `Summary:` line (22 for UCHL1); a reproduction that is not reported cannot be checked.

A new gene starts with the exploration, in the researcher's words for the tissues:

```bash
uv run SKILL/scripts/expression.py check
uv run SKILL/scripts/expression.py explore --gene MYGENE --tissue "<tissue in their words>" --tissue "<another>" --out OUT
```

Then show `OUT/exploration/summary.md` and create the study from the researcher's choices
(Workflow step 1). Every command works unchanged in PowerShell, bash and zsh.

## Utility Scripts

All in `scripts/`; `uv run SKILL/scripts/<script> <subcommand> --help` lists the options. Every
stage and fetch records itself in `manifest.json` when it succeeds, so there is no separate
recording command to run.

| Command | Does |
|---|---|
| `expression.py check [--json]` | Checks `uv`, dependency skills, terms of use, `ALPHAGENOME_API_KEY` in the project's `.env`; names what is missing and hands it to `alphagenome-skills-setup`; installs nothing |
| `expression.py explore --gene G --out OUT --tissue T…` | Step 0: read-only exploration → `OUT/exploration/summary.md` (Portuguese), `exploration.json`, `raw/`, `atlas/`, `manifest.json` |
| `study.py init --from-exploration OUT --study-dir D …` | Step 1: `study.yaml` from the chosen window, CURIEs, labels and scorers; copies the exploration captures into `raw/` with their manifest steps |
| `expression.py fetch <kind> --study-dir D [--dry-run]` | Steps 2-13 captures: `gtf`, `regulatory`, `sequence`, `ols`, `ols-term`, `tracks [--catalogue]`, `avi`, `scores`, `links`, `interpret` |
| `expression.py locus` / `cohorts` / `variants` / `master` | Steps 4, 8, 11: window and Figure 1; cohorts and track registry; SNV table; Tier 1 table with two QC gates |
| `expression.py contrast [--propose]` | Step 12: the triage's contrast (plan generated from the groups when absent) → tissue contrast, top variants; `--propose` only writes the plan for a maintainer |
| `expression.py report` | Step 13: `report.md` (study, decisions, results, skills used, limits, citations) |
| `expression.py dashboard [--port N] [--export-html P] [--edit] [--legacy]` | Step 14: analysis dashboard (a free port if N is busy; prints the URL), a static HTML export, the marimo editor for maintainers (`--edit`, changes the skill's `assets/dashboard.py`) or the previous workbench (`--legacy`) |
| `expression.py views --study-dir D [--print NAME \| --save NAME choices…]` | The dashboard's saved views (`analyses/*.yaml`): list (who saved each, name for `--view`), one in full, or save the agent's analysis as a view the researcher opens on the page |
| `expression.py summarize --study-dir D [--view NAME \| choices…]` | The dashboard's view as text for the same choices (groups, question, filters): headline, measurements, noise, filter steps, top variants. Writes nothing |
| `expression.py figures --study-dir D [--view NAME \| choices… --top-k K --show kept\|all --highlight V…]` | The dashboard's charts for a saved view (or given choices) as PNG in `figures/dashboard/<view>/`, with `choices.json`; `report` shows each view drawn |
| `expression.py regions [--only NAME]` | Prints the regulatory search region, the window and `max_window` (the fetch commands use them already) |
| `expression.py reproduce --reference uchl1 [--until S] [--keep DIR]` | Re-runs the stages on the bundled reference and compares with `expected.json` |
| `study.py validate / pending [--json] / decide / require / window` | Study file, researcher decisions, open questions, `--print-max-window` |
| `manifest.py show` | Prints the manifest (`manifest.py add` remains for a step run outside the skill) |
| `models.py validate --dir D/raw` | Checks every raw file against the fields the dependencies write |

## Workflow

**Decisions (every stage).** When a stage exits 1 with `N decision(s) needed: run study.py
pending.`, run `study.py pending --study-dir D` (or `--json`) and first do what its "For the agent" lines
say. Ask one question at a time (or questions of the same kind together), in the researcher's language:
the title, the question, why it matters, the names with their identifiers, each option's label and
what happens, and which one is recommended and why. With more than about six options, show the first
ones and offer the rest on request. Never show codes or keys. Record the answer with the chosen
option's code: `study.py decide --study-dir D --key K --value CODE --step N --note "…"` (one `--key`
per key when several share an answer), and run the stage again. The recommendation is the stage's
suggestion, not an answer. For a question you ask yourself (the choices of step 1), pass what you
offered with `--option` (repeated) and `--recommended`. When the researcher has delegated ("take the
recommended options"), record the recommendation and say so in the note.

**Final answer.** List, in the researcher's language: the decisions they made (question, options,
recommended, chosen, reason, identifiers); the rules the skill applied (`decided_by: skill` in
`decisions.yaml`: what was left out and why, summarising many of a kind, e.g. "5 positions left out:
chr…"); and the choices you made that no stage raised (e.g. a
search text or a scorer list), so the researcher can tell the three apart.

0. **Exploration.** `expression.py check`; anything missing goes to `alphagenome-skills-setup`.
   Ask the researcher only for the gene and the tissues or cell types in their own words, then
   `expression.py explore`. Show `summary.md`: windows with size, SNVs and cost; candidate terms
   with tracks per scorer, source and assay; served descendants of terms without tracks; the
   gene's predicted activity per tissue, with a warning where it is below the floor. Do not ask
   for the study's choices before this page is shown, because they depend on it.
   To change the groups after the results (a redefinition), go back to the exploration:
   `exploration/activity.md` lists every catalogue term with the gene's activity, kind of sample and assays
   (`explore --out OUT --activity-only` writes it for an older exploration). Then record the new groups with
   `decide` and recapture (`references/workflow.md`, *Changing the groups after the results*). Never rank
   tissues with your own code.
1. **Study.** Ask the researcher, with the words of the exploration page, which window (by its name
   on the page), which terms form the target and the control groups (with a short name for each group),
   and which kinds of prediction to compare (gene expression is always included). `study.py init --from-exploration`,
   then `decide` `cohort.target` and `cohort.control` with the terms offered (`--option`,
   `--recommended`), and `study.py validate`.
2. **Gene structure** (`fetch gtf`) and 3. **Regulatory Build and sequence** (`fetch regulatory`,
   `fetch sequence`): already copied from the exploration; run them only for a study without one.
4. **Locus and window.** `expression.py locus`; questions `locus.transcript`, `locus.feature_*`,
   `locus.no_feature_chain` (the rest of the window follows its definition, by rule). Show Figure 1.
5. **Terms.** Also from the exploration; `fetch ols --query "…"` for a tissue added later.
6. **Descendants** (`fetch ols-term --curie C`) for a CURIE without tracks or without RNA_SEQ tracks.
7. **Track catalogue.** `fetch tracks` (every study CURIE × scorer; zero-track pairs are recorded).
8. **Cohorts.** `expression.py cohorts`; questions `cohort.unserved:*`, `cohort.no_rna_seq:*`. A
   `replace:` answer also updates `study.yaml`; then `fetch tracks` for the new CURIE and run `cohorts`
   again. A prediction type without measurements in one group is left out by rule.
9. **Saturation and AVI** (`fetch avi`) and 10. **Scores per tissue** (`fetch scores`), together.
11. **Variants and Tier 1.** `expression.py variants` (positions the Atlas does not match are left out
    by rule), then `master` (two QC gates).
12. **Contrast.** `expression.py contrast` (the triage; no question).
13. **Interpretation and report.** `fetch links`; optionally `fetch interpret` on the top variants
    with a target CURIE the researcher chose; `expression.py report`; read `report.md` back.
14. **Dashboard.** Read `references/dashboard.md`. `expression.py dashboard` and give the researcher the
    link (tell them to reload the page if you restart it), with the short presentation of that file: the
    three choices, the row of how to show, where filters and *Avançado* are. When they open *Avançado* or
    ask about a field, explain it with the field table there, including what the page does without saying
    (empty = all, ignored filters, display limits). When they ask about what they see, run
    `expression.py summarize` with the same choices and explain from its numbers; do not describe results
    you have not computed; with a view they saved (*Visões salvas*), use `--view <name>` so you read exactly
    their choices. When you build a new analysis they want to see on the page, offer to save it with
    `expression.py views --save "<name>" --note "…"` (same choice options) and tell them to open it in
    *Visões salvas* (*Atualizar lista*, *Abrir*). For each saved view they want in the report, run
    `expression.py figures --view <name>`,
    then `expression.py report` again, so `report.md` shows each view's figures with its choices and note.
    The page saves only when they click *Salvar*.

Copy this checklist into the conversation and tick it as you go:

```
expression-analysis progress:
- [ ] Step 0: check passes (anything missing set up by alphagenome-skills-setup); explore run; summary.md shown
- [ ] Step 1: window, target, control, labels, scorers chosen by the researcher; init --from-exploration; cohort decisions recorded
- [ ] Step 4: locus run; locus.* answered by the researcher; Figure 1 shown
- [ ] Steps 6-8: fetch tracks; cohorts run; cohort.* answered (study.yaml updated on replace)
- [ ] Steps 9-10: fetch avi and fetch scores, same window and pull
- [ ] Step 11: variants and master run; positions left out by rule noted; QC gates passed
- [ ] Step 12: contrast run
- [ ] Step 13: atlas links (and interpretations) fetched; report.md written and read back
- [ ] Step 14: dashboard opened and presented (choices, how to show, filters, Avançado, saved views); questions answered from `summarize --view`; the saved views wanted in the report drawn with `figures --view` and `report` run again
- [ ] Final answer: researcher decisions, the skill's rules and agent choices in separate lists
```

## Scope

The skill answers one question: which variants of the window the model predicts to change the gene's
expression in the target tissues and not in the control. Questions beyond it, such as which
transcription factor or motif a variant touches, the mechanism, or what the literature says, are out of
scope. When the researcher asks one:

- Say first that it is outside the skill: nothing in the study checks it, and the answer is the agent's,
  not the model's.
- Work outside the study folder (e.g. `OUT/extra/<topic>/`), never in `D/`: the stages, the manifest and
  `report.md` must describe only what the skill ran.
- Record what was done there (`README.md`: question, sources, commands, coordinates and strand checked,
  what is uncertain), cite every source you actually opened, and never generalise beyond it.
- Never propose experiments, assays or validations; describe what the predictions and sources say.

## Rate Limiting

`regulatory_features.py` paces its own Ensembl requests (15/s, file lock in
`~/.cache/alphagenome-gene-skills/`). Every other request comes from a dependency script with its
own limits (Ensembl and OLS through `polite-http`; the Atlas client's `--workers`). The fetch
commands give `polite-http` its own lock folder and reset a reservation left far in the future by
another clock, which otherwise makes OLS and Ensembl calls sleep silently for hours. Run the fetches
one at a time rather than spreading them over parallel subagents. Model-API calls of step 13 take
about 20 s each; run them sequentially.

## Common Mistakes

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
- **Reading a predicted effect in a tissue where the gene is silent as tissue-specific.** Low
  baseline activity inflates predicted fold changes; check the exploration's activity table.
- **Explaining the dashboard from memory.** Its numbers depend on the researcher's current choices; run
  `summarize` with the same choices and quote it, and use the wording of `references/dashboard.md`.
- **Presenting a filter as a clinical threshold.** A minimum % change, a multiple of the noise, AVI
  Phred 20 or a large Δ rank predicted effects; none is a pathogenicity cut-off.
- **Writing your own scripts in the study folder.** Ad-hoc code in `D/` is not in the manifest and is
  mistaken for the skill's output; use the commands here, and keep out-of-scope work outside `D/` (*Scope*).
- **Choosing for the researcher, or editing `decisions.yaml` by hand.** Record answers with
  `study.py decide`, which also clears the pending question and keeps the timestamps.
- **Showing codes instead of the question.** `keep_full_unit`, `drop_curie` or `abort` mean nothing to
  the researcher; ask with the labels and effects, in their language, and keep the identifiers.
- **Installing a dependency, writing a `.licenses/` file or handling a key from this skill.** That is
  `alphagenome-skills-setup`'s job; doing it here repeats its questions.
- **Printing or committing a `.env`, or putting the key elsewhere.** The key lives only in the project's
  `.env`, listed in `.gitignore`; the fetch commands load it and it never reaches the screen, the manifest or git.

## Outputs

```
OUT/exploration/                     # step 0, before any study
  summary.md                         # Portuguese options page for the researcher
  activity.md                        # every catalogue term with the gene's predicted activity (choosing groups)
  exploration.json                   # the same options, read by study.py init --from-exploration
  raw/                               # captures with the study's file names (copied into the study)
  atlas/                             # atlas_tracks.tsv (sources, assays), scores_RNA_SEQ_ACTIVE_baseline.tsv
  manifest.json                      # every command, version and output
<study_dir>/                         # chosen by the researcher, never inside the skill
  study.yaml                         # the request (gene, cohorts, scorers, window policy, top_n)
  decisions.yaml                     # the researcher's answers and the skill's rules (references/formats.md)
  contrasts.yaml                     # the triage's contrast plan (step 12)
  manifest.json                      # provenance: every command, output, version and decision key
  raw/                               # dependency outputs, unedited
  processed/                         # this skill's outputs (locus, registry, master, contrast, top variants)
  figures/                           # fig1_genomic_architecture.html, dashboard.html (export), variants/ (step 13),
                                     # dashboard/<view>/ (step 14: each saved view drawn as PNG + choices.json, shown in report.md)
  analyses/                          # the dashboard's saved views (<name>.yaml), saved only when the researcher asks
  report.md                          # Portuguese report
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
  (UBERON, CL, EFO, CLO); cite them as their skills ask. `report.md` lists the skills actually used,
  read from `manifest.json`.
