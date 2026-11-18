# Changelog

Changes to the `expression-analysis` skill that matter to the people who use it and to the studies it
produces. Format based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/) as described in [Versioning](#versioning) below.

Each study records the skill version it ran with (`skill_version` in `manifest.json` and in each
step's `source_versions`; also shown in `report.md`). To know what a study's numbers rest on, find that version here.

## [Unreleased]

### Added

- **Term registry** (`assets/terms.yaml`, `scripts/terms.py`): one entry per quantity the skill shows, with the
  API variable it comes from, its name in Portuguese and English, definition, formula, unit, what it may be
  compared with, status (official, this skill's convention, undecided) and source.
  - `references/glossary.md` is generated from it; the agent reads it before explaining a number.
  - `expression.py terms` checks the registry, the glossary and the code and documentation for deprecated
    names; `reproduce` runs the same check.
  - `expression.py check` notes when the installed `alphagenome` skill differs from the version the
    definitions were checked against (`skills-lock.json`).

### Changed

- **`references/interpretation.md` rewritten** for the agent's explanations: how to explain a number (the
  glossary's name; what it is, relative to what, unit, status; values from a command; "not established" for
  an open point; levels next to a large percentage), each AlphaGenome score with its formula and unit, the
  exploration's expression level step by step, the dashboard's quantities with their formulas (percent,
  group value, between-track disagreement and its 1.4826 × median scale, REF reconstruction), the three
  thresholds, very low levels, the model's stated limits (tissue specificity, no uncertainty estimate), and
  answers to the questions researchers ask. The fixed triage is one section. SKILL.md: explaining a number
  follows it (step 0, step 14, final answer, checklist).
- Statements aligned with the open points: comparability of levels and of effects across assays, the
  quantile background and the Atlas model are marked not established where no source settles them; the
  exploration page no longer says that effects in a low-level target are unreliable without saying why.
- **Names on the pages** (wording only; no number changes). The exploration pages, the dashboard, `summarize`
  and the report now use the registry's names:
  - *nível de expressão previsto* for the gene's level (before: *atividade*);
  - *trilha* for an Atlas output of one sample, assay and source (before: *medição*);
  - *discordância entre trilhas* and *discordância típica* (before: *ruído* and σ). The filter steps read
    "alvo − controle ≥ N × a discordância típica"; the pages say the two tracks of a pair differ in assay or
    source and are not replicates;
  - *abaixo de 5% do máximo do gene* (before: *gene não expresso*, *praticamente desligado*);
  - *variante de sondagem*, *tipo de predição*, *pseudocontagem da pontuação RNA_SEQ*.
  The option `--min-sigma` and the columns `diff_sigma` and `sigma` keep their names.
- The exploration no longer says that a tissue without expression is a good control "para mostrar ausência
  de efeito": a term of very low level can dominate a group value.

## [0.4.0] - 2026-11-17

The researcher's decisions become part of the study, and the agent and the open dashboard can see and move
each other's view. This version changes the questions asked and the study files (`decisions.yaml` gains
`said`, `delegated` and `previous`; `study.py init --from-exploration` and `decide` have new required
options), and the order of the dashboard's lists. A study made with 0.3.0 still opens; its numbers do not
change, but lists ordered by a difference come in another order.

### Added

- **Link between the agent and the open dashboard** (`expression.py view`, `references/dashboard.md`).
  The page writes what the researcher has on screen to `.dashboard/researcher_view.json` and follows
  `.dashboard/agent_view.json`: `view` without choices prints the screen; with choices, the page moves to
  that view without reloading, on top of what is on screen (`--reset`: on top of the opening view), with
  a line showing the agent's `--message`. `summarize` and `figures` take `--from-dashboard`. The link can
  also send context tissues (`--shown`), highlights and the stretch of the window (`--window`,
  `--full-window`), and switch selections off (`--no-control`, `--all-domains`, …). Live dashboard only.
- **Each answer carries the researcher's words.** `study.py decide` and `init --from-exploration` require
  `--said` (the researcher's words the answer comes from) and take `--delegated` when the researcher
  explicitly left the choice to the agent; `decisions.yaml` stores `said` and `delegated`, and the report's
  *Escolhas do pesquisador* table shows them (*Nas palavras do pesquisador*), listing delegated answers apart.
- **The researcher talks freely; the workflow keeps its order** (SKILL.md). Answers given early are kept,
  matched to the real options and shown, and only what is missing is asked; requests ahead of their step
  are listed and done at that step; hypotheses and context answer nothing by themselves; delegation must
  name the question (a general "go ahead" is not delegation).
  - Step 0 asks for the tissues of both groups even when the request names some; tissues named only as
    context or in a hypothesis are not groups. Tissues that arrive after the exploration are explored
    before the study is created.
  - No default by silence; picking "the recommended one" after seeing the options is an answer, not a
    delegation; one recommendation per question.
  - `decide` refuses `--recommended` for `cohort.target`/`cohort.control` (a group has no single
    recommended term), and keeps the options' labels when a question is answered again with codes only.
  - Tested with four simulated researchers (all answers upfront, "do everything", hypotheses, change of
    mind): none of the 20 required answers was filled by the agent.
- **The researcher's step-1 answers are required.** `study.py init --from-exploration` refuses to create a
  study without `--window`, `--target`, `--control`, `--scorers` and `--question`, lists what is missing,
  and records `study.window`, `study.scorers` and `analysis.question` in `decisions.yaml` with the options
  offered. The agent asks the missing ones together and never fills one from the recommendation.
- **Score and group value on the dashboard** (*Avançado*: *Pontuação*, *Valor de cada grupo*;
  `summarize`/`figures` `--score raw|quantile`, `--group-agg mean|median`). Quantiles switch off the %
  filters, with a note in the filter steps.
- **`report.md` always has figures.** Without a view drawn by the researcher, `report` draws the opening
  view (the study's groups and question, no filter) under *Visão inicial do estudo*; Figure 1 is included
  as an image (`figures/fig1_genomic_architecture.svg`).

- **Decision history.** Deciding a key again with another value keeps the earlier answer in `previous`
  (value, label, date, note, options) in `decisions.yaml`; `report.md` lists them under *Redefinições*.
- **Scope section in SKILL.md.** Questions beyond the screen (transcription factors, motifs, mechanism,
  literature) are flagged as out of scope, worked outside the study folder with a record of what was done,
  and never end in proposed experiments.

- **`exploration/activity.md`**: every term of the Atlas catalogue with the gene's predicted activity (and
  relative to the most active), kind of sample, assays, sources and number of measurements, written by
  `explore`; `explore --out OUT --activity-only` writes it for an existing exploration without new captures.
  `summary.md` points to it.
- **Changing the groups after the results** is documented as a redefinition that follows the study's
  definition. The agent chooses from `activity.md` with the researcher, records the change with `decide`
  and recaptures; it no longer ranks tissues with its own code.
- **`summarize` prints the reference expression (REF) of each chosen tissue**: its % of the highest REF and
  what the current floor does to it, and the value that marks a tissue as not expressed (5% of the highest
  REF, as in the exploration). The floor can now be chosen from the numbers.

### Changed

- **'Ordenar por' follows the question.** The change in the target and the target − control differences
  are ordered in the direction asked (largest rise first with an increase, largest fall with a
  repression), not by absolute value: a variant whose control changes more than the target no longer
  comes first. `summarize` lists its first variants in this same order.
- **The blue box names the first variants in the order of 'Ordenar por'**, the same as the table's first
  rows (before: the ones whose least-changing target tissue changed most, whatever the order chosen). It
  also names the tissues that left a group for lack of measurements with the chosen assays and sources,
  and words the floor warning by the floor chosen.
- With every region ticked, *domínios* is no longer listed as a filter step.
- **Redefining a group needs no manual step.** `decide cohort.target|cohort.control` writes the new terms
  (and `--group-label`) into `study.yaml`, and a contrast plan the skill generated is generated again from
  the current groups.
- **The dashboard opens on the researcher's question** (`analysis.question`), and `summarize`/`figures`
  use it when `--question` is not given.
- `expression.py dashboard` prints `Dashboard READY` only after marimo answers, or `Dashboard FAILED`.

### Fixed

- **Statements that the sources contradict** (terminology review of 18/11/2026; wording only, no number
  changes):
  - the exploration's "atividade" is no longer described as the gene's activity "sem variante": it is the
    `RNA_SEQ_ACTIVE` score of one probe SNV, the larger of its two alleles; the pages say that this SNV
    was scored, that the 5% cut is this skill's convention and not absence of the transcript, and that
    tracks of different assays and sources are not directly comparable;
  - two tracks of one term are no longer called replicates (they differ in assay or source);
  - quantiles are ranks among common variants, not "against the genome", and a large quantile is an
    unusual score, not a strong effect;
  - `summarize` no longer presents 5% of the study's largest REF as the exploration's rule.
- **Documentation brought up to the code**: `README.md` (key in the project's `.env`; `check` installs
  nothing), `references/workflow.md` step 14 (`view`, `--from-dashboard`), `references/formats.md` (the
  two files of `.dashboard/`), SKILL.md (commands, common mistakes), `references/dependency_outputs.md`
  (the baseline file is copied into a study).
- **Where each score is defined**: `references/interpretation.md` now sends the agent to the `alphagenome`
  dependency skill's scorer table before it explains a score.
- The variants table showed `None` in the control column when no control tissue was chosen; the column is
  now left out.
- **Terms named exactly as an Atlas sample are found.** The exploration recommends the catalogue term whose
  sample name matches the researcher's text written any way ("MCF-7" = "MCF7", "T-cell" = "T cell"), even
  when the ontology search ranks it low, misses it or no longer knows it. Such a term enters the cohorts
  with the catalogue's name, as the rule `cohort.term_without_ols:{curie}`.
- A failed lookup of a term's descendants is noted in the exploration (`descendants_error`) instead of
  stopping the other tissues; an OLS record with `in_subset: null` is read as no subsets.
- `study.py init --from-exploration` leaves out raw files of the exploration without a manifest step (the
  error reply of a failed search) and says so, instead of refusing to create the study; `cohorts` skips
  unreadable OLS files of terms the study does not use.
- The *Links do Atlas* table in `report.md` had blank lines between rows and did not render as a table.

- **The noise filter follows the question** (*Alvo − controle acima do ruído*, `--min-sigma`): with
  repression the target must be below the control, with increase above, with either both ways. Before, a
  variant that raised the target could pass a repression question.
- A redefinition of the groups ends with `figures` and `report` run again on the redefined analysis
  (`references/workflow.md`).

- **Keys**: the AlphaGenome key is read only from the project's `.env`, not from `~/.env` or environment
  variables, so each study runs with the key of its own project.
  - `fetch` passes that file to every script that needs it, the interpretation step included.
  - `check` warns when the project's `.env` is not ignored by git: the AlphaGenome terms make the key
    personal.
- **Installation belongs to `alphagenome-skills-setup`**. `check` names what is missing and hands it to
  that skill, or says how to install it. It no longer prints install commands for each dependency, and
  SKILL.md no longer has its own install and terms procedure. This avoids two skills asking the
  researcher the same questions.

## [0.3.0] - 2026-11-16

First release after 0.1.0. It includes the work of 0.2.0, which was never released on its own.

### Added

- **Exploration before the study** (`expression.py explore`). Before any study is fixed, it gathers:
  - the gene structure and the Ensembl Regulatory Build;
  - candidate windows with their size and cost;
  - ontology terms for the tissues named by the researcher, with Atlas tracks per tissue;
  - the gene's predicted activity per tissue (`RNA_SEQ_ACTIVE` of one baseline variant).

  The result is a Portuguese options page, `exploration/summary.md`, with `exploration.json`. Then
  `study.py init --from-exploration` creates the study from the researcher's choices and copies the
  captures, the baseline included, with their provenance.
- **Captures by step** (`expression.py fetch <kind>`). Each step runs its dependency skill with pinned
  package versions, writes the file the stages read, and records the exact command in `manifest.json`.
  This replaces hand-typed commands.
- **Analysis dashboard** (`assets/dashboard.py`, computation in `scripts/core/analysis.py`). The
  researcher chooses on the page:
  - what they look for (decrease, increase, any change);
  - the target and control tissues;
  - how results are shown (order, how many, whether variants that do not pass appear in grey);
  - filters in %, all off at first: minimum change in every target tissue, maximum change in the
    control, target − control above the noise, regions;
  - technical settings: assays and sources, expression floor, prediction type.

  The page shows:
  - the change per tissue in %;
  - the noise between measurements of the same tissue;
  - the gene's predicted reference expression.

  Its charts are the window lollipop, target × control, the variants table, per-measurement bars and
  variants × tissues, each with PNG/SVG download. The page never calls the network and saves nothing.
- **`expression.py summarize`**: the dashboard's view as text for given choices, so the agent explains
  what the researcher sees from computed numbers.
- **`expression.py figures`**: the dashboard's charts for the researcher's choices, as PNG in
  `figures/dashboard/`, with `choices.json`.
- **Report section *Análise do pesquisador (dashboard)***. It shows:
  - the choices in the page's words, with CURIEs;
  - the result;
  - the filter steps;
  - figures D1–D4 with captions.

  The report warns when these figures are older than the variant table.
- **Questions in plain language.** Each question asked of the researcher has a title, why it matters,
  the identifiers it is about (CURIEs, transcript and region IDs, coordinates, counts), and options with
  a label and what happens. The agent asks in the researcher's language and records a code.
  - `study.py pending --json` prints the questions in a machine-readable form.
  - `decide` refuses answers that were not offered.
- **Rules recorded by the skill.** Checks the researcher cannot judge are applied by fixed rules and
  recorded in `decisions.yaml` with `decided_by: skill` and the rule in words. They are evaluated again on
  every run. They cover:
  - positions where the Atlas and the reference sequence disagree;
  - prediction types without data in one group;
  - the window's definition;
  - the triage's contrast plan.

  Above 1% of the window, excluded positions become a question.
- **Window question `exclude_feature`.** A neighbouring regulatory region can be left out of the window.
  The window questions come one at a time.
- **`expression.py check --json`**: missing prerequisites as questions, with the skill that resolves
  each. Routes to the planned `alphagenome-skills-setup` skill when it is installed.
- **`references/dashboard.md`**: what to say when presenting the page, every *Avançado* field, and what
  the page does without saying.
- **`references/formats.md`**: the fields of the study files (questions, decisions, exploration,
  contrast plan).
- Repository: `.gitattributes`, LF line endings, and study results and evals kept out of git.

### Changed

- **Researchers are asked fewer questions.** Asked: the transcript, the window's extension, and tissues
  without predictions. Recorded as rules instead:
  - `cohort.scorer_without_cohort:*`;
  - `variants.ref_mismatch:*` and `variants.missing_in_atlas:*`;
  - `locus.exon1_overlap:*`, and `locus.multiple_tss_features` when every region touches the start of
    the gene.
- **The contrast plan is no longer asked** (step 12). Without `contrasts.yaml`, `contrast` generates it
  from the groups. `contrast.review:*` became notes in the report.
- **Tissue replacements** are ranked by RNA-seq data and marked tissue or cell type, with what each
  predicts. Removing a tissue is recommended only while its group keeps RNA-seq.
- **`report.md`**:
  - *Decisões* separates the researcher's choices (labels and codes) from the rules applied by the skill;
  - *Resultados* is now *Resultados da triagem fixa*, and the triage labels (*Inconsistent*, *Mixed*,
    the quantile filter) are described as heuristics;
  - the exploration page no longer shows decision keys.
- **Workbench**: the previous triage workbench is still available as `dashboard --legacy`;
  `dashboard --edit` opens the marimo editor, for maintainers.
- **Dependencies are pinned in every script** (PEP 723). `polars` is fixed at 1.31.0, which runs under
  emulated x86-64 Python on Windows ARM64; newer versions stop with `unknown feature flag: 'sse3'`.
  Consoles use UTF-8.
- **`SKILL.md`, `workflow.md` and `interpretation.md`** were rewritten for the exploration, the
  dashboard, plain-language questions and rules, and the final answer. The final answer separates the
  researcher's decisions, the skill's rules and the agent's choices.

### Fixed

- **Percent thresholds for decreases.** A −30% threshold is now exactly −30% (log 0.7). Before, it used
  the same log as +30%.
- **Filters that could not apply** (no control tissue, no noise estimate) removed every variant. They are
  now listed as *ignorado*, with the reason.
- **The control filter** checks every control tissue in both directions, not only the mean.
- **The master QC** counts positions excluded by rule.
- **`decide`**:
  - accepts `yes`/`no`-like codes, placeholder codes such as `fixed:<value>`, and Ensembl IDs with or
    without version;
  - questions a stage no longer raises leave `pending`.
- **A researcher's approval of an earlier contrast plan** is no longer carried to a newly generated one.
- **Dashboard**:
  - *Quantas mostrar* limits what is drawn instead of acting as a filter;
  - variants that do not pass are grey in every chart;
  - highlights and crop no longer move the window unexpectedly.
- **Ctrl+C** on `dashboard` ends cleanly (exit 130) and records the session.

### Known limitations

- `report.md` and the dashboard are in Portuguese. The questions follow the researcher's language.
- `study.py`, `manifest.py`, `models.py` and `reference.py` are shared with the splicing skill
  (`cnpem/plp1-splicing-analysis`), which still has the 0.1.0 versions.
- `reproduce` covers the pipeline stages, not the dashboard, `summarize` or `figures`.

## [0.2.0] - 2026-11-15 (not released)

Development version on a working branch. Its changes are listed under 0.3.0.

## [0.1.0] - 2026-11-13

First version:
- a 14-step workflow from gene to report;
- dependency outputs validated against schemas;
- study, decisions and manifest files;
- locus and regulatory window, cohorts against the Atlas catalogue, saturated SNV table with QC gates,
  Tier 2 contrasts with direction and scope labels;
- a Portuguese report and a marimo triage workbench;
- an offline UCHL1 reference study checked by `reproduce`.

## Versioning

- **Version number**: `metadata.version` in `skills/expression-analysis/SKILL.md`. It is what studies
  record.
- **While the version is 0.x**:
  - the middle number goes up when a change can alter a study's numbers, the questions asked, or the
    format of study files (0.3.0 → 0.4.0);
  - the last number goes up for fixes and documentation that do not change results (0.3.0 → 0.3.1);
  - 1.0.0 marks a skill considered stable for use outside the team.
- **One version per merge to `main`.** A pull request that changes the skill updates `metadata.version`
  and adds its entry here. Small changes collect under *Unreleased* until the next version.
- **Every released version** gets a git tag `vX.Y.Z` on the merge commit and a GitHub release whose
  notes are its entry here.
- **Changes to the files shared with the splicing skill** (`study.py`, `manifest.py`, `models.py`,
  `reference.py`) go to both repositories.

[Unreleased]: https://github.com/cnpem/uchl1-expression-analysis/compare/v0.4.0...HEAD
[0.4.0]: https://github.com/cnpem/uchl1-expression-analysis/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/cnpem/uchl1-expression-analysis/compare/v0.1.0...v0.3.0
[0.1.0]: https://github.com/cnpem/uchl1-expression-analysis/releases/tag/v0.1.0
