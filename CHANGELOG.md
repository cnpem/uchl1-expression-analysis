# Changelog

Changes to the `expression-analysis` skill that matter to the people who use it and to the studies it
produces. Format based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/) as described in [Versioning](#versioning) below.

Each study records the skill version it ran with (`skill_version` in `manifest.json` and in each
step's `source_versions`; also shown in `report.md`). To know what a study's numbers rest on, find that version here.

## [Unreleased]

### Changed

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

[Unreleased]: https://github.com/cnpem/uchl1-expression-analysis/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/cnpem/uchl1-expression-analysis/compare/v0.1.0...v0.3.0
[0.1.0]: https://github.com/cnpem/uchl1-expression-analysis/releases/tag/v0.1.0
