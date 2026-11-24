# Commands, prerequisites and services: expression-analysis

Every command of the skill, what `check` tests, the dependency skills, the services reached and their
limits, the reference reproduction and the citations. The step-by-step commands and questions are in
`workflow.md`.

## Contents

- [What the skill does](#what-the-skill-does)
- [Quick start](#quick-start)
- [Utility scripts](#utility-scripts)
- [Prerequisites: what `check` tests](#prerequisites-what-check-tests)
- [Dependencies](#dependencies)
- [Services and rate limiting](#services-and-rate-limiting)
- [Citation](#citation)

## What the skill does

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
created.

**Dependency skills fetch, this skill analyses.** Gene structure, ontology terms, track catalogue and
scores come from the DeepMind science skills and the K-Dense `alphagenome` skill. `expression.py fetch
<kind>` runs each of them with pinned package versions, writes the output unedited under the name the next
stage reads, and records the exact command in `manifest.json`. The stages read only the study folder, so a
study can be re-run offline and every number traces back to a recorded command. You choose the parameters;
the script runs and records them, so no shell syntax has to be adapted between PowerShell, bash and zsh.

## Quick start

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

## Utility scripts

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
| `expression.py view --study-dir D [choices…] [--reset] [--message "…"]` | Link with the open dashboard: without choices, prints what the researcher has on screen; with them, the page moves to that view without reloading (`dashboard.md`) |
| `expression.py summarize --study-dir D [choices… \| --from-dashboard]` | The dashboard's view as text for the same choices (groups, question, filters): headline, tracks, between-track disagreement, filter steps, top variants. Writes nothing |
| `expression.py figures --study-dir D [choices… \| --from-dashboard] [--top-k K] [--show kept\|all] [--highlight V…]` | The dashboard's charts for those choices as PNG in `figures/dashboard/`, with `choices.json`; `report` then shows them |
| `expression.py terms [--write-glossary]` | Term registry (`assets/terms.yaml`): checks that every entry has a source, that `glossary.md` is current and that no deprecated name is left in the code or the documentation; `--write-glossary` regenerates the glossary (maintainers) |
| `expression.py regions [--only NAME]` | Prints the regulatory search region, the window and `max_window` (the fetch commands use them already) |
| `expression.py reproduce --reference uchl1 [--until S] [--keep DIR]` | Re-runs the stages on the bundled reference and compares with `expected.json` |
| `study.py validate / pending [--json] / decide / require / window` | Study file, researcher decisions, open questions, `--print-max-window` |
| `manifest.py show` | Prints the manifest (`manifest.py add` remains for a step run outside the skill) |
| `models.py validate --dir D/raw` | Checks every raw file against the fields the dependencies write |

## Prerequisites: what `check` tests

The **`alphagenome-skills-setup`** skill installs and configures everything this skill needs, with the
researcher's approval for each item; this skill installs nothing and asks nothing about installation.
Run `expression.py check`. When it reports something missing, hand the setup to that skill. If it is not
installed, `check` prints how to install it; ask the researcher first (the repository is private for now).
Do not install dependency skills, accept terms of use or handle keys from this skill: doing it twice,
from two skills, asks the researcher the same questions twice.

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

Cited by the `name` of each skill. Versions tested are in `dependency_outputs.md`; the pinned package
lists are in `scripts/deps.py`.

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

## Services and rate limiting

**Services reached** (all through the fetch commands; nothing else goes to the network): the AlphaGenome Atlas
and model API (Google), Ensembl REST (`rest.ensembl.org`), EBI OLS (`www.ebi.ac.uk/ols4`), and PyPI and GitHub
when `uv` installs pinned packages. `reproduce` reaches none of them. When Ensembl does not answer, the
exploration goes on without the Regulatory Build and the reference sequence and says so in `summary.md`; the study
waits until a later `explore` fills them in (Workflow step 0).

`regulatory_features.py` paces its own Ensembl requests (15/s, file lock in
`~/.cache/alphagenome-gene-skills/`). Every other request comes from a dependency script with its
own limits (Ensembl and OLS through `polite-http`; the Atlas client's `--workers`). The fetch
commands give `polite-http` its own lock folder and reset a reservation left far in the future by
another clock, which otherwise makes OLS and Ensembl calls sleep silently for hours. Run the fetches
one at a time rather than spreading them over parallel subagents. Model-API calls of step 13 take
about 20 s each; run them sequentially.

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
