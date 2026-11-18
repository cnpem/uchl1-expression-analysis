# Changelog: alphagenome-skills-setup

Changes to the `alphagenome-skills-setup` skill. The `expression-analysis` skill keeps its own
changelog at the repository root ([`CHANGELOG.md`](../../CHANGELOG.md)); each skill has its own version.

Format based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow [Semantic Versioning](https://semver.org/), with the rules of the analysis skills' changelog (`cnpem/uchl1-expression-analysis`).

## [Unreleased]

## [0.2.0] - 2026-11-18

Skills go where the agent app reads them. In a test with Antigravity, the analysis skill sat in a
`.skills/` folder next to this one; 0.1.0 counted it as installed, the app never listed it, and the agent
ran its own analysis with the dependency skills, choosing the window and the tissues without the researcher.

### Added

- **The agent app is asked first and recorded** (`setup.agent:<id>`, `agents` in `.alphagenome-setup/state.json`).
  `plan` offers nothing else before the answer; later commands prepare the recorded apps (`--agents id[,id]`
  overrides; an older state infers them from earlier installs).
- **`setup.place:<name>`**: a copy of an analysis skill, or of this setup skill, found outside the app's folder
  (next to this skill, in `.skills/`, in the project root or in another app's folder) is copied into it with
  the researcher's approval; the old copy is left as it is and the placement is recorded.
- **Notices for services that cannot be reached** in `plan` and `verify`, with what each one affects (e.g.
  Ensembl REST: the Regulatory Build and the reference sequence). `assets/targets.json` names each service.
- `references/keys.md`: step-by-step guide to creating the AlphaGenome account and key (and the optional NCBI key).

### Changed

- **A skill counts as installed only in a folder the chosen app reads**: its project folder or its home
  folder, from the agent table of npx skills 1.9.0 (Claude: `.claude/skills`; Antigravity, Codex, Cursor,
  Gemini CLI, OpenCode: `.agents/skills`; `assets/targets.json`). The folder this skill sits in no longer
  counts. `setup.skill:<name>` installs only for the apps that lack it, and checks the folder after npx.
- **`verify` checks placement** (each analysis skill and this one in every app's folder) before `check`
  and `reproduce`, and the setup is finished only after `verify` and `report`. The report shows the apps,
  their folders, copies outside them, and says when verify has not run.
- `--agents` takes comma-separated ids (`--agents antigravity,claude-code`) or is repeated; with several ids
  after it, argparse took the subcommand for an id.

- **Lives next to `expression-analysis`** in `cnpem/uchl1-expression-analysis` (`skills/alphagenome-skills-setup`),
  moved with its history from `romulorcruz/alphagenome-skills-setup`. Install it with
  `npx skills add cnpem/uchl1-expression-analysis --skill alphagenome-skills-setup --copy`.
- **Keys** live only in the project's `.env`; `~/.env` and environment variables are not used. A key found there is reported as unused.
  - `install setup.key:<VAR>` adds `.env` to `.gitignore` first, then opens the page where the key is created (AlphaGenome: https://alphagenome.google/api) and the project's `.env` in the text editor.
  - `diagnose`, `key` and the report say where the key is, and warn when git does not ignore `.env`.
- **Terms of use**: `terms <skill> --open` explains in plain words what each set of terms covers and opens the pages (`assets/terms.json`; AlphaGenome, EMBL-EBI OLS, Ensembl).


## [0.1.0] - 2026-11-16

### Added

- `setup.py`, standard library only:
  - `diagnose`, `plan`, `terms`, `install`, `key`, `verify` and `report`;
  - each install needs `--approved` naming the same item.
- Installing the skills with `npx skills add --copy -y --json` for the agents chosen. Each skill's source and hash are recorded in `.alphagenome-setup/state.json` and the report.
- Terms of use read from each dependency skill's SKILL.md and recorded in `.licenses/` in the format those skills expect.
- The AlphaGenome key: `~/.env` with an empty line opened in the system's text editor. Only the key's presence is checked.
- `assets/targets.json`: the analysis skills (expression-analysis, splicing-analysis), their dependencies, keys, checks and reference studies.
- Platform notes and troubleshooting.

### Tested

- Windows 11 on ARM64, with a workspace holding only expression-analysis 0.3.0:
  - the 8 dependency skills installed one by one;
  - 3 terms of use recorded;
  - key detected;
  - `check` and `reproduce --reference uchl1` passed.
- Not yet tested on macOS or Linux, or on splicing-analysis.
