# Changelog: alphagenome-skills-setup

Changes to the `alphagenome-skills-setup` skill. The `expression-analysis` skill keeps its own
changelog at the repository root ([`CHANGELOG.md`](../../CHANGELOG.md)); each skill has its own version.

Format based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow [Semantic Versioning](https://semver.org/), with the rules of the analysis skills' changelog (`cnpem/uchl1-expression-analysis`).

## [Unreleased]

### Changed

- **Lives next to `expression-analysis`** in `cnpem/uchl1-expression-analysis` (`skills/alphagenome-skills-setup`),
  moved with its history from `romulorcruz/alphagenome-skills-setup`. Install it with
  `npx skills add cnpem/uchl1-expression-analysis --skill alphagenome-skills-setup --copy`.
- **Keys** live only in the project's `.env`; `~/.env` and environment variables are not used. A key found there is reported as unused.
  - `install setup.key:<VAR>` adds `.env` to `.gitignore` first, then opens the page where the key is created (AlphaGenome: https://alphagenome.google/api) and the project's `.env` in the text editor.
  - `diagnose`, `key` and the report say where the key is, and warn when git does not ignore `.env`.
- **Terms of use**: `terms <skill> --open` explains in plain words what each set of terms covers and opens the pages (`assets/terms.json`; AlphaGenome, EMBL-EBI OLS, Ensembl).

### Added

- `references/keys.md`: step-by-step guide to creating the AlphaGenome account and key (and the optional NCBI key).

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
