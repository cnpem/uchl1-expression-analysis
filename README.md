# expression-analysis

Agent skill for comparing AlphaGenome-predicted regulatory and expression effects across every SNV in a human gene's regulatory window (GRCh38) between target and control tissues or cell types (UBERON/CL), with an interactive marimo dashboard in which the researcher chooses the groups, the question and the filters.

> **Research use only.** All outputs are predicted molecular impacts used to prioritise experiments, never measurements, clinical classifications or diagnoses.

This repository holds two skills:

| Skill | What it is for |
|---|---|
| [`expression-analysis`](skills/expression-analysis/SKILL.md) | the analysis itself |
| [`alphagenome-skills-setup`](skills/alphagenome-skills-setup/SKILL.md) | prepares the computer for it (and for `splicing-analysis`), for researchers who never used a terminal |

## Installation

### First time on this computer: the setup skill

Install the setup skill and ask the agent, for example, "prepare my computer to use expression-analysis":

```bash
npx skills add cnpem/uchl1-expression-analysis --skill alphagenome-skills-setup --copy
```

It checks what is there and installs each missing piece (uv, Node.js, `expression-analysis`, its dependency
skills, the terms of use, the AlphaGenome key in the project's `.env`) only after the researcher approves
that item, then runs `check` and the offline `reproduce`. Details in
[`skills/alphagenome-skills-setup/SKILL.md`](skills/alphagenome-skills-setup/SKILL.md).

### The analysis skill only

```bash
npx skills add cnpem/uchl1-expression-analysis --skill expression-analysis
```

Or copy `skills/expression-analysis` into the folder your agent app reads: `.claude/skills/` for Claude Code,
`.agents/skills/` for Antigravity, Codex, Cursor, Gemini CLI and OpenCode (or the app's folder in your home).
A copy kept anywhere else (for example a `.skills/` folder, or `.agents/skills/` when the agent is Claude Code)
is not listed among the agent's skills, and the agent then answers with the dependency skills alone, choosing
windows and tissues without the researcher. `expression.py check` warns when the skill is in no agent's folder,
and the setup skill checks the folder of each agent chosen.

## Prerequisites & Dependencies

Run the prerequisite check (it names what is missing and never installs anything; the `alphagenome-skills-setup` skill installs it, one approved item at a time):

```bash
uv run skills/expression-analysis/scripts/expression.py check
```

External dependency skills (installed separately by the researcher as needed):
- From `google-deepmind/science-skills`: `uv`, `credentials`, `alphagenome-variant-impact-score`, `alphagenome-single-variant-analysis`, `alphagenome-atlas-website-links`, `embl-ebi-ols`, `ensembl-database`
- From `K-Dense-AI/scientific-agent-skills`: `alphagenome`

Set `ALPHAGENOME_API_KEY` in the project's `.env` (the only place the key is read from; keep `.env` out of git) and create the required `.licenses/` acknowledgement files only after reviewing each dependency's terms of use. The setup skill guides both.

## Quick Start

Reproduce the bundled UCHL1 reference study offline:

```bash
uv run skills/expression-analysis/scripts/expression.py reproduce --reference uchl1 --keep /tmp/uchl1-demo
uv run skills/expression-analysis/scripts/expression.py dashboard --study-dir /tmp/uchl1-demo
```

See [`skills/expression-analysis/SKILL.md`](skills/expression-analysis/SKILL.md) and [`skills/expression-analysis/references/workflow.md`](skills/expression-analysis/references/workflow.md) for the complete 14-step workflow.

## What the numbers mean

Every quantity the pages show is defined once in [`assets/terms.yaml`](skills/expression-analysis/assets/terms.yaml), with its name in Portuguese and English, formula, unit and source; [`references/glossary.md`](skills/expression-analysis/references/glossary.md) is generated from it (`expression.py terms --write-glossary`), and [`references/interpretation.md`](skills/expression-analysis/references/interpretation.md) explains how each score is computed and what the model does not establish. `expression.py terms` checks that the code and the documentation use only the registry's names.

## Versions

The skill version is `metadata.version` in [`SKILL.md`](skills/expression-analysis/SKILL.md); each study records it in its manifest and report. What changed in each version, and the versioning rules, are in [`CHANGELOG.md`](CHANGELOG.md). Released versions are tagged `vX.Y.Z`.

`alphagenome-skills-setup` has its own version (`metadata.version` in its [`SKILL.md`](skills/alphagenome-skills-setup/SKILL.md)) and
[changelog](skills/alphagenome-skills-setup/CHANGELOG.md); its tags are `setup-vX.Y.Z`.

## Notes on Reference Data & License

- Code and skill documentation in this repository are licensed under the MIT License ([`LICENSE`](LICENSE)).
- `skills/expression-analysis/assets/reference/uchl1/raw/` contains bundled reference outputs (AlphaGenome Atlas, Ensembl, EMBL-EBI OLS4) used for offline reproduction. Verify AlphaGenome redistribution terms before sharing publicly.
