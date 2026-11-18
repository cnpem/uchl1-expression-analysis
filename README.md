# expression-analysis

Agent skill for comparing AlphaGenome-predicted regulatory and expression effects across every SNV in a human gene's regulatory window (GRCh38) between target and control tissues or cell types (UBERON/CL), with an interactive marimo dashboard in which the researcher chooses the groups, the question and the filters.

> **Research use only.** All outputs are predicted molecular impacts used to prioritise experiments, never measurements, clinical classifications or diagnoses.

## Installation

Install the skill with the Skills CLI:

```bash
npx skills add cnpem/uchl1-expression-analysis --skill expression-analysis
```

Or copy `skills/expression-analysis` into `.claude/skills/`, `.agents/skills/`, `~/.claude/skills/`, or `~/.agents/skills/`.

## Prerequisites & Dependencies

Run the prerequisite check (it names what is missing and never installs anything; the `alphagenome-skills-setup` skill installs it, one approved item at a time):

```bash
uv run skills/expression-analysis/scripts/expression.py check
```

External dependency skills (installed separately by the researcher as needed):
- From `google-deepmind/science-skills`: `uv`, `credentials`, `alphagenome-variant-impact-score`, `alphagenome-single-variant-analysis`, `alphagenome-atlas-website-links`, `embl-ebi-ols`, `ensembl-database`
- From `K-Dense-AI/scientific-agent-skills`: `alphagenome`

Set `ALPHAGENOME_API_KEY` in the project's `.env` (the only place the key is read from; keep `.env` out of git) and create the required `.licenses/` acknowledgement files only after reviewing each dependency's terms of use.

## Quick Start

Reproduce the bundled UCHL1 reference study offline:

```bash
uv run skills/expression-analysis/scripts/expression.py reproduce --reference uchl1 --keep /tmp/uchl1-demo
uv run skills/expression-analysis/scripts/expression.py dashboard --study-dir /tmp/uchl1-demo
```

See [`skills/expression-analysis/SKILL.md`](skills/expression-analysis/SKILL.md) and [`skills/expression-analysis/references/workflow.md`](skills/expression-analysis/references/workflow.md) for the complete 14-step workflow.

## Versions

The skill version is `metadata.version` in [`SKILL.md`](skills/expression-analysis/SKILL.md); each study records it in its manifest and report. What changed in each version, and the versioning rules, are in [`CHANGELOG.md`](CHANGELOG.md). Released versions are tagged `vX.Y.Z`.

## Notes on Reference Data & License

- Code and skill documentation in this repository are licensed under the MIT License ([`LICENSE`](LICENSE)).
- `skills/expression-analysis/assets/reference/uchl1/raw/` contains bundled reference outputs (AlphaGenome Atlas, Ensembl, EMBL-EBI OLS4) used for offline reproduction. Verify AlphaGenome redistribution terms before sharing publicly.
