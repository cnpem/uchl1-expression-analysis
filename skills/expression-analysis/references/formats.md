# Study file formats (expression-analysis)

The fields of the files this skill writes and reads around the researcher's choices. The raw files
written by the dependency skills are in `dependency_outputs.md`. Read this when you need a field the
commands do not print; edit none of these files by hand except `study.yaml` (with the researcher) and,
for a maintainer, `contrasts.yaml`.

## Contents

- [Questions: `processed/pending_decisions.json`](#questions-processedpending_decisionsjson)
- [Prerequisites: `expression.py check --json`](#prerequisites-expressionpy-check---json)
- [Answers and rules: `decisions.yaml`](#answers-and-rules-decisionsyaml)
- [Exploration: `exploration/exploration.json`](#exploration-explorationexplorationjson)
- [Contrast plan: `contrasts.yaml`](#contrast-plan-contrastsyaml)

## Questions: `processed/pending_decisions.json`

A list; `study.py pending --json` prints it. One item per open question:

| Field | Meaning |
|---|---|
| `key` | `<area>.<topic>[:<subject>]`, e.g. `locus.exon1_overlap:ENSR…`; what `decide --key` takes |
| `title` | a short name for the question |
| `question` | the question in plain English, with names and identifiers |
| `why` | why the answer matters for the conclusions |
| `identifiers` | what the question is about: CURIEs, transcript and region IDs, coordinates, counts |
| `options` | list of `{code, label, effect}`: the code is recorded, the label and effect are shown |
| `recommended` | the code of the suggested option, or null |
| `agent_note` | an instruction for the agent (what to run before asking, how to present a long list); never shown to the researcher |
| `workflow_step` | the step that asked; what `decide --step` takes |
| `asked_at` | when it was first raised (ISO-8601, UTC) |

The texts are written once, in English; present them in the researcher's language. Older studies may
list bare codes in `options`; they are read as `{code, label: code}`.

## Prerequisites: `expression.py check --json`

`{ready, setup_skill, setup_skill_installed, items}`. Each item has the fields of a question above
(`key` under `setup.`: `setup.uv`, `setup.skill:{name}`, `setup.licence:{skill}`, `setup.key:{VAR}`)
plus `handled_by`, the skill that resolves it (`uv`, `credentials`, or the setup skill). Nothing is
recorded in a study: the setup happens before one exists. Exit code 1 while items remain.

## Answers and rules: `decisions.yaml`

A mapping from key to entry:

| Field | Always | Meaning |
|---|---|---|
| `value` | yes | the recorded answer: an option's code, or a list (e.g. CURIEs for `cohort.target`) |
| `workflow_step` | yes | the step that asked |
| `asked_at`, `answered_at` | yes | ISO-8601 times |
| `note` | yes | the researcher's reason, in their words (may be empty for a rule) |
| `options`, `recommended` | together, when the question offered options | what was offered (`{code, label}`; bare codes in older studies) and the recommended code |
| `value_label` | when the value is one of the options | the label the researcher saw |
| `decided_by` | only for rules | `skill`: applied by a fixed rule, without asking (absent = the researcher) |
| `rule` | only for rules | the rule and what it affected, in words |
| `previous` | only after a redefinition | earlier answers to the same key, oldest first: `value`, `value_label`, `answered_at`, `note`, `options`, `recommended` as they were; `decide` appends the old entry when the value changes |

Rules in use: `locus.multiple_tss_features` (`union`), `locus.exon1_overlap:{feature}`
(`keep_full_unit`), `cohort.scorer_without_cohort:{scorer}` (`drop_scorer`), `variants.ref_mismatch:{pos}` and
`variants.missing_in_atlas:{pos}` (`exclude_position`; above 1% of the window only after the researcher's
`variants.exclusion_limit` = `exclude_all`), `contrast.plan` (`generated`, with the plan's sha256 in
`note`; `existing` when the `contrasts.yaml` in the study was written or edited outside the skill). They are evaluated again on every run and removed when they no
longer apply; a researcher's answer under the same key is never overwritten.

## Exploration: `exploration/exploration.json`

What `explore` found, read by `study.py init --from-exploration`:

| Field | Meaning |
|---|---|
| `gene`, `transcript`, `transcript_note`, `chromosome`, `strand`, `tss`, `gene_span`, `n_transcripts` | the gene and the transcript the exploration used (`transcript_note` says when the study will ask `locus.transcript`) |
| `search_region`, `flank_bp` | where the Regulatory Build was searched |
| `scorers` | the prediction types looked up |
| `features` | regulatory regions of the search region: `id`, `biotype`, `start`, `end`, `distance_to_tss` |
| `windows` | the candidate windows: `id` (what `init --window` takes), `label`, `region`, `start`, `end`, `bp`, `snvs`, `policy` (the `study.yaml` window policy), `why`, `study_questions` (the step 4 questions it will raise), `features_included`, `features_excluded`, `needs_max_window`, `rows_per_scorer` (Atlas rows, the cost) |
| `recommended_window` | the id of the suggested window |
| `terms` | one per tissue text: `text`, `ols_file`, `candidates` (`curie`, `label`, `ontology`, `definition`, `per_scorer` track counts, `biosample_names`), `recommended` CURIE, `related_in_catalogue` (catalogue terms whose name contains the text), `descendants` (more specific terms with tracks, for a term without RNA-seq) |
| `baseline_variant` | the single SNV whose `RNA_SEQ_ACTIVE` score gives the gene's activity per tissue |
| `activity` | `max`, `floor`, `median_relative`, `top` (most active terms) and `by_curie`: `name`, `activity`, `relative` (to the maximum), `tracks`, `expressed` (relative ≥ floor) |
| `recommended_tracks_per_scorer` | tracks of the recommended terms, per prediction type |

## Contrast plan: `contrasts.yaml`

Generated by `contrast` from the groups when absent (the fixed triage behind the report); a maintainer
may write one before running `contrast`, starting from `contrast --propose`.

```yaml
study_id: <study_id>
contrasts:
- name: rna_seq               # used in column names and in priority
  scorer: RNA_SEQ
  strand: unstranded          # sense | antisense | unstranded | any (relative to the gene)
  track_filter: {}            # optional allow-lists: assay_title, data_source, biosample_type
  target: [<CURIE>, ...]      # the tissues of each group that this scorer measures
  control: [<CURIE>, ...]
  tissue_aggregation: mean    # mean | median of a tissue's tracks
priority:
  direction_contrast: rna_seq # the contrast that orders the top list
  scope_contrast: rna_seq     # optional; default the direction contrast
```
