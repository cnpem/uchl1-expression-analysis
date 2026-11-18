# Interpreting the results (expression-analysis)

Read this before explaining any number to the researcher or reporting results. Every number comes from
AlphaGenome Atlas predictions: they prioritise experiments, they are not measurements, and they are for
research use only (never a diagnosis or a clinical finding).

## Contents

- [How to explain a number](#how-to-explain-a-number)
- [Where each score is defined](#where-each-score-is-defined)
- [The exploration's expression level (step 0)](#the-explorations-expression-level-step-0)
- [The dashboard's numbers](#the-dashboards-numbers)
- [Very low expression levels](#very-low-expression-levels)
- [Limits of the model to say out loud](#limits-of-the-model-to-say-out-loud)
- [Questions researchers ask](#questions-researchers-ask)
- [The fixed triage (report and reproduce)](#the-fixed-triage-report-and-reproduce)
- [Scope and limits](#scope-and-limits)
- [Legacy workbench](#legacy-workbench-expressionpy-dashboard---legacy)
- [Reporting](#reporting)

## How to explain a number

1. **Name it as the glossary does.** `references/glossary.md` (generated from `assets/terms.yaml`) has one
   entry per quantity: its name on the pages, what it is, formula, unit, what it may be compared with,
   status and source. Use the Portuguese name written there; do not coin another one.
2. **Say four things about it:** what it is; relative to what (a reference allele, a maximum, a group, the
   common variants); its unit, or that it has no official unit; and its status (stated by a source, this
   skill's convention, or not established).
3. **Take the value from a command, never from memory:** `summarize` (or `summarize --from-dashboard` for
   what is on screen), `report.md`, `view`, or the exploration pages. Do not describe values you have not
   computed.
4. **For an entry marked undecided, say that it is not established** and what would settle it. Do not fill
   the gap with an inference presented as fact; if you add your own reading, label it as yours.
5. **Give the level next to a large percentage** in a tissue where the reference level is low (see
   *Very low expression levels*).
6. **Keep the researcher's question in view:** a difference in the target that is not larger than the
   control's, or than the between-track disagreement, is not evidence that the effect is the target's own.
7. **Write names as the pages do** (`references/glossary.md`, *Writing conventions*): the gene symbol in
   italics and the protein in roman type; a variant by its Atlas identifier with its genomic HGVS beside it
   (`chr17:7676154:G>C`, `NC_000017.11:g.7676154G>C`), never calling the Atlas identifier HGVS; coordinates
   without thousands separators; Portuguese numbers with a decimal comma.

## Where each score is defined

The formula, the sign and the window of every AlphaGenome score are defined by the `alphagenome`
dependency skill, not here. Before explaining what a score is, read its `references/model-api.md`
(table *Recommended scorers*, column *Meaning of `raw_score`*) and its `references/interpretation.md`
(*Raw score and quantile: always both*, *Tissue matching*, *What the model cannot see*). The official
documentation of AlphaGenome has the same definitions (*Variant scoring*, tab *Active-allele scorers*;
*Exploring model metadata* for the track units). Do not define a score from memory.

| Score | What it is | Formula | Unit |
|---|---|---|---|
| `RNA_SEQ` | change of the gene's predicted RNA-seq signal, alternative against reference allele, in one track | `ln((mean_exons(ALT) + 0.001) / (mean_exons(REF) + 0.001))` | none (natural log of a ratio) |
| `RNA_SEQ_ACTIVE` | the gene's predicted level: the larger of the two alleles' mean predicted coverage over the gene's exons | `max(mean_exons(ALT), mean_exons(REF))` | no official unit; the track unit is "normalized read signal"; not TPM |
| `CAGE`, `DNASE`, `ATAC` | change of the predicted signal summed over 501 bp centred on the variant; belongs to no gene | `log2((sum_501bp(ALT) + 1) / (sum_501bp(REF) + 1))` | none (log2 of a ratio) |
| quantile score | the rank of a raw score among common variants, for the same prediction type and the same track (common variants, frequency above 1% in gnomAD v3: about 300,000 per the AlphaGenome FAQ, 348,126 SNPs of chromosome 22 per the preprint; whether the Atlas uses exactly that background is not established) | rank mapped to −1..1 (signed scores) or 0..1 | none |

The mask of `RNA_SEQ` and `RNA_SEQ_ACTIVE` is the union of the gene's annotated exons (GENCODE v46) that
fall in the scored interval. The Atlas holds these scores precomputed for every possible SNV. Per the
`alphagenome` skill's `references/atlas.md`, they come from a distilled AlphaGenome model with the AlphaGenome 1
scorer definitions; this was not checked against the Atlas report (status: not established).

## The exploration's expression level (step 0)

What the exploration computes, in order:

1. **One probe SNV.** The middle of exon 1 of the chosen transcript, with the first base other than the
   reference. It is arbitrary: the Atlas only answers per variant, so one variant is needed to read the
   gene's level. It is the only variant scored before the study.
2. **One Atlas request:** `RNA_SEQ_ACTIVE` of that SNV for the gene, on every RNA-seq track of the
   catalogue (`exploration/atlas/scores_RNA_SEQ_ACTIVE_baseline.tsv`, column `raw_score`).
3. **A mean per term:** the plain mean of the tracks of each ontology term. `activity.md` names each term by
   the Atlas catalogue's sample name, `summary.md` by the ontology's label, so one CURIE can show two names
   (for example, a stem-cell term whose catalogue sample is *neuronal* and whose ontology label is *neural*). A term's tracks can be of
   different assays (polyA plus, total RNA-seq) and sources (ENCODE, GTEx); the columns *Ensaios* and
   *Fontes* of `activity.md` say which went into each mean.
4. **A share of the gene's maximum:** each term's mean divided by the largest term mean in the whole
   catalogue (not only the tissues the researcher named).
5. **The 5% cut:** a term whose share is below 0.05 is flagged *abaixo de 5% do máximo do gene*.

How to read it:

- **It is a predicted level, not a change and not a measurement.** It is the probe SNV's larger allele: equal
  to the reference level when that SNV does not raise expression, and larger otherwise. In one track, when the
  probe's `RNA_SEQ` score L is positive, `level − REF = (level + 0.001) × (1 − e^(−L))`: L = 0.08 puts REF about
  7.7% below the level. The exploration does not fetch L, and its value per term is a mean over tracks, so the
  difference is unknown at step 0; once a study scores the window, `summarize` shows the reference level (REF),
  reconstructed per track (it is clipped at 0, and is imprecise when the level is near 0.001).
- **It has no official unit and is not TPM** (TPM is normalised across all genes of a sample; this score is
  computed one gene at a time). Its exact unit is not established. The track unit is "normalized read signal" (per-base read
  coverage normalised to reads per million, per the preprint's methods). Do not compare it with TPM from
  GTEx or the Human Protein Atlas.
- **Compare terms only within one kind of track.** No source says that levels of different assays or
  sources are comparable (status: not established). ENCODE polyA plus, ENCODE total and GTEx are separate
  kinds of track, and for one gene they can differ several-fold for the same biosample. A share of the
  maximum that mixes kinds, or a maximum that is a track of another kind, says little about the biology.
- **The 5% cut is this skill's convention.** No source gives a level above which a gene counts as
  expressed. A flagged term does not mean the transcript is absent; it means the gene's predicted level
  there is under 5% of its highest term.
- **It is the model's prediction for a training biosample**, with the limits below; the model's weak side
  is exactly which tissues express a gene more than others.

## The dashboard's numbers

`references/dashboard.md` describes every part of the page; these are the quantities behind it.

- **Change (%)** (*mudança prevista (%)*): `100 × (e^score − 1)` for `RNA_SEQ`, `100 × (2^score − 1)` for
  `CAGE`, `DNASE`, `ATAC`. For a tissue with several tracks the score is first averaged over the tracks, and
  the % is computed on that mean of logs (a geometric mean of the ratios).
- **Group value** (*valor do grupo*): the mean, or the median, of the group's tissue values, each tissue
  weighing the same whatever its number of tracks.
- **Target − control difference** (*diferença alvo − controle*): target group value minus control group
  value, in the score's own scale (log).
- **Between-track disagreement** (*discordância entre trilhas*): for each variant, the absolute difference
  between the scores of two tracks of the same tissue. The two tracks of a pair always differ in assay or
  source; they are not replicates, and the value says nothing about the model's error. It is pooled over
  every tissue with two or more tracks among the chosen assays and sources, in five bins of the pair's mean
  absolute score.
- **Typical disagreement** (*discordância típica*): `1.4826 × median` of that difference in each bin (a
  robust scale: the constant makes it equal to the standard deviation when the differences are normally
  distributed around zero; a few very large differences do not move it). The filter *alvo − controle acima da discordância* and the order *Diferença em unidades de
  discordância* divide the target − control difference by the typical disagreement of the bin of
  `(|target| + |control|) / 2`. With no tissue holding two tracks, there is none and the filter is ignored.
  Say how many tissues provide pairs: in a study where only a few control tissues do, it measures the
  disagreement between assay kinds in those tissues.
- **REF** (*nível de expressão previsto com o alelo de referência*): the reference allele's level in each
  track, reconstructed from the probe SNV's two scores in the study:
  `REF = ACTIVE` when its `RNA_SEQ` score is ≤ 0, otherwise `(ACTIVE + 0.001) × e^(−score) − 0.001`.
  Same scale as the exploration's level; a tissue's REF is the mean over its tracks. This is algebra on the
  two scorer formulas, not a value the Atlas returns.
- **Three different thresholds.** Do not mix them up:

  | Threshold | Of what | Where |
  |---|---|---|
  | *abaixo de 5% do máximo do gene* | share of the gene's catalogue maximum of the exploration's level | exploration pages |
  | *piso do nível REF* | an absolute REF value chosen by the researcher (off by default) | dashboard *Avançado*, `--ref-floor` |
  | 0.01 | an absolute REF value (10 × the `RNA_SEQ` pseudocount), used by the blue box when no floor is set | dashboard headline |

  `summarize` also prints each tissue's REF as a share of the largest REF **among the study's tissues**;
  that share is not the exploration's rule, whose maximum is over the whole catalogue.
- **The blue box warns about low REF only for target tissues.** A control tissue with a very low level is
  not flagged there: check the REF table of `summarize` before explaining a control value.
- **Quantile mode** (*Pontuação: quantil do Atlas*): the % filters are switched off. A large quantile means
  the score is unusual among common variants, not that the effect is large: in a promoter window most scores
  can sit above 0.99 while being changes of a few percent. The per-track bars stay in the raw scale.

## Very low expression levels

Where REF is close to the `RNA_SEQ` pseudocount (0.001), what can be said is arithmetic:

- A large % can be a minimal absolute change. The alternative allele's level is
  `ALT ≈ (REF + 0.001) × (1 + %/100) − 0.001` (exact for a tissue with one track). Give REF and ALT when you
  quote a large % in such a tissue: "+120% from 0.012 to 0.026" is not "the variant turns the gene on".
- A group value can be dominated by such tissues. When a control value is much larger than the target's,
  look at its per-tissue columns and REF before concluding that the effect is shared or stronger outside
  the target; re-run `summarize` with the floor applied to the control (`--ref-floor N --floor-applies
  control`) and report both.
- **What the sources say, and what they do not.** The paper reports lower accuracy across tissues for lowly
  and for narrowly expressed genes (how well the model predicts which tissues express a gene more). It does not
  say whether a variant's predicted change is less reliable in a tissue where the gene's level is low: that is
  not established. The pseudocount alone shrinks a log fold change near 0.001, it does not inflate it. Say that
  a large % there is fragile because it is a small absolute change, and give the model's limit on tissue
  specificity, without claiming a known behaviour of variant effects at low levels.
- Call it *aumento previsto a partir de um nível muito baixo*, with the levels; it is not evidence that a
  tissue expresses the gene.

## Limits of the model to say out loud

From the AlphaGenome paper (Avsec et al., Nature 2026) and its documentation, as gathered in the skill's
terminology review:

- **Tissue specificity is the hard case.** The authors write that capturing cell-type-specific expression
  deviations remains challenging; accuracy across tissues is lower for lowly and for narrowly expressed
  genes. The skill's question, target against control, is that case: say so with any claim that an effect is
  the target's own.
- **No uncertainty estimate.** The model returns no confidence for a prediction; neither the
  between-track disagreement nor a quantile replaces one.
- **Typical failures** named by the authors: magnitudes predicted too low, and errors in which tissues
  express a gene.
- **Window and variant type:** only the chosen window is mutated, SNVs only; effects of distal elements
  (beyond about 100 kb) are underestimated.
- **The Atlas uses a distilled model** with a 1 Mb input and the AlphaGenome 1 scorers, per the `alphagenome`
  skill; not checked against the Atlas report.

## Questions researchers ask

- *"What is this level, and in what unit?"* It is the gene's predicted expression level (`RNA_SEQ_ACTIVE` of
  the probe SNV, larger allele), averaged over a term's tracks. No official unit; the track unit is
  normalised read signal; not TPM. "N% do máximo" is a share of the gene's highest term in the catalogue.
- *"Is it the same as expression, activation, REF?"* The level and REF are the same kind of quantity (the
  level is ≥ REF). A change (%) is a different quantity. "Activation" is not an AlphaGenome quantity.
- *"Can I compare a GTEx tissue with an ENCODE cell line, or polyA with total?"* Not their levels: no source
  says they are comparable, and they differ several-fold for the same biosample. Not their effects either, as a
  fact: the AlphaGenome documentation says raw scores of different tracks are on different scales and points
  to the quantile score, and the `RNA_SEQ` pseudocount makes a track's score depend on its own level. Whether
  effects can be compared across assay kinds is not established; rank with quantiles, and say so.
- *"Where does the 5% come from? How many TPM is that?"* It is this skill's convention, relative to the
  gene's own maximum; there is no TPM equivalent and no source gives an expression threshold.
- *"RNA_SEQ or RNA_SEQ_ACTIVE?"* `RNA_SEQ` is the change between alleles (log); `RNA_SEQ_ACTIVE` is the level
  of the larger allele (linear). Both are means over the gene's exons.
- *"+X% in the target: X% of what?"* Of the gene's predicted expression with the reference allele, averaged in log scale over
  the group's tissues and tracks.
- *"What is the quantile?"* The score's rank among common variants for the same prediction type and track;
  how unusual, not how large; not a p-value and not a %.
- *"What is the typical disagreement?"* `1.4826 × median` of the difference between two tracks of one tissue,
  by size of change; tracks of different assays or sources, not replicates. Older pages called it by another
  name; the glossary lists the names no longer used.
- *"Why this floor value?"* A floor suggested as 5% of the largest REF among the study's tissues (`summarize` prints it
  next to its REF table) (a third maximum, not
  the exploration's catalogue maximum); once set it is an absolute REF value, and any value is the
  researcher's choice, in REF units (no official unit). 0.01 is ten times the pseudocount.
- *"A tissue at 0.1% of the maximum has +300%: does the variant turn the gene on there?"* No: give REF and
  ALT; it remains far below the gene's other tissues.
- *"Are CAGE, DNase and ATAC % comparable with RNA-seq %?"* No: 501 bp local sums in log2 against a gene's
  exon mean in natural log. Use quantiles to rank across them, and say what each measures.
- *"Is the effect specific to my tissue?"* Look at the target − control difference against the typical
  disagreement, the scatter, the per-tissue columns and REF; then say the model's limit on tissue specificity.

## The fixed triage (report and reproduce)

The analysis the researcher sees is the dashboard's. The tier files, the labels and the legacy workbench are
the pipeline's fixed triage, kept for `report.md` (*Resultados da triagem fixa*) and for `reproduce`; their
labels (**Inconsistent**, **Mixed**, the quantile cut) are heuristics, so do not present them as the study's
result.

### Tiers

| Tier | File | Content |
| :--- | :--- | :--- |
| 1 | `processed/master.parquet` | One row per SNV: coordinates, region, AVI (`avi_raw`, `avi_phred`, `top_feature_key`, `fi_*`) and per-track `{prefix}_raw` / `{prefix}_quantile`. No aggregation. |
| 2 | `processed/tissue_contrast.parquet` | Tier 1 plus, per contrast of the approved plan, target and control mean, max and min raw effect and `{contrast}_target_vs_control_diff`, the top tissue per cohort and the priority labels. |

A tissue value is the mean (or median, per the plan) of the plan's tracks for that tissue; cohort
statistics are over tissues. `top_variants.tsv/.json` hold the first `report.top_n` variants
(consistent direction first, then the largest |difference|; a variant without a difference, e.g.
no score in one cohort, comes last).

### Labels

Computed from tissue-level **quantiles** of the plan's `direction_contrast`, which put tissues and assays on
one rank scale; the difference column is in raw scores.

| Label | Meaning |
| :--- | :--- |
| `predicted_direction` | **Increase** / **Decrease** when every target tissue has that sign; otherwise **Inconsistent (within target)**. |
| `target_min_abs_quantile` | The weakest target effect: smallest \|quantile\| across target tissues. A large value means that, in every target tissue, the score is unusual among common variants; the size of the effect is the raw score, not the quantile. |
| `specificity_scope` | Control tissues vote against the target direction: **Shared** (no control goes the other way), **Opposite in control** (most controls go the other way), **Mixed**; **No scope** when the direction is inconsistent, or when no control tissue has a score for the variant (no control data, so no vote). |

Tissue-specific effects are the ones with a consistent direction, a high `target_min_abs_quantile`
and scope **Opposite in control** (or a large `{contrast}_target_vs_control_diff`). **Shared** means
the target change is not specific to the target cohort. Chosen thresholds are yours; none is built
into the labels.

### Direction

`Increase` / `Decrease` refer to the predicted ALT vs REF effect on the direction contrast's assay
in the target tissues (for gene expression, higher or lower predicted RNA), read from the quantiles. The sign
of `{contrast}_target_vs_control_diff` is target minus control in raw scores, so a row labelled `Increase`
can carry a negative difference: read both.

## Scope and limits

- **Raw scores** are comparable within a prediction type, not across them; use quantiles to rank across.
- **Track sources** may differ between groups (e.g. GTEx regions vs ENCODE tissues, polyA vs total
  RNA-seq). The report notes when a contrast of the triage mixes them (`plan_notes`); in the dashboard
  the researcher chooses which assays and sources enter.
- **Window scope:** only the analysed window is mutated; distal enhancers and introns are not.
- **Transcript:** annotations follow the selected transcript; other transcripts that start inside the
  window are listed in `processed/transcripts_in_window.tsv`.
- **Dropped cohorts:** CURIEs without Atlas tracks of their own and prediction types without tracks in a
  group are dropped by decision; `processed/ontology_resolution.json` lists them.
- **`top_track_*` columns** of the master table are null when the AVI region export has no track
  columns (the case of the bundled reference). The workbench does not use them: Figures 3 and 4
  read the per-track `{prefix}_raw` / `{prefix}_quantile` columns through the track registry.
- **AVI Phred:** per the `alphagenome` skill, 10 = top 10% and 20 = top 1% of all SNVs genome-wide; it
  summarises the overall predicted impact and does not say in which tissue. Its scale and the modality
  names were not checked against a primary source (status: not established).

## Legacy workbench (`expression.py dashboard --legacy`)

Figure 1 is the window (footprint, Regulatory Build, transcript). Figure 2 is one lollipop per SNV
along the window (height = chosen metric; click to pick variants; the menu saves SVG or 4x PNG).
Figure 3 is variants by tissues for one contrast. Figure 4 is per-tissue bars for the picked
variants (mean and SD when several).

Defaults when it opens: every genomic domain, direction and scope selected, |quantile| >= 0.99 in
every target tissue, Top-K = all, Y metric = difference of the direction contrast. The `workbench`
stage of `reproduce` counts the variants shown by these defaults (`default_filter`; 713 for the
bundled reference) and the domain filter options in order of position (enhancer, promoter, 5' UTR,
start codon, CDS). Figure 3 shows at most 60 rows by default.

## Reporting

`expression.py report` writes `report.md` (Estudo, Decisões, Resultados, Skills utilizadas,
Limitações, Citações). Include `raw/atlas_links.md` (Atlas deep links, built from
`processed/top_variants.json`) and `raw/interpret_*.txt` (per-variant interpretation) in the study
before running it when you want them in the report; both are optional. When you summarise the report
for the researcher, follow *How to explain a number*: the glossary's names, the four points, the levels
next to large percentages, and the model's limits on tissue specificity.
