# Interpreting the results (expression-analysis)

Read this before reporting results. All numbers come from AlphaGenome Atlas predictions: they
prioritise experiments, they are not measurements, and they are for research use only (never a
diagnosis or a clinical finding).

The analysis the researcher sees is the dashboard's (`references/dashboard.md`): % change per tissue,
target and control groups chosen on the page, the noise between measurements of the same tissue, and
filters that start off. The tier files, the labels and the workbench below are the pipeline's fixed
triage, kept for `report.md` and for `reproduce`; their labels (**Inconsistent**, **Mixed**, the
quantile cut) are heuristics, so do not present them as the study's result.

## Tiers

| Tier | File | Content |
| :--- | :--- | :--- |
| 1 | `processed/master.parquet` | One row per SNV: coordinates, region, AVI (`avi_raw`, `avi_phred`, `top_feature_key`, `fi_*`) and per-track `{prefix}_raw` / `{prefix}_quantile`. No aggregation. |
| 2 | `processed/tissue_contrast.parquet` | Tier 1 plus, per contrast of the approved plan, target and control mean, max and min raw effect and `{contrast}_target_vs_control_diff`, the top tissue per cohort and the priority labels. |

A tissue value is the mean (or median, per the plan) of the plan's tracks for that tissue; cohort
statistics are over tissues. `top_variants.tsv/.json` hold the first `report.top_n` variants
(consistent direction first, then the largest |difference|; a variant without a difference, e.g.
no score in one cohort, comes last).

## Labels (legacy triage)

Computed from tissue-level **quantiles** of the plan's `direction_contrast`; the quantile scale
[-1, 1] is what makes tissues and assays comparable.

| Label | Meaning |
| :--- | :--- |
| `predicted_direction` | **Increase** / **Decrease** when every target tissue has that sign; otherwise **Inconsistent (within target)**. |
| `target_min_abs_quantile` | The weakest target effect: smallest \|quantile\| across target tissues. Large values mean a strong effect in every target tissue. |
| `specificity_scope` | Control tissues vote against the target direction: **Shared** (no control goes the other way), **Opposite in control** (most controls go the other way), **Mixed**; **No scope** when the direction is inconsistent, or when no control tissue has a score for the variant (no control data, so no vote). |

Tissue-specific effects are the ones with a consistent direction, a high `target_min_abs_quantile`
and scope **Opposite in control** (or a large `{contrast}_target_vs_control_diff`). **Shared** means
the target change is not specific to the target cohort. Chosen thresholds are yours; none is built
into the labels.

## Direction

`Increase` / `Decrease` refer to the predicted ALT vs REF effect on the direction contrast's assay
in the target tissues (for gene expression, higher or lower predicted RNA). The sign of
`{contrast}_target_vs_control_diff` is target minus control, so it can differ from the direction of
the target alone: read both.

## Scope and limits

- **Raw scores** are comparable within a scorer, not across scorers; use quantiles to compare assays.
- **Track sources** may differ between cohorts (e.g. GTEx regions vs ENCODE tissues, polyA vs total
  RNA-seq). The plan review points (`contrast.review:*`) record whether the researcher accepted this.
- **Window scope:** only the analysed window is mutated; distal enhancers and introns are not.
- **Transcript:** annotations follow the selected transcript; other transcripts that start inside the
  window are listed in `processed/transcripts_in_window.tsv`.
- **Dropped cohorts:** CURIEs without Atlas tracks of their own and scorers without tracks in a
  cohort are dropped by decision; `processed/ontology_resolution.json` lists them.
- **`top_track_*` columns** of the master table are null when the AVI region export has no track
  columns (the case of the bundled reference). The workbench does not use them: Figures 3 and 4
  read the per-track `{prefix}_raw` / `{prefix}_quantile` columns through the track registry.
- **AVI Phred:** 10 = top 10%, 20 = top 1% of all SNVs genome-wide; it summarises the overall
  predicted impact and does not say in which tissue.

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
before running it when you want them in the report; both are optional.
