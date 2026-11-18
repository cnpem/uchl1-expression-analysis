# Glossary (generated from `assets/terms.yaml`; do not edit by hand)

One name per quantity. Read this before explaining a number to the researcher: say what it is, relative to
what, in what unit, and its status. For an entry marked *undecided*, say that it is not established; do not
infer. The Portuguese name is the one the pages show; use it as written. For the full definition of a score,
read the document named in *Source* (for the AlphaGenome scores, the `alphagenome` dependency skill).

Regenerate after changing the registry: `expression.py terms --write-glossary`. `expression.py terms --check`
fails when a deprecated name is still in the code or the documentation.

## Names

| Name on the pages (Portuguese) | In the documentation (English) | Comes from | Status |
|---|---|---|---|
| trilha | track | Atlas track metadata (`var` of the AnnData; one row of `tracks`) | official (stated by a source) |
| termo (tecido, tipo celular ou linhagem) | biosample term | `ontology_curie` of a track (UBERON, CL, EFO, CLO) | official (stated by a source) |
| ensaio | assay type | `assay_title` of a track | official (stated by a source) |
| tipo de predição | scorer (prediction type) | `requested_scorers` of the Atlas query (RNA_SEQ, RNA_SEQ_ACTIVE, CAGE, DNASE, ATAC) | official (stated by a source) |
| mudança prevista da expressão do gene | predicted expression change (natural-log fold change) | scorer `RNA_SEQ`, raw_score | official (stated by a source) |
| mudança prevista (%) | predicted percent change | core.analysis.pct | this skill's convention |
| nível de expressão previsto | predicted expression level (active allele score) | scorer `RNA_SEQ_ACTIVE`, raw_score | official (stated by a source) |
| variante de sondagem | probe SNV | exploration.json: baseline_variant | this skill's convention |
| nível de expressão previsto com o alelo de referência (REF) | reference expression level (REF) | core.analysis.reference_expression | this skill's convention |
| fração do máximo do gene | share of the gene's maximum | exploration.json: activity.by_curie[*].relative | this skill's convention |
| abaixo de 5% do máximo do gene | below 5% of the gene's maximum | exploration.json: activity.floor (default 0.05) | this skill's convention |
| piso do nível REF | REF-level floor | Profile.ref_floor | this skill's convention |
| pseudocontagem da pontuação RNA_SEQ (0,001) | RNA_SEQ scorer pseudocount (0.001) | scorer `RNA_SEQ`, constant in the formula | official (stated by a source) |
| quantil do Atlas | quantile score | scorer `any`, quantile_score | official (stated by a source) |
| mudança prevista do sinal local | predicted local signal change (log2 fold change, 501 bp) | scorer `CAGE, DNASE, ATAC`, raw_score | official (stated by a source) |
| valor do grupo | group value | Profile.group_agg; columns target_mean, control_mean | this skill's convention |
| diferença alvo − controle | target − control difference | column diff | this skill's convention |
| discordância entre trilhas | between-track disagreement | core.analysis.noise_model | this skill's convention |
| discordância típica | typical disagreement (robust scale) | noise model column sigma; flag --min-sigma; column diff_sigma | this skill's convention |
| aumento previsto a partir de um nível muito baixo | predicted increase from a very low level | n/a | this skill's convention |
| AVI Phred | AVI Phred | scorer `AVI_SCORE`, avi_phred | undecided (no source consulted settles it) |

## Definitions

### trilha (track)

- **What it is:** The predicted signal of one biosample, assay and data source (and strand).
- **Unit:** normalized read signal (RNA-seq, CAGE); normalized insertion signal (DNase, ATAC)
- **Comparable with:** A variant's two alleles within the same track.
- **Comes from:** Atlas track metadata (`var` of the AnnData; one row of `tracks`)
- **Status:** official (stated by a source)
- **Source:** AlphaGenome docs, exploring_model_metadata.md, output types table (stated)
- **Say with it:** A track is a prediction, not an experiment. Two tracks of one term differ in assay or source; they are not replicates.
- **Names no longer used for it:** "medição"; "medições"; "measurements of the same tissue"; "tissue with one measurement"; "two measurements"

### termo (tecido, tipo celular ou linhagem) (biosample term)

- **What it is:** The ontology term a track belongs to. It can be a tissue, a primary cell, a cell line or an in-vitro differentiated cell.
- **Unit:** none
- **Comparable with:** n/a
- **Comes from:** `ontology_curie` of a track (UBERON, CL, EFO, CLO)
- **Status:** official (stated by a source)
- **Source:** AlphaGenome docs, exploring_model_metadata.md (stated)
- **Say with it:** The pages say "tecido" for short; a term may be a cell line.

### ensaio (assay type)

- **What it is:** The kind of experiment a track was trained on, such as polyA plus RNA-seq or total RNA-seq.
- **Unit:** none
- **Comparable with:** n/a
- **Comes from:** `assay_title` of a track
- **Status:** official (stated by a source)
- **Source:** AlphaGenome docs, exploring_model_metadata.md (stated)
- **Say with it:** Tracks of different assay types or sources are different kinds of track. Whether their levels can be compared is not established (see the open points); do not compare them.

### tipo de predição (scorer (prediction type))

- **What it is:** The rule that turns the predicted tracks of the two alleles into one score per variant.
- **Unit:** none
- **Comparable with:** n/a
- **Comes from:** `requested_scorers` of the Atlas query (RNA_SEQ, RNA_SEQ_ACTIVE, CAGE, DNASE, ATAC)
- **Status:** official (stated by a source)
- **Source:** K-Dense alphagenome skill, references/model-api.md, Recommended scorers (stated)
- **Names no longer used for it:** "tipo de sinal"; "tipos de sinal"

### mudança prevista da expressão do gene (predicted expression change (natural-log fold change))

- **What it is:** How much the alternative allele changes the gene's predicted RNA-seq signal, compared with the reference allele, in one track.
- **Formula:** `ln((mean_exons(ALT) + 0.001) / (mean_exons(REF) + 0.001))`
- **Unit:** none (natural log of a ratio)
- **Comparable with:** Variants within the same track; not comparable with CAGE, DNase or ATAC scores.
- **Comes from:** scorer `RNA_SEQ`, raw_score
- **Status:** official (stated by a source)
- **Source:** K-Dense alphagenome skill, references/model-api.md, Recommended scorers; AlphaGenome docs, variant_scoring.md (stated)

### mudança prevista (%) (predicted percent change)

- **What it is:** The score written as a percent change of the signal it measures.
- **Formula:** `100 * (e^score - 1) for RNA_SEQ; 100 * (2^score - 1) for CAGE, DNASE and ATAC`
- **Unit:** %
- **Comparable with:** Within one prediction type.
- **Comes from:** core.analysis.pct
- **Status:** this skill's convention
- **Source:** arithmetic on the scorer formula (scripts/core/analysis.py) (inferred)
- **Say with it:** For a tissue or a group, the percent is computed on the mean of the log scores, so each tissue weighs the same.

### nível de expressão previsto (predicted expression level (active allele score))

- **What it is:** The larger of the two alleles' mean predicted RNA-seq coverage over the gene's exons, in one track.
- **Formula:** `max(mean_exons(ALT), mean_exons(REF))`
- **Unit:** no official unit (normalized read signal; not TPM)
- **Comparable with:** The same gene across variants in one track, and across terms within one kind of track.
- **Comes from:** scorer `RNA_SEQ_ACTIVE`, raw_score
- **Status:** official (stated by a source)
- **Source:** K-Dense alphagenome skill, references/model-api.md, Recommended scorers; AlphaGenome docs, variant_scoring.md, Active-allele scorers (stated)
- **Say with it:** It is a level, not a change, and exists only per variant. It equals the reference level when the variant does not raise expression, and is larger otherwise.
- **Names no longer used for it:** "atividade do gene"; "atividade prevista"; "atividade basal"; "atividade da variante"; "mais ativo"; "gene activity"; "baseline activity"; "predicted activity"; "the gene's activity"; "relative activity"; "most active"

### variante de sondagem (probe SNV)

- **What it is:** The single SNV in the middle of exon 1 that the exploration scores only to read the gene's level in every track.
- **Unit:** none
- **Comparable with:** n/a
- **Comes from:** exploration.json: baseline_variant
- **Status:** this skill's convention
- **Source:** scripts/explore.py, baseline_variant (stated)
- **Say with it:** It is arbitrary. The exploration does not measure how much it changes expression.
- **Names no longer used for it:** "variante de referência"

### nível de expressão previsto com o alelo de referência (REF) (reference expression level (REF))

- **What it is:** The predicted expression level of the reference allele in a track, reconstructed from the two scores of the probe SNV.
- **Formula:** `REF = ACTIVE when LFC <= 0; otherwise (ACTIVE + 0.001) * exp(-LFC) - 0.001`
- **Unit:** no official unit (same scale as the expression level)
- **Comparable with:** As the expression level.
- **Comes from:** core.analysis.reference_expression
- **Status:** this skill's convention
- **Source:** algebra on the two scorer formulas (scripts/core/analysis.py) (inferred)
- **Say with it:** REF also names the reference base of a variant; say "nível REF" when the level is meant.

### fração do máximo do gene (share of the gene's maximum)

- **What it is:** A term's mean expression level divided by the largest term mean of the gene in the catalogue.
- **Unit:** none
- **Comparable with:** Terms whose tracks are of the same kind; the catalogue-wide maximum mixes kinds of track.
- **Comes from:** exploration.json: activity.by_curie[*].relative
- **Status:** this skill's convention
- **Source:** scripts/explore.py, baseline_activity (stated)
- **Say with it:** The mean of a term may mix assays and sources. Read the assay and source columns before comparing two terms.

### abaixo de 5% do máximo do gene (below 5% of the gene's maximum)

- **What it is:** The exploration's flag for a term whose share of the maximum is under the cut.
- **Unit:** none
- **Comparable with:** n/a
- **Comes from:** exploration.json: activity.floor (default 0.05)
- **Status:** this skill's convention
- **Source:** scripts/explore.py, DEFAULT_FLOOR (stated)
- **Say with it:** Not calibrated, and no source gives an "expressed" threshold. It does not mean the transcript is absent.
- **Names no longer used for it:** "não expresso"; "praticamente desligado"; "gene desligado"; "sem expressão"; "quase desligado"; "not expressed"; "gene is silent"; "practically off"; "almost off"

### piso do nível REF (REF-level floor)

- **What it is:** A reference level chosen by the researcher below which a tissue is warned about, or left out of the control. A value may be suggested as 5% of the largest REF among the study's tissues; once set, it is an absolute level.
- **Unit:** no official unit (same scale as the expression level)
- **Comparable with:** n/a
- **Comes from:** Profile.ref_floor
- **Status:** this skill's convention
- **Source:** scripts/core/analysis.py, expression_flags (stated)
- **Say with it:** It is an absolute level, unlike the exploration's cut, which is a share of the gene's maximum over the whole catalogue. A suggestion computed from the study's tissues uses a third maximum, that of the study.

### pseudocontagem da pontuação RNA_SEQ (0,001) (RNA_SEQ scorer pseudocount (0.001))

- **What it is:** The constant added to both exon means in the RNA_SEQ log fold change.
- **Unit:** same scale as the expression level
- **Comparable with:** n/a
- **Comes from:** scorer `RNA_SEQ`, constant in the formula
- **Status:** official (stated by a source)
- **Source:** K-Dense alphagenome skill, references/model-api.md, Recommended scorers (stated)
- **Names no longer used for it:** "pseudocount do modelo"

### quantil do Atlas (quantile score)

- **What it is:** The rank of a raw score among common variants, for the same prediction type and the same track.
- **Unit:** none (-1 to 1 for signed scores, 0 to 1 otherwise)
- **Comparable with:** Says how unusual a score is; the size of the effect is the raw score.
- **Comes from:** scorer `any`, quantile_score
- **Status:** official (stated by a source)
- **Source:** AlphaGenome docs, faqs.md; K-Dense alphagenome skill, references/atlas.md, Raw, quantile, Phred (stated)
- **Say with it:** It is not a percentage and not a p-value. For RNA_SEQ_ACTIVE it ranks the level, not an effect.
- **Names no longer used for it:** "ranks against the genome"

### mudança prevista do sinal local (predicted local signal change (log2 fold change, 501 bp))

- **What it is:** How much the alternative allele changes the predicted signal summed over 501 bp centred on the variant.
- **Formula:** `log2((sum_501bp(ALT) + 1) / (sum_501bp(REF) + 1))`
- **Unit:** none (log2 of a ratio)
- **Comparable with:** Variants within the same track; it belongs to no gene and is not comparable with RNA_SEQ scores.
- **Comes from:** scorer `CAGE, DNASE, ATAC`, raw_score
- **Status:** official (stated by a source)
- **Source:** K-Dense alphagenome skill, references/model-api.md, Recommended scorers (stated)

### valor do grupo (group value)

- **What it is:** The mean (or median) over a group's tissues of each tissue's score; a tissue's score is the mean of its tracks.
- **Unit:** that of the score
- **Comparable with:** Target against control, within one prediction type.
- **Comes from:** Profile.group_agg; columns target_mean, control_mean
- **Status:** this skill's convention
- **Source:** scripts/core/analysis.py, compute (stated)

### diferença alvo − controle (target − control difference)

- **What it is:** The target group value minus the control group value.
- **Unit:** that of the score
- **Comparable with:** Within one prediction type.
- **Comes from:** column diff
- **Status:** this skill's convention
- **Source:** scripts/core/analysis.py, compute (stated)

### discordância entre trilhas (between-track disagreement)

- **What it is:** The absolute difference between the scores of two tracks of the same term, for the same variant.
- **Unit:** that of the score
- **Comparable with:** n/a
- **Comes from:** core.analysis.noise_model
- **Status:** this skill's convention
- **Source:** scripts/core/analysis.py, _pair_differences (stated)
- **Say with it:** The two tracks of a pair differ in assay or source; this is not replicate noise and says nothing about the model's error.
- **Names no longer used for it:** "ruído"; "ruídos"; "noise between measurements"; "noise of the study"

### discordância típica (typical disagreement (robust scale))

- **What it is:** The usual size of the between-track disagreement, by size of change.
- **Formula:** `1.4826 * median(|a - b|) within each of five bins of the mean absolute score`
- **Unit:** that of the score
- **Comparable with:** n/a
- **Comes from:** noise model column sigma; flag --min-sigma; column diff_sigma
- **Status:** this skill's convention
- **Source:** scripts/core/analysis.py, noise_model (stated)
- **Say with it:** The filter and the ordering divide the target − control difference by it. The code keeps the names sigma, --min-sigma and diff_sigma.
- **Names no longer used for it:** "σ"

### aumento previsto a partir de um nível muito baixo (predicted increase from a very low level)

- **What it is:** A positive expression change in a tissue whose reference level is close to the scorer pseudocount.
- **Unit:** n/a
- **Comparable with:** n/a
- **Comes from:** n/a
- **Status:** this skill's convention
- **Source:** this skill (inferred)
- **Say with it:** Give the reference and the alternative levels with it; a large percentage there can be a minimal absolute change.
- **Names no longer used for it:** "ativação ectópica"; "ectopic activation"

### AVI Phred (AVI Phred)

- **What it is:** The Atlas overall variant-impact score on a Phred-like scale.
- **Unit:** Phred-like
- **Comparable with:** Variants genome-wide; it does not say in which tissue.
- **Comes from:** scorer `AVI_SCORE`, avi_phred
- **Status:** undecided (no source consulted settles it)
- **Source:** K-Dense alphagenome skill, references/atlas.md, AlphaGenome Variant Impact (AVI); the Atlas report was not read (not found)
- **Say with it:** The scale and the meaning of each modality were not checked against a primary source.

## Writing conventions

How names and numbers are written on the pages; write them the same way in your answers.

| What | Rule | Example | Source |
|---|---|---|---|
| Human gene | The HGNC-approved symbol, in italics. Ask the researcher when they give a disease, a protein or an alias, and confirm the symbol before the exploration. | *TP53*; for phenylketonuria (PKU) the gene is *PAH* | HGNC guidelines for human gene nomenclature |
| Protein | The same letters as the gene, in roman type, or the protein's name; the skill's numbers are about the gene, not the protein. | p53 (protein of *TP53*); PAH (protein of *PAH*) | HGNC guidelines for human gene nomenclature |
| Mouse or rat gene | Only the first letter in upper case, in italics; the protein in upper case, roman. The skill analyses human genes only. | *Uchl1* | MGI and RGD nomenclature guidelines |
| Variant | The Atlas identifier, which every command and link uses, with the genomic HGVS description on the GRCh38 RefSeq chromosome next to it where a variant is listed. Do not call the Atlas identifier HGVS. | chr17:7676154:G>C (NC_000017.11:g.7676154G>C) | HGVS nomenclature; NCBI GRCh38.p14 assembly report for the chromosome accessions |
| Transcript, gene and region identifiers, ontology terms | As issued (Ensembl, Regulatory Build, CURIE), in roman or as code, never translated. | ENST00000284440.9, ENSR4_837Q9H, CL:0000100 | this skill |
| Genomic coordinates | No thousands separator. | chr17:7676001-7676800 | this skill |
| Numbers in Portuguese text | Decimal comma, and a dot every three digits for counts. | 1.234 variantes; nível 12,5; +5,2% | this skill |
| Species | The scientific name in italics. | *Homo sapiens* | binomial nomenclature |

## Not decided by any source consulted

State these as open when the researcher asks; a rewording does not settle them.

- The exact unit of the expression level (RNA_SEQ_ACTIVE) and its numeric relation to RPKM (it is not TPM, which is normalised across genes).
- Whether expression levels of different assay types and sources can be compared.
- Any threshold above which a gene counts as expressed.
- Whether a low level makes predicted changes unreliable (the pseudocount alone shrinks them).
- Whether the Atlas computes exactly the AlphaGenome 1 scorers with the distilled model the alphagenome skill describes (not checked against the Atlas report).
- Whether one reference level applies to every variant of a window.
- Whether the Atlas's quantile background is exactly the one the AlphaGenome documentation describes (common gnomAD variants), and what a mean of quantiles over tracks means.
- The AVI Phred scale and the AVI modality names.
