# The analysis dashboard (`assets/dashboard.py`)

How to present the dashboard to the researcher and how to answer questions about it. Read it before
opening the dashboard (Workflow step 14). The page is in Portuguese; quote its labels as written.

## Contents

- [What the page is for](#what-the-page-is-for)
- [Layout, top to bottom](#layout-top-to-bottom)
- [Explaining a view: `summarize`](#explaining-a-view-summarize)
- [What each quantity means](#what-each-quantity-means)
- [Questions researchers ask](#questions-researchers-ask)
- [Limits to say out loud](#limits-to-say-out-loud)

## What the page is for

It answers one question about the study's window: **which variants change the predicted signal in the
target tissues, and is the change their own or shared with the control?** It opens on the study's own
groups with every filter off, so the researcher sees the data before choosing anything (asking for
thresholds first makes them choose blind). Every choice is a control on the page; nothing is fixed in
the code. The page saves no analysis: the only files it touches are its marimo session folder
(`D/.workbench/`) and one `14-dashboard` step in `manifest.json` when the session ends (the export adds
the HTML). The dashboard never calls the network: a tissue that was not captured needs a new capture
(Workflow steps 7-10), which is your job, not the page's.

## Layout, top to bottom

| Part | What it shows | What to tell the researcher |
|---|---|---|
| Title, *Como usar* | gene, window, number of variants (every single-base change) | the change is a percent of what the scorer measures |
| Three choices | *O que você procura?* (diminui / aumenta / qualquer mudança), *Onde* (target tissues), *Comparado com* (control tissues) | these are the only choices needed to start |
| How to show (row below) | *Ordenar por* (order of the table, which variants *Quantas mostrar* keeps, height in the window chart), *Quantas mostrar* (0 = all; otherwise only the first N in that order are drawn and listed), *Variantes nos gráficos* (only those that pass, or all with the rest in grey) | change what is drawn and in which order; the headline and the filter steps still count every variant that passes |
| *De onde vêm as medições* | for each chosen tissue, its measurements (GTEx or ENCODE, tissue or cell type, polyA or total RNA-seq) and how well they agree in this gene | tissues whose measurements agree poorly are unstable comparisons |
| Headline (blue box) | how many variants pass, the ones that change the target most (% in target and in control), a warning when the gene is predicted to be almost off in a target tissue | read it first; it states the result of the current choices |
| *Ao longo da janela* | one stem per variant at its position; height = change in the target (or what *Ordenar por* chose); colour, point size and *Uma alternativa por posição* are next to it | click to highlight; wheel zooms, drag pans, double-click resets; *Recortar* limits the window for the other parts; *Salvar PNG/SVG* below |
| *O efeito é próprio do alvo?* | each variant as a point: change in target (x) vs control (y) | on the dashed diagonal = shared; near the horizontal zero line = own to the target; opposite quadrants = opposite directions |
| *Variantes* | table in % per group and per tissue | tick rows to highlight |
| *Variantes destacadas* | per-measurement bars for the highlighted variants, Atlas link | bars of one tissue that disagree = unstable tissue |
| *Filtrar resultados* (closed) | minimum change in every target tissue, maximum change in the control, target − control above the noise (× σ), regions; the count after each step | all start off; the noise of the study is the reference. A filter that cannot apply (no control tissue, no tissue with two measurements) is listed as *ignorado* with the reason instead of removing everything |
| *Avançado* (closed) | group names, context tissues, which measurements enter, the expression floor, the captured tissues table, the noise table, the prediction type | technical refinements, not for a first look |
| *Comparar todos os tecidos* (closed) | variant × tissue matrix | red = increase, blue = decrease |

## Explaining a view: `summarize`

When the researcher asks about what they see ("why is this variant first?", "is the effect specific?"),
run `expression.py summarize --study-dir D` with the same choices they set. It runs the dashboard's
computation and prints the headline, the measurements of each tissue, the noise, the filter steps and
the top variants. Explain from those numbers; never describe values you have not computed. If they want
a view kept, write its choices in your answer: the dashboard saves nothing.

You cannot read the page's state, so ask for the choices that differ from the opening view (or have
the researcher read them out) and map them:

| On the page | `summarize` option |
|---|---|
| *O que você procura?* Diminui / Aumenta / Qualquer mudança | `--question repression / increase / either` |
| *Onde* (target tissues) | `--target CURIE` once per tissue (the *Tecidos capturados* table in *Avançado* gives the CURIEs) |
| *Comparado com* (control tissues) | `--control CURIE` once per tissue |
| *Mudança mínima em cada tecido do alvo (%)* | `--min-change-pct N` |
| *Exigir pouca mudança no controle* on, *Mudança máxima no controle (%)* | `--max-control-pct N` (omit when off) |
| *Alvo − controle acima do ruído (× σ)* | `--min-sigma N` |
| *Regiões da janela* | `--domain NAME` once per region kept |
| *Ensaios que entram* / *Fontes que entram* | `--assay`, `--source` (repeatable) |
| *Piso (REF)* on, *O que fazer com o piso* | `--ref-floor N --floor-applies warn/control/target/both` |
| *Uma alternativa por posição* | `--one-alt` |
| *Ordenar por* | `--y-metric` |
| *Quantas mostrar* | nothing: it only limits what is drawn; `--top N` lists the first N |
| *Nome do grupo alvo / controle* | `--target-label`, `--control-label` |
| *Tipo de predição* | `--scorer` |

Crop, highlights, colours, the point size, *Quantas mostrar* and *Variantes nos gráficos* change only what
is drawn, not the numbers.

## What each quantity means

- **Change (%)**: RNA_SEQ is a natural-log fold change of the whole gene's predicted expression;
  CAGE, DNASE and ATAC are log2 ratios of the signal around the variant. The page shows both as a
  percent change of that quantity.
- **Target change**: the mean over the target tissues. The *minimum change* filter applies to every
  target tissue, so a variant that moves only one of them does not pass.
- **Control change**: the mean over the control tissues (the per-tissue columns show if one tissue
  dominates, e.g. a control mean pulled by one tissue).
- **Noise**: how much two measurements of the same tissue differ, by size of change (p95, p99, σ).
  A target − control difference below it cannot be told from the disagreement between measurements.
  It is pooled over every tissue with two or more measurements (with the chosen assays and sources), so
  it is the study's typical disagreement, not that of one tissue; a tissue with one measurement has none
  of its own and uses the pooled value. With no such tissue the σ filter is ignored.
- **REF**: the gene's predicted expression in the tissue without the variant, reconstructed from the
  exploration's baseline variant (copied into `raw/` by `study.py init --from-exploration`). It is empty
  when the study has no baseline (an older exploration, or a study not created from one) or the tissue
  has no RNA-seq track. Near 0.001, the model's pseudocount, the gene is practically off and predicted
  changes there are unstable. The headline warns below 0.01 (ten times the pseudocount) or below the
  floor the researcher set in *Avançado*; the floor can also drop control tissues below it, or skip the
  analysis when a target tissue is below it.
- **AVI Phred / modality**: the Atlas overall impact score and its main contributor; *Splicing* means
  the variant probably acts on splicing, not transcription.

## Questions researchers ask

- *"Is this variant specific to my tissue?"* Look at the scatter: near the horizontal zero line and
  far from the centre. Then check the per-tissue columns and the noise; a difference within the noise is
  not evidence of specificity.
- *"Why does the control change so much?"* Often one control tissue drives the mean, or the gene is
  barely expressed there (low REF), where predicted fold changes are larger and unstable.
- *"Can I trust a change in a tissue where the gene is off?"* Only as a possible ectopic activation,
  and with low confidence; say so.
- *"What threshold should I use?"* There is no validated one. Start from the study's noise and from the
  size of change that matters for their question; show how the count changes.

## Limits to say out loud

Predicted effects, not measurements; research use only. SNVs only (no indels or rearrangements) inside
the chosen window. The groups compare the tracks the Atlas has; a tissue with one measurement has no
noise of its own. Quantiles are ranks against the genome, not effect sizes, and are not used here. The
static export (`--export-html`) shows the opening view, not the researcher's choices.
