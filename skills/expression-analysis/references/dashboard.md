# The analysis dashboard (`assets/dashboard.py`)

How to present the dashboard to the researcher and how to answer questions about it. Read it before
opening the dashboard (Workflow step 14). The page is in Portuguese; quote its labels as written.

## Contents

- [What the page is for](#what-the-page-is-for)
- [Presenting it: what to say first](#presenting-it-what-to-say-first)
- [Layout, top to bottom](#layout-top-to-bottom)
- [Avançado, field by field](#avançado-field-by-field)
- [What the page does without saying](#what-the-page-does-without-saying)
- [Link with the page: `view`](#link-with-the-page-view)
- [Explaining a view: `summarize`](#explaining-a-view-summarize)
- [Figures for the report: `figures`](#figures-for-the-report-figures)
- [What each quantity means](#what-each-quantity-means)
- [Questions researchers ask](#questions-researchers-ask)
- [Limits to say out loud](#limits-to-say-out-loud)

## What the page is for

It answers one question about the study's window: **which variants change the predicted signal in the
target tissues, and is the change their own or shared with the control?** It opens on the study's own
groups with every filter off, so the researcher sees the data before choosing anything (asking for
thresholds first makes them choose blind). Every choice is a control on the page; nothing is fixed in
the code. The page saves no analysis: the only files it touches are its marimo session folder
(`D/.workbench/`), the two files of the link with you (`D/.dashboard/`, next section but one) and one
`14-dashboard` step in `manifest.json` when the session ends (the export adds the HTML). The dashboard never calls the network: a tissue that was not captured needs a new capture
(Workflow steps 7-10), which is your job, not the page's.

## Presenting it: what to say first

Say it in the researcher's language, in a few short lines, when you give the link:

1. The page opens with the study's groups and no filter: every single-base change of the window is
   there, and the blue box gives the result of the current choices.
2. Three choices drive everything: what they look for (decrease, increase, any change), the target
   tissues and the control tissues. Change one and the whole page follows.
3. The row below changes only how results are shown: the order (*Ordenar por*), how many are drawn
   (*Quantas mostrar*, 0 = all) and whether the variants that do not pass appear in grey.
4. Thresholds are in *Filtrar resultados*, all off; technical settings are in *Avançado* (offer to
   explain them when they open it, with the next section).
5. Clicking a variant (in a chart or the table) highlights it everywhere and shows its tracks.
6. You can see what they have on screen and move the page to a view for them; when they settle on a
   view, you draw it for the report (`figures --from-dashboard`).

## Layout, top to bottom

| Part | What it shows | What to tell the researcher |
|---|---|---|
| Title, *Como usar* | gene, window, number of variants (every single-base change) | the change is a percent of what the scorer measures |
| Three choices | *O que você procura?* (diminui / aumenta / qualquer mudança), *Onde* (target tissues), *Comparado com* (control tissues) | these are the only choices needed to start |
| How to show (row below) | *Ordenar por* (order of the table, which variants *Quantas mostrar* keeps, height in the window chart), *Quantas mostrar* (0 = all; otherwise only the first N in that order are drawn and listed), *Variantes nos gráficos* (only those that pass, or all with the rest in grey) | change what is drawn and in which order; the headline and the filter steps still count every variant that passes |
| *De onde vêm as trilhas* | for each chosen tissue, its tracks (GTEx or ENCODE, tissue or cell type, polyA or total RNA-seq) and how well they agree in this gene | tissues whose tracks agree poorly are unstable comparisons |
| Headline (blue box) | how many variants pass, the first three in the order of *Ordenar por* (% in target and in control; the same as the table's first rows), a warning when a tissue left a group for lack of tracks with the chosen assays and sources, and one when the reference level is very low in a target tissue (REF below 0.01, or below the floor chosen) | read it first; it states the result of the current choices |
| *Ao longo da janela* | one stem per variant at its position; height = change in the target (or what *Ordenar por* chose); colour, point size and *Uma alternativa por posição* are next to it | click to highlight; wheel zooms, drag pans, double-click resets; *Recortar* limits the window for the other parts; *Salvar PNG/SVG* below |
| *O efeito é próprio do alvo?* | each variant as a point: change in target (x) vs control (y) | on the dashed diagonal = shared; near the horizontal zero line = own to the target; opposite quadrants = opposite directions |
| *Variantes* | table in % per group and per tissue | tick rows to highlight |
| *Variantes destacadas* | per-track bars for the highlighted variants, Atlas link | bars of one tissue that disagree = unstable tissue |
| *Filtrar resultados* (closed) | minimum change in every target tissue, maximum change in the control, target − control above the between-track disagreement (× the typical one), regions; the count after each step | all start off; the between-track disagreement of the study is the reference. A filter that cannot apply (no control tissue, no tissue with two tracks) is listed as *ignorado* with the reason instead of removing everything |
| *Avançado* (closed) | group names, context tissues, which tracks enter, the expression floor, the captured tissues table, the between-track disagreement table, the prediction type | technical refinements, not for a first look |
| *Comparar todos os tecidos* (closed) | variant × tissue matrix | red = increase, blue = decrease |

## Avançado, field by field

These are closed by default because they change the comparison itself, not just its display. When
the researcher opens *Avançado* or asks about one, explain the field with this, and say what its
current value does:

| Field | What it does | When to change it |
|---|---|---|
| *Nome do grupo alvo / controle* | the names used in titles, legends, the blue box and the report | to use the researcher's own names; nothing is recomputed |
| *Tecidos de contexto (só exibidos)* | tissues shown in the per-tissue table, bars and matrix, but not in the target or control means | to compare with a tissue without making it part of a group |
| *Ensaios que entram* / *Fontes que entram* (empty = all) | which tracks of each tissue count: polyA or total RNA-seq; ENCODE or GTEx | when *De onde vêm as trilhas* says a tissue's tracks agree poorly, or to compare like with like (only GTEx, only polyA); it also changes the between-track disagreement |
| *Usar piso de expressão de referência*, *Piso (REF)*, *O que fazer com o piso* | REF is the gene's predicted expression level in the tissue with the reference allele; below the floor the level is treated as too low to compare. *Só avisar* only warns; *Tirar controles abaixo do piso* removes those control tissues; *Não analisar se o alvo estiver abaixo* lets no variant pass while a target tissue is below it (the filter step *piso de expressão no alvo* shows it) | to avoid comparing against a tissue where the level is very low (a large % there can be a minimal absolute change); keep *Só avisar* when the question is an increase from a very low level |
| *Tecidos capturados* (table) | every tissue with predictions in the study: group at capture, number of tracks, REF min and max, where the tracks come from | to find a CURIE, or a tissue to add to a group; a tissue not listed needs a new capture (your job) |
| *Discordância entre trilhas* (table) | how much two tracks of the same tissue disagree, by size of change: p95, p99 and the typical value (column *típica*). The two tracks differ in assay or source; they are not replicates | the reference for the thresholds; the *típica* column is the unit of the *acima da discordância* filter |
| *Tipo de predição* | RNA-seq (expression of the whole gene), CAGE (transcription start), DNase / ATAC (open chromatin around the variant) | to ask about the mechanism near the variant rather than the gene's expression; groups may lose tissues without that prediction |
| *Pontuação* | the model's value (a change, shown in %) or the Atlas quantile of the same prediction (the variant's rank among common variants, −1 to 1) | quantiles compare kinds of prediction and tracks on one scale; they are not a %, so the % filters are switched off (the filter steps say so) and the question still keeps every target tissue in the direction asked |
| *Valor de cada grupo* | the mean of the group's tissues, or the median | the median weighs less a tissue that disagrees with the rest (for example, one cell line far from the others) |

## What the page does without saying

Mention these when they matter for what the researcher is looking at:

- An empty selection means all (regions, assays, sources).
- *Mudança mínima* applies to **every** target tissue in the direction asked; *mudança máxima no
  controle* applies to every control tissue in both directions.
- *Alvo − controle acima da discordância entre trilhas* follows the question: with *Diminui*, the target must be
  that many times the typical disagreement **below** the control; with *Aumenta*, above; with *Qualquer mudança*, either way. A variant that
  raises the target never passes a repression question through this filter.
- A filter that cannot apply is listed as *ignorado* in the filter steps, with the reason (no control
  tissue, no tissue with two tracks to estimate the disagreement).
- *Quantas mostrar* and *Variantes nos gráficos* do not change the counts: the blue box and the filter
  steps always count every variant that passes.
- *Ordenar por* follows the question for the change in the target and for the two differences: with
  *Aumenta*, the largest rise (or the target furthest above the control) comes first; with *Diminui*, the
  largest fall; with *Qualquer mudança*, the size. A variant whose control changes more than the target goes
  to the end of the list, not to the top. "no sentido procurado, no tecido do alvo que menos muda" ranks
  variants that move all target tissues; "diferença em unidades de discordância" ranks those whose target − control
  difference is largest compared with the disagreement between tracks. The change in the control and
  the AVI Phred are ordered by size.
- A tissue with no track under the chosen assays and sources leaves its group (its column leaves the
  table and the group value is computed without it); the blue box names it. With fewer tissues to agree,
  more variants may pass.
- With every region ticked, *domínios* is not a filter and is not listed in the filter steps.
- *Uma alternativa por posição* keeps, at each position, the base change that comes first in the current
  order.
- *Recortar* limits the window for the scatter and the table too; *Janela inteira* undoes it.
- A tissue with one track has no disagreement of its own; the pooled value of the study is used.
- The dashboard saves no analysis. Reloading the page returns to the last view you sent with `view`, or
  to the opening view when you sent none; the researcher's own changes since then are lost.
- Highlights and the stretch of the window are not part of the choices: a view you send leaves them as they
  are unless it carries `--highlight`, `--clear-highlight`, `--window`, `--full-window` or `--reset`.

## Link with the page: `view`

While the dashboard is open, the page and you share two files in `D/.dashboard/`:

- `researcher_view.json`, written by the page on every change: the choices on screen, the highlights, the
  stretch of the window shown and how many variants pass.
- `agent_view.json`, written by `expression.py view`: the page follows it and moves its controls to the view
  you sent, without reloading. A line *Visão enviada pelo agente* shows your `--message`.

**Reading.** `expression.py view --study-dir D` prints what the researcher has on screen. Use it before
answering anything about the page ("what am I seeing?", "what did I change?"): compare with the opening view
(the study's groups, the decided question, no filter) and say what differs. It is the current state only,
not a history: you cannot tell the order of the changes or what was tried and undone.

**Sending.** Say first what you will change and why, then send it:

```bash
uv run SKILL/scripts/expression.py view --study-dir D --min-sigma 1 --ref-floor 0.78 --floor-applies control --y-metric diff_sigma --message "<what this view is, in the researcher's language>"
```

- Only the options you pass change; the rest stays as the researcher left it. `--reset` starts from the
  study's opening view instead, clears the highlights and shows the whole window.
- The options are those of `summarize` and `figures` (table below), plus `--shown CURIE` (context
  tissues), `--top-k`, `--show kept|all`, `--color-by`, `--highlight VARIANT`, `--window START END`.
- To switch something off: `--no-max-control`, `--no-floor`, `--clear-highlight`, `--no-control`,
  `--no-shown`, `--all-domains`, `--all-assays`, `--all-sources`, `--full-window`.
- Highlights do not move the window; send `--window` when the highlighted variants should be in view.
- A value the page cannot take (a prediction type or a source the study did not capture, an unknown name)
  is ignored and the control keeps its value: read the view back to confirm what was applied.

**What the link does not reach.** The point size, the order and number of rows of the matrix, the page of
the table and which blocks are open stay with the researcher.

**The page must be open.** Each browser tab is its own session. With no tab open, the view you send waits
and is applied when the page is next opened (the page then opens on it, not on the study's opening view);
`view` says when it has never seen an open page. With two tabs open, both follow you and
`researcher_view.json` holds the state of the one that wrote last. The static export does not take part.

## Explaining a view: `summarize`

When the researcher asks about what they see ("why is this variant first?", "is the effect specific?"),
run `expression.py summarize --study-dir D --from-dashboard`: it takes the choices on screen. It runs the dashboard's
computation and prints the headline, the tracks of each tissue, the between-track disagreement, the filter steps and
the top variants. Explain from those numbers; never describe values you have not computed. If they want
a view kept, draw it for the report with `figures` (next section) and write its choices in your answer:
the dashboard saves no analysis.

With the page closed (or for a view that is not on screen), pass the choices yourself:

| On the page | `summarize` option |
|---|---|
| *O que você procura?* Diminui / Aumenta / Qualquer mudança | `--question repression / increase / either` |
| *Onde* (target tissues) | `--target CURIE` once per tissue (the *Tecidos capturados* table in *Avançado* gives the CURIEs) |
| *Comparado com* (control tissues) | `--control CURIE` once per tissue |
| *Mudança mínima em cada tecido do alvo (%)* | `--min-change-pct N` |
| *Exigir pouca mudança no controle* on, *Mudança máxima no controle (%)* | `--max-control-pct N` (omit when off) |
| *Alvo − controle acima da discordância entre trilhas (× a típica)* | `--min-sigma N` |
| *Regiões da janela* | `--domain NAME` once per region kept |
| *Ensaios que entram* / *Fontes que entram* | `--assay`, `--source` (repeatable) |
| *Piso (REF)* on, *O que fazer com o piso* | `--ref-floor N --floor-applies warn/control/target/both` |
| *Uma alternativa por posição* | `--one-alt` |
| *Ordenar por* | `--y-metric` |
| *Quantas mostrar* | nothing: it only limits what is drawn; `--top N` lists the first N |
| *Nome do grupo alvo / controle* | `--target-label`, `--control-label` |
| *Tecidos de contexto (só exibidos)* | `--shown CURIE` once per tissue |
| *Tipo de predição* | `--scorer` |
| *Pontuação* | `--score raw\|quantile` |
| *Valor de cada grupo* | `--group-agg mean\|median` |

Crop, highlights, colours, the point size, *Quantas mostrar* and *Variantes nos gráficos* change only what
is drawn, not the numbers.

## Figures for the report: `figures`

When the researcher settles on a view (or asks for figures for the report), run `expression.py figures
--study-dir D --from-dashboard` (the choices, the way of showing and the highlights on screen; the figures
always draw the whole window, not the stretch shown), or pass the same choices as `summarize` plus how it is
shown: `--top-k` (*Quantas mostrar*),
`--show kept|all` (*Variantes nos gráficos*), `--color-by` (*Cor no gráfico da janela*) and `--highlight
VARIANT` for each highlighted variant (without it, the first five that pass are highlighted). It draws the
page's charts as PNG in `figures/dashboard/` (window, target × control, highlighted variants, variants ×
tissues) and writes `figures/dashboard/choices.json`. Then run `expression.py report` again: `report.md`
gets a section *Análise do pesquisador (dashboard)* with the choices, the blue box, the filter steps and
the figures with captions. Read that section back to the researcher. One view at a time: a new `figures`
run replaces the previous one; the report warns when the figures are older than the variant table.

Before the researcher draws a view, `report` draws the **opening view** itself (the study's groups and
question, no filter, the first 60 variants, all shown) into the same folder, so `report.md` always has
figures, under *Visão inicial do estudo (dashboard)*, plus Figure 1 as an image. A `figures` run with the
researcher's choices replaces it.

**Starting it.** `expression.py dashboard` prints `Dashboard READY: open <url>` only when marimo answers, or
`Dashboard FAILED …` if it exits first. When it runs in the background, wait for one of those lines and give
the link only after READY; on FAILED, say so and show the error.

## What each quantity means

Every quantity on the page has an entry in `references/glossary.md` (name, definition, unit, status). The
scores themselves (formula, sign, window) are defined by the `alphagenome` dependency skill: see *Where each
score is defined* in `references/interpretation.md` before explaining one.

- **Change (%)**: RNA_SEQ is a natural-log fold change of the whole gene's predicted expression;
  CAGE, DNASE and ATAC are log2 ratios of the signal around the variant. The page shows both as a
  percent change of that quantity.
- **Target change**: the mean over the target tissues. The *minimum change* filter applies to every
  target tissue, so a variant that moves only one of them does not pass.
- **Control change**: the mean over the control tissues (the per-tissue columns show if one tissue
  dominates, e.g. a control mean pulled by one tissue).
- **Between-track disagreement** (*discordância entre trilhas*): how much two tracks of the same tissue differ,
  by size of change (p95, p99 and the typical value, 1.4826 × the median difference). The two tracks of a pair
  differ in assay or source: they are not replicates, and the value says nothing about the model's error. A
  target − control difference below it cannot be told from the disagreement between kinds of track. It is
  pooled over every tissue with two or more tracks (with the chosen assays and sources), so it is the study's
  typical disagreement, not that of one tissue; a tissue with one track has none of its own and uses the
  pooled value. With no such tissue the filter is ignored.
- **REF**: the gene's predicted expression level in the tissue with the reference allele, reconstructed from
  the exploration's probe SNV (copied into `raw/` by `study.py init --from-exploration`). It is empty
  when the study has no baseline (an older exploration, or a study not created from one) or the tissue
  has no RNA-seq track. Near 0.001, the pseudocount of the RNA_SEQ score, the level is almost nil, and a large
  % change there can be a minimal absolute change. The headline warns below 0.01 (ten times the pseudocount) or below the
  floor the researcher set in *Avançado*; the floor can also drop control tissues below it, or skip the
  analysis when a target tissue is below it.
- **AVI Phred / modality**: the Atlas overall impact score and its main contributor; *Splicing* means
  the variant probably acts on splicing, not transcription.

## Questions researchers ask

- *"Is this variant specific to my tissue?"* Look at the scatter: near the horizontal zero line and
  far from the centre. Then check the per-tissue columns and the between-track disagreement; a difference within it is
  not evidence of specificity.
- *"Why does the control change so much?"* Often one control tissue drives the mean, or the gene's
  level is very low there (low REF), where a large % can be a minimal absolute change.
- *"Can I trust a change in a tissue where the level is almost nil?"* Only as a predicted increase from a very
  low level, with low confidence; give the REF and ALT levels and say so.
- *"What threshold should I use?"* There is no validated one. Start from the study's between-track disagreement and from the
  size of change that matters for their question; show how the count changes.

## Limits to say out loud

Predicted effects, not measurements; research use only. SNVs only (no indels or rearrangements) inside
the chosen window. The groups compare the tracks the Atlas has; a tissue with one track has no
disagreement of its own. Quantiles are ranks among common variants (per prediction type and per track), not effect sizes. The static export
(`--export-html`) shows the opening view, not the researcher's choices, and has no link with you. As the code
is written (`scripts/viz/variant_canvas.py`), the two interactive charts load vega-embed from `esm.sh` in the
browser, so they may not draw on a computer without internet access; this was not tested.
