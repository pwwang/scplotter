# Visualize gene set enrichment and over-representation analysis results

Gene set enrichment analysis identifies biological pathways, gene
ontologies, or functional categories that are statistically
over-represented among a list of genes of interest (e.g., differentially
expressed genes from a single-cell RNA-seq experiment). Rather than
interpreting individual genes in isolation, enrichment analysis places
gene-level results into a broader biological context, revealing which
processes, functions, or diseases are perturbed.

`EnrichmentPlot` generates publication-quality visualizations for
enrichment results across eight distinct plot types, each suited to a
different analytical perspective:

- **bar** — Horizontal bar chart of the top enriched terms, ordered by
  significance. Best for a quick overview or when showing a small number
  of terms.

- **dot** — Dot plot where x-axis shows a continuous metric (default:
  `GeneRatio`), dot size reflects gene count, and dot color reflects
  significance. Ideal for comparing terms along two dimensions
  simultaneously.

- **lollipop** — Lollipop chart combining dot and bar aesthetics.
  Similar to the dot plot but with stems emphasizing the ranking.

- **comparison** — Side-by-side dot plot comparing enrichment across
  groups (e.g., cell types, conditions). Requires `group_by`.

- **network** — Network visualization where nodes are enriched terms and
  edges represent overlapping gene sets. Reveals functional modules and
  redundant terms.

- **enrichmap** — Enrichment map similar to the network plot but
  optimized for large term sets (default `top_term = 100`). Nodes are
  terms and edges represent gene overlap.

- **wordcloud** — Word cloud where term size reflects significance. Can
  display either enrichment terms (`word_type = "term"`) or individual
  gene symbols (`word_type = "feature"`).

- **heatmap** — Heatmap of enrichment significance across groups
  (`group_by` is mapped to columns). Useful for comparing enrichment
  patterns across multiple conditions or cell types.

The function auto-detects the input data format (clusterProfiler or
enrichR) and delegates visualization to the appropriate plotthis
plotting function.

## Usage

``` r
EnrichmentPlot(
  data,
  top_term = NULL,
  plot_type = c("bar", "dot", "lollipop", "network", "enrichmap", "wordcloud",
    "comparison", "heatmap"),
  x_by = NULL,
  size_by = NULL,
  fill_cutoff_name = NULL,
  fill_name = NULL,
  values_fill = 0,
  character_width = 50,
  expand = NULL,
  word_type = c("term", "feature"),
  split_by = NULL,
  split_by_sep = "_",
  facet_by = NULL,
  facet_scales = NULL,
  group_by = NULL,
  group_by_sep = "_",
  metric = "p.adjust",
  cutoff = NULL,
  palette = "Spectral",
  xlab = NULL,
  ylab = NULL,
  ...
)
```

## Arguments

- data:

  A data frame with enrichment results. Must be the output of a
  clusterProfiler function (`enrichGO`, `enrichKEGG`, `enrichPathway`,
  `enrichWP`, etc.) or an enrichR result processed through
  [`plotthis::prepare_enrichr_result()`](https://pwwang.github.io/plotthis/reference/prepare_enrichr_result.html).
  The function auto-detects the format based on column names.

- top_term:

  Integer. Number of top terms to display per group/facet combination,
  ranked by the coloring metric: smallest values for p-value columns,
  largest values otherwise. Default: `6` for all plot types except
  `"enrichmap"` which defaults to `100`. Note that terms are not
  filtered globally — the top terms are selected independently within
  each combination of `split_by`, `group_by`, and `facet_by` levels.

- plot_type:

  Character. The type of plot to generate. One of: `"bar"`, `"dot"`,
  `"lollipop"`, `"network"`, `"enrichmap"`, `"wordcloud"`,
  `"comparison"`, or `"heatmap"`. See the Description section for
  guidance on choosing a plot type. Default: `"bar"`.

- x_by:

  Character. Column name(s) to use for the x-axis. Works only for
  `"dot"` and `"lollipop"` plot types. Default: `NULL` (defaults to
  `"GeneRatio"` internally).

- size_by:

  Character. Column name(s) to map to point size. Works only for
  `"comparison"`, `"dot"`, and `"lollipop"` plot types. Default: `NULL`
  (defaults to `"GeneRatio"` for comparison, `"Count"` for dot and
  lollipop).

- fill_cutoff_name:

  Character. Legend label for terms that exceed the `cutoff` (shown in
  gray). Applies to `"comparison"`, `"dot"`, and `"lollipop"` plot
  types. Default: `NULL` (defaults to `"Non-significant"` when `cutoff`
  is set).

- fill_name:

  Character. Legend title for the fill color scale (the coloring
  metric). Applies to `"comparison"`, `"dot"`, and `"lollipop"` plot
  types. Default: `NULL` (auto-generated from `metric`, e.g.
  `"-log10(p.adjust)"` or `"zScore"`).

- values_fill:

  Numeric. The fill value for missing entries in the heatmap matrix.
  Used only for `"heatmap"` plot type. Default: `0`.

- character_width:

  Integer. Maximum character width for term descriptions before
  line-wrapping. Applies to all plot types; for `"heatmap"` the wrapping
  is deferred to the Heatmap function. Default: `50`.

- expand:

  Numeric vector of length 1, 2, or 4. Axis expansion factors passed to
  [`plotthis::BarPlot()`](https://pwwang.github.io/plotthis/reference/barplot.html).
  Used only for `"bar"` plot type. Default: `NULL` (defaults to
  `c(0.1, 0.6, 0, 0.6)`).

- word_type:

  Character. What to display in the wordcloud. One of `"term"`
  (enrichment term descriptions) or `"feature"` (gene symbols from the
  enriched gene list). Used only for `"wordcloud"` plot type. Default:
  `"term"`.

- split_by:

  Character vector. Column name(s) in `data` to split the data and
  generate separate plots for each unique value. Multiple columns are
  concatenated with `split_by_sep`. Default: `NULL`.

- split_by_sep:

  Character. Separator used when concatenating multiple `split_by`
  columns. Default: `"_"`.

- facet_by:

  Character vector. Column name(s) in `data` to use for faceting
  (generating sub-panels within each plot). Default: `NULL`.

- facet_scales:

  Character. Facet scale behavior — `"fixed"` (same scales), `"free"`,
  `"free_x"`, or `"free_y"`. Default: `NULL` (defaults to `"free_y"` for
  bar, dot, lollipop, and comparison plots).

- group_by:

  Character vector. Column name(s) in `data` to group terms. Behavior
  depends on `plot_type`:

  - `"comparison"` — Groups are shown as x-axis categories in a dot plot
    comparing enrichment across groups. Required for this type.

  - `"heatmap"` — Groups are used as the columns of the heatmap (mapped
    to `columns_by` in
    [`plotthis::Heatmap()`](https://pwwang.github.io/plotthis/reference/Heatmap.html)).

  - All other types — `group_by` is not supported and will raise an
    error. Use `facet_by` or `split_by` instead.

  Multiple columns are concatenated with `group_by_sep`. Default:
  `NULL`.

- group_by_sep:

  Character. Separator used when concatenating multiple `group_by`
  columns. Used only for `"comparison"` plot type. Default: `"_"`.

- metric:

  Character. The column name in `data` to use for coloring terms, or two
  column names separated by `":"` to use one column for coloring and
  another for significance. `"p.adjust"` (default), `"pvalue"`, and
  `"qvalue"` are transformed with \\-log\_{10}\\ so that more
  significant terms have higher values; other columns are used as-is.
  With `"zScore:p.adjust"`, cells are colored by `zScore` while `cutoff`
  is applied to `p.adjust`. A single column is shorthand for
  `"column:column"`. See the *Metric transformation* section.

- cutoff:

  Numeric. A significance threshold applied to the significance metric
  (the column after `":"` in `metric`, or `metric` itself when no
  modifier is given). Default: `NULL` (no marking). When `metric`
  supplies a separate significance metric, a term failing the cutoff in
  every group is dropped before `top_term` selection. A term that passes
  in at least one group is kept in all of them, so its cells are not
  left empty. Dot, lollipop, and comparison plots are the exception:
  they show insignificant terms in gray, so nothing is dropped. The
  behavior on the plot itself depends on `plot_type`:

  - `"bar"` — Adds a vertical dashed line at the transformed cutoff
    (e.g., \\-log\_{10}(0.05)\\). No line is drawn when the significance
    metric differs from the coloring metric, since it would sit on a
    different scale; terms are still filtered.

  - `"dot"`, `"lollipop"`, `"comparison"` — Terms that fail the cutoff
    are colored gray with the legend label from `fill_cutoff_name`,
    while fill color still follows the coloring metric.

  - `"heatmap"` — Adds asterisk (`*`) labels to cells whose significance
    metric passes the cutoff, plus a legend entry labelled
    `"<significance metric> < <cutoff>"` (legend title from
    `label_name`, default `"Significance"`).

  - `"network"`, `"enrichmap"`, `"wordcloud"` — Nothing is marked on the
    plot; terms are still filtered when a significance metric is given.

  Without a separate significance metric, this parameter only marks
  terms — it does not filter them. Use `top_term` to control how many
  terms are shown.

- palette:

  Character. Color palette name for the fill scale. See
  [`plotthis::show_palettes()`](https://pwwang.github.io/plotthis/reference/show_palettes.html)
  for available palettes. Default: `"Spectral"`.

- xlab:

  Character. Custom x-axis label. Default: `NULL` (auto-generated based
  on plot type and `x_by`/`metric`).

- ylab:

  Character. Custom y-axis label. Default: `NULL` (auto-generated based
  on plot type).

- ...:

  Additional arguments passed to the underlying plotthis plotting
  function, determined by `plot_type`:

  `"bar"`

  :   [`plotthis::BarPlot()`](https://pwwang.github.io/plotthis/reference/barplot.html)

  `"dot"`

  :   [`plotthis::DotPlot()`](https://pwwang.github.io/plotthis/reference/dotplot.html)

  `"lollipop"`

  :   [`plotthis::LollipopPlot()`](https://pwwang.github.io/plotthis/reference/dotplot.html)

  `"network"`

  :   [`plotthis::EnrichNetwork()`](https://pwwang.github.io/plotthis/reference/enrichmap1.html)

  `"enrichmap"`

  :   [`plotthis::EnrichMap()`](https://pwwang.github.io/plotthis/reference/enrichmap1.html)

  `"wordcloud"`

  :   [`plotthis::WordCloudPlot()`](https://pwwang.github.io/plotthis/reference/WordCloudPlot.html)

  `"comparison"`

  :   [`plotthis::DotPlot()`](https://pwwang.github.io/plotthis/reference/dotplot.html)

  `"heatmap"`

  :   [`plotthis::Heatmap()`](https://pwwang.github.io/plotthis/reference/Heatmap.html)

## Value

A ggplot object (or a `patchwork` object when `split_by` generates
multiple plots and `combine = TRUE`), or a list of ggplot objects if
`combine = FALSE`. The specific return type depends on the underlying
plotthis function dispatched by `plot_type`.

## Note

- The function auto-detects clusterProfiler vs enrichR input format. For
  enrichR input, it must contain `P.value` and `Adjusted.P.value`
  columns. If your enrichR results came from a different pipeline,
  pre-process them with
  [`plotthis::prepare_enrichr_result()`](https://pwwang.github.io/plotthis/reference/prepare_enrichr_result.html).

- The `cutoff` parameter only **marks** terms — it does not filter them.
  To reduce the number of displayed terms, use `top_term`.

- `GeneRatio` strings (e.g., `"38/225"`) and `BgRatio` strings are
  automatically converted to numeric values by dividing the numerator by
  the denominator.

- When using `group_by` with `plot_type = "comparison"`, `size_by`
  defaults to `"GeneRatio"` and each group's terms are shown
  side-by-side. For `plot_type = "heatmap"`, `group_by` becomes the
  heatmap columns.

- `group_by` is not supported for `"bar"`, `"dot"`, `"lollipop"`,
  `"network"`, `"enrichmap"`, and `"wordcloud"` — use `facet_by` or
  `split_by` to separate groups for those types.

- For `plot_type = "wordcloud"` with `word_type = "feature"`, individual
  gene symbols are extracted from the `geneID` column. Gene-level
  significance scores are aggregated by summing \\-log\_{10}(p)\\
  values.

## Input data formats

The function auto-detects the input format by checking for
characteristic column names:

- **clusterProfiler** (`enrichGO`, `enrichKEGG`, `enrichPathway`, etc.)
  — recognized by the presence of `pvalue`, `p.adjust`, and `qvalue`
  columns.

- **enrichR** (web-based enrichment tool) — recognized by the presence
  of `P.value` and `Adjusted.P.value` columns. enrichR results are
  automatically converted to clusterProfiler-compatible format via
  plotthis's `prepare_enrichr_result()`.

If neither format is detected, the function stops with an error.

## Metric transformation

When a metric column is a p-value column (`pvalue`, `p.adjust`, or
`qvalue`), the function applies a \\-log\_{10}\\ transformation so that
more significant terms have higher values on the plot. The transformed
metric is stored internally as `.metric`. Other columns are used as-is.

`metric` accepts `"color_column:significance_column"`. The first column
drives the fill color and the `top_term` selection; the second drives
`cutoff`, which both drops the terms failing it and marks the surviving
ones. Both columns are transformed independently, so
`metric = "zScore:p.adjust"` colors cells by raw `zScore`, keeps only
the terms with `p.adjust < cutoff`, and marks them. A single column is
shorthand for `"column:column"`, so `metric = "p.adjust"` marks cells by
`p.adjust` without filtering.

For "bar", "dot", "lollipop", "comparison", and "heatmap" plot types,
`GeneRatio` (stored as strings like `"38/225"`) and `BgRatio` are
automatically converted to numeric ratios.

## Term ordering and selection

For each unique combination of `split_by`, `group_by`, and `facet_by`
levels, the function selects the `top_term` terms by the coloring metric
(the first column of `metric`): the smallest values when it is a p-value
column (`pvalue`, `p.adjust`, `qvalue`), and the largest values
otherwise (e.g. `zScore`). This ensures that each facet or split shows
its own most relevant terms rather than the globally most relevant ones.
When a significance metric is given (`"coloring:sig"`) together with
`cutoff`, terms that fail the cutoff in every group are dropped first
(except for dot, lollipop, and comparison plots, which show them in
gray); a term passing in any group is kept in all groups. The default
`top_term` is 6 for most plot types and 100 for "enrichmap" (which
benefits from showing more terms to reveal the network structure of gene
set relationships).

## See also

[`plotthis::BarPlot()`](https://pwwang.github.io/plotthis/reference/barplot.html),
[`plotthis::DotPlot()`](https://pwwang.github.io/plotthis/reference/dotplot.html),
[`plotthis::LollipopPlot()`](https://pwwang.github.io/plotthis/reference/dotplot.html),
[`plotthis::EnrichNetwork()`](https://pwwang.github.io/plotthis/reference/enrichmap1.html),
[`plotthis::EnrichMap()`](https://pwwang.github.io/plotthis/reference/enrichmap1.html),
[`plotthis::WordCloudPlot()`](https://pwwang.github.io/plotthis/reference/WordCloudPlot.html),
[`plotthis::Heatmap()`](https://pwwang.github.io/plotthis/reference/Heatmap.html),
[`plotthis::show_palettes()`](https://pwwang.github.io/plotthis/reference/show_palettes.html)

## Examples

``` r
# \donttest{
set.seed(8525)
data(enrich_example, package = "plotthis")
enrich_example$Group <- sample(LETTERS[1:3], nrow(enrich_example), replace = TRUE)
data(enrich_multidb_example, package = "plotthis")

EnrichmentPlot(enrich_example)

EnrichmentPlot(enrich_example, cutoff = 0.05)

EnrichmentPlot(enrich_example, palette = "Paired")


enrich_example$Description <- enrich_example$ID
EnrichmentPlot(enrich_example, plot_type = "heatmap", group_by = "Group",
 show_row_names = TRUE, show_column_names = TRUE, cutoff = 0.05)


# Color by zScore, and mark the cells whose p.adjust passes the cutoff
EnrichmentPlot(enrich_example, plot_type = "heatmap", group_by = "Group",
 metric = "zScore:p.adjust", cutoff = 0.05)

EnrichmentPlot(enrich_example, plot_type = "dot",
 metric = "zScore:p.adjust", cutoff = 1e-10)


# Multiple databases#'
EnrichmentPlot(enrich_multidb_example, facet_by = "Database", facet_nrow = 2)


enrich_example$Group <- sample(c("A", "B"), nrow(enrich_example), replace = TRUE)
EnrichmentPlot(enrich_example, plot_type = "comparison", group_by = "Group")

EnrichmentPlot(enrich_example, plot_type = "dot", top_term = 10)

EnrichmentPlot(enrich_example, plot_type = "lollipop", top_term = 10)

EnrichmentPlot(enrich_example, plot_type = "network")

EnrichmentPlot(enrich_example, plot_type = "enrichmap")

EnrichmentPlot(enrich_example, plot_type = "wordcloud")

# Wordcloud with feature
EnrichmentPlot(enrich_example, plot_type = "wordcloud", word_type = "feature")

# }
```
