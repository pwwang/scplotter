# Visualize differential expression markers

Visualize differential expression (DE) results — typically the output of
[`Seurat::FindMarkers()`](https://satijalab.org/seurat/reference/FindMarkers.html)
or
[`Seurat::FindAllMarkers()`](https://satijalab.org/seurat/reference/FindAllMarkers.html)
— across a variety of plot types. `MarkersPlot()` bridges the gap
between DE testing and visualization by providing a unified interface
for both **summary-level DE visualizations** (volcano, jitter, heatmap,
and dot plots of fold changes and significance) and **expression-level
visualizations** (violin, box, bar, ridge, heatmap, and dot plots of
actual expression values from a Seurat object).

The function handles two broad categories of plots:

- **DE summary plots** (no `object` required): visualize the DE
  statistics themselves — log2 fold change, percentage difference,
  p-values, and adjusted p-values — across groups or comparisons.

  - `"volcano"` / `"volcano_log2fc"` — Volcano plot with log2 fold
    change on the x-axis and \\-log\_{10}(p)\\ on the y-axis. Genes
    passing the `cutoff` are highlighted and top genes are labeled.
    Ideal for overview of effect size vs. significance.

  - `"volcano_pct"` — Volcano plot with percentage-point difference
    (`pct.1 - pct.2`) on the x-axis. Useful when the biological question
    is about detection rate rather than expression magnitude.

  - `"jitter"` / `"jitter_log2fc"` — Jitter plot of log2 fold changes
    across groups (defined by `each`). Dot size encodes
    \\-log\_{10}(p)\\. Reveals distribution of effect sizes per cluster
    or condition.

  - `"jitter_pct"` — Jitter plot of percentage-point differences across
    groups.

  - `"heatmap_log2fc"` — Heatmap of log2 fold changes (genes × groups).
    Cells can be marked for significance via `cutoff` and `sig_mark`.

  - `"heatmap_pct"` — Heatmap of percentage-point differences (genes ×
    groups). Same significance-marking support.

  - `"dot_log2fc"` — Dot plot of log2 fold changes (genes × groups). Dot
    size encodes \\-log\_{10}(p)\\.

  - `"dot_pct"` — Dot plot of percentage-point differences (genes ×
    groups). Dot size encodes \\-log\_{10}(p)\\.

- **Expression plots** (`object` required): visualize the actual
  expression values of the selected marker genes in the context of the
  original Seurat object. These are useful for validating DE results by
  inspecting the underlying expression distributions.

  - `"heatmap"` — Expression heatmap of selected marker genes.

  - `"violin"` — Violin plots of expression per gene.

  - `"box"` — Box plots of expression per gene.

  - `"bar"` — Bar plots of mean expression per gene.

  - `"ridge"` — Ridge plots of expression distribution per gene.

  - `"dot"` — Dot plot of expression (fraction expressing × mean
    expression) per gene.

## Usage

``` r
MarkersPlot(
  markers,
  object = NULL,
  plot_type = c("volcano", "volcano_log2fc", "volcano_pct", "jitter", "jitter_log2fc",
    "jitter_pct", "heatmap_log2fc", "heatmap_pct", "dot_log2fc", "dot_pct", "heatmap",
    "violin", "box", "bar", "ridge", "dot"),
  subset_by = NULL,
  each = NULL,
  subset_as_facet = FALSE,
  facet_each = FALSE,
  comparison_by = NULL,
  p_adjust = TRUE,
  cutoff = NULL,
  show_labels = FALSE,
  sig_mark = "*",
  order_by = "desc(abs(avg_log2FC))",
  select = ifelse(plot_type %in% c("volcano", "volcano_log2fc", "volcano_pct",
    "jitter", "jitter_log2fc", "jitter_pct"), 5, 10),
  ...
)
```

## Arguments

- markers:

  A data frame of differential expression results, typically the output
  of
  [`Seurat::FindMarkers()`](https://satijalab.org/seurat/reference/FindMarkers.html)
  or
  [`Seurat::FindAllMarkers()`](https://satijalab.org/seurat/reference/FindAllMarkers.html).
  Must contain columns `"gene"` (or gene symbols as rownames),
  `"p_val"`, and `"avg_log2FC"`. For percentage-based plots
  (`volcano_pct`, `jitter_pct`, `heatmap_pct`, `dot_pct`), columns
  `"pct.1"` and `"pct.2"` are also required.

- object:

  A Seurat object. Required for expression-based plot types:
  `"heatmap"`, `"violin"`, `"box"`, `"bar"`, `"ridge"`, and `"dot"`. Not
  used for DE summary plot types. Default: `NULL`.

- plot_type:

  The type of plot to generate. One of `"volcano"`, `"volcano_log2fc"`,
  `"volcano_pct"`, `"jitter"`, `"jitter_log2fc"`, `"jitter_pct"`,
  `"heatmap_log2fc"`, `"heatmap_pct"`, `"dot_log2fc"`, `"dot_pct"`,
  `"heatmap"`, `"violin"`, `"box"`, `"bar"`, `"ridge"`, or `"dot"`. See
  Description for details on each type.

- subset_by:

  Deprecated. Use `each` instead.

- each:

  A column name in `markers` indicating the grouping from which each
  marker was identified (e.g., the `cluster` column from
  `FindAllMarkers()`). Supports the `"marker_column:metadata_column"`
  syntax for linking to Seurat object metadata (see **Metadata column
  mapping** section). For jitter and DE heatmap/dot plot types, `each`
  is required and defines the x-axis or column groups. For expression
  plot types, `each` controls faceting or splitting. Default: `NULL`.

- subset_as_facet:

  Deprecated. Use `facet_each` instead.

- facet_each:

  Logical. If `TRUE`, facet the plot by `each` groups instead of
  splitting into separate plots. Most useful for expression plot types.
  For volcano plots, controls whether faceting or split_by dispatch is
  used. Default: `FALSE`.

- comparison_by:

  A column name in `markers` indicating the comparison (e.g., `"g1:g2"`
  for a pairwise comparison, or a single group name for one-vs-rest).
  Required for expression-based plot types (`"heatmap"`, `"violin"`,
  `"box"`, `"bar"`, `"ridge"`, `"dot"`). Supports the
  `"marker_column:metadata_column"` syntax (see **Metadata column
  mapping** section). If the comparison values contain a colon (e.g.,
  `"G2M:G1"`), the two groups on either side of the colon are used to
  subset the object. If only a single group is present, a one-vs-other
  comparison is assumed. Default: `NULL`.

- p_adjust:

  Logical. If `TRUE` (default), use adjusted p-value (`p_val_adj`
  column) for significance calculations and y-axis transformations. If
  `FALSE`, use raw p-value (`p_val` column).

- cutoff:

  Numeric. The p-value (or adjusted p-value, depending on `p_adjust`)
  threshold for labeling significance. For volcano plots, sets
  `y_cutoff`. For heatmap-based DE plots (`heatmap_log2fc`,
  `heatmap_pct`), controls which cells receive significance marks.
  Default: `NULL` (no cutoff; defaults to `0.05` for volcano plots).

- show_labels:

  Logical. For `heatmap_log2fc` and `heatmap_pct` plot types only. If
  `TRUE`, display numeric values in heatmap cells. When combined with
  `cutoff`, both values and significance marks are shown. Default:
  `FALSE`.

- sig_mark:

  Character. The symbol or compound mark used to annotate statistically
  significant cells in `heatmap_log2fc` and `heatmap_pct` plots. Must be
  a valid ComplexHeatmap mark: single characters (`"-"`, `"|"`, `"+"`,
  `"/"`, `"\\"`, `"x"`, `"o"`) or compound marks (`"[*]"`, `"<*>"`,
  `"(*)"`, `"{*}"`). Note that `"*"` conflicts with `show_labels = TRUE`
  because both use the label layer — use a compound mark instead.
  Default: `"*"`.

- order_by:

  A string expression to order markers by (evaluated with
  [`dplyr::arrange()`](https://dplyr.tidyverse.org/reference/arrange.html)).
  Can reference columns in `markers` as well as columns from the object
  metadata (when `object` is provided and `each` enables merging). Only
  the first value of merged metadata columns is used. Example:
  `"desc(avg_log2FC)"`. The ordering affects which markers are selected
  when `select` is numeric. Default: `desc(abs(avg_log2FC))`.

- select:

  How to select markers for labeling or display. See **Marker selection
  and filtering** section for full details.

  - Numeric: Top N markers per `each` group (default: `5` for
    volcano/jitter types, `10` for others).

  - Character expression: Filter condition for
    [`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html).

  - Character vector: Multiple filter expressions; those containing the
    `each` column name filter the overall data, others filter within
    remaining data.

- ...:

  Additional arguments passed to the underlying plotting function,
  depending on `plot_type`:

  For `volcano`, `volcano_log2fc`, `volcano_pct`

  :   Passed to
      [`plotthis::VolcanoPlot()`](https://pwwang.github.io/plotthis/reference/VolcanoPlot.html).
      Common arguments: `x_cutoff`, `x_cutoff_name`, `label_by`,
      `color_by`, `nlabel`, `flip_negative`.

  For `jitter`, `jitter_log2fc`, `jitter_pct`

  :   Passed to
      [`plotthis::JitterPlot()`](https://pwwang.github.io/plotthis/reference/jitterplot.html).
      Common arguments: `add_hline`, `shape`, `size_by`, `nlabel`.

  For `heatmap_log2fc`, `heatmap_pct`, `dot_log2fc`, `dot_pct`

  :   Passed to
      [`plotthis::Heatmap()`](https://pwwang.github.io/plotthis/reference/Heatmap.html).
      Common arguments: `show_row_names`, `show_column_names`,
      `values_fill`, `palette`, `cluster_rows`, `cluster_columns`,
      `add_reticle`.

  For `heatmap`, `violin`, `box`, `bar`, `ridge`, `dot`

  :   Passed to
      [`FeatureStatPlot`](https://pwwang.github.io/scplotter/reference/FeatureStatPlot.md).
      Common arguments: `name`, `palette`, `ncol`, `nrow`, `stack`,
      `columns_split_by`.

## Value

A ggplot object (from
[`plotthis::VolcanoPlot()`](https://pwwang.github.io/plotthis/reference/VolcanoPlot.html)
or
[`plotthis::JitterPlot()`](https://pwwang.github.io/plotthis/reference/jitterplot.html)),
a Heatmap object (from
[`plotthis::Heatmap()`](https://pwwang.github.io/plotthis/reference/Heatmap.html)),
or a ggplot/patchwork object (from
[`FeatureStatPlot`](https://pwwang.github.io/scplotter/reference/FeatureStatPlot.md)).
When `split_by` or faceting generates multiple plots and
`combine = TRUE` (default), a combined patchwork object is returned;
when `combine = FALSE`, a list of individual plots is returned.

## Note

- `each` is required for jitter plots (`"jitter"`, `"jitter_log2fc"`,
  `"jitter_pct"`) and DE heatmap/dot plots (`"heatmap_log2fc"`,
  `"heatmap_pct"`, `"dot_log2fc"`, `"dot_pct"`). Without it, there is no
  grouping axis.

- `comparison_by` is required for expression-based plot types
  (`"heatmap"`, `"violin"`, `"box"`, `"bar"`, `"ridge"`, `"dot"`) — it
  tells the function which comparison groups to extract from the object.

- When `object` is provided and `each` maps to a metadata column, the
  markers data frame is left-joined with the object metadata. Only the
  first row per group is kept for non-key columns, which is sufficient
  for most annotation purposes but can cause issues if per-cell metadata
  is needed.

- For expression-based heatmap and dot plots, when `each_2` is available
  (i.e., the metadata column is mapped), genes are automatically grouped
  by `each` via `columns_split_by`, and `group_by` is set to `NULL`.

- The function calculates \\-log\_{10}(p)\\ (or
  \\-log\_{10}(p\_{adj})\\) internally and stores it in a temporary
  `neg_log10_p` column. This column is available for use in `order_by`.

- When the comparison involves only a single group (one-vs-rest), cells
  not in the comparison group are labeled `"Other"` in the object
  metadata.

## Metadata column mapping

Both `each` and `comparison_by` support a special
`"marker_column:metadata_column"` syntax for linking columns in the
markers data frame to columns in the Seurat object's metadata.

- The part before the colon refers to a column in `markers`.

- The part after the colon refers to a column in `object@meta.data`.

- If only one name is provided (no colon), it is used for both the
  markers column and the metadata column (if a matching metadata column
  exists).

- Example: `each = "cluster:RNA_snn_res.0.8"` maps the `cluster` column
  in the DE results to the `RNA_snn_res.0.8` column in the Seurat
  metadata. For expression-based plots, this allows the function to
  extract the relevant expression values from the Seurat object for the
  specified groups. You can also specify `each = "cluster:NULL"` to use
  "cluster" to select markers in each cluster but don't split the
  expression plots by cluster.

When the markers data frame and object metadata are merged via `each`,
only the first value of each non-key column within each group is
retained — this is by design to avoid duplication.

## Marker selection and filtering

The `select` argument supports three modes:

- **Numeric** — Select the top `N` markers (ordered by `order_by`)
  within each group defined by `each`. For volcano and jitter plots, all
  markers are plotted but only the top `N` per group are labeled. For
  other plot types, only the selected markers are shown.

- **Single expression** — A filter expression string evaluated by
  [`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html).
  For example, `"p_val_adj < 0.05 & avg_log2FC > 1"`. All markers
  matching the condition are retained across all groups.

- **Multiple expressions** (character vector) — Each element is
  evaluated independently. Expressions that mention the `each` column
  filter the overall data (removing groups); other expressions filter
  within the remaining data. For example,
  `select = c("cluster %in% c('0', '1')", "p_val_adj < 0.05")` first
  restricts to clusters 0 and 1, then keeps only significant markers. A
  numeric string like `"5"` among the expressions is treated as a top-N
  selection.

Default `select`: `5` for volcano and jitter plot types, `10` for all
other plot types.

## Significance marking in heatmaps

For `heatmap_log2fc` and `heatmap_pct`, the `cutoff` and `sig_mark`
arguments control how statistically significant cells are annotated in
the heatmap:

- When `cutoff` is set and `show_labels = FALSE`, cells with p-value (or
  adjusted p-value) below the cutoff are marked with `sig_mark` using
  ComplexHeatmap's mark system. Valid `sig_mark` values include `"-"`,
  `"|"`, `"+"`, `"/"`, `"\\"`, `"x"`, `"o"`, and compound marks like
  `"[*]"`, `"<*>"`, `"(*)"`, `"{*}"`.

- When `cutoff` is set and `show_labels = TRUE`, both numeric values and
  significance marks are displayed (`cell_type = "label+mark"`). Note
  that `sig_mark = "*"` does not work with `show_labels = TRUE` — use
  compound marks instead.

- When `cutoff = NULL` and `show_labels = TRUE`, all cells are labeled
  with their numeric values.

## See also

[`plotthis::VolcanoPlot()`](https://pwwang.github.io/plotthis/reference/VolcanoPlot.html),
[`plotthis::JitterPlot()`](https://pwwang.github.io/plotthis/reference/jitterplot.html),
[`plotthis::Heatmap()`](https://pwwang.github.io/plotthis/reference/Heatmap.html),
[`FeatureStatPlot`](https://pwwang.github.io/scplotter/reference/FeatureStatPlot.md),
[`Seurat::FindMarkers()`](https://satijalab.org/seurat/reference/FindMarkers.html),
[`Seurat::FindAllMarkers()`](https://satijalab.org/seurat/reference/FindAllMarkers.html)

## Examples

``` r
# \donttest{
data(pancreas_sub)
markers <- Seurat::FindMarkers(pancreas_sub,
 group.by = "Phase", ident.1 = "G2M", ident.2 = "G1")
#> For a (much!) faster implementation of the Wilcoxon Rank Sum Test,
#> (default method for FindMarkers) please install the presto package
#> --------------------------------------------
#> install.packages('devtools')
#> devtools::install_github('immunogenomics/presto')
#> --------------------------------------------
#> After installation of presto, Seurat will automatically use the more 
#> efficient implementation (no further action necessary).
#> This message will be shown once per session
allmarkers <- Seurat::FindAllMarkers(pancreas_sub)  # seurat_clusters
#> Calculating cluster 0
#> Calculating cluster 1
#> Calculating cluster 2
#> Calculating cluster 3
#> Calculating cluster 4
#> Calculating cluster 5
#> Calculating cluster 6

MarkersPlot(markers)
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf

MarkersPlot(markers, x_cutoff = 2)
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf

MarkersPlot(allmarkers, each = "cluster", ncol = 2, facet_each = TRUE)

MarkersPlot(markers, plot_type = "volcano_pct", flip_negative = TRUE)
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf


MarkersPlot(allmarkers, plot_type = "jitter", each = "cluster")

MarkersPlot(allmarkers, plot_type = "jitter_pct",
    each = "cluster", add_hline = 0, shape = 16)


MarkersPlot(allmarkers, plot_type = "heatmap_log2fc", each = "cluster",
    order_by = "desc(avg_log2FC)")

MarkersPlot(allmarkers, plot_type = "heatmap_log2fc", each = "cluster",
    label = scales::label_number(accuracy = 0.01), select = 3,
    cutoff = 0.05, show_labels = TRUE, sig_mark = '{}')

MarkersPlot(allmarkers, plot_type = "heatmap_pct", each = "cluster",
    cutoff = 0.05, select = 3)


MarkersPlot(allmarkers, plot_type = "dot_log2fc", each = "cluster",
    add_reticle = TRUE, select = 3)


MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "heatmap",
   columns_split_by = "CellType", layer = "data",
   comparison_by = "cluster:seurat_clusters")
#> Warning: Layer counts isn't present in the assay object; returning NULL


# Suppose we did a DE between G2M and G1 phases in each cluster and
# stored the results in a new column "comparison"
allmarkers$comparison <- "G1:G2M"
MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "heatmap",
   comparison_by = "comparison:Phase", each = "cluster:seurat_clusters",
   layer = "data")
#> Warning: Layer counts isn't present in the assay object; returning NULL

# Select by each cluster, but don't split by seurat clusters in the plot
MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "heatmap",
   comparison_by = "comparison:Phase", each = "cluster:NULL", select = 3,
   layer = "data")
#> Warning: Layer counts isn't present in the assay object; returning NULL

MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "dot", select = 2,
   comparison_by = "Phase", each = "cluster:seurat_clusters", layer = "data")
#> Warning: Layer counts isn't present in the assay object; returning NULL


MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "violin", select = 2,
   comparison_by = "Phase", each = "cluster:seurat_clusters", layer = "data")
#> Warning: Layer counts isn't present in the assay object; returning NULL


# select markers with a custom condition, e.g.,
# significant markers in cluster 0, 1, and 2 with pct.2 - pct.1 > 0.6
# Note that other clusters are still included in the plot
MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "violin", each = "cluster",
   select = c('cluster %in% c("1", "2", "0") & pct.2 - pct.1 > 0.6'),
   comparison_by = "cluster:seurat_clusters", cutoff = 0.05, layer = "data")
#> Warning: [MarkersPlot] `each` 'cluster' is ignored, since it is not found in the object's metadata. Set `each` to 'cluster:<object_metadata_column>' to make it work.
#> Warning: Layer counts isn't present in the assay object; returning NULL


# To exclude other clusters, you can separate the filtering conditions into
# multiple expressions
MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "violin", each = "cluster",
   select = c('cluster %in% c("1", "2", "0")', 'pct.2 - pct.1 > 0.6'),
   comparison_by = "cluster:seurat_clusters", cutoff = 0.05, layer = "data")
#> Warning: [MarkersPlot] `each` 'cluster' is ignored, since it is not found in the object's metadata. Set `each` to 'cluster:<object_metadata_column>' to make it work.
#> Warning: Layer counts isn't present in the assay object; returning NULL


MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "box", select = 3,
   comparison_by = "Phase", each = "cluster:seurat_clusters", layer = "data")
#> Warning: Layer counts isn't present in the assay object; returning NULL


MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "ridge", select = 2,
   comparison_by = "Phase", each = "cluster:seurat_clusters", layer = "data",
   ncol = 4)
#> Warning: Layer counts isn't present in the assay object; returning NULL
#> Picking joint bandwidth of 0.464
#> Picking joint bandwidth of 0.303
#> Picking joint bandwidth of 0.177
#> Picking joint bandwidth of 0.142
#> Picking joint bandwidth of 0.326
#> Picking joint bandwidth of 0.464
#> Picking joint bandwidth of 0.359
#> Picking joint bandwidth of 0.303
#> Picking joint bandwidth of 0.222
#> Picking joint bandwidth of 0.32
#> Picking joint bandwidth of 0.33
#> Picking joint bandwidth of 0.433
#> Picking joint bandwidth of 0.881
#> Picking joint bandwidth of 0.375
#> Picking joint bandwidth of 0.202
#> Picking joint bandwidth of 0.319
#> Picking joint bandwidth of 0.586
#> Picking joint bandwidth of 0.0593
#> Picking joint bandwidth of 0.0355
#> Picking joint bandwidth of 0.102
#> Picking joint bandwidth of 0.384
#> Picking joint bandwidth of 0.427
#> Picking joint bandwidth of 0.323
#> Picking joint bandwidth of 0.133
#> Picking joint bandwidth of 0.368
#> Picking joint bandwidth of 0.512
#> Picking joint bandwidth of 0.423
#> Picking joint bandwidth of 0.512
#> Picking joint bandwidth of 0.234
#> Picking joint bandwidth of 0.515
#> Picking joint bandwidth of 0.363
#> Picking joint bandwidth of 0.426
#> Picking joint bandwidth of 0.125
#> Picking joint bandwidth of 0.131
#> Picking joint bandwidth of 0.187
#> Picking joint bandwidth of 0.327
#> Picking joint bandwidth of 0.497
#> Picking joint bandwidth of 0.0635
#> Picking joint bandwidth of 0.327
#> Picking joint bandwidth of 0.327
#> Picking joint bandwidth of 0.0789
#> Picking joint bandwidth of 0.116
#> Picking joint bandwidth of 0.102
#> Picking joint bandwidth of 0.0982
#> Picking joint bandwidth of 0.0543
#> Picking joint bandwidth of 0.392
#> Picking joint bandwidth of 0.115
#> Picking joint bandwidth of 0.0379
#> Picking joint bandwidth of 0.229
#> Picking joint bandwidth of 0.193
#> Picking joint bandwidth of 0.167
#> Picking joint bandwidth of 0.325
#> Picking joint bandwidth of 0.241
#> Picking joint bandwidth of 0.173
#> Picking joint bandwidth of 0.468
#> Picking joint bandwidth of 0.358
#> Picking joint bandwidth of 0.203
#> Picking joint bandwidth of 0.499
#> Picking joint bandwidth of 0.0658
#> Picking joint bandwidth of 0.0904
#> Picking joint bandwidth of 0.267
#> Picking joint bandwidth of 0.28
#> Picking joint bandwidth of 0.12
#> Picking joint bandwidth of 0.146
#> Picking joint bandwidth of 0.0428
#> Picking joint bandwidth of 0.131
#> Picking joint bandwidth of 0.179
#> Picking joint bandwidth of 0.364
#> Picking joint bandwidth of 0.0428
#> Picking joint bandwidth of 0.25
#> Picking joint bandwidth of 0.464
#> Picking joint bandwidth of 0.303
#> Picking joint bandwidth of 0.177
#> Picking joint bandwidth of 0.142
#> Picking joint bandwidth of 0.326
#> Picking joint bandwidth of 0.464
#> Picking joint bandwidth of 0.359
#> Picking joint bandwidth of 0.303
#> Picking joint bandwidth of 0.222
#> Picking joint bandwidth of 0.32
#> Picking joint bandwidth of 0.33
#> Picking joint bandwidth of 0.433
#> Picking joint bandwidth of 0.881
#> Picking joint bandwidth of 0.375
#> Picking joint bandwidth of 0.202
#> Picking joint bandwidth of 0.319
#> Picking joint bandwidth of 0.586
#> Picking joint bandwidth of 0.0593
#> Picking joint bandwidth of 0.0355
#> Picking joint bandwidth of 0.102
#> Picking joint bandwidth of 0.384
#> Picking joint bandwidth of 0.427
#> Picking joint bandwidth of 0.323
#> Picking joint bandwidth of 0.133
#> Picking joint bandwidth of 0.368
#> Picking joint bandwidth of 0.512
#> Picking joint bandwidth of 0.423
#> Picking joint bandwidth of 0.512
#> Picking joint bandwidth of 0.234
#> Picking joint bandwidth of 0.515
#> Picking joint bandwidth of 0.363
#> Picking joint bandwidth of 0.426
#> Picking joint bandwidth of 0.125
#> Picking joint bandwidth of 0.131
#> Picking joint bandwidth of 0.187
#> Picking joint bandwidth of 0.327
#> Picking joint bandwidth of 0.497
#> Picking joint bandwidth of 0.0635
#> Picking joint bandwidth of 0.327
#> Picking joint bandwidth of 0.327
#> Picking joint bandwidth of 0.0789
#> Picking joint bandwidth of 0.116
#> Picking joint bandwidth of 0.102
#> Picking joint bandwidth of 0.0982
#> Picking joint bandwidth of 0.0543
#> Picking joint bandwidth of 0.392
#> Picking joint bandwidth of 0.115
#> Picking joint bandwidth of 0.0379
#> Picking joint bandwidth of 0.229
#> Picking joint bandwidth of 0.193
#> Picking joint bandwidth of 0.167
#> Picking joint bandwidth of 0.325
#> Picking joint bandwidth of 0.241
#> Picking joint bandwidth of 0.173
#> Picking joint bandwidth of 0.468
#> Picking joint bandwidth of 0.358
#> Picking joint bandwidth of 0.203
#> Picking joint bandwidth of 0.499
#> Picking joint bandwidth of 0.0658
#> Picking joint bandwidth of 0.0904
#> Picking joint bandwidth of 0.267
#> Picking joint bandwidth of 0.28
#> Picking joint bandwidth of 0.12
#> Picking joint bandwidth of 0.146
#> Picking joint bandwidth of 0.0428
#> Picking joint bandwidth of 0.131
#> Picking joint bandwidth of 0.179
#> Picking joint bandwidth of 0.364
#> Picking joint bandwidth of 0.0428
#> Picking joint bandwidth of 0.25

# }
```
