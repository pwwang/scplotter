# Visualize differential expression markers

Visualize differential expression (DE) results — typically the output of
[`Seurat::FindMarkers()`](https://satijalab.org/seurat/reference/FindMarkers.html)
or
[`Seurat::FindAllMarkers()`](https://satijalab.org/seurat/reference/FindAllMarkers.html)
— across a variety of plot types. You can also compose the DE results
from other tools into a data frame with the required columns and use
this function to visualize them.

`MarkersPlot()` bridges the gap between DE testing and visualization by
providing a unified interface for both **summary-level DE
visualizations** (volcano, jitter, heatmap, and dot plots of fold
changes and significance) and **expression-level visualizations**
(violin, box, bar, ridge, heatmap, and dot plots of actual expression
values from a Seurat object).

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
  group_by = NULL,
  each = NULL,
  facet_each = FALSE,
  p_adjust = TRUE,
  cutoff = NULL,
  show_labels = FALSE,
  sig_mark = "*",
  order_by = "desc(abs(avg_log2FC))",
  select = ifelse(plot_type %in% c("volcano", "volcano_log2fc", "volcano_pct",
    "jitter", "jitter_log2fc", "jitter_pct"), 5, 10),
  flatten_markers = FALSE,
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

- group_by:

  Used only for expression-based plot types (ignored for DE summary plot
  types). A column in the Seurat object's metadata to group cells by,
  e.g., a condition column — useful when the DEs were calculated between
  conditions (such as cell cycle phases) and you want to compare the
  expression of the markers across those conditions. A single value is
  passed directly to
  [`FeatureStatPlot`](https://pwwang.github.io/scplotter/reference/FeatureStatPlot.md):
  for `heatmap` and `dot` plots it is applied as the column annotation
  (`ident`), and it only takes effect when `each` includes a metadata
  column mapping; for `violin`, `box`, `bar`, and `ridge` plots it is
  passed as `group_by`. The `"marker_column:metadata_column"` syntax
  (see **Metadata column mapping**) restricts the object to only the
  cells involved in the comparisons: for example, if a `comparison`
  column in the markers data frame holds `"G1:G2M"`, passing
  `group_by = "comparison:Phase"` keeps only G1 and G2M cells in the
  plot, with the `Phase` column re-factored to these two levels in the
  order they first appear in the `comparison` column. Without the
  restriction, e.g., `group_by = "Phase"`, all phase cells (G1, G2M,
  and S) are included in the plot. Default: `NULL`.

- each:

  A column name in `markers` indicating the grouping from which each
  marker was identified (e.g., the `cluster` column from
  `FindAllMarkers()`). Required for jitter and DE heatmap/dot plot
  types, where it defines the x-axis or column groups. For volcano plot
  types, it splits the plot by group (or facets it, with
  `facet_each = TRUE`). For expression plot types, `each` is used to
  select the markers within each group; a plain column name does not
  split the plot — use the `"marker_column:metadata_column"` syntax (see
  **Metadata column mapping**) to also split the plot by the mapped
  metadata column. Default: `NULL`.

- facet_each:

  Logical. Only for volcano plot types: if `TRUE`, facet the volcano
  plot by the `each` groups instead of splitting it into separate
  subplots. Ignored for other plot types. Default: `FALSE`.

- p_adjust:

  Logical. If `TRUE` (default), use adjusted p-value (`p_val_adj`
  column) for significance calculations and y-axis transformations. If
  `FALSE`, use raw p-value (`p_val` column).

- cutoff:

  Numeric. The p-value (or adjusted p-value, depending on `p_adjust`)
  threshold for labeling significance. For volcano plots, sets
  `y_cutoff`. For DE heatmap plots (`heatmap_log2fc`, `heatmap_pct`),
  controls which cells receive significance marks. For expression plot
  types with a numeric `select`, only markers with a p-value below
  `cutoff` are eligible for selection. Ignored by DE dot plots
  (`dot_log2fc`, `dot_pct`). Default: `NULL` (no cutoff; defaults to
  `0.05` for volcano plots).

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

  A string of one or more comma-separated expressions used to order the
  markers (evaluated with
  [`dplyr::arrange()`](https://dplyr.tidyverse.org/reference/arrange.html)).
  Can reference columns in `markers` as well as metadata columns merged
  in via a colon-form `each` (see **Metadata column mapping**). Only the
  first value of each merged metadata column is kept. Example:
  `"desc(avg_log2FC)"` or `"desc(avg_log2FC), desc(pct.1)"`. The
  ordering determines which markers are selected when `select` is
  numeric. For jitter plots, it is also passed to
  [`plotthis::JitterPlot()`](https://pwwang.github.io/plotthis/reference/jitterplot.html).
  Default: `"desc(abs(avg_log2FC))"`.

- select:

  How to select markers for display or labeling. See **Marker selection
  and filtering** section for full details.

  - Numeric: Top N markers per `each` group, or overall when `each` is
    `NULL` (default: `5` for volcano/jitter types, `10` for others).

  - Single expression: Filter condition for
    [`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html).

  - Character vector of multiple expressions (DE heatmap/dot plot types
    only): expressions mentioning the `each` column name filter the
    overall data, others filter within the remaining data.

- flatten_markers:

  Logical. Only for the expression `heatmap` and `dot` plot types. When
  `each` is used to select markers per group, the markers are by default
  provided to
  [`FeatureStatPlot`](https://pwwang.github.io/scplotter/reference/FeatureStatPlot.md)
  as a named list (one entry per group), which splits the feature rows
  of the plot by group. With `flatten_markers = TRUE`, the selected
  markers are collapsed into a single vector so the plot shows one
  unsplit block of features — useful e.g. to mimic
  [`Seurat::DoHeatmap()`](https://satijalab.org/seurat/reference/DoHeatmap.html)
  on globally selected markers. Default: `FALSE`.

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
      `layer`, `cell_type`. Note that `group_by`, `ident`, and
      `columns_split_by` are set by `MarkersPlot()` from the `group_by`
      and `each` arguments.

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

- `plot_type` determines which underlying plotting function is called
  and also what to be plotted. `volcano`, `volcano_log2fc`,
  `volcano_pct` `jitter`, `jitter_log2fc`, `jitter_pct`,
  `heatmap_log2fc`, `heatmap_pct`, `dot_log2fc`, and `dot_pct` are DE
  summary plots that visualize the DE statistics themselves, while
  `heatmap`, `violin`, `box`, `bar`, `ridge`, and `dot` are
  expression-based plots that visualize the actual expression values of
  the selected marker genes in the context of the original Seurat
  object.

- `each` is required for jitter plots (`"jitter"`, `"jitter_log2fc"`,
  `"jitter_pct"`) and DE heatmap/dot plots (`"heatmap_log2fc"`,
  `"heatmap_pct"`, `"dot_log2fc"`, `"dot_pct"`). Its role depends on the
  plot type:

  - Volcano plot types: the plot is split by the `each` groups (faceted
    when `facet_each = TRUE`).

  - Jitter plot types: the x-axis grouping.

  - DE heatmap/dot plot types: the columns of the heatmap/dot plot.

  - Expression plot types: used to select the markers within each group;
    it does not split or facet the plot. Pass
    `"marker_column:metadata_column"` (e.g.,
    `"cluster:seurat_clusters"`) to also split the plot by the mapped
    metadata column (via `columns_split_by` for heatmap/dot, or `ident`
    for violin/box/bar).

- When `each` uses the `"marker_column:metadata_column"` form, the
  markers data frame is left-joined with the object metadata. Only the
  first row per group is kept for non-key columns, which is sufficient
  for most annotation purposes but can cause issues if per-cell metadata
  is needed.

- The function calculates \\-log\_{10}(p)\\ (or
  \\-log\_{10}(p\_{adj})\\) internally and stores it in a temporary
  `neg_log10_p` column. This column is available for use in `order_by`.

## Metadata column mapping

Both `each` and `group_by` accept a `"marker_column:metadata_column"`
syntax that links a column in the markers data frame to a column in the
Seurat object's metadata.

- The part before the colon must be a column in `markers` (e.g.,
  `cluster`); the part after the colon must be a column in
  `object@meta.data` (e.g., `seurat_clusters`).

- This syntax requires `object` to be provided; otherwise an error is
  raised.

- Every value in the marker column must exist in the metadata column,
  otherwise an error is raised.

- For `each`, the metadata is merged into the markers data frame
  (keeping the first row of each metadata group for non-key columns), so
  metadata columns become available for arguments like `order_by`. On
  name conflicts, the merged columns get a `.meta` suffix.

- For `group_by`, the object is subset to the cells whose metadata
  values occur in the marker column, and the metadata column is
  re-factored with those values in the order they first appear in the
  marker column. Values separated by a colon (e.g., `"G1:G2M"`) are
  split into individual groups.

## Marker selection and filtering

How `select` picks the markers depends on the plot type and the value
provided:

- **Numeric** — Select the top `N` markers (ordered by `order_by`)
  within each group defined by `each`, or overall when `each` is `NULL`.
  Jitter plots label the top `N` markers per group (a numeric `select`
  is required). Volcano plots ignore `select` — labeling is controlled
  via `...` (e.g., `nlabel`). For expression plot types, a numeric
  `select` only keeps markers with a p-value below `cutoff` (when set)
  before the top-N selection.

- **Single expression** — A filter expression string evaluated by
  [`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html).
  For example, `"p_val_adj < 0.05 & avg_log2FC > 1"`. All markers
  matching the condition are retained across all groups.

- **Multiple expressions** (character vector) — Only for DE heatmap/dot
  plot types (`"heatmap_log2fc"`, `"heatmap_pct"`, `"dot_log2fc"`,
  `"dot_pct"`). Each element is evaluated independently: expressions
  that mention the `each` column filter the overall data (removing
  groups); other expressions filter within the remaining data. For
  example, `select = c("cluster %in% c('0', '1')", "p_val_adj < 0.05")`
  first restricts to clusters 0 and 1, then keeps only significant
  markers. A numeric string like `"5"` among the expressions is treated
  as a top-N selection.

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

MarkersPlot(allmarkers, plot_type = "jitter_pct", order_by = "desc(abs(pct.1 - pct.2))",
    each = "cluster", add_hline = 0, shape = 16)


MarkersPlot(allmarkers, plot_type = "heatmap_log2fc", each = "cluster",
    order_by = "desc(avg_log2FC)", select = 3)

MarkersPlot(allmarkers, plot_type = "heatmap_log2fc", each = "cluster",
    label = scales::label_number(accuracy = 0.01), select = 3,
    cutoff = 0.05, show_labels = TRUE, sig_mark = '{}')

MarkersPlot(allmarkers, plot_type = "heatmap_pct", each = "cluster",
    cutoff = 0.05, select = 3)


MarkersPlot(allmarkers, plot_type = "dot_log2fc", each = "cluster",
    add_reticle = TRUE, select = 3)


topmarkers <- allmarkers[order(allmarkers$avg_log2FC, decreasing = TRUE), ]
# Mimic Seurat's DoHeatmap()
MarkersPlot(topmarkers[1:20, ], object = pancreas_sub, plot_type = "heatmap",
   layer = "data", cell_type = "bars", flatten_markers = TRUE, cluster_rows = FALSE,
   show_column_names = "inplace", each = "cluster:seurat_clusters")
#> Warning: Layer counts isn't present in the assay object; returning NULL


# Select top 3 markers per cluster
MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "heatmap",
   order_by = "desc(avg_log2FC)", select = 3,
   layer = "data", cell_type = "bars",
   show_column_names = "inplace", each = "cluster:seurat_clusters")
#> Warning: Layer counts isn't present in the assay object; returning NULL

# Suppose we did a DE between G2M and G1 phases in each cluster and
# stored the results in a new column "comparison"
allmarkers$comparison <- "G1:G2M"
MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "heatmap",
   group_by = "comparison:Phase", each = "cluster:seurat_clusters",
   order_by = "desc(avg_log2FC)", select = 3, layer = "data")
#> Warning: Layer counts isn't present in the assay object; returning NULL


MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "dot", select = 2,
   flatten_markers = TRUE, order_by = "desc(avg_log2FC)",
   group_by = "Phase", each = "cluster:seurat_clusters", layer = "data")
#> Warning: Layer counts isn't present in the assay object; returning NULL


MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "violin", select = 2,
   position_dodge_preserve = "single", add_bg = TRUE, add_box = TRUE,
   group_by = "comparison:Phase", each = "cluster:seurat_clusters", layer = "data")
#> Warning: Layer counts isn't present in the assay object; returning NULL


# select markers with a custom condition, e.g.,
# significant markers in cluster 0, 1, and 2 with pct.2 - pct.1 > 0.6
# Note that other clusters are still included in the plot
MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "violin",
  select = c('cluster %in% c("1", "2", "0") & pct.2 - pct.1 > 0.6'),
  each = "cluster:seurat_clusters", cutoff = 0.05, layer = "data")
#> Warning: Layer counts isn't present in the assay object; returning NULL


MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "box", select = 3,
  group_by = "Phase", each = "cluster:seurat_clusters", layer = "data")
#> Warning: Layer counts isn't present in the assay object; returning NULL


MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "ridge", select = 2,
   group_by = "Phase", each = "cluster:seurat_clusters", layer = "data",
   ncol = 4)
#> Warning: Layer counts isn't present in the assay object; returning NULL
#> Picking joint bandwidth of 0.441
#> Picking joint bandwidth of 0.226
#> Picking joint bandwidth of 0.157
#> Picking joint bandwidth of 0.144
#> Picking joint bandwidth of 0.232
#> Picking joint bandwidth of 0.441
#> Picking joint bandwidth of 0.293
#> Picking joint bandwidth of 0.334
#> Picking joint bandwidth of 0.229
#> Picking joint bandwidth of 0.345
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
#> Picking joint bandwidth of 0.283
#> Picking joint bandwidth of 0.337
#> Picking joint bandwidth of 0.242
#> Picking joint bandwidth of 0.129
#> Picking joint bandwidth of 0.27
#> Picking joint bandwidth of 0.498
#> Picking joint bandwidth of 0.358
#> Picking joint bandwidth of 0.498
#> Picking joint bandwidth of 0.223
#> Picking joint bandwidth of 0.498
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
#> Picking joint bandwidth of 0.0624
#> Picking joint bandwidth of 0.255
#> Picking joint bandwidth of 0.102
#> Picking joint bandwidth of 0.0994
#> Picking joint bandwidth of 0.0465
#> Picking joint bandwidth of 0.393
#> Picking joint bandwidth of 0.126
#> Picking joint bandwidth of 0.216
#> Picking joint bandwidth of 0.203
#> Picking joint bandwidth of 0.294
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
#> Picking joint bandwidth of 0.525
#> Picking joint bandwidth of 0.16
#> Picking joint bandwidth of 0.195
#> Picking joint bandwidth of 0.465
#> Picking joint bandwidth of 0.413
#> Picking joint bandwidth of 0.457
#> Picking joint bandwidth of 0.481
#> Picking joint bandwidth of 0.574
#> Picking joint bandwidth of 0.413
#> Picking joint bandwidth of 0.206
#> Picking joint bandwidth of 0.441
#> Picking joint bandwidth of 0.226
#> Picking joint bandwidth of 0.157
#> Picking joint bandwidth of 0.144
#> Picking joint bandwidth of 0.232
#> Picking joint bandwidth of 0.441
#> Picking joint bandwidth of 0.293
#> Picking joint bandwidth of 0.334
#> Picking joint bandwidth of 0.229
#> Picking joint bandwidth of 0.345
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
#> Picking joint bandwidth of 0.283
#> Picking joint bandwidth of 0.337
#> Picking joint bandwidth of 0.242
#> Picking joint bandwidth of 0.129
#> Picking joint bandwidth of 0.27
#> Picking joint bandwidth of 0.498
#> Picking joint bandwidth of 0.358
#> Picking joint bandwidth of 0.498
#> Picking joint bandwidth of 0.223
#> Picking joint bandwidth of 0.498
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
#> Picking joint bandwidth of 0.0624
#> Picking joint bandwidth of 0.255
#> Picking joint bandwidth of 0.102
#> Picking joint bandwidth of 0.0994
#> Picking joint bandwidth of 0.0465
#> Picking joint bandwidth of 0.393
#> Picking joint bandwidth of 0.126
#> Picking joint bandwidth of 0.216
#> Picking joint bandwidth of 0.203
#> Picking joint bandwidth of 0.294
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
#> Picking joint bandwidth of 0.525
#> Picking joint bandwidth of 0.16
#> Picking joint bandwidth of 0.195
#> Picking joint bandwidth of 0.465
#> Picking joint bandwidth of 0.413
#> Picking joint bandwidth of 0.457
#> Picking joint bandwidth of 0.481
#> Picking joint bandwidth of 0.574
#> Picking joint bandwidth of 0.413
#> Picking joint bandwidth of 0.206

# }
```
