# Cell statistics plot

This function creates a plot to visualize the statistics of cells in a
Seurat object, a Giotto object, a path to an .h5ad file or an opened
`H5File` by `hdf5r` package. It can create various types of plots,
including bar plots, circos plots, pie charts, pies (heatmap with
cell_type = 'pie'), ring/donut plots, trend plots area plots,
sankey/alluvial plots, heatmaps, radar plots, spider plots, violin
plots, and box plots. The function allows for grouping, splitting, and
faceting the data based on metadata columns. It also supports
calculating fractions of cells based on specified groupings.#'

## Usage

``` r
CellStatPlot(
  object,
  ident = NULL,
  group_by = NULL,
  group_by_sep = "_",
  spat_unit = NULL,
  feat_type = NULL,
  split_by = NULL,
  split_by_sep = "_",
  facet_by = NULL,
  rows_by = NULL,
  columns_split_by = NULL,
  frac = c("none", "group", "ident", "cluster", "all"),
  rows_name = NULL,
  name = NULL,
  plot_type = c("bar", "circos", "pie", "pies", "ring", "donut", "trend", "area",
    "sankey", "alluvial", "heatmap", "radar", "spider", "violin", "box"),
  swap = FALSE,
  ylab = NULL,
  ...
)
```

## Arguments

- object:

  A Seurat object, a Giotto object, a path to h5ad file or an opened
  `H5File` (from `hdf5r` package) object a data frame (for internal use)
  containing cell metadata.

- ident:

  The column with the cell identities. i.e. clusters. Default: NULL If
  NULL, the active identity of the Seurat object and the name "Identity"
  will be used. For 'pies', this will be used as the `pie_group_by`. For
  'heatmap' plot, this will be used as the rows_by of the heatmap.

- group_by:

  The column name in the meta data to group the cells. Default: NULL
  This should work as the columns of the plot_type: heatmap. For
  violin/box plot, at most 2 `group_by` columns are allowed and they
  will not be concatenated. The first one is used to break down the
  values in groups, and the second one works as the `group_by` argument
  in plotthis::ViolinPlot/plotthis::BoxPlot.

- group_by_sep:

  The separator to use when combining multiple columns in `group_by`.
  Default: "\_" For 'sankey'/'heatmap' plot, multiple columns will not
  be combined, and each of them will be used as a node.

- spat_unit:

  The spatial unit to use for the plot. Only applied to Giotto objects.

- feat_type:

  feature type of the features (e.g. "rna", "dna", "protein"), only
  applied to Giotto objects.

- split_by:

  The column name in the meta data to split the cells. Default: NULL
  Each split will be plotted in a separate plot.

- split_by_sep:

  The separator to use when combining multiple columns in `split_by`.
  Default: "\_"

- facet_by:

  The column name in the meta data to facet the plots. Default: NULL Not
  available for 'circos', 'sankey', and 'heatmap' plots.

- rows_by:

  The column names in the data used as the rows of the 'pies' (heatmap
  with cell_type = 'pie'). Default: NULL. Only available for 'pies'
  plot. The values don't matter, and they only indicate the cells
  overlapping with the columns and distributed in different `ident`
  values.

- columns_split_by:

  The column name in the meta data to split the columns of the
  'pies'/'heatmap' plot. Default: NULL

- frac:

  The way of calculating the fraction. Default is "none". Possible
  values are "group", "ident", "cluster", "all", "none". Note that the
  fractions are calculated in each split and facet group if `split_by`
  and `facet_by` are specified.

  - group: calculate the fraction in each group. The total fraction of
    the cells of idents in each group will be 1. When `group-by` is not
    specified, it will be the same as `all`.

  - ident: calculate the fraction in each ident. The total fraction of
    the cells of groups in each ident will be 1. Only works when
    `group-by` is specified.

  - cluster: alias of `ident`.

  - all: calculate the fraction against all cells.

  - none: do not calculate the fraction, use the number of cells
    instead.

- rows_name:

  The name of the rows in the 'pies'/'heatmap' plot. Default is NULL.

- name:

  The name of the 'pies'/'heatmap' plot, shown as the name of the main
  legend. Default is NULL.

- plot_type:

  The type of plot to use. Default is "bar". Possible values are "bar",
  "circos", "pie", "pies", "ring"/"donut", "trend", "area", "heatmap",
  "sankey"/"alluvial", "radar" and "spider". 'pie' vs 'pies': 'pie' plot
  will plot a single pie chart for each group, while 'pies' plot will
  plot multiple pie charts for each group and split. 'pies' basically is
  a heatmap with 'cell_type = "pie"'.

- swap:

  Whether to swap the cluster and group, that is, using group as the
  x-axis and cluster to fill the plot. For circos plot, when transposed,
  the arrows will be drawn from the idents (by `ident`) to the the
  groups (by `group_by`). Only works when `group_by` is specified. For
  'pies' plot, this will swap the group_by and ident. For 'heatmap'
  plot, this will swap the group_by and columns_split_by. Note that this
  is different from `flip`. `flip` is used to transpose the plot. But
  `swap` is used to swap the x and group_by during plotting.

- ylab:

  The y-axis label. Default is NULL.

- ...:

  Other arguments passed to the specific plot function.

  - For `bar` plot, see
    [`plotthis::BarPlot()`](https://pwwang.github.io/plotthis/reference/barplot.html).

  - For `circos` plot, see
    [`plotthis::CircosPlot()`](https://pwwang.github.io/plotthis/reference/chordplot.html).

  - For `pie` chart, see
    [`plotthis::PieChart()`](https://pwwang.github.io/plotthis/reference/PieChart.html).

  - For `pies` plot, see
    [`plotthis::Heatmap()`](https://pwwang.github.io/plotthis/reference/Heatmap.html).

  - For `heatmap` plot, see
    [`plotthis::Heatmap()`](https://pwwang.github.io/plotthis/reference/Heatmap.html).

  - For `ring`/`donut` plot, see
    [`plotthis::RingPlot()`](https://pwwang.github.io/plotthis/reference/RingPlot.html).

  - For `trend` plot, see
    [`plotthis::TrendPlot()`](https://pwwang.github.io/plotthis/reference/TrendPlot.html).

  - For `area` plot, see
    [`plotthis::AreaPlot()`](https://pwwang.github.io/plotthis/reference/AreaPlot.html).

  - For `sankey`/`alluvial` plot, see
    [`plotthis::SankeyPlot()`](https://pwwang.github.io/plotthis/reference/sankeyplot.html).

  - For `radar` plot, see
    [`plotthis::RadarPlot()`](https://pwwang.github.io/plotthis/reference/radarplot.html).

  - For `spider` plot, see
    [`plotthis::SpiderPlot()`](https://pwwang.github.io/plotthis/reference/radarplot.html).

  - For `violin` plot, see
    [`plotthis::ViolinPlot()`](https://pwwang.github.io/plotthis/reference/boxviolinplot.html).

  - For `box` plot, see
    [`plotthis::BoxPlot()`](https://pwwang.github.io/plotthis/reference/boxviolinplot.html).

## Value

A ggplot object or a list if `combine` is FALSE

## Details

See:

- <https://pwwang.github.io/scplotter/articles/Giotto_Xenium.html>

for examples of using this function with a Giotto object.

And see:

- <https://pwwang.github.io/scplotter/articles/Working_with_anndata_h5ad_files.html>

for examples of using this function with .h5ad files.

## Examples

``` r
# \donttest{
# library(patchwork)
data(ifnb_sub)

# Bar plot
p1 <- CellStatPlot(ifnb_sub)
p2 <- CellStatPlot(ifnb_sub, group_by = "stim")
p3 <- CellStatPlot(ifnb_sub, group_by = "stim", palette = "Set2")
p4 <- CellStatPlot(ifnb_sub, group_by = "stim", position = "stack")
(p1 | p2) / (p3 | p4)


# Fraction of cells
p1 <- CellStatPlot(ifnb_sub, group_by = "stim", frac = "group",
                   ident = "seurat_annotations", x_text_angle = 60)
p2 <- CellStatPlot(ifnb_sub, group_by = "stim", frac = "group")
p3 <- CellStatPlot(ifnb_sub, group_by = "stim", frac = "ident",
                   position = "stack", alpha = .6)
p4 <- CellStatPlot(ifnb_sub, group_by = "stim", frac = "group",
                   swap = TRUE, position = "stack")
(p1 | p2) / (p3 | p4)


# Splitting/Facetting the plot
CellStatPlot(ifnb_sub, split_by = "stim")

CellStatPlot(ifnb_sub, facet_by = "stim")

CellStatPlot(ifnb_sub, facet_by = "stim", facet_nrow = 2)


# Circos plot
CellStatPlot(ifnb_sub, group_by = "stim", plot_type = "circos")

CellStatPlot(ifnb_sub, group_by = "stim", ident = "seurat_annotations",
             plot_type = "circos", labels_rot = TRUE)


# Pie plot
CellStatPlot(ifnb_sub, plot_type = "pie")

CellStatPlot(ifnb_sub, plot_type = "pie", split_by = "stim")


# Ring plot
CellStatPlot(ifnb_sub, plot_type = "ring", group_by = "stim",
             palette = "Spectral")
#> 'frac' is forced to 'group' for 'ring' plot.


# Trend plot
CellStatPlot(ifnb_sub, plot_type = "trend", frac = "group",
             x_text_angle = 90, group_by = c("stim", "seurat_annotations"))
#> Multiple columns are provided in 'group_by'. They will be concatenated into one column.


# Sankey plot
CellStatPlot(ifnb_sub, plot_type = "sankey", group_by = c("seurat_clusters", "stim"),
             links_alpha = .6)
#> Missing alluvia for some stratum combinations.

CellStatPlot(ifnb_sub, plot_type = "sankey", links_alpha = .6,
             group_by = c("stim", "seurat_annotations", "orig.ident"))
#> Missing alluvia for some stratum combinations.

CellStatPlot(ifnb_sub, plot_type = "sankey", links_alpha = .6,
             group_by = c("seurat_clusters", "stim", "seurat_annotations", "orig.ident"))
#> Missing alluvia for some stratum combinations.


# Area plot
CellStatPlot(ifnb_sub, plot_type = "area", frac = "group", x_text_angle = 90,
             group_by = "seurat_annotations", split_by = "stim")


# Pies
# Simulate some sets of cells (e.g. clones)
ifnb_sub$r1 <- ifelse(ifnb_sub$seurat_clusters %in% c("0", "1", "2"), 1, 0)
ifnb_sub$r2 <- sample(c(1, 0), ncol(ifnb_sub), prob = c(0.5, 0.5), replace = TRUE)
ifnb_sub$r3 <- sample(c(1, 0), ncol(ifnb_sub), prob = c(0.7, 0.3), replace = TRUE)
CellStatPlot(ifnb_sub, plot_type = "pies", group_by = "stim", rows_name = "Clones",
   rows_by = c("r1", "r2", "r3"), show_row_names = TRUE, add_reticle = TRUE,
   show_column_names = TRUE, column_names_side = "top", cluster_columns = FALSE,
   row_names_side = "right", pie_size = "identity", pie_values = "sum")


# Heatmap
CellStatPlot(ifnb_sub, plot_type = "heatmap", group_by = "stim", palette = "Blues")

CellStatPlot(ifnb_sub, plot_type = "heatmap", group_by = "stim",
   frac = "group", columns_split_by = "seurat_annotations", swap = TRUE)


# Radar plot/Spider plot
pr <- CellStatPlot(ifnb_sub, plot_type = "radar", group_by = "stim")
ps <- CellStatPlot(ifnb_sub, plot_type = "spider", group_by = "stim")
pr | ps


# Box/Violin plot
ifnb_sub$group <- sample(paste0("g", 1:10), nrow(ifnb_sub), replace = TRUE)
CellStatPlot(ifnb_sub, group_by = c("group", "stim"), frac = "group",
   plot_type = "violin", add_box = TRUE, ident = "seurat_annotations",
   x_text_angle = 60, comparisons = TRUE, aspect.ratio = 0.8)
#> Warning: Computation failed in `stat_pwc()`.
#> Caused by error in `mutate()`:
#> ℹ In argument: `data = map(.data$data, .f, ...)`.
#> Caused by error in `map()`:
#> ℹ In index: 13.
#> Caused by error in `utils::combn()`:
#> ! n < m
#> Warning: Computation failed in `stat_pwc()`.
#> Caused by error in `mutate()`:
#> ℹ In argument: `data = map(.data$data, .f, ...)`.
#> Caused by error in `map()`:
#> ℹ In index: 13.
#> Caused by error in `utils::combn()`:
#> ! n < m

# }
```
