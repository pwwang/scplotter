# Cell-Cell Communication Plot

Plot the cell-cell communication. See also:

- The review:
  <https://www.sciencedirect.com/science/article/pii/S2452310021000081>

- The LIANA package:
  <https://liana-py.readthedocs.io/en/latest/notebooks/basic_usage.html#Tileplot>

- The CCPlotR package: <https://github.com/Sarah145/CCPlotR>

## Usage

``` r
CCCPlot(
  data,
  plot_type = c("dot", "network", "chord", "circos", "heatmap", "sankey", "alluvial",
    "box", "violin", "ridge"),
  method = c("aggregation", "interaction"),
  magnitude = waiver(),
  specificity = waiver(),
  magnitude_agg = length,
  magnitude_name = "No. of interactions",
  meta_specificity = "sumlog",
  split_by = NULL,
  x_text_angle = 90,
  link_curvature = 0.2,
  link_alpha = 0.6,
  facet_by = NULL,
  show_row_names = TRUE,
  show_column_names = TRUE,
  ...
)
```

## Arguments

- data:

  A data frame with the cell-cell communication data. A typical data
  frame should have the following columns:

  - `source` The source cell type.

  - `target` The target cell type.

  - `ligand` The ligand gene.

  - `receptor` The receptor gene.

  - `ligand_means` The mean expression of the ligand gene per cell type.

  - `receptor_means` The mean expression of the receptor gene per cell
    type.

  - `ligand_props` The proportion of cells that express the entity.

  - `receptor_props` The proportion of cells that express the entity.

  - `<magnitude>` The magnitude of the communication.

  - `<specificity>` The specificity of the communication. Depends on the
    `plot_type`, some columns are optional. But the `source`, `target`,
    `ligand`, `receptor` and `<magnitude>` are required.

- plot_type:

  The type of plot to use. Default is "dot". Possible values are
  "network", "chord", "circos", "heatmap", "sankey", "alluvial", "dot",
  "box", "violin" and "ridge". For "box", "violin" and "ridge", the
  `method` should be "interaction".

  - network: A network plot with the source and target cells as the
    nodes and the communication as the edges.

  - chord: A chord plot with the source and target cells as the nodes
    and the communication as the chords.

  - circos: Alias of "chord".

  - heatmap: A heatmap plot with the source and target cells as the rows
    and columns.

  - sankey: A sankey plot with the source and target cells as the nodes
    and the communication as the flows.

  - alluvial: Alias of "sankey".

  - dot: A dot plot with the source and target cells as the nodes and
    the communication as the dots.

  - box: Box plots for source cell types. Each x is a target cell type
    and the values will be the interaction strengths of the
    ligand-receptor pairs.

  - violin: Violin plots for source cell types. Each x is a target cell
    type and the values will be the interaction strengths of the
    ligand-receptor pairs.

  - ridge: Ridge plots for source cell types. Each row is a target cell
    type and the values will be the interaction strengths of the
    ligand-receptor pairs.

- method:

  The method to determine the plot entities.

  - aggregation: Aggregate the ligand-receptor pairs interactions for
    each source-target pair. Only the source / target pairs will be
    plotted.

  - interaction: Plot the ligand-receptor pairs interactions directly.
    The ligand-receptor pairs will also be plotted.

- magnitude:

  The column name in the data to use as the magnitude of the
  communication. By default, the second last column will be used. See
  `li.mt.show_methods()` for the available methods in `LIANA`. or
  <https://liana-py.readthedocs.io/en/latest/notebooks/basic_usage.html#Tileplot>

- specificity:

  The column name in the data to use as the specificity of the
  communication. By default, the last column will be used. If the method
  doesn't have a specificity, set it to NULL.

- magnitude_agg:

  A function to aggregate the magnitude of the communication. Default is
  `length`.

- magnitude_name:

  The name of the magnitude in the plot. Default is "No. of
  interactions".

- meta_specificity:

  The method to calculate the specificity when there are multiple
  ligand-receptor pairs interactions. Default is "sumlog". It should be
  one of the methods in the `metap` package. Current available methods
  are:

  - `invchisq`: Combine p values using the inverse chi squared method

  - `invt`: Combine p values using the inverse t method

  - `logitp`: Combine p values using the logit method

  - `meanp`: Combine p values by the mean p method

  - `meanz`: Combine p values using the mean z method

  - `sumlog`: Combine p-values by the sum of logs (Fisher's) method

  - `sump`: Combine p-values using the sum of p (Edgington's) method

  - `two2one`: Convert two-sided p-values to one-sided

  - `votep`: Combine p-values by the vote counting method

  - `wilkinsonp`: Combine p-values using Wilkinson's method

- split_by:

  A character vector of column names to split the plots. Default is
  NULL.

- x_text_angle:

  The angle of the x-axis text. Default is 90. Only used when
  `plot_type` is "dot".

- link_curvature:

  The curvature of the links. Default is 0.2. Only used when `plot_type`
  is "network".

- link_alpha:

  The transparency of the links. Default is 0.6. Only used when
  `plot_type` is "network".

- facet_by:

  A character vector of column names to facet the plots. Default is
  NULL. It should always be NULL.

- show_row_names:

  Whether to show the row names in the heatmap. Default is TRUE. Only
  used when `plot_type` is "heatmap".

- show_column_names:

  Whether to show the column names in the heatmap. Default is TRUE. Only
  used when `plot_type` is "heatmap".

- ...:

  Other arguments passed to the specific plot function.

  - For `Network`, see
    [`plotthis::Network()`](https://pwwang.github.io/plotthis/reference/Network.html).

  - For `ChordPlot`, see
    [`plotthis::ChordPlot()`](https://pwwang.github.io/plotthis/reference/chordplot.html).

  - For `Heatmap`, see
    [`plotthis::Heatmap()`](https://pwwang.github.io/plotthis/reference/Heatmap.html).

  - For `SankeyPlot`, see
    [`plotthis::SankeyPlot()`](https://pwwang.github.io/plotthis/reference/sankeyplot.html).

  - For `DotPlot`, see
    [`plotthis::DotPlot()`](https://pwwang.github.io/plotthis/reference/dotplot.html).

## Value

A ggplot object or a list if `combine` is FALSE

## Examples

``` r
# \donttest{
set.seed(8525)
data(cellphonedb_res)
CCCPlot(data = cellphonedb_res, plot_type = "network", legend.position = "none",
  theme = "theme_blank", theme_args = list(add_coord = FALSE))

CCCPlot(cellphonedb_res, plot_type = "chord")

CCCPlot(cellphonedb_res, plot_type = "heatmap")

CCCPlot(cellphonedb_res, plot_type = "dot",
  magnitude_agg = mean, magnitude_name = "Average Interaction Strength")

CCCPlot(cellphonedb_res, plot_type = "sankey")
#> Missing alluvia for some stratum combinations.


cellphonedb_res_sub <- cellphonedb_res[
  cellphonedb_res$source %in% c("Dendritic", "CD14+ Monocyte"),]
CCCPlot(cellphonedb_res_sub, plot_type = "dot", method = "interaction")
#> Multiple columns are provided in 'y'. They will be concatenated into one column.

CCCPlot(cellphonedb_res_sub, plot_type = "network", method = "interaction",
  node_size_by = 1)

CCCPlot(cellphonedb_res_sub, plot_type = "heatmap", method = "interaction",
  palette = "Reds")

CCCPlot(cellphonedb_res_sub, plot_type = "box", method = "interaction")

CCCPlot(cellphonedb_res_sub, plot_type = "violin", method = "interaction",
  add_box = TRUE)

CCCPlot(cellphonedb_res_sub, plot_type = "ridge", method = "interaction")
#> Picking joint bandwidth of 0.311
#> Picking joint bandwidth of 0.285
#> Picking joint bandwidth of 0.311
#> Picking joint bandwidth of 0.285

# }
```
