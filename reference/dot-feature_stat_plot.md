# Feature statistic plot when given a composed data frame

Feature statistic plot when given a composed data frame

## Usage

``` r
.feature_stat_plot(
  data,
  features,
  plot_type,
  should_shrink,
  should_pivot,
  downsample,
  graph = NULL,
  bg_cutoff = 0,
  dims = 1:2,
  rows_name = "Features",
  ident = NULL,
  agg = mean,
  group_by = NULL,
  pos_only = c("no", "any", "all"),
  split_by = NULL,
  facet_by = NULL,
  xlab = NULL,
  ylab = NULL,
  x_text_angle = NULL,
  ...
)
```

## Arguments

- data:

  A data frame containing the feature data and metadata.

- features:

  A character vector of feature names

- plot_type:

  Type of the plot. It can be "violin", "box", "bar", "ridge", "dim",
  "cor", "heatmap" or "dot"

- downsample:

  A numeric the number of cells in each identity group to downsample to
  for violin, box, or ridge plots. If n \> 1, it is treated as the
  number of cells to downsample to. If 0 \< n \<= 1, it is treated as
  the fraction of cells to downsample to.

- graph:

  Specify the graph name to add edges between cell neighbors to the
  plot, only used when `plot_type` is "dim".

- bg_cutoff:

  Background cutoff for the dim plot, only used when `plot_type` is
  "dim".

- dims:

  Dimensions to plot, only used when `plot_type` is "dim".

- rows_name:

  The name of the rows in the heatmap, only used when `plot_type` is
  "heatmap".

- ident:

  The column name in the meta data to identify the cells.

- agg:

  The aggregation function to use for the bar plot.

- group_by:

  The column name in the meta data to group the cells.

- pos_only:

  Whether to only include cells with positive feature values.

  - "no": Do not filter cells based on feature values. (default)

  - "any": Include cells with positive values for any of the features.

  - "all": Include cells with positive values for all of the features.
    If you have named features (i.e. a named list), `pos_only` will be
    applied to all flattened features.

- split_by:

  Column name in the meta data to split the cells to different plots. If
  TRUE, the cells are split by the features.

- facet_by:

  Column name in the meta data to facet the plots. Should be always
  NULL.

- xlab:

  The x-axis label.

- ylab:

  The y-axis label.

- x_text_angle:

  The angle of the x-axis text. Only used when `plot_type` is "violin",
  "bar", or "box".

- ...:

  Arguments passed on to
  [`FeatureStatPlot`](https://pwwang.github.io/scplotter/reference/FeatureStatPlot.md)

  `object`

  :   A seurat object, a giotto object, a path to an .h5ad file or an
      opened `H5File` by `hdf5r` package.

  `spat_unit`

  :   The spatial unit to use for the plot. Only applied to Giotto
      objects.

  `feat_type`

  :   feature type of the features (e.g. "rna", "dna", "protein"), only
      applied to Giotto objects.

  `reduction`

  :   Name of the reduction to plot (for example, "umap"), only used
      when `plot_type` is "dim" or you can to use the reduction as
      feature.

  `assay`

  :   The assay to use for the feature data.

  `layer`

  :   The layer to use for the feature data.

## Value

A ggplot object or a list if `combine` is FALSE
