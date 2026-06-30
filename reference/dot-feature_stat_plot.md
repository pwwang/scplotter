# Internal dispatcher for FeatureStatPlot — plot from a composed data frame

This is the internal workhorse of `FeatureStatPlot`. After the S3
methods (`FeatureStatPlot.Seurat`, `FeatureStatPlot.giotto`, etc.)
extract expression data and metadata into a unified data frame, this
function handles the common logic: pivoting wide feature columns to long
format, optionally downsampling cells within identity groups, and
dispatching to the appropriate plotthis plotting function based on
`plot_type`.

Named feature lists are handled by converting the names to a `rows_data`
data frame with `rows_split_by`, enabling automatic row grouping in
heatmap and dot plots.

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

  A data frame containing the feature expression columns, metadata
  columns (cell identities, groupings), and optionally dimension
  reduction coordinates.

- features:

  A character vector or a named list of character vectors specifying the
  features to plot. Features can be gene names, gene signature scores,
  or any column present in the expression matrix or metadata. Named
  lists (e.g., `list(Beta = c("Ins1", "Ins2"))`) enable automatic row
  grouping in heatmap and dot plots.

- plot_type:

  Character. The type of plot to generate. One of: `"violin"`, `"box"`,
  `"bar"`, `"ridge"`, `"dim"`, `"cor"`, `"heatmap"`, or `"dot"`. See the
  Description section for guidance on choosing a plot type. Default:
  `"violin"`.

- should_shrink:

  Logical. If `TRUE`, the data frame is subset to only the columns
  needed for the plot (features, ident, group_by, split_by, and
  optionally dims). This prevents unnecessary data from being passed
  through the plotting pipeline.

- should_pivot:

  Logical. If `TRUE`, the wide-format feature columns are pivoted to
  long format with `.features` and `.value` columns, and features are
  used for faceting.

- downsample:

  Numeric. Number or fraction of cells to downsample to per identity
  group. Used for `"violin"`, `"box"`, and `"ridge"` plot types to
  reduce overplotting in large datasets:

  - If `downsample > 1`: Exact number of cells per group.

  - If `0 < downsample <= 1`: Fraction of cells per group.

  Downsampling is performed within each identity group (via
  `dplyr::slice_sample(by = ident)`). Default: `NULL` (no downsampling).

- graph:

  Character or `TRUE`. A graph (nearest-neighbor network) name to
  overlay edges between connected cells on the dim plot. For Seurat
  objects, the name should exist in `object@graphs`. For Giotto objects,
  the name should exist in the nearest network list. If `TRUE`, the
  first available graph is used. Only used when `plot_type = "dim"`.
  Default: `NULL` (no edges).

- bg_cutoff:

  Numeric. Expression cutoff for the background in dim plots. Cells with
  expression below this value are shown in the background color
  (typically gray). Set to `-Inf` to color all cells. Only used when
  `plot_type = "dim"`. Default: `0`.

- dims:

  Integer vector of length 2. The dimensions (columns) of the reduction
  to plot on the x and y axes. Only used when `plot_type = "dim"`.
  Default: `1:2`.

- rows_name:

  Character. The column name used to identify feature rows in the
  heatmap and as the key when merging `rows_data`. Only used when
  `plot_type = "heatmap"` or `"dot"`. Default: `"Features"`.

- ident:

  Character. The metadata column name identifying cell groups (e.g.,
  `"CellType"`, `"SubCellType"`, `"seurat_clusters"`). Used as the
  x-axis for violin, box, bar, and ridge plots; as the heatmap columns
  for heatmap and dot plots; and as the grouping variable for
  correlation plots. For Seurat objects, defaults to the active identity
  (`"Identity"`). Required for non-dim plots on Giotto and H5File
  objects. Default: `NULL`.

- agg:

  Function. The aggregation function applied when `plot_type = "bar"`.
  Common choices: `mean` (default), `median`, `sum`. Applied within each
  group defined by `ident`, `group_by`, and `split_by`.

- group_by:

  Character. A metadata column name to further subdivide cells within
  each identity group (e.g., coloring by treatment within cell type).
  Works with `"violin"`, `"box"`, `"bar"`, and `"ridge"` plot types.
  Default: `NULL`.

- pos_only:

  Character. Whether to restrict to cells with positive feature values:

  - `"no"` — Include all cells (default).

  - `"any"` — Include cells where at least one feature is \> 0.

  - `"all"` — Include only cells where all features are \> 0.

  For named feature lists, filtering is applied to all flattened
  features.

- split_by:

  Character vector or `TRUE`. Metadata column name(s) to split the data
  into separate plots (one per unique value). If `TRUE`, splits by the
  features themselves, creating one plot per feature. Multiple columns
  are concatenated. Default: `NULL`.

- facet_by:

  Character vector. Metadata column name(s) to facet the data within a
  plot. Note: for `"violin"`, `"box"`, `"bar"`, and `"ridge"` plot
  types, `facet_by` should always be `NULL` because the plot is already
  faceted by features. Works normally with `"dim"`, `"heatmap"`,
  `"dot"`, and `"cor"` plot types. Default: `NULL`.

- xlab:

  Character. Custom x-axis label. Default: `NULL` (auto-generated based
  on plot type).

- ylab:

  Character. Custom y-axis label. Default: `NULL` (auto-generated based
  on plot type).

- x_text_angle:

  Numeric. Angle (in degrees) for x-axis text labels. Used for
  `"violin"`, `"box"`, and `"bar"` plot types. Default: `NULL` (defaults
  to `45`).

- ...:

  Arguments passed on to
  [`FeatureStatPlot`](https://pwwang.github.io/scplotter/reference/FeatureStatPlot.md)

  `object`

  :   An object containing single-cell expression data. Supported types:
      a Seurat object, a Giotto object, a character path to an `.h5ad`
      file, or an opened `H5File` object from the hdf5r package.

  `spat_unit`

  :   Character. The spatial unit to extract data from. Only applicable
      to Giotto objects. Default: `NULL` (auto-detected by Giotto).

  `feat_type`

  :   Character. The feature type to extract (e.g., `"rna"`, `"dna"`,
      `"protein"`). Only applicable to Giotto objects. Default: `NULL`
      (auto-detected by Giotto).

  `reduction`

  :   Character. Name of the dimensionality reduction to use (e.g.,
      `"umap"`, `"tsne"`, `"pca"`). Required when `plot_type = "dim"`;
      optional for other types where the reduction coordinates can be
      used as feature values. For Seurat objects, defaults to the
      default reduction if available. Default: `NULL`.

  `assay`

  :   Character. The assay name to extract feature data from (e.g.,
      `"RNA"`, `"SCT"`, `"integrated"`). For Seurat objects, passed to
      [`SeuratObject::GetAssayData()`](https://satijalab.github.io/seurat-object/reference/AssayData.html).
      For H5File objects, determines whether to read from `X` (when
      `assay = "RNA"`) or `layers/<assay>`. Not applicable to Giotto
      objects. Default: `NULL`.

  `layer`

  :   Character. The layer name within the assay to extract data from
      (e.g., `"data"`, `"counts"`, `"scale.data"`). For Seurat objects,
      passed to
      [`SeuratObject::GetAssayData()`](https://satijalab.github.io/seurat-object/reference/AssayData.html).
      For Giotto objects, passed to
      [`GiottoClass::getExpression()`](https://giotto-suite.github.io/GiottoClass/reference/getExpression.html).
      Default: `NULL`.

## Value

A ggplot object (or a `patchwork` object when `split_by` generates
multiple plots and `combine = TRUE`), or a list of ggplot objects if
`combine = FALSE`.
