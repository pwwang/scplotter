# ClonalVolumePlot

ClonalVolumePlot

## Usage

``` r
ClonalVolumePlot(
  data,
  clone_call = "aa",
  chain = "both",
  scale = FALSE,
  plot_type = c("bar", "box", "violin"),
  x = "Sample",
  group_by = NULL,
  facet_by = NULL,
  split_by = NULL,
  order = list(),
  ylab = NULL,
  ...
)
```

## Arguments

- data:

  The product of
  [scRepertoire::combineTCR](https://www.borch.dev/uploads/scRepertoire/reference/combineTCR.html),
  [scRepertoire::combineTCR](https://www.borch.dev/uploads/scRepertoire/reference/combineTCR.html),
  or
  [scRepertoire::combineExpression](https://www.borch.dev/uploads/scRepertoire/reference/combineExpression.html).

- clone_call:

  How to call the clone - VDJC gene (gene), CDR3 nucleotide (nt), CDR3
  amino acid (aa), VDJC gene + CDR3 nucleotide (strict) or a custom
  variable in the data

- chain:

  indicate if both or a specific chain should be used - e.g. "both",
  "TRA", "TRG", "IGH", "IGL"

- scale:

  Whether to use clone proportion or clone size for the plot.

- plot_type:

  The type of plot to use. Default is "bar". Possible values are "bar",
  "box", and "violin". When "box" or "violin" is used, the data will be
  broken down by the Sample and plotted for each group.

- x:

  The column name in the meta data to use as the x-axis. Default:
  "Sample"

- group_by:

  The column name in the meta data to group the cells. Default: NULL

- facet_by:

  The column name in the meta data to facet the plots. Default: NULL

- split_by:

  The column name in the meta data to split the plots. Default: NULL

- order:

  The order of the x-axis items or groups. Default is an empty list. It
  should be a list of values. The names are the column names, and the
  values are the order.

- ylab:

  The y-axis label.

- ...:

  Other arguments passed to the specific plot function.

  - For `bar` plot, see
    [`plotthis::BarPlot()`](https://pwwang.github.io/plotthis/reference/barplot.html).

  - For `box` plot, see
    [`plotthis::BoxPlot()`](https://pwwang.github.io/plotthis/reference/boxviolinplot.html).

  - For `violin` plot, see
    [`plotthis::ViolinPlot()`](https://pwwang.github.io/plotthis/reference/boxviolinplot.html).

## Value

A ggplot object or a list if `combine` is FALSE

## Examples

``` r
# \donttest{
set.seed(8525)
data(contig_list, package = "scRepertoire")
data <- scRepertoire::combineTCR(contig_list)
data <- scRepertoire::addVariable(data,
    variable.name = "Type",
    variables = sample(c("B", "L"), 8, replace = TRUE)
)
data <- scRepertoire::addVariable(data,
    variable.name = "Sex",
    variables = sample(c("M", "F"), 8, replace = TRUE)
)

ClonalVolumePlot(data)

ClonalVolumePlot(data, x = "Type")

ClonalVolumePlot(data, x = "Type", order = list(Type = c("L", "B")))

ClonalVolumePlot(data, x = c("Type", "Sex"), scale = TRUE, fill_by = "Type")
#> Multiple columns are provided in 'x'. They will be concatenated into one column.

ClonalVolumePlot(data, x = "Type", group_by = "Sex", position = "stack")

ClonalVolumePlot(data,
    plot_type = "box", x = "Type", comparisons = TRUE,
    group_by = "Sex"
)

ClonalVolumePlot(data, plot_type = "violin", x = "Type", add_box = TRUE)


# on a Seurat object
data(scRep_example, package = "scRepertoire")
data(contig_list, package = "scRepertoire")
combined <- scRepertoire::combineTCR(contig_list,
    samples = c("P17B", "P17L", "P18B", "P18L", "P19B", "P19L", "P20B", "P20L")
)
sobj <- scRepertoire::combineExpression(combined, scRep_example)
ClonalVolumePlot(sobj)
#> Warning: The 'Sample' column is not found in the meta data, 'orig.indent' will be used instead.

ClonalVolumePlot(sobj, x = "seurat_clusters")
#> Warning: The 'Sample' column is not found in the meta data, 'orig.indent' will be used instead.

ClonalVolumePlot(sobj, group_by = "seurat_clusters")
#> Warning: The 'Sample' column is not found in the meta data, 'orig.indent' will be used instead.

ClonalVolumePlot(sobj, x = "seurat_clusters", plot_type = "box")
#> Warning: The 'Sample' column is not found in the meta data, 'orig.indent' will be used instead.

# }
```
