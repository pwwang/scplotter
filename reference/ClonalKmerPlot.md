# ClonalKmerPlot

Explore the k-mer frequency of CDR3 sequences.

## Usage

``` r
ClonalKmerPlot(
  data,
  chain = "TRB",
  clone_call = "aa",
  k = 3,
  top = 25,
  group_by = "Sample",
  group_by_sep = "_",
  facet_by = NULL,
  split_by = NULL,
  plot_type = c("bar", "line", "heatmap"),
  theme_args = list(),
  aspect.ratio = NULL,
  facet_ncol = NULL,
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

- chain:

  The chain to be analyzed. Default is "TRB".

- clone_call:

  The column name of the clone call. Default is "aa".

- k:

  The length of the k-mer. Default is 3.

- top:

  The number of top k-mers to display. Default is 25.

- group_by:

  The variable to group the data by. Default is "Sample".

- group_by_sep:

  The separator to use when combining groupings. Default is "\_".

- facet_by:

  A character vector of column names to facet the plots. Default is
  NULL.

- split_by:

  A character vector of column names to split the plots. Default is
  NULL.

- plot_type:

  The type of plot to generate. Default is "bar".

  - "bar": Bar plot.

  - "line": Line plot.

  - "heatmap": Heatmap.

- theme_args:

  A list of arguments to be passed to the
  [ggplot2::theme](https://ggplot2.tidyverse.org/reference/theme.html)
  function.

- aspect.ratio:

  The aspect ratio of the plot. Default is NULL.

- facet_ncol:

  The number of columns in the facet grid. Default is NULL.

- ...:

  Other arguments passed to the specific plot function.

  - For "bar",
    [`plotthis::BarPlot()`](https://pwwang.github.io/plotthis/reference/barplot.html).

  - For "line",
    [`plotthis::LinePlot()`](https://pwwang.github.io/plotthis/reference/LinePlot.html).

  - For "heatmap",
    [`plotthis::Heatmap()`](https://pwwang.github.io/plotthis/reference/Heatmap.html).

## Value

A ggplot object or a list if `combine` is FALSE

## Examples

``` r
# \donttest{
set.seed(8525)
data(contig_list, package = "scRepertoire")
data <- scRepertoire::combineTCR(contig_list,
    samples = c("P17B", "P17L", "P18B", "P18L", "P19B","P19L", "P20B", "P20L"))
data <- scRepertoire::addVariable(data,
    variable.name = "Type",
    variables = rep(c("B", "L"), 4)
)
data <- scRepertoire::addVariable(data,
    variable.name = "Subject",
    variables = rep(c("P17", "P18", "P19", "P20"), each = 2)
)

ClonalKmerPlot(data)

ClonalKmerPlot(data, group_by = "Type")

ClonalKmerPlot(data, group_by = "Type", plot_type = "line")

ClonalKmerPlot(data, group_by = "Type", plot_type = "heatmap")

# }
```
