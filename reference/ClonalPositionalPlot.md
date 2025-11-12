# ClonalPositionalPlot

Visualize the positional entropy, property or amino acid frequency of
CDR3 sequences.

## Usage

``` r
ClonalPositionalPlot(
  data,
  chain = "TRB",
  aa_length = 20,
  group_by = "Sample",
  group_by_sep = "_",
  split_by = NULL,
  method = c("AA", "shannon", "inv.simpson", "norm.entropy", "Atchley", "Kidera",
    "stScales", "tScales", "VHSE"),
  plot_type = c("bar", "line", "heatmap", "box", "violin"),
  theme_args = list(),
  xlab = NULL,
  ylab = NULL,
  facet_by = NULL,
  facet_ncol = NULL,
  facet_nrow = NULL,
  aspect.ratio = NULL,
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

- aa_length:

  The length of the amino acid sequence. Default is 20.

- group_by:

  The variable to group the data by. Default is "Sample".

- group_by_sep:

  The separator to use when combining groupings. Default is "\_".

- split_by:

  The variable to split the data by. Default is NULL.

- method:

  The method to calculate the positional entropy. Default is "AA".

  - "AA": Amino acid frequency.

  - "shannon": Shannon entropy.

  - "inv.simpson": Inverse Simpson index.

  - "norm.entropy": Normalized entropy.

  - "Atchley": Atchley factors.

  - "Kidera": Kidera factors.

  - "stScales": stScales factors.

  - "tScales": tScales factors.

  - "VHSE": Vectors of Hydrophobic, Steric, and Electronic properties.
    See also
    [scRepertoire::percentAA](https://www.borch.dev/uploads/scRepertoire/reference/percentAA.html),
    [scRepertoire::positionalEntropy](https://www.borch.dev/uploads/scRepertoire/reference/positionalEntropy.html)
    and
    [scRepertoire::positionalProperty](https://www.borch.dev/uploads/scRepertoire/reference/positionalProperty.html).

- plot_type:

  The type of plot to generate. Default is "bar".

  - "bar": Bar plot.

  - "line": Line plot.

  - "heatmap": Heatmap.

  - "box": Box plot.

  - "violin": Violin plot.

- theme_args:

  A list of arguments to be passed to the
  [ggplot2::theme](https://ggplot2.tidyverse.org/reference/theme.html)
  function.

- xlab:

  The x-axis label. Default is NULL.

- ylab:

  The y-axis label. Default is NULL.

- facet_by:

  A character vector of column names to facet the plots. Default is
  NULL.

- facet_ncol:

  The number of columns in the facet grid. Default is NULL.

- facet_nrow:

  The number of rows in the facet grid. Default is NULL.

- aspect.ratio:

  The aspect ratio of the plot. Default is NULL.

- ...:

  Other arguments passed to the specific plot function.

  - For "bar",
    [`plotthis::BarPlot()`](https://pwwang.github.io/plotthis/reference/barplot.html).

  - For "line",
    [`plotthis::LinePlot()`](https://pwwang.github.io/plotthis/reference/LinePlot.html).

  - For "heatmap",
    [`plotthis::Heatmap()`](https://pwwang.github.io/plotthis/reference/Heatmap.html).

  - For "box",
    [`plotthis::BoxPlot()`](https://pwwang.github.io/plotthis/reference/boxviolinplot.html).

  - For "violin",
    [`plotthis::ViolinPlot()`](https://pwwang.github.io/plotthis/reference/boxviolinplot.html).

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

ClonalPositionalPlot(data)

ClonalPositionalPlot(data, method = "shannon")

ClonalPositionalPlot(data, method = "norm.entropy", plot_type = "heatmap")

ClonalPositionalPlot(data, method = "Atchley", group_by = "Type", plot_type = "bar")

ClonalPositionalPlot(data, method = "Atchley", plot_type = "line")

# }
```
