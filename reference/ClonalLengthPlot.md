# Clonal CDR3 Length Plot

Visualizes the distribution of CDR3 sequence lengths across the immune
repertoire. CDR3 length is a key feature of T-cell and B-cell receptor
diversity — different clones have different CDR3 lengths, and shifts in
length distribution can indicate clonal selection, antigen-specific
expansion, or repertoire bias.

`ClonalLengthPlot` computes CDR3 length data via
[`scRepertoire::clonalLength()`](https://www.borch.dev/uploads/scRepertoire/reference/clonalLength.html)
and visualizes the distribution as bar, box, violin, or density plots.
Length is measured in amino acids (when `clone_call = "aa"`) or
nucleotides (when `clone_call = "nt"`).

## Usage

``` r
ClonalLengthPlot(
  data,
  clone_call = "aa",
  chain = "both",
  plot_type = c("bar", "box", "violin", "density"),
  x_nbreaks = 10,
  group_by = "Sample",
  order = NULL,
  xlab = "Length",
  ylab = NULL,
  position = "dodge",
  facet_by = NULL,
  split_by = NULL,
  ...
)
```

## Arguments

- data:

  The product of
  [`scRepertoire::combineTCR()`](https://www.borch.dev/uploads/scRepertoire/reference/combineTCR.html),
  [`scRepertoire::combineBCR()`](https://www.borch.dev/uploads/scRepertoire/reference/combineBCR.html),
  or
  [`scRepertoire::combineExpression()`](https://www.borch.dev/uploads/scRepertoire/reference/combineExpression.html).

- clone_call:

  How to define a clone. Only `"nt"` (CDR3 nucleotide length) or `"aa"`
  (CDR3 amino acid length, default) are supported.

- chain:

  Which chain(s) to use: `"both"` (default), `"TRA"`, `"TRB"`, `"TRD"`,
  `"TRG"`, `"IGH"`, or `"IGL"`.

- plot_type:

  The visualization type. One of `"bar"` (default), `"box"`, `"violin"`,
  or `"density"`.

  - `"bar"` — Bar chart of clone counts at each CDR3 length. Empty
    length bins (zero clones) are padded with zeros to maintain a
    continuous x-axis.

  - `"box"` — Box plot of per-group length distributions.

  - `"violin"` — Violin plot of per-group length distributions.

  - `"density"` — Kernel density estimate of the length distribution,
    using raw (unaggregated) data.

- x_nbreaks:

  Number of x-axis breaks for the bar plot. Default is `10`. Breaks are
  computed as quantiles of the length range and rounded to the nearest
  10.

- group_by:

  Metadata column used to group (color) the data. Default is `"Sample"`.

- order:

  A named list controlling the order of factor levels. List names are
  column names; list values are the desired order. Default is `NULL`.

- xlab:

  X-axis label. Default is `"Length"`.

- ylab:

  Y-axis label. Default is `NULL`, which auto-generates
  `"Number of CDR3 (AA)"` or `"Number of CDR3 (NT)"` based on
  `clone_call`.

- position:

  Bar position for the bar plot. One of `"dodge"` (default), `"stack"`,
  or `"fill"`.

- facet_by:

  Metadata column used to facet the plot into separate panels. Default
  is `NULL`.

- split_by:

  Metadata column used to split the data into separate plots. Default is
  `NULL`.

- ...:

  Additional arguments passed to the underlying plotthis function:

  - `"bar"` —
    [`plotthis::BarPlot()`](https://pwwang.github.io/plotthis/reference/barplot.html)
    (`palette`, `alpha`, `position_dodge_preserve`, ...)

  - `"box"` —
    [`plotthis::BoxPlot()`](https://pwwang.github.io/plotthis/reference/boxviolinplot.html)
    (`comparisons`, `alpha`, `palette`, ...)

  - `"violin"` —
    [`plotthis::ViolinPlot()`](https://pwwang.github.io/plotthis/reference/boxviolinplot.html)
    (`add_box`, `comparisons`, `palette`, ...)

  - `"density"` —
    [`plotthis::DensityPlot()`](https://pwwang.github.io/plotthis/reference/densityhistoplot.html)
    (`palette`, `alpha`, `bw`, ...)

## Value

A `ggplot` object, or a list of `ggplot` objects if `combine = FALSE` is
passed via `...`.

## Examples

``` r
# \donttest{
set.seed(8525)
data(contig_list, package = "scRepertoire")
data <- scRepertoire::combineTCR(contig_list)
data <- scRepertoire::addVariable(data, variable.name = "Type",
 variables = factor(sample(c("B", "L"), 8, replace = TRUE), levels = c("L", "B")))
data <- scRepertoire::addVariable(data, variable.name = "Sex",
 variables = factor(sample(c("M", "F"), 8, replace = TRUE), levels = c("M", "F")))

ClonalLengthPlot(data)

ClonalLengthPlot(data, plot_type = "box")

ClonalLengthPlot(data, clone_call = "nt", plot_type = "violin", chain = "TRB",
 group_by = "Type", comparisons = TRUE)
#> Warning: [Box/Violin/BeeswarmPlot] Some pairwise comparisons may fail due to insufficient data points or variability. Adjusting data to ensure valid comparisons.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.

ClonalLengthPlot(data, plot_type = "density", chain = "TRA")

# }
```
