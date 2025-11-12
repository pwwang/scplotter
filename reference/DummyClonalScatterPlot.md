# DummyClonalScatterPlot

Function to plot the scatter plot of the clonal data for a dummy group
pair.

## Usage

``` r
DummyClonalScatterPlot(df, title, group_by, scatter_cor, size_by, ...)
```

## Arguments

- df:

  The data frame with the clonal data. The data frame should have the
  columns: group_by, 'count', 'fraction'.

## Value

A ggplot object
