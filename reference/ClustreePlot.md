# Clustree plot

This function generates a clustree plot from a data frame or a Seurat
object.

## Usage

``` r
ClustreePlot(object, ...)
```

## Arguments

- object:

  The data frame or Seurat object

- ...:

  Other arguments passed to
  [`plotthis::ClustreePlot()`](https://pwwang.github.io/plotthis/reference/ClustreePlot.html)

## Value

A ggplot object or a list if `combine` is FALSE

## Examples

``` r
data(ifnb_sub)
ClustreePlot(ifnb_sub, prefix = "RNA_snn_res.")
#> Registered S3 method overwritten by 'gglogger':
#>   method from   
#>   +.gg   ggplot2
```
