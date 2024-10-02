# scplotter <a href="https://pwwang.github.io/scplotter/"><img src="man/figures/logo.png" align="right" height="139" alt="scplotter website" /></a>

`scplotter` is an R package that is built upon [`plotthis`][1]. It provides a set of functions to visualize single-cell sequencing data in an easy and efficient way.

## Installation

```r
remotes::install_github("pwwang/scplotter")
# or
devtools::install_github("pwwang/scplotter")
# or using conda
$ conda install pwwang::r-scplotter
```

## Credits

`scplotter` is greatly inspired by the [`SCP`][2] package, where the visualization was detached from the analysis and implemented in the [`plotthis`][1] package. The `scplotter` package is built upon the `plotthis` package and provides a set of even higher-level functions to visualize single-cell sequencing data.

## Gallery

[`CellDimPlot`][3]

![CellDimPlot](./man/figures/celldimplot.png)

[`CellStatPlot`][4]

![CellStatPlot](./man/figures/cellstatplot.png)

[`ClustreePlot`][5]

![ClustreePlot](./man/figures/clustreeplot.png)

[`FeatureStatPlot`][6]

![FeatureStatPlot](./man/figures/featurestatplot.png)

[`EnrichmentPlot`][7]

![EnrichmentPlot](./man/figures/enrichmentplot.png)

[`GSEASummaryPlot`][8] | [`GSEAPlot`][8]

![GSEAPlot](./man/figures/gseaplot.png)

[`VolcanoPlot`][9]

![VolcanoPlot](./man/figures/volcanoplot.png)

[1]: https://github.com/pwwang/plotthis
[2]: https://zhanghao-njmu.github.io/SCP/index.html
[3]: https://pwwang.github.io/scplotter/reference/CellDimPlot.html
[4]: https://pwwang.github.io/scplotter/reference/CellStatPlot.html
[5]: https://pwwang.github.io/scplotter/reference/ClustreePlot.html
[6]: https://pwwang.github.io/scplotter/reference/FeatureStatPlot.html
[7]: https://pwwang.github.io/scplotter/reference/EnrichmentPlot.html
[8]: https://pwwang.github.io/plotthis/reference/gsea.html
[9]: https://pwwang.github.io/plotthis/reference/VolcanoPlot.html
