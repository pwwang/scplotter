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

`scplotter` draws significant inspiration from the [`SCP`][2] package, which separates visualization from analysis and implements it in the [`plotthis`][1] package. Building on `plotthis`, `scplotter` offers advanced functions for visualizing single-cell sequencing data. Special thanks to the [`scRepertoire`][21] package for its APIs that facilitate the analysis of single-cell TCR/BCR sequencing data.

## Gallery

### scRNA-seq

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

[`CCCPlot`][10]

![CCCPlot](./man/figures/cccplot.png)

### scTCR-seq/scBCR-seq

[`ClonalVolumePlot`][11] | [`ClonalAbundancePlot`][12] | [`ClonalResidencyPlot`][13] | [`ClonalCompositionPlot`][14] | [`ClonalOverlapPlot`][15] | [`ClonalGeneUsagePlot`][16]

![clonalstat](./man/figures/clonalstat.png)

[`ClonalRarefactionPlot`][17] | [`ClonalGeneUsagePlot`][18] | [`ClonalDiversityPlot`][19] | [`ClonalPositionalPlot`][20]

![clonaldiv](./man/figures/clonaldiv.png)

[1]: https://github.com/pwwang/plotthis
[2]: https://zhanghao-njmu.github.io/SCP/index.html
[3]: https://pwwang.github.io/scplotter/reference/CellDimPlot.html
[4]: https://pwwang.github.io/scplotter/reference/CellStatPlot.html
[5]: https://pwwang.github.io/scplotter/reference/ClustreePlot.html
[6]: https://pwwang.github.io/scplotter/reference/FeatureStatPlot.html
[7]: https://pwwang.github.io/scplotter/reference/EnrichmentPlot.html
[8]: https://pwwang.github.io/plotthis/reference/gsea.html
[9]: https://pwwang.github.io/plotthis/reference/VolcanoPlot.html
[10]: https://pwwang.github.io/plotthis/reference/CCCPlot.html
[11]: https://pwwang.github.io/scplotter/reference/ClonalVolumePlot.html
[12]: https://pwwang.github.io/scplotter/reference/ClonalAbundancePlot.html
[13]: https://pwwang.github.io/scplotter/reference/ClonalResidencyPlot.html
[14]: https://pwwang.github.io/scplotter/reference/ClonalCompositionPlot.html
[15]: https://pwwang.github.io/scplotter/reference/ClonalOverlapPlot.html
[16]: https://pwwang.github.io/scplotter/reference/ClonalGeneUsagePlot.html
[17]: https://pwwang.github.io/scplotter/reference/ClonalRarefactionPlot.html
[18]: https://pwwang.github.io/scplotter/reference/ClonalGeneUsagePlot.html
[19]: https://pwwang.github.io/scplotter/reference/ClonalDiversityPlot.html
[20]: https://pwwang.github.io/scplotter/reference/ClonalPositionalPlot.html
[21]: https://github.com/ncborcherding/scRepertoire
