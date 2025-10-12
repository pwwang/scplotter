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

## Gallery

### scRNA-seq

[`CellDimPlot`][3]

![](./man/figures/celldimplot.png)

[`CellStatPlot`][4]

![](./man/figures/cellstatplot.png)

[`ClustreePlot`][5] | [`CellVelocityPlot`][2]

![](./man/figures/clustreeplot.png)

[`FeatureStatPlot`][6]

![](./man/figures/featurestatplot.png)

[`EnrichmentPlot`][7]

![](./man/figures/enrichmentplot.png)

[`GSEASummaryPlot`][8] | [`GSEAPlot`][8]

![](./man/figures/gseaplot.png)

[`MarkersPlot`][9]

![](./man/figures/volcanoplot.png)

[`CCCPlot`][10] (Cell-Cell Communication Plot)

![](./man/figures/cccplot.png)

### scTCR-seq/scBCR-seq

[`ClonalVolumePlot`][11] | [`ClonalAbundancePlot`][12] | [`ClonalResidencyPlot`][13] | [`ClonalDynamicsPlot`][22] | [`ClonalCompositionPlot`][14] | [`ClonalOverlapPlot`][15] | [`ClonalGeneUsagePlot`][16]

![](./man/figures/clonalstat.png)

[`ClonalRarefactionPlot`][17] | [`ClonalKmerPlot`][18] | [`ClonalDiversityPlot`][19] | [`ClonalPositionalPlot`][20] | [`ClonalLengthPlot`][29] | [`ClonalStatPlot`][30]

![](./man/figures/clonaldiv.png)

## Spatial data

[`SpatDimPlot`][23] | [`SpatFeaturePlot`][24]

![](./man/figures/spatialplot.png)

## Visualization with LLMs

```r
provider <- tidyprompt::llm_provider_openai(api_key = Sys.getenv("OPENAI_API_KEY"))
chat <- SCPlotterChat$new(provider = provider)
chat$ask("Generate a cell-cell communication plot for the cellphonedb_res data.")
# Tool identified:  CCCPlot
# Data object identified:  cellphonedb_res
# Running tool:  CCCPlot
```

![](./man/figures/scplotter-chat1.png)

```r
chat$ask("Do a heatmap instead")
# Tool identified:  CCCPlot
# Data object identified:  cellphonedb_res
# Running tool:  CCCPlot
```

![](./man/figures/scplotter-chat2.png)

## Credits

`scplotter` is built upon the following fantastic packages:

- [`plotthis`][1] for the core plotting functions.
- [`tidyprompt`][21] for the LLM interface.
- [`Seurat`][25] for the Seurat object support.
- [`LIANA`][26] for the cell-cell communication analysis.
- [`scRepertoire`][27] for the TCR/BCR repertoire analysis.
- [`Giotto`][28] for the spatial data analysis.


[1]: https://github.com/pwwang/plotthis
[2]: https://pwwang.github.io/scplotter/reference/CellVelocityPlot.html
[3]: https://pwwang.github.io/scplotter/reference/CellDimPlot.html
[4]: https://pwwang.github.io/scplotter/reference/CellStatPlot.html
[5]: https://pwwang.github.io/scplotter/reference/ClustreePlot.html
[6]: https://pwwang.github.io/scplotter/reference/FeatureStatPlot.html
[7]: https://pwwang.github.io/scplotter/reference/EnrichmentPlot.html
[8]: https://pwwang.github.io/plotthis/reference/gsea.html
[9]: https://pwwang.github.io/scplotter/reference/MarkersPlot.html
[10]: https://pwwang.github.io/scplotter/reference/CCCPlot.html
[11]: https://pwwang.github.io/scplotter/reference/ClonalVolumePlot.html
[12]: https://pwwang.github.io/scplotter/reference/ClonalAbundancePlot.html
[13]: https://pwwang.github.io/scplotter/reference/ClonalResidencyPlot.html
[14]: https://pwwang.github.io/scplotter/reference/ClonalCompositionPlot.html
[15]: https://pwwang.github.io/scplotter/reference/ClonalOverlapPlot.html
[16]: https://pwwang.github.io/scplotter/reference/ClonalGeneUsagePlot.html
[17]: https://pwwang.github.io/scplotter/reference/ClonalRarefactionPlot.html
[18]: https://pwwang.github.io/scplotter/reference/ClonalKmerPlot.html
[19]: https://pwwang.github.io/scplotter/reference/ClonalDiversityPlot.html
[20]: https://pwwang.github.io/scplotter/reference/ClonalPositionalPlot.html
[21]: https://github.com/tjarkvandemerwe/tidyprompt
[22]: https://pwwang.github.io/scplotter/reference/ClonalDynamicsPlot.html
[23]: https://pwwang.github.io/scplotter/reference/SpatDimPlot.html
[24]: https://pwwang.github.io/scplotter/reference/SpatFeaturePlot.html
[25]: https://satijalab.org/seurat/
[26]: https://github.com/saezlab/liana-py
[27]: https://github.com/BorchLab/scRepertoire
[28]: https://drieslab.github.io/Giotto_website/
[29]: https://pwwang.github.io/scplotter/reference/ClonalLengthPlot.html
[30]: https://pwwang.github.io/scplotter/reference/ClonalStatPlot.html
