# Package index

## scRNA-seq

Functions for plotting single cell RNA-seq data

- [`CellDimPlot()`](https://pwwang.github.io/scplotter/reference/CellDimPlot.md)
  : Cell Dimension Reduction Plot
- [`CellStatPlot()`](https://pwwang.github.io/scplotter/reference/CellStatPlot.md)
  : Cell statistics plot
- [`CellVelocityPlot()`](https://pwwang.github.io/scplotter/reference/CellVelocityPlot.md)
  : Cell Velocity Plot
- [`FeatureStatPlot()`](https://pwwang.github.io/scplotter/reference/FeatureStatPlot.md)
  : Visualize feature expression and statistics across cell groups
- [`ClustreePlot()`](https://pwwang.github.io/scplotter/reference/ClustreePlot.md)
  : Visualize cluster stability across clustering resolutions
- [`EnrichmentPlot()`](https://pwwang.github.io/scplotter/reference/EnrichmentPlot.md)
  : Visualize gene set enrichment and over-representation analysis
  results
- [`CCCPlot()`](https://pwwang.github.io/scplotter/reference/CCCPlot.md)
  : Visualize Cell-Cell Communication (CCC) Interactions
- [`MarkersPlot()`](https://pwwang.github.io/scplotter/reference/MarkersPlot.md)
  : Visualize differential expression markers

## scTCR-seq/scBCR-seq

Functions for plotting single cell TCR-seq/BCR-seq data

- [`ClonalVolumePlot()`](https://pwwang.github.io/scplotter/reference/ClonalVolumePlot.md)
  : Clonal Volume Plot
- [`ClonalAbundancePlot()`](https://pwwang.github.io/scplotter/reference/ClonalAbundancePlot.md)
  : Clonal Abundance Plot
- [`ClonalLengthPlot()`](https://pwwang.github.io/scplotter/reference/ClonalLengthPlot.md)
  : Clonal CDR3 Length Plot
- [`ClonalResidencyPlot()`](https://pwwang.github.io/scplotter/reference/ClonalResidencyPlot.md)
  : Clonal Residency Plot
- [`ClonalCompositionPlot()`](https://pwwang.github.io/scplotter/reference/ClonalCompositionPlot.md)
  : Clonal Composition Plot
- [`ClonalOverlapPlot()`](https://pwwang.github.io/scplotter/reference/ClonalOverlapPlot.md)
  : Clonal Overlap Plot
- [`ClonalStatPlot()`](https://pwwang.github.io/scplotter/reference/ClonalStatPlot.md)
  : Visualize clone abundance, frequency, and dynamics across groups
- [`ClonalDiversityPlot()`](https://pwwang.github.io/scplotter/reference/ClonalDiversityPlot.md)
  : Clonal Diversity Plot
- [`ClonalGeneUsagePlot()`](https://pwwang.github.io/scplotter/reference/ClonalGeneUsagePlot.md)
  : Visualize TCR/BCR gene segment usage
- [`ClonalPositionalPlot()`](https://pwwang.github.io/scplotter/reference/ClonalPositionalPlot.md)
  : Visualize positional properties of CDR3 sequences
- [`ClonalKmerPlot()`](https://pwwang.github.io/scplotter/reference/ClonalKmerPlot.md)
  : Visualize CDR3 k-mer (motif) frequency
- [`ClonalRarefactionPlot()`](https://pwwang.github.io/scplotter/reference/ClonalRarefactionPlot.md)
  : Clonal Rarefaction Plot

## Spatial data visualization

Functions for plotting spatial data

- [`SpatDimPlot()`](https://pwwang.github.io/scplotter/reference/SpatDimPlot.md)
  : Visualize categorical groups on spatial coordinates
- [`SpatFeaturePlot()`](https://pwwang.github.io/scplotter/reference/SpatFeaturePlot.md)
  : Visualize feature expression on spatial coordinates

## Visualizing using LLMs

Functions/Classes for visualizing using LLMs

- [`SCPlotterChat`](https://pwwang.github.io/scplotter/reference/SCPlotterChat.md)
  : LLM-Powered Chat Interface for Single-Cell Visualization

## data

Built-in data sets

- [`ifnb_sub`](https://pwwang.github.io/scplotter/reference/ifnb_sub.md)
  : A subsetted version of 'ifnb' datasets
- [`pancreas_sub`](https://pwwang.github.io/scplotter/reference/pancreas_sub.md)
  : A subsetted version of mouse 'pancreas' datasets
- [`cellphonedb_res`](https://pwwang.github.io/scplotter/reference/cellphonedb_res.md)
  : A toy example of CellPhoneDB output from LIANA

## Clone selectors

Functions for selecting clones

- [`top()`](https://pwwang.github.io/scplotter/reference/clone_selectors.md)
  [`sel()`](https://pwwang.github.io/scplotter/reference/clone_selectors.md)
  [`uniq()`](https://pwwang.github.io/scplotter/reference/clone_selectors.md)
  [`shared()`](https://pwwang.github.io/scplotter/reference/clone_selectors.md)
  [`gt()`](https://pwwang.github.io/scplotter/reference/clone_selectors.md)
  [`ge()`](https://pwwang.github.io/scplotter/reference/clone_selectors.md)
  [`lt()`](https://pwwang.github.io/scplotter/reference/clone_selectors.md)
  [`le()`](https://pwwang.github.io/scplotter/reference/clone_selectors.md)
  [`eq()`](https://pwwang.github.io/scplotter/reference/clone_selectors.md)
  [`ne()`](https://pwwang.github.io/scplotter/reference/clone_selectors.md)
  [`and()`](https://pwwang.github.io/scplotter/reference/clone_selectors.md)
  [`or()`](https://pwwang.github.io/scplotter/reference/clone_selectors.md)
  : Programmatic clone selection for TCR/BCR repertoire analysis

## Re-exports

Functions that are re-exported from plotthis

- [`reexports`](https://pwwang.github.io/scplotter/reference/reexports.md)
  [`GSEASummaryPlot`](https://pwwang.github.io/scplotter/reference/reexports.md)
  [`GSEAPlot`](https://pwwang.github.io/scplotter/reference/reexports.md)
  : Objects exported from other packages
