# scplotter - AI Coding Agent Instructions

## Project Overview

`scplotter` is an R package for publication-quality visualizations of
single-cell sequencing data. It’s built as a **wrapper layer** around
[`plotthis`](https://github.com/pwwang/plotthis), extracting data from
Seurat/Giotto objects and .h5ad files to pass to plotthis functions.

**Key Architecture Principle**: This package does NOT implement core
plotting logic. It transforms single-cell data structures into formats
compatible with `plotthis`, which handles the actual visualization.

## Core Components

### Plot Function Categories

1.  **scRNA-seq**: `CellDimPlot`, `CellStatPlot`, `FeatureStatPlot`,
    `EnrichmentPlot`, `MarkersPlot`, `CCCPlot` (cell-cell communication)
2.  **scTCR/BCR-seq**: `Clonal*Plot` family (15+ functions for TCR/BCR
    repertoire analysis via `scRepertoire`)
3.  **Spatial**: `SpatDimPlot`, `SpatFeaturePlot` with multi-platform
    support (Visium, VisiumHD, Xenium, CosMx, CODEX, etc.)
4.  **LLM Integration**: `SCPlotterChat` R6 class using `tidyprompt` for
    natural language plot generation

### Data Flow Pattern

All plot functions follow this structure:

``` r
Plot*() -> UseMethod() -> Plot*.Seurat / Plot*.giotto / Plot*.H5File -> plotthis::PlotFunction()
```

Functions extract embeddings, metadata, graphs from objects, then
delegate to `plotthis`. See `R/celldimplot.R` (lines 1-100) for the
canonical pattern.

## Critical Development Patterns

### 1. Object Support via S3 Methods

Functions must support multiple object types through S3 dispatch: -
`Seurat` objects (via `SeuratObject` package) - `giotto` objects (via
`GiottoClass`) - `.h5ad` files (paths or
[`hdf5r::H5File`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.md)
objects)

Example: `SpatPlot.Seurat`, `SpatPlot.giotto`, `SpatPlot.H5File` in
`R/spatialplot.R`

### 2. Import Strategy

- Use `@importFrom` for specific functions (never `@import` for entire
  packages except R6)
- Import from `plotthis` heavily:
  `@importFrom plotthis DimPlot ChordPlot Heatmap ...`
- Import helper functions via
  [`getFromNamespace()`](https://rdrr.io/r/utils/getFromNamespace.html):
  `check_columns <- getFromNamespace("check_columns", "plotthis")` (see
  `R/utils.R`)
- Tidyverse imports: `@importFrom dplyr filter mutate group_by`,
  `@importFrom tidyr pivot_wider`

### 3. Clone Selector DSL (scRepertoire)

Unique feature: expression-based clone selection in `Clonal*Plot`
functions - Users write: `clones = "top(3)"` or
`clones = "shared() & gt(5)"` - Implementation:
[`rlang::parse_expr()`](https://rlang.r-lib.org/reference/parse_expr.html)
evaluates selectors as functions - Selector functions:
[`top()`](https://pwwang.github.io/scplotter/reference/clone_selectors.md),
[`shared()`](https://pwwang.github.io/scplotter/reference/clone_selectors.md),
[`uniq()`](https://pwwang.github.io/scplotter/reference/clone_selectors.md),
[`gt()`](https://pwwang.github.io/scplotter/reference/clone_selectors.md),
[`ge()`](https://pwwang.github.io/scplotter/reference/clone_selectors.md),
[`lt()`](https://pwwang.github.io/scplotter/reference/clone_selectors.md),
[`le()`](https://pwwang.github.io/scplotter/reference/clone_selectors.md),
[`and()`](https://pwwang.github.io/scplotter/reference/clone_selectors.md),
[`or()`](https://pwwang.github.io/scplotter/reference/clone_selectors.md),
[`sel()`](https://pwwang.github.io/scplotter/reference/clone_selectors.md) -
See `R/clonalutils.R` (lines 360+) and test files in
`tests/testthat/test-clone_selectors_*.R`

### 4. H5AD File Support

- Read .h5ad files using `hdf5r` package
- Helpers:
  [`h5group_to_dataframe()`](https://pwwang.github.io/scplotter/reference/h5group_to_dataframe.md),
  [`h5group_to_matrix()`](https://pwwang.github.io/scplotter/reference/h5group_to_matrix.md)
  in `R/utils.R`
- Handle sparse matrices (CSR format) and categorical data from AnnData
- See vignette: `vignettes/Working_with_anndata_h5ad_files.Rmd`

### 5. Spatial Data Complexity

Spatial functions (`R/spatialplot.R`) handle multiple coordinate
systems: - **Visium**: `imagerow`/`imagecol` coordinates,
scale_factor=1 - **VisiumHD**: Flexible binning, requires `scale_factor`
parameter - **FOV-based** (Xenium, CosMx): `x`/`y` coordinates, optional
molecule points (`nmols`) - **Giotto**: Uses `spat_unit`, `feat_type`,
`spat_loc_name` parameters

Check object type first: `class(object@images[[first_image]])` →
dispatch to appropriate handler

## Build & Test Workflow

### Development Commands (Makefile)

``` fish
make readme        # Build README from README.Rmd
make docs          # Generate documentation + pkgdown site
make install       # Build and install package locally
make test          # Run testthat tests
make notebooks     # Convert notebooks to HTML (use EXECUTE=true to run)
```

### Testing

- Framework: `testthat` (see `tests/testthat/test-*.R`)
- Focus on clone selector logic, not plot outputs
- Run via
  [`devtools::test()`](https://devtools.r-lib.org/reference/test.html)
  or `make test`

### Documentation

- Roxygen2 for function docs (roxygen2 7.3.3)
- pkgdown site config: `_pkgdown.yml` (Bootstrap 5, litera theme)
- Spatial examples in `vignettes/articles/` (16 platform-specific
  tutorials)

### CI/CD (.github/workflows/main.yml)

- Runs on Ubuntu (R 4.4.1)
- Installs system dep: `libglpk40` for igraph
- Builds pkgdown site, deploys to gh-pages
- Requires `OPENAI_API_KEY` secret for LLM vignette

## Package Dependencies

**Core**: `plotthis` (main plotting engine), `scRepertoire` (\>= 2.0.8,
\< 2.3.2), `Seurat` (\>= 5.0.0), `SeuratObject`

**Spatial**: `GiottoClass`, `GiottoData`, `hdf5r`, `terra`

**TCR/BCR**: `scRepertoire`, `iNEXT` (rarefaction), `metap` (\>= 1.11,
p-value combination)

**LLM**: `tidyprompt` (from GitHub: KennispuntTwente/tidyprompt),
`callr`

**Remotes**: `plotthis`, `GiottoClass`, `GiottoData`, `tidyprompt`
installed from GitHub via `Remotes:` field

## Common Pitfalls

1.  **Don’t implement plots**: Call `plotthis::*Plot()`, don’t reinvent
    visualization logic
2.  **Version constraints**: `scRepertoire` locked to 2.0.8-2.3.1 (API
    changes), `ggVennDiagram >= 1.5.0`
3.  **Lazy data**: `LazyData: true` with `xz` compression for included
    datasets (`data/*.rda`)
4.  **Clone selector scope**: Selectors evaluate in caller environment -
    use
    [`rlang::caller_env()`](https://rlang.r-lib.org/reference/stack.html)
    carefully
5.  **Spatial coords**: Always check platform-specific coordinate naming
    (x/y vs imagerow/imagecol)

## When Adding New Plot Functions

1.  Create S3 generic + methods for Seurat/giotto/H5File
2.  Extract necessary data (embeddings, metadata, matrices)
3.  Transform to format expected by `plotthis::*Plot()`
4.  Add `@importFrom` statements for all dependencies
5.  Write examples using included datasets: `pancreas_sub`, `ifnb_sub`,
    `cellphonedb_res`
6.  Add vignette if introducing new data type or analysis workflow

## LLM Chat Integration

`SCPlotterChat` (R6 class in `R/chat.R`) enables: - Natural language
plot requests: `chat$ask("Plot cell-cell communication as heatmap")` -
Auto-detects datasets from `.GlobalEnv` and package data - Maintains
conversation history for context - Tool discovery from package namespace

Uses `tidyprompt` for provider abstraction (OpenAI, Anthropic, etc.)
