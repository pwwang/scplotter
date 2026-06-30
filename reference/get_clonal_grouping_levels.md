# Extract grouping levels from clonal data for later restoration

When clonal data is transformed with multiple groupings (e.g., via
`merge_clonal_groupings`), the original factor levels of the grouping
variables can be lost. This function captures those levels before
transformation so they can be restored later, preserving the intended
ordering of groups in plots and analyses.

The function handles multiple data formats: Seurat objects, single data
frames, and lists of data frames (the standard scRepertoire combined
TCR/BCR format). For list-format data, if all samples share the same
value for a grouping variable, the per-sample values are collected;
otherwise the levels from the first sample are used.

## Usage

``` r
get_clonal_grouping_levels(data, groupings, order = NULL)
```

## Arguments

- data:

  The product of
  [`scRepertoire::combineTCR()`](https://www.borch.dev/uploads/scRepertoire/reference/combineTCR.html),
  [`scRepertoire::combineBCR()`](https://www.borch.dev/uploads/scRepertoire/reference/combineBCR.html),
  or
  [`scRepertoire::combineExpression()`](https://www.borch.dev/uploads/scRepertoire/reference/combineExpression.html).
  May also be a Seurat object or a single data frame.

- groupings:

  A character vector of column names in the metadata to group the cells
  by.

- order:

  A named list specifying the order of the levels for each grouping
  variable. Each element name should match a grouping column, and the
  element value should be a character vector of level names in the
  desired order. Default: `NULL` (use the order present in the data).

## Value

A named list where each element corresponds to a grouping variable and
contains a character vector of its factor levels.
