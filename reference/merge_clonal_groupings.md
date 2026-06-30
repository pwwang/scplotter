# Merge multiple grouping columns into a single composite grouping

Because
[`scRepertoire::clonalQuant()`](https://www.borch.dev/uploads/scRepertoire/reference/clonalQuant.html)
and related functions do not support multiple grouping variables, this
function merges all grouping columns into a single composite `.group`
column (using
[`tidyr::unite()`](https://tidyr.tidyverse.org/reference/unite.html)
with a separator). The original groupings can later be restored by
splitting on the separator.

The function handles both Seurat objects (operating on `@meta.data`) and
the standard scRepertoire list-of-data-frames format. When no groupings
are provided, a placeholder `.group` column set to an empty string is
created, effectively treating all cells as a single group.

## Usage

``` r
merge_clonal_groupings(data, groupings, sep = " // ")
```

## Arguments

- data:

  The product of
  [`scRepertoire::combineTCR()`](https://www.borch.dev/uploads/scRepertoire/reference/combineTCR.html),
  [`scRepertoire::combineBCR()`](https://www.borch.dev/uploads/scRepertoire/reference/combineBCR.html),
  or
  [`scRepertoire::combineExpression()`](https://www.borch.dev/uploads/scRepertoire/reference/combineExpression.html).
  May also be a Seurat object.

- groupings:

  A character vector of column names to merge into the composite
  `.group` column. If empty, `.group` is set to an empty string.

- sep:

  The separator string used to join grouping values. Default: `" // "`.

## Value

The input data with an added `.group` column containing the concatenated
grouping values (or an empty string if no groupings are provided). For
Seurat objects, the column is added to `@meta.data`. For list-format
data, each data frame gains the column.
