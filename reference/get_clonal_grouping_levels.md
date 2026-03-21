# Get the grouping levels for clonal data for later to restore

Data transform with multiple groupings may lose the original grouping
levels. This function is to get the grouping levels for later to
restore.

## Usage

``` r
get_clonal_grouping_levels(data, groupings, order = NULL)
```

## Arguments

- data:

  The product of
  [scRepertoire::combineTCR](https://www.borch.dev/uploads/scRepertoire/reference/combineTCR.html),
  [scRepertoire::combineTCR](https://www.borch.dev/uploads/scRepertoire/reference/combineTCR.html),
  or
  [scRepertoire::combineExpression](https://www.borch.dev/uploads/scRepertoire/reference/combineExpression.html).

- groupings:

  The column names in the meta data to group the cells.

- order:

  A list specifying the order of the levels for each grouping variable.

## Value

A list of the grouping levels for each grouping variable.
