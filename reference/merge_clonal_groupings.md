# merge_clonal_groupings

Merge the multiple clonal groupings into a single grouping.

## Usage

``` r
merge_clonal_groupings(data, groupings, sep = " // ")
```

## Arguments

- data:

  The product of
  [scRepertoire::combineTCR](https://www.borch.dev/uploads/scRepertoire/reference/combineTCR.html),
  [scRepertoire::combineTCR](https://www.borch.dev/uploads/scRepertoire/reference/combineTCR.html),
  or
  [scRepertoire::combineExpression](https://www.borch.dev/uploads/scRepertoire/reference/combineExpression.html).

- groupings:

  A list of the clonal groupings. Each element is a column in the data.

## Value

The data with the combined groupings (`.group`)

## Details

Because
[scRepertoire::clonalQuant](https://www.borch.dev/uploads/scRepertoire/reference/clonalQuant.html)
and families don't support mutliple groupings, this is trying to merge
the multiple groupings into a single grouping. And then later restore
the original groupings.
