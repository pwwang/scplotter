# clonal_size_data

Function to get the clonal size data for all group_by values.

## Usage

``` r
clonal_size_data(data, clone_call, chain, groupings)
```

## Arguments

- data:

  The product of
  [scRepertoire::combineTCR](https://www.borch.dev/uploads/scRepertoire/reference/combineTCR.html),
  [scRepertoire::combineTCR](https://www.borch.dev/uploads/scRepertoire/reference/combineTCR.html),
  or
  [scRepertoire::combineExpression](https://www.borch.dev/uploads/scRepertoire/reference/combineExpression.html).

- clone_call:

  How to call the clone - VDJC gene (gene), CDR3 nucleotide (nt), CDR3
  amino acid (aa), VDJC gene + CDR3 nucleotide (strict) or a custom
  variable in the data

- chain:

  indicate if both or a specific chain should be used - e.g. "both",
  "TRA", "TRG", "IGH", "IGL"

- groupings:

  The column names in the meta data to group the cells.

## Value

A data frame with the clonal size data.
