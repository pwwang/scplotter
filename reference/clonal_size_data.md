# Prepare clonal abundance data for visualization

This is the central data preparation function used by all clonal
visualization functions in scplotter. It transforms raw scRepertoire
combined TCR/BCR data (a list of data frames or Seurat object) into a
standardized long-format data frame suitable for plotting.

The function performs the following steps:

1.  Extracts grouping levels via
    [`get_clonal_grouping_levels()`](https://pwwang.github.io/scplotter/reference/get_clonal_grouping_levels.md).

2.  Merges multiple grouping columns into a single composite `.group`
    column via
    [`merge_clonal_groupings()`](https://pwwang.github.io/scplotter/reference/merge_clonal_groupings.md).

3.  Reorganizes list-format data by splitting on `.group` and combining
    cells that share the same group across samples.

4.  Calls scRepertoire's internal `.data.wrangle()` and `.theCall()` to
    compute clone-level counts per group.

5.  Pivots the resulting wide-format table (clones × groups) to long
    format with columns `CloneID`, `.group`, `count`, and `fraction`.

6.  Separates the composite `.group` column back into the original
    grouping columns.

7.  Restores original factor levels for each grouping variable.

The output data frame contains one row per clone per group, with columns
for each grouping variable, `CloneID`, `count` (number of cells), and
`fraction` (proportion of cells within the group).

## Usage

``` r
clonal_size_data(data, clone_call, chain, groupings, order)
```

## Arguments

- data:

  The product of
  [`scRepertoire::combineTCR()`](https://www.borch.dev/uploads/scRepertoire/reference/combineTCR.html),
  [`scRepertoire::combineBCR()`](https://www.borch.dev/uploads/scRepertoire/reference/combineBCR.html),
  or
  [`scRepertoire::combineExpression()`](https://www.borch.dev/uploads/scRepertoire/reference/combineExpression.html).
  A named list of data frames where each element represents a sample, or
  a Seurat object with TCR/BCR information in the metadata.

- clone_call:

  How to identify a clone. One of `"gene"` (VDJC gene segment), `"nt"`
  (CDR3 nucleotide sequence), `"aa"` (CDR3 amino acid sequence),
  `"strict"` (VDJC gene + CDR3 nucleotide), or a custom column name
  present in the data.

- chain:

  Which TCR/BCR chain(s) to include. One of `"both"` (both chains
  combined), `"TRA"`, `"TRB"`, `"TRG"`, `"TRD"`, `"IGH"`, `"IGL"`, or
  `"IGK"`.

- groupings:

  A character vector of column names in the metadata to use for grouping
  cells. These define the categories across which clone abundances are
  computed.

- order:

  A named list specifying the order of the levels for each grouping
  variable. Default: `NULL` (uses the order in the data).

## Value

A data frame in long format with columns: the grouping variables,
`CloneID` (clone identifier), `count` (number of cells per clone per
group), and `fraction` (proportion of cells per clone within each
group). Only rows with `count > 0` are retained.
