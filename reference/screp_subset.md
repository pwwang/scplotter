# Subset an scRepertoire object with a filter expression

Filters an scRepertoire object (either a Seurat object or a named list
of data frames) using a character-specified filter expression. This is
useful for restricting analysis to a subset of cells (e.g., cells with
short CDR3 sequences, or cells from specific samples) without manually
manipulating the underlying data structure.

For list-format data (the standard scRepertoire combined TCR/BCR
format), a `Sample` column is temporarily added to each data frame to
enable sample-based filtering. Empty data frames after filtering are
removed from the result.

## Usage

``` r
screp_subset(screp, subset)
```

## Arguments

- screp:

  An scRepertoire object — either a Seurat object (requires tidyseurat
  for dplyr filter support) or a named list of data frames (the standard
  output of
  [`scRepertoire::combineTCR()`](https://www.borch.dev/uploads/scRepertoire/reference/combineTCR.html)
  or
  [`scRepertoire::combineBCR()`](https://www.borch.dev/uploads/scRepertoire/reference/combineBCR.html)).

- subset:

  A character string containing a filter expression (e.g.,
  `"nchar(CTaa) < 20"` or `"Sample %in% c('P17B', 'P17L')"`). The
  expression is parsed and evaluated within the data context using
  [`rlang::parse_expr()`](https://rlang.r-lib.org/reference/parse_expr.html).

## Value

The filtered scRepertoire object in the same format as the input. For
list-format data, samples with zero rows after filtering are dropped.

## Examples

``` r
# \donttest{
data(contig_list, package = "scRepertoire")
screp <- scRepertoire::combineTCR(contig_list,
    samples = c("P17B", "P17L", "P18B", "P18L", "P19B", "P19L", "P20B", "P20L")
)

head(scplotter:::screp_subset(screp, "nchar(CTaa) < 20")[[1]])
#>                   barcode sample              TCR1      cdr3_aa1
#> 1 P17B_AAACGGGAGAGCCCAA-1   P17B TRAV20.TRAJ8.TRAC CAVRGEGFQKLVF
#> 2 P17B_AAAGTAGAGGCTACGA-1   P17B              <NA>          <NA>
#> 3 P17B_AAATGCCGTCGACTGC-1   P17B              <NA>          <NA>
#> 4 P17B_AACCATGCAACGATGG-1   P17B              <NA>          <NA>
#> 5 P17B_AACCGCGCACTCGACG-1   P17B              <NA>          <NA>
#> 6 P17B_AACTCAGGTCTAACGT-1   P17B              <NA>          <NA>
#>                                  cdr3_nt1                       TCR2
#> 1 TGTGCTGTGCGAGGAGAAGGCTTTCAGAAACTTGTATTT                       <NA>
#> 2                                    <NA> TRBV3-1.None.TRBJ1-1.TRBC1
#> 3                                    <NA> TRBV3-1.None.TRBJ1-1.TRBC1
#> 4                                    <NA>  TRBV27.None.TRBJ2-1.TRBC2
#> 5                                    <NA> TRBV4-3.None.TRBJ1-1.TRBC1
#> 6                                    <NA> TRBV6-2.None.TRBJ2-1.TRBC2
#>           cdr3_aa2                                         cdr3_nt2
#> 1             <NA>                                             <NA>
#> 2   CAAGQGVMNTEAFF       TGTGCCGCGGGGCAGGGGGTCATGAACACTGAAGCTTTCTTT
#> 3   CAAGQGVMNTEAFF       TGTGCCGCGGGGCAGGGGGTCATGAACACTGAAGCTTTCTTT
#> 4 CASSLGSGGTGNEQFF TGTGCCAGCAGTTTAGGGTCGGGGGGGACGGGGAATGAGCAGTTCTTC
#> 5    CASSQDSFTEAFF          TGCGCCAGCAGCCAAGACAGTTTCACTGAAGCTTTCTTT
#> 6 CASSWSKTSGRDEQFF TGTGCCAGCAGTTGGAGTAAGACTAGCGGGAGGGATGAGCAGTTCTTC
#>                          CTgene
#> 1          TRAV20.TRAJ8.TRAC_NA
#> 2 NA_TRBV3-1.None.TRBJ1-1.TRBC1
#> 3 NA_TRBV3-1.None.TRBJ1-1.TRBC1
#> 4  NA_TRBV27.None.TRBJ2-1.TRBC2
#> 5 NA_TRBV4-3.None.TRBJ1-1.TRBC1
#> 6 NA_TRBV6-2.None.TRBJ2-1.TRBC2
#>                                                  CTnt                CTaa
#> 1          TGTGCTGTGCGAGGAGAAGGCTTTCAGAAACTTGTATTT_NA    CAVRGEGFQKLVF_NA
#> 2       NA_TGTGCCGCGGGGCAGGGGGTCATGAACACTGAAGCTTTCTTT   NA_CAAGQGVMNTEAFF
#> 3       NA_TGTGCCGCGGGGCAGGGGGTCATGAACACTGAAGCTTTCTTT   NA_CAAGQGVMNTEAFF
#> 4 NA_TGTGCCAGCAGTTTAGGGTCGGGGGGGACGGGGAATGAGCAGTTCTTC NA_CASSLGSGGTGNEQFF
#> 5          NA_TGCGCCAGCAGCCAAGACAGTTTCACTGAAGCTTTCTTT    NA_CASSQDSFTEAFF
#> 6 NA_TGTGCCAGCAGTTGGAGTAAGACTAGCGGGAGGGATGAGCAGTTCTTC NA_CASSWSKTSGRDEQFF
#>                                                                            CTstrict
#> 1                   TRAV20.TRAJ8.TRAC;TGTGCTGTGCGAGGAGAAGGCTTTCAGAAACTTGTATTT_NA;NA
#> 2       NA;NA_TRBV3-1.None.TRBJ1-1.TRBC1;TGTGCCGCGGGGCAGGGGGTCATGAACACTGAAGCTTTCTTT
#> 3       NA;NA_TRBV3-1.None.TRBJ1-1.TRBC1;TGTGCCGCGGGGCAGGGGGTCATGAACACTGAAGCTTTCTTT
#> 4  NA;NA_TRBV27.None.TRBJ2-1.TRBC2;TGTGCCAGCAGTTTAGGGTCGGGGGGGACGGGGAATGAGCAGTTCTTC
#> 5          NA;NA_TRBV4-3.None.TRBJ1-1.TRBC1;TGCGCCAGCAGCCAAGACAGTTTCACTGAAGCTTTCTTT
#> 6 NA;NA_TRBV6-2.None.TRBJ2-1.TRBC2;TGTGCCAGCAGTTGGAGTAAGACTAGCGGGAGGGATGAGCAGTTCTTC
#>   Sample
#> 1   P17B
#> 2   P17B
#> 3   P17B
#> 4   P17B
#> 5   P17B
#> 6   P17B
names(scplotter:::screp_subset(screp, "Sample %in% c('P17B', 'P17L')"))
#> [1] "P17B" "P17L"
# }
```
