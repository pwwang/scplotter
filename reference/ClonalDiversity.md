# Calculate the clonal diversities.

Calculate the clonal diversities.

## Usage

``` r
ClonalDiversity(
  input.data,
  cloneCall = "gene",
  chain = "both",
  method = c("shannon", "gini.coeff", "inv.simpson", "norm.entropy", "gini.simpson",
    "chao1", "ACE", "d50", "dXX"),
  d = 50,
  group_by = NULL,
  n_boots = 0
)
```
