# Helper functions to select clones based on various criteria

These helper functions allow for the selection of clones based on
various criteria such as size, group comparison, and existence in
specific groups.

## Usage

``` r
top(
  n,
  group_by = NULL,
  data = NULL,
  order = NULL,
  id = NULL,
  in_form = NULL,
  within = NULL,
  output_within = NULL,
  output = NULL
)

sel(
  expr,
  group_by = NULL,
  data = NULL,
  id = NULL,
  in_form = NULL,
  top = NULL,
  order = NULL,
  within = NULL,
  output_within = NULL,
  output = NULL
)

uniq(
  group1,
  group2,
  ...,
  group_by = NULL,
  data = NULL,
  id = NULL,
  in_form = NULL,
  top = NULL,
  order = NULL,
  within = NULL,
  output_within = NULL,
  output = NULL
)

shared(
  group1,
  group2,
  ...,
  group_by = NULL,
  data = NULL,
  id = NULL,
  in_form = NULL,
  top = NULL,
  order = NULL,
  within = NULL,
  output_within = NULL,
  output = NULL
)

gt(
  group1,
  group2,
  include_zeros = TRUE,
  group_by = NULL,
  data = NULL,
  id = NULL,
  in_form = NULL,
  top = NULL,
  order = NULL,
  within = NULL,
  output_within = NULL,
  output = NULL
)

ge(
  group1,
  group2,
  include_zeros = TRUE,
  group_by = NULL,
  data = NULL,
  id = NULL,
  in_form = NULL,
  top = NULL,
  order = NULL,
  within = NULL,
  output_within = NULL,
  output = NULL
)

lt(
  group1,
  group2,
  include_zeros = TRUE,
  group_by = NULL,
  data = NULL,
  id = NULL,
  in_form = NULL,
  top = NULL,
  order = NULL,
  within = NULL,
  output_within = NULL,
  output = NULL
)

le(
  group1,
  group2,
  include_zeros = TRUE,
  group_by = NULL,
  data = NULL,
  id = NULL,
  in_form = NULL,
  top = NULL,
  order = NULL,
  within = NULL,
  output_within = NULL,
  output = NULL
)

eq(
  group1,
  group2,
  group_by = NULL,
  data = NULL,
  id = NULL,
  top = NULL,
  order = NULL,
  in_form = NULL,
  within = NULL,
  output_within = NULL,
  output = NULL
)

ne(
  group1,
  group2,
  include_zeros = TRUE,
  group_by = NULL,
  data = NULL,
  id = NULL,
  in_form = NULL,
  top = NULL,
  order = NULL,
  within = NULL,
  output_within = NULL,
  output = NULL
)

and(x, y, ...)

or(x, y, ...)
```

## Arguments

- n:

  The number of top clones to select or the threshold size.

- group_by:

  The column names in the meta data to group the cells. By default, it
  is assumed `facet_by` and `split_by` to be in the parent frame if used
  in scplotter functions. When used in dplyr verbs, it should be a
  character vector of the grouping columns, where the first column is
  used to extract the values (count) for `group1`, `group2`, and `...`
  and the rest are used to keep the groupings.

- data:

  The data frame containing clone information. Default is NULL. If NULL,
  when used in scplotter functions, it will get data from parent.frame.
  A typical `data` should have a column named `CloneID` and other
  columns for the groupings. Supposingly it should be a grouped data
  frame with the grouping columns. Under each grouping column, the value
  should be the size of the clone. By default, the data is assumed to be
  in the parent frame. When used in dplyr verbs, it should be the parent
  data frame passed to the dplyr verb.

- order:

  The order of the clones to select. It can be an expression to order
  the clones by a specific column. Only used in `top()`.

- id:

  The column name that contains the clone ID. Default is "CTaa".

- in_form:

  The format of the input data. It can be "long" or "wide". If "long",
  the data should be in a long format with a column for the clone IDs
  and a column for the size. If "wide", the data should be in a wide
  format with columns for the clone IDs and the size for each group.
  When used in dplyr verbs, it should be "long" by default. If used in
  scplotter functions, it should be "wide" by default.#'

- within:

  An expression passed to subset the data before applying the selection
  criteria. Only works for `long` format. It is useful when you want to
  select clones based on the criteria within a specific subset of the
  data. Note that this subsetting is only applied to determine the
  selection of clones, not to the final output. So if a cell belongs to
  a clone that is selected based on the subsetted data, it will be
  included in the final output, even if it does not meet the `within`
  criteria. If you want the clones returned to also meet the `within`
  criteria, you can set `output_within` to TRUE, which will return the
  clones that meet both the selection criteria and the `within`
  criteria.

- output_within:

  An expression passed to subset the data after applying the selection
  criteria. Can work with both `long` and `wide` format. It is useful
  when you want to return clones that meet both the selection criteria
  and this criteria. If set to TRUE (only works when `within` is
  specified), the `within` criteria will be applied to filter the final
  output to include only the clones that meet both the selection
  criteria and the `within` criteria. If FALSE or NULL (default), the
  `within` criteria will only be applied to determine the selection of
  clones, not to the final output.

- output:

  There are three options for the output: "id" (or "ids"), "logical" (or
  "bool", "boolean", "indicator"), and "data" (or "df", "data.frame").

  - "id" (or "ids"): return a vector with the same length as the input
    data, with the selected clones' CTaa values (clone IDs) and NA for
    others. It is useful for adding a new column to the data frame.

  - "logical" (or "bool", "boolean", "indicator"): return a logical
    vector indicating whether each clone is selected or not. Same as
    `id` but with TRUE for selected clones and FALSE for others.

  - "data" (or "df", "data.frame"): return a subset of the data frame
    with only the selected clones. This is useful for filtering the data
    frame to only include the clones that meet the criteria. It is used
    internally in some other scplotter functions, such as
    `ClonalStatPlot`, to select clones. The groupings are also applied,
    and defaulting to `facet_by` and `split_by` in the parent frame. By
    default, it is NULL, which will return "id" when used in dplyr verbs
    and "data" when used in scplotter functions.

- expr:

  The expression (in characters) to filter the clones (e.g. "group1 \>
  group2" to select clones where group1 is larger than group2).

- top:

  The number of top clones to select based on the expression. If
  specified, it will select the top N clones that meet the criteria
  defined by `expr` and ordered by `order`. If `order` is not specified,
  it will select the top N clones based on the order they appear in the
  data after filtering by `expr`.

- group1:

  The first group to compare.

- group2:

  The second group to compare.

- ...:

  Additional vectors to compare in logical operations (and/or).

- include_zeros:

  Whether to include clones with zero size in the comparison. If TRUE,
  in a comparison (s1 \> s2) for a clone to be selected, both s1 and s2
  must be greater than 0. If FALSE, only the first group must be greater
  than the second group.

- x:

  The first vector to compare in logical operations (and/or).

- y:

  The second vector to compare in logical operations (and/or).

## Value

A vector of CTaas or a data frame with the selected clones based on the
criteria.

## Details

These helper functions are designed to be used in a dplyr pipeline or
used internally in other scplotter functions to select clones based on
various criteria.

- When used in a dplyr pipeline, they will return a vector with the same
  length as the input data, with the selected clones' CTaa values (clone
  IDs) and NA for others. It is useful for adding a new column to the
  data frame. For the functions that need `group1`, `group2`, and/or
  `...`, `group_by` should be provided to specify the grouping columns.
  Then `group1`, `group2`, and `...` can be the values in the grouping
  column. To include more grouping columns, just use
  `c(grouping1, grouping2, ...)`, where `grouping1` is used for values
  of `group1`, `group2` and `...`; `grouping2` and so on will be kept as
  the groupings where the clones are selected in each combination of the
  grouping values.

- When used in a scplotter function, they will return a subset of the
  data frame with only the selected clones. This is useful for filtering
  the data frame to only include the clones that meet the criteria. It
  is used internally in some other scplotter functions, such as
  `ClonalStatPlot`, to select clones. The groupings are also applied,
  and defaulting to `facet_by` and `split_by` in the parent frame.

- When used independently, you should pass the arguments explicitly,
  such as `group_by` and `output`, to control the behavior and the
  output of the function.

## Examples

``` r
set.seed(8525)
data <- data.frame(
    CTaa = c("AA1", "AA2", "AA3", "AA4", "AA5", "AA6", "AA7", "AA8", "AA9", "AA10"),
    group1 = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
    group2 = c(7, 3, 8, 2, 1, 5, 9, 4, 6, 0),
    groups = c("A", "A", "A", "A", "B", "B", "B", "B", "B", "B")
)
data <- data[order(data$group1 + data$group2, decreasing = TRUE), ]
top(3)
#>   CTaa group1 group2 groups
#> 1  AA7      6      9      B
#> 2  AA9      8      6      B
#> 3  AA8      7      4      B
top(3, group_by = "groups")
#> # A tibble: 6 × 4
#> # Groups:   groups [2]
#>   CTaa  group1 group2 groups
#>   <chr>  <dbl>  <dbl> <chr> 
#> 1 AA7        6      9 B     
#> 2 AA9        8      6 B     
#> 3 AA8        7      4 B     
#> 4 AA3        2      8 A     
#> 5 AA1        0      7 A     
#> 6 AA4        3      2 A     
sel(group1 == 0 | group2 == 0)
#>   CTaa group1 group2 groups
#> 1 AA10      9      0      B
#> 2  AA1      0      7      A
uniq(group1, group2)
#>   CTaa group1 group2 groups
#> 1 AA10      9      0      B
shared(group1, group2)
#>   CTaa group1 group2 groups
#> 1  AA7      6      9      B
#> 2  AA9      8      6      B
#> 3  AA8      7      4      B
#> 4  AA3      2      8      A
#> 5  AA6      5      5      B
#> 6  AA4      3      2      A
#> 7  AA5      4      1      B
#> 8  AA2      1      3      A
gt(group1, group2)
#>   CTaa group1 group2 groups
#> 1  AA9      8      6      B
#> 2  AA8      7      4      B
#> 3 AA10      9      0      B
#> 4  AA4      3      2      A
#> 5  AA5      4      1      B
lt(group1, group2)
#>   CTaa group1 group2 groups
#> 1  AA7      6      9      B
#> 2  AA3      2      8      A
#> 3  AA1      0      7      A
#> 4  AA2      1      3      A
le(group1, group2)
#>   CTaa group1 group2 groups
#> 1  AA7      6      9      B
#> 2  AA3      2      8      A
#> 3  AA6      5      5      B
#> 4  AA1      0      7      A
#> 5  AA2      1      3      A
lt(group1, group2, include_zeros = FALSE)
#>   CTaa group1 group2 groups
#> 1  AA7      6      9      B
#> 2  AA3      2      8      A
#> 3  AA2      1      3      A
eq(group1, group2)
#>   CTaa group1 group2 groups
#> 1  AA6      5      5      B
ne(group1, group2)
#>   CTaa group1 group2 groups
#> 1  AA7      6      9      B
#> 2  AA9      8      6      B
#> 3  AA8      7      4      B
#> 4  AA3      2      8      A
#> 5 AA10      9      0      B
#> 6  AA1      0      7      A
#> 7  AA4      3      2      A
#> 8  AA5      4      1      B
#> 9  AA2      1      3      A

# Use them in a dplyr pipeline
data <- tidyr::pivot_longer(data,
    cols = c("group1", "group2"),
    names_to = "group", values_to = "value"
)
data <- tidyr::uncount(data, !!rlang::sym("value"))
data$subset <- sample(c("S1", "S2"), nrow(data), replace = TRUE)
# Take a glimpse of the data
data[sample(1:nrow(data), 10), ]
#> # A tibble: 10 × 4
#>    CTaa  groups group  subset
#>    <chr> <chr>  <chr>  <chr> 
#>  1 AA8   B      group2 S1    
#>  2 AA8   B      group1 S1    
#>  3 AA5   B      group1 S1    
#>  4 AA10  B      group1 S1    
#>  5 AA7   B      group1 S2    
#>  6 AA7   B      group1 S1    
#>  7 AA1   A      group2 S2    
#>  8 AA1   A      group2 S1    
#>  9 AA5   B      group1 S1    
#> 10 AA9   B      group1 S2    

unique(dplyr::mutate(data, Top3 = top(3))$Top3)
#> [1] "AA7" "AA9" "AA8" NA   
# Note that AA9 also reported in S2, even though the comparison is only applied within S1
dplyr::distinct(
    dplyr::mutate(data, Top3 = top(3, within = subset == "S1")),
    CTaa, subset, Top3
)
#> # A tibble: 19 × 3
#>    CTaa  subset Top3 
#>    <chr> <chr>  <chr>
#>  1 AA7   S2     AA7  
#>  2 AA7   S1     AA7  
#>  3 AA9   S1     AA9  
#>  4 AA9   S2     AA9  
#>  5 AA8   S2     NA   
#>  6 AA8   S1     NA   
#>  7 AA3   S2     NA   
#>  8 AA3   S1     NA   
#>  9 AA6   S1     NA   
#> 10 AA6   S2     NA   
#> 11 AA10  S2     AA10 
#> 12 AA10  S1     AA10 
#> 13 AA1   S1     NA   
#> 14 AA1   S2     NA   
#> 15 AA4   S1     NA   
#> 16 AA4   S2     NA   
#> 17 AA5   S1     NA   
#> 18 AA2   S1     NA   
#> 19 AA2   S2     NA   
# Note that AA9 is now excluded
dplyr::distinct(
    dplyr::mutate(data, Top3 = top(3, within = subset == "S1", output_within = TRUE)),
    CTaa, subset, Top3
)
#> # A tibble: 19 × 3
#>    CTaa  subset Top3 
#>    <chr> <chr>  <chr>
#>  1 AA7   S2     NA   
#>  2 AA7   S1     AA7  
#>  3 AA9   S1     AA9  
#>  4 AA9   S2     NA   
#>  5 AA8   S2     NA   
#>  6 AA8   S1     NA   
#>  7 AA3   S2     NA   
#>  8 AA3   S1     NA   
#>  9 AA6   S1     NA   
#> 10 AA6   S2     NA   
#> 11 AA10  S2     NA   
#> 12 AA10  S1     AA10 
#> 13 AA1   S1     NA   
#> 14 AA1   S2     NA   
#> 15 AA4   S1     NA   
#> 16 AA4   S2     NA   
#> 17 AA5   S1     NA   
#> 18 AA2   S1     NA   
#> 19 AA2   S2     NA   
# We can also exclude S1 clones even when the comparison is applied within S1
dplyr::distinct(
    dplyr::mutate(data, Top3 = top(3, within = subset == "S1", output_within = subset == "S2")),
    CTaa, subset, Top3
)
#> # A tibble: 19 × 3
#>    CTaa  subset Top3 
#>    <chr> <chr>  <chr>
#>  1 AA7   S2     AA7  
#>  2 AA7   S1     NA   
#>  3 AA9   S1     NA   
#>  4 AA9   S2     AA9  
#>  5 AA8   S2     NA   
#>  6 AA8   S1     NA   
#>  7 AA3   S2     NA   
#>  8 AA3   S1     NA   
#>  9 AA6   S1     NA   
#> 10 AA6   S2     NA   
#> 11 AA10  S2     AA10 
#> 12 AA10  S1     NA   
#> 13 AA1   S1     NA   
#> 14 AA1   S2     NA   
#> 15 AA4   S1     NA   
#> 16 AA4   S2     NA   
#> 17 AA5   S1     NA   
#> 18 AA2   S1     NA   
#> 19 AA2   S2     NA   
unique(dplyr::mutate(data, Top3 = top(3, group_by = "groups"))$Top3)
#> [1] "AA7" "AA9" "AA8" "AA3" NA    "AA1" "AA4"
unique(dplyr::mutate(data, Unique = sel(group1 == 0 | group2 == 0, group_by = "group"))$Unique)
#> [1] NA     "AA10" "AA1" 
unique(dplyr::mutate(data, UniqueInG1 = uniq(group1, group2, group_by = "group"))$UniqueInG1)
#> [1] NA     "AA10"
unique(dplyr::mutate(data, Shared = shared(group1, group2, group_by = "group"))$Shared)
#> [1] "AA7" "AA9" "AA8" "AA3" "AA6" NA    "AA4" "AA5" "AA2"
unique(dplyr::mutate(data, Greater = gt(group1, group2, group_by = "group"))$Greater)
#> [1] NA     "AA9"  "AA8"  "AA10" "AA4"  "AA5" 
unique(dplyr::mutate(data, Less = lt(group1, group2, group_by = "group"))$Less)
#> [1] "AA7" NA    "AA3" "AA1" "AA2"
unique(dplyr::mutate(data, LessEqual = le(group1, group2, group_by = "group"))$LessEqual)
#> [1] "AA7" NA    "AA3" "AA6" "AA1" "AA2"
unique(dplyr::mutate(data, GreaterEqual = ge(group1, group2, group_by = "group"))$GreaterEqual)
#> [1] NA     "AA9"  "AA8"  "AA6"  "AA10" "AA4"  "AA5" 
unique(dplyr::mutate(data, Equal = eq(group1, group2, group_by = "group"))$Equal)
#> [1] NA    "AA6"
unique(dplyr::mutate(data, NotEqual = ne(group1, group2, group_by = "group"))$NotEqual)
#>  [1] "AA7"  "AA9"  "AA8"  "AA3"  NA     "AA10" "AA1"  "AA4"  "AA5"  "AA2" 
# Compond expressions
unique(
    dplyr::mutate(data,
        Top3OrEqual = or(top(3), eq(group1, group2, group_by = "group"))
    )$Top3OrEqual
)
#> [1] "AA7" "AA9" "AA8" NA    "AA6"

unique(
    dplyr::mutate(data,
        SharedAndGreater = and(
            shared(group1, group2, group_by = "group"),
            gt(group1, group2, group_by = "group")
        )
    )$SharedAndGreater
)
#> [1] NA    "AA9" "AA8" "AA4" "AA5"
dplyr::mutate(data,
    SharedAndGreater = and(
        shared(group1, group2, group_by = "group", output = "logical"),
        gt(group1, group2, group_by = "group", output = "logical")
    )
)$SharedAndGreater
#>  [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#> [13] FALSE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
#> [25]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
#> [37]  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#> [49] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#> [61] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#> [73] FALSE FALSE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
#> [85]  TRUE  TRUE FALSE FALSE FALSE FALSE
```
