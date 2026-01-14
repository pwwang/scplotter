# Utility functions for clone selectors

These utility functions are used internally by the clone selectors to
handle different environments and return the results in the desired
format.

## Usage

``` r
.get_envtype()

.get_data(envtype)

.return_what(out, id, return_ids)

.top_long(n, groups, data, order, id, within, return_ids)

.top_wide(n, groups, data, order, id, return_ids)

.sel_long(expr, groups, data, id, within, return_ids)

.sel_wide(expr, groups, data, id, return_ids)

.bquote(x)
```

## Arguments

- envtype:

  The type of environment to use. It can be "tidy" for dplyr verbs,
  "scplotter" for scplotter functions, or "unknown" if the context
  cannot be determined.

- out:

  The output data frame to be processed.

- id:

  The column name that contains the clone ID. Default is "CTaa".

- return_ids:

  If TRUE, the function returns a vector with the same length as the
  data, with CTaa values for selected clones and NA for others. If
  FALSE, it returns a subset data frame with only the selected clones.
  Default is NULL, which will be determined based on the data. If the
  function is used in a context of dplyr verbs, it defaults to TRUE.
  Otherwise, it defaults to FALSE.

- n:

  The number of top clones to select or the threshold size.

- groups:

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
  the clones by a specific column. Only used in
  [`top()`](https://pwwang.github.io/scplotter/reference/clone_selectors.md).

- within:

  An expression passed to subset the data before applying the selection
  criteria. Note that this subsetting is only applied to determine the
  selection of clones, not to the final output. So if a cell belongs to
  a clone that is selected based on the subsetted data, it will be
  included in the final output, even if it does not meet the `within`
  criteria.

- expr:

  The expression (in characters) to filter the clones (e.g. "group1 \>
  group2" to select clones where group1 is larger than group2).

- x:

  The name to be backtick-quoted
