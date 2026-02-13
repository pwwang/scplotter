# Utility functions for clone selectors

These utility functions are used internally by the clone selectors to
handle different environments and return the results in the desired
format.

## Usage

``` r
.get_envtype()

.get_data(envtype)

.return_what(out, id, output)

.top_long(n, groups, data, order, id, within, output_within, output)

.top_wide(n, groups, data, order, id, output_within, output)

.sel_long(expr, groups, data, id, within, output_within, output)

.sel_wide(expr, groups, data, id, output_within, output)

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

- expr:

  The expression (in characters) to filter the clones (e.g. "group1 \>
  group2" to select clones where group1 is larger than group2).

- x:

  The name to be backtick-quoted
