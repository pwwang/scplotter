# A toy example of CellPhoneDB output from LIANA

The dataset is generated using python package LIANA

## Format

A `data.frame` with 10 rows and 6 columns

## Source

<https://liana-py.readthedocs.io/en/latest/notebooks/basic_usage.html#Tileplot>

## Examples

``` r
if (FALSE) { # \dontrun{
# Python code
# # import liana
# import liana as li
# # needed for visualization and toy data
# import scanpy as sc
#
# from liana.method import cellphonedb
# adata = sc.datasets.pbmc68k_reduced()
# cellphonedb(adata,
#             groupby='bulk_labels',
#             # NOTE by default the resource uses HUMAN gene symbols
#             resource_name='consensus',
#             expr_prop=0.1,
#             verbose=True, key_added='cpdb_res')
# cellphonedb_res = adata.uns['cpdb_res']
# cellphonedb_res = cellphonedb_res[cellphonedb_res['cellphone_pvals'] < 0.05]
} # }
```
