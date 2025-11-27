# Subset function with automatic layer repair for Seurat v5 objects

This function handles corrupted layers in Seurat v5 assays that can
cause "incorrect number of dimensions" errors during subsetting.

## Usage

``` r
subset_seurat(object, ...)
```

## Arguments

- object:

  Seurat object to subset

- ...:

  Arguments passed to [`subset()`](https://rdrr.io/r/base/subset.html)

## Value

Subsetted Seurat object with repaired layers if needed

## Details

Strategy:

1.  Evaluate subset normally

2.  Check for corrupted layers proactively

3.  Repair any corrupted layers before subsetting

4.  Perform the subset operation

5.  If subset still fails, re-throw the error

Example: subset_seurat(obj, EFS %in% c("EFS_L", "EFS_S"), nCount_RNA \>
1000)
