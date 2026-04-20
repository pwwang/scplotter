# Get default dimension reduction for Seurat object

Get default dimension reduction for Seurat object

## Usage

``` r
default_dimreduc(object)
```

## Arguments

- object:

  Seurat object to query

## Value

Name of the default dimension reduction, or NULL if not set

## Details

`DefaultDimReduc` was just introduced in Seurat v5.4.0. This function
will first check if we have `DefaultDimReduc` value in the misc slot
(set by our fallback setter), and if not, it will try to call the
official `DefaultDimReduc` getter function. If that also fails, it will
return NULL.
