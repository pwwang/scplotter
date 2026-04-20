# Set default dimension reduction for Seurat object

Set default dimension reduction for Seurat object

## Usage

``` r
default_dimreduc(object) <- value
```

## Arguments

- object:

  Seurat object to modify

- value:

  Name of the default dimension reduction to set

## Value

Modified Seurat object with updated default dimension reduction

## Details

`DefaultDimReduc<-` was just introduced in Seurat v5.4.0 See:
https://github.com/satijalab/seurat-object/pull/268 and
https://github.com/satijalab/seurat-object/releases/tag/v5.4.0 This
function is a fallback for older Seurat versions that do not have this
setter function. It will simply save the default dimension reduction
name in the Seurat object's misc slot under `DefaultDimReduc`. When the
setter function is available, it will use that instead.
