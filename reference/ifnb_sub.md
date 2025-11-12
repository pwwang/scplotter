# A subsetted version of 'ifnb' datasets

Human PBMC control/IFNB-stimulated dataset

## Format

A `Seurat` object.

## Source

<https://www.nature.com/articles/nbt.4042>

## Details

After processing, this object is further slimmed by:

- removing the other layers than "data" of the RNA assay

- reducing the dimension of reductions (pca and umap) to 2

## Examples

``` r
if (FALSE) { # \dontrun{
if (interactive()) {
  if (!require("SeuratData", quietly = TRUE)) {
    devtools::install_github("satijalab/seurat-data")
  }
  library(SeuratData)
  library(Seurat)
  suppressWarnings(InstallData("ifnb"))
  data("ifnb")
  set.seed(11)
  cells_sub <- unlist(lapply(split(colnames(ifnb), ifnb$stim), function(x) sample(x, size = 1000)))
  ifnb_sub <- subset(ifnb, cells = cells_sub)
  ifnb_sub <- ifnb_sub[rowSums(ifnb_sub@assays$RNA@counts) > 0, ]
  ifnb_sub <- UpdateSeuratObject(ifnb_sub)
  ifnb_sub <- NormalizeData(ifnb_sub)
  ifnb_sub <- FindVariableFeatures(ifnb_sub)
  ifnb_sub <- ScaleData(ifnb_sub)
  ifnb_sub <- RunPCA(ifnb_sub)
  ifnb_sub <- FindNeighbors(ifnb_sub)
  ifnb_sub <- FindClusters(
      ifnb_sub,
      resolution = c(setdiff(seq(0.2, 1.2, 0.2), 0.4), 0.4)
  )
  ifnb_sub <- RunUMAP(ifnb_sub, dims = 1:30)
  # usethis::use_data(ifnb_sub, compress = "xz")
}
} # }
```
