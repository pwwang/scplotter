#' A subsetted version of 'ifnb' datasets
#'
#' Human PBMC control/IFNB-stimulated dataset
#'
#' @format A \code{Seurat} object.
#' @concept data
#' @source \url{https://www.nature.com/articles/nbt.4042}
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   if (!require("SeuratData", quietly = TRUE)) {
#'     devtools::install_github("satijalab/seurat-data")
#'   }
#'   library(SeuratData)
#'   library(Seurat)
#'   suppressWarnings(InstallData("ifnb"))
#'   data("ifnb")
#'   set.seed(11)
#'   cells_sub <- unlist(lapply(split(colnames(ifnb), ifnb$stim), function(x) sample(x, size = 1000)))
#'   ifnb_sub <- subset(ifnb, cells = cells_sub)
#'   ifnb_sub <- ifnb_sub[rowSums(ifnb_sub@assays$RNA@counts) > 0, ]
#'   ifnb_sub <- UpdateSeuratObject(ifnb_sub)
#'   ifnb_sub <- NormalizeData(ifnb_sub)
#'   ifnb_sub <- FindVariableFeatures(ifnb_sub)
#'   ifnb_sub <- ScaleData(ifnb_sub)
#'   ifnb_sub <- RunPCA(ifnb_sub)
#'   ifnb_sub <- FindNeighbors(ifnb_sub)
#'   ifnb_sub <- FindClusters(
#'       ifnb_sub,
#'       resolution = c(setdiff(seq(0.2, 1.2, 0.2), 0.4), 0.4)
#'   )
#'   ifnb_sub <- RunUMAP(ifnb_sub, dims = 1:30)
#'   # usethis::use_data(ifnb_sub, compress = "xz")
#' }
#' }
#' @name ifnb_sub
NULL

#' A subsetted version of mouse 'pancreas' datasets
#'
#' Mouse pancreatic endocrinogenesis dataset from \href{https://doi.org/10.1242/dev.173849}{Bastidas-Ponce et al. (2019)}. A total of 1000 cells were downsampled to form the \code{pancreas_sub} dataset.
#'
#' @format A \code{Seurat} object.
#' @concept data
#' @source \url{https://scvelo.readthedocs.io/scvelo.datasets.pancreas/} \url{https://github.com/theislab/scvelo_notebooks/raw/master/data/Pancreas/endocrinogenesis_day15.h5ad}
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   library(Seurat)
#'   library(reticulate)
#'   library(slingshot)
#'   check_Python("scvelo")
#'   scv <- import("scvelo")
#'   adata <- scv$datasets$pancreas()
#'   pancreas <- adata_to_srt(adata)
#'   set.seed(11)
#'   pancreas_sub <- subset(pancreas, cells = sample(colnames(pancreas), size = 1000))
#'   pancreas_sub <- pancreas_sub[rowSums(pancreas_sub@assays$RNA@counts) > 0, ]
#'   pancreas_sub[["CellType"]] <- pancreas_sub[["clusters_coarse"]]
#'   pancreas_sub[["SubCellType"]] <- pancreas_sub[["clusters"]]
#'   pancreas_sub[["clusters_coarse"]] <- pancreas_sub[["clusters"]] <- NULL
#'   pancreas_sub[["Phase"]] <- ifelse(pancreas_sub$S_score > pancreas_sub$G2M_score, "S", "G2M")
#'   pancreas_sub[["Phase"]][
#'       apply(pancreas_sub[[]][, c("S_score", "G2M_score")], 1, max) < 0, ] <- "G1"
#'   pancreas_sub[["Phase", drop = TRUE]] <- factor(pancreas_sub[["Phase", drop = TRUE]],
#'       levels = c("G1", "S", "G2M"))
#'   pancreas_sub[["PCA"]] <- pancreas_sub[["X_pca"]]
#'   pancreas_sub[["UMAP"]] <- pancreas_sub[["X_umap"]]
#'   pancreas_sub[["X_umap"]] <- pancreas_sub[["X_pca"]] <- NULL
#'   VariableFeatures(pancreas_sub) <- rownames(pancreas_sub[["RNA"]])[
#'        which(pancreas_sub[["RNA"]]@meta.features$highly_variable_genes == "True")]
#'   pancreas_sub <- NormalizeData(pancreas_sub)
#'   pancreas_sub <- FindVariableFeatures(pancreas_sub)
#'   pancreas_sub <- ScaleData(pancreas_sub)
#'   pancreas_sub <- RunPCA(pancreas_sub)
#'   pancreas_sub <- FindNeighbors(pancreas_sub)
#'   pancreas_sub <- FindClusters(pancreas_sub, resolution = 0.5)
#'   pancreas_sub <- RunUMAP(pancreas_sub, dims = 1:30)
#'   # run slingshot
#'   reduction <- DefaultDimReduc(pancreas_sub)
#'   sl <- slingshot(
#'       data = as.data.frame(pancreas_sub[[reduction]]@cell.embeddings[, 1:2]),
#'       clusterLabels = as.character(pancreas_sub$SubCellType)
#'   )
#'   df <- as.data.frame(slingPseudotime(sl))
#'   pancreas_sub <- AddMetaData(pancreas_sub, metadata = df)
#'   pancreas_sub <- AddMetaData(pancreas_sub, metadata = slingBranchID(sl), col.name = "BranchID")
#'   # usethis::use_data(pancreas_sub, compress = "xz")
#' }
#' }
#' @name pancreas_sub
NULL

#' A toy example of CellPhoneDB output from LIANA
#'
#' The dataset is generated using python package LIANA
#'
#' @format A \code{data.frame} with 10 rows and 6 columns
#' @concept data
#' @source \url{https://liana-py.readthedocs.io/en/latest/notebooks/basic_usage.html#Tileplot}
#' @examples
#' \dontrun{
#' # Python code
#' # # import liana
#' # import liana as li
#' # # needed for visualization and toy data
#' # import scanpy as sc
#' #
#' # from liana.method import cellphonedb
#' # adata = sc.datasets.pbmc68k_reduced()
#' # cellphonedb(adata,
#' #             groupby='bulk_labels',
#' #             # NOTE by default the resource uses HUMAN gene symbols
#' #             resource_name='consensus',
#' #             expr_prop=0.1,
#' #             verbose=True, key_added='cpdb_res')
#' # cellphonedb_res = adata.uns['cpdb_res']
#' # cellphonedb_res = cellphonedb_res[cellphonedb_res['cellphone_pvals'] < 0.05]
#' }
#' @name cellphonedb_res
NULL
