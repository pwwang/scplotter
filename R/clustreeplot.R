
#' Clustree plot
#'
#' @description This function generates a clustree plot from a data frame or a Seurat object.
#' @param object The data frame or Seurat object
#' @param ... Other arguments passed to [plotthis::ClustreePlot()]
#' @return A ggplot object or a list if `combine` is FALSE
#' @export
#' @examples
#' data(ifnb_sub)
#' ClustreePlot(ifnb_sub, prefix = "RNA_snn_res.")
ClustreePlot <- function(object, ...) {
    plotthis::ClustreePlot(object@meta.data, ...)
}
