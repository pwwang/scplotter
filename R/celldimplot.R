#' CellDimPlot
#'
#' @description Dimension reduction plot
#' @param object A seurat object
#' @param reduction Name of the reduction to plot (for example, "umap").
#' @param graph Specify the graph name to add edges between cell neighbors to the plot.
#' @param group_by A character vector of column name(s) to group the data. Default is NULL.
#' @param ... Other arguments passed to [plotthis::DimPlot()].
#' @return A ggplot object or a list if `combine` is FALSE
#' @export
#' @importFrom SeuratObject DefaultDimReduc Embeddings Graphs Reductions Idents
#' @importFrom plotthis DimPlot
#' @examples
#' \donttest{
#' data(pancreas_sub)
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP")
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             theme = "theme_blank")
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             palette = "seurat", theme = "theme_blank")
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             theme = ggplot2::theme_classic, theme_args = list(base_size = 16))
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             raster = TRUE, raster_dpi = 30)
#'
#' # Highlight cells
#' CellDimPlot(pancreas_sub,
#'   group_by = "SubCellType", reduction = "UMAP",
#'   highlight = 'SubCellType == "Epsilon"'
#' )
#' CellDimPlot(pancreas_sub,
#'   group_by = "SubCellType", split_by = "Phase", reduction = "UMAP",
#'   highlight = TRUE, theme = "theme_blank", legend.position = "none"
#' )
#' CellDimPlot(pancreas_sub,
#'   group_by = "SubCellType", facet_by = "Phase", reduction = "UMAP",
#'   highlight = TRUE, theme = "theme_blank", legend.position = "none"
#' )
#'
#' # Add group labels
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             label = TRUE)
#' CellDimPlot(pancreas_sub,
#'   group_by = "SubCellType", reduction = "UMAP",
#'   label = TRUE, label_fg = "orange", label_bg = "red", label_size = 5
#' )
#' CellDimPlot(pancreas_sub,
#'   group_by = "SubCellType", reduction = "UMAP",
#'   label = TRUE, label_insitu = TRUE
#' )
#' CellDimPlot(pancreas_sub,
#'   group_by = "SubCellType", reduction = "UMAP",
#'   label = TRUE, label_insitu = TRUE, label_repel = TRUE,
#'   label_segment_color = "red"
#' )
#'
#' # Add various shape of marks
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             add_mark = TRUE)
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             add_mark = TRUE, mark_expand = grid::unit(1, "mm"))
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             add_mark = TRUE, mark_alpha = 0.3)
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             add_mark = TRUE, mark_linetype = 2)
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             add_mark = TRUE, mark_type = "ellipse")
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             add_mark = TRUE, mark_type = "rect")
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             add_mark = TRUE, mark_type = "circle")
#'
#' # Add a density layer
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             add_density = TRUE)
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             add_density = TRUE, density_filled = TRUE)
#' CellDimPlot(pancreas_sub,
#'   group_by = "SubCellType", reduction = "UMAP",
#'   add_density = TRUE, density_filled = TRUE, density_filled_palette = "Blues",
#'   highlight = TRUE
#' )
#'
#' # Add statistical charts
#' CellDimPlot(pancreas_sub,
#'   group_by = "CellType", reduction = "UMAP", stat_by = "Phase")
#' CellDimPlot(pancreas_sub,
#'   group_by = "CellType", reduction = "UMAP", stat_by = "Phase",
#'   stat_plot_type = "ring", stat_plot_label = TRUE, stat_plot_size = 0.15)
#' CellDimPlot(pancreas_sub,
#'   group_by = "CellType", reduction = "UMAP", stat_by = "Phase",
#'   stat_plot_type = "bar", stat_type = "count")
#' CellDimPlot(pancreas_sub,
#'   group_by = "CellType", reduction = "UMAP", stat_by = "Phase",
#'   stat_plot_type = "line", stat_type = "count", stat_args = list(point_size = 1))
#'
#' # Chane the plot type from point to the hexagonal bin
#' CellDimPlot(pancreas_sub, group_by = "CellType", reduction = "UMAP",
#'             hex = TRUE)
#' CellDimPlot(pancreas_sub, group_by = "CellType", reduction = "UMAP",
#'             hex = TRUE, hex_bins = 20)
#' CellDimPlot(pancreas_sub, group_by = "CellType", reduction = "UMAP",
#'             hex = TRUE, hex_count = FALSE)
#'
#' # Show neighbors graphs on the plot
#' CellDimPlot(pancreas_sub, group_by = "CellType", reduction = "UMAP",
#'             graph = "RNA_nn")
#' CellDimPlot(pancreas_sub, group_by = "CellType", reduction = "UMAP",
#'             graph = "RNA_snn", edge_color = "grey80")
#'
#' # Show lineages on the plot based on the pseudotime
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             lineages = paste0("Lineage", 1:3))
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             lineages = paste0("Lineage", 1:3), lineages_whiskers = TRUE)
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             lineages = paste0("Lineage", 1:3), lineages_span = 0.1)
#' }
CellDimPlot <- function(object, reduction = NULL, graph = NULL, group_by = NULL, ...) {
    reduction <- reduction %||% DefaultDimReduc(object)

    if (!reduction %in% Reductions(object)) {
        stop("The object does not have reduction:", reduction)
    }

    if (!is.null(graph)) {
        if (!graph %in% Graphs(object)) {
            stop("The object does not have graph:", graph)
        }
        graph <- object@graphs[[graph]]
    }

    data <- cbind(Embeddings(object, reduction = reduction), object@meta.data)
    if (is.null(group_by)) {
        group_by <- "Identity"
        data[[group_by]] <- Idents(object)
    }

    DimPlot(data, graph = graph, group_by = group_by, ...)
}
