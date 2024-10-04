#' FeatureStatPlot
#' @description Statistic plot for features
#' @param object A seurat object
#' @param features A character vector of feature names
#' @param plot_type Type of the plot. It can be "violin", "box", "bar", "ridge", "dim" or "heatmap".
#' @param reduction Name of the reduction to plot (for example, "umap"), only used when `plot_type` is "dim" or you can to use the reduction as feature.
#' @param dims Dimensions to plot, only used when `plot_type` is "dim".
#' @param rows_name The name of the rows in the heatmap, only used when `plot_type` is "heatmap".
#' @param graph Specify the graph name to add edges between cell neighbors to the plot, only used when `plot_type` is "dim".
#' @param bg_cutoff Background cutoff for the dim plot, only used when `plot_type` is "dim".
#' @param ident The column name in the meta data to identify the cells.
#' @param assay The assay to use for the feature data.
#' @param layer The layer to use for the feature data.
#' @param agg The aggregation function to use for the bar plot.
#' @param group_by The column name in the meta data to group the cells.
#' @param split_by Column name in the meta data to split the cells to different plots.
#'   If TRUE, the cells are split by the features.
#' @param facet_by Column name in the meta data to facet the plots. Should be always NULL.
#' @param xlab The x-axis label.
#' @param ylab The y-axis label.
#' @param x_text_angle The angle of the x-axis text. Only used when `plot_type` is "violin", "bar", or "box".
#' @param ... Other arguments passed to the plot functions.
#'  * For `plot_type` "violin", the arguments are passed to \code{\link{plotthis::ViolinPlot}}.
#'  * For `plot_type` "box", the arguments are passed to \code{\link{plotthis::BoxPlot}}.
#'  * For `plot_type` "bar", the arguments are passed to \code{\link{plotthis::BarPlot}}.
#'  * For `plot_type` "ridge", the arguments are passed to \code{\link{plotthis::RidgePlot}}.
#'  * For `plot_type` "dim", the arguments are passed to \code{\link{plotthis::FeatureDimPlot}}.
#'  * For `plot_type` "heatmap", the arguments are passed to \code{\link{plotthis::Heatmap}}.
#' @return A ggplot object or a list if `combine` is FALSE
#' @export
#' @importFrom rlang %||% syms
#' @importFrom dplyr group_by summarise
#' @importFrom tidyr pivot_longer
#' @importFrom SeuratObject GetAssayData Embeddings DefaultDimReduc Graphs
#' @importFrom plotthis ViolinPlot BoxPlot BarPlot DotPlot RidgePlot FeatureDimPlot Heatmap
#' @examples
#' data(pancreas_sub)
#'
#' FeatureStatPlot(pancreas_sub, features = c("G2M_score", "nCount_RNA"),
#'    ident = "SubCellType", facet_scales = "free_y")
#' FeatureStatPlot(pancreas_sub, features = c("G2M_score", "nCount_RNA"),
#'    ident = "SubCellType", plot_type = "box", facet_scales = "free_y")
#' FeatureStatPlot(pancreas_sub, features = c("G2M_score", "nCount_RNA"),
#'    ident = "SubCellType", plot_type = "bar", facet_scales = "free_y")
#' FeatureStatPlot(pancreas_sub, features = c("G2M_score", "nCount_RNA"),
#'    ident = "SubCellType", plot_type = "ridge", flip = TRUE, facet_scales = "free_y")
#' FeatureStatPlot(pancreas_sub, features = c("G2M_score", "nCount_RNA"),
#'    ident = "SubCellType", facet_scales = "free_y", add_point = TRUE)
#' FeatureStatPlot(pancreas_sub, features = c("G2M_score", "nCount_RNA"),
#'    ident = "SubCellType", facet_scales = "free_y", add_trend = TRUE)
#' FeatureStatPlot(pancreas_sub, features = c("G2M_score", "nCount_RNA"),
#'    ident = "SubCellType", facet_scales = "free_y", add_stat = mean)
#' FeatureStatPlot(pancreas_sub, features = c("G2M_score", "nCount_RNA"),
#'    ident = "SubCellType", facet_scales = "free_y", group_by = "Phase")
#' FeatureStatPlot(pancreas_sub, features = c("G2M_score"),
#'    ident = "SubCellType", group_by = "Phase", comparisons = TRUE)
#' FeatureStatPlot(pancreas_sub, features = c("Rbp4", "Pyy"), ident = "SubCellType",
#'    add_bg = TRUE, add_box = TRUE, stack = TRUE)
#' FeatureStatPlot(pancreas_sub, features = c(
#'        "Sox9", "Anxa2", "Bicc1", # Ductal
#'        "Neurog3", "Hes6", # EPs
#'        "Fev", "Neurod1", # Pre-endocrine
#'        "Rbp4", "Pyy", # Endocrine
#'        "Ins1", "Gcg", "Sst", "Ghrl" # Beta, Alpha, Delta, Epsilon
#'    ), ident = "SubCellType", add_bg = TRUE, stack = TRUE,
#'    legend.position = "top", legend.direction = "horizontal")
#' FeatureStatPlot(pancreas_sub, plot_type = "box", features = c(
#'       "Sox9", "Anxa2", "Bicc1", # Ductal
#'       "Neurog3", "Hes6", # EPs
#'       "Fev", "Neurod1", # Pre-endocrine
#'       "Rbp4", "Pyy", # Endocrine
#'       "Ins1", "Gcg", "Sst", "Ghrl" # Beta, Alpha, Delta, Epsilon
#'    ), ident = "SubCellType", add_bg = TRUE, stack = TRUE, flip = TRUE,
#'    legend.position = "top", legend.direction = "horizontal")
#' # Use splitting instead of facetting
#' FeatureStatPlot(pancreas_sub, features = c("Neurog3", "Rbp4", "Ins1"),
#'    ident = "CellType", split_by = TRUE)
#'
#' FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "G2M_score", reduction = "UMAP")
#' FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "G2M_score", reduction = "UMAP",
#'    bg_cutoff = -Inf)
#' FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "G2M_score", reduction = "UMAP",
#'    theme = "theme_blank")
#' FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "G2M_score", reduction = "UMAP",
#'    theme = ggplot2::theme_classic, theme_args = list(base_size = 16))
#'
#' # Label and highlight cell points
#' FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "Rbp4", reduction = "UMAP",
#'    highlight = 'SubCellType == "Delta"')
#' FeatureStatPlot(pancreas_sub, plot_type = "dim",
#'    features = "Rbp4", split_by = "Phase", reduction = "UMAP",
#'    highlight = TRUE, theme = "theme_blank")
#'
#' # Add a density layer
#' FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "Rbp4", reduction = "UMAP",
#'    add_density = TRUE)
#' FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "Rbp4", reduction = "UMAP",
#'    add_density = TRUE, density_filled = TRUE)
#'
#' # Change the plot type from point to the hexagonal bin
#' FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "Rbp4", reduction = "UMAP",
#'    hex = TRUE)
#' FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "Rbp4", reduction = "UMAP",
#'    hex = TRUE, hex_bins = 20)
#'
#' # Show lineages on the plot based on the pseudotime
#' FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "Lineage3", reduction = "UMAP",
#'    lineages = "Lineage3")
#' FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "Lineage3", reduction = "UMAP",
#'    lineages = "Lineage3", lineages_whiskers = TRUE)
#' FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "Lineage3", reduction = "UMAP",
#'    lineages = "Lineage3", lineages_span = 0.1)
#'
#' FeatureStatPlot(pancreas_sub, plot_type = "dim",
#'   features = c("Sox9", "Anxa2", "Bicc1"), reduction = "UMAP",
#'   theme = "theme_blank",
#'   theme_args = list(plot.subtitle = ggplot2::element_text(size = 10),
#'      strip.text = ggplot2::element_text(size = 8))
#' )
#'
#' # Plot multiple features with different scales
#' endocrine_markers <- c("Ins1", "Gcg", "Sst", "Ghrl")
#' FeatureStatPlot(pancreas_sub, endocrine_markers, reduction = "UMAP", plot_type = "dim")
#' FeatureStatPlot(pancreas_sub, endocrine_markers, reduction = "UMAP", lower_quantile = 0,
#'    upper_quantile = 0.8, plot_type = "dim")
#' FeatureStatPlot(pancreas_sub, endocrine_markers, reduction = "UMAP",
#'    lower_cutoff = 1, upper_cutoff = 4, plot_type = "dim")
#' FeatureStatPlot(pancreas_sub, endocrine_markers, reduction = "UMAP", bg_cutoff = 2,
#'    lower_cutoff = 2, upper_cutoff = 4, plot_type = "dim")
#' FeatureStatPlot(pancreas_sub, c("Sst", "Ghrl"), split_by = "Phase", reduction = "UMAP",
#'    plot_type = "dim")
#'
#' # Heatmap
#' features <- c(
#'    "Sox9", "Anxa2", "Bicc1", # Ductal
#'    "Neurog3", "Hes6", # EPs
#'    "Fev", "Neurod1", # Pre-endocrine
#'    "Rbp4", "Pyy", # Endocrine
#'    "Ins1", "Gcg", "Sst", "Ghrl" # Beta, Alpha, Delta, Epsilon
#' )
#' rows_data <- data.frame(
#'    features = features,
#'    group = c(
#'        "Ductal", "Ductal", "Ductal", "EPs", "EPs", "Pre-endocrine",
#'        "Pre-endocrine", "Endocrine", "Endocrine", "Beta", "Alpha", "Delta", "Epsilon"),
#'    TF = c(TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE,
#'        TRUE, TRUE, TRUE),
#'    CSPA = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE,
#'        FALSE, FALSE, FALSE)
#' )
#' FeatureStatPlot(pancreas_sub, features = features, ident = "SubCellType",
#'    plot_type = "heatmap", name = "Expression Level")
#' FeatureStatPlot(pancreas_sub, features = features, ident = "Phase",
#'    plot_type = "heatmap", name = "Expression Level", columns_split_by = "SubCellType")
#' FeatureStatPlot(pancreas_sub, features = features, ident = "SubCellType",
#'    plot_type = "heatmap", cell_type = "bars", name = "Expression Level")
#' FeatureStatPlot(pancreas_sub, features = features, ident = "SubCellType", cell_type = "dot",
#'    plot_type = "heatmap", name = "Expression Level", dot_size = mean, add_bg = TRUE,
#'    rows_data = rows_data, show_row_names = TRUE, rows_split_by = "group", cluster_rows = FALSE,
#'    column_annotation = c("Phase", "G2M_score"),
#'    column_annotation_type = list(Phase = "pie", G2M_score = "violin"),
#'    column_annotation_params = list(G2M_score = list(show_legend = FALSE)),
#'    row_annotation = c("TF", "CSPA"),
#'    row_annotation_side = "right",
#'    row_annotation_type = list(TF = "simple", CSPA = "simple"))
#' FeatureStatPlot(pancreas_sub, features = features, ident = "SubCellType", cell_type = "dot",
#'    plot_type = "heatmap", name = "Expression Level", dot_size = mean, add_bg = TRUE,
#'    rows_data = rows_data, show_column_names = TRUE, rows_split_by = "group",
#'    cluster_rows = FALSE, flip = TRUE, palette = "YlOrRd",
#'    column_annotation = c("Phase", "G2M_score"),
#'    column_annotation_type = list(Phase = "pie", G2M_score = "violin"),
#'    column_annotation_params = list(G2M_score = list(show_legend = FALSE)),
#'    row_annotation = c("TF", "CSPA"),
#'    row_annotation_side = "right",
#'    row_annotation_type = list(TF = "simple", CSPA = "simple"))
#' FeatureStatPlot(pancreas_sub, features = features, ident = "SubCellType", cell_type = "violin",
#'    plot_type = "heatmap", name = "Expression Level", show_row_names = TRUE,
#'    cluster_columns = FALSE, rows_split_by = "group", rows_data = rows_data)
#' FeatureStatPlot(pancreas_sub, features = features, ident = "SubCellType", cell_type = "dot",
#'    plot_type = "heatmap", dot_size = mean, palette = "viridis", add_reticle = TRUE,
#'    rows_data = rows_data, name = "Expression Level", show_row_names = TRUE,
#'    rows_split_by = "group")
#' named_features <- list(
#'    Ductal = c("Sox9", "Anxa2", "Bicc1"),
#'    EPs = c("Neurog3", "Hes6"),
#'    `Pre-endocrine` = c("Fev", "Neurod1"),
#'    Endocrine = c("Rbp4", "Pyy"),
#'    Beta = "Ins1", Alpha = "Gcg", Delta = "Sst", Epsilon = "Ghrl"
#' )
#' FeatureStatPlot(pancreas_sub, features = named_features, ident = "SubCellType",
#'    plot_type = "heatmap", name = "Expression Level", show_row_names = TRUE)
FeatureStatPlot <- function(
    object, features, plot_type = c("violin", "box", "bar", "ridge", "dim", "heatmap", "cor"),
    reduction = NULL, graph = NULL, bg_cutoff = 0, dims = 1:2, rows_name = "Features",
    ident = "seurat_clusters", assay = NULL, layer = NULL, agg = mean, group_by = NULL,
    split_by = NULL, facet_by = NULL, xlab = NULL, ylab = NULL, x_text_angle = NULL, ...
) {
    if (!is.null(facet_by)) {
        stop("Cannot facet plots because the plots are facetted by the 'features'.")
    }
    plot_type <- match.arg(plot_type)

    reduction <- reduction %||% DefaultDimReduc(object)
    # dim plot may use expression for highlighting cells
    # Heatmap may use other variables as annotations, but shrinking only includes minimal columns
    should_shrink <- !plot_type %in% c("dim", "heatmap")
    should_pivot <- !plot_type %in% c("dim", "heatmap", "cor")

    unlisted_features <- unname(unlist(features))
    assay_data <- GetAssayData(object, assay = assay, layer = layer)
    assay_feature <- intersect(unlisted_features, rownames(assay_data))
    assay_data <- t(as.matrix(assay_data[assay_feature, , drop = FALSE]))
    data <- cbind(Embeddings(object, reduction = reduction), object@meta.data, assay_data)

    if (should_shrink) {
        dims <- if (is.null(dims)) NULL else colnames(data)[dims]
        data <- data[, c(dims, ident, unlisted_features, group_by, if (isTRUE(split_by)) NULL else split_by), drop = FALSE]
    }
    if (should_pivot) {
        data <- pivot_longer(data, cols = unlisted_features, names_to = ".features", values_to = ".value")

    if (isTRUE(split_by)) {
        split_by <- ".features"
        facet_by <- NULL
    } else {
        facet_by <- ".features"
        }
    }

    if (plot_type == "violin") {
        ViolinPlot(
            data, x = ident, y = ".value", group_by = group_by, split_by = split_by, facet_by = facet_by,
            xlab = xlab %||% "", ylab = ylab %||% "", x_text_angle = x_text_angle %||% 45, ...)
    } else if (plot_type == "box") {
        BoxPlot(
            data, x = ident, y = ".value", group_by = group_by, split_by = split_by, facet_by = facet_by,
            xlab = xlab %||% "", ylab = ylab %||% "", x_text_angle = x_text_angle %||% 45, ...)
    } else if (plot_type == "bar") {
        data <- data %>% group_by(!!!syms(unique(c(ident, split_by, ".features")))) %>%
            summarise(.value = agg(.value))
        BarPlot(
            data, x = ident, y = ".value", group_by = group_by, split_by = split_by, facet_by = facet_by,
            xlab = xlab %||% "", ylab = ylab %||% "", x_text_angle = x_text_angle %||% 45, ...)
    } else if (plot_type == "ridge") {
        RidgePlot(
            data, x = ".value", group_by = ident, split_by = split_by, facet_by = facet_by,
            xlab = xlab %||% "", ylab = ylab %||% "", ...)
    } else if (plot_type == "dim") {
        if (!is.null(graph)) {
            if (!graph %in% Graphs(object)) {
                stop("The object does not have graph:", graph)
            }
            graph <- object@graphs[[graph]]
        }
        FeatureDimPlot(
            data, dims = dims, features = unlisted_features, graph = graph, bg_cutoff = bg_cutoff,
            split_by = split_by, xlab = xlab, ylab = ylab, ...)
    } else if (plot_type == "heatmap") {
        Heatmap(data, rows = features, columns_by = ident, rows_name = rows_name, split_by = split_by, ...)
    }
}
