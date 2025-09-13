#' Feature statistic plot when given a composed data frame
#'
#' @inheritParams FeatureStatPlot
#' @inheritDotParams FeatureStatPlot
#' @param data A data frame containing the feature data and metadata.
#' @return A ggplot object or a list if `combine` is FALSE
#' @keywords internal
#' @importFrom rlang %||% syms sym
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr summarise %>%
.feature_stat_plot <- function(
    data, features, plot_type, should_shrink, should_pivot, downsample,
    graph = NULL, bg_cutoff = 0, dims = 1:2, rows_name = "Features",
    ident = NULL, agg = mean, group_by = NULL,
    split_by = NULL, facet_by = NULL, xlab = NULL, ylab = NULL, x_text_angle = NULL,
    ...
) {
    unlisted_features <- unname(unlist(features))
    if (should_shrink) {
        dims <- if (is.null(dims)) NULL else colnames(data)[dims]
        selected_columns <- c(dims, ident, unlisted_features, group_by, if (!isTRUE(split_by)) split_by)
        nonexisting_columns <- setdiff(selected_columns, colnames(data))
        if (length(nonexisting_columns) > 0) {
            stop("[FeatureStatPlot] The following columns are not found in the object: ", paste(nonexisting_columns, collapse = ", "))
        }
        data <- data[, selected_columns, drop = FALSE]
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
    downsample_data <- function() {
        if (!is.null(downsample) && !is.null(ident)) {
            if (downsample > 0 && downsample <= 1) {
                slice_sample(data, prop = downsample, by = !!sym(ident))
            } else {  # downsample > 1
                slice_sample(data, n = downsample, by = !!sym(ident))
            }
        } else {
            data
        }
    }

    if (plot_type == "violin") {
        data <- downsample_data()
        ViolinPlot(
            data, x = ident, y = ".value", group_by = group_by, split_by = split_by, facet_by = facet_by,
            xlab = xlab %||% "", ylab = ylab %||% "", x_text_angle = x_text_angle %||% 45, ...)
    } else if (plot_type == "box") {
        data <- downsample_data()
        BoxPlot(
            data, x = ident, y = ".value", group_by = group_by, split_by = split_by, facet_by = facet_by,
            xlab = xlab %||% "", ylab = ylab %||% "", x_text_angle = x_text_angle %||% 45, ...)
    } else if (plot_type == "bar") {
        data <- data %>% dplyr::group_by(!!!syms(unique(c(ident, group_by, split_by, ".features")))) %>%
            summarise(.value = agg(!!sym(".value")))
        BarPlot(
            data, x = ident, y = ".value", group_by = group_by, split_by = split_by, facet_by = facet_by,
            xlab = xlab %||% "", ylab = ylab %||% "", x_text_angle = x_text_angle %||% 45, ...)
    } else if (plot_type == "ridge") {
        data <- downsample_data()
        RidgePlot(
            data, x = ".value", group_by = ident, split_by = split_by, facet_by = facet_by,
            xlab = xlab %||% "", ylab = ylab %||% "", ...)
    } else if (plot_type == "dim") {
        FeatureDimPlot(
            data, dims = dims, features = unlisted_features, graph = graph, bg_cutoff = bg_cutoff,
            split_by = split_by, facet_by = facet_by, xlab = xlab, ylab = ylab, ...)
    } else if (plot_type == "heatmap") {
        args <- rlang::dots_list(...)
        args$data <- data
        args$rows_by <- unlisted_features
        args$columns_by <- ident
        args$rows_name <- rows_name
        args$split_by <- split_by
        args$values_fill <- args$values_fill %||% 0
        args$show_row_names <- args$show_row_names %||% TRUE
        args$show_column_names <- args$show_column_names %||% TRUE
        if (!is.null(names(features))) {
            # handle named features.
            # names will be used as rows_split_by
            if (!is.null(args$rows_split_by)) {
                stop("[FeatureStatPlot] Cannot use both `rows_split_by` and named features.")
            }
            feature_groups <- lapply(names(features), function(x) rep(x, length(features[[x]])))
            feature_groups <- do.call(c, feature_groups)
            fdata <- data.frame(
                Feautures = unlisted_features,
                FeatureGroups = feature_groups
            )
            args$rows_split_name <- args$rows_split_name %||% "Feature Groups"
            names(fdata) <- c(rows_name, args$rows_split_name)
            args$rows_split_by <- args$rows_split_name
            if (!is.null(args$rows_data)) {
                args$rows_data <- merge(fdata, args$rows_data, by = rows_name, all.x = TRUE)
            } else {
                args$rows_data <- fdata
            }
        }
        do.call(Heatmap, args)
    } else if (plot_type == "dot") {
        args <- rlang::dots_list(...)
        args$data <- data
        args$rows_by <- unlisted_features
        args$columns_by <- ident
        args$rows_name <- rows_name
        args$split_by <- split_by
        args$cell_type <- "dot"
        args$name <- args$name %||% "Expression Level"
        args$dot_size <- args$dot_size %||% function(x) sum(x > 0) / length(x)
        args$dot_size_name <- args$dot_size_name %||% "Percent Expressed"
        args$row_name_annotation <- args$row_name_annotation %||% TRUE
        args$column_name_annotation <- args$column_name_annotation %||% TRUE
        args$show_row_names <- args$show_row_names %||% TRUE
        args$show_column_names <- args$show_column_names %||% TRUE
        args$add_reticle <- args$add_reticle %||% TRUE
        args$cluster_rows <- args$cluster_rows %||% FALSE
        args$cluster_columns <- args$cluster_columns %||% FALSE
        args$row_names_side <- args$row_names_side %||% "left"
        if (!is.null(names(features))) {
            # handle named features.
            # names will be used as rows_split_by
            if (!is.null(args$rows_split_by)) {
                stop("[FeatureStatPlot] Cannot use both `rows_split_by` and named features.")
            }
            fdata <- data.frame(
                Feautures = unlisted_features,
                FeatureGroups = unlist(sapply(names(features), function(x) rep(x, length(features[[x]]))))
            )
            args$rows_split_name <- args$rows_split_name %||% "Feature Groups"
            names(fdata) <- c(rows_name, args$rows_split_name)
            args$rows_split_by <- args$rows_split_name
            if (!is.null(args$rows_data)) {
                args$rows_data <- merge(fdata, args$rows_data, by = rows_name, all.x = TRUE)
            } else {
                args$rows_data <- fdata
            }
        }
        do.call(Heatmap, args)
    } else if (plot_type == "cor") {
        if (length(unlisted_features) < 2) {
            stop("The number of features should be at least 2 for correlation plot.")
        }
        if (length(unlisted_features) == 2) {
            CorPlot(data, x = unlisted_features[1], y = unlisted_features[2], group_by = ident, split_by = split_by,
                facet_by = facet_by, xlab = xlab, ylab = ylab, ...)
        } else {
            CorPairsPlot(
                data, columns = unlisted_features, group_by = ident, split_by = split_by, facet_by = facet_by, ...)
        }
    }
}

#' Feature statistic plot
#'
#' @description This function creates various types of feature statistic plots for a Seurat object, a Giotto object,
#' a path to an .h5ad file or an opened `H5File` by `hdf5r` package.
#' It allows for plotting features such as gene expression, scores, or other metadata across different groups or conditions.
#' The function supports multiple plot types including violin, box, bar, ridge, dimension reduction, correlation, heatmap, and dot plots.
#' It can also handle multiple features and supports faceting, splitting, and grouping by metadata columns.
#' @param object A seurat object, a giotto object, a path to an .h5ad file or an opened `H5File` by `hdf5r` package.
#' @param features A character vector of feature names
#' @param plot_type Type of the plot. It can be "violin", "box", "bar", "ridge", "dim", "cor", "heatmap" or "dot"
#' @param spat_unit The spatial unit to use for the plot. Only applied to Giotto objects.
#' @param feat_type feature type of the features (e.g. "rna", "dna", "protein"), only applied to Giotto objects.
#' @param reduction Name of the reduction to plot (for example, "umap"), only used when `plot_type` is "dim" or you can to use the reduction as feature.
#' @param dims Dimensions to plot, only used when `plot_type` is "dim".
#' @param rows_name The name of the rows in the heatmap, only used when `plot_type` is "heatmap".
#' @param graph Specify the graph name to add edges between cell neighbors to the plot, only used when `plot_type` is "dim".
#' @param bg_cutoff Background cutoff for the dim plot, only used when `plot_type` is "dim".
#' @param ident The column name in the meta data to identify the cells.
#' @param assay The assay to use for the feature data.
#' @param layer The layer to use for the feature data.
#' @param agg The aggregation function to use for the bar plot.
#' @param downsample A numeric the number of cells in each identity group to downsample to for violin, box, or ridge plots.
#' If n > 1, it is treated as the number of cells to downsample to.
#' If 0 < n <= 1, it is treated as the fraction of cells to downsample to.
#' @param group_by The column name in the meta data to group the cells.
#' @param split_by Column name in the meta data to split the cells to different plots.
#'   If TRUE, the cells are split by the features.
#' @param facet_by Column name in the meta data to facet the plots. Should be always NULL.
#' @param xlab The x-axis label.
#' @param ylab The y-axis label.
#' @param x_text_angle The angle of the x-axis text. Only used when `plot_type` is "violin", "bar", or "box".
#' @param ... Other arguments passed to the plot functions.
#'  * For `plot_type` "violin", the arguments are passed to [plotthis::ViolinPlot()].
#'  * For `plot_type` "box", the arguments are passed to [plotthis::BoxPlot()].
#'  * For `plot_type` "bar", the arguments are passed to [plotthis::BarPlot()].
#'  * For `plot_type` "ridge", the arguments are passed to [plotthis::RidgePlot()].
#'  * For `plot_type` "dim", the arguments are passed to [plotthis::FeatureDimPlot()].
#'  * For `plot_type` "heatmap", the arguments are passed to [plotthis::Heatmap()].
#'  * For `plot_type` "cor" with 2 features, the arguments are passed to [plotthis::CorPlot()].
#'  * For `plot_type` "cor" with more than 2 features, the arguments are passed to [plotthis::CorPairsPlot()].
#'  * For `plot_type` "dot", the arguments are passed to [plotthis::Heatmap()] with `cell_type` set to "dot".
#' @return A ggplot object or a list if `combine` is FALSE
#' @export
#' @importFrom rlang %||%
#' @importFrom SeuratObject GetAssayData Embeddings DefaultDimReduc Graphs Reductions Idents
#' @importFrom plotthis ViolinPlot BoxPlot BarPlot DotPlot RidgePlot FeatureDimPlot Heatmap CorPlot CorPairsPlot
#' @details
#' See:
#' * <https://pwwang.github.io/scplotter/articles/Giotto_Visium.html>
#' * <https://pwwang.github.io/scplotter/articles/Giotto_Xenium.html>
#'
#' for examples of using this function with Giotto objects.
#'
#' And see:
#' * <https://pwwang.github.io/scplotter/articles/Working_with_anndata_h5ad_files.html>
#'
#' for examples of using this function with .h5ad files.
#' @examples
#' \donttest{
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
#' FeatureStatPlot(
#'    subset(pancreas_sub,
#'        subset = SubCellType %in% c("Ductal", "Ngn3 low EP", "Ngn3 high EP")),
#'    features = c("G2M_score"),
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
#' FeatureStatPlot(pancreas_sub, features = c("G2M_score", "nCount_RNA"),
#'    ident = "SubCellType", plot_type = "dim", facet_by = "Phase", split_by = TRUE, ncol = 1)
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
#'    Features = features,  # 'rows_name' default is "Features"
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
#'    plot_type = "heatmap", name = "Expression Level", dot_size = function(x) sum(x > 0) / length(x),
#'    dot_size_name = "Percent Expressed", add_bg = TRUE, rows_data = rows_data,
#'    show_row_names = TRUE, rows_split_by = "group", cluster_rows = FALSE,
#'    column_annotation = c("Phase", "G2M_score"),
#'    column_annotation_type = list(Phase = "pie", G2M_score = "violin"),
#'    column_annotation_params = list(G2M_score = list(show_legend = FALSE)),
#'    row_annotation = c("TF", "CSPA"),
#'    row_annotation_side = "right",
#'    row_annotation_type = list(TF = "simple", CSPA = "simple"))
#' FeatureStatPlot(pancreas_sub, features = features, ident = "SubCellType", cell_type = "dot",
#'    plot_type = "heatmap", name = "Expression Level", dot_size = function(x) sum(x > 0) / length(x),
#'    dot_size_name = "Percent Expressed", add_bg = TRUE,
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
#'    plot_type = "heatmap", dot_size = function(x) sum(x > 0) / length(x),
#'    dot_size_name = "Percent Expressed", palette = "viridis", add_reticle = TRUE,
#'    rows_data = rows_data, name = "Expression Level", show_row_names = TRUE,
#'    rows_split_by = "group")
#' # Visualize the markers for each sub-cell type (the markers can overlap)
#' # Say: markers <- Seurat::FindAllMarkers(pancreas_sub, ident = "SubCellType")
#' markers <- data.frame(
#'     avg_log2FC = c(
#'          3.44, 2.93, 2.72, 2.63, 2.13, 1.97, 2.96, 1.92, 5.22, 3.91, 3.64, 4.52,
#'          3.45, 2.45, 1.75, 2.08, 9.10, 4.45, 3.61, 6.30, 4.96, 3.49, 3.91, 3.90,
#'          10.58, 5.84, 4.73, 3.34, 7.22, 4.52, 10.10, 4.25),
#'     cluster = factor(rep(
#'          c("Ductal", "Ngn3 low EP", "Ngn3 high EP", "Pre-endocrine", "Beta",
#'            "Alpha", "Delta", "Epsilon"), each = 4),
#'          levels = levels(pancreas_sub$SubCellType)),
#'     gene = c(
#'          "Cyr61", "Adamts1", "Anxa2", "Bicc1", "1700011H14Rik", "Gsta3", "8430408G22Rik",
#'          "Anxa2", "Ppp1r14a", "Btbd17", "Neurog3", "Gadd45a", "Fev", "Runx1t1", "Hmgn3",
#'          "Cryba2", "Ins2", "Ppp1r1a", "Gng12", "Sytl4", "Irx1", "Tmem27", "Peg10", "Irx2",
#'          "Sst", "Ptprz1", "Arg1", "Frzb", "Irs4", "Mboat4", "Ghrl", "Arg1"
#'     )
#' )
#' FeatureStatPlot(pancreas_sub,
#'   features = unique(markers$gene), ident = "SubCellType", cell_type = "bars",
#'   plot_type = "heatmap", rows_data = markers, rows_name = "gene", rows_split_by = "cluster",
#'   show_row_names = TRUE, show_column_names = TRUE, name = "Expression Level",
#'   cluster_rows = FALSE, cluster_columns = FALSE, rows_split_palette = "Paired")
#'
#' # Use plot_type = "dot" to as a shortcut for heatmap with cell_type = "dot"
#' FeatureStatPlot(pancreas_sub, features = features, ident = "SubCellType", plot_type = "dot")
#'
#' named_features <- list(
#'    Ductal = c("Sox9", "Anxa2", "Bicc1"),
#'    EPs = c("Neurog3", "Hes6"),
#'    `Pre-endocrine` = c("Fev", "Neurod1"),
#'    Endocrine = c("Rbp4", "Pyy"),
#'    Beta = "Ins1", Alpha = "Gcg", Delta = "Sst", Epsilon = "Ghrl"
#' )
#' FeatureStatPlot(pancreas_sub, features = named_features, ident = "SubCellType",
#'    plot_type = "heatmap", name = "Expression Level", show_row_names = TRUE)
#'
#' # Correlation plot
#' FeatureStatPlot(pancreas_sub, features = c("Pyy", "Rbp4"), plot_type = "cor",
#'    anno_items = c("eq", "r2", "spearman"))
#' FeatureStatPlot(pancreas_sub, features = c("Ins1", "Gcg", "Sst", "Ghrl"),
#'    plot_type = "cor")
#' }
FeatureStatPlot <- function(
    object, features, plot_type = c("violin", "box", "bar", "ridge", "dim", "cor", "heatmap", "dot"),
    spat_unit = NULL, feat_type = NULL, downsample = NULL,
    reduction = NULL, graph = NULL, bg_cutoff = 0, dims = 1:2, rows_name = "Features",
    ident = NULL, assay = NULL, layer = NULL, agg = mean, group_by = NULL,
    split_by = NULL, facet_by = NULL, xlab = NULL, ylab = NULL, x_text_angle = NULL, ...
) {
    UseMethod("FeatureStatPlot")
}

#' @export
FeatureStatPlot.giotto <- function(
    object, features, plot_type = c("violin", "box", "bar", "ridge", "dim", "cor", "heatmap", "dot"),
    spat_unit = NULL, feat_type = NULL, downsample = NULL,
    reduction = NULL, graph = NULL, bg_cutoff = 0, dims = 1:2, rows_name = "Features",
    ident = NULL, assay = NULL, layer = NULL, agg = mean, group_by = NULL,
    split_by = NULL, facet_by = NULL, xlab = NULL, ylab = NULL, x_text_angle = NULL, ...
) {
    plot_type <- match.arg(plot_type)
    if (!is.null(facet_by) && plot_type != "dim") {
        stop("Cannot facet plots because the plots are facetted by the 'features'.")
    }
    stopifnot("[FeatureStatPlot] 'assay' should not be used for Giotto object." = is.null(assay))
    if (plot_type != "dim" && is.null(ident)) {
        stop("[FeatureStatPlot] 'ident' should be provided for non-dim plot for giotto object.")
    }
    # dim plot may use expression for highlighting cells
    # Heatmap may use other variables as annotations, but shrinking only includes minimal columns
    should_shrink <- !plot_type %in% c("dim", "heatmap", "dot")
    should_pivot <- !plot_type %in% c("dim", "heatmap", "dot", "cor")

    spat_unit <- GiottoClass::set_default_spat_unit(
        gobject = object,
        spat_unit = spat_unit
    )

    feat_type <- GiottoClass::set_default_feat_type(
        gobject = object,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    unlisted_features <- unname(unlist(features))
    assay_data <- GiottoClass::getExpression(
        object, values = layer,
        spat_unit = spat_unit, feat_type = feat_type,
        output = "matrix"
    )

    assay_feature <- intersect(unlisted_features, rownames(assay_data))
    assay_data <- t(as.matrix(assay_data[assay_feature, , drop = FALSE]))

    if (is.null(reduction) && plot_type != "dim") {
        data <- cbind(
            GiottoClass::getCellMetadata(
                gobject = object,
                spat_unit = spat_unit,
                feat_type = feat_type,
                output = "data.table"
            ),
            assay_data
        )
    } else {
        data <- cbind(
            GiottoClass::getDimReduction(
                gobject = object,
                reduction = "cells",
                reduction_method = reduction,
                spat_unit = spat_unit,
                feat_type = feat_type,
                output = "matrix"
            ),
            GiottoClass::getCellMetadata(
                gobject = object,
                spat_unit = spat_unit,
                feat_type = feat_type,
                output = "data.table"
            ),
            assay_data
        )
    }

    if (plot_type == "dim" && !is.null(graph)) {
        # spat_unit	feat_type	nn_type	name
        # cell	rna	sNN	sNN.pca
        graph_info <- GiottoClass::list_nearest_networks(
            gobject = object,
            spat_unit = spat_unit,
            feat_type = feat_type
        )
        if (is.null(graph_info)) {
            stop(
                "[FeatureStatPlot] The object does not have any nearest networks with given ",
                "spat_unit (", spat_unit, ") and feat_type (", feat_type, ")."
            )
        }
        if (isTRUE(graph)) { graph <- graph_info$name[1] }
        graph_info <- graph_info[, graph_info$name == graph, drop = FALSE]
        if (nrow(graph_info) == 0) {
            stop(
                "[FeatureStatPlot] The object does not have network: ", graph,
                " with given spat_unit (", spat_unit, ") and feat_type (", feat_type, ")."
            )
        }
        graph_info <- graph_info[1, , drop = FALSE]
        graph <- igraph::as_adjacency_matrix(GiottoClass::getNearestNetwork(
            object,
            spat_unit = spat_unit,
            feat_type = feat_type,
            nn_type = graph_info$nn_type,
            name = graph_info$name,
            output = "igraph"
        ))
    }

    .feature_stat_plot(
        data = data, features = features, plot_type = plot_type,
        should_shrink = should_shrink, should_pivot = should_pivot,
        graph = graph, bg_cutoff = bg_cutoff, downsample = downsample,
        dims = dims, rows_name = rows_name, ident = ident, agg = agg,
        group_by = group_by, split_by = split_by, facet_by = facet_by,
        xlab = xlab, ylab = ylab, x_text_angle = x_text_angle, ...
    )
}

#' @export
FeatureStatPlot.Seurat <- function(
    object, features, plot_type = c("violin", "box", "bar", "ridge", "dim", "cor", "heatmap", "dot"),
    spat_unit = NULL, feat_type = NULL, downsample = NULL,
    reduction = NULL, graph = NULL, bg_cutoff = 0, dims = 1:2, rows_name = "Features",
    ident = NULL, assay = NULL, layer = NULL, agg = mean, group_by = NULL,
    split_by = NULL, facet_by = NULL, xlab = NULL, ylab = NULL, x_text_angle = NULL, ...
) {
    plot_type <- match.arg(plot_type)
    if (!is.null(facet_by) && plot_type != "dim") {
        stop("Cannot facet plots because the plots are facetted by the 'features'.")
    }
    stopifnot("[FeatureStatPlot] 'spat_unit' and 'feat_type' should not be used for Seurat object." = is.null(spat_unit) && is.null(feat_type))

    reduction <- reduction %||% (
        if(is.null(Reductions(object))) NULL else DefaultDimReduc(object)
    )
    # dim plot may use expression for highlighting cells
    # Heatmap may use other variables as annotations, but shrinking only includes minimal columns
    should_shrink <- !plot_type %in% c("dim", "heatmap", "dot")
    should_pivot <- !plot_type %in% c("dim", "heatmap", "dot", "cor")

    unlisted_features <- unname(unlist(features))
    assay_data <- GetAssayData(object, assay = assay, layer = layer)
    assay_feature <- intersect(unlisted_features, rownames(assay_data))
    assay_data <- t(as.matrix(assay_data[assay_feature, , drop = FALSE]))
    if (is.null(reduction)) {
        data <- cbind(object@meta.data, assay_data)
    } else {
        emb <- Embeddings(object, reduction = reduction)
        data <- cbind(
            emb,
            object@meta.data[rownames(emb), , drop = FALSE],
            assay_data[rownames(emb), , drop = FALSE]
        )
    }

    if (is.null(ident)) {
        ident <- "Identity"
        data$Identity <- Idents(object)
    }

    if (plot_type == "dim" && !is.null(graph)) {
        if (!graph %in% Graphs(object)) {
            stop("The object does not have graph:", graph)
        }
        graph <- object@graphs[[graph]]
    }

    .feature_stat_plot(
        data = data, features = features, plot_type = plot_type,
        should_shrink = should_shrink, should_pivot = should_pivot,
        graph = graph, bg_cutoff = bg_cutoff, downsample = downsample,
        dims = dims, rows_name = rows_name, ident = ident, agg = agg,
        group_by = group_by, split_by = split_by, facet_by = facet_by,
        xlab = xlab, ylab = ylab, x_text_angle = x_text_angle, ...
    )
}

#' @export
FeatureStatPlot.character <- function(
    object, features, plot_type = c("violin", "box", "bar", "ridge", "dim", "cor", "heatmap", "dot"),
    spat_unit = NULL, feat_type = NULL, downsample = NULL,
    reduction = NULL, graph = NULL, bg_cutoff = 0, dims = 1:2, rows_name = "Features",
    ident = NULL, assay = NULL, layer = NULL, agg = mean, group_by = NULL,
    split_by = NULL, facet_by = NULL, xlab = NULL, ylab = NULL, x_text_angle = NULL, ...
) {
    if (!endsWith(object, ".h5ad")) {
        stop("[FeatureStatPlot] Currently only supports .h5ad files when called with a string/path.")
    }

    object <- hdf5r::H5File$new(object, mode = "r")
    on.exit(object$close_all())

    FeatureStatPlot.H5File(
        object, features = features, plot_type = plot_type,
        spat_unit = spat_unit, feat_type = feat_type, downsample = downsample,
        reduction = reduction, graph = graph, bg_cutoff = bg_cutoff,
        dims = dims, rows_name = rows_name, ident = ident,
        assay = assay, layer = layer, agg = agg,
        group_by = group_by, split_by = split_by, facet_by = facet_by,
        xlab = xlab, ylab = ylab, x_text_angle = x_text_angle, ...
    )
}

#' @export
FeatureStatPlot.H5File <- function(
    object, features, plot_type = c("violin", "box", "bar", "ridge", "dim", "cor", "heatmap", "dot"),
    spat_unit = NULL, feat_type = NULL, downsample = NULL,
    reduction = NULL, graph = NULL, bg_cutoff = 0, dims = 1:2, rows_name = "Features",
    ident = NULL, assay = NULL, layer = NULL, agg = mean, group_by = NULL,
    split_by = NULL, facet_by = NULL, xlab = NULL, ylab = NULL, x_text_angle = NULL, ...
) {
    plot_type <- match.arg(plot_type)
    if (!is.null(facet_by) && plot_type != "dim") {
        stop("Cannot facet plots because the plots are facetted by the 'features'.")
    }
    stopifnot("[FeatureStatPlot] 'spat_unit' and 'feat_type' should not be used for H5File object." = is.null(spat_unit) && is.null(feat_type))

    # dim plot may use expression for highlighting cells
    # Heatmap may use other variables as annotations, but shrinking only includes minimal columns
    should_shrink <- !plot_type %in% c("dim", "heatmap", "dot")
    should_pivot <- !plot_type %in% c("dim", "heatmap", "dot", "cor")
    if (plot_type != "dim" && is.null(ident)) {
        stop("[FeatureStatPlot] 'ident' should be provided for non-dim plot for H5File object.")
    }

    unlisted_features <- unname(unlist(features))
    assay <- assay %||% "RNA"
    assay_key <- ifelse(identical(assay, "RNA"), "X", paste0("layers/", assay))
    assay_data <- h5group_to_matrix(object[[assay_key]])

    feature_key <- hdf5r::h5attr(object[["var"]], "_index")  # index
    assay_features <- object[["var"]][[feature_key]]$read()
    rownames(assay_data) <- assay_features
    assay_features <- intersect(unlisted_features, assay_features)
    assay_data <- t(as.matrix(assay_data[assay_features, , drop = FALSE]))

    if (is.null(reduction) && plot_type != "dim") {
        data <- cbind(h5group_to_dataframe(object[["obs"]]), assay_data)
    } else {
        reductions <- names(object[['obsm']])
        if (is.null(reductions) || length(reductions) == 0) {
            stop("[FeatureStatPlot] The object does not have any reductions.")
        }

        reduction <- reduction %||% reductions[1]
        if (!startsWith(reduction, "X_") && !reduction %in% reductions) {
            reduction <- paste0("X_", reduction)
        }
        if (!reduction %in% reductions) {
            stop("[FeatureStatPlot] The object does not have reduction:", reduction)
        }

        data <- cbind(
            t(object[["obsm"]][[reduction]]$read()),
            h5group_to_dataframe(object[["obs"]]),
            assay_data
        )
    }

    if (plot_type == "dim" && !is.null(graph)) {
        ggrp <- object[['obsp']][['connectivities']]
        if (is.null(ggrp)) {
            stop("[FeatureStatPlot] The object does not have any graph/connectivities.")
        }

        graph <- h5group_to_matrix(ggrp)
        graph <- igraph::graph_from_adjacency_matrix(graph, mode = "undirected", weighted = TRUE)
        graph <- igraph::as_adjacency_matrix(graph, attr = "weight", sparse = TRUE)
        colnames(graph) <- rownames(graph) <- object[['obs']][['index']]$read()
        rownames(data) <- rownames(graph)
    }

    .feature_stat_plot(
        data = data, features = features, plot_type = plot_type,
        should_shrink = should_shrink, should_pivot = should_pivot,
        graph = graph, bg_cutoff = bg_cutoff, downsample = downsample,
        dims = dims, rows_name = rows_name, ident = ident, agg = agg,
        group_by = group_by, split_by = split_by, facet_by = facet_by,
        xlab = xlab, ylab = ylab, x_text_angle = x_text_angle, ...
    )
}
