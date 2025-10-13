#' Visualize Markers
#'
#' @description
#' Plot markers, typically identified by [Seurat::FindMarkers()] or [Seurat::FindAllMarkers()].
#' @param markers A data frame of markers, typically identified by [Seurat::FindMarkers()] or [Seurat::FindAllMarkers()].
#' @param object A Seurat object. Required for some plot types, see `plot_type`.
#' @param plot_type Type of plot to generate. Options include:
#' * `volcano`/`volcano_log2fc`: Volcano plot with log2 fold change on x-axis and -log10(p-value) on y-axis.
#'   If `p_adjust` is TRUE, -log10(adjusted p-value) is used instead.
#' * `volcano_pct`: Volcano plot with difference in percentage of cells expressing the gene between two groups on x-axis and -log10(p-value) on y-axis.
#'   If `p_adjust` is TRUE, -log10(adjusted p-value) is used instead.
#' * `jitter`/`jitter_log2fc`: Jitter plot of log2 fold change for each gene.
#'   The x-axis is the groups defined by `subset_by`, and the y-axis is log2 fold change.
#'   The size of the dots represents -log10(p-value) or -log10(adjusted p-value) if `p_adjust` is TRUE.
#' * `jitter_pct`: Jitter plot of difference in percentage of cells expressing the gene between two groups for each gene.
#'   The x-axis is the groups defined by `subset_by`, and the y-axis is the difference in percentage of cells expressing the gene between two groups.
#'   The size of the dots represents -log10(p-value) or -log10(adjusted p-value) if `p_adjust` is TRUE.
#' * `heatmap_log2fc`: Heatmap of log2 fold change for each gene across groups defined by `subset_by`. By specifying `cutoff`,
#'   The heatmap cells will be labeled with "*" for p-value < cutoff, if `p_adjust` is TRUE, adjusted p-value < cutoff.
#' * `heatmap_pct`: Heatmap of difference in percentage of cells expressing the gene between two groups for each gene across groups defined by `subset_by`.
#'   By specifying `cutoff`, The heatmap cells will be labeled with "*" for p-value < cutoff, if `p_adjust` is TRUE, adjusted p-value < cutoff.
#' * `dot_log2fc`: Dot plot of log2 fold change for each gene across groups defined by `subset_by`.
#'   The size of the dots represents -log10(p-value) or -log10(adjusted p-value) if `p_adjust` is TRUE.
#' * `dot_pct`: Dot plot of difference in percentage of cells expressing the gene between two groups for each gene across groups defined by `subset_by`.
#'   The size of the dots represents -log10(p-value) or -log10(adjusted p-value) if `p_adjust` is TRUE.
#' * `heatmap`: Heatmap of expression values for each gene across groups defined by `subset_by`. Requires `object`.
#' * `violin`: Violin plot of expression values for each gene across groups defined by `subset_by`. Requires `object`.
#' * `box`: Box plot of expression values for each gene across groups defined by `subset_by`. Requires `object`.
#' * `bar`: Bar plot of average expression values for each gene across groups defined by `subset_by`. Requires `object`.
#' * `ridge`: Ridge plot of expression values for each gene across groups defined by `subset_by`. Requires `object`.
#' * `dot`: Dot plot of expression values for each gene across groups defined by `subset_by`. Requires `object`.
#' @param subset_by A column in markers indicating where the markers are identified from, e.g., cluster or condition.
#' If object is provided, you can provide the corresponding metadata column in `subset_by` to merge the markers with the object's metadata, using `:` as separator. For example, if the markers are identified from different clusters (identified by `FindAllMarkers`), and the object's ident column is "RNA_snn_res.0.8", you can set `subset_by = "cluster:RNA_snn_res.0.8"`.
#' Note that for other columns to be merged, only the first value of each group defined by `subset_by` will be used.
#' For some plots, this is used to split the markers into multiple plots, e.g., one plot for each cluster.
#' See `plot_type` for details.
#' @param subset_as_facet Logical, whether to facet the plots by `subset_by` if applicable.
#' @param comparison_by The metadata column in markers indicating the comparion.
#' When visualizing the expression values, this column should also be in the object's metadata.
#' The values of this column should either a single value, indicating the comparison is between this group and all other cells (in the subset),
#' Similar as `subset_by`, you can provide the corresponding object metadata column by using `:` as separator.
#' If `NULL`, all markers are treated as from one comparison.
#' @param p_adjust Logical, whether to use adjusted p-value for plots that involve p-values.
#' Default is TRUE.
#' @param cutoff Numeric, p-value or adjusted p-value cutoff to label significance in heatmap plots.
#' Default is NULL, no cutoff.
#' @param order_by A string of expression to order the markers within each group defined by `subset_by`.
#' In addition to the columns in `markers`, you can also use the columns from the object's metadata if `object` is provided.
#' The object's metadata will be merged with `markers` by `subset_by`. Be carefull that only the first value of other columns will be used.
#' @param select Number of top markers (ordered by `order_by`) to select for each group defined by `subset_by` or a string of expression to filter markers.
#' It will be evaluated by [dplyr::filter()].
#' For `volcano`, `volcano_log2fc`, `volcano_pct`, `jitter`, `jitter_log2fc`, and `jitter_pct` plots, the selected markers will be labeled in the plot.
#' FOr  other plot types, only the selected markers will be plotted.
#' Default is 5 for `volcano`, `volcano_log2fc`, `volcano_pct`, `jitter`, `jitter_log2fc`, `jitter_pct`, and 10 for other plot types.
#' @param ... Additional arguments passed to specific plotting functions.
#' See Details.
#' @return A ggplot object or a list of ggplot objects if not merged.
#' @details
#' Additional arguments passed to specific plotting functions:
#' * For `heatmap_log2fc`, `heatmap_pct`, `dot_log2fc`, and `dot_pct`, the arguments will be passed to [plotthis::Heatmap()]
#' * For `violin`, "box", "bar", "ridge" and "dot", the arguments will be passed to [scplotter::FeatureStatPlot()]
#' * For `volcano`, `volcano_log2fc`, `volcano_pct`, the arguments will be passed to [plotthis::VolcanoPlot()]
#' * For `jitter`, `jitter_log2fc`, `jitter_pct`, the arguments will be passed to [plotthis::JitterPlot()]
#' @examples
#' \donttest{
#' data(pancreas_sub)
#' markers <- Seurat::FindMarkers(pancreas_sub,
#'  group.by = "Phase", ident.1 = "G2M", ident.2 = "G1")
#' allmarkers <- Seurat::FindAllMarkers(pancreas_sub)  # seurat_clusters
#'
#' MarkersPlot(markers)
#' MarkersPlot(markers, x_cutoff = 2)
#' MarkersPlot(allmarkers, p_adjust = FALSE,
#'     subset_by = "cluster", ncol = 2, subset_as_facet = TRUE)
#' MarkersPlot(markers, plot_type = "volcano_pct", flip_negative = TRUE)
#'
#' MarkersPlot(allmarkers, plot_type = "jitter", subset_by = "cluster")
#' MarkersPlot(allmarkers, plot_type = "jitter_pct",
#'     subset_by = "cluster", add_hline = 0, shape = 16)
#'
#' MarkersPlot(allmarkers, plot_type = "heatmap_log2fc", subset_by = "cluster")
#' MarkersPlot(allmarkers, plot_type = "heatmap_pct", subset_by = "cluster",
#'     cutoff = 0.05)
#'
#' MarkersPlot(allmarkers, plot_type = "dot_log2fc", subset_by = "cluster",
#'     add_reticle = TRUE)
#'
#' MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "heatmap",
#'    columns_split_by = "CellType",
#'    comparison_by = "cluster:seurat_clusters")
#'
#' # Suppose we did a DE between g1 and g2 in each cluster
#' allmarkers$comparison <- "g1:g2"
#' MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "heatmap",
#'    comparison_by = "Phase", subset_by = "cluster:seurat_clusters")
#' MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "dot",
#'    comparison_by = "Phase", subset_by = "cluster:seurat_clusters")
#'
#' MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "violin", select = 3,
#'    comparison_by = "Phase", subset_by = "cluster:seurat_clusters")
#'
#' MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "box", select = 3,
#'    comparison_by = "Phase", subset_by = "cluster:seurat_clusters")
#'
#' MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "ridge", select = 2,
#'    comparison_by = "Phase", subset_by = "cluster:seurat_clusters",
#'    ncol = 2)
#' }
#' @export
MarkersPlot <- function(
    markers, object = NULL,
    plot_type = c(
        # Plots that don't need object
        "volcano", "volcano_log2fc", "volcano_pct",
        "jitter", "jitter_log2fc", "jitter_pct",
        "heatmap_log2fc", "heatmap_pct", "dot_log2fc", "dot_pct",
        # Plots that need object, basically expression values, but
        # markers are selected from the markers data frame
        "heatmap", "violin", "box", "bar", "ridge", "dot"
    ),
    subset_by = NULL,
    subset_as_facet = FALSE,
    comparison_by = NULL,
    p_adjust = TRUE,
    cutoff = NULL,
    order_by = NULL,
    select = ifelse(plot_type %in% c(
        "volcano", "volcano_log2fc", "volcano_pct",
        "jitter", "jitter_log2fc", "jitter_pct"
    ), 5, 10),
    ...
) {
    plot_type <- match.arg(plot_type)

    # check if object is provided for plot types that need it
    plot_types_need_object <- c(
        "heatmap", "violin", "box", "bar", "ridge", "dot"
    )
    if (plot_type %in% plot_types_need_object & is.null(object)) {
        stop("[MarkersPlot] `object` is required for plot_type '", plot_type, "'")
    }

    # Add the gene column if missing (likely from FindMarkers)
    if (!"gene" %in% colnames(markers)) {
        markers$gene <- rownames(markers)
    }
    # Check subset_by
    check_columns <- utils::getFromNamespace("check_columns", "plotthis")
    if (!is.null(subset_by)) {
        subset_by_1 <- strsplit(subset_by, ":")[[1]][1]
        subset_by_2 <- strsplit(subset_by, ":")[[1]][2]  # NA if not provided
        subset_by_2 <- if (!is.na(subset_by_2)) subset_by_2 else NULL
    } else {
        subset_by_1 <- NULL
        subset_by_2 <- NULL
    }
    subset_by_1 <- check_columns(markers, subset_by_1)
    if (!is.null(subset_by_1) && !is.null(subset_by_2) && !is.null(object)) {
        if (!subset_by_2 %in% colnames(object@meta.data)) {
            stop("[MarkersPlot] `subset_by` '", subset_by_2, "' is not found in the object's metadata.")
        }
        # check if subset_by values are consistent between markers and object
        sub_markers <- unique(markers[[subset_by_1]])
        sub_object <- unique(object@meta.data[[subset_by_2]])
        nonexisting_sub <- setdiff(sub_markers, sub_object)
        if (length(nonexisting_sub) > 0) {
            stop(
                "[MarkersPlot] The following values in `subset_by` '", subset_by_1,
                "' are not found in the object's metadata (", subset_by_2, "): ",
                paste(nonexisting_sub, collapse = ", "))
        }
        meta <- dplyr::summarise(object@meta.data, dplyr::across(dplyr::everything(), ~ .[1]), .by = !!rlang::sym(subset_by_2))
        markers <- dplyr::left_join(markers, meta, by = stats::setNames(subset_by_2, subset_by_1), suffix = c("", ".meta"))
    }

    # Check comparison_by
    if (!is.null(comparison_by)) {
        comparison_by_1 <- strsplit(comparison_by, ":")[[1]][1]
        comparison_by_2 <- strsplit(comparison_by, ":")[[1]][2]  # NA if not provided
        comparison_by_2 <- if (!is.na(comparison_by_2)) comparison_by_2 else NULL
    } else {
        comparison_by_1 <- NULL
        comparison_by_2 <- NULL
    }
    comparison_by_1 <- check_columns(markers, comparison_by_1)
    if (!is.null(comparison_by_2) && !is.null(object)) {
        if (!comparison_by_2 %in% colnames(object@meta.data)) {
            stop("[MarkersPlot] `comparison_by` '", comparison_by_2, "' is not found in the object's metadata.")
        }
        # check if comparison_by values are consistent between markers and object
        comp_markers <- unique(unlist(strsplit(unique(as.character(markers[[comparison_by_1]])), ":")))
        comp_object <- unique(object@meta.data[[comparison_by_2]])
        nonexisting_comp <- setdiff(comp_markers, comp_object)
        if (length(nonexisting_comp) > 0) {
            stop(
                "[MarkersPlot] The following values in `comparison_by` '", comparison_by_1,
                "' are not found in the object's metadata (", comparison_by_2, "): ",
                paste(nonexisting_comp, collapse = ", "))
        }
    }

    pcol <- ifelse(p_adjust, "p_val_adj", "p_val")

    # calculate pct.1 - pct.2 if needed
    if (plot_type %in% c("volcano", "volcano_log2fc", "volcano_pct", "jitter_pct", "heatmap_pct", "dot_pct")) {
        if (!all(c("pct.1", "pct.2") %in% colnames(markers))) {
            stop("[MarkersPlot] `markers` must contain 'pct.1' and 'pct.2' columns for plot_type '", plot_type, "'")
        }
        markers <- dplyr::mutate(markers, pct_diff = !!sym("pct.1") - !!sym("pct.2"))
    }

    # calcualte -log10(p-value) or -log10(adjusted p-value)
    markers <- dplyr::mutate(markers, neg_log10_p = -log10(!!rlang::sym(pcol)))

    # order markers by order_by
    if (!is.null(order_by)) {
        markers <- dplyr::arrange(markers, !!rlang::parse_expr(order_by))
    }

    if (plot_type %in% c("volcano", "volcano_log2fc", "volcano_pct")) {
        args <- list(
            markers,
            x = ifelse(plot_type == "volcano_pct", "pct_diff", "avg_log2FC"),
            y = pcol,
            y_cutoff = cutoff %||% 0.05,
            y_cutoff_name = paste0(pcol, " = ", cutoff %||% 0.05),
            label_by = "gene",
            ...
        )
        args$color_by <- args$color_by %||% ifelse(plot_type == "volcano_pct", "avg_log2FC", "pct_diff")
        if (!is.null(subset_by_1) && subset_as_facet) {
            args$facet_by <- subset_by_1
        } else if (!is.null(subset_by_1)) {
            args$split_by <- subset_by_1
        }
        do.call(plotthis::VolcanoPlot, args)
    } else if (plot_type %in% c("jitter", "jitter_log2fc", "jitter_pct")) {
        if (is.null(subset_by_1)) {
            stop("[MarkersPlot] `subset_by` is required for plot_type '", plot_type, "'. Consider using volcano plot if you don't have groups.")
        }
        if (!is.numeric(select)) {
            stop("[MarkersPlot] `select` must be numeric for plot_type '", plot_type, "', to label top N markers in each group.")
        }
        args <- list(
            markers,
            x = subset_by_1,
            y = ifelse(plot_type == "jitter_pct", "pct_diff", "avg_log2FC"),
            size_by = "neg_log10_p",
            size_name = paste0("-log10(", pcol, ")"),
            label_by = "gene",
            nlabel = select,
            ...
        )
        if (!is.null(args$order_by)) {
            args$order_by <- order_by
        }
        do.call(plotthis::JitterPlot, args)
    } else if (plot_type %in% c("heatmap_log2fc", "heatmap_pct", "dot_log2fc", "dot_pct")) {
        if (is.null(subset_by)) {
            stop("[MarkersPlot] `subset_by` is required for plot_type '", plot_type, "'")
        }
        y <- ifelse(endsWith(plot_type, "_pct"), "pct_diff", "avg_log2FC")
        y_max <- max(markers[[y]], na.rm = TRUE)
        y_min <- min(markers[[y]], na.rm = TRUE)
        if (y_max > 0 && y_min < 0) {
            # center the color bar at 0
            y_max <- max(abs(y_max), abs(y_min))
            y_min <- -y_max
        }
        if (is.numeric(select)) {
            genes <- dplyr::slice_head(markers, n = select, by = !!rlang::sym(subset_by_1))$gene
        } else {
            genes <- dplyr::filter(markers, !!rlang::parse_expr(select))$gene
        }
        genes <- unique(genes)
        markers <- dplyr::filter(markers, !!sym("gene") %in% genes)
        if (!is.factor(markers$gene)) {
            markers$gene <- factor(markers$gene, levels = unique(markers$gene))
        }
        if (!is.factor(markers[[subset_by_1]])) {
            markers[[subset_by_1]] <- factor(markers[[subset_by_1]], levels = unique(markers[[subset_by_1]]))
        }
        args <- list(
            data = markers,
            values_by = y,
            rows_by = "gene",
            columns_by = subset_by_1,
            in_form = "long",
            upper_cutoff = y_max,
            lower_cutoff = y_min,
            ...
        )
        args$show_row_names <- args$show_row_names %||% TRUE
        args$show_column_names <- args$show_column_names %||% TRUE

        # add label if cutoff is provided for heatmap
        genes <- levels(markers$gene)
        groups <- levels(markers[[subset_by_1]])
        if (!is.null(cutoff) && startsWith(plot_type, "heatmap_")) {
            sig_mat <- tidyr::pivot_wider(
                markers,
                id_cols = !!sym("gene"),
                names_from = subset_by_1,
                values_from = !!rlang::sym(pcol),
                values_fill = 1
            )
            sig_mat <- as.data.frame(sig_mat)
            rownames(sig_mat) <- sig_mat$gene
            sig_mat$gene <- NULL
            sig_mat <- sig_mat[genes, groups, drop = FALSE]
            sig_mat <- as.matrix(sig_mat)

            args$cell_type <- "label"
            args$label <- function(x, i, j) {
                pval <- ComplexHeatmap::pindex(sig_mat, i, j)
                ifelse(pval < cutoff, "*", NA)
            }
        }

        # set dot_size fot dot plot
        if (startsWith(plot_type, "dot_")) {
            ds_mat <- tidyr::pivot_wider(
                markers,
                id_cols = !!sym("gene"),
                names_from = subset_by_1,
                values_from = "neg_log10_p"
            )
            ds_mat <- as.data.frame(ds_mat)
            rownames(ds_mat) <- ds_mat$gene
            ds_mat$gene <- NULL
            ds_mat <- ds_mat[genes, groups, drop = FALSE]
            ds_mat <- as.matrix(ds_mat)

            args$cell_type <- "dot"
            args$dot_size <- function(x, i, j) {
                ComplexHeatmap::pindex(ds_mat, i, j)
            }
            args$dot_size_name <- paste0("-log10(", pcol, ")")
        }
        do.call(plotthis::Heatmap, args)
    } else {  # if (plot_type %in% c("heatmap", "violin", "box", "bar", "ridge", "dot")) {
        if (is.null(comparison_by)) {
            stop("[MarkersPlot] `comparison_by` is required for plot_type '", plot_type, "'")
        }
        if (is.null(comparison_by_2)) {
            comparison_by_2 <- comparison_by_1
            if (!comparison_by_2 %in% colnames(object@meta.data)) {
                stop("[MarkersPlot] `comparison_by` '", comparison_by_2, "' is not found in the object's metadata.")
            }
        }
        if (!is.null(subset_by_1) && is.null(subset_by_2)) {
            if (subset_by_1 %in% colnames(object@meta.data)) {
                subset_by_2 <- subset_by_1
            } else {
                warning(
                    "[MarkersPlot] `subset_by` '", subset_by_1, "' is ignored, ",
                    "since it is not found in the object's metadata. ",
                    "Set `subset_by` to '", subset_by_1, ":<object_metadata_column>' to make it work.",
                    immediate. = TRUE)
            }
        }

        if (is.numeric(select)) {
            if (!is.null(subset_by)) {
                genes <- dplyr::slice_head(markers, n = select, by = !!rlang::sym(subset_by_1))$gene
            } else {
                genes <- dplyr::slice_head(markers, n = select)$gene
            }
        } else {
            genes <- dplyr::filter(markers, !!rlang::parse_expr(select))$gene
        }
        genes <- unique(genes)
        # subset the object to only include the selected genes
        object <- Seurat::DietSeurat(object, features = genes)
        # subset the object to only include the comparison groups
        comp_groups <- unique(unlist(strsplit(unique(as.character(markers[[comparison_by_1]])), ":")))
        object <- subset(object, subset = !!rlang::sym(comparison_by_2) %in% comp_groups)

        args <- list(
            object,
            features = genes,
            group_by = comparison_by_2,
            plot_type = plot_type,
            facet_by = if (subset_as_facet) subset_by_2 else NULL,
            split_by = if (!subset_as_facet) subset_by_2 else NULL,
            ident = comparison_by_2,
            ...
        )
        if (plot_type %in% c("heatmap", "dot")) {
            args$name <- args$name %||% "Expression"
            if (!is.null(subset_by_2)) {
                args$facet_by <- NULL
                args$split_by <- NULL
                args$columns_split_by <- subset_by_2
                args$group_by <- NULL
            }
        } else if (plot_type %in% c("violin", "box", "bar")) {
            args$stack <- args$stack %||% TRUE
            if (!is.null(subset_by_2)) {
                args$ident <- subset_by_2
                args$group_by <- comparison_by_2
                args$split_by <- NULL
                args$facet_by <- NULL
            }
        }
        do.call(scplotter::FeatureStatPlot, args)
    }
}
