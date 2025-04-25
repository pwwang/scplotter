#' Cell-Cell Communication Plot
#'
#' @description Plot the cell-cell communication.
#'  See also:
#'  * The review: \url{https://www.sciencedirect.com/science/article/pii/S2452310021000081}
#'  * The LIANA package: \url{https://liana-py.readthedocs.io/en/latest/notebooks/basic_usage.html#Tileplot}
#'  * The CCPlotR package: \url{https://github.com/Sarah145/CCPlotR}
#'
#' @param data A data frame with the cell-cell communication data.
#'  A typical data frame should have the following columns:
#'  * `source` The source cell type.
#'  * `target` The target cell type.
#'  * `ligand` The ligand gene.
#'  * `receptor` The receptor gene.
#'  * `ligand_means` The mean expression of the ligand gene per cell type.
#'  * `receptor_means` The mean expression of the receptor gene per cell type.
#'  * `ligand_props` The proportion of cells that express the entity.
#'  * `receptor_props` The proportion of cells that express the entity.
#'  * `<magnitude>` The magnitude of the communication.
#'  * `<specificity>` The specificity of the communication.
#'  Depends on the `plot_type`, some columns are optional. But the `source`, `target`,
#'  `ligand`, `receptor` and `<magnitude>` are required.
#' @param plot_type The type of plot to use. Default is "dot".
#'  Possible values are "network", "chord", "circos", "heatmap", "sankey", "alluvial", "dot",
#'  "box", "violin" and "ridge".
#'  For "box", "violin" and "ridge", the `method` should be "interaction".
#'  * network: A network plot with the source and target cells as the nodes and the communication as the edges.
#'  * chord: A chord plot with the source and target cells as the nodes and the communication as the chords.
#'  * circos: Alias of "chord".
#'  * heatmap: A heatmap plot with the source and target cells as the rows and columns.
#'  * sankey: A sankey plot with the source and target cells as the nodes and the communication as the flows.
#'  * alluvial: Alias of "sankey".
#'  * dot: A dot plot with the source and target cells as the nodes and the communication as the dots.
#'  * box: Box plots for source cell types. Each x is a target cell type and the values will be
#'    the interaction strengths of the ligand-receptor pairs.
#'  * violin: Violin plots for source cell types. Each x is a target cell type and the values will be
#'    the interaction strengths of the ligand-receptor pairs.
#'  * ridge: Ridge plots for source cell types. Each row is a target cell type and the values will be
#'    the interaction strengths of the ligand-receptor pairs.
#' @param method The method to determine the plot entities.
#'  * aggregation: Aggregate the ligand-receptor pairs interactions for each source-target pair.
#'    Only the source / target pairs will be plotted.
#'  * interaction: Plot the ligand-receptor pairs interactions directly.
#'    The ligand-receptor pairs will also be plotted.
#' @param magnitude The column name in the data to use as the magnitude of the communication.
#'  By default, the second last column will be used.
#'  See `li.mt.show_methods()` for the available methods in `LIANA`.
#'  or \url{https://liana-py.readthedocs.io/en/latest/notebooks/basic_usage.html#Tileplot}
#' @param specificity The column name in the data to use as the specificity of the communication.
#'  By default, the last column will be used.
#'  If the method doesn't have a specificity, set it to NULL.
#' @param magnitude_agg A function to aggregate the magnitude of the communication.
#'  Default is `length`.
#' @param magnitude_name The name of the magnitude in the plot.
#'  Default is "No. of interactions".
#' @param meta_specificity The method to calculate the specificity when there are multiple
#'  ligand-receptor pairs interactions. Default is "sumlog".
#'  It should be one of the methods in the `metap` package.
#' @param split_by A character vector of column names to split the plots. Default is NULL.
#' @param x_text_angle The angle of the x-axis text. Default is 90.
#'  Only used when `plot_type` is "dot".
#' @param link_curvature The curvature of the links. Default is 0.2.
#'  Only used when `plot_type` is "network".
#' @param link_alpha The transparency of the links. Default is 0.6.
#'  Only used when `plot_type` is "network".
#' @param facet_by A character vector of column names to facet the plots. Default is NULL.
#'  It should always be NULL.
#' @param show_row_names Whether to show the row names in the heatmap. Default is TRUE.
#'  Only used when `plot_type` is "heatmap".
#' @param show_column_names Whether to show the column names in the heatmap. Default is TRUE.
#'  Only used when `plot_type` is "heatmap".
#' @param ... Other arguments passed to the specific plot function.
#'  * For `Network`, see [plotthis::Network()].
#'  * For `ChordPlot`, see [plotthis::ChordPlot()].
#'  * For `Heatmap`, see [plotthis::Heatmap()].
#'  * For `SankeyPlot`, see [plotthis::SankeyPlot()].
#'  * For `DotPlot`, see [plotthis::DotPlot()].
#' @return A ggplot object or a list if `combine` is FALSE
#' @importFrom utils getFromNamespace
#' @importFrom rlang syms sym
#' @importFrom dplyr group_by summarise distinct filter n
#' @importFrom tidyr replace_na pivot_wider
#' @importFrom ggplot2 waiver
#' @importFrom plotthis Network ChordPlot Heatmap SankeyPlot DotPlot BoxPlot ViolinPlot RidgePlot
#' @export
#' @examples
#' \donttest{
#' set.seed(8525)
#' data(cellphonedb_res)
#' CCCPlot(data = cellphonedb_res, plot_type = "network", legend.position = "none",
#'   theme = "theme_blank", theme_args = list(add_coord = FALSE))
#' CCCPlot(cellphonedb_res, plot_type = "chord")
#' CCCPlot(cellphonedb_res, plot_type = "heatmap")
#' CCCPlot(cellphonedb_res, plot_type = "dot",
#'   magnitude_agg = mean, magnitude_name = "Average Interaction Strength")
#' CCCPlot(cellphonedb_res, plot_type = "sankey")
#'
#' cellphonedb_res_sub <- cellphonedb_res[
#'   cellphonedb_res$source %in% c("Dendritic", "CD14+ Monocyte"),]
#' CCCPlot(cellphonedb_res_sub, plot_type = "dot", method = "interaction")
#' CCCPlot(cellphonedb_res_sub, plot_type = "network", method = "interaction",
#'   node_size_by = 1)
#' CCCPlot(cellphonedb_res_sub, plot_type = "heatmap", method = "interaction")
#' CCCPlot(cellphonedb_res_sub, plot_type = "box", method = "interaction")
#' CCCPlot(cellphonedb_res_sub, plot_type = "violin", method = "interaction",
#'   add_box = TRUE)
#' CCCPlot(cellphonedb_res_sub, plot_type = "ridge", method = "interaction")
#' }
CCCPlot <- function(
    data,
    plot_type = c("dot", "network", "chord", "circos", "heatmap", "sankey", "alluvial",
        "box", "violin", "ridge"),
    method = c("aggregation", "interaction"),
    magnitude = waiver(),
    specificity = waiver(),
    magnitude_agg = length,
    magnitude_name = "No. of interactions",
    meta_specificity = "sumlog",
    split_by = NULL,
    x_text_angle = 90,
    link_curvature = 0.2,
    link_alpha = 0.6,
    facet_by = NULL,
    show_row_names = TRUE,
    show_column_names = TRUE,
    ...
) {
    stopifnot("[CCCPlot] 'facet_by' is not supported." = is.null(facet_by))

    plot_type <- match.arg(plot_type)
    method <- match.arg(method)
    source_col <- check_columns(data, "source", force_factor = TRUE)
    target_col <- check_columns(data, "target", force_factor = TRUE)
    ligand_col <- check_columns(data, "ligand", force_factor = TRUE)
    receptor_col <- check_columns(data, "receptor", force_factor = TRUE)

    stopifnot("[CCCPlot] Columns 'source', 'target', 'ligand', and 'receptor' are required." =
        !is.null(source_col) && !is.null(target_col) && !is.null(ligand_col) && !is.null(receptor_col))

    if (inherits(magnitude, "waiver")) {
        magnitude <- names(data)[ncol(data) - 1]
        if (!is.numeric(data[[magnitude]])) {
            stop("[CCCPlot] The column '", magnitude, "' is not numeric, specify it explicitly.")
        }
    }
    if (inherits(specificity, "waiver")) {
        specificity <- names(data)[ncol(data)]
        if (!is.numeric(data[[specificity]])) {
            stop("[CCCPlot] The column '", specificity, "' is not numeric, specify it explicitly.")
        }
    }
    stopifnot("[CCCPlot] At least one of 'magnitude' and 'specificity' is required." =
        !inherits(magnitude, "waiver") || !inherits(specificity, "waiver"))

    if (method == "aggregation") {
        links <- data %>% group_by(!!!syms(c(source_col, target_col, split_by)))
        if (is.null(magnitude)) {
            magnitiude <- "mag_score"
            links[[magnitiude]] <- NA
        }

        metap_fn <- getFromNamespace(meta_specificity, "metap")
        links <- suppressWarnings({ links %>%
            filter(!is.na(!!sym(specificity))) %>%
            summarise(
                interactionStrength = magnitude_agg(!!sym(magnitude)),
                .specificity = if (is.null(specificity)) {
                    NA
                } else if (n() == 1) {
                    !!sym(specificity)
                } else {
                    metap_fn(!!sym(specificity))$p
                },
                .groups = "drop") %>%
            replace_na(list(.specificity = 0))
        })

        if (plot_type == "network") {
            Network(links, from = source_col, to = target_col, node_fill_by = "name", split_by = split_by,
                link_curvature = link_curvature, link_weight_name = magnitude_name, link_alpha = link_alpha,
                node_fill_name = "Source/Target", link_weight_by = "interactionStrength", ...)
        } else if (plot_type %in% c("chord", "circos")) {
            ChordPlot(links, y = "interactionStrength", from = source_col, to = target_col,
                split_by = split_by, ...)
        } else if (plot_type == "heatmap") {
            sources <- if (is.factor(links[[source_col]])) {
                levels(links[[source_col]])
            } else {
                unique(links[[source_col]])
            }
            links <- pivot_wider(links, names_from = source_col, values_from = "interactionStrength",
                values_fill = 0)
            Heatmap(links, rows = sources, columns_by = target_col, rows_name = "source", split_by = split_by,
                name = magnitude_name, show_row_names = show_row_names, show_column_names = show_column_names,
                ...)
        } else if (plot_type %in% c("sankey", "alluvial")) {
            SankeyPlot(links, y = "interactionStrength", x = c(source_col, target_col), split_by = split_by,
                links_fill_by = source_col, flow = TRUE, xlab = "", ylab = "Strength", ...)
        } else if (plot_type == "dot") {
            if (!is.null(specificity)) {
                DotPlot(links, x = source_col, y = target_col, size_by = "interactionStrength",
                    fill_by = ".specificity", fill_name = paste0(meta_specificity, "(", specificity, ")"),
                    size_name = magnitude_name, x_text_angle = x_text_angle, split_by = split_by, ...)
            } else {
                DotPlot(links, x = source_col, y = target_col, size_by = "interactionStrength",
                    size_name = magnitude_name, x_text_angle = x_text_angle, split_by = split_by, ...)
            }
        } else {
            stop("[CCCPlot] Plot type '", plot_type, "' is not supported for method 'aggregation' yet.")
        }
    } else if (method == "interaction") {
        stopifnot("[CCCPlot] 'magnitude' is required when 'method' is 'interaction'." =
            !is.null(magnitude))

        if (plot_type == "dot") {
            if (!is.null(specificity)) {
                data[[specificity]] <- -log10(data[[specificity]])
            }
            data[[source_col]] <- paste0("source: ", data[[source_col]])
            DotPlot(data, x = target_col, y = c(ligand_col, receptor_col), y_sep = " -> ",
                fill_by = specificity, fill_name = paste0("-log10(", specificity, ")"),
                size_by = magnitude, x_text_angle = x_text_angle, split_by = split_by,
                facet_by = source_col, ...)
        } else if (plot_type == "network") {
            data$source_target <- paste0(data[[source_col]], " -> ", data[[target_col]])
            Network(data, from = "ligand", to = "receptor",
                link_weight_by = magnitude, link_alpha = link_alpha, link_color_by = "source_target",
                link_color_name = "source -> target", split_by = split_by, ...)
        } else if (plot_type == "heatmap") {
            data$ligand_receptor <- paste0(data[[ligand_col]], " -> ", data[[receptor_col]])
            all_lrs <- unique(data$ligand_receptor)
            data <- pivot_wider(data, names_from = "ligand_receptor", values_from = magnitude,
                values_fill = 0)
            Heatmap(data, rows = all_lrs, rows_name = "Ligand -> Receptor", split_by = split_by,
                name = magnitude, columns_by = target_col, columns_split_by = source_col,
                show_row_names = show_row_names, show_column_names = show_column_names,
                ...)
        } else if (plot_type == "box") {
            data[[source_col]] <- paste0("source: ", data[[source_col]])
            BoxPlot(data, x = target_col, y = magnitude, facet_by = source_col,
                xlab = "Target", ylab = "Interaction Strength", split_by = split_by,
                x_text_angle = x_text_angle, ...)
        } else if (plot_type == "violin") {
            data[[source_col]] <- paste0("source: ", data[[source_col]])
            ViolinPlot(data, x = target_col, y = magnitude, facet_by = source_col,
                xlab = "Target", ylab = "Interaction Strength", split_by = split_by,
                x_text_angle = x_text_angle, ...)
        } else if (plot_type == "ridge") {
            data[[source_col]] <- paste0("source: ", data[[source_col]])
            RidgePlot(data, x = magnitude, group_by = target_col, facet_by = source_col,
                xlab = "Interaction Strength", ylab = "", split_by = split_by,
                ...)
        } else {
            stop("[CCCPlot] Plot type '", plot_type, "' is not supported for method 'interaction' yet.")
        }
    }
}
