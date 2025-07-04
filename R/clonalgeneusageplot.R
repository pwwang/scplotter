#' ClonalGeneUsagePlot
#'
#' @param data The product of [scRepertoire::combineTCR], [scRepertoire::combineTCR], or
#'  [scRepertoire::combineExpression].
#' @param genes The prefix of genes to be plotted. Default is "TRBV".
#'  If two sets of genes are provided (e.g. c("TRBV", "TRBJ")), the second dimension will be the
#'  second set of genes instead of the group_by variable.
#' @param scale Whether to use the proportion that is scaled to the group or the count.
#' @param top The number of top genes/genepairs to be plotted.
#' @param plot_type The type of plot to be generated. Default is "bar".
#'  Options are "bar", "heatmap", "circos" (aka "chord").
#' @param group_by The variable to group the data by. Default is "Sample".
#' @param facet_by A character vector of column names to facet the plots. Default is NULL.
#'  Should not be specified manually.
#' @param facet_ncol The number of columns in the facet grid. Default is 1.
#' @param split_by A character vector of column names to split the plots. Default is NULL.
#' @param aspect.ratio The aspect ratio of the plot. Only available for "bar" plot. Default is 2/top.
#' @param theme_args A list of arguments to be passed to the [ggplot2::theme] function.
#' @param ylab The y-axis label. Default is NULL.
#' @param row_annotation A list of row annotations to be added to the heatmap. Default is NULL.
#' @param row_annotation_type A list of row annotation types.
#' @param row_annotation_side The side of the row annotation. Default is "right".
#' @param row_annotation_agg A list of row annotation aggregation functions.
#' @param show_row_names Whether to show row names in the heatmap. Default is TRUE.
#' @param show_column_names Whether to show column names in the heatmap. Default is TRUE.
#' @param ... Other arguments passed to the specific plot function.
#'  * For "bar", [plotthis::BarPlot()].
#'  * For "heatmap", [plotthis::Heatmap()].
#'  * For "circos", [plotthis::ChordPlot()].
#'  * For "chord", [plotthis::ChordPlot()].
#' @return A ggplot object or a list if `combine` is FALSE
#' @importFrom rlang := syms
#' @importFrom dplyr mutate rename slice_max summarise pull filter all_of ungroup
#' @importFrom tidyr separate
#' @importFrom ggplot2 unit element_blank
#' @importFrom scRepertoire vizGenes
#' @importFrom plotthis BarPlot Heatmap ChordPlot SankeyPlot
#' @export
#' @examples
#' \donttest{
#' set.seed(8525)
#' data(contig_list, package = "scRepertoire")
#' data <- scRepertoire::combineTCR(contig_list,
#'     samples = c("P17B", "P17L", "P18B", "P18L", "P19B","P19L", "P20B", "P20L"))
#' data <- scRepertoire::addVariable(data,
#'     variable.name = "Type",
#'     variables = rep(c("B", "L"), 4)
#' )
#' data <- scRepertoire::addVariable(data,
#'     variable.name = "Subject",
#'     variables = rep(c("P17", "P18", "P19", "P20"), each = 2)
#' )
#'
#' ClonalGeneUsagePlot(data)
#' ClonalGeneUsagePlot(data, genes = c("TRBJ", "TRBV"))
#' ClonalGeneUsagePlot(data, top = 40, plot_type = "heatmap")
#' ClonalGeneUsagePlot(data, genes = c("TRBV", "TRBJ"), plot_type = "heatmap")
#' ClonalGeneUsagePlot(data, genes = "TRBV", group_by = "Type", plot_type = "chord")
#' ClonalGeneUsagePlot(data, genes = c("TRBV", "TRBJ"), group_by = "Type", plot_type = "chord")
#' ClonalGeneUsagePlot(data, genes = c("TRBV", "TRBJ"), plot_type = "alluvial",
#'      facet_scales = "free_y")
#' ClonalGeneUsagePlot(data, genes = c("TRBV", "TRBJ"), plot_type = "alluvial",
#'      group_by = NULL)
#' }
ClonalGeneUsagePlot <- function(
    data, genes = "TRBV", scale = TRUE, top = 20,
    plot_type = c("bar", "heatmap", "circos", "chord", "alluvial", "sankey"),
    group_by = "Sample", facet_by = NULL, facet_ncol = 1, split_by = NULL,
    aspect.ratio = 2 / top, theme_args = list(), ylab = NULL,
    show_row_names = TRUE, show_column_names = TRUE,
    row_annotation = NULL, row_annotation_type = list(), row_annotation_side = "right",
    row_annotation_agg = list(), ...
) {
    plot_type <- match.arg(plot_type)
    if (plot_type == "alluvial") { plot_type <- "sankey" }
    if (!is.null(facet_by)) {
        stop("'facet_by' should not be specified in ClonalGeneUsagePlot.")
    }

    if (length(genes) > 1) {
        genes2 <- genes[2]
        genes <- genes[1]
    } else {
        genes2 <- NULL
    }

    all_groupings <- unique(c(group_by, split_by))
    data <- merge_clonal_groupings(data, all_groupings)
    data <- vizGenes(data, x.axis = genes, y.axis = genes2, group.by = ".group", scale = scale, exportTable = TRUE)
    if (is.null(genes2)) {
        stopifnot("[ClonalGeneUsagePlot] 'genes' must be of length 2." = plot_type != "sankey")
        axis1 <- genes
        axis2 <- group_by

        data <- separate(data, "y.axis", into = all_groupings, sep = " // ") %>%
            rename(!!sym(axis1) := "x.axis")

        selected_genes <- data %>% dplyr::group_by(!!sym(genes)) %>%
            summarise(total = sum(!!sym("count")), .groups = "drop") %>%
            slice_max(order_by = !!sym("total"), n = top) %>%
            pull(axis1)%>%
            as.character()

        selected_genes_levels <- if (is.factor(data[[axis1]])) {
            levels(data[[axis1]])
        } else {
            unique(data[[axis1]])
        }
        selected_genes_levels <- selected_genes_levels[selected_genes_levels %in% selected_genes]

        data <- data %>% filter(!!sym(axis1) %in% selected_genes) %>%
            mutate(!!sym(axis1) := factor(!!sym(axis1), levels = selected_genes_levels))
    } else {
        axis1 <- genes
        axis2 <- genes2

        if (!is.null(all_groupings)) {
            data <- separate(data, "element.names", into = all_groupings, sep = " // ")
        }
        data <- data %>%
            rename(!!sym(genes) := "x.axis", !!sym(genes2) := "y.axis") %>%
            unite("GenePairs", c(genes, genes2), sep = " // ", remove = FALSE)

        genepairs <- data %>% dplyr::group_by(!!!syms(c(genes, genes2))) %>%
            summarise(total = sum(!!sym("count")), .groups = "drop") %>%
            slice_max(order_by = !!sym("total"), n = top) %>%
            unite("GenePairs", c(genes, genes2), sep = " // ") %>%
            pull("GenePairs")

        data <- data %>% filter(!!sym("GenePairs") %in% genepairs)
        data$GenePairs <- NULL
    }

    if (plot_type == "bar") {
        # theme_args$panel.spacing <- theme_args$panel.spacing %||% unit(-0.1, "lines")
        theme_args$panel.grid.major.y <- theme_args$panel.grid.major.y %||% element_blank()
        ylab <- ylab %||% ifelse(scale, "Gene Usage Fraction", "Gene Usage Count")

        BarPlot(data, x = axis1, y = ifelse(scale, "proportion", "count"), facet_by = axis2, facet_ncol = facet_ncol,
            split_by = split_by, legend.position = "none", x_text_angle = 90, aspect.ratio = aspect.ratio,
            facet_args = list(strip.position = "right"), theme_args = theme_args, ylab = ylab, ...)
    } else if (plot_type == "heatmap") {
        row_annotation_type[["Total Usage"]] <- row_annotation_type[["Total Usage"]] %||% "lines"
        row_annotation_agg[["Total Usage"]] <- row_annotation_agg[["Total Usage"]] %||% function(x) ifelse(length(x) > 1, x[1], 0)
        if (is.null(genes2)) {
            data <- data %>% dplyr::group_by(!!!syms(c(axis1, split_by))) %>%
                mutate(.total = sum(!!sym(ifelse(scale, "proportion", "count")))) %>%
                ungroup()
            row_annotation <- row_annotation %||% list(`Total Usage` = ".total")
        }

        Heatmap(
            data,
            values_by = ifelse(scale, "proportion", "count"), values_fill = 0,
            rows_by = axis1, columns_by = axis2, split_by = split_by,
            rows_name = axis1, name = ifelse(scale, "Gene Usage Fraction", "Gene Usage Count"),
            row_annotation = row_annotation, row_annotation_type = row_annotation_type, row_annotation_agg = row_annotation_agg,
            row_annotation_side = row_annotation_side,
            show_row_names = show_row_names, show_column_names = show_column_names,
            ...)
    } else if (plot_type %in% c("circos", "chord")) {
        if (is.null(genes2)) {
            ChordPlot(data, from = axis1, to = axis2, y = ifelse(scale, "proportion", "count"),
                split_by = split_by, theme_args = theme_args, ...)
        } else {
            if (!is.null(split_by)) {
                stop("[ClonalGeneUsagePlot] 'split_by' should not be specified when 'genes' has length 2, since the plot will be split by 'group_by'.")
            }
            ChordPlot(data, from = axis1, to = axis2, y = ifelse(scale, "proportion", "count"),
                split_by = group_by, theme_args = theme_args, ...)
        }
    } else {  # alluvial / sankey
        SankeyPlot(data, x = c(axis1, axis2), y = ifelse(scale, "proportion", "count"),
            links_fill_by = axis1, facet_by = group_by,
            split_by = split_by, theme_args = theme_args, ...)
    }
}
