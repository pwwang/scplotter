#' Cell statistics plot
#'
#' @description Plot the statistics of the cells.
#'
#' @param object A Seurat object
#' @param ident The column with the cell identities. i.e. clusters. Default: seurat_clusters
#'  For 'pies', this will be used as the `pie_group_by`.
#'  For 'heatmap' plot, this will be used as the rows of the heatmap.
#' @param group_by The column name in the meta data to group the cells. Default: NULL
#'  This should work as the columns of the plot_type: heatmap.
#' @param group_by_sep The separator to use when combining multiple columns in `group_by`. Default: "_"
#'  For 'sankey'/'heatmap' plot, multiple columns will not be combined, and each of them will be used as a node.
#' @param split_by The column name in the meta data to split the cells. Default: NULL
#'  Each split will be plotted in a separate plot.
#' @param split_by_sep The separator to use when combining multiple columns in `split_by`. Default: "_"
#' @param columns_split_by The column name in the meta data to split the columns of the 'pies'/'heatmap' plot. Default: NULL
#' @param rows The column names in the data used as the rows of the 'pies' (heatmap with cell_type = 'pie').
#'  Default: NULL. Only available for 'pies' plot.
#'  The values don't matter, and they only indicate the cells overlapping with the columns and distributed in different `ident` values.
#' @param facet_by The column name in the meta data to facet the plots. Default: NULL
#'  Not available for 'circos', 'sankey', and 'heatmap' plots.
#' @param plot_type The type of plot to use. Default is "bar".
#'   Possible values are "bar", "circos", "pie", "pies", "ring"/"donut", "trend", "area", "heatmap", and "sankey"/"alluvial".
#'   'pie' vs 'pies': 'pie' plot will plot a single pie chart for each group, while 'pies' plot will plot multiple pie charts for each group and split.
#'   'pies' basically is a heatmap with 'cell_type = "pie"'.
#' @param frac The way of calculating the fraction. Default is "none".
#'   Possible values are "group", "ident", "cluster", "all", "none".
#'   Note that the fractions are calculated in each split and facet group if `split_by` and `facet_by` are specified.
#'   * group: calculate the fraction in each group.
#'       The total fraction of the cells of idents in each group will be 1.
#'       When `group-by` is not specified, it will be the same as `all`.
#'   * ident: calculate the fraction in each ident.
#'       The total fraction of the cells of groups in each ident will be 1.
#'       Only works when `group-by` is specified.
#'   * cluster: alias of `ident`.
#'   * all: calculate the fraction against all cells.
#'   * none: do not calculate the fraction, use the number of cells instead.
#' @param swap Whether to swap the cluster and group, that is,
#'   using group as the x-axis and cluster to fill the plot.
#'   For circos plot, when transposed, the arrows will be drawn from the idents (by `ident`) to the
#'   the groups (by `group_by`).
#'   Only works when `group_by` is specified.
#'   For 'pies' plot, this will swap the group_by and ident.
#'   For 'heatmap' plot, this will swap the group_by and columns_split_by.
#'   Note that this is different from `flip`. `flip` is used to transpose the plot. But `swap` is used to swap the x and group_by during plotting.
#' @param rows_name The name of the rows in the 'pies'/'heatmap' plot. Default is NULL.
#' @param name The name of the 'pies'/'heatmap' plot, shown as the name of the main legend. Default is NULL.
#' @param ylab The y-axis label. Default is NULL.
#' @param ... Other arguments passed to the specific plot function.
#'   * For `bar` plot, see \code{\link{plotthis::BarPlot}}.
#'   * For `circos` plot, see \code{\link{plotthis::CircosPlot}}.
#'   * For `pie` chart, see \code{\link{plotthis::PieChart}}.
#'   * For `pies` plot, see \code{\link{plotthis::Heatmap}}.
#'   * For `heatmap` plot, see \code{\link{plotthis::Heatmap}}.
#'   * For `ring`/`donut` plot, see \code{\link{plotthis::RingPlot}}.
#'   * For `trend` plot, see \code{\link{plotthis::TrendPlot}}.
#'   * For `area` plot, see \code{\link{plotthis::AreaPlot}}.
#'   * For `sankey`/`alluvial` plot, see \code{\link{plotthis::SankeyPlot}}.
#'
#' @return A ggplot object or a list if `combine` is FALSE
#' @importFrom rlang sym syms
#' @importFrom dplyr %>% summarise mutate ungroup n
#' @importFrom tidyr drop_na pivot_wider
#' @importFrom plotthis BarPlot CircosPlot PieChart RingPlot TrendPlot AreaPlot SankeyPlot Heatmap
#' @export
#' @examples
#' library(patchwork)
#' data(ifnb_sub)
#'
#' # Bar plot
#' p1 <- CellStatPlot(ifnb_sub)
#' p2 <- CellStatPlot(ifnb_sub, group_by = "stim")
#' p3 <- CellStatPlot(ifnb_sub, group_by = "stim", palette = "Set2")
#' p4 <- CellStatPlot(ifnb_sub, group_by = "stim", position = "stack")
#' (p1 | p2) / (p3 | p4)
#'
#' # Fraction of cells
#' p1 <- CellStatPlot(ifnb_sub, group_by = "stim", frac = "group",
#'                    ident = "seurat_annotations", x_text_angle = 60)
#' p2 <- CellStatPlot(ifnb_sub, group_by = "stim", frac = "group")
#' p3 <- CellStatPlot(ifnb_sub, group_by = "stim", frac = "ident",
#'                    position = "stack", alpha = .6)
#' p4 <- CellStatPlot(ifnb_sub, group_by = "stim", frac = "group",
#'                    swap = TRUE, position = "stack")
#' (p1 | p2) / (p3 | p4)
#'
#' # Splitting/Facetting the plot
#' CellStatPlot(ifnb_sub, split_by = "stim")
#' CellStatPlot(ifnb_sub, facet_by = "stim")
#' CellStatPlot(ifnb_sub, facet_by = "stim", facet_nrow = 2)
#'
#' # Circos plot
#' CellStatPlot(ifnb_sub, group_by = "stim", plot_type = "circos")
#' CellStatPlot(ifnb_sub, group_by = "stim", ident = "seurat_annotations",
#'              plot_type = "circos", labels_rot = TRUE)
#'
#' # Pie plot
#' CellStatPlot(ifnb_sub, plot_type = "pie")
#' CellStatPlot(ifnb_sub, plot_type = "pie", split_by = "stim")
#'
#' # Ring plot
#' CellStatPlot(ifnb_sub, plot_type = "ring", group_by = "stim",
#'              palette = "Spectral")
#'
#' # Trend plot
#' CellStatPlot(ifnb_sub, plot_type = "trend", frac = "group",
#'              x_text_angle = 90, group_by = c("stim", "seurat_annotations"))
#'
#' # Sankey plot
#' CellStatPlot(ifnb_sub, plot_type = "sankey", group_by = "stim", alpha = .6)
#' CellStatPlot(ifnb_sub, plot_type = "sankey", alpha = .6,
#'              group_by = c("stim", "seurat_annotations", "orig.ident"))
#'
#' # Area plot
#' CellStatPlot(ifnb_sub, plot_type = "area", frac = "group", x_text_angle = 90,
#'              group_by = "seurat_annotations", split_by = "stim")
#'
#' # Pies
#' ifnb_sub$r1 <- sample(c(1, NA), ncol(ifnb_sub), replace = TRUE)
#' ifnb_sub$r2 <- sample(c(1, NA), ncol(ifnb_sub), replace = TRUE)
#' ifnb_sub$r3 <- sample(c(1, NA), ncol(ifnb_sub), replace = TRUE)
#' CellStatPlot(ifnb_sub, plot_type = "pies", group_by = "stim",
#'    rows = c("r1", "r2", "r3"), show_row_names = TRUE,
#'    show_column_names = TRUE, column_names_side = "top", cluster_columns = FALSE,
#'    row_names_side = "left", pie_size = sqrt)
#'
#' # Heatmap
#' CellStatPlot(ifnb_sub, plot_type = "heatmap", group_by = "stim", palette = "Blues")
#' CellStatPlot(ifnb_sub, plot_type = "heatmap", group_by = "stim",
#'    frac = "group", columns_split_by = "seurat_annotations", swap = TRUE)
CellStatPlot <- function(
    object, ident = "seurat_clusters", group_by = NULL, group_by_sep = "_",
    split_by = NULL, split_by_sep = "_", facet_by = NULL, rows = NULL, columns_split_by = NULL,
    frac = c("none", "group", "ident", "cluster", "all"), rows_name = NULL, name = NULL,
    plot_type = c("bar", "circos", "pie", "pies", "ring", "donut", "trend", "area", "sankey", "alluvial", "heatmap"),
    swap = FALSE, ylab = NULL, ...
) {
    data <- object@meta.data

    plot_type <- match.arg(plot_type)
    if (plot_type == "donut") plot_type <- "ring"
    if (plot_type == "alluvial") plot_type <- "sankey"
    if (isFALSE(swap) && plot_type %in% c("sankey", "heatmap")) {
        group_by <- plotthis:::check_columns(data, group_by, force_factor = TRUE,
            allow_multi = TRUE)
    } else {
        group_by <- plotthis:::check_columns(data, group_by, force_factor = TRUE,
            allow_multi = TRUE, concat_multi = TRUE, concat_sep = group_by_sep)
    }
    split_by <- plotthis:::check_columns(data, split_by, force_factor = TRUE,
        allow_multi = TRUE, concat_multi = TRUE, concat_sep = split_by_sep)
    facet_by <- plotthis:::check_columns(data, facet_by, force_factor = TRUE,
        allow_multi = TRUE)
    rows <- plotthis:::check_columns(data, rows, force_factor = TRUE,
        allow_multi = TRUE)

    frac <- match.arg(frac)
    if (frac == "cluster") frac <- "ident"
    if (identical(plot_type, "ring") && !identical(frac, "group")) {
        message("'frac' is forced to 'group' for 'ring' plot.")
        frac <- "group"
    }
    if (is.null(group_by) && frac == "ident") {
        stop("Cannot calculate the fraction by 'ident' without specifying 'group_by'.")
    }

    if (isTRUE(swap) && is.null(group_by)) {
        stop("Cannot swap the 'ident' and 'group_by' without specifying 'group_by'.")
    }

    data <- data %>% drop_na(!!sym(ident))
    if (nrow(data) == 0) {
        stop("No cells found with ident:'", ident, "'")
    }

    if (plot_type != "pies") {
        # Heatmap itself will handle the data for cell_type = pie
        if (is.null(group_by)) {
            data <- data %>%
                dplyr::group_by(!!!syms(unique(c(split_by, facet_by, columns_split_by, ident)))) %>%
                summarise(.n = n(), .groups = "drop")

            if (frac != "none") {
                data <- data %>%
                    dplyr::group_by(!!!syms(unique(c(split_by, facet_by, columns_split_by)))) %>%
                    mutate(.frac = .n / sum(.n))
            } else {
                data <- data %>% mutate(.frac = 1)  # not used
            }
        } else if (length(group_by) > 1) {
            # calculate the frac for each group. we don't want them to be concatenated.
            data <- do.call(rbind, lapply(group_by, function(g) {
                dat <- data %>%
                    dplyr::group_by(!!!syms(unique(c(split_by, facet_by, columns_split_by, ident, g)))) %>%
                    summarise(.n = n(), .groups = "drop")
                dat <- dat[!is.na(dat[[g]]), , drop = FALSE]

                if (frac == "group") {
                    dat <- dat %>%
                        dplyr::group_by(!!!syms(unique(c(split_by, facet_by, g, columns_split_by)))) %>%
                        mutate(.frac = .n / sum(.n)) %>%
                        ungroup()
                } else if (frac == "ident") {
                    dat <- dat %>%
                        dplyr::group_by(!!!syms(unique(c(split_by, facet_by, columns_split_by, ident)))) %>%
                        mutate(.frac = .n / sum(.n)) %>%
                        ungroup()
                } else if (frac == "all") {
                    dat <- dat %>%
                        dplyr::group_by(!!!syms(unique(c(split_by, facet_by, columns_split_by)))) %>%
                        mutate(.frac = .n / sum(.n)) %>%
                        ungroup()
                } else {
                    dat <- dat %>% mutate(.frac = 1)  # not used
                }

                dat[, setdiff(group_by, g)] <- NA
                dat
            }))
        } else {
            data <- data %>%
                dplyr::group_by(!!!syms(unique(c(split_by, facet_by, group_by, columns_split_by, ident)))) %>%
                summarise(.n = n(), .groups = "drop")

            if (frac == "group") {
                data <- data %>%
                    dplyr::group_by(!!!syms(unique(c(split_by, facet_by, group_by, columns_split_by)))) %>%
                    mutate(.frac = .n / sum(.n)) %>%
                    ungroup()
            } else if (frac == "ident") {
                data <- data %>%
                    dplyr::group_by(!!!syms(unique(c(split_by, facet_by, columns_split_by, ident)))) %>%
                    mutate(.frac = .n / sum(.n)) %>%
                    ungroup()
            } else if (frac == "all") {
                data <- data %>%
                    dplyr::group_by(!!!syms(unique(c(split_by, facet_by, columns_split_by)))) %>%
                    mutate(.frac = .n / sum(.n)) %>%
                    ungroup()
            } else {
                data <- data %>% mutate(.frac = 1)  # not used
            }
        }
    }

    ylab <- ylab %||% paste0(ifelse(identical(frac, "none"), "Number", "Fraction"), " of cells")
    if (plot_type == "bar") {
        if (is.null(group_by)) {
            BarPlot(
                data, x = ident, y = ifelse(identical(frac, "none"), ".n", ".frac"),
                split_by = split_by, facet_by = facet_by, ylab = ylab, ...)
        } else {
            BarPlot(
                data,
                x = ifelse(swap, group_by, ident),
                y = ifelse(identical(frac, "none"), ".n", ".frac"),
                group_by = ifelse(swap, ident, group_by),
                split_by = split_by, facet_by = facet_by, ylab = ylab, ...)
        }
    } else if (plot_type == "circos") {
        if (is.null(group_by)) {
            stop("Cannot create a circos plot without specifying 'group_by'.")
        }
        CircosPlot(
            data,
            from = ifelse(swap, ident, group_by),
            to = ifelse(swap, group_by, ident),
            split_by = split_by, facet_by = facet_by, ...)
    } else if (plot_type == "pie") {
        if (!is.null(group_by)) {
            stop("Cannot create a pie plot with 'group_by'. You may want to use 'ring' plot instead.")
        }
        PieChart(
            data, x = ident, y = ifelse(identical(frac, "none"), ".n", ".frac"),
            split_by = split_by, facet_by = facet_by, ...)
    } else if (plot_type == "ring") {
        if (!identical(frac, "group")) {
            stop("'frac' must be 'group' for 'ring' plot.")
        }
        RingPlot(
            data,
            x = ifelse(swap, ident, group_by),
            y = ifelse(identical(frac, "none"), ".n", ".frac"),
            group_by = ifelse(swap, group_by, ident),
            split_by = split_by, facet_by = facet_by, ...)
    } else if (plot_type == "trend") {
        if (is.null(group_by)) {
            stop("Cannot create a trend plot without specifying 'group_by'.")
        }
        TrendPlot(
            data,
            x = ifelse(swap, ident, group_by),
            y = ifelse(identical(frac, "none"), ".n", ".frac"),
            group_by = ifelse(swap, group_by, ident),
            ylab = ylab, split_by = split_by, facet_by = facet_by, ...)
    } else if (plot_type == "area") {
        if (is.null(group_by)) {
            stop("Cannot create an area plot without specifying 'group_by'.")
        }
        AreaPlot(
            data,
            x = ifelse(swap, ident, group_by),
            y = ifelse(identical(frac, "none"), ".n", ".frac"),
            group_by = ifelse(swap, group_by, ident),
            ylab = ylab, split_by = split_by, facet_by = facet_by, ...)
    } else if (plot_type == "sankey") {
        if (is.null(group_by)) {
            stop("Cannot create a sankey plot without specifying 'group_by'.")
        }
        if (frac == "ident") {
            stop("Cannot calculate the fraction by 'ident' for 'trend' plot.")
        }
        if (isTRUE(swap)) {
            stop("'swap = TRUE' is not supported for 'sankey' plot.")
        }
        SankeyPlot(
            data,
            nodes_by = c(ident, group_by),
            # links_by = ident,
            y = ifelse(identical(frac, "none"), ".n", ".frac"),
            split_by = split_by, facet_by = facet_by, ...)
    } else if (plot_type == "pies") {
        if (is.null(group_by)) {
            stop("Cannot create a heatmap (cell_type = 'pie') plot without specifying 'group_by'.")
        }
        if (is.null(rows)) {
            stop("Cannot create a heatmap (cell_type = 'pie') plot without specifying 'rows'.")
        }
        if (!is.null(facet_by)) {
            stop("Cannot create a heatmap (cell_type = 'pie') plot with 'facet_by'.")
        }

        rows_name <- rows_name %||% "rows"
        name <- name %||% "value"

        Heatmap(
            data, rows = rows, cell_type = "pie", rows_name = rows_name, name = name,
            columns_by = if (swap) ident else group_by,
            pie_group_by = if (swap) group_by else ident,
            columns_split_by = columns_split_by, split_by = split_by, ...)
    } else if (plot_type == "heatmap") {
        if (is.null(group_by)) {
            stop("Cannot create a heatmap plot without specifying 'group_by', which should work as the columns of the heatmap.")
        }
        if (!is.null(facet_by)) {
            stop("Cannot create a heatmap plot with 'facet_by'.")
        }
        if (swap && is.null(columns_split_by)) {
            stop("Cannot swap between 'group_by' and 'columns_split_by' without specifying 'columns_split_by'.")
        }
        idents <- if (is.factor(data[[ident]])) levels(data[[ident]]) else unique(data[[ident]])
        data <- pivot_wider(data, id_cols = unique(c(split_by, group_by, columns_split_by)),
            names_from = ident, values_fill = 0,
            values_from = if (identical(frac, "none")) ".n" else ".frac")

        rows_name <- rows_name %||% ident
        name <- name %||% ifelse(identical(frac, "none"), "Number of cells", "Fraction of cells")

        Heatmap(
            data, rows = idents, rows_name = rows_name, name = name,
            columns_by = if (swap) columns_split_by else group_by,
            columns_split_by = if (swap) group_by else columns_split_by,
            split_by = split_by, ...)
    }
}
