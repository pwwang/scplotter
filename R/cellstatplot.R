#' Cell statistics plot
#'
#' @description This function creates a plot to visualize the statistics of cells in a Seurat object, a Giotto object,
#' a path to an .h5ad file or an opened `H5File` by `hdf5r` package.
#' It can create various types of plots, including bar plots, circos plots, pie charts, pies (heatmap with cell_type = 'pie'), ring/donut plots, trend plots
#' area plots, sankey/alluvial plots, heatmaps, radar plots, spider plots, violin plots, and box plots.
#' The function allows for grouping, splitting, and faceting the data based on metadata columns.
#' It also supports calculating fractions of cells based on specified groupings.#'
#' @param object A Seurat object, a Giotto object, a path to h5ad file or an opened `H5File` (from `hdf5r` package) object a data frame (for internal use) containing cell metadata.
#' @param ident The column with the cell identities. i.e. clusters. Default: NULL
#'  If NULL, the active identity of the Seurat object and the name "Identity" will be used.
#'  For 'pies', this will be used as the `pie_group_by`.
#'  For 'heatmap' plot, this will be used as the rows of the heatmap.
#' @param group_by The column name in the meta data to group the cells. Default: NULL
#'  This should work as the columns of the plot_type: heatmap.
#'  For violin/box plot, at most 2 `group_by` columns are allowed and they will not be concatenated.
#'  The first one is used to break down the values in groups, and the second one works as the `group_by`
#'  argument in [plotthis::ViolinPlot]/[plotthis::BoxPlot].
#' @param group_by_sep The separator to use when combining multiple columns in `group_by`. Default: "_"
#'  For 'sankey'/'heatmap' plot, multiple columns will not be combined, and each of them will be used as a node.
#' @param spat_unit The spatial unit to use for the plot. Only applied to Giotto objects.
#' @param feat_type feature type of the features (e.g. "rna", "dna", "protein"), only applied to Giotto objects.
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
#'   Possible values are "bar", "circos", "pie", "pies", "ring"/"donut", "trend", "area", "heatmap", "sankey"/"alluvial", "radar" and "spider".
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
#'   * For `bar` plot, see [plotthis::BarPlot()].
#'   * For `circos` plot, see [plotthis::CircosPlot()].
#'   * For `pie` chart, see [plotthis::PieChart()].
#'   * For `pies` plot, see [plotthis::Heatmap()].
#'   * For `heatmap` plot, see [plotthis::Heatmap()].
#'   * For `ring`/`donut` plot, see [plotthis::RingPlot()].
#'   * For `trend` plot, see [plotthis::TrendPlot()].
#'   * For `area` plot, see [plotthis::AreaPlot()].
#'   * For `sankey`/`alluvial` plot, see [plotthis::SankeyPlot()].
#'   * For `radar` plot, see [plotthis::RadarPlot()].
#'   * For `spider` plot, see [plotthis::SpiderPlot()].
#'   * For `violin` plot, see [plotthis::ViolinPlot()].
#'   * For `box` plot, see [plotthis::BoxPlot()].
#'
#' @return A ggplot object or a list if `combine` is FALSE
#' @importFrom rlang sym syms
#' @importFrom SeuratObject Idents
#' @importFrom dplyr %>% summarise mutate ungroup n
#' @importFrom tidyr drop_na pivot_wider pivot_longer
#' @importFrom plotthis BarPlot CircosPlot PieChart RingPlot TrendPlot AreaPlot SankeyPlot Heatmap RadarPlot SpiderPlot ViolinPlot BoxPlot
#' @details See:
#' * <https://pwwang.github.io/scplotter/articles/Giotto_Xenium.html>
#'
#' for examples of using this function with a Giotto object.
#'
#' And see:
#' * <https://pwwang.github.io/scplotter/articles/Working_with_anndata_h5ad_files.html>
#'
#' for examples of using this function with .h5ad files.
#' @export
#' @examples
#' \donttest{
#' # library(patchwork)
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
#' CellStatPlot(ifnb_sub, plot_type = "sankey", group_by = c("seurat_clusters", "stim"),
#'              links_alpha = .6)
#' CellStatPlot(ifnb_sub, plot_type = "sankey", links_alpha = .6,
#'              group_by = c("stim", "seurat_annotations", "orig.ident"))
#' CellStatPlot(ifnb_sub, plot_type = "sankey", links_alpha = .6,
#'              group_by = c("seurat_clusters", "stim", "seurat_annotations", "orig.ident"))
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
#'
#' # Radar plot/Spider plot
#' pr <- CellStatPlot(ifnb_sub, plot_type = "radar", group_by = "stim")
#' ps <- CellStatPlot(ifnb_sub, plot_type = "spider", group_by = "stim")
#' pr | ps
#'
#' # Box/Violin plot
#' ifnb_sub$group <- sample(paste0("g", 1:10), nrow(ifnb_sub), replace = TRUE)
#' CellStatPlot(ifnb_sub, group_by = c("group", "stim"), frac = "group",
#'    plot_type = "violin", add_box = TRUE, ident = "seurat_annotations",
#'    x_text_angle = 60, comparisons = TRUE, aspect.ratio = 0.8)
#' }
CellStatPlot <- function(
    object, ident = NULL, group_by = NULL, group_by_sep = "_", spat_unit = NULL, feat_type = NULL,
    split_by = NULL, split_by_sep = "_", facet_by = NULL, rows = NULL, columns_split_by = NULL,
    frac = c("none", "group", "ident", "cluster", "all"), rows_name = NULL, name = NULL,
    plot_type = c("bar", "circos", "pie", "pies", "ring", "donut", "trend", "area", "sankey", "alluvial", "heatmap", "radar", "spider", "violin", "box"),
    swap = FALSE, ylab = NULL, ...
) {
    UseMethod("CellStatPlot")
}

#' @export
CellStatPlot.giotto <- function(
    object, ident = NULL, group_by = NULL, group_by_sep = "_", spat_unit = NULL, feat_type = NULL,
    split_by = NULL, split_by_sep = "_", facet_by = NULL, rows = NULL, columns_split_by = NULL,
    frac = c("none", "group", "ident", "cluster", "all"), rows_name = NULL, name = NULL,
    plot_type = c("bar", "circos", "pie", "pies", "ring", "donut", "trend", "area", "sankey", "alluvial", "heatmap", "radar", "spider", "violin", "box"),
    swap = FALSE, ylab = NULL, ...
) {
    stopifnot("[CellStatPlot] 'ident' is required for Giotto object." = !is.null(ident))

    spat_unit <- GiottoClass::set_default_spat_unit(
        gobject = object,
        spat_unit = spat_unit
    )

    feat_type <- GiottoClass::set_default_feat_type(
        gobject = object,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    data <- GiottoClass::getCellMetadata(
        gobject = object,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "data.table"
    )

    CellStatPlot.data.frame(
        data, ident = ident, group_by = group_by, group_by_sep = group_by_sep,
        split_by = split_by, split_by_sep = split_by_sep, facet_by = facet_by,
        rows = rows, columns_split_by = columns_split_by, frac = frac,
        rows_name = rows_name, name = name, plot_type = plot_type,
        swap = swap, ylab = ylab, ...
    )
}

#' @export
CellStatPlot.Seurat <- function(
    object, ident = NULL, group_by = NULL, group_by_sep = "_", spat_unit = NULL, feat_type = NULL,
    split_by = NULL, split_by_sep = "_", facet_by = NULL, rows = NULL, columns_split_by = NULL,
    frac = c("none", "group", "ident", "cluster", "all"), rows_name = NULL, name = NULL,
    plot_type = c("bar", "circos", "pie", "pies", "ring", "donut", "trend", "area", "sankey", "alluvial", "heatmap", "radar", "spider", "violin", "box"),
    swap = FALSE, ylab = NULL, ...
) {
    data <- object@meta.data
    if (is.null(ident)) {
        ident <- "Identity"
        data[[ident]] <- Idents(object)
    }

    CellStatPlot.data.frame(
        data, ident = ident, group_by = group_by, group_by_sep = group_by_sep,
        split_by = split_by, split_by_sep = split_by_sep, facet_by = facet_by,
        rows = rows, columns_split_by = columns_split_by, frac = frac,
        rows_name = rows_name, name = name, plot_type = plot_type,
        swap = swap, ylab = ylab, ...
    )
}

#' @export
CellStatPlot.character <- function(
    object, ident = NULL, group_by = NULL, group_by_sep = "_", spat_unit = NULL, feat_type = NULL,
    split_by = NULL, split_by_sep = "_", facet_by = NULL, rows = NULL, columns_split_by = NULL,
    frac = c("none", "group", "ident", "cluster", "all"), rows_name = NULL, name = NULL,
    plot_type = c("bar", "circos", "pie", "pies", "ring", "donut", "trend", "area", "sankey", "alluvial", "heatmap", "radar", "spider", "violin", "box"),
    swap = FALSE, ylab = NULL, ...
) {
    if (!endsWith(object, ".h5ad")) {
        stop("[CellStatPlot] Currently only supports .h5ad files when called with a string/path.")
    }

    object <- hdf5r::H5File$new(object, mode = "r")
    on.exit(object$close_all())

    CellStatPlot.H5File(
        object, ident = ident, group_by = group_by, group_by_sep = group_by_sep,
        spat_unit = spat_unit, feat_type = feat_type,
        split_by = split_by, split_by_sep = split_by_sep, facet_by = facet_by,
        rows = rows, columns_split_by = columns_split_by, frac = frac,
        rows_name = rows_name, name = name, plot_type = plot_type,
        swap = swap, ylab = ylab, ...
    )
}

#' @export
CellStatPlot.H5File <- function(
    object, ident = NULL, group_by = NULL, group_by_sep = "_", spat_unit = NULL, feat_type = NULL,
    split_by = NULL, split_by_sep = "_", facet_by = NULL, rows = NULL, columns_split_by = NULL,
    frac = c("none", "group", "ident", "cluster", "all"), rows_name = NULL, name = NULL,
    plot_type = c("bar", "circos", "pie", "pies", "ring", "donut", "trend", "area", "sankey", "alluvial", "heatmap", "radar", "spider", "violin", "box"),
    swap = FALSE, ylab = NULL, ...
) {
    stopifnot("[CellStatPlot] 'ident' is required for anndata (h5ad) object." = !is.null(ident))

    object <- h5group_to_dataframe(object[["obs"]])

    CellStatPlot.data.frame(
        object, ident = ident, group_by = group_by, group_by_sep = group_by_sep,
        spat_unit = spat_unit, feat_type = feat_type,
        split_by = split_by, split_by_sep = split_by_sep, facet_by = facet_by,
        rows = rows, columns_split_by = columns_split_by, frac = frac,
        rows_name = rows_name, name = name, plot_type = plot_type,
        swap = swap, ylab = ylab, ...
    )
}

#' @export
CellStatPlot.data.frame <- function(
    object, ident = NULL, group_by = NULL, group_by_sep = "_", spat_unit = NULL, feat_type = NULL,
    split_by = NULL, split_by_sep = "_", facet_by = NULL, rows = NULL, columns_split_by = NULL,
    frac = c("none", "group", "ident", "cluster", "all"), rows_name = NULL, name = NULL,
    plot_type = c("bar", "circos", "pie", "pies", "ring", "donut", "trend", "area", "sankey", "alluvial", "heatmap", "radar", "spider", "violin", "box"),
    swap = FALSE, ylab = NULL, ...
) {
    plot_type <- match.arg(plot_type)
    if (plot_type == "donut") plot_type <- "ring"
    if (plot_type == "alluvial") plot_type <- "sankey"
    if (isFALSE(swap) && plot_type %in% c("sankey", "heatmap", "violin", "box")) {
        group_by <- check_columns(object, group_by, force_factor = TRUE,
            allow_multi = TRUE)
    } else {
        group_by <- check_columns(object, group_by, force_factor = TRUE,
            allow_multi = TRUE, concat_multi = TRUE, concat_sep = group_by_sep)
    }
    split_by <- check_columns(object, split_by, force_factor = TRUE,
        allow_multi = TRUE, concat_multi = TRUE, concat_sep = split_by_sep)
    facet_by <- check_columns(object, facet_by, force_factor = TRUE,
        allow_multi = TRUE)
    rows <- check_columns(object, rows, force_factor = TRUE,
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

    object <- object %>% drop_na(!!sym(ident))
    if (nrow(object) == 0) {
        stop("No cells found with ident:'", ident, "'")
    }

    if (plot_type != "pies") {
        # Heatmap itself will handle the data for cell_type = pie
        if (is.null(group_by)) {
            object <- object %>%
                dplyr::group_by(!!!syms(unique(c(split_by, facet_by, columns_split_by, ident)))) %>%
                summarise(.n = n(), .groups = "drop")

            if (frac != "none") {
                object <- object %>%
                    dplyr::group_by(!!!syms(unique(c(split_by, facet_by, columns_split_by)))) %>%
                    mutate(.frac = !!sym(".n") / sum(!!sym(".n")))
            } else {
                object <- object %>% mutate(.frac = 1)  # not used
            }
        } else if (length(group_by) > 1 && !plot_type %in% c("sankey", "violin", "box")) {
            # calculate the frac for each group. we don't want them to be concatenated.
            object <- do.call(rbind, lapply(group_by, function(g) {
                dat <- object %>%
                    dplyr::group_by(!!!syms(unique(c(split_by, facet_by, columns_split_by, ident, g)))) %>%
                    summarise(.n = n(), .groups = "drop")
                dat <- dat[!is.na(dat[[g]]), , drop = FALSE]

                if (frac == "group") {
                    dat <- dat %>%
                        dplyr::group_by(!!!syms(unique(c(split_by, facet_by, g, columns_split_by)))) %>%
                        mutate(.frac = !!sym(".n") / sum(!!sym(".n"))) %>%
                        ungroup()
                } else if (frac == "ident") {
                    dat <- dat %>%
                        dplyr::group_by(!!!syms(unique(c(split_by, facet_by, columns_split_by, ident)))) %>%
                        mutate(.frac = !!sym(".n") / sum(!!sym(".n"))) %>%
                        ungroup()
                } else if (frac == "all") {
                    dat <- dat %>%
                        dplyr::group_by(!!!syms(unique(c(split_by, facet_by, columns_split_by)))) %>%
                        mutate(.frac = !!sym(".n") / sum(!!sym(".n"))) %>%
                        ungroup()
                } else {
                    dat <- dat %>% mutate(.frac = 1)  # not used
                }

                dat[, setdiff(group_by, g)] <- NA
                dat
            }))
        } else {
            object <- object %>%
                dplyr::group_by(!!!syms(unique(c(split_by, facet_by, group_by, columns_split_by, ident)))) %>%
                summarise(.n = n(), .groups = "drop")

            if (frac == "group") {
                object <- object %>%
                    dplyr::group_by(!!!syms(unique(c(split_by, facet_by, group_by, columns_split_by)))) %>%
                    mutate(.frac = !!sym(".n") / sum(!!sym(".n"))) %>%
                    ungroup()
            } else if (frac == "ident") {
                object <- object %>%
                    dplyr::group_by(!!!syms(unique(c(split_by, facet_by, columns_split_by, ident)))) %>%
                    mutate(.frac = !!sym(".n") / sum(!!sym(".n"))) %>%
                    ungroup()
            } else if (frac == "all") {
                object <- object %>%
                    dplyr::group_by(!!!syms(unique(c(split_by, facet_by, columns_split_by)))) %>%
                    mutate(.frac = !!sym(".n") / sum(!!sym(".n"))) %>%
                    ungroup()
            } else {
                object <- object %>% mutate(.frac = 1)  # not used
            }
        }
    }

    ylab <- ylab %||% paste0(ifelse(identical(frac, "none"), "Number", "Fraction"), " of cells")
    if (plot_type == "bar") {
        if (is.null(group_by)) {
            BarPlot(
                object, x = ident, y = ifelse(identical(frac, "none"), ".n", ".frac"),
                split_by = split_by, facet_by = facet_by, ylab = ylab, ...)
        } else {
            BarPlot(
                object,
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
            object,
            from = ifelse(swap, ident, group_by),
            to = ifelse(swap, group_by, ident),
            split_by = split_by, facet_by = facet_by, ...)
    } else if (plot_type == "pie") {
        if (!is.null(group_by)) {
            stop("Cannot create a pie plot with 'group_by'. You may want to use 'ring' plot instead.")
        }
        PieChart(
            object, x = ident, y = ifelse(identical(frac, "none"), ".n", ".frac"),
            split_by = split_by, facet_by = facet_by, ...)
    } else if (plot_type == "ring") {
        if (!identical(frac, "group")) {
            stop("'frac' must be 'group' for 'ring' plot.")
        }
        RingPlot(
            object,
            x = ifelse(swap, ident, group_by),
            y = ifelse(identical(frac, "none"), ".n", ".frac"),
            group_by = ifelse(swap, group_by, ident),
            split_by = split_by, facet_by = facet_by, ...)
    } else if (plot_type == "trend") {
        if (is.null(group_by)) {
            stop("Cannot create a trend plot without specifying 'group_by'.")
        }
        TrendPlot(
            object,
            x = ifelse(swap, ident, group_by),
            y = ifelse(identical(frac, "none"), ".n", ".frac"),
            group_by = ifelse(swap, group_by, ident),
            ylab = ylab, split_by = split_by, facet_by = facet_by, ...)
    } else if (plot_type == "area") {
        if (is.null(group_by)) {
            stop("Cannot create an area plot without specifying 'group_by'.")
        }
        AreaPlot(
            object,
            x = ifelse(swap, ident, group_by),
            y = ifelse(identical(frac, "none"), ".n", ".frac"),
            group_by = ifelse(swap, group_by, ident),
            ylab = ylab, split_by = split_by, facet_by = facet_by, ...)
    } else if (plot_type == "sankey") {
        if (is.null(group_by)) {
            stop("Cannot create a sankey plot without specifying 'group_by'.")
        }
        if (frac == "ident") {
            stop("Cannot calculate the fraction by 'ident' for 'sankey' plot.")
        }
        if (isTRUE(swap)) {
            stop("'swap = TRUE' is not supported for 'sankey' plot.")
        }
        SankeyPlot(
            object, x = group_by, links_fill_by = ident, xlab = "", ylab = "", flow = TRUE,
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
            object, rows = rows, cell_type = "pie", rows_name = rows_name, name = name,
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
        idents <- if (is.factor(object[[ident]])) levels(object[[ident]]) else unique(object[[ident]])
        object <- pivot_wider(object, id_cols = unique(c(split_by, group_by, columns_split_by)),
            names_from = ident, values_fill = 0,
            values_from = if (identical(frac, "none")) ".n" else ".frac")

        rows_name <- rows_name %||% ident
        name <- name %||% ifelse(identical(frac, "none"), "Number of cells", "Fraction of cells")

        Heatmap(
            object, rows = idents, rows_name = rows_name, name = name,
            columns_by = if (swap) columns_split_by else group_by,
            columns_split_by = if (swap) group_by else columns_split_by,
            split_by = split_by, ...)
    } else if (plot_type %in% c("violin", "box")) {
        if (is.null(group_by)) {
            stop("Cannot create a 'violin'/'box' plot without specifying 'group_by'.")
        }
        if (length(group_by) > 2) {
            stop("Cannot create a 'violin'/'box' plot with more than 2 'group_by'.")
        }
        fn <- ifelse(plot_type == "violin", ViolinPlot, BoxPlot)
        if (length(group_by) == 1) {
            fn(
                object,
                x = ifelse(swap, group_by, ident),
                y = ifelse(identical(frac, "none"), ".n", ".frac"),
                in_form = "long",
                ylab = ylab, split_by = split_by, facet_by = facet_by, ...)
        } else {
            fn(
                object,
                x = ifelse(swap, group_by[2], ident),
                y = ifelse(identical(frac, "none"), ".n", ".frac"),
                group_by = ifelse(swap, ident, group_by[2]),
                ylab = ylab, split_by = split_by, facet_by = facet_by, ...)
        }
    } else {  # radar/spider plot
        if (is.null(group_by)) {
            stop("Cannot create a '", plot_type, "' plot without specifying 'group_by'.")
        }
        if (plot_type == "radar") {
            RadarPlot(
                object,
                x = ifelse(swap, group_by, ident),
                y = ifelse(identical(frac, "none"), ".n", ".frac"),
                group_by = ifelse(swap, ident, group_by),
                split_by = split_by, facet_by = facet_by, ...)
        } else if (plot_type == "spider") {
            SpiderPlot(
                object,
                x = ifelse(swap, group_by, ident),
                y = ifelse(identical(frac, "none"), ".n", ".frac"),
                group_by = ifelse(swap, ident, group_by),
                split_by = split_by, facet_by = facet_by, ...)
        }
    }
}
