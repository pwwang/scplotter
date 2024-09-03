
#' Clustree plot
#'
#' @description This function generates a clustree plot from a data frame or a Seurat object.
#'
#' @rdname ClustreePlot
#'
#' @param x The data frame or Seurat object
#' @param prefix The prefix of cluster columns
#' @param edge_palette The color palette to use for edges. Default is "YlOrRd".
#' @param edge_palcolor Custom colors for edge_palette. Default is NULL.
#' @inheritParams validate_common_arguments
#' @param ... Additional arguments to be passed to the clustree::clustree() function
#' @return A ggplot object or a list with the plot and the height and width of the plot if guess_size is TRUE
#' @export
#' @examples
#' data(ifnb_sub)
#' ClustreePlot(ifnb_sub, prefix = "RNA_snn_res.")
ClustreePlot <- function(x, prefix, ...) UseMethod("ClustreePlot")

#' @rdname ClustreePlot
#' @export
ClustreePlot.default <- function(x, prefix, ...) {
    stop("ClustreePlot() is not implemented for objects of class ", class(x))
}

#' @inheritParams ClustreePlot
#' @importFrom dplyr %>% select starts_with
#' @importFrom ggplot2 scale_color_manual labs coord_cartesian element_line element_blank
#' @rdname ClustreePlot
#' @export
ClustreePlot.data.frame <- function(
    x, prefix,
    palette = "Paired", palcolor = NULL, edge_palette = "YlOrRd", edge_palcolor = NULL,
    aspect.ratio = 1, legend.position = "right", legend.direction = "vertical",
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
    theme = "theme_scp", theme_args = list(), res = 100, guess_size = FALSE, seed = 8525,
    ...) {
    set.seed(seed)
    # ! Unknown guide: edge_colourbar
    suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(ggraph))))

    x <- x %>% select(starts_with(prefix))
    x <- x[complete.cases(x), , drop = FALSE]

    if (ncol(x) == 0 || nrow(x) == 0) {
        stop("No data found with prefix '", prefix, "'")
    }
    resolutions <- substring(colnames(x), nchar(prefix) + 1)
    nres <- length(resolutions)
    if (all(startsWith(resolutions, "_"))) {
        prefix <- paste0(prefix, "_")
    } else if (all(startsWith(resolutions, "."))) {
        prefix <- paste0(prefix, ".")
    }
    resolutions <- as.numeric(substring(colnames(x), nchar(prefix) + 1))
    # make sure res.0.50 to be access by res.0.5
    colnames(x) <- paste0(prefix, resolutions)
    nres <- length(resolutions)
    max_clusters <- length(unique(x[[paste0(prefix, max(resolutions))]]))
    clustree_args <- list(...)
    clustree_args$x <- x
    clustree_args$prefix <- prefix
    clustree_args$node_alpha <- clustree_args$node_alpha %||% 0.85
    clustree_args$edge_width <- clustree_args$edge_width %||% 0.9
    clustree_args$show_axis <- clustree_args$show_axis %||% TRUE
    clustree_args$layout <- clustree_args$layout %||% "sugiyama"

    p <- suppressMessages(
        do.call(clustree::clustree, clustree_args) +
        scale_color_manual(
            values = palette_scp(n = nres, palette = palette, palcolor = palcolor),
            guide = "none") +
        scale_edge_color_gradientn(
            name = "count",
            n.breaks = 5,
            colors = palette_scp(palette = edge_palette, palcolor = edge_palcolor),
            na.value = "grey80") +
        do.call(theme, theme_args) +
        coord_cartesian(clip = "off") +
        ggplot2::theme(
            aspect.ratio = aspect.ratio,
            legend.position = legend.position,
            legend.direction = legend.direction,
            panel.grid.major.y = element_line(colour = "grey80", linetype = 2),
            axis.line.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            panel.border = element_blank()) +
        labs(title = title, subtitle = subtitle, x = xlab %||% "", y = ylab %||% sub("\\.|_$", "", prefix))
    )
    if (isFALSE(guess_size)) {
        p
    } else {
        list(
            plot = p,
            height = 1.15 * (nres) * res,
            width = ifelse(max_clusters > 20, 11, ifelse(max_clusters > 15, 9, 7)) * res
        )
    }
}

#' @rdname ClustreePlot
#' @export
ClustreePlot.Seurat <- function(x, prefix, ...) {
    ClustreePlot(x@meta.data, prefix, ...)
}
