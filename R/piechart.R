#' Pie chart
#'
#' @description Pie chart showing the composition of the data.
#' @inheritParams atomic_plots_common_arguments
#' @param label Which column to use as the label. NULL means no label.
#' @keywords internal
#' @rdname piechart
#' @importFrom rlang sym
#' @importFrom gglogger ggplot
#' @importFrom dplyr lead if_else mutate
#' @importFrom tidyr complete replace_na
#' @importFrom ggplot2 coord_polar geom_col
PieChartAtomic <- function(
    data, x, y, label = "y",
    theme = "theme_scp", theme_args = list(), palette = "Paired", palcolor = NULL,
    alpha = 1, facet = FALSE, facet_by = NULL, facet_scales = "fixed",
    aspect.ratio = 1, legend.position = "right", legend.direction = "vertical",
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL, keep_empty = FALSE,
    nrow = NULL, ncol = NULL, byrow = TRUE, ...
) {
    base_size <- theme_args$base_size %||% 12
    text_size_scale <- base_size / 12

    x <- check_columns(data, x, force_factor = TRUE)
    y <- check_columns(data, y)
    concated_facet_by <- check_columns(data, facet_by, force_factor = TRUE, allow_multi = TRUE, concat_multi = TRUE)
    # if keep_empty is TRUE, fill the empty levels with 0
    if (keep_empty) {
        if (isTRUE(facet) && !is.null(concated_facet_by)) {
            data <- do.call(rbind, lapply(split(data, data[[concated_facet_by]]), function(d) {
                d <- d %>% complete(!!sym(x))
                d[[y]] <- d[[y]] %>% replace_na(0)
                d
            }))
        } else {
            data <- data %>% complete(!!sym(x))
            data[[y]] <- data[[y]] %>% replace_na(0)
        }
    }
    # order the data by the levels of x
    data <- data[order(data[[x]]), , drop = FALSE]

    if (isTRUE(facet) && !is.null(concated_facet_by)) {
        pos_df <- do.call(rbind, lapply(split(data, data[[concated_facet_by]]), function(d) {
            d %>% mutate(
                csum = rev(cumsum(rev(!!sym(y)))),
                pos = !!sym(y) / 2 + lead(csum, 1),
                pos = if_else(is.na(pos), !!sym(y) / 2, pos)
            )
        }))
    } else {
        pos_df <- data %>% mutate(
            csum = rev(cumsum(rev(!!sym(y)))),
            pos = !!sym(y) / 2 + lead(csum, 1),
            pos = if_else(is.na(pos), !!sym(y) / 2, pos)
        )
    }
    if (!is.null(label)) {
        pos_df[[label]] <- data[[label]]
    }
    colors <- palette_scp(
        levels(data[[x]]),
        palette = palette, palcolor = palcolor, keep_names = TRUE
    )

    p <- ggplot(data, aes(x = "", y = !!sym(y), fill = !!sym(x))) +
        geom_col(width = 1, alpha = alpha) +
        scale_fill_manual(name = x, drop = !keep_empty, values = colors) +
        labs(title = title, subtitle = subtitle, x = xlab, y = ylab) +
        coord_polar(theta = "y") +
        do.call(theme, theme_args) +
        ggplot2::theme(
            aspect.ratio = aspect.ratio,
            legend.position = legend.position,
            legend.direction = legend.direction,
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            panel.border = element_blank()
        )

    if (!is.null(label)) {
        p <- p + geom_label_repel(
            data = pos_df,
            aes(y = pos, label = !!sym(label)),
            nudge_x = 0.1,
            color = "grey20",
            fill = "#fcfcfc",
            size = text_size_scale * 3
        )
    }

    height <- 4.5
    width <- 4.5
    if (legend.position %in% c("right", "left")) {
        width <- width + 1
    } else if (legend.direction == "horizontal") {
        height <- height + 1
    } else {
        width <- width + 2
    }

    attr(p, "height") <- height
    attr(p, "width") <- width
    facet_plot(p, facet, facet_scales, facet_by, nrow, ncol, byrow)
}

#' @export
#' @rdname piechart
#' @inheritParams PieChartAtomic
#' @inheritParams validate_common_arguments
#' @return A ggplot object.
PieChart <- function(
    data, x, y, label = "y", split_by = NULL, split_by_sep = "_",
    facet = FALSE, facet_scales = "fixed",
    theme = "theme_scp", theme_args = list(), palette = "Paired", palcolor = NULL,
    alpha = 1, aspect.ratio = 1,
    legend.position = "right", legend.direction = "vertical",
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL, keep_empty = FALSE,
    combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, seed = NULL, ...
) {
    validate_common_arguments(split_by, seed, facet)

    x <- check_columns(data, x, force_factor = TRUE)
    y <- check_columns(data, y)
    split_by <- check_columns(data, split_by, force_factor = TRUE, allow_multi = TRUE, concat_multi = TRUE, concat_sep = split_by_sep)

    if (is.null(split_by) || isTRUE(facet)) {
        datas <- list(data)
    } else {
        datas <- split(data, data[[split_by]])
        datas <- datas[levels(data[[split_by]])]
    }

    plots <- lapply(
        datas, PieChartAtomic,
        x = x, y = y, label = label, facet_by = split_by,
        facet = facet, facet_scales = facet_scales,
        theme = theme, theme_args = theme_args, palette = palette, palcolor = palcolor,
        alpha = alpha, aspect.ratio = aspect.ratio,
        legend.position = legend.position, legend.direction = legend.direction,
        title = title, subtitle = subtitle, xlab = xlab, ylab = ylab, keep_empty = keep_empty
    )

    combine_plots(plots, combine = combine, nrow = nrow, ncol = ncol, byrow = byrow)
}
