#' Bar plot
#'
#' @description Provides a base function to create bar plots.
#'
#' @rdname barplot
#' @inheritParams base_plots_common_arguments
#' @param fill_by_x A logical value indicating whether to fill the bars by the x-axis values.
#'   If FALSE, the bars will be filled a single color (the first color in the palette).
#' @param width A numeric value specifying the width of the bars.
#' @return A ggplot object.
#' @keywords internal
#' @importFrom rlang sym %||%
#' @importFrom gglogger ggplot
#' @importFrom ggplot2 aes geom_bar scale_fill_manual labs
# ' @examples
# ' data <- data.frame(
# '    x = c("A", "B", "C", "D"),
# '    y = c(10, 8, 16, 4)
# ' )
# ' BarPlotSingle(data, x = "x", y = "y")
# ' BarPlotSingle(data, x = "x", y = "y", fill_by_x = F)
BarPlotSingle <- function(
    data, x, y,
    theme = "theme_scp", theme_args = list(), palette = "Paired", palcolor = NULL,
    alpha = 1, x_text_angle = 0, aspect.ratio = 1,
    legend.position = "right", legend.direction = "vertical",
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL, keep_empty = FALSE,
    expand = 0, fill_by_x = TRUE, width = 0.9, ...
) {
    expand <- expand %||% 0
    expand <- expand_expansion(expand)

    x <- check_columns(data, x, force_factor = TRUE)
    y <- check_columns(data, y)

    if (isTRUE(fill_by_x)) {
        p <- ggplot(data, aes(x = !!sym(x), y = !!sym(y), fill = !!sym(x)))
        colors <- palette_scp(
            levels(data[[x]]),
            palette = palette, palcolor = palcolor, keep_names = TRUE
        )
        guide = "legend"
    } else {
        p <- ggplot(data, aes(x = !!sym(x), y = !!sym(y), fill = "fill"))
        colors <- palette_scp("fill", palette = palette, palcolor = palcolor, keep_names = TRUE)
        guide = "none"
    }
    just <- calc_just(x_text_angle)

    p <- p + geom_col(alpha = alpha, width = width) +
        scale_fill_manual(name = x, values = colors, guide = guide) +
        labs(title = title, subtitle = subtitle, x = xlab %||% x, y = ylab %||% y) +
        scale_x_discrete(drop = !keep_empty, expand = expand$x) +
        scale_y_continuous(expand = expand$y) +
        do.call(theme, theme_args) +
        ggplot2::theme(
            aspect.ratio = aspect.ratio,
            legend.position = legend.position,
            legend.direction = legend.direction,
            panel.grid.major = element_line(colour = "grey80", linetype = 2),
            axis.text.x = element_text(angle = x_text_angle, hjust = just$h, vjust = just$v)
        )

    height <- 4.5
    width <- .5 + length(levels(data[[x]])) * .8
    if (guide == "legend") {
        if (legend.position %in% c("right", "left")) {
            width <- width + 1
        } else if (legend.direction == "horizontal") {
            height <- height + 1
        } else {
            width <- width + 2
        }
    }
    attr(p, "height") <- height
    attr(p, "width") <- width
    p
}

#' @rdname barplot
#' @inheritParams BarPlotSingle
#' @param position A character string indicating the position of the bars.
#'  If "auto", the position will be "stack" if group_by has more than 5 levels, otherwise "dodge".
#'  "fill" is also a valid option. Only works when group_by is not NULL.
#' @param position_dodge_preserve Should dodging preserve the "total" width of all elements at a position, or the width of a "single" element?
#' @param add_bg A logical value indicating whether to add a background to the plot.
#' @param bg_palette A character string indicating the palette to use for the background.
#' @param bg_palcolor A character string indicating the color to use for the background.
#' @param bg_alpha A numeric value indicating the alpha of the background.
#' @return A ggplot object.
#' @keywords internal
#' @importFrom rlang sym %||%
#' @importFrom gglogger ggplot
#' @importFrom ggplot2 aes geom_bar scale_fill_manual labs position_dodge2
BarPlotGrouped <- function(
    data, x, y, group_by, group_by_sep = "_",
    theme = "theme_scp", theme_args = list(), palette = "Paired", palcolor = NULL,
    add_bg = FALSE, bg_palette = "Paired", bg_palcolor = NULL, bg_alpha = 0.2,
    alpha = 1, x_text_angle = 0, aspect.ratio = 1,
    position = "auto", position_dodge_preserve = "total",
    legend.position = "right", legend.direction = "vertical",
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL, keep_empty = FALSE,
    expand = c(0.05, 0.2, 0, 0.2), width = 0.9, ...
) {
    expand <- expand %||% c(0.05, 0.2, 0, 0.2)
    expand <- expand_expansion(expand)
    x <- check_columns(data, x, force_factor = TRUE)
    y <- check_columns(data, y)
    group_by <- check_columns(data, group_by, force_factor = TRUE, allow_multi = TRUE, concat_multi = TRUE, concat_sep = group_by_sep)

    p <- ggplot(data, aes(x = !!sym(x), y = !!sym(y), fill = !!sym(group_by)))
    if (isTRUE(add_bg)) {
        p <- p + bg_layer(data, x, bg_palette, bg_palcolor, bg_alpha, keep_empty)
    }

    colors <- palette_scp(
        levels(data[[group_by]]),
        palette = palette, palcolor = palcolor, keep_names = TRUE
    )
    just <- calc_just(x_text_angle)
    if (position == "auto") {
        position <- if (length(colors) <= 5) position_dodge2(preserve = position_dodge_preserve) else "stack"
    } else if (position == "dodge") {
        position <- position_dodge2(preserve = position_dodge_preserve)
    }

    p <- p + geom_col(alpha = alpha, position = position, width = width) +
        scale_fill_manual(name = group_by, values = colors, drop = !keep_empty) +
        labs(title = title, subtitle = subtitle, x = xlab %||% x, y = ylab %||% y) +
        scale_x_discrete(drop = !keep_empty, expand = expand$x) +
        scale_y_continuous(expand = expand$y) +
        do.call(theme, theme_args) +
        ggplot2::theme(
            aspect.ratio = aspect.ratio,
            legend.position = legend.position,
            legend.direction = legend.direction,
            panel.grid.major = element_line(colour = "grey80", linetype = 2),
            axis.text.x = element_text(angle = x_text_angle, hjust = just$h, vjust = just$v)
        )

    height <- 4.5
    width <- .5 + length(levels(data[[x]])) * length(unique(data[[group_by]])) * .5
    if (legend.position %in% c("right", "left")) {
        width <- width + 1
    } else if (legend.direction == "horizontal") {
        height <- height + 1
    } else {
        width <- width + 2
    }
    attr(p, "height") <- height
    attr(p, "width") <- width
    p
}

#' @rdname barplot
#' @inheritParams BarPlotGrouped
#' @param fill_by_x_if_no_group A logical value indicating whether to fill the bars by the x-axis values if there is no group_by.
#' @param facet A logical value indicating whether to facet the plot.
#' @param facet_by A column name to facet the plot by.
#' @param facet_scales A character string indicating the scales of the facets.
#' @param nrow An integer value indicating the number of rows using facet_wrap.
#' @param ncol An integer value indicating the number of columns using facet_wrap.
#' @param byrow A logical value indicating whether to fill the facets by row.
#' @return A ggplot object.
#' @keywords internal
BarPlotAtomic <- function(
    data, x, y, group_by = NULL, fill_by_x_if_no_group = TRUE,
    theme = "theme_scp", theme_args = list(), palette = "Paired", palcolor = NULL,
    alpha = 1, x_text_angle = 0, aspect.ratio = 1,
    position = "auto", position_dodge_preserve = "total",
    legend.position = "right", legend.direction = "vertical",
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL, keep_empty = FALSE,
    expand = NULL, width = 0.9, facet = TRUE, facet_by = NULL, facet_scales = "fixed",
    nrow = NULL, ncol = NULL, byrow = TRUE, ...
) {
    if (is.null(group_by)) {
        p <- BarPlotSingle(
            data, x, y,
            theme = theme, theme_args = theme_args, palette = palette, palcolor = palcolor,
            alpha = alpha, x_text_angle = x_text_angle, aspect.ratio = aspect.ratio,
            legend.position = legend.position, legend.direction = legend.direction,
            title = title, subtitle = subtitle, xlab = xlab, ylab = ylab, keep_empty = keep_empty,
            expand = expand, fill_by_x = fill_by_x_if_no_group, width = width, ...
        )
    } else {
        p <- BarPlotGrouped(
            data, x, y, group_by,
            theme = theme, theme_args = theme_args, palette = palette, palcolor = palcolor,
            alpha = alpha, x_text_angle = x_text_angle, aspect.ratio = aspect.ratio,
            position = position, position_dodge_preserve = position_dodge_preserve,
            legend.position = legend.position, legend.direction = legend.direction,
            title = title, subtitle = subtitle, xlab = xlab, ylab = ylab, keep_empty = keep_empty,
            expand = expand, width = width, ...
        )
    }

    facet_plot(p, facet, facet_scales, facet_by, nrow, ncol, byrow)
}

#' @rdname barplot
#' @inheritParams BarPlotAtomic
#' @inheritParams validate_common_arguments
#' @return A ggplot object
#' @export
BarPlot <- function(
    data, x, y, group_by = NULL, group_by_sep = "_", split_by = NULL, split_by_sep = "_",
    facet = FALSE, facet_scales = "fixed", fill_by_x_if_no_group = TRUE,
    add_bg = FALSE, bg_palette = "Paired", bg_palcolor = NULL, bg_alpha = 0.2,
    theme = "theme_scp", theme_args = list(), palette = "Paired", palcolor = NULL,
    alpha = 1, x_text_angle = 0, aspect.ratio = 1,
    position = "auto", position_dodge_preserve = "total",
    legend.position = "right", legend.direction = "vertical",
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL, keep_empty = FALSE,
    expand = NULL, width = 0.9, combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, seed = 8525, ...
) {
    validate_common_arguments(split_by, seed, facet)

    x <- check_columns(data, x, force_factor = TRUE)
    y <- check_columns(data, y)
    group_by <- check_columns(data, group_by, force_factor = TRUE, allow_multi = TRUE, concat_multi = TRUE, concat_sep = group_by_sep)
    split_by <- check_columns(data, split_by, force_factor = TRUE, allow_multi = TRUE, concat_multi = TRUE, concat_sep = split_by_sep)

    if (is.null(split_by) || isTRUE(facet)) {
        datas <- list(data)
    } else {
        datas <- split(data, data[[split_by]])
        datas <- datas[levels(data[[split_by]])]
    }

    plots <- lapply(
        datas, BarPlotAtomic,
        x = x, y = y, group_by = group_by, fill_by_x_if_no_group = fill_by_x_if_no_group,
        theme = theme, theme_args = theme_args, palette = palette, palcolor = palcolor,
        add_bg = add_bg, bg_palette = bg_palette, bg_palcolor = bg_palcolor, bg_alpha = bg_alpha,
        alpha = alpha, x_text_angle = x_text_angle, aspect.ratio = aspect.ratio,
        position = position, position_dodge_preserve = position_dodge_preserve,
        legend.position = legend.position, legend.direction = legend.direction,
        title = title, subtitle = subtitle, xlab = xlab, ylab = ylab, keep_empty = keep_empty,
        expand = expand, width = width, facet = facet, facet_scales = facet_scales,
        facet_by = split_by, nrow = nrow, ncol = ncol, byrow = byrow, ...
    )

    combine_plots(plots, combine = combine, nrow = nrow, ncol = ncol, byrow = byrow)
}
