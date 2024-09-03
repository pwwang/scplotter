#' Plot with a single value per group.
#'
#' @description These functions create plots with a single value per group.
#'   BarPlotBaseSingle: Bar plot.
#'   PieChartBaseSingle: Pie chart.
#'   LinePlotBaseSingle: Line plot.
#'
#' @rdname base_single_value_plots
#' @inheritParams base_plots_common_arguments
#' @param fill_by_x A logical value indicating whether to fill the bars by the x-axis values.
#'   If FALSE, the bars will be filled a single color (the first color in the palette).
#' @return A ggplot object.
#' @export
#' @importFrom rlang sym %||%
#' @importFrom gglogger ggplot
#' @importFrom ggplot2 aes geom_bar scale_fill_manual labs
#' @examples
#' data <- data.frame(
#'    x = c("A", "B", "C", "D"),
#'    y = c(10, 8, 16, 4)
#' )
#' BarPlotBaseSingle(data, x = "x", y = "y")
#' BarPlotBaseSingle(data, x = "x", y = "y", fill_by_x = F)
BarPlotBaseSingle <- function(
    data, x, y,
    theme = "theme_scp", theme_args = list(), palette = "Paired", palcolor = NULL,
    alpha = 1, x_text_angle = 0, aspect.ratio = 1,
    legend.position = "right", legend.direction = "vertical",
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL, keep_empty = FALSE,
    expand_x = 0, expand_y = 0, fill_by_x = TRUE, ...
) {
    if (length(expand_x) == 1) {
        expand_x <- c(expand_x, expand_x)
    }
    if (length(expand_y) == 1) {
        expand_y <- c(expand_y, expand_y)
    }
    if (isTRUE(fill_by_x)) {
        p <- ggplot(data, aes(x = !!sym(x), y = !!sym(y), fill = !!sym(x)))
        colors <- palette_scp(
            if (is.factor(data[[x]])) levels(data[[x]]) else unique(data[[x]]),
            palette = palette, palcolor = palcolor, keep_names = TRUE
        )
        guide = "legend"
    } else {
        p <- ggplot(data, aes(x = !!sym(x), y = !!sym(y), fill = "fill"))
        colors <- palette_scp("fill", palette = palette, palcolor = palcolor, keep_names = TRUE)
        guide = "none"
    }
    just <- calc_just(x_text_angle)

    p + geom_col(alpha = alpha) +
        scale_fill_manual(name = x, values = colors, guide = guide) +
        labs(title = title, subtitle = subtitle, x = xlab %||% x, y = ylab %||% y) +
        scale_x_discrete(drop = !keep_empty, expand = c(expand_x[1], expand_x[2])) +
        scale_y_continuous(expand = c(expand_y[1], expand_y[2])) +
        do.call(theme, theme_args) +
        ggplot2::theme(
            aspect.ratio = aspect.ratio,
            legend.position = legend.position,
            legend.direction = legend.direction,
            panel.grid.major = element_line(colour = "grey80", linetype = 2),
            axis.text.x = element_text(angle = x_text_angle, hjust = just$h, vjust = just$v)
        )
}

#' @rdname base_single_value_plots
#' @inheritParams base_plots_common_arguments
#' @param label Which column to use as the label. NULL means no label.
#' @export
#' @importFrom rlang sym
#' @importFrom gglogger ggplot
#' @importFrom dplyr lead if_else
#' @importFrom tidyr complete replace_na
#' @importFrom ggplot2 coord_polar geom_col
#' @examples
#' data <- data.frame(
#'    x = c("A", "B", "C", "D"),
#'    y = c(10, 8, 16, 4)
#' )
#' PieChartBaseSingle(data, x = "x", y = "y")
#' PieChartBaseSingle(data, x = "x", y = "y", label = NULL)
PieChartBaseSingle <- function(
    data, x, y, label = "y",
    theme = "theme_scp", theme_args = list(), palette = "Paired", palcolor = NULL,
    alpha = 1, aspect.ratio = 1,
    legend.position = "right", legend.direction = "vertical",
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL, keep_empty = FALSE, ...
) {
    base_size <- theme_args$base_size %||% 12
    text_size_scale <- base_size / 12
    if (!is.factor(data[[x]])) {
        data[[x]] <- factor(data[[x]], levels = unique(data[[x]]))
    }
    # if keep_empty is TRUE, fill the empty levels with 0
    if (keep_empty) {
        data <- data %>% complete(!!sym(x))
        data[[y]] <- data[[y]] %>% replace_na(0)
    }
    # order the data by the levels of x
    data <- data[order(data[[x]]), , drop = FALSE]

    pos_df <- data %>% mutate(
        csum = rev(cumsum(rev(!!sym(y)))),
        pos = !!sym(y) / 2 + lead(csum, 1),
        pos = if_else(is.na(pos), !!sym(y) / 2, pos)
    )
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

    p
}

#' @rdname base_single_value_plots
#' @inheritParams base_plots_common_arguments
#' @param fill_point_by_x A logical value indicating whether to color the points by the x-axis values.
#'   If FALSE, the lines will be colored a single color (the first color in the palette).
#' @param color_line_by_x A logical value indicating whether to color the lines by the x-axis values.
#'  If FALSE, the lines will be colored a single color (the first color in the palette).
#' @param line_type The type of line to draw.
#' @param line_width The width of the line.
#' @param line_alpha The alpha value of the line.
#' @param point_alpha The alpha value of the points.
#' @param point_size The size of the points.
#' @param add_bg A logical value indicating whether to add a background to the plot.
#' @param bg_palette The palette to use for the background.
#' @param bg_palcolor The color to use for the background.
#' @param bg_alpha The alpha value of the background.
#' @param add_errorbars A logical value indicating whether to add error bars to the plot.
#' @param errorbar_color The color to use for the error bars.
#'   If "line", the error bars will be colored the same as the lines.
#' @param errorbar_alpha The alpha value of the error bars.
#' @param errorbar_linewidth The line width of the error bars.
#' @param errorbar_width The width of the error bars.
#' @param errorbar_min The column in the data frame containing the lower bound of the error bars.
#' @param errorbar_max The column in the data frame containing the upper bound of the error bars.
#' @param errorbar_sd The column in the data frame containing the standard deviation of the error bars.
#'   If errorbar_min and errorbar_max are not provided, this column will be used to calculate the error bars.
#'   errorbar_min = y - errorbar_sd, errorbar_max = y + errorbar_sd.
#'   If errorbar_min and errorbar_max are provided, this column will be ignored.
#'
#' @export
#' @importFrom rlang sym
#' @importFrom gglogger ggplot
#' @importFrom ggplot2 geom_line scale_color_manual labs geom_rect geom_errorbar
#' @examples
#' data <- data.frame(
#'    x = c("A", "B", "C", "D"),
#'    y = c(10, 8, 16, 4)
#' )
#' LinePlotBaseSingle(data, x = "x", y = "y")
#' LinePlotBaseSingle(data, x = "x", y = "y", fill_point_by_x = FALSE)
#' LinePlotBaseSingle(data, x = "x", y = "y", add_bg = TRUE)
#'
#' data$sd <- c(1, 2, 3, 4)
#' LinePlotBaseSingle(data, x = "x", y = "y", add_errorbars = TRUE, errorbar_sd = "sd")
LinePlotBaseSingle <- function(
    data, x, y, fill_point_by_x = TRUE, color_line_by_x = TRUE,
    add_bg = FALSE, bg_palette = "Paired", bg_palcolor = NULL, bg_alpha = 0.2,
    add_errorbars = FALSE, errorbar_width = 0.1, errorbar_alpha = 1,
    errorbar_color = "grey30", errorbar_linewidth = .75, errorbar_min = NULL, errorbar_max = NULL, errorbar_sd = NULL,
    point_alpha = 1, point_size = 5,
    line_type = "solid", line_width = 1, line_alpha = .8,
    theme = "theme_scp", theme_args = list(), palette = "Paired", palcolor = NULL,
    x_text_angle = 0, aspect.ratio = 1,
    legend.position = "right", legend.direction = "vertical",
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL, keep_empty = FALSE, ...
) {
    if (!is.factor(data[[x]])) {
        data[[x]] <- factor(data[[x]], levels = unique(data[[x]]))
    }
    if (isTRUE(add_bg)) {
        bg_color <- palette_scp(
            levels(data[[x]]),
            palette = bg_palette, palcolor = bg_palcolor, keep_names = TRUE
        )
        bg_data <- na.omit(unique(data[, x, drop = FALSE]))
        bg_data$x <- as.numeric(bg_data[[x]])
        bg_data$xmin <- ifelse(bg_data$x == min(bg_data$x), -Inf, bg_data$x - 0.5)
        bg_data$xmax <- ifelse(bg_data$x == max(bg_data$x), Inf, bg_data$x + 0.5)
        bg_data$ymin <- -Inf
        bg_data$ymax <- Inf
        bg_data$fill <- bg_color[as.character(data$x)]
    }
    if (isTRUE(add_errorbars)) {
        if (is.null(errorbar_sd) && (is.null(errorbar_min) || is.null(errorbar_max))) {
            stop("If 'errorbar_min' and 'errorbar_max' are not provided, 'errorbar_sd' must be provided.")
        }
        if (is.null(errorbar_min) || is.null(errorbar_max)) {
            data$errorbar_min <- data[[y]] - data[[errorbar_sd]]
            data$errorbar_max <- data[[y]] + data[[errorbar_sd]]
            errorbar_min <- "errorbar_min"
            errorbar_max <- "errorbar_max"
        }
    }

    p <- ggplot(data, aes(x = !!sym(x), y = !!sym(y)))
    if (isTRUE(add_bg)) {
        p <- p + geom_rect(
            data = bg_data,
            xmin = bg_data$xmin, xmax = bg_data$xmax, ymin = bg_data$ymin, ymax = bg_data$ymax,
            fill = bg_data$fill, alpha = bg_alpha, inherit.aes = FALSE)
    }
    colors <- palette_scp(
        levels(data[[x]]),
        palette = palette, palcolor = palcolor, keep_names = TRUE
    )
    if (isTRUE(color_line_by_x)) {
        p <- p + geom_line(
            aes(color = !!sym(x), group = 1),
            alpha = line_alpha, linetype = line_type, linewidth = line_width) +
            scale_color_manual(name = x, values = colors, guide = "legend", drop = !keep_empty)
    } else {
        p <- p + geom_line(
            aes(group = 1), color = colors[[1]],
            alpha = line_alpha, linetype = line_type, linewidth = line_width)
    }
    if (isTRUE(add_errorbars)) {
        if (errorbar_color == "line" && isTRUE(color_line_by_x)) {
            p <- p + geom_errorbar(
                aes(ymin = !!sym(errorbar_min), ymax = !!sym(errorbar_max), color = !!sym(x)),
                alpha = errorbar_alpha, width = errorbar_width, linewidth = errorbar_linewidth)
        } else if (errorbar_color == "line") {
            p <- p + geom_errorbar(
                aes(ymin = !!sym(errorbar_min), ymax = !!sym(errorbar_max)), color = colors[[1]],
                alpha = errorbar_alpha, width = errorbar_width, linewidth = errorbar_linewidth)
        } else {
            p <- p + geom_errorbar(
                aes(ymin = !!sym(errorbar_min), ymax = !!sym(errorbar_max)), color = errorbar_color,
                alpha = errorbar_alpha, width = errorbar_width, linewidth = errorbar_linewidth)
        }
    }
    if (isTRUE(fill_point_by_x)) {
        p <- p + geom_point(
            aes(fill = !!sym(x)),
            color = "grey20", alpha = point_alpha, size = point_size, shape = 21) +
            scale_fill_manual(name = x, values = colors, guide = "legend", drop = !keep_empty)
    } else {
        p <- p + geom_point(
            fill = colors[[1]],
            color = "grey20", alpha = point_alpha, size = point_size, shape = 21)
    }
    just <- calc_just(x_text_angle)
    p <- p + scale_x_discrete(drop = !keep_empty) +
        labs(title = title, subtitle = subtitle, x = xlab %||% x, y = ylab %||% y) +
        do.call(theme, theme_args) +
        ggplot2::theme(
            aspect.ratio = aspect.ratio,
            legend.position = legend.position,
            legend.direction = legend.direction,
            panel.grid.major = element_line(colour = "grey80", linetype = 2),
            axis.text.x = element_text(angle = x_text_angle, hjust = just$h, vjust = just$v)
        )
    p
}
