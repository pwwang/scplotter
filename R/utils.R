#' Adjust_network_layout
#'
#' @keywords internal
#' @importFrom igraph degree neighbors V
adjust_network_layout <- function(graph, layout, width, height = 2, scale = 100, iter = 100) {
    w <- width / 2
    layout[, 1] <- layout[, 1] / diff(range(layout[, 1])) * scale
    layout[, 2] <- layout[, 2] / diff(range(layout[, 2])) * scale

    adjusted <- c()
    # for (i in seq_len(iter)) {
    for (v in order(degree(graph), decreasing = TRUE)) {
        adjusted <- c(adjusted, v)
        neighbors <- as.numeric(neighbors(graph, V(graph)[v]))
        neighbors <- setdiff(neighbors, adjusted)
        x <- layout[v, 1]
        y <- layout[v, 2]
        r <- w[v]
        for (neighbor in neighbors) {
            nx <- layout[neighbor, 1]
            ny <- layout[neighbor, 2]
            ndist <- sqrt((nx - x)^2 + (ny - y)^2)
            nr <- w[neighbor]
            expect <- r + nr
            if (ndist < expect) {
                dx <- (x - nx) * (expect - ndist) / ndist
                dy <- (y - ny) * (expect - ndist) / ndist
                layout[neighbor, 1] <- nx - dx
                layout[neighbor, 2] <- ny - dy
                adjusted <- c(adjusted, neighbor)
            }
        }
    }
    # }

    for (i in seq_len(iter)) {
        dist_matrix <- as.matrix(dist(layout))
        nearest_neighbors <- apply(dist_matrix, 2, function(x) which(x == min(x[x > 0])), simplify = FALSE)
        # nearest_neighbors <- apply(dist_matrix, 2, function(x) {
        #   head(order(x), 3)[-1]
        # }, simplify = FALSE)
        for (v in sample(seq_len(nrow(layout)))) {
            neighbors <- unique(nearest_neighbors[[v]])
            x <- layout[v, 1]
            y <- layout[v, 2]
            r <- w[v]
            for (neighbor in neighbors) {
                nx <- layout[neighbor, 1]
                ny <- layout[neighbor, 2]
                nr <- w[neighbor]
                if (abs(nx - x) < (r + nr) && abs(ny - y) < height) {
                    dx <- r + nr - (nx - x)
                    dy <- height - (ny - y)
                    if (sample(c(1, 0), 1) == 1) {
                        dx <- 0
                    } else {
                        dy <- 0
                    }
                    layout[neighbor, 1] <- nx - dx
                    layout[neighbor, 2] <- ny - dy
                }
            }
        }
    }
    return(layout)
}

#' Convert RGBA to RGB
#' @keywords internal
rgba_to_rgb <- function(RGBA, BackGround = c(1, 1, 1)) {
    A <- RGBA[[length(RGBA)]]
    RGB <- RGBA[[-length(RGBA)]] * A + BackGround * (1 - A)
    return(RGB)
}

#' Blend two colors
#' @keywords internal
blend_to_color <- function(C1, C2, mode = "blend") {
    c1 <- C1[[1]]
    c1a <- C1[[2]]
    c2 <- C2[[1]]
    c2a <- C2[[2]]
    A <- 1 - (1 - c1a) * (1 - c2a)
    if (A < 1.0e-6) {
        return(list(c(0, 0, 0), 1))
    }
    if (mode == "blend") {
        out <- (c1 * c1a + c2 * c2a * (1 - c1a)) / A
        A <- 1
    }
    if (mode == "average") {
        out <- (c1 + c2) / 2
        out[out > 1] <- 1
    }
    if (mode == "screen") {
        out <- 1 - (1 - c1) * (1 - c2)
    }
    if (mode == "multiply") {
        out <- c1 * c2
    }
    return(list(out, A))
}

#' Blend a list of colors
#' @keywords internal
blend_rgblist <- function(Clist, mode = "blend", RGB_BackGround = c(1, 1, 1)) {
    N <- length(Clist)
    ClistUse <- Clist
    while (N != 1) {
        temp <- ClistUse
        ClistUse <- list()
        for (C in temp[1:(length(temp) - 1)]) {
            c1 <- C[[1]]
            a1 <- C[[2]]
            c2 <- temp[[length(temp)]][[1]]
            a2 <- temp[[length(temp)]][[2]]
            ClistUse <- append(ClistUse, list(blend_to_color(C1 = list(c1, a1 * (1 - 1 / N)), C2 = list(c2, a2 * 1 / N), mode = mode)))
        }
        N <- length(ClistUse)
    }
    Result <- list(ClistUse[[1]][[1]], ClistUse[[1]][[2]])
    Result <- rgba_to_rgb(Result, BackGround = RGB_BackGround)
    return(Result)
}

#' Blend colors
#'
#' This function blends a list of colors using the specified blend mode.
#'
#' @param colors Color vectors.
#' @param mode Blend mode. One of "blend", "average", "screen", or "multiply".
#'
#' @examples
#' blend <- c("red", "green", blend_colors(c("red", "green"), mode = "blend"))
#' average <- c("red", "green", blend_colors(c("red", "green"), mode = "average"))
#' screen <- c("red", "green", blend_colors(c("red", "green"), mode = "screen"))
#' multiply <- c("red", "green", blend_colors(c("red", "green"), mode = "multiply"))
#' show_palettes(list("blend" = blend, "average" = average, "screen" = screen, "multiply" = multiply))
#'
#' @keywords internal
#' @return The blended color.
blend_colors <- function(colors, mode = c("blend", "average", "screen", "multiply")) {
    mode <- match.arg(mode)
    colors <- colors[!is.na(colors)]
    if (length(colors) == 0) {
        return(NA)
    }
    if (length(colors) == 1) {
        return(colors)
    }
    rgb <- as.list(as.data.frame(col2rgb(colors) / 255))
    Clist <- lapply(rgb, function(x) {
        list(x, 1)
    })
    blend_color <- blend_rgblist(Clist, mode = mode)
    blend_color <- rgb(blend_color[1], blend_color[2], blend_color[3])
    return(blend_color)
}

#' Calculate hjust and vjust based on angle
#' @keywords internal
calc_just <- function(angle) {
    angle <- angle %% 360
    if (angle < 0) {
        angle <- angle + 360
    }
    if (angle < 10) {
        h <- 0.5
        v <- 1
    } else if (angle < 90) {
        h <- 1
        v <- 1
    } else if (angle < 135) {
        h <- 1
        v <- 0.5
    } else if (angle < 180) {
        h <- 1
        v <- 0.5
    } else if (angle < 225) {
        h <- 0
        v <- 0
    } else if (angle < 270) {
        h <- 0
        v <- 0
    } else if (angle < 315) {
        h <- 0
        v <- 0.5
    } else if (angle < 360) {
        h <- 0
        v <- 1
    } else {
        h <- 0.5
        v <- 1
    }
    list(h = h, v = v)
}

#' Common arguments for base plots
#'
#' @name base_plots_common_arguments
#' @param data A data frame.
#' @param x A character string specifying the column name of the data frame to plot for the x-axis.
#' @param y A character string specifying the column name of the data frame to plot for the y-axis.
#' @param group_by Columns to group the data for plotting
#'   For those plotting functions that do not support multiple groups,
#'   They will be concatenated into one column, using \code{\link{group_by_sep}} as the separator
#' @param theme A character string or a theme class (i.e. ggplot2::theme_classic) specifying the theme to use.
#'   Default is "theme_scp".
#' @param theme_args A list of arguments to pass to the theme function.
#' @param palette A character string specifying the palette to use.
#' @param palcolor A character string specifying the color to use in the palette.
#' @param alpha A numeric value specifying the transparency of the plot.
#' @param x_text_angle A numeric value specifying the angle of the x-axis text.
#' @param aspect.ratio A numeric value specifying the aspect ratio of the plot.
#' @param title A character string specifying the title of the plot.
#' @param subtitle A character string specifying the subtitle of the plot.
#' @param xlab A character string specifying the x-axis label.
#' @param ylab A character string specifying the y-axis label.
#' @param legend.position A character string specifying the position of the legend.
#' @param legend.direction A character string specifying the direction of the legend.
#' @param keep_empty A logical value indicating whether to keep empty groups.
#'   If FALSE, empty groups will be removed.
#' @param expand The values to expand the x and y axes. It is like CSS padding.
#'   When a single value is provided, it is used for both axes on both sides.
#'   When two values are provided, the first value is used for the top/bottom side and the second value is used for the left/right side.
#'   When three values are provided, the first value is used for the top side, the second value is used for the left/right side, and the third value is used for the bottom side.
#'   When four values are provided, the values are used for the top, right, bottom, and left sides, respectively.
#'   The values will be applied as 'mult' to the 'expansion' function.
#'   See also \url{https://ggplot2.tidyverse.org/reference/expansion.html}
#' @param ... Additional arguments.
#' @keywords internal
NULL

#' Common arguments for atomic plots
#'
#' @name atomic_plots_common_arguments
#' @inheritParams base_plots_common_arguments
#' @param facet A logical value indicating whether to facet the plot.
#' @param facet_by A character string specifying the column name of the data frame to facet the plot.
#' @param facet_scales A character string specifying the scales of the facets.
#' @param nrow A numeric value specifying the number of rows in the facet.
#' @param ncol A numeric value specifying the number of columns in the facet.
#' @param byrow A logical value indicating whether to fill the plots by row.
#' @keywords internal
NULL

#' Common params for plotting functions and validate the ones related to plotting
#'
#' @keywords internal
#' @inheritParams base_plots_common_arguments
#' @param plot_type The type of plot.
#' @param split_by The column(s) to split data by and plot separately or facet by.
#'   When `facet` is TRUE, the columns specified here will be used to facet the plot
#'   When `facet` is FALSE, the columns specified here will be used to split the data
#'   and generate multiple plots and combine them into one
#'   If multiple columns are specified
#'   When `facet` is TRUE, up to 2 columns are allowed. If one column,
#'       `ggplot2::facet_wrap` is used, the number of rows and columns is determined
#'       by `nrow` and `ncol`
#'       If two columns, `ggplot2::facet_grid` is used. The first column is used for rows
#'      and the second column is used for columns
#'   When `facet` is FALSE, a warning will be issued and the columns will be
#'       concatenated into one column, using `split_by_sep` as the separator
#' @param split_by_sep The separator for multiple split_by columns. See \code{\link{split_by}}
#' @param group_by_sep The separator for multiple group_by columns. See \code{\link{group_by}}
#' @param facet Whether to facet the plot. Default is FALSE.
#'   Otherwise, the data will be split by \code{\link{split_by}} and generate multiple plots
#'   and combine them into one using \code{patchwork::wrap_plots}
#' @param facet_scales Whether to scale the axes of facets. Default is "fixed"
#'   Other options are "free", "free_x", "free_y". See \code{\link{ggplot2::facet_wrap}}
#' @param combine Whether to combine the plots into one when facet is FALSE. Default is TRUE.
#' @param nrow The number of rows in facet_wrap/wrap_plots. Default is NULL.
#' @param ncol The number of columns in facet_wrap/wrap_plots. Default is NULL.
#' @param byrow Whether to fill the plots by row. Default is TRUE.
#' @param seed The random seed to use. Default is 8525.
#' @param ... Additional arguments for future expansion.
validate_common_arguments <- function(
    split_by,
    seed,
    facet,
    plot_type = NULL,
    split_by_sep = "_",
    group_by = NULL,
    group_by_sep = "_",
    facet_scales = "fixed",
    theme = "theme_scp",
    theme_args = list(),
    palette = NULL,
    palcolor = NULL,
    keep_empty = FALSE,
    alpha = 1,
    x_text_angle = 0,
    aspect.ratio = 1,
    legend.position = "right",
    legend.direction = "vertical",
    title = NULL,
    subtitle = NULL,
    xlab = NULL,
    ylab = NULL,
    combine = TRUE,
    nrow = NULL,
    ncol = NULL,
    byrow = TRUE,
    ...) {
    set.seed(seed)

    if (length(split_by) > 2 && isTRUE(facet)) {
        stop("When 'facet' is TRUE, up to 2 columns are allowed in 'split_by'")
    }

    if (is.null(split_by) && isTRUE(facet)) {
        stop("When 'facet' is TRUE, 'split_by' must be specified")
    }
}

#' Facetting a plot
#' @keywords internal
#' @param plot The plot to facet or a list list(plot, height, width) if guess_size is TRUE
#' @param facet Whether to facet the plot
#' @param facet_by The column(s) to split data by and plot separately or facet by
#' @param facet_scales Whether to scale the axes of facets.
#' @param nrow The number of rows in facet_wrap
#' @param ncol The number of columns in facet_wrap
#' @param byrow Whether to fill the plots by row
#' @param recalc_size Whether to re-calculate the size of the plot
#' @return The faceted plot. If guess_size is TRUE, attr(p, "height") and attr(p, "width") will be set
#' @importFrom ggplot2 facet_wrap facet_grid ggplot_build
#' @keywords internal
facet_plot <- function(plot, facet, facet_scales, facet_by, nrow, ncol, byrow, recalc_size = TRUE) {
    if (isFALSE(facet) || is.null(facet_by)) {
        return(plot)
    }

    if (recalc_size) {
        p <- facet_plot(plot, facet, facet_scales, facet_by, nrow, ncol, byrow, recalc_size = FALSE)
        d <- wrap_dims(length(unique(ggplot_build(p)$data[[1]]$PANEL)))
        attr(p, "height") <- d[1] * attr(plot, "height")
        attr(p, "width") <- d[2] * attr(plot, "width")
        return(p)
    }

    if (length(facet_by) == 1) {
        plot <- plot + ggplot2::facet_wrap(facets = facet_by, scales = facet_scales, nrow = nrow, ncol = ncol, dir = if (byrow) "h" else "v")
    } else {
        plot <- plot + ggplot2::facet_grid(rows = facet_by[[1]], cols = facet_by[[2]], scales = facet_scales)
    }

    return(plot)
}

#' Combine plots into one
#' @keywords internal
#' @param plots A list of plots
#' @param combine Whether to combine the plots into one
#' @param nrow The number of rows in the combined plot
#' @param ncol The number of columns in the combined plot
#' @param byrow Whether to fill the plots by row
#' @param recalc_size Whether to re-calculate the size of the combined plot
#' @return The faceted plot. If guess_size is TRUE, attr(p, "height") and attr(p, "width") will be set
#' @importFrom patchwork wrap_plots
#' @importFrom rlang %||%
#' @importFrom ggplot2 wrap_dims
combine_plots <- function(plots, combine, nrow, ncol, byrow, recalc_size = TRUE) {
    if (isFALSE(combine)) {
        return(plots)
    }

    if (recalc_size) {
        d <- wrap_dims(length(plots))
        nrow <- nrow %||% d[1]
        ncol <- ncol %||% d[2]
        p <- combine_plots(plots, TRUE, nrow, ncol, byrow, recalc_size = FALSE)
        attr(p, "height") <- nrow * max(sapply(plots, function(x) attr(x, "height")))
        attr(p, "width") <- ncol * max(sapply(plots, function(x) attr(x, "width")))
        return(p)
    }

    if (length(plots) == 1) {
        return(plots[[1]])
    }

    wrap_plots(plots, nrow = nrow, ncol = ncol, byrow = byrow)
}

#' @keywords internal
#' @importFrom grid is.grob grobWidth grobHeight
#' @importFrom gtable is.gtable gtable_add_rows gtable_add_cols gtable_add_grob
add_grob <- function(gtable, grob, position = c("top", "bottom", "left", "right", "none"), space = NULL, clip = "on") {
    position <- match.arg(position)
    if (position == "none" || is.null(grob)) {
        return(gtable)
    }

    if (is.null(space)) {
        if (is.gtable(grob)) {
            if (position %in% c("top", "bottom")) {
                space <- sum(grob$heights)
            } else {
                space <- sum(grob$widths)
            }
        } else if (is.grob(grob)) {
            if (position %in% c("top", "bottom")) {
                space <- grobHeight(grob)
            } else {
                space <- grobWidth(grob)
            }
        }
    }

    if (position == "top") {
        gtable <- gtable_add_rows(gtable, space, 0)
        gtable <- gtable_add_grob(gtable, grob, t = 1, l = mean(gtable$layout[grepl(pattern = "panel", x = gtable$layout$name), "l"]), clip = clip)
    }
    if (position == "bottom") {
        gtable <- gtable_add_rows(gtable, space, -1)
        gtable <- gtable_add_grob(gtable, grob, t = dim(gtable)[1], l = mean(gtable$layout[grepl(pattern = "panel", x = gtable$layout$name), "l"]), clip = clip)
    }
    if (position == "left") {
        gtable <- gtable_add_cols(gtable, space, 0)
        gtable <- gtable_add_grob(gtable, grob, t = mean(gtable$layout[grep("panel", gtable$layout$name), "t"]), l = 1, clip = clip)
    }
    if (position == "right") {
        gtable <- gtable_add_cols(gtable, space, -1)
        gtable <- gtable_add_grob(gtable, grob, t = mean(gtable$layout[grep("panel", gtable$layout$name), "t"]), l = dim(gtable)[2], clip = clip)
    }
    return(gtable)
}

#' Expand the plot area with CSS-like padding
#' @keywords internal
#' @param expand A numeric vector of length 1, 2, 3, or 4
#' @return A list with x and y values for expand
expand_expansion <- function(expand) {
    if (length(expand) == 1) {
        expand <- rep(expand, 4)
    } else if (length(expand) == 2) {
        expand <- c(expand[1], expand[2], expand[1], expand[2])
    } else if (length(expand) == 3) {
        expand <- c(expand[1], expand[2], expand[3], expand[2])
    } else if (length(expand) > 4) {
        stop("'expand' must have 1, 2, 3, or 4 values")
    }

    return(list(
        x = c(expand[4], 0, expand[2], 0),
        y = c(expand[3], 0, expand[1], 0)
    ))
}

#' Check the columns if columns found in the data
#' @keywords internal
#' @param df A data frame
#' @param columns A character vector of column names
#' @param force_factor Whether to force the columns to be factors
#' @param allow_multi Whether to allow multiple columns
#' @param concat_multi Whether to concatenate multiple columns
#' @param concat_sep The separator to use for concatenation
#' @return A character string of the valid column
#' @importFrom tidyr unite
#' @importFrom rlang syms
check_columns <- function(
    df,
    columns,
    force_factor = FALSE,
    allow_multi = FALSE,
    concat_multi = FALSE,
    concat_sep = "_"
) {
    if (is.null(columns)) {
        return(NULL)
    }
    param_name = deparse(substitute(columns))
    if (isFALSE(allow_multi)) {
        if (length(columns) > 1) {
            stop(paste0("Only one column is allowed in '", param_name, "'"))
        }
        if (!columns %in% colnames(df)) {
            stop(paste0("'", columns, "' is not in the data."))
        }
    } else {
        notfound <- setdiff(columns, colnames(df))
        if (length(notfound) > 0) {
            stop(paste0("'", paste0(notfound, collapse = ", "), "' is not in the data."))
        }
        if (isTRUE(concat_multi) && length(columns) > 1) {
            warning(
                paste0("Multiple columns are provided in '", param_name, "'. They will be concatenated into one column."),
                immediate. = TRUE
            )
            new_col <- paste(columns, collapse = concat_sep)
            df <- unite(df, new_col, !!!syms(columns), sep = concat_sep)
            columns <- new_col
        }
    }
    if (isTRUE(force_factor)) {
        p <- parent.frame()
        df_name <- deparse(substitute(df))
        for (col in columns) {
            p[[df_name]][[col]] <- factor(df[[col]], levels = unique(df[[col]]))
        }
    }
    return(columns)
}


#' Get a ggplot layer for background
#' @keywords internal
#' @param data A data frame
#' @param x A character string specifying the column name of the data frame to plot for the x-axis
#' @param palette A character string specifying the palette to use
#' @param palcolor A character string specifying the color to use in the palette
#' @param alpha A numeric value specifying the transparency of the plot
#' @param keep_empty A logical value indicating whether to keep empty groups
#' @return A ggplot layer for background
#' @importFrom ggplot2 geom_rect
bg_layer <- function(data, x, palette, palcolor, alpha, keep_empty) {
    f <- data[[x]]
    if (isFALSE(keep_empty)) {
        f <- droplevels(f)
    }
    bg_color <- palette_scp(levels(f), palette = palette, palcolor = palcolor, keep_names = TRUE)

    bg_data <- data.frame(x = factor(levels(f), levels = levels(f)))
    bg_data$x <- as.numeric(bg_data[[x]])
    bg_data$xmin <- ifelse(bg_data$x == min(bg_data$x), -Inf, bg_data$x - 0.5)
    bg_data$xmax <- ifelse(bg_data$x == max(bg_data$x), Inf, bg_data$x + 0.5)
    bg_data$ymin <- -Inf
    bg_data$ymax <- Inf
    bg_data$fill <- bg_color[levels(f)]

    geom_rect(
        data = bg_data,
        xmin = bg_data$xmin, xmax = bg_data$xmax, ymin = bg_data$ymin, ymax = bg_data$ymax,
        fill = bg_data$fill, alpha = alpha, inherit.aes = FALSE
    )
}
