#' A ggplot2 theme and palettes for scplotter
#' Borrowed from the `theme_scp` function in the `SCP` pipeline
#'
#' @seealso \url{https://github.com/zhanghao-njmu/SCP}
#'
#' @param aspect.ratio The aspect ratio of the plot
#' @param base_size The base size of the text
#' @param font_family The font family of the text
#' @param ... Other arguments for `theme()`
#'
#' @return A ggplot2 theme
#' @importFrom ggplot2 element_text element_rect margin theme element_blank unit
#' @export
theme_scp <- function(aspect.ratio = NULL, base_size = 12, font_family = NULL, ...) {
    text_size_scale <- base_size / 12
    args1 <- list(
        aspect.ratio = aspect.ratio,
        text = element_text(size = 12 * text_size_scale, family = font_family, color = "black"),
        plot.title = element_text(size = 14 * text_size_scale, family = font_family, colour = "black", vjust = 1),
        plot.subtitle = element_text(size = 13 * text_size_scale, family = font_family, hjust = 0, margin = margin(b = 3)),
        plot.background = element_rect(fill = "white", color = "white"),
        plot.margin = margin(10, 10, 10, 10),
        axis.line = element_blank(),
        axis.title = element_text(size = 13 * text_size_scale, family = font_family, colour = "black"),
        axis.text = element_text(size = 12 * text_size_scale, family = font_family, colour = "black"),
        strip.text = element_text(size = 12.5 * text_size_scale, family = font_family, colour = "black", hjust = 0.5, margin = margin(3, 3, 3, 3)),
        strip.background = element_rect(fill = "transparent", linetype = 0),
        strip.switch.pad.grid = unit(-1, "pt"),
        strip.switch.pad.wrap = unit(-1, "pt"),
        strip.placement = "outside",
        legend.title = element_text(size = 12 * text_size_scale, family = font_family, colour = "black", hjust = 0),
        legend.text = element_text(size = 11 * text_size_scale, family = font_family, colour = "black"),
        legend.key = element_rect(fill = "transparent", color = "transparent"),
        legend.key.size = unit(10, "pt"),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = "white"),
        panel.border = element_rect(fill = "transparent", colour = "black", linewidth = 1)
    )
    args2 <- as.list(match.call())[-1]
    call.envir <- parent.frame(1)
    args2 <- lapply(args2, function(arg) {
        if (is.symbol(arg)) {
            eval(arg, envir = call.envir)
        } else if (is.call(arg)) {
            eval(arg, envir = call.envir)
        } else {
            arg
        }
    })
    for (n in names(args2)) {
        args1[[n]] <- args2[[n]]
    }
    args <- args1[names(args1) %in% formalArgs(theme)]
    do.call(what = theme, args = args)
}

#' Blank theme
#'
#' This function creates a theme with all elements blank except for axis lines and labels.
#' It can optionally add coordinate axes in the plot.
#'
#' @param add_coord Whether to add coordinate arrows. Default is \code{TRUE}.
#' @param xlen_npc The length of the x-axis arrow in "npc".
#' @param ylen_npc The length of the y-axis arrow in "npc".
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param lab_size Label size.
#' @param ... Arguments passed to the \code{\link[ggplot2]{theme}}.
#'
#' @examples
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(x = wt, y = mpg, colour = factor(cyl))) +
#'     geom_point()
#' p + theme_blank()
#' p + theme_blank(xlab = "x-axis", ylab = "y-axis", lab_size = 16)
#' @importFrom ggplot2 theme element_blank margin annotation_custom coord_cartesian
#' @importFrom grid grobTree gList linesGrob textGrob arrow gpar
#' @export
theme_blank <- function(add_coord = TRUE, xlen_npc = 0.15, ylen_npc = 0.15, xlab = "", ylab = "", lab_size = 12, ...) {
    args1 <- list(
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.background = element_blank(),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.margin = margin(0, 0, 0, 0),
        legend.key.size = unit(10, "pt"),
        plot.margin = margin(lab_size + 2, lab_size + 2, lab_size + 2, lab_size + 2, unit = "points")
    )
    args2 <- as.list(match.call())[-1]
    call.envir <- parent.frame(1)
    args2 <- lapply(args2, function(arg) {
        if (is.symbol(arg)) {
            eval(arg, envir = call.envir)
        } else if (is.call(arg)) {
            eval(arg, envir = call.envir)
        } else {
            arg
        }
    })
    for (n in names(args2)) {
        args1[[n]] <- args2[[n]]
    }
    args <- args1[names(args1) %in% formalArgs(theme)]
    out <- do.call(
        what = theme,
        args = args
    )
    if (isTRUE(add_coord)) {
        g <- grobTree(gList(
            linesGrob(x = unit(c(0, xlen_npc), "npc"), y = unit(c(0, 0), "npc"), arrow = arrow(length = unit(0.02, "npc")), gp = gpar(lwd = 2)),
            textGrob(label = xlab, x = unit(0, "npc"), y = unit(0, "npc"), vjust = 4 / 3, hjust = 0, gp = gpar(fontsize = lab_size)),
            linesGrob(x = unit(c(0, 0), "npc"), y = unit(c(0, ylen_npc), "npc"), arrow = arrow(length = unit(0.02, "npc")), gp = gpar(lwd = 2)),
            textGrob(label = ylab, x = unit(0, "npc"), y = unit(0, "npc"), vjust = -2 / 3, hjust = 0, rot = 90, gp = gpar(fontsize = lab_size))
        ))
        return(list(
            list(annotation_custom(g)),
            list(theme_scp() + out),
            list(coord_cartesian(clip = "off"))
        ))
    } else {
        return(list(
            list(theme_scp() + out)
        ))
    }
}

#' Color palettes collected in scplotter.
#'
#' @param x A vector of character/factor or numeric values. If missing, numeric values 1:n will be used as x.
#' @param n The number of colors to return for numeric values.
#' @param palette Palette name. All available palette names can be queried with \code{show_palettes()}.
#' @param palcolor Custom colors used to create a color palette.
#' @param type Type of \code{x}. Can be one of "auto", "discrete" or "continuous". The default is "auto", which automatically detects if \code{x} is a numeric value.
#' @param keep_names Whether to keep the names of the color vector.
#' @param matched If \code{TRUE}, will return a color vector of the same length as \code{x}.
#' @param reverse Whether to invert the colors.
#' @param NA_keep Whether to keep the color assignment to NA in \code{x}.
#' @param NA_color Color assigned to NA if NA_keep is \code{TRUE}.
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom stats setNames
#' @export
#'
palette_scp <- function(
    x, n = 100, palette = "Paired", palcolor = NULL, type = "auto", keep_names = FALSE,
    matched = FALSE, reverse = FALSE, NA_keep = FALSE, NA_color = "grey80") {
    palette_list <- scplotter::palette_list
    if (missing(x)) {
        x <- 1:n
        type <- "continuous"
    }
    if (!palette %in% names(palette_list)) {
        stop("The palette name (", palette, ") is invalid! You can check the available palette names with 'show_palettes()'. Or pass palette colors via the 'palcolor' parameter.")
    }
    if (is.list(palcolor)) {
        palcolor <- unlist(palcolor)
    }
    if (all(palcolor == "")) {
        palcolor <- palette_list[[palette]]
    }
    if (is.null(palcolor) || length(palcolor) == 0) {
        palcolor <- palette_list[[palette]]
    }
    if (!is.null(names(palcolor))) {
        if (all(x %in% names(palcolor))) {
            palcolor <- palcolor[intersect(names(palcolor), x)]
        }
    }
    pal_n <- length(palcolor)

    if (!type %in% c("auto", "discrete", "continuous")) {
        stop("'type' must be one of 'auto','discrete' and 'continuous'.")
    }
    if (type == "auto") {
        if (is.numeric(x)) {
            type <- "continuous"
        } else {
            type <- "discrete"
        }
    }

    if (type == "discrete") {
        if (!is.factor(x)) {
            x <- factor(x, levels = unique(x))
        }
        n_x <- nlevels(x)
        if (isTRUE(attr(palcolor, "type") == "continuous")) {
            color <- colorRampPalette(palcolor)(n_x)
        } else {
            color <- ifelse(rep(n_x, n_x) <= pal_n,
                palcolor[1:n_x],
                colorRampPalette(palcolor)(n_x)
            )
        }
        names(color) <- levels(x)
        if (any(is.na(x))) {
            color <- c(color, setNames(NA_color, "NA"))
        }
        if (isTRUE(matched)) {
            color <- color[x]
            color[is.na(color)] <- NA_color
        }
    } else if (type == "continuous") {
        if (!is.numeric(x) && all(!is.na(x))) {
            stop("'x' must be type of numeric when use continuous color palettes.")
        }
        if (all(is.na(x))) {
            values <- as.factor(rep(0, n))
        } else if (length(unique(na.omit(as.numeric(x)))) == 1) {
            values <- as.factor(rep(unique(na.omit(as.numeric(x))), n))
        } else {
            if (isTRUE(matched)) {
                values <- cut(x, breaks = seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = n + 1), include.lowest = TRUE)
            } else {
                values <- cut(1:100, breaks = seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = n + 1), include.lowest = TRUE)
            }
        }

        n_x <- nlevels(values)
        color <- ifelse(rep(n_x, n_x) <= pal_n,
            palcolor[1:n_x],
            colorRampPalette(palcolor)(n_x)
        )
        names(color) <- levels(values)
        if (any(is.na(x))) {
            color <- c(color, setNames(NA_color, "NA"))
        }
        if (isTRUE(matched)) {
            if (all(is.na(x))) {
                color <- NA_color
            } else if (length(unique(na.omit(x))) == 1) {
                color <- color[as.character(unique(na.omit(x)))]
                color[is.na(color)] <- NA_color
            } else {
                color <- color[as.character(values)]
                color[is.na(color)] <- NA_color
            }
        }
    }

    if (isTRUE(reverse)) {
        color <- rev(color)
    }
    if (!isTRUE(NA_keep)) {
        color <- color[names(color) != "NA"]
    }
    if (!isTRUE(keep_names)) {
        names(color) <- NULL
    }
    return(color)
}

#' Show the color palettes
#'
#' This function displays color palettes using ggplot2.
#'
#' @param palettes A list of color palettes. If `NULL`, uses default palettes.
#' @param type A character vector specifying the type of palettes to include. Default is "discrete".
#' @param index A numeric vector specifying the indices of the palettes to include. Default is `NULL`.
#' @param palette_names A character vector specifying the names of the SCP palettes to include. Default is `NULL`.
#' @param return_names A logical value indicating whether to return the names of the selected palettes. Default is `TRUE`.
#' @param return_palettes A logical value indicating whether to return the colors of selected palettes. Default is `FALSE`.
#'
#' @seealso \code{\link{palette_list}}
#' @seealso \href{../articles/all-available-palettes.html}{All available palettes}
#'
#' @examples
#' show_palettes(palettes = list(c("red", "blue", "green"), c("yellow", "purple", "orange")))
#' all_palettes <- show_palettes(return_palettes = TRUE)
#' names(all_palettes)
#' all_palettes[["simspec"]]
#' show_palettes(index = 1:10)
#' show_palettes(type = "discrete", index = 1:10)
#' show_palettes(type = "continuous", index = 1:10)
#' show_palettes(palette_names = c("Paired", "nejm", "simspec", "Spectral", "jet"), return_palettes = TRUE)
#'
#' @importFrom gglogger ggplot
#' @importFrom ggplot2 geom_col scale_fill_manual scale_x_continuous element_blank aes element_text margin element_rect unit theme
#' @export
show_palettes <- function(palettes = NULL, type = c("discrete", "continuous"), index = NULL, palette_names = NULL, return_names = TRUE, return_palettes = FALSE) {
    if (!is.null(palettes)) {
        palette_list <- palettes
    } else {
        palette_list <- scplotter::palette_list[unlist(lapply(palette_list, function(x) isTRUE(attr(x, "type") %in% type)))]
    }
    index <- index[index %in% seq_along(palette_list)]
    if (!is.null(index)) {
        palette_list <- palette_list[index]
    }
    if (is.null(names(palette_list))) {
        names(palette_list) <- seq_along(palette_list)
    }
    if (is.null(palette_names)) {
        palette_names <- palette_names %||% names(palette_list)
    }
    if (any(!palette_names %in% names(palette_list))) {
        stop(paste("Can not find the palettes: ", paste0(palette_names[!palette_names %in% names(palette_list)], collapse = ",")))
    }
    palette_list <- palette_list[palette_names]

    df <- data.frame(palette = rep(names(palette_list), sapply(palette_list, length)), color = unlist(palette_list))
    df$palette <- factor(df$palette, levels = rev(unique(df$palette)))
    df$color_order <- factor(seq_len(nrow(df)), levels = seq_len(nrow(df)))
    df$proportion <- as.numeric(1 / table(df$palette)[df$palette])
    p <- ggplot(data = df, aes(y = palette, x = proportion, fill = color_order)) +
        geom_col(show.legend = FALSE) +
        scale_fill_manual(values = df[["color"]]) +
        scale_x_continuous(expand = c(0, 0), trans = "reverse") +
        theme_scp(
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_blank(),
            panel.border = element_blank()
        )
    print(p)

    if (isTRUE(return_palettes)) {
        return(palette_list)
    }
    if (isTRUE(return_names)) {
        return(palette_names)
    }
}
