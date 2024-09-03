
#' Cell statistics plot
#'
#' @description Plot the statistics of the cells.
#'
#' @rdname CellStatPlot
#'
#' @param x The data frame or Seurat object
#' @param ident The column with the cell identities. i.e. clusters. Default: seurat_clusters
#' @param plot_type The type of plot to use. Default is "bar".
#'   Possible values are "bar", "circos", "pie", "ring"/"donut", "trend", "area", and "sankey"/"alluvial".
#' @param frac The way of calculating the fraction. Default is "none".
#'   Possible values are "group", "ident", "cluster", "all", "none".
#'   - group: calculate the fraction in each group.
#'       The total fraction of the cells of idents in each group will be 1.
#'       When `group-by` is not specified, it will be the same as `all`.
#'   - ident: calculate the fraction in each ident.
#'       The total fraction of the cells of groups in each ident will be 1.
#'       Only works when `group-by` is specified.
#'   - cluster: alias of `ident`.
#'   - all: calculate the fraction against all cells.
#'   - none: do not calculate the fraction, use the number of cells instead.
#' @param flip Whether to flip the cluster and group, that is,
#'   using group as the x-axis and cluster to fill the plot.
#'   For circos plot, when transposed, the arrows will be drawn from the idents (by `ident`) to the
#'   the groups (by `group_by`).
#'   Only works when `group_by` is specified.
#' @param bar_position The position of the bars. Default is "auto".
#'   - stack: Use `position_stack()`.
#'   - fill: Use `position_fill()`.
#'   - dodge: Use `position_dodge()`.
#'   - auto: Use `stack` when there are more than 5 groups, otherwise use `dodge`.
#' @param pie_label Whether to add labels to the pie plot. Default is TRUE.
#' @param circos_labels_rot Whether to rotate the labels in the circos plot.
#'   In case the labels are too long. Default is FALSE.
#' @inheritParams validate_common_arguments
#' @param ... Additional arguments
#'
#' @return A ggplot object or a list with the plot and the height and width of the plot if guess_size is TRUE
#' @export
#' @examples
#' library(patchwork)
#' data(ifnb_sub)
#'
#' # Bar plot
#' p1 <- CellStatPlot(ifnb_sub)
#' p2 <- CellStatPlot(ifnb_sub, group_by = "stim")
#' p3 <- CellStatPlot(ifnb_sub, group_by = "stim", palette = "Set2")
#' p4 <- CellStatPlot(ifnb_sub, group_by = "stim", bar_position = "dodge")
#' (p1 | p2) / (p3 | p4)
#'
#' # Fraction of cells
#' p1 <- CellStatPlot(ifnb_sub, group_by = "stim", frac = "group",
#'                    ident = "seurat_annotations", x_text_angle = 60)
#' p2 <- CellStatPlot(ifnb_sub, group_by = "stim", frac = "group")
#' p3 <- CellStatPlot(ifnb_sub, group_by = "stim", frac = "ident",
#'                    bar_position = "stack", alpha = .6)
#' p4 <- CellStatPlot(ifnb_sub, group_by = "stim", frac = "ident",
#'                    flip = TRUE, bar_position = "stack")
#' (p1 | p2) / (p3 | p4)
#'
#' # Splitting/Facetting the plot
#' CellStatPlot(ifnb_sub, split_by = "stim")
#' CellStatPlot(ifnb_sub, split_by = "stim", facet = TRUE)
#' CellStatPlot(ifnb_sub, split_by = "stim", facet = TRUE, nrow = 2)
#'
#' # Circos plot
#' CellStatPlot(ifnb_sub, group_by = "stim", plot_type = "circos")
#' CellStatPlot(ifnb_sub, group_by = "stim", ident = "seurat_annotations",
#'              plot_type = "circos", circos_labels_rot = TRUE)
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
#'              group_by = c("stim", "seurat_annotations"))
#'
#' # Sankey plot
#' CellStatPlot(ifnb_sub, plot_type = "sankey", group_by = "stim", alpha = .6)
#' CellStatPlot(ifnb_sub, plot_type = "sankey", alpha = .6,
#'              group_by = c("stim", "seurat_annotations", "orig.ident"))
#'
#' # Area plot
#' CellStatPlot(ifnb_sub, plot_type = "area", frac = "group",
#'              group_by = "seurat_annotations", split_by = "stim")
#'
CellStatPlot <- function(x, ...) { UseMethod("CellStatPlot") }

#' Default method for `CellStatPlot`
#' @rdname CellStatPlot
#' @export
CellStatPlot.default <- function(x, ...) {
    stop("CellStatPlot() is not implemented for objects of class '", class(x), "'")
}

#' CellStatPlot method for data frames
#'
#' @rdname CellStatPlot
#' @importFrom rlang sym syms
#' @importFrom dplyr %>% group_by summarise mutate
#' @importFrom tidyr drop_na
#' @importFrom patchwork wrap_plots
#' @export
CellStatPlot.data.frame <- function(
    x, ident = "seurat_clusters",
    plot_type = c("bar", "circos", "pie", "ring", "donut", "trend", "area", "sankey", "alluvial"),
    frac = c("none", "group", "ident", "cluster", "all"), flip = FALSE,
    bar_position = c("stack", "fill", "dodge", "auto"), pie_label = TRUE,
    group_by = NULL, group_by_sep = "_", split_by = NULL, split_by_sep = "_",
    palette = "Paired", palcolor = NULL, alpha = 1, keep_empty = FALSE,
    x_text_angle = 90, aspect.ratio = 1, legend.position = "right", legend.direction = "vertical",
    theme = "theme_scp", theme_args = list(), facet = FALSE, facet_scales = "fixed",
    res = 100, guess_size = FALSE, seed = 8525, nrow = NULL, ncol = NULL, byrow = TRUE,
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
    circos_labels_rot = FALSE, combine = TRUE,
    ...) {
    validate_common_arguments(
        split_by = split_by, seed = seed, res = res,
        guess_size = guess_size, facet = facet
    )

    if (!ident %in% colnames(x)) {
        stop("Column '", ident, "' not found in the data frame.")
    }
    if (!is.null(group_by)) {
        non_existing <- setdiff(group_by, colnames(x))
        if (length(non_existing) > 0) {
            stop("'group_by' columns '", paste(non_existing, collapse = "', '"), "' not found in the data frame.")
        }
    }
    if (!is.null(split_by)) {
        non_existing <- setdiff(split_by, colnames(x))
        if (length(non_existing) > 0) {
            stop("'split_by' columns '", paste(non_existing, collapse = "', '"), "' not found in the data frame.")
        }
        if (length(split_by) > 1) {
            x <- concat_cols(x, split_by, sep = split_by_sep)
            split_by <- attr(x, "new_col")
        }
    }
    if (isTRUE(flip) && is.null(group_by)) {
        stop("Cannot flip the axes without specifying 'group_by'.")
    }

    plot_type <- match.arg(plot_type)
    if (plot_type == "donut") plot_type <- "ring"
    if (plot_type == "alluvial") plot_type <- "sankey"

    bar_position <- match.arg(bar_position)

    frac <- match.arg(frac)
    if (frac == "cluster") frac <- "ident"
    if (is.null(group_by) && frac == "ident") {
        stop("Cannot calculate the fraction by 'ident' without specifying 'group_by'.")
    }

    cells <- x %>% drop_na(!!sym(ident))
    if (nrow(cells) == 0) {
        stop("No data found with column '", ident, "'")
    }

    if (!is.factor(cells[[ident]])) {
        cells[[ident]] <- factor(cells[[ident]], levels = unique(cells[[ident]]))
    }
    if (!is.null(group_by)) {
        for (g in group_by) {
            if (!is.factor(cells[[g]])) {
                cells[[g]] <- factor(cells[[g]], levels = unique(cells[[g]]))
            }
        }
    }
    handle_split_df <- function(df) {
        if (!is.null(group_by)) {
            df <- df %>%
                group_by(!!sym(ident), !!!syms(group_by), .drop = FALSE) %>%
                summarise(.n = n(), .groups = "drop")
            if (frac == "group") {
                df <- df %>% group_by(!!!syms(group_by), .drop = FALSE) %>% mutate(.frac = .n / sum(.n))
            } else if (frac == "ident") {
                df <- df %>% group_by(!!sym(ident), .drop = FALSE) %>% mutate(.frac = .n / sum(.n))
            } else if (frac == "all") {
                df <- df %>% mutate(.frac = .n / sum(.n))
            } else {
                df <- df %>% mutate(.frac = 1)
            }
        } else {
            df <- df %>%
                group_by(!!sym(ident), .drop = FALSE) %>%
                summarise(.n = n(), .groups = "drop")
            if (frac != "none") {
                df <- df %>% mutate(.frac = .n / sum(.n))
            } else {
                df <- df %>% mutate(.frac = 1)
            }
        }
        df
    }
    if (!is.null(split_by)) {
        cells <- split(cells, cells[[split_by]])
        split_values <- names(cells)
        cells <- lapply(split_values, function(name) {
            df <- handle_split_df(cells[[name]])
            if (isTRUE(facet)) df[[split_by]] <- name
            df
        })
        names(cells) <- split_values
        if (isTRUE(facet)) {
            cells <- list(do.call(rbind, cells))
        }
    } else {
        cells <- list(handle_split_df(cells))
    }

    plots <- lapply(cells, function(df) {
        class(df) <- c(plot_type, class(df))
        suppressWarnings(CellStatPlotAtomic(
            df, ident = ident, frac = frac, flip = flip, bar_position = bar_position, pie_label = pie_label, # 6
            group_by = group_by, group_by_sep = group_by_sep, split_by = split_by, split_by_sep = split_by_sep, # 4
            palette = palette, palcolor = palcolor, alpha = alpha, keep_empty = keep_empty, # 4
            x_text_angle = x_text_angle, aspect.ratio = aspect.ratio, legend.position = legend.position, # 3
            legend.direction = legend.direction, theme = theme, theme_args = theme_args, # 3
            facet = facet, facet_scales = facet_scales, res = res, guess_size = guess_size, seed = seed, # 5
            title = title, subtitle = title, xlab = xlab, ylab = ylab, # 4
            nrow = nrow, ncol = ncol, byrow = byrow, circos_labels_rot = circos_labels_rot, ...)) # 4
    })

    combine_plots(plots, combine = combine, nrow = nrow, ncol = ncol, byrow = byrow)
}

#' CellStatPlot method for Seurat objects
#'
#' @rdname CellStatPlot
#' @export
CellStatPlot.Seurat <- function(x, ident = "seurat_clusters", ...) {
    if (!ident %in% colnames(x@meta.data)) {
        x@meta.data[[ident]] <- Seurat::Idents(x)
    }
    CellStatPlot(x@meta.data, ident = ident, ...)
}

#' CellStatPlotAtomic
#'
#' @inheritParams CellStatPlot
#'
#' @return A ggplot object or a list with the plot and the height and width of the plot if guess_size is TRUE
#' @keywords internal
CellStatPlotAtomic <- function(x, ...) UseMethod("CellStatPlotAtomic")

#' Default method for `CellStatPlotAtomic`
#' @keywords internal
CellStatPlotAtomic.default <- function(x, ...) {
    stop("CellStatPlotAtomic() is not implemented for objects of class '", class(x), "'")
}


#' CellStatPlotAtomic method for bar plot
#'
#' @inheritParams CellStatPlot
#' @keywords internal
#' @importFrom rlang syms
#' @importFrom dplyr %>% mutate
#' @importFrom ggplot2 geom_col scale_fill_manual labs facet_grid geom_rect
CellStatPlotAtomic.bar <- function(
    x, ident = "seurat_clusters",
    frac = "none", flip = FALSE, bar_position = "auto", pie_label = TRUE,
    group_by = NULL, group_by_sep = "_", split_by = NULL, split_by_sep = "_",
    palette = "Paired", palcolor = NULL, alpha = 1, keep_empty = FALSE,
    x_text_angle = 90, aspect.ratio = 1, legend.position = "right", legend.direction = "vertical",
    theme = "theme_scp", theme_args = list(), facet = FALSE, facet_scales = "fixed",
    res = 100, guess_size = FALSE, seed = 8525, nrow = NULL, ncol = NULL, byrow = TRUE,
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
    circos_labels_rot = FALSE,
    ...) {
    if (length(group_by) > 1) {
        x <- concat_cols(x, group_by, sep = group_by_sep)
        group_by <- attr(x, "new_col")
    }
    ngroups <- ifelse(is.null(group_by), 1, length(levels(x[[group_by]])))
    nidents <- length(unique(x[[ident]]))
    bar_position <- ifelse(bar_position == "auto", ifelse(ngroups > 5, "stack", "dodge"), bar_position)

    # theme_args$add_coord <- theme_args$add_coord %||% FALSE
    if (is.null(group_by)) {
        x_by <- ident
        p <- ggplot(x, aes(
            x = !!sym(x_by),
            y = if (frac == "none") .n else .frac,
            fill = !!sym(x_by)
        ))
        ncolors <- nidents
        guide <- "legend"
        nx <- nidents
    } else {
        x_by <- ifelse(isTRUE(flip), group_by, ident)
        p <- ggplot(x, aes(
            x = !!sym(x_by),
            y = if (frac == "none") .n else .frac,
            fill = !!sym(ifelse(is.null(group_by) || isTRUE(flip), ident, group_by))
        ))
        ncolors <- ifelse(isTRUE(flip), nidents, ngroups)
        guide <- "legend"
        nx <- ifelse(isTRUE(flip), ngroups, nidents)
    }

    if (bar_position == "dodge" && ngroups > 1) {
        bg_color <- palette_scp(levels(x[[x_by]]), nx, palcolor = rep(c("transparent", "grey85"), nx), keep_names = TRUE)
        bg_data <- na.omit(unique(x[, x_by, drop = FALSE]))
        bg_data$x <- as.numeric(bg_data[[x_by]])
        bg_data$xmin <- ifelse(bg_data$x == min(bg_data$x), -Inf, bg_data$x - 0.5)
        bg_data$xmax <- ifelse(bg_data$x == max(bg_data$x), Inf, bg_data$x + 0.5)
        bg_data$ymin <- -Inf
        bg_data$ymax <- Inf
        bg_data$fill <- bg_color[as.character(bg_data$x)]
        p <- p + geom_rect(
            data = bg_data,
            xmin = bg_data$xmin,
            xmax = bg_data$xmax,
            ymin = bg_data$ymin,
            ymax = bg_data$ymax,
            fill = bg_data$fill,
            inherit.aes = FALSE,
            alpha = 0.2
        ) +
        geom_col(width = 0.8, color = "black", position = bar_position, alpha = alpha) +
        scale_y_continuous(expand = c(0, 0, 0.02, 0))
    } else {
        p <- p +
            geom_col(width = 0.9, color = "black", position = bar_position, alpha = alpha) +
            scale_x_discrete(drop = !keep_empty, expand = c(0, 0)) +
            scale_y_continuous(expand = c(0, 0))
    }

    p <- p +
        scale_fill_manual(
            values = palette_scp(n = ncolors, palette = palette, palcolor = palcolor),
            guide = guide) +
        labs(title = title, subtitle = subtitle, x = xlab %||% "", y = ylab %||% ifelse(frac == "none", "Number of cells", "Fraction of cells")) +
        do.call(theme, theme_args) +
        ggplot2::theme(
            aspect.ratio = aspect.ratio,
            legend.position = legend.position,
            legend.direction = legend.direction,
            panel.grid.major = element_line(colour = "grey80", linetype = 2),
            axis.text.x = element_text(angle = x_text_angle, hjust = 1)
        )

    if (isTRUE(guess_size)) {
        p <- list(
            plot = p,
            height = 4 * res,
            width = 0.4 * nx * res + 120
        )
    }
    facet_plot(p, facet, facet_scales, split_by, nrow, ncol, byrow)
}

#' CellStatPlotAtomic method for circos plot
#'
#' @inheritParams CellStatPlot
#' @keywords internal
#' @import gridGraphics
#' @importFrom circlize chordDiagram circos.clear circos.par circos.text circos.track circos.axis get.cell.meta.data get.all.sector.index mm_h CELL_META
#' @importFrom patchwork wrap_elements
#' @importFrom rlang sym syms
#' @importFrom dplyr %>% select mutate
CellStatPlotAtomic.circos <- function(
    x, ident = "seurat_clusters",
    frac = "none", flip = FALSE, bar_position = "auto", pie_label = TRUE,
    group_by = NULL, group_by_sep = "_", split_by = NULL, split_by_sep = "_",
    palette = "Paired", palcolor = NULL, alpha = 1, keep_empty = FALSE,
    x_text_angle = 90, aspect.ratio = 1, legend.position = "right", legend.direction = "vertical",
    theme = "theme_scp", theme_args = list(), facet = FALSE, facet_scales = "fixed",
    res = 100, guess_size = FALSE, seed = 8525, nrow = NULL, ncol = NULL, byrow = TRUE,
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
    circos_labels_rot = FALSE,
    ...) {
    if (is.null(group_by)) {
        stop("Cannot create a circos plot without specifying 'group_by'.")
    }
    if (length(group_by) > 1) {
        x <- concat_cols(x, group_by, sep = group_by_sep)
        group_by <- attr(x, "new_col")
    }
    if (isTRUE(facet)) {
        stop("Cannot facet a circos plot. Set 'facet' to FALSE to use patchwork to combine the plots.")
    }
    if (isTRUE(flip)) {
        circos_df <- x %>% select(from = !!sym(ident), to = !!sym(group_by), value = .n)
    } else {
        circos_df <- x %>% select(from = !!sym(group_by), to = !!sym(ident), value = .n)
    }
    circos_df$from <- as.character(circos_df$from)
    circos_df$to <- as.character(circos_df$to)

    groups <- unique(circos_df$from)
    idents <- unique(circos_df$to)
    grid_cols <- palette_scp(n = length(idents))
    names(grid_cols) <- idents
    gcols <- rep("grey50", length(groups))
    names(gcols) <- groups
    grid_cols <- c(grid_cols, gcols)
    link_cols <- grid_cols[circos_df$to]

    circos.clear()
    circos.par(track.margin = c(0.01, 0.02))

    if (!isTRUE(circos_labels_rot)) {
        p <- ~ {
            chordDiagram(
                circos_df,
                grid.col = grid_cols,
                col = link_cols,
                direction = 1,
                annotationTrack = "grid",
                direction.type = c("diffHeight", "arrows"),
                link.arr.type = "big.arrow",
                link.arr.length = .04,
                preAllocateTracks = list(
                    list(track.height = mm_h(1)),
                    list(track.height = mm_h(.1))
                )
            )
            circos.track(track.index = 1, panel.fun = function(x, y) {
                circos.text(
                    CELL_META$xcenter, CELL_META$ylim[1] + 5.5,
                    CELL_META$sector.index,
                    niceFacing = TRUE,
                    adj = c(0.5, 0.5))
            }, bg.border = NA)
            circos.track(track.index = 2, panel.fun = function(x, y) {
                for (si in get.all.sector.index()) {
                    start.degree <- get.cell.meta.data("cell.start.degree", sector.index = si)
                    end.degree <- get.cell.meta.data("cell.end.degree", sector.index = si)
                    if (abs(end.degree - start.degree) > 2) {
                        # otherwise: patchwork wrap_elements 'x' and 'units' must have length > 0
                        circos.axis(h = "top", labels.cex = 0.7, labels.niceFacing = TRUE, sector.index = si)
                    }
                }
            }, bg.border = NA) # here set bg.border to NA is important
        }
    } else {
        allnames <- unique(c(groups, idents))
        p <- ~ {
            chordDiagram(
                circos_df,
                grid.col = grid_cols,
                col = link_cols,
                direction = 1,
                annotationTrack = "grid",
                direction.type = c("diffHeight", "arrows"),
                link.arr.type = "big.arrow",
                link.arr.length = .04,
                preAllocateTracks = list(
                    list(track.height = max(strwidth(allnames))),
                    list(track.height = mm_h(.1))
                )
            )
            circos.track(track.index = 1, panel.fun = function(x, y) {
                circos.text(
                    CELL_META$xcenter, CELL_META$ylim[1] +.15,
                    CELL_META$sector.index,
                    facing = "clockwise",
                    niceFacing = TRUE,
                    adj = c(0, 0.5))
            }, bg.border = NA)
            circos.track(track.index = 2, panel.fun = function(x, y) {
                for (si in get.all.sector.index()) {
                    start.degree <- get.cell.meta.data("cell.start.degree", sector.index = si)
                    end.degree <- get.cell.meta.data("cell.end.degree", sector.index = si)
                    if (abs(end.degree - start.degree) > 2) {
                        # otherwise: patchwork wrap_elements 'x' and 'units' must have length > 0
                        circos.axis(h = "top", labels.cex = 0.7, labels.niceFacing = TRUE, sector.index = si)
                    }
                }
            }, bg.border = NA) # here set bg.border to NA is important
        }
    }
    p <- wrap_elements(full = p)

    if (isTRUE(guess_size)) {
        base_size <- 7
        if (isTRUE(circos_labels_rot)) {
            maxchar <- max(c(nchar(ident), nchar(group_by)))
            if (maxchar < 16) {
                base_size <- base_size + 2
            } else if (maxchar < 32) {
                base_size <- base_size + 4
            } else {
                base_size <- base_size + 6
            }
        }

        p <- list(
            plot = p,
            height = base_size * res,
            width = base_size * res
        )
    }
    facet_plot(p, facet, facet_scales, split_by, nrow, ncol, byrow)
}


#' CellStatPlotAtomic method for pie plot
#'
#' @inheritParams CellStatPlot
#' @keywords internal
#' @importFrom dplyr %>% mutate if_else lead
#' @importFrom ggplot2 geom_col scale_fill_manual labs facet_grid wrap_dims ggplot_build coord_polar
#' @importFrom ggrepel geom_label_repel
CellStatPlotAtomic.pie <- function(
    x, ident = "seurat_clusters",
    frac = "none", flip = FALSE, bar_position = "auto", pie_label = TRUE,
    group_by = NULL, group_by_sep = "_", split_by = NULL, split_by_sep = "_",
    palette = "Paired", palcolor = NULL, alpha = 1, keep_empty = FALSE,
    x_text_angle = 90, aspect.ratio = 1, legend.position = "right", legend.direction = "vertical",
    theme = "theme_scp", theme_args = list(), facet = FALSE, facet_scales = "fixed",
    res = 100, guess_size = FALSE, seed = 8525, nrow = NULL, ncol = NULL, byrow = TRUE,
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
    circos_labels_rot = FALSE,
    ...) {
    if (!is.null(group_by)) {
        stop("Cannot create a pie plot with 'group_by'. You may want to use 'ring' plot instead.")
    }
    if (isTRUE(facet)) {
        stop("Cannot facet a pie plot. Set 'facet' to FALSE to use patchwork to combine the plots.")
    }
    base_size <- theme_args$base_size %||% 12
    text_size_scale <- base_size / 12

    pos_df <- x %>% mutate(
        csum = rev(cumsum(rev(.n))),
        pos = .n / 2 + lead(csum, 1),
        pos = if_else(is.na(pos), .n / 2, pos)
    )
    p = ggplot(x, aes(x="", y = .n, fill=!!sym(ident))) +
        geom_col(width = 1, color = "grey50", alpha = alpha) +
        coord_polar("y", start=0) +
        scale_fill_manual(values = palette_scp(n = length(unique(x[[ident]])), palette = palette, palcolor = palcolor)) +
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
        ) +
        labs(title = title, subtitle = subtitle, x = xlab %||% "", y = ylab %||% "")

    if (pie_label) {
        p <- p + geom_label_repel(
            data = pos_df,
            if (frac != "none")
                aes(y = pos, label = sprintf("%.1f%%", .frac * 100))
            else
                aes(y = pos, label = .n),
            nudge_x = 0.1,
            color = "grey20",
            fill = "#fcfcfc",
            size = text_size_scale * 3
        )
    }

    if (isTRUE(guess_size)) {
        d <- wrap_dims(length(unique(ggplot_build(p)$data[[1]]$PANEL)))
        height <- 6 * res
        if (d[2] > 1) {
            width <- 7 * res
        } else {
            width <- 8 * res
        }

        p <- list(plot = p, height = height, width = width)
    }
    facet_plot(p, facet, facet_scales, split_by, nrow, ncol, byrow)
}


#' CellStatPlotAtomic method for ring plot
#'
#' @inheritParams CellStatPlot
#' @keywords internal
#' @importFrom dplyr %>% mutate if_else lead summarise group_by first n ungroup
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 geom_col scale_fill_manual labs facet_grid wrap_dims ggplot_build coord_polar geom_label scale_x_discrete
CellStatPlotAtomic.ring <- function(
    x, ident = "seurat_clusters",
    frac = "none", flip = FALSE, bar_position = "auto", pie_label = TRUE,
    group_by = NULL, group_by_sep = "_", split_by = NULL, split_by_sep = "_",
    palette = "Paired", palcolor = NULL, alpha = 1, keep_empty = FALSE,
    x_text_angle = 90, aspect.ratio = 1, legend.position = "right", legend.direction = "vertical",
    theme = "theme_scp", theme_args = list(), facet = FALSE, facet_scales = "fixed",
    res = 100, guess_size = FALSE, seed = 8525, nrow = NULL, ncol = NULL, byrow = TRUE,
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
    circos_labels_rot = FALSE,
    ...) {
    if (isTRUE(facet)) {
        stop("Cannot facet a pie plot. Set 'facet' to FALSE to use patchwork to combine the plots.")
    }
    if (length(group_by) > 1) {
        x <- concat_cols(x, group_by, sep = group_by_sep)
        group_by <- attr(x, "new_col")
    }
    base_size <- theme_args$base_size %||% 12
    text_size_scale <- base_size / 12

    if (is.null(group_by)) {
        x$group <- "group"
        x$group <- factor(x$group, levels = unique(x$group))
        group_by <- "group"
        label <- FALSE
    } else {
        label <- TRUE
    }
    if (frac != "group") {
        warning("'frac' is forced to 'group' for 'ring' plot.", immediate. = TRUE)
        # recalculate the .frac
        if (!".n" %in% colnames(x)) {
            x <- x %>% group_by(!!sym(ident), !!sym(group_by)) %>%
                summarise(.n = n(), .groups = "drop")
        }
        x <- x %>% group_by(!!sym(group_by)) %>%
            mutate(.frac = .n / sum(.n)) %>%
            ungroup()
    }

    groups <- levels(x[[group_by]])
    label_df <- data.frame(groups = groups)

    p = ggplot(x, aes(x = !!sym(group_by), y = .frac, fill = !!sym(ident))) +
        geom_col(width = 0.9, color = "grey50", alpha = alpha) +
        coord_polar("y", start=0) +
        scale_fill_manual(values = palette_scp(n = length(unique(x[[ident]])), palette = palette, palcolor = palcolor)) +
        scale_x_discrete(limits = c(" ", groups)) +
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
        ) +
        labs(title = title, subtitle = subtitle, x = xlab %||% "", y = ylab %||% "")

    if (label) {
        p <- p + geom_label(
            aes(label = groups, x = groups, y = 0),
            data = label_df,
            inherit.aes = FALSE,
            color = "grey20",
            size = text_size_scale * 3
        )
    }

    if (isTRUE(guess_size)) {
        d <- wrap_dims(length(unique(ggplot_build(p)$data[[1]]$PANEL)))
        height <- 6 * res
        if (d[2] > 1) {
            width <- 7 * res
        } else {
            width <- 8 * res
        }

        p <- list(plot = p, height = height, width = width)
    }
    facet_plot(p, facet, facet_scales, split_by, nrow, ncol, byrow)
}


#' CellStatPlotAtomic method for sankey plot
#'
#' @inheritParams CellStatPlot
#' @keywords internal
#' @importFrom rlang sym syms
#' @importFrom ggnewscale new_scale_fill
#' @importFrom dplyr %>% select group_by summarise all_of
#' @importFrom gglogger ggplot
#' @importFrom ggplot2 aes geom_text scale_fill_manual guides guide_legend geom_col labs
#' @importFrom ggalluvial geom_stratum geom_alluvium to_lodes_form StatStratum
CellStatPlotAtomic.sankey <- function(
    x, ident = "seurat_clusters",
    frac = "none", flip = FALSE, bar_position = "auto", pie_label = TRUE,
    group_by = NULL, group_by_sep = "_", split_by = NULL, split_by_sep = "_",
    palette = "Paired", palcolor = NULL, alpha = 1, keep_empty = FALSE,
    x_text_angle = 90, aspect.ratio = 1, legend.position = "right", legend.direction = "vertical",
    theme = "theme_scp", theme_args = list(), facet = FALSE, facet_scales = "fixed",
    res = 100, guess_size = FALSE, seed = 8525, nrow = NULL, ncol = NULL, byrow = TRUE,
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
    circos_labels_rot = FALSE,
    ...) {
    if (is.null(group_by)) {
        stop("Cannot create a sankey plot without specifying 'group_by'.")
    }
    if (isTRUE(facet)) {
        stop("Cannot facet a sankey/alluvial plot. Set 'facet' to FALSE to use patchwork to combine the plots.")
    }
    if (!ident %in% group_by) {
        group_by <- c(ident, group_by)
    }
    idents <- levels(x[[ident]])
    x <- to_lodes_form(x, key = "Group", value = "GroupValue", id = "Ident", axes = group_by)
    colors <- palette_scp(unique(x$GroupValue), palette = palette, palcolor = palcolor, keep_names = TRUE)
    just <- calc_just(x_text_angle)

    p <- ggplot(x)
    for (i in seq_along(group_by)) {
        if (group_by[i] == ident) next
        gvalues <- levels(x$GroupValue)
        gvalues <- gvalues[gvalues %in% x[x$Group == group_by[i], "GroupValue", drop = TRUE]]
        p <- p +
            geom_col(aes(x = Group, fill = GroupValue, y = 0), width = 0, color = "grey20") +
            scale_fill_manual(name = group_by[i], values = colors, breaks = gvalues,
                              guide = guide_legend(order = i + 1)) +
            new_scale_fill()
    }
    p <- p +
        geom_alluvium(aes(x = Group, stratum = GroupValue, alluvium = Ident, y = .n, fill = GroupValue), alpha = alpha) +
        geom_stratum(aes(x = Group, stratum = GroupValue, alluvium = Ident, y = .n, fill = GroupValue),
                     width = 0.25, color = "grey50") +
        scale_fill_manual(name = ident, values = colors, breaks = idents,
                          guide = guide_legend(order = 1)) +
        scale_x_discrete(drop = !keep_empty, expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        do.call(theme, theme_args) +
        labs(title = title, subtitle = subtitle, x = xlab %||% "", y = xlab %||% "Number of cells") +
        ggplot2::theme(
            aspect.ratio = aspect.ratio,
            legend.position = legend.position,
            legend.direction = legend.direction,
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_text(angle = x_text_angle, hjust = just$h, vjust = just$v)
        )

    if (isTRUE(guess_size)) {
        base_size <- 7
        if (isTRUE(circos_labels_rot)) {
            maxchar <- max(c(nchar(ident), nchar(group_by)))
            if (maxchar < 16) {
                base_size <- base_size + 2
            } else if (maxchar < 32) {
                base_size <- base_size + 4
            } else {
                base_size <- base_size + 6
            }
        }

        p <- list(plot = p, height = base_size * res, width = base_size * res)
    }
    facet_plot(p, facet, facet_scales, split_by, nrow, ncol, byrow)
}

#' CellStatPlotAtomic method for trend plot
#'
#' @inheritParams CellStatPlot
#' @keywords internal
#' @importFrom ggplot2 geom_area geom_col scale_fill_manual labs facet_grid position_stack
CellStatPlotAtomic.trend <- function(
    x, ident = "seurat_clusters",
    frac = "none", flip = FALSE, bar_position = "auto", pie_label = TRUE,
    group_by = NULL, group_by_sep = "_", split_by = NULL, split_by_sep = "_",
    palette = "Paired", palcolor = NULL, alpha = 1, keep_empty = FALSE,
    x_text_angle = 90, aspect.ratio = 1, legend.position = "right", legend.direction = "vertical",
    theme = "theme_scp", theme_args = list(), facet = FALSE, facet_scales = "fixed",
    res = 100, guess_size = FALSE, seed = 8525, nrow = NULL, ncol = NULL, byrow = TRUE,
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
    circos_labels_rot = FALSE,
    ...) {
    if (is.null(group_by)) {
        stop("Cannot create a trend plot without 'group_by'.")
    }
    if (length(group_by) > 1) {
        x <- concat_cols(x, group_by, sep = group_by_sep)
        group_by <- attr(x, "new_col")
    }
    if (frac == "ident") {
        stop("Cannot calculate the fraction by 'ident' for 'trend' plot.")
    }

    idents <- levels(x[[ident]])
    ngroups <- length(levels(x[[group_by]]))

    nr <- nrow(x)
    dat_area <- x[rep(seq_len(nr), each = 2), , drop = FALSE]
    dat_area[[group_by]] <- as.numeric(dat_area[[group_by]])
    dat_area[seq(1, nr * 2, 2), group_by] <- dat_area[seq(1, nr * 2, 2), group_by] - 0.2
    dat_area[seq(2, nr * 2, 2), group_by] <- dat_area[seq(2, nr * 2, 2), group_by] + 0.2

    position <- position_stack(vjust = 0.5)
    just <- calc_just(x_text_angle)
    p <- ggplot(x, aes(x = !!sym(group_by), y = if (frac == "none") .n else .frac, fill = !!sym(ident))) +
        geom_area(
            data = dat_area, mapping = aes(x = !!sym(group_by), fill = !!sym(ident)),
            alpha = alpha / 2, color = "grey50", position = position
        ) +
        geom_col(aes(fill = !!sym(ident)),
            width = 0.4,
            color = "black",
            alpha = alpha,
            position = position
        ) +
        scale_x_discrete(drop = !keep_empty, expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0), labels = if (frac == "none") scales::number else scales::percent) +
        scale_fill_manual(values = palette_scp(idents, palette = palette, palcolor = palcolor, keep_names = TRUE)) +
        labs(title = title, subtitle = subtitle, x = xlab %||% "", y = ylab %||% ifelse(frac == "none", "Number of cells", "Fraction of cells")) +
        do.call(theme, theme_args) +
        ggplot2::theme(
            aspect.ratio = aspect.ratio,
            legend.position = legend.position,
            legend.direction = legend.direction,
            panel.grid.major = element_line(colour = "grey80", linetype = 2),
            axis.text.x = element_text(angle = x_text_angle, hjust = just$h, vjust = just$v)
        )

    if (isTRUE(guess_size)) {
        p <- list(plot = p, height = 6 * res, width = ngroups * 0.5 * res)
    }
    facet_plot(p, facet, facet_scales, split_by, nrow, ncol, byrow)
}

#' CellStatPlotAtomic method for area plot
#'
#' @inheritParams CellStatPlot
#' @keywords internal
#' @importFrom ggplot2 geom_area geom_col scale_fill_manual labs facet_grid position_stack
CellStatPlotAtomic.area <- function(
    x, ident = "seurat_clusters",
    frac = "none", flip = FALSE, bar_position = "auto", pie_label = TRUE,
    group_by = NULL, group_by_sep = "_", split_by = NULL, split_by_sep = "_",
    palette = "Paired", palcolor = NULL, alpha = 1, keep_empty = FALSE,
    x_text_angle = 90, aspect.ratio = 1, legend.position = "right", legend.direction = "vertical",
    theme = "theme_scp", theme_args = list(), facet = FALSE, facet_scales = "fixed",
    res = 100, guess_size = FALSE, seed = 8525, nrow = NULL, ncol = NULL, byrow = TRUE,
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
    circos_labels_rot = FALSE,
    ...) {
    if (is.null(group_by)) {
        stop("Cannot create an area plot without 'group_by'.")
    }
    if (length(group_by) > 1) {
        x <- concat_cols(x, group_by, sep = group_by_sep)
        group_by <- attr(x, "new_col")
    }
    if (frac == "ident") {
        stop("Cannot calculate the fraction by 'ident' for 'area' plot.")
    }

    idents <- levels(x[[ident]])
    groups <- levels(x[[group_by]])

    just <- calc_just(x_text_angle)
    p <- ggplot(x, aes(x = as.numeric(!!sym(group_by)), y = if (frac == "none") .n else .frac, fill = !!sym(ident))) +
        geom_area(alpha = alpha, color = "grey50", position = position_stack(vjust = 0.5)) +
        scale_x_discrete(drop = !keep_empty, expand = c(0, 0), breaks = groups, limits = groups) +
        scale_y_continuous(expand = c(0, 0), labels = if (frac == "none") scales::number else scales::percent) +
        scale_fill_manual(values = palette_scp(idents, palette = palette, palcolor = palcolor, keep_names = TRUE)) +
        labs(title = title, subtitle = subtitle, x = xlab %||% "", y = ylab %||% ifelse(frac == "none", "Number of cells", "Fraction of cells")) +
        do.call(theme, theme_args) +
        ggplot2::theme(
            aspect.ratio = aspect.ratio,
            legend.position = legend.position,
            legend.direction = legend.direction,
            panel.grid.major = element_line(colour = "grey80", linetype = 2),
            axis.text.x = element_text(angle = x_text_angle, hjust = just$h, vjust = just$v)
        )

    if (isTRUE(guess_size)) {
        p <- list(plot = p, height = 6 * res, width = length(groups) * 0.5 * res)
    }
    facet_plot(p, facet, facet_scales, split_by, nrow, ncol, byrow)
}
