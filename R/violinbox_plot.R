# #' Violin, box or violin-box plot
# #'
# #' @inheritParams validate_common_arguments
# #' @param x Data frame with the data to plot.
# #' @param y A character string specifying the column name of the data frame to plot
# #'   for the y-axis.
# #' @param mode A character string specifying the mode of the plot.
# #'   It can be:
# #'   - violin: Only violin plot. Alias: v.
# #'   - box: Only box plot. Alias: b.
# #'   - violin+box: Violin and box plots together. Alias: v+b.
# #'   - jitter+box: Jitter and box plots together. Alias: j+b.
# #'   - jitter+violin: Jitter and violin plots together. Alias: j+v.
# #'   - jitter+violin+box: Jitter, violin and box plots together. Alias: j+v+b.
# #'   - violin-box: Violin plot on the left and box plot on the right. Alias: v-b.
# #'   - box-violin: Box plot on the left and violin plot on the right. Alias: b-v.
# #'   - jitter-box: Jitter plot on the left and box plot on the right. Alias: j-b.
# #'   - jitter-violin: Jitter plot on the left and violin plot on the right. Alias: j-v.
# #'   - box-jitter: Box plot on the left and jitter plot on the right. Alias: b-j.
# #'   - violin-jitter: Violin plot on the left and jitter plot on the right. Alias: v-j.
# #'   - jitter-violinbox: Jitter plot on the left, violin and box plot on the right. Alias: j-vb.
# #'   - jitterbox-violin: Jitter plot on the left, box and half violin plot on the right. Alias: jb-v.
# #'   - jitterviolin-box: Jitter and violin plot on the left, box plot on the right. Alias: jv-b.
# #'   - violinbox-jitter: Violin and box plot on the left, half jitter plot on the right. Alias: vb-j.
# #'   - violin-jitterbox: Violin plot on the left, jitter and box plot on the right. Alias: v-jb.
# #'   - box-jitterviolin: Box plot on the left, jitter and violin plot on the right. Alias: b-jv.
# #' @param dodge_by A character string specifying the variable to dodge the plot for side-by-side comparison.
# #'   It can also be a list of character strings. If so 'dodge_by_sep' will be used to concatenate the values.
# #' @param dodge_by_sep A character string specifying the separator to concatenate the values of 'dodge_by'.
# #' @param dodge_position_preserve A character string specifying the dodge position to preserve.
# #'   It can be "single" or "total".
# #' @param return_layer A logical value indicating whether to return the plot as a layer.
# #' @param ... Additional arguments
# #' @return A ggplot object if 'return_layer' is FALSE, otherwise a ggplot layer.
# #'   If 'guess_size' is TRUE, it returns a list with the ggplot object and the height and width of the plot.
# #'
# #' @importFrom gglogger ggplot
# #' @importFrom ggplot2 geom_boxplot geom_jitter geom_violin
# ViolinBoxPlot <- function(
#     x, y,
#     group_by = NULL, group_by_sep = "_", split_by = NULL, split_by_sep = "_", seed = 8525,
#     facet = FALSE, facet_scales = "fixed", guess_size = FALSE, res = 100,
#     theme = "theme_scp", theme_args = list(), palette = "Paired", palcolor = NULL,
#     keep_empty = FALSE, alpha = 1, x_text_angle = 0, aspect.ratio = 1, title = NULL, subtitle = NULL,
#     xlab = NULL, ylab = NULL, legend.position = "right", legend.direction = "vertical",
#     combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, dodge_position_preserve = "total",
#     mode = c(
#         "v+b", "jitter+box", "violin", "v", "box", "b", "violin+box", "j+b", "jitter+violin", "j+v",
#         "jitter+violin+box", "j+v+b", "violin-box", "v-b", "box-violin", "b-v", "jitter-box", "j-b",
#         "jitter-violin", "j-v", "box-jitter", "b-j", "violin-jitter", "v-j", "jitter-violinbox", "j-vb",
#         "jitterbox-violin", "jb-v", "jitterviolin-box", "jv-b", "violinbox-jitter", "vb-j",
#         "violin-jitterbox", "v-jb", "box-jitterviolin", "b-jv"
#     ),
#     dodge_by = NULL, dodge_by_sep = "_", return_layer = FALSE,
#     ...
# ) {
#     validate_common_arguments(
#         split_by = split_by, seed = seed, res = res,
#         guess_size = guess_size, facet = facet
#     )

#     mode <- match.arg(mode)
#     # convert the fullnames to shortcuts
#     mode <- sub("violin", "v", mode, fixed = TRUE)
#     mode <- sub("box", "b", mode, fixed = TRUE)
#     mode <- sub("jitter", "j", mode, fixed = TRUE)

#     if (guess_size && return_layer) {
#         stop("guess_size and return_layer cannot be TRUE at the same time.")
#     }

#     if (!y %in% colnames(x)) {
#         stop(paste0("'", y, "' is not in the x column names."))
#     }

#     if (is.null(group_by)) {
#         stop("'group_by' must be provided to create the ViolinBoxPlot.")
#     }
#     if (length(group_by) > 1) {
#         x <- concat_cols(x, group_by, group_by_sep)
#         group_by <- attr(x, "new_col")
#     }
#     if (!is.null(group_by) && !is.factor(x[[group_by]])) {
#         x[[group_by]] <- factor(x[[group_by]], levels = unique(x[[group_by]]))
#     }
#     if (length(dodge_by) > 1) {
#         x <- concat_cols(x, dodge_by, dodge_by_sep)
#         dodge_by <- attr(x, "new_col")
#     }
#     if (!is.null(dodge_by) && !is.factor(x[[dodge_by]])) {
#         x[[dodge_by]] <- factor(x[[dodge_by]], levels = unique(x[[dodge_by]]))
#     }

#     # split data
#     if (is.null(split_by) || isTRUE(facet)) {
#         xs <- list(x)
#     } else {
#         if (length(split_by) > 1) {
#             x <- concat_cols(x, split_by, split_by_sep)
#             split_by <- attr(x, "new_col")
#         }
#         if (!is.factor(x[[split_by]])) {
#             x[[split_by]] <- factor(x[[split_by]], levels = unique(x[[split_by]]))
#         }
#         xs <- split(x, x[[split_by]])
#         xs <- xs[levels(x[[split_by]])]
#     }
#     rm(x)

#     # plot
#     plots <- lapply(xs, function(x) {
#         ViolinBoxPlotAtomic(
#             x, y,
#             group_by = group_by, group_by_sep = group_by_sep, split_by = split_by, split_by_sep = split_by_sep,
#             seed = seed, facet = facet, facet_scales = facet_scales, guess_size = guess_size, res = res,
#             theme = theme, theme_args = theme_args, palette = palette, palcolor = palcolor,
#             keep_empty = keep_empty, alpha = alpha, x_text_angle = x_text_angle, aspect.ratio = aspect.ratio,
#             title = title, subtitle = subtitle, xlab = xlab, ylab = ylab, legend.position = legend.position,
#             legend.direction = legend.direction, combine = FALSE, nrow = nrow, ncol = ncol, byrow = byrow,
#             mode = mode, dodge_by = dodge_by, dodge_by_sep = dodge_by_sep, return_layer = return_layer,
#             dodge_position_preserve = dodge_position_preserve,
#             ...
#         )
#     })

#     combine_plots(plots, combine = combine, nrow = nrow, ncol = ncol, byrow = byrow)
# }

# #' ViolinBoxPlotAtomic
# #'
# #' @inheritParams ViolinBoxPlot
# #' @importFrom rlang sym
# #' @importFrom gglogger ggplot
# #' @importFrom ggplot2 position_dodge2
# #' @keywords internal
# ViolinBoxPlotAtomic <- function(
#     x, y,
#     group_by = NULL, group_by_sep = "_", split_by = NULL, split_by_sep = "_", seed = 8525,
#     facet = FALSE, facet_scales = "fixed", guess_size = FALSE, res = 100,
#     theme = "theme_scp", theme_args = list(), palette = "Paired", palcolor = NULL,
#     keep_empty = FALSE, alpha = 1, x_text_angle = 0, aspect.ratio = 1, title = NULL, subtitle = NULL,
#     xlab = NULL, ylab = NULL, legend.position = "right", legend.direction = "vertical",
#     combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, dodge_position_preserve = "total",
#     mode = "v+b", dodge_by = NULL, dodge_by_sep = "_", return_layer = FALSE,
#     ...
# ) {
#     ngroups <- length(unique(x[[group_by]]))
#     if (!is.null(dodge_by)) {
#         ndodges <- length(unique(x[[dodge_by]]))
#         position <- position_dodge2(preserve = dodge_position_preserve)
#         p <- ggplot(x, aes(x = !!sym(group_by), y = !!sym(y), fill = !!sym(dodge_by))) +
#             scale_fill_manual(
#                 name = dodge_by,
#                 values = palette_scp(levels(x[[dodge_by]]), palette = palette, palcolor = palcolor, keep_names = TRUE)
#             )
#     } else {
#         ndodges <- 1
#         position <- "identity"
#         p <- ggplot(x, aes(x = !!sym(group_by), y = !!sym(y), fill = !!sym(group_by))) +
#             scale_fill_manual(
#                 name = group_by,
#                 values = palette_scp(levels(x[[group_by]]), palette = palette, palcolor = palcolor, keep_names = TRUE)
#             )
#     }

#     half <- grepl("-", mode, fixed = TRUE) && nchar(mode) > 1
#     full_violin <- FALSE
#     fill_box <- FALSE
#     fill_jitter <- FALSE
#     if (grepl("v", mode, fixed = TRUE) && !half) {
#         p <- p + geom_violin(alpha = alpha, position = position)
#         full_violin <- TRUE
#     }
#     if (grepl("b", mode, fixed = TRUE) && !half) {
#         if (isTRUE(full_violin)) {
#             p <- p + geom_boxplot(alpha = alpha, position = position, width = 0.1, fill = "white")
#         } else {
#             p <- p + geom_boxplot(alpha = alpha, position = position, width = 0.8)
#         }
#         fill_box <- TRUE
#     }
#     if (grepl("j", mode, fixed = TRUE) && !half) {
#         p <- p + geom_jitter(alpha = alpha, position = position)
#         fill_jitter <- TRUE
#     }
#     p <- p +
#         labs(title = title, subtitle = subtitle, x = xlab, y = ylab) +
#         do.call(theme, theme_args) +
#         ggplot2::theme(
#             aspect.ratio = aspect.ratio,
#             legend.position = legend.position,
#             legend.direction = legend.direction
#         )

#     if (isTRUE(facet)) {
#         p <- facet_plot(p, facet, facet_scales, split_by, nrow, ncol, byrow)
#     }
#     if (isTRUE(guess_size)) {
#         list(
#             plot = p,
#             height = 5 * res,
#             width = (2 + ngroups * 0.8 * ndodges * 0.6) * res
#         )
#     } else {
#         p
#     }
# }