#' Dimension reduction plot
#'
#' @description Dimension reduction plot for single cell data.
#'
#' @rdname CellDimPlot
#' @inheritParams validate_common_arguments
#'
#' @param x Seurat object
#' @param group_by Name of one or more meta.data columns to group (color) cells by (for example, seurat_clusters).
#' @param reduction Name of the reduction to plot (for example, "umap").
#' @param bg_color Color value for background(NA) points.
#' @param pt_size Point size.
#' @param pt_alpha Point transparency.
#' @param cells_highlight A vector of cell names to highlight or an expression passed to \code{\link{subset}} to select cells to highlight.
#' @param cols_highlight Color used to highlight the cells.
#' @param sizes_highlight Size of highlighted cell points.
#' @param alpha_highlight Transparency of highlighted cell points.
#' @param stroke_highlight Border width of highlighted cell points.
#' @param dims Dimensions to plot, must be a two-length numeric vector specifying x- and y-dimensions
#' @param show_na Whether to assign a color from the color palette to NA group. If \code{FALSE}, cell points with NA level will colored by \code{bg_color}.
#' @param show_stat Whether to show statistical information on the plot.
#' @param label Whether to label the cell groups.
#' @param label_insitu Whether to place the raw labels (group names) in the center of the cells with the corresponding group. Default is \code{FALSE}, which using numbers instead of raw labels.
#' @param label_size Size of labels.
#' @param label_fg Foreground color of label.
#' @param label_bg Background color of label.
#' @param label_bg_r Background ratio of label.
#' @param label_repel Logical value indicating whether the label is repel away from the center points.
#' @param label_repulsion Force of repulsion between overlapping text labels. Defaults to 20.
#' @param label_point_size Size of the center points.
#' @param label_point_color Color of the center points.
#' @param label_segment_color Color of the line segment for labels.
#' @param add_mark Whether to add a mark layer on the plot.
#' @param mark_type Type of mark to add to the plot. Default is "hull".
#' @param mark_expand Expansion of the mark.
#' @param mark_alpha Transparency of the mark.
#' @param mark_linetype Line type of the mark.
#' @param add_density Whether to add a density layer on the plot.
#' @param density_color Color of the density contours lines.
#' @param density_filled Whether to add filled contour bands instead of contour lines.
#' @param density_filled_palette Color palette used to fill contour bands.
#' @param density_filled_palcolor Custom colors used to fill contour bands.
#' @param lineages Lineages/pseudotime to add to the plot. If specified, curves will be fitted using \code{\link[stats]{loess}} method.
#' @param lineages_trim Trim the leading and the trailing data in the lineages.
#' @param lineages_span The parameter α which controls the degree of smoothing in \code{\link[stats]{loess}} method.
#' @param lineages_palette Color palette used for lineages.
#' @param lineages_palcolor Custom colors used for lineages.
#' @param lineages_arrow Set arrows of the lineages. See \code{\link[grid]{arrow}}.
#' @param lineages_linewidth Width of fitted curve lines for lineages.
#' @param lineages_line_bg Background color of curve lines for lineages.
#' @param lineages_line_bg_stroke Border width of curve lines background.
#' @param lineages_whiskers Whether to add whiskers for lineages.
#' @param lineages_whiskers_linewidth Width of whiskers for lineages.
#' @param lineages_whiskers_alpha Transparency of whiskers for lineages.
#' @param stat_by The name of a metadata column to stat.
#' @param stat_frac The way of calculating the fraction. Default is "none".
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
#' @param stat_pie_label Whether to add labels in the pie plot.
#' @param stat_plot_type Set the statistical plot type.
#' @param stat_plot_size Set the statistical plot size. Defaults to 0.1
#' @param stat_plot_palette Color palette used in statistical plot.
#' @param stat_palcolor Custom colors used in statistical plot
#' @param stat_plot_position Position adjustment in statistical plot.
#' @param stat_plot_alpha Transparency of the statistical plot.
#' @param stat_plot_label Whether to add labels in the statistical plot.
#' @param stat_plot_label_size Label size in the statistical plot.
#' @param graph Specify the graph name to add edges between cell neighbors to the plot.
#' @param edge_size Size of edges.
#' @param edge_alpha Transparency of edges.
#' @param edge_color Color of edges.
#' @param paga Specify the calculated paga results to add a PAGA graph layer to the plot.
#' @param paga_type PAGA plot type. "connectivities" or "connectivities_tree".
#' @param paga_node_size Size of the nodes in PAGA plot.
#' @param paga_edge_threshold Threshold of edge connectivities in PAGA plot.
#' @param paga_edge_size Size of edges in PAGA plot.
#' @param paga_edge_color Color of edges in PAGA plot.
#' @param paga_edge_alpha Transparency of edges in PAGA plot.
#' @param paga_show_transition Whether to show transitions between edges.
#' @param paga_transition_threshold Threshold of transition edges in PAGA plot.
#' @param paga_transition_size Size of transition edges in PAGA plot.
#' @param paga_transition_color Color of transition edges in PAGA plot.
#' @param paga_transition_alpha Transparency of transition edges in PAGA plot.
#' @param velocity Specify the calculated RNA velocity mode to add a velocity layer to the plot.
#' @param velocity_plot_type Set the velocity plot type.
#' @param velocity_n_neighbors Set the number of neighbors used in velocity plot.
#' @param velocity_density Set the density value used in velocity plot.
#' @param velocity_smooth Set the smooth value used in velocity plot.
#' @param velocity_scale Set the scale value used in velocity plot.
#' @param velocity_min_mass Set the min_mass value used in velocity plot.
#' @param velocity_cutoff_perc Set the cutoff_perc value used in velocity plot.
#' @param velocity_arrow_color Color of arrows in velocity plot.
#' @param velocity_arrow_angle Angle of arrows in velocity plot.
#' @param streamline_L Typical length of a streamline in x and y units
#' @param streamline_minL Minimum length of segments to show.
#' @param streamline_res Resolution parameter (higher numbers increases the resolution).
#' @param streamline_n Number of points to draw.
#' @param streamline_width Size of streamline.
#' @param streamline_alpha Transparency of streamline.
#' @param streamline_color Color of streamline.
#' @param streamline_palette Color palette used for streamline.
#' @param streamline_palcolor Custom colors used for streamline.
#' @param streamline_bg_color Background color of streamline.
#' @param streamline_bg_stroke Border width of streamline background.
#' @param hex Whether to chane the plot type from point to the hexagonal bin.
#' @param hex_count Whether show cell counts in each hexagonal bin.
#' @param hex_bins Number of hexagonal bins.
#' @param hex_binwidth Hexagonal bin width.
#' @param hex_linewidth Border width of hexagonal bins.
#' @param raster Convert points to raster format, default is NULL which automatically rasterizes if plotting more than 100,000 cells
#' @param raster_dpi Pixel resolution for rasterized plots, passed to geom_scattermore(). Default is c(512, 512).
#' @param ... Additional arguments
#' @return A ggplot object or a list with the plot and the height and width of the plot if guess_size is TRUE
#' @export
#'
#' @importFrom Seurat Reductions Embeddings Key SplitObject GetAssay
#' @importFrom SeuratObject DefaultDimReduc
#' @importFrom rlang enexpr as_label
#'
#' @examples
#' data("pancreas_sub")
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP")
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             theme = "theme_blank")
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             palette = "seurat", theme = "theme_blank")
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             theme = ggplot2::theme_classic, theme_args = list(base_size = 16))
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             raster = TRUE, raster_dpi = 30)
#'
#' # Highlight cells
#' CellDimPlot(pancreas_sub,
#'   group_by = "SubCellType", reduction = "UMAP",
#'   cells_highlight = SubCellType == "Epsilon"
#' )
#' CellDimPlot(pancreas_sub,
#'   group_by = "SubCellType", split_by = "Phase", reduction = "UMAP",
#'   cells_highlight = TRUE, theme = "theme_blank", legend.position = "none"
#' )
#' CellDimPlot(pancreas_sub,
#'   group_by = "SubCellType", split_by = "Phase", reduction = "UMAP",
#'   cells_highlight = TRUE, theme = "theme_blank", legend.position = "none",
#'   facet = TRUE
#' )
#'
#' # Add group labels
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             label = TRUE)
#' CellDimPlot(pancreas_sub,
#'   group_by = "SubCellType", reduction = "UMAP",
#'   label = TRUE, label_fg = "orange", label_bg = "red", label_size = 5
#' )
#' CellDimPlot(pancreas_sub,
#'   group_by = "SubCellType", reduction = "UMAP",
#'   label = TRUE, label_insitu = TRUE
#' )
#' CellDimPlot(pancreas_sub,
#'   group_by = "SubCellType", reduction = "UMAP",
#'   label = TRUE, label_insitu = TRUE, label_repel = TRUE,
#'   label_segment_color = "red"
#' )
#'
#' # Add various shape of marks
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             add_mark = TRUE)
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             add_mark = TRUE, mark_expand = unit(1, "mm"))
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             add_mark = TRUE, mark_alpha = 0.3)
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             add_mark = TRUE, mark_linetype = 2)
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             add_mark = TRUE, mark_type = "ellipse")
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             add_mark = TRUE, mark_type = "rect")
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             add_mark = TRUE, mark_type = "circle")
#'
#' # Add a density layer
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             add_density = TRUE)
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             add_density = TRUE, density_filled = TRUE)
#' CellDimPlot(pancreas_sub,
#'   group_by = "SubCellType", reduction = "UMAP",
#'   add_density = TRUE, density_filled = TRUE, density_filled_palette = "Blues",
#'   cells_highlight = TRUE
#' )
#'
#' # Add statistical charts
#' CellDimPlot(pancreas_sub,
#'   group_by = "CellType", reduction = "UMAP", stat_by = "Phase")
#' CellDimPlot(pancreas_sub,
#'   group_by = "CellType", reduction = "UMAP", stat_by = "Phase",
#'   stat_plot_type = "ring", stat_plot_label = TRUE, stat_plot_size = 0.15)
#' CellDimPlot(pancreas_sub,
#'   group_by = "CellType", reduction = "UMAP", stat_by = "Phase",
#'   stat_plot_type = "bar", stat_type = "count", stat_plot_position = "dodge")
#'
#' # Chane the plot type from point to the hexagonal bin
#' CellDimPlot(pancreas_sub, group_by = "CellType", reduction = "UMAP",
#'             hex = TRUE)
#' CellDimPlot(pancreas_sub, group_by = "CellType", reduction = "UMAP",
#'             hex = TRUE, hex_bins = 20)
#' CellDimPlot(pancreas_sub, group_by = "CellType", reduction = "UMAP",
#'             hex = TRUE, hex_count = FALSE)
#'
#' # Show neighbors graphs on the plot
#' pancreas_sub <- Standard_SCP(pancreas_sub)
#' CellDimPlot(pancreas_sub, group_by = "CellType", reduction = "UMAP",
#'             graph = "RNA_nn")
#' CellDimPlot(pancreas_sub, group_by = "CellType", reduction = "UMAP",
#'             graph = "RNA_snn", edge_color = "grey80")
#'
#' # Show lineages on the plot based on the pseudotime
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             lineages = paste0("Lineage", 1:3))
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             lineages = paste0("Lineage", 1:3), lineages_whiskers = TRUE)
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             lineages = paste0("Lineage", 1:3), lineages_span = 0.1)
CellDimPlot <- function(
    x,
    # common arguments
    group_by = NULL, group_by_sep = "_", split_by = NULL, split_by_sep = "_",
    seed = 8525, facet = FALSE, guess_size = FALSE, res = 100, facet_scales = "fixed",
    theme = "theme_scp", theme_args = list(), palette = "Paired", palcolor = NULL,
    aspect.ratio = 1, legend.position = "right", legend.direction = "vertical",
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL, nrow = NULL, ncol = NULL, byrow = TRUE,
    # specific arguments
    reduction = NULL, bg_color = "grey80", pt_size = NULL, pt_alpha = 1, add_mark = FALSE, mark_type = c("ellipse", "hull", "rect", "circle"), stat_pie_label = FALSE,
    mark_expand = unit(3, "mm"), mark_alpha = 0.1, mark_linetype = 1,
    cells_highlight = NULL, cols_highlight = "black", sizes_highlight = 1, alpha_highlight = 1, stroke_highlight = 0.5,
    dims = c(1, 2), show_na = FALSE, show_stat = !identical(theme, "theme_blank") && !facet,
    label = FALSE, label_insitu = FALSE, label_size = 4, label_fg = "white", label_bg = "black", label_bg_r = 0.1,
    label_repel = FALSE, label_repulsion = 20, label_point_size = 1, label_point_color = "black", label_segment_color = "black",
    add_density = FALSE, density_color = "grey80", density_filled = FALSE, density_filled_palette = "Greys", density_filled_palcolor = NULL,
    lineages = NULL, lineages_trim = c(0.01, 0.99), lineages_span = 0.75, lineages_palette = "Dark2", lineages_palcolor = NULL,
    lineages_arrow = arrow(length = unit(0.1, "inches")), lineages_linewidth = 1, lineages_line_bg = "white", lineages_line_bg_stroke = 0.5,
    lineages_whiskers = FALSE, lineages_whiskers_linewidth = 0.5, lineages_whiskers_alpha = 0.5,
    stat_by = NULL, stat_frac = "none", stat_plot_type = "pie", stat_plot_size = 0.1, stat_plot_palette = "Set1", stat_palcolor = NULL,
    stat_plot_position = c("stack", "dodge"), stat_plot_alpha = 1, stat_plot_label = FALSE, stat_plot_label_size = 3,
    graph = NULL, edge_size = c(0.05, 0.5), edge_alpha = 0.1, edge_color = "grey40",
    paga = NULL, paga_type = "connectivities", paga_node_size = 4, paga_edge_threshold = 0.01, paga_edge_size = c(0.2, 1),
    paga_edge_color = "grey40", paga_edge_alpha = 0.5, paga_show_transition = FALSE, paga_transition_threshold = 0.01,
    paga_transition_size = c(0.2, 1), paga_transition_color = "black", paga_transition_alpha = 1,
    velocity = NULL, velocity_plot_type = "raw", velocity_n_neighbors = ceiling(ncol(GetAssay(x)) / 50),
    velocity_density = 1, velocity_smooth = 0.5, velocity_scale = 1, velocity_min_mass = 1, velocity_cutoff_perc = 5,
    velocity_arrow_color = "black", velocity_arrow_angle = 20,
    streamline_L = 5, streamline_minL = 1, streamline_res = 1, streamline_n = 15,
    streamline_width = c(0, 0.8), streamline_alpha = 1, streamline_color = NULL, streamline_palette = "RdYlBu", streamline_palcolor = NULL,
    streamline_bg_color = "white", streamline_bg_stroke = 0.5,
    hex = FALSE, hex_count = TRUE, hex_bins = 50, hex_binwidth = NULL, hex_linewidth = 0.5,
    raster = NULL, raster_dpi = c(512, 512), combine = TRUE,
    ...) {
    set.seed(seed)
    # Check arguments ----------------------------------------------------------
    validate_common_arguments(
        split_by = split_by, seed = seed, res = res,
        guess_size = guess_size, facet = facet
    )
    if ("group.by" %in% names(list(...))) {
        stop("Please use 'group_by' instead of 'group.by'")
    }
    if ("split.by" %in% names(list(...))) {
        stop("Please use 'split_by' instead of 'split.by'")
    }
    if ("stat_by" %in% names(list(...))) {
        stop("Please use 'stat_by' instead of 'stat_by'")
    }

    mark_type <- match.arg(mark_type)

    if (!is.null(lineages)) {
        non_exist <- setdiff(lineages, colnames(x@meta.data))
        if (length(non_exist) > 0) {
            stop(paste0("Lineages ", paste(non_exist, collapse = ", "), " not found in the meta.data of seurat object."))
        }
    }

    if (!is.null(stat_by) && !stat_by %in% colnames(x@meta.data)) {
        stop(paste0("Stat by ", stat_by, " not found in the meta.data of seurat object."))
    }

    if (is.null(group_by)) {
        group_by <- ".group_by"
        x@meta.data$.group_by <- factor("")
    } else {
        if (length(group_by) > 1) {
            x@meta.data <- concat_cols(x@meta.data, group_by, group_by_sep)
            group_by <- attr(x@meta.data, "new_col")
        }
    }

    if (isTRUE(show_na) && any(is.na(x@meta.data[[group_by]]))) {
        raw_levels <- unique(c(levels(x@meta.data[[group_by]]), "NA"))
        x@meta.data[[group_by]] <- as.character(x@meta.data[[groupo_by]])
        x@meta.data[[group_by]][is.na(x@meta.data[[group_by]])] <- "NA"
        x@meta.data[[group_by]] <- factor(x@meta.data[[group_by]], levels = raw_levels)
    }
    if (!is.null(graph) && !graph %in% names(x@graphs)) {
        stop("Graph ", graph, " is not exist in the seurat object.")
    }
    if (!is.null(graph)) {
        graph <- x@graphs[[graph]]
    }
    if (is.null(reduction)) {
        reduction <- DefaultDimReduc(x)
    }
    if (!reduction %in% names(x@reductions)) {
        stop(paste0(reduction, " is not in the x reduction names."))
    }
    if (!is.null(x = raster_dpi)) {
        if (!is.numeric(x = raster_dpi)) {
            stop("'raster_dpi' must be a two-length numeric vector")
        }
        if (length(raster_dpi) == 1) {
            raster_dpi <- rep(raster_dpi, 2)
        }
    }
    nocells_msg <- "No cells in 'cells_highlight' found in the seurat object."
    tryCatch({
        if (is.null(cells_highlight)) {
            cells_highlight <- NULL
        } else if (isTRUE(cells_highlight)) {
            cells_highlight <- TRUE
        } else {
            cells_names <- colnames(GetAssay(x))
            if (!any(cells_highlight %in% cells_names)) {
                stop(nocells_msg)
            }
            if (!all(cells_highlight %in% cells_names)) {
                warning("Some cells in 'cells_highlight' not found in seurat object.", immediate. = TRUE)
            }
            cells_highlight <- intersect(cells_highlight, cells_names)
            rm(cells_names)
        }
    }, error = function(e) {
        if (identical(e$message, nocells_msg)) {
            stop(nocells_msg)
        }
        cells_highlight <- as_label(enexpr(cells_highlight))
        x_highlight <- eval(parse(text = paste0("subset(x, subset = ", cells_highlight, ")")))
        cells_highlight <<- colnames(GetAssay(x_highlight))
        rm(x_highlight)
    })

    # Split data --------------------------------------------------------------
    if (is.null(split_by) || isTRUE(facet)) {
        xs <- list(x)
    } else if (isFALSE(facet)){
        if (length(split_by) > 1) {
            x@meta.data <- concat_cols(x@meta.data, split_by, split_by_sep)
            split_by <- attr(x@meta.data, "new_col")
        } else if (!is.factor(x@meta.data[[split_by]])) {
            x@meta.data[[split_by]] <- factor(
                x@meta.data[[split_by]],
                levels = unique(x@meta.data[[split_by]])
            )
        }
        xs <- SplitObject(x, split.by = split_by)
        xs <- xs[levels(x@meta.data[[split_by]])]
    }
    rm(x)

    # Plot for each split data ------------------------------------------------
    plots <- lapply(xs, function(s) {
        CellDimPlotAtomic(
            s,
            # common arguments
            group_by = group_by, split_by = split_by, seed = seed, facet = facet, guess_size = guess_size, res = res, facet_scales = facet_scales,
            theme = theme, theme_args = theme_args, palette = palette, palcolor = palcolor, aspect.ratio = aspect.ratio,
            legend.position = legend.position, legend.direction = legend.direction, nrow = nrow, ncol = ncol, byrow = byrow,
            title = title, subtitle = subtitle, xlab = xlab, ylab = ylab,
            # specific arguments
            reduction = reduction, bg_color = bg_color, pt_size = pt_size, pt_alpha = pt_alpha, mark_type = mark_type,
            cells_highlight = cells_highlight, cols_highlight = cols_highlight, sizes_highlight = sizes_highlight, alpha_highlight = alpha_highlight,
            stroke_highlight = stroke_highlight, add_mark = add_mark, mark_expand = mark_expand, mark_alpha = mark_alpha, mark_linetype = mark_linetype,
            dims = dims, show_na = show_na, show_stat = show_stat, label = label, label_insitu = label_insitu, label_size = label_size,
            label_fg = label_fg, label_bg = label_bg, label_bg_r = label_bg_r, label_repel = label_repel, label_repulsion = label_repulsion,
            label_point_size = label_point_size, label_point_color = label_point_color, label_segment_color = label_segment_color, stat_pie_label = stat_pie_label,
            add_density = add_density, density_color = density_color, density_filled = density_filled, density_filled_palette = density_filled_palette,
            density_filled_palcolor = density_filled_palcolor, lineages = lineages, lineages_trim = lineages_trim, lineages_span = lineages_span,
            lineages_palette = lineages_palette, lineages_palcolor = lineages_palcolor, lineages_arrow = lineages_arrow, lineages_linewidth = lineages_linewidth,
            lineages_line_bg = lineages_line_bg, lineages_line_bg_stroke = lineages_line_bg_stroke, lineages_whiskers = lineages_whiskers,
            lineages_whiskers_linewidth = lineages_whiskers_linewidth, lineages_whiskers_alpha = lineages_whiskers_alpha, stat_by = stat_by,
            stat_frac = stat_frac, stat_plot_type = stat_plot_type, stat_plot_size = stat_plot_size, stat_plot_palette = stat_plot_palette,
            stat_palcolor = stat_palcolor, stat_plot_position = stat_plot_position, stat_plot_alpha = stat_plot_alpha, stat_plot_label = stat_plot_label,
            stat_plot_label_size = stat_plot_label_size, graph = graph, edge_size = edge_size, edge_alpha = edge_alpha, edge_color = edge_color,
            paga = paga, paga_type = paga_type, paga_node_size = paga_node_size, paga_edge_threshold = paga_edge_threshold,
            paga_edge_size = paga_edge_size, paga_edge_color = paga_edge_color, paga_edge_alpha = paga_edge_alpha, paga_show_transition = paga_show_transition,
            paga_transition_threshold = paga_transition_threshold, paga_transition_size = paga_transition_size, paga_transition_color = paga_transition_color,
            paga_transition_alpha = paga_transition_alpha, velocity = velocity, velocity_plot_type = velocity_plot_type, velocity_n_neighbors = velocity_n_neighbors,
            velocity_density = velocity_density, velocity_smooth = velocity_smooth, velocity_scale = velocity_scale, velocity_min_mass = velocity_min_mass,
            velocity_cutoff_perc = velocity_cutoff_perc, velocity_arrow_color = velocity_arrow_color, velocity_arrow_angle = velocity_arrow_angle,
            streamline_L = streamline_L, streamline_minL = streamline_minL, streamline_res = streamline_res, streamline_n = streamline_n,
            streamline_width = streamline_width, streamline_alpha = streamline_alpha, streamline_color = streamline_color, streamline_palette = streamline_palette,
            streamline_palcolor = streamline_palcolor, streamline_bg_color = streamline_bg_color, streamline_bg_stroke = streamline_bg_stroke,
            hex = hex, hex_count = hex_count, hex_bins = hex_bins, hex_binwidth = hex_binwidth, hex_linewidth = hex_linewidth,
            raster = raster, raster_dpi = raster_dpi,
            ...
        )
    })

    combine_plots(plots, combine = combine, nrow = nrow, ncol = ncol, byrow = byrow)
}

#' CellDimPlotAtomic
#'
#' @inheritParams CellDimPlot
#' @keywords internal
#' @importFrom gglogger ggplot
#' @importFrom cowplot as_grob get_plot_component
#' @importFrom ggplot2 aes geom_point geom_density_2d stat_density_2d geom_segment labs scale_x_continuous scale_y_continuous scale_size_continuous facet_grid scale_color_manual scale_fill_manual guides guide_legend geom_hex geom_path theme_void annotation_custom scale_linewidth_continuous after_stat layer_scales
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggnewscale new_scale_color new_scale_fill new_scale
#' @importFrom gtable gtable_add_cols gtable_add_grob
#' @importFrom ggforce geom_mark_hull geom_mark_ellipse geom_mark_circle geom_mark_rect
#' @importFrom grid arrow unit
#' @importFrom patchwork wrap_plots
#' @importFrom stats median loess aggregate
CellDimPlotAtomic <- function(
    x,
    # common arguments
    group_by = NULL, group_by_sep = "_", split_by = NULL, split_by_sep = "_",
    seed = 8525, facet = FALSE, guess_size = FALSE, res = 100, facet_scales = "fixed",
    theme = "theme_scp", theme_args = list(), palette = "Paired", palcolor = NULL,
    aspect.ratio = 1, legend.position = "right", legend.direction = "vertical",
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL, nrow = NULL, ncol = NULL, byrow = TRUE,
    # specific arguments
    reduction = NULL, bg_color = "grey80", pt_size = NULL, pt_alpha = 1, add_mark = FALSE, mark_type = c("ellipse", "hull", "rect", "circle"),
    cells_highlight = NULL, cols_highlight = "black", sizes_highlight = 1, alpha_highlight = 1, stroke_highlight = 0.5,
    dims = c(1, 2), show_na = FALSE, show_stat = !identical(theme, "theme_blank") && !facet,
    mark_expand = unit(3, "mm"), mark_alpha = 0.1, mark_linetype = 1, stat_pie_label = FALSE,
    label = FALSE, label_insitu = FALSE, label_size = 4, label_fg = "white", label_bg = "black", label_bg_r = 0.1,
    label_repel = FALSE, label_repulsion = 20, label_point_size = 1, label_point_color = "black", label_segment_color = "black",
    add_density = FALSE, density_color = "grey80", density_filled = FALSE, density_filled_palette = "Greys", density_filled_palcolor = NULL,
    lineages = NULL, lineages_trim = c(0.01, 0.99), lineages_span = 0.75, lineages_palette = "Dark2", lineages_palcolor = NULL,
    lineages_arrow = arrow(length = unit(0.1, "inches")), lineages_linewidth = 1, lineages_line_bg = "white", lineages_line_bg_stroke = 0.5,
    lineages_whiskers = FALSE, lineages_whiskers_linewidth = 0.5, lineages_whiskers_alpha = 0.5,
    stat_by = NULL, stat_frac = "none", stat_plot_type = "pie", stat_plot_size = 0.15, stat_plot_palette = "Set1", stat_palcolor = NULL,
    stat_plot_position = c("stack", "dodge"), stat_plot_alpha = 1, stat_plot_label = FALSE, stat_plot_label_size = 3,
    graph = NULL, edge_size = c(0.05, 0.5), edge_alpha = 0.1, edge_color = "grey40",
    paga = NULL, paga_type = "connectivities", paga_node_size = 4, paga_edge_threshold = 0.01, paga_edge_size = c(0.2, 1),
    paga_edge_color = "grey40", paga_edge_alpha = 0.5, paga_show_transition = FALSE, paga_transition_threshold = 0.01,
    paga_transition_size = c(0.2, 1), paga_transition_color = "black", paga_transition_alpha = 1,
    velocity = NULL, velocity_plot_type = "raw", velocity_n_neighbors = ceiling(ncol(x@assays[[1]]) / 50),
    velocity_density = 1, velocity_smooth = 0.5, velocity_scale = 1, velocity_min_mass = 1, velocity_cutoff_perc = 5,
    velocity_arrow_color = "black", velocity_arrow_angle = 20,
    streamline_L = 5, streamline_minL = 1, streamline_res = 1, streamline_n = 15,
    streamline_width = c(0, 0.8), streamline_alpha = 1, streamline_color = NULL, streamline_palette = "RdYlBu", streamline_palcolor = NULL,
    streamline_bg_color = "white", streamline_bg_stroke = 0.5,
    hex = FALSE, hex_count = TRUE, hex_bins = 50, hex_binwidth = NULL, hex_linewidth = 0.5,
    raster = NULL, raster_dpi = c(512, 512),
    ...) {
    reduction_key <- x@reductions[[reduction]]@key
    # keep split_by for facetting
    dat_meta <- x@meta.data[, unique(c(group_by, split_by)), drop = FALSE]
    dat_dim <- x@reductions[[reduction]]@cell.embeddings
    colnames(dat_dim) <- paste0(reduction_key, seq_len(ncol(dat_dim)))
    rownames(dat_dim) <- rownames(dat_dim) %||% colnames(x@assays[[1]])
    dat <- cbind(dat_dim, dat_meta[row.names(dat_dim), , drop = FALSE])

    if (is.null(pt_size)) {
        pt_size <- min(3000 / nrow(dat), 0.5)
    }
    raster <- raster %||% (nrow(dat) > 1e5)
    if (isTRUE(raster)) {
        requireNamespace("scattermore", quietly = TRUE)
    }
    if (!is.null(stat_by)) {
        if (isTRUE(facet)) {
            stop('Cannot facet a CellDimPlot with "stat_by"')
        }
        subplots <- CellStatPlot(
            x,
            ident = stat_by, split_by = group_by, pie_label = stat_pie_label,
            frac = stat_frac, plot_type = stat_plot_type, position = stat_plot_position,
            palette = stat_plot_palette, palcolor = stat_palcolor, alpha = stat_plot_alpha,
            label = stat_plot_label, label_size = stat_plot_label_size,
            legend.position = "bottom", legend.direction = legend.direction,
            theme = theme, theme_args = theme_args, combine = FALSE
        )
    }
    if (!is.null(lineages)) {
        if (isTRUE(facet)) {
            stop('Cannot facet a CellDimPlot with "lineages"')
        }
        lineages_layers <- LineagePlot(
            x, lineages = lineages, reduction = reduction, dims = dims,
            trim = lineages_trim, span = lineages_span,
            palette = lineages_palette, palcolor = lineages_palcolor, lineages_arrow = lineages_arrow,
            linewidth = lineages_linewidth, line_bg = lineages_line_bg, line_bg_stroke = lineages_line_bg_stroke,
            whiskers = lineages_whiskers, whiskers_linewidth = lineages_whiskers_linewidth, whiskers_alpha = lineages_whiskers_alpha,
            aspect.ratio = aspect.ratio, title = title, subtitle = subtitle, xlab = xlab, ylab = ylab,
            legend.position = "bottom", legend.direction = legend.direction,
            theme = theme, theme_args = theme_args,
            return_layer = TRUE
        )
        lineages_layers <- lineages_layers[!names(lineages_layers) %in% c("lab_layer", "theme_layer")]
    }
    if (!is.null(paga)) {
        if (facet) {
            stop('Cannot facet a CellDimPlot with "paga"')
        }
        stop("PAGA plot is not supported in CellDimPlot yet.")
        # paga_layers <- PAGAPlot(x,
        #     paga = paga, type = paga_type, reduction = reduction, dims = dims,
        #     node_palette = palette, node_palcolor = palcolor, node_size = paga_node_size,
        #     edge_threshold = paga_edge_threshold, edge_size = paga_edge_size, edge_color = paga_edge_color, edge_alpha = paga_edge_alpha,
        #     transition_threshold = paga_transition_threshold, transition_size = paga_transition_size, transition_color = paga_transition_color, transition_alpha = paga_transition_alpha, show_transition = paga_show_transition,
        #     aspect.ratio = aspect.ratio, title = title, subtitle = subtitle, xlab = xlab, ylab = ylab,
        #     legend.position = "bottom", legend.direction = legend.direction,
        #     theme = theme, theme_args = theme_args,
        #     return_layer = TRUE
        # )
        # paga_layers <- paga_layers[!names(paga_layers) %in% c("lab_layer", "theme_layer")]
    }
    if (!is.null(velocity)) {
        if (facet) {
            stop('Cannot facet a CellDimPlot with "velocity"')
        }
        stop("Velocity plot is not supported in CellDimPlot yet.")
        # velocity_layers <- VelocityPlot(x,
        #     reduction = reduction, dims = dims, velocity = velocity, plot_type = velocity_plot_type, group_by = group_by, group_palette = palette, group_palcolor = palcolor,
        #     n_neighbors = velocity_n_neighbors, density = velocity_density, smooth = velocity_smooth, scale = velocity_scale, min_mass = velocity_min_mass, cutoff_perc = velocity_cutoff_perc,
        #     arrow_color = velocity_arrow_color, arrow_angle = velocity_arrow_angle,
        #     streamline_L = streamline_L, streamline_minL = streamline_minL, streamline_res = streamline_res, streamline_n = streamline_n,
        #     streamline_width = streamline_width, streamline_alpha = streamline_alpha, streamline_color = streamline_color, streamline_palette = streamline_palette, streamline_palcolor = streamline_palcolor,
        #     streamline_bg_color = streamline_bg_color, streamline_bg_stroke = streamline_bg_stroke,
        #     aspect.ratio = aspect.ratio, title = title, subtitle = subtitle, xlab = xlab, ylab = ylab,
        #     legend.position = "bottom", legend.direction = legend.direction,
        #     theme = theme_void, theme_args = theme_args,
        #     return_layer = TRUE
        # )
        # velocity_layers <- velocity_layers[!names(velocity_layers) %in% c("lab_layer", "theme_layer")]
    }
    xlab <- xlab %||% paste0(reduction_key, dims[1])
    ylab <- ylab %||% paste0(reduction_key, dims[2])
    if (identical(theme, "theme_blank")) {
        theme_args[["xlab"]] <- xlab
        theme_args[["ylab"]] <- ylab
    }
    colors <- palette_scp(levels(dat[[group_by]]), palette = palette, palcolor = palcolor, NA_keep = TRUE, keep_names = TRUE)
    legend_list <- list()
    labels_tb <- table(dat[[group_by]])
    labels_tb <- labels_tb[labels_tb != 0]
    cells.highlight_use <- cells_highlight
    if (isTRUE(cells.highlight_use)) {
        cells.highlight_use <- rownames(dat)[!is.na(dat[[group_by]])]
    }
    if (isTRUE(label_insitu)) {
        if (isTRUE(show_stat)) {
            label_use <- paste0(names(labels_tb), "(", labels_tb, ")")
        } else {
            label_use <- paste0(names(labels_tb))
        }
    } else {
        if (isTRUE(label)) {
            if (isTRUE(show_stat)) {
                label_use <- paste0(seq_along(labels_tb), ": ", names(labels_tb), "(", labels_tb, ")")
            } else {
                label_use <- paste0(seq_along(labels_tb), ": ", names(labels_tb))
            }
        } else {
            if (isTRUE(show_stat)) {
                label_use <- paste0(names(labels_tb), "(", labels_tb, ")")
            } else {
                label_use <- paste0(names(labels_tb))
            }
        }
    }
    dat$x <- dat[[paste0(reduction_key, dims[1])]]
    dat$y <- dat[[paste0(reduction_key, dims[2])]]
    if (isTRUE(facet) && !is.null(split_by)) {
        split_values <- unique(dat[[split_by]])
        # Add cells from other split to show them as grey
        if (length(split_values) > 1) {
            for (sp in split_values) {
                rest_cells <- dat[dat[[split_by]] != sp, , drop = FALSE]
                rest_cells[[split_by]] <- sp
                rest_cells[[group_by]] <- "NA"
                dat <- rbind(dat, rest_cells)
            }
            rm(rest_cells)
        }
    }
    dat$group_by <- dat[[group_by]]

    dat <- dat[order(dat$group_by, decreasing = FALSE, na.last = FALSE), , drop = FALSE]
    naindex <- which(is.na(dat$group_by))
    naindex <- ifelse(length(naindex) > 0, max(naindex), 1)
    dat <- dat[c(1:naindex, sample((min(naindex + 1, nrow(dat))):nrow(dat))), , drop = FALSE]

    sp <- ifelse(isFALSE(facet) && !is.null(split_by), as.character(dat[[split_by]])[1], "")
    if (isTRUE(show_stat)) {
        subtitle_use <- subtitle %||% paste0(sp, " nCells:", sum(!is.na(dat$group_by)))
    } else {
        subtitle_use <- subtitle %||% sp
    }

    if (isTRUE(add_mark)) {
        mark_fun <- switch(mark_type,
            "ellipse" = "geom_mark_ellipse",
            "hull" = "geom_mark_hull",
            "rect" = "geom_mark_rect",
            "circle" = "geom_mark_circle"
        )
        mark <- list(
            do.call(
                mark_fun,
                list(
                    data = dat[!is.na(dat$group_by), , drop = FALSE],
                    mapping = aes(x = x, y = y, color = group_by, fill = group_by),
                    expand = mark_expand, alpha = mark_alpha, linetype = mark_linetype, show.legend = FALSE, inherit.aes = FALSE
                ),
            ),
            scale_fill_manual(values = colors[names(labels_tb)]),
            scale_color_manual(values = colors[names(labels_tb)]),
            new_scale_fill(),
            new_scale_color()
        )
    } else {
        mark <- NULL
    }

    if (!is.null(graph)) {
        if (isTRUE(facet)) {
            stop('Cannot facet a CellDimPlot with "graph"')
        }
        net_mat <- as.matrix(graph)[rownames(dat), rownames(dat)]
        net_mat[net_mat == 0] <- NA
        net_mat[upper.tri(net_mat)] <- NA
        net_df <- reshape2::melt(net_mat, na.rm = TRUE, stringsAsFactors = FALSE)
        net_df[, "value"] <- as.numeric(net_df[, "value"])
        net_df[, "Var1"] <- as.character(net_df[, "Var1"])
        net_df[, "Var2"] <- as.character(net_df[, "Var2"])
        net_df[, "x"] <- dat[net_df[, "Var1"], "x"]
        net_df[, "y"] <- dat[net_df[, "Var1"], "y"]
        net_df[, "xend"] <- dat[net_df[, "Var2"], "x"]
        net_df[, "yend"] <- dat[net_df[, "Var2"], "y"]
        net <- list(
            geom_segment(
                data = net_df, mapping = aes(x = x, y = y, xend = xend, yend = yend, linewidth = value),
                color = edge_color, alpha = edge_alpha, show.legend = FALSE
            ),
            scale_linewidth_continuous(range = edge_size)
        )
    } else {
        net <- NULL
    }

    if (isTRUE(add_density)) {
        if (isTRUE(density_filled)) {
            filled_color <- palette_scp(palette = density_filled_palette, palcolor = density_filled_palcolor)
            density <- list(
                stat_density_2d(
                    geom = "raster", aes(x = x, y = y, fill = after_stat(density)),
                    contour = FALSE, inherit.aes = FALSE, show.legend = FALSE
                ),
                scale_fill_gradientn(name = "Density", colours = filled_color),
                new_scale_fill()
            )
        } else {
            density <- geom_density_2d(aes(x = x, y = y),
                color = density_color, inherit.aes = FALSE, show.legend = FALSE
            )
        }
    } else {
        density <- NULL
    }

    p <- ggplot(dat) +
        mark +
        net +
        density +
        labs(title = title, subtitle = subtitle_use, x = xlab, y = ylab) +
        scale_x_continuous(limits = c(min(dat_dim[, paste0(reduction_key, dims[1])], na.rm = TRUE), max(dat_dim[, paste0(reduction_key, dims[1])], na.rm = TRUE))) +
        scale_y_continuous(limits = c(min(dat_dim[, paste0(reduction_key, dims[2])], na.rm = TRUE), max(dat_dim[, paste0(reduction_key, dims[2])], na.rm = TRUE))) +
        do.call(theme, theme_args) +
        ggplot2::theme(
            aspect.ratio = aspect.ratio,
            legend.position = legend.position,
            legend.direction = legend.direction
        )

    if (isTRUE(raster)) {
        p <- p + scattermore::geom_scattermore(
            data = dat[is.na(dat$group_by), , drop = FALSE],
            mapping = aes(x = x, y = y), color = bg_color,
            pointsize = ceiling(pt_size), alpha = pt_alpha, pixels = raster_dpi
        ) + scattermore::geom_scattermore(
            data = dat[!is.na(dat$group_by), , drop = FALSE],
            mapping = aes(x = x, y = y, color = group_by, fill = group_by),
            pointsize = ceiling(pt_size), alpha = pt_alpha, pixels = raster_dpi
        )
    } else if (isTRUE(hex)) {
        requireNamespace("hexbin", quietly = TRUE)
        if (isTRUE(hex_count)) {
            p <- p + geom_hex(
                mapping = aes(x = x, y = y, fill = group_by, color = group_by, alpha = after_stat(count)),
                linewidth = hex_linewidth, bins = hex_bins, binwidth = hex_binwidth
            )
        } else {
            p <- p + geom_hex(
                mapping = aes(x = x, y = y, fill = group_by, color = group_by),
                linewidth = hex_linewidth, bins = hex_bins, binwidth = hex_binwidth
            )
        }
    } else {
        p <- p + geom_point(
            mapping = aes(x = x, y = y, color = group_by, fill = group_by),
            size = pt_size, alpha = pt_alpha
        )
    }

    if (!is.null(cells.highlight_use) && !isTRUE(hex)) {
        cell_df <- subset(p$data, rownames(p$data) %in% cells.highlight_use)
        if (nrow(cell_df) > 0) {
            if (isTRUE(raster)) {
                p <- p + scattermore::geom_scattermore(
                    data = cell_df, aes(x = x, y = y), color = cols_highlight,
                    pointsize = floor(sizes_highlight) + stroke_highlight, alpha = alpha_highlight, pixels = raster_dpi
                ) +
                    scattermore::geom_scattermore(
                        data = cell_df, aes(x = x, y = y, color = group_by),
                        pointsize = floor(sizes_highlight), alpha = alpha_highlight, pixels = raster_dpi
                    )
            } else {
                p <- p + geom_point(
                    data = cell_df, aes(x = x, y = y), color = cols_highlight,
                    size = sizes_highlight + stroke_highlight, alpha = alpha_highlight
                ) +
                    geom_point(
                        data = cell_df, aes(x = x, y = y, color = group_by),
                        size = sizes_highlight, alpha = alpha_highlight
                    )
            }
        }
    }
    p <- p + scale_color_manual(
        name = paste0(group_by, ":"),
        values = colors[names(labels_tb)],
        labels = label_use,
        na.value = bg_color,
        guide = guide_legend(
            title.hjust = 0,
            order = 1,
            override.aes = list(size = 4, alpha = 1)
        )
    ) + scale_fill_manual(
        name = paste0(group_by, ":"),
        values = colors[names(labels_tb)],
        labels = label_use,
        na.value = bg_color,
        guide = guide_legend(
            title.hjust = 0,
            order = 1
        )
    )
    p_base <- p

    if (!is.null(stat_by)) {
        coor_df <- aggregate(p$data[, c("x", "y")], by = list(p$data$group_by), FUN = median)
        colnames(coor_df)[1] <- "group"
        coor_df <- coor_df[!is.na(coor_df$group) & coor_df$group != "NA", , drop = FALSE]
        x_range <- diff(layer_scales(p)$x$range$range)
        y_range <- diff(layer_scales(p)$y$range$range)
        stat_plot <- subplots[levels(dat$group_by)]

        stat_plot_list <- list()
        for (i in seq_len(nrow(coor_df))) {
            stat_plot_list[[i]] <- annotation_custom(
                as_grob(stat_plot[[coor_df[i, "group"]]] + theme_void() + ggplot2::theme(legend.position = "none")),
                xmin = coor_df[i, "x"] - x_range * stat_plot_size / 2, ymin = coor_df[i, "y"] - y_range * stat_plot_size / 2,
                xmax = coor_df[i, "x"] + x_range * stat_plot_size / 2, ymax = coor_df[i, "y"] + y_range * stat_plot_size / 2
            )
        }
        p <- p + stat_plot_list
        leg <- get_plot_component(
            stat_plot[[coor_df[i, "group"]]] + ggplot2::theme(legend.position = "bottom"),
            "guide-box-bottom"
        )
        if (!inherits(leg, "zeroGrob")) {
            legend_list$stat_by <- get_plot_component(
                stat_plot[[coor_df[i, "group"]]] + ggplot2::theme(legend.position = "bottom"),
                "guide-box-bottom"
            )
        }
    }
    if (!is.null(lineages)) {
        lineages_layers <- c(list(new_scale_color()), lineages_layers)
        suppressMessages({
            legend_list[["lineages"]] <- get_plot_component(
                ggplot() + lineages_layers +
                theme_scp(
                    legend.position = "bottom",
                    legend.direction = legend.direction
                ),
                "guide-box-bottom"
            )
        })
        p <- suppressWarnings({
            p + lineages_layers + ggplot2::theme(legend.position = "none")
        })
        if (is.null(legend_list[["lineages"]])) {
            legend_list["lineages"] <- list(NULL)
        }
    }
    # if (!is.null(paga)) {
    #     paga_layers <- c(list(new_scale_color()), paga_layers)
    #     if (group_by != paga$groups) {
    #         suppressMessages({
    #             legend_list[["paga"]] <- get_legend(ggplot() +
    #                 paga_layers +
    #                 theme_scp(
    #                     legend.position = "bottom",
    #                     legend.direction = legend.direction
    #                 ))
    #         })
    #     }
    #     p <- suppressWarnings({
    #         p + paga_layers + ggplot2::theme(legend.position = "none")
    #     })
    #     if (is.null(legend_list[["paga"]])) {
    #         legend_list["paga"] <- list(NULL)
    #     }
    # }
    # if (!is.null(velocity)) {
    #     velocity_layers <- c(list(new_scale_color()), list(new_scale("size")), velocity_layers)
    #     if (velocity_plot_type != "raw") {
    #         suppressMessages({
    #             legend_list[["velocity"]] <- get_legend(ggplot() +
    #                 velocity_layers +
    #                 theme_scp(
    #                     legend.position = "bottom",
    #                     legend.direction = legend.direction
    #                 ))
    #         })
    #     }
    #     p <- suppressWarnings({
    #         p + velocity_layers + ggplot2::theme(legend.position = "none")
    #     })
    #     if (is.null(legend_list[["velocity"]])) {
    #         legend_list["velocity"] <- list(NULL)
    #     }
    # }
    if (isTRUE(label)) {
        label_df <- aggregate(p$data[, c("x", "y")], by = list(p$data[["group_by"]]), FUN = median)
        colnames(label_df)[1] <- "label"
        label_df <- label_df[!is.na(label_df[, "label"]), , drop = FALSE]
        if (!isTRUE(label_insitu)) {
            label_df[, "label"] <- seq_len(nrow(label_df))
        }
        if (isTRUE(label_repel)) {
            p <- p + geom_point(
                data = label_df, mapping = aes(x = x, y = y),
                color = label_point_color, size = label_point_size
            ) + geom_text_repel(
                data = label_df, aes(x = x, y = y, label = label),
                fontface = "bold", min.segment.length = 0, segment.color = label_segment_color,
                point.size = label_point_size, max.overlaps = 100, force = label_repulsion,
                color = label_fg, bg.color = label_bg, bg.r = label_bg_r, size = label_size, inherit.aes = FALSE
            )
        } else {
            p <- p + geom_text_repel(
                data = label_df, aes(x = x, y = y, label = label),
                fontface = "bold", min.segment.length = 0, segment.color = label_segment_color,
                point.size = NA, max.overlaps = 100, force = 0,
                color = label_fg, bg.color = label_bg, bg.r = label_bg_r, size = label_size, inherit.aes = FALSE
            )
        }
    }
    if (length(legend_list) > 0) {
        legend_list <- legend_list[!sapply(legend_list, is.null)]
        legend_base <- get_plot_component(
            p_base + theme_scp(
                legend.position = "bottom",
                legend.direction = legend.direction
            ),
            "guide-box-bottom"
        )

        if (legend.direction == "vertical") {
            legend <- do.call(cbind, c(list(base = legend_base), legend_list))
        } else {
            legend <- do.call(rbind, c(list(base = legend_base), legend_list))
        }
        gtable <- as_grob(p + ggplot2::theme(legend.position = "none"))
        gtable <- add_grob(gtable, legend, legend.position)
        p <- wrap_plots(gtable)
    }

    facet_plot(p, facet, facet_scales, split_by, nrow, ncol, byrow)
}
