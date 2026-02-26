#' Arguments for spatial plot functions
#' @keywords internal
#' @name spatialplot_args
#' @param object A Seurat object or a Giotto object.
#' @param fov The name of the field of view (FOV) to plot, only works for Seurat objects.
#' @param boundaries The name of the boundaries to plot, only works for Seurat objects.
#' @param image The name of the image to plot. Possible values are:
#' * NULL: For Seurat objects with Visium data, the first image will be used.
#' For Giotto objects and Seurat objects with other spatial data, no image will be plotted.
#' * image name(s): the name of the image(s) to plot.
#' * color name: a color to use as a background for the plot.
#' * TRUE: For Seurat objects with Visium data, the first image will be used.
#' For Seurat objects with other spatial data, an error will be raised.
#' For Giotto objects with FOV, all (non-overlapping) images will be plotted; otherwise first image will be used.
#' * FALSE: disable image plotting.
#' @param masks Logical, whether to plot masks. Not supported for Seurat or Giotto objects for now.
#' @param shapes Plot shapes on the spatial plot.
#' Not supported for Seurat objects with Visium or SlideSeq data.
#' For Seurat objects with FOV, when TRUE, `boundaries` will be
#' used as boundaries for the shapes. Otherwise itself will be used as boundaries.
#' For Giotto objects or Seurat objects with FOV, this defaults to TRUE when `shapes_fill_by` is provided.
#' Set to FALSE to disable shapes plotting.
#' @param points Logical, whether to plot points. If TRUE, the points will be plotted using the coordinates from the object.
#' Defaults to TRUE.
#' @param ext The extent of the plot. If NULL, the extent will be calculated from the data.
#' If a numeric vector of length 4, it should be in the format c(xmin, xmax, ymin, ymax).
#' It can also be an object created by [terra::ext()].
#' @param crop Whether to crop the plot to the extent of the available data. Similar to `crop` argument in
#' [Seurat::SpatialDimPlot()]. Defaults to TRUE.
#' @param group_by The name of the metadata column to group the points by. Should be a character or factor column.
#' A special value "molecules" can be used to plot molecules in the FOV.
#' @param features A character vector of feature names to plot. If provided, the points will be colored by the features.
#' For `SpatDimPlot()`, this will be used to plot the molecules in the FOV.
#' For `SpatFeaturePlot()`, the plots will be faceted by the features.
#' @param layer The layer to use for the feature expression data.
#' Applicable for both Seurat and Giotto objects.
#' Defaults to "data" for Seurat objects, and "normalized" for Giotto objects.
#' For Giotto objects, it can also be "scaled", "raw", "counts", or "custom".
#' For Seurat objects, it can be "data", "scale.data", or "counts".
#' @param scale_factor Internal use only. The scale factor to use for the image, which will be extracted from the object.
#' @param layers A character vector of layers to plot. Possible values are:
#' * "image": plot the image as a background, which should be the first layer if provided.
#' * "masks": plot the masks
#' * "shapes": plot the shapes
#' * "points": plot the points
#' The order of the layers matters, as the first layer will be plotted first.
#' You can also use it to disable some layers by excluding them from the vector.
#' @param flip_y Internal use mostly, unless you want to flip the y-axis of the plot.
#' @param padding The padding to add to the extent of the plot, only available when `crop = TRUE` and `ext = NULL`.
#' For Seurat objects with FOV, this defaults to 0. In other cases, this defaults to 0 when image is plotted, and 0.05 otherwise.
#' @param image_scale Choose the scale factor ("lowres"/"hires") to apply in order to matchthe plot with the specified `image`.
#' Similar to `image.scale` argument in [Seurat::SpatialDimPlot()].
#' @param x Internal use only, the name of the x coordinate column in the data. Used to adopt different data types.
#' @param y Internal use only, the name of the y coordinate column in the data. Used to adopt different data types.
#' @param nmols Max number of each molecule specified in `features` for dim plot
#' Similar to `nmols` argument in [Seurat::ImageDimPlot()]. It also applied to Giotto objects.
#' @param shapes_fill_by The name of the variable to fill the shapes by.
#' It can also be a color name, in which case the shapes will be filled with that color.
#' When this is provided, the `shapes` argument will be set to TRUE by default.
#' @param graph The name of the graph to use for the spatial plot. Currently only supported for Giotto objects.
#' The graph data is obtained using [GiottoClass::getSpatialNetwork()].
#' When TRUE, the default graph will be used.
#' When given as a character, it should be the name of the graph to use.
#' If there is ":" in the name, the first part will be used as spat_unit, and the second part as the graph name.
#' @param shape The shape of the points, alias of `points_shape`.
#' See \url{https://ggplot2.tidyverse.org/reference/aes_linetype_size_shape.html} for more details.
#' @param legend.position The position of the legend. Defaults to "right".
#' @param legend.direction The direction of the legend. Defaults to "vertical".
#' @param theme The theme to use for the plot. Defaults to `"theme_box"`.
#' It can be the name of the theme (e.g. "ggplot2::theme_bw") or the function itself.
#' There are three themes that can be passed without namespace: "theme_box", "theme_this" and "theme_blank", which
#' are actually aliases of [plotthis::theme_box()], [plotthis::theme_this()] and [ggplot2::theme_void()] (without braces).
#' @param theme_args A list of arguments to pass to the theme function.
#' @param title The title of the plot. If NULL, no title will be added.
#' @param subtitle The subtitle of the plot. If NULL, no subtitle will be added.
#' @param xlab The label for the x-axis. If NULL, no label will be added.
#' @param ylab The label for the y-axis. If NULL, no label will be added.
#' @param facet_scales The scales to use for the facets. Defaults to "free". Can be "free", "fixed", "free_x", "free_y".
#' @param facet_nrow The number of rows to use for the facets. Defaults to NULL, which means the number of rows will be calculated automatically.
#' @param facet_ncol The number of columns to use for the facets. Defaults to NULL, which means the number of columns will be calculated automatically.
#' @param facet_byrow Logical, whether to facet by row. Defaults to FALSE.
#' @param feat_type feature type of the features (e.g. "rna", "dna", "protein"), only applied to Giotto objects.
#' @param use_overlap use polygon and feature coordinates overlap results, only applied to Giotto objects.
#' @param shapes_feat_type feature type of the features to use for shapes (e.g. "rna", "dna", "protein"), only applied to Giotto objects.
#' @param shapes_alpha The alpha value to use for the shapes.
#' When "points" are plotted, this defaults to 0.5; otherwise it defaults to 1.
#' @param spat_unit The spatial unit to use for the plot. Only applied to Giotto objects.
#' @param spat_loc_name The name of the spatial location to use for the plot. Only applied to Giotto objects.
#' @param spat_enr_names The names of the spatial enrichment results to use for the plot. Only applied to Giotto objects.
#' @param ... Additional arguments that will be passed to the spatial plot function.
#' With the `image_` prefix, these arguments will be used to plot the image ([plotthis::SpatImagePlot()]).
#' With the `masks_` prefix, these arguments will be used to plot the masks ([plotthis::SpatMasksPlot()]).
#' With the `shapes_` prefix, these arguments will be used to plot the shapes ([plotthis::SpatShapesPlot()]).
#' With the `points_` prefix, these arguments will be used to plot the points ([plotthis::SpatPointsPlot()]).
#' If the prefix is not provided, the arguments will be used as points arguments, but
#' with lower priority than the `points_` prefixed arguments.
NULL

#' Draw a rectangle as the image (background) of a spatial plot and fill with a given color
#' @keywords internal
#' @param args A list of arguments (the dot arguments) passed to the spatial plot function.
#' @param fill The color to fill the rectangle.
#' @return A ggplot2 layer object that can be added to a ggplot.
.rect_bg_image <- function(args, fill) {
    # use as a color to plot a background rect
    image_args <- args[startsWith(names(args), "image_")]
    names(image_args) <- sub("^image_", "", names(image_args))
    image_args$data = data.frame(x = 0, y = 0)  # dummy data
    image_args$xmin = -Inf
    image_args$xmax = Inf
    image_args$ymin = -Inf
    image_args$ymax = Inf
    image_args$fill <- fill
    image_args$color <- NA
    do_call(ggplot2::geom_rect, image_args)
}

#' Get points arguments from the dot-dot-dot arguments in the spatial plot function
#' @keywords internal
#' @param args A list of arguments (the dot arguments) passed to the spatial plot function.
#' @param ... Additional arguments that will be parsed.
#' @return A list of points arguments.
.points_args <- function(args, ...) {
    args <- c(args, rlang::dots_list(..., .named = TRUE))
    if (length(args) == 0) {
        return(list())
    }

    # The arguments without prefices are assumed to be the points arguments
    points_args <- args[
        !startsWith(names(args), "image_") & !startsWith(names(args), "masks_") &
        !startsWith(names(args), "shapes_") & !startsWith(names(args), "points_")
    ]
    # The arguments with points_ prefix are also considered as points arguments
    # but with higher priority
    points_args2 <- args[startsWith(names(args), "points_")]
    names(points_args2) <- sub("^points_", "", names(points_args2))
    conflicts <- intersect(names(points_args), names(points_args2))
    if (length(conflicts) > 0) {
        warning(paste0(
            "[SpatPlot] The following arguments are defined as points_* arguments: ",
            paste(conflicts, collapse = ", "),
            ", and also passed directly (without 'points_' prefix). The points_* arguments will be used."
        ))
    }
    out <- c(points_args2, points_args)
    dup_names <- duplicated(names(out))
    if (any(dup_names)) {
        out[-which(dup_names)]
    } else {
        out
    }
}

#' Process points layer for Seurat spatial plots
#' @keywords internal
#' @inheritParams spatialplot_args
#' @param swap_xy Logical, whether to swap x and y coordinates in the points data.
#' Seurat objects, loaded with Visium, Xenium or SlideSeq data, have x and y coordinates
#' in the opposite order. But when loaded with FOV data, the coordinates are in the correct order.
#' @return A list containing the ggplot2 layer object and the facet_by variable if applicable.
.seurat_points_layer <- function(
    object, fov = NULL, boundaries = NULL, x = "x", y = "y", swap_xy = TRUE,
    image, args, crop, points_data, ext_unscaled, scale_factor, group_by, shape,
    features, layer, legend.position, legend.direction, flip_y, ext
) {
    # The arguments passed as points_* are collected in args
    if (shape != 16) {
        points_args <- .points_args(args, x = x, y = y, shape = shape)
    } else {
        # allow shape to be overridden by points_shape without warning
        points_args <- .points_args(args, x = x, y = y)
    }
    get_cells <- if (is.null(fov)) rownames else function(x) x$cell

    if (crop) {
        # attach metadata for highlighting selection
        points_args$data <- object@meta.data[get_cells(points_data), , drop = FALSE]
        # Handle missing columns for VisiumV2 vs SlideSeq difference
        if (ncol(points_args$data) > 0) {
            points_args$data <- points_args$data[, colnames(points_args$data)[!is.na(colnames(points_args$data))], drop = FALSE]
        }
        points_args$data <- cbind(points_data, points_args$data)
        # check if there is duplicated columns
        dup_cols <- duplicated(colnames(points_args$data))
        if (any(dup_cols)) {
            points_args$data <- points_args$data[, !dup_cols, drop = FALSE]
        }
        # x and y are already handled in points_data, we don't need to handle the swapping here
        points_args$data[[points_args$x]] <- points_args$data[[points_args$x]] * scale_factor
        points_args$data[[points_args$y]] <- points_args$data[[points_args$y]] * scale_factor
        points_args$ext <- ext %||% (ext_unscaled * scale_factor)
    } else {
        points_args$data <- if (is.null(fov)) {
            Seurat::GetTissueCoordinates(object, image = if(isFALSE(image)) NULL else image)
        } else {
            ggplot2::fortify(object[[fov]][[boundaries]])
        }
        # remove the NA columns (NA as column names, which will error later)
        points_args$data <- points_args$data[, colnames(points_args$data)[!is.na(colnames(points_args$data))], drop = FALSE]
        if (swap_xy) {
            points_args$data$.tmp <- points_args$data[[points_args$x]]
            points_args$data[[points_args$x]] <- points_args$data[[points_args$y]] * scale_factor
            points_args$data[[points_args$y]] <- points_args$data$.tmp * scale_factor
            points_args$data$.tmp <- NULL
        } else {
            points_args$data[[points_args$x]] <- points_args$data[[points_args$x]] * scale_factor
            points_args$data[[points_args$y]] <- points_args$data[[points_args$y]] * scale_factor
        }
        points_args$data <- cbind(
            object@meta.data[rownames(points_args$data), , drop = FALSE],
            points_args$data
        )
    }

    facet_by <- NULL
    if (!is.null(group_by)) {
        if (nrow(object@meta.data) == nrow(points_args$data)) {
            points_args$data[[group_by]] <- object@meta.data[[group_by]]
        } else {
            # data is subset by fov
            points_args$data[[group_by]] <- object@meta.data[get_cells(points_args$data), group_by, drop = TRUE]
        }
        points_args$color_by <- group_by
    } else if (!is.null(features)) {
        cells_by_image <- utils::getFromNamespace("CellsByImage", "Seurat")
        cells <- unique(cells_by_image(object, images = if(isFALSE(image)) NULL else image, unlist = TRUE))
        featdata <- Seurat::FetchData(
            object = object,
            vars = features,
            cells = cells,
            layer = layer,
            clean = FALSE
        )
        features <- colnames(featdata)
        points_args$data[, features] <- featdata
        points_args$color_by <- features
        if (length(features) == 1) {
            points_args$color_name <- points_args$color_name %||% features
        } else {
            points_args$color_name <- points_args$color_name %||% "feature"
            facet_by <- ".facet_var"
        }
    }
    points_args$legend.position <- legend.position
    points_args$legend.direction <- legend.direction
    points_args$flip_y <- flip_y
    points_args$return_layer <- TRUE
    player <- do_call(SpatPointsPlot, points_args)

    list(player = player, facet_by = facet_by)
}

#' Process points layer with molecules for Seurat spatial plots
#' @keywords internal
#' @inheritParams .seurat_points_layer
#' @return A list containing the ggplot2 layer object and the facet_by variable if applicable.
.seurat_points_layer_molecules <- function(
    object, fov, boundaries, x = "x", y = "y", swap_xy = TRUE,
    image, args, nmols, crop, points_data, ext_unscaled, scale_factor, group_by, shape,
    features, layer, legend.position, legend.direction, flip_y, ext
) {

    points_args <- .points_args(args, x = x, y = y, shape = shape)

    points_args$data <- Seurat::FetchData(object[[fov]], vars = features, nmols = nmols)
    if (swap_xy) {
        points_args$data$.tmp <- points_args$data[[points_args$x]]
        points_args$data[[points_args$x]] <- points_args$data[[points_args$y]] * scale_factor
        points_args$data[[points_args$y]] <- points_args$data$.tmp * scale_factor
        points_args$data$.tmp <- NULL
    } else {
        points_args$data[[points_args$x]] <- points_args$data[[points_args$x]] * scale_factor
        points_args$data[[points_args$y]] <- points_args$data[[points_args$y]] * scale_factor
    }
    points_args$data$molecules <- points_args$data$molecule
    points_args$data$molecule <- NULL

    if (crop) {
        points_args$ext <- ext %||% (ext_unscaled * scale_factor)
    }

    facet_by <- NULL

    points_args$color_by <- group_by  # molecules
    points_args$color_name <- points_args$color_name %||% "Molecules"

    points_args$legend.position <- legend.position
    points_args$legend.direction <- legend.direction
    points_args$flip_y <- flip_y
    points_args$return_layer <- TRUE
    player <- do_call(SpatPointsPlot, points_args)

    list(player = player, facet_by = facet_by)
}

#' Add new scale layers
#' @keywords internal
#' @param scales The scales have been already applied in the previous layers
#' (so that we need to add new scales only for the new layers).
#' @return A list of ggnewscale layers for the specified scales.
.new_scale_layers <- function(scales) {
    layers <- list()
    for (scale in scales) {
        if (scale == "fill") {
            layers <- c(layers, list(ggnewscale::new_scale_fill()))
        } else if (scale == "color") {
            layers <- c(layers, list(ggnewscale::new_scale_color()))
        } else {
            stop(paste("Unknown scale:", scale))
        }
    }
    layers
}
