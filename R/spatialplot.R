
#' Plot features for spatial data
#'
#' The features can include  expression, dimension reduction components, metadata, etc
#'
#' @inheritParams spatialplot_args
#' @return A ggplot object
#' @keywords internal
SpatPlot <- function(
    object, fov = NULL, boundaries = NULL, image = NULL, masks = NULL, shapes = NULL, points = NULL,
    ext = NULL, crop = TRUE, group_by = NULL, features = NULL, layer = NULL, scale_factor = NULL,
    layers = NULL, flip_y = NULL, padding = NULL, image_scale = NULL,
    x = NULL, y = NULL, nmols = 1000, shapes_fill_by = NULL, graph = NULL, shape = 16,
    legend.position = "right", legend.direction = "vertical", theme = "theme_box", theme_args = list(),
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
    facet_scales = "fixed", facet_nrow = NULL, facet_ncol = NULL, facet_byrow = TRUE,
    feat_type = "rna", use_overlap = FALSE, shapes_feat_type = "cell", shapes_alpha = NULL,
    spat_unit = NULL, spat_loc_name = NULL, spat_enr_names = NULL,
    ...
) {
    UseMethod("SpatPlot", object)
}

#' @keywords internal
SpatPlot.Seurat <- function(
    object, fov = NULL, boundaries = NULL, image = NULL, masks = NULL, shapes = NULL, points = NULL,
    ext = NULL, crop = TRUE, group_by = NULL, features = NULL, layer = NULL, scale_factor = NULL,
    layers = NULL, flip_y = NULL, padding = NULL, image_scale = NULL,
    x = NULL, y = NULL, nmols = 1000, shapes_fill_by = NULL, graph = NULL, shape = 16,
    legend.position = "right", legend.direction = "vertical", theme = "theme_box", theme_args = list(),
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
    facet_scales = "fixed", facet_nrow = NULL, facet_ncol = NULL, facet_byrow = TRUE,
    feat_type = "rna", use_overlap = FALSE, shapes_feat_type = "cell", shapes_alpha = NULL,
    spat_unit = NULL, spat_loc_name = NULL, spat_enr_names = NULL,
    ...
) {
    first_image <- Seurat::Images(object)[1]
    if (is.null(first_image)) {
        stop("[SpatPlot] No images found in the Seurat object. Is this an object with spatial data? ")
    }
    stype <- class(object@images[[first_image]])

    if ("VisiumV1" %in% stype) {
        if (!is.null(fov) || !is.null(boundaries)) {
            stop("[SpatPlot] 'fov' and 'boundaries' are not supported for Seurat objects with Visium data.")
        }
        if (isTRUE(image) || is.null(image)) {
            image <- first_image
        }
        SpatPlot.Seurat.Visium(
            object, fov = fov, boundaries = boundaries, image = image, masks = masks, shapes = shapes, points = points,
            ext = ext, crop = crop, group_by = group_by, features = features, layer = layer, scale_factor = 1,
            layers = layers, flip_y = flip_y, padding = padding, image_scale = image_scale,
            x = "imagerow", y = "imagecol", nmols = nmols, shapes_fill_by = shapes_fill_by, graph = graph, shape = shape,
            legend.position = legend.position, legend.direction = legend.direction, theme = theme, theme_args = theme_args,
            title = title, subtitle = subtitle, xlab = xlab, ylab = ylab,
            facet_scales = facet_scales, facet_nrow = facet_nrow, facet_ncol = facet_ncol, facet_byrow = facet_byrow,
            feat_type = feat_type, use_overlap = use_overlap, shapes_feat_type = shapes_feat_type, shapes_alpha = shapes_alpha,
            spat_unit = spat_unit, spat_loc_name = spat_loc_name, spat_enr_names = spat_enr_names,
            ...
        )
    } else if ("VisiumV2" %in% stype) {
        if (!is.null(fov) || !is.null(boundaries)) {
            stop("[SpatPlot] 'fov' and 'boundaries' are not supported for Seurat objects with Visium data.")
        }
        if (isTRUE(image) || is.null(image)) {
            image <- first_image
        }
        SpatPlot.Seurat.Visium(
            object, fov = fov, boundaries = boundaries, image = image, masks = masks, shapes = shapes, points = points,
            ext = ext, crop = crop, group_by = group_by, features = features, layer = layer, scale_factor = scale_factor,
            layers = layers, flip_y = flip_y, padding = padding %||% 0.05, image_scale = image_scale,
            x = x, y = y, nmols = nmols, shapes_fill_by = shapes_fill_by, graph = graph, shape = shape,
            legend.position = legend.position, legend.direction = legend.direction, theme = theme, theme_args = theme_args,
            title = title, subtitle = subtitle, xlab = xlab, ylab = ylab,
            facet_scales = facet_scales, facet_nrow = facet_nrow, facet_ncol = facet_ncol, facet_byrow = facet_byrow,
            feat_type = feat_type, use_overlap = use_overlap, shapes_feat_type = shapes_feat_type, shapes_alpha = shapes_alpha,
            spat_unit = spat_unit, spat_loc_name = spat_loc_name, spat_enr_names = spat_enr_names,
            ...
        )
    } else if ("FOV" %in% stype) {
        assay <- object@images[[first_image]]@assay
        if (identical(assay, "Nanostring")) {
            x <- x %||% "x"
            y <- y %||% "y"
        } else {
            x <- x %||% "y"
            y <- y %||% "x"
        }

        SpatPlot.Seurat.FOV(
            object, fov = fov, boundaries = boundaries, image = image, masks = masks, shapes = shapes, points = points,
            ext = ext, crop = crop, group_by = group_by, features = features, layer = layer, scale_factor = scale_factor,
            layers = layers, flip_y = flip_y, padding = padding, image_scale = image_scale,
            x = x, y = y, nmols = nmols, shapes_fill_by = shapes_fill_by, graph = graph, shape = shape,
            legend.position = legend.position, legend.direction = legend.direction, theme = theme, theme_args = theme_args,
            title = title, subtitle = subtitle, xlab = xlab, ylab = ylab,
            facet_scales = facet_scales, facet_nrow = facet_nrow, facet_ncol = facet_ncol, facet_byrow = facet_byrow,
            feat_type = feat_type, use_overlap = use_overlap, shapes_feat_type = shapes_feat_type, shapes_alpha = shapes_alpha,
            spat_unit = spat_unit, spat_loc_name = spat_loc_name, spat_enr_names = spat_enr_names,
            ...
        )
    } else if ("SlideSeq" %in% stype) {
        if (!is.null(fov) || !is.null(boundaries)) {
            stop("[SpatPlot] 'fov' and 'boundaries' are not supported for Seurat objects with SlideSeq data.")
        }
        SpatPlot.Seurat.SlideSeq(
            object, fov = fov, boundaries = boundaries, image = image, masks = masks, shapes = shapes, points = points,
            ext = ext, crop = crop, group_by = group_by, features = features, layer = layer, scale_factor = scale_factor,
            layers = layers, flip_y = flip_y, padding = padding, image_scale = image_scale,
            x = x, y = y, nmols = nmols, shapes_fill_by = shapes_fill_by, graph = graph, shape = shape,
            legend.position = legend.position, legend.direction = legend.direction, theme = theme, theme_args = theme_args,
            title = title, subtitle = subtitle, xlab = xlab, ylab = ylab,
            facet_scales = facet_scales, facet_nrow = facet_nrow, facet_ncol = facet_ncol, facet_byrow = facet_byrow,
            feat_type = feat_type, use_overlap = use_overlap, shapes_feat_type = shapes_feat_type, shapes_alpha = shapes_alpha,
            spat_unit = spat_unit, spat_loc_name = spat_loc_name, spat_enr_names = spat_enr_names,
            ...
        )
    } else {
        stop("[SpatPlot] Unsupported image type. Supported types are 'VisiumV1', 'VisiumV2', 'FOV' and 'SlideSeq'.")
    }
}

#' @keywords internal
#' @importFrom plotthis SpatImagePlot SpatMasksPlot SpatShapesPlot SpatPointsPlot
SpatPlot.Seurat.Visium <- function(
    object, fov = NULL, boundaries = NULL, image = NULL, masks = NULL, shapes = NULL, points = NULL,
    ext = NULL, crop = TRUE, group_by = NULL, features = NULL, layer = NULL, scale_factor = NULL,
    layers = NULL, flip_y = NULL, padding = NULL, image_scale = NULL,
    x = NULL, y = NULL, nmols = 1000, shapes_fill_by = NULL, graph = NULL, shape = 16,
    legend.position = "right", legend.direction = "vertical", theme = "theme_box", theme_args = list(),
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
    facet_scales = "fixed", facet_nrow = NULL, facet_ncol = NULL, facet_byrow = TRUE,
    feat_type = "rna", use_overlap = FALSE, shapes_feat_type = "cell", shapes_alpha = NULL,
    spat_unit = NULL, spat_loc_name = NULL, spat_enr_names = NULL,
    ...
) {
    # because we don't have molecules in Visium data so they should be both not provided
    stopifnot("[SpatPlot] Either 'group_by' or 'features' should be provided, not both." =
        is.null(group_by) || is.null(features))
    if (is.character(theme) && theme %in% c("theme_box", "theme_this", "theme_blank")) {
        theme <- utils::getFromNamespace(theme, "plotthis")
    }

    x <- x %||% "x"
    y <- y %||% "y"
    flip_y <- flip_y %||% TRUE
    layer <- layer %||% "data"
    points <- points %||% TRUE

    layers <- intersect(
        layers %||% c("image", "masks", "shapes", "points"),
        c(
            if (!is.null(image) && !isFALSE(image)) "image",
            if (!is.null(masks) && !isFALSE(masks)) "masks",
            if (!is.null(shapes) && !isFALSE(shapes)) "shapes",
            if (!is.null(points) && !isFALSE(points)) "points"
        )
    )

    stopifnot('Either "image", "masks", "shapes", or "points" must be provided.' =
        any(layers %in% c("image", "masks", "shapes", "points")))
    if ("image" %in% layers && layers[1] != "image") {
        stop('If "image" is provided, it must be the first element in "layers".')
    }

    players <- list()
    scales_used <- c()
    args <- rlang::dots_list(...)
    the_scale_factor <- 1
    facet_by <- NULL
    ext_unscaled <- NULL
    if (crop) {
        points_data <- Seurat::GetTissueCoordinates(object)
        points_data <- points_data[, colnames(points_data)[!is.na(colnames(points_data))], drop = FALSE]
        points_data$.tmp <- points_data[[x]]
        points_data[[x]] <- points_data[[y]]
        points_data[[y]] <- points_data$.tmp
        points_data$.tmp <- NULL
        if (is.null(padding)) {
            padding <- if ("image" %in% layers) 0 else 0.05
        }
        delta_x <- diff(range(points_data[[x]], na.rm = TRUE)) * padding
        delta_y <- diff(range(points_data[[y]], na.rm = TRUE)) * padding
        ext_unscaled <- c(
            min(points_data[[x]], na.rm = TRUE) - delta_x,
            max(points_data[[x]], na.rm = TRUE) + delta_x,
            min(points_data[[y]], na.rm = TRUE) - delta_y,
            max(points_data[[y]], na.rm = TRUE) + delta_y
        )
    }

    for (element in layers) {
        if (element == "image") {
            if (image %in% names(object@images)) {
                image_obj <- object@images[[image]]
                img <- terra::rast(image_obj@image)
                image_scale <- image_scale %||% which.min(image_obj@scale.factors)
                the_scale_factor <- scale_factor %||% image_obj@scale.factors[[image_scale]]
                # img <- terra::crop(img, ext, extend = TRUE)
                image_args <- args[startsWith(names(args), "image_")]
                names(image_args) <- sub("^image_", "", names(image_args))
                image_args$data <- terra::flip(img, direction = "vertical")
                image_args$flip_y <- flip_y
                image_args$return_layer <- TRUE
                image_args$ext <- if (!is.null(ext_unscaled)) ext_unscaled * the_scale_factor
                player <- do.call(SpatImagePlot, image_args)
                scales_reused <- intersect(scales_used, attr(player, "scales"))
                players <- c(players, list(player))
                scales_used <- unique(c(scales_used, attr(player, "scales")))
            } else {
                # use as a color to plot a background rect
                player <- .rect_bg_image(args, image)
                scales_reused <- intersect(scales_used, "fill")
                players <- c(players, .new_scale_layers(scales_reused), list(player))
                scales_used <- unique(c(scales_used, attr(player, "scales")))
            }
        }
        if (element == "points") {
            points_layer <- .seurat_points_layer(
                object = object, image = image, args = args, crop = crop, points_data = points_data,
                x = x, y = y, shape = shape,
                ext_unscaled = ext_unscaled, scale_factor = the_scale_factor, group_by = group_by,
                features = features, layer = layer, legend.position = legend.position,
                legend.direction = legend.direction, flip_y = flip_y, ext = ext
            )
            player <- points_layer$player
            facet_by <- points_layer$facet_by
            scales_reused <- intersect(scales_used, attr(player, "scales"))
            players <- c(players, .new_scale_layers(scales_reused), list(player))
            scales_used <- unique(c(scales_used, attr(player, "scales")))
        }
        if (element == "masks") {
            stop("[SpatPlot] 'masks' is not supported for Seurat Visium objects.")
        }
        if (element == "shapes") {
            stop("[SpatPlot] 'shapes' is not supported for Seurat Visium objects.")
        }
    }

    if (!is.null(ext_unscaled)) {
        ext <- ext %||% (ext_unscaled * the_scale_factor)
    }
    p <- utils::getFromNamespace(".wrap_spatial_layers", "plotthis")(
        layers = players,
        ext = ext,
        flip_y = flip_y,
        legend.position = legend.position,
        legend.direction = legend.direction,
        title = title,
        subtitle = subtitle,
        xlab = xlab,
        ylab = ylab,
        theme = theme,
        theme_args = theme_args
    )

    if (!is.null(facet_by)) {
        p <- utils::getFromNamespace("facet_plot", "plotthis")(
            p, facet_by, facet_scales, facet_nrow, facet_ncol, facet_byrow,
            legend.position = legend.position, legend.direction = legend.direction
        )
    }

    p
}

#' @keywords internal
SpatPlot.Seurat.SlideSeq <- function(
    object, fov = NULL, boundaries = NULL, image = NULL, masks = NULL, shapes = NULL, points = NULL,
    ext = NULL, crop = TRUE, group_by = NULL, features = NULL, layer = NULL, scale_factor = NULL,
    layers = NULL, flip_y = NULL, padding = NULL, image_scale = NULL,
    x = NULL, y = NULL, nmols = 1000, shapes_fill_by = NULL, graph = NULL, shape = 16,
    legend.position = "right", legend.direction = "vertical", theme = "theme_box", theme_args = list(),
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
    facet_scales = "fixed", facet_nrow = NULL, facet_ncol = NULL, facet_byrow = TRUE,
    feat_type = "rna", use_overlap = FALSE, shapes_feat_type = "cell", shapes_alpha = NULL,
    spat_unit = NULL, spat_loc_name = NULL, spat_enr_names = NULL,
    ...
) {
    stopifnot("[SpatPlot] Either 'group_by' or 'features' should be provided, not both." =
        is.null(group_by) || is.null(features))
    if (is.character(theme) && theme %in% c("theme_box", "theme_this", "theme_blank")) {
        theme <- utils::getFromNamespace(theme, "plotthis")
    }

    flip_y <- flip_y %||% TRUE
    layer <- layer %||% "data"
    points <- points %||% TRUE

    layers <- intersect(
        layers %||% c("image", "masks", "shapes", "points"),
        c(
            if (!is.null(image) && !isFALSE(image)) "image",
            if (!is.null(masks) && !isFALSE(masks)) "masks",
            if (!is.null(shapes) && !isFALSE(shapes)) "shapes",
            if (!is.null(points) && !isFALSE(points)) "points"
        )
    )
    stopifnot('Either "image", "masks", "shapes", or "points" must be provided.' =
        any(layers %in% c("image", "masks", "shapes", "points")))

    players <- list()
    scales_used <- c()
    args <- rlang::dots_list(...)
    scale_factor <- 1
    facet_by <- NULL
    ext_unscaled <- NULL
    if (crop) {
        points_data <- Seurat::GetTissueCoordinates(object)
        points_data <- points_data[, colnames(points_data)[!is.na(colnames(points_data))], drop = FALSE]
        points_data$.y <- points_data$x
        points_data$x <- points_data$y
        points_data$y <- points_data$.y
        points_data$.y <- NULL
        if (is.null(padding)) {
            padding <- if ("image" %in% layers) 0 else 0.05
        }
        delta_x <- diff(range(points_data$x, na.rm = TRUE)) * padding
        delta_y <- diff(range(points_data$y, na.rm = TRUE)) * padding
        ext_unscaled <- c(
            min(points_data$x, na.rm = TRUE) - delta_x,
            max(points_data$x, na.rm = TRUE) + delta_x,
            min(points_data$y, na.rm = TRUE) - delta_y,
            max(points_data$y, na.rm = TRUE) + delta_y
        )
    }
    for (element in layers) {
        if (element == "image") {
            # use as a color to plot a background rect
            player <- .rect_bg_image(args, image)
            scales_reused <- intersect(scales_used, "fill")
            players <- c(players, .new_scale_layers(scales_reused), list(player))
            scales_used <- unique(c(scales_used, attr(player, "scales")))
        }
        if (element == "points") {
            points_layer <- .seurat_points_layer(
                object = object, image = image, args = args, crop = crop, points_data = points_data,
                ext_unscaled = ext_unscaled, scale_factor = scale_factor, group_by = group_by, shape = shape,
                features = features, layer = layer, legend.position = legend.position,
                legend.direction = legend.direction, flip_y = flip_y, ext = ext
            )
            player <- points_layer$player
            facet_by <- points_layer$facet_by
            scales_reused <- intersect(scales_used, attr(player, "scales"))
            players <- c(players, .new_scale_layers(scales_reused), list(player))
            scales_used <- unique(c(scales_used, attr(player, "scales")))
        }
        if (element == "masks") {
            stop("[SpatPlot] 'masks' is not supported for Seurat SlideSeq objects.")
        }
        if (element == "shapes") {
            stop("[SpatPlot] 'shapes' is not supported for Seurat SlideSeq objects.")
        }
    }

    if (!is.null(ext_unscaled)) {
        ext <- ext %||% (ext_unscaled * scale_factor)
    }

    p <- utils::getFromNamespace(".wrap_spatial_layers", "plotthis")(
        layers = players,
        ext = ext,
        flip_y = flip_y,
        legend.position = legend.position,
        legend.direction = legend.direction,
        title = title,
        subtitle = subtitle,
        xlab = xlab,
        ylab = ylab,
        theme = theme,
        theme_args = theme_args
    )

    if (!is.null(facet_by)) {
        p <- utils::getFromNamespace("facet_plot", "plotthis")(
            p, facet_by, facet_scales, facet_nrow, facet_ncol, facet_byrow,
            legend.position = legend.position, legend.direction = legend.direction
        )
    }

    p
}

#' @keywords internal
#' @importFrom plotthis SpatImagePlot SpatMasksPlot SpatShapesPlot SpatPointsPlot
SpatPlot.Seurat.FOV <- function(
    object, fov = NULL, boundaries = NULL, image = NULL, masks = NULL, shapes = NULL, points = NULL,
    ext = NULL, crop = TRUE, group_by = NULL, features = NULL, layer = NULL, scale_factor = NULL,
    layers = NULL, flip_y = NULL, padding = NULL, image_scale = NULL,
    x = NULL, y = NULL, nmols = 1000, shapes_fill_by = NULL, graph = NULL, shape = 16,
    legend.position = "right", legend.direction = "vertical", theme = "theme_box", theme_args = list(),
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
    facet_scales = "fixed", facet_nrow = NULL, facet_ncol = NULL, facet_byrow = TRUE,
    feat_type = "rna", use_overlap = FALSE, shapes_feat_type = "cell", shapes_alpha = NULL,
    spat_unit = NULL, spat_loc_name = NULL, spat_enr_names = NULL,
    ...
) {
    if (!identical(group_by, "molecules")) {
        stopifnot("[SpatPlot] Either 'group_by' or 'features' should be provided, not both." =
            is.null(group_by) || is.null(features))
    }
    if (is.character(theme) && theme %in% c("theme_box", "theme_this", "theme_blank")) {
        theme <- utils::getFromNamespace(theme, "plotthis")
    }

    x <- x %||% "y"
    y <- y %||% "x"
    flip_y <- flip_y %||% FALSE
    layer <- layer %||% "data"
    points <- points %||% TRUE
    shapes <- shapes %||% (if(!is.null(shapes_fill_by)) TRUE)

    layers <- intersect(
        layers %||% c("image", "masks", "shapes", "points"),
        c(
            if (!is.null(image) && !isFALSE(image)) "image",
            if (!is.null(masks) && !isFALSE(masks)) "masks",
            if (!is.null(shapes) && !isFALSE(shapes)) "shapes",
            if (!is.null(points) && !isFALSE(points)) "points"
        )
    )
    stopifnot('Either "image", "masks", "shapes", or "points" must be provided.' =
        any(layers %in% c("image", "masks", "shapes", "points")))
    if ("image" %in% layers && layers[1] != "image") {
        stop('If "image" is provided, it must be the first element in "layers".')
    }

    fov <- fov %||% SeuratObject::DefaultFOV(object)
    boundaries <- boundaries %||% SeuratObject::DefaultBoundary(object[[fov]])

    players <- list()
    scales_used <- c()
    args <- rlang::dots_list(...)
    the_scale_factor <- 1
    facet_by <- NULL
    ext_unscaled <- NULL
    if (crop) {
        # x	y	cell
        # <dbl>	<dbl>	<chr>
        # 6152.972	1477.041	caijolig-1
        # 6158.783	1455.548	caikcgdg-1
        points_data <- ggplot2::fortify(object[[fov]][[boundaries]])

        points_data <- points_data[, colnames(points_data)[!is.na(colnames(points_data))], drop = FALSE]
        # We don't need to swap x and y here
        # points_data$.tmp <- points_data[[x]]
        # points_data[[x]] <- points_data[[y]]
        # points_data[[y]] <- points_data$.tmp
        # points_data$.tmp <- NULL
        padding <- padding %||% 0
        delta_x <- diff(range(points_data[[x]], na.rm = TRUE)) * padding
        delta_y <- diff(range(points_data[[y]], na.rm = TRUE)) * padding
        ext_unscaled <- c(
            min(points_data[[x]], na.rm = TRUE) - delta_x,
            max(points_data[[x]], na.rm = TRUE) + delta_x,

            min(points_data[[y]], na.rm = TRUE) - delta_y,
            max(points_data[[y]], na.rm = TRUE) + delta_y
        )
    }
    for (element in layers) {
        if (element == "image") {
            # use as a color to plot a background rect
            player <- .rect_bg_image(args, image)
            scales_reused <- intersect(scales_used, "fill")
            players <- c(players, .new_scale_layers(scales_reused), list(player))
            scales_used <- unique(c(scales_used, attr(player, "scales")))
        }
        if (element == "points") {
            if (identical(group_by, "molecules")) {
                points_layer <- .seurat_points_layer_molecules(
                    object = object, image = image, args = args, crop = crop,
                    points_data = points_data, nmols = nmols, swap_xy = FALSE,
                    fov = fov, boundaries = boundaries, x = x, y = y, shape = shape,
                    ext_unscaled = ext_unscaled, scale_factor = the_scale_factor, group_by = group_by,
                    features = features, layer = layer, legend.position = legend.position,
                    legend.direction = legend.direction, flip_y = flip_y, ext = ext
                )
            } else {
                points_layer <- .seurat_points_layer(
                    object = object, image = if (is.character(image) && image %in% names(object@images)) image,
                    args = args, crop = crop, points_data = points_data, swap_xy = FALSE,
                    fov = fov, boundaries = boundaries, x = x, y = y, shape = shape,
                    ext_unscaled = ext_unscaled, scale_factor = the_scale_factor, group_by = group_by,
                    features = features, layer = layer, legend.position = legend.position,
                    legend.direction = legend.direction, flip_y = flip_y, ext = ext
                )
            }
            player <- points_layer$player
            facet_by <- points_layer$facet_by
            scales_reused <- intersect(scales_used, attr(player, "scales"))
            players <- c(players, .new_scale_layers(scales_reused), list(player))
            scales_used <- unique(c(scales_used, attr(player, "scales")))
        }
        if (element == "masks") {
            stop("[SpatPlot] 'masks' is not supported for Seurat FOV objects.")
        }
        if (element == "shapes") {
            if (isTRUE(shapes)) {
                if ("points" %in% layers) {
                    warning(
                        "[SpatPlot] 'shapes' is set to TRUE, meaning the same boundaries as points will be used. ",
                        "You may want to provide a different boundaries for shapes. Otherwise the shapes is plotted as points."
                    )
                }
                shapes <- boundaries
            }
            shapes_args <- args[startsWith(names(args), "shapes_")]
            names(shapes_args) <- sub("^shapes_", "", names(shapes_args))
            shapes_args$data <- ggplot2::fortify(object[[fov]][[shapes]])
            if (!is.null(shapes_fill_by)) {
                shapes_args$data[[shapes_fill_by]] <- object@meta.data[shapes_args$data$cell, shapes_fill_by, drop = TRUE]
            }
            shapes_args$flip_y <- flip_y
            shapes_args$return_layer <- TRUE
            shapes_args$ext <- ext
            shapes_args$x <- "x"
            shapes_args$y <- "y"
            shapes_args$group <- "cell"
            shapes_args$fill_by <- shapes_fill_by
            shapes_args$alpha <- ifelse("points" %in% layers, 0.5, 1)
            shapes_args$legend.direction <- legend.direction
            shapes_args$legend.position <- legend.position
            player <- do.call(SpatShapesPlot, shapes_args)
            scales_reused <- intersect(scales_used, attr(player, "scales"))
            players <- c(players, .new_scale_layers(scales_reused), list(player))
            scales_used <- unique(c(scales_used, attr(player, "scales")))
        }
    }

    if (!is.null(ext_unscaled)) {
        ext <- ext %||% (ext_unscaled * the_scale_factor)
    }

    p <- utils::getFromNamespace(".wrap_spatial_layers", "plotthis")(
        layers = players,
        ext = ext,
        flip_y = flip_y,
        legend.position = legend.position,
        legend.direction = legend.direction,
        title = title,
        subtitle = subtitle,
        xlab = xlab,
        ylab = ylab,
        theme = theme,
        theme_args = theme_args
    )

    if (!is.null(facet_by)) {
        p <- utils::getFromNamespace("facet_plot", "plotthis")(
            p, facet_by, facet_scales, facet_nrow, facet_ncol, facet_byrow,
            legend.position = legend.position, legend.direction = legend.direction
        )
    }

    p
}

#' @keywords internal
SpatPlot.giotto <- function(
    object, fov = NULL, boundaries = NULL, image = NULL, masks = NULL, shapes = NULL, points = NULL,
    ext = NULL, crop = TRUE, group_by = NULL, features = NULL, layer = NULL, scale_factor = NULL,
    layers = NULL, flip_y = NULL, padding = NULL, image_scale = NULL,
    x = NULL, y = NULL, nmols = 1000, shapes_fill_by = NULL, graph = NULL, shape = 16,
    legend.position = "right", legend.direction = "vertical", theme = "theme_box", theme_args = list(),
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
    facet_scales = "fixed", facet_nrow = NULL, facet_ncol = NULL, facet_byrow = TRUE,
    feat_type = "rna", use_overlap = FALSE, shapes_feat_type = "cell", shapes_alpha = NULL,
    spat_unit = NULL, spat_loc_name = NULL, spat_enr_names = NULL,
    ...
) {
    if (!identical(group_by, "molecules")) {
        stopifnot(
            "[SpatPlot] Either 'group_by' or 'features' should be provided, not both." =
            is.null(group_by) || is.null(features))
    }

    x <- x %||% "x"
    y <- y %||% "y"
    flip_y <- flip_y %||% FALSE
    layer <- layer %||% "normalized"

    if (is.character(theme) && theme %in% c("theme_box", "theme_this", "theme_blank")) {
        theme <- utils::getFromNamespace(theme, "plotthis")
    }

    spat_unit <- GiottoClass::set_default_spat_unit(
        gobject = object,
        spat_unit = spat_unit
    )

    feat_type <- GiottoClass::set_default_feat_type(
        gobject = object,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    if (is.null(spat_loc_name)) {
        if (!is.null(methods::slot(object, "spatial_locs"))) {
            spat_loc_name <- GiottoClass::list_spatial_locations_names(object, spat_unit = spat_unit)[[1]]
        } else {
            stop("[SpatPlot] No spatial locations found in the Giotto object. ")
        }
    }

    if (isTRUE(image)) {
        if ("fov" %in% colnames(GiottoClass::pDataDT(object))) {
            image <- c()
            image_exts <- list()
            for (im in names(object@images)) {
                image_ext <- terra::ext(object@images[[im]])
                if (length(image_exts) == 0 || !any(sapply(image_exts, function(ie) ie == image_ext))) {
                    image <- c(image, im)
                    image_exts <- c(image_exts, list(image_ext))
                }
            }
        } else {
            image <- names(object@images)[1]
        }
    }
    # if (is.null(points) && !is.null(graph)) {
    #     points <- TRUE
    # }
    points <- points %||% TRUE
    shapes <- shapes %||% (if(!is.null(shapes_fill_by)) TRUE)
    layers <- intersect(
        layers %||% c("image", "masks", "shapes", "points"),
        c(
            if (!is.null(image) && !isFALSE(image)) "image",
            if (!is.null(masks) && !isFALSE(masks)) "masks",
            if (!is.null(shapes) && !isFALSE(shapes)) "shapes",
            if (!is.null(points) && !isFALSE(points)) "points"
        )
    )
    stopifnot('Either "image", "masks", "shapes", or "points" must be provided.' =
        any(layers %in% c("image", "masks", "shapes", "points")))

    if (isTRUE(crop) && is.null(ext)) {
        # set the ext to the range of the spatial locations
        ext <- GiottoClass::ext(GiottoClass::getSpatialLocations(
            gobject = object,
            spat_unit = spat_unit
        ))
        if (is.null(padding)) {
            padding <- if ("image" %in% layers) 0 else 0.05
        }
        if (padding > 0) {
            delta_x <- diff(ext[c(1, 2)]) * padding
            delta_y <- diff(ext[c(3, 4)]) * padding
            ext <- terra::ext(c(
                ext[1] - delta_x, ext[2] + delta_x,
                ext[3] - delta_y, ext[4] + delta_y
            ))
        }
    }

    players <- list()
    scales_used <- c()
    facet_by <- NULL
    args <- rlang::dots_list(...)
    shapes_alpha <- shapes_alpha %||% ifelse("image" %in% layers, 0.5, 1)
    for (element in layers) {
        if (element == "image") {
            if (length(image) > 1 || image %in% names(object@images)) {
                image_obj <- GiottoClass::getGiottoImage(gobject = object, name = image)
                if (!is.list(image_obj)) {
                    image_obj <- list(image_obj)
                }
                for (imobj in image_obj) {
                    if (length(image_obj) == 1) {
                        ext <- ext %||% as.vector(terra::ext(imobj))
                    }
                    image_args <- args[startsWith(names(args), "image_")]
                    names(image_args) <- sub("^image_", "", names(image_args))
                    if ("raster_object" %in% methods::slotNames(imobj)) {
                        image_args$data <- imobj@raster_object
                    } else if ("mg_object" %in% methods::slotNames(imobj)) {
                        image_args$data <- terra::rast(as.numeric(imobj@mg_object[[1]]), extent = ext)
                    } else {
                        stop("[SpatPlot] No raster object found in the Giotto image.")
                    }
                    image_args$flip_y <- flip_y
                    image_args$return_layer <- TRUE
                    image_args$ext <- ext
                    image_colors <- attr(imobj, "colors")
                    if (!is.null(image_colors)) { image_args$palcolor <- image_colors }
                    player <- do.call(SpatImagePlot, image_args)
                    scales_reused <- intersect(scales_used, attr(player, "scales"))
                    players <- c(players, .new_scale_layers(scales_reused), list(player))
                    scales_used <- unique(c(scales_used, attr(player, "scales")))
                }
            } else {
                # use as a color to plot a background rect
                player <- .rect_bg_image(args, image)
                scales_reused <- intersect(scales_used, "fill")
                players <- c(players, .new_scale_layers(scales_reused), list(player))
                scales_used <- unique(c(scales_used, attr(player, "scales")))
            }
        } else if (element == "masks") {
            stop("[SpatPlot] 'masks' is not supported for Giotto objects yet.")
        } else if (element == "shapes") {
            # If shapes is given, that means this is a feature plot
            # prepare shapes
            avail_shapes_names <- GiottoClass::list_spatial_info_names(gobject = object)
            if (identical(shapes_feat_type, "cell") && !"cell" %in% avail_shapes_names) {
                shapes_feat_type <- spat_unit
                warning(
                    "[SpatPlot] 'cell' not discovered in shape names. Defaulting to spat_unit (",
                    spat_unit, ")."
                )
            }
            shapes_feat_type <- GiottoClass::set_default_spat_unit(
                gobject = object,
                spat_unit = shapes_feat_type
            )
            feat_type <- GiottoClass::set_default_feat_type(
                gobject = object,
                spat_unit = shapes_feat_type,
                feat_type = feat_type
            )
            shapes_combo <- GiottoClass::combineCellData(
                gobject = object,
                spat_loc_name = spat_loc_name,
                feat_type = feat_type,
                include_poly_info = TRUE,
                poly_info = shapes_feat_type
            )
            shapes_dt <- data.table::rbindlist(shapes_combo, fill = TRUE)

            shapes_args <- args[startsWith(names(args), "shapes_")]
            names(shapes_args) <- sub("^shapes_", "", names(shapes_args))

            shapes_args$x <- "x"
            shapes_args$y <- "y"
            shapes_args$group <- "cell_ID"

            if (!is.null(shapes_fill_by)) {
                # It should be a color name
                if (length(shapes_fill_by) > 1 || shapes_fill_by %in% GiottoClass::featIDs(object, feat_type = shapes_feat_type)) {
                    shapes_dt <- merge(
                        shapes_dt,
                        suppressMessages(GiottoClass::spatValues(
                            object,
                            feats = shapes_fill_by,
                            spat_unit = shapes_feat_type,
                            feat_type = feat_type,
                            spat_enr_name = spat_enr_names
                        )),
                        by = "cell_ID",
                        all.x = TRUE,
                        suffixes = c(".x", "")
                    )
                }
                shapes_args$data <- shapes_dt
                shapes_args$fill_by <- shapes_fill_by
                shapes_args$legend.position <- shapes_args$legend.position %||% "right"
                shapes_args$legend.direction <- shapes_args$legend.direction %||% "vertical"
            } else {
                shapes_args$data <- shapes_dt
            }
            shapes_args$flip_y <- flip_y
            shapes_args$return_layer <- TRUE
            shapes_args$alpha <- shapes_alpha
            shapes_args$ext <- ext %||% c(range(shapes_dt$x, na.rm = TRUE), range(shapes_dt$y, na.rm = TRUE))
            player <- do.call(SpatShapesPlot, shapes_args)
            scales_reused <- intersect(scales_used, attr(player, "scales"))
            players <- c(players, .new_scale_layers(scales_reused), list(player))
            scales_used <- unique(c(scales_used, attr(player, "scales")))

        } else if (element == "points") {
            if (shape != 16) {
                points_args <- .points_args(args, shape = shape)
            } else {
                points_args <- .points_args(args)
            }

            # prepare the data
            if (identical(group_by, "molecules")) {
                # prepare the feature metadata
                if (isTRUE(use_overlap) && !is.null(shapes) && !isFALSE(shapes)) {
                    points_args$data <- tryCatch({
                        suppressWarnings(GiottoClass::combineFeatureOverlapData(
                            gobject = object,
                            feat_type = feat_type,
                            sel_feats = NULL,
                            poly_info = shapes_feat_type
                        ))
                    }, error = function(e) {
                        warning(
                            "[SpatPlot] Error in combineFeatureOverlapData: ",
                            e$message, ". Trying use_overlap = FALSE."
                        )
                        NULL
                    })
                }

                if (is.null(object@feat_info)) {
                    # No feature info, we can only do non-insitu plots
                    points_args$data <- GiottoClass::getSpatialLocations(
                        gobject = object,
                        spat_unit = spat_unit,
                        name = spat_loc_name,
                        output = "data.table"
                    )
                } else {
                    if (is.null(points_args$data)) {
                        points_args$data <- GiottoClass::combineFeatureData(
                            gobject = object,
                            spat_unit = spat_unit,
                            feat_type = feat_type,
                            sel_feats = NULL
                        )
                    }

                    points_args$data <- data.table::rbindlist(points_args$data, fill = TRUE)
                }
            } else {
                # do we need in-situ data?
                # points_args$data <- GiottoClass::combineCellData(
                #     gobject = object,
                #     feat_type = feat_type,
                #     spat_loc_name = spat_loc_name,
                #     spat_enr_name = spat_enr_names,
                #     include_spat_locs = FALSE
                # )[[feat_type]]

                # expression or just locations or group_by
                points_args$data <- GiottoClass::getSpatialLocations(
                    gobject = object,
                    spat_unit = spat_unit,
                    name = spat_loc_name,
                    output = "data.table"
                )

                # also add cell metadata (in case features are cell metadata)
                points_args$data <- merge(
                    points_args$data,
                    GiottoClass::getCellMetadata(
                        gobject = object,
                        spat_unit = spat_unit,
                        feat_type = feat_type,
                        output = "data.table"
                    ),
                    by = "cell_ID",
                    all.x = TRUE,
                    suffixes = c("", ".y")
                )
            }

            if (!is.null(graph)) {
                if (isTRUE(graph)) {
                    points_args$graph <- GiottoClass::getSpatialNetwork(
                        gobject = object,
                        spat_unit = spat_unit,
                        output = "networkDT",
                        verbose = FALSE
                    )
                } else if (is.character(graph) && grepl(":", graph)) {
                    graph_splits <- strsplit(graph, ":")[[1]]
                    points_args$graph <- GiottoClass::getSpatialNetwork(
                        gobject = object,
                        spat_unit = graph_splits[1],
                        name = graph_splits[2],
                        output = "networkDT",
                        verbose = FALSE
                    )
                } else if (is.character(graph)) {
                    points_args$graph <- GiottoClass::getSpatialNetwork(
                        gobject = object,
                        spat_unit = spat_unit,
                        name = graph,
                        output = "networkDT",
                        verbose = FALSE
                    )
                } else {
                    stop("[SpatPlot] 'graph' must be TRUE/NULL or the graph name or 'spat_unit:graph_name'.")
                }
                points_args$graph_x <- "sdimx_begin"
                points_args$graph_y <- "sdimy_begin"
                points_args$graph_xend <- "sdimx_end"
                points_args$graph_yend <- "sdimy_end"
                points_args$graph_value <- "weight"
            }

            if (identical(group_by, "molecules")) {
                if (is.null(features)) {
                    stop("[SpatPlot] 'features' must be provided as molecules when 'group_by' is 'molecules'.")
                }
                points_args$data <- data.table::rbindlist(
                    lapply(split(points_args$data, points_args$data$feat_ID), function(x) {
                        if (x$feat_ID[1] %in% features) {
                            if (nrow(x) > nmols) {
                                x[sample(nrow(x), nmols), ]
                            } else {
                                x
                            }
                        }
                    }), fill = TRUE
                )
                points_args$color_by <- "feat_ID"
                points_args$color_name <- points_args$color_name %||% "Molecules"
                points_args$legend.position <- points_args$legend.position %||% "right"
                points_args$legend.direction <- points_args$legend.direction %||% "vertical"
            } else if (!is.null(features)) {
                expr_feats <- intersect(features, GiottoClass::featIDs(object, feat_type = feat_type))
                if (length(expr_feats) > 0) {
                    expr_df <- t(
                        as.data.frame(
                            GiottoClass::getExpression(
                                gobject = object,
                                spat_unit = spat_unit,
                                feat_type = feat_type,
                                values = layer,
                                output = "matrix"
                            )
                        )
                    )
                    points_args$data[, expr_feats] <- expr_df[points_args$data$cell_ID, expr_feats, drop = FALSE]
                }
                points_args$color_by <- features
                if (length(features) == 1) {
                    points_args$color_name <- points_args$color_name %||% features
                } else {
                    points_args$color_name <- points_args$color_name %||% "Expression"
                    facet_by <- ".facet_var"
                }
                points_args$legend.position <- points_args$legend.position %||% "right"
                points_args$legend.direction <- points_args$legend.direction %||% "vertical"
            } else if (!is.null(group_by)) {
                if (group_by %in% names(points_args$data)) {
                    points_args$data[[group_by]] <- as.factor(points_args$data[[group_by]])
                }
                points_args$color_by <- group_by
                points_args$color_name <- points_args$color_name %||% group_by
                points_args$legend.position <- points_args$legend.position %||% "right"
                points_args$legend.direction <- points_args$legend.direction %||% "vertical"
            }
            points_args$flip_y <- flip_y
            points_args$return_layer <- TRUE
            points_args$ext <- ext
            player <- do.call(SpatPointsPlot, points_args)
            scales_reused <- intersect(scales_used, attr(player, "scales"))
            players <- c(players, .new_scale_layers(scales_reused), list(player))
            scales_used <- unique(c(scales_used, attr(player, "scales")))
        }
    }

    p <- utils::getFromNamespace(".wrap_spatial_layers", "plotthis")(
        layers = players,
        ext = ext,
        flip_y = flip_y,
        legend.position = legend.position,
        legend.direction = legend.direction,
        title = title,
        subtitle = subtitle,
        xlab = xlab,
        ylab = ylab,
        theme = theme,
        theme_args = theme_args
    )

    if (!is.null(facet_by)) {
        p <- utils::getFromNamespace("facet_plot", "plotthis")(
            p, facet_by, facet_scales, facet_nrow, facet_ncol, facet_byrow,
            legend.position = legend.position, legend.direction = legend.direction
        )
    }

    p
}

#' Plot features for spatial data
#'
#' The features can include  expression, dimension reduction components, metadata, etc
#'
#' @inheritParams SpatPlot
#' @return A ggplot object
#' @export
#' @details
#' See:
#' * <https://pwwang.github.io/scplotter/articles/Knowing_your_spatial_data_and_visualization.html> for more details, and
#' * <https://pwwang.github.io/scplotter/articles/Knowing_your_spatial_data_and_visualization.html#examples> for examples.
SpatFeaturePlot <- function(
    object, fov = NULL, boundaries = NULL, image = NULL, masks = NULL, shapes = NULL, points = NULL,
    ext = NULL, crop = TRUE, group_by = NULL, features = NULL, layer = NULL, scale_factor = NULL,
    layers = NULL, flip_y = NULL, padding = NULL, image_scale = NULL,
    x = NULL, y = NULL, nmols = 1000, shapes_fill_by = NULL, graph = NULL, shape = 16,
    legend.position = "right", legend.direction = "vertical", theme = "theme_box", theme_args = list(),
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
    facet_scales = "fixed", facet_nrow = NULL, facet_ncol = NULL, facet_byrow = TRUE,
    feat_type = "rna", use_overlap = FALSE, shapes_feat_type = "cell", shapes_alpha = NULL,
    spat_unit = NULL, spat_loc_name = NULL, spat_enr_names = NULL,
    ...
) {
    UseMethod("SpatFeaturePlot", object)
}

#' @export
SpatFeaturePlot.Seurat <- function(
    object, fov = NULL, boundaries = NULL, image = NULL, masks = NULL, shapes = NULL, points = NULL,
    ext = NULL, crop = TRUE, group_by = NULL, features = NULL, layer = NULL, scale_factor = NULL,
    layers = NULL, flip_y = NULL, padding = NULL, image_scale = NULL,
    x = NULL, y = NULL, nmols = 1000, shapes_fill_by = NULL, graph = NULL, shape = 16,
    legend.position = "right", legend.direction = "vertical", theme = "theme_box", theme_args = list(),
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
    facet_scales = "fixed", facet_nrow = NULL, facet_ncol = NULL, facet_byrow = TRUE,
    feat_type = "rna", use_overlap = FALSE, shapes_feat_type = "cell", shapes_alpha = NULL,
    spat_unit = NULL, spat_loc_name = NULL, spat_enr_names = NULL,
    ...
) {
    SpatPlot(
        object, fov = fov, boundaries = boundaries, image = image, masks = masks, shapes = shapes,
        points = points, ext = ext, crop = crop, group_by = group_by, features = features,
        layer = layer, scale_factor = scale_factor, layers = layers, flip_y = flip_y,
        padding = padding, image_scale = image_scale,
        x = x, y = y, nmols = nmols, shapes_fill_by = shapes_fill_by, graph = graph, shape = shape,
        legend.position = legend.position, legend.direction = legend.direction, theme = theme, theme_args = theme_args,
        title = title, subtitle = subtitle, xlab = xlab, ylab = ylab,
        facet_scales = facet_scales, facet_nrow = facet_nrow,
        facet_ncol = facet_ncol, facet_byrow = facet_byrow,
        feat_type = feat_type, use_overlap = use_overlap,
        shapes_feat_type = shapes_feat_type, shapes_alpha = shapes_alpha,
        spat_unit = spat_unit, spat_loc_name = spat_loc_name,
        spat_enr_names = spat_enr_names,
        ...
    )
}

#' @export
SpatFeaturePlot.giotto <- function(
    object, fov = NULL, boundaries = NULL, image = NULL, masks = NULL, shapes = NULL, points = NULL,
    ext = NULL, crop = TRUE, group_by = NULL, features = NULL, layer = NULL, scale_factor = NULL,
    layers = NULL, flip_y = NULL, padding = NULL, image_scale = NULL,
    x = NULL, y = NULL, nmols = 1000, shapes_fill_by = NULL, graph = NULL, shape = 16,
    legend.position = "right", legend.direction = "vertical", theme = "theme_box", theme_args = list(),
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
    facet_scales = "fixed", facet_nrow = NULL, facet_ncol = NULL, facet_byrow = TRUE,
    feat_type = "rna", use_overlap = FALSE, shapes_feat_type = "cell", shapes_alpha = NULL,
    spat_unit = NULL, spat_loc_name = NULL, spat_enr_names = NULL,
    ...
) {
    stopifnot(
        "[SpatFeaturePlot] 'group_by' is not supported for Giotto objects. Use 'SpatDimPlot' instead." =
        is.null(group_by)
    )
    SpatPlot.giotto(
        object, fov = fov, boundaries = boundaries, image = image, masks = masks, shapes = shapes,
        points = points, ext = ext, crop = crop, group_by = group_by, features = features,
        layer = layer, scale_factor = scale_factor, layers = layers, flip_y = flip_y,
        padding = padding, image_scale = image_scale,
        x = x, y = y, nmols = nmols, shapes_fill_by = shapes_fill_by, graph = graph, shape = shape,
        legend.position = legend.position, legend.direction = legend.direction,
        theme = theme, theme_args = theme_args,
        title = title, subtitle = subtitle, xlab = xlab, ylab = ylab,
        facet_scales = facet_scales, facet_nrow = facet_nrow,
        facet_ncol = facet_ncol, facet_byrow = facet_byrow,
        feat_type = feat_type, use_overlap = use_overlap,
        shapes_feat_type = shapes_feat_type, shapes_alpha = shapes_alpha,
        spat_unit = spat_unit, spat_loc_name = spat_loc_name,
        spat_enr_names = spat_enr_names,
        ...
    )
}

#' Plot categories for spatial data
#'
#' @inheritParams SpatPlot
#' @return A ggplot object
#' @export
#' @details
#' See:
#' * <https://pwwang.github.io/scplotter/articles/Knowing_your_spatial_data_and_visualization.html> for more details, and
#' * <https://pwwang.github.io/scplotter/articles/Knowing_your_spatial_data_and_visualization.html#examples> for examples.
SpatDimPlot <- function(
    object, fov = NULL, boundaries = NULL, image = NULL, masks = NULL, shapes = NULL, points = NULL,
    ext = NULL, crop = TRUE, group_by = NULL, features = NULL, layer = NULL, scale_factor = NULL,
    layers = NULL, flip_y = NULL, padding = NULL, image_scale = NULL,
    x = NULL, y = NULL, nmols = 1000, shapes_fill_by = NULL, graph = NULL, shape = 16,
    legend.position = "right", legend.direction = "vertical", theme = "theme_box", theme_args = list(),
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
    facet_scales = "fixed", facet_nrow = NULL, facet_ncol = NULL, facet_byrow = TRUE,
    feat_type = "rna", use_overlap = FALSE, shapes_feat_type = "cell", shapes_alpha = NULL,
    spat_unit = NULL, spat_loc_name = NULL, spat_enr_names = NULL,
    ...
) {
    UseMethod("SpatDimPlot", object)
}

#' @export
SpatDimPlot.Seurat <- function(
    object, fov = NULL, boundaries = NULL, image = NULL, masks = NULL, shapes = NULL, points = NULL,
    ext = NULL, crop = TRUE, group_by = NULL, features = NULL, layer = NULL, scale_factor = NULL,
    layers = NULL, flip_y = NULL, padding = NULL, image_scale = NULL,
    x = NULL, y = NULL, nmols = 1000, shapes_fill_by = NULL, graph = NULL, shape = 16,
    legend.position = "right", legend.direction = "vertical", theme = "theme_box", theme_args = list(),
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
    facet_scales = "fixed", facet_nrow = NULL, facet_ncol = NULL, facet_byrow = TRUE,
    feat_type = "rna", use_overlap = FALSE, shapes_feat_type = "cell", shapes_alpha = NULL,
    spat_unit = NULL, spat_loc_name = NULL, spat_enr_names = NULL,
    ...
) {
    if (
        "images" %in% methods::slotNames(object) &&
        length(object@images) > 0 &&
        "FOV" %in% class(object@images[[1]]) &&
        !is.null(features)
    ) {
        group_by <- group_by %||% "molecules"
    } else if (is.null(group_by)) {
        group_by <- "Identity"
        object@meta.data$Identity <- Seurat::Idents(object)
    }
    SpatPlot.Seurat(
        object, fov = fov, boundaries = boundaries, image = image, masks = masks, shapes = shapes,
        points = points, ext = ext, crop = crop, group_by = group_by, features = features,
        layer = layer, scale_factor = scale_factor, layers = layers, flip_y = flip_y,
        padding = padding, image_scale = image_scale,
        x = x, y = y, nmols = nmols, shapes_fill_by = shapes_fill_by, graph = graph, shape = shape,
        legend.position = legend.position, legend.direction = legend.direction,
        theme = theme, theme_args = theme_args,
        title = title, subtitle = subtitle, xlab = xlab, ylab = ylab,
        facet_scales = facet_scales, facet_nrow = facet_nrow,
        facet_ncol = facet_ncol, facet_byrow = facet_byrow,
        feat_type = feat_type, use_overlap = use_overlap,
        shapes_feat_type = shapes_feat_type, shapes_alpha = shapes_alpha,
        spat_unit = spat_unit, spat_loc_name = spat_loc_name,
        spat_enr_names = spat_enr_names,
        ...
    )
}

#' @export
SpatDimPlot.giotto <- function(
    object, fov = NULL, boundaries = NULL, image = NULL, masks = NULL, shapes = NULL, points = NULL,
    ext = NULL, crop = TRUE, group_by = NULL, features = NULL, layer = NULL, scale_factor = NULL,
    layers = NULL, flip_y = NULL, padding = NULL, image_scale = NULL,
    x = NULL, y = NULL, nmols = 1000, shapes_fill_by = NULL, graph = NULL, shape = 16,
    legend.position = "right", legend.direction = "vertical", theme = "theme_box", theme_args = list(),
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
    facet_scales = "fixed", facet_nrow = NULL, facet_ncol = NULL, facet_byrow = TRUE,
    feat_type = "rna", use_overlap = FALSE, shapes_feat_type = "cell", shapes_alpha = NULL,
    spat_unit = NULL, spat_loc_name = NULL, spat_enr_names = NULL,
    ...
) {
    if (!is.null(features)) {
        group_by <- group_by %||% "molecules"
    }
    SpatPlot.giotto(
        object, fov = fov, boundaries = boundaries, image = image, masks = masks, shapes = shapes,
        points = points, ext = ext, crop = crop, group_by = group_by, features = features,
        layer = layer, scale_factor = scale_factor, layers = layers, flip_y = flip_y,
        padding = padding, image_scale = image_scale,
        x = x, y = y, nmols = nmols, shapes_fill_by = shapes_fill_by, graph = graph, shape = shape,
        legend.position = legend.position, legend.direction = legend.direction,
        theme = theme, theme_args = theme_args,
        title = title, subtitle = subtitle, xlab = xlab, ylab = ylab,
        facet_scales = facet_scales, facet_nrow = facet_nrow,
        facet_ncol = facet_ncol, facet_byrow = facet_byrow,
        feat_type = feat_type, use_overlap = use_overlap,
        shapes_feat_type = shapes_feat_type, shapes_alpha = shapes_alpha,
        spat_unit = spat_unit, spat_loc_name = spat_loc_name,
        spat_enr_names = spat_enr_names,
        ...
    )
}
