#' Convert a data.frame to a terra SpatVector polygon
#' @keywords internal
.df_to_polygon_spatvector <- function(
    df, id, x, y,
    attrs = NULL
) {
    shapes <- split(df[, c(x, y), drop = FALSE], df[[id]])
    sv <- terra::vect(lapply(shapes, as.matrix), type = "polygons")
    if (!is.null(attrs)) {
        for (attr in attrs) {
            # only get the first value of each polygon
            sv[[attr]] <- sapply(split(df[[attr]], df[[id]]), function(x) x[1])
        }
    }
    sv
}

#' Process points layer for Seurat spatial plots
#' @keywords internal
.seurat_points_layer <- function(
    object, points_x = "x", points_y = "y", image, args, crop, points_data, ext_unscaled,
    scale_factor, group_by, features, layer, legend.position, legend.direction,
    label, label_size, label_fg, label_bg, label_bg_r, label_repel, label_repulsion,
    label_pt_size, label_pt_color, label_segment_color, label_insitu,
    highlight, highlight_alpha, highlight_size, highlight_color, highlight_stroke,
    palette, palette_reverse, palcolor, flip_y, ext
) {
    points_args <- args[startsWith(names(args), "points_")]
    names(points_args) <- sub("^points_", "", names(points_args))

    if (crop) {
        # attach metadata for highlighting selection
        points_args$data <- object@meta.data[rownames(points_data), , drop = FALSE]
        # Handle missing columns for VisiumV2 vs SlideSeq difference
        if (ncol(points_args$data) > 0) {
            points_args$data <- points_args$data[, colnames(points_args$data)[!is.na(colnames(points_args$data))], drop = FALSE]
        }
        points_args$data <- cbind(points_args$data, points_data)
        points_args$data[[points_x]] <- points_args$data[[points_x]] * scale_factor
        points_args$data[[points_y]] <- points_args$data[[points_y]] * scale_factor
        points_args$ext <- ext %||% (ext_unscaled * scale_factor)
    } else {
        points_args$data <- Seurat::GetTissueCoordinates(object, image = if(isFALSE(image)) NULL else image)
        points_args$data <- points_args$data[, colnames(points_args$data)[!is.na(colnames(points_args$data))], drop = FALSE]
        points_args$data$.tmp <- points_args$data[[points_x]]
        points_args$data[[points_x]] <- points_args$data[[points_y]] * scale_factor
        points_args$data[[points_y]] <- points_args$data$.tmp * scale_factor
        points_args$data$.tmp <- NULL
        points_args$data <- cbind(
            object@meta.data[rownames(points_args$data), , drop = FALSE],
            points_args$data
        )
    }

    facet_by <- NULL

    if (!is.null(group_by)) {
        points_args$data[[group_by]] <- object@meta.data[[group_by]]
        points_args$color_by <- group_by
        points_args$label <- points_args$label %||% label
        points_args$label_size <- points_args$label_size %||% label_size
        points_args$label_fg <- points_args$label_fg %||% label_fg
        points_args$label_bg <- points_args$label_bg %||% label_bg
        points_args$label_bg_r <- points_args$label_bg_r %||% label_bg_r
        points_args$label_repel <- points_args$label_repel %||% label_repel
        points_args$label_repulsion <- points_args$label_repulsion %||% label_repulsion
        points_args$label_pt_size <- points_args$label_pt_size %||% label_pt_size
        points_args$label_pt_color <- points_args$label_pt_color %||% label_pt_color
        points_args$label_segment_color <- points_args$label_segment_color %||% label_segment_color
        points_args$label_insitu <- points_args$label_insitu %||% label_insitu
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

    points_args$x <- points_x
    points_args$y <- points_y
    points_args$highlight <- points_args$highlight %||% highlight
    points_args$highlight_alpha <- points_args$highlight_alpha %||% highlight_alpha
    points_args$highlight_size <- points_args$highlight_size %||% highlight_size
    points_args$highlight_color <- points_args$highlight_color %||% highlight_color
    points_args$highlight_stroke <- points_args$highlight_stroke %||% highlight_stroke
    points_args$palette <- points_args$palette %||% palette
    points_args$palette_reverse <- points_args$palette_reverse %||% palette_reverse
    points_args$palcolor <- points_args$palcolor %||% palcolor
    points_args$legend.position <- legend.position
    points_args$legend.direction <- legend.direction
    points_args$flip_y <- points_args$flip_y
    points_args$return_layer <- TRUE
    player <- do.call(SpatialPointsPlot, points_args)

    list(player = player, facet_by = facet_by)
}

#' Add new scale layers
#' @keywords internal
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

#' Plot features for spatial data
#'
#' The features can include  expression, dimension reduction components, metadata, etc
#'
#' @param object A Seurat object or a Giotto object.
#' @param image The name of the image to plot. If NULL, the first image in the object will be used.
#' @param ... Additional arguments passed to the plotting functions.
#' @return A ggplot object
#' @keywords internal
#' @rdname SpatialPlot
SpatPlot <- function(object, ...) {
    UseMethod("SpatPlot", object)
}

#' @keywords internal
#' @rdname SpatialPlot
SpatPlot.Seurat <- function(object, image = NULL, ...) {
    first_image <- Seurat::Images(object)[1]
    if (is.null(first_image)) {
        stop("[SpatPlot] No images found in the Seurat object. Is this an object with spatial data? ")
    }
    image <- image %||% first_image
    stype <- class(object@images[[first_image]])
    if ("VisiumV1" %in% stype) {
        SpatPlot.Seurat.Visium(object, image = image, points_x = "imagerow", points_y = "imagecol", scale_factor = 1, ...)
    } else if ("VisiumV2" %in% stype) {
        SpatPlot.Seurat.Visium(object, image = image, ...)
    } else if ("FOV" %in% stype) {
        SpatPlot.Seurat.FOV(object, image = image, ...)
    } else if ("SlideSeq" %in% stype) {
        SpatPlot.Seurat.SlideSeq(object, image = NULL, ...)
    } else {
        stop("[SpatPlot] Unsupported image type. Supported types are 'VisiumV1', 'VisiumV2', and 'SlideSeq'.")
    }
}

#' @keywords internal
#' @rdname SpatialPlot
#' @importFrom plotthis SpatialImagePlot SpatialMasksPlot SpatialShapesPlot SpatialPointsPlot
SpatPlot.Seurat.Visium <- function(
    object, image = NULL, masks = NULL, shapes = NULL, points = NULL, ext = NULL, points_x = "x", points_y = "y",
    image_scale = NULL, crop = TRUE, group_by = NULL, features = NULL, layer = "data", scale_factor = NULL,
    layers = NULL, flip_y = TRUE, theme = "theme_box", theme_args = list(),
    label = FALSE, label_size = 4, label_fg = "white", label_bg = "black", label_bg_r = 0.1,
    label_repel = FALSE, label_repulsion = 20, label_pt_size = 1, label_pt_color = "black",
    label_segment_color = "black", label_insitu = FALSE,
    palette = NULL, palette_reverse = FALSE, palcolor = NULL,
    highlight = NULL, highlight_alpha = 1, highlight_size = 1, highlight_color = "black", highlight_stroke = 0.8,
    legend.position = "right", legend.direction = "vertical",
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
    facet_scales = "fixed", facet_nrow = NULL, facet_ncol = NULL, facet_byrow = TRUE,
    ...

) {
    stopifnot("[SpatPlot] Either 'group_by' or 'features' should be provided, not both." = is.null(group_by) || is.null(features))
    if (is.character(theme) && theme %in% c("theme_box", "theme_this", "theme_blank")) {
        theme <- utils::getFromNamespace(theme, "plotthis")
    }

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
    stopifnot('Either "image", "masks", "shapes", or "points" must be provided.' = any(layers %in% c("image", "masks", "shapes", "points")))
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
        points_data$.tmp <- points_data[[points_x]]
        points_data[[points_x]] <- points_data[[points_y]]
        points_data[[points_y]] <- points_data$.tmp
        points_data$.tmp <- NULL
        padding <- 0.05
        delta_x <- diff(range(points_data[[points_x]], na.rm = TRUE)) * padding
        delta_y <- diff(range(points_data[[points_y]], na.rm = TRUE)) * padding
        ext_unscaled <- c(
            min(points_data[[points_x]], na.rm = TRUE) - delta_x,
            max(points_data[[points_x]], na.rm = TRUE) + delta_x,
            min(points_data[[points_y]], na.rm = TRUE) - delta_y,
            max(points_data[[points_y]], na.rm = TRUE) + delta_y
        )
    }
    for (element in layers) {
        if (element == "image") {
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
            player <- do.call(SpatialImagePlot, image_args)
            scales_reused <- intersect(scales_used, attr(player, "scales"))
            players <- c(players, list(player))
            scales_used <- unique(c(scales_used, attr(player, "scales")))
        }
        if (element == "points") {
            points_layer <- .seurat_points_layer(
                object = object, image = image, args = args, crop = crop, points_data = points_data,
                points_x = points_x, points_y = points_y,
                ext_unscaled = ext_unscaled, scale_factor = the_scale_factor, group_by = group_by,
                features = features, layer = layer, legend.position = legend.position,
                legend.direction = legend.direction, label = label, label_size = label_size,
                label_fg = label_fg, label_bg = label_bg, label_bg_r = label_bg_r,
                label_repel = label_repel, label_repulsion = label_repulsion,
                label_pt_size = label_pt_size, label_pt_color = label_pt_color,
                label_segment_color = label_segment_color, label_insitu = label_insitu,
                highlight = highlight, highlight_alpha = highlight_alpha,
                highlight_size = highlight_size, highlight_color = highlight_color,
                highlight_stroke = highlight_stroke, palette = palette,
                palette_reverse = palette_reverse, palcolor = palcolor,
                flip_y = flip_y, ext = ext
            )
            player <- points_layer$player
            facet_by <- points_layer$facet_by
            scales_reused <- intersect(scales_used, attr(player, "scales"))
            if ("fill" %in% scales_reused) {
                players <- c(players, list(ggnewscale::new_scale_fill()))
            }
            players <- c(players, list(player))
            scales_used <- unique(c(scales_used, attr(player, "scales")))
        }
        if (element == "masks") {
            stop("[SpatPlot] 'masks' is not supported for Seurat Visium V2 objects.")
        }
        if (element == "shapes") {
            stop("[SpatPlot] 'shapes' is not supported for Seurat Visium V2 objects.")
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
#' @rdname SpatialPlot
SpatPlot.Seurat.SlideSeq <- function(
    object, image = NULL, masks = NULL, shapes = NULL, points = NULL, ext = NULL,
    image_scale = NULL, crop = TRUE, group_by = NULL, features = NULL, layer = "data",
    layers = NULL, flip_y = TRUE, theme = "theme_box", theme_args = list(),
    label = FALSE, label_size = 4, label_fg = "white", label_bg = "black", label_bg_r = 0.1,
    label_repel = FALSE, label_repulsion = 20, label_pt_size = 1, label_pt_color = "black",
    label_segment_color = "black", label_insitu = FALSE,
    palette = NULL, palette_reverse = FALSE, palcolor = NULL,
    highlight = NULL, highlight_alpha = 1, highlight_size = 1, highlight_color = "black", highlight_stroke = 0.8,
    legend.position = "right", legend.direction = "vertical",
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
    facet_scales = "fixed", facet_nrow = NULL, facet_ncol = NULL, facet_byrow = TRUE,
    ...

) {
    stopifnot("[SpatPlot] Either 'group_by' or 'features' should be provided, not both." = is.null(group_by) || is.null(features))
    if (is.character(theme) && theme %in% c("theme_box", "theme_this", "theme_blank")) {
        theme <- utils::getFromNamespace(theme, "plotthis")
    }

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
    stopifnot('Either "image", "masks", "shapes", or "points" must be provided.' = any(layers %in% c("image", "masks", "shapes", "points")))

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
        padding <- 0.05
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
            stop("[SpatPlot] 'image' is not supported for Seurat SlideSeq objects.")
        }
        if (element == "points") {
            points_layer <- .seurat_points_layer(
                object = object, image = image, args = args, crop = crop, points_data = points_data,
                ext_unscaled = ext_unscaled, scale_factor = scale_factor, group_by = group_by,
                features = features, layer = layer, legend.position = legend.position,
                legend.direction = legend.direction, label = label, label_size = label_size,
                label_fg = label_fg, label_bg = label_bg, label_bg_r = label_bg_r,
                label_repel = label_repel, label_repulsion = label_repulsion,
                label_pt_size = label_pt_size, label_pt_color = label_pt_color,
                label_segment_color = label_segment_color, label_insitu = label_insitu,
                highlight = highlight, highlight_alpha = highlight_alpha,
                highlight_size = highlight_size, highlight_color = highlight_color,
                highlight_stroke = highlight_stroke, palette = palette,
                palette_reverse = palette_reverse, palcolor = palcolor,
                flip_y = flip_y, ext = ext
            )
            player <- points_layer$player
            facet_by <- points_layer$facet_by
            scales_reused <- intersect(scales_used, attr(player, "scales"))
            if ("fill" %in% scales_reused) {
                players <- c(players, list(ggnewscale::new_scale_fill()))
            }
            players <- c(players, list(player))
            scales_used <- unique(c(scales_used, attr(player, "scales")))
        }
        if (element == "masks") {
            stop("[SpatPlot] 'masks' is not supported for Seurat Visium V2 objects.")
        }
        if (element == "shapes") {
            stop("[SpatPlot] 'shapes' is not supported for Seurat Visium V2 objects.")
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
#' @rdname SpatialPlot
#' @importFrom plotthis SpatialImagePlot SpatialMasksPlot SpatialShapesPlot SpatialPointsPlot
SpatPlot.Seurat.FOV <- function(
    object, image = NULL, masks = NULL, shapes = NULL, points = NULL, ext = NULL, points_x = "x", points_y = "y",
    image_scale = NULL, crop = TRUE, group_by = NULL, features = NULL, layer = "data", scale_factor = NULL,
    layers = NULL, flip_y = TRUE, theme = "theme_box", theme_args = list(),
    label = FALSE, label_size = 4, label_fg = "white", label_bg = "black", label_bg_r = 0.1,
    label_repel = FALSE, label_repulsion = 20, label_pt_size = 1, label_pt_color = "black",
    label_segment_color = "black", label_insitu = FALSE,
    palette = NULL, palette_reverse = FALSE, palcolor = NULL,
    highlight = NULL, highlight_alpha = 1, highlight_size = 1, highlight_color = "black", highlight_stroke = 0.8,
    legend.position = "right", legend.direction = "vertical",
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
    facet_scales = "fixed", facet_nrow = NULL, facet_ncol = NULL, facet_byrow = TRUE,
    ...

) {
    stopifnot("[SpatPlot] Either 'group_by' or 'features' should be provided, not both." = is.null(group_by) || is.null(features))
    if (is.character(theme) && theme %in% c("theme_box", "theme_this", "theme_blank")) {
        theme <- utils::getFromNamespace(theme, "plotthis")
    }

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
    stopifnot('Either "image", "masks", "shapes", or "points" must be provided.' = any(layers %in% c("image", "masks", "shapes", "points")))
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
        points_data$.tmp <- points_data[[points_x]]
        points_data[[points_x]] <- points_data[[points_y]]
        points_data[[points_y]] <- points_data$.tmp
        points_data$.tmp <- NULL
        padding <- 0.05
        delta_x <- diff(range(points_data[[points_x]], na.rm = TRUE)) * padding
        delta_y <- diff(range(points_data[[points_y]], na.rm = TRUE)) * padding
        ext_unscaled <- c(
            min(points_data[[points_x]], na.rm = TRUE) - delta_x,
            max(points_data[[points_x]], na.rm = TRUE) + delta_x,
            min(points_data[[points_y]], na.rm = TRUE) - delta_y,
            max(points_data[[points_y]], na.rm = TRUE) + delta_y
        )
    }
    for (element in layers) {
        if (element == "image") {
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
            player <- do.call(SpatialImagePlot, image_args)
            scales_reused <- intersect(scales_used, attr(player, "scales"))
            players <- c(players, list(player))
            scales_used <- unique(c(scales_used, attr(player, "scales")))
        }
        if (element == "points") {
            points_layer <- .seurat_points_layer(
                object = object, image = image, args = args, crop = crop, points_data = points_data,
                points_x = points_x, points_y = points_y,
                ext_unscaled = ext_unscaled, scale_factor = the_scale_factor, group_by = group_by,
                features = features, layer = layer, legend.position = legend.position,
                legend.direction = legend.direction, label = label, label_size = label_size,
                label_fg = label_fg, label_bg = label_bg, label_bg_r = label_bg_r,
                label_repel = label_repel, label_repulsion = label_repulsion,
                label_pt_size = label_pt_size, label_pt_color = label_pt_color,
                label_segment_color = label_segment_color, label_insitu = label_insitu,
                highlight = highlight, highlight_alpha = highlight_alpha,
                highlight_size = highlight_size, highlight_color = highlight_color,
                highlight_stroke = highlight_stroke, palette = palette,
                palette_reverse = palette_reverse, palcolor = palcolor,
                flip_y = flip_y, ext = ext
            )
            player <- points_layer$player
            facet_by <- points_layer$facet_by
            scales_reused <- intersect(scales_used, attr(player, "scales"))
            if ("fill" %in% scales_reused) {
                players <- c(players, list(ggnewscale::new_scale_fill()))
            }
            players <- c(players, list(player))
            scales_used <- unique(c(scales_used, attr(player, "scales")))
        }
        if (element == "masks") {
            stop("[SpatPlot] 'masks' is not supported for Seurat Visium V2 objects.")
        }
        if (element == "shapes") {
            stop("[SpatPlot] 'shapes' is not supported for Seurat Visium V2 objects.")
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
#' @rdname SpatialPlot
SpatPlot.giotto <- function(
    object,
    image = NULL,
    masks = NULL,
    shapes = NULL,
    points = NULL,
    ext = NULL,
    layer = c("normalized", "scaled", "custom", "raw"),
    layers = NULL,
    flip_y = FALSE,
    group_by = "shapes",
    feat_type = "rna",
    features = NULL,
    nmols = 1000,
    use_overlap = FALSE,
    spat_unit = NULL,
    spat_loc_name = NULL,
    spat_enr_names = NULL,
    shapes_feat_type = "cell",
    shapes_fill_by = NULL,
    shapes_alpha = NULL,
    theme = "theme_box", theme_args = list(),
    legend.position = "right", legend.direction = "vertical",
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
    facet_scales = "fixed", facet_nrow = NULL, facet_ncol = NULL, facet_byrow = TRUE,
    ...) {
    if (!identical(group_by, "molecules") && !identical(group_by, "shapes")) {
        stopifnot(
            "[SpatPlot] Either 'group_by' or 'features' should be provided, not both." =
            is.null(group_by) || is.null(features))
    }
    layer <- match.arg(layer)
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
        image <- names(object@images)[1]
    }
    points <- points %||% TRUE
    shapes <- shapes %||% (if(identical(group_by, "shapes") || !is.null(shapes_fill_by)) TRUE)
    layers <- intersect(
        layers %||% c("image", "masks", "shapes", "points"),
        c(
            if (!is.null(image) && !isFALSE(image)) "image",
            if (!is.null(masks) && !isFALSE(masks)) "masks",
            if (!is.null(shapes) && !isFALSE(shapes)) "shapes",
            if (!is.null(points) && !isFALSE(points)) "points"
        )
    )
    stopifnot('Either "image", "masks", "shapes", or "points" must be provided.' = any(layers %in% c("image", "masks", "shapes", "points")))

    players <- list()
    scales_used <- c()
    facet_by <- NULL
    args <- rlang::dots_list(...)
    shapes_alpha <- shapes_alpha %||% ifelse("image" %in% layers, 0.5, 1)
    for (element in layers) {
        if (element == "image") {
            image_obj <- GiottoClass::getGiottoImage(gobject = object, name = image)
            ext <- ext %||% as.vector(terra::ext(image_obj))
            image_args <- args[startsWith(names(args), "image_")]
            names(image_args) <- sub("^image_", "", names(image_args))
            image_args$data <- image_obj@raster_object
            image_args$flip_y <- flip_y
            image_args$return_layer <- TRUE
            image_args$ext <- ext
            image_colors <- attr(image_obj, "colors")
            if (!is.null(image_colors)) { image_args$palcolor <- image_colors }
            player <- do.call(SpatialImagePlot, image_args)
            scales_reused <- intersect(scales_used, attr(player, "scales"))
            players <- c(players, .new_scale_layers(scales_reused), list(player))
            scales_used <- unique(c(scales_used, attr(player, "scales")))
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

            if (!is.null(shapes_fill_by)) {
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
                shapes_args$data <- .df_to_polygon_spatvector(shapes_dt, "cell_ID", "x", "y", attrs = shapes_fill_by)
                shapes_args$fill_by <- shapes_fill_by
                shapes_args$legend.position <- shapes_args$legend.position %||% "right"
                shapes_args$legend.direction <- shapes_args$legend.direction %||% "vertical"
            } else {
                shapes_args$data <- .df_to_polygon_spatvector(shapes_dt, "cell_ID", "x", "y")
            }
            shapes_args$flip_y <- flip_y
            shapes_args$return_layer <- TRUE
            shapes_args$alpha <- shapes_alpha
            shapes_args$ext <- ext %||% as.vector(terra::ext(shapes_args$data))
            player <- do.call(SpatialShapesPlot, shapes_args)
            scales_reused <- intersect(scales_used, attr(player, "scales"))
            players <- c(players, .new_scale_layers(scales_reused), list(player))
            scales_used <- unique(c(scales_used, attr(player, "scales")))

        } else if (element == "points") {
            points_args <- args[startsWith(names(args), "points_")]
            names(points_args) <- sub("^points_", "", names(points_args))

            # prepare the data
            if (identical(group_by, "shapes") || identical(group_by, "molecules")) {
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

                if (is.null(points_args$data)) {
                    points_args$data <- GiottoClass::combineFeatureData(
                        gobject = object,
                        spat_unit = spat_unit,
                        feat_type = feat_type,
                        sel_feats = NULL
                    )
                }

                points_args$data <- data.table::rbindlist(points_args$data, fill = TRUE)
            } else if (!is.null(group_by)) {
                # prepare the cell metadata
            } else {  # expression or just locations
                points_args$data <- GiottoClass::getSpatialLocations(
                    gobject = object,
                    spat_unit = spat_unit,
                    name = spat_loc_name,
                    output = "data.table"
                )
            }

            if (identical(group_by, "shapes")) {
                points_args$label <- FALSE
            } else if (identical(group_by, "molecules")) {
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
                points_args$label <- FALSE
                points_args$color_by <- "feat_ID"
                points_args$color_name <- points_args$color_name %||% "Molecules"
                points_args$legend.position <- points_args$legend.position %||% "right"
                points_args$legend.direction <- points_args$legend.direction %||% "vertical"
            } else if (!is.null(features)) {
                expr_df <- t(
                    as.data.frame(
                        GiottoClass::getExpression(
                            gobject = object,
                            spat_unit = spat_unit,
                            feat_type = feat_type,
                            values = layer,
                            output = "matrix"
                        )
                    )[features, , drop = FALSE]
                )
                points_args$data[, features] <- expr_df[points_args$data$cell_ID, features, drop = FALSE]
                points_args$color_by <- features
                if (length(features) == 1) {
                    points_args$color_name <- points_args$color_name %||% features
                } else {
                    points_args$color_name <- points_args$color_name %||% "Expression"
                    facet_by <- ".facet_var"
                }
                points_args$legend.position <- points_args$legend.position %||% "right"
                points_args$legend.direction <- points_args$legend.direction %||% "vertical"
            }
            points_args$flip_y <- flip_y
            points_args$return_layer <- TRUE
            points_args$ext <- ext
            player <- do.call(SpatialPointsPlot, points_args)
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
#' @param object A Seurat object or a Giotto object.
#' @return A ggplot object
#' @export
#' @rdname SpatialPlot
SpatFeaturePlot <- function(object, ...) {
    UseMethod("SpatFeaturePlot", object)
}

#' @export
#' @rdname SpatialPlot
SpatFeaturePlot.Seurat <- function(object, image = NULL, ...) {
    SpatPlot(object, image = image, ...)
}

#' @export
#' @rdname SpatialPlot
SpatFeaturePlot.giotto <- function(object, image = NULL, group_by = NULL, ...) {
    stopifnot(
        "[SpatFeaturePlot] 'group_by' is not supported for Giotto objects. Use 'SpatDimPlot' instead." =
        is.null(group_by)
    )
    SpatPlot(object, image = image, group_by = NULL, ...)
}

#' Plot categories for spatial data
#'
#' @param object A Seurat object or a Giotto object.
#' @return A ggplot object
#' @export
#' @rdname SpatialPlot
SpatDimPlot <- function(object, ...) {
    UseMethod("SpatDimPlot", object)
}

#' @export
#' @rdname SpatialPlot
SpatDimPlot.Seurat <- function(object, image = NULL, group_by = NULL, ...) {
    if (is.null(group_by)) {
        group_by <- "Identity"
        object@meta.data$Identity <- Seurat::Idents(object)
    }
    SpatPlot.Seurat(object, image = image, group_by = group_by, ...)
}

#' @export
#' @rdname SpatialPlot
SpatDimPlot.giotto <- function(object, image = NULL, group_by = NULL, ...) {
    group_by <- group_by %||% "shapes"
    SpatPlot.giotto(object, image = image, group_by = group_by, ...)
}
