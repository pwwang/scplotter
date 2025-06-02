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
SpatialPlot <- function(object, ...) {
    UseMethod("SpatialPlot", object)
}

#' @keywords internal
#' @rdname SpatialPlot
SpatialPlot.Seurat <- function(object, image = NULL, ...) {
    first_image <- Seurat::Images(object)[1]
    if (is.null(first_image)) {
        stop("[SpatialPlot] No images found in the Seurat object. Is this an object with spatial data? ")
    }
    image <- image %||% first_image
    stype <- class(object@images[[first_image]])
    if ("VisiumV2" %in% stype) {
        SpatialPlot.Seurat.VisiumV2(object, image = image, ...)
    } else if ("SlideSeq" %in% stype) {
        SpatialPlot.Seurat.SlideSeq(object, image = NULL, ...)
    }
}

#' @keywords internal
#' @rdname SpatialPlot
#' @importFrom plotthis SpatialImagePlot SpatialMasksPlot SpatialShapesPlot SpatialPointsPlot
SpatialPlot.Seurat.VisiumV2 <- function(
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
    stopifnot("[SpatialPlot] Either 'group_by' or 'features' should be provided, not both." = is.null(group_by) || is.null(features))
    if (theme %in% c("theme_box", "theme_this", "theme_blank")) {
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
        if (element == "image" && !is.null(image)) {
            image_obj <- object@images[[image]]
            img <- terra::rast(image_obj@image)
            image_scale <- image_scale %||% which.min(image_obj@scale.factors)
            scale_factor <- image_obj@scale.factors[[image_scale]]
            # img <- terra::crop(img, ext, extend = TRUE)
            image_args <- args[startsWith(names(args), "image_")]
            names(image_args) <- sub("^image_", "", names(image_args))
            image_args$data <- terra::flip(img, direction = "vertical")
            image_args$flip_y <- flip_y
            image_args$return_layer <- TRUE
            image_args$ext <- if (!is.null(ext_unscaled)) ext_unscaled * scale_factor
            player <- do.call(SpatialImagePlot, image_args)
            scales_reused <- intersect(scales_used, attr(player, "scales"))
            players <- c(players, list(player))
            scales_used <- unique(c(scales_used, attr(player, "scales")))
        }
        if (element == "points" && !is.null(points)) {
            points_args <- args[startsWith(names(args), "points_")]
            names(points_args) <- sub("^points_", "", names(points_args))
            if (crop) {
                # attach metadata for highlighting selection
                points_args$data <- object@meta.data[rownames(points_data), , drop = FALSE]
                points_args$data <- points_args$data[, colnames(points_args$data)[!is.na(colnames(points_args$data))], drop = FALSE]
                points_args$data <- cbind(points_args$data, points_data)
                points_args$data$x <- points_args$data$x * scale_factor
                points_args$data$y <- points_args$data$y * scale_factor
                points_args$ext <- ext %||% (ext_unscaled * scale_factor)
            } else {
                points_args$data <- Seurat::GetTissueCoordinates(object, image = if(isFALSE(image)) NULL else image)
                points_args$data$.y <- points_args$data$x * scale_factor
                points_args$data$x <- points_args$data$y * scale_factor
                points_args$data$y <- points_args$data$.y
                points_args$data$.y <- NULL
                points_args$data <- cbind(
                    object@meta.data[rownames(points_args$data), , drop = FALSE],
                    points_args$data
                )
            }
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
                points_args$color_by <- colnames(featdata)
                points_args$legend.position <- legend.position
                points_args$legend.direction <- legend.direction
                if (length(features) == 1) {
                    points_args$color_name <- points_args$color_name %||% features
                } else {
                    points_args$color_name <- points_args$color_name %||% "feature"
                    facet_by <- ".facet_var"
                }
            }
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
            points_args$flip_y <- flip_y
            points_args$return_layer <- TRUE
            player <- do.call(SpatialPointsPlot, points_args)
            scales_reused <- intersect(scales_used, attr(player, "scales"))
            if ("fill" %in% scales_reused) {
                players <- c(players, list(ggnewscale::new_scale_fill()))
            }
            players <- c(players, list(player))
            scales_used <- unique(c(scales_used, attr(player, "scales")))
        }
        if (element == "masks" && !is.null(masks)) {
            stop("[SpatialPlot] 'masks' is not supported for Seurat Visium V2 objects.")
        }
        if (element == "shapes" && !is.null(shapes)) {
            stop("[SpatialPlot] 'shapes' is not supported for Seurat Visium V2 objects.")
        }
    }

    if (!is.null(ext_unscaled)) {
        ext <- ext %||% ext_unscaled * scale_factor
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
SpatialPlot.Seurat.SlideSeq <- function(
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
    stopifnot("[SpatialPlot] Either 'group_by' or 'features' should be provided, not both." = is.null(group_by) || is.null(features))
    if (theme %in% c("theme_box", "theme_this", "theme_blank")) {
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
        if (element == "image" && !is.null(image)) {
            stop("[SpatialPlot] 'image' is not supported for Seurat SlideSeq objects.")
        }
        if (element == "points" && !is.null(points)) {
            points_args <- args[startsWith(names(args), "points_")]
            names(points_args) <- sub("^points_", "", names(points_args))
            if (crop) {
                # attach metadata for highlighting selection
                points_args$data <- object@meta.data[rownames(points_data), , drop = FALSE]
                points_args$data <- cbind(points_args$data, points_data)
                points_args$data$x <- points_args$data$x * scale_factor
                points_args$data$y <- points_args$data$y * scale_factor
                points_args$ext <- ext %||% (ext_unscaled * scale_factor)
            } else {
                points_args$data <- Seurat::GetTissueCoordinates(object, image = if(isFALSE(image)) NULL else image)
                points_args$data <- points_args$data[, colnames(points_args$data)[!is.na(colnames(points_args$data))], drop = FALSE]
                points_args$data$.y <- points_args$data$x * scale_factor
                points_args$data$x <- points_args$data$y * scale_factor
                points_args$data$y <- points_args$data$.y
                points_args$data$.y <- NULL
                points_args$data <- cbind(
                    object@meta.data[rownames(points_args$data), , drop = FALSE],
                    points_args$data
                )
            }
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
                points_args$color_by <- colnames(featdata)
                points_args$legend.position <- legend.position
                points_args$legend.direction <- legend.direction
                if (length(features) == 1) {
                    points_args$color_name <- points_args$color_name %||% features
                } else {
                    points_args$color_name <- points_args$color_name %||% "feature"
                    facet_by <- ".facet_var"
                }
            }
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
            points_args$flip_y <- flip_y
            points_args$return_layer <- TRUE
            player <- do.call(SpatialPointsPlot, points_args)
            scales_reused <- intersect(scales_used, attr(player, "scales"))
            if ("fill" %in% scales_reused) {
                players <- c(players, list(ggnewscale::new_scale_fill()))
            }
            players <- c(players, list(player))
            scales_used <- unique(c(scales_used, attr(player, "scales")))
        }
        if (element == "masks" && !is.null(masks)) {
            stop("[SpatialPlot] 'masks' is not supported for Seurat Visium V2 objects.")
        }
        if (element == "shapes" && !is.null(shapes)) {
            stop("[SpatialPlot] 'shapes' is not supported for Seurat Visium V2 objects.")
        }
    }

    if (!is.null(ext_unscaled)) {
        ext <- ext %||% ext_unscaled * scale_factor
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
SpatialFeaturePlot <- function(object, ...) {
    UseMethod("SpatialFeaturePlot", object)
}

#' @export
#' @rdname SpatialPlot
SpatialFeaturePlot.Seurat <- function(object, image = NULL, ...) {
    SpatialPlot(object, image = image, ...)
}

#' Plot categories for spatial data
#'
#' @param object A Seurat object or a Giotto object.
#' @return A ggplot object
#' @export
#' @rdname SpatialPlot
SpatialDimPlot <- function(object, ...) {
    UseMethod("SpatialDimPlot", object)
}

#' @export
#' @rdname SpatialPlot
SpatialDimPlot.Seurat <- function(object, image = NULL, group_by = NULL, ...) {
    if (is.null(group_by)) {
        group_by <- "Identity"
        object@meta.data$Identity <- Seurat::Idents(object)
    }
    SpatialPlot.Seurat(object, image = image, group_by = group_by, ...)
}
