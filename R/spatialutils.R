#' Shared argument definitions for spatial plot functions
#'
#' @description
#' This documentation block defines the full set of parameters shared across all
#' spatial plotting functions in \pkg{scplotter}: \code{\link{SpatFeaturePlot}},
#' \code{\link{SpatDimPlot}}, and the internal \code{SpatPlot()} workhorse. Each
#' function inherits these parameter definitions via \code{@inheritParams}.
#' Individual function documentation pages may override or supplement parameter
#' descriptions when behavior differs between functions.
#'
#' @section Spatial data technologies:
#' \pkg{scplotter} supports several spatial transcriptomics technologies through
#' both Seurat and Giotto objects:
#' \itemize{
#'   \item \strong{Visium (10x Genomics)} — Spot-based spatial transcriptomics.
#'     Coordinates are stored with \code{imagerow}/\code{imagecol} axes. The
#'     tissue image is always available and plotted as the background layer.
#'     Supported via \code{SpatPlot.Seurat.Visium}.
#'   \item \strong{Slide-seq} — Bead-based spatial transcriptomics with
#'     randomized barcode positions. Similar to Visium in coordinate handling
#'     but without a registered tissue image. Supported via
#'     \code{SpatPlot.Seurat.SlideSeq}.
#'   \item \strong{FOV-based (Xenium, CosMx, MERSCOPE, Nanostring)} —
#'     Single-molecule resolution technologies with cell segmentation
#'     boundaries. Supports molecule-level visualization, cell boundary shapes,
#'     and registered images. Supported via \code{SpatPlot.Seurat.FOV} and
#'     \code{SpatPlot.giotto}.
#' }
#'
#' @section Layer system:
#' Spatial plots in \pkg{scplotter} are composed of ordered layers. The
#' \code{layers} argument controls which layers are drawn and in what order:
#' \enumerate{
#'   \item \strong{"image"} — The tissue image or a colored background rectangle.
#'     Must be the first layer when included. Use \code{image = "colorname"} to
#'     fill the background with a solid color instead of an actual image.
#'   \item \strong{"masks"} — Cell segmentation masks (not yet supported for
#'     Seurat or Giotto objects).
#'   \item \strong{"shapes"} — Cell boundary polygons, filled by metadata or
#'     feature expression when \code{shapes_fill_by} is set.
#'   \item \strong{"points"} — Individual points representing cells, spots, or
#'     molecules, colored by \code{group_by}, \code{features}, or
#'     \code{shapes_fill_by}.
#' }
#' Layer-specific styling is controlled via prefixed \code{...} arguments:
#' \code{image_*} for the image layer, \code{masks_*} for masks,
#' \code{shapes_*} for shapes, and \code{points_*} for points.
#'
#' @section Coordinate systems:
#' Different spatial technologies use different coordinate conventions:
#' \itemize{
#'   \item \strong{Visium/Slide-seq} — Coordinates are swapped (x and y are
#'     reversed) relative to the tissue image, and the y-axis is flipped
#'     (\code{flip_y = TRUE}) so the image origin aligns with the plot origin.
#'   \item \strong{FOV-based (Seurat)} — Coordinates are in image space with
#'     no swapping needed. \code{flip_y = FALSE} by default.
#'   \item \strong{FOV-based (Giotto)} — Coordinates are stored as \code{x}
#'     and \code{y} directly. \code{flip_y = FALSE} by default.
#' }
#' These are handled automatically; the \code{x}, \code{y}, and \code{flip_y}
#' arguments are primarily for internal use.
#'
#' @param object A Seurat object (with spatial data loaded via
#'   \pkg{SeuratObject}) or a Giotto object (created with \pkg{Giotto}).
#'   The spatial technology is auto-detected from the object's image class.
#' @param fov The name of the field of view (FOV) to plot. For Seurat FOV-based
#'   objects (Xenium, CosMx, etc.), defaults to
#'   \code{\link[SeuratObject:DefaultFOV]{SeuratObject::DefaultFOV()}}. Not
#'   applicable to Visium or Slide-seq objects.
#' @param boundaries The name of the segmentation boundaries within the FOV
#'   to use for cell outlines. For Seurat FOV-based objects, defaults to
#'   \code{\link[SeuratObject:DefaultBoundary]{SeuratObject::DefaultBoundary()}}.
#'   Not applicable to Visium or Slide-seq objects.
#' @param image Controls the image/background layer. Possible values:
#'   \itemize{
#'     \item \code{NULL} — Default behavior: for Visium, the first image is
#'       used; for Giotto and FOV objects, no image is plotted.
#'     \item A character string naming an image in the object — that specific
#'       image is plotted.
#'     \item A color name (e.g., \code{"white"}, \code{"lightgray"}) — fills
#'       the background with a solid color rectangle.
#'     \item \code{TRUE} — For Visium: uses the first image. For Giotto FOV:
#'       plots all non-overlapping images. For Seurat FOV: raises an error
#'       (no single default image).
#'     \item \code{FALSE} — Disables the image layer entirely.
#'   }
#' @param masks Logical. Whether to plot cell segmentation masks. Currently
#'   not supported — setting this to \code{TRUE} will produce an error for
#'   all object types.
#' @param shapes Controls the cell boundary (shapes) layer. Possible values:
#'   \itemize{
#'     \item \code{TRUE} — For Seurat FOV objects, uses \code{boundaries} as
#'     the shape boundaries. For other object types, uses the default
#'     boundaries. Not supported for Visium or Slide-seq objects.
#'     \item A character string — The name of a specific set of boundaries
#'     within the FOV (Seurat) or spatial info name (Giotto) to use as shapes.
#'     \item \code{FALSE} — Disables the shapes layer.
#'   }
#'   Defaults to \code{TRUE} when \code{shapes_fill_by} is provided, and
#'   \code{NULL} otherwise.
#' @param points Logical. Whether to plot the points layer (cells, spots, or
#'   molecules as points on the spatial coordinates). Default: \code{TRUE}.
#' @param ext The spatial extent (bounding box) of the plot. If \code{NULL},
#'   the extent is calculated automatically from the data or, when
#'   \code{crop = TRUE}, from the tissue coordinates. Can be a numeric vector
#'   in the format \code{c(xmin, xmax, ymin, ymax)} or a
#'   \code{\link[terra:ext]{terra::SpatExtent}} object.
#' @param crop Logical. Whether to crop the plot to the extent of the
#'   tissue/spots. When \code{TRUE} (default), the plot is automatically
#'   zoomed to the data extent with optional \code{padding}. Analogous to the
#'   \code{crop} argument in
#'   \code{\link[Seurat:SpatialDimPlot]{Seurat::SpatialDimPlot()}}.
#' @param group_by A metadata column name used to color the points. Must be a
#'   character or factor column in the object's metadata. For
#'   \code{SpatDimPlot()}, if \code{NULL} and the object has FOV data with
#'   features, defaults to \code{"molecules"}; otherwise defaults to
#'   \code{"Identity"} (the active cluster identities). For
#'   \code{SpatFeaturePlot()}, \code{group_by} is ignored — use
#'   \code{\link{SpatDimPlot}} for categorical grouping. The special value
#'   \code{"molecules"} enables molecule-level visualization in FOV-based
#'   data.
#' @param features A character vector of feature names to visualize. For
#'   \code{SpatFeaturePlot()}, each feature is plotted as a separate facet
#'   (or combined theme when a single feature is given), with expression
#'   values coloring the points. For \code{SpatDimPlot()}, features are
#'   treated as molecule names to plot at single-molecule resolution (FOV
#'   data only). Can include gene names, metadata column names, or dimension
#'   reduction components.
#' @param layer The assay layer from which to extract feature expression
#'   values. For Seurat objects, one of \code{"data"} (default),
#'   \code{"scale.data"}, or \code{"counts"}. For Giotto objects, one of
#'   \code{"normalized"} (default), \code{"scaled"}, \code{"raw"},
#'   \code{"counts"}, or \code{"custom"}.
#' @param scale_factor Internal use. The image scale factor extracted from
#'   the object, used to map between pixel and tissue coordinate spaces.
#'   Automatically determined from the object's image data.
#' @param layers A character vector specifying which layers to include and in
#'   what order. Possible values are \code{"image"}, \code{"masks"},
#'   \code{"shapes"}, and \code{"points"}. Order matters — the first element
#'   is drawn first (bottom). Omit a layer name to disable it. Default:
#'   \code{c("image", "masks", "shapes", "points")} intersected with which
#'   layers are non-\code{FALSE}/non-\code{NULL}.
#' @param flip_y Logical. Whether to flip the y-axis. This is primarily for
#'   internal coordinate system alignment — Visium/Slide-seq objects default
#'   to \code{TRUE}, FOV-based objects to \code{FALSE}. In most cases you do
#'   not need to set this manually.
#' @param padding Numeric. Extra space added around the data extent when
#'   \code{crop = TRUE}. Expressed as a fraction of the data range (e.g.,
#'   \code{0.05} adds 5\% padding on each side). For Seurat FOV objects,
#'   defaults to \code{0}. For other object types, defaults to \code{0} when
#'   an image layer is present and \code{0.05} otherwise.
#' @param image_scale The image scale factor name to use, typically
#'   \code{"lowres"} or \code{"hires"}. Controls which resolution of the
#'   stored image is rendered. Analogous to the \code{image.scale} argument
#'   in \code{\link[Seurat:SpatialDimPlot]{Seurat::SpatialDimPlot()}}.
#' @param x Internal use. The name of the x-coordinate column in the spatial
#'   data. Auto-detected based on the spatial technology (\code{"imagerow"}
#'   for Visium V1, varies for other types).
#' @param y Internal use. The name of the y-coordinate column in the spatial
#'   data. Auto-detected based on the spatial technology (\code{"imagecol"}
#'   for Visium V1, varies for other types).
#' @param nmols Integer. Maximum number of molecules to plot per feature when
#'   \code{group_by = "molecules"}. Analogous to the \code{nmols} argument in
#'   \code{\link[Seurat:ImageDimPlot]{Seurat::ImageDimPlot()}}. Default:
#'   \code{1000}.
#' @param shapes_fill_by A column name in the metadata (or a feature/gene
#'   name) used to fill the cell boundary shapes. If a single color string is
#'   provided, all shapes are filled with that color. When set, \code{shapes}
#'   defaults to \code{TRUE}.
#' @param graph The name of a spatial network graph to overlay on the plot.
#'   Currently supported only for Giotto objects. Possible values:
#'   \itemize{
#'     \item \code{TRUE} — Use the default spatial network.
#'     \item A character string — The graph name. If the name contains
#'       \code{":"}, the part before the colon is used as \code{spat_unit}
#'       and the part after as the graph name.
#'     \item \code{NULL} (default) — No graph overlay.
#'   }
#'   The graph data is retrieved via
#'   \code{\link[GiottoClass:getSpatialNetwork]{GiottoClass::getSpatialNetwork()}}.
#' @param shape Numeric. The point shape (ggplot2 shape aesthetic). Default:
#'   \code{16} (filled circle). See
#'   \url{https://ggplot2.tidyverse.org/reference/aes_linetype_size_shape.html}
#'   for the full shape palette.
#' @param legend.position Character. Legend position. One of
#'   \code{"right"} (default), \code{"left"}, \code{"top"},
#'   \code{"bottom"}, or \code{"none"}.
#' @param legend.direction Character. Legend direction. One of
#'   \code{"vertical"} (default) or \code{"horizontal"}.
#' @param theme A theme function or a character string naming one. Default:
#'   \code{"theme_box"}. Built-in aliases (usable without namespace):
#'   \code{"theme_box"} (\code{\link[plotthis:theme_box]{plotthis::theme_box()}}),
#'   \code{"theme_this"} (\code{\link[plotthis:theme_this]{plotthis::theme_this()}}),
#'   \code{"theme_blank"} (\code{\link[ggplot2:theme_void]{ggplot2::theme_void()}}).
#'   Any \pkg{ggplot2} theme can be used with its fully qualified name (e.g.,
#'   \code{"ggplot2::theme_bw"}).
#' @param theme_args A named list of additional arguments passed to the theme
#'   function. Default: \code{list()}.
#' @param title Character. Plot title. Default: \code{NULL} (no title).
#' @param subtitle Character. Plot subtitle. Default: \code{NULL}.
#' @param xlab Character. x-axis label. Default: \code{NULL}.
#' @param ylab Character. y-axis label. Default: \code{NULL}.
#' @param facet_scales Character. Whether facet scales are \code{"fixed"}
#'   (default), \code{"free"}, \code{"free_x"}, or \code{"free_y"}. Passed
#'   to \code{\link[ggplot2:facet_wrap]{ggplot2::facet_wrap()}} when
#'   multiple features are plotted.
#' @param facet_nrow Integer. Number of facet rows. Default: \code{NULL}
#'   (auto-calculated).
#' @param facet_ncol Integer. Number of facet columns. Default: \code{NULL}
#'   (auto-calculated).
#' @param facet_byrow Logical. Whether to fill facets by row. Default:
#'   \code{TRUE}.
#' @param feat_type Character. The feature type (modality) to query for
#'   expression values in Giotto objects. Common values: \code{"rna"}
#'   (default), \code{"dna"}, \code{"protein"}. Ignored for Seurat objects.
#' @param use_overlap Logical. For Giotto FOV objects, whether to use
#'   pre-computed polygon-feature overlap results (from
#'   \code{GiottoClass::combineFeatureOverlapData()}) instead of cell-level
#'   expression. Default: \code{FALSE}.
#' @param shapes_feat_type Character. The feature type to use when extracting
#'   metadata for shape filling in Giotto objects. Default: \code{"cell"}.
#' @param shapes_alpha Numeric. Transparency (alpha) value for the shapes
#'   layer, between 0 and 1. When points are also plotted, defaults to
#'   \code{0.5} so points remain visible on top of shapes; otherwise defaults
#'   to \code{1}.
#' @param spat_unit Character. The spatial unit to query in a Giotto object
#'   (e.g., \code{"cell"}, \code{"subcellular"}). Auto-detected if
#'   \code{NULL}. Ignored for Seurat objects.
#' @param spat_loc_name Character. The spatial locations name to query in a
#'   Giotto object. Auto-detected from available spatial locations if
#'   \code{NULL}. Ignored for Seurat objects.
#' @param spat_enr_names Character. Spatial enrichment results names in a
#'   Giotto object (for enrichment-based feature extraction). Ignored for
#'   Seurat objects.
#' @param ... Additional arguments passed to the underlying layer functions.
#'   Arguments are dispatched by prefix:
#'   \describe{
#'     \item{\code{image_*}}{Arguments passed to
#'       \code{\link[plotthis:SpatImagePlot]{plotthis::SpatImagePlot()}}
#'       (e.g., \code{image_alpha}, \code{image_interpolation}).}
#'     \item{\code{masks_*}}{Arguments passed to
#'       \code{\link[plotthis:SpatMasksPlot]{plotthis::SpatMasksPlot()}}.}
#'     \item{\code{shapes_*}}{Arguments passed to
#'       \code{\link[plotthis:SpatShapesPlot]{plotthis::SpatShapesPlot()}}
#'       (e.g., \code{shapes_color}, \code{shapes_linewidth}).}
#'     \item{\code{points_*}}{Arguments passed to
#'       \code{\link[plotthis:SpatPointsPlot]{plotthis::SpatPointsPlot()}}
#'       (e.g., \code{points_size}, \code{points_alpha}).}
#'     \item{No prefix}{Arguments without a recognized prefix are treated as
#'       points arguments, but with lower priority than \code{points_*}
#'       arguments.}
#'   }
#' @keywords internal
#' @name spatialplot_args
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
