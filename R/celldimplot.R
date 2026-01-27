#' Cell Dimension Reduction Plot
#'
#' @description This function creates a dimension reduction plot for a Seurat object
#' a Giotto object, a path to an .h5ad file or an opened `H5File` by `hdf5r` package.
#' It allows for various customizations such as grouping by metadata,
#' adding edges between cell neighbors, highlighting specific cells, and more.
#' This function is a wrapper around [plotthis::DimPlot()], which provides a
#' flexible way to visualize cell clusters in reduced dimensions. This function
#' extracts the necessary data from the Seurat or Giotto object and passes it to
#' [plotthis::DimPlot()].
#' @param object A seurat object, a giotto object, a path to an .h5ad file or an opened `H5File` by `hdf5r` package.
#' @param reduction Name of the reduction to plot (for example, "umap").
#' @param graph Specify the graph name to add edges between cell neighbors to the plot.
#' @param velocity The name of velocity reduction to plot cell velocities.
#' It is typically `"stochastic_<reduction>"`, `"deterministic_<reduction>"`, or `"dynamical_<reduction>"`.
#' @param group_by A character vector of column name(s) to group the data. Default is NULL.
#' @param ident Alias for `group_by`.
#' @param spat_unit The spatial unit to use for the plot. Only applied to Giotto objects.
#' @param feat_type feature type of the features (e.g. "rna", "dna", "protein"), only applied to Giotto objects.
#' @param ... Other arguments passed to [plotthis::DimPlot()].
#' @return A ggplot object or a list if `combine` is FALSE
#' @export
#' @seealso [scplotter::CellStatPlot()] [scplotter::CellVelocityPlot()]
#' @importFrom SeuratObject DefaultDimReduc Embeddings Graphs Reductions Idents
#' @importFrom plotthis DimPlot
#' @details
#' See
#' * <https://pwwang.github.io/scplotter/articles/Giotto_CODEX.html>
#' * <https://pwwang.github.io/scplotter/articles/Giotto_seqFISH.html>
#' * <https://pwwang.github.io/scplotter/articles/Giotto_SlideSeq.html>
#' * <https://pwwang.github.io/scplotter/articles/Giotto_Spatial_CITE-Seq.html>
#' * <https://pwwang.github.io/scplotter/articles/Giotto_Visium.html>
#' * <https://pwwang.github.io/scplotter/articles/Giotto_VisiumHD.html>
#' * <https://pwwang.github.io/scplotter/articles/Giotto_Xenium.html>
#'
#' for examples of using this function with a Giotto object.
#'
#' And see:
#' * <https://pwwang.github.io/scplotter/articles/Working_with_anndata_h5ad_files.html>
#'
#' for examples of using this function with .h5ad files.
#' @examples
#' \donttest{
#' data(pancreas_sub)
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
#'   highlight = 'SubCellType == "Epsilon"'
#' )
#' CellDimPlot(pancreas_sub,
#'   group_by = "SubCellType", split_by = "Phase", reduction = "UMAP",
#'   highlight = TRUE, theme = "theme_blank", legend.position = "none"
#' )
#' CellDimPlot(pancreas_sub,
#'   group_by = "SubCellType", facet_by = "Phase", reduction = "UMAP",
#'   highlight = TRUE, theme = "theme_blank", legend.position = "none"
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
#'             add_mark = TRUE, mark_expand = grid::unit(1, "mm"))
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
#'   highlight = TRUE
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
#'   stat_plot_type = "bar", stat_type = "count")
#' CellDimPlot(pancreas_sub,
#'   group_by = "CellType", reduction = "UMAP", stat_by = "Phase",
#'   stat_plot_type = "line", stat_type = "count", stat_args = list(point_size = 1))
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
#'
#' # Velocity
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "PCA",
#'   velocity = "stochastic_PCA")
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "PCA",
#'   velocity = "stochastic_PCA", velocity_plot_type = "grid", pt_alpha = 0.5)
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "PCA",
#'   velocity = "stochastic_PCA", velocity_plot_type = "stream", pt_alpha = 0.5)
#' }
CellDimPlot <- function(
    object, reduction = NULL, graph = NULL, group_by = NULL, ident = NULL,
    spat_unit = NULL, feat_type = NULL, velocity = NULL, ...
) {
    UseMethod("CellDimPlot")
}

#' @export
CellDimPlot.giotto <- function(
    object, reduction = NULL, graph = NULL, group_by = NULL, ident = NULL,
    spat_unit = NULL, feat_type = NULL, velocity = NULL, ...
) {
    stopifnot("[CellDimPlot] Either 'group_by' or 'ident' must be provided, not both." =
        is.null(group_by) || is.null(ident) || !identical(group_by, ident))
    group_by <- group_by %||% ident
    stopifnot("[CellDimPlot] 'group_by' is required for giotto objects." = !is.null(group_by))

    spat_unit <- GiottoClass::set_default_spat_unit(
        gobject = object,
        spat_unit = spat_unit
    )

    feat_type <- GiottoClass::set_default_feat_type(
        gobject = object,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    if (!is.null(graph)) {
        # spat_unit	feat_type	nn_type	name
        # cell	rna	sNN	sNN.pca
        graph_info <- GiottoClass::list_nearest_networks(
            gobject = object,
            spat_unit = spat_unit,
            feat_type = feat_type
        )
        if (is.null(graph_info)) {
            stop(
                "[CellDimPlot] The object does not have any nearest networks with given ",
                "spat_unit (", spat_unit, ") and feat_type (", feat_type, ")."
            )
        }
        if (isTRUE(graph)) { graph <- graph_info$name[1] }
        graph_info <- graph_info[, graph_info$name == graph, drop = FALSE]
        if (nrow(graph_info) == 0) {
            stop(
                "[CellDimPlot] The object does not have network: ", graph,
                " with given spat_unit (", spat_unit, ") and feat_type (", feat_type, ")."
            )
        }
        graph_info <- graph_info[1, , drop = FALSE]
        graph <- igraph::as_adjacency_matrix(GiottoClass::getNearestNetwork(
            object,
            spat_unit = spat_unit,
            feat_type = feat_type,
            nn_type = graph_info$nn_type,
            name = graph_info$name,
            output = "igraph"
        ))
    }

    data <- cbind(
        GiottoClass::getDimReduction(
            gobject = object,
            reduction = "cells",
            reduction_method = reduction,
            spat_unit = spat_unit,
            feat_type = feat_type,
            output = "matrix"
        ),
        GiottoClass::getCellMetadata(
            gobject = object,
            spat_unit = spat_unit,
            feat_type = feat_type,
            output = "data.table"
        )
    )

    if (!is.null(velocity)) {
        velocity <- GiottoClass::getDimReduction(
            gobject = object,
            reduction = "cells",
            reduction_method = velocity,
            spat_unit = spat_unit,
            feat_type = feat_type,
            output = "matrix"
        )[, 1:2, drop = FALSE]
    }

    rownames(data) <- data$cell_ID

    DimPlot(data, graph = graph, group_by = group_by, velocity = velocity, ...)
}

#' @export
CellDimPlot.Seurat <- function(
    object, reduction = NULL, graph = NULL, group_by = NULL, ident = NULL,
    spat_unit = NULL, feat_type = NULL, velocity = NULL,
    ...
) {
    stopifnot("[CellDimPlot] 'spat_unit' and 'feat_type' are not used for Seurat objects." =
        is.null(spat_unit) && is.null(feat_type))
    stopifnot("[CellDimPlot] Either 'group_by' or 'ident' must be provided, not both." =
        is.null(group_by) || is.null(ident) || !identical(group_by, ident))
    group_by <- group_by %||% ident

    reduction <- reduction %||% DefaultDimReduc(object)

    if (!reduction %in% Reductions(object)) {
        stop("The object does not have reduction:", reduction)
    }

    if (!is.null(graph)) {
        if (!graph %in% Graphs(object)) {
            stop("The object does not have graph:", graph)
        }
        graph <- object@graphs[[graph]]
    }

    emb <- Embeddings(object, reduction = reduction)
    data <- cbind(emb, object@meta.data[rownames(emb), , drop = FALSE])
    if (is.null(group_by)) {
        group_by <- "Identity"
        data[[group_by]] <- Idents(object)
    }

    if (!is.null(velocity)) {
        velocity <- Embeddings(object, reduction = velocity)[, 1:2, drop = FALSE]
    }

    DimPlot(data, graph = graph, group_by = group_by, velocity = velocity, ...)
}

#' @export
CellDimPlot.character <- function(
    object, reduction = NULL, graph = NULL, group_by = NULL, ident = NULL,
    spat_unit = NULL, feat_type = NULL, velocity = NULL,
    ...
) {
    stopifnot("[CellDimPlot] Currently only supports .h5ad files when called with a string/path." =
        endsWith(object, ".h5ad"))
    stopifnot("[CellDimPlot] Either 'group_by' or 'ident' must be provided, not both." =
        is.null(group_by) || is.null(ident) || !identical(group_by, ident))
    group_by <- group_by %||% ident

    object <- hdf5r::H5File$new(object, mode = "r")
    on.exit(object$close_all())

    CellDimPlot.H5File(
        object = object, reduction = reduction, graph = graph,
        group_by = group_by, spat_unit = spat_unit, feat_type = feat_type,
        velocity = velocity, ...
    )
}

#' @export
CellDimPlot.H5File <- function(
    object, reduction = NULL, graph = NULL, group_by = NULL, ident = NULL,
    spat_unit = NULL, feat_type = NULL, velocity = NULL,
    ...
) {
    stopifnot("[CellDimPlot] 'spat_unit' and 'feat_type' are not used for anndata." =
        is.null(spat_unit) && is.null(feat_type))
    stopifnot("[CellDimPlot] Either 'group_by' or 'ident' must be provided, not both." =
        is.null(group_by) || is.null(ident) || !identical(group_by, ident))
    group_by <- group_by %||% ident
    stopifnot("[CellDimPlot] 'group_by' is required for anndata (h5ad) objects." = !is.null(group_by))

    reductions <- names(object[['obsm']])
    if (is.null(reductions) || length(reductions) == 0) {
        stop("[CellDimPlot] The object does not have any reductions.")
    }

    reduction <- reduction %||% reductions[1]
    if (!startsWith(reduction, "X_") && !reduction %in% reductions) {
        reduction <- paste0("X_", reduction)
    }

    if (!reduction %in% reductions) {
        stop("[CellDimPlot] The object does not have reduction:", reduction)
    }

    if (!is.null(graph)) {
        ggrp <- object[['obsp']][['connectivities']]
        if (is.null(ggrp)) {
            stop("[CellDimPlot] The object does not have any graph/connectivities.")
        }

        graph <- h5group_to_matrix(ggrp)
        graph <- igraph::graph_from_adjacency_matrix(graph, mode = "undirected", weighted = TRUE)
        graph <- igraph::as_adjacency_matrix(graph, attr = "weight", sparse = TRUE)
        colnames(graph) <- rownames(graph) <- object[['obs']][['index']]$read()
    }

    data <- cbind(t(object[["obsm"]][[reduction]]$read()), h5group_to_dataframe(object[['obs']]))
    rownames(data) <- object[['obs']][['index']]$read()

    if (!is.null(velocity)) {
        velocity <- t(object[["obsm"]][[velocity]]$read())
    }

    DimPlot(data, graph = graph, group_by = group_by, velocity = velocity, ...)
}


#' Cell Velocity Plot
#'
#' @description This function creates a cell velocity plot for a Seurat object,
#' a Giotto object, a path to an .h5ad file or an opened `H5File` by `hdf5r` package.
#' It allows for various customizations such as grouping by metadata,
#' adding edges between cell neighbors, highlighting specific cells, and more.
#' This function is a wrapper around [plotthis::VelocityPlot()], which provides a
#' flexible way to visualize cell velocities in reduced dimensions. This function
#' extracts the cell embeddings and velocity embeddings from the Seurat or Giotto object
#' and passes them to [plotthis::VelocityPlot()].
#' @param object A seurat object, a giotto object, a path to an .h5ad file or an opened `H5File` by `hdf5r` package.
#' @param reduction Name of the reduction to plot (for example, "umap").
#' @param v_reduction Name of the velocity reduction to plot (for example, "stochastic_umap").
#' It should be the same as the reduction used to calculate the velocity.
#' @param spat_unit The spatial unit to use for the plot. Only applied to Giotto objects.
#' @param feat_type feature type of the features (e.g. "rna", "dna", "protein"), only applied to Giotto objects.
#' @param group_by A character vector of metadata column name(s) to group (color) the data. Default is NULL.
#' @param ident Alias for `group_by`.
#' @param ... Other arguments passed to [plotthis::VelocityPlot()].
#' @return A ggplot object
#' @export
#' @seealso [scplotter::CellDimPlot()]
#' @importFrom SeuratObject DefaultDimReduc Embeddings Reductions
#' @importFrom plotthis VelocityPlot
#' @details See:
#' * <https://pwwang.github.io/scplotter/articles/Working_with_anndata_h5ad_files.html>
#'
#' for examples of using this function with .h5ad files.
#' @examples
#' \donttest{
#' data(pancreas_sub)
#'
#' CellVelocityPlot(pancreas_sub, reduction = "PCA", v_reduction = "stochastic_PCA")
#' CellVelocityPlot(pancreas_sub, reduction = "PCA", v_reduction = "stochastic_PCA",
#'  plot_type = "grid")
#' CellVelocityPlot(pancreas_sub, reduction = "PCA", v_reduction = "stochastic_PCA",
#'  plot_type = "stream")
#' CellVelocityPlot(pancreas_sub, reduction = "PCA", v_reduction = "stochastic_PCA",
#'  group_by = "SubCellType")
#' }
CellVelocityPlot <- function(
    object, reduction, v_reduction, spat_unit = NULL, feat_type = NULL, group_by = NULL, ident = NULL, ...
) {
    UseMethod("CellVelocityPlot")
}

#' @export
CellVelocityPlot.giotto <- function(
    object, reduction, v_reduction, spat_unit = NULL, feat_type = NULL, group_by = NULL, ident = NULL, ...
) {
    stopifnot("[CellVelocityPlot] Either 'group_by' or 'ident' must be provided, not both." =
        is.null(group_by) || is.null(ident) || !identical(group_by, ident))
    group_by <- group_by %||% ident

    spat_unit <- GiottoClass::set_default_spat_unit(
        gobject = object,
        spat_unit = spat_unit
    )

    feat_type <- GiottoClass::set_default_feat_type(
        gobject = object,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    reduc_info <- GiottoClass::list_dim_reductions(
        gobject = object,
        data_type = "cells",
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    if (!reduction %in% reduc_info$name) {
        stop("[CellVelocityPlot] The object does not have reduction:", reduction)
    }
    if (!v_reduction %in% reduc_info$name) {
        stop("[CellVelocityPlot] The object does not have velocity reduction:", v_reduction)
    }

    if (!is.null(group_by)) {
        metadata <- GiottoClass::getCellMetadata(
            gobject = object,
            spat_unit = spat_unit,
            feat_type = feat_type,
            output = "data.table"
        )
        if (!group_by %in% colnames(metadata)) {
            stop("[CellVelocityPlot] The object does not have metadata column:", group_by)
        }
        group_by <- metadata[[group_by]]
    }

    VelocityPlot(
        embedding = GiottoClass::getDimReduction(
            gobject = object,
            reduction = "cells",
            reduction_method = reduction,
            spat_unit = spat_unit,
            feat_type = feat_type,
            output = "matrix"
        )[, 1:2, drop = FALSE],
        v_embedding = GiottoClass::getDimReduction(
            gobject = object,
            reduction = "cells",
            reduction_method = v_reduction,
            spat_unit = spat_unit,
            feat_type = feat_type,
            output = "matrix"
        )[, 1:2, drop = FALSE],
        group_by = group_by,
        ...
    )
}

#' @export
CellVelocityPlot.Seurat <- function(
    object, reduction, v_reduction, spat_unit = NULL, feat_type = NULL, group_by = NULL, ident = NULL, ...
) {
    stopifnot("[CellVelocityPlot] 'spat_unit' and 'feat_type' are not used for Seurat objects." =
        is.null(spat_unit) && is.null(feat_type))
    stopifnot("[CellVelocityPlot] Either 'group_by' or 'ident' must be provided, not both." =
        is.null(group_by) || is.null(ident) || !identical(group_by, ident))
    group_by <- group_by %||% ident

    if (!reduction %in% Reductions(object)) {
        stop("[CellVelocityPlot] The object does not have reduction:", reduction)
    }
    if (!v_reduction %in% Reductions(object)) {
        stop("[CellVelocityPlot] The object does not have velocity reduction:", v_reduction)
    }

    if (!is.null(group_by) && !group_by %in% colnames(object@meta.data)) {
        stop("[CellVelocityPlot] The object does not have metadata column:", group_by)
    } else if (!is.null(group_by)) {
        group_by <- object@meta.data[[group_by]]
    }

    VelocityPlot(
        embedding = Embeddings(object, reduction = reduction)[, 1:2, drop = FALSE],
        v_embedding = Embeddings(object, reduction = v_reduction)[, 1:2, drop = FALSE],
        group_by = group_by,
        ...
    )
}


#' @export
CellVelocityPlot.character <- function(
    object, reduction, v_reduction, spat_unit = NULL, feat_type = NULL, group_by = NULL, ident = NULL, ...
) {
    stopifnot("[CellVelocityPlot] Currently only supports .h5ad files when called with a string/path." =
        endsWith(object, ".h5ad"))
    stopifnot("[CellVelocityPlot] Either 'group_by' or 'ident' must be provided, not both." =
        is.null(group_by) || is.null(ident) || !identical(group_by, ident))
    group_by <- group_by %||% ident

    object <- hdf5r::H5File$new(object, mode = "r")
    on.exit(object$close_all())

    CellVelocityPlot.H5File(
        object = object, reduction = reduction, v_reduction = v_reduction,
        spat_unit = spat_unit, feat_type = feat_type, group_by = group_by, ...
    )
}

#' @export
CellVelocityPlot.H5File <- function(
    object, reduction, v_reduction, spat_unit = NULL, feat_type = NULL, group_by = NULL, ident = NULL, ...
) {
    stopifnot("[CellVelocityPlot] 'spat_unit' and 'feat_type' are not used for anndata." =
        is.null(spat_unit) && is.null(feat_type))
    stopifnot("[CellVelocityPlot] Either 'group_by' or 'ident' must be provided, not both." =
        is.null(group_by) || is.null(ident) || !identical(group_by, ident))
    group_by <- group_by %||% ident

    reduc_info <- names(object[['obsm']])
    if (is.null(reduc_info) || length(reduc_info) == 0) {
        stop("[CellVelocityPlot] The object does not have any reductions.")
    }

    if (!startsWith(reduction, "X_") && !reduction %in% reduc_info) {
        reduction <- paste0("X_", reduction)
    }

    if (!reduction %in% reduc_info) {
        stop("[CellVelocityPlot] The object does not have reduction:", reduction)
    }
    if (!v_reduction %in% reduc_info) {
        stop("[CellVelocityPlot] The object does not have velocity reduction:", v_reduction)
    }

    if (!is.null(group_by)) {
        metakeys <- names(object[['obs']])
        if (!group_by %in% metakeys) {
            stop("[CellVelocityPlot] The object does not have metadata column:", group_by)
        }
        if (inherits(object[['obs']][[group_by]], "H5Group") &&
            setequal(names(object[['obs']][[group_by]]), c("categories", "codes"))) {
            codes <- object[['obs']][[group_by]][['codes']]$read()
            categories <- object[['obs']][[group_by]][['categories']]$read()
            group_by <- factor(categories[codes + 1], levels = categories)
        } else if (inherits(object[['obs']][[group_by]], "H5Group")) {
            stop("[CellVelocityPlot] Complex data structure detected in metadata (obs) column:", group_by)
        } else {
            # Read the metadata column directly
            group_by <- object[['obs']][[group_by]]$read()
        }
    }

    VelocityPlot(
        embedding = t(object[["obsm"]][[reduction]]$read())[, 1:2, drop = FALSE],
        v_embedding = t(object[["obsm"]][[v_reduction]]$read())[, 1:2, drop = FALSE],
        group_by = group_by,
        ...
    )
}
