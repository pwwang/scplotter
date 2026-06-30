#' Cell Dimension Reduction Plot
#'
#' @description
#' Visualizes single-cell data in reduced dimension space (e.g., UMAP, t-SNE,
#' PCA). This is the primary function for exploring cell clustering, cell
#' identity, and spatial relationships in transcriptomics datasets. It creates
#' scatter plots where each point represents a cell, positioned by its
#' coordinates in the reduced dimension space and colored by metadata variables
#' such as cell type, sample condition, or cluster assignment.
#'
#' `CellDimPlot` serves as a unified interface across multiple single-cell data
#' containers:
#' \itemize{
#'   \item **Seurat objects** — Extracts embeddings from `Reductions()` and
#'         metadata from `@meta.data`. The default reduction is auto-detected
#'         via `default_dimreduc()`.
#'   \item **Giotto objects** — Extracts spatial dimension reductions and cell
#'         metadata using `spat_unit` and `feat_type` to identify the correct
#'         spatial unit and feature type.
#'   \item **h5ad files** (.h5ad or opened `H5File`) — Reads from `obsm` for
#'         embeddings and `obs` for metadata. Reduction names are automatically
#'         prefixed with `"X_"` when needed (e.g., `"umap"` → `"X_umap"`).
#' }
#'
#' Beyond basic cluster visualization, `CellDimPlot` supports a rich set of
#' visual overlays and analytical enhancements:
#' \itemize{
#'   \item **Cluster highlighting** — Emphasize cells matching a logical
#'         expression while dimming others (`highlight`).
#'   \item **Group labels** — Add text labels at group centroids (`label`,
#'         `label_insitu`).
#'   \item **Group marks** — Draw boundary shapes around groups: ellipse,
#'         rectangle, or circle (`add_mark`, `mark_type`).
#'   \item **Density contours** — Overlay 2D density estimates (`add_density`).
#'   \item **Neighbor graphs** — Draw edges between neighboring cells from
#'         k-NN or shared-nearest-neighbor graphs (`graph`).
#'   \item **Lineage trajectories** — Overlay pseudotime lineage curves
#'         (`lineages`).
#'   \item **Velocity arrows** — Overlay RNA velocity vectors on the embedding
#'         (`velocity`). For dedicated velocity visualization with grid or
#'         stream plots, see \code{\link{CellVelocityPlot}}.
#'   \item **Statistical charts** — Embed small bar, ring, or line charts at
#'         group positions showing composition of a second variable (`stat_by`,
#'         `stat_plot_type`).
#'   \item **Hexagonal binning** — Replace scatter points with binned hexagons
#'         for large datasets (`hex`).
#'   \item **3D visualization** — Plot three dimensions by specifying
#'         `dims = 1:3`.
#'   \item **Rasterization** — Render points as a raster image for performance
#'         with large cell counts (`raster`).
#' }
#'
#' @section Data extraction:
#' The function extracts cell embeddings and metadata from the input object,
#' combines them into a single data frame, and passes the result to
#' \code{\link[plotthis:DimPlot]{plotthis::DimPlot()}} along with any additional
#' arguments provided via `...`. The extraction logic varies by object type:
#'
#' **Seurat objects:** `Embeddings(object, reduction)` provides the coordinates;
#' `object@meta.data` provides the metadata. When `group_by` is `NULL`,
#' `Idents(object)` is used as the default grouping variable (column name
#' `"Identity"`). The `graph` parameter references `object@graphs`.
#'
#' **Giotto objects:** `getDimReduction()` with `output = "matrix"` provides
#' coordinates; `getCellMetadata()` with `output = "data.table"` provides
#' metadata. Both `spat_unit` and `feat_type` are resolved to defaults if not
#' specified. The `graph` parameter references nearest-neighbor networks
#' retrieved via `getNearestNetwork()` and converted with \pkg{igraph}.
#'
#' **h5ad files:** `obsm[[reduction]]` provides the embedding matrix;
#' `obs` provides metadata (decoded via `h5group_to_dataframe()` for
#' categorical variables). The `graph` parameter references
#' `obsp[["connectivities"]]`, converted via \pkg{igraph} to an adjacency
#' matrix.
#'
#' @section Velocity overlay:
#' When `velocity` is specified, velocity vectors are overlaid on the dimension
#' reduction plot. The `velocity` parameter names a second reduction
#' (e.g., `"stochastic_UMAP"`) whose first two dimensions encode the velocity
#' arrows. This is a lightweight overlay — for dedicated velocity visualization
#' with grid or stream plots, use \code{\link{CellVelocityPlot}} which delegates to
#' \code{\link[plotthis:VelocityPlot]{plotthis::VelocityPlot()}}.
#'
#' @param object A Seurat object, a Giotto object, a path to an `.h5ad` file,
#'   or an opened `H5File` from the \pkg{hdf5r} package.
#' @param reduction Name of the dimension reduction to plot. Typical values are
#'   `"umap"`, `"tsne"`, or `"pca"`. For Seurat objects, the default reduction
#'   is auto-detected via `default_dimreduc()`. For h5ad files, the prefix
#'   `"X_"` is added automatically when needed (i.e., `"umap"` is treated as
#'   `"X_umap"` if `"umap"` alone is not found). For Giotto objects, this
#'   parameter is required unless a default reduction has been set.
#' @param graph Name of the neighbor graph used to draw edges between
#'   neighboring cells on the plot. For Seurat objects, this is a graph name
#'   in `Graphs(object)` (e.g., `"RNA_nn"`, `"RNA_snn"`). For Giotto objects,
#'   this is a nearest-neighbor network name; set `graph = TRUE` to use the
#'   first available network. For h5ad files, connectivities are read from
#'   `obsp[["connectivities"]]`. Setting `graph` to `NULL` (the default)
#'   suppresses edge drawing.
#' @param velocity Name of a velocity reduction whose first two dimensions
#'   encode RNA velocity vectors to overlay on the plot. Typical values are
#'   `"stochastic_<reduction>"`, `"deterministic_<reduction>"`, or
#'   `"dynamical_<reduction>"` (e.g., `"stochastic_UMAP"`). For dedicated
#'   velocity visualization (grid or stream plots), use \code{\link{CellVelocityPlot}}
#'   instead. Default is `NULL` (no velocity overlay).
#' @param group_by Character vector of metadata column name(s) used to color
#'   the cells. Can be a single column (e.g., `"CellType"`) or multiple
#'   columns for combined grouping. Default is `NULL`, which falls back to
#'   the active identity for Seurat objects and is required for other object
#'   types.
#' @param ident Alias for `group_by`, provided for compatibility with Seurat's
#'   naming convention. When both `group_by` and `ident` are specified, they
#'   must be identical.
#' @param spat_unit Spatial unit name for Giotto objects (e.g., `"cell"`).
#'   Ignored for Seurat and h5ad inputs. If `NULL`, the default spatial unit
#'   is auto-detected via `GiottoClass::set_default_spat_unit()`.
#' @param feat_type Feature type name for Giotto objects (e.g., `"rna"`,
#'   `"dna"`, `"protein"`). Ignored for Seurat and h5ad inputs. If `NULL`, the
#'   default feature type is auto-detected via
#'   `GiottoClass::set_default_feat_type()`.
#' @param ... Additional arguments passed to
#'   \code{\link[plotthis:DimPlot]{plotthis::DimPlot()}}. Key parameters
#'   include:
#'   \itemize{
#'     \item **Appearance:** `pt_size`, `pt_alpha`, `palette`, `palcolor`,
#'           `theme`, `theme_args`, `legend.position`, `legend.direction`
#'     \item **Highlighting:** `highlight` — a logical expression as a string
#'           (e.g., `'CellType == "Beta"'`) to emphasize matching cells
#'     \item **Labels:** `label`, `label_insitu`, `label_repel`,
#'           `label_fg`, `label_bg`, `label_size`, `label_segment_color`
#'     \item **Group marks:** `add_mark`, `mark_type` (`"ellipse"`,
#'           `"rect"`, `"circle"`), `mark_expand`, `mark_alpha`,
#'           `mark_linetype`
#'     \item **Density:** `add_density`, `density_filled`,
#'           `density_filled_palette`
#'     \item **Lineages:** `lineages`, `lineages_whiskers`, `lineages_span`
#'     \item **Statistical charts:** `stat_by`, `stat_plot_type`
#'           (`"ring"`, `"bar"`, `"line"`), `stat_type` (`"percent"`,
#'           `"count"`), `stat_plot_label`, `stat_plot_size`, `stat_args`
#'     \item **Hexagonal binning:** `hex`, `hex_bins`, `hex_count`
#'     \item **Rendering:** `raster`, `raster_dpi`
#'     \item **Dimensionality:** `dims` — which dimensions to plot
#'           (default `1:2`; set to `1:3` for 3D)
#'     \item **Layout:** `split_by`, `facet_by`, `combine`, `nrow`, `ncol`
#'   }
#' @return A `ggplot` object, or a list of `ggplot` objects if `combine = FALSE`
#'   is passed via `...`.
#' @note
#' **Default reduction for Seurat objects:** When `reduction = NULL`,
#' `CellDimPlot` calls `default_dimreduc(object)` to determine the default
#' reduction. For Seurat >= 5.4.0, this uses the official
#' `DefaultDimReduc` getter; for older versions, it falls back to
#' `object@misc$DefaultDimReduc`.
#'
#' **Giotto spatial units and feature types:** The `spat_unit` and `feat_type`
#' parameters are required to locate the correct data within Giotto's
#' hierarchical spatial data structure. When `NULL`, Giotto's own default
#' resolution is used, which is typically sufficient for standard analyses.
#'
#' **Performance with large datasets:** For datasets with many cells (e.g.,
#' over 50,000), consider using `raster = TRUE` to render points as a raster
#' image, or `hex = TRUE` to use hexagonal binning. Both options significantly
#' reduce rendering time and output file size.
#'
#' **Factor ordering:** The order of groups in the legend follows the factor
#' levels of the `group_by` column. Set factor levels on your metadata column
#' before plotting to control legend order.
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[plotthis:DimPlot]{plotthis::DimPlot()}} — The
#'         underlying plotting engine
#'   \item \code{\link{CellVelocityPlot}} — Dedicated RNA velocity visualization
#'   \item \code{\link{CellStatPlot}} — Statistical summaries and comparisons
#' }
#'
#' @export
#' @importFrom SeuratObject Embeddings Graphs Reductions Idents
#' @importFrom plotthis DimPlot
#' @details
#' See the following vignettes for examples with Giotto objects:
#' * <https://pwwang.github.io/scplotter/articles/Giotto_CODEX.html>
#' * <https://pwwang.github.io/scplotter/articles/Giotto_seqFISH.html>
#' * <https://pwwang.github.io/scplotter/articles/Giotto_SlideSeq.html>
#' * <https://pwwang.github.io/scplotter/articles/Giotto_Spatial_CITE-Seq.html>
#' * <https://pwwang.github.io/scplotter/articles/Giotto_Visium.html>
#' * <https://pwwang.github.io/scplotter/articles/Giotto_VisiumHD.html>
#' * <https://pwwang.github.io/scplotter/articles/Giotto_Xenium.html>
#'
#' And for examples with h5ad files:
#' * <https://pwwang.github.io/scplotter/articles/Working_with_anndata_h5ad_files.html>
#' @examples
#' \donttest{
#' set.seed(8525)
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
#'
#' # 3D plot
#' pancreas_sub@reductions$UMAP@cell.embeddings <- cbind(
#'  pancreas_sub@reductions$UMAP@cell.embeddings,
#'  umap_3 = rnorm(1000)  # fake the 3rd dimension
#' )
#' CellDimPlot(pancreas_sub, group_by = "SubCellType", reduction = "UMAP",
#'             dims = 1:3, label = TRUE)
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
        is.null(group_by) || is.null(ident) || identical(group_by, ident))
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
        is.null(group_by) || is.null(ident) || identical(group_by, ident))
    group_by <- group_by %||% ident

    reduction <- reduction %||% default_dimreduc(object)

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
    dots <- list(...)
    dims <- dots$dims %||% 1:2
    # Check if dims exist
    if (max(dims) > ncol(emb)) {
        stop("[CellDimPlot] The reduction does not have enough dimensions for the specified 'dims' [",
            paste(dims, collapse = ", "), "].")
    }

    data <- cbind(emb, object@meta.data[rownames(emb), setdiff(colnames(object@meta.data), colnames(emb)), drop = FALSE])
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
        is.null(group_by) || is.null(ident) || identical(group_by, ident))
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
        is.null(group_by) || is.null(ident) || identical(group_by, ident))
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
#' @description
#' Visualizes RNA velocity on a reduced-dimension embedding. RNA velocity
#' infers the future transcriptional state of individual cells by modeling the
#' ratio of unspliced (nascent) to spliced (mature) mRNA transcripts. On a
#' dimension reduction plot, velocity is displayed as arrows (or grid/stream
#' fields) showing the predicted direction and magnitude of transcriptional
#' change for each cell — effectively revealing the "flow" of cells through
#' differentiation, development, or other state transitions.
#'
#' `CellVelocityPlot` serves as a unified interface across multiple single-cell
#' data containers:
#' \itemize{
#'   \item **Seurat objects** — Extracts embeddings from both the main
#'         reduction (`reduction`) and the velocity reduction
#'         (`v_reduction`) via `Embeddings()`; metadata for grouping via
#'         `@meta.data`.
#'   \item **Giotto objects** — Extracts dimension reductions via
#'         `getDimReduction()` using `spat_unit` and `feat_type` to identify
#'         the correct spatial unit and feature type.
#'   \item **h5ad files** (.h5ad or opened `H5File`) — Reads from `obsm` for
#'         both the main and velocity embeddings; `obs` for metadata. Reduction
#'         names are automatically prefixed with `"X_"` when needed.
#' }
#'
#' @section RNA velocity background:
#' RNA velocity requires two dimension reductions stored in the object:
#' \enumerate{
#'   \item **Main reduction** (`reduction`): A standard dimension reduction
#'         (e.g., UMAP, PCA, t-SNE) computed on the cells' current expression
#'         profiles. This defines the positions of cells in the plot.
#'   \item **Velocity reduction** (`v_reduction`): A dimension reduction
#'         computed on velocity vectors derived from spliced/unspliced RNA
#'         ratios. The first two dimensions of this reduction encode the
#'         direction and magnitude of predicted transcriptional change.
#' }
#' Velocity reductions typically follow the naming convention
#' `"<model>_<reduction>"`, where `<model>` is one of:
#' \itemize{
#'   \item `"stochastic"` — Stochastic model of transcriptional dynamics
#'   \item `"deterministic"` — Deterministic model assuming constant
#'         transcription rates
#'   \item `"dynamical"` — Dynamical model accounting for time-dependent
#'         transcription rates
#' }
#' For example, `"stochastic_UMAP"` or `"dynamical_PCA"`.
#'
#' @section Relationship to CellDimPlot:
#' While \code{\link{CellDimPlot}} can overlay velocity arrows on a dimension reduction
#' plot via its `velocity` parameter, `CellVelocityPlot` is the dedicated
#' velocity visualization function. The key differences are:
#' \itemize{
#'   \item `CellVelocityPlot` delegates to
#'         \code{\link[plotthis:VelocityPlot]{plotthis::VelocityPlot()}},
#'         which supports multiple plot types: `"arrow"` (default), `"grid"`,
#'         and `"stream"`.
#'   \item `CellVelocityPlot` requires both `reduction` and `v_reduction` as
#'         explicit arguments (they have no defaults).
#'   \item `CellVelocityPlot` does not support the full set of `DimPlot`
#'         overlays (marks, density, lineages, etc.) — it focuses solely on
#'         velocity visualization.
#' }
#'
#' @section Velocity plot types:
#' Passed via `...` to \pkg{plotthis}'s `VelocityPlot()`, the `plot_type`
#' parameter controls how velocities are rendered:
#' \itemize{
#'   \item **`"arrow"`** (default) — An arrow is drawn from each cell's
#'         current position to its predicted future position. Best for small
#'         to medium datasets where individual cell trajectories are of
#'         interest.
#'   \item **`"grid"`** — The embedding space is divided into a grid, and
#'         average velocity vectors are computed per grid cell. Reduces
#'         visual clutter for large datasets and shows regional trends.
#'   \item **`"stream"`** — Streamlines are computed from the velocity
#'         field, showing continuous flow trajectories. Best for visualizing
#'         differentiation paths and developmental trajectories.
#' }
#'
#' @param object A Seurat object, a Giotto object, a path to an `.h5ad` file,
#'   or an opened `H5File` from the \pkg{hdf5r} package.
#' @param reduction Name of the main dimension reduction that defines cell
#'   positions in the plot (e.g., `"umap"`, `"pca"`, `"tsne"`). This
#'   reduction must already exist in the object. Unlike \code{\link{CellDimPlot}},
#'   this parameter is required — there is no default.
#' @param v_reduction Name of the velocity reduction that encodes RNA velocity
#'   vectors (e.g., `"stochastic_UMAP"`, `"dynamical_PCA"`). Only the first
#'   two dimensions are used. This reduction must already exist in the object
#'   and is typically generated by RNA velocity analysis tools such as
#'   \pkg{velocyto.R} or \pkg{scVelo} (Python, exported to h5ad).
#' @param spat_unit Spatial unit name for Giotto objects (e.g., `"cell"`).
#'   Ignored for Seurat and h5ad inputs. If `NULL`, the default spatial unit
#'   is auto-detected via `GiottoClass::set_default_spat_unit()`.
#' @param feat_type Feature type name for Giotto objects (e.g., `"rna"`,
#'   `"dna"`, `"protein"`). Ignored for Seurat and h5ad inputs. If `NULL`, the
#'   default feature type is auto-detected via
#'   `GiottoClass::set_default_feat_type()`.
#' @param group_by Metadata column name used to color cells by group (e.g.,
#'   `"CellType"`, `"Phase"`). Unlike \code{\link{CellDimPlot}}, this parameter is
#'   optional — when `NULL`, all cells are plotted in a uniform color, which
#'   can be useful for focusing on the velocity flow patterns without visual
#'   distraction from cluster colors.
#' @param ident Alias for `group_by`, provided for compatibility with Seurat's
#'   naming convention. When both `group_by` and `ident` are specified, they
#'   must be identical.
#' @param ... Additional arguments passed to
#'   \code{\link[plotthis:VelocityPlot]{plotthis::VelocityPlot()}}. Key
#'   parameters include:
#'   \itemize{
#'     \item **Plot type:** `plot_type` — `"arrow"` (default), `"grid"`,
#'           or `"stream"` (see the \emph{Velocity plot types} section)
#'     \item **Appearance:** `pt_size`, `pt_alpha`, `arrow_length`,
#'           `arrow_thickness`, `palette`, `palcolor`, `theme`, `theme_args`,
#'           `legend.position`
#'     \item **Grid parameters (for `plot_type = "grid"`):** `grid_n`,
#'           `grid_color`, `grid_alpha`, `grid_lwd`
#'     \item **Stream parameters (for `plot_type = "stream"`):**
#'           `stream_n`, `stream_color`, `stream_alpha`, `stream_lwd`
#'     \item **Layout:** `split_by`, `combine`, `nrow`, `ncol`
#'   }
#' @return A `ggplot` object.
#' @note
#' **Velocity computation prerequisites:** This function only *visualizes*
#' pre-computed RNA velocity. Velocity must be calculated separately using
#' tools such as \pkg{velocyto.R} (R) or \pkg{scVelo} (Python). The computed
#' velocity reduction must already exist in the object before calling
#' `CellVelocityPlot`.
#'
#' **Velocity overlay in CellDimPlot:** For a quick velocity overlay on
#' a standard dimension reduction plot (with all `DimPlot` features like
#' highlighting, density contours, and group marks), use
#' `CellDimPlot(velocity = "stochastic_UMAP")` instead. `CellVelocityPlot`
#' is best when you need grid or stream visualizations specifically.
#'
#' **Required parameters:** Both `reduction` and `v_reduction` are required
#' and have no defaults — unlike `CellDimPlot()` where `reduction` can be
#' auto-detected for Seurat objects. This is because velocity visualization
#' intrinsically requires two distinct reductions, and auto-detection of the
#' velocity reduction is not reliable across analysis pipelines.
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[plotthis:VelocityPlot]{plotthis::VelocityPlot()}} —
#'         The underlying plotting engine supporting arrow, grid, and stream
#'         plot types
#'   \item \code{\link{CellDimPlot}} — Standard dimension reduction plot; can overlay
#'         velocity arrows via the `velocity` parameter
#'   \item \code{\link{CellStatPlot}} — Statistical summaries and comparisons
#' }
#'
#' @export
#' @importFrom SeuratObject Embeddings Reductions
#' @importFrom plotthis VelocityPlot
#' @details See:
#' * <https://pwwang.github.io/scplotter/articles/Working_with_anndata_h5ad_files.html>
#'
#' for examples of using this function with .h5ad files.
#' @examples
#' \donttest{
#' data(pancreas_sub)
#'
#' # Arrow plot (default) — each cell shows its predicted direction
#' CellVelocityPlot(pancreas_sub, reduction = "PCA", v_reduction = "stochastic_PCA")
#'
#' # Grid plot — velocity vectors averaged over grid cells
#' CellVelocityPlot(pancreas_sub, reduction = "PCA", v_reduction = "stochastic_PCA",
#'   plot_type = "grid")
#'
#' # Stream plot — continuous flow trajectories
#' CellVelocityPlot(pancreas_sub, reduction = "PCA", v_reduction = "stochastic_PCA",
#'   plot_type = "stream")
#'
#' # With group coloring
#' CellVelocityPlot(pancreas_sub, reduction = "PCA", v_reduction = "stochastic_PCA",
#'   group_by = "SubCellType")
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
        is.null(group_by) || is.null(ident) || identical(group_by, ident))
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
        is.null(group_by) || is.null(ident) || identical(group_by, ident))
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
        is.null(group_by) || is.null(ident) || identical(group_by, ident))
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
        is.null(group_by) || is.null(ident) || identical(group_by, ident))
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
