#' Internal dispatcher for FeatureStatPlot — plot from a composed data frame
#'
#' @description
#' This is the internal workhorse of \code{FeatureStatPlot}. After the S3
#' methods (\code{FeatureStatPlot.Seurat}, \code{FeatureStatPlot.giotto}, etc.)
#' extract expression data and metadata into a unified data frame, this function
#' handles the common logic: pivoting wide feature columns to long format,
#' optionally downsampling cells within identity groups, and dispatching to the
#' appropriate \pkg{plotthis} plotting function based on \code{plot_type}.
#'
#' Named feature lists are handled by converting the names to a
#' \code{rows_data} data frame with \code{rows_split_by}, enabling automatic
#' row grouping in heatmap and dot plots.
#' @inheritParams FeatureStatPlot
#' @inheritDotParams FeatureStatPlot
#' @param data A data frame containing the feature expression columns, metadata
#'   columns (cell identities, groupings), and optionally dimension reduction
#'   coordinates.
#' @param should_shrink Logical. If \code{TRUE}, the data frame is subset to
#'   only the columns needed for the plot (features, ident, group_by,
#'   split_by, and optionally dims). This prevents unnecessary data from
#'   being passed through the plotting pipeline.
#' @param should_pivot Logical. If \code{TRUE}, the wide-format feature
#'   columns are pivoted to long format with \code{.features} and
#'   \code{.value} columns, and features are used for faceting.
#' @return A ggplot object (or a \code{patchwork} object when \code{split_by}
#'   generates multiple plots and \code{combine = TRUE}), or a list of ggplot
#'   objects if \code{combine = FALSE}.
#' @keywords internal
#' @importFrom rlang %||% syms sym
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr summarise %>%
.feature_stat_plot <- function(
    data, features, plot_type, should_shrink, should_pivot, downsample,
    graph = NULL, bg_cutoff = 0, dims = 1:2, rows_name = "Features",
    ident = NULL, agg = mean, group_by = NULL, pos_only = c("no", "any", "all"),
    split_by = NULL, facet_by = NULL, xlab = NULL, ylab = NULL, x_text_angle = NULL,
    ...
) {
    pos_only <- match.arg(pos_only)
    unlisted_features <- unname(unlist(features))
    if (pos_only == "all") {
        data <- data[rowSums(data[, unlisted_features, drop = FALSE] > 0) == length(unlisted_features), , drop = FALSE]
    } else if (pos_only == "any") {
        data <- data[rowSums(data[, unlisted_features, drop = FALSE] > 0) > 0, , drop = FALSE]
    }
    if (should_shrink) {
        dims <- if (is.null(dims)) NULL else colnames(data)[dims]
        selected_columns <- c(dims, ident, unlisted_features, group_by, if (!isTRUE(split_by)) split_by)
        nonexisting_columns <- setdiff(selected_columns, colnames(data))
        if (length(nonexisting_columns) > 0) {
            stop("[FeatureStatPlot] The following columns are not found in the object: ", paste(nonexisting_columns, collapse = ", "))
        }
        data <- data[, selected_columns, drop = FALSE]
    }
    if (should_pivot) {
        data <- pivot_longer(data, cols = unlisted_features, names_to = ".features", values_to = ".value")

        if (isTRUE(split_by)) {
            split_by <- ".features"
            facet_by <- NULL
        } else {
            facet_by <- ".features"
        }
    }
    downsample_data <- function() {
        if (!is.null(downsample) && !is.null(ident)) {
            if (downsample > 0 && downsample <= 1) {
                slice_sample(data, prop = downsample, by = !!sym(ident))
            } else {  # downsample > 1
                slice_sample(data, n = downsample, by = !!sym(ident))
            }
        } else {
            data
        }
    }

    if (plot_type == "violin") {
        data <- downsample_data()
        ViolinPlot(
            data, x = ident, y = ".value", group_by = group_by, split_by = split_by, facet_by = facet_by,
            xlab = xlab %||% "", ylab = ylab %||% "", x_text_angle = x_text_angle %||% 45, ...)
    } else if (plot_type == "box") {
        data <- downsample_data()
        BoxPlot(
            data, x = ident, y = ".value", group_by = group_by, split_by = split_by, facet_by = facet_by,
            xlab = xlab %||% "", ylab = ylab %||% "", x_text_angle = x_text_angle %||% 45, ...)
    } else if (plot_type == "bar") {
        data <- data %>% dplyr::group_by(!!!syms(unique(c(ident, group_by, split_by, ".features")))) %>%
            summarise(.value = agg(!!sym(".value")))
        BarPlot(
            data, x = ident, y = ".value", group_by = group_by, split_by = split_by, facet_by = facet_by,
            xlab = xlab %||% "", ylab = ylab %||% "", x_text_angle = x_text_angle %||% 45, ...)
    } else if (plot_type == "ridge") {
        data <- downsample_data()
        RidgePlot(
            data, x = ".value", group_by = ident, split_by = split_by, facet_by = facet_by,
            xlab = xlab %||% "", ylab = ylab %||% "", ...)
    } else if (plot_type == "dim") {
        FeatureDimPlot(
            data, dims = dims, features = unlisted_features, graph = graph, bg_cutoff = bg_cutoff,
            split_by = split_by, facet_by = facet_by, xlab = xlab, ylab = ylab, ...)
    } else if (plot_type %in% c("heatmap", "dot")) {
        args <- rlang::dots_list(...)
        args$data <- data
        args$rows_by <- unlisted_features
        args$columns_by <- ident
        args$rows_name <- rows_name
        args$split_by <- split_by
        args$show_row_names <- args$show_row_names %||% TRUE
        if (plot_type == "heatmap") {
            args$show_column_names <- args$show_column_names %||% !identical(args$cell_type, "bars")
            args$values_fill <- args$values_fill %||% 0
        } else {  # dot
            args$show_column_names <- args$show_column_names %||% TRUE
            args$cell_type <- "dot"
            args$name <- args$name %||% "Expression Level"
            args$dot_size <- args$dot_size %||% function(x) sum(x > 0, na.rm = TRUE) / length(x)
            args$dot_size_name <- args$dot_size_name %||% "Percent Expressed"
            args$add_reticle <- args$add_reticle %||% TRUE
            args$cluster_rows <- args$cluster_rows %||% FALSE
            args$cluster_columns <- args$cluster_columns %||% FALSE
            args$row_names_side <- args$row_names_side %||% "left"
        }
        if (!is.null(names(features))) {
            # handle named features.
            # names will be used as rows_split_by
            if (!is.null(args$rows_split_by)) {
                stop("[FeatureStatPlot] Cannot use both `rows_split_by` and named features.")
            }
            feature_groups <- lapply(names(features), function(x) rep(x, length(features[[x]])))
            fdata <- data.frame(
                Feautures = unlisted_features,
                FeatureGroups = do_call(c, feature_groups)
            )
            args$rows_split_name <- args$rows_split_name %||% "Feature Groups"
            names(fdata) <- c(rows_name, args$rows_split_name)
            args$rows_split_by <- args$rows_split_name
            if (!is.null(args$rows_data)) {
                args$rows_data <- merge(fdata, args$rows_data, by = rows_name, all.x = TRUE)
            } else {
                args$rows_data <- fdata
            }
        }
        do_call(Heatmap, args)
    } else if (plot_type == "cor") {
        if (length(unlisted_features) < 2) {
            stop("[FeatureStatPlot] The number of features should be at least 2 for correlation plot.")
        }
        if (length(unlisted_features) == 2) {
            CorPlot(data, x = unlisted_features[1], y = unlisted_features[2], group_by = ident, split_by = split_by,
                facet_by = facet_by, xlab = xlab, ylab = ylab, ...)
        } else {
            CorPairsPlot(
                data, columns = unlisted_features, group_by = ident, split_by = split_by, facet_by = facet_by, ...)
        }
    }
}

#' Visualize feature expression and statistics across cell groups
#'
#' @description
#' A central question in single-cell analysis is how features — genes, gene
#' signatures, module scores, or other molecular measurements — vary across cell
#' types, conditions, or experimental groups. \code{FeatureStatPlot} answers this
#' question by providing eight complementary visualization types, each suited to
#' a different analytical perspective:
#' \itemize{
#'   \item \strong{violin} — Violin plot showing the full distribution of feature
#'     values per identity group. Best for comparing expression distributions
#'     and detecting bimodality or outliers.
#'   \item \strong{box} — Box plot summarizing feature values with quartiles and
#'     outliers. A compact alternative to the violin plot.
#'   \item \strong{bar} — Bar chart of aggregated feature values (default: mean)
#'     per group. Useful for summary-level comparisons with error bars.
#'   \item \strong{ridge} — Ridge (joy) plot showing density curves per group.
#'     Effective when comparing many groups or when distribution shape matters.
#'   \item \strong{dim} — Dimensionality reduction plot (UMAP, t-SNE, PCA) with
#'     cells colored by feature expression. Reveals spatial patterns of gene
#'     expression in the reduced space.
#'   \item \strong{cor} — Correlation plot between two features (scatter with
#'     fitted line and annotations) or among multiple features (pairs plot).
#'     Reveals co-expression relationships.
#'   \item \strong{heatmap} — Heatmap of feature expression across identity
#'     groups. Supports rich annotations (row/column metadata, bar charts,
#'     pie charts, violin plots) and flexible clustering. The go-to choice
#'     for visualizing many features across many groups.
#'   \item \strong{dot} — Dot plot (a shortcut for heatmap with
#'     \code{cell_type = "dot"}) where dot size reflects the fraction of
#'     expressing cells and dot color reflects mean expression. A compact,
#'     publication-ready format for marker gene visualization.
#' }
#'
#' The function is an S3 generic with methods for \pkg{Seurat} objects,
#' \pkg{Giotto} objects, AnnData (.h5ad) file paths, and \code{H5File} objects
#' (via \pkg{hdf5r}). Each method extracts the relevant expression matrix and
#' metadata, then delegates to the internal \code{.feature_stat_plot()} which
#' dispatches to the appropriate \pkg{plotthis} plotting function.
#'
#' @section Input objects:
#' \code{FeatureStatPlot} supports four input types, each handled by its own
#' S3 method:
#' \itemize{
#'   \item \strong{Seurat} (\code{FeatureStatPlot.Seurat}) — Uses
#'     \code{\link[SeuratObject:GetAssayData]{SeuratObject::GetAssayData()}} to
#'     extract expression data and \code{@meta.data} for cell metadata.
#'     Reductions are accessed via \code{\link[SeuratObject:Embeddings]{SeuratObject::Embeddings()}}.
#'   \item \strong{Giotto} (\code{FeatureStatPlot.giotto}) — Uses
#'     \code{GiottoClass::getExpression()} and
#'     \code{GiottoClass::getCellMetadata()}. Requires \code{spat_unit} and
#'     \code{feat_type} parameters.
#'   \item \strong{.h5ad path} (\code{FeatureStatPlot.character}) — Opens the
#'     file via \pkg{hdf5r}, then delegates to \code{FeatureStatPlot.H5File}.
#'     Currently only \code{.h5ad} (AnnData) files are supported.
#'   \item \strong{H5File} (\code{FeatureStatPlot.H5File}) — Reads expression
#'     from \code{X} (or \code{layers/<assay>}), cell metadata from \code{obs},
#'     and reductions from \code{obsm}.
#' }
#'
#' @section Feature specification:
#' Features can be provided as:
#' \itemize{
#'   \item \strong{A character vector} — e.g., \code{c("Sox9", "Neurog3", "Ins1")}.
#'     All features are treated equally.
#'   \item \strong{A named list} — e.g.,
#'     \code{list(Ductal = c("Sox9", "Anxa2"), Endocrine = c("Ins1", "Gcg"))}.
#'     The names are used as row group annotations in heatmap and dot plots
#'     (via \code{rows_split_by}), automatically creating a \code{rows_data}
#'     data frame that maps each feature to its group. This is especially
#'     useful for organizing marker genes by cell type.
#' }
#' When \code{pos_only} is \code{"any"} or \code{"all"}, cells are filtered
#' based on whether they have positive values for the specified features.
#' For named feature lists, the filter applies to the flattened set of all
#' features.
#'
#' @section Faceting and splitting behavior:
#' For most plot types (\code{"violin"}, \code{"box"}, \code{"bar"},
#' \code{"ridge"}), features are automatically faceted — each feature appears
#' in its own panel. The \code{facet_by} parameter is therefore restricted
#' for these types. Use \code{split_by} (or \code{split_by = TRUE} to split
#' by feature) for creating separate plots per category. For \code{"dim"},
#' \code{"heatmap"}, \code{"dot"}, and \code{"cor"} plots, features are
#' incorporated into a single visualization and \code{split_by}/\code{facet_by}
#' behave normally.
#'
#' @param object An object containing single-cell expression data. Supported
#'   types: a \pkg{Seurat} object, a \pkg{Giotto} object, a character path to
#'   an \code{.h5ad} file, or an opened \code{H5File} object from the
#'   \pkg{hdf5r} package.
#' @param features A character vector or a named list of character vectors
#'   specifying the features to plot. Features can be gene names, gene
#'   signature scores, or any column present in the expression matrix or
#'   metadata. Named lists (e.g., \code{list(Beta = c("Ins1", "Ins2"))})
#'   enable automatic row grouping in heatmap and dot plots.
#' @param plot_type Character. The type of plot to generate. One of:
#'   \code{"violin"}, \code{"box"}, \code{"bar"}, \code{"ridge"},
#'   \code{"dim"}, \code{"cor"}, \code{"heatmap"}, or \code{"dot"}.
#'   See the Description section for guidance on choosing a plot type.
#'   Default: \code{"violin"}.
#' @param spat_unit Character. The spatial unit to extract data from.
#'   Only applicable to \pkg{Giotto} objects. Default: \code{NULL}
#'   (auto-detected by Giotto).
#' @param feat_type Character. The feature type to extract (e.g.,
#'   \code{"rna"}, \code{"dna"}, \code{"protein"}). Only applicable to
#'   \pkg{Giotto} objects. Default: \code{NULL} (auto-detected by Giotto).
#' @param reduction Character. Name of the dimensionality reduction to use
#'   (e.g., \code{"umap"}, \code{"tsne"}, \code{"pca"}). Required when
#'   \code{plot_type = "dim"}; optional for other types where the reduction
#'   coordinates can be used as feature values. For Seurat objects, defaults
#'   to the default reduction if available. Default: \code{NULL}.
#' @param dims Integer vector of length 2. The dimensions (columns) of the
#'   reduction to plot on the x and y axes. Only used when
#'   \code{plot_type = "dim"}. Default: \code{1:2}.
#' @param rows_name Character. The column name used to identify feature rows
#'   in the heatmap and as the key when merging \code{rows_data}. Only used
#'   when \code{plot_type = "heatmap"} or \code{"dot"}. Default:
#'   \code{"Features"}.
#' @param graph Character or \code{TRUE}. A graph (nearest-neighbor network)
#'   name to overlay edges between connected cells on the dim plot. For
#'   Seurat objects, the name should exist in \code{object@graphs}. For
#'   Giotto objects, the name should exist in the nearest network list. If
#'   \code{TRUE}, the first available graph is used. Only used when
#'   \code{plot_type = "dim"}. Default: \code{NULL} (no edges).
#' @param bg_cutoff Numeric. Expression cutoff for the background in dim
#'   plots. Cells with expression below this value are shown in the
#'   background color (typically gray). Set to \code{-Inf} to color all
#'   cells. Only used when \code{plot_type = "dim"}. Default: \code{0}.
#' @param pos_only Character. Whether to restrict to cells with positive
#'   feature values:
#'   \itemize{
#'     \item \code{"no"} — Include all cells (default).
#'     \item \code{"any"} — Include cells where at least one feature is > 0.
#'     \item \code{"all"} — Include only cells where all features are > 0.
#'   }
#'   For named feature lists, filtering is applied to all flattened features.
#' @param ident Character. The metadata column name identifying cell groups
#'   (e.g., \code{"CellType"}, \code{"SubCellType"}, \code{"seurat_clusters"}).
#'   Used as the x-axis for violin, box, bar, and ridge plots; as the
#'   heatmap columns for heatmap and dot plots; and as the grouping variable
#'   for correlation plots. For Seurat objects, defaults to the active
#'   identity (\code{"Identity"}). Required for non-dim plots on Giotto and
#'   H5File objects. Default: \code{NULL}.
#' @param assay Character. The assay name to extract feature data from
#'   (e.g., \code{"RNA"}, \code{"SCT"}, \code{"integrated"}). For Seurat
#'   objects, passed to \code{\link[SeuratObject:GetAssayData]{SeuratObject::GetAssayData()}}.
#'   For H5File objects, determines whether to read from \code{X} (when
#'   \code{assay = "RNA"}) or \code{layers/<assay>}. Not applicable to
#'   Giotto objects. Default: \code{NULL}.
#' @param layer Character. The layer name within the assay to extract data
#'   from (e.g., \code{"data"}, \code{"counts"}, \code{"scale.data"}). For
#'   Seurat objects, passed to \code{\link[SeuratObject:GetAssayData]{SeuratObject::GetAssayData()}}.
#'   For Giotto objects, passed to \code{GiottoClass::getExpression()}.
#'   Default: \code{NULL}.
#' @param agg Function. The aggregation function applied when
#'   \code{plot_type = "bar"}. Common choices: \code{mean} (default),
#'   \code{median}, \code{sum}. Applied within each group defined by
#'   \code{ident}, \code{group_by}, and \code{split_by}.
#' @param downsample Numeric. Number or fraction of cells to downsample to
#'   per identity group. Used for \code{"violin"}, \code{"box"}, and
#'   \code{"ridge"} plot types to reduce overplotting in large datasets:
#'   \itemize{
#'     \item If \code{downsample > 1}: Exact number of cells per group.
#'     \item If \code{0 < downsample <= 1}: Fraction of cells per group.
#'   }
#'   Downsampling is performed within each identity group
#'   (via \code{dplyr::slice_sample(by = ident)}). Default: \code{NULL}
#'   (no downsampling).
#' @param group_by Character. A metadata column name to further subdivide
#'   cells within each identity group (e.g., coloring by treatment within
#'   cell type). Works with \code{"violin"}, \code{"box"}, \code{"bar"},
#'   and \code{"ridge"} plot types. Default: \code{NULL}.
#' @param split_by Character vector or \code{TRUE}. Metadata column name(s)
#'   to split the data into separate plots (one per unique value).
#'   If \code{TRUE}, splits by the features themselves, creating one plot
#'   per feature. Multiple columns are concatenated. Default: \code{NULL}.
#' @param facet_by Character vector. Metadata column name(s) to facet the
#'   data within a plot. Note: for \code{"violin"}, \code{"box"},
#'   \code{"bar"}, and \code{"ridge"} plot types, \code{facet_by} should
#'   always be \code{NULL} because the plot is already faceted by features.
#'   Works normally with \code{"dim"}, \code{"heatmap"}, \code{"dot"}, and
#'   \code{"cor"} plot types. Default: \code{NULL}.
#' @param xlab Character. Custom x-axis label. Default: \code{NULL}
#'   (auto-generated based on plot type).
#' @param ylab Character. Custom y-axis label. Default: \code{NULL}
#'   (auto-generated based on plot type).
#' @param x_text_angle Numeric. Angle (in degrees) for x-axis text labels.
#'   Used for \code{"violin"}, \code{"box"}, and \code{"bar"} plot types.
#'   Default: \code{NULL} (defaults to \code{45}).
#' @param ... Additional arguments passed to the underlying \pkg{plotthis}
#'   plotting function, determined by \code{plot_type}:
#'   \describe{
#'     \item{\code{"violin"}}{\code{\link[plotthis:ViolinPlot]{plotthis::ViolinPlot()}}}
#'     \item{\code{"box"}}{\code{\link[plotthis:BoxPlot]{plotthis::BoxPlot()}}}
#'     \item{\code{"bar"}}{\code{\link[plotthis:BarPlot]{plotthis::BarPlot()}}}
#'     \item{\code{"ridge"}}{\code{\link[plotthis:RidgePlot]{plotthis::RidgePlot()}}}
#'     \item{\code{"dim"}}{\code{\link[plotthis:FeatureDimPlot]{plotthis::FeatureDimPlot()}}}
#'     \item{\code{"heatmap"}}{\code{\link[plotthis:Heatmap]{plotthis::Heatmap()}}}
#'     \item{\code{"dot"}}{\code{\link[plotthis:Heatmap]{plotthis::Heatmap()}} with \code{cell_type = "dot"}}
#'     \item{\code{"cor"}}{\code{\link[plotthis:CorPlot]{plotthis::CorPlot()}} (2 features) or
#'       \code{\link[plotthis:CorPairsPlot]{plotthis::CorPairsPlot()}} (3+ features)}
#'   }
#'   Common arguments include \code{palette}, \code{flip}, \code{add_bg},
#'   \code{add_point}, \code{add_box}, \code{stack}, \code{comparisons},
#'   \code{theme}, \code{legend.position}, and many more — see the
#'   documentation of the specific \pkg{plotthis} function for details.
#' @return A ggplot object (or a \code{patchwork} object when
#'   \code{split_by} generates multiple plots and \code{combine = TRUE}), or
#'   a list of ggplot objects if \code{combine = FALSE}. The specific return
#'   type depends on the underlying \pkg{plotthis} function dispatched by
#'   \code{plot_type}.
#' @note
#' \itemize{
#'   \item For \code{"violin"}, \code{"box"}, \code{"bar"}, and
#'     \code{"ridge"} plot types, the data is automatically pivoted from wide
#'     to long format and features are faceted. Do not use \code{facet_by}
#'     with these types — use \code{split_by} instead.
#'   \item For \code{"heatmap"} and \code{"dot"} plot types, named feature
#'     lists are automatically converted to a \code{rows_data} data frame
#'     with \code{rows_split_by} set. You can also provide your own
#'     \code{rows_data} for richer annotations (e.g., TF status, CSPA
#'     membership).
#'   \item The \code{dot} plot type is a convenience shortcut for
#'     \code{plot_type = "heatmap"} with \code{cell_type = "dot"}. It
#'     pre-configures sensible defaults: dot size = fraction of expressing
#'     cells, dot color = mean expression, reticle added, no clustering.
#'   \item For \code{"dim"} plots with \code{graph}, edges are drawn between
#'     connected cells. This helps reveal whether feature expression follows
#'     the neighborhood structure (e.g., continuous gradients vs scattered
#'     expression).
#'   \item When \code{ident} is not specified for Seurat objects, the active
#'     identity (\code{Idents(object)}) is used automatically.
#'   \item Giotto objects require \code{spat_unit} and \code{feat_type} to
#'     locate the correct expression matrix and metadata.
#' }
#' @seealso
#' \code{\link[plotthis:ViolinPlot]{plotthis::ViolinPlot()}},
#' \code{\link[plotthis:BoxPlot]{plotthis::BoxPlot()}},
#' \code{\link[plotthis:BarPlot]{plotthis::BarPlot()}},
#' \code{\link[plotthis:RidgePlot]{plotthis::RidgePlot()}},
#' \code{\link[plotthis:FeatureDimPlot]{plotthis::FeatureDimPlot()}},
#' \code{\link[plotthis:Heatmap]{plotthis::Heatmap()}},
#' \code{\link[plotthis:CorPlot]{plotthis::CorPlot()}},
#' \code{\link[plotthis:CorPairsPlot]{plotthis::CorPairsPlot()}},
#' \code{\link{CellDimPlot}}, \code{\link{CellStatPlot}}
#' @export
#' @importFrom rlang %||%
#' @importFrom SeuratObject GetAssayData Embeddings Graphs Reductions Idents
#' @importFrom plotthis ViolinPlot BoxPlot BarPlot DotPlot RidgePlot FeatureDimPlot Heatmap CorPlot CorPairsPlot
#' @details
#' See the vignettes for examples with non-Seurat objects:
#' \itemize{
#'   \item \href{https://pwwang.github.io/scplotter/articles/Giotto_Visium.html}{Giotto Visium}
#'   \item \href{https://pwwang.github.io/scplotter/articles/Giotto_Xenium.html}{Giotto Xenium}
#'   \item \href{https://pwwang.github.io/scplotter/articles/Working_with_anndata_h5ad_files.html}{Working with AnnData .h5ad files}
#' }
#' @examples
#' \donttest{
#' data(pancreas_sub)
#'
#' FeatureStatPlot(pancreas_sub, features = c("G2M_score", "nCount_RNA"),
#'    ident = "SubCellType", facet_scales = "free_y")
#' FeatureStatPlot(pancreas_sub, features = c("G2M_score", "nCount_RNA"),
#'    ident = "SubCellType", plot_type = "box", facet_scales = "free_y")
#' FeatureStatPlot(pancreas_sub, features = c("G2M_score", "nCount_RNA"),
#'    ident = "SubCellType", plot_type = "bar", facet_scales = "free_y")
#' FeatureStatPlot(pancreas_sub, features = c("G2M_score", "nCount_RNA"),
#'    ident = "SubCellType", plot_type = "ridge", flip = TRUE, facet_scales = "free_y")
#' FeatureStatPlot(pancreas_sub, features = c("G2M_score", "nCount_RNA"),
#'    ident = "SubCellType", facet_scales = "free_y", add_point = TRUE)
#' FeatureStatPlot(pancreas_sub, features = c("G2M_score", "nCount_RNA"),
#'    ident = "SubCellType", facet_scales = "free_y", add_trend = TRUE)
#' FeatureStatPlot(pancreas_sub, features = c("G2M_score", "nCount_RNA"),
#'    ident = "SubCellType", facet_scales = "free_y", add_stat = mean)
#' FeatureStatPlot(pancreas_sub, features = c("G2M_score", "nCount_RNA"),
#'    ident = "SubCellType", facet_scales = "free_y", group_by = "Phase")
#' if (requireNamespace("ggpubr", quietly = TRUE)) {
#'   # https://github.com/kassambara/ggpubr/issues/751
#'   library(ggpubr)
#'   FeatureStatPlot(
#'      subset(pancreas_sub,
#'          subset = SubCellType %in% c("Ductal", "Ngn3 low EP", "Ngn3 high EP")),
#'      features = c("G2M_score"),
#'      ident = "SubCellType", group_by = "Phase", comparisons = TRUE)
#' }
#' FeatureStatPlot(pancreas_sub, features = c("Rbp4", "Pyy"), ident = "SubCellType",
#'    add_bg = TRUE, add_box = TRUE, stack = TRUE)
#' # Use `pos_only` to include only cells with positive expression of all features
#' FeatureStatPlot(pancreas_sub, features = c("Rbp4", "Pyy"), ident = "SubCellType",
#'    add_bg = TRUE, add_box = TRUE, stack = TRUE, pos_only = "all")
#' FeatureStatPlot(pancreas_sub, features = c(
#'        "Sox9", "Anxa2", "Bicc1", # Ductal
#'        "Neurog3", "Hes6", # EPs
#'        "Fev", "Neurod1", # Pre-endocrine
#'        "Rbp4", "Pyy", # Endocrine
#'        "Ins1", "Gcg", "Sst", "Ghrl" # Beta, Alpha, Delta, Epsilon
#'    ), ident = "SubCellType", add_bg = TRUE, stack = TRUE,
#'    legend.position = "top", legend.direction = "horizontal")
#' FeatureStatPlot(pancreas_sub, plot_type = "box", features = c(
#'       "Sox9", "Anxa2", "Bicc1", # Ductal
#'       "Neurog3", "Hes6", # EPs
#'       "Fev", "Neurod1", # Pre-endocrine
#'       "Rbp4", "Pyy", # Endocrine
#'       "Ins1", "Gcg", "Sst", "Ghrl" # Beta, Alpha, Delta, Epsilon
#'    ), ident = "SubCellType", add_bg = TRUE, stack = TRUE, flip = TRUE,
#'    legend.position = "top", legend.direction = "horizontal")
#' # Use splitting instead of facetting
#' FeatureStatPlot(pancreas_sub, features = c("Neurog3", "Rbp4", "Ins1"),
#'    ident = "CellType", split_by = TRUE)
#'
#' FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "G2M_score", reduction = "UMAP")
#' FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "G2M_score", reduction = "UMAP",
#'    bg_cutoff = -Inf)
#' FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "G2M_score", reduction = "UMAP",
#'    theme = "theme_blank")
#' FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "G2M_score", reduction = "UMAP",
#'    theme = ggplot2::theme_classic, theme_args = list(base_size = 16))
#'
#' # Label and highlight cell points
#' FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "Rbp4", reduction = "UMAP",
#'    highlight = 'SubCellType == "Delta"')
#' FeatureStatPlot(pancreas_sub, plot_type = "dim",
#'    features = "Rbp4", split_by = "Phase", reduction = "UMAP",
#'    highlight = TRUE, theme = "theme_blank")
#'
#' # Add a density layer
#' FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "Rbp4", reduction = "UMAP",
#'    add_density = TRUE)
#' FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "Rbp4", reduction = "UMAP",
#'    add_density = TRUE, density_filled = TRUE)
#'
#' # Change the plot type from point to the hexagonal bin
#' FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "Rbp4", reduction = "UMAP",
#'    hex = TRUE)
#' FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "Rbp4", reduction = "UMAP",
#'    hex = TRUE, hex_bins = 20)
#'
#' # Show lineages on the plot based on the pseudotime
#' FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "Lineage3", reduction = "UMAP",
#'    lineages = "Lineage3")
#' FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "Lineage3", reduction = "UMAP",
#'    lineages = "Lineage3", lineages_whiskers = TRUE)
#' FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "Lineage3", reduction = "UMAP",
#'    lineages = "Lineage3", lineages_span = 0.1)
#'
#' FeatureStatPlot(pancreas_sub, plot_type = "dim",
#'   features = c("Sox9", "Anxa2", "Bicc1"), reduction = "UMAP",
#'   theme = "theme_blank",
#'   theme_args = list(plot.subtitle = ggplot2::element_text(size = 10),
#'      strip.text = ggplot2::element_text(size = 8))
#' )
#'
#' # Plot multiple features with different scales
#' endocrine_markers <- c("Ins1", "Gcg", "Sst", "Ghrl")
#' FeatureStatPlot(pancreas_sub, endocrine_markers, reduction = "UMAP", plot_type = "dim")
#' FeatureStatPlot(pancreas_sub, endocrine_markers, reduction = "UMAP", lower_quantile = 0,
#'    upper_quantile = 0.8, plot_type = "dim")
#' FeatureStatPlot(pancreas_sub, endocrine_markers, reduction = "UMAP",
#'    lower_cutoff = 1, upper_cutoff = 4, plot_type = "dim")
#' FeatureStatPlot(pancreas_sub, endocrine_markers, reduction = "UMAP", bg_cutoff = 2,
#'    lower_cutoff = 2, upper_cutoff = 4, plot_type = "dim")
#' FeatureStatPlot(pancreas_sub, c("Sst", "Ghrl"), split_by = "Phase", reduction = "UMAP",
#'    plot_type = "dim")
#' FeatureStatPlot(pancreas_sub, features = c("G2M_score", "nCount_RNA"),
#'    ident = "SubCellType", plot_type = "dim", facet_by = "Phase", split_by = TRUE, ncol = 1)
#'
#' # Heatmap
#' features <- c(
#'    "Sox9", "Anxa2", "Bicc1", # Ductal
#'    "Neurog3", "Hes6", # EPs
#'    "Fev", "Neurod1", # Pre-endocrine
#'    "Rbp4", "Pyy", # Endocrine
#'    "Ins1", "Gcg", "Sst", "Ghrl" # Beta, Alpha, Delta, Epsilon
#' )
#' rows_data <- data.frame(
#'    Features = features,  # 'rows_name' default is "Features"
#'    group = c(
#'        "Ductal", "Ductal", "Ductal", "EPs", "EPs", "Pre-endocrine",
#'        "Pre-endocrine", "Endocrine", "Endocrine", "Beta", "Alpha", "Delta", "Epsilon"),
#'    TF = c(TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE,
#'        TRUE, TRUE, TRUE),
#'    CSPA = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE,
#'        FALSE, FALSE, FALSE)
#' )
#' FeatureStatPlot(pancreas_sub, features = features, ident = "SubCellType",
#'    plot_type = "heatmap", name = "Expression Level")
#' FeatureStatPlot(pancreas_sub, features = features, ident = "Phase",
#'    plot_type = "heatmap", name = "Expression Level", columns_split_by = "SubCellType")
#' FeatureStatPlot(pancreas_sub, features = features, ident = "SubCellType",
#'    plot_type = "heatmap", cell_type = "bars", name = "Expression Level")
#' FeatureStatPlot(pancreas_sub, features = features, ident = "SubCellType", cell_type = "dot",
#'    plot_type = "heatmap", name = "Expression Level", dot_size = function(x) sum(x > 0) / length(x),
#'    dot_size_name = "Percent Expressed", add_bg = TRUE, rows_data = rows_data,
#'    show_row_names = TRUE, rows_split_by = "group", cluster_rows = FALSE,
#'    column_annotation = c("Phase", "G2M_score"),
#'    column_annotation_type = list(Phase = "pie", G2M_score = "violin"),
#'    column_annotation_params = list(G2M_score = list(show_legend = FALSE)),
#'    row_annotation = c("TF", "CSPA"),
#'    row_annotation_side = "right",
#'    row_annotation_type = list(TF = "simple", CSPA = "simple"))
#' FeatureStatPlot(pancreas_sub, features = features, ident = "SubCellType", cell_type = "dot",
#'    plot_type = "heatmap", name = "Expression Level", dot_size = function(x) sum(x > 0) / length(x),
#'    dot_size_name = "Percent Expressed", add_bg = TRUE,
#'    rows_data = rows_data, show_column_names = TRUE, rows_split_by = "group",
#'    cluster_rows = FALSE, flip = TRUE, palette = "YlOrRd",
#'    column_annotation = c("Phase", "G2M_score"),
#'    column_annotation_type = list(Phase = "pie", G2M_score = "violin"),
#'    column_annotation_params = list(G2M_score = list(show_legend = FALSE)),
#'    row_annotation = c("TF", "CSPA"),
#'    row_annotation_side = "right",
#'    row_annotation_type = list(TF = "simple", CSPA = "simple"))
#' FeatureStatPlot(pancreas_sub, features = features, ident = "SubCellType", cell_type = "violin",
#'    plot_type = "heatmap", name = "Expression Level", show_row_names = TRUE,
#'    cluster_columns = FALSE, rows_split_by = "group", rows_data = rows_data)
#' FeatureStatPlot(pancreas_sub, features = features, ident = "SubCellType", cell_type = "dot",
#'    plot_type = "heatmap", dot_size = function(x) sum(x > 0) / length(x),
#'    dot_size_name = "Percent Expressed", palette = "viridis", add_reticle = TRUE,
#'    rows_data = rows_data, name = "Expression Level", show_row_names = TRUE,
#'    rows_split_by = "group")
#' # Visualize the markers for each sub-cell type (the markers can overlap)
#' # Say: markers <- Seurat::FindAllMarkers(pancreas_sub, ident = "SubCellType")
#' markers <- data.frame(
#'     avg_log2FC = c(
#'          3.44, 2.93, 2.72, 2.63, 2.13, 1.97, 2.96, 1.92, 5.22, 3.91, 3.64, 4.52,
#'          3.45, 2.45, 1.75, 2.08, 9.10, 4.45, 3.61, 6.30, 4.96, 3.49, 3.91, 3.90,
#'          10.58, 5.84, 4.73, 3.34, 7.22, 4.52, 10.10, 4.25),
#'     cluster = factor(rep(
#'          c("Ductal", "Ngn3 low EP", "Ngn3 high EP", "Pre-endocrine", "Beta",
#'            "Alpha", "Delta", "Epsilon"), each = 4),
#'          levels = levels(pancreas_sub$SubCellType)),
#'     gene = c(
#'          "Cyr61", "Adamts1", "Anxa2", "Bicc1", "1700011H14Rik", "Gsta3", "8430408G22Rik",
#'          "Anxa2", "Ppp1r14a", "Btbd17", "Neurog3", "Gadd45a", "Fev", "Runx1t1", "Hmgn3",
#'          "Cryba2", "Ins2", "Ppp1r1a", "Gng12", "Sytl4", "Irx1", "Tmem27", "Peg10", "Irx2",
#'          "Sst", "Ptprz1", "Arg1", "Frzb", "Irs4", "Mboat4", "Ghrl", "Arg1"
#'     )
#' )
#' FeatureStatPlot(pancreas_sub,
#'   features = unique(markers$gene), ident = "SubCellType", cell_type = "bars",
#'   plot_type = "heatmap", rows_data = markers, rows_name = "gene", rows_split_by = "cluster",
#'   show_row_names = TRUE, show_column_names = TRUE, name = "Expression Level",
#'   cluster_rows = FALSE, cluster_columns = FALSE,
#'   row_annotation_palette = list(.row = "Paired"))
#'
#' # Use plot_type = "dot" to as a shortcut for heatmap with cell_type = "dot"
#' FeatureStatPlot(pancreas_sub, features = features, ident = "SubCellType", plot_type = "dot")
#'
#' named_features <- list(
#'    Ductal = c("Sox9", "Anxa2", "Bicc1"),
#'    EPs = c("Neurog3", "Hes6"),
#'    `Pre-endocrine` = c("Fev", "Neurod1"),
#'    Endocrine = c("Rbp4", "Pyy"),
#'    Beta = "Ins1", Alpha = "Gcg", Delta = "Sst", Epsilon = "Ghrl"
#' )
#' FeatureStatPlot(pancreas_sub, features = named_features, ident = "SubCellType",
#'    plot_type = "heatmap", name = "Expression Level", show_row_names = TRUE)
#'
#' # Correlation plot
#' FeatureStatPlot(pancreas_sub, features = c("Pyy", "Rbp4"), plot_type = "cor",
#'    anno_items = c("eq", "r2", "spearman"))
#' FeatureStatPlot(pancreas_sub, features = c("Ins1", "Gcg", "Sst", "Ghrl"),
#'    plot_type = "cor")
#' }
FeatureStatPlot <- function(
    object, features, plot_type = c("violin", "box", "bar", "ridge", "dim", "cor", "heatmap", "dot"),
    spat_unit = NULL, feat_type = NULL, downsample = NULL, pos_only = c("no", "any", "all"),
    reduction = NULL, graph = NULL, bg_cutoff = 0, dims = 1:2, rows_name = "Features",
    ident = NULL, assay = NULL, layer = NULL, agg = mean, group_by = NULL,
    split_by = NULL, facet_by = NULL, xlab = NULL, ylab = NULL, x_text_angle = NULL, ...
) {
    UseMethod("FeatureStatPlot")
}

#' @export
FeatureStatPlot.giotto <- function(
    object, features, plot_type = c("violin", "box", "bar", "ridge", "dim", "cor", "heatmap", "dot"),
    spat_unit = NULL, feat_type = NULL, downsample = NULL, pos_only = c("no", "any", "all"),
    reduction = NULL, graph = NULL, bg_cutoff = 0, dims = 1:2, rows_name = "Features",
    ident = NULL, assay = NULL, layer = NULL, agg = mean, group_by = NULL,
    split_by = NULL, facet_by = NULL, xlab = NULL, ylab = NULL, x_text_angle = NULL, ...
) {
    plot_type <- match.arg(plot_type)
    if (!is.null(facet_by) && plot_type != "dim") {
        stop("Cannot facet plots because the plots are facetted by the 'features'.")
    }
    stopifnot("[FeatureStatPlot] 'assay' should not be used for Giotto object." = is.null(assay))
    if (plot_type != "dim" && is.null(ident)) {
        stop("[FeatureStatPlot] 'ident' should be provided for non-dim plot for giotto object.")
    }
    # dim plot may use expression for highlighting cells
    # Heatmap may use other variables as annotations, but shrinking only includes minimal columns
    should_shrink <- !plot_type %in% c("dim", "heatmap", "dot")
    should_pivot <- !plot_type %in% c("dim", "heatmap", "dot", "cor")

    spat_unit <- GiottoClass::set_default_spat_unit(
        gobject = object,
        spat_unit = spat_unit
    )

    feat_type <- GiottoClass::set_default_feat_type(
        gobject = object,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    unlisted_features <- unname(unlist(features))
    assay_data <- GiottoClass::getExpression(
        object, values = layer,
        spat_unit = spat_unit, feat_type = feat_type,
        output = "matrix"
    )

    assay_feature <- intersect(unlisted_features, rownames(assay_data))
    assay_data <- t(as.matrix(assay_data[assay_feature, , drop = FALSE]))

    if (is.null(reduction) && plot_type != "dim") {
        data <- cbind(
            GiottoClass::getCellMetadata(
                gobject = object,
                spat_unit = spat_unit,
                feat_type = feat_type,
                output = "data.table"
            ),
            assay_data
        )
    } else {
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
            ),
            assay_data
        )
    }

    if (plot_type == "dim" && !is.null(graph)) {
        # spat_unit	feat_type	nn_type	name
        # cell	rna	sNN	sNN.pca
        graph_info <- GiottoClass::list_nearest_networks(
            gobject = object,
            spat_unit = spat_unit,
            feat_type = feat_type
        )
        if (is.null(graph_info)) {
            stop(
                "[FeatureStatPlot] The object does not have any nearest networks with given ",
                "spat_unit (", spat_unit, ") and feat_type (", feat_type, ")."
            )
        }
        if (isTRUE(graph)) { graph <- graph_info$name[1] }
        graph_info <- graph_info[, graph_info$name == graph, drop = FALSE]
        if (nrow(graph_info) == 0) {
            stop(
                "[FeatureStatPlot] The object does not have network: ", graph,
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

    .feature_stat_plot(
        data = data, features = features, plot_type = plot_type, pos_only = pos_only,
        should_shrink = should_shrink, should_pivot = should_pivot,
        graph = graph, bg_cutoff = bg_cutoff, downsample = downsample,
        dims = dims, rows_name = rows_name, ident = ident, agg = agg,
        group_by = group_by, split_by = split_by, facet_by = facet_by,
        xlab = xlab, ylab = ylab, x_text_angle = x_text_angle, ...
    )
}

#' @export
FeatureStatPlot.Seurat <- function(
    object, features, plot_type = c("violin", "box", "bar", "ridge", "dim", "cor", "heatmap", "dot"),
    spat_unit = NULL, feat_type = NULL, downsample = NULL, pos_only = c("no", "any", "all"),
    reduction = NULL, graph = NULL, bg_cutoff = 0, dims = 1:2, rows_name = "Features",
    ident = NULL, assay = NULL, layer = NULL, agg = mean, group_by = NULL,
    split_by = NULL, facet_by = NULL, xlab = NULL, ylab = NULL, x_text_angle = NULL, ...
) {
    plot_type <- match.arg(plot_type)
    if (!is.null(facet_by) && plot_type != "dim") {
        stop("Cannot facet plots because the plots are facetted by the 'features'.")
    }
    stopifnot("[FeatureStatPlot] 'spat_unit' and 'feat_type' should not be used for Seurat object." = is.null(spat_unit) && is.null(feat_type))

    reduction <- reduction %||% (
        if(is.null(Reductions(object))) NULL else default_dimreduc(object)
    )
    # dim plot may use expression for highlighting cells
    # Heatmap may use other variables as annotations, but shrinking only includes minimal columns
    should_shrink <- !plot_type %in% c("dim", "heatmap", "dot")
    should_pivot <- !plot_type %in% c("dim", "heatmap", "dot", "cor")

    unlisted_features <- unname(unlist(features))
    assay_data <- GetAssayData(object, assay = assay, layer = layer)
    assay_feature <- intersect(unlisted_features, rownames(assay_data))
    assay_data <- t(as.matrix(assay_data[assay_feature, , drop = FALSE]))
    if (is.null(reduction)) {
        data <- cbind(object@meta.data, assay_data)
    } else {
        emb <- Embeddings(object, reduction = reduction)
        data <- cbind(
            emb,
            object@meta.data[rownames(emb), , drop = FALSE],
            assay_data[rownames(emb), , drop = FALSE]
        )
    }

    if (is.null(ident)) {
        ident <- "Identity"
        data$Identity <- Idents(object)
    }

    if (plot_type == "dim" && !is.null(graph)) {
        if (!graph %in% Graphs(object)) {
            stop("The object does not have graph:", graph)
        }
        graph <- object@graphs[[graph]]
    }

    .feature_stat_plot(
        data = data, features = features, plot_type = plot_type, pos_only = pos_only,
        should_shrink = should_shrink, should_pivot = should_pivot,
        graph = graph, bg_cutoff = bg_cutoff, downsample = downsample,
        dims = dims, rows_name = rows_name, ident = ident, agg = agg,
        group_by = group_by, split_by = split_by, facet_by = facet_by,
        xlab = xlab, ylab = ylab, x_text_angle = x_text_angle, ...
    )
}

#' @export
FeatureStatPlot.character <- function(
    object, features, plot_type = c("violin", "box", "bar", "ridge", "dim", "cor", "heatmap", "dot"),
    spat_unit = NULL, feat_type = NULL, downsample = NULL, pos_only = c("no", "any", "all"),
    reduction = NULL, graph = NULL, bg_cutoff = 0, dims = 1:2, rows_name = "Features",
    ident = NULL, assay = NULL, layer = NULL, agg = mean, group_by = NULL,
    split_by = NULL, facet_by = NULL, xlab = NULL, ylab = NULL, x_text_angle = NULL, ...
) {
    if (!endsWith(object, ".h5ad")) {
        stop("[FeatureStatPlot] Currently only supports .h5ad files when called with a string/path.")
    }

    object <- hdf5r::H5File$new(object, mode = "r")
    on.exit(object$close_all())

    FeatureStatPlot.H5File(
        object, features = features, plot_type = plot_type,
        spat_unit = spat_unit, feat_type = feat_type, downsample = downsample,
        reduction = reduction, graph = graph, bg_cutoff = bg_cutoff,
        dims = dims, rows_name = rows_name, ident = ident,
        assay = assay, layer = layer, agg = agg,
        group_by = group_by, split_by = split_by, facet_by = facet_by,
        xlab = xlab, ylab = ylab, x_text_angle = x_text_angle, ...
    )
}

#' @export
FeatureStatPlot.H5File <- function(
    object, features, plot_type = c("violin", "box", "bar", "ridge", "dim", "cor", "heatmap", "dot"),
    spat_unit = NULL, feat_type = NULL, downsample = NULL, pos_only = c("no", "any", "all"),
    reduction = NULL, graph = NULL, bg_cutoff = 0, dims = 1:2, rows_name = "Features",
    ident = NULL, assay = NULL, layer = NULL, agg = mean, group_by = NULL,
    split_by = NULL, facet_by = NULL, xlab = NULL, ylab = NULL, x_text_angle = NULL, ...
) {
    plot_type <- match.arg(plot_type)
    if (!is.null(facet_by) && plot_type != "dim") {
        stop("Cannot facet plots because the plots are facetted by the 'features'.")
    }
    stopifnot("[FeatureStatPlot] 'spat_unit' and 'feat_type' should not be used for H5File object." = is.null(spat_unit) && is.null(feat_type))

    # dim plot may use expression for highlighting cells
    # Heatmap may use other variables as annotations, but shrinking only includes minimal columns
    should_shrink <- !plot_type %in% c("dim", "heatmap", "dot")
    should_pivot <- !plot_type %in% c("dim", "heatmap", "dot", "cor")
    if (plot_type != "dim" && is.null(ident)) {
        stop("[FeatureStatPlot] 'ident' should be provided for non-dim plot for H5File object.")
    }

    unlisted_features <- unname(unlist(features))
    assay <- assay %||% "RNA"
    assay_key <- ifelse(identical(assay, "RNA"), "X", paste0("layers/", assay))
    assay_data <- h5group_to_matrix(object[[assay_key]])

    feature_key <- hdf5r::h5attr(object[["var"]], "_index")  # index
    assay_features <- object[["var"]][[feature_key]]$read()
    rownames(assay_data) <- assay_features
    assay_features <- intersect(unlisted_features, assay_features)
    assay_data <- t(as.matrix(assay_data[assay_features, , drop = FALSE]))

    if (is.null(reduction) && plot_type != "dim") {
        data <- cbind(h5group_to_dataframe(object[["obs"]]), assay_data)
    } else {
        reductions <- names(object[['obsm']])
        if (is.null(reductions) || length(reductions) == 0) {
            stop("[FeatureStatPlot] The object does not have any reductions.")
        }

        reduction <- reduction %||% reductions[1]
        if (!startsWith(reduction, "X_") && !reduction %in% reductions) {
            reduction <- paste0("X_", reduction)
        }
        if (!reduction %in% reductions) {
            stop("[FeatureStatPlot] The object does not have reduction:", reduction)
        }

        data <- cbind(
            t(object[["obsm"]][[reduction]]$read()),
            h5group_to_dataframe(object[["obs"]]),
            assay_data
        )
    }

    if (plot_type == "dim" && !is.null(graph)) {
        ggrp <- object[['obsp']][['connectivities']]
        if (is.null(ggrp)) {
            stop("[FeatureStatPlot] The object does not have any graph/connectivities.")
        }

        graph <- h5group_to_matrix(ggrp)
        graph <- igraph::graph_from_adjacency_matrix(graph, mode = "undirected", weighted = TRUE)
        graph <- igraph::as_adjacency_matrix(graph, attr = "weight", sparse = TRUE)
        colnames(graph) <- rownames(graph) <- object[['obs']][['index']]$read()
        rownames(data) <- rownames(graph)
    }

    .feature_stat_plot(
        data = data, features = features, plot_type = plot_type, pos_only = pos_only,
        should_shrink = should_shrink, should_pivot = should_pivot,
        graph = graph, bg_cutoff = bg_cutoff, downsample = downsample,
        dims = dims, rows_name = rows_name, ident = ident, agg = agg,
        group_by = group_by, split_by = split_by, facet_by = facet_by,
        xlab = xlab, ylab = ylab, x_text_angle = x_text_angle, ...
    )
}
