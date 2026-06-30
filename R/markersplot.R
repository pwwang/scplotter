#' Visualize differential expression markers
#'
#' @description
#' Visualize differential expression (DE) results — typically the output of
#' \code{\link[Seurat:FindMarkers]{Seurat::FindMarkers()}} or
#' \code{\link[Seurat:FindAllMarkers]{Seurat::FindAllMarkers()}} — across a
#' variety of plot types. \code{MarkersPlot()} bridges the gap between DE
#' testing and visualization by providing a unified interface for both
#' \strong{summary-level DE visualizations} (volcano, jitter, heatmap, and dot
#' plots of fold changes and significance) and \strong{expression-level
#' visualizations} (violin, box, bar, ridge, heatmap, and dot plots of actual
#' expression values from a Seurat object).
#'
#' The function handles two broad categories of plots:
#' \itemize{
#'   \item \strong{DE summary plots} (no \code{object} required): visualize the
#'     DE statistics themselves — log2 fold change, percentage difference,
#'     p-values, and adjusted p-values — across groups or comparisons.
#'     \itemize{
#'       \item \code{"volcano"} / \code{"volcano_log2fc"} — Volcano plot with
#'         log2 fold change on the x-axis and \eqn{-log_{10}(p)} on the y-axis.
#'         Genes passing the \code{cutoff} are highlighted and top genes are
#'         labeled. Ideal for overview of effect size vs. significance.
#'       \item \code{"volcano_pct"} — Volcano plot with percentage-point
#'         difference (\code{pct.1 - pct.2}) on the x-axis. Useful when the
#'         biological question is about detection rate rather than expression
#'         magnitude.
#'       \item \code{"jitter"} / \code{"jitter_log2fc"} — Jitter plot of log2
#'         fold changes across groups (defined by \code{subset_by}). Dot size
#'         encodes \eqn{-log_{10}(p)}. Reveals distribution of effect sizes
#'         per cluster or condition.
#'       \item \code{"jitter_pct"} — Jitter plot of percentage-point
#'         differences across groups.
#'       \item \code{"heatmap_log2fc"} — Heatmap of log2 fold changes (genes
#'         × groups). Cells can be marked for significance via \code{cutoff}
#'         and \code{sig_mark}.
#'       \item \code{"heatmap_pct"} — Heatmap of percentage-point differences
#'         (genes × groups). Same significance-marking support.
#'       \item \code{"dot_log2fc"} — Dot plot of log2 fold changes (genes ×
#'         groups). Dot size encodes \eqn{-log_{10}(p)}.
#'       \item \code{"dot_pct"} — Dot plot of percentage-point differences
#'         (genes × groups). Dot size encodes \eqn{-log_{10}(p)}.
#'     }
#'   \item \strong{Expression plots} (\code{object} required): visualize the
#'     actual expression values of the selected marker genes in the context of
#'     the original Seurat object. These are useful for validating DE results
#'     by inspecting the underlying expression distributions.
#'     \itemize{
#'       \item \code{"heatmap"} — Expression heatmap of selected marker genes.
#'       \item \code{"violin"} — Violin plots of expression per gene.
#'       \item \code{"box"} — Box plots of expression per gene.
#'       \item \code{"bar"} — Bar plots of mean expression per gene.
#'       \item \code{"ridge"} — Ridge plots of expression distribution per gene.
#'       \item \code{"dot"} — Dot plot of expression (fraction expressing ×
#'         mean expression) per gene.
#'     }
#' }
#'
#' @section Metadata column mapping:
#' Both \code{subset_by} and \code{comparison_by} support a special
#' \code{"marker_column:metadata_column"} syntax for linking columns in the
#' markers data frame to columns in the Seurat object's metadata.
#' \itemize{
#'   \item The part before the colon refers to a column in \code{markers}.
#'   \item The part after the colon refers to a column in
#'     \code{object@meta.data}.
#'   \item If only one name is provided (no colon), it is used for both the
#'     markers column and the metadata column (if a matching metadata column
#'     exists).
#'   \item Example: \code{subset_by = "cluster:RNA_snn_res.0.8"} maps the
#'     \code{cluster} column in the DE results to the
#'     \code{RNA_snn_res.0.8} column in the Seurat metadata.
#' }
#' When the markers data frame and object metadata are merged via
#' \code{subset_by}, only the first value of each non-key column within each
#' group is retained — this is by design to avoid duplication.
#'
#' @section Marker selection and filtering:
#' The \code{select} argument supports three modes:
#' \itemize{
#'   \item \strong{Numeric} — Select the top \code{N} markers (ordered by
#'     \code{order_by}) within each group defined by \code{subset_by}. For
#'     volcano and jitter plots, all markers are plotted but only the top
#'     \code{N} per group are labeled. For other plot types, only the selected
#'     markers are shown.
#'   \item \strong{Single expression} — A filter expression string evaluated by
#'     \code{\link[dplyr:filter]{dplyr::filter()}}. For example,
#'     \code{"p_val_adj < 0.05 & avg_log2FC > 1"}. All markers matching the
#'     condition are retained across all groups.
#'   \item \strong{Multiple expressions} (character vector) — Each element is
#'     evaluated independently. Expressions that mention the \code{subset_by}
#'     column filter the overall data (removing groups); other expressions
#'     filter within the remaining data. For example,
#'     \code{select = c("cluster \%in\% c('0', '1')", "p_val_adj < 0.05")}
#'     first restricts to clusters 0 and 1, then keeps only significant
#'     markers. A numeric string like \code{"5"} among the expressions is
#'     treated as a top-N selection.
#' }
#'
#' Default \code{select}: \code{5} for volcano and jitter plot types,
#' \code{10} for all other plot types.
#'
#' @section Significance marking in heatmaps:
#' For \code{heatmap_log2fc} and \code{heatmap_pct}, the \code{cutoff} and
#' \code{sig_mark} arguments control how statistically significant cells are
#' annotated in the heatmap:
#' \itemize{
#'   \item When \code{cutoff} is set and \code{show_labels = FALSE}, cells
#'     with p-value (or adjusted p-value) below the cutoff are marked with
#'     \code{sig_mark} using ComplexHeatmap's mark system. Valid \code{sig_mark}
#'     values include \code{"-"}, \code{"|"}, \code{"+"}, \code{"/"},
#'     \code{"\\\\"}, \code{"x"}, \code{"o"}, and compound marks like
#'     \code{"[*]"}, \code{"<*>"}, \code{"(*)"}, \code{"{*}"}.
#'   \item When \code{cutoff} is set and \code{show_labels = TRUE}, both
#'     numeric values and significance marks are displayed
#'     (\code{cell_type = "label+mark"}). Note that \code{sig_mark = "*"}
#'     does not work with \code{show_labels = TRUE} — use compound marks
#'     instead.
#'   \item When \code{cutoff = NULL} and \code{show_labels = TRUE}, all cells
#'     are labeled with their numeric values.
#' }
#'
#' @param markers A data frame of differential expression results, typically
#'   the output of \code{\link[Seurat:FindMarkers]{Seurat::FindMarkers()}} or
#'   \code{\link[Seurat:FindAllMarkers]{Seurat::FindAllMarkers()}}. Must
#'   contain columns \code{"gene"} (or gene symbols as rownames),
#'   \code{"p_val"}, and \code{"avg_log2FC"}. For percentage-based plots
#'   (\code{volcano_pct}, \code{jitter_pct}, \code{heatmap_pct},
#'   \code{dot_pct}), columns \code{"pct.1"} and \code{"pct.2"} are also
#'   required.
#' @param object A Seurat object. Required for expression-based plot types:
#'   \code{"heatmap"}, \code{"violin"}, \code{"box"}, \code{"bar"},
#'   \code{"ridge"}, and \code{"dot"}. Not used for DE summary plot types.
#'   Default: \code{NULL}.
#' @param plot_type The type of plot to generate. One of \code{"volcano"},
#'   \code{"volcano_log2fc"}, \code{"volcano_pct"}, \code{"jitter"},
#'   \code{"jitter_log2fc"}, \code{"jitter_pct"}, \code{"heatmap_log2fc"},
#'   \code{"heatmap_pct"}, \code{"dot_log2fc"}, \code{"dot_pct"},
#'   \code{"heatmap"}, \code{"violin"}, \code{"box"}, \code{"bar"},
#'   \code{"ridge"}, or \code{"dot"}. See Description for details on each type.
#' @param subset_by A column name in \code{markers} indicating the grouping
#'   from which each marker was identified (e.g., the \code{cluster} column
#'   from \code{FindAllMarkers()}). Supports the
#'   \code{"marker_column:metadata_column"} syntax for linking to Seurat
#'   object metadata (see \strong{Metadata column mapping} section). For
#'   jitter and DE heatmap/dot plot types, \code{subset_by} is required and
#'   defines the x-axis or column groups. For expression plot types,
#'   \code{subset_by} controls faceting or splitting. Default: \code{NULL}.
#' @param subset_as_facet Logical. If \code{TRUE}, facet the plot by
#'   \code{subset_by} groups instead of splitting into separate plots. Most
#'   useful for expression plot types. For volcano plots, controls whether
#'   faceting or split_by dispatch is used. Default: \code{FALSE}.
#' @param comparison_by A column name in \code{markers} indicating the
#'   comparison (e.g., \code{"g1:g2"} for a pairwise comparison, or a single
#'   group name for one-vs-rest). Required for expression-based plot types
#'   (\code{"heatmap"}, \code{"violin"}, \code{"box"}, \code{"bar"},
#'   \code{"ridge"}, \code{"dot"}). Supports the
#'   \code{"marker_column:metadata_column"} syntax (see \strong{Metadata
#'   column mapping} section). If the comparison values contain a colon
#'   (e.g., \code{"G2M:G1"}), the two groups on either side of the colon
#'   are used to subset the object. If only a single group is present, a
#'   one-vs-other comparison is assumed. Default: \code{NULL}.
#' @param p_adjust Logical. If \code{TRUE} (default), use adjusted p-value
#'   (\code{p_val_adj} column) for significance calculations and y-axis
#'   transformations. If \code{FALSE}, use raw p-value (\code{p_val} column).
#' @param cutoff Numeric. The p-value (or adjusted p-value, depending on
#'   \code{p_adjust}) threshold for labeling significance. For volcano plots,
#'   sets \code{y_cutoff}. For heatmap-based DE plots
#'   (\code{heatmap_log2fc}, \code{heatmap_pct}), controls which cells
#'   receive significance marks. Default: \code{NULL} (no cutoff; defaults
#'   to \code{0.05} for volcano plots).
#' @param show_labels Logical. For \code{heatmap_log2fc} and
#'   \code{heatmap_pct} plot types only. If \code{TRUE}, display numeric
#'   values in heatmap cells. When combined with \code{cutoff}, both values
#'   and significance marks are shown. Default: \code{FALSE}.
#' @param sig_mark Character. The symbol or compound mark used to annotate
#'   statistically significant cells in \code{heatmap_log2fc} and
#'   \code{heatmap_pct} plots. Must be a valid ComplexHeatmap mark: single
#'   characters (\code{"-"}, \code{"|"}, \code{"+"}, \code{"/"},
#'   \code{"\\\\"}, \code{"x"}, \code{"o"}) or compound marks
#'   (\code{"[*]"}, \code{"<*>"}, \code{"(*)"}, \code{"{*}"}). Note that
#'   \code{"*"} conflicts with \code{show_labels = TRUE} because both use
#'   the label layer — use a compound mark instead. Default: \code{"*"}.
#' @param order_by A string expression to order markers by (evaluated with
#'   \code{\link[dplyr:arrange]{dplyr::arrange()}}). Can reference columns
#'   in \code{markers} as well as columns from the object metadata (when
#'   \code{object} is provided and \code{subset_by} enables merging). Only
#'   the first value of merged metadata columns is used. Example:
#'   \code{"desc(avg_log2FC)"}. The ordering affects which markers are
#'   selected when \code{select} is numeric. Default: \code{NULL}.
#' @param select How to select markers for labeling or display. See
#'   \strong{Marker selection and filtering} section for full details.
#'   \itemize{
#'     \item Numeric: Top N markers per \code{subset_by} group (default:
#'       \code{5} for volcano/jitter types, \code{10} for others).
#'     \item Character expression: Filter condition for
#'       \code{\link[dplyr:filter]{dplyr::filter()}}.
#'     \item Character vector: Multiple filter expressions; those containing
#'       the \code{subset_by} column name filter the overall data, others
#'       filter within remaining data.
#'   }
#' @param ... Additional arguments passed to the underlying plotting
#'   function, depending on \code{plot_type}:
#'   \describe{
#'     \item{For \code{volcano}, \code{volcano_log2fc}, \code{volcano_pct}}{
#'       Passed to \code{\link[plotthis:VolcanoPlot]{plotthis::VolcanoPlot()}}.
#'       Common arguments: \code{x_cutoff}, \code{x_cutoff_name},
#'       \code{label_by}, \code{color_by}, \code{nlabel}, \code{flip_negative}.
#'     }
#'     \item{For \code{jitter}, \code{jitter_log2fc}, \code{jitter_pct}}{
#'       Passed to \code{\link[plotthis:JitterPlot]{plotthis::JitterPlot()}}.
#'       Common arguments: \code{add_hline}, \code{shape}, \code{size_by},
#'       \code{nlabel}.
#'     }
#'     \item{For \code{heatmap_log2fc}, \code{heatmap_pct}, \code{dot_log2fc},
#'       \code{dot_pct}}{
#'       Passed to \code{\link[plotthis:Heatmap]{plotthis::Heatmap()}}.
#'       Common arguments: \code{show_row_names}, \code{show_column_names},
#'       \code{values_fill}, \code{palette}, \code{cluster_rows},
#'       \code{cluster_columns}, \code{add_reticle}.
#'     }
#'     \item{For \code{heatmap}, \code{violin}, \code{box}, \code{bar},
#'       \code{ridge}, \code{dot}}{
#'       Passed to \code{\link{FeatureStatPlot}}. Common arguments:
#'       \code{name}, \code{palette}, \code{ncol}, \code{nrow},
#'       \code{stack}, \code{columns_split_by}.
#'     }
#'   }
#' @return A ggplot object (from \code{\link[plotthis:VolcanoPlot]{plotthis::VolcanoPlot()}}
#'   or \code{\link[plotthis:JitterPlot]{plotthis::JitterPlot()}}), a
#'   Heatmap object (from \code{\link[plotthis:Heatmap]{plotthis::Heatmap()}}),
#'   or a ggplot/patchwork object (from \code{\link{FeatureStatPlot}}). When
#'   \code{split_by} or faceting generates multiple plots and
#'   \code{combine = TRUE} (default), a combined patchwork object is
#'   returned; when \code{combine = FALSE}, a list of individual plots is
#'   returned.
#' @note
#' \itemize{
#'   \item \code{subset_by} is required for jitter plots
#'     (\code{"jitter"}, \code{"jitter_log2fc"}, \code{"jitter_pct"}) and
#'     DE heatmap/dot plots (\code{"heatmap_log2fc"}, \code{"heatmap_pct"},
#'     \code{"dot_log2fc"}, \code{"dot_pct"}). Without it, there is no
#'     grouping axis.
#'   \item \code{comparison_by} is required for expression-based plot types
#'     (\code{"heatmap"}, \code{"violin"}, \code{"box"}, \code{"bar"},
#'     \code{"ridge"}, \code{"dot"}) — it tells the function which
#'     comparison groups to extract from the object.
#'   \item When \code{object} is provided and \code{subset_by} maps to a
#'     metadata column, the markers data frame is left-joined with the object
#'     metadata. Only the first row per group is kept for non-key columns,
#'     which is sufficient for most annotation purposes but can cause issues
#'     if per-cell metadata is needed.
#'   \item For expression-based heatmap and dot plots, when
#'     \code{subset_by_2} is available (i.e., the metadata column is mapped),
#'     genes are automatically grouped by \code{subset_by} via
#'     \code{columns_split_by}, and \code{group_by} is set to \code{NULL}.
#'   \item The function calculates \eqn{-log_{10}(p)} (or
#'     \eqn{-log_{10}(p_{adj})}) internally and stores it in a temporary
#'     \code{neg_log10_p} column. This column is available for use in
#'     \code{order_by}.
#'   \item When the comparison involves only a single group (one-vs-rest),
#'     cells not in the comparison group are labeled \code{"Other"} in the
#'     object metadata.
#' }
#' @seealso
#' \code{\link[plotthis:VolcanoPlot]{plotthis::VolcanoPlot()}},
#' \code{\link[plotthis:JitterPlot]{plotthis::JitterPlot()}},
#' \code{\link[plotthis:Heatmap]{plotthis::Heatmap()}},
#' \code{\link{FeatureStatPlot}},
#' \code{\link[Seurat:FindMarkers]{Seurat::FindMarkers()}},
#' \code{\link[Seurat:FindAllMarkers]{Seurat::FindAllMarkers()}}
#' @examples
#' \donttest{
#' data(pancreas_sub)
#' markers <- Seurat::FindMarkers(pancreas_sub,
#'  group.by = "Phase", ident.1 = "G2M", ident.2 = "G1")
#' allmarkers <- Seurat::FindAllMarkers(pancreas_sub)  # seurat_clusters
#'
#' MarkersPlot(markers)
#' MarkersPlot(markers, x_cutoff = 2)
#' MarkersPlot(allmarkers,
#'     subset_by = "cluster", ncol = 2, subset_as_facet = TRUE)
#' MarkersPlot(markers, plot_type = "volcano_pct", flip_negative = TRUE)
#'
#' MarkersPlot(allmarkers, plot_type = "jitter", subset_by = "cluster")
#' MarkersPlot(allmarkers, plot_type = "jitter_pct",
#'     subset_by = "cluster", add_hline = 0, shape = 16)
#'
#' MarkersPlot(allmarkers, plot_type = "heatmap_log2fc", subset_by = "cluster")
#' MarkersPlot(allmarkers, plot_type = "heatmap_log2fc", subset_by = "cluster",
#'     label = scales::label_number(accuracy = 0.01),
#'     cutoff = 0.05, show_labels = TRUE, sig_mark = '{}')
#' MarkersPlot(allmarkers, plot_type = "heatmap_pct", subset_by = "cluster",
#'     cutoff = 0.05)
#'
#' MarkersPlot(allmarkers, plot_type = "dot_log2fc", subset_by = "cluster",
#'     add_reticle = TRUE)
#'
#' MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "heatmap",
#'    columns_split_by = "CellType",
#'    comparison_by = "cluster:seurat_clusters")
#'
#' # Suppose we did a DE between g1 and g2 in each cluster
#' allmarkers$comparison <- "g1:g2"
#' MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "heatmap",
#'    comparison_by = "Phase", subset_by = "cluster:seurat_clusters")
#' MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "dot",
#'    comparison_by = "Phase", subset_by = "cluster:seurat_clusters")
#'
#' MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "violin", select = 3,
#'    comparison_by = "Phase", subset_by = "cluster:seurat_clusters")
#'
#' # select markers with a custom condition, e.g.,
#' # significant markers in cluster 0, 1, and 2 with pct.2 - pct.1 > 0.6
#' # Note that other clusters are still included in the plot
#' MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "violin", subset_by = "cluster",
#'    select = c('cluster %in% c("1", "2", "0") & pct.2 - pct.1 > 0.6'),
#'    comparison_by = "cluster:seurat_clusters",
#'    cutoff = 0.05)
#'
#' # To exclude other clusters, you can separate the filtering conditions into
#' # multiple expressions
#' MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "violin", subset_by = "cluster",
#'    select = c('cluster %in% c("1", "2", "0")', 'pct.2 - pct.1 > 0.6'),
#'    comparison_by = "cluster:seurat_clusters",
#'    cutoff = 0.05)
#'
#' MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "box", select = 3,
#'    comparison_by = "Phase", subset_by = "cluster:seurat_clusters")
#'
#' MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "ridge", select = 2,
#'    comparison_by = "Phase", subset_by = "cluster:seurat_clusters",
#'    ncol = 2)
#' }
#' @export
MarkersPlot <- function(
    markers, object = NULL,
    plot_type = c(
        # Plots that don't need object
        "volcano", "volcano_log2fc", "volcano_pct",
        "jitter", "jitter_log2fc", "jitter_pct",
        "heatmap_log2fc", "heatmap_pct", "dot_log2fc", "dot_pct",
        # Plots that need object, basically expression values, but
        # markers are selected from the markers data frame
        "heatmap", "violin", "box", "bar", "ridge", "dot"
    ),
    subset_by = NULL,
    subset_as_facet = FALSE,
    comparison_by = NULL,
    p_adjust = TRUE,
    cutoff = NULL,
    show_labels = FALSE,
    sig_mark = "*",
    order_by = NULL,
    select = ifelse(plot_type %in% c(
        "volcano", "volcano_log2fc", "volcano_pct",
        "jitter", "jitter_log2fc", "jitter_pct"
    ), 5, 10),
    ...
) {
    plot_type <- match.arg(plot_type)

    # check if object is provided for plot types that need it
    plot_types_need_object <- c(
        "heatmap", "violin", "box", "bar", "ridge", "dot"
    )
    if (plot_type %in% plot_types_need_object & is.null(object)) {
        stop("[MarkersPlot] `object` is required for plot_type '", plot_type, "'")
    }

    # Add the gene column if missing (likely from FindMarkers)
    if (!"gene" %in% colnames(markers)) {
        markers$gene <- rownames(markers)
    }
    # Check subset_by
    check_columns <- utils::getFromNamespace("check_columns", "plotthis")
    if (!is.null(subset_by)) {
        subset_by_1 <- strsplit(subset_by, ":")[[1]][1]
        subset_by_2 <- strsplit(subset_by, ":")[[1]][2]  # NA if not provided
        subset_by_2 <- if (!is.na(subset_by_2)) subset_by_2 else NULL
        if (is.null(subset_by_2) && !is.null(object) && subset_by_1 %in% colnames(object@meta.data)) {
            subset_by_2 <- subset_by_1
        }
    } else {
        subset_by_1 <- NULL
        subset_by_2 <- NULL
    }
    subset_by_1 <- check_columns(markers, subset_by_1)
    if (!is.null(subset_by_1) && !is.null(subset_by_2) && !is.null(object)) {
        if (!subset_by_2 %in% colnames(object@meta.data)) {
            stop("[MarkersPlot] `subset_by` '", subset_by_2, "' is not found in the object's metadata.")
        }
        # check if subset_by values are consistent between markers and object
        sub_markers <- unique(markers[[subset_by_1]])
        sub_object <- unique(object@meta.data[[subset_by_2]])
        nonexisting_sub <- setdiff(sub_markers, sub_object)
        if (length(nonexisting_sub) > 0) {
            stop(
                "[MarkersPlot] The following values in `subset_by` '", subset_by_1,
                "' are not found in the object's metadata (", subset_by_2, "): ",
                paste(nonexisting_sub, collapse = ", "))
        }
        meta <- dplyr::summarise(object@meta.data, dplyr::across(dplyr::everything(), ~ .[1]), .by = !!rlang::sym(subset_by_2))
        markers <- dplyr::left_join(markers, meta, by = stats::setNames(subset_by_2, subset_by_1), suffix = c("", ".meta"))
    }

    # Check comparison_by
    if (!is.null(comparison_by)) {
        comparison_by_1 <- strsplit(comparison_by, ":")[[1]][1]
        comparison_by_2 <- strsplit(comparison_by, ":")[[1]][2]  # NA if not provided
        comparison_by_2 <- if (!is.na(comparison_by_2)) comparison_by_2 else NULL
        if (is.null(comparison_by_2) && !is.null(object) && comparison_by_1 %in% colnames(object@meta.data)) {
            comparison_by_2 <- comparison_by_1
        }
    } else {
        comparison_by_1 <- NULL
        comparison_by_2 <- NULL
    }
    comparison_by_1 <- check_columns(markers, comparison_by_1)
    if (!is.null(comparison_by_2) && !is.null(object)) {
        if (!comparison_by_2 %in% colnames(object@meta.data)) {
            stop("[MarkersPlot] `comparison_by` '", comparison_by_2, "' is not found in the object's metadata.")
        }
        # check if comparison_by values are consistent between markers and object
        comp_markers <- unique(unlist(strsplit(unique(as.character(markers[[comparison_by_1]])), ":")))
        comp_object <- unique(object@meta.data[[comparison_by_2]])
        nonexisting_comp <- setdiff(comp_markers, comp_object)
        if (length(nonexisting_comp) > 0) {
            stop(
                "[MarkersPlot] The following values in `comparison_by` '", comparison_by_1,
                "' are not found in the object's metadata (", comparison_by_2, "): ",
                paste(nonexisting_comp, collapse = ", "))
        }
    }

    pcol <- ifelse(p_adjust, "p_val_adj", "p_val")

    # calculate pct.1 - pct.2 if needed
    if (plot_type %in% c("volcano", "volcano_log2fc", "volcano_pct", "jitter_pct", "heatmap_pct", "dot_pct")) {
        if (!all(c("pct.1", "pct.2") %in% colnames(markers))) {
            stop("[MarkersPlot] `markers` must contain 'pct.1' and 'pct.2' columns for plot_type '", plot_type, "'")
        }
        markers <- dplyr::mutate(markers, pct_diff = !!sym("pct.1") - !!sym("pct.2"))
    }

    # calcualte -log10(p-value) or -log10(adjusted p-value)
    markers <- dplyr::mutate(markers, neg_log10_p = -log10(!!rlang::sym(pcol)))

    # order markers by order_by
    if (!is.null(order_by)) {
        markers <- dplyr::arrange(markers, !!rlang::parse_expr(order_by))
    }

    if (plot_type %in% c("volcano", "volcano_log2fc", "volcano_pct")) {
        args <- list(
            markers,
            x = ifelse(plot_type == "volcano_pct", "pct_diff", "avg_log2FC"),
            y = pcol,
            y_cutoff = cutoff %||% 0.05,
            y_cutoff_name = paste0(pcol, " = ", cutoff %||% 0.05),
            label_by = "gene",
            ...
        )
        args$color_by <- args$color_by %||% ifelse(plot_type == "volcano_pct", "avg_log2FC", "pct_diff")
        if (!is.null(subset_by_1) && subset_as_facet) {
            args$facet_by <- subset_by_1
        } else if (!is.null(subset_by_1)) {
            args$split_by <- subset_by_1
        }
        do_call(plotthis::VolcanoPlot, args)
    } else if (plot_type %in% c("jitter", "jitter_log2fc", "jitter_pct")) {
        if (is.null(subset_by_1)) {
            stop("[MarkersPlot] `subset_by` is required for plot_type '", plot_type, "'. Consider using volcano plot if you don't have groups.")
        }
        if (!is.numeric(select)) {
            stop("[MarkersPlot] `select` must be numeric for plot_type '", plot_type, "', to label top N markers in each group.")
        }
        args <- list(
            markers,
            x = subset_by_1,
            y = ifelse(plot_type == "jitter_pct", "pct_diff", "avg_log2FC"),
            size_by = "neg_log10_p",
            size_name = paste0("-log10(", pcol, ")"),
            label_by = "gene",
            nlabel = select,
            ...
        )
        if (!is.null(order_by)) {
            args$order_by <- order_by
        }
        do_call(plotthis::JitterPlot, args)
    } else if (plot_type %in% c("heatmap_log2fc", "heatmap_pct", "dot_log2fc", "dot_pct")) {
        if (is.null(subset_by)) {
            stop("[MarkersPlot] `subset_by` is required for plot_type '", plot_type, "'")
        }
        y <- ifelse(endsWith(plot_type, "_pct"), "pct_diff", "avg_log2FC")
        y_max <- max(markers[[y]], na.rm = TRUE)
        y_min <- min(markers[[y]], na.rm = TRUE)
        if (y_max > 0 && y_min < 0) {
            # center the color bar at 0
            y_max <- max(abs(y_max), abs(y_min))
            y_min <- -y_max
        }
        if (is.numeric(select)) {
            genes <- dplyr::slice_head(markers, n = select, by = !!rlang::sym(subset_by_1))$gene
        } else if (length(select) == 1) {
            genes <- dplyr::filter(markers, !!rlang::parse_expr(select))$gene
        } else {
            # The expressions in select with the entire `subset_by` word in it are
            # supposed to be the ones to filter the data
            select_sb <- grepl(paste0("\\b", subset_by_1, "\\b"), select)
            if (any(select_sb)) {
                markers <- dplyr::filter(markers, !!!rlang::parse_exprs(select[select_sb]))
                if (all(select_sb)) {
                    genes <- markers$gene
                } else {
                    select_non_sb <- select[!select_sb]
                    if (length(select_non_sb) == 1 && grepl("^\\d+$", select_non_sb)) {
                        genes <- dplyr::slice_head(markers, n = as.numeric(select_non_sb), by = !!rlang::sym(subset_by_1))$gene
                    } else {
                        genes <- dplyr::filter(markers, !!!rlang::parse_exprs(select[!select_sb]))$gene
                    }
                }
            } else {
                genes <- dplyr::filter(markers, !!!rlang::parse_exprs(select))$gene
            }
        }
        genes <- unique(genes)
        markers <- dplyr::filter(markers, !!sym("gene") %in% genes)
        if (!is.factor(markers$gene)) {
            markers$gene <- factor(markers$gene, levels = unique(markers$gene))
        }
        if (!is.factor(markers[[subset_by_1]])) {
            markers[[subset_by_1]] <- factor(markers[[subset_by_1]], levels = unique(markers[[subset_by_1]]))
        }
        args <- list(
            data = markers,
            values_by = y,
            rows_by = "gene",
            columns_by = subset_by_1,
            in_form = "long",
            upper_cutoff = y_max,
            lower_cutoff = y_min,
            ...
        )
        args$show_row_names <- args$show_row_names %||% TRUE
        args$show_column_names <- args$show_column_names %||% TRUE
        args$values_fill <- args$values_fill %||% 0

        # add label if cutoff is provided for heatmap
        genes <- levels(markers$gene)
        groups <- levels(markers[[subset_by_1]])
        if (!is.null(cutoff) && startsWith(plot_type, "heatmap_")) {
            sig_mat <- tidyr::pivot_wider(
                markers,
                id_cols = !!sym("gene"),
                names_from = subset_by_1,
                values_from = !!rlang::sym(pcol),
                values_fill = 1
            )
            sig_mat <- as.data.frame(sig_mat)
            rownames(sig_mat) <- sig_mat$gene
            sig_mat$gene <- NULL
            # There might be some groups failed to run DE analysis
            # So we just exclude them
            failed_groups <- setdiff(groups, colnames(sig_mat))
            if (length(failed_groups) > 1) {
                warning(
                    "[MarkersPlot] The following groups in `subset_by` '", subset_by_1,
                    "' are not found in the markers data frame and will be ignored: ",
                    paste(failed_groups, collapse = ", "),
                    immediate. = TRUE
                )
            }
            groups <- intersect(groups, colnames(sig_mat))
            sig_mat <- sig_mat[genes, groups, drop = FALSE]
            sig_mat <- as.matrix(sig_mat)

            # Using Heatmap' mark
            if (
                sig_mark %in% c('-', '|', '+', '/', '\\', 'x', 'o') |
                startsWith(sig_mark, "[") && endsWith(sig_mark, "]") |
                startsWith(sig_mark, "<") && endsWith(sig_mark, ">") |
                startsWith(sig_mark, "(") && endsWith(sig_mark, ")") |
                startsWith(sig_mark, "{") && endsWith(sig_mark, "}")
            ) {
                pname <- ifelse(p_adjust, "padj", "p")
                if (show_labels) {
                    args$cell_type <- "label+mark"
                    args$mark <- function(x, i, j) {
                        pval <- ComplexHeatmap::pindex(sig_mat, i, j)
                        if (pval < cutoff) list(sig_mark, legend = paste0(pname, " < ", cutoff))
                        else NA
                    }
                } else {
                    args$cell_type <- "mark"
                    args$mark <- function(x, i, j) {
                        pval <- ComplexHeatmap::pindex(sig_mat, i, j)
                        if (pval < cutoff) list(sig_mark, legend = paste0(pname, " < ", cutoff))
                        else NA
                    }
                }
            } else if (show_labels && !(isFALSE(sig_mark) || is.null(sig_mark) || sig_mark == "")) {
                stop(
                    "[MarkersPlot] `cutoff` is provided and `show_labels` is TRUE, ",
                    "`sig_mark` must be a valid `mark` for Heatmap ",
                    "(e.g., '-', '|', '+', '/', '\\', 'x', 'o', '[]', '<>', '()', '{}', or a compound mark like '[*]')"
                )
            } else if (show_labels) {
                args$cell_type <- "label"
            } else {  # arbitrary sig_mark
                args$cell_type <- "label"
                args$label <- function(x, i, j) {
                    pval <- ComplexHeatmap::pindex(sig_mat, i, j)
                    ifelse(pval < cutoff, sig_mark, NA)
                }
            }
        }

        # set dot_size fot dot plot
        if (startsWith(plot_type, "dot_")) {
            ds_mat <- tidyr::pivot_wider(
                markers,
                id_cols = !!sym("gene"),
                names_from = subset_by_1,
                values_from = "neg_log10_p"
            )
            ds_mat <- as.data.frame(ds_mat)
            rownames(ds_mat) <- ds_mat$gene
            ds_mat$gene <- NULL
            # There might be some groups failed to run DE analysis
            # So we just exclude them
            failed_groups <- setdiff(groups, colnames(ds_mat))
            if (length(failed_groups) > 1) {
                warning(
                    "[MarkersPlot] The following groups in `subset_by` '", subset_by_1,
                    "' are not found in the markers data frame and will be ignored: ",
                    paste(failed_groups, collapse = ", "),
                    immediate. = TRUE
                )
            }
            groups <- intersect(groups, colnames(ds_mat))
            ds_mat <- ds_mat[genes, groups, drop = FALSE]
            ds_mat <- as.matrix(ds_mat)

            args$cell_type <- "dot"
            args$dot_size <- function(x, i, j) {
                ComplexHeatmap::pindex(ds_mat, i, j)
            }
            args$dot_size_name <- paste0("-log10(", pcol, ")")
        }
        do_call(plotthis::Heatmap, args)
    } else {  # if (plot_type %in% c("heatmap", "violin", "box", "bar", "ridge", "dot")) {
        if (is.null(comparison_by)) {
            stop("[MarkersPlot] `comparison_by` is required for plot_type '", plot_type, "'")
        }
        if (is.null(comparison_by_2)) {
            comparison_by_2 <- comparison_by_1
            if (!comparison_by_2 %in% colnames(object@meta.data)) {
                stop("[MarkersPlot] `comparison_by` '", comparison_by_2, "' is not found in the object's metadata.")
            }
        }
        if (!is.null(subset_by_1) && is.null(subset_by_2)) {
            if (subset_by_1 %in% colnames(object@meta.data)) {
                subset_by_2 <- subset_by_1
            } else if (is.numeric(select)) {
                warning(
                    "[MarkersPlot] `subset_by` '", subset_by_1, "' is only used to select markers, ",
                    "but not in plotting, since it is not found in the object's metadata. ",
                    "Set `subset_by` to '", subset_by_1, ":<object_metadata_column>' to make it work.",
                    immediate. = TRUE)
            } else {
                warning(
                    "[MarkersPlot] `subset_by` '", subset_by_1, "' is ignored, ",
                    "since it is not found in the object's metadata. ",
                    "Set `subset_by` to '", subset_by_1, ":<object_metadata_column>' to make it work.",
                    immediate. = TRUE)

            }
        }

        if (is.numeric(select)) {
            if (!is.null(subset_by)) {
                genes <- dplyr::slice_head(markers, n = select, by = !!rlang::sym(subset_by_1))
                if (plot_type %in% c("heatmap", "dot")) {
                    genes <- dplyr::summarise(genes, gene = list(!!sym("gene")), .by = !!rlang::sym(subset_by_1))
                    genes <- stats::setNames(genes$gene, genes[[subset_by_1]])
                } else {
                    genes <- genes$gene
                }
            } else {
                genes <- dplyr::slice_head(markers, n = select)$gene
            }
        } else if (is.null(subset_by) || length(select) == 1) {
            genes <- dplyr::filter(markers, !!!rlang::parse_exprs(select))$gene
        } else {
            # The expressions in select with the entire `subset_by` word in it are
            # supposed to be the ones to filter the data
            select_sb <- grepl(paste0("\\b", subset_by_1, "\\b"), select)
            if (any(select_sb)) {
                markers <- dplyr::filter(markers, !!!rlang::parse_exprs(select[select_sb]))
                if (all(select_sb)) {
                    genes <- markers$gene
                } else {
                    select_non_sb <- select[!select_sb]
                    if (length(select_non_sb) == 1 && grepl("^\\d+$", select_non_sb)) {
                        if (is.null(subset_by)) {
                            genes <- dplyr::slice_head(markers, n = as.numeric(select_non_sb))$gene
                        } else {
                            genes <- dplyr::slice_head(markers, n = as.numeric(select_non_sb), by = !!rlang::sym(subset_by_1))$gene
                        }
                    } else {
                        genes <- dplyr::filter(markers, !!!rlang::parse_exprs(select_non_sb))$gene
                    }
                }
            } else {
                genes <- dplyr::filter(markers, !!!rlang::parse_exprs(select))$gene
            }
        }

        # subset the object to only include the comparison groups
        comp_groups <- unique(unlist(strsplit(unique(as.character(markers[[comparison_by_1]])), ":")))
        if (length(comp_groups) > 1) {
            object <- subset_seurat(object, subset = !!rlang::sym(comparison_by_2) %in% comp_groups)
        } else {
            object@meta.data[[comparison_by_2]] <- ifelse(object@meta.data[[comparison_by_2]] == comp_groups, comp_groups, "Other")
            object@meta.data[[comparison_by_2]] <- factor(object@meta.data[[comparison_by_2]], levels = c(comp_groups, "Other"))
        }
        if (!is.list(genes)) {
            unigenes <- unique(genes)
        } else {
            unigenes <- unique(unlist(genes))
        }
        # subset the object to only include the selected genes
        object <- tryCatch({
            # In case the features do not exist in some assays
            subset_seurat(object, features = unigenes)
        }, error = function(e) {
            object
        })

        args <- list(
            object,
            features = genes,
            group_by = comparison_by_2,
            plot_type = plot_type,
            facet_by = if (subset_as_facet) subset_by_2 else NULL,
            split_by = if (!subset_as_facet) subset_by_2 else NULL,
            ident = comparison_by_2,
            ...
        )
        if (plot_type %in% c("heatmap", "dot")) {
            args$name <- args$name %||% "Expression"
            args$cluster_columns <- args$cluster_columns %||% FALSE
            if (!is.null(subset_by_2)) {
                args$facet_by <- NULL
                args$split_by <- NULL
                args$columns_split_by <- subset_by_2
                args$group_by <- NULL
            }
        } else if (plot_type %in% c("violin", "box", "bar")) {
            args$stack <- args$stack %||% TRUE
            if (!is.null(subset_by_2)) {
                args$ident <- subset_by_2
                args$group_by <- comparison_by_2
                args$split_by <- NULL
                args$facet_by <- NULL
            }
        }
        do_call(scplotter::FeatureStatPlot, args)
    }
}
