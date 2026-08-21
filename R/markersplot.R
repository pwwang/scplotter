#' Visualize differential expression markers
#'
#' @description
#' Visualize differential expression (DE) results — typically the output of
#' \code{\link[Seurat:FindMarkers]{Seurat::FindMarkers()}} or
#' \code{\link[Seurat:FindAllMarkers]{Seurat::FindAllMarkers()}} — across a
#' variety of plot types. You can also compose the DE results from other
#' tools into a data frame with the required columns and use this function to visualize them.
#'
#' \code{MarkersPlot()} bridges the gap between DE
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
#'         fold changes across groups (defined by \code{each}). Dot size
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
#' Both \code{each} and \code{group_by} accept a
#' \code{"marker_column:metadata_column"} syntax that links a column in the
#' markers data frame to a column in the Seurat object's metadata.
#' \itemize{
#'   \item The part before the colon must be a column in \code{markers}
#'     (e.g., \code{cluster}) or be empty; the part after the colon must be
#'     a column in \code{object@meta.data} (e.g., \code{seurat_clusters}).
#'   \item This syntax requires \code{object} to be provided; otherwise an
#'     error is raised.
#'   \item When the marker part is non-empty, every value in the marker
#'     column must exist in the metadata column, otherwise an error is
#'     raised.
#'   \item For \code{each} with a non-empty marker part, the metadata is
#'     merged into the markers data frame (keeping the first row of each
#'     metadata group for non-key columns), so metadata columns become
#'     available for arguments like \code{order_by}. On name conflicts, the
#'     merged columns get a \code{.meta} suffix. With an empty marker part
#'     (\code{":metadata_column"}), no merging or per-group selection
#'     happens; the metadata column is used only to split/annotate the
#'     expression plot (\code{columns_split_by} for heatmap/dot, \code{ident}
#'     for violin/box/bar).
#'   \item For \code{group_by}, the object is subset to the cells whose
#'     metadata values occur in the marker column, and the metadata column
#'     is re-factored with those values in the order they first appear in
#'     the marker column. Values separated by a colon (e.g.,
#'     \code{"G1:G2M"}) are split into individual groups.
#' }
#'
#' @section Marker selection and filtering:
#' How \code{select} picks the markers depends on the plot type and the
#' value provided:
#' \itemize{
#'   \item \strong{Numeric} — Select the top \code{N} markers (ordered by
#'     \code{order_by}) within each group defined by \code{each}, or overall
#'     when \code{each} is \code{NULL}. Jitter plots label the top
#'     \code{N} markers per group (a numeric \code{select} is required).
#'     Volcano plots ignore \code{select} — labeling is controlled via
#'     \code{...} (e.g., \code{nlabel}). For expression plot types, a
#'     numeric \code{select} only keeps markers with a p-value below
#'     \code{cutoff} (when set) before the top-N selection.
#'   \item \strong{Single expression} — A filter expression string evaluated by
#'     \code{\link[dplyr:filter]{dplyr::filter()}}. For example,
#'     \code{"p_val_adj < 0.05 & avg_log2FC > 1"}. All markers matching the
#'     condition are retained across all groups.
#'   \item \strong{Multiple expressions} (character vector) — Only for DE
#'     heatmap/dot plot types (\code{"heatmap_log2fc"},
#'     \code{"heatmap_pct"}, \code{"dot_log2fc"}, \code{"dot_pct"}). Each
#'     element is evaluated independently: expressions that mention the
#'     \code{each} column filter the overall data (removing groups); other
#'     expressions filter within the remaining data. For example,
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
#' @param group_by Used only for expression-based plot types (ignored for DE
#'   summary plot types). A column in the Seurat object's metadata to group
#'   cells by, e.g., a condition column — useful when the DEs were calculated
#'   between conditions (such as cell cycle phases) and you want to compare
#'   the expression of the markers across those conditions. A single value is
#'   passed directly to \code{\link{FeatureStatPlot}}: for \code{heatmap} and
#'   \code{dot} plots it is applied as the column annotation (\code{ident}),
#'   and it only takes effect when \code{each} includes a metadata column
#'   mapping; for \code{violin}, \code{box}, \code{bar}, and \code{ridge}
#'   plots it is passed as \code{group_by}.
#'   The \code{"marker_column:metadata_column"} syntax (see \strong{Metadata
#'   column mapping}) restricts the object to only the cells involved in the
#'   comparisons: for example, if a \code{comparison} column in the markers
#'   data frame holds \code{"G1:G2M"}, passing
#'   \code{group_by = "comparison:Phase"} keeps only G1 and G2M cells in the
#'   plot, with the \code{Phase} column re-factored to these two levels in
#'   the order they first appear in the \code{comparison} column. Without
#'   the restriction, e.g., \code{group_by = "Phase"}, all phase cells
#'   (G1, G2M, and S) are included in the plot. Default: \code{NULL}.
#' @param each A column name in \code{markers} indicating the grouping
#'   from which each marker was identified (e.g., the \code{cluster} column
#'   from \code{FindAllMarkers()}). Required for jitter and DE heatmap/dot
#'   plot types, where it defines the x-axis or column groups. For volcano
#'   plot types, it splits the plot by group (or facets it, with
#'   \code{facet_each = TRUE}). For expression plot types, \code{each} is
#'   used to select the markers within each group; a plain column name does
#'   not split the plot — use the \code{"marker_column:metadata_column"}
#'   syntax (see \strong{Metadata column mapping}) to also split the plot by
#'   the mapped metadata column. Alternatively, pass
#'   \code{":metadata_column"} with an empty marker part to split the
#'   expression plot by the metadata column directly, without selecting
#'   markers per group (markers are selected overall) and without merging
#'   metadata. Default: \code{NULL}.
#' @param facet_each Logical. Only for volcano plot types: if \code{TRUE},
#'   facet the volcano plot by the \code{each} groups instead of splitting
#'   it into separate subplots. Ignored for other plot types. Default:
#'   \code{FALSE}.
#' @param p_adjust Logical. If \code{TRUE} (default), use adjusted p-value
#'   (\code{p_val_adj} column) for significance calculations and y-axis
#'   transformations. If \code{FALSE}, use raw p-value (\code{p_val} column).
#' @param cutoff Numeric. The p-value (or adjusted p-value, depending on
#'   \code{p_adjust}) threshold for labeling significance. For volcano plots,
#'   sets \code{y_cutoff}. For DE heatmap plots (\code{heatmap_log2fc},
#'   \code{heatmap_pct}), controls which cells receive significance marks.
#'   For expression plot types with a numeric \code{select}, only markers
#'   with a p-value below \code{cutoff} are eligible for selection. Ignored
#'   by DE dot plots (\code{dot_log2fc}, \code{dot_pct}). Default:
#'   \code{NULL} (no cutoff; defaults to \code{0.05} for volcano plots).
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
#' @param order_by A string of one or more comma-separated expressions used
#'   to order the markers (evaluated with
#'   \code{\link[dplyr:arrange]{dplyr::arrange()}}). Can reference columns
#'   in \code{markers} as well as metadata columns merged in via a
#'   colon-form \code{each} (see \strong{Metadata column mapping}). Only
#'   the first value of each merged metadata column is kept. Example:
#'   \code{"desc(avg_log2FC)"} or \code{"desc(avg_log2FC), desc(pct.1)"}.
#'   The ordering determines which markers are selected when \code{select}
#'   is numeric. For jitter plots, it is also passed to
#'   \code{\link[plotthis:JitterPlot]{plotthis::JitterPlot()}}. Default:
#'   \code{"desc(abs(avg_log2FC))"}.
#' @param select How to select markers for display or labeling. See
#'   \strong{Marker selection and filtering} section for full details.
#'   \itemize{
#'     \item Numeric: Top N markers per \code{each} group, or overall when
#'       \code{each} is \code{NULL} (default: \code{5} for volcano/jitter
#'       types, \code{10} for others).
#'     \item Single expression: Filter condition for
#'       \code{\link[dplyr:filter]{dplyr::filter()}}.
#'     \item Character vector of multiple expressions (DE heatmap/dot plot
#'       types only): expressions mentioning the \code{each} column name
#'       filter the overall data, others filter within the remaining data.
#'   }
#' @param flatten_markers Logical. Only for the expression \code{heatmap} and
#'   \code{dot} plot types. When \code{each} is used to select markers per
#'   group, the markers are by default provided to
#'   \code{\link{FeatureStatPlot}} as a named list (one entry per group),
#'   which splits the feature rows of the plot by group. With
#'   \code{flatten_markers = TRUE}, the selected markers are collapsed into a
#'   single vector so the plot shows one unsplit block of features — useful
#'   e.g. to mimic \code{\link[Seurat:DoHeatmap]{Seurat::DoHeatmap()}} on
#'   globally selected markers. Default: \code{FALSE}.
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
#'       \code{stack}, \code{layer}, \code{cell_type}. Note that
#'       \code{group_by}, \code{ident}, and \code{columns_split_by} are set
#'       by \code{MarkersPlot()} from the \code{group_by} and \code{each}
#'       arguments.
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
#'   \item \code{plot_type} determines which underlying plotting function is called and
#'     also what to be plotted. `volcano`, `volcano_log2fc`, `volcano_pct`
#'     `jitter`, `jitter_log2fc`, `jitter_pct`, `heatmap_log2fc`, `heatmap_pct`, `dot_log2fc`, and `dot_pct`
#'     are DE summary plots that visualize the DE statistics themselves,
#'     while `heatmap`, `violin`, `box`, `bar`, `ridge`, and `dot` are expression-based plots
#'     that visualize the actual expression values of the selected marker genes in the context of the original Seurat object.
#'   \item \code{each} is required for jitter plots
#'     (\code{"jitter"}, \code{"jitter_log2fc"}, \code{"jitter_pct"}) and
#'     DE heatmap/dot plots (\code{"heatmap_log2fc"}, \code{"heatmap_pct"},
#'     \code{"dot_log2fc"}, \code{"dot_pct"}). Its role depends on the plot
#'     type:
#'     \itemize{
#'       \item Volcano plot types: the plot is split by the \code{each}
#'         groups (faceted when \code{facet_each = TRUE}).
#'       \item Jitter plot types: the x-axis grouping.
#'       \item DE heatmap/dot plot types: the columns of the heatmap/dot plot.
#'       \item Expression plot types: used to select the markers within each
#'         group; it does not split or facet the plot. Pass
#'         \code{"marker_column:metadata_column"} (e.g.,
#'         \code{"cluster:seurat_clusters"}) to also split the plot by the
#'         mapped metadata column (via \code{columns_split_by} for
#'         heatmap/dot, or \code{ident} for violin/box/bar), or
#'         \code{":metadata_column"} (e.g., \code{":seurat_clusters"}) to
#'         split the plot by the metadata column without per-group marker
#'         selection.
#'     }
#'   \item When \code{each} uses the
#'     \code{"marker_column:metadata_column"} form with a non-empty marker
#'     column, the markers data frame is left-joined with the object
#'     metadata. Only the first row per group is kept for non-key columns,
#'     which is sufficient for most annotation purposes but can cause issues
#'     if per-cell metadata is needed. The \code{":metadata_column"} form
#'     (empty marker part) skips the join entirely.
#'   \item The function calculates \eqn{-log_{10}(p)} (or
#'     \eqn{-log_{10}(p_{adj})}) internally and stores it in a temporary
#'     \code{neg_log10_p} column. This column is available for use in
#'     \code{order_by}.
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
#' MarkersPlot(allmarkers, each = "cluster", ncol = 2, facet_each = TRUE)
#' MarkersPlot(markers, plot_type = "volcano_pct", flip_negative = TRUE)
#'
#' MarkersPlot(allmarkers, plot_type = "jitter", each = "cluster")
#' MarkersPlot(allmarkers, plot_type = "jitter_pct", order_by = "desc(abs(pct.1 - pct.2))",
#'     each = "cluster", add_hline = 0, shape = 16)
#'
#' MarkersPlot(allmarkers, plot_type = "heatmap_log2fc", each = "cluster",
#'     order_by = "desc(avg_log2FC)", select = 3)
#' MarkersPlot(allmarkers, plot_type = "heatmap_log2fc", each = "cluster",
#'     label = scales::label_number(accuracy = 0.01), select = 3,
#'     cutoff = 0.05, show_labels = TRUE, sig_mark = '{}')
#' MarkersPlot(allmarkers, plot_type = "heatmap_pct", each = "cluster",
#'     cutoff = 0.05, select = 3)
#'
#' MarkersPlot(allmarkers, plot_type = "dot_log2fc", each = "cluster",
#'     add_reticle = TRUE, select = 3)
#'
#' topmarkers <- allmarkers[order(allmarkers$avg_log2FC, decreasing = TRUE), ]
#' # Mimic Seurat's DoHeatmap()
#' MarkersPlot(topmarkers[1:20, ], object = pancreas_sub, plot_type = "heatmap",
#'    layer = "data", cell_type = "bars", flatten_markers = TRUE, cluster_rows = FALSE,
#'    show_column_names = "inplace", each = "cluster:seurat_clusters")
#'
#' # Select top 3 markers per cluster
#' MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "heatmap",
#'    order_by = "desc(avg_log2FC)", select = 3,
#'    layer = "data", cell_type = "bars",
#'    show_column_names = "inplace", each = "cluster:seurat_clusters")
#' # Suppose we did a DE between G2M and G1 phases in each cluster and
#' # stored the results in a new column "comparison"
#' allmarkers$comparison <- "G1:G2M"
#' MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "heatmap",
#'    group_by = "comparison:Phase", each = "cluster:seurat_clusters",
#'    order_by = "desc(avg_log2FC)", select = 3, layer = "data")
#'
#' MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "dot", select = 2,
#'    flatten_markers = TRUE, order_by = "desc(avg_log2FC)",
#'    group_by = "Phase", each = "cluster:seurat_clusters", layer = "data")
#'
#' MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "violin", select = 2,
#'    position_dodge_preserve = "single", add_bg = TRUE, add_box = TRUE,
#'    group_by = "comparison:Phase", each = "cluster:seurat_clusters", layer = "data")
#'
#' # select markers with a custom condition, e.g.,
#' # significant markers in cluster 0, 1, and 2 with pct.2 - pct.1 > 0.6
#' # Note that other clusters are still included in the plot
#' MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "violin",
#'   select = c('cluster %in% c("1", "2", "0") & pct.2 - pct.1 > 0.6'),
#'   each = "cluster:seurat_clusters", cutoff = 0.05, layer = "data")
#'
#' MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "box", select = 3,
#'   group_by = "Phase", each = "cluster:seurat_clusters", layer = "data")
#'
#' MarkersPlot(allmarkers, object = pancreas_sub, plot_type = "ridge", select = 2,
#'    group_by = "Phase", each = "cluster:seurat_clusters", layer = "data",
#'    ncol = 4)
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
    group_by = NULL,
    each = NULL,
    facet_each = FALSE,
    p_adjust = TRUE,
    cutoff = NULL,
    show_labels = FALSE,
    sig_mark = "*",
    order_by = "desc(abs(avg_log2FC))",
    select = ifelse(plot_type %in% c(
        "volcano", "volcano_log2fc", "volcano_pct",
        "jitter", "jitter_log2fc", "jitter_pct"
    ), 5, 10),
    flatten_markers = FALSE,
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
    # Check each
    check_columns <- utils::getFromNamespace("check_columns", "plotthis")
    if (!is.null(each)) {
        if (grepl(":", each)) {
            if (length(strsplit(each, ":")[[1]]) != 2) {
                stop("[MarkersPlot] `each` must be in the format 'marker_column:metadata_column' or a single column name.")
            }
            each_1 <- strsplit(each, ":")[[1]][1]
            each_2 <- strsplit(each, ":")[[1]][2]
            if (!each_2 %in% colnames(object@meta.data)) {
                stop("[MarkersPlot] `each` '", each_2, "' is not found in the object's metadata.")
            }
            if (nchar(each_1) == 0) {
                each_1 <- NULL
            } else {
                # check if each values are consistent between markers and object
                sub_markers <- unique(markers[[each_1]])
                sub_object <- unique(object@meta.data[[each_2]])
                nonexisting_sub <- setdiff(sub_markers, sub_object)
                if (length(nonexisting_sub) > 0) {
                    stop(
                        "[MarkersPlot] The following values in `each` '", each_1,
                        "' are not found in the object's metadata (", each_2, "): ",
                        paste(nonexisting_sub, collapse = ", "))
                }
                # Get the first row of each group in the metadata to avoid duplication
                # User has to make sure that the metadata columns are consistent within each group
                meta <- dplyr::summarise(object@meta.data, dplyr::across(dplyr::everything(), ~ .[1]), .by = !!rlang::sym(each_2))
                markers <- dplyr::left_join(markers, meta, by = stats::setNames(each_2, each_1), suffix = c("", ".meta"))
            }
        } else {
            each_1 <- each
            each_2 <- NULL
        }

        each_1 <- check_columns(markers, each_1)
    } else {
        each_1 <- NULL
        each_2 <- NULL
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
        markers <- dplyr::arrange(markers, !!!rlang::parse_exprs(order_by))
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
        if (!is.null(each_1) && facet_each) {
            args$facet_by <- each_1
        } else if (!is.null(each_1)) {
            args$split_by <- each_1
        }
        do_call(plotthis::VolcanoPlot, args)
    } else if (plot_type %in% c("jitter", "jitter_log2fc", "jitter_pct")) {
        if (is.null(each_1)) {
            stop("[MarkersPlot] `each` is required for plot_type '", plot_type, "'. Consider using volcano plot if you don't have groups.")
        }
        if (!is.numeric(select)) {
            stop("[MarkersPlot] `select` must be numeric for plot_type '", plot_type, "', to label top N markers in each group.")
        }
        args <- list(
            markers,
            x = each_1,
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
        if (is.null(each)) {
            stop("[MarkersPlot] `each` is required for plot_type '", plot_type, "'")
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
            genes <- dplyr::slice_head(markers, n = select, by = !!rlang::sym(each_1))$gene
        } else if (length(select) == 1) {
            genes <- dplyr::filter(markers, !!rlang::parse_expr(select))$gene
        } else {
            # The expressions in select with the entire `each` word in it are
            # supposed to be the ones to filter the data
            select_sb <- grepl(paste0("\\b", each_1, "\\b"), select)
            if (any(select_sb)) {
                markers <- dplyr::filter(markers, !!!rlang::parse_exprs(select[select_sb]))
                if (all(select_sb)) {
                    genes <- markers$gene
                } else {
                    select_non_sb <- select[!select_sb]
                    if (length(select_non_sb) == 1 && grepl("^\\d+$", select_non_sb)) {
                        genes <- dplyr::slice_head(markers, n = as.numeric(select_non_sb), by = !!rlang::sym(each_1))$gene
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
        if (!is.factor(markers[[each_1]])) {
            markers[[each_1]] <- factor(markers[[each_1]], levels = unique(markers[[each_1]]))
        }
        args <- list(
            data = markers,
            values_by = y,
            rows_by = "gene",
            columns_by = each_1,
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
        groups <- levels(markers[[each_1]])
        if (!is.null(cutoff) && startsWith(plot_type, "heatmap_")) {
            sig_mat <- tidyr::pivot_wider(
                markers,
                id_cols = !!sym("gene"),
                names_from = each_1,
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
                    "[MarkersPlot] The following groups in `each` '", each_1,
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
                names_from = each_1,
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
                    "[MarkersPlot] The following groups in `each` '", each_1,
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

        if (is.numeric(select)) {
            if (!is.null(cutoff)) {
                markers <- dplyr::filter(markers, !!rlang::sym(pcol) < cutoff)
            }
            if (!is.null(each_1)) {
                # e.g. each cluster
                genes <- dplyr::slice_head(markers, n = select, by = !!rlang::sym(each_1))

                if (plot_type %in% c("heatmap", "dot") && !flatten_markers) {
                    genes <- dplyr::summarise(genes, gene = list(!!sym("gene")), .by = !!rlang::sym(each_1))
                    # keep the order of the group name
                    genes <- dplyr::arrange(genes, !!rlang::sym(each_1))
                    # convert to a named list, with each group name as the list name
                    genes <- stats::setNames(genes$gene, genes[[each_1]])
                } else {
                    genes <- genes$gene
                }
            } else {
                # generally, select top N markers overall
                genes <- dplyr::slice_head(markers, n = select)$gene
            }
        } else {
            filtered <- dplyr::filter(markers, !!rlang::parse_expr(select))
            if (!is.null(each_1)) {
                # e.g. each cluster
                if (plot_type %in% c("heatmap", "dot") && !flatten_markers) {
                    genes <- dplyr::summarise(filtered, gene = list(!!sym("gene")), .by = !!rlang::sym(each_1))
                    # keep the order of the group name
                    genes <- dplyr::arrange(genes, !!rlang::sym(each_1))
                    # convert to a named list, with each group name as the list name
                    genes <- stats::setNames(genes$gene, genes[[each_1]])
                } else {
                    genes <- filtered$gene
                }
            } else {
                genes <- filtered$gene
            }
        }

        if (!is.null(group_by)) {
            if (grepl(":", group_by)) {
                if (length(strsplit(group_by, ":")[[1]]) != 2) {
                    stop("[MarkersPlot] `group_by` must be in the format 'marker_column:metadata_column' or a single column name.")
                }
                group_by_1 <- strsplit(group_by, ":")[[1]][1]
                group_by_2 <- strsplit(group_by, ":")[[1]][2]
                group_by_1 <- check_columns(markers, group_by_1)
            } else {
                group_by_1 <- NULL
                group_by_2 <- group_by
            }
            if (!group_by_2 %in% colnames(object@meta.data)) {
                stop("[MarkersPlot] `group_by` '", group_by_2, "' is not found in the object's metadata.")
            }
        } else {
            group_by_1 <- NULL
            group_by_2 <- NULL
        }

        if (!is.null(group_by_1) && !is.null(group_by_2)) {
            # check if group_by values are consistent between markers and object
            groups <- unique(unlist(strsplit(unique(as.character(markers[[group_by_1]])), ":")))
            if (!all(groups %in% unique(as.character(object@meta.data[[group_by_2]])))) {
                stop("[MarkersPlot] The following values in `group_by` '", group_by_1, "' are not found in the object's metadata (", group_by_2, "): ",
                    paste(setdiff(groups, unique(as.character(object@meta.data[[group_by_2]]))), collapse = ", ")
                )
            }
            object <- subset_seurat(object, subset = !!rlang::sym(group_by_2) %in% groups)
            object@meta.data[[group_by_2]] <- factor(object@meta.data[[group_by_2]], levels = groups)
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
            plot_type = plot_type,
            ...
        )
        if (plot_type %in% c("heatmap", "dot")) {
            args$name <- args$name %||% "Expression"
            args$cluster_columns <- args$cluster_columns %||% FALSE
            if (!is.null(each_2)) {
                args$facet_by <- NULL
                args$split_by <- NULL
                args$group_by <- NULL
                if (!is.null(each_2) && !is.null(group_by_2)) {
                    args$columns_split_by <- each_2
                    args$ident <- group_by_2
                } else if (!is.null(each_2)) {
                    args$ident <- each_2
                } else if (!is.null(group_by_2)) {
                    args$ident <- group_by_2
                } else {
                    args$ident <- NULL
                }
            }
        } else if (plot_type %in% c("violin", "box", "bar")) {
            args$stack <- args$stack %||% TRUE
            args$group_by <- group_by_2
            args$split_by <- NULL
            args$facet_by <- NULL
            args$ident <- each_2 %||% "orig.ident"
        } else {
            args$group_by <- group_by_2
        }
        do_call(scplotter::FeatureStatPlot, args)
    }
}
