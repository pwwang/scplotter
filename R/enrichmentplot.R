#' Visualize gene set enrichment and over-representation analysis results
#'
#' @description
#' Gene set enrichment analysis identifies biological pathways, gene ontologies,
#' or functional categories that are statistically over-represented among a list
#' of genes of interest (e.g., differentially expressed genes from a single-cell
#' RNA-seq experiment). Rather than interpreting individual genes in isolation,
#' enrichment analysis places gene-level results into a broader biological context,
#' revealing which processes, functions, or diseases are perturbed.
#'
#' \code{EnrichmentPlot} generates publication-quality visualizations for enrichment
#' results across eight distinct plot types, each suited to a different analytical
#' perspective:
#' \itemize{
#'   \item \strong{bar} — Horizontal bar chart of the top enriched terms, ordered by
#'     significance. Best for a quick overview or when showing a small number of terms.
#'   \item \strong{dot} — Dot plot where x-axis shows a continuous metric (default:
#'     \code{GeneRatio}), dot size reflects gene count, and dot color reflects
#'     significance. Ideal for comparing terms along two dimensions simultaneously.
#'   \item \strong{lollipop} — Lollipop chart combining dot and bar aesthetics.
#'     Similar to the dot plot but with stems emphasizing the ranking.
#'   \item \strong{comparison} — Side-by-side dot plot comparing enrichment across
#'     groups (e.g., cell types, conditions). Requires \code{group_by}.
#'   \item \strong{network} — Network visualization where nodes are enriched terms
#'     and edges represent overlapping gene sets. Reveals functional modules and
#'     redundant terms.
#'   \item \strong{enrichmap} — Enrichment map similar to the network plot but
#'     optimized for large term sets (default \code{top_term = 100}). Nodes are
#'     terms and edges represent gene overlap.
#'   \item \strong{wordcloud} — Word cloud where term size reflects significance.
#'     Can display either enrichment terms (\code{word_type = "term"}) or
#'     individual gene symbols (\code{word_type = "feature"}).
#'   \item \strong{heatmap} — Heatmap of enrichment significance across groups
#'     (\code{group_by} is mapped to columns). Useful for comparing enrichment
#'     patterns across multiple conditions or cell types.
#' }
#'
#' The function auto-detects the input data format (\pkg{clusterProfiler} or
#' \pkg{enrichR}) and delegates visualization to the appropriate \pkg{plotthis}
#' plotting function.
#'
#' @section Input data formats:
#' The function auto-detects the input format by checking for characteristic
#' column names:
#' \itemize{
#'   \item \strong{clusterProfiler} (\code{enrichGO}, \code{enrichKEGG},
#'     \code{enrichPathway}, etc.) — recognized by the presence of
#'     \code{pvalue}, \code{p.adjust}, and \code{qvalue} columns.
#'   \item \strong{enrichR} (web-based enrichment tool) — recognized by the
#'     presence of \code{P.value} and \code{Adjusted.P.value} columns.
#'     enrichR results are automatically converted to clusterProfiler-compatible
#'     format via \pkg{plotthis}'s \code{prepare_enrichr_result()}.
#' }
#' If neither format is detected, the function stops with an error.
#'
#' @section Metric transformation:
#' When a metric column is a p-value column (\code{pvalue}, \code{p.adjust},
#' or \code{qvalue}), the function applies a \eqn{-log_{10}} transformation
#' so that more significant terms have higher values on the plot. The transformed
#' metric is stored internally as \code{.metric}. Other columns are used as-is.
#'
#' \code{metric} accepts \code{"color_column:significance_column"}. The first
#' column drives the fill color and the \code{top_term} selection; the second
#' drives \code{cutoff}, which both drops the terms failing it and marks the
#' surviving ones. Both columns are transformed independently, so
#' \code{metric = "zScore:p.adjust"} colors cells by raw \code{zScore}, keeps
#' only the terms with \code{p.adjust < cutoff}, and marks them. A single column
#' is shorthand for \code{"column:column"}, so \code{metric = "p.adjust"} marks
#' cells by \code{p.adjust} without filtering.
#'
#' For "bar", "dot", "lollipop", "comparison", and "heatmap" plot types,
#' \code{GeneRatio} (stored as strings like \code{"38/225"}) and \code{BgRatio}
#' are automatically converted to numeric ratios.
#'
#' @section Term ordering and selection:
#' For each unique combination of \code{split_by}, \code{group_by}, and
#' \code{facet_by} levels, the function selects the \code{top_term} terms by
#' the coloring metric (the first column of \code{metric}): the smallest values
#' when it is a p-value column (\code{pvalue}, \code{p.adjust}, \code{qvalue}),
#' and the largest values otherwise (e.g. \code{zScore}). This ensures that each
#' facet or split shows its own most relevant terms rather than the globally
#' most relevant ones. When a significance metric is given
#' (\code{"coloring:sig"}) together with \code{cutoff}, terms that fail the
#' cutoff in every group are dropped first (except for dot, lollipop, and
#' comparison plots, which show them in gray); a term passing in any group is
#' kept in all groups. The default \code{top_term} is 6 for most plot types
#' and 100 for "enrichmap" (which benefits from showing more terms to reveal
#' the network structure of gene set relationships).
#'
#' @param data A data frame with enrichment results. Must be the output of a
#'   \pkg{clusterProfiler} function (\code{enrichGO}, \code{enrichKEGG},
#'   \code{enrichPathway}, \code{enrichWP}, etc.) or an \pkg{enrichR} result
#'   processed through \code{plotthis::prepare_enrichr_result()}. The function
#'   auto-detects the format based on column names.
#' @param top_term Integer. Number of top terms to display per group/facet
#'   combination, ranked by the coloring metric: smallest values for p-value
#'   columns, largest values otherwise. Default: \code{6} for all plot types
#'   except \code{"enrichmap"} which defaults to \code{100}. Note that terms are
#'   not filtered globally — the top terms are selected independently within each
#'   combination of \code{split_by}, \code{group_by}, and \code{facet_by} levels.
#' @param plot_type Character. The type of plot to generate. One of:
#'   \code{"bar"}, \code{"dot"}, \code{"lollipop"}, \code{"network"},
#'   \code{"enrichmap"}, \code{"wordcloud"}, \code{"comparison"}, or
#'   \code{"heatmap"}. See the Description section for guidance on choosing a
#'   plot type. Default: \code{"bar"}.
#' @param metric Character. The column name in \code{data} to use for coloring
#'   terms, or two column names separated by \code{":"} to use one column for
#'   coloring and another for significance. \code{"p.adjust"} (default),
#'   \code{"pvalue"}, and \code{"qvalue"} are transformed with
#'   \eqn{-log_{10}} so that more significant terms have higher values; other
#'   columns are used as-is. With \code{"zScore:p.adjust"}, cells are colored by
#'   \code{zScore} while \code{cutoff} is applied to \code{p.adjust}. A single
#'   column is shorthand for \code{"column:column"}. See the
#'   \emph{Metric transformation} section.
#' @param cutoff Numeric. A significance threshold applied to the significance
#'   metric (the column after \code{":"} in \code{metric}, or \code{metric}
#'   itself when no modifier is given). Default: \code{NULL} (no marking).
#'   When \code{metric} supplies a separate significance metric, a term failing
#'   the cutoff in every group is dropped before \code{top_term} selection. A
#'   term that passes in at least one group is kept in all of them, so its
#'   cells are not left empty. Dot, lollipop, and comparison plots are the
#'   exception: they show insignificant terms in gray, so nothing is dropped.
#'   The behavior on the plot itself depends on \code{plot_type}:
#'   \itemize{
#'     \item \code{"bar"} — Adds a vertical dashed line at the transformed cutoff
#'       (e.g., \eqn{-log_{10}(0.05)}). No line is drawn when the significance
#'       metric differs from the coloring metric, since it would sit on a
#'       different scale; terms are still filtered.
#'     \item \code{"dot"}, \code{"lollipop"}, \code{"comparison"} — Terms that
#'       fail the cutoff are colored gray with the legend label from
#'       \code{fill_cutoff_name}, while fill color still follows the coloring
#'       metric.
#'     \item \code{"heatmap"} — Adds asterisk (\code{*}) labels to cells whose
#'       significance metric passes the cutoff, plus a legend entry labelled
#'       \code{"<significance metric> < <cutoff>"} (legend title from
#'       \code{label_name}, default \code{"Significance"}).
#'     \item \code{"network"}, \code{"enrichmap"}, \code{"wordcloud"} — Nothing
#'       is marked on the plot; terms are still filtered when a significance
#'       metric is given.
#'   }
#'   Without a separate significance metric, this parameter only marks terms —
#'   it does not filter them. Use
#'   \code{top_term} to control how many terms are shown.
#' @param x_by Character. Column name(s) to use for the x-axis. Works only for
#'   \code{"dot"} and \code{"lollipop"} plot types. Default: \code{NULL}
#'   (defaults to \code{"GeneRatio"} internally).
#' @param size_by Character. Column name(s) to map to point size. Works only for
#'   \code{"comparison"}, \code{"dot"}, and \code{"lollipop"} plot types.
#'   Default: \code{NULL} (defaults to \code{"GeneRatio"} for comparison,
#'   \code{"Count"} for dot and lollipop).
#' @param fill_cutoff_name Character. Legend label for terms that exceed the
#'   \code{cutoff} (shown in gray). Applies to \code{"comparison"},
#'   \code{"dot"}, and \code{"lollipop"} plot types. Default: \code{NULL}
#'   (defaults to \code{"Non-significant"} when \code{cutoff} is set).
#' @param fill_name Character. Legend title for the fill color scale (the
#'   coloring metric). Applies to \code{"comparison"}, \code{"dot"}, and
#'   \code{"lollipop"} plot types. Default: \code{NULL} (auto-generated from
#'   \code{metric}, e.g. \code{"-log10(p.adjust)"} or \code{"zScore"}).
#' @param values_fill Numeric. The fill value for missing entries in the heatmap
#'   matrix. Used only for \code{"heatmap"} plot type. Default: \code{0}.
#' @param character_width Integer. Maximum character width for term descriptions
#'   before line-wrapping. Applies to all plot types; for \code{"heatmap"} the
#'   wrapping is deferred to the Heatmap function. Default: \code{50}.
#' @param expand Numeric vector of length 1, 2, or 4. Axis expansion factors
#'   passed to \code{\link[plotthis:BarPlot]{plotthis::BarPlot()}}. Used only
#'   for \code{"bar"} plot type. Default: \code{NULL} (defaults to
#'   \code{c(0.1, 0.6, 0, 0.6)}).
#' @param word_type Character. What to display in the wordcloud. One of
#'   \code{"term"} (enrichment term descriptions) or \code{"feature"} (gene
#'   symbols from the enriched gene list). Used only for \code{"wordcloud"}
#'   plot type. Default: \code{"term"}.
#' @param split_by Character vector. Column name(s) in \code{data} to split the
#'   data and generate separate plots for each unique value. Multiple columns
#'   are concatenated with \code{split_by_sep}. Default: \code{NULL}.
#' @param split_by_sep Character. Separator used when concatenating multiple
#'   \code{split_by} columns. Default: \code{"_"}.
#' @param facet_by Character vector. Column name(s) in \code{data} to use for
#'   faceting (generating sub-panels within each plot). Default: \code{NULL}.
#' @param facet_scales Character. Facet scale behavior — \code{"fixed"}
#'   (same scales), \code{"free"}, \code{"free_x"}, or \code{"free_y"}.
#'   Default: \code{NULL} (defaults to \code{"free_y"} for bar, dot, lollipop,
#'   and comparison plots).
#' @param group_by Character vector. Column name(s) in \code{data} to group
#'   terms. Behavior depends on \code{plot_type}:
#'   \itemize{
#'     \item \code{"comparison"} — Groups are shown as x-axis categories in a
#'       dot plot comparing enrichment across groups. Required for this type.
#'     \item \code{"heatmap"} — Groups are used as the columns of the heatmap
#'       (mapped to \code{columns_by} in \code{\link[plotthis:Heatmap]{plotthis::Heatmap()}}).
#'     \item All other types — \code{group_by} is not supported and will raise
#'       an error. Use \code{facet_by} or \code{split_by} instead.
#'   }
#'   Multiple columns are concatenated with \code{group_by_sep}.
#'   Default: \code{NULL}.
#' @param group_by_sep Character. Separator used when concatenating multiple
#'   \code{group_by} columns. Used only for \code{"comparison"} plot type.
#'   Default: \code{"_"}.
#' @param palette Character. Color palette name for the fill scale. See
#'   \code{\link[plotthis:show_palettes]{plotthis::show_palettes()}} for
#'   available palettes. Default: \code{"Spectral"}.
#' @param xlab Character. Custom x-axis label. Default: \code{NULL}
#'   (auto-generated based on plot type and \code{x_by}/\code{metric}).
#' @param ylab Character. Custom y-axis label. Default: \code{NULL}
#'   (auto-generated based on plot type).
#' @param ... Additional arguments passed to the underlying \pkg{plotthis}
#'   plotting function, determined by \code{plot_type}:
#'   \describe{
#'     \item{\code{"bar"}}{\code{\link[plotthis:BarPlot]{plotthis::BarPlot()}}}
#'     \item{\code{"dot"}}{\code{\link[plotthis:DotPlot]{plotthis::DotPlot()}}}
#'     \item{\code{"lollipop"}}{\code{\link[plotthis:LollipopPlot]{plotthis::LollipopPlot()}}}
#'     \item{\code{"network"}}{\code{\link[plotthis:EnrichNetwork]{plotthis::EnrichNetwork()}}}
#'     \item{\code{"enrichmap"}}{\code{\link[plotthis:EnrichMap]{plotthis::EnrichMap()}}}
#'     \item{\code{"wordcloud"}}{\code{\link[plotthis:WordCloudPlot]{plotthis::WordCloudPlot()}}}
#'     \item{\code{"comparison"}}{\code{\link[plotthis:DotPlot]{plotthis::DotPlot()}}}
#'     \item{\code{"heatmap"}}{\code{\link[plotthis:Heatmap]{plotthis::Heatmap()}}}
#'   }
#' @return A ggplot object (or a \code{patchwork} object when \code{split_by}
#'   generates multiple plots and \code{combine = TRUE}), or a list of ggplot
#'   objects if \code{combine = FALSE}. The specific return type depends on the
#'   underlying \pkg{plotthis} function dispatched by \code{plot_type}.
#' @note
#' \itemize{
#'   \item The function auto-detects \pkg{clusterProfiler} vs \pkg{enrichR} input
#'     format. For \pkg{enrichR} input, it must contain \code{P.value} and
#'     \code{Adjusted.P.value} columns. If your enrichR results came from a
#'     different pipeline, pre-process them with
#'     \code{plotthis::prepare_enrichr_result()}.
#'   \item The \code{cutoff} parameter only \strong{marks} terms — it does not
#'     filter them. To reduce the number of displayed terms, use \code{top_term}.
#'   \item \code{GeneRatio} strings (e.g., \code{"38/225"}) and \code{BgRatio}
#'     strings are automatically converted to numeric values by dividing the
#'     numerator by the denominator.
#'   \item When using \code{group_by} with \code{plot_type = "comparison"},
#'     \code{size_by} defaults to \code{"GeneRatio"} and each group's terms are
#'     shown side-by-side. For \code{plot_type = "heatmap"}, \code{group_by}
#'     becomes the heatmap columns.
#'   \item \code{group_by} is not supported for \code{"bar"}, \code{"dot"},
#'     \code{"lollipop"}, \code{"network"}, \code{"enrichmap"}, and
#'     \code{"wordcloud"} — use \code{facet_by} or \code{split_by} to separate
#'     groups for those types.
#'   \item For \code{plot_type = "wordcloud"} with \code{word_type = "feature"},
#'     individual gene symbols are extracted from the \code{geneID} column.
#'     Gene-level significance scores are aggregated by summing
#'     \eqn{-log_{10}(p)} values.
#' }
#' @seealso
#' \code{\link[plotthis:BarPlot]{plotthis::BarPlot()}},
#' \code{\link[plotthis:DotPlot]{plotthis::DotPlot()}},
#' \code{\link[plotthis:LollipopPlot]{plotthis::LollipopPlot()}},
#' \code{\link[plotthis:EnrichNetwork]{plotthis::EnrichNetwork()}},
#' \code{\link[plotthis:EnrichMap]{plotthis::EnrichMap()}},
#' \code{\link[plotthis:WordCloudPlot]{plotthis::WordCloudPlot()}},
#' \code{\link[plotthis:Heatmap]{plotthis::Heatmap()}},
#' \code{\link[plotthis:show_palettes]{plotthis::show_palettes()}}
#' @importFrom rlang sym syms
#' @importFrom stringr str_wrap
#' @importFrom dplyr %>% group_by slice_min ungroup
#' @importFrom plotthis BarPlot DotPlot LollipopPlot EnrichNetwork EnrichMap WordCloudPlot Heatmap
#' @export
#' @examples
#' \donttest{
#' set.seed(8525)
#' data(enrich_example, package = "plotthis")
#' enrich_example$Group <- sample(LETTERS[1:3], nrow(enrich_example), replace = TRUE)
#' data(enrich_multidb_example, package = "plotthis")
#'
#' EnrichmentPlot(enrich_example)
#' EnrichmentPlot(enrich_example, cutoff = 0.05)
#' EnrichmentPlot(enrich_example, palette = "Paired")
#'
#' enrich_example$Description <- enrich_example$ID
#' EnrichmentPlot(enrich_example, plot_type = "heatmap", group_by = "Group",
#'  show_row_names = TRUE, show_column_names = TRUE, cutoff = 0.05)
#'
#' # Color by zScore, and mark the cells whose p.adjust passes the cutoff
#' EnrichmentPlot(enrich_example, plot_type = "heatmap", group_by = "Group",
#'  metric = "zScore:p.adjust", cutoff = 0.05)
#' EnrichmentPlot(enrich_example, plot_type = "dot",
#'  metric = "zScore:p.adjust", cutoff = 1e-10)
#'
#' # Multiple databases#'
#' EnrichmentPlot(enrich_multidb_example, facet_by = "Database", facet_nrow = 2)
#'
#' enrich_example$Group <- sample(c("A", "B"), nrow(enrich_example), replace = TRUE)
#' EnrichmentPlot(enrich_example, plot_type = "comparison", group_by = "Group")
#' EnrichmentPlot(enrich_example, plot_type = "dot", top_term = 10)
#' EnrichmentPlot(enrich_example, plot_type = "lollipop", top_term = 10)
#' EnrichmentPlot(enrich_example, plot_type = "network")
#' EnrichmentPlot(enrich_example, plot_type = "enrichmap")
#' EnrichmentPlot(enrich_example, plot_type = "wordcloud")
#' # Wordcloud with feature
#' EnrichmentPlot(enrich_example, plot_type = "wordcloud", word_type = "feature")
#' }
EnrichmentPlot <- function(
    data, top_term = NULL,
    plot_type = c("bar", "dot", "lollipop", "network", "enrichmap", "wordcloud", "comparison", "heatmap"),
    x_by = NULL, size_by = NULL, fill_cutoff_name = NULL, fill_name = NULL, values_fill = 0,
    character_width = 50, expand = NULL, word_type = c("term", "feature"),
    split_by = NULL, split_by_sep = "_", facet_by = NULL, facet_scales = NULL,
    group_by = NULL, group_by_sep = "_", metric = "p.adjust", cutoff = NULL,
    palette = "Spectral", xlab = NULL, ylab = NULL,
    ...
) {
    if (all(c("pvalue", "p.adjust", "qvalue") %in% colnames(data))) {
        in_form <- "clusterProfiler"
    } else if (all(c("P.value", "Adjusted.P.value") %in% colnames(data))) {
        in_form <- "enrichr"
    } else {
        stop("[EnrichmentPlot] Cannot infer the input data format. The data should be the output of either 'clusterProfiler' or 'enrichr'.")
    }

    if (in_form == "enrichr") {
        prepare_enrichr_result <- utils::getFromNamespace("prepare_enrichr_result", "plotthis")
        data <- tryCatch(
            prepare_enrichr_result(data, n_input = NULL),
            error = function(e) {
                n_input <- length(unique(unlist(strsplit(unlist(data$Genes), ";", fixed = TRUE)))) * 2
                prepare_enrichr_result(data, n_input = n_input)
            }
        )
    }

    descr_col <- "Description"
    plot_type <- match.arg(plot_type)
    split_by <- check_columns(
        data, split_by, force_factor = TRUE, allow_multi = TRUE,
        concat_multi = TRUE, concat_sep = split_by_sep
    )
    group_by <- check_columns(
        data, group_by, force_factor = TRUE, allow_multi = TRUE,
        concat_multi = TRUE, concat_sep = group_by_sep
    )
    facet_by <- check_columns(
        data, facet_by, force_factor = TRUE, allow_multi = TRUE
    )
    descr_col <- check_columns(data, descr_col, force_factor = TRUE)

    p_metrics <- c("pvalue", "p.adjust", "qvalue")
    metric <- strsplit(metric, ":", fixed = TRUE)[[1]]
    if (length(metric) > 2) {
        stop("[EnrichmentPlot] 'metric' must be 'column' or 'column:column'.")
    }
    fill_metric <- metric[1]
    has_sig_metric <- length(metric) == 2
    sig_metric <- if (has_sig_metric) metric[2] else fill_metric
    if (!sig_metric %in% colnames(data)) {
        stop("[EnrichmentPlot] Significance metric '", sig_metric, "' not found in data.")
    }
    metric_label <- if (fill_metric %in% p_metrics) {
        paste0("-log10(", fill_metric, ")")
    } else {
        fill_metric
    }

    # A separate significance metric makes the cutoff a filter, not just a mark.
    # Dot/lollipop/comparison can show insignificant terms in gray, so keep them.
    if (has_sig_metric && !is.null(cutoff) &&
            !plot_type %in% c("dot", "lollipop", "comparison")) {
        # Keep a term in every group as long as it passes the cutoff in one of
        # them, so its cells in the other groups are not left empty.
        sig_terms <- unique(data[[descr_col]][data[[sig_metric]] < cutoff])
        data <- data[data[[descr_col]] %in% sig_terms, , drop = FALSE]
        if (nrow(data) == 0) {
            stop(
                "[EnrichmentPlot] No terms pass the cutoff: ",
                sig_metric, " < ", cutoff, "."
            )
        }
    }

    top_term <- top_term %||% ifelse(plot_type == "enrichmap", 100 , 6)
    # Keep the top terms for each split_by, group_by, and facet_by.
    # Smallest values for p-value metrics; largest otherwise.
    slice_terms <- if (fill_metric %in% p_metrics) slice_min else slice_max
    if (!is.null(split_by) || !is.null(group_by) || !is.null(facet_by)) {
        data <- data %>%
            group_by(!!!syms(unique(c(split_by, group_by, facet_by)))) %>%
            slice_terms(!!sym(fill_metric), n = top_term, with_ties = FALSE) %>%
            ungroup()
    } else {
        data <- data %>% slice_terms(!!sym(fill_metric), n = top_term, with_ties = FALSE)
    }

    # preprocessing
    if (plot_type %in% c("bar", "comparison", "dot", "lollipop", "heatmap")) {
        if (fill_metric %in% p_metrics) {
            data$.metric <- -log10(data[[fill_metric]])
        } else {
            data$.metric <- data[[fill_metric]]
        }
        if (sig_metric %in% p_metrics) {
            data$.sig_metric <- -log10(data[[sig_metric]])
        } else {
            data$.sig_metric <- data[[sig_metric]]
        }
        if (!is.null(cutoff) && sig_metric %in% p_metrics) {
            sig_cutoff <- -log10(cutoff)
        } else {
            sig_cutoff <- cutoff
        }
        if (plot_type != "heatmap") {
            # we lost order?
            data[[descr_col]] <- str_wrap(data[[descr_col]], width = character_width)
        }
        # Convert GeneRatio from something like "38/225" to 0.169
        data$GeneRatio <- as.numeric(sapply(strsplit(data$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))
        if (!is.null(data$BgRatio) && !all(is.na(data$BgRatio))) {
            data$BgRatio <- as.numeric(sapply(strsplit(data$BgRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))
        }
    }

    if (plot_type == "heatmap") {
        if (!is.null(facet_by)) {
            stop('[EnrichmentPlot] "heatmap" plot does not support "facet_by". Use "split_by" to split the heatmap.')
        }
        if (is.null(group_by)) {
            stop('[EnrichmentPlot] "group_by" is required for "heatmap" plot.')
        }
        data[[fill_metric]] <- NULL
        data[[sig_metric]] <- NULL
        hm_args <- list(
            data = data, in_form = "long", values_by = ".metric", name = metric_label,
            rows_by = descr_col, columns_by = group_by, values_fill = values_fill
        )
        dots <- list(...)
        if (is.null(cutoff)) {
            do_call(Heatmap, c(hm_args, dots))
        } else {
            legend <- paste0(sig_metric, " < ", cutoff)
            hm_args$cell_type <- "label"
            if (sig_metric == fill_metric) {
                hm_args$label <- function(x) {
                    if (isTRUE(x > sig_cutoff)) list("*", legend = legend) else NA
                }
            } else {
                # plotthis names the heatmap columns "<split level> // <group>" (and the rows
                # "<row split> // <description>") through `unite(sep = " // ")`, so the
                # significance matrix has to use the same keys. Keying it by `group_by` alone
                # paints a significant cell into the split panel where that term has no data
                # at all. `i`/`j` index the displayed (clustered) matrix, so they can't be
                # used here, and `pindex()` does not undo clustering either.
                split_keys <- function(split_by, column_by, concat_sep) {
                    split_by <- unique(split_by)
                    if (length(split_by) > 1) {
                        # plotthis concatenates multiple split columns first
                        split_by <- paste(split_by, collapse = concat_sep)
                    }
                    unique(c(split_by, column_by))
                }
                col_keys <- split_keys(
                    dots$columns_split_by, group_by,
                    dots$columns_split_by_sep %||% "_"
                )
                row_keys <- split_keys(
                    dots$rows_split_by, descr_col,
                    dots$rows_split_by_sep %||% "_"
                )
                sig_mat <- data %>%
                    select(all_of(c(col_keys, row_keys, ".sig_metric"))) %>%
                    tidyr::unite(".col", all_of(col_keys), sep = " // ") %>%
                    tidyr::unite(".row", all_of(row_keys), sep = " // ") %>%
                    # duplicated cells are aggregated by plotthis with `cell_agg` ("nanmean")
                    summarise(
                        .sig_metric = if (all(is.na(.sig_metric))) NA_real_ else
                            mean(.sig_metric, na.rm = TRUE),
                        .by = c(".col", ".row")
                    ) %>%
                    pivot_wider(
                        id_cols = ".row", names_from = ".col",
                        values_from = ".sig_metric", values_fill = NA
                    ) %>%
                    as.data.frame()
                rownames(sig_mat) <- sig_mat[[".row"]]
                sig_mat[[".row"]] <- NULL
                sig_mat <- as.matrix(sig_mat)
                hm_args$label <- function(x, i, j, rname, cname) {
                    if (!rname %in% rownames(sig_mat) || !cname %in% colnames(sig_mat)) {
                        return(NA)
                    }
                    if (isTRUE(sig_mat[rname, cname] > sig_cutoff)) {
                        list("*", legend = legend)
                    } else {
                        NA
                    }
                }
            }
            dots$label_name <- dots$label_name %||% "Significance"
            do_call(Heatmap, c(hm_args, dots))
        }
    } else if (plot_type == "bar") {
        if (!is.null(group_by)) {
            stop("'group_by' is not supported for Enrichment bar plot. Use 'facet_by'/'split_by' to split the plots.")
        }
        expand <- expand %||% c(0.1, 0.6, 0, 0.6)
        if (!is.null(cutoff) && sig_metric == fill_metric) {
            add_line <- sig_cutoff
            line_name <- paste0(sig_metric, " = ", cutoff)
        } else {
            if (!is.null(cutoff)) {
                warning(
                    "[EnrichmentPlot] No cutoff line is drawn for bar plots when the significance metric ",
                    "differs from the fill metric; terms failing the cutoff are still dropped.", immediate. = TRUE
                )
            }
            add_line <- NULL
            line_name <- NULL
        }
        BarPlot(data, x = descr_col, y = ".metric", flip = TRUE, label = "Count",
            add_line = add_line, line_name = line_name,
            fill_by = "Count", palette = palette, expand = expand,
            split_by = split_by, facet_by = facet_by, facet_scales = facet_scales %||% "free_y",
            ylab = xlab %||% metric_label, xlab = ylab %||% "", ...)
    } else if (plot_type == "comparison") {
        if (is.null(group_by)) {
            stop("'group_by' is required for Enrichment comparison plot.")
        }
        size_by <- size_by %||% "GeneRatio"
        if (!is.null(cutoff)) {
            fill_cutoff <- sig_cutoff
            fill_cutoff_name <- fill_cutoff_name %||% "Non-significant"
        } else {
            fill_cutoff <- NULL
            fill_cutoff_name <- fill_cutoff_name
        }
        DotPlot(data, x = group_by, y = descr_col, fill_by = ".metric",
            fill_cutoff = fill_cutoff, fill_cutoff_by = ".sig_metric",
            fill_cutoff_name = fill_cutoff_name,
            size_by = size_by, palette = palette, fill_name = fill_name %||% metric_label,
            split_by = split_by, facet_by = facet_by, facet_scales = facet_scales %||% "free_y",
            xlab = xlab %||% "", ylab = ylab %||% group_by, ...)
    } else if (plot_type == "dot") {
        if (!is.null(group_by)) {
            stop("'group_by' is not supported for Enrichment dot plot. Use 'facet_by'/'split_by' to split the plots.")
        }
        x_by <- x_by %||% "GeneRatio"
        size_by <- size_by %||% "Count"
        if (!is.null(cutoff)) {
            fill_cutoff <- sig_cutoff
            fill_cutoff_name <- fill_cutoff_name %||% "Non-significant"
        } else {
            fill_cutoff <- NULL
            fill_cutoff_name <- fill_cutoff_name
        }
        DotPlot(data, x = x_by, y = descr_col, fill_by = ".metric",
            fill_cutoff = fill_cutoff, fill_cutoff_by = ".sig_metric",
            fill_cutoff_name = fill_cutoff_name,
            size_by = size_by, palette = palette, fill_name = fill_name %||% metric_label,
            split_by = split_by, facet_by = facet_by, facet_scales = facet_scales %||% "free_y",
            xlab = xlab %||% x_by, ylab = ylab %||% "", ...)
    } else if (plot_type == "lollipop") {
        if (!is.null(group_by)) {
            stop("'group_by' is not supported for Enrichment lollipop plot. Use 'facet_by'/'split_by' to split the plots.")
        }
        x_by <- x_by %||% "GeneRatio"
        size_by <- size_by %||% "Count"
        if (!is.null(cutoff)) {
            fill_cutoff <- sig_cutoff
            fill_cutoff_name <- fill_cutoff_name %||% "Non-significant"
        } else {
            fill_cutoff <- NULL
            fill_cutoff_name <- fill_cutoff_name
        }
        LollipopPlot(data, x = x_by, y = descr_col, fill_by = ".metric",
            fill_cutoff = fill_cutoff, fill_cutoff_by = ".sig_metric",
            fill_cutoff_name = fill_cutoff_name,
            size_by = size_by, fill_name = fill_name %||% metric_label,
            palette = palette, split_by = split_by, facet_by = facet_by, facet_scales = facet_scales %||% "free_y",
            xlab = xlab %||% x_by, ylab = ylab %||% "", ...)
    } else if (plot_type == "network") {
        if (!is.null(group_by)) {
            stop("'group_by' is not supported for Enrichment network plot. Use 'facet_by'/'split_by' to split the plots.")
        }
        EnrichNetwork(data, in_form = "clusterProfiler", character_width = character_width, split_by = split_by,
            facet_by = facet_by, facet_scales = facet_scales, metric = fill_metric, palette = palette, xlab = xlab, ylab = ylab, ...)
    } else if (plot_type == "enrichmap") {
        if (!is.null(group_by)) {
            stop("'group_by' is not supported for Enrichment enrichmap plot. Use 'facet_by'/'split_by' to split the plots.")
        }
        EnrichMap(data, in_form = "clusterProfiler", character_width = character_width, split_by = split_by,
            facet_by = facet_by, facet_scales = facet_scales, metric = fill_metric, palette = palette, xlab = xlab, ylab = ylab, ...)
    } else if (plot_type == "wordcloud") {
        if (!is.null(group_by)) {
            stop("'group_by' is not supported for Enrichment wordcloud plot. Use 'facet_by'/'split_by' to split the plots.")
        }
        word_type <- match.arg(word_type)
        data$.metric <- -log10(data[[fill_metric]])
        if (word_type == "feature") {
            data$Genes <- lapply(data$geneID, gsub, pattern = "/", replacement = " ", fixed = TRUE)
            WordCloudPlot(data, sentence_by = "Genes", score_name = paste0("sum(-log10(", fill_metric, "))"),
                split_by = split_by, facet_by = facet_by, score_by = ".metric", score_agg = sum,
                facet_scales = facet_scales, palette = palette, xlab = xlab, ylab = ylab, ...)
        } else {
            WordCloudPlot(data, sentence_by = descr_col, score_name = paste0("sum(-log10(", fill_metric, "))"),
                split_by = split_by, facet_by = facet_by, score_by = ".metric", score_agg = sum,
                facet_scales = facet_scales, palette = palette, xlab = xlab, ylab = ylab, ...)
        }
    }
}
