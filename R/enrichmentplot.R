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
#' When the \code{metric} is a p-value column (\code{pvalue}, \code{p.adjust},
#' or \code{qvalue}), the function applies a \eqn{-log_{10}} transformation
#' so that more significant terms have higher values on the plot. The transformed
#' metric is stored internally as \code{.metric}. When \code{cutoff} is specified,
#' it is also transformed (e.g., \code{p.adjust = 0.05} becomes a line at
#' \eqn{-log_{10}(0.05) = 1.3}).
#'
#' For "bar", "dot", "lollipop", "comparison", and "heatmap" plot types,
#' \code{GeneRatio} (stored as strings like \code{"38/225"}) and \code{BgRatio}
#' are automatically converted to numeric ratios.
#'
#' @section Term ordering and selection:
#' For each unique combination of \code{split_by}, \code{group_by}, and
#' \code{facet_by} levels, the function selects the \code{top_term} terms with
#' the smallest metric values (i.e., most significant). This ensures that each
#' facet or split shows its own most relevant terms rather than the globally
#' most significant ones. The default \code{top_term} is 6 for most plot types
#' and 100 for "enrichmap" (which benefits from showing more terms to reveal
#' the network structure of gene set relationships).
#'
#' @param data A data frame with enrichment results. Must be the output of a
#'   \pkg{clusterProfiler} function (\code{enrichGO}, \code{enrichKEGG},
#'   \code{enrichPathway}, \code{enrichWP}, etc.) or an \pkg{enrichR} result
#'   processed through \code{plotthis::prepare_enrichr_result()}. The function
#'   auto-detects the format based on column names.
#' @param top_term Integer. Number of top terms (by significance) to display per
#'   group/facet combination. Default: \code{6} for all plot types except
#'   \code{"enrichmap"} which defaults to \code{100}. Note that terms are not
#'   filtered globally — the top terms are selected independently within each
#'   combination of \code{split_by}, \code{group_by}, and \code{facet_by} levels.
#' @param plot_type Character. The type of plot to generate. One of:
#'   \code{"bar"}, \code{"dot"}, \code{"lollipop"}, \code{"network"},
#'   \code{"enrichmap"}, \code{"wordcloud"}, \code{"comparison"}, or
#'   \code{"heatmap"}. See the Description section for guidance on choosing a
#'   plot type. Default: \code{"bar"}.
#' @param metric Character. The column name in \code{data} to use as the
#'   significance metric for ordering and coloring terms. Common choices are
#'   \code{"p.adjust"} (default), \code{"pvalue"}, or \code{"qvalue"}. When the
#'   metric is a p-value column, a \eqn{-log_{10}} transformation is applied
#'   automatically so that more significant terms have higher values.
#' @param cutoff Numeric. A significance threshold to mark on the plot.
#'   Default: \code{NULL} (no marking). The behavior depends on \code{plot_type}:
#'   \itemize{
#'     \item \code{"bar"} — Adds a vertical dashed line at the transformed cutoff
#'       (e.g., \eqn{-log_{10}(0.05)}).
#'     \item \code{"dot"}, \code{"lollipop"}, \code{"comparison"} — Terms above
#'       the cutoff are colored gray with the legend label from
#'       \code{fill_cutoff_name}.
#'     \item \code{"heatmap"} — Adds asterisk (\code{*}) labels to cells where
#'       the metric exceeds the cutoff.
#'     \item \code{"network"}, \code{"enrichmap"}, \code{"wordcloud"} — No effect.
#'   }
#'   This parameter only marks terms — it does not filter them. Use
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
#'   significance metric). Applies to \code{"comparison"}, \code{"dot"}, and
#'   \code{"lollipop"} plot types. Default: \code{NULL} (auto-generated as
#'   \code{"-log10(metric)"}).
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

    top_term <- top_term %||% ifelse(plot_type == "enrichmap", 100 , 6)
    # if (!is.null(cutoff)) {
    #     data <- data[data[[metric]] < cutoff, , drop = FALSE]
    # }

    # Keep the top terms for each split_by, group_by, and facet_by
    if (!is.null(split_by) || !is.null(group_by) || !is.null(facet_by)) {
        data <- data %>%
            group_by(!!!syms(unique(c(split_by, group_by, facet_by)))) %>%
            slice_min(!!sym(metric), n = top_term, with_ties = FALSE) %>%
            ungroup()
    } else {
        data <- data %>% slice_min(!!sym(metric), n = top_term, with_ties = FALSE)
    }

    # preprocessing
    if (plot_type %in% c("bar", "comparison", "dot", "lollipop", "heatmap")) {
        if (metric %in% c("pvalue", "p.adjust", "qvalue")) {
            data$.metric <- -log10(data[[metric]])
        } else {
            data$.metric <- data[[metric]]
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
        data[[metric]] <- NULL
        if (metric %in% c("pvalue", "p.adjust", "qvalue")) {
            metric <- paste0("-log10(", metric, ")")
            if (!is.null(cutoff)) {
                cutoff <- -log10(cutoff)
            }
        }
        if (is.null(cutoff)) {
            Heatmap(data, in_form = "long", values_by = ".metric", name = metric,
                rows_by = descr_col, columns_by = group_by, values_fill = values_fill, ...)
        } else {
            Heatmap(data, in_form = "long", values_by = ".metric", name = metric,
                rows_by = descr_col, columns_by = group_by, values_fill = values_fill,
                cell_type = "label", label = function(x) {
                    ifelse(x > cutoff, '*', NA)
                }, ...)
        }
    } else if (plot_type == "bar") {
        if (!is.null(group_by)) {
            stop("'group_by' is not supported for Enrichment bar plot. Use 'facet_by'/'split_by' to split the plots.")
        }
        expand <- expand %||% c(0.1, 0.6, 0, 0.6)
        if (!is.null(cutoff)) {
            add_line <- -log10(cutoff)
            line_name <- paste0(metric, " = ", cutoff)
        } else {
            add_line <- NULL
            line_name <- NULL
        }
        BarPlot(data, x = descr_col, y = ".metric", flip = TRUE, label = "Count",
            add_line = add_line, line_name = line_name,
            fill_by = "Count", palette = palette, expand = expand,
            split_by = split_by, facet_by = facet_by, facet_scales = facet_scales %||% "free_y",
            ylab = xlab %||% paste0("-log10(", metric, ")"), xlab = ylab %||% "", ...)
    } else if (plot_type == "comparison") {
        if (is.null(group_by)) {
            stop("'group_by' is required for Enrichment comparison plot.")
        }
        size_by <- size_by %||% "GeneRatio"
        if (!is.null(cutoff)) {
            fill_cutoff <- -log10(cutoff)
            fill_cutoff_name <- fill_cutoff_name %||% "Non-significant"
        } else {
            fill_cutoff <- NULL
            fill_cutoff_name <- fill_cutoff_name
        }
        DotPlot(data, x = group_by, y = descr_col, fill_by = ".metric",
            fill_cutoff = fill_cutoff, fill_cutoff_name = fill_cutoff_name,
            size_by = size_by, palette = palette, fill_name = fill_name %||% paste0("-log10(", metric, ")"),
            split_by = split_by, facet_by = facet_by, facet_scales = facet_scales %||% "free_y",
            xlab = xlab %||% "", ylab = ylab %||% group_by, ...)
    } else if (plot_type == "dot") {
        if (!is.null(group_by)) {
            stop("'group_by' is not supported for Enrichment dot plot. Use 'facet_by'/'split_by' to split the plots.")
        }
        x_by <- x_by %||% "GeneRatio"
        size_by <- size_by %||% "Count"
        if (!is.null(cutoff)) {
            fill_cutoff <- -log10(cutoff)
            fill_cutoff_name <- fill_cutoff_name %||% "Non-significant"
        } else {
            fill_cutoff <- NULL
            fill_cutoff_name <- fill_cutoff_name
        }
        DotPlot(data, x = x_by, y = descr_col, fill_by = ".metric",
            fill_cutoff = fill_cutoff, fill_cutoff_name = fill_cutoff_name,
            size_by = size_by, palette = palette, fill_name = fill_name %||% paste0("-log10(", metric, ")"),
            split_by = split_by, facet_by = facet_by, facet_scales = facet_scales %||% "free_y",
            xlab = xlab %||% x_by, ylab = ylab %||% "", ...)
    } else if (plot_type == "lollipop") {
        if (!is.null(group_by)) {
            stop("'group_by' is not supported for Enrichment lollipop plot. Use 'facet_by'/'split_by' to split the plots.")
        }
        x_by <- x_by %||% "GeneRatio"
        size_by <- size_by %||% "Count"
        if (!is.null(cutoff)) {
            fill_cutoff <- -log10(cutoff)
            fill_cutoff_name <- fill_cutoff_name %||% "Non-significant"
        } else {
            fill_cutoff <- NULL
            fill_cutoff_name <- fill_cutoff_name
        }
        LollipopPlot(data, x = x_by, y = descr_col, fill_by = ".metric",
            fill_cutoff = fill_cutoff, fill_cutoff_name = fill_cutoff_name,
            size_by = size_by, fill_name = fill_name %||% paste0("-log10(", metric, ")"),
            palette = palette, split_by = split_by, facet_by = facet_by, facet_scales = facet_scales %||% "free_y",
            xlab = xlab %||% x_by, ylab = ylab %||% "", ...)
    } else if (plot_type == "network") {
        if (!is.null(group_by)) {
            stop("'group_by' is not supported for Enrichment network plot. Use 'facet_by'/'split_by' to split the plots.")
        }
        EnrichNetwork(data, in_form = "clusterProfiler", character_width = character_width, split_by = split_by,
            facet_by = facet_by, facet_scales = facet_scales, metric = metric, palette = palette, xlab = xlab, ylab = ylab, ...)
    } else if (plot_type == "enrichmap") {
        if (!is.null(group_by)) {
            stop("'group_by' is not supported for Enrichment enrichmap plot. Use 'facet_by'/'split_by' to split the plots.")
        }
        EnrichMap(data, in_form = "clusterProfiler", character_width = character_width, split_by = split_by,
            facet_by = facet_by, facet_scales = facet_scales, metric = metric, palette = palette, xlab = xlab, ylab = ylab, ...)
    } else if (plot_type == "wordcloud") {
        if (!is.null(group_by)) {
            stop("'group_by' is not supported for Enrichment wordcloud plot. Use 'facet_by'/'split_by' to split the plots.")
        }
        word_type <- match.arg(word_type)
        data$.metric <- -log10(data[[metric]])
        if (word_type == "feature") {
            data$Genes <- lapply(data$geneID, gsub, pattern = "/", replacement = " ", fixed = TRUE)
            WordCloudPlot(data, sentence_by = "Genes", score_name = paste0("sum(-log10(", metric, "))"),
                split_by = split_by, facet_by = facet_by, score_by = ".metric", score_agg = sum,
                facet_scales = facet_scales, palette = palette, xlab = xlab, ylab = ylab, ...)
        } else {
            WordCloudPlot(data, sentence_by = descr_col, score_name = paste0("sum(-log10(", metric, "))"),
                split_by = split_by, facet_by = facet_by, score_by = ".metric", score_agg = sum,
                facet_scales = facet_scales, palette = palette, xlab = xlab, ylab = ylab, ...)
        }
    }
}
