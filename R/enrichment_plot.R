#' Enrichment analysis result plot
#'
#' @description This function generates various types of plots for enrichment (over-representation) analysis.
#'
#' @inheritParams validate_common_arguments
#'
#' @param x A data frame with Enrichment results
#' @param enrichtool The enrichment tool used to generate the enrichment results. Default is NULL.
#'   Will automatically detect the enrichment tool from the enrichment results.
#'   Currently supported tools are "enrichr" and "clusterprofiler".
#' @param db The database to use for enrichment plot. Default: all Database in the enrichment results.
#'   If Database is not specified in the enrichment results, the default value is this value.
#' @param plot_type The type of plot to generate. Options are: "bar", "dot", "lollipop", "network", "enrichmap", "wordcloud", "comparison". Default is "bar".
#' @param group_by The grouping variable(s) for the plot. Default is NULL.
#'   Only work when \code{plot_type} is "bar" or "comparison".
#'   For bar plot, bars will be colored by this variable.
#' @param pvalue_cutoff The p-value cutoff. Default is NULL. Only work when \code{padjust_cutoff} is NULL.
#' @param padjust_cutoff The p-adjusted cutoff. Default is 0.05.
#' @param top_term The number of top terms to display. Default is 6, or 100 if 'plot_type' is "enrichmap".
#' @param bar_showcutoff Whether to show the p-value or p-adjusted cutoff line. Default is TRUE.
#' @param compare_only_sig Whether to compare only significant terms. Default is FALSE.
#' @param top_word The number of top words to display for wordcloud. Default is 100.
#' @param word_type The type of words to display in wordcloud. Options are "term" and "feature". Default is "term".
#' @param word_size The size range for words in wordcloud. Default is c(2, 8).
#' @param words_excluded Words to be excluded from the wordcloud. The default value is NULL, which means that the built-in words (scplotter::words_excluded) will be used.
#' @param network_layout The layout algorithm to use for network plot. Options are "fr", "kk","random", "circle", "tree", "grid", or other algorithm from 'igraph' package. Default is "fr".
#' @param network_labelsize The label size for network plot. Default is 5.
#' @param network_blendmode The blend mode for network plot. Default is "blend".
#' @param network_layoutadjust Whether to adjust the layout of the network plot to avoid overlapping words. Default is TRUE.
#' @param network_adjscale The scale for adjusting network plot layout. Default is 60.
#' @param network_adjiter The number of iterations for adjusting network plot layout. Default is 100.
#' @param enrichmap_layout The layout algorithm to use for enrichmap plot. Options are "fr", "kk","random", "circle", "tree", "grid", or other algorithm from 'igraph' package. Default is "fr".
#' @param enrichmap_cluster The clustering algorithm to use for enrichmap plot. Options are "walktrap", "fast_greedy", or other algorithm from 'igraph' package. Default is "fast_greedy".
#' @param enrichmap_label  The label type for enrichmap plot. Options are "term" and "feature". Default is "term".
#' @param enrichmap_labelsize The label size for enrichmap plot. Default is 5.
#' @param enrlichmap_nlabel The number of labels to display for each cluster in enrichmap plot. Default is 4.
#' @param enrichmap_show_keyword Whether to show the keyword of terms or features in enrichmap plot. Default is FALSE.
#' @param enrichmap_mark The mark shape for enrichmap plot. Options are "ellipse" and "hull". Default is "ellipse".
#' @param enrichmap_expand The expansion factor for enrichmap plot. Default is c(0.5, 0.5).
#' @param character_width  The maximum width of character of descriptions. Default is 50.
#' @param lineheight The line height for y-axis labels. Default is 0.5.
#' @return A ggplot object or a list with the plot and the height and width of the plot if guess_size is TRUE
#' @export
#' @importFrom patchwork wrap_plots
#' @examples
#' data(enrichr_example)
#'
#' EnrichmentPlot(enrichr_example,
#'     db = "KEGG", plot_type = "bar",
#'     theme_args = list(base_size = 10)
#' )
#'
#' # Coloring by group
#' EnrichmentPlot(enrichr_example,
#'     db = "KEGG", plot_type = "bar",
#'     group_by = "Group", theme_args = list(base_size = 10)
#' )
#'
#' # Using a different palette
#' EnrichmentPlot(enrichr_example,
#'     db = "KEGG", plot_type = "bar", palette = "Paired",
#'     group_by = "Group", theme_args = list(base_size = 10)
#' )
#'
#' # Multiple databases
#' data(enrichr_multidb_example)
#'
#' EnrichmentPlot(enrichr_multidb_example, plot_type = "bar", facet = TRUE, nrow = 2)
#'
#' # Using `split_by` and combining using `patchwork::wrap_plots`
#' EnrichmentPlot(enrichr_multidb_example, plot_type = "bar", facet = FALSE, nrow = 2)
#'
#' # Comparison plot
#' EnrichmentPlot(enrichr_example,
#'     db = "KEGG", plot_type = "comparison",
#'     group_by = "Group"
#' )
#'
#' # Dot plot
#' EnrichmentPlot(enrichr_example,
#'     db = "KEGG", plot_type = "dot",
#'     theme_args = list(base_size = 10)
#' )
#'
#' # Lollipop plot
#' EnrichmentPlot(enrichr_example,
#'     db = "KEGG", plot_type = "lollipop",
#'     theme_args = list(base_size = 10)
#' )
#'
#' # Network
#' EnrichmentPlot(enrichr_example, db = "KEGG", plot_type = "network")
#'
#' # Enrichmap
#' EnrichmentPlot(enrichr_example, db = "KEGG", plot_type = "enrichmap")
#'
#' # Wordcloud
#' EnrichmentPlot(enrichr_example, db = "KEGG", plot_type = "wordcloud")
#'
#' # Wordcloud with feature
#' EnrichmentPlot(enrichr_example,
#'     db = "KEGG", plot_type = "wordcloud",
#'     word_type = "feature"
#' )
EnrichmentPlot <- function(
    x,
    enrichtool = NULL, db = NULL,
    plot_type = c("bar", "dot", "lollipop", "network", "enrichmap", "wordcloud", "comparison"),
    group_by = NULL, group_by_sep = "_", split_by = "Database", split_by_sep = "_",
    pvalue_cutoff = NULL, padjust_cutoff = 0.05, bar_showcutoff = TRUE,
    top_term = ifelse(plot_type == "enrichmap", 100, 10), compare_only_sig = FALSE,
    top_word = 100, word_type = c("term", "feature"), word_size = c(2, 8), words_excluded = scplotter::words_excluded,
    network_layout = "fr", network_labelsize = 5, network_blendmode = "blend",
    network_layoutadjust = TRUE, network_adjscale = 60, network_adjiter = 100,
    enrichmap_layout = "fr", enrichmap_cluster = "fast_greedy", enrichmap_label = "term",
    enrichmap_labelsize = 5, enrichmap_nlabel = 4, enrichmap_show_keyword = FALSE, enrlichmap_nlabel = 4,
    enrichmap_mark = "ellipse", enrichmap_expand = c(0.5, 0.5), x_text_angle = 0,
    character_width = 50, lineheight = 0.8, palette = "Spectral", palcolor = NULL,
    aspect.ratio = 1, legend.position = "right", legend.direction = "vertical",
    theme = "theme_scp", theme_args = list(), facet = FALSE, facet_scales = "fixed",
    nrow = NULL, ncol = NULL, byrow = TRUE, seed = 8525,
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
    res = 100, guess_size = FALSE, combine = TRUE,
    ...) {
    validate_common_arguments(
        split_by = split_by, seed = seed, res = res,
        guess_size = guess_size, facet = facet
    )

    plot_type <- match.arg(plot_type)
    word_type <- match.arg(word_type)
    enrichmap_label <- match.arg(enrichmap_label)
    enrichmap_mark <- match.arg(enrichmap_mark)

    enrichtool <- enrichtool %||% guess_enrichment_tool(x)
    class(x) <- c(paste0("enrichment_", enrichtool), class(x))
    x <- normalize_enrichment(x, db)

    if (is.null(db) && is.null(x$Database)) {
        stop("The 'db' argument is required when the 'Database' column is not present in the enrichment results.")
    }
    if (is.null(x$Database) && length(db) > 1) {
        stop("The 'db' argument should have only one value when the 'Database' column is not present in the enrichment results.")
    }
    db <- db %||% unique(x$Database)
    x$Database <- x$Database %||% db[1]

    if (is.null(pvalue_cutoff) && is.null(padjust_cutoff)) {
        stop("One of 'pvalue_cutoff' or 'padjust_cutoff' must be specified")
    }

    if (!is.factor(x$Database)) {
        x$Database <- factor(x$Database, levels = unique(x$Database))
    }

    metric <- ifelse(is.null(padjust_cutoff), "pvalue", "p.adjust")
    metric_value <- ifelse(is.null(padjust_cutoff), pvalue_cutoff, padjust_cutoff)
    x <- x[x[[metric]] < metric_value, , drop = FALSE]
    if (nrow(x) == 0) {
        stop(
            "\nNo significant terms enriched using the threshold:\n",
            paste0("- pvalue_cutoff = ", pvalue_cutoff), "\n",
            paste0("- padjust_cutoff = ", padjust_cutoff)
        )
    }
    enrichment_dbs <- unique(x$Database)
    x <- x[x$Database %in% db, , drop = FALSE]
    if (nrow(x) == 0) {
        stop(
            "\nHave you specified the correct 'db' value? \n",
            "- Available databases: ", paste(enrichment_dbs, collapse = ", "),
            "\n- Your 'db' value: ", db
        )
    }
    x <- x[order(x[[metric]], decreasing = FALSE), , drop = FALSE]
    x$Term <- factor(x$Term, levels = unique(x$Term))


    if (length(setdiff(split_by, colnames(x))) > 0) {
        stop(
            "The 'split_by' argument must be one or more of the columns in the enrichment results. ",
            "Not found: ", paste(setdiff(split_by, colnames(x)), collapse = ", ")
        )
    }
    if (!is.null(group_by) && (length(group_by) > 1 || !group_by %in% colnames(x))) {
        stop(
            "The 'group_by' argument must be one of the columns in the enrichment results. ",
            "Not found: ", group_by
        )
    }

    if (!is.null(split_by) && isFALSE(facet)) {
        df_list <- split(x, formula(paste0("~", split_by, collapse = "+")))
        df_list <- df_list[lapply(df_list, nrow) > 0]
    } else {
        df_list <- list(x)
    }

    plots <- lapply(df_list, function(df) {
        class(df) <- c(plot_type, class(df))
        suppressWarnings(EnrichmentPlotAtomic(
            df, metric = metric, metric_value = metric_value, enrichtool = enrichtool, db = db, plot_type = plot_type,
            group_by = group_by, group_by_sep = group_by_sep, split_by = split_by, split_by_sep = split_by_sep,
            pvalue_cutoff = pvalue_cutoff, padjust_cutoff = padjust_cutoff, top_term = top_term,
            bar_showcutoff = bar_showcutoff, compare_only_sig = compare_only_sig, keep_empty = keep_empty,
            top_word = top_word, word_type = word_type, word_size = word_size, words_excluded = words_excluded,
            network_layout = network_layout, network_labelsize = network_labelsize, network_blendmode = network_blendmode,
            network_layoutadjust = network_layoutadjust, network_adjscale = network_adjscale, network_adjiter = network_adjiter,
            enrichmap_layout = enrichmap_layout, enrichmap_cluster = enrichmap_cluster, enrichmap_label = enrichmap_label, enrlichmap_nlabel = enrlichmap_nlabel,
            enrichmap_labelsize = enrichmap_labelsize, enrichmap_nlabel = enrichmap_nlabel, enrichmap_show_keyword = enrichmap_show_keyword,
            enrichmap_mark = enrichmap_mark, enrichmap_expand = enrichmap_expand,
            x_text_angle = x_text_angle, character_width = character_width, lineheight = lineheight,
            palette = palette, palcolor = palcolor, aspect.ratio = aspect.ratio,
            legend.position = legend.position, legend.direction = legend.direction,
            theme = theme, theme_args = theme_args,
            title = title, subtitle = subtitle, xlab = xlab, ylab = ylab,
            facet = facet, facet_scales = facet_scales,
            nrow = nrow, ncol = ncol, byrow = byrow,
            seed = seed, res = res, guess_size = guess_size,
            ...
        ))
    })

    combine_plots(plots, combine = combine, nrow = nrow, ncol = ncol, byrow = byrow)
}

#' Guess_enrichment_tool
#'
#' @description Guess the enrichment tool from the enrichment results.
#'
#' @param x A data frame with Enrichment results
#' @keywords internal
#' @return The enrichment tool used to generate the enrichment results.
guess_enrichment_tool <- function(x) {
    if (!inherits(x, "data.frame")) {
        stop("The enrichment results must be a data frame.")
    }
    enrichr_cols <- c("Term", "Overlap", "P.value", "Adjusted.P.value", "Genes")
    if (length(setdiff(enrichr_cols, colnames(x))) == 0) {
        return("enrichr")
    }

    clusterprofiler_cols <- c("ID", "Description", "Overlap", "BgRatio", "pvalue", "p.adjust", "geneID")
    if (length(setdiff(clusterprofiler_cols, colnames(x))) == 0) {
        return("clusterprofiler")
    }

    return("unknown")
}

#' Normalize_enrichment
#'
#' @description Normalize the enrichment results.
#'
#' @param x A data frame with Enrichment results
#' @keywords internal
#' @return The normalized enrichment results.
normalize_enrichment <- function(x, db, ...) UseMethod("normalize_enrichment")

#' Default method for `normalize_enrichment`
#' @keywords internal
normalize_enrichment.default <- function(x, db, ...) {
    return(x)
}

#' Normalize_enrichment for "enrichment_enrichr"
#'
#' @keywords internal
normalize_enrichment.enrichment_enrichr <- function(x, db, ...) {
    x$pvalue <- x$P.value
    x$p.adjust <- x$Adjusted.P.value
    x$Overlap <- x$Overlap
    x$EnrichmentScore <- x$Combined.Score
    x$GeneRatio <- sapply(x$Overlap, function(x) {
        sp <- strsplit(x, "/")[[1]]
        as.numeric(sp[1]) / as.numeric(sp[2])
    })
    x$Count <- sapply(x$Overlap, function(x) {
        sp <- strsplit(x, "/")[[1]]
        as.numeric(sp[1])
    })
    if (!is.null(x$Database)) {
        dfs <- split(x, x$Database)
        x <- do.call(rbind, lapply(dfs, function(df) {
            df$ID <- paste0(df$Database[1], "_", seq_len(nrow(df)))
            df
        }))
    } else {
        x$ID <- paste0(db, "_", seq_len(nrow(x)))
    }
    x
}

#' Normalize_enrichment for "enrichment_clusterprofiler"
#'
#' @keywords internal
normalize_enrichment.enrichment_clusterprofiler <- function(x, db, ...) {
    x$Term <- x$Description
    x$Overlap <- x$GeneRatio
    x$Genes <- sapply(x$geneID, function(x) {
        gsub("/", ";", x, fixed = TRUE)
    })
    x$BgRatio <- sapply(x$BgRatio, function(x) {
        sp <- strsplit(x, "/")[[1]]
        as.numeric(sp[1]) / as.numeric(sp[2])
    })
    x$GeneRatio <- sapply(x$GeneRatio, function(x) {
        sp <- strsplit(x, "/")[[1]]
        as.numeric(sp[1]) / as.numeric(sp[2])
    })
    x$Odds.Ratio <- x$GeneRatio / x$BgRatio
    x
}

#' EnrichmentPlotAtomic
#'
#' @keywords internal
EnrichmentPlotAtomic <- function(x, ...) UseMethod("EnrichmentPlotAtomic")

#' EnrichmentPlotAtomic for default
#'
#' @keywords internal
EnrichmentPlotAtomic.default <- function(x, ...) {
    stop("Don't know how to generate the plot for ", class(x))
}

#' EnrichmentPlotAtomic for "bar" plot
#'
#' @inheritParams EnrichmentPlot
#' @param metric The metric to use for coloring the bars. Default is "p.adjust".
#' @param metric_value The metric value to use for coloring the bars. Default is 0.05.
#' @importFrom dplyr %>% slice_min
#' @importFrom rlang sym
#' @importFrom gglogger ggplot
#' @importFrom ggplot2 aes geom_bar geom_text geom_hline annotate labs scale_fill_manual scale_y_continuous facet_wrap coord_flip expansion element_line wrap_dims element_text ggplot_build
#' @importFrom stringr str_wrap
#' @keywords internal
EnrichmentPlotAtomic.bar <- function(
    x,
    enrichtool = NULL, db = NULL, metric = "p.adjust", metric_value = 0.05,
    group_by = NULL, group_by_sep = "_", split_by = "Database", split_by_sep = "_",
    pvalue_cutoff = NULL, padjust_cutoff = 0.05, bar_showcutoff = TRUE,
    top_term = 10, compare_only_sig = FALSE, keep_empty = FALSE,
    top_word = 100, word_type = c("term", "feature"), word_size = c(2, 8), words_excluded = scplotter::words_excluded,
    network_layout = "fr", network_labelsize = 5, network_blendmode = "blend",
    network_layoutadjust = TRUE, network_adjscale = 60, network_adjiter = 100,
    enrichmap_layout = "fr", enrichmap_cluster = "fast_greedy", enrichmap_label = "term",
    enrichmap_labelsize = 5, enrichmap_nlabel = 4, enrichmap_show_keyword = FALSE,
    enrichmap_mark = "ellipse", enrichmap_expand = c(0.5, 0.5), x_text_angle = 0,
    character_width = 50, lineheight = 0.8, palette = "Spectral", palcolor = NULL,
    aspect.ratio = 1, legend.position = "right", legend.direction = "vertical",
    theme = "theme_scp", theme_args = list(), facet = FALSE, facet_scales = "fixed",
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
    nrow = NULL, ncol = NULL, byrow = TRUE, seed = 8525, res = 100, guess_size = FALSE,
    ...) {
    base_size <- theme_args$base_size %||% 12
    text_size_scale <- base_size / 12
    top_term <- min(top_term, max(unlist(lapply(split(x, x$Database), nrow))))
    x <- x %>% slice_min(
        n = top_term,
        by = Database,
        order_by = !!sym(metric),
        with_ties = FALSE
    )
    x$metric <- -log10(x[[metric]])
    x$Term <- str_wrap(x$Term, width = character_width)
    x$Term <- factor(x$Term, levels = unique(rev(x$Term)))
    if (is.null(group_by)) {
        x$group_by <- "1"
        x$group_by <- factor(x$group_by, levels = unique(x$group_by))
        group_by <- "group_by"
        guide <- "none"
    } else if (!is.factor(x[[group_by]])) {
        x[[group_by]] <- factor(x[[group_by]], levels = unique(x[[group_by]]))
        guide <- "legend"
    }

    p <- ggplot(
        x,
        aes(x = Term, y = metric, fill = !!sym(group_by), label = Overlap)
    ) +
        geom_bar(width = 0.9, stat = "identity", color = "black") +
        geom_text(hjust = -0.5, size = text_size_scale * 11, size.unit = "pt")

    if (bar_showcutoff) {
        p <- p + geom_hline(yintercept = -log10(metric_value), linetype = 2, color = "#848484") +
            # add metric value text
            annotate(
                "text",
                x = Inf, y = -log10(metric_value), label = paste0(metric, " = ", metric_value),
                size = text_size_scale * 10, size.unit = "pt",
                color = "#848484", hjust = 1, vjust = -1
            )
    }

    just <- calc_just(x_text_angle)
    p <- p +
        scale_fill_manual(
            values = palette_scp(levels(x[[group_by]]), palette = palette, palcolor = palcolor),
            na.value = "grey80",
            guide = guide
        ) +
        scale_y_continuous(limits = c(0, 1.3 * max(x$metric, na.rm = TRUE)), expand = expansion(0, 0)) +
        coord_flip(clip = "off") +
        labs(title = title, subtitle = subtitle, x = xlab %||% "", y = ylab %||% paste0("-log10(", metric, ")")) +
        do.call(theme, theme_args) +
        ggplot2::theme(
            aspect.ratio = aspect.ratio,
            legend.position = legend.position,
            legend.direction = legend.direction,
            panel.grid.major = element_line(colour = "grey80", linetype = 2),
            axis.text.x = element_text(angle = x_text_angle, hjust = just$h, vjust = just$v),
            axis.text.y = element_text(
                lineheight = lineheight, hjust = 1,
                face = ifelse(any(grepl("\n", levels(x$Term))), "italic", "plain")
            )
        )

    if (bar_showcutoff) {
        p <- p + theme(plot.margin = margin(t = 24, r = 10, b = 10, l = 10))
    }

    if (isTRUE(guess_size)) {
        p <- list(
            plot = p,
            height = (top_term * .5 + .8) * res,
            width = (5 + (character_width * .04)) * res
        )
    }
    facet_plot(p, facet, facet_scales, split_by, nrow, ncol, byrow)
}

#' EnrichmentPlotAtomic for "comparison" plot
#'
#' @inheritParams EnrichmentPlot
#' @param metric The metric to use for coloring the bars. Default is "p.adjust".
#' @param metric_value The metric value to use for coloring the bars. Default is 0.05.
#' @importFrom rlang sym
#' @importFrom gglogger ggplot
#' @importFrom ggplot2 aes geom_point scale_size_area scale_fill_gradientn labs scale_color_manual guides guide_none facet_wrap element_line element_text guide_legend guide_colorbar wrap_dims
#' @importFrom stringr str_wrap
#' @keywords internal
EnrichmentPlotAtomic.comparison <- function(
    x,
    enrichtool = NULL, db = NULL, metric = "p.adjust", metric_value = 0.05,
    group_by = NULL, group_by_sep = "_", split_by = "Database", split_by_sep = "_",
    pvalue_cutoff = NULL, padjust_cutoff = 0.05, bar_showcutoff = TRUE,
    top_term = 10, compare_only_sig = FALSE,
    top_word = 100, word_type = c("term", "feature"), word_size = c(2, 8), words_excluded = scplotter::words_excluded,
    network_layout = "fr", network_labelsize = 5, network_blendmode = "blend",
    network_layoutadjust = TRUE, network_adjscale = 60, network_adjiter = 100,
    enrichmap_layout = "fr", enrichmap_cluster = "fast_greedy", enrichmap_label = "term",
    enrichmap_labelsize = 5, enrichmap_nlabel = 4, enrichmap_show_keyword = FALSE,
    enrichmap_mark = "ellipse", enrichmap_expand = c(0.5, 0.5), x_text_angle = 0,
    character_width = 50, lineheight = 0.8, palette = "Spectral", palcolor = NULL,
    aspect.ratio = 1, legend.position = "right", legend.direction = "vertical",
    theme = "theme_scp", theme_args = list(), facet = FALSE, facet_scales = "fixed",
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
    nrow = NULL, ncol = NULL, byrow = TRUE, seed = 8525, res = 100, guess_size = FALSE,
    ...) {
    if (is.null(group_by)) {
        stop("The 'group_by' argument is required for the 'comparison' plot (providing grouping on x-axis).")
    }
    df_groups <- split(x, list(x$Database, x[[group_by]]))
    df_groups <- lapply(df_groups, function(group) {
        group[head(seq_len(nrow(group)), top_term), , drop = FALSE]
    })
    enrichment_sub <- do.call(rbind, df_groups)

    enrichment_sub$Term <- str_wrap(enrichment_sub$Term, width = character_width)
    enrichment_sub$Term <- factor(enrichment_sub$Term, levels = unique(rev(enrichment_sub$Term)))
    if (isTRUE(compare_only_sig)) {
        enrichment_sub <- enrichment_sub[enrichment_sub[[metric]] < metric_value, , drop = FALSE]
    }
    ngroups <- length(unique(enrichment_sub[[group_by]]))
    nterms <- length(unique(enrichment_sub$Term))
    just <- calc_just(x_text_angle)
    p <- ggplot(enrichment_sub, aes(x = !!sym(group_by), y = Term)) +
        geom_point(aes(size = GeneRatio, fill = !!sym(metric), color = ""), shape = 21) +
        scale_size_area(name = "GeneRatio", max_size = 6, n.breaks = 4) +
        guides(size = guide_legend(override.aes = list(fill = "grey30", shape = 21), order = 1)) +
        scale_fill_gradientn(
            name = metric,
            limits = c(0, min(metric_value, 1)),
            n.breaks = 3,
            colors = palette_scp(palette = palette, palcolor = palcolor, reverse = TRUE),
            na.value = "grey80",
            guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0, order = 2)
        ) +
        labs(title = title, subtitle = subtitle, x = xlab %||% "", y = ylab %||% "") +
        scale_color_manual(values = NA, na.value = "black") +
        guides(colour = if (isTRUE(compare_only_sig)) guide_none() else guide_legend("Non-signif", override.aes = list(colour = "black", fill = "grey80", size = 3))) +
        do.call(theme, theme_args) +
        ggplot2::theme(
            aspect.ratio = aspect.ratio,
            legend.position = legend.position,
            legend.direction = legend.direction,
            panel.grid.major = element_line(colour = "grey80", linetype = 2),
            axis.text.x = element_text(angle = x_text_angle, hjust = just$h, vjust = just$v),
            axis.text.y = element_text(
                lineheight = lineheight, hjust = 1,
                face = ifelse(any(grepl("\n", levels(enrichment_sub$Term))), "italic", "plain")
            )
        )

    if (isTRUE(guess_size)) {
        p <- list(
            plot = p,
            height = (nterms * .3 + 1.4) * res,
            width = (1.5 * ngroups + (character_width * .15) + 1) * res
        )
    }
    facet_plot(p, facet, facet_scales, split_by, nrow, ncol, byrow)
}

#' EnrichmentPlotAtomic for "dot" plot
#'
#' @inheritParams EnrichmentPlot
#' @param metric The metric to use for coloring the dots. Default is "p.adjust".
#' @param metric_value The metric value to use for coloring the dots. Default is 0.05.
#' @importFrom gglogger ggplot
#' @importFrom ggplot2 aes geom_point labs scale_size scale_fill_gradientn scale_y_continuous coord_flip facet_wrap element_line element_text expansion guide_legend guide_colorbar wrap_dims
#' @importFrom stringr str_wrap
#' @importFrom scales breaks_extended
#' @keywords internal
EnrichmentPlotAtomic.dot <- function(
    x,
    enrichtool = NULL, db = NULL, metric = "p.adjust", metric_value = 0.05,
    group_by = NULL, group_by_sep = "_", split_by = "Database", split_by_sep = "_",
    pvalue_cutoff = NULL, padjust_cutoff = 0.05, bar_showcutoff = TRUE,
    top_term = 10, compare_only_sig = FALSE,
    top_word = 100, word_type = c("term", "feature"), word_size = c(2, 8), words_excluded = scplotter::words_excluded,
    network_layout = "fr", network_labelsize = 5, network_blendmode = "blend",
    network_layoutadjust = TRUE, network_adjscale = 60, network_adjiter = 100,
    enrichmap_layout = "fr", enrichmap_cluster = "fast_greedy", enrichmap_label = "term",
    enrichmap_labelsize = 5, enrichmap_nlabel = 4, enrichmap_show_keyword = FALSE,
    enrichmap_mark = "ellipse", enrichmap_expand = c(0.5, 0.5), x_text_angle = 0,
    character_width = 50, lineheight = 0.8, palette = "Spectral", palcolor = NULL,
    aspect.ratio = 1, legend.position = "right", legend.direction = "vertical",
    theme = "theme_scp", theme_args = list(), facet = FALSE, facet_scales = "fixed",
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
    nrow = NULL, ncol = NULL, byrow = TRUE, seed = 8525, res = 100, guess_size = FALSE,
    ...) {
    df_groups <- split(x, list(x$Database))
    df_groups <- lapply(df_groups, function(group) {
        group[head(seq_len(nrow(group)), top_term), , drop = FALSE]
    })
    df <- do.call(rbind, df_groups)

    df <- df[order(df$GeneRatio, decreasing = TRUE), ]
    df$metric <- -log10(df[[metric]])
    df$Term <- str_wrap(df$Term, width = character_width)
    df$Term <- factor(df$Term, levels = unique(rev(df$Term)))

    nterms <- length(unique(df$Term))
    p <- ggplot(df, aes(x = Term, y = GeneRatio)) +
        geom_point(aes(fill = metric, size = Count), color = "black", shape = 21) +
        labs(title = title, subtitle = subtitle, x = xlab %||% "", y = ylab %||% "GeneRatio") +
        scale_size(name = "Count", range = c(3, 6), breaks_extended(n = 4)) +
        guides(size = guide_legend(override.aes = list(fill = "grey30", shape = 21), order = 1)) +
        scale_fill_gradientn(
            name = paste0("-log10(", metric, ")"),
            n.breaks = 3,
            colors = palette_scp(palette = palette, palcolor = palcolor),
            na.value = "grey80",
            guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0)
        ) +
        scale_y_continuous(limits = c(0, 1.3 * max(df$GeneRatio, na.rm = TRUE)), expand = expansion(0, 0)) +
        coord_flip() +
        do.call(theme, theme_args) +
        ggplot2::theme(
            aspect.ratio = aspect.ratio,
            legend.position = legend.position,
            legend.direction = legend.direction,
            panel.grid.major = element_line(colour = "grey80", linetype = 2),
            axis.text.y = element_text(
                lineheight = lineheight, hjust = 1,
                face = ifelse(any(grepl("\n", levels(df$Term))), "italic", "plain")
            )
        )

    if (isTRUE(guess_size)) {
        p <- list(
            plot = p,
            height = (nterms * .3 + 1.4) * res,
            width = (4 + (character_width * .08) + 1) * res
        )
    }
    facet_plot(p, facet, facet_scales, split_by, nrow, ncol, byrow)
}

#' EnrichmentPlotAtomic for "lollipop" plot
#' @inheritParams EnrichmentPlot
#' @param metric The metric to use for coloring the lollipops. Default is "p.adjust".
#' @param metric_value The metric value to use for coloring the lollipops. Default is 0.05.
#' @importFrom gglogger ggplot
#' @importFrom ggplot2 aes geom_blank geom_segment geom_point scale_size scale_fill_gradientn scale_y_continuous coord_flip facet_wrap element_line element_text expansion guide_legend guide_colorbar wrap_dims
#' @importFrom stringr str_wrap
#' @importFrom scales breaks_extended
#' @keywords internal
EnrichmentPlotAtomic.lollipop <- function(
    x,
    enrichtool = NULL, db = NULL, metric = "p.adjust", metric_value = 0.05,
    group_by = NULL, group_by_sep = "_", split_by = "Database", split_by_sep = "_",
    pvalue_cutoff = NULL, padjust_cutoff = 0.05, bar_showcutoff = TRUE,
    top_term = 10, compare_only_sig = FALSE,
    top_word = 100, word_type = c("term", "feature"), word_size = c(2, 8), words_excluded = scplotter::words_excluded,
    network_layout = "fr", network_labelsize = 5, network_blendmode = "blend",
    network_layoutadjust = TRUE, network_adjscale = 60, network_adjiter = 100,
    enrichmap_layout = "fr", enrichmap_cluster = "fast_greedy", enrichmap_label = "term",
    enrichmap_labelsize = 5, enrichmap_nlabel = 4, enrichmap_show_keyword = FALSE,
    enrichmap_mark = "ellipse", enrichmap_expand = c(0.5, 0.5), x_text_angle = 0,
    character_width = 50, lineheight = 0.8, palette = "Spectral", palcolor = NULL,
    aspect.ratio = 1, legend.position = "right", legend.direction = "vertical",
    theme = "theme_scp", theme_args = list(), facet = FALSE, facet_scales = "fixed",
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
    nrow = NULL, ncol = NULL, byrow = TRUE, seed = 8525, res = 100, guess_size = FALSE,
    ...) {
    df_groups <- split(x, list(x$Database))
    df_groups <- lapply(df_groups, function(group) {
        group[head(seq_len(nrow(group)), top_term), , drop = FALSE]
    })
    df <- do.call(rbind, df_groups)

    df$metric <- -log10(df[[metric]])
    df$Term <- str_wrap(df$Term, width = character_width)
    df$Term <- factor(df$Term, levels = unique(df[order(df[["Odds.Ratio"]]), "Term"]))
    nterms <- length(unique(df$Term))
    p <- ggplot(df, aes(x = Term, y = Odds.Ratio, fill = metric)) +
        geom_blank() +
        geom_segment(
            aes(y = 0, xend = Term, yend = Odds.Ratio),
            color = "black", linewidth = 2
        ) +
        geom_segment(
            aes(y = 0, xend = Term, yend = Odds.Ratio, color = metric),
            linewidth = 1
        ) +
        geom_point(aes(size = GeneRatio), shape = 21, color = "black") +
        scale_size(name = "GeneRatio", range = c(3, 6), breaks_extended(n = 4)) +
        guides(size = guide_legend(override.aes = list(fill = "grey30", shape = 21), order = 1)) +
        scale_y_continuous(limits = c(0, 1.2 * max(df[["Odds.Ratio"]], na.rm = TRUE)), expand = expansion(0, 0)) +
        labs(title = title, subtitle = subtitle, x = xlab %||% "", y = ylab %||% "Odds Ratio") +
        scale_fill_gradientn(
            name = paste0("-log10(", metric, ")"),
            n.breaks = 3,
            colors = palette_scp(palette = palette, palcolor = palcolor),
            na.value = "grey80",
            guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0),
            aesthetics = c("color", "fill")
        ) +
        coord_flip() +
        do.call(theme, theme_args) +
        ggplot2::theme(
            aspect.ratio = aspect.ratio,
            legend.position = legend.position,
            legend.direction = legend.direction,
            panel.grid.major = element_line(colour = "grey80", linetype = 2),
            axis.text.y = element_text(
                lineheight = lineheight, hjust = 1,
                face = ifelse(any(grepl("\n", levels(df$Term))), "italic", "plain")
            )
        )

    if (isTRUE(guess_size)) {
        p <- list(
            plot = p,
            height = (nterms * .3 + 1.4) * res,
            width = (4 + (character_width * .08) + 1) * res
        )
    }
    facet_plot(p, facet, facet_scales, split_by, nrow, ncol, byrow)
}

#' EnrichmentPlotAtomic for "network" plot
#'
#' @inheritParams EnrichmentPlot
#' @param metric The metric to use for coloring the nodes. Default is "p.adjust".
#' @param metric_value The metric value to use for coloring the nodes. Default is 0.05.
#' @importFrom tidyr unnest
#' @importFrom stringr str_wrap
#' @importFrom igraph graph_from_data_frame as_data_frame layout_in_circle layout_as_tree layout_on_grid layout_with_fr
#' @importFrom gglogger ggplot
#' @importFrom ggplot2 aes geom_segment geom_point geom_text scale_color_manual scale_fill_manual element_text element_line scale_color_identity scale_fill_identity .pt draw_key_point
#' @importFrom ggrepel geom_text_repel
#' @importFrom grid grobTree gpar
#' @keywords internal
EnrichmentPlotAtomic.network <- function(
    x,
    enrichtool = NULL, db = NULL, metric = "p.adjust", metric_value = 0.05,
    group_by = NULL, group_by_sep = "_", split_by = "Database", split_by_sep = "_",
    pvalue_cutoff = NULL, padjust_cutoff = 0.05, bar_showcutoff = TRUE,
    top_term = 10, compare_only_sig = FALSE,
    top_word = 100, word_type = c("term", "feature"), word_size = c(2, 8), words_excluded = scplotter::words_excluded,
    network_layout = "fr", network_labelsize = 5, network_blendmode = "blend",
    network_layoutadjust = TRUE, network_adjscale = 60, network_adjiter = 100,
    enrichmap_layout = "fr", enrichmap_cluster = "fast_greedy", enrichmap_label = "term",
    enrichmap_labelsize = 5, enrichmap_nlabel = 4, enrichmap_show_keyword = FALSE,
    enrichmap_mark = "ellipse", enrichmap_expand = c(0.5, 0.5), x_text_angle = 0,
    character_width = 50, lineheight = 0.8, palette = "Spectral", palcolor = NULL,
    aspect.ratio = 1, legend.position = "right", legend.direction = "vertical",
    theme = "theme_scp", theme_args = list(), facet = FALSE, facet_scales = "fixed",
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
    nrow = NULL, ncol = NULL, byrow = TRUE, seed = 8525, res = 100, guess_size = FALSE,
    ...) {
    n_dbs <- length(unique(x$Database))
    if (n_dbs > 1) {
        stop(
            "Can't create network plot with multiple databases. ",
            "Try using 'split_by = \"Database\"' to split the data."
        )
    }
    if (!is.null(split_by) && split_by != "Database" && isTRUE(facet)) {
        stop(
            "Can't facet the network plot. \n",
            "Try 'facet = FALSE' to split the data and separate the plots."
        )
    }

    df_groups <- split(x, list(x$Database))
    df_groups <- lapply(df_groups, function(group) {
        group[head(seq_len(nrow(group)), top_term), , drop = FALSE]
    })
    df <- do.call(rbind, df_groups)

    df$metric <- -log10(df[[metric]])
    df$Term <- str_wrap(df$Term, width = character_width)
    df$Term <- factor(df$Term, levels = unique(df$Term))
    df$Genes <- strsplit(df$Genes, ";")
    df_unnest <- unnest(df, cols = "Genes")

    nodes <- rbind(
        data.frame("ID" = df$Term, class = "term", metric = df$metric),
        data.frame("ID" = unique(df_unnest$Genes), class = "gene", metric = 0)
    )
    nodes$Database <- df$Database[1]
    edges <- as.data.frame(df_unnest[, c("Term", "Genes")])
    colnames(edges) <- c("from", "to")
    edges$weight <- 1
    graph <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
    if (network_layout %in% c("circle", "tree", "grid")) {
        layout <- switch(network_layout,
            "circle" = layout_in_circle(graph),
            "tree" = layout_as_tree(graph),
            "grid" = layout_on_grid(graph)
        )
    } else {
        layout <- do.call(paste0("layout_with_", network_layout), list(graph))
    }
    df_graph <- as_data_frame(graph, what = "both")

    df_nodes <- df_graph$vertices
    if (isTRUE(network_layoutadjust)) {
        width <- nchar(df_nodes$name)
        width[df_nodes$class == "term"] <- 8
        layout <- adjust_network_layout(
            graph = graph, layout = layout, width = width, height = 2,
            scale = network_adjscale, iter = network_adjiter
        )
    }
    df_nodes$dim1 <- layout[, 1]
    df_nodes$dim2 <- layout[, 2]

    df_edges <- df_graph$edges
    df_edges$from_dim1 <- df_nodes[df_edges$from, "dim1"]
    df_edges$from_dim2 <- df_nodes[df_edges$from, "dim2"]
    df_edges$to_dim1 <- df_nodes[df_edges$to, "dim1"]
    df_edges$to_dim2 <- df_nodes[df_edges$to, "dim2"]

    colors <- palette_scp(levels(df$Term), palette = palette, palcolor = palcolor, keep_names = TRUE)
    df_edges$color <- colors[df_edges$from]
    node_colors <- aggregate(
        df_unnest$Term,
        by = list(df_unnest$Genes),
        FUN = function(x) blend_colors(colors = colors[x], mode = network_blendmode)
    )
    colors <- c(colors, setNames(node_colors[, 2], node_colors[, 1]))
    label_colors <- ifelse(colSums(col2rgb(colors)) > 255 * 2, "black", "white")
    df_nodes$color <- colors[df_nodes$name]
    df_nodes$label_color <- label_colors[df_nodes$name]
    df_nodes$label <- NA
    df_nodes[levels(df$Term), "label"] <- seq_len(nlevels(df$Term))

    draw_key_cust <- function(data, params, size) {
        data_text <- data
        data_text$label <- which(levels(df$Term) %in% names(colors)[colors == data_text$fill])
        data_text$colour <- "black"
        data_text$alpha <- 1
        data_text$size <- 11 / .pt
        grobTree(
            draw_key_point(data, list(color = "white", shape = 21)),
            ggrepel:::shadowtextGrob(label = data_text$label, bg.colour = "black", bg.r = 0.1, gp = gpar(col = "white", fontface = "bold"))
        )
    }

    p <- ggplot() +
        geom_segment(data = df_edges, aes(x = from_dim1, y = from_dim2, xend = to_dim1, yend = to_dim2, color = color), alpha = 1, lineend = "round", show.legend = FALSE) +
        geom_label(data = df_nodes[df_nodes$class == "gene", ], aes(x = dim1, y = dim2, label = name, fill = color, color = label_color), size = 3, show.legend = FALSE) +
        geom_point(data = df_nodes[df_nodes$class == "term", ], aes(x = dim1, y = dim2), size = 8, color = "black", fill = "black", stroke = 1, shape = 21, show.legend = FALSE) +
        geom_point(data = df_nodes[df_nodes$class == "term", ], aes(x = dim1, y = dim2, fill = color), size = 7, color = "white", stroke = 1, shape = 21, key_glyph = draw_key_cust) +
        geom_text_repel(
            data = df_nodes[df_nodes$class == "term", ], aes(x = dim1, y = dim2, label = label),
            fontface = "bold", min.segment.length = 0, segment.color = "black",
            point.size = NA, max.overlaps = 100, force = 0, color = "white", bg.color = "black", bg.r = 0.1, size = network_labelsize
        ) +
        scale_color_identity(guide = "none") +
        scale_fill_identity(
            name = "Term:", guide = "legend",
            labels = levels(df$Term),
            breaks = colors[levels(df$Term)]
        ) +
        guides(color = guide_legend(override.aes = list(color = "transparent"))) +
        labs(title = title, subtitle = subtitle, x = xlab %||% "", y = ylab %||% "") +
        do.call(theme, theme_args) +
        ggplot2::theme(
            aspect.ratio = aspect.ratio,
            legend.position = legend.position,
            legend.direction = legend.direction
        )

    if (isTRUE(guess_size)) {
        p <- list(
            plot = p,
            height = 8 * res,
            width = 12 * res
        )
    }
    facet_plot(p, facet, facet_scales, split_by, nrow, ncol, byrow)
}

#' EnrichmentPlotAtomic for "enrichmap" plot
#'
#' @inheritParams EnrichmentPlot
#' @param metric The metric to use for coloring the nodes. Default is "p.adjust".
#' @param metric_value The metric value to use for coloring the nodes. Default is 0.05.
#' @importFrom rlang sym
#' @importFrom dplyr %>% group_by reframe slice_head distinct arrange mutate filter ungroup n
#' @importFrom tidyr unnest
#' @importFrom stringr str_wrap
#' @importFrom gglogger ggplot
#' @importFrom ggplot2 aes geom_point geom_label scale_color_manual scale_fill_manual element_text element_line scale_linewidth
#' @importFrom ggrepel geom_label_repel
#' @importFrom igraph as_data_frame cluster_fast_greedy cluster_leading_eigen cluster_edge_betweenness cluster_walktrap
#' @importFrom scales breaks_extended
#' @importFrom stringr str_wrap
#' @keywords internal
EnrichmentPlotAtomic.enrichmap <- function(
    x,
    enrichtool = NULL, db = NULL, metric = "p.adjust", metric_value = 0.05,
    group_by = NULL, group_by_sep = "_", split_by = "Database", split_by_sep = "_",
    pvalue_cutoff = NULL, padjust_cutoff = 0.05, bar_showcutoff = TRUE,
    top_term = 10, compare_only_sig = FALSE,
    top_word = 100, word_type = c("term", "feature"), word_size = c(2, 8), words_excluded = scplotter::words_excluded,
    network_layout = "fr", network_labelsize = 5, network_blendmode = "blend",
    network_layoutadjust = TRUE, network_adjscale = 60, network_adjiter = 100,
    enrichmap_layout = "fr", enrichmap_cluster = "fast_greedy", enrichmap_label = "term", enrlichmap_nlabel = 4,
    enrichmap_labelsize = 5, enrichmap_nlabel = 4, enrichmap_show_keyword = FALSE,
    enrichmap_mark = "ellipse", enrichmap_expand = c(0.5, 0.5), x_text_angle = 0,
    character_width = 50, lineheight = 0.8, palette = "Spectral", palcolor = NULL,
    aspect.ratio = 1, legend.position = "right", legend.direction = "vertical",
    theme = "theme_scp", theme_args = list(), facet = FALSE, facet_scales = "fixed",
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
    nrow = NULL, ncol = NULL, byrow = TRUE, seed = 8525, res = 100, guess_size = FALSE,
    ...) {
    n_dbs <- length(unique(x$Database))
    if (n_dbs > 1) {
        stop(
            "Can't create network plot with multiple databases. ",
            "Try using 'split_by = \"Database\"' to split the data."
        )
    }

    df_groups <- split(x, list(x$Database))
    df_groups <- lapply(df_groups, function(group) {
        group[head(seq_len(nrow(group)), top_term), , drop = FALSE]
    })
    df <- do.call(rbind, df_groups)

    df$metric <- -log10(df[[metric]])
    df$Term <- factor(df$Term, levels = unique(df$Term))
    df$Genes <- strsplit(df$Genes, ";")
    rownames(df) <- df$ID

    nodes <- df
    edges <- as.data.frame(t(combn(nodes$ID, 2)))
    colnames(edges) <- c("from", "to")
    edges$weight <- mapply(function(x, y) length(intersect(df[[x, "Genes"]], df[[y, "Genes"]])), edges$from, edges$to)
    edges <- edges[edges$weight > 0, , drop = FALSE]
    nodes <- nodes[c("ID", setdiff(colnames(nodes), "ID"))]
    graph <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
    if (enrichmap_layout %in% c("circle", "tree", "grid")) {
        layout <- switch(enrichmap_layout,
            "circle" = layout_in_circle(graph),
            "tree" = layout_as_tree(graph),
            "grid" = layout_on_grid(graph)
        )
    } else {
        layout <- do.call(paste0("layout_with_", enrichmap_layout), list(graph))
    }
    clusters <- do.call(paste0("cluster_", enrichmap_cluster), list(graph))
    df_graph <- as_data_frame(graph, what = "both")

    df_nodes <- df_graph$vertices
    df_nodes$dim1 <- layout[, 1]
    df_nodes$dim2 <- layout[, 2]
    df_nodes$clusters <- factor(paste0("C", clusters$membership), paste0("C", unique(sort(clusters$membership))))

    if (isTRUE(enrichmap_show_keyword)) {
        df_keyword1 <- df_nodes %>%
            mutate(keyword = strsplit(tolower(as.character(Term)), "\\s|\\n", perl = TRUE)) %>%
            unnest(cols = "keyword") %>%
            group_by(keyword, Database, clusters) %>%
            reframe(
                keyword = capitalize(keyword),
                score = sum(-(log10(!!sym(metric)))),
                count = n(),
                Database = Database,
                .groups = "keep"
            ) %>%
            filter(!grepl(pattern = "\\[.*\\]", x = keyword)) %>%
            filter(nchar(keyword) >= 1) %>%
            filter(!tolower(keyword) %in% tolower(words_excluded)) %>%
            distinct() %>%
            group_by(Database, clusters) %>%
            arrange(desc(score)) %>%
            slice_head(n = enrlichmap_nlabel) %>%
            reframe(keyword = paste0(keyword, collapse = " ")) %>%
            as.data.frame()
        rownames(df_keyword1) <- as.character(df_keyword1$clusters)
        df_keyword1$keyword <- str_wrap(df_keyword1$keyword, width = character_width)
        df_keyword1$label <- paste0(df_keyword1$clusters, ":\n", df_keyword1$keyword)
    } else {
        if (enrichmap_label == "term") {
            df_nodes$Term <- str_wrap(df_nodes$Term, width = character_width)
        }
        df_keyword1 <- df_nodes %>%
            group_by(Database, clusters) %>%
            arrange(desc(metric)) %>%
            reframe(keyword = Term) %>%
            distinct() %>%
            group_by(Database, clusters) %>%
            slice_head(n = enrlichmap_nlabel) %>%
            reframe(keyword = paste0(keyword, collapse = "\n")) %>%
            as.data.frame()

        rownames(df_keyword1) <- as.character(df_keyword1$clusters)
        df_keyword1$label <- paste0(df_keyword1$clusters, ":\n", df_keyword1$keyword)
    }

    df_keyword2 <- df_nodes %>%
        mutate(keyword = Genes) %>%
        unnest(cols = "keyword") %>%
        group_by(keyword, Database, clusters) %>%
        reframe(
            keyword = keyword,
            score = sum(-(log10(!!sym(metric)))),
            count = n(),
            Database = Database,
            .groups = "keep"
        ) %>%
        distinct() %>%
        group_by(Database, clusters) %>%
        arrange(desc(score)) %>%
        slice_head(n = enrlichmap_nlabel) %>%
        reframe(keyword = paste0(keyword, collapse = " ")) %>%
        as.data.frame()

    rownames(df_keyword2) <- as.character(df_keyword2$clusters)
    df_keyword2$keyword <- str_wrap(df_keyword2$keyword, width = character_width)
    df_keyword2$label <- paste0(df_keyword2$clusters, ":\n", df_keyword2$keyword)

    df_nodes$keyword1 <- df_keyword1[as.character(df_nodes$clusters), "keyword"]
    df_nodes$keyword2 <- df_keyword2[as.character(df_nodes$clusters), "keyword"]

    df_edges <- df_graph$edges
    df_edges$from_dim1 <- df_nodes[df_edges$from, "dim1"]
    df_edges$from_dim2 <- df_nodes[df_edges$from, "dim2"]
    df_edges$to_dim1 <- df_nodes[df_edges$to, "dim1"]
    df_edges$to_dim2 <- df_nodes[df_edges$to, "dim2"]

    if (enrichmap_mark == "hull") {
        check_R("concaveman")
    }
    mark_layer <- do.call(
        getFromNamespace(switch(enrichmap_mark,
            "ellipse" = "geom_mark_ellipse",
            "hull" = "geom_mark_hull"
        ), ns = "ggforce"),
        list(
            data = df_nodes, aes(
                x = dim1, y = dim2, color = clusters, fill = clusters,
                label = clusters, description = if (enrichmap_label == "term") keyword1 else keyword2
            ),
            expand = unit(3, "mm"),
            alpha = 0.1,
            label.margin = margin(1, 1, 1, 1, "mm"),
            label.fontsize = enrichmap_labelsize * 2,
            label.fill = "grey95",
            label.minwidth = unit(character_width, "in"),
            label.buffer = unit(0, "mm"),
            con.size = 1,
            con.cap = 0
        )
    )

    p <- ggplot() +
        mark_layer +
        geom_segment(data = df_edges, aes(x = from_dim1, y = from_dim2, xend = to_dim1, yend = to_dim2, linewidth = weight), alpha = 0.1, lineend = "round") +
        geom_point(data = df_nodes, aes(x = dim1, y = dim2, size = Count, fill = clusters), color = "black", shape = 21) +
        labs(title = title, subtitle = subtitle, x = xlab %||% "", y = ylab %||% "") +
        scale_size(name = "Count", range = c(2, 6), breaks_extended(n = 4)) +
        guides(size = guide_legend(override.aes = list(fill = "grey30", shape = 21), order = 1)) +
        scale_linewidth(name = "Intersection", range = c(0.3, 3), breaks_extended(n = 4)) +
        guides(linewidth = guide_legend(override.aes = list(alpha = 1, color = "grey"), order = 2)) +
        scale_fill_manual(
            name = switch(enrichmap_label,
                "term" = "Feature:",
                "feature" = "Term:"
            ),
            values = palette_scp(levels(df_nodes$clusters), palette = palette, palcolor = palcolor),
            labels = if (enrichmap_label == "term") df_keyword2[levels(df_nodes$clusters), "label"] else df_keyword1[levels(df_nodes$clusters), "label"],
            na.value = "grey80",
            aesthetics = c("colour", "fill")
        ) +
        guides(fill = guide_legend(override.aes = list(alpha = 1, color = "black", shape = NA), byrow = TRUE, order = 3)) +
        guides(color = guide_none()) +
        scale_x_continuous(expand = expansion(c(enrichmap_expand[1], enrichmap_expand[1]), 0)) +
        scale_y_continuous(expand = expansion(c(enrichmap_expand[2], enrichmap_expand[2]), 0)) +
        # facet_grid(facet, scales = "free") +
        do.call(theme, theme_args) +
        ggplot2::theme(
            aspect.ratio = aspect.ratio,
            legend.position = legend.position,
            legend.direction = legend.direction
        )

    if (isTRUE(guess_size)) {
        p <- list(
            plot = p,
            height = 8 * res,
            width = 12 * res
        )
    }
    facet_plot(p, facet, facet_scales, split_by, nrow, ncol, byrow)
}

#' EnrichmentPlotAtomic for "wordcloud" plot
#'
#' @inheritParams EnrichmentPlot
#' @param metric The metric to use for coloring the nodes. Default is "p.adjust".
#' @param metric_value The metric value to use for coloring the nodes. Default is 0.05.
#' @importFrom dplyr %>% reframe filter mutate group_by distinct
#' @importFrom rlang sym
#' @importFrom gglogger ggplot
#' @importFrom ggplot2 aes scale_color_gradientn scale_size guides guide_legend coord_flip element_line element_text
#' @importFrom scales rescale
#' @importFrom tidyr unnest
#' @importFrom ggwordcloud geom_text_wordcloud
#' @keywords internal
EnrichmentPlotAtomic.wordcloud <- function(
    x,
    enrichtool = NULL, db = NULL, metric = "p.adjust", metric_value = 0.05,
    group_by = NULL, group_by_sep = "_", split_by = "Database", split_by_sep = "_",
    pvalue_cutoff = NULL, padjust_cutoff = 0.05, bar_showcutoff = TRUE,
    top_term = 10, compare_only_sig = FALSE,
    top_word = 100, word_type = c("term", "feature"), word_size = c(2, 8), words_excluded = scplotter::words_excluded,
    network_layout = "fr", network_labelsize = 5, network_blendmode = "blend",
    network_layoutadjust = TRUE, network_adjscale = 60, network_adjiter = 100,
    enrichmap_layout = "fr", enrichmap_cluster = "fast_greedy", enrichmap_label = "term",
    enrichmap_labelsize = 5, enrichmap_nlabel = 4, enrichmap_show_keyword = FALSE,
    enrichmap_mark = "ellipse", enrichmap_expand = c(0.5, 0.5), x_text_angle = 0,
    character_width = 50, lineheight = 0.8, palette = "Spectral", palcolor = NULL,
    aspect.ratio = 1, legend.position = "right", legend.direction = "vertical",
    theme = "theme_scp", theme_args = list(), facet = FALSE, facet_scales = "fixed",
    title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
    nrow = NULL, ncol = NULL, byrow = TRUE, seed = 8525, res = 100, guess_size = FALSE,
    ...) {
    if (word_type == "term") {
        df_groups <- split(x, list(x$Database))
        df_groups <- df_groups[sapply(df_groups, nrow) > 0]
        for (i in seq_along(df_groups)) {
            df_sub <- df_groups[[i]] %>%
                mutate(keyword = strsplit(tolower(as.character(Term)), " ")) %>%
                unnest(cols = "keyword") %>%
                group_by(keyword, Database) %>%
                reframe(
                    keyword = keyword,
                    score = sum(-(log10(!!sym(metric)))),
                    count = n(),
                    Database = Database,
                    .groups = "keep"
                ) %>%
                filter(!grepl(pattern = "\\[.*\\]", x = keyword)) %>%
                filter(nchar(keyword) >= 2) %>%
                filter(!tolower(keyword) %in% tolower(words_excluded)) %>%
                distinct() %>%
                mutate(angle = 90 * sample(c(0, 1), n(), replace = TRUE, prob = c(60, 40))) %>%
                as.data.frame()

            df_sub <- df_sub[head(order(df_sub[["score"]], decreasing = TRUE), top_word), , drop = FALSE]
            df_groups[[i]] <- df_sub
        }
        df <- do.call(rbind, df_groups)
    } else {
        df <- x %>%
            mutate(keyword = strsplit(as.character(Genes), ";")) %>%
            unnest(cols = "keyword") %>%
            group_by(keyword, Database) %>%
            reframe(
                keyword = keyword,
                score = sum(-(log10(!!sym(metric)))),
                count = n(),
                Database = Database,
                .groups = "keep"
            ) %>%
            distinct() %>%
            mutate(angle = 90 * sample(c(0, 1), n(), replace = TRUE, prob = c(60, 40))) %>%
            as.data.frame()
        df <- df[head(order(df[["score"]], decreasing = TRUE), top_word), , drop = FALSE]
    }

    colors <- palette_scp(df[["score"]], type = "continuous", palette = palette, palcolor = palcolor, matched = FALSE)
    colors_value <- seq(min(df[["score"]], na.rm = TRUE), quantile(df[["score"]], 0.99, na.rm = TRUE) + 0.001, length.out = 100)
    p <- ggplot(df, aes(label = keyword, size = count, color = score, angle = angle)) +
        geom_text_wordcloud(rm_outside = TRUE, eccentricity = 1, shape = "square", show.legend = TRUE, grid_margin = 3) +
        scale_color_gradientn(
            name = "Score:", colours = colors, values = rescale(colors_value),
            guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0)
        ) +
        scale_size(name = "Count", range = word_size, breaks = ceiling(seq(min(df[["count"]], na.rm = TRUE), max(df[["count"]], na.rm = TRUE), length.out = 3))) +
        guides(size = guide_legend(override.aes = list(colour = "black", label = "G"), order = 1)) +
        coord_flip() +
        do.call(theme, theme_args) +
        ggplot2::theme(
            aspect.ratio = aspect.ratio,
            legend.position = legend.position,
            legend.direction = legend.direction
        )

    if (isTRUE(guess_size)) {
        p <- list(
            plot = p,
            height = 6 * res,
            width = 8 * res
        )
    }
    facet_plot(p, facet, facet_scales, split_by, nrow, ncol, byrow)
}
