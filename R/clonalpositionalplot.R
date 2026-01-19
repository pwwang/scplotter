#' ClonalPositionalPlot
#'
#' Visualize the positional entropy, property or amino acid frequency of CDR3 sequences.
#'
#' @param data The product of [scRepertoire::combineTCR], [scRepertoire::combineTCR], or
#'  [scRepertoire::combineExpression].
#' @param chain The chain to be analyzed. Default is "TRB".
#' @param aa_length The length of the amino acid sequence. Default is 20.
#' @param group_by The variable to group the data by. Default is "Sample".
#' @param group_by_sep The separator to use when combining groupings. Default is "_".
#' @param split_by The variable to split the data by. Default is NULL.
#' @param method The method to calculate the positional entropy. Default is "AA".
#'  * "AA": Amino acid frequency.
#'  * "shannon": Shannon entropy.
#'  * "inv.simpson": Inverse Simpson index.
#'  * "norm.entropy": Normalized entropy.
#'  * "Atchley": Atchley factors.
#'  * "Kidera": Kidera factors.
#'  * "stScales": stScales factors.
#'  * "tScales": tScales factors.
#'  * "VHSE": Vectors of Hydrophobic, Steric, and Electronic properties.
#'  See also [scRepertoire::percentAA], [scRepertoire::positionalEntropy] and
#'  [scRepertoire::positionalProperty].
#' @param plot_type The type of plot to generate. Default is "bar".
#'  * "bar": Bar plot.
#'  * "line": Line plot.
#'  * "heatmap": Heatmap.
#'  * "box": Box plot.
#'  * "violin": Violin plot.
#' @param theme_args A list of arguments to be passed to the [ggplot2::theme] function.
#' @param xlab The x-axis label. Default is NULL.
#' @param ylab The y-axis label. Default is NULL.
#' @param facet_by A character vector of column names to facet the plots. Default is NULL.
#' @param facet_ncol The number of columns in the facet grid. Default is NULL.
#' @param facet_nrow The number of rows in the facet grid. Default is NULL.
#' @param aspect.ratio The aspect ratio of the plot. Default is NULL.
#' @param ... Other arguments passed to the specific plot function.
#'  * For "bar", [plotthis::BarPlot()].
#'  * For "line", [plotthis::LinePlot()].
#'  * For "heatmap", [plotthis::Heatmap()].
#'  * For "box", [plotthis::BoxPlot()].
#'  * For "violin", [plotthis::ViolinPlot()].
#' @return A ggplot object or a list if `combine` is FALSE
#' @export
#' @importFrom ggplot2 element_blank
#' @importFrom scRepertoire percentAA positionalEntropy positionalProperty
#' @importFrom plotthis BarPlot LinePlot Heatmap BoxPlot ViolinPlot
#' @examples
#' \donttest{
#' set.seed(8525)
#' data(contig_list, package = "scRepertoire")
#' data <- scRepertoire::combineTCR(contig_list,
#'    samples = c("P17B", "P17L", "P18B", "P18L", "P19B","P19L", "P20B", "P20L"))
#' data <- scRepertoire::addVariable(data,
#'   variable.name = "Type",
#'   variables = rep(c("B", "L"), 4)
#' )
#'
#' ClonalPositionalPlot(data)
#' ClonalPositionalPlot(data, method = "shannon")
#' ClonalPositionalPlot(data, method = "norm.entropy", plot_type = "heatmap")
#' ClonalPositionalPlot(data, method = "Atchley", group_by = "Type", plot_type = "bar")
#' ClonalPositionalPlot(data, method = "Atchley", plot_type = "line")
#' }
ClonalPositionalPlot <- function (
    data, chain = "TRB", aa_length = 20, group_by = "Sample", group_by_sep = "_", split_by = NULL,
    method = c("AA", "shannon", "inv.simpson", "norm.entropy", "Atchley",
        "Kidera", "stScales", "tScales", "VHSE"),
    plot_type = c("bar", "line", "heatmap", "box", "violin"), theme_args = list(),
    xlab = NULL, ylab = NULL, facet_by = NULL, facet_ncol = NULL, facet_nrow = NULL,
    aspect.ratio = NULL,
    ...
) {
    method <- match.arg(method)
    plot_type <- match.arg(plot_type)
    if (plot_type %in% c("box", "violin")) {
        if (is.null(group_by) || identical(group_by, "Sample")) {
            stop("'group_by' must be provided for box/violin ClonalPositionalPlot")
        }
        all_groupings <- unique(c("Sample", group_by, facet_by, split_by))
    } else {
        all_groupings <- unique(c(group_by, facet_by, split_by))
    }
    grouping_levels <- sapply(all_groupings, function(g) {
        dg <- if (inherits(data, "Seurat")) {
            data@meta.data[[g]]
        } else {
            data[[g]]
        }
        if (is.null(dg)) return(NULL)
        if (!is.factor(dg)) dg <- factor(dg)
        levels(dg)
    })
    grouping_levels <- grouping_levels[!sapply(grouping_levels, is.null)]

    data <- merge_clonal_groupings(data, all_groupings)

    if (method == "AA") {
        if (!is.null(facet_by)) {
            stop("'facet_by' should not be specified for AA bar plot in ClonalPositionalPlot.")
        }
        data <- percentAA(data, chain = chain, aa.length = aa_length, group.by = ".group",
            exportTable = TRUE)
        data <- separate(data, "group", into = all_groupings, sep = " // ")
        for (gl in names(grouping_levels)) {
            if (!is.null(data[[gl]])) {
                data[[gl]] <- factor(data[[gl]], levels = grouping_levels[[gl]])
            }
        }

        if (plot_type == "bar") {
            theme_args$panel.grid.major.y <- theme_args$panel.grid.major.y %||% element_blank()

            BarPlot(data, x = "variable", y = "value", group_by = "AA", position = "stack",
                xlab = xlab %||% "Position", ylab = ylab %||% "Amino Acid Frequency",
                split_by = split_by, facet_by = group_by, facet_ncol = facet_ncol %||% 1, facet_nrow = facet_nrow,
                x_text_angle = 90, facet_args = list(strip.position = "right"),
                aspect.ratio = aspect.ratio %||% (2 / aa_length), theme_args = theme_args, ...
            )
        } else if (plot_type == "heatmap") {
            data <- data %>% unite(".group", !!!syms(group_by), sep = group_by_sep)
            allgroups <- unique(data$.group)
            data <- data %>%
                pivot_wider(names_from = ".group", values_from = "value") %>%
                rename(Position = "variable")

            Heatmap(data, columns_by = "Position", rows_by = allgroups, rows_name = paste(group_by, collapse = group_by_sep),
                cell_type = "pie", pie_group_by = "AA", cluster_rows = FALSE, cluster_columns = FALSE,
                pie_values = "sum", ...)
        } else {
            stop("Only 'bar' and 'heatmap' plot types are supported for AA in ClonalPositionalPlot.")
        }
    } else if (method %in% c("shannon", "inv.simpson", "norm.entropy")) {
        data <- positionalEntropy(data, chain = chain, aa.length = aa_length, group.by = ".group",
            method = method, exportTable = TRUE) %>%
            separate("Var1", into = all_groupings, sep = " // ") %>%
            rename(Position = "Var2") %>%
            unite(".group", !!!syms(group_by), sep = group_by_sep)

        for (gl in names(grouping_levels)) {
            if (!is.null(data[[gl]])) {
                data[[gl]] <- factor(data[[gl]], levels = grouping_levels[[gl]])
            }
        }

        group_by <- paste(group_by, sep = group_by_sep)
        data <- rename(data, !!sym(group_by) := ".group")

        if (plot_type == "bar") {
            if (!is.null(facet_by)) {
                stop("'facet_by' should not be specified for entropy bar plot in ClonalPositionalPlot.")
            }
            theme_args$panel.grid.major.y <- theme_args$panel.grid.major.y %||% element_blank()

            BarPlot(data, x = "Position", y = "value",
                xlab = xlab %||% "Position", ylab = ylab %||% method, split_by = split_by,
                facet_by = group_by, facet_ncol = facet_ncol %||% 1, facet_nrow = facet_nrow,
                x_text_angle = 90, legend.position = "none", facet_args = list(strip.position = "right"),
                aspect.ratio = aspect.ratio %||% (2 / aa_length), theme_args = theme_args, ...
            )
        } else if (plot_type == "line") {
            LinePlot(data, x = "Position", y = "value", group_by = group_by, pt_size = 2,
                xlab = xlab %||% "Position", ylab = ylab %||% method, split_by = split_by,
                facet_by = facet_by, facet_ncol = facet_ncol, facet_nrow = facet_nrow, x_text_angle = 90,
                facet_args = list(strip.position = "right"), aspect.ratio = aspect.ratio %||% (6 / aa_length),
                theme_args = theme_args, ...
            )
        } else if (plot_type == "heatmap") {
            allgroups <- unique(data[[group_by]])
            data <- data %>% pivot_wider(names_from = group_by, values_from = "value")

            Heatmap(data, columns_by = "Position", rows_by = allgroups, rows_name = group_by,
                name = method, cluster_columns = FALSE, show_column_names = TRUE, show_row_names = TRUE,
                ...)
        } else if (plot_type == "box") {
            BoxPlot(data, x = "Position", y = "value", xlab = xlab %||% "Position",
                ylab = ylab %||% method, split_by = split_by, group_by = group_by, facet_ncol = facet_ncol,
                facet_nrow = facet_nrow, x_text_angle = 90, theme_args = theme_args,
                aspect.ratio = aspect.ratio %||% (10 / aa_length), ...
            )
        } else if (plot_type == "violin") {
            ViolinPlot(data, x = "Position", y = "value", xlab = xlab %||% "Position",
                ylab = ylab %||% method, split_by = split_by, group_by = group_by, facet_ncol = facet_ncol,
                facet_nrow = facet_nrow, x_text_angle = 90, theme_args = theme_args,
                aspect.ratio = aspect.ratio %||% (10 / aa_length), ...
            )
        }
    } else {
        # https://github.com/ncborcherding/scRepertoire/issues/420
        data <- positionalProperty(data, chain = chain, aa.length = aa_length, group.by = ".group",
            method = method)$data %>%
            separate("group", into = all_groupings, sep = " // ") %>%
            rename(Position = "position") %>%
            unite(".group", !!!syms(group_by), sep = group_by_sep)

        for (gl in names(grouping_levels)) {
            if (!is.null(data[[gl]])) {
                data[[gl]] <- factor(data[[gl]], levels = grouping_levels[[gl]])
            }
        }

        group_by <- paste(group_by, sep = group_by_sep)
        data <- rename(data, !!sym(group_by) := ".group")

        n_properties <- length(unique(data$property))

        if (plot_type == "bar") {
            if (!is.null(facet_by)) {
                stop("'facet_by' should not be specified for property bar plot in ClonalPositionalPlot.")
            }
            theme_args$panel.grid.major.y <- theme_args$panel.grid.major.y %||% element_blank()

            BarPlot(data, x = "Position", y = "mean",
                xlab = xlab %||% "Position", ylab = ylab %||% "Mean Values", split_by = split_by,
                facet_by = c("property", group_by), facet_ncol = facet_ncol, facet_nrow = facet_nrow %||% n_properties,
                x_text_angle = 90, legend.position = "none",
                aspect.ratio = aspect.ratio %||% (4 / aa_length), theme_args = theme_args, ...
            )
        } else if (plot_type == "line") {
            theme_args$panel.grid.major.y <- theme_args$panel.grid.major.y %||% element_blank()
            LinePlot(data, x = "Position", y = "mean", group_by = group_by, pt_size = 2,
                xlab = xlab %||% "Position", ylab = ylab %||% "Mean Values", split_by = split_by,
                facet_by = "property", facet_ncol = facet_ncol %||% 1, facet_nrow = facet_nrow, x_text_angle = 90,
                facet_args = list(strip.position = "right"), aspect.ratio = aspect.ratio %||% (6 / aa_length),
                theme_args = theme_args, ...
            )
        } else {
            stop("Only 'bar' and 'line' plot types are supported for property in ClonalPositionalPlot.")
        }
    }
}

#' ClonalKmerPlot
#'
#' Explore the k-mer frequency of CDR3 sequences.
#' @param data The product of [scRepertoire::combineTCR], [scRepertoire::combineTCR], or
#'  [scRepertoire::combineExpression].
#' @param chain The chain to be analyzed. Default is "TRB".
#' @param clone_call The column name of the clone call. Default is "aa".
#' @param k The length of the k-mer. Default is 3.
#' @param top The number of top k-mers to display. Default is 25.
#' @param group_by The variable to group the data by. Default is "Sample".
#' @param group_by_sep The separator to use when combining groupings. Default is "_".
#' @param facet_by A character vector of column names to facet the plots. Default is NULL.
#' @param split_by A character vector of column names to split the plots. Default is NULL.
#' @param plot_type The type of plot to generate. Default is "bar".
#'  * "bar": Bar plot.
#'  * "line": Line plot.
#'  * "heatmap": Heatmap.
#' @param theme_args A list of arguments to be passed to the [ggplot2::theme] function.
#' @param aspect.ratio The aspect ratio of the plot. Default is NULL.
#' @param facet_ncol The number of columns in the facet grid. Default is NULL.
#' @param ... Other arguments passed to the specific plot function.
#'  * For "bar", [plotthis::BarPlot()].
#'  * For "line", [plotthis::LinePlot()].
#'  * For "heatmap", [plotthis::Heatmap()].
#' @return A ggplot object or a list if `combine` is FALSE
#' @export
#' @importFrom tidyr pivot_longer separate unite
#' @importFrom dplyr %>% rename
#' @importFrom scRepertoire percentKmer
#' @importFrom plotthis BarPlot Heatmap
#' @examples
#' \donttest{
#' set.seed(8525)
#' data(contig_list, package = "scRepertoire")
#' data <- scRepertoire::combineTCR(contig_list,
#'     samples = c("P17B", "P17L", "P18B", "P18L", "P19B","P19L", "P20B", "P20L"))
#' data <- scRepertoire::addVariable(data,
#'     variable.name = "Type",
#'     variables = rep(c("B", "L"), 4)
#' )
#' data <- scRepertoire::addVariable(data,
#'     variable.name = "Subject",
#'     variables = rep(c("P17", "P18", "P19", "P20"), each = 2)
#' )
#'
#' ClonalKmerPlot(data)
#' ClonalKmerPlot(data, group_by = "Type")
#' ClonalKmerPlot(data, group_by = "Type", plot_type = "line")
#' ClonalKmerPlot(data, group_by = "Type", plot_type = "heatmap")
#' }
ClonalKmerPlot <- function (
    data, chain = "TRB", clone_call = "aa", k = 3, top = 25, group_by = "Sample",
    group_by_sep = "_", facet_by = NULL, split_by = NULL,
    plot_type = c("bar", "line", "heatmap"), theme_args = list(), aspect.ratio = NULL,
    facet_ncol = NULL, ...
) {
    plot_type <- match.arg(plot_type)
    all_groupings <- unique(c(group_by, split_by))
    grouping_levels <- sapply(all_groupings, function(g) {
        dg <- if (inherits(data, "Seurat")) {
            data@meta.data[[g]]
        } else {
            data[[g]]
        }
        if (is.null(dg)) return(NULL)
        if (!is.factor(dg)) dg <- factor(dg)
        levels(dg)
    })
    grouping_levels <- grouping_levels[!sapply(grouping_levels, is.null)]

    data <- merge_clonal_groupings(data, all_groupings)
    data <- percentKmer(data, chain = chain, cloneCall = clone_call, motif.length = k,
        top.motifs = top, group.by = ".group", exportTable = TRUE)
    data <- as.data.frame(data)
    motifs <- colnames(data)
    data$.group <- rownames(data)
    data <- data %>%
        separate(".group", into = all_groupings, sep = " // ") %>%
        unite(".group", !!!syms(group_by), sep = group_by_sep) %>%
        rename(!!sym(paste(group_by, sep = group_by_sep)) := ".group")

    for (gl in names(grouping_levels)) {
        if (!is.null(data[[gl]])) {
            data[[gl]] <- factor(data[[gl]], levels = grouping_levels[[gl]])
        }
    }

    group_by <- paste(group_by, sep = group_by_sep)

    if (plot_type == "bar") {
        if (!is.null(facet_by)) {
            stop("'facet_by' should not be specified in bar ClonalKmerPlot.")
        }

        data <- data %>% pivot_longer(cols = motifs, names_to = "Motifs", values_to = "Frequency")
        theme_args$panel.grid.major.y <- theme_args$panel.grid.major.y %||% element_blank()

        BarPlot(data, x = "Motifs", y = "Frequency", facet_by = group_by,
            xlab = "Motifs", ylab = "Frequency", split_by = split_by,
            facet_ncol = facet_ncol %||% 1, x_text_angle = 90, facet_args = list(strip.position = "right"),
            aspect.ratio = aspect.ratio %||% (4 / length(motifs)), legend.position = "none",
            theme_args = theme_args, ...
        )
    } else if (plot_type == "line") {
        data <- data %>% pivot_longer(cols = motifs, names_to = "Motifs", values_to = "Frequency")
        theme_args$panel.grid.major.y <- theme_args$panel.grid.major.y %||% element_blank()

        LinePlot(data, x = "Motifs", y = "Frequency", group_by = group_by, pt_size = 2,
            xlab = "Motifs", ylab = "Frequency", split_by = split_by, facet_by = facet_by,
            facet_ncol = facet_ncol, x_text_angle = 90,
            aspect.ratio = aspect.ratio %||% (8 / length(motifs)), theme_args = theme_args, ...
        )
    } else if (plot_type == "heatmap") {
        args <- rlang::dots_list(...)
        args$data <- data
        args$columns_by <- group_by
        args$rows_by <- motifs
        args$rows_name <- "Motifs"
        args$name <- "Frequency"
        args$show_row_names <- args$show_row_names %||% TRUE
        args$show_column_names <- args$show_column_names %||% TRUE
        do.call(Heatmap, args)
    }
}
