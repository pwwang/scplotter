#' ClonalVolumePlot
#'
#' @param data The product of [scRepertoire::combineTCR], [scRepertoire::combineTCR], or
#'  [scRepertoire::combineExpression].
#' @param clone_call How to call the clone - VDJC gene (gene), CDR3 nucleotide (nt),
#'  CDR3 amino acid (aa), VDJC gene + CDR3 nucleotide (strict) or a custom variable
#'  in the data
#' @param chain indicate if both or a specific chain should be used - e.g. "both",
#'  "TRA", "TRG", "IGH", "IGL"
#' @param scale Whether to use clone proportion or clone size for the plot.
#' @param plot_type The type of plot to use. Default is "bar".
#'  Possible values are "bar", "box", and "violin".
#'  When "box" or "violin" is used, the data will be broken down by the Sample and plotted
#'  for each group.
#' @param x The column name in the meta data to use as the x-axis. Default: "Sample"
#' @param ylab The y-axis label.
#' @param group_by The column name in the meta data to group the cells. Default: NULL
#' @param facet_by The column name in the meta data to facet the plots. Default: NULL
#' @param split_by The column name in the meta data to split the plots. Default: NULL
#' @param order The order of the x-axis items or groups. Default is an empty list.
#'  It should be a list of values. The names are the column names, and the values are the order.
#' @param ... Other arguments passed to the specific plot function.
#'  * For `bar` plot, see [plotthis::BarPlot()].
#'  * For `box` plot, see [plotthis::BoxPlot()].
#'  * For `violin` plot, see [plotthis::ViolinPlot()].
#' @return A ggplot object or a list if `combine` is FALSE
#' @importFrom tidyr separate
#' @importFrom plotthis BarPlot BoxPlot ViolinPlot
#' @importFrom scRepertoire clonalQuant
#' @export
#' @examples
#' \donttest{
#' set.seed(8525)
#' data(contig_list, package = "scRepertoire")
#' data <- scRepertoire::combineTCR(contig_list)
#' data <- scRepertoire::addVariable(data,
#'     variable.name = "Type",
#'     variables = sample(c("B", "L"), 8, replace = TRUE)
#' )
#' data <- scRepertoire::addVariable(data,
#'     variable.name = "Sex",
#'     variables = sample(c("M", "F"), 8, replace = TRUE)
#' )
#'
#' ClonalVolumePlot(data)
#' ClonalVolumePlot(data, x = "Type")
#' ClonalVolumePlot(data, x = "Type", order = list(Type = c("L", "B")))
#' ClonalVolumePlot(data, x = c("Type", "Sex"), scale = TRUE)
#' ClonalVolumePlot(data, x = "Type", group_by = "Sex", position = "stack")
#' ClonalVolumePlot(data,
#'     plot_type = "box", x = "Type", comparisons = TRUE,
#'     group_by = "Sex"
#' )
#' ClonalVolumePlot(data, plot_type = "violin", x = "Type", add_box = TRUE)
#'
#' # on a Seurat object
#' data(scRep_example, package = "scRepertoire")
#' data(contig_list, package = "scRepertoire")
#' combined <- scRepertoire::combineTCR(contig_list,
#'     samples = c("P17B", "P17L", "P18B", "P18L", "P19B", "P19L", "P20B", "P20L")
#' )
#' sobj <- scRepertoire::combineExpression(combined, scRep_example)
#' ClonalVolumePlot(sobj)
#' ClonalVolumePlot(sobj, x = "seurat_clusters")
#' ClonalVolumePlot(sobj, group_by = "seurat_clusters")
#' ClonalVolumePlot(sobj, x = "seurat_clusters", plot_type = "box")
#' }
ClonalVolumePlot <- function(
    data, clone_call = "aa", chain = "both", scale = FALSE,
    plot_type = c("bar", "box", "violin"), x = "Sample", group_by = NULL,
    facet_by = NULL, split_by = NULL, order = list(), ylab = NULL, ...) {
    plot_type <- match.arg(plot_type)
    if (plot_type %in% c("box", "violin")) {
        all_groupings <- unique(c("Sample", x, group_by, facet_by, split_by))
    } else {
        all_groupings <- unique(c(x, group_by, facet_by, split_by))
    }
    data <- merge_clonal_groupings(data, all_groupings)
    data <- clonalQuant(data,
        cloneCall = clone_call, chain = chain, scale = scale,
        group.by = ".group", exportTable = TRUE
    )
    # restore the groups
    data <- separate(data, ".group", into = all_groupings, sep = " // ")
    for (group in all_groupings) {
        if (!is.null(order[[group]])) {
            data[[group]] <- factor(data[[group]], levels = order[[group]])
        }
    }
    if (scale) {
        data$scaled <- data$scaled / 100
    }

    if (plot_type == "bar") {
        BarPlot(data,
            x = x, y = ifelse(scale, "scaled", "contigs"),
            group_by = group_by, facet_by = facet_by, split_by = split_by,
            ylab = ylab %||% ifelse(scale, "Fraction of Unique Clones", "Number of Unique Clones"),
            ...
        )
    } else if (plot_type == "box") {
        BoxPlot(data,
            x = x, y = ifelse(scale, "scaled", "contigs"), group_by = group_by,
            facet_by = facet_by, split_by = split_by,
            ylab = ylab %||% ifelse(scale, "Fraction of Unique Clones", "Number of Unique Clones"),
            ...
        )
    } else {
        ViolinPlot(data,
            x = x, y = ifelse(scale, "scaled", "contigs"),
            group_by = group_by, facet_by = facet_by, split_by = split_by,
            ylab = ylab %||% ifelse(scale, "Fraction of Unique Clones", "Number of Unique Clones"),
            ...
        )
    }
}


#' ClonalAbundancePlot
#'
#' @description Plot the count or density of the clones at different abundance levels.
#'
#' @param data The product of [scRepertoire::combineTCR], [scRepertoire::combineTCR], or
#'  [scRepertoire::combineExpression].
#' @param clone_call How to call the clone - VDJC gene (gene), CDR3 nucleotide (nt),
#'  CDR3 amino acid (aa), VDJC gene + CDR3 nucleotide (strict) or a custom variable
#'  in the data
#' @param chain indicate if both or a specific chain should be used - e.g. "both",
#'  "TRA", "TRG", "IGH", "IGL"
#' @param plot_type The type of plot to use. Default is "trend".
#'  Possible values are "trend", "histogram" and "density".
#' @param trend_skip_zero Whether to skip the zero values in the trend line. Default is TRUE.
#' @param binwidth The binwidth for the histogram plot. Default is 0.1.
#' @param group_by The column name in the meta data to group the cells. Default: "Sample"
#' @param group_by_sep The separator to use when combining the group_by columns. Default: "_"
#' @param facet_by The column name in the meta data to facet the plots. Default: NULL
#' @param split_by The column name in the meta data to split the plots. Default: NULL
#' @param order The order of the x-axis items or groups. Default is an empty list.
#'  It should be a list of values. The names are the column names, and the values are the order.
#' @param bw The smoothing bandwidth to be used for density plots. Default is 0.5.
#' @param xlab The x-axis label. Default is "Abundance".
#' @param ylab The y-axis label. Default is "Number of Clones" for trend and histogram, and
#'  "Density of Clones" for density.
#' @param xtrans The transformation to apply to the x-axis. Default is "log10".
#' @param ytrans The transformation to apply to the y-axis. Default is "identity".
#' @param theme_args The theme arguments to be passed to the plot function.
#' @param ... Other arguments passed to the specific plot function.
#'  * For `trend` plot, see [plotthis::Histogram()].
#'  * For `histogram` plot, see [plotthis::Histogram()].
#'  * For `density` plot, see [plotthis::DensityPlot()].
#' @return A ggplot object or a list if `combine` is FALSE
#' @importFrom tidyr full_seq complete
#' @importFrom dplyr group_modify
#' @importFrom plotthis Histogram DensityPlot
#' @importFrom scRepertoire clonalAbundance
#' @export
#' @examples
#' \donttest{
#' set.seed(8525)
#' data(contig_list, package = "scRepertoire")
#' data <- scRepertoire::combineTCR(contig_list)
#' data <- scRepertoire::addVariable(data,
#'     variable.name = "Type",
#'     variables = sample(c("B", "L"), 8, replace = TRUE)
#' )
#' data <- scRepertoire::addVariable(data,
#'     variable.name = "Sex",
#'     variables = sample(c("M", "F"), 8, replace = TRUE)
#' )
#'
#' ClonalAbundancePlot(data)
#' ClonalAbundancePlot(data, ytrans = "log10")
#' ClonalAbundancePlot(data, plot_type = "histogram")
#' ClonalAbundancePlot(data, plot_type = "histogram", add_trend = TRUE, trend_skip_zero = TRUE)
#' ClonalAbundancePlot(data, plot_type = "density")
#' }
ClonalAbundancePlot <- function(
    data, clone_call = "aa", chain = "both", xtrans = "log10", ytrans = "identity",
    plot_type = c("trend", "histogram", "density"), binwidth = 0.1, trend_skip_zero = TRUE,
    bw = 0.5, group_by = "Sample", group_by_sep = "_", facet_by = NULL, split_by = NULL,
    order = list(), xlab = "Abundance", ylab = NULL, theme_args = list(), ...) {
    plot_type <- match.arg(plot_type)

    all_groupings <- unique(c(group_by, facet_by, split_by))
    data <- merge_clonal_groupings(data, all_groupings)
    if (length(all_groupings) > 0) {
        data <- clonalAbundance(data,
            cloneCall = clone_call, chain = chain,
            group.by = ".group", exportTable = TRUE
        )

        # restore the groups
        data <- separate(data, ".group", into = all_groupings, sep = " // ")
        for (group in all_groupings) {
            if (!is.null(order[[group]])) {
                data[[group]] <- factor(data[[group]], levels = order[[group]])
            }
        }
    } else {
        data <- clonalAbundance(data,
            cloneCall = clone_call, chain = chain,
            group.by = group_by, exportTable = TRUE
        )
    }

    if (plot_type == "trend") {
        theme_args$panel.grid.major.y <- theme_args$panel.grid.major.y %||%
            element_line(color = "grey", linetype = 2)
        Histogram(data,
            x = "Abundance", group_by = group_by, group_by_sep = group_by_sep,
            facet_by = facet_by, split_by = split_by, xtrans = xtrans, xlab = xlab, binwidth = binwidth,
            ytrans = ytrans, ylab = ylab %||% "Number of Clones", use_trend = TRUE,
            trend_skip_zero = trend_skip_zero, theme_args = theme_args, ...
        )
    } else if (plot_type == "histogram") {
        Histogram(data,
            x = "Abundance", group_by = group_by, group_by_sep = group_by_sep,
            facet_by = facet_by, split_by = split_by, xtrans = xtrans, xlab = xlab, binwidth = binwidth,
            trend_skip_zero = trend_skip_zero, ytrans = ytrans, theme_args = theme_args,
            ylab = ylab %||% "Number of Clones", ...
        )
    } else if (plot_type == "density") {
        DensityPlot(data,
            x = "Abundance", group_by = group_by, group_by_sep = group_by_sep,
            facet_by = facet_by, split_by = split_by, xtrans = xtrans, xlab = xlab,
            ytrans = ytrans, ylab = ylab %||% "Density of Clones", theme_args = theme_args,
            bw = bw, ...
        )
    }
}

#' ClonalLengthPlot
#'
#' @description Plot the length distribution of the CDR3 sequences
#' @param data The product of [scRepertoire::combineTCR], [scRepertoire::combineTCR], or
#'  [scRepertoire::combineExpression].
#' @param clone_call How to call the clone - only "nt" or "aa" is supported.
#' @param chain indicate if both or a specific chain should be used - e.g. "both",
#'  "TRA", "TRB", "TRD", "TRG", "IGH", or "IGL" to specify a specific chain.
#' @param plot_type The type of plot to use. Default is "bar".
#'  Possible values are "box", "violin" and "density".
#' @param x_nbreaks The number of breaks for the x-axis. Default is 10.
#' @param group_by The column name in the meta data to group the cells. Default: "Sample"
#' @param facet_by The column name in the meta data to facet the plots. Default: NULL
#' @param split_by The column name in the meta data to split the plots. Default: NULL
#' @param order The order of the groups. Default is an empty list.
#'  It should be a list of values. The names are the column names, and the values are the order.
#' @param xlab The x-axis label.
#' @param ylab The y-axis label.
#' @param position The position of the bars for bar plot on the x-axis. Default is "dodge".
#' @param ... Other arguments passed to the specific plot function.
#' * For `bar` plot, see [plotthis::BarPlot()].
#' * For `box` plot, see [plotthis::BoxPlot()].
#' * For `violin` plot, see [plotthis::ViolinPlot()].
#' * For `density` plot, see [plotthis::DensityPlot()].
#' @return A ggplot object or a list if `combine` is FALSE
#' @importFrom stats quantile
#' @importFrom rlang syms
#' @importFrom dplyr summarise n
#' @importFrom tidyr separate
#' @importFrom ggplot2 element_blank scale_x_discrete element_line
#' @importFrom plotthis BarPlot
#' @importFrom scRepertoire clonalLength
#' @export
#' @examples
#' \donttest{
#' set.seed(8525)
#' data(contig_list, package = "scRepertoire")
#' data <- scRepertoire::combineTCR(contig_list)
#' data <- scRepertoire::addVariable(data, variable.name = "Type",
#'  variables = sample(c("B", "L"), 8, replace = TRUE))
#' data <- scRepertoire::addVariable(data, variable.name = "Sex",
#'  variables = sample(c("M", "F"), 8, replace = TRUE))
#'
#' ClonalLengthPlot(data)
#' ClonalLengthPlot(data, plot_type = "box")
#' ClonalLengthPlot(data, clone_call = "nt", plot_type = "violin", chain = "TRB",
#'  group_by = "Type", comparisons = TRUE)
#' ClonalLengthPlot(data, plot_type = "density", chain = "TRA")
#' }
ClonalLengthPlot <- function(
    data, clone_call = "aa", chain = "both", plot_type = c("bar", "box", "violin", "density"),
    x_nbreaks = 10, group_by = "Sample", order = list(), xlab = "Length", ylab = NULL,
    position = "dodge", facet_by = NULL, split_by = NULL, ...) {
    plot_type <- match.arg(plot_type)

    if (plot_type %in% c("box", "violin")) {
        all_groupings <- unique(c("Sample", group_by, facet_by, split_by))
    } else {
        all_groupings <- unique(c(group_by, facet_by, split_by))
    }
    data <- merge_clonal_groupings(data, all_groupings)

    if (identical(all_groupings, "Sample")) {
        data <- clonalLength(data, cloneCall = clone_call, chain = chain, exportTable = TRUE)
        data$Sample <- data$values
    } else {
        data <- clonalLength(data,
            cloneCall = clone_call, chain = chain, group.by = ".group",
            exportTable = TRUE
        )
        data <- separate(data, ".group", into = all_groupings, sep = " // ")
    }

    # restore the groups
    for (group in all_groupings) {
        if (!is.null(order[[group]])) {
            data[[group]] <- factor(data[[group]], levels = order[[group]])
        }
    }

    if (plot_type == "density") {
        rawdata <- data
    }
    data <- data %>%
        dplyr::group_by(!!!syms(c("length", all_groupings))) %>%
        summarise(.n = n(), .groups = "drop")

    default_ylab <- ifelse(clone_call == "aa", "Number of CDR3 (AA)", "Number of CDR3 (NT)")
    if (plot_type == "bar") {
        args <- list(...)
        args$data <- data
        args$x <- "length"
        args$y <- ".n"
        args$group_by <- group_by
        args$facet_by <- facet_by
        args$split_by <- split_by
        args$position <- position
        args$xlab <- "Length"
        args$ylab <- ylab %||% default_ylab
        args$theme_args <- args$theme_args %||% list()
        args$theme_args$panel.grid.major.x <- element_blank()
        breaks <- quantile(data$length, probs = seq(0, 1, length.out = x_nbreaks), type = 3)
        if (is.null(args$split_by)) {
            suppressMessages({
                do.call(BarPlot, args) + scale_x_discrete(breaks = breaks)
            })
        } else {
            suppressMessages({
                do.call(BarPlot, args) & scale_x_discrete(breaks = breaks)
            })
        }
    } else if (plot_type == "box") {
        args <- list(...)
        args$data <- data
        args$x <- "length"
        args$y <- ".n"
        args$group_by <- if (identical(group_by, "Sample")) NULL else group_by
        args$facet_by <- facet_by
        args$split_by <- split_by
        args$xlab <- "Length"
        args$ylab <- ylab %||% default_ylab
        args$theme_args <- args$theme_args %||% list()
        if (clone_call == "nt" || chain == "both") {
            args$theme_args$panel.grid.major.x <- args$theme_args$panel.grid.major.x %||% element_blank()
        }
        args$theme_args$panel.grid.major.y <- args$theme_args$panel.grid.major.y %||%
            element_line(color = "grey", linetype = 2)
        args$legend.position <- args$legend.position %||% ifelse(is.null(args$split_by), "none", "right")
        do.call(BoxPlot, args)
    } else if (plot_type == "violin") {
        args <- list(...)
        args$data <- data
        args$x <- "length"
        args$y <- ".n"
        args$group_by <- if (identical(group_by, "Sample")) NULL else group_by
        args$facet_by <- facet_by
        args$split_by <- split_by
        args$xlab <- "Length"
        args$ylab <- ylab %||% default_ylab
        args$theme_args <- args$theme_args %||% list()
        if (clone_call == "nt" || chain == "both") {
            args$theme_args$panel.grid.major.x <- args$theme_args$panel.grid.major.x %||% element_blank()
        }
        args$theme_args$panel.grid.major.y <- args$theme_args$panel.grid.major.y %||%
            element_line(color = "grey", linetype = 2)
        args$legend.position <- args$legend.position %||% ifelse(is.null(args$split_by), "none", "right")
        do.call(ViolinPlot, args)
    } else if (plot_type == "density") {
        DensityPlot(
            rawdata,
            x = "length",
            group_by = group_by,
            facet_by = facet_by,
            split_by = split_by,
            xlab = xlab %||% "Length",
            ylab = ylab %||% ifelse(clone_call == "aa", "Density of CDR3 (AA)", "Density of CDR3 (NT)"),
            ...
        )
    }
}

#' DummyClonalScatterPlot
#'
#' @description Function to plot the scatter plot of the clonal data for a dummy group pair.
#' @param df The data frame with the clonal data.
#'  The data frame should have the columns: group_by, 'count', 'fraction'.
#' @return A ggplot object
#' @importFrom rlang sym
#' @importFrom stats cor.test
#' @importFrom circlize colorRamp2
#' @importFrom dplyr case_when mutate distinct summarise group_by first
#' @importFrom ggplot2 geom_segment aes element_blank theme scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 guides geom_point scale_fill_manual guide_legend
#' @importFrom ggnewscale new_scale_fill
#' @importFrom plotthis ScatterPlot palette_this
#' @keywords internal
DummyClonalScatterPlot <- function(df, title, group_by, scatter_cor, size_by, ...) {
    if (nrow(df) == 0) {
        stop("No data found for the group_by: ", group_by,
            ". Did you specify the correct 'group_by'/'groups' parameters?")
    }
    pair <- unique(as.character(df[[group_by]]))
    if (length(pair) != 2) {
        stop("The group_by should have exactly 2 unique values.")
    }

    exponent <- function(x) { floor(log10(abs(x))) }

    mantissa <- function(x) {
        mant <- log10(abs(x))
        10 ^ (mant - floor(mant))
    }

    df <- pivot_wider(df, names_from = group_by, values_from = c("count", "fraction"), values_fill = 0)
    suf1 <- paste0("count_", pair[1])
    suf2 <- paste0("count_", pair[2])
    df <- df[df[[suf1]] > 0 | df[[suf2]] > 0, , drop = FALSE]
    N <- nrow(df)
    dual <- which(df[[suf1]] > 0 & df[[suf2]] > 0)
    if (length(dual) <= 2) {
        test <- list(estimate = NA, p.value = NA)
    } else {
        test <- cor.test(log(df[[suf1]][dual]), log(df[[suf2]][dual]), method = scatter_cor)
    }
    sum_counts1 <- sum(df[[suf1]])
    sum_counts2 <- sum(df[[suf2]])

    counts1_norm <- jitter(1 + df[[suf1]], amount = 0.25) / sum_counts1
    counts2_norm <- jitter(1 + df[[suf2]], amount = 0.25) / sum_counts2

    # Avoid some points always overlaying each other
    oo <- sample(length(counts1_norm))
    plotdata <- data.frame(x = counts1_norm[oo], y = counts2_norm[oo])
    names(plotdata) <- pair

    n_singlet1 <- sum(df[[suf1]] == 1 & df[[suf2]] == 0)
    n_singlet2 <- sum(df[[suf1]] == 0 & df[[suf2]] == 1)
    n_expanded1 <- sum(df[[suf1]] > 1 & df[[suf2]] == 0)
    n_expanded2 <- sum(df[[suf1]] == 0 & df[[suf2]] > 1)
    n_dual1 <- sum(df[[suf1]] > 0 & df[[suf2]] > 0 & df[[suf1]] > df[[suf2]])
    n_dual2 <- sum(df[[suf1]] > 0 & df[[suf2]] > 0 & df[[suf1]] < df[[suf2]])
    n_dual <- sum(df[[suf1]] > 0 & df[[suf2]] > 0 & df[[suf1]] == df[[suf2]])

    labels <- c(
        paste0(pair[1], " Singlet (", n_singlet1, ")"),
        paste0(pair[1], " Expanded (", n_expanded1, ")"),
        paste0(pair[2], " Singlet (", n_singlet2, ")"),
        paste0(pair[2], " Expanded (", n_expanded2, ")"),
        paste0("Dual (", pair[2], " < ", pair[1], ") (", n_dual1, ")"),
        paste0("Dual (", pair[2], " > ", pair[1], ") (", n_dual2, ")"),
        paste0("Dual (Equal) (", n_dual, ")")
    )

    plotdata <- mutate(
        plotdata,
        Type = case_when(
            df[[suf1]][oo] == 1 & df[[suf2]][oo] == 0 ~ 1,
            df[[suf1]][oo] > 1 & df[[suf2]][oo] == 0 ~ 3,
            df[[suf1]][oo] == 0 & df[[suf2]][oo] > 1 ~ 9,
            df[[suf1]][oo] == 0 & df[[suf2]][oo] == 1 ~ 11,
            df[[suf1]][oo] > df[[suf2]][oo] ~ 5,
            df[[suf1]][oo] < df[[suf2]][oo] ~ 7,
            TRUE ~ 6
        ),
        TypeName = case_when(
            Type == 1 ~ labels[1],
            Type == 3 ~ labels[2],
            Type == 9 ~ labels[4],
            Type == 11 ~ labels[3],
            Type == 5 ~ labels[5],
            Type == 7 ~ labels[6],
            TRUE ~ labels[7]
        ),
        Max_Size = pmax(df[[suf1]][oo], df[[suf2]][oo]),
        Total_Size = df[[suf1]][oo] + df[[suf2]][oo],
        # Make sure color similar in each category
        NumType = !!sym("Type") + scales::rescale(!!sym("Total_Size"), to = c(0, 1))
    )

    xbreaks <- c(
        1 / sum_counts1,
        0.001 + 1 / sum_counts1,
        0.01 + 1 / sum_counts1,
        0.1 + 1 / sum_counts1
    )
    ybreaks <- c(
        1 / sum_counts2,
        0.001 + 1 / sum_counts2,
        0.01 + 1 / sum_counts2,
        0.1 + 1 / sum_counts2
    )

    minx <- min(plotdata[[pair[1]]], na.rm = TRUE)
    miny <- min(plotdata[[pair[1]]], na.rm = TRUE)
    maxx <- max(plotdata[[pair[2]]], na.rm = TRUE)
    maxy <- max(plotdata[[pair[2]]], na.rm = TRUE)

    n_formatted <- formatC(length(oo), format = "f", big.mark = ",", digits = 0)
    r_formatted <- format(test$estimate, digits = 2, scientific = F)
    if (is.na(test$p.value)) {
        subtitle <- bquote(
            italic(N) == .(N) ~ ~ italic(N)[dual] == .(length(dual)) ~ ~ italic(r) ==
            .(r_formatted) ~ ~ italic(p) == "NA"
        )
    } else if (test$p.value < 1e-4) {
        P_mant <- format(mantissa(test$p.value), digits = 2)
        P_exp <- exponent(test$p.value)
        subtitle <- bquote(
            italic(N) == .(N) ~ ~ italic(N)[dual] == .(length(dual)) ~ ~ italic(r) ==
            .(r_formatted) ~ ~ italic(p) == .(P_mant) %*% 10^.(P_exp)
        )
    } else {
        P_formatted <- format(test$p.value, digits = 2)
        subtitle <- bquote(
            italic(N) == .(N) ~ ~ italic(N)[dual] == .(length(dual)) ~ ~ italic(r) ==
            .(r_formatted) ~ ~ italic(p) == .(P_formatted)
        )
    }

    colfun <- colorRamp2(
        seq(min(plotdata$NumType, na.rm = TRUE), max(plotdata$NumType, na.rm = TRUE), length.out = 100),
        palette_this(palette = "Spectral", palcolor = NULL)
    )
    label_df <- plotdata %>% group_by(!!sym("TypeName")) %>%
        summarise(
            x = first(!!sym(pair[1])),
            y = first(!!sym(pair[2])),
            Type = mean(!!sym("Type"), na.rm = TRUE), .groups = "drop")
    label_df$TypeName <- factor(label_df$TypeName, levels = labels)
    label_df <- label_df[order(label_df$TypeName), , drop = FALSE]

    p <- ScatterPlot(
            plotdata, x = pair[1], y = pair[2], color_by = "NumType",
            title = if (identical(title, "...")) NULL else title, subtitle = subtitle,
            border_color = TRUE,
            size_by = ifelse(size_by == "max", "Max_Size", "Total_Size"),
            size_name = ifelse(size_by == "max", "Max Size", "Total Size"),
            aspect.ratio = 1, ...) +
        guides(color = "none") +
        new_scale_fill() +
        geom_point(data = label_df, aes(x = !!sym("x"), y = !!sym("y"), fill = !!sym("TypeName")), shape = 21, alpha = 0) +
        scale_fill_manual(
            values = colfun(label_df$Type),
            breaks = label_df$TypeName,
            guide = guide_legend(
                title = "Category",
                order = 9,
                override.aes = list(alpha = 1, color = "grey90", size = 3)
            )) +
        theme(panel.grid.major = element_blank()) +
        geom_segment(
            data = data.frame(
                # diagnal, horizontal, vertical, horizontal short, vertical short
                x = c(1.5 / sum_counts1, minx, 1.5 / sum_counts1, minx, 2.5 / sum_counts1),
                xend = c(maxx, maxx, 1.5 / sum_counts1, 1.5 / sum_counts1, 2.5 / sum_counts1),
                y = c(1.5 / sum_counts2, 1.5 / sum_counts2, miny, 2.5 / sum_counts2, miny),
                yend = c(maxy, 1.5 / sum_counts2, maxy, 2.5 / sum_counts2, 1.5 / sum_counts2)
            ),
            aes(x = !!sym("x"), y = !!sym("y"), xend = !!sym("xend"), yend = !!sym("yend")), color = "gray90"
        )

    # remove the x, y scales
    p$scales$scales[[4]] <- NULL
    p$scales$scales[[4]] <- NULL  # list shrinks after removing one element
    p + scale_x_continuous(trans = "log2", limits = c(minx, maxx), breaks = xbreaks,
            labels = c("0", "0.001", "0.01", "0.1")) +
        scale_y_continuous(trans = "log2", limits = c(miny, maxy), breaks = ybreaks,
            labels = c("0", "0.001", "0.01", "0.1"))
}

#' ClonalResidencyPlot
#'
#' @description Plot the residency of the clones in different samples.
#' @param data The product of [scRepertoire::combineTCR], [scRepertoire::combineTCR], or
#'  [scRepertoire::combineExpression].
#' @param clone_call How to call the clone - VDJC gene (gene), CDR3 nucleotide (nt),
#'  CDR3 amino acid (aa), VDJC gene + CDR3 nucleotide (strict) or a custom variable
#'  in the data
#' @param chain indicate if both or a specific chain should be used - e.g. "both",
#'  "TRA", "TRG", "IGH", "IGL"
#' @param plot_type The type of plot to use. Default is "scatter".
#'  Possible values are "scatter", "venn", and "upset".
#' @param group_by The column name in the meta data to group the cells. Default: "Sample"
#' @param groups The groups to compare. Default is NULL.
#'  If NULL, all the groups in `group_by` will be compared.
#'  Note that for "scatter" plot, only two groups can be compared.
#' @param facet_by The column name in the meta data to facet the plots. Default: NULL
#' @param split_by The column name in the meta data to split the plots. Default: NULL
#' @param split_by_sep The separator used to concatenate the split_by when multiple columns are used.
#' @param scatter_cor The correlation method to use for the scatter plot. Default is "pearson".
#' @param scatter_size_by The size of the points in the scatter plot. Default is "max".
#'  Possible values are "max" and "total".
#'  * "max" - The max size of the clone in the two groups.
#'  * "total" - The total size of the clone in the two groups.
#' @param order The order of the x-axis items or groups. Default is an empty list.
#'  It should be a list of values. The names are the column names, and the values are the order.
#' @param combine Whether to combine the plots into a single plot. Default is TRUE.
#' @param nrow The number of rows in the combined plot. Default is NULL.
#' @param ncol The number of columns in the combined plot. Default is NULL.
#' @param byrow Whether to fill the combined plot by row. Default is TRUE.
#' @param ... Other arguments passed to the specific plot function.
#'  * For `scatter` plot, see [plotthis::ScatterPlot()].
#'  * For `venn` plot, see [plotthis::VennDiagram()].
#'  * For `upset` plot, see [plotthis::UpsetPlot()].
#' @return A ggplot object or a list if `combine` is FALSE
#' @export
#' @importFrom utils combn
#' @importFrom dplyr distinct rename_with select across starts_with ends_with
#' @importFrom tidyr pivot_longer separate pivot_wider unite
#' @importFrom scRepertoire clonalScatter
#' @importFrom plotthis ScatterPlot VennDiagram UpsetPlot
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
#' ClonalResidencyPlot(data, groups = c("P18B", "P18L"))
#' ClonalResidencyPlot(data, group_by = "Type", split_by = "Subject")
#' ClonalResidencyPlot(data, plot_type = "venn", groups = c("B", "L"), group_by = "Type",
#'  split_by = "Subject")
#' ClonalResidencyPlot(data, plot_type = "upset", groups = c("P18B", "P18L"))
#' }
ClonalResidencyPlot <- function(
    data, clone_call = "aa", chain = "both", plot_type = c("scatter", "venn", "upset"),
    group_by = "Sample", groups = NULL, facet_by = NULL, split_by = NULL, split_by_sep = "_",
    scatter_cor = "pearson", scatter_size_by = c("max", "total"),
    order = list(), combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, ...
) {

    stopifnot("ClonalResidencyPlot supports only a single group_by column" = length(group_by) == 1)
    stopifnot("'facet_by' is not supported for 'ClonalResidencyPlot'" = is.null(facet_by))

    plot_type <- match.arg(plot_type)
    scatter_size_by <- match.arg(scatter_size_by)

    all_groupings <- unique(c(group_by, facet_by, split_by))
    data <- clonal_size_data(data, clone_call, chain, all_groupings)

    # restore the groups
    for (group in all_groupings) {
        if (!is.null(order[[group]])) {
            data[[group]] <- factor(data[[group]], levels = order[[group]])
        }
    }

    groups <- groups %||% unique(data[[group_by]])
    if (plot_type == "scatter") {
        if (!is.null(split_by)) {
            data <- unite(data, ".split", split_by, sep = split_by_sep, remove = FALSE)
        } else {
            data$.split <- "..."
        }
        data <- split(data, data$.split)
        spdata <- list()
        if (!is.list(groups)) {
            groups <- as.list(as.data.frame(combn(groups, 2, simplify = TRUE)))
        }
        for (spname in names(data)) {
            sdata <- data[[spname]]
            for (gp in groups) {
                gname <- paste(gp, collapse = " - ")
                gname <- ifelse(identical(spname, "..."), gname, paste(spname, gname, sep = ": "))
                spdata[[gname]] <- sdata[sdata[[group_by]] %in% gp, , drop = FALSE]
            }
        }
        if (is.null(nrow) && is.null(ncol) && isTRUE(byrow) && length(groups) > 1 && length(data) > 1) {
            nrow <- length(data)
        }
        data <- spdata
        rm(spdata)
    } else {
        groups <- unique(unlist(groups))
    }

    if (plot_type == "scatter") {
        plots <- lapply(names(data), function(nm) {
            DummyClonalScatterPlot(data[[nm]], nm, group_by, scatter_cor, scatter_size_by, ...)
        })

        combine_plots(plots, combine = combine, nrow = nrow, ncol = ncol, byrow = byrow)
    } else if (plot_type == "venn") {
        if (length(groups) > 4) {
            stop("Too many groups for venn plot. Please use 'upset' plot instead.")
        }
        data <- data[data[[group_by]] %in% groups, , drop = FALSE]
        data <- data[data$count > 0, , drop = FALSE]
        if (!is.null(split_by)) {
            data <- unite(data, ".split", split_by, sep = split_by_sep, remove = FALSE)
        }
        # Calculate the # singlets for each group
        label_fun <- function(df) {
            label <- c()
            for (i in 1:nrow(df)) {
                if (grepl("/", df$id[i], fixed = TRUE)) {
                    label <- c(label, df$count[i])
                } else {
                    indicator <- data$CloneID %in% df$item[[i]] &
                        data[[group_by]] == df$name[i] &
                        data$count == 1
                    if (!is.null(split_by)) {
                        indicator <- indicator & data$.split == df$.split[i]
                    }
                    ns <- data[indicator, , drop = FALSE]
                    label <- c(label, paste0(df$count[i], "\n(singlets: ", nrow(ns), ")"))
                }
            }
            label
        }

        VennDiagram(data,
            in_form = "long", id_by = "CloneID", group_by = group_by, label = label_fun,
            split_by = split_by, split_by_sep = split_by_sep, ...
        )
    } else if (plot_type == "upset") {
        data$fraction <- NULL
        data <- data[data[[group_by]] %in% groups, , drop = FALSE]
        data <- data[data$count > 0, , drop = FALSE]
        data <- data %>%
            pivot_wider(names_from = !!sym("group_by"), names_prefix = "count_", values_from = !!sym("count"), values_fill = 0) %>%
            mutate(across(
                .cols = starts_with("count_"),
                .names = "{.col} Singlet",
                .fns = ~ .x == 1
            )) %>%
            mutate(across(
                .cols = starts_with("count_") & !ends_with(" Singlet"),
                .names = "{.col} Expanded",
                .fns = ~ .x > 1
            )) %>%
            dplyr::select(!starts_with("count_") | ends_with(" Singlet") | ends_with(" Expanded")) %>%
            rename_with(
                function(x) substring(x, 7),
                .cols = starts_with("count_")
            )

        UpsetPlot(data, in_form = "wide", id_by = "CloneID", split_by = split_by, split_by_sep = split_by_sep, ...)
    }
}

#' ClonalCompositionPlot
#'
#' @description Plot the composition of the clones in different samples/groups.
#' @param data The product of [scRepertoire::combineTCR], [scRepertoire::combineTCR], or
#'  [scRepertoire::combineExpression].
#' @param clone_call How to call the clone - VDJC gene (gene), CDR3 nucleotide (nt),
#'  CDR3 amino acid (aa), VDJC gene + CDR3 nucleotide (strict) or a custom variable
#'  in the data
#' @param chain indicate if both or a specific chain should be used - e.g. "both",
#'  "TRA", "TRG", "IGH", "IGL"
#' @param method The method of plot to use. Default is "homeostasis".
#'  Possible values are "homeostasis", "homeo", "rel", "top", and "rare".
#'  * "homeostasis" - Plot the homeostasis/relative abundance of the clones. The `clone_split` will
#'    be the fraction of the clones in each sample.
#'  * "homeo" - Same as "homeostasis".
#'  * "rel" - Same as "homeostasis".
#'  * "top" - Plot the top clones. The `clone_split` will be indexes to cut the clones.
#'  * "rare" - Plot the rare clones. The `clone_split` will be the clone sizes.
#' @param clone_split The split for the clones. Default is NULL.
#'  * For "homeostasis", "homeo", "rel" - Default is `list(Rare = 1e-04, Small = 0.001, Medium = 0.01, Large = 0.1, Hyperexpanded = 1)`.
#'  * For "top" - Default is `c(10, 100, 1000, 10000, 30000, 100000)`.
#'  * For "rare" - Default is `c(1, 3, 10, 30, 100)`.
#' @param scale Whether to scale the values on the y-axis. Default is TRUE.
#'  * TRUE: The values of each group (on the x-axis) will be scaled to 1.
#'  * FALSE: No scaling.
#'  * "sample"/"Sample": The values in each sample will be scaled to 1.
#' @param facet_by The column name in the meta data to facet the plots. Default: NULL
#' @param group_by The column name in the meta data to group the cells. Default: NULL
#' @param split_by The column name in the meta data to split the plots. Default: NULL
#' @param xlab The x-axis label. Default is NULL.
#' @param ylab The y-axis label. Default is NULL.
#' @param plot_type The type of plot to use. Default is "bar".
#'  Possible values are "bar", "ring", "box", and "violin".
#' @param order The order of the x-axis items or groups. Default is an empty list.
#'  It should be a list of values. The names are the column names, and the values are the order.
#' @param ... Other arguments passed to the specific plot function.
#'  * For `bar` plot, see [plotthis::BarPlot()].
#'  * For `ring` plot, see [plotthis::RingPlot()].
#'  * For `box` plot, see [plotthis::BoxPlot()].
#'  * For `violin` plot, see [plotthis::ViolinPlot()].
#' @return A ggplot object or a list if `combine` is FALSE
#' @importFrom dplyr %>% mutate summarise ungroup n
#' @importFrom tidyr pivot_longer separate
#' @importFrom plotthis BarPlot RingPlot BoxPlot ViolinPlot
#' @importFrom scRepertoire clonalHomeostasis clonalProportion
#' @export
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
#' ClonalCompositionPlot(data)
#' ClonalCompositionPlot(data, method = "top")
#' ClonalCompositionPlot(data, plot_type = "ring")
#' ClonalCompositionPlot(data, group_by = "Type", plot_type = "box", comparison = TRUE)
#' ClonalCompositionPlot(data, group_by = "Type", plot_type = "violin", add_box = TRUE,
#'  add_bg = TRUE)
#' ClonalCompositionPlot(data, method = "rare")
#' }
ClonalCompositionPlot <- function(
    data, clone_call = "aa", chain = "both", method = c("homeostasis", "homeo", "rel", "top", "rare"),
    clone_split = NULL, scale = TRUE, facet_by = NULL, group_by = NULL, split_by = NULL,
    xlab = NULL, ylab = NULL, plot_type = c("bar", "ring", "box", "violin"), order = list(), ...
) {
    plot_type <- match.arg(plot_type)
    method <- match.arg(method)
    if (plot_type %in% c("box", "violin") && is.null(group_by)) {
        stop("'group_by' must be provided for box/violin ClonalCompositionPlot")
    }

    if (is.null(clone_split)) {
        if (method %in% c("homeostasis", "homeo", "rel")) {
            clone_split <- list(Rare = 1e-04, Small = 0.001, Medium = 0.01, Large = 0.1, Hyperexpanded = 1)
        } else if (method == "top") {
            # clone indexes
            clone_split <- c(10, 100, 1000, 10000, 30000, 100000)
        } else {  # rare
            # clone sizes
            clone_split <- c(1, 3, 10, 30, 100)
        }
    }
    all_groupings <- unique(c("Sample", group_by, facet_by, split_by))
    if (method == "homeostasis" || method == "homeo" || method == "rel" || method == "top") {
        data <- merge_clonal_groupings(data, all_groupings)
        if (method == "top") {
            data <- clonalProportion(data, cloneCall = clone_call, chain = chain,
                clonalSplit = clone_split, group.by = ".group", exportTable = TRUE)
        } else {
            data <- clonalHomeostasis(data, cloneCall = clone_call, chain = chain,
                cloneSize = clone_split, group.by = ".group", exportTable = TRUE)
        }
        data <- as.data.frame(data)
        data$.group <- rownames(data)
        # restore the groups
        data <- separate(data, ".group", into = all_groupings, sep = " // ")
        data <- pivot_longer(data, cols = -all_groupings, names_to = ".names", values_to = ".values")
        # Sample Type  .names        .values
        name_levels <- unique(data$.names)
    } else {  # rare
        data <- clonal_size_data(data, clone_call, chain, all_groupings)
        # CloneID Sample Type count fraction
        clone_split <- sort(clone_split)
        clone_split <- c(-Inf, clone_split, Inf)
        labels <- sapply(1:(length(clone_split) - 1), function(i) {
            if (clone_split[i] == -Inf && clone_split[i + 1] == 1) {
                "1"
            } else if (clone_split[i] == -Inf) {
                paste0("[1:", clone_split[i + 1], "]")
            } else if (i == length(clone_split) - 1) {
                paste0("[", clone_split[i], ":MAX]")
            } else {
                paste0("[", clone_split[i] + 1, ":", clone_split[i + 1], "]")
            }
        })
        data$.names <- cut(data$count, breaks = clone_split, labels = labels)
        name_levels <- levels(data$.names)
        data <- data %>%
            dplyr::group_by(!!!syms(setdiff(colnames(data), c("CloneID", "count", "fraction")))) %>%
            summarise(.values = n(), .groups = "drop")
    }

    for (group in all_groupings) {
        if (!is.null(order[[group]])) {
            data[[group]] <- factor(data[[group]], levels = order[[group]])
        }
    }

    if (isTRUE(scale)) {
        data <- data %>%
            dplyr::group_by(!!!syms(c(group_by, facet_by, split_by))) %>%
            mutate(.values = !!sym(".values") / sum(!!sym(".values")))
        ylab <- ylab %||% "Relative Abundance"
    } else if (identical(scale, "sample") || identical(scale, "Sample")) {
        data <- data %>%
            dplyr::group_by(!!!syms(c("Sample", group_by, facet_by, split_by))) %>%
            mutate(.values = !!sym(".values") / sum(!!sym(".values")))
        ylab <- ylab %||% "Relative Abundance"
    } else {
        ylab <- ylab %||% "Abundance"
    }
    if (method == "homeostasis" || method == "homeo" || method == "rel") {
        if (plot_type %in% c("bar", "ring") && !is.null(group_by)) {
            data <- data %>%
                dplyr::group_by(!!!syms(setdiff(colnames(data), c("Sample", ".values")))) %>%
                summarise(.values = sum(!!sym(".values")), .groups = "drop")
        }
        data$.names <- factor(data$.names, levels = name_levels)

        if (plot_type == "bar") {
            BarPlot(data, x = group_by %||% "Sample", y = ".values", group_by = ".names",
                group_name = "Clonal Group", position = "stack", facet_by = facet_by, split_by = split_by,
                xlab = xlab %||% group_by %||% "Sample", ylab = ylab %||% "Relative Abundance", ...)
        } else if (plot_type == "ring") {
            RingPlot(data, x = group_by %||% "Sample", y = ".values", group_by = ".names",
                group_name = "Clonal Group", facet_by = facet_by, split_by = split_by,
                xlab = xlab %||% group_by %||% "Sample", ylab = ylab %||% "Relative Abundance", ...)
        } else if (plot_type == "box") {
            BoxPlot(data, x = ".names", y = ".values", group_by = group_by,
                facet_by = facet_by, split_by = split_by,
                xlab = xlab %||% "Clonal Group", ylab = ylab %||% "Relative Abundance", ...)
        } else {  # violin
            ViolinPlot(data, x = ".names", y = ".values", group_by = group_by,
                facet_by = facet_by, split_by = split_by,
                xlab = xlab %||% "Clonal Group", ylab = ylab %||% "Relative Abundance", ...)
        }
    } else {
        group_name <- ifelse(method == "top", "Clonal Indices", "Clonal Sizes")

        if (plot_type %in% c("bar", "ring") && !is.null(group_by)) {
            data <- data %>%
                dplyr::group_by(!!!syms(setdiff(colnames(data), c("Sample", ".values")))) %>%
                summarise(.values = sum(!!sym(".values")), .groups = "drop")
        }
        data$.names <- factor(data$.names, levels = name_levels)

        if (plot_type == "bar") {
            BarPlot(data, x = group_by %||% "Sample", y = ".values", group_by = ".names",
                group_name = group_name, position = "stack", facet_by = facet_by, split_by = split_by,
                xlab = xlab %||% group_by %||% "Sample", ylab = ylab, ...)
        } else if (plot_type == "ring") {
            RingPlot(data, x = group_by %||% "Sample", y = ".values", group_by = ".names",
                group_name = group_name, facet_by = facet_by, split_by = split_by,
                xlab = xlab %||% group_by %||% "Sample", ylab = ylab, ...)
        } else if (plot_type == "box") {
            BoxPlot(data, x = ".names", y = ".values", group_by = group_by,
                facet_by = facet_by, split_by = split_by,
                xlab = xlab %||% group_name, ylab = ylab, ...)
        } else {  # violin
            ViolinPlot(data, x = ".names", y = ".values", group_by = group_by,
                facet_by = facet_by, split_by = split_by,
                xlab = xlab %||% group_name, ylab = ylab, ...)
        }
    }
}

#' ClonalOverlapPlot
#'
#' @description Plot the overlap of the clones in different samples/groups.
#' @param data The product of [scRepertoire::combineTCR], [scRepertoire::combineTCR], or
#'  [scRepertoire::combineExpression].
#' @param clone_call How to call the clone - VDJC gene (gene), CDR3 nucleotide (nt),
#'  CDR3 amino acid (aa), VDJC gene + CDR3 nucleotide (strict) or a custom variable
#'  in the data
#' @param chain indicate if both or a specific chain should be used - e.g. "both",
#'  "TRA", "TRG", "IGH", "IGL"
#' @param group_by The column name in the meta data to group the cells. Default: "Sample"
#' @param group_by_sep The separator used to concatenate the group_by when multiple columns are used.
#' @param full Whether to plot the full heatmap, or just a triangle. Default is TRUE.
#' @param split_by The column name in the meta data to split the plots. Default: NULL
#' @param order The order of the groups. Default is an empty list.
#'  It should be a list of values. The names are the column names, and the values are the order.
#' @param method The method to calculate the overlap. Default is "raw".
#'  * "overlap" - overlap coefficient
#'  * "morisita" - Morisita’s overlap index
#'  * "jaccard" - Jaccard index
#'  * "cosine" - cosine similarity
#'  * "raw" - exact number of overlapping clones
#'  See also [scRepertoire::clonalOverlap].
#' @param palette The color palette to use. Default is "Blues".
#' @param label_accuracy The accuracy of the labels. Default is NULL.
#'  If NULL, it will be 1 for "raw" and 0.01 for other methods.
#' @param label_cutoff The cutoff for the labels to show. Default is 1e-3.
#' @param cluster_rows Whether to cluster the rows. Default is FALSE.
#' @param cluster_columns Whether to cluster the columns. Default is FALSE.
#' @param show_row_names Whether to show the row names. Default is TRUE.
#' @param show_column_names Whether to show the column names. Default is TRUE.
#' @param ... Other arguments passed to the specific plot function [plotthis::Heatmap()].
#' @return A ComplexHeatmap object or a list if `combine` is FALSE
#' @importFrom stats as.dist
#' @importFrom rlang syms
#' @importFrom dplyr %>% filter
#' @importFrom tidyr separate pivot_longer unite
#' @importFrom scRepertoire clonalOverlap
#' @export
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
#' ClonalOverlapPlot(data)
#' ClonalOverlapPlot(data, clone_call = "strict", label_cutoff = 0,
#'   label_accuracy = 0.001, method = "morisita", full = FALSE)
#' ClonalOverlapPlot(data, group_by = c("Subject", "Type"))
#' ClonalOverlapPlot(data, group_by = "Type", split_by = "Subject")
#' }
ClonalOverlapPlot <- function(
    data, clone_call = "aa", chain = "both", group_by = "Sample", group_by_sep = "_", full = TRUE,
    split_by = NULL, order = list(), method = c("raw", "overlap", "morisita", "jaccard", "cosine"),
    palette = "Blues", label_accuracy = NULL, label_cutoff = 1e-3, cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = TRUE, show_column_names = TRUE, ...
) {
    method <- match.arg(method)

    all_groupings <- unique(c(group_by, split_by))
    data <- merge_clonal_groupings(data, all_groupings)
    data <- clonalOverlap(data, cloneCall = clone_call, chain = chain, group.by = ".group",
        method = method, exportTable = TRUE)
    if (isTRUE(full)) {
        data[lower.tri(data)] <- t(data)[lower.tri(data)]
    }
    #          B // P17 B // P18 B // P19
    # B // P17       NA        0    0.117
    # B // P18       NA       NA    0.001
    # B // P19       NA       NA       NA
    data$.group <- rownames(data)
    data <- data %>%
        separate(".group", into = all_groupings, sep = " // ") %>%
        pivot_longer(cols = -all_groupings, names_to = ".names", values_to = ".values") %>%
        separate(".names", into = paste(".names", all_groupings, sep = "_"), sep = " // ")

    if (!is.null(split_by)) {
        data <- data %>%
            unite(".split", split_by, sep = " // ", remove = FALSE) %>%
            unite(".names.split", !!!syms(paste(".names", split_by, sep = "_")), sep = " // ", remove = FALSE) %>%
            filter(!!sym(".split") == !!sym(".names.split"))

        data <- data[, setdiff(colnames(data), c(".split", ".names.split", paste(".names", split_by, sep = "_"))), drop = FALSE]
    }

    columns_by <- paste(group_by, collapse = group_by_sep)
    data <- data %>%
        unite("rows", !!!syms(paste(".names", group_by, sep = "_")), sep = group_by_sep) %>%
        unite(!!sym(columns_by), !!!syms(group_by), sep = group_by_sep)

    name <- switch(method,
        overlap = "Overlap Coefficient",
        morisita = "Morisita's Overlap Index",
        jaccard = "Jaccard Index",
        cosine = "Cosine Similarity",
        "# Overlap Clones"
    )
    label_accuracy <- label_accuracy %||% ifelse(method == "raw", 1, 0.01)

    clustering_distance <- function(m) {
        values <- m[upper.tri(m)]
        values <- 1 - scales::rescale(values, to = c(0, 1))
        m[upper.tri(m)] <- values
        m[lower.tri(m)] <- m[upper.tri(m)]
        diag(m) <- 0
        as.dist(m)
    }

    Heatmap(
        data, rows_by = "rows", columns_by = columns_by, values_by = ".values", split_by = split_by,
        clustering_distance_rows = function(m) { clustering_distance(t(m)) }, values_fill = 0,
        clustering_distance_columns = clustering_distance,
        # same rows_name as columns_name (inferred from columns_by) will cause an error
        rows_name = paste0(" ", columns_by),
        name = name, palette = palette, label = function(x) {
            ifelse(x > label_cutoff, scales::number(x, accuracy = label_accuracy), NA)
        },
        cluster_rows = cluster_rows, cluster_columns = cluster_columns, cell_type = "label",
        show_row_names = show_row_names, show_column_names = show_column_names, ...)
}
