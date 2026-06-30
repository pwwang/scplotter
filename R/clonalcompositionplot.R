#' Clonal Volume Plot
#'
#' @description
#' Visualizes the number (or fraction) of unique T-cell or B-cell clones across
#' samples and metadata groups. Clonal volume — the count of distinct clonotypes
#' detected in a sample — is a fundamental measure of immune repertoire diversity.
#' Higher clonal volume indicates a more diverse repertoire, while lower volume
#' may reflect clonal expansion in response to antigen stimulation.
#'
#' `ClonalVolumePlot` computes clonal counts via
#' \code{\link[scRepertoire:clonalQuant]{scRepertoire::clonalQuant()}} and
#' visualizes them as bar, box, or violin plots. It accepts both
#' \pkg{scRepertoire} combined TCR/BCR data and Seurat objects with clonal
#' information integrated via \code{scRepertoire::combineExpression()}.
#'
#' @param data The product of \code{\link[scRepertoire:combineTCR]{scRepertoire::combineTCR()}},
#'   \code{\link[scRepertoire:combineBCR]{scRepertoire::combineBCR()}}, or
#'   \code{\link[scRepertoire:combineExpression]{scRepertoire::combineExpression()}}.
#' @param clone_call How to define a clone. One of:
#'   \itemize{
#'     \item `"gene"` — V(D)JC gene combination
#'     \item `"nt"` — CDR3 nucleotide sequence
#'     \item `"aa"` — CDR3 amino acid sequence (default)
#'     \item `"strict"` — V(D)JC gene + CDR3 nucleotide
#'   }
#'   Or a custom variable name in the data.
#' @param chain Which chain(s) to use: `"both"` (default), `"TRA"`, `"TRB"`,
#'   `"TRD"`, `"TRG"`, `"IGH"`, or `"IGL"`.
#' @param scale Logical; if `TRUE`, values are scaled to clone proportion
#'   (fraction of unique clones) instead of absolute clone counts. Default
#'   is `FALSE`.
#' @param plot_type The visualization type. One of `"bar"` (default),
#'   `"box"`, or `"violin"`. For `"box"` and `"violin"`, the data is broken
#'   down by Sample and grouped by `group_by` to show per-sample
#'   distributions.
#' @param x The metadata column used as the x-axis. Default is `"Sample"`.
#' @param ylab Y-axis label. Default is `NULL`, which auto-generates
#'   `"Number of Unique Clones"` or `"Fraction of Unique Clones"` depending
#'   on `scale`.
#' @param group_by Metadata column used to group (color) the data. Default
#'   is `NULL`.
#' @param facet_by Metadata column used to facet the plot into separate
#'   panels. Default is `NULL`.
#' @param split_by Metadata column used to split the data into separate
#'   plots. Default is `NULL`.
#' @param order A named list controlling the order of factor levels. List
#'   names are column names; list values are the desired order. Default is
#'   `NULL` (use existing factor levels or alphabetical order).
#' @param ... Additional arguments passed to the underlying \pkg{plotthis}
#'   function:
#'   \itemize{
#'     \item `"bar"` — \code{\link[plotthis:BarPlot]{plotthis::BarPlot()}}
#'           (`position`, `palette`, `fill_by`, ...)
#'     \item `"box"` — \code{\link[plotthis:BoxPlot]{plotthis::BoxPlot()}}
#'           (`comparisons`, `alpha`, `palette`, ...)
#'     \item `"violin"` — \code{\link[plotthis:ViolinPlot]{plotthis::ViolinPlot()}}
#'           (`add_box`, `comparisons`, `palette`, ...)
#'   }
#' @return A `ggplot` object, or a list of `ggplot` objects if
#'   `combine = FALSE` is passed via `...`.
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
#'     variables = factor(sample(c("B", "L"), 8, replace = TRUE), levels = c("B", "L"))
#' )
#' data <- scRepertoire::addVariable(data,
#'     variable.name = "Sex",
#'     variables = factor(sample(c("M", "F"), 8, replace = TRUE), levels = c("M", "F"))
#' )
#'
#' ClonalVolumePlot(data)
#' ClonalVolumePlot(data, x = "Type")
#' ClonalVolumePlot(data, x = "Type", order = list(Type = c("L", "B")))
#' ClonalVolumePlot(data, x = c("Type", "Sex"), scale = TRUE, fill_by = "Type")
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
    facet_by = NULL, split_by = NULL, order = NULL, ylab = NULL, ...) {
    plot_type <- match.arg(plot_type)
    if (plot_type %in% c("box", "violin")) {
        all_groupings <- unique(c("Sample", x, group_by, facet_by, split_by))
    } else {
        all_groupings <- unique(c(x, group_by, facet_by, split_by))
    }
    grouping_levels <- get_clonal_grouping_levels(data, all_groupings, order)
    data <- merge_clonal_groupings(data, all_groupings)
    data <- clonalQuant(data,
        cloneCall = clone_call, chain = chain, scale = scale,
        group.by = ".group", exportTable = TRUE
    )
    # restore the groups
    data <- separate(data, ".group", into = all_groupings, sep = " // ")
    for (gl in names(grouping_levels)) {
        if (!is.null(data[[gl]])) {
            data[[gl]] <- factor(data[[gl]], levels = grouping_levels[[gl]])
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


#' Clonal Abundance Plot
#'
#' @description
#' Visualizes the distribution of clonal abundances — how many clones are
#' present at each abundance level (frequency) in the repertoire. Clonal
#' abundance distributions typically follow a power-law pattern: a small
#' number of highly expanded clones and a large number of rare clones.
#' This function helps characterize repertoire structure by showing whether
#' the immune response is dominated by a few large clones (clonal expansion)
#' or evenly distributed across many clones (high diversity).
#'
#' `ClonalAbundancePlot` computes clonal abundance data via
#' \code{\link[scRepertoire:clonalAbundance]{scRepertoire::clonalAbundance()}}
#' and visualizes it as trend lines, histograms, or density curves.
#'
#' @param data The product of \code{\link[scRepertoire:combineTCR]{scRepertoire::combineTCR()}},
#'   \code{\link[scRepertoire:combineBCR]{scRepertoire::combineBCR()}}, or
#'   \code{\link[scRepertoire:combineExpression]{scRepertoire::combineExpression()}}.
#' @param clone_call How to define a clone. One of `"gene"`, `"nt"`,
#'   `"aa"` (default), `"strict"`, or a custom variable name in the data.
#' @param chain Which chain(s) to use: `"both"` (default), `"TRA"`, `"TRB"`,
#'   `"TRD"`, `"TRG"`, `"IGH"`, or `"IGL"`.
#' @param plot_type The visualization type. One of:
#'   \itemize{
#'     \item `"trend"` (default) — Smoothed trend line showing the number
#'           of clones at each abundance level. The x-axis is transformed
#'           by `xtrans` (default log10), and a LOESS trend is fitted.
#'     \item `"histogram"` — Binned histogram of clonal abundances.
#'           Optionally overlay a trend line with `add_trend = TRUE`.
#'     \item `"density"` — Kernel density estimate of the abundance
#'           distribution.
#'   }
#' @param trend_skip_zero Logical; if `TRUE` (default), zero-abundance bins
#'   are excluded from the trend line fit. Improves fit quality when many
#'   abundance bins have zero clones.
#' @param binwidth The histogram bin width (in log10-transformed abundance
#'   units). Default is `0.1`. Only used for `"trend"` and `"histogram"`
#'   plot types.
#' @param group_by Metadata column used to group (color) the data. Default
#'   is `"Sample"`.
#' @param group_by_sep Separator used when concatenating multiple `group_by`
#'   columns. Default is `"_"`.
#' @param facet_by Metadata column used to facet the plot into separate
#'   panels. Default is `NULL`.
#' @param split_by Metadata column used to split the data into separate
#'   plots. Default is `NULL`.
#' @param order A named list controlling the order of factor levels. List
#'   names are column names; list values are the desired order. Default is
#'   `NULL`.
#' @param bw Smoothing bandwidth for density plots. Higher values produce
#'   smoother curves. Default is `0.5`.
#' @param xlab X-axis label. Default is `"Abundance"`.
#' @param ylab Y-axis label. Default is `NULL`, which auto-generates
#'   `"Number of Clones"` (trend/histogram) or `"Density of Clones"`
#'   (density).
#' @param xtrans Transformation applied to the x-axis. Default is `"log10"`,
#'   which spreads low-abundance clones for better visibility. Use
#'   `"identity"` for linear scale.
#' @param ytrans Transformation applied to the y-axis. Default is
#'   `"identity"`. Use `"log10"` to better visualize distributions spanning
#'   multiple orders of magnitude.
#' @param theme_args A list of theme elements passed to the underlying
#'   \pkg{plotthis} function. Default is an empty list.
#' @param ... Additional arguments passed to the underlying \pkg{plotthis}
#'   function:
#'   \itemize{
#'     \item `"trend"` — \code{\link[plotthis:Histogram]{plotthis::Histogram()}}
#'           (with `use_trend = TRUE`; `palette`, `alpha`, ...)
#'     \item `"histogram"` — \code{\link[plotthis:Histogram]{plotthis::Histogram()}}
#'           (`add_trend`, `palette`, `alpha`, ...)
#'     \item `"density"` — \code{\link[plotthis:DensityPlot]{plotthis::DensityPlot()}}
#'           (`palette`, `alpha`, ...)
#'   }
#' @return A `ggplot` object, or a list of `ggplot` objects if
#'   `combine = FALSE` is passed via `...`.
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
#'     variables = factor(sample(c("B", "L"), 8, replace = TRUE), levels = c("L", "B"))
#' )
#' data <- scRepertoire::addVariable(data,
#'     variable.name = "Sex",
#'     variables = factor(sample(c("M", "F"), 8, replace = TRUE), levels = c("M", "F"))
#' )
#'
#' ClonalAbundancePlot(data)
#' ClonalAbundancePlot(data, ytrans = "log10")
#' ClonalAbundancePlot(data, plot_type = "histogram")
#' ClonalAbundancePlot(data, plot_type = "histogram", add_trend = TRUE, trend_skip_zero = TRUE)
#' ClonalAbundancePlot(data, plot_type = "density")
#' ClonalAbundancePlot(data, group_by = "Type")
#' }
ClonalAbundancePlot <- function(
    data, clone_call = "aa", chain = "both", xtrans = "log10", ytrans = "identity",
    plot_type = c("trend", "histogram", "density"), binwidth = 0.1, trend_skip_zero = TRUE,
    bw = 0.5, group_by = "Sample", group_by_sep = "_", facet_by = NULL, split_by = NULL,
    order = NULL, xlab = "Abundance", ylab = NULL, theme_args = list(), ...) {
    plot_type <- match.arg(plot_type)

    all_groupings <- unique(c(group_by, facet_by, split_by))
    grouping_levels <- get_clonal_grouping_levels(data, all_groupings, order)
    data <- merge_clonal_groupings(data, all_groupings)
    if (length(all_groupings) > 0) {
        data <- clonalAbundance(data,
            cloneCall = clone_call, chain = chain,
            group.by = ".group", exportTable = TRUE
        )

        # restore the groups
        data <- separate(data, ".group", into = all_groupings, sep = " // ")
        for (gl in names(grouping_levels)) {
            if (!is.null(data[[gl]])) {
                data[[gl]] <- factor(data[[gl]], levels = grouping_levels[[gl]])
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

#' Clonal CDR3 Length Plot
#'
#' @description
#' Visualizes the distribution of CDR3 sequence lengths across the immune
#' repertoire. CDR3 length is a key feature of T-cell and B-cell receptor
#' diversity — different clones have different CDR3 lengths, and shifts in
#' length distribution can indicate clonal selection, antigen-specific
#' expansion, or repertoire bias.
#'
#' `ClonalLengthPlot` computes CDR3 length data via
#' \code{\link[scRepertoire:clonalLength]{scRepertoire::clonalLength()}} and
#' visualizes the distribution as bar, box, violin, or density plots. Length
#' is measured in amino acids (when `clone_call = "aa"`) or nucleotides
#' (when `clone_call = "nt"`).
#'
#' @param data The product of \code{\link[scRepertoire:combineTCR]{scRepertoire::combineTCR()}},
#'   \code{\link[scRepertoire:combineBCR]{scRepertoire::combineBCR()}}, or
#'   \code{\link[scRepertoire:combineExpression]{scRepertoire::combineExpression()}}.
#' @param clone_call How to define a clone. Only `"nt"` (CDR3 nucleotide
#'   length) or `"aa"` (CDR3 amino acid length, default) are supported.
#' @param chain Which chain(s) to use: `"both"` (default), `"TRA"`, `"TRB"`,
#'   `"TRD"`, `"TRG"`, `"IGH"`, or `"IGL"`.
#' @param plot_type The visualization type. One of `"bar"` (default),
#'   `"box"`, `"violin"`, or `"density"`.
#'   \itemize{
#'     \item `"bar"` — Bar chart of clone counts at each CDR3 length.
#'           Empty length bins (zero clones) are padded with zeros to
#'           maintain a continuous x-axis.
#'     \item `"box"` — Box plot of per-group length distributions.
#'     \item `"violin"` — Violin plot of per-group length distributions.
#'     \item `"density"` — Kernel density estimate of the length
#'           distribution, using raw (unaggregated) data.
#'   }
#' @param x_nbreaks Number of x-axis breaks for the bar plot. Default is
#'   `10`. Breaks are computed as quantiles of the length range and rounded
#'   to the nearest 10.
#' @param group_by Metadata column used to group (color) the data. Default
#'   is `"Sample"`.
#' @param facet_by Metadata column used to facet the plot into separate
#'   panels. Default is `NULL`.
#' @param split_by Metadata column used to split the data into separate
#'   plots. Default is `NULL`.
#' @param order A named list controlling the order of factor levels. List
#'   names are column names; list values are the desired order. Default is
#'   `NULL`.
#' @param xlab X-axis label. Default is `"Length"`.
#' @param ylab Y-axis label. Default is `NULL`, which auto-generates
#'   `"Number of CDR3 (AA)"` or `"Number of CDR3 (NT)"` based on
#'   `clone_call`.
#' @param position Bar position for the bar plot. One of `"dodge"`
#'   (default), `"stack"`, or `"fill"`.
#' @param ... Additional arguments passed to the underlying \pkg{plotthis}
#'   function:
#'   \itemize{
#'     \item `"bar"` — \code{\link[plotthis:BarPlot]{plotthis::BarPlot()}}
#'           (`palette`, `alpha`, `position_dodge_preserve`, ...)
#'     \item `"box"` — \code{\link[plotthis:BoxPlot]{plotthis::BoxPlot()}}
#'           (`comparisons`, `alpha`, `palette`, ...)
#'     \item `"violin"` — \code{\link[plotthis:ViolinPlot]{plotthis::ViolinPlot()}}
#'           (`add_box`, `comparisons`, `palette`, ...)
#'     \item `"density"` — \code{\link[plotthis:DensityPlot]{plotthis::DensityPlot()}}
#'           (`palette`, `alpha`, `bw`, ...)
#'   }
#' @return A `ggplot` object, or a list of `ggplot` objects if
#'   `combine = FALSE` is passed via `...`.
#' @importFrom stats quantile
#' @importFrom rlang syms
#' @importFrom dplyr summarise n n_distinct
#' @importFrom tidyr separate
#' @importFrom SeuratObject Idents<-
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
#'  variables = factor(sample(c("B", "L"), 8, replace = TRUE), levels = c("L", "B")))
#' data <- scRepertoire::addVariable(data, variable.name = "Sex",
#'  variables = factor(sample(c("M", "F"), 8, replace = TRUE), levels = c("M", "F")))
#'
#' ClonalLengthPlot(data)
#' ClonalLengthPlot(data, plot_type = "box")
#' ClonalLengthPlot(data, clone_call = "nt", plot_type = "violin", chain = "TRB",
#'  group_by = "Type", comparisons = TRUE)
#' ClonalLengthPlot(data, plot_type = "density", chain = "TRA")
#' }
ClonalLengthPlot <- function(
    data, clone_call = "aa", chain = "both", plot_type = c("bar", "box", "violin", "density"),
    x_nbreaks = 10, group_by = "Sample", order = NULL, xlab = "Length", ylab = NULL,
    position = "dodge", facet_by = NULL, split_by = NULL, ...) {
    plot_type <- match.arg(plot_type)

    if (plot_type %in% c("box", "violin")) {
        all_groupings <- unique(c("Sample", group_by, facet_by, split_by))
    } else {
        all_groupings <- unique(c(group_by, facet_by, split_by))
    }
    grouping_levels <- get_clonal_grouping_levels(data, all_groupings, order)
    data <- merge_clonal_groupings(data, all_groupings)

    if (identical(all_groupings, "Sample")) {
        if (inherits(data, "Seurat")) {
            Idents(data) <- "Sample"
        }
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
    for (gl in names(grouping_levels)) {
        if (!is.null(data[[gl]])) {
            data[[gl]] <- factor(data[[gl]], levels = grouping_levels[[gl]])
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
        args$position_dodge_preserve <- args$position_dodge_preserve %||% "single"
        args$xlab <- "Length"
        args$ylab <- ylab %||% default_ylab
        args$theme_args <- args$theme_args %||% list()
        args$theme_args$panel.grid.major.x <- element_blank()
        minx <- min(data$length, na.rm = TRUE)
        maxx <- max(data$length, na.rm = TRUE)
        # Make both minx and maxx to the nearest 10s
        minx <- floor(minx / 10) * 10
        maxx <- ceiling(maxx / 10) * 10
        breaks <- unname(
            quantile(
                minx:maxx,
                probs = seq(0, 1, length.out = x_nbreaks + 1),
                type = 3
            )
        )
        missing_breaks <- setdiff(breaks, data$length)
        missing_data <- data.frame(
            length = missing_breaks,
            .n = 0
        )
        for (col in setdiff(colnames(data), c("length", ".n"))) {
            missing_data[[col]] <- data[[col]][1]
        }
        args$data <- rbind(data, missing_data)
        args$data <- args$data[order(args$data$length), , drop = FALSE]
        if (is.null(args$split_by)) {
            suppressMessages({
                do_call(BarPlot, args) + scale_x_discrete(breaks = breaks)
            })
        } else {
            suppressMessages({
                do_call(BarPlot, args) & scale_x_discrete(breaks = breaks)
            })
        }
    } else if (plot_type == "box" || plot_type == "violin") {
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
        args$legend.position <- args$legend.position %||% ifelse(
            is.null(args$split_by),
            ifelse(!is.null(args$group_by) && n_distinct(data[[group_by]]) < 10, "right", "none"),
            "right"
        )
        fn <- if (plot_type == "box") BoxPlot else ViolinPlot
        do_call(fn, args)
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
DummyClonalScatterPlot <- function(df, group_by, scatter_cor, size_by, ...) {
    if (nrow(df) == 0) {
        stop("[ClonalResidencyPlot] No data found for the group_by: ", group_by,
            ". Did you specify the correct 'group_by'/'groups' parameters?")
    }
    pair <- rev(levels(df[[group_by]]))

    if (length(pair) != 2) {
        stop("[ClonalResidencyPlot] The group_by should have exactly 2 unique values, got: ", length(pair))
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
    sum_counts_max <- max(sum_counts1, sum_counts2)

    counts1_norm <- jitter(1 + df[[suf1]], amount = 0.25) / sum_counts_max
    counts2_norm <- jitter(1 + df[[suf2]], amount = 0.25) / sum_counts_max

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
        1 / sum_counts_max,
        0.001 + 1 / sum_counts_max,
        0.01 + 1 / sum_counts_max,
        0.1 + 1 / sum_counts_max
    )
    ybreaks <- c(
        1 / sum_counts_max,
        0.001 + 1 / sum_counts_max,
        0.01 + 1 / sum_counts_max,
        0.1 + 1 / sum_counts_max
    )
    minx <- min(plotdata[[pair[1]]], na.rm = TRUE)
    maxx <- max(plotdata[[pair[1]]], na.rm = TRUE)
    miny <- min(plotdata[[pair[2]]], na.rm = TRUE)
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

    args <- rlang::dots_list(...)
    args$palette <- args$palette %||% "turbo"
    args$data <- plotdata
    args$x <- pair[1]
    args$y <- pair[2]
    args$color_by <- "NumType"
    args$title <- if (identical(args$title, "...")) NULL else args$title
    args$subtitle <- subtitle
    args$border_color <- TRUE
    args$size_by <- ifelse(size_by == "max", "Max_Size", "Total_Size")
    args$size_name <- ifelse(size_by == "max", "Max Size", "Total Size")
    args$aspect.ratio <- 1

    colfun <- colorRamp2(
        seq(min(plotdata$NumType, na.rm = TRUE), max(plotdata$NumType, na.rm = TRUE), length.out = 100),
        palette_this(palette = args$palette, palcolor = args$palcolor)
    )
    label_df <- plotdata %>% group_by(!!sym("TypeName")) %>%
        summarise(
            x = first(!!sym(pair[1])),
            y = first(!!sym(pair[2])),
            Type = mean(!!sym("Type"), na.rm = TRUE), .groups = "drop")
    label_df$TypeName <- factor(label_df$TypeName, levels = labels)
    label_df <- label_df[order(label_df$TypeName), , drop = FALSE]

    diag_df <- data.frame(
        # diagnal, horizontal, vertical, horizontal short, vertical short
        x = c(1.5 / sum_counts_max, minx, 1.5 / sum_counts_max, minx, 2.5 / sum_counts_max),
        xend = c(maxx, maxx, 1.5 / sum_counts_max, 1.5 / sum_counts_max, 2.5 / sum_counts_max),
        y = c(1.5 / sum_counts_max, 1.5 / sum_counts_max, miny, 2.5 / sum_counts_max, miny),
        yend = c(maxy, 1.5 / sum_counts_max, maxy, 2.5 / sum_counts_max, 1.5 / sum_counts_max)
    )
    p <- do_call(ScatterPlot, args) +
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
            data = diag_df,
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

#' Clonal Residency Plot
#'
#' @description
#' Visualizes the sharing (residency) of T-cell or B-cell clones across
#' different samples or metadata groups. Clonal residency analysis reveals
#' how clonotypes are distributed — whether a clone is private to one
#' condition or shared across multiple conditions — which is critical for
#' understanding immune responses, tracking antigen-specific clones, and
#' identifying public vs. private repertoires.
#'
#' `ClonalResidencyPlot` supports three visualization modes:
#' \itemize{
#'   \item **Scatter plot** — Compares clone sizes between two groups on
#'         log-transformed axes. Points are colored by clonal category:
#'         singletons (unique to one group), expanded clones, and dual
#'         clones (shared between groups). Correlation statistics are
#'         displayed in the subtitle.
#'   \item **Venn diagram** — Shows the overlap of clone sets between up
#'         to 4 groups. When `with_class = TRUE`, labels include singlet
#'         counts.
#'   \item **UpSet plot** — Shows intersection sizes for any number of
#'         groups. When `with_class = TRUE`, clone classes (singlet,
#'         expanded) are displayed as separate intersections.
#' }
#'
#' @section The groups parameter:
#' The `groups` parameter controls which groups are compared and how they
#' are displayed:
#' \itemize{
#'   \item **`NULL` (default):** All levels of the `group_by` column are
#'         used. For scatter plots, all pairwise combinations are plotted.
#'   \item **Character vector:** Specifies a subset of groups to include.
#'         For scatter plots, the specified pairs are plotted individually;
#'         use `":"` notation for explicit pairings (e.g., `c("L:B", "Y:X")`
#'         compares L vs B and Y vs X).
#'   \item **Named vector/list:** The names are used as display labels and
#'         the values match groups in the data. For example,
#'         `c(B = "P17B", L = "P17L")` labels groups as "B" and "L".
#'         For scatter plots, use `c("L:B" = "group1:group2")` where the
#'         first group in the pair is on the y-axis and the second is on
#'         the x-axis.
#' }
#'
#' @param data The product of \code{\link[scRepertoire:combineTCR]{scRepertoire::combineTCR()}},
#'   \code{\link[scRepertoire:combineBCR]{scRepertoire::combineBCR()}}, or
#'   \code{\link[scRepertoire:combineExpression]{scRepertoire::combineExpression()}}.
#' @param clone_call How to define a clone. One of `"gene"`, `"nt"`,
#'   `"aa"` (default), `"strict"`, or a custom variable name in the data.
#' @param chain Which chain(s) to use: `"both"` (default), `"TRA"`, `"TRB"`,
#'   `"TRD"`, `"TRG"`, `"IGH"`, or `"IGL"`.
#' @param plot_type The visualization type. One of `"scatter"` (default),
#'   `"venn"`, or `"upset"`.
#' @param group_by Metadata column used to define the groups being compared.
#'   Default is `"Sample"`. Multiple columns are concatenated using
#'   `group_by_sep`.
#' @param group_by_sep Separator used when concatenating multiple `group_by`
#'   columns. Default is `"_"`.
#' @param groups The groups to compare. See the \emph{The groups parameter}
#'   section for detailed usage. Default is `NULL` (all groups).
#' @param facet_by Not supported for `ClonalResidencyPlot`. Must be `NULL`.
#' @param split_by Metadata column used to split the data into separate
#'   plots. Default is `NULL`.
#' @param split_by_sep Separator used when concatenating multiple `split_by`
#'   columns. Default is `"_"`.
#' @param scatter_cor Correlation method for scatter plots. One of
#'   `"pearson"` (default), `"spearman"`, or `"kendall"`. Correlation is
#'   computed on log-transformed clone sizes of dual (shared) clones.
#' @param scatter_size_by How point sizes are determined in scatter plots.
#'   \itemize{
#'     \item `"max"` (default) — Size reflects the larger clone size
#'           between the two groups.
#'     \item `"total"` — Size reflects the sum of clone sizes across both
#'           groups.
#'   }
#' @param with_class Logical; if `TRUE` (default), clonal class information
#'   (singlet vs. expanded) is included in Venn and UpSet plot labels. Only
#'   applicable for `"venn"` and `"upset"` plot types.
#' @param order A named list controlling the order of factor levels. List
#'   names are column names; list values are the desired order. Default is
#'   `NULL`.
#' @param combine Logical; if `TRUE` (default), multiple plots are combined
#'   into a single layout using \code{\link[plotthis:combine_plots]{plotthis::combine_plots()}}.
#' @param nrow Number of rows in the combined plot layout. Default is `NULL`
#'   (auto-determined).
#' @param ncol Number of columns in the combined plot layout. Default is
#'   `NULL` (auto-determined).
#' @param byrow Logical; if `TRUE` (default), the combined layout is filled
#'   row by row.
#' @param ... Additional arguments passed to the underlying \pkg{plotthis}
#'   function:
#'   \itemize{
#'     \item `"scatter"` — \code{\link[plotthis:ScatterPlot]{plotthis::ScatterPlot()}}
#'           (`palette`, `palcolor`, `title`, ...)
#'     \item `"venn"` — \code{\link[plotthis:VennDiagram]{plotthis::VennDiagram()}}
#'           (`palette`, `alpha`, ...)
#'     \item `"upset"` — \code{\link[plotthis:UpsetPlot]{plotthis::UpsetPlot()}}
#'           (`palette`, ...)
#'   }
#' @return A `ggplot` object, or a list of `ggplot` objects if
#'   `combine = FALSE`.
#' @note
#' **Scatter plot group limits:** Each scatter plot compares exactly two
#' groups. When more than two groups exist and `groups` is not specified
#' with `":"` notation, all pairwise combinations are generated
#' automatically.
#'
#' **Venn diagram group limits:** Venn diagrams support at most 4 groups.
#' For more groups, use `plot_type = "upset"` instead.
#'
#' @export
#' @importFrom utils combn
#' @importFrom dplyr distinct rename_with select across starts_with ends_with %>% filter
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
#'     variables = factor(rep(c("B", "L", "X", "Y"), 2), levels = c("Y", "B", "L", "X"))
#' )
#' data <- scRepertoire::addVariable(data,
#'     variable.name = "Subject",
#'     variables = rep(c("P17", "P18", "P19", "P20"), each = 2)
#' )
#'
#' ClonalResidencyPlot(data, groups = c("P18B", "P18L"))
#' ClonalResidencyPlot(data, group_by = "Type", groups = c("L", "B"),
#'  split_by = "Subject")
#' ClonalResidencyPlot(data, group_by = "Type", groups = c("L:B", "Y:X"),
#'  split_by = "Subject")
#' ClonalResidencyPlot(data, plot_type = "venn", groups = c("B", "L"),
#'  group_by = "Type", split_by = "Subject")
#' ClonalResidencyPlot(data, groups = c("P17_B", "P17_L", "P18_X", "P18_Y"),
#'    plot_type = "venn", with_class = TRUE, palette = "Blues",
#'    group_by = c("Subject", "Type"))
#' ClonalResidencyPlot(data, plot_type = "upset", groups = c("P18B", "P18L"))
#' ClonalResidencyPlot(data, plot_type = "upset", with_class = FALSE)
#' }
ClonalResidencyPlot <- function(
    data, clone_call = "aa", chain = "both", plot_type = c("scatter", "venn", "upset"),
    group_by = "Sample", group_by_sep = "_", groups = NULL, facet_by = NULL, split_by = NULL, split_by_sep = "_",
    scatter_cor = "pearson", scatter_size_by = c("max", "total"), with_class = TRUE,
    order = NULL, combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, ...
) {

    stopifnot("[ClonalResidencyPlot] 'facet_by' is not supported" = is.null(facet_by))

    plot_type <- match.arg(plot_type)
    scatter_size_by <- match.arg(scatter_size_by)

    all_groupings <- unique(c(group_by, facet_by, split_by))
    data <- clonal_size_data(data, clone_call, chain, all_groupings, order)

    check_columns <- utils::getFromNamespace("check_columns", "plotthis")
    orig_group_by <- group_by
    group_by <- check_columns(
        data, group_by, force_factor = TRUE,
        allow_multi = TRUE, concat_multi = TRUE, concat_sep = group_by_sep
    )
    all_groupings <- unique(c(group_by, facet_by, split_by))
    if (length(orig_group_by) > 1) {
        for (og in orig_group_by) {
            data[[og]] <- NULL
        }
    }

    # restore the groups
    for (group in all_groupings) {
        if (!is.null(order[[group]])) {
            data[[group]] <- factor(data[[group]], levels = order[[group]])
        }
    }

    # handle the groups
    # the order of groups should override the order in the data, even specified by `order`
    if (is.null(groups)) {
        groups <- levels(data[[group_by]])
    } else {
        groups <- unlist(groups)
        if (is.null(names(groups))) {
            names(groups) <- groups
        }

        if (
            identical(plot_type, "scatter") &&
            any(grepl(":", c(names(groups), groups), fixed = TRUE)) &&
            !all(grepl(":", c(names(groups), groups), fixed = TRUE))) {
            stop("[ClonalResidencyPlot] All 'groups' should contain ':' for scatter plot if any of them does.")
        }

        if (identical(plot_type, "scatter") && any(grepl(":", groups, fixed = TRUE))) {
            recoded_groups <- list()
            # recode the groups to the format of "group1:group2"
            for (old_group in names(groups)) {
                new_group <- groups[old_group]
                old_groups <- strsplit(old_group, ":", fixed = TRUE)[[1]]
                new_groups <- strsplit(new_group, ":", fixed = TRUE)[[1]]
                if (length(old_groups) != 2 || length(new_groups) != 2) {
                    stop("[ClonalResidencyPlot] The name and value of 'groups' parameter should contain exactly two values separated by ':' for scatter plot.")
                }
                old_groups <- trimws(old_groups)
                new_groups <- trimws(new_groups)
                for (i in 1:2) {
                    if (!old_groups[i] %in% data[[group_by]]) {
                        stop(paste0("[ClonalResidencyPlot] '", old_groups[i], "' in the values of 'groups' parameter are not present in the data."))
                    }
                    if (is.null(recoded_groups[[new_groups[i]]])) {
                        recoded_groups[[new_groups[i]]] <- old_groups[i]
                    } else {
                        if (!identical(recoded_groups[[new_groups[i]]], old_groups[i])) {
                            stop(paste0("[ClonalResidencyPlot] The name '", new_groups[i], "' in 'groups' parameter is mapped to different groups in the data: '", paste(recoded_groups[[new_groups[i]]], collapse = ", "), "' and '", old_groups[i], "'."))
                        }
                    }
                }
            }
            data[[group_by]] <- forcats::fct_recode(data[[group_by]], !!!recoded_groups)
            data[[group_by]] <- factor(data[[group_by]], levels = rev(unique(names(recoded_groups))))
        } else {
            # Handle case where names are single-column names (without : or even if they contain : but not in scatter plot)
            non_existing_groups <- unique(setdiff(names(groups), data[[group_by]]))
            if (length(non_existing_groups) > 0) {
                stop(paste0("[ClonalResidencyPlot] '", paste(non_existing_groups, collapse = ", "), "' in the names of 'groups' parameter are not present in the data."))
            }
            data[[group_by]] <- forcats::fct_recode(data[[group_by]], !!!stats::setNames(names(groups), groups))
            data[[group_by]] <- factor(data[[group_by]], levels = unname(groups))
        }
    }

    if (plot_type == "scatter") {
        if (!is.null(split_by)) {
            # data <- unite(data, ".split", split_by, sep = split_by_sep, remove = FALSE)
            # keep the order of the splits in the data
            split_by <- check_columns(
                data, split_by,
                force_factor = TRUE, allow_multi = TRUE, concat_multi = TRUE, concat_sep = split_by_sep
            )
        } else {
            data$.split <- "..."
            split_by <- ".split"
        }
        data <- split(data, data[[split_by]])
        if (any(grepl(":", groups, fixed = TRUE))) {
            groups <- lapply(groups, function(g) {
                gs <- strsplit(g, ":", fixed = TRUE)[[1]]
                if (length(gs) != 2) {
                    stop("[ClonalResidencyPlot] The 'groups' parameter should contain exactly two values separated by ':' for scatter plot.")
                }
                gs <- trimws(gs)
                gs
            })
        } else {
            groups <- as.list(as.data.frame(combn(groups, 2, simplify = TRUE)))
        }

        plots <- list()
        for (nm in names(data)) {
            for (g in groups) {
                d <- data[[nm]] %>% filter(!!sym(group_by) %in% g)
                if (length(unique(as.character(d[[group_by]]))) != 2) {
                    warning(paste0("[ClonalResidencyPlot] Not both groups '", paste(g, collapse = ", "), "' is not present in the data for '", nm, "'. Skipping."))
                    next
                }
                d[[group_by]] <- factor(d[[group_by]], levels = g)
                plots[[length(plots) + 1]] <- DummyClonalScatterPlot(
                    df = d,
                    title = nm,
                    group_by = group_by,
                    scatter_cor = scatter_cor,
                    size_by = scatter_size_by,
                    ...
                )
            }
        }

        combine_plots(plots, combine = combine, nrow = nrow, ncol = ncol, byrow = byrow)
    } else if (plot_type == "venn") {
        if (length(groups) > 4) {
            stop("[ClonalResidencyPlot] Too many groups for venn plot. Please use 'upset' plot instead.")
        }
        data <- data[data[[group_by]] %in% groups, , drop = FALSE]
        data <- data[data$count > 0, , drop = FALSE]
        if (!is.null(split_by)) {
            # # data <- unite(data, ".split", split_by, sep = split_by_sep, remove = FALSE)
            split_by <- check_columns(
                data, split_by,
                force_factor = TRUE, allow_multi = TRUE, concat_multi = TRUE, concat_sep = split_by_sep
            )
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
                        indicator <- indicator & data[[split_by]] == df$.split[i]
                    }
                    if (with_class) {
                        ns <- data[indicator, , drop = FALSE]
                        label <- c(label, paste0(df$count[i], "\n(singlets: ", nrow(ns), ")"))
                    } else {
                        label <- c(label, df$count[i])
                    }
                }
            }
            label
        }

        VennDiagram(data,
            in_form = "long", id_by = "CloneID", group_by = group_by, label = label_fun,
            split_by = split_by, ...
        )
    } else if (plot_type == "upset") {
        # keep the order of group_by
        columns <- if (is.factor(data[[group_by]])) {
            levels(data[[group_by]])
        } else {
            unique(data[[group_by]])
        }
        columns <- paste0("count_", columns)

        data$fraction <- NULL
        data <- data[data[[group_by]] %in% groups, , drop = FALSE]
        data[[group_by]] <- factor(data[[group_by]], levels = rev(groups))
        data <- data[data$count > 0, , drop = FALSE]
        data <- data %>%
            pivot_wider(
                names_from = !!sym("group_by"),
                names_prefix = "count_",
                values_from = !!sym("count"),
                values_fill = 0
            )
        data <- data[, intersect(unique(c("CloneID", split_by, columns)), colnames(data)), drop = FALSE]

        if (with_class) {
            data <- data %>%
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
                dplyr::select(!starts_with("count_") | ends_with(" Singlet") | ends_with(" Expanded"))
        } else {
            data <- data %>%
                mutate(across(
                    .cols = starts_with("count_"),
                    .names = "{.col}",
                    .fns = ~ .x > 0
                ))
        }

        data <- data %>%
            rename_with(
                function(x) substring(x, 7),
                .cols = starts_with("count_")
            )

        UpsetPlot(data, in_form = "wide", id_by = "CloneID", split_by = split_by, split_by_sep = split_by_sep, ...)
    }
}

#' Clonal Composition Plot
#'
#' @description
#' Visualizes the composition of the immune repertoire by categorizing clones
#' into abundance groups (Rare, Small, Medium, Large, Hyperexpanded) and
#' plotting their relative proportions across samples or metadata groups.
#' This reveals the overall structure of the repertoire — whether it is
#' dominated by a few large clones (clonal expansion) or composed of many
#' small clones (high diversity).
#'
#' `ClonalCompositionPlot` supports three analysis methods:
#' \itemize{
#'   \item **Homeostasis** (`"homeostasis"`, `"homeo"`, `"rel"`) —
#'         Clones are binned by their frequency (fraction of the total
#'         repertoire) into categories such as Rare, Small, Medium, Large,
#'         and Hyperexpanded. Uses
#'         \code{\link[scRepertoire:clonalHomeostasis]{scRepertoire::clonalHomeostasis()}}.
#'   \item **Top clones** (`"top"`) — Clones are ranked and binned by
#'         their rank index (e.g., top 10, top 100, etc.). Uses
#'         \code{\link[scRepertoire:clonalProportion]{scRepertoire::clonalProportion()}}.
#'   \item **Rare clones** (`"rare"`) — Clones are binned by their
#'         absolute size (clone count). Uses clone size thresholds directly.
#' }
#'
#' @param data The product of \code{\link[scRepertoire:combineTCR]{scRepertoire::combineTCR()}},
#'   \code{\link[scRepertoire:combineBCR]{scRepertoire::combineBCR()}}, or
#'   \code{\link[scRepertoire:combineExpression]{scRepertoire::combineExpression()}}.
#' @param clone_call How to define a clone. One of `"gene"`, `"nt"`,
#'   `"aa"` (default), `"strict"`, or a custom variable name in the data.
#' @param chain Which chain(s) to use: `"both"` (default), `"TRA"`, `"TRB"`,
#'   `"TRD"`, `"TRG"`, `"IGH"`, or `"IGL"`.
#' @param method The clonal categorization method. One of:
#'   \itemize{
#'     \item `"homeostasis"` (default), `"homeo"`, `"rel"` — Frequency-based
#'           binning using `clone_split` as abundance thresholds.
#'     \item `"top"` — Rank-based binning using `clone_split` as rank
#'           cutoffs.
#'     \item `"rare"` — Size-based binning using `clone_split` as clone
#'           size thresholds.
#'   }
#' @param clone_split Threshold values defining the clonal categories.
#'   Default is `NULL`, which picks sensible defaults per method:
#'   \itemize{
#'     \item For `"homeostasis"`/`"homeo"`/`"rel"` — a named list of
#'           frequency thresholds:
#'           `list(Rare = 1e-04, Small = 0.001, Medium = 0.01, Large = 0.1, Hyperexpanded = 1)`
#'     \item For `"top"` — rank cutoffs:
#'           `c(10, 100, 1000, 10000, 30000, 100000)`
#'     \item For `"rare"` — clone size thresholds:
#'           `c(1, 3, 10, 30, 100)`
#'   }
#' @param scale How to normalize the values. One of:
#'   \itemize{
#'     \item `TRUE` (default) — Values within each x-axis group sum to 1
#'           (group-wise proportion).
#'     \item `FALSE` — Raw values (clone counts) are used.
#'     \item `"sample"` or `"Sample"` — Values within each sample sum to
#'           1 (sample-wise proportion).
#'   }
#' @param facet_by Metadata column used to facet the plot into separate
#'   panels. Default is `NULL`.
#' @param group_by Metadata column used to group (color) the data. Default
#'   is `NULL`. Required for `"box"` and `"violin"` plot types.
#' @param split_by Metadata column used to split the data into separate
#'   plots. Default is `NULL`.
#' @param xlab X-axis label. Default is `NULL`, which uses the `group_by`
#'   column name or `"Sample"`.
#' @param ylab Y-axis label. Default is `NULL`, which auto-generates
#'   `"Abundance"` or `"Relative Abundance"` depending on `scale`.
#' @param plot_type The visualization type. One of:
#'   \itemize{
#'     \item `"bar"` (default) — Stacked bar chart of clonal categories
#'           across groups. Best for comparing composition across
#'           categories.
#'     \item `"ring"` — Ring (donut) chart alternative to stacked bars.
#'     \item `"box"` — Box plot showing the distribution of each clonal
#'           category's abundance across samples. Requires `group_by`.
#'     \item `"violin"` — Violin plot alternative to box plot. Requires
#'           `group_by`.
#'   }
#' @param order A named list controlling the order of factor levels. List
#'   names are column names; list values are the desired order. Default is
#'   `NULL`.
#' @param ... Additional arguments passed to the underlying \pkg{plotthis}
#'   function:
#'   \itemize{
#'     \item `"bar"` — \code{\link[plotthis:BarPlot]{plotthis::BarPlot()}}
#'           (`palette`, `position`, `alpha`, ...)
#'     \item `"ring"` — \code{\link[plotthis:RingPlot]{plotthis::RingPlot()}}
#'           (`palette`, `alpha`, ...)
#'     \item `"box"` — \code{\link[plotthis:BoxPlot]{plotthis::BoxPlot()}}
#'           (`comparisons`, `add_bg`, `palette`, ...)
#'     \item `"violin"` — \code{\link[plotthis:ViolinPlot]{plotthis::ViolinPlot()}}
#'           (`add_box`, `add_bg`, `comparisons`, `palette`, ...)
#'   }
#' @return A `ggplot` object, or a list of `ggplot` objects if
#'   `combine = FALSE` is passed via `...`.
#' @note
#' **group_by for box/violin:** The `group_by` parameter is required when
#' `plot_type` is `"box"` or `"violin"`. These plot types show per-sample
#' distributions, with `group_by` determining the coloring.
#'
#' **Bar/ring aggregation:** When `group_by` is specified for bar or ring
#' plots, data is aggregated across samples within each group (Sample
#' values are summed) before plotting, to show group-level composition.
#'
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
#'     variables = factor(rep(c("B", "L"), 4), levels = c("L", "B"))
#' )
#' data <- scRepertoire::addVariable(data,
#'     variable.name = "Subject",
#'     variables = rep(c("P17", "P18", "P19", "P20"), each = 2)
#' )
#'
#' ClonalCompositionPlot(data)
#' ClonalCompositionPlot(data, method = "top")
#' ClonalCompositionPlot(data, plot_type = "ring")
#' ClonalCompositionPlot(data, group_by = "Type", plot_type = "box", comparison = TRUE,
#'  clone_split = list(Small = 0.001, Medium = 0.01, Large = 0.1, Hyperexpanded = 1))
#' ClonalCompositionPlot(data, group_by = "Type", plot_type = "violin", add_box = TRUE,
#'  add_bg = TRUE)
#' ClonalCompositionPlot(data, method = "rare")
#' }
ClonalCompositionPlot <- function(
    data, clone_call = "aa", chain = "both", method = c("homeostasis", "homeo", "rel", "top", "rare"),
    clone_split = NULL, scale = TRUE, facet_by = NULL, group_by = NULL, split_by = NULL,
    xlab = NULL, ylab = NULL, plot_type = c("bar", "ring", "box", "violin"), order = NULL, ...
) {
    plot_type <- match.arg(plot_type)
    method <- match.arg(method)
    if (plot_type %in% c("box", "violin") && is.null(group_by)) {
        stop("[ClonalCompositionPlot] 'group_by' must be provided for box/violin ClonalCompositionPlot")
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
        grouping_levels <- get_clonal_grouping_levels(data, all_groupings, order)
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

        for (gl in names(grouping_levels)) {
            if (!is.null(data[[gl]])) {
                data[[gl]] <- factor(data[[gl]], levels = grouping_levels[[gl]])
            }
        }
        # Sample Type  .names        .values
        name_levels <- unique(data$.names)
    } else {  # rare
        data <- clonal_size_data(data, clone_call, chain, all_groupings, order)

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

#' Clonal Overlap Plot
#'
#' @description
#' Visualizes the overlap (sharing) of T-cell or B-cell clonotypes between
#' samples or metadata groups as a heatmap. Each cell in the heatmap
#' quantifies the degree of clonal sharing between two groups, using one of
#' several similarity or overlap metrics. This is a key analysis for
#' identifying public clones shared across individuals, tracking
#' antigen-specific clones across time points or tissues, and comparing
#' repertoire similarity between conditions.
#'
#' `ClonalOverlapPlot` computes pairwise overlap via
#' \code{\link[scRepertoire:clonalOverlap]{scRepertoire::clonalOverlap()}}
#' and visualizes the resulting matrix as a labeled heatmap using
#' \code{\link[plotthis:Heatmap]{plotthis::Heatmap()}}.
#'
#' @param data The product of \code{\link[scRepertoire:combineTCR]{scRepertoire::combineTCR()}},
#'   \code{\link[scRepertoire:combineBCR]{scRepertoire::combineBCR()}}, or
#'   \code{\link[scRepertoire:combineExpression]{scRepertoire::combineExpression()}}.
#' @param clone_call How to define a clone. One of `"gene"`, `"nt"`,
#'   `"aa"` (default), `"strict"`, or a custom variable name in the data.
#' @param chain Which chain(s) to use: `"both"` (default), `"TRA"`, `"TRB"`,
#'   `"TRD"`, `"TRG"`, `"IGH"`, or `"IGL"`.
#' @param group_by Metadata column(s) used to define the groups being
#'   compared. Default is `"Sample"`. Multiple columns are concatenated
#'   using `group_by_sep` to form compound group labels.
#' @param group_by_sep Separator used when concatenating multiple `group_by`
#'   columns. Default is `"_"`.
#' @param full Logical; if `TRUE` (default), the full symmetric heatmap is
#'   displayed. If `FALSE`, only the upper triangle is shown (lower triangle
#'   values are mirrored from the upper triangle).
#' @param split_by Metadata column used to split the data into separate
#'   heatmaps. When specified, overlaps are only calculated within each
#'   split group (not across splits). Default is `NULL`.
#' @param order A named list controlling the order of factor levels. List
#'   names are column names; list values are the desired order. Default is
#'   `NULL`.
#' @param method The overlap or similarity metric. One of:
#'   \itemize{
#'     \item `"raw"` (default) — Absolute number of overlapping clones
#'           between two groups.
#'     \item `"overlap"` — Overlap coefficient: size of intersection
#'           divided by the size of the smaller set.
#'     \item `"morisita"` — Morisita’s overlap index, accounting for
#'           clone size (abundance) in addition to presence/absence.
#'     \item `"jaccard"` — Jaccard similarity index: size of intersection
#'           divided by size of union.
#'     \item `"cosine"` — Cosine similarity between clone abundance
#'           vectors.
#'   }
#' @param palette Color palette for the heatmap. Default is `"Blues"`.
#' @param label_accuracy Numeric; the number of decimal places shown in
#'   cell labels. Default is `NULL`, which uses `1` for `"raw"` and
#'   `0.01` for other methods.
#' @param label_cutoff Numeric; values below this threshold are not labeled
#'   in the heatmap cells. Default is `1e-3`. Set to `0` to show all
#'   labels.
#' @param cluster_rows Logical; if `TRUE`, rows are hierarchically
#'   clustered. Default is `FALSE`.
#' @param cluster_columns Logical; if `TRUE`, columns are hierarchically
#'   clustered. Default is `FALSE`. Clustering distance is computed as
#'   `1 - rescaled_overlap_value`.
#' @param show_row_names Logical; if `TRUE` (default), row names are
#'   displayed.
#' @param show_column_names Logical; if `TRUE` (default), column names are
#'   displayed.
#' @param ... Additional arguments passed to
#'   \code{\link[plotthis:Heatmap]{plotthis::Heatmap()}} (e.g., `name`,
#'   `cell_type`, `width`, `height`).
#' @return A `ComplexHeatmap` object, or a list if `combine = FALSE` is
#'   passed via `...`.
#' @note
#' **Clustering:** When `cluster_rows` or `cluster_columns` is `TRUE`, the
#' clustering distance is `1 - rescaled(values)` for the upper triangle,
#' ensuring that groups with high overlap are placed close together.
#'
#' **Split behavior:** When `split_by` is specified, overlap is calculated
#' independently within each split group. Groups from different splits are
#' never compared against each other.
#'
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
#'     variables = factor(rep(c("B", "L"), 4), levels = c("L", "B"))
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
    split_by = NULL, order = NULL, method = c("raw", "overlap", "morisita", "jaccard", "cosine"),
    palette = "Blues", label_accuracy = NULL, label_cutoff = 1e-3, cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = TRUE, show_column_names = TRUE, ...
) {
    method <- match.arg(method)

    all_groupings <- unique(c(group_by, split_by))
    grouping_levels <- get_clonal_grouping_levels(data, all_groupings, order)
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

    for (gl in names(grouping_levels)) {
        if (!is.null(data[[gl]])) {
            data[[gl]] <- factor(data[[gl]], levels = grouping_levels[[gl]])
        }
    }

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
