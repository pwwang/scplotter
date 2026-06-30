#' Calculate the clonal diversities.
#'
#' @keywords internal
#' @importFrom utils getFromNamespace
#' @importFrom dplyr slice_sample
ClonalDiversity <- function(
    input.data, cloneCall = "gene", chain = "both",
    method = c("shannon", "gini.coeff", "inv.simpson", "norm.entropy", "gini.simpson", "chao1", "ACE", "d50", "dXX"),
    d = 50, group_by = NULL, n_boots = 0) {
    method <- match.arg(method)
    if (method == "gini.coeff") {
        div_fn <- function(dat) {
            dat <- sort(dat)
            n <- length(dat)
            1 / n * (n + 1 - 2 * sum((n + 1 - 1:n) * dat) / sum(dat))
        }
    } else if (method == "d50" || method == "dXX") {
        if (method == "d50") {
            d <- 50
        }
        div_fn <- function(dat) {
            dat <- sort(dat, decreasing = TRUE)
            dat <- 100 * cumsum(dat) / sum(dat)
            which(dat > d)[1]
        }
    } else {
        div_fn <- getFromNamespace(paste0(".", gsub(".", "", method, fixed = TRUE)), "scRepertoire")
    }

    is_seurat_object <- getFromNamespace("is_seurat_object", "scRepertoire")
    is_se_object <- getFromNamespace("is_se_object", "scRepertoire")
    .data.wrangle <- getFromNamespace(".data.wrangle", "scRepertoire")
    .theCall <- getFromNamespace(".theCall", "scRepertoire")
    .groupList <- getFromNamespace(".groupList", "scRepertoire")
    .short.check <- getFromNamespace(".short.check", "scRepertoire")

    sco <- is_seurat_object(input.data) | is_se_object(input.data)
    input.data <- .data.wrangle(input.data, group_by, .theCall(input.data,
        cloneCall,
        check.df = FALSE
    ), chain)
    cloneCall <- .theCall(input.data, cloneCall)
    mat <- NULL
    sample <- c()
    if (!is.null(group_by) & !sco) {
        input.data <- .groupList(input.data, group_by)
    }
    min <- .short.check(input.data, cloneCall)
    for (i in seq_along(input.data)) {
        data <- as.data.frame(table(input.data[[i]][, cloneCall]))
        mat_a <- NULL
        sample <- c()
        if (n_boots == 0) {
            sample <- div_fn(data$Freq)
            mat_a <- rbind(mat_a, sample)
            mat_a[is.na(mat_a)] <- 0
            mat <- rbind(mat, mat_a)
            mat <- as.data.frame(mat)
        } else {
            for (j in seq(seq_len(n_boots))) {
                x <- slice_sample(data, n = min)
                sample <- div_fn(x$Freq)
                mat_a <- rbind(mat_a, sample)
            }
            mat_a[is.na(mat_a)] <- 0
            mat_b <- colMeans(mat_a)
            mat_b <- as.data.frame(t(mat_b))
            mat <- rbind(mat, mat_b)
        }
    }
    if (is.null(group_by)) {
        group_by <- "Group"
    }
    colnames(mat) <- method
    mat[, group_by] <- names(input.data)

    mat
}

#' Clonal Diversity Plot
#'
#' @description
#' Visualizes clonal diversity metrics across samples or metadata groups.
#' Clonal diversity quantifies the richness and evenness of the immune
#' repertoire — how many distinct clonotypes are present and how evenly
#' cells are distributed among them. High diversity indicates a broad,
#' well-distributed repertoire; low diversity may indicate clonal expansion
#' (oligoclonality) in response to antigen stimulation or disease.
#'
#' `ClonalDiversityPlot` computes diversity scores using a custom
#' implementation that wraps several \pkg{scRepertoire} methods and adds
#' three \pkg{scplotter}-specific metrics (Gini coefficient, D50, DXX).
#' Results are visualized as bar, box, or violin plots.
#'
#' @section Diversity metrics:
#' The `method` parameter selects the diversity metric:
#'
#' **Richness and evenness metrics:**
#' \itemize{
#'   \item `"shannon"` (default) — Shannon entropy index. Higher values
#'         indicate greater diversity. Sensitive to both richness and
#'         evenness.
#'   \item `"inv.simpson"` — Inverse Simpson index. The effective number
#'         of equally abundant clones. Less sensitive to rare clones than
#'         Shannon.
#'   \item `"norm.entropy"` — Normalized entropy (Pielou's evenness).
#'         Shannon entropy divided by the log of richness; ranges from 0
#'         to 1.
#'   \item `"gini.simpson"` — Gini-Simpson index. The probability that
#'         two randomly selected cells belong to different clones.
#' }
#'
#' **Richness estimators (account for unobserved clones):**
#' \itemize{
#'   \item `"chao1"` — Chao1 richness estimator. Estimates the total
#'         number of clones including those not yet observed, based on
#'         the number of singletons and doubletons.
#'   \item `"ACE"` — Abundance-based Coverage Estimator. Estimates
#'         richness with a correction for sample coverage.
#' }
#'
#' **\pkg{scplotter}-specific metrics:**
#' \itemize{
#'   \item `"gini.coeff"` — Gini coefficient. Measures inequality in
#'         clone size distribution. `0` indicates perfect equality (all
#'         clones the same size); `1` indicates perfect inequality (one
#'         clone dominates).
#'   \item `"d50"` — The number of top clones that together account for
#'         50% of the total repertoire.
#'   \item `"dXX"` — The number of top clones that together account for
#'         `XX`% of the total repertoire. Use the `d` parameter to set
#'         the percentage.
#' }
#'
#' @param data The product of \code{\link[scRepertoire:combineTCR]{scRepertoire::combineTCR()}},
#'   \code{\link[scRepertoire:combineBCR]{scRepertoire::combineBCR()}}, or
#'   \code{\link[scRepertoire:combineExpression]{scRepertoire::combineExpression()}}.
#' @param clone_call How to define a clone. One of `"gene"` (default),
#'   `"nt"`, `"aa"`, `"strict"`, or a custom variable name in the data.
#' @param chain Which chain(s) to use: `"both"` (default), `"TRA"`, `"TRB"`,
#'   `"TRD"`, `"TRG"`, `"IGH"`, or `"IGL"`.
#' @param method The diversity metric to compute. One of `"shannon"`
#'   (default), `"inv.simpson"`, `"norm.entropy"`, `"gini.simpson"`,
#'   `"chao1"`, `"ACE"`, `"gini.coeff"`, `"d50"`, or `"dXX"`. See the
#'   \emph{Diversity metrics} section for details on each metric.
#' @param d The percentage threshold for the `"dXX"` method. For example,
#'   `d = 90` computes the number of clones accounting for 90% of the
#'   repertoire. Default is `50`.
#' @param plot_type The visualization type. One of `"bar"` (default),
#'   `"box"`, or `"violin"`. For `"box"` and `"violin"`, `group_by` is
#'   required to provide the x-axis grouping.
#' @param position Bar position adjustment for `"bar"` plot type. One of
#'   `"dodge"` (default), `"stack"`, or `"fill"`.
#' @param order A named list controlling the order of factor levels. List
#'   names are column names; list values are the desired order. Default is
#'   `NULL`.
#' @param group_by Metadata column used to group (color) the data. Default
#'   is `NULL`. Required for `"box"` and `"violin"` plot types.
#' @param facet_by Metadata column used to facet the plot into separate
#'   panels. Default is `NULL`.
#' @param split_by Metadata column used to split the data into separate
#'   plots. Default is `NULL`.
#' @param xlab X-axis label. Default is `NULL`, which uses the `group_by`
#'   column name or `"Sample"`.
#' @param ylab Y-axis label. Default is `NULL`, which auto-generates the
#'   full metric name (e.g., `"Shannon Index"`, `"Gini Coefficient"`).
#' @param ... Additional arguments passed to the underlying \pkg{plotthis}
#'   function:
#'   \itemize{
#'     \item `"bar"` — \code{\link[plotthis:BarPlot]{plotthis::BarPlot()}}
#'           (`palette`, `alpha`, `fill_by`, ...)
#'     \item `"box"` — \code{\link[plotthis:BoxPlot]{plotthis::BoxPlot()}}
#'           (`comparisons`, `alpha`, `palette`, ...)
#'     \item `"violin"` — \code{\link[plotthis:ViolinPlot]{plotthis::ViolinPlot()}}
#'           (`add_box`, `comparisons`, `palette`, ...)
#'   }
#' @return A `ggplot` object, or a list of `ggplot` objects if
#'   `combine = FALSE` is passed via `...`.
#' @note
#' **Bootstrap support:** The underlying `ClonalDiversity()` function
#' supports bootstrap resampling (`n_boots`). This is not exposed in
#' `ClonalDiversityPlot` directly but is used internally.
#'
#' **group_by required for box/violin:** The `group_by` parameter is
#' required when `plot_type` is `"box"` or `"violin"`. These types show
#' per-sample distributions grouped by the `group_by` variable.
#'
#' @export
#' @importFrom tidyr separate
#' @importFrom scRepertoire clonalDiversity
#' @importFrom plotthis BarPlot BoxPlot ViolinPlot
#' @examples
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
#' ClonalDiversityPlot(data)
#' ClonalDiversityPlot(data, group_by = "Type")
#' ClonalDiversityPlot(data, group_by = "Type", plot_type = "box")
#' ClonalDiversityPlot(data, group_by = "Type", plot_type = "violin")
#' ClonalDiversityPlot(data, group_by = "Type", plot_type = "violin",
#'   method = "gini.coeff", add_box = TRUE)
#' ClonalDiversityPlot(data, group_by = "Type", plot_type = "violin",
#'   method = "inv.simpson", add_box = TRUE)
#' ClonalDiversityPlot(data, group_by = "Type", plot_type = "violin",
#'   method = "d50", add_box = TRUE)
ClonalDiversityPlot <- function(
    data, clone_call = "gene", chain = "both",
    method = c("shannon", "gini.coeff", "inv.simpson", "norm.entropy", "gini.simpson", "chao1", "ACE", "d50", "dXX"),
    d = 50, plot_type = c("bar", "box", "violin"), position = "dodge", order = NULL,
    group_by = NULL, facet_by = NULL, split_by = NULL,
    xlab = NULL, ylab = NULL,
    ...) {
    method <- match.arg(method)
    plot_type <- match.arg(plot_type)

    if (plot_type %in% c("box", "violin") && is.null(group_by)) {
        stop("'group_by' must be provided for box/violin ClonalDiversityPlot")
    }
    all_groupings <- unique(c("Sample", group_by, facet_by, split_by))
    method_name <- switch(method,
        shannon = "Shannon Index",
        gini.coeff = "Gini Coefficient",
        inv.simpson = "Inverse Simpson Index",
        norm.entropy = "Normalized Entropy",
        gini.simpson = "Gini-Simpson Index",
        chao1 = "Chao1 Index",
        ACE = "ACE Index"
    )

    grouping_levels <- get_clonal_grouping_levels(data, all_groupings, order)
    data <- merge_clonal_groupings(data, all_groupings)
    data <- ClonalDiversity(data,
        cloneCall = clone_call, chain = chain, method = method, d = d, group_by = ".group"
    )
    data <- separate(data, ".group", into = all_groupings, sep = " // ")
    for (gl in names(grouping_levels)) {
        if (!is.null(data[[gl]])) {
            data[[gl]] <- factor(data[[gl]], levels = grouping_levels[[gl]])
        }
    }

    if (plot_type == "bar") {
        x <- group_by %||% "Sample"
        group_by <- if(is.null(group_by)) NULL else "Sample"
        BarPlot(data,
            x = x, y = method, group_by = group_by, position = position,
            xlab = xlab %||% group_by, ylab = ylab %||% method_name,
            split_by = split_by, facet_by = facet_by, ...
        )
    } else if (plot_type == "box") {
        BoxPlot(data,
            x = group_by, y = method,
            xlab = xlab %||% group_by, ylab = ylab %||% method_name,
            split_by = split_by, facet_by = facet_by, ...
        )
    } else if (plot_type == "violin") {
        ViolinPlot(data,
            x = group_by, y = method,
            xlab = xlab %||% group_by, ylab = ylab %||% method_name,
            split_by = split_by, facet_by = facet_by, ...
        )
    }
}

#' Clonal Rarefaction Plot
#'
#' @description
#' Visualizes clonal rarefaction curves — estimates of clone richness as a
#' function of sampling depth. Rarefaction addresses a fundamental challenge
#' in immune repertoire analysis: the number of clones observed depends on
#' how many cells are sequenced. By repeatedly subsampling (bootstrapping)
#' the data at varying depths, rarefaction curves reveal whether the
#' repertoire has been sampled to saturation or whether additional
#' sequencing would uncover many more clones.
#'
#' `ClonalRarefactionPlot` extracts clone count data from the repertoire,
#' optionally groups it by metadata columns, and generates rarefaction
#' curves via \code{\link[plotthis:RarefactionPlot]{plotthis::RarefactionPlot()}}.
#' When `split_by` is specified, separate plots are generated for each split
#' group and combined into a multi-panel layout.
#'
#' @section Hill numbers (the q parameter):
#' The `q` parameter selects the diversity order (Hill number) used for
#' rarefaction:
#' \itemize{
#'   \item **`q = 0`** — Species richness (clone count). Counts the number
#'         of distinct clonotypes regardless of their size. Most sensitive
#'         to rare clones.
#'   \item **`q = 1`** — Shannon entropy (exponential). Weighs clones
#'         proportionally to their abundance. Balances rare and dominant
#'         clones.
#'   \item **`q = 2`** — Simpson index (inverse). Weighs dominant clones
#'         more heavily. Least sensitive to rare clones.
#' }
#' Higher values of `q` increasingly emphasize abundant clones over rare
#' ones.
#'
#' @param data The product of \code{\link[scRepertoire:combineTCR]{scRepertoire::combineTCR()}},
#'   \code{\link[scRepertoire:combineBCR]{scRepertoire::combineBCR()}}, or
#'   \code{\link[scRepertoire:combineExpression]{scRepertoire::combineExpression()}}.
#' @param clone_call How to define a clone. One of `"gene"`, `"nt"`,
#'   `"aa"` (default), `"strict"`, or a custom variable name in the data.
#' @param chain Which chain(s) to use: `"both"` (default), `"TRA"`, `"TRB"`,
#'   `"TRD"`, `"TRG"`, `"IGH"`, or `"IGL"`.
#' @param group_by Metadata column(s) used to define the curves (each unique
#'   group produces one rarefaction curve). Multiple columns are concatenated
#'   using `group_by_sep`. Default is `"Sample"`.
#' @param group_by_sep Separator used when concatenating multiple `group_by`
#'   columns. Default is `"_"`.
#' @param order A named list controlling the order of factor levels. List
#'   names are column names; list values are the desired order. Default is
#'   `NULL`.
#' @param n_boots Number of bootstrap iterations for estimating confidence
#'   intervals. Higher values produce smoother confidence bands but increase
#'   computation time. Default is `20`.
#' @param q The diversity order (Hill number). `0` for species richness,
#'   `1` for Shannon entropy, `2` for Simpson index. Default is `0`. See
#'   the \emph{Hill numbers} section for details.
#' @param facet_by Not supported for `ClonalRarefactionPlot`. Use
#'   `split_by` or `group_by` instead. Must be `NULL`.
#' @param split_by Metadata column used to split the data into separate
#'   rarefaction plots. When specified, an independent rarefaction is
#'   performed for each split group, and all plots are combined. Default
#'   is `NULL`.
#' @param split_by_sep Separator used when concatenating multiple
#'   `split_by` columns. Default is `"_"`.
#' @param palette Color palette for distinguishing curves from different
#'   groups. Default is `"Paired"`.
#' @param combine Logical; if `TRUE` (default), multiple plots (from
#'   `split_by`) are combined into a single layout.
#' @param nrow Number of rows in the combined plot layout. Default is
#'   `NULL` (auto-determined).
#' @param ncol Number of columns in the combined plot layout. Default is
#'   `NULL` (auto-determined).
#' @param byrow Logical; if `TRUE` (default), the combined layout is filled
#'   row by row.
#' @param ... Additional arguments passed to
#'   \code{\link[plotthis:RarefactionPlot]{plotthis::RarefactionPlot()}}.
#'   Key parameters include:
#'   \itemize{
#'     \item `type` — Plot type: `1` (line only), `2` (line with
#'           confidence band), or `3` (confidence band only).
#'     \item `title` — Plot title.
#'     \item `xlab`, `ylab` — Axis labels.
#'   }
#' @return A `ggplot` object, or a list of `ggplot` objects if
#'   `combine = FALSE`.
#' @note
#' **Bootstrap iterations:** The `n_boots` parameter controls the number of
#' resampling iterations. Higher values give more stable estimates but
#' increase computation time linearly. For exploratory analysis, `n_boots =
#' 20` is typically sufficient; for publication-quality figures, consider
#' using `n_boots = 100` or more.
#'
#' **facet_by not supported:** Unlike many other \pkg{scplotter} functions,
#' `ClonalRarefactionPlot` does not support `facet_by`. Use `split_by` for
#' separate plots or `group_by` to show multiple curves on the same axes.
#'
#' @export
#' @importFrom plotthis RarefactionPlot
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
#' ClonalRarefactionPlot(data, type = 1, q = 0, n_boots = 2)
#' ClonalRarefactionPlot(data, type = 2, q = 0, n_boots = 2)
#' ClonalRarefactionPlot(data, type = 3, q = 0, n_boots = 2)
#' ClonalRarefactionPlot(data, q = 1, n_boots = 2)
#' ClonalRarefactionPlot(data, q = 1, n_boots = 2, group_by = "Type")
#' ClonalRarefactionPlot(data, n_boots = 2, split_by = "Type")
#' }
ClonalRarefactionPlot <- function(
    data, clone_call = "aa", chain = "both",
    group_by = "Sample", group_by_sep = "_", order = NULL,
    n_boots = 20, q = 0, facet_by = NULL, split_by = NULL, split_by_sep = "_",
    palette = "Paired", combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, ...
) {
    if (!is.null(facet_by)) {
        stop("'facet_by' is not supported in ClonalRarefactionPlot.")
    }

    all_groupings <- unique(c(group_by, split_by))
    grouping_levels <- get_clonal_grouping_levels(data, all_groupings, order)
    data <- merge_clonal_groupings(data, all_groupings)

    .data.wrangle <- getFromNamespace(".data.wrangle", "scRepertoire")
    .theCall <- getFromNamespace(".theCall", "scRepertoire")
    .groupList <- getFromNamespace(".groupList", "scRepertoire")
    is_se_object <- getFromNamespace("is_se_object", "scRepertoire")
    is_seurat_object <- getFromNamespace("is_seurat_object", "scRepertoire")

    data <- .data.wrangle(data, ".group", .theCall(data, clone_call, check.df = FALSE), chain)
    cloneCall <- .theCall(data, clone_call)

    if (!is_seurat_object(data) && !is_se_object(data)) {
        data <- .groupList(data, group.by = ".group")
        if (length(grouping_levels) == 1) {
            data <- data[order(match(names(data), grouping_levels[[1]]))]
        } else if (length(grouping_levels) > 1) {
            # all combinations of grouping levels, concat with ' // '
            combs <- expand.grid(grouping_levels, stringsAsFactors = FALSE)
            combs <- apply(combs, 1, paste, collapse = " // ")
            data <- data[match(combs, names(data))]
        }
    }

    if (is.null(split_by)) {
        matlist <- lapply(data, function(x) { table(x[, cloneCall]) })
        names(matlist) <- gsub(" // ", group_by_sep, names(data), fixed = TRUE)
        RarefactionPlot(matlist, q = q, nboot = n_boots, palette = palette,
            group_name = paste(group_by, collapse = group_by_sep), ...)
    } else {
        datas <- list()
        for (nm in names(data)) {
            nms <- strsplit(nm, " // ", fixed = TRUE)[[1]]
            names(nms) <- all_groupings
            gname <- paste(nms[group_by], collapse = group_by_sep)
            sname <- paste(nms[split_by], collapse = split_by_sep)
            d <- list(table(data[[nm]][, cloneCall]))
            names(d) <- gname
            datas[[sname]] <- c(datas[[sname]], d)
        }

        plots <- lapply(names(datas), function(nm) {
            RarefactionPlot(datas[[nm]], q = q, nboot = n_boots, palette = palette,
                group_name = paste(group_by, collapse = group_by_sep), title = nm, ...)
        })

        combine_plots(plots, combine = combine, nrow = nrow, ncol = ncol, byrow = byrow)
    }
}
