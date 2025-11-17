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

#' ClonalDiversityPlot
#'
#' Plot the clonal diversities of the samples/groups.
#'
#' @param data The product of [scRepertoire::combineTCR], [scRepertoire::combineTCR], or
#'  [scRepertoire::combineExpression].
#' @param clone_call How to call the clone - VDJC gene (gene), CDR3 nucleotide (nt),
#'  CDR3 amino acid (aa), VDJC gene + CDR3 nucleotide (strict) or a custom variable
#'  in the data
#' @param chain indicate if both or a specific chain should be used - e.g. "both",
#'  "TRA", "TRG", "IGH", "IGL"
#' @param method The method to calculate the diversity. Options are "shannon" (default),
#'  "inv.simpson", "norm.entropy", "gini.simpson", "chao1", "ACE", "gini.coeff", "d50" and "dXX".
#'  See [scRepertoire::clonalDiversity] for details.
#'  The last 3 methods are supported by `scplotter` only:
#'  * "gini.coeff" - The Gini Coefficient. A measure of inequality in the distribution of clones.
#'    0 indicates perfect equality, 1 indicates perfect inequality.
#'  * "d50" - The number of clones that make up `50%` of the total number of clones.
#'  * "dXX" - The number of clones that make up `XX%` of the total number of clones.
#' @param d The percentage for the "dXX" method. Default is 50.
#' @param plot_type The type of plot. Options are "bar", "box" and "violin".
#' @param position The position adjustment for the bars. Default is "dodge".
#' @param group_by A character vector of column names to group the samples. Default is NULL.
#' @param facet_by A character vector of column names to facet the plots. Default is NULL.
#' @param split_by A character vector of column names to split the plots. Default is NULL.
#' @param xlab The x-axis label. Default is NULL.
#' @param ylab The y-axis label. Default is NULL.
#' @param ... Other arguments passed to the specific plot function.
#'  * For "bar", [plotthis::BarPlot()].
#'  * For "box", [plotthis::BoxPlot()].
#'  * For "violin", [plotthis::ViolinPlot()].
#' @return A ggplot object or a list if `combine` is FALSE
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
#'     variables = rep(c("B", "L"), 4)
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
    d = 50, plot_type = c("bar", "box", "violin"), position = "dodge",
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

    data <- merge_clonal_groupings(data, all_groupings)
    data <- ClonalDiversity(data,
        cloneCall = clone_call, chain = chain, method = method, d = d, group_by = ".group"
    )
    data <- separate(data, ".group", into = all_groupings, sep = " // ")

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

#' ClonalRarefactionPlot
#'
#' Plot the rarefaction curves
#'
#' @param data The product of [scRepertoire::combineTCR], [scRepertoire::combineTCR], or
#'  [scRepertoire::combineExpression].
#' @param clone_call How to call the clone - VDJC gene (gene), CDR3 nucleotide (nt),
#'  CDR3 amino acid (aa), VDJC gene + CDR3 nucleotide (strict) or a custom variable
#' @param chain indicate if both or a specific chain should be used - e.g. "both",
#'  "TRA", "TRG", "IGH", "IGL"
#' @param group_by A character vector of column names to group the samples. Default is "Sample".
#' @param group_by_sep The separator for the group_by column. Default is "_".
#' @param n_boots The number of bootstrap samples. Default is 20.
#' @param q The hill number. Default is 0.
#'  * 0 - Species richness
#'  * 1 - Shannon entropy
#'  * 2 - Simpson index#'
#' @param facet_by A character vector of column names to facet the plots. Default is NULL.
#' @param split_by A character vector of column names to split the plots. Default is NULL.
#' @param split_by_sep The separator for the split_by column. Default is "_".
#' @param palette The color palette to use. Default is "Paired".
#' @param combine Whether to combine the plots into a single plot. Default is TRUE.
#' @param nrow The number of rows in the combined plot. Default is NULL.
#' @param ncol The number of columns in the combined plot. Default is NULL.
#' @param byrow Whether to fill the combined plot by row. Default is TRUE.
#' @param ... Other arguments passed to [plotthis::RarefactionPlot()].
#' @return A ggplot object or a list if `combine` is FALSE
#' @importFrom plotthis RarefactionPlot
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
#' ClonalRarefactionPlot(data, type = 1, q = 0, n_boots = 2)
#' ClonalRarefactionPlot(data, type = 2, q = 0, n_boots = 2)
#' ClonalRarefactionPlot(data, type = 3, q = 0, n_boots = 2)
#' ClonalRarefactionPlot(data, q = 1, n_boots = 2)
#' ClonalRarefactionPlot(data, q = 1, n_boots = 2, group_by = "Type")
#' ClonalRarefactionPlot(data, n_boots = 2, split_by = "Type")
#' }
ClonalRarefactionPlot <- function(
    data, clone_call = "aa", chain = "both",
    group_by = "Sample", group_by_sep = "_",
    n_boots = 20, q = 0, facet_by = NULL, split_by = NULL, split_by_sep = "_",
    palette = "Paired", combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, ...
) {
    if (!is.null(facet_by)) {
        stop("'facet_by' is not supported in ClonalRarefactionPlot.")
    }

    all_groupings <- unique(c(group_by, split_by))
    data <- merge_clonal_groupings(data, all_groupings)

    .data.wrangle <- getFromNamespace(".data.wrangle", "scRepertoire")
    .theCall <- getFromNamespace(".theCall", "scRepertoire")
    .groupList <- getFromNamespace(".groupList", "scRepertoire")
    is_se_object <- getFromNamespace("is_se_object", "scRepertoire")
    is_seurat_object <- getFromNamespace("is_seurat_object", "scRepertoire")

    data <- .data.wrangle(data, ".group", .theCall(data, clone_call, check.df = FALSE),
        chain)
    cloneCall <- .theCall(data, clone_call)

    if (!is_seurat_object(data) && !is_se_object(data)) {
        data <- .groupList(data, group.by = ".group")
    }

    if (is.null(split_by)) {
        matlist <- lapply(data, function(x) { table(x[, cloneCall]) })
        RarefactionPlot(matlist, q = q, nboot = n_boots, palette = palette,
            group_name = paste(group_by, sep = group_by_sep), ...)
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
                group_name = paste(group_by, sep = group_by_sep), title = nm, ...)
        })

        combine_plots(plots, combine = combine, nrow = nrow, ncol = ncol, byrow = byrow)
    }
}
