#' Extract grouping levels from clonal data for later restoration
#'
#' @description
#' When clonal data is transformed with multiple groupings (e.g., via
#' `merge_clonal_groupings`), the original factor levels of the grouping
#' variables can be lost. This function captures those levels before
#' transformation so they can be restored later, preserving the intended
#' ordering of groups in plots and analyses.
#'
#' The function handles multiple data formats: Seurat objects, single data
#' frames, and lists of data frames (the standard scRepertoire combined
#' TCR/BCR format). For list-format data, if all samples share the same
#' value for a grouping variable, the per-sample values are collected;
#' otherwise the levels from the first sample are used.
#' @param data The product of \code{\link[scRepertoire:combineTCR]{scRepertoire::combineTCR()}},
#'   \code{\link[scRepertoire:combineBCR]{scRepertoire::combineBCR()}}, or
#'   \code{\link[scRepertoire:combineExpression]{scRepertoire::combineExpression()}}.
#'   May also be a Seurat object or a single data frame.
#' @param groupings A character vector of column names in the metadata to group the cells by.
#' @param order A named list specifying the order of the levels for each grouping variable.
#'   Each element name should match a grouping column, and the element value should be a
#'   character vector of level names in the desired order. Default: \code{NULL} (use the
#'   order present in the data).
#' @return A named list where each element corresponds to a grouping variable and contains
#'   a character vector of its factor levels.
#' @keywords internal
get_clonal_grouping_levels <- function(data, groupings, order = NULL) {
    grouping_levels <- sapply(groupings, function(g) {
        if (!is.null(order) && !is.null(order[[g]])) {
            return(order[[g]])
        }
        dg <- if (inherits(data, "Seurat")) {
            data@meta.data[[g]]
        } else if (is.data.frame(data)) {
            data[[g]]
        } else if (is.factor(data[[1]][[g]])) {
            data[[1]][[g]]
        } else {
            if (all(sapply(data, function(d) length(unique(d[[g]])) == 1))) {
                # if all data frames have the same value for this grouping variable, just use the value in the first data frame
                sapply(data, function(d) d[[g]][1])
            } else {
                data[[1]][[g]]
            }
        }
        if (is.null(dg)) {
            return(NULL)
        }
        if (!is.factor(dg)) dg <- factor(dg)
        levels(dg)
    }, simplify = FALSE)

    grouping_levels[!sapply(grouping_levels, is.null)]
}

#' Prepare clonal abundance data for visualization
#'
#' @description
#' This is the central data preparation function used by all clonal visualization
#' functions in \pkg{scplotter}. It transforms raw scRepertoire combined TCR/BCR
#' data (a list of data frames or Seurat object) into a standardized long-format
#' data frame suitable for plotting.
#'
#' The function performs the following steps:
#' \enumerate{
#'   \item Extracts grouping levels via \code{get_clonal_grouping_levels()}.
#'   \item Merges multiple grouping columns into a single composite \code{.group}
#'     column via \code{merge_clonal_groupings()}.
#'   \item Reorganizes list-format data by splitting on \code{.group} and
#'     combining cells that share the same group across samples.
#'   \item Calls scRepertoire's internal \code{.data.wrangle()} and \code{.theCall()}
#'     to compute clone-level counts per group.
#'   \item Pivots the resulting wide-format table (clones × groups) to long format
#'     with columns \code{CloneID}, \code{.group}, \code{count}, and \code{fraction}.
#'   \item Separates the composite \code{.group} column back into the original
#'     grouping columns.
#'   \item Restores original factor levels for each grouping variable.
#' }
#'
#' The output data frame contains one row per clone per group, with columns for
#' each grouping variable, \code{CloneID}, \code{count} (number of cells), and
#' \code{fraction} (proportion of cells within the group).
#' @param data The product of \code{\link[scRepertoire:combineTCR]{scRepertoire::combineTCR()}},
#'   \code{\link[scRepertoire:combineBCR]{scRepertoire::combineBCR()}}, or
#'   \code{\link[scRepertoire:combineExpression]{scRepertoire::combineExpression()}}.
#'   A named list of data frames where each element represents a sample, or a
#'   Seurat object with TCR/BCR information in the metadata.
#' @param clone_call How to identify a clone. One of \code{"gene"} (VDJC gene
#'   segment), \code{"nt"} (CDR3 nucleotide sequence), \code{"aa"} (CDR3 amino
#'   acid sequence), \code{"strict"} (VDJC gene + CDR3 nucleotide), or a custom
#'   column name present in the data.
#' @param chain Which TCR/BCR chain(s) to include. One of \code{"both"} (both
#'   chains combined), \code{"TRA"}, \code{"TRB"}, \code{"TRG"}, \code{"TRD"},
#'   \code{"IGH"}, \code{"IGL"}, or \code{"IGK"}.
#' @param groupings A character vector of column names in the metadata to use
#'   for grouping cells. These define the categories across which clone
#'   abundances are computed.
#' @param order A named list specifying the order of the levels for each
#'   grouping variable. Default: \code{NULL} (uses the order in the data).
#' @return A data frame in long format with columns: the grouping variables,
#'   \code{CloneID} (clone identifier), \code{count} (number of cells per clone
#'   per group), and \code{fraction} (proportion of cells per clone within each
#'   group). Only rows with \code{count > 0} are retained.
#' @importFrom rlang sym
#' @importFrom dplyr distinct filter rename
#' @importFrom tidyr pivot_longer
#' @keywords internal
clonal_size_data <- function(data, clone_call, chain, groupings, order) {
    if (length(groupings > 0)) {
        grouping_levels <- get_clonal_grouping_levels(data, groupings, order)
        data <- merge_clonal_groupings(data, groupings)
    } else {
        if (inherits(data, "Seurat")) {
            data$.group <- "All"
        } else {
            data <- lapply(data, function(d) {
                d$.group <- "All"
                d
            })
        }

        groupings <- ".group"
        grouping_levels <- list()
    }

    if (inherits(data, "Seurat")) {
        # all_gvalues <- unique(data@meta.data$.group)
    } else {
        # clonalScatter only returns data for each sample
        # need to re-organize the data to get the data for each group
        # all_gvalues <- unique(unlist(sapply(data, function(x) x$.group)))
        newdata <- list()
        for (d in data) {
            nd <- split(d, d$.group)
            for (gn in names(nd)) {
                if (is.null(newdata[[gn]])) {
                    newdata[[gn]] <- nd[[gn]]
                } else {
                    newdata[[gn]] <- rbind(newdata[[gn]], nd[[gn]])
                }
            }
        }
        data <- newdata
    }

    .data.wrangle <- utils::getFromNamespace(".data.wrangle", "scRepertoire")
    .theCall <- utils::getFromNamespace(".theCall", "scRepertoire")

    data <- .data.wrangle(data, ".group", .theCall(data, clone_call, check.df = FALSE), chain)
    all_gvalues <- names(data)
    if (length(all_gvalues) < 2) {
        gv_pairs <- list(all_gvalues[1])
    } else {
        gv_pairs <- as.list(as.data.frame(combn(all_gvalues, 2, simplify = TRUE)))
    }
    clone_call <- .theCall(data, clone_call)

    df <- do_call(rbind, lapply(gv_pairs, function(gv) {
        if (length(gv) == 1) {
            mat <- as.data.frame(table(data[[gv[1]]][, clone_call]))
            x_col <- colnames(mat)[2] <- paste0(gv[1], ".count")
            mat[is.na(mat)] <- 0
            mat[, paste0(gv[1], ".fraction")] <- mat[, x_col] / sum(mat[, x_col])
        } else {
            x_df <- as.data.frame(table(data[[gv[1]]][, clone_call]))
            x_col <- colnames(x_df)[2] <- paste0(gv[1], ".count")
            y_df <- as.data.frame(table(data[[gv[2]]][, clone_call]))
            y_col <- colnames(y_df)[2] <- paste0(gv[2], ".count")

            mat <- merge(x_df, y_df, by = "Var1", all = TRUE)
            mat[is.na(mat)] <- 0
            mat[, paste0(gv[1], ".fraction")] <- mat[, x_col] / sum(mat[, x_col])
            mat[, paste0(gv[2], ".fraction")] <- mat[, y_col] / sum(mat[, y_col])
        }
        pivot_longer(
            mat,
            cols = -"Var1",
            names_to = c(".group", ".value"),
            names_pattern = "(.+)\\.(.+)"
        )
    })) %>%
        distinct(!!sym("Var1"), !!sym(".group"), .keep_all = TRUE) %>%
        separate(".group", into = groupings, sep = " // ") %>%
        filter(!!sym("count") > 0) %>%
        rename(CloneID = "Var1")

    for (gl in names(grouping_levels)) {
        if (!is.null(df[[gl]])) {
            df[[gl]] <- factor(df[[gl]], levels = grouping_levels[[gl]])
        }
    }

    df
}

#' Merge multiple grouping columns into a single composite grouping
#'
#' @description
#' Because \code{\link[scRepertoire:clonalQuant]{scRepertoire::clonalQuant()}} and related
#' functions do not support multiple grouping variables, this function merges all
#' grouping columns into a single composite \code{.group} column (using
#' \code{\link[tidyr:unite]{tidyr::unite()}} with a separator). The original groupings
#' can later be restored by splitting on the separator.
#'
#' The function handles both Seurat objects (operating on \code{@meta.data}) and
#' the standard scRepertoire list-of-data-frames format. When no groupings are
#' provided, a placeholder \code{.group} column set to an empty string is created,
#' effectively treating all cells as a single group.
#' @param data The product of \code{\link[scRepertoire:combineTCR]{scRepertoire::combineTCR()}},
#'   \code{\link[scRepertoire:combineBCR]{scRepertoire::combineBCR()}}, or
#'   \code{\link[scRepertoire:combineExpression]{scRepertoire::combineExpression()}}.
#'   May also be a Seurat object.
#' @param groupings A character vector of column names to merge into the composite
#'   \code{.group} column. If empty, \code{.group} is set to an empty string.
#' @param sep The separator string used to join grouping values. Default: \code{" // "}.
#' @return The input data with an added \code{.group} column containing the
#'   concatenated grouping values (or an empty string if no groupings are provided).
#'   For Seurat objects, the column is added to \code{@meta.data}. For list-format
#'   data, each data frame gains the column.
#' @importFrom rlang syms
#' @importFrom tidyr unite
#' @keywords internal
merge_clonal_groupings <- function(data, groupings, sep = " // ") {
    if (inherits(data, "Seurat")) {
        if (!"Sample" %in% colnames(data@meta.data)) {
            warning("The 'Sample' column is not found in the meta data, 'orig.indent' will be used instead.")
            data$Sample <- data$orig.ident
        }
        # samples <- unique(data$Sample)
        if (length(groupings) == 0) {
            data@meta.data$.group <- ""
        } else {
            # combine Sample and group_by so that the count/frequency is calculated for each
            # combined group
            data@meta.data <- unite(data@meta.data, ".group", !!!syms(groupings), sep = sep, remove = FALSE)
        }
    } else {
        samples <- names(data)
        data <- lapply(samples, function(s) {
            if (nrow(data[[s]]) == 0) {
                data[[s]]$Sample <- character()
            } else {
                data[[s]]$Sample <- s
            }
            data[[s]]
        })
        names(data) <- samples
        # combine Sample and group_by so that the count/frequency is calculated for each
        # combined group
        data <- lapply(data, function(d) {
            if (length(groupings) == 0) {
                d$.group <- ""
            } else {
                d <- unite(d, ".group", !!!syms(groupings), sep = sep, remove = TRUE)
            }
            d
        })
    }

    data
}

#' Subset an scRepertoire object with a filter expression
#'
#' @description
#' Filters an scRepertoire object (either a Seurat object or a named list of
#' data frames) using a character-specified filter expression. This is useful
#' for restricting analysis to a subset of cells (e.g., cells with short CDR3
#' sequences, or cells from specific samples) without manually manipulating the
#' underlying data structure.
#'
#' For list-format data (the standard scRepertoire combined TCR/BCR format), a
#' \code{Sample} column is temporarily added to each data frame to enable
#' sample-based filtering. Empty data frames after filtering are removed from
#' the result.
#' @param screp An scRepertoire object — either a Seurat object (requires
#'   \pkg{tidyseurat} for dplyr filter support) or a named list of data frames
#'   (the standard output of \code{\link[scRepertoire:combineTCR]{scRepertoire::combineTCR()}}
#'   or \code{\link[scRepertoire:combineBCR]{scRepertoire::combineBCR()}}).
#' @param subset A character string containing a filter expression (e.g.,
#'   \code{"nchar(CTaa) < 20"} or \code{"Sample \%in\% c('P17B', 'P17L')"}).
#'   The expression is parsed and evaluated within the data context using
#'   \code{\link[rlang:parse_expr]{rlang::parse_expr()}}.
#' @return The filtered scRepertoire object in the same format as the input.
#'   For list-format data, samples with zero rows after filtering are dropped.
#' @keywords internal
#' @importFrom rlang parse_expr
#' @importFrom dplyr filter
#' @examples
#' \donttest{
#' data(contig_list, package = "scRepertoire")
#' screp <- scRepertoire::combineTCR(contig_list,
#'     samples = c("P17B", "P17L", "P18B", "P18L", "P19B", "P19L", "P20B", "P20L")
#' )
#'
#' head(scplotter:::screp_subset(screp, "nchar(CTaa) < 20")[[1]])
#' names(scplotter:::screp_subset(screp, "Sample %in% c('P17B', 'P17L')"))
#' }
screp_subset <- function(screp, subset) {
    if (inherits(screp, "Seurat")) {
        # You need tidyseurat to work with dplyr verbs on Seurat objects
        if (is.null(utils::getS3method("filter", "Seurat", optional = TRUE, envir = asNamespace("dplyr")))) {
            stop("'tidyseurat' package is required to use 'screp_subset' on Seurat objects.")
        }
        dplyr::filter(screp, !!parse_expr(subset))
    } else {
        screp <- sapply(names(screp), function(x) {
            y <- screp[[x]]
            if (nrow(y) == 0) {
                y$Sample <- character()
            } else {
                y$Sample <- x
            }
            filter(y, !!parse_expr(subset))
        }, simplify = FALSE, USE.NAMES = TRUE)
        screp[unlist(sapply(screp, nrow)) > 0]
    }
}
