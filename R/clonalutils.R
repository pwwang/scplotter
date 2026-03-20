#' Get the grouping levels for clonal data for later to restore
#'
#' @description Data transform with multiple groupings may lose the original grouping levels.
#' This function is to get the grouping levels for later to restore.
#' @param data The product of [scRepertoire::combineTCR], [scRepertoire::combineTCR], or [scRepertoire::combineExpression].
#' @param groupings The column names in the meta data to group the cells.
#' @param order A list specifying the order of the levels for each grouping variable.
#' @return A list of the grouping levels for each grouping variable.
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

#' clonal_size_data
#'
#' @description Function to get the clonal size data for all group_by values.
#'
#' @param data The product of [scRepertoire::combineTCR], [scRepertoire::combineTCR], or
#'  [scRepertoire::combineExpression].
#' @param clone_call How to call the clone - VDJC gene (gene), CDR3 nucleotide (nt),
#'  CDR3 amino acid (aa), VDJC gene + CDR3 nucleotide (strict) or a custom variable
#'  in the data
#' @param chain indicate if both or a specific chain should be used - e.g. "both",
#'  "TRA", "TRG", "IGH", "IGL"
#' @param groupings The column names in the meta data to group the cells.
#' @param order A list specifying the order of the levels for each grouping variable. Default is NULL, which will use the order in the data.
#' @return A data frame with the clonal size data.
#' @keywords internal
#' @importFrom rlang sym
#' @importFrom dplyr distinct filter rename
#' @importFrom tidyr pivot_longer
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

#' merge_clonal_groupings
#'
#' @description Merge the multiple clonal groupings into a single grouping.
#' @details Because [scRepertoire::clonalQuant] and families don't support mutliple groupings,
#'  this is trying to merge the multiple groupings into a single grouping. And then
#'  later restore the original groupings.
#' @param data The product of [scRepertoire::combineTCR], [scRepertoire::combineTCR], or
#'  [scRepertoire::combineExpression].
#' @param groupings A list of the clonal groupings. Each element is a column in the data.
#' @return The data with the combined groupings (`.group`)
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

#' Subset scRepertorie object
#'
#' @param screp The scRepertorie object. It is either a Seurat object or a list of data.frames
#' @param subset The subset expression (in characters)
#' @return The subsetted scRepertorie object
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
