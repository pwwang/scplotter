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
#' @return A data frame with the clonal size data.
#' @keywords internal
#' @importFrom rlang sym
#' @importFrom dplyr distinct filter rename
#' @importFrom tidyr pivot_longer
clonal_size_data <- function(data, clone_call, chain, groupings) {
    if (length(groupings > 0)) {
        grouping_levels <- sapply(groupings, function(g) {
            dg <- if (inherits(data, "Seurat")) {
                data@meta.data[[g]]
            } else {
                data[[g]]
            }
            if (is.null(dg)) {
                return(NULL)
            }
            if (!is.factor(dg)) dg <- factor(dg)
            levels(dg)
        })
        grouping_levels <- grouping_levels[!sapply(grouping_levels, is.null)]

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

    df <- do.call(rbind, lapply(gv_pairs, function(gv) {
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

#' Utility functions for clone selectors
#'
#' @rdname clone_selector_utils
#' @inheritParams CloneSelectors
#' @description These utility functions are used internally by the clone selectors to handle
#' different environments and return the results in the desired format.
#' @importFrom rlang caller_env env_has env_get sym
#' @importFrom dplyr filter left_join mutate summarise arrange
#' @importFrom tidyr pivot_wider
#' @importFrom rlang parse_expr
#' @importFrom dplyr across
#' @param envtype The type of environment to use. It can be "tidy" for dplyr verbs, "scplotter" for scplotter functions, or "unknown" if the context cannot be determined.
#' @param out The output data frame to be processed.
#' @param x The name to be backtick-quoted
#' @keywords internal
.get_envtype <- function() {
    # Returns:
    # - "tidy" if the function is called in a tidyverse context (e.g. dplyr verbs)
    # - "scplotter" if the function is called in a scplotter function context
    # - "unknown" if the context cannot be determined
    env <- caller_env(n = 2)
    if (env_has(env, ".data") && env_has(env, ".__tidyeval_data_mask__.")) {
        "tidy"
    } else if (env_has(env, "data") && is.data.frame(env$data)) {
        "scplotter"
    } else {
        "unknown"
    }
}

#' @rdname clone_selector_utils
#' @keywords internal
.get_data <- function(envtype) {
    # Get the data from the environment if not provided
    if (envtype == "tidy") {
        return(across(everything()))
    } else if (envtype == "scplotter") {
        env <- caller_env(n = 2)
        if (env_has(env, "data") && is.data.frame(env$data)) {
            return(env_get(env, "data"))
        }
    }
    stop("No data not found in caller environment. Please provide it to the 'data' argument.")
}

#' @rdname clone_selector_utils
#' @keywords internal
.return_what <- function(out, id, output) {
    if (identical(output, "id") || identical(output, "ids")) {
        if (!id %in% colnames(out)) {
            stop(paste0("The id column '", id, "' is not found in the data."))
        }
        ifelse(out$.indicator, as.character(out[[id]]), NA_character_)
    } else if (identical(output, "logical") || identical(output, "bool") || identical(output, "boolean") || identical(output, "indicator")) {
        out$.indicator
    } else { # data, df, data.frame
        out %>%
            filter(!!sym(".indicator")) %>%
            select(-!!sym(".indicator"))
    }
}

#' @rdname clone_selector_utils
#' @keywords internal
.top_long <- function(n, group_by, data, order, id, within, output_within, output) {
    if (!is.null(within)) {
        data2 <- data %>% filter(!!parse_expr(within))
    } else {
        data2 <- data
    }
    # data from clonal_size_data()
    if ("count" %in% colnames(data2) && "fraction" %in% colnames(data2)) {
        data2 <- data2 %>%
            dplyr::group_by(!!!syms(unique(c(id, group_by))), .drop = FALSE) %>%
            summarise(.n = sum(!!sym("count")), .groups = "drop")
    } else {
        data2 <- data2 %>%
            dplyr::group_by(!!!syms(unique(c(id, group_by))), .drop = FALSE) %>%
            summarise(.n = n(), .groups = "drop")
    }

    order <- order %||% "-.n"
    data2 <- data2 %>% arrange(!!parse_expr(order))

    if (!is.null(group_by)) {
        data2 <- data2 %>% dplyr::group_by(!!!syms(group_by))
    }
    # to gain the ids
    out <- data2 %>% mutate(.indicator = row_number() <= n)
    # to attach indicator to the original data
    out <- data %>% left_join(out, by = unique(c(id, group_by)))
    if (isTRUE(output_within)) {
        output_within <- within
    }
    if (!is.null(output_within) && !isFALSE(output_within)) {
        out <- out %>% mutate(.indicator = !!sym(".indicator") & !!parse_expr(output_within))
    }
    .return_what(out, id, output)
}

#' @rdname clone_selector_utils
#' @keywords internal
.top_wide <- function(n, group_by, data, order, id, output_within, output) {
    if (!is.null(group_by)) {
        data <- data %>% dplyr::group_by(!!!syms(group_by))
    }

    if (!is.null(order)) {
        data <- data %>% arrange(!!parse_expr(order))
    }

    out <- data %>% mutate(.indicator = row_number() <= n)
    if (!is.null(output_within) && !isFALSE(output_within)) {
        out <- out %>% mutate(.indicator = !!sym(".indicator") & !!parse_expr(output_within))
    }
    .return_what(out, id, output)
}


#' @rdname clone_selector_utils
#' @keywords internal
.sel_long <- function(expr, group_by, data, id, top, order, within, output_within, output) {
    if (!is.null(within)) {
        data2 <- data %>% filter(!!parse_expr(within))
    } else {
        data2 <- data
    }
    # data from clonal_size_data()
    if ("count" %in% colnames(data2) && "fraction" %in% colnames(data2)) {
        data2 <- data2 %>%
            dplyr::group_by(!!!syms(unique(c(id, group_by))), .drop = FALSE) %>%
            summarise(.n = sum(!!sym("count")), .groups = "drop")
    } else {
        data2 <- data2 %>%
            dplyr::group_by(!!!syms(unique(c(id, group_by))), .drop = FALSE) %>%
            summarise(.n = n(), .groups = "drop")
    }
    if (!is.null(group_by)) {
        data2 <- tidyr::pivot_wider(
            data2,
            # id_cols = unique(c(id, group_by[-1])),
            names_from = !!sym(group_by[1]),
            values_from = !!sym(".n"),
            values_fill = 0
        )
        if (length(group_by) > 1) {
            data2 <- dplyr::group_by(data2, !!!syms(group_by[-1]))
        }
    }

    out <- data2 %>% mutate(.indicator = !!parse_expr(expr))
    if (!is.null(top)) {
        if (is.null(order)) {
            out <- out %>% mutate(
                .indicator.index = ifelse(!!sym(".indicator"), cumsum(!!sym(".indicator")), Inf),
                .indicator = !!sym(".indicator") & !!sym(".indicator.index") <= top
            ) %>% select(-".indicator.index")
        } else {
            out2 <- out %>% arrange(!!parse_expr(order)) %>% mutate(
                .indicator.index = ifelse(!!sym(".indicator"), cumsum(!!sym(".indicator")), Inf),
                .indicator = !!sym(".indicator") & !!sym(".indicator.index") <= top
            )
            out <- out %>%
                left_join(
                    select(out2, !!!syms(unique(c(id, group_by[-1], ".indicator")))),
                    by = unique(c(id, group_by[-1])),
                    suffix = c(".x", "")
                ) %>%
                select(-".indicator.x")
        }
    }
    out <- data %>% left_join(out, by = unique(c(id, group_by[-1])))
    if (isTRUE(output_within)) {
        output_within <- within
    }
    if (!is.null(output_within) && !isFALSE(output_within)) {
        out <- out %>% mutate(.indicator = !!sym(".indicator") & !!parse_expr(output_within))
    }
    .return_what(out, id, output)
}

#' @rdname clone_selector_utils
#' @keywords internal
.sel_wide <- function(expr, group_by, data, id, top, order, output_within, output) {
    if (!is.null(group_by)) {
        data <- data %>% dplyr::group_by(!!!syms(group_by))
    }

    out <- data %>% mutate(.indicator = !!parse_expr(expr))
    if (!is.null(top)) {
        if (is.null(order)) {
            out <- out %>% mutate(
                .indicator.index = ifelse(!!sym(".indicator"), cumsum(!!sym(".indicator")), Inf),
                .indicator = !!sym(".indicator") & !!sym(".indicator.index") <= top
            ) %>% select(-".indicator.index")
        } else {
            out2 <- out %>% arrange(!!parse_expr(order)) %>% mutate(
                .indicator.index = ifelse(!!sym(".indicator"), cumsum(!!sym(".indicator")), Inf),
                .indicator = !!sym(".indicator") & !!sym(".indicator.index") <= top
            )
            out <- out %>%
                left_join(out2[, unique(c(id, group_by, ".indicator"))], by = unique(c(id, group_by)), suffix = c(".x", "")) %>%
                select(-".indicator.x")
        }
    }
    if (!is.null(output_within) && !isFALSE(output_within)) {
        out <- out %>% mutate(.indicator = !!sym(".indicator") & !!parse_expr(output_within))
    }
    .return_what(out, id, output)
}

#' @rdname clone_selector_utils
#' @keywords internal
.bquote <- function(x) {
    if (is.character(x)) {
        if (!suppressWarnings(is.na(as.numeric(x)))) {
            return(x) # if x is a number, return it as is
        }
        if (grepl("`", x)) {
            return(x) # already quoted
        } else {
            return(paste0("`", x, "`")) # backtick-quote the name
        }
    } else {
        return(x)
    }
}

#' @rdname clone_selector_utils
#' @keywords internal
.to_chr <- function(val, expr, allow_bool = FALSE) {
    ok <- tryCatch(
        is.character(val) || is.null(val) || (allow_bool && (isTRUE(val) || isFALSE(val))),
        error = function(e) FALSE
    )
    if (!ok) expr_name(expr) else val
}

#' @rdname clone_selector_utils
#' @keywords internal
.selector_setup <- function(
    orig_data, data, in_form, output, group_by, id, within, output_within,
    envtype, pframe, check_group_by = FALSE, factor_group = FALSE
) {
    if (!is.null(orig_data)) {
        stopifnot("`in_form` must be provided when `data` is provided explictly." = !is.null(in_form))
        in_form <- match.arg(in_form, c("wide", "long"))
        stopifnot("`output` must be provided when `data` is provided explictly." = !is.null(output))
    } else {
        in_form <- in_form %||% ifelse(envtype == "tidy", "long", "wide")
        output <- output %||% ifelse(envtype == "tidy", "id", "data")
    }
    if (envtype == "scplotter" && is.null(group_by)) {
        group_by <- unique(c(pframe$facet_by, pframe$split_by))
    }
    id <- id %||% ifelse(envtype == "tidy", "CTaa", "CloneID")
    if (check_group_by) {
        stopifnot("`group_by` must be provided when `in_form` is 'long'." = !is.null(group_by) || in_form == "wide")
    }
    stopifnot("`within` can only be used with `long` format data." = is.null(within) || in_form == "long")
    stopifnot("`output_within = TRUE` can only be used when `within` is specified." = !isTRUE(output_within) || !is.null(within))
    if (factor_group && in_form == "long" && length(group_by) == 1) {
        data[[group_by]] <- as.factor(data[[group_by]])
    }
    list(data = data, in_form = in_form, output = output, group_by = group_by, id = id)
}

#' Helper functions to select clones based on various criteria
#'
#' @name CloneSelectors
#' @description These helper functions allow for the selection of clones based on various criteria such as size, group comparison, and existence in specific groups.
#' @details These helper functions are designed to be used in a dplyr pipeline or used internally in other scplotter
#' functions to select clones based on various criteria.
#' * When used in a dplyr pipeline, they will return a vector with the same length as the input data, with the selected
#' clones' CTaa values (clone IDs) and NA for others. It is useful for adding a new column to the data frame. For the
#' functions that need `group1`, `group2`, and/or `...`, `group_by` should be provided to specify the grouping columns.
#' Then `group1`, `group2`, and `...` can be the values in the grouping column. To include more grouping columns, just use
#' `c(grouping1, grouping2, ...)`, where `grouping1` is used for values of `group1`, `group2` and `...`; `grouping2` and
#' so on will be kept as the groupings where the clones are selected in each combination of the grouping values.
#' * When used in a scplotter function, they will return a subset of the data frame with only the selected clones.
#' This is useful for filtering the data frame to only include the clones that meet the criteria. It is used internally in
#' some other scplotter functions, such as `ClonalStatPlot`, to select clones. The groupings are also applied, and defaulting
#' to `facet_by` and `split_by` in the parent frame.
#' * When used independently, you should pass the arguments explicitly, such as `group_by` and `output`, to control the
#' behavior and the output of the function.
#'
#' @param n The number of top clones to select or the threshold size.
#' @param expr The expression (in characters) to filter the clones
#' (e.g. "group1 > group2" to select clones where group1 is larger than group2).
#' @param group1 The first group to compare.
#' @param group2 The second group to compare.
#' @param ... More groups to compare.
#' @param include_zeros Whether to include clones with zero size in the comparison.
#' If TRUE, in a comparison (s1 > s2) for a clone to be selected, both s1 and s2 must be greater than 0.
#' If FALSE, only the first group must be greater than the second group.
#' @param group_by The column names in the meta data to group the cells.
#' By default, it is assumed `facet_by` and `split_by` to be in the parent frame if used in scplotter functions.
#' When used in dplyr verbs, it should be a character vector of the grouping columns, where the first column is used to
#' extract the values (count) for `group1`, `group2`, and `...` and the rest are used to keep the groupings.
#' @param top The number of top clones to select based on the expression.
#' If specified, it will select the top N clones that meet the criteria defined by `expr` and ordered by `order`.
#' If `order` is not specified, it will select the top N clones based on the order they appear in the data after filtering by `expr`.
#' @param order The order of the clones to select. It can be an expression to order the clones by a specific column. Only used in `top()`.
#' @param in_form The format of the input data. It can be "long" or "wide".
#' If "long", the data should be in a long format with a column for the clone IDs and a column for the size.
#' If "wide", the data should be in a wide format with columns for the clone IDs and the size for each group.
#' When used in dplyr verbs, it should be "long" by default.
#' If used in scplotter functions, it should be "wide" by default.#'
#' @param data The data frame containing clone information. Default is NULL. If NULL,
#' when used in scplotter functions, it will get data from parent.frame.
#' A typical `data` should have a column named `CloneID` and other columns for the groupings.
#' Supposingly it should be a grouped data frame with the grouping columns.
#' Under each grouping column, the value should be the size of the clone.
#' By default, the data is assumed to be in the parent frame.
#' When used in dplyr verbs, it should be the parent data frame passed to the dplyr verb.
#' @param id The column name that contains the clone ID. Default is "CTaa".
#' @param within An expression passed to subset the data before applying the selection criteria.
#' Only works for `long` format. It is useful when you want to select clones based on the criteria within a specific subset of the data.
#' Note that this subsetting is only applied to determine the selection of clones, not to the final output.
#' So if a cell belongs to a clone that is selected based on the subsetted data, it will be included in the final output, even
#' if it does not meet the `within` criteria. If you want the clones returned to also meet the `within` criteria,
#' you can set `output_within` to TRUE, which will return the clones that meet both the selection criteria and the `within` criteria.
#' @param output_within An expression passed to subset the data after applying the selection criteria.
#' Can work with both `long` and `wide` format.
#' It is useful when you want to return clones that meet both the selection criteria and this criteria.
#' If set to TRUE (only works when `within` is specified), the `within` criteria will be applied to filter the final output to include only the clones that meet both the selection criteria and the `within` criteria.
#' If FALSE or NULL (default), the `within` criteria will only be applied to determine the selection of clones, not to the final output.
#' @param output There are three options for the output: "id" (or "ids"), "logical" (or "bool", "boolean", "indicator"), and "data" (or "df", "data.frame").
#' * "id" (or "ids"): return a vector with the same length as the input data, with the selected clones' CTaa values (clone IDs) and NA for others. It is useful for adding a new column to the data frame.
#' * "logical" (or "bool", "boolean", "indicator"): return a logical vector indicating whether each clone is selected or not. Same as `id` but with TRUE for selected clones and FALSE for others.
#' * "data" (or "df", "data.frame"): return a subset of the data frame with only the selected clones. This is useful for filtering the data frame to only include the clones that meet the criteria. It is used internally in some other scplotter functions, such as `ClonalStatPlot`, to select clones. The groupings are also applied, and defaulting to `facet_by` and `split_by` in the parent frame.
#' By default, it is NULL, which will return "id" when used in dplyr verbs and "data" when used in scplotter functions.
#' @param x The first vector to compare in logical operations (and/or).
#' @param y The second vector to compare in logical operations (and/or).
#' @param ... Additional vectors to compare in logical operations (and/or).
#' @return A vector of CTaas or a data frame with the selected clones based on the criteria.
#' @importFrom rlang parse_expr syms sym caller_env env_has env_get enquo expr_name
#' @importFrom dplyr summarise filter reframe pull mutate select across everything row_number left_join
#' @rdname clone_selectors
#' @examples
#' set.seed(8525)
#' data <- data.frame(
#'     CTaa = c("AA1", "AA2", "AA3", "AA4", "AA5", "AA6", "AA7", "AA8", "AA9", "AA10"),
#'     group1 = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
#'     group2 = c(7, 3, 8, 2, 1, 5, 9, 4, 6, 0),
#'     groups = c("A", "A", "A", "A", "B", "B", "B", "B", "B", "B")
#' )
#' data <- data[order(data$group1 + data$group2, decreasing = TRUE), ]
#' top(3)
#' top(3, group_by = "groups")
#' sel(group1 == 0 | group2 == 0)
#' uniq(group1, group2)
#' shared(group1, group2)
#' gt(group1, group2)
#' lt(group1, group2)
#' le(group1, group2)
#' lt(group1, group2, include_zeros = FALSE)
#' eq(group1, group2)
#' ne(group1, group2)
#'
#' # Use them in a dplyr pipeline
#' data <- tidyr::pivot_longer(data,
#'     cols = c("group1", "group2"),
#'     names_to = "group", values_to = "value"
#' )
#' data <- tidyr::uncount(data, !!rlang::sym("value"))
#' data$subset <- sample(c("S1", "S2"), nrow(data), replace = TRUE)
#' # Take a glimpse of the data
#' data[sample(1:nrow(data), 10), ]
#'
#' unique(dplyr::mutate(data, Top3 = top(3))$Top3)
#' # Note that AA9 also reported in S2, even though the comparison is only applied within S1
#' dplyr::distinct(
#'     dplyr::mutate(data, Top3 = top(3, within = subset == "S1")),
#'     CTaa, subset, Top3
#' )
#' # Note that AA9 is now excluded
#' dplyr::distinct(
#'     dplyr::mutate(data, Top3 = top(3, within = subset == "S1", output_within = TRUE)),
#'     CTaa, subset, Top3
#' )
#' # We can also exclude S1 clones even when the comparison is applied within S1
#' dplyr::distinct(
#'     dplyr::mutate(data, Top3 = top(3, within = subset == "S1", output_within = subset == "S2")),
#'     CTaa, subset, Top3
#' )
#' unique(dplyr::mutate(data, Top3 = top(3, group_by = "groups"))$Top3)
#' unique(dplyr::mutate(data, Unique = sel(group1 == 0 | group2 == 0, group_by = "group"))$Unique)
#' unique(dplyr::mutate(data, UniqueInG1 = uniq(group1, group2, group_by = "group"))$UniqueInG1)
#' unique(dplyr::mutate(data, Shared = shared(group1, group2, group_by = "group"))$Shared)
#' unique(dplyr::mutate(data, Greater = gt(group1, group2, group_by = "group"))$Greater)
#' unique(dplyr::mutate(data, Less = lt(group1, group2, group_by = "group"))$Less)
#' unique(dplyr::mutate(data, LessEqual = le(group1, group2, group_by = "group"))$LessEqual)
#' unique(dplyr::mutate(data, GreaterEqual = ge(group1, group2, group_by = "group"))$GreaterEqual)
#' unique(dplyr::mutate(data, Equal = eq(group1, group2, group_by = "group"))$Equal)
#' unique(dplyr::mutate(data, NotEqual = ne(group1, group2, group_by = "group"))$NotEqual)
#' # Compond expressions
#' unique(
#'     dplyr::mutate(data,
#'         Top3OrEqual = or(top(3), eq(group1, group2, group_by = "group"))
#'     )$Top3OrEqual
#' )
#'
#' unique(
#'     dplyr::mutate(data,
#'         SharedAndGreater = and(
#'             shared(group1, group2, group_by = "group"),
#'             gt(group1, group2, group_by = "group")
#'         )
#'     )$SharedAndGreater
#' )
#' dplyr::mutate(data,
#'     SharedAndGreater = and(
#'         shared(group1, group2, group_by = "group", output = "logical"),
#'         gt(group1, group2, group_by = "group", output = "logical")
#'     )
#' )$SharedAndGreater
NULL

#' @rdname clone_selectors
#' @export
top <- function(n, group_by = NULL, data = NULL, order = NULL, id = NULL, in_form = NULL, within = NULL, output_within = NULL, output = NULL) {
    group_by <- .to_chr(group_by, substitute(group_by))
    order <- .to_chr(order, substitute(order))
    id <- .to_chr(id, substitute(id))
    within <- .to_chr(within, substitute(within))
    output_within <- .to_chr(output_within, substitute(output_within), allow_bool = TRUE)
    envtype <- .get_envtype()
    orig_data <- data
    if (is.null(data)) data <- .get_data(envtype)
    ctx <- .selector_setup(orig_data, data, in_form, output, group_by, id, within, output_within, envtype, parent.frame())
    in_form <- ctx$in_form; output <- ctx$output; group_by <- ctx$group_by; id <- ctx$id
    if (in_form == "long") {
        stopifnot("`id` must be provided when `in_form` is 'long'." = !is.null(id))
        .top_long(n, group_by, data, order, id, within, output_within, output)
    } else {
        .top_wide(n, group_by, data, order, id, output_within, output)
    }
}


#' @rdname clone_selectors
#' @export
sel <- function(
    expr,
    group_by = NULL,
    data = NULL,
    id = NULL,
    in_form = NULL,
    top = NULL,
    order = NULL,
    within = NULL,
    output_within = NULL,
    output = NULL
) {
    expr <- .to_chr(expr, substitute(expr))
    group_by <- .to_chr(group_by, substitute(group_by))
    id <- .to_chr(id, substitute(id))
    order <- .to_chr(order, substitute(order))
    within <- .to_chr(within, substitute(within))
    output_within <- .to_chr(output_within, substitute(output_within), allow_bool = TRUE)
    envtype <- .get_envtype()
    orig_data <- data
    if (is.null(data)) data <- .get_data(envtype)
    ctx <- .selector_setup(orig_data, data, in_form, output, group_by, id, within, output_within, envtype, parent.frame())
    in_form <- ctx$in_form; output <- ctx$output; group_by <- ctx$group_by; id <- ctx$id
    if (in_form == "long") {
        stopifnot("`id` must be provided when `in_form` is 'long'." = !is.null(id))
        .sel_long(expr, group_by, data, id, top, order, within, output_within, output)
    } else {
        .sel_wide(expr, group_by, data, id, top, order, output_within, output)
    }
}

#' @rdname clone_selectors
#' @export
uniq <- function(
    group1,
    group2,
    ...,
    group_by = NULL,
    data = NULL,
    id = NULL,
    in_form = NULL,
    top = NULL,
    order = NULL,
    within = NULL,
    output_within = NULL,
    output = NULL
) {
    group1 <- .to_chr(group1, substitute(group1))
    group2 <- .to_chr(group2, substitute(group2))
    group_by <- .to_chr(group_by, substitute(group_by))
    id <- .to_chr(id, substitute(id))
    order <- .to_chr(order, substitute(order))
    within <- .to_chr(within, substitute(within))
    output_within <- .to_chr(output_within, substitute(output_within), allow_bool = TRUE)
    other_groups <- as.character(substitute(list(...)))[-1]
    envtype <- .get_envtype()
    orig_data <- data
    if (is.null(data)) data <- .get_data(envtype)
    expr <- paste0(.bquote(group1), " > 0 & ", .bquote(group2), " == 0")
    for (g in other_groups) {
        expr <- paste0(expr, " & ", .bquote(g), " == 0")
    }
    ctx <- .selector_setup(orig_data, data, in_form, output, group_by, id, within, output_within, envtype, parent.frame(), check_group_by = TRUE, factor_group = TRUE)
    data <- ctx$data; in_form <- ctx$in_form; output <- ctx$output; group_by <- ctx$group_by; id <- ctx$id
    order <- order %||% paste0("desc(", paste(sapply(c(group1, group2, other_groups), .bquote), collapse = " + "), ")")
    sel(expr, group_by, data, id = id, in_form = in_form, top = top, order = order, within = within, output_within = output_within, output = output)
}

#' @rdname clone_selectors
#' @export
shared <- function(
    group1,
    group2,
    ...,
    group_by = NULL,
    data = NULL,
    id = NULL,
    in_form = NULL,
    top = NULL,
    order = NULL,
    within = NULL,
    output_within = NULL,
    output = NULL
) {
    group1 <- .to_chr(group1, substitute(group1))
    group2 <- .to_chr(group2, substitute(group2))
    group_by <- .to_chr(group_by, substitute(group_by))
    id <- .to_chr(id, substitute(id))
    order <- .to_chr(order, substitute(order))
    within <- .to_chr(within, substitute(within))
    output_within <- .to_chr(output_within, substitute(output_within), allow_bool = TRUE)
    other_groups <- as.character(substitute(list(...)))[-1]
    envtype <- .get_envtype()
    orig_data <- data
    if (is.null(data)) data <- .get_data(envtype)
    expr <- paste0(.bquote(group1), " > 0 & ", .bquote(group2), " > 0")
    for (g in other_groups) {
        expr <- paste0(expr, " & ", .bquote(g), " > 0")
    }
    ctx <- .selector_setup(orig_data, data, in_form, output, group_by, id, within, output_within, envtype, parent.frame(), check_group_by = TRUE, factor_group = TRUE)
    data <- ctx$data; in_form <- ctx$in_form; output <- ctx$output; group_by <- ctx$group_by; id <- ctx$id
    order <- order %||% paste0("desc(", paste(sapply(c(group1, group2, other_groups), .bquote), collapse = " + "), ")")
    sel(expr, group_by, data, id = id, in_form = in_form, top = top, order = order, within = within, output_within = output_within, output = output)
}

#' @rdname clone_selectors
#' @export
gt <- function(
    group1,
    group2,
    include_zeros = TRUE,
    group_by = NULL,
    data = NULL,
    id = NULL,
    in_form = NULL,
    top = NULL,
    order = NULL,
    within = NULL,
    output_within = NULL,
    output = NULL
) {
    group1 <- .to_chr(group1, substitute(group1))
    group2 <- .to_chr(group2, substitute(group2))
    group_by <- .to_chr(group_by, substitute(group_by))
    id <- .to_chr(id, substitute(id))
    order <- .to_chr(order, substitute(order))
    within <- .to_chr(within, substitute(within))
    output_within <- .to_chr(output_within, substitute(output_within), allow_bool = TRUE)
    envtype <- .get_envtype()
    orig_data <- data
    if (is.null(data)) data <- .get_data(envtype)
    expr <- paste0(.bquote(group1), " > ", .bquote(group2))
    if (!include_zeros) {
        expr <- paste0(.bquote(group2), " > 0 & ", expr)
    }
    ctx <- .selector_setup(orig_data, data, in_form, output, group_by, id, within, output_within, envtype, parent.frame(), check_group_by = TRUE, factor_group = TRUE)
    data <- ctx$data; in_form <- ctx$in_form; output <- ctx$output; group_by <- ctx$group_by; id <- ctx$id
    order <- order %||% paste0("desc(", .bquote(group1), " - ", .bquote(group2), ")")
    sel(expr, group_by, data, id = id, in_form = in_form, top = top, order = order, within = within, output_within = output_within, output = output)
}

#' @rdname clone_selectors
#' @export
ge <- function(
    group1,
    group2,
    include_zeros = TRUE,
    group_by = NULL,
    data = NULL,
    id = NULL,
    in_form = NULL,
    top = NULL,
    order = NULL,
    within = NULL,
    output_within = NULL,
    output = NULL
) {
    group1 <- .to_chr(group1, substitute(group1))
    group2 <- .to_chr(group2, substitute(group2))
    group_by <- .to_chr(group_by, substitute(group_by))
    id <- .to_chr(id, substitute(id))
    order <- .to_chr(order, substitute(order))
    within <- .to_chr(within, substitute(within))
    output_within <- .to_chr(output_within, substitute(output_within), allow_bool = TRUE)
    envtype <- .get_envtype()
    orig_data <- data
    if (is.null(data)) data <- .get_data(envtype)
    expr <- paste0(.bquote(group1), " >= ", .bquote(group2))
    if (!include_zeros) {
        expr <- paste0(.bquote(group2), " > 0 & ", expr)
    }
    ctx <- .selector_setup(orig_data, data, in_form, output, group_by, id, within, output_within, envtype, parent.frame(), check_group_by = TRUE, factor_group = TRUE)
    data <- ctx$data; in_form <- ctx$in_form; output <- ctx$output; group_by <- ctx$group_by; id <- ctx$id
    order <- order %||% paste0("desc(", .bquote(group1), " - ", .bquote(group2), ")")
    sel(expr, group_by, data, id = id, in_form = in_form, top = top, order = order, within = within, output_within = output_within, output = output)
}

#' @rdname clone_selectors
#' @export
lt <- function(
    group1,
    group2,
    include_zeros = TRUE,
    group_by = NULL,
    data = NULL,
    id = NULL,
    in_form = NULL,
    top = NULL,
    order = NULL,
    within = NULL,
    output_within = NULL,
    output = NULL
) {
    group1 <- .to_chr(group1, substitute(group1))
    group2 <- .to_chr(group2, substitute(group2))
    group_by <- .to_chr(group_by, substitute(group_by))
    id <- .to_chr(id, substitute(id))
    order <- .to_chr(order, substitute(order))
    within <- .to_chr(within, substitute(within))
    output_within <- .to_chr(output_within, substitute(output_within), allow_bool = TRUE)
    envtype <- .get_envtype()
    orig_data <- data
    if (is.null(data)) data <- .get_data(envtype)
    expr <- paste0(.bquote(group1), " < ", .bquote(group2))
    if (!include_zeros) {
        expr <- paste0(.bquote(group1), " > 0 & ", expr)
    }
    ctx <- .selector_setup(orig_data, data, in_form, output, group_by, id, within, output_within, envtype, parent.frame(), check_group_by = TRUE, factor_group = TRUE)
    data <- ctx$data; in_form <- ctx$in_form; output <- ctx$output; group_by <- ctx$group_by; id <- ctx$id
    order <- order %||% paste0(.bquote(group1), " - ", .bquote(group2))
    sel(expr, group_by, data, id = id, in_form = in_form, top = top, order = order, within = within, output_within = output_within, output = output)
}

#' @rdname clone_selectors
#' @export
le <- function(
    group1,
    group2,
    include_zeros = TRUE,
    group_by = NULL,
    data = NULL,
    id = NULL,
    in_form = NULL,
    top = NULL,
    order = NULL,
    within = NULL,
    output_within = NULL,
    output = NULL
) {
    group1 <- .to_chr(group1, substitute(group1))
    group2 <- .to_chr(group2, substitute(group2))
    group_by <- .to_chr(group_by, substitute(group_by))
    id <- .to_chr(id, substitute(id))
    order <- .to_chr(order, substitute(order))
    within <- .to_chr(within, substitute(within))
    output_within <- .to_chr(output_within, substitute(output_within), allow_bool = TRUE)
    envtype <- .get_envtype()
    orig_data <- data
    if (is.null(data)) data <- .get_data(envtype)
    expr <- paste0(.bquote(group1), " <= ", .bquote(group2))
    if (!include_zeros) {
        expr <- paste0(.bquote(group1), " > 0 & ", expr)
    }
    ctx <- .selector_setup(orig_data, data, in_form, output, group_by, id, within, output_within, envtype, parent.frame(), check_group_by = TRUE, factor_group = TRUE)
    data <- ctx$data; in_form <- ctx$in_form; output <- ctx$output; group_by <- ctx$group_by; id <- ctx$id
    order <- order %||% paste0(.bquote(group1), " - ", .bquote(group2))
    sel(expr, group_by, data, id = id, in_form = in_form, top = top, order = order, within = within, output_within = output_within, output = output)
}

#' @rdname clone_selectors
#' @export
eq <- function(
    group1,
    group2,
    group_by = NULL,
    data = NULL,
    id = NULL,
    top = NULL,
    order = NULL,
    in_form = NULL,
    within = NULL,
    output_within = NULL,
    output = NULL
) {
    group1 <- .to_chr(group1, substitute(group1))
    group2 <- .to_chr(group2, substitute(group2))
    group_by <- .to_chr(group_by, substitute(group_by))
    id <- .to_chr(id, substitute(id))
    order <- .to_chr(order, substitute(order))
    within <- .to_chr(within, substitute(within))
    output_within <- .to_chr(output_within, substitute(output_within), allow_bool = TRUE)
    envtype <- .get_envtype()
    orig_data <- data
    if (is.null(data)) data <- .get_data(envtype)
    expr <- paste0(.bquote(group1), " == ", .bquote(group2))
    ctx <- .selector_setup(orig_data, data, in_form, output, group_by, id, within, output_within, envtype, parent.frame(), check_group_by = TRUE, factor_group = TRUE)
    data <- ctx$data; in_form <- ctx$in_form; output <- ctx$output; group_by <- ctx$group_by; id <- ctx$id
    order <- order %||% paste0("desc(", .bquote(group1), " + ", .bquote(group2), ")")
    sel(expr, group_by, data, id = id, in_form = in_form, top = top, order = order, within = within, output_within = output_within, output = output)
}

#' @rdname clone_selectors
#' @export
ne <- function(
    group1,
    group2,
    include_zeros = TRUE,
    group_by = NULL,
    data = NULL,
    id = NULL,
    in_form = NULL,
    top = NULL,
    order = NULL,
    within = NULL,
    output_within = NULL,
    output = NULL
) {
    group1 <- .to_chr(group1, substitute(group1))
    group2 <- .to_chr(group2, substitute(group2))
    group_by <- .to_chr(group_by, substitute(group_by))
    id <- .to_chr(id, substitute(id))
    order <- .to_chr(order, substitute(order))
    within <- .to_chr(within, substitute(within))
    output_within <- .to_chr(output_within, substitute(output_within), allow_bool = TRUE)
    envtype <- .get_envtype()
    orig_data <- data
    if (is.null(data)) data <- .get_data(envtype)
    expr <- paste0(.bquote(group1), " != ", .bquote(group2))
    if (!include_zeros) {
        expr <- paste0(.bquote(group1), " > 0 & ", .bquote(group2), " > 0 & ", expr)
    }
    ctx <- .selector_setup(orig_data, data, in_form, output, group_by, id, within, output_within, envtype, parent.frame(), check_group_by = TRUE, factor_group = TRUE)
    data <- ctx$data; in_form <- ctx$in_form; output <- ctx$output; group_by <- ctx$group_by; id <- ctx$id
    order <- order %||% paste0("desc(abs(", .bquote(group1), " - ", .bquote(group2), "))")
    sel(expr, group_by, data, id = id, in_form = in_form, top = top, order = order, within = within, output_within = output_within, output = output)
}

#' @rdname clone_selectors
#' @export
and <- function(x, y, ...) {
    type_of <- function(o) {
        if (is.character(o) || anyNA(o)) {
            "id"
        } else if (is.logical(o)) {
            "logical"
        } else {
            stop("[and] Unsupported output type. The output of the elements must be either 'id' or 'logical'.")
        }
    }
    length_x <- length(x)
    if (length(y) != length_x) {
        stop("[and] The length of all elements must be the same. The first element has length ", length_x, " but the second element has length ", length(y), ".")
    }
    type_x <- type_of(x)
    if (type_of(y) != type_x) {
        stop("[and] The output types of all elements must be the same. The first element has output type '", type_x, "' but the second element has output type '", type_of(y), "'.")
    }
    if (type_x == "logical") {
        mask <- x & y
    } else {
        mask <- is.na(x) | is.na(y)
    }
    for (arg in list(...)) {
        if (length(arg) != length_x) {
            stop("[and] The length of all elements must be the same. The first element has length ", length_x, " but one of the other elements has length ", length(arg), ".")
        }
        if (type_of(arg) != type_x) {
            stop("[and] The output types of all elements must be the same. The first element has output type '", type_x, "' but one of the other elements has output type '", type_of(arg), "'.")
        }
        if (type_x == "logical") {
            mask <- mask & arg
        } else {
            mask <- mask | is.na(arg)
        }
    }

    if (type_x == "logical") {
        return(mask)
    } else {
        x[mask] <- NA
        return(x)
    }
}

#' @rdname clone_selectors
#' @export
or <- function(x, y, ...) {
    type_of <- function(o) {
        if (is.character(o) || anyNA(o)) {
            "id"
        } else if (is.logical(o)) {
            "logical"
        } else {
            stop("[or] Unsupported output type. The output of the elements must be either 'id' or 'logical'.")
        }
    }
    length_x <- length(x)
    if (length(y) != length_x) {
        stop("[or] The length of all elements must be the same. The first element has length ", length_x, " but the second element has length ", length(y), ".")
    }
    type_x <- type_of(x)
    if (type_of(y) != type_x) {
        stop("[or] The output types of all elements must be the same. The first element has output type '", type_x, "' but the second element has output type '", type_of(y), "'.")
    }
    if (type_x == "logical") {
        out <- x | y
    } else {
        out <- dplyr::coalesce(x, y)
    }
    for (arg in list(...)) {
        if (length(arg) != length_x) {
            stop("[or] The length of all elements must be the same. The first element has length ", length_x, " but one of the other elements has length ", length(arg), ".")
        }
        if (type_of(arg) != type_x) {
            stop("[or] The output types of all elements must be the same. The first element has output type '", type_x, "' but one of the other elements has output type '", type_of(arg), "'.")
        }
        if (type_x == "logical") {
            out <- out | arg
        } else {
            out <- dplyr::coalesce(out, arg)
        }
    }

    return(out)
}
