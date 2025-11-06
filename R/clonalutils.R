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
    data <- merge_clonal_groupings(data, groupings)

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
    gv_pairs <- as.list(as.data.frame(combn(all_gvalues, 2, simplify = TRUE)))
    clone_call <- .theCall(data, clone_call)

    do.call(rbind, lapply(gv_pairs, function(gv) {
        x_df <- as.data.frame(table(data[[gv[1]]][, clone_call]))
        x_col <- colnames(x_df)[2] <- paste0(gv[1], ".count")
        y_df <- as.data.frame(table(data[[gv[2]]][, clone_call]))
        y_col <- colnames(y_df)[2] <- paste0(gv[2], ".count")

        mat <- merge(x_df, y_df, by = "Var1", all = TRUE)
        mat[is.na(mat)] <- 0
        mat[, paste0(gv[1], ".fraction")] <- mat[, x_col]/sum(mat[,  x_col])
        mat[, paste0(gv[2], ".fraction")] <- mat[, y_col]/sum(mat[,  y_col])

        mat <- pivot_longer(
            mat,
            cols = -"Var1",
            names_to = c(".group", ".value"),
            names_sep = "\\."
        )
        mat
    })) %>%
        distinct(!!sym("Var1"), !!sym(".group"), .keep_all = TRUE) %>%
        separate(".group", into = groupings, sep = " // ") %>%
        filter(!!sym("count") > 0) %>%
        rename(CloneID = "Var1")
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
#' @importFrom dplyr filter left_join mutate group_by summarise arrange
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
.return_what <- function(out, id, return_ids) {
    if (return_ids) {
        if (!id %in% colnames(out)) {
            stop(paste0("The id column '", id, "' is not found in the data."))
        }
        ifelse(out$.indicator, out[[id]], NA)
    } else {
        out %>% filter(!!sym(".indicator")) %>% select(-!!sym(".indicator"))
    }
}

#' @rdname clone_selector_utils
#' @keywords internal
.top_long <- function(n, groups, data, order, id, return_ids) {
    data2 <- data %>%
        group_by(!!!syms(unique(c(id, groups)))) %>%
        summarise(.n = n(), .groups = "drop")

    order <- order %||% "-.n"
    data2 <- data2 %>% arrange(!!parse_expr(order))

    if (!is.null(groups)) {
        data2 <- data2 %>% group_by(!!!syms(groups))
    }

    out <- data2 %>% mutate(.indicator = row_number() <= n)
    out <- data %>% left_join(out, by = unique(c(id, groups)))
    .return_what(out, id, return_ids)
}

#' @rdname clone_selector_utils
#' @keywords internal
.top_wide <- function(n, groups, data, order, id, return_ids) {
    if (!is.null(groups)) {
        data <- data %>% group_by(!!!syms(groups))
    }

    if (!is.null(order)) {
        data <- data %>% arrange(!!parse_expr(order))
    }

    out <- data %>% mutate(.indicator = row_number() <= n)
    .return_what(out, id, return_ids)
}


#' @rdname clone_selector_utils
#' @keywords internal
.sel_long <- function(expr, groups, data, id, return_ids) {
    data2 <- data %>%
        group_by(!!!syms(unique(c(id, groups)))) %>%
        summarise(.n = n(), .groups = "drop")

    if (!is.null(groups)) {
        data2 <- tidyr::pivot_wider(
            data2,
            names_from = !!sym(groups[1]),
            values_from = !!sym(".n"),
            values_fill = 0
        )
        if (length(groups) > 1) {
            data2 <- dplyr::group_by(data2, !!!syms(groups[-1]))
        }
    }

    out <- data2 %>% mutate(.indicator = !!parse_expr(expr))
    out <- data %>% left_join(out, by = unique(c(id, groups[-1])))
    .return_what(out, id, return_ids)
}

#' @rdname clone_selector_utils
#' @keywords internal
.sel_wide <- function(expr, groups, data, id, return_ids) {
    if (!is.null(groups)) {
        data <- data %>% group_by(!!!syms(groups))
    }

    out <- data %>% mutate(.indicator = !!parse_expr(expr))
    .return_what(out, id, return_ids)
}

#' @rdname clone_selector_utils
#' @keywords internal
.bquote <- function(x) {
    if (is.character(x)) {
        if (!suppressWarnings(is.na(as.numeric(x)))) {
            return(x)  # if x is a number, return it as is
        }
        if (grepl("`", x)) {
            return(x)  # already quoted
        } else {
            return(paste0("`", x, "`"))  # backtick-quote the name
        }
    } else {
        return(x)
    }
}

#' Helper functions to select clones based on various criteria
#'
#' @name CloneSelectors
#' @description These helper functions allow for the selection of clones based on various criteria such as size, group comparison, and existence in specific groups.
#' @details These helper functions are designed to be used in a dplyr pipeline or used internally in other scplotter
#' functions to select clones based on various criteria.
#' * When used in a dplyr pipeline, they will return a vector with the same length as the input data, with the selected
#' clones' CTaa values (clone IDs) and NA for others. It is useful for adding a new column to the data frame. For the
#' functions that need `group1`, `group2`, and/or `...`, `groups` should be provided to specify the grouping columns.
#' Then `group1`, `group2`, and `...` can be the values in the grouping column. To include more grouping columns, just use
#' `c(grouping1, grouping2, ...)`, where `grouping1` is used for values of `group1`, `group2` and `...`; `grouping2` and
#' so on will be kept as the groupings where the clones are selected in each combination of the grouping values.
#' * When used in a scplotter function, they will return a subset of the data frame with only the selected clones.
#' This is useful for filtering the data frame to only include the clones that meet the criteria. It is used internally in
#' some other scplotter functions, such as `ClonalStatPlot`, to select clones. The groupings are also applied, and defaulting
#' to `facet_by` and `split_by` in the parent frame.
#' * When used independently, you should pass the arguments explicitly, such as `groups` and `return_ids`, to control the
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
#' @param groups The column names in the meta data to group the cells.
#' By default, it is assumed `facet_by` and `split_by` to be in the parent frame if used in scplotter functions.
#' When used in dplyr verbs, it should be a character vector of the grouping columns, where the first column is used to
#' extract the values (count) for `group1`, `group2`, and `...` and the rest are used to keep the groupings.
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
#' @param return_ids If TRUE, the function returns a vector with the same length as the data, with CTaa values for selected clones and NA for others.
#' If FALSE, it returns a subset data frame with only the selected clones.
#' Default is NULL, which will be determined based on the data. If the function is used in a context of dplyr verbs, it defaults to TRUE.
#' Otherwise, it defaults to FALSE.
#' @param x The first vector to compare in logical operations (and/or).
#' @param y The second vector to compare in logical operations (and/or).
#' @return A vector of CTaas or a data frame with the selected clones based on the criteria.
#' @importFrom rlang parse_expr syms sym caller_env env_has env_get enquo expr_name
#' @importFrom dplyr group_by summarise filter reframe pull mutate select across everything row_number left_join
#' @rdname clone_selectors
#' @examples
#' data <- data.frame(
#'     CTaa = c("AA1", "AA2", "AA3", "AA4", "AA5", "AA6", "AA7", "AA8", "AA9", "AA10"),
#'     group1 = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
#'     group2 = c(7, 3, 8, 2, 1, 5, 9, 4, 6, 0),
#'     groups = c("A", "A", "A", "A", "B", "B", "B", "B", "B", "B")
#' )
#' data <- data[order(data$group1 + data$group2, decreasing = TRUE), ]
#' top(3)
#' top(3, groups = "groups")
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
#' data <- tidyr::pivot_longer(data, cols = c("group1", "group2"),
#'     names_to = "group", values_to = "value")
#' data <- tidyr::uncount(data, !!rlang::sym("value"))
#'
#' unique(dplyr::mutate(data, Top3 = top(3))$Top3)
#' unique(dplyr::mutate(data, Top3 = top(3, groups = "groups"))$Top3)
#' unique(dplyr::mutate(data, Unique = sel(group1 == 0 | group2 == 0, groups = "group"))$Unique)
#' unique(dplyr::mutate(data, UniqueInG1 = uniq(group1, group2, groups = "group"))$UniqueInG1)
#' unique(dplyr::mutate(data, Shared = shared(group1, group2, groups = "group"))$Shared)
#' unique(dplyr::mutate(data, Greater = gt(group1, group2, groups = "group"))$Greater)
#' unique(dplyr::mutate(data, Less = lt(group1, group2, groups = "group"))$Less)
#' unique(dplyr::mutate(data, LessEqual = le(group1, group2, groups = "group"))$LessEqual)
#' unique(dplyr::mutate(data, GreaterEqual = ge(group1, group2, groups = "group"))$GreaterEqual)
#' unique(dplyr::mutate(data, Equal = eq(group1, group2, groups = "group"))$Equal)
#' unique(dplyr::mutate(data, NotEqual = ne(group1, group2, groups = "group"))$NotEqual)
#' # Compond expressions
#' unique(
#'   dplyr::mutate(data,
#'      Top3OrEqual = or(top(3), eq(group1, group2, groups = "group")))$Top3OrEqual
#' )
#'
#' unique(
#'   dplyr::mutate(data,
#'      SharedAndGreater = and(
#'         shared(group1, group2, groups = "group"),
#'         gt(group1, group2, groups = "group")
#'      ))$SharedAndGreater
#' )
NULL

#' @rdname clone_selectors
#' @export
top <- function(n, groups = NULL, data = NULL, order = NULL, id = NULL, in_form = NULL, return_ids = NULL) {
    if (!tryCatch(is.character(groups) || is.null(groups), error = function(e) FALSE)) groups <- expr_name(substitute(groups))
    if (!tryCatch(is.character(order) || is.null(order), error = function(e) FALSE)) order <- expr_name(substitute(order))
    if (!tryCatch(is.character(id) || is.null(id), error = function(e) FALSE)) id <- expr_name(substitute(id))
    envtype <- .get_envtype()
    if (!is.null(data)) {
        stopifnot("`in_form` must be provided when `data` is provided explictly." = !is.null(in_form))
        in_form <- match.arg(in_form, c("wide", "long"))
        stopifnot("`return_ids` must be provided when `data` is provided explictly." = !is.null(return_ids))
    } else {
        data <- .get_data(envtype)  # envtype is ensured to be "tidy" or "scplotter"
        in_form <- in_form %||% ifelse(envtype == "tidy", "long", "wide")
        return_ids <- return_ids %||% (envtype == "tidy")
    }
    if (envtype == "scplotter" && is.null(groups)) {
        env <- caller_env()
        groups <- unique(c(env$facet_by, env$split_by))
    }
    id <- id %||% ifelse(envtype == "tidy", "CTaa", "CloneID")
    if (in_form == "long") {
        stopifnot("`id` must be provided when `in_form` is 'long'." = !is.null(id))
        .top_long(n, groups, data, order, id, return_ids)
    } else {
        .top_wide(n, groups, data, order, id, return_ids)
    }
}


#' @rdname clone_selectors
#' @export
sel <- function(expr, groups = NULL, data = NULL, id = NULL, in_form = NULL, return_ids = NULL) {
    if (!tryCatch(is.character(expr) || is.null(expr), error = function(e) FALSE)) expr <- expr_name(substitute(expr))
    if (!tryCatch(is.character(groups) || is.null(groups), error = function(e) FALSE)) groups <- expr_name(substitute(groups))
    if (!tryCatch(is.character(id) || is.null(id), error = function(e) FALSE)) id <- expr_name(substitute(id))
    envtype <- .get_envtype()
    if (!is.null(data)) {
        stopifnot("`in_form` must be provided when `data` is provided explictly." = !is.null(in_form))
        in_form <- match.arg(in_form, c("wide", "long"))
        stopifnot("`return_ids` must be provided when `data` is provided explictly." = !is.null(return_ids))
    } else {
        data <- .get_data(envtype)  # envtype is ensured to be "tidy" or "scplotter"
        in_form <- in_form %||% ifelse(envtype == "tidy", "long", "wide")
        return_ids <- return_ids %||% (envtype == "tidy")
    }
    if (envtype == "scplotter" && is.null(groups)) {
        env <- caller_env()
        groups <- unique(c(env$facet_by, env$split_by))
    }
    id <- id %||% ifelse(envtype == "tidy", "CTaa", "CloneID")
    if (in_form == "long") {
        stopifnot("`id` must be provided when `in_form` is 'long'." = !is.null(id))
        .sel_long(expr, groups, data, id, return_ids)
    } else {
        .sel_wide(expr, groups, data, id, return_ids)
    }
}

#' @rdname clone_selectors
#' @export
uniq <- function(group1, group2, ..., groups = NULL, data = NULL, id = NULL, in_form = NULL, return_ids = NULL) {
    if (!tryCatch(is.character(group1) || is.null(group1), error = function(e) FALSE)) group1 <- expr_name(substitute(group1))
    if (!tryCatch(is.character(group2) || is.null(group2), error = function(e) FALSE)) group2 <- expr_name(substitute(group2))
    if (!tryCatch(is.character(groups) || is.null(groups), error = function(e) FALSE)) groups <- expr_name(substitute(groups))
    if (!tryCatch(is.character(id) || is.null(id), error = function(e) FALSE)) id <- expr_name(substitute(id))
    other_groups <- as.character(substitute(list(...)))[-1]
    envtype <- .get_envtype()
    if (!is.null(data)) {
        stopifnot("`in_form` must be provided when `data` is provided explictly." = !is.null(in_form))
        in_form <- match.arg(in_form, c("wide", "long"))
        stopifnot("`return_ids` must be provided when `data` is provided explictly." = !is.null(return_ids))
    } else {
        data <- .get_data(envtype)  # envtype is ensured to be "tidy" or "scplotter"
        in_form <- in_form %||% ifelse(envtype == "tidy", "long", "wide")
        return_ids <- return_ids %||% (envtype == "tidy")
    }
    if (envtype == "scplotter" && is.null(groups)) {
        env <- caller_env()
        groups <- unique(c(env$facet_by, env$split_by))
    }
    id <- id %||% ifelse(envtype == "tidy", "CTaa", "CloneID")
    expr <- paste0(.bquote(group1), " > 0 & ", .bquote(group2), " == 0")
    for (g in other_groups) {
        expr <- paste0(expr, " & ", .bquote(g), " == 0")
    }
    stopifnot("`groups` must be provided when `in_form` is 'long'." = !is.null(groups) || in_form == "wide")
    sel(expr, groups, data, id = id, in_form = in_form, return_ids = return_ids)
}

#' @rdname clone_selectors
#' @export
shared <- function(group1, group2, ..., groups = NULL, data = NULL, id = NULL, in_form = NULL, return_ids = NULL) {
    if (!tryCatch(is.character(group1) || is.null(group1), error = function(e) FALSE)) group1 <- expr_name(substitute(group1))
    if (!tryCatch(is.character(group2) || is.null(group2), error = function(e) FALSE)) group2 <- expr_name(substitute(group2))
    if (!tryCatch(is.character(groups) || is.null(groups), error = function(e) FALSE)) groups <- expr_name(substitute(groups))
    if (!tryCatch(is.character(id) || is.null(id), error = function(e) FALSE)) id <- expr_name(substitute(id))
    other_groups <- as.character(substitute(list(...)))[-1]
    envtype <- .get_envtype()
    if (!is.null(data)) {
        stopifnot("`in_form` must be provided when `data` is provided explictly." = !is.null(in_form))
        in_form <- match.arg(in_form, c("wide", "long"))
        stopifnot("`return_ids` must be provided when `data` is provided explictly." = !is.null(return_ids))
    } else {
        data <- .get_data(envtype)  # envtype is ensured to be "tidy" or "scplotter"
        in_form <- in_form %||% ifelse(envtype == "tidy", "long", "wide")
        return_ids <- return_ids %||% (envtype == "tidy")
    }
    if (envtype == "scplotter" && is.null(groups)) {
        env <- caller_env()
        groups <- unique(c(env$facet_by, env$split_by))
    }
    id <- id %||% ifelse(envtype == "tidy", "CTaa", "CloneID")
    expr <- paste0(.bquote(group1), " > 0 & ", .bquote(group2), " > 0")
    for (g in other_groups) {
        expr <- paste0(expr, " & ", .bquote(g), " > 0")
    }
    stopifnot("`groups` must be provided when `in_form` is 'long'." = !is.null(groups) || in_form == "wide")
    sel(expr, groups, data, id = id, in_form = in_form, return_ids = return_ids)
}

#' @rdname clone_selectors
#' @export
gt <- function(group1, group2, include_zeros = TRUE, groups = NULL, data = NULL, id = NULL, in_form = NULL, return_ids = NULL) {
    if (!tryCatch(is.character(group1) || is.null(group1), error = function(e) FALSE)) group1 <- expr_name(substitute(group1))
    if (!tryCatch(is.character(group2) || is.null(group2), error = function(e) FALSE)) group2 <- expr_name(substitute(group2))
    if (!tryCatch(is.character(groups) || is.null(groups), error = function(e) FALSE)) groups <- expr_name(substitute(groups))
    if (!tryCatch(is.character(id) || is.null(id), error = function(e) FALSE)) id <- expr_name(substitute(id))
    envtype <- .get_envtype()
    if (!is.null(data)) {
        stopifnot("`in_form` must be provided when `data` is provided explictly." = !is.null(in_form))
        in_form <- match.arg(in_form, c("wide", "long"))
        stopifnot("`return_ids` must be provided when `data` is provided explictly." = !is.null(return_ids))
    } else {
        data <- .get_data(envtype)  # envtype is ensured to be "tidy" or "scplotter"
        in_form <- in_form %||% ifelse(envtype == "tidy", "long", "wide")
        return_ids <- return_ids %||% (envtype == "tidy")
    }
    if (envtype == "scplotter" && is.null(groups)) {
        env <- caller_env()
        groups <- unique(c(env$facet_by, env$split_by))
    }
    id <- id %||% ifelse(envtype == "tidy", "CTaa", "CloneID")
    expr <- paste0(.bquote(group1), " > ", .bquote(group2))
    if (!include_zeros) {
        expr <- paste0(.bquote(group2), " > 0 & ", expr)
    }
    stopifnot("`groups` must be provided when `in_form` is 'long'." = !is.null(groups) || in_form == "wide")
    sel(expr, groups, data, id = id, in_form = in_form, return_ids = return_ids)
}

#' @rdname clone_selectors
#' @export
ge <- function(group1, group2, include_zeros = TRUE, groups = NULL, data = NULL, id = NULL, in_form = NULL, return_ids = NULL) {
    if (!tryCatch(is.character(group1) || is.null(group1), error = function(e) FALSE)) group1 <- expr_name(substitute(group1))
    if (!tryCatch(is.character(group2) || is.null(group2), error = function(e) FALSE)) group2 <- expr_name(substitute(group2))
    if (!tryCatch(is.character(groups) || is.null(groups), error = function(e) FALSE)) groups <- expr_name(substitute(groups))
    if (!tryCatch(is.character(id) || is.null(id), error = function(e) FALSE)) id <- expr_name(substitute(id))
    envtype <- .get_envtype()
    if (!is.null(data)) {
        stopifnot("`in_form` must be provided when `data` is provided explictly." = !is.null(in_form))
        in_form <- match.arg(in_form, c("wide", "long"))
        stopifnot("`return_ids` must be provided when `data` is provided explictly." = !is.null(return_ids))
    } else {
        data <- .get_data(envtype)  # envtype is ensured to be "tidy" or "scplotter"
        in_form <- in_form %||% ifelse(envtype == "tidy", "long", "wide")
        return_ids <- return_ids %||% (envtype == "tidy")
    }
    if (envtype == "scplotter" && is.null(groups)) {
        env <- caller_env()
        groups <- unique(c(env$facet_by, env$split_by))
    }
    id <- id %||% ifelse(envtype == "tidy", "CTaa", "CloneID")
    expr <- paste0(.bquote(group1), " >= ", .bquote(group2))
    if (!include_zeros) {
        expr <- paste0(.bquote(group2), " > 0 & ", expr)
    }
    stopifnot("`groups` must be provided when `in_form` is 'long'." = !is.null(groups) || in_form == "wide")
    sel(expr, groups, data, id = id, in_form = in_form, return_ids = return_ids)
}

#' @rdname clone_selectors
#' @export
lt <- function(group1, group2, include_zeros = TRUE, groups = NULL, data = NULL, id = NULL, in_form = NULL, return_ids = NULL) {
    if (!tryCatch(is.character(group1) || is.null(group1), error = function(e) FALSE)) group1 <- expr_name(substitute(group1))
    if (!tryCatch(is.character(group2) || is.null(group2), error = function(e) FALSE)) group2 <- expr_name(substitute(group2))
    if (!tryCatch(is.character(groups) || is.null(groups), error = function(e) FALSE)) groups <- expr_name(substitute(groups))
    if (!tryCatch(is.character(id) || is.null(id), error = function(e) FALSE)) id <- expr_name(substitute(id))
    envtype <- .get_envtype()
    if (!is.null(data)) {
        stopifnot("`in_form` must be provided when `data` is provided explictly." = !is.null(in_form))
        in_form <- match.arg(in_form, c("wide", "long"))
        stopifnot("`return_ids` must be provided when `data` is provided explictly." = !is.null(return_ids))
    } else {
        data <- .get_data(envtype)  # envtype is ensured to be "tidy" or "scplotter"
        in_form <- in_form %||% ifelse(envtype == "tidy", "long", "wide")
        return_ids <- return_ids %||% (envtype == "tidy")
    }
    if (envtype == "scplotter" && is.null(groups)) {
        env <- caller_env()
        groups <- unique(c(env$facet_by, env$split_by))
    }
    id <- id %||% ifelse(envtype == "tidy", "CTaa", "CloneID")
    expr <- paste0(.bquote(group1), " < ", .bquote(group2))
    if (!include_zeros) {
        expr <- paste0(.bquote(group1), " > 0 & ", expr)
    }
    stopifnot("`groups` must be provided when `in_form` is 'long'." = !is.null(groups) || in_form == "wide")
    sel(expr, groups, data, id = id, in_form = in_form, return_ids = return_ids)
}

#' @rdname clone_selectors
#' @export
le <- function(group1, group2, include_zeros = TRUE, groups = NULL, data = NULL, id = NULL, in_form = NULL, return_ids = NULL) {
    if (!tryCatch(is.character(group1) || is.null(group1), error = function(e) FALSE)) group1 <- expr_name(substitute(group1))
    if (!tryCatch(is.character(group2) || is.null(group2), error = function(e) FALSE)) group2 <- expr_name(substitute(group2))
    if (!tryCatch(is.character(groups) || is.null(groups), error = function(e) FALSE)) groups <- expr_name(substitute(groups))
    if (!tryCatch(is.character(id) || is.null(id), error = function(e) FALSE)) id <- expr_name(substitute(id))
    envtype <- .get_envtype()
    if (!is.null(data)) {
        stopifnot("`in_form` must be provided when `data` is provided explictly." = !is.null(in_form))
        in_form <- match.arg(in_form, c("wide", "long"))
        stopifnot("`return_ids` must be provided when `data` is provided explictly." = !is.null(return_ids))
    } else {
        data <- .get_data(envtype)  # envtype is ensured to be "tidy" or "scplotter"
        in_form <- in_form %||% ifelse(envtype == "tidy", "long", "wide")
        return_ids <- return_ids %||% (envtype == "tidy")
    }
    if (envtype == "scplotter" && is.null(groups)) {
        env <- caller_env()
        groups <- unique(c(env$facet_by, env$split_by))
    }
    id <- id %||% ifelse(envtype == "tidy", "CTaa", "CloneID")
    expr <- paste0(.bquote(group1), " <= ", .bquote(group2))
    if (!include_zeros) {
        expr <- paste0(.bquote(group1), " > 0 & ", expr)
    }
    stopifnot("`groups` must be provided when `in_form` is 'long'." = !is.null(groups) || in_form == "wide")
    sel(expr, groups, data, id = id, in_form = in_form, return_ids = return_ids)
}

#' @rdname clone_selectors
#' @export
eq <- function(group1, group2, groups = NULL, data = NULL, id = NULL, in_form = NULL, return_ids = NULL) {
    if (!tryCatch(is.character(group1) || is.null(group1), error = function(e) FALSE)) group1 <- expr_name(substitute(group1))
    if (!tryCatch(is.character(group2) || is.null(group2), error = function(e) FALSE)) group2 <- expr_name(substitute(group2))
    if (!tryCatch(is.character(groups) || is.null(groups), error = function(e) FALSE)) groups <- expr_name(substitute(groups))
    if (!tryCatch(is.character(id) || is.null(id), error = function(e) FALSE)) id <- expr_name(substitute(id))
    envtype <- .get_envtype()
    if (!is.null(data)) {
        stopifnot("`in_form` must be provided when `data` is provided explictly." = !is.null(in_form))
        in_form <- match.arg(in_form, c("wide", "long"))
        stopifnot("`return_ids` must be provided when `data` is provided explictly." = !is.null(return_ids))
    } else {
        data <- .get_data(envtype)  # envtype is ensured to be "tidy" or "scplotter"
        in_form <- in_form %||% ifelse(envtype == "tidy", "long", "wide")
        return_ids <- return_ids %||% (envtype == "tidy")
    }
    if (envtype == "scplotter" && is.null(groups)) {
        env <- caller_env()
        groups <- unique(c(env$facet_by, env$split_by))
    }
    id <- id %||% ifelse(envtype == "tidy", "CTaa", "CloneID")
    expr <- paste0(.bquote(group1), " == ", .bquote(group2))
    stopifnot("`groups` must be provided when `in_form` is 'long'." = !is.null(groups) || in_form == "wide")
    sel(expr, groups, data, id = id, in_form = in_form, return_ids = return_ids)
}

#' @rdname clone_selectors
#' @export
ne <- function(group1, group2, include_zeros = TRUE, groups = NULL, data = NULL, id = NULL, in_form = NULL, return_ids = NULL) {
    if (!tryCatch(is.character(group1) || is.null(group1), error = function(e) FALSE)) group1 <- expr_name(substitute(group1))
    if (!tryCatch(is.character(group2) || is.null(group2), error = function(e) FALSE)) group2 <- expr_name(substitute(group2))
    if (!tryCatch(is.character(groups) || is.null(groups), error = function(e) FALSE)) groups <- expr_name(substitute(groups))
    if (!tryCatch(is.character(id) || is.null(id), error = function(e) FALSE)) id <- expr_name(substitute(id))
    envtype <- .get_envtype()
    if (!is.null(data)) {
        stopifnot("`in_form` must be provided when `data` is provided explictly." = !is.null(in_form))
        in_form <- match.arg(in_form, c("wide", "long"))
        stopifnot("`return_ids` must be provided when `data` is provided explictly." = !is.null(return_ids))
    } else {
        data <- .get_data(envtype)  # envtype is ensured to be "tidy" or "scplotter"
        in_form <- in_form %||% ifelse(envtype == "tidy", "long", "wide")
        return_ids <- return_ids %||% (envtype == "tidy")
    }
    if (envtype == "scplotter" && is.null(groups)) {
        env <- caller_env()
        groups <- unique(c(env$facet_by, env$split_by))
    }
    id <- id %||% ifelse(envtype == "tidy", "CTaa", "CloneID")
    expr <- paste0(.bquote(group1), " != ", .bquote(group2))
    if (!include_zeros) {
        expr <- paste0(.bquote(group1), " > 0 & ", .bquote(group2), " > 0 & ", expr)
    }
    stopifnot("`groups` must be provided when `in_form` is 'long'." = !is.null(groups) || in_form == "wide")
    sel(expr, groups, data, id = id, in_form = in_form, return_ids = return_ids)
}

#' @rdname clone_selectors
#' @export
and <- function(x, y) {
    mask_x <- is.na(x)
    mask_y <- is.na(y)
    x[mask_x | mask_y] <- NA
    x
}

#' @rdname clone_selectors
#' @export
or <- function(x, y) {
    dplyr::coalesce(x, y)
}
