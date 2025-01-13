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
        all_gvalues <- unique(data@meta.data$.group)
    } else {
        # clonalScatter only returns data for each sample
        # need to re-organize the data to get the data for each group
        all_gvalues <- unique(sapply(data, function(x) x$.group[1]))
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

    gv_pairs <- as.list(as.data.frame(combn(all_gvalues, 2, simplify = TRUE)))
    do.call(rbind, lapply(gv_pairs, function(gv) {
        d <- clonalScatter(data,
            cloneCall = clone_call, chain = chain,
            x.axis = gv[1], y.axis = gv[2], exportTable = TRUE
        )
        d$class <- NULL
        d$sum <- NULL
        d$size <- NULL
        names(d)[2:3] <- paste(names(d)[2:3], "count", sep = ".")
        d <- pivot_longer(
            d,
            cols = -"Var1",
            names_to = c(".group", ".value"),
            names_sep = "\\."
        )
        d
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
#' @importFrom scRepertoire addVariable
#' @keywords internal
merge_clonal_groupings <- function(data, groupings, sep = " // ") {
    if (inherits(data, "Seurat")) {
        if (!"Sample" %in% colnames(data@meta.data)) {
            warning("The 'Sample' column is not found in the meta data, 'orig.indent' will be used instead.")
            data$Sample <- data$orig.ident
        }
        samples <- unique(data$Sample)
        # combine Sample and group_by so that the count/frequency is calculated for each
        # combined group
        data@meta.data <- unite(data@meta.data, ".group", !!!syms(groupings), sep = sep, remove = FALSE)
    } else {
        samples <- names(data)
        data <- addVariable(data, variable.name = "Sample", variables = samples)
        # combine Sample and group_by so that the count/frequency is calculated for each
        # combined group
        vs <- sapply(samples, function(s) {
            paste0(sapply(groupings, function(group) {
                data[[s]][, group, drop = TRUE][1]
            }), collapse = sep)
        })
        data <- addVariable(data, variable.name = ".group", variables = vs)
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
#'    samples = c("P17B", "P17L", "P18B", "P18L", "P19B","P19L", "P20B", "P20L"))
#'
#' head(scplotter:::screp_subset(screp, "nchar(CTaa) < 20")[[1]])
#' names(scplotter:::screp_subset(screp, "Sample %in% c('P17B', 'P17L')"))
#' }
screp_subset <- function(screp, subset) {
    if (inherits(screp, "Seurat")) {
        eval(parse(text = paste('base::subset(screp, subset = "', subset, '")')))
    } else {
        screp <- sapply(names(screp), function(x) {
            y <- screp[[x]]
            y$Sample <- x
            filter(y, !!parse_expr(subset))
        }, simplify = FALSE, USE.NAMES = TRUE)
        screp[sapply(screp, nrow) > 0]
    }
}

#' Helper functions to select clones based on various criteria
#'
#' @description These helper functions allow for the selection of clones based on various criteria such as size, group comparison, and existence in specific groups.
#'
#' @param n The number of top clones to select or the threshold size.
#' @param expr The expression (in characters) to filter the clones
#' (e.g. "group1 > group2" to select clones where group1 is larger than group2).
#' @param group1 The first group to compare.
#' @param group2 The second group to compare.
#' @param ... More groups to compare.
#' @param include_eq Whether to include equal-sized clones.
#' @param groups The column names in the meta data to group the cells.
#' By default, it is assumed `facet_by` and `split_by` to be in the parent frame.
#' @param data The data frame containing clone information. Default is NULL. If NULL, it will get data from parent.frame.
#' A typical `data` should have a column named `CloneID` and other columns for the groupings.
#' Under each grouping column, the value should be the size of the clone.
#' By default, the data is assumed to be in the parent frame.
#' @return A vector of selected clones.
#' @importFrom rlang parse_expr syms sym
#' @importFrom dplyr group_by summarise filter reframe pull
#' @export
#' @rdname clone_selectors
#' @examples
#' data <- data.frame(
#'    CloneID = 1:10,
#'    group1 = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
#'    group2 = c(7, 3, 8, 2, 1, 5, 9, 4, 6, 0),
#'    groups = c("A", "A", "A", "A", "B", "B", "B", "B", "B", "B")
#' )
#' data <- data[order(data$group1 + data$group2, decreasing = TRUE), ]
#' scplotter:::top(3)
#' scplotter:::top(3, groups = "groups")
#' scplotter:::select(group1 == 0 | group2 == 0)
#' scplotter:::uniq(group1, group2)
#' scplotter:::shared(group1, group2)
#' scplotter:::larger(group1, group2)
#' scplotter:::smaller(group1, group2)
#' scplotter:::smaller(group1, group2, include_eq = TRUE)
#' scplotter:::smaller(group1, group2, shared = TRUE)
#' scplotter:::eq(group1, group2)
top <- function(n, groups = NULL, data = NULL) {
    data <- data %||% parent.frame()$data
    groups <- groups %||% c(parent.frame()$split_by, parent.frame()$facet_by)
    if (is.null(groups)) {
        return(data$CloneID[1:n])
    }
    data %>%
        group_by(!!!syms(groups)) %>%
        reframe(CloneID = (!!sym("CloneID"))[1:n])
}

#' @rdname clone_selectors
#' @keywords internal
select <- function(expr, groups = NULL, data = NULL) {
    data <- data %||% parent.frame()$data
    groups <- groups %||% c(parent.frame()$split_by, parent.frame()$facet_by)
    is_char <- tryCatch({
        is.character(expr)
    }, error = function(e) {
        FALSE
    })
    if (!is_char) {
        expr <- deparse(substitute(expr))
    }
    if (is.null(groups) || isFALSE(groups)) {
        return(data %>% filter(!!parse_expr(expr)) %>% pull("CloneID"))
    }
    data %>%
        group_by(!!!syms(groups)) %>%
        reframe(CloneID = select(expr, FALSE, data), .groups = "drop")
}

#' @rdname clone_selectors
#' @keywords internal
uniq <- function(group1, group2, ..., groups = NULL, data = NULL) {
    data <- data %||% parent.frame()$data
    groups <- groups %||% c(parent.frame()$split_by, parent.frame()$facet_by)
    group1 <- as.character(substitute(group1))
    group2 <- as.character(substitute(group2))
    other_groups <- as.character(substitute(list(...)))[-1]
    expr <- paste0("`", group1, "` > 0 & `", group2, "` == 0")
    for (g in other_groups) {
        expr <- paste0(expr, " & `", g, "` == 0")
    }
    return(select(expr, groups, data))
}

#' @rdname clone_selectors
#' @keywords internal
shared <- function(group1, group2, ..., groups = NULL, data = NULL) {
    data <- data %||% parent.frame()$data
    groups <- groups %||% c(parent.frame()$split_by, parent.frame()$facet_by)
    group1 <- as.character(substitute(group1))
    group2 <- as.character(substitute(group2))
    other_groups <- as.character(substitute(list(...)))[-1]
    expr <- paste0("`", group1, "` > 0 & `", group2, "` > 0")
    for (g in other_groups) {
        expr <- paste0(expr, " & `", g, "` > 0")
    }
    return(select(expr, groups, data))
}

#' @rdname clone_selectors
#' @keywords internal
larger <- function(group1, group2, include_eq = FALSE, shared = FALSE, groups = NULL, data = NULL) {
    data <- data %||% parent.frame()$data
    groups <- groups %||% c(parent.frame()$split_by, parent.frame()$facet_by)
    group1 <- as.character(substitute(group1))
    group2 <- as.character(substitute(group2))
    expr <- ifelse(include_eq, ">=", ">")
    expr <- paste0("`", group1, "` ", expr, " `", group2, "`")
    if (shared) {
        expr <- paste0("`", group1, "` > 0 & `", group2, "` > 0 & ", expr)
    }
    return(select(expr, groups, data))
}

#' @rdname clone_selectors
#' @keywords internal
smaller <- function(group1, group2, include_eq = FALSE, shared = FALSE, groups = NULL, data = NULL) {
    data <- data %||% parent.frame()$data
    groups <- groups %||% c(parent.frame()$split_by, parent.frame()$facet_by)
    group1 <- as.character(substitute(group1))
    group2 <- as.character(substitute(group2))
    expr <- ifelse(include_eq, "<=", "<")
    expr <- paste0("`", group1, "` ", expr, " `", group2, "`")
    if (shared) {
        expr <- paste0("`", group1, "` > 0 & `", group2, "` > 0 & ", expr)
    }
    return(select(expr, groups, data))
}

#' @rdname clone_selectors
#' @keywords internal
eq <- function(group1, group2, groups = NULL, shared = FALSE, data = NULL) {
    data <- data %||% parent.frame()$data
    groups <- groups %||% c(parent.frame()$split_by, parent.frame()$facet_by)
    group1 <- as.character(substitute(group1))
    group2 <- as.character(substitute(group2))
    expr <- paste0("`", group1, "` == `", group2, "`")
    if (shared) {
        expr <- paste0("`", group1, "` > 0 & `", group2, "` > 0 & ", expr)
    }
    return(select(expr, groups, data))
}
