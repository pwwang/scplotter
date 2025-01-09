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

