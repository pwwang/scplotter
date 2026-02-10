
#' ClonalStatPlot
#'
#' @description Visualize the statistics of the clones.
#' @param data The product of [scRepertoire::combineTCR], [scRepertoire::combineTCR], or
#' [scRepertoire::combineExpression].
#' @param clones The specific clones to track. This argument must be provided.
#' If multiple character values are provided, they will be treated as clone IDs.
#' If a single character value is provided with parentheses, it will be evaluated as an expression to select the clones.
#' The clones will be selected per subgrouping/facetting/splitting group. For example, if you have
#' `top(3)` will select the top 3 clones in each facetting/splitting group.
#' You can change this behavior by passing the `groups` argument explicitly.
#' For example `top(3, groups = "Sample")` will select the top 3 clones in each sample.
#' For expression, see also [`clone_selectors`](clone_selectors.html).
#' This can also be a named list of expressions, which need to be quoted. Then basic unit for visualization will be the
#' the clone groups defined by the names of the list, instead of single clones.
#' @param top The number of top clones to select. Default is 10.
#' A shortcut for `top(10)` if `clones` is not provided.
#' If `clones` is provided, this will be a limit for the number of clones selected (based on the `orderby` expression).
#' If `clones` is a list, this will be applied to each clone group.
#' @param clone_call How to call the clone - VDJC gene (gene), CDR3 nucleotide (nt),
#'  CDR3 amino acid (aa), VDJC gene + CDR3 nucleotide (strict) or a custom variable
#'  in the data
#' @param chain indicate if both or a specific chain should be used - e.g. "both",
#'  "TRA", "TRG", "IGH", "IGL"
#' @param relabel Whether to relabel the clones. Default is FALSE.
#' The clone ids, especially using CDR3 sequences, can be long and hard to read.
#' If TRUE, the clones will be relabeled as "clone1", "clone2", etc.
#' Only works for visualizations for single clones.
#' @param plot_type The type of plot to use. Default is "bar".
#' Possible values are:
#' * "bar" - bar plot showing the total size of the selected clones in each group.
#' * "box" - box plot showing the distribution of the clone sizes in each group.
#' * "violin" - violin plot showing the distribution of the clone sizes in each group.
#' * "heatmap" - heatmap showing the clone sizes in each group.
#' * "pies" - heatmap with pie charts showing the clone sizes and subgroup compositions in each group. Requires `subgroup_by` to be provided.
#' * "sankey" - sankey plot showing the dynamics of the clones between groups. The clone groups will be defined by the `clones` argument. The flows will be colored by the clone groups.
#' * "alluvial" - same as "sankey".
#' * "trend" - line plot showing the trend of the clone sizes in each group. The clone groups will be defined by the `clones` argument. The lines will be colored by the clone groups.
#' * "col" - same as "col-rel".
#' * "col-rel" - column plot showing the relative size of the clones in each group.
#' * "col-abs" - column plot showing the absolute size of the clones in each group.
#' Note that for "col-rel" and "col-abs", the plot will be faceted by the groups, so "facet_by" is not supported. Please use "split_by" instead if you want to split the plot by another variable.
#' @param group_by The column name in the meta data to group the cells. Default: "Sample"
#' @param groups The groups to include in the plot. Default is NULL.
#'  If NULL, all the groups in `group_by` will be included.
#' @param subgroup_by The column name in the meta data to subgroup the nodes (group nodes on each `x`). Default: NULL.
#' This argument is only supported for "sankey"/"alluvial" plot.
#' If NULL, the nodes will be grouped/colored by the clones
#' @param subgroups The subgroups to include in the plot. Default is NULL.
#' If NULL, all the subgroups in `subgroup_by` will be included.
#' If a vector, the subgroups will be included in the order of the vector for all `groups`.
#' If a list, the subgroups will be used for each `groups`, with `groups` as the names.
#' @param within_subgroup Whether to select the clones within each subgroup.
#' @param subgroups The subgroups to include in the plot. Default is NULL.
#' @param facet_by The column name in the meta data to facet the plots. Default: NULL.
#'  This argument is not supported and will raise an error if provided.
#' @param split_by The column name in the meta data to split the plots. Default: NULL
#' @param orderby An expression to order the clones by. Default is NULL.
#' Note that the clones will be ordered by the value of this expression in descending order.
#' @param y The y-axis variable to use for the plot. Default is NULL.
#' * For `bar` plot, Either "TotalSize" or "Count" can be used, representing the total size (# cells) of the selected clones or the number of selected clones, respectively.
#' @param xlab The x-axis label. Default is NULL.
#' @param ylab The y-axis label. Default is NULL.
#' @param ... Other arguments passed to the specific plot function.
#' * For `bar` plot, see [plotthis::BarPlot()].
#' * For `trend` plot, see [plotthis::TrendPlot()].
#' * For `sankey` plot, see [plotthis::SankeyPlot()].
#' @return A ggplot object or a list if `combine` is FALSE
#' @importFrom rlang parse_expr as_label enexpr
#' @importFrom dplyr %>% summarise arrange desc mutate ungroup across slice_head n
#' @importFrom tidyr pivot_wider
#' @importFrom plotthis SankeyPlot TrendPlot
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
#' # add a fake variable (e.g. cell type from scRNA-seq)
#' data <- lapply(data, function(x) {
#'     x$CellType <- sample(c("CD4", "CD8", "B", "NK"), nrow(x), replace = TRUE)
#'     # x <- x[x$CTaa == "CAVRKTTGTASKLTF_CASSLFGDKGETQYF", , drop = F]
#'     return(x)
#' })
#'
#' # showing the top 10 clones in P17B and P17L
#' ClonalStatPlot(data, group_by = "Sample", groups = c("P17B", "P17L"))
#' # showing the top 10 clones in P17B and P17L, with the clones relabeled
#' ClonalStatPlot(data, group_by = "Sample", groups = c("P17B", "P17L"), relabel = TRUE)
#' # showing the top 2 clones in groups B and L, with subgroups in each group
#' ClonalStatPlot(data, group_by = "Type", subgroup_by = "Sample", top = 2,
#'     subgroups = c("P17B", "P17L", "P18B", "P18L", "P19B","P19L"), relabel = TRUE)
#' # showing selected clones in P17B and P17L
#' ClonalStatPlot(data, group_by = "Sample", groups = c("P17B", "P17L"),
#'     clones = c("CVVSDNTGGFKTIF_CASSVRRERANTGELFF", "NA_CASSVRRERANTGELFF"), relabel = TRUE)
#' # facetting is supported
#' ClonalStatPlot(data, group_by = "Subject", groups = c("P17", "P19"),
#'     facet_by = "Type", relabel = TRUE)
#' # as well as splitting
#' ClonalStatPlot(data, group_by = "Subject", groups = c("P17", "P19"),
#'     split_by = "Type", relabel = TRUE)
#' # showing shared clones between P17B and P17L (top 10 clones that are present in both samples)
#' ClonalStatPlot(data, group_by = "Sample", groups = c("P17B", "P17L"),
#'      clones = "shared(P17B, P17L)", relabel = TRUE, top = 10)
#' # showing shared clones but with a different order
#' ClonalStatPlot(data, group_by = "Sample", groups = c("P17B", "P17L"), top = 10,
#'      clones = "shared(P17B, P17L)", relabel = TRUE, orderby = "P17B")
#' # showing clones larger than 10 in P17L and ordered by the clone size in P17L descendingly
#' ClonalStatPlot(data, group_by = "Sample", groups = c("P17B", "P17L"),
#'      clones = "sel(P17L > 10)", relabel = TRUE, top = 5, orderby = "P17L")
#' # using trend plot
#' ClonalStatPlot(data, group_by = "Sample", groups = c("P17B", "P17L"),
#'     clones = sel(P17L > 10 & P17B > 0), relabel = TRUE, orderby = "P17L",
#'     plot_type = "trend")
#' # using heatmap
#' ClonalStatPlot(data, group_by = "Sample", groups = c("P17B", "P17L"),
#'     clones = sel(P17L > 10 & P17B > 0), relabel = TRUE, orderby = "P17L",
#'     plot_type = "heatmap")
#' # using heatmap with subgroups
#' ClonalStatPlot(data, group_by = "Sample", groups = c("P17B", "P17L"),
#'     clones = list(
#'         ExpandedClonesInP17L = "sel(P17L > 20)",
#'         ExpandedClonesInP17B = "sel(P17B > 20)"
#'     ), subgroup_by = "CellType", pie_size = sqrt,
#'     plot_type = "pies", show_row_names = TRUE, show_column_names = TRUE)
#' # using clone groups and showing dynamics using sankey plot
#' ClonalStatPlot(data, group_by = "Sample", groups = c("P17B", "P17L"),
#'     clones = list(
#'       "Hyper-expanded clones in P17B" = "sel(P17B > 10)",
#'       "Hyper-expanded clones in P17L" = "sel(P17L > 10)"
#'     ), plot_type = "sankey")
#' # col-rel/col-abs
#' ClonalStatPlot(data, plot_type = "col-rel")
#' ClonalStatPlot(data, plot_type = "col-abs", facet_scale = "free")
#' ClonalStatPlot(data, plot_type = "col", groups = c("P17B", "P17L"),
#'     top = 20, facet_ncol = 1, legend.position = "right")
#' ClonalStatPlot(data, plot_type = "col", groups = c("P17B", "P17L"),
#'     top = 20, facet_ncol = 1, legend.position = "right",
#'     relabel = TRUE, fill_by = "CloneGroups", fill_name = "Clones")
#' }
ClonalStatPlot <- function(
    data, clones = NULL, top = 10, orderby = NULL, clone_call = "aa", chain = "both",
    plot_type = c("bar", "box", "violin", "heatmap", "pies", "sankey", "alluvial", "trend", "col", "col-rel", "col-abs"),
    group_by = "Sample", groups = NULL, subgroup_by = NULL, subgroups = NULL,
    within_subgroup = match.arg(plot_type) != "pies", relabel = FALSE,
    facet_by = NULL, split_by = NULL, y = NULL, xlab = NULL, ylab = NULL, ...
) {
    plot_type <- match.arg(plot_type)
    if (plot_type == "alluvial") plot_type <- "sankey"
    if (plot_type == "circos") plot_type <- "chord"
    if (plot_type == "col") plot_type <- "col-rel"
    stopifnot("Only a single group_by is supported for 'ClonalStatPlot'" = length(unique(group_by)) <= 1)
    stopifnot("'subgroup_by' is not supported for 'ClonalStatPlot' with plot_type = 'trend'" = is.null(subgroup_by) || plot_type != "trend")
    stopifnot("'subgroup_by' is not supported for 'ClonalStatPlot' with plot_type = 'col-rel' or 'col-abs'" = is.null(subgroup_by) || !(plot_type %in% c("col-rel", "col-abs")))

    if (!is.null(groups)) {
        data <- screp_subset(data, paste0('`', group_by, '`', ' %in% c(', paste0('"', groups, '"', collapse = ', '), ')'))
    }

    all_groupings <- unique(c(group_by, subgroup_by, facet_by, split_by))
    data <- clonal_size_data(data, clone_call, chain, all_groupings)

    # data$fraction <- NULL
    if (!is.null(group_by)) {
        groups <- groups %||% unique(data[[group_by]])
        nonexist_groups <- setdiff(groups, unique(data[[group_by]]))
        if (length(nonexist_groups) > 0) {
            stop(paste("The following groups do not exist in the data:", paste(nonexist_groups, collapse = ", ")))
        }
    }
    if (length(groups) < 2 && !plot_type %in% c("col-rel", "col-abs")) {
        stop("At least 2 groups are required for ClonalStatPlot")
    }
    if (identical(plot_type, "chord") && length(groups) > 2) {
        stop("'chord'/'circos' plot only supports up to 2 groups. Please use 'sankey' plot instead.")
    }

    topn <- top
    top <- getFromNamespace("top", "scplotter")
    orderby <- tryCatch({
        orderby
    }, error = function(e) {
        # if orderby is not a valid expression, use the string representation
        as_label(enexpr(orderby))
    })
    if (plot_type %in% c("col-rel", "col-abs")) {
        orderby <- orderby %||% "desc(count)"
    } else if (!is.null(groups)) {
        orderby <- orderby %||% paste0("desc(", paste0("`", groups, "`", collapse = "+"), ")")
    }

    clones <- tryCatch({
        if (!is.list(clones)) {
            clones <- list(.selected.clones = clones)
        } else {
            clones
        }
    }, error = function(e) {
        # if orderby is not a valid expression, use the string representation
        list(.selected.clones = as_label(enexpr(clones)))
    })

    by_clones <- identical(names(clones), ".selected.clones")
    if (by_clones && identical(clones$.selected.clones, "list(...)")) {
        stop("[CloneStatPlot] When 'clones' is provided as a list of expressions; the expressions must be quoted.")
    }
    if (!plot_type %in% c("col-rel", "col-abs")) {
        data <- data %>%
            pivot_wider(names_from = group_by, values_from = "count", values_fill = 0) %>%
            dplyr::group_by(!!!syms(unique(c("CloneID", subgroup_by, facet_by, split_by)))) %>%
            summarise(across(groups, sum), .groups = "drop") %>%
            # keep the grouping columns for the selectors
            dplyr::group_by(!!!syms(unique(c(subgroup_by, facet_by, split_by))))
    }
    data <- data %>% arrange(!!parse_expr(orderby))

    fulldata <- data
    # CloneID                                    Sample CellType count
    # <fct>                                      <chr>  <chr>    <dbl>
    # 1 CAAADTGTASKLTF_CASSYSPQGGYEQYF           P17B   CD8          1
    # 2 CAAADTGTASKLTF_CASSYSPQGGYEQYF           P17L   B            1

    if (!is.null(subgroup_by) && !within_subgroup) {
        # if subgroup_by is provided, but selection is not within subgroup,
        # we need to merge the subgroup_by values
        data <- data %>%
            dplyr::group_by(!!!syms(unique(c("CloneID", facet_by, split_by)))) %>%
            summarise(across(groups, sum), .groups = "drop") %>%
            dplyr::group_by(!!!syms(unique(c(facet_by, split_by)))) %>%
            arrange(!!parse_expr(orderby))
    }

    selected_data <- NULL
    for (nc in names(clones)) {
        clones_expr <- clones[[nc]]
        if (is.null(clones_expr)) {
            selected <- data
        } else if (length(clones_expr) == 1 && grepl("(", clones_expr, fixed = TRUE) && grepl(")", clones_expr, fixed = TRUE)) {
            selected <- eval(parse(text = clones_expr))
        } else {
            selected <- filter(data, !!sym("CloneID") %in% clones_expr)
        }
        if (!is.null(topn)) {
            groupings <- unique(intersect(all_groupings, colnames(selected)))
            if (length(groupings) == 0) {
                selected <- slice_head(selected, n = topn)
            } else {
                selected <- selected %>%
                    group_by(!!!syms(groupings)) %>%
                    slice_head(n = topn) %>%
                    ungroup()
            }
        }
        if (nrow(selected) == 0) {
            stop("No clones selected in the data with clone selector: ", ifelse(by_clones, "default", nc), ".")
        }
        selected$CloneGroups <- if (by_clones) selected$CloneID else nc
        if (!is.null(subgroup_by) && !within_subgroup) {
            # attach the subgroup_by to the selected data
            # a clone/clone group can have multiple subgroup_by values
            # split the values in each group for the subgroup_by
            #
            # For example (data):
            # CloneID                          P17B  P17L  CloneGroups
            # <fct>                            <dbl> <dbl> <fct>
            # CAVSNTGNQFYF_CASSLGLAGAGVDTQYF   1     11    ExpandedClonesInP17L
            #
            # To (fulldata, subgroup_by = "CellType"):
            # CloneID                          CellType  P17B  P17L  CloneGroups
            # <fct>                            <fct>    <dbl> <dbl> <fct>
            # CAVSNTGNQFYF_CASSLGLAGAGVDTQYF   CD4      1     1     ExpandedClonesInP17L
            # CAVSNTGNQFYF_CASSLGLAGAGVDTQYF   CD8      0     5     ExpandedClonesInP17L
            # CAVSNTGNQFYF_CASSLGLAGAGVDTQYF   NK       0     5     ExpandedClonesInP17L

            selected <- dplyr::left_join(
                fulldata,
                selected[, unique(c("CloneID", facet_by, split_by, "CloneGroups")), drop = FALSE],
                by = unique(c("CloneID", facet_by, split_by))
            )
            selected <- selected[!is.na(selected$CloneGroups), , drop = FALSE]
        }
        if (is.null(selected_data)) {
            selected_data <- selected
        } else {
            selected_data <- rbind(selected_data, selected)
        }
    }

    if (!is.null(selected_data)) {
        data <- ungroup(selected_data)
        rm(selected_data)
    }
    data$CloneGroups <- factor(data$CloneGroups, levels = unique(data$CloneGroups))

    if (by_clones) {
        data$CloneGroups <- factor(data$CloneGroups, levels = unique(data$CloneGroups))
        if (relabel) {
            data$CloneGroups <- paste0("clone", as.numeric(data$CloneGroups))
        }
    }

    clone_groups_name <- if (by_clones) "Clones" else "Clone Groups"

    if (plot_type %in% c("col-rel", "col-abs")) {
        if (!is.null(facet_by)) {
            stop("'facet_by' is not supported for 'col-rel'/'col-abs' plot. Please use 'split_by' instead.")
        }
        args <- rlang::dots_list(...)
        args$data <- data
        args$x <- "CloneGroups"
        args$y <- ifelse(plot_type == "col-rel", "fraction", "count")
        args$split_by <- split_by
        args$facet_by <- group_by
        args$facet_scales <- args$facet_scale %||% "free_x"
        args$fill_by <- args$fill_by %||% FALSE
        args$xlab <- xlab %||% clone_groups_name
        args$ylab <- ylab %||% ifelse(plot_type == "col-rel", "Relative Abundance", "Size")
        args$legend.position <- args$legend.position %||% "none"
        args$aspect.ratio <- args$aspect.ratio %||% 0.6
        args$theme_args <- args$theme_args %||% list()
        args$theme_args$panel.grid.major.x <- args$theme_args$panel.grid.major.x %||% ggplot2::element_blank()
        args$theme_args$axis.text.x <- args$theme_args$axis.text.x %||% ggplot2::element_blank()
        axis.text.x <- args$theme_args$axis.text.x
        args$theme_args$axis.text.x <- NULL
        p <- do.call(BarPlot, args) & ggplot2::theme(axis.text.x = axis.text.x)
        attr(p, "height") <- 0.6 * attr(p, "width")
        p
    } else if (identical(plot_type, "bar")) {
        data <- data %>%
            pivot_longer(cols = groups, names_to = group_by, values_to = "Size") %>%
            dplyr::filter(!!sym("Size") > 0) %>%
            dplyr::group_by(!!!syms(unique(c("CloneGroups", group_by, subgroup_by, facet_by, split_by)))) %>%
            summarise(TotalSize = sum(!!sym("Size")), Count = n(), .groups = "drop")

        BarPlot(data, x = group_by, y = y %||% "TotalSize", group_by = "CloneGroups", group_name = clone_groups_name,
            facet_by = facet_by, split_by = split_by, xlab = xlab, ylab = ylab, ...)
    } else if (identical(plot_type, "box")) {
        data <- data %>% pivot_longer(cols = groups, names_to = group_by, values_to = "Size") %>%
            dplyr::filter(!!sym("Size") > 0)

        BoxPlot(data, x = group_by, y = "Size", group_by = "CloneGroups", group_name = clone_groups_name,
            facet_by = facet_by, split_by = split_by, xlab = xlab, ylab = ylab, ...)
    } else if (identical(plot_type, "violin")) {
        data <- data %>% pivot_longer(cols = groups, names_to = group_by, values_to = "Size") %>%
            dplyr::filter(!!sym("Size") > 0)

        ViolinPlot(data, x = group_by, y = "Size", group_by = "CloneGroups", group_name = clone_groups_name,
           split_by = split_by, xlab = xlab, ylab = ylab, ...)
    } else if (identical(plot_type, "heatmap")) {
        if (!is.null(facet_by)) {
            stop("'facet_by' is not supported for 'heatmap' plot. Please use 'split_by' instead.")
        }
        data[[clone_groups_name]] <- data$CloneGroups
        data$CloneGroups <- NULL
        Heatmap(data, rows_by = groups, columns_by = clone_groups_name, rows_name = group_by,
            split_by = split_by, name = "Clone Size", ...)
    } else if (identical(plot_type, "pies")) {
        if (!is.null(facet_by)) {
            stop("'facet_by' is not supported for 'pies' plot. Please use 'split_by' instead.")
        }
        if (is.null(subgroup_by)) {
            stop("'subgroup_by' is required for 'pies' plot. Please provide it.")
        }
        data[[clone_groups_name]] <- data$CloneGroups
        data$CloneGroups <- NULL
        args <- rlang::dots_list(...)
        args$data <- data
        args$rows_by <- groups
        args$columns_by <- clone_groups_name
        args$rows_name <- group_by
        args$split_by <- split_by
        args$cell_type <- "pie"
        args$pie_group_by <- subgroup_by
        args$pie_values <- "sum"
        args$pie_size <- args$pie_size %||% sqrt
        args$pie_size_name <- args$pie_size_name %||% "Size"
        args$show_row_names <- args$show_row_names %||% TRUE
        args$show_column_names <- args$show_column_names %||% TRUE
        args$add_reticle <- args$add_reticle %||% TRUE
        do.call(Heatmap, args)
    } else if (identical(plot_type, "sankey")) {
        SankeyPlot(data, x = groups, links_name = clone_groups_name,
            links_fill_by = "CloneGroups", flow = TRUE, xlab = xlab %||% group_by, ylab = ylab,
            facet_by = facet_by, split_by = split_by, ...)
    } else {
        data <- data %>% pivot_longer(cols = groups, names_to = group_by, values_to = "Size") %>%
            dplyr::filter(!!sym("Size") > 0) %>%
            dplyr::group_by(!!!syms(unique(c("CloneGroups", group_by, subgroup_by, facet_by, split_by)))) %>%
            summarise(TotalSize = sum(!!sym("Size")), .groups = "drop")
        TrendPlot(data, x = group_by, y = "TotalSize", group_by = "CloneGroups", group_name = clone_groups_name,
            facet_by = facet_by, split_by = split_by, xlab = xlab, ylab = ylab, ...)
    }
}

#' ClonalDynamicsPlot
#'
#' @description This function is deprecated. Please use [scplotter::ClonalStatPlot()] instead.
#' @param data The product of [scRepertoire::combineTCR], [scRepertoire::combineTCR], or
#' [scRepertoire::combineExpression].
#' @param ... Other arguments.
#' @export
ClonalDynamicsPlot <- function(data, ...) {
    stop("ClonalDynamicsPlot is deprecated. Please use ClonalStatPlot with plot_type 'sankey' instead.")
}
