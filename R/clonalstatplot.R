
#' ClonalStatPlot
#'
#' @description Visualize the statistics of the clones.
#' @param data The product of [scRepertoire::combineTCR], [scRepertoire::combineTCR], or
#' [scRepertoire::combineExpression].
#' @param clones The specific clones to track. This argument must be provided.
#' If multiple character values are provided, they will be treated as clone IDs.
#' If a single character value is provided with parentheses, it will be evaluated as an expression to select the clones.
#' The clones will be selected per facetting/splitting group. For example, if you have
#' `top(3)` will select the top 3 clones in each facetting/splitting group.
#' You can change this behavior by passing the `group_by` argument explicitly.
#' For example `top(3, group_by = "Sample")` will select the top 3 clones in each sample.
#' For expression, see also [`clone_selectors`](clone_selectors.html).
#' This can also be a named list of expressions, which need to be quoted. Then basic unit for visualization will be the
#' the clone groups defined by the names of the list, instead of single clones.
#' @param clone_call How to call the clone - VDJC gene (gene), CDR3 nucleotide (nt),
#'  CDR3 amino acid (aa), VDJC gene + CDR3 nucleotide (strict) or a custom variable
#'  in the data
#' @param chain indicate if both or a specific chain should be used - e.g. "both",
#'  "TRA", "TRG", "IGH", "IGL"
#' @param values_by The variable to use for the values of the clones.
#' Default is "count", which represents the number of cells in each clone.
#' "fraction" can also be used to represent the fraction of cells in each clone out of the total cells in the group.
#' "n" can be used to represent the number of cells in each clone, same as "count".
#' @param relabel Whether to relabel the clones. Default is FALSE.
#' The clone ids, especially using CDR3 sequences, can be long and hard to read.
#' If TRUE, the clones will be relabeled as "clone1", "clone2", etc. (ordered by the descending clone sizes)
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
#' * "col" - column plot showing the size of the clones in each group.
#' Note that for "col", the plot will be faceted by the groups, so "facet_by" is not supported. Please use "split_by" instead if you want to split the plot by another variable.
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
#' @importFrom dplyr %>% summarise arrange desc mutate ungroup across slice_head n filter
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
#'     return(x)
#' })
#' # showing the top 10 clones (by default)
#' ClonalStatPlot(data, group_by = "Sample", title = "Top 10 clones")
#' # showing the top 10 clones in P17B and in P17L, with the clones relabeled
#' ClonalStatPlot(data, clones = "top(10, group_by = 'Sample')", group_by = "Sample",
#'     groups = c("P17B", "P17L"), relabel = TRUE, values_by = "fraction",
#'     title = "Top 10 clones in P17B and in P17L (relabelled)")
#' # showing the top 10 clones in each sample using violin plots
#' ClonalStatPlot(data, group_by = "Sample",
#'     plot_type = "violin", clones = "top(10, group_by = 'Sample')",
#'     subgroup_by = "CellType", subgroups = c("CD4", "CD8"), add_box = TRUE,
#'     comparison = TRUE, title = "Violin plots showing top 10 clones in each sample")
#' # showing selected clones in P17B and P17L
#' ClonalStatPlot(data, group_by = "Sample", groups = c("P17B", "P17L"),
#'     clones = c("CVVSDNTGGFKTIF_CASSVRRERANTGELFF", "NA_CASSVRRERANTGELFF"),
#'     title = "Selected clones in P17B and P17L")
#' # facetting is supported, note that selection of clones is done within each facet
#' ClonalStatPlot(data, group_by = "Subject", groups = c("P17", "P19"),
#'     facet_by = "Type", relabel = TRUE,
#'     title = "Top 10 clones in Type B and L for P17 and P19")
#' # as well as splitting
#' ClonalStatPlot(data, group_by = "Subject", groups = c("P17", "P19"),
#'     split_by = "Type", relabel = TRUE)
#' # showing top 10 shared clones between P17B and P17L
#' ClonalStatPlot(data, group_by = "Sample", groups = c("P17B", "P17L"),
#'     clones = "shared(P17B, P17L, group_by = 'Sample', top = 10)", relabel = TRUE,
#'     title = "Shared clones between P17B and P17L")
#' # showing clones larger than 10 in P17L and ordered by the clone size in P17L descendingly
#' ClonalStatPlot(data, group_by = "Sample", groups = c("P17B", "P17L"),
#'      clones = "sel(P17B > 10, group_by = 'Sample', top = 5, order = desc(P17B))",
#'      relabel = TRUE, position = "stack", title = "Top 5 clones larger than 10 in P17B")
#' # using trend plot
#' ClonalStatPlot(data, group_by = "Sample", groups = c("P17B", "P17L"),
#'     clones = "sel(P17L > 10 & P17B > 0, group_by = 'Sample')", relabel = TRUE,
#'     plot_type = "trend", title = "Clones larger than 10 in P17L and existing in P17B")
#' # using heatmap
#' ClonalStatPlot(data, group_by = "Sample", groups = c("P17B", "P17L"),
#'     clones = "sel(P17L > 10 & P17B > 0, group_by = 'Sample')", relabel = TRUE,
#'     plot_type = "heatmap", show_row_names = TRUE, show_column_names = TRUE,
#'     title = "Clones larger than 10 in P17L and existing in P17B (heatmap)")
#' # using heatmap with subgroups for groups of clones
#' ClonalStatPlot(data, group_by = "Sample", groups = c("P17B", "P17L"),
#'     clones = list(
#'         ExpandedClonesInP17L = "sel(P17L > 20, group_by = 'Sample')",
#'         ExpandedClonesInP17B = "sel(P17B > 20, group_by = 'Sample')"
#'     ), subgroup_by = "CellType", pie_size = sqrt,
#'     plot_type = "pies", show_row_names = TRUE, show_column_names = TRUE,
#'     title = "Clones larger than 20 in P17L and P17B (pies with subgroups by CellType)")
#' # using clone groups and showing dynamics using sankey plot
#' ClonalStatPlot(data, group_by = "Sample", groups = c("P17B", "P17L"),
#'     clones = list(
#'       "Hyper-expanded clones in P17B" = "sel(P17B > 10, group_by = 'Sample')",
#'       "Hyper-expanded clones in P17L" = "sel(P17L > 10, group_by = 'Sample')"
#'     ), plot_type = "sankey", title = "Hyper-expanded clones in P17B and P17L")
#' # col plot
#' ClonalStatPlot(data, clones = "top(5, group_by = 'Sample')", plot_type = "col",
#'     title = "Top 5 clones in each sample (col plot)")
#' ClonalStatPlot(data, clones = "top(5, group_by = 'Sample')", plot_type = "col",
#'     values_by = "fraction", facet_scale = "free",
#'     title = "Top 5 clones in each sample (col plot, showing fraction)")
#' ClonalStatPlot(data, plot_type = "col", groups = c("P17B", "P17L"),
#'     facet_ncol = 1, legend.position = "right",
#'     relabel = TRUE, fill_by = ".Clones", fill_name = "Clones")
#' # showing top 10 shared clones between P17B and P17L
#' ClonalStatPlot(data, group_by = "Sample", groups = c("P17B", "P17L"),
#'     clones = "shared(P17B, P17L, group_by = 'Sample', top = 10)", relabel = TRUE,
#'     title = "Shared clones between P17B and P17L")
#' # showing clones larger than 10 in P17L and ordered by the clone size in P17L descendingly
#' ClonalStatPlot(data, group_by = "Sample", groups = c("P17B", "P17L"),
#'      clones = "sel(P17B > 10, group_by = 'Sample', top = 5, order = desc(P17B))",
#'      relabel = TRUE, position = "stack", title = "Top 5 clones larger than 10 in P17B")
#' # using trend plot
#' ClonalStatPlot(data, group_by = "Sample", groups = c("P17B", "P17L"),
#'     clones = "sel(P17L > 10 & P17B > 0, group_by = 'Sample')", relabel = TRUE,
#'     plot_type = "trend", title = "Clones larger than 10 in P17L and existing in P17B")
#' # using heatmap
#' ClonalStatPlot(data, group_by = "Sample", groups = c("P17B", "P17L"),
#'     clones = "sel(P17L > 10 & P17B > 0, group_by = 'Sample')", relabel = TRUE,
#'     plot_type = "heatmap", show_row_names = TRUE, show_column_names = TRUE,
#'     title = "Clones larger than 10 in P17L and existing in P17B (heatmap)")
#' # chord plot
#' ClonalStatPlot(data, group_by = "Sample", groups = c("P17B", "P17L"),
#'     clones = "sel(P17L > 10 & P17B > 0, group_by = 'Sample')",
#'     plot_type = "chord", labels_rot = TRUE,
#'     title = "Clones larger than 10 in P17L and existing in P17B (chord plot)")
#' # using heatmap with subgroups for groups of clones
#' ClonalStatPlot(data, group_by = "Sample", groups = c("P17B", "P17L"),
#'     clones = list(
#'         ExpandedClonesInP17L = "sel(P17L > 20, group_by = 'Sample')",
#'         ExpandedClonesInP17B = "sel(P17B > 20, group_by = 'Sample')"
#'     ), subgroup_by = "CellType", pie_size = sqrt,
#'     plot_type = "pies", show_row_names = TRUE, show_column_names = TRUE,
#'     title = "Clones larger than 20 in P17L and P17B (pies with subgroups by CellType)")
#' # using clone groups and showing dynamics using sankey plot
#' ClonalStatPlot(data, group_by = "Sample", groups = c("P17B", "P17L"),
#'     clones = list(
#'       "Hyper-expanded clones in P17B" = "sel(P17B > 10, group_by = 'Sample')",
#'       "Hyper-expanded clones in P17L" = "sel(P17L > 10, group_by = 'Sample')"
#'     ), plot_type = "sankey", title = "Hyper-expanded clones in P17B and P17L")
#' # col plot
#' ClonalStatPlot(data, clones = "top(5, group_by = 'Sample')", plot_type = "col",
#'     title = "Top 5 clones in each sample (col plot)")
#' ClonalStatPlot(data, clones = "top(5, group_by = 'Sample')", plot_type = "col",
#'     values_by = "fraction", facet_scale = "free",
#'     title = "Top 5 clones in each sample (col plot, showing fraction)")
#' ClonalStatPlot(data, plot_type = "col", groups = c("P17B", "P17L"),
#'     facet_ncol = 1, legend.position = "right",
#'     relabel = TRUE, fill_by = ".Clones", fill_name = "Clones")
#' }
ClonalStatPlot <- function(
    data, clones = "top(10)", clone_call = "aa", chain = "both", values_by = c("count", "fraction", "n"),
    plot_type = c("bar", "box", "violin", "heatmap", "pies", "circos", "chord", "sankey", "alluvial", "trend", "col"),
    group_by = "Sample", groups = NULL, subgroup_by = NULL, subgroups = NULL,
    within_subgroup = match.arg(plot_type) != "pies", relabel = plot_type %in% c("col", "chord", "circos"),
    facet_by = NULL, split_by = NULL, y = NULL, xlab = NULL, ylab = NULL, ...
) {
    plot_type <- match.arg(plot_type)
    if (plot_type == "alluvial") plot_type <- "sankey"
    if (plot_type == "circos") plot_type <- "chord"
    values_by <- match.arg(values_by)
    stopifnot("Only a single group_by is supported for 'ClonalStatPlot'" = length(unique(group_by)) <= 1)
    stopifnot("'subgroup_by' is not supported for 'ClonalStatPlot' with plot_type = 'bar'" = is.null(subgroup_by) || plot_type != "bar")
    stopifnot("'subgroup_by' is not supported for 'ClonalStatPlot' with plot_type = 'trend'" = is.null(subgroup_by) || plot_type != "trend")
    stopifnot("'subgroup_by' is not supported for 'ClonalStatPlot' with plot_type = 'col'" = is.null(subgroup_by) || plot_type != "col")

    all_groupings <- unique(c(group_by, subgroup_by, facet_by, split_by))
    data <- clonal_size_data(data, clone_call, chain, all_groupings)
    # Selectors need CTaa in tidy environment
    data$CTaa <- data$CloneID
    data$CloneID <- NULL
    all_groupings <- setdiff(all_groupings, group_by)
    if (!within_subgroup) {
        all_groupings <- setdiff(all_groupings, subgroup_by)
    }
    data <- dplyr::group_by(data, !!!syms(all_groupings))

    # data$fraction <- NULL
    if (!is.null(group_by)) {
        check_columns <- utils::getFromNamespace("check_columns", "plotthis")
        group_by <- check_columns(data, group_by, force_factor = TRUE)
        groups <- groups %||% levels(data[[group_by]])
        nonexist_groups <- setdiff(groups, levels(data[[group_by]]))
        if (length(nonexist_groups) > 0) {
            stop(paste("The following groups do not exist in the data:", paste(nonexist_groups, collapse = ", ")))
        }
        data <- data %>% filter(!!sym(group_by) %in% groups)
    }
    if (!is.null(subgroup_by)) {
        check_columns <- utils::getFromNamespace("check_columns", "plotthis")
        subgroup_by <- check_columns(data, subgroup_by, force_factor = TRUE)
        subgroups <- subgroups %||% levels(data[[subgroup_by]])
        nonexist_subgroups <- setdiff(subgroups, levels(data[[subgroup_by]]))
        if (length(nonexist_subgroups) > 0) {
            stop(paste("The following subgroups do not exist in the data:", paste(nonexist_subgroups, collapse = ", ")))
        }
        data <- data %>% filter(!!sym(subgroup_by) %in% subgroups)
    }
    if (length(groups) < 2 && !plot_type %in% c("bar", "col")) {
        stop("At least 2 groups are required for ClonalStatPlot")
    }
    if (identical(plot_type, "chord") && length(groups) > 2) {
        stop("'chord'/'circos' plot only supports up to 2 groups. Please use 'sankey' plot instead.")
    }
    # select clones
    if (!is.list(clones)) {
        clones <- list(.Clones = clones)
    }
    by_clones <- TRUE
    cg_data <- NULL
    if (!identical(names(clones), ".Clones")) {
        # clone groups
        orig_data <- data
        by_clones <- FALSE
    }

    for (clname in names(clones)) {
        clones_expr <- clones[[clname]]
        if (is.null(clones_expr)) {
            data[[clname]] <- data$CTaa
        } else if (length(clones_expr) == 1 && grepl("(", clones_expr, fixed = TRUE) && grepl(")", clones_expr, fixed = TRUE)) {
            # if it's an expression, validate it can be parsed
            tryCatch(
                parse(text = clones_expr),
                error = function(e) {
                    stop(
                        "Invalid clone selector expression: '", clones_expr, "'.\n",
                        "See ?clone_selectors for available selector functions.\n",
                        "Parse error: ", e$message,
                        call. = FALSE
                    )
                }
            )
            if (clname == ".Clones") {
                # A single unnamed clone selector expression
                data <- data %>%
                    mutate(!!sym(clname) := !!parse_expr(clones_expr)) %>%
                    filter(!is.na(!!sym(clname)))

                if (nrow(data) == 0) {
                    stop("[ClonalStatPlot] No clones selected for selector `", clones_expr, "`.")
                }
            } else {
                tmp_data <- orig_data %>%
                    mutate(
                        .CloneGroups = !!parse_expr(clones_expr),
                        .CloneGroups = ifelse(is.na(!!sym(".CloneGroups")), NA_character_, clname)
                    ) %>%
                    filter(!is.na(!!sym(".CloneGroups")))

                if (nrow(tmp_data) == 0) {
                    stop("[ClonalStatPlot] No clones selected for clone group `", clname, "`.")
                } else {
                    cg_data <- rbind(cg_data, tmp_data)
                }
            }
        } else {
            if (clname == ".Clones") {
                data <- data %>%
                    mutate(!!sym(clname) := ifelse(!!sym("CTaa") %in% clones_expr, as.character(!!sym("CTaa")), NA_character_)) %>%
                    filter(!is.na(!!sym(clname)))

                if (nrow(data) == 0) {
                    stop("[ClonalStatPlot] No clones selected for selector '", paste(clones_expr[1:3], collapse = ", "), ", ...'.")
                }
            } else {
                tmp_data <- orig_data %>%
                    mutate(.CloneGroups = ifelse(!!sym("CTaa") %in% clones_expr, clname, NA_character_)) %>%
                    filter(!is.na(!!sym(".CloneGroups")))

                if (nrow(tmp_data) == 0) {
                    stop("[ClonalStatPlot] No clones selected for clone group `", clname, "`.")
                } else {
                    cg_data <- rbind(cg_data, tmp_data)
                }
            }
        }
    }
    if (!is.null(cg_data)) {
        data <- cg_data
        rm(cg_data)
    }
    x <- ifelse(by_clones, ".Clones", ".CloneGroups")
    ylab <- ifelse(values_by == "count", "Clone Size", ifelse(values_by == "fraction", "Relative Abundance", "Number of Clones"))

    data <- data %>%
        dplyr::group_by(!!!syms(unique(c(x, group_by, subgroup_by, facet_by, split_by)))) %>%
        summarise(count = sum(!!sym("count")), n = n(), fraction = sum(!!sym("fraction")), .groups = "drop") %>%
        arrange(!!sym(group_by), desc(!!sym("count")))

    if (by_clones && relabel) {
        clone_ids <- factor(data[[x]], levels = unique(data[[x]]))
        data[[x]] <- paste0("clone", as.numeric(clone_ids))
    }

    if (identical(plot_type, "chord")) {
        ChordPlot(data, from = x, to = group_by, y = values_by,
            facet_by = facet_by, split_by = split_by, ...)
    } else if (identical(plot_type, "col")) {
        if (!is.null(facet_by)) {
            stop("'facet_by' is not supported for 'col' plot. Please use 'split_by' instead.")
        }
        args <- rlang::dots_list(...)
        args$data <- data
        args$x <- x
        args$y <- values_by
        args$split_by <- split_by
        args$facet_by <- group_by
        args$facet_scales <- args$facet_scale %||% "free_x"
        args$fill_by <- args$fill_by %||% FALSE
        args$xlab <- xlab %||% ifelse(by_clones, "Clones", "Clone Groups")
        args$ylab <- ylab
        args$x_text_angle <- args$x_text_angle %||% 90
        args$legend.position <- args$legend.position %||% "none"
        args$aspect.ratio <- args$aspect.ratio %||% 0.6
        args$theme_args <- args$theme_args %||% list()
        args$theme_args$panel.grid.major.x <- args$theme_args$panel.grid.major.x %||% ggplot2::element_blank()
        p <- do_call(BarPlot, args)
        attr(p, "height") <- 0.6 * attr(p, "width")
        p
    } else if (identical(plot_type, "bar")) {
        BarPlot(data, x = group_by, y = values_by, group_by = x, group_name = xlab,
            facet_by = facet_by, split_by = split_by, xlab = xlab, ylab = ylab, ...)
    } else if (identical(plot_type, "box")) {
        BoxPlot(data, x = group_by, y = values_by, group_by = subgroup_by,
            facet_by = facet_by, split_by = split_by, xlab = xlab, ylab = ylab, ...)
    } else if (identical(plot_type, "violin")) {
        ViolinPlot(data, x = group_by, y = values_by, group_by = subgroup_by,
            facet_by = facet_by, split_by = split_by, xlab = xlab, ylab = ylab, ...)
    } else if (identical(plot_type, "heatmap")) {
        if (!is.null(facet_by)) {
            stop("'facet_by' is not supported for 'heatmap' plot. Please use 'split_by' instead.")
        }
        xlab <- xlab %||% ifelse(by_clones, "Clones", "Clone Groups")
        data[[xlab]] <- data[[x]]
        data[[x]] <- NULL
        args <- rlang::dots_list(...)
        args$data <- data
        args$in_form <- "long"
        args$values_by <- values_by
        args$rows_by <- group_by
        args$columns_by <- xlab
        args$name <- ylab
        args$split_by <- split_by
        args$show_row_names <- args$show_row_names %||% TRUE
        args$show_column_names <- args$show_column_names %||% TRUE

        do_call(Heatmap, args)
    } else if (identical(plot_type, "pies")) {
        if (!is.null(facet_by)) {
            stop("'facet_by' is not supported for 'pies' plot. Please use 'split_by' instead.")
        }
        if (is.null(subgroup_by)) {
            stop("'subgroup_by' is required for 'pies' plot. Please provide it.")
        }
        xlab <- xlab %||% ifelse(by_clones, "Clones", "Clone Groups")
        data[[xlab]] <- data[[x]]
        data[[x]] <- NULL
        args <- rlang::dots_list(...)
        args$data <- data
        args$name <- ylab
        args$rows_by <- group_by
        args$columns_by <- xlab
        args$values_by <- values_by
        args$in_form = "long"
        args$split_by <- split_by
        args$cell_type <- "pie"
        args$pie_group_by <- subgroup_by
        args$pie_values <- "sum"
        args$pie_size <- args$pie_size %||% sqrt
        args$pie_size_name <- args$pie_size_name %||% "Size"
        args$show_row_names <- args$show_row_names %||% TRUE
        args$show_column_names <- args$show_column_names %||% TRUE
        args$add_reticle <- args$add_reticle %||% TRUE

        do_call(Heatmap, args)
    } else if (identical(plot_type, "sankey")) {
        SankeyPlot(data, x = group_by, y = values_by, links_name = xlab, in_form = "long",
            alluvium = x, links_fill_by = x, flow = TRUE, xlab = xlab, ylab = ylab,
            facet_by = facet_by, split_by = split_by, ...)
    } else {
        TrendPlot(data, x = group_by, y = values_by, group_by = x,
            group_name = ifelse(by_clones, "Clones", "Clone Groups"),
            facet_by = facet_by, split_by = split_by, xlab = xlab, ylab = ylab, ...)
    }
}
