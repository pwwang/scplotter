
#' Visualize clone abundance, frequency, and dynamics across groups
#'
#' @description
#' ClonalStatPlot provides a unified interface for visualizing the abundance, frequency,
#' and dynamics of T cell and B cell clones across experimental groups. It is the most
#' versatile clone visualization function in scplotter, offering multiple plot types
#' for different analytical purposes.
#'
#' The function operates on the output of \code{\link[scRepertoire:combineTCR]{scRepertoire::combineTCR()}},
#' \code{\link[scRepertoire:combineBCR]{scRepertoire::combineBCR()}}, or
#' \code{\link[scRepertoire:combineExpression]{scRepertoire::combineExpression()}}. Clones are
#' identified by their CDR3 amino acid sequence, nucleotide sequence, V(D)J gene usage,
#' or a combination thereof (via \code{clone_call}). The function then computes clone-level
#' statistics (size, fraction, or count of clones) within each group and renders them
#' using one of ten supported plot types.
#'
#' A defining feature of ClonalStatPlot is its flexible clone selection system. Clones
#' can be specified directly by their IDs, or selected programmatically using expression
#' selectors such as \code{top()}, \code{sel()}, \code{shared()}, \code{uniq()}, and
#' comparison operators (\code{gt()}, \code{lt()}, \code{eq()}, etc.). These selectors
#' evaluate within the context of each faceting/splitting group, enabling per-group
#' selection of the most expanded clones, clones shared between conditions, or clones
#' meeting custom abundance thresholds. See the \strong{Clone selection} section below
#' and \code{\link{CloneSelectors}} for full details.
#'
#' Clones can also be aggregated into named groups (by passing a named list to
#' \code{clones}), where each group is defined by its own selection expression. In
#' this mode, the visualization unit becomes the clone group rather than individual
#' clones, enabling comparisons such as "hyper-expanded clones in condition A" vs.
#' "hyper-expanded clones in condition B."
#'
#' @section Clone selection:
#'
#' The \code{clones} argument accepts three forms:
#' \describe{
#'   \item{Character vector of clone IDs}{Directly specifies which clones to track.
#'     Clone IDs are matched against the column identified by \code{clone_call}
#'     (e.g., CDR3 amino acid sequences when \code{clone_call = "aa"}).}
#'   \item{Selection expression (single string with parentheses)}{A string containing
#'     a clone selector function call, e.g. \code{"top(10)"}, \code{"shared(P17B, P17L, group_by = 'Sample')"},
#'     or \code{"sel(P17L > 10 & P17B > 0, group_by = 'Sample')"}. The expression is
#'     parsed and evaluated within the data context. Available selectors include:
#'     \itemize{
#'       \item \code{top(n, ...)} — select the \code{n} largest clones (by total count)
#'       \item \code{sel(expr, ...)} — select clones matching a logical expression
#'       \item \code{shared(g1, g2, ...)} — select clones present in all specified groups
#'       \item \code{uniq(g1, g2, ...)} — select clones unique to group 1
#'       \item \code{gt(g1, g2)}, \code{lt(g1, g2)}, \code{eq(g1, g2)}, etc. — comparison-based selection
#'     }
#'     All selectors accept \code{group_by}, \code{top}, \code{order}, \code{within},
#'     and \code{output_within} arguments. See \code{\link{CloneSelectors}} for
#'     complete documentation.
#'   }
#'   \item{Named list of expressions}{Defines clone groups. Each element is a selection
#'     expression (as above), and the element name becomes the group label. For example:
#'     \code{list(ExpandedInA = "sel(A > 20, group_by = 'Sample')", ExpandedInB = "sel(B > 20, group_by = 'Sample')")}.
#'     In this mode, the visualization aggregates clones within each group rather than
#'     showing individual clones.}
#' }
#'
#' By default, clone selection operates within each faceting/splitting group (i.e.,
#' \code{top(3)} selects the top 3 clones per facet). Pass \code{group_by} explicitly
#' within the selector expression to change this behavior.
#'
#' @section Plot types:
#'
#' ClonalStatPlot supports ten plot types, each suited to different analytical questions:
#'
#' \describe{
#'   \item{\code{"bar"} (default)}{Stacked or grouped bar plot showing the total
#'     abundance (size or fraction) of each selected clone across groups. Best for
#'     comparing the composition of the top clones between conditions. Requires at
#'     least 1 group.}
#'   \item{\code{"box"}}{Box plot showing the distribution of individual clone sizes
#'     within each group. Useful for assessing whether clone size distributions differ
#'     between conditions. Optionally colored by \code{subgroup_by}.}
#'   \item{\code{"violin"}}{Violin plot alternative to box plot, showing the full
#'     density distribution of clone sizes. Supports \code{subgroup_by} for split
#'     violins.}
#'   \item{\code{"heatmap"}}{Heatmap where rows are clones (or clone groups) and
#'     columns are groups from \code{group_by}. Cell color encodes clone abundance.
#'     When \code{subgroup_by} is provided, rows are split by group and colored by
#'     subgroup. \code{facet_by} is not supported; use \code{split_by} instead.}
#'   \item{\code{"pies"}}{Heatmap variant where each cell contains a pie chart showing
#'     the composition of the clone (or clone group) across \code{subgroup_by} levels.
#'     The pie size reflects total abundance. \code{subgroup_by} is required.
#'     \code{within_subgroup} defaults to \code{FALSE} for this plot type.}
#'   \item{\code{"chord"} / \code{"circos"}}{Chord diagram showing clone flow between
#'     exactly 2 groups. Clones are represented as arcs, with ribbons indicating shared
#'     clones. For more than 2 groups, use \code{"sankey"} instead.}
#'   \item{\code{"sankey"} / \code{"alluvial"}}{Sankey (alluvial) diagram showing
#'     clone dynamics across groups. Flows are colored by clone groups (when using
#'     clone groups) or by individual clones. Best for tracking clone expansion,
#'     contraction, or sharing across multiple time points or conditions.}
#'   \item{\code{"trend"}}{Line plot showing the abundance trajectory of each clone
#'     (or clone group) across groups. Lines are colored by clone identity. Best for
#'     longitudinal data or dose-response experiments where group order is meaningful.}
#'   \item{\code{"col"}}{Column plot where each clone gets its own column, faceted
#'     by \code{group_by}. Unlike \code{"bar"}, this places clones on the x-axis.
#'     \code{facet_by} is not supported; use \code{split_by} instead. Clones are
#'     auto-relabeled by default.}
#' }
#'
#' @section Value types:
#'
#' The \code{values_by} parameter controls what is plotted on the y-axis:
#' \describe{
#'   \item{\code{"count"}}{The sum of cell counts for each clone within the group
#'     (i.e., clone size). This is the default.}
#'   \item{\code{"fraction"}}{The fraction of cells belonging to each clone, calculated
#'     as the clone's cell count divided by the total cells in the group. Suitable
#'     when group sizes differ and proportions are more meaningful than absolute counts.}
#'   \item{\code{"n"}}{The number of distinct clones (not cells) meeting the selection
#'     criteria. Shorthand for \code{"count"} and produces the same result.}
#' }
#'
#' @param data The product of \code{\link[scRepertoire:combineTCR]{scRepertoire::combineTCR()}},
#'   \code{\link[scRepertoire:combineBCR]{scRepertoire::combineBCR()}}, or
#'   \code{\link[scRepertoire:combineExpression]{scRepertoire::combineExpression()}}.
#'   A list of data frames where each element represents a sample, with columns for
#'   clone identifiers (CTaa, CTnt, CTgene, etc.) and cell-level metadata.
#' @param clones Which clones to track and visualize. Default: \code{"top(10)"}.
#'   Accepts three forms: (1) a character vector of clone IDs, (2) a single selection
#'   expression string (e.g. \code{"top(10)"}, \code{"sel(P17B > 5, group_by = 'Sample')"}),
#'   or (3) a named list of selection expressions to define clone groups (e.g.
#'   \code{list(Expanded = "sel(A > 20, group_by = 'Sample')")}). See the
#'   \strong{Clone selection} section for details. When a single unnamed expression
#'   is used, individual clones are visualized. When a named list is used, clone
#'   groups become the visualization unit.
#' @param clone_call How to identify a clone. One of \code{"gene"} (VDJC gene segment),
#'   \code{"nt"} (CDR3 nucleotide sequence), \code{"aa"} (CDR3 amino acid sequence,
#'   default), \code{"strict"} (VDJC gene + CDR3 nucleotide), or a custom column name
#'   present in the data.
#' @param chain Which TCR/BCR chain(s) to include. One of \code{"both"} (default,
#'   both chains combined), \code{"TRA"}, \code{"TRB"}, \code{"TRG"}, \code{"TRD"},
#'   \code{"IGH"}, \code{"IGL"}, or \code{"IGK"}. When \code{"both"}, dual-chain
#'   data (e.g., TRA and TRB) is combined.
#' @param values_by The metric to plot on the y-axis. One of \code{"count"} (default,
#'   number of cells per clone), \code{"fraction"} (proportion of cells per clone
#'   within each group), or \code{"n"} (equivalent to \code{"count"}). See the
#'   \strong{Value types} section.
#' @param relabel Whether to relabel clone IDs as "clone1", "clone2", etc., ordered
#'   by descending clone size. Default: \code{TRUE} for \code{"col"}, \code{"chord"},
#'   and \code{"circos"} plot types; \code{FALSE} otherwise. Useful when clone IDs
#'   are long CDR3 sequences. Only applies when visualizing individual clones
#'   (not clone groups).
#' @param plot_type The type of plot to generate. One of \code{"bar"} (default),
#'   \code{"box"}, \code{"violin"}, \code{"heatmap"}, \code{"pies"}, \code{"chord"}
#'   (or \code{"circos"}), \code{"sankey"} (or \code{"alluvial"}), \code{"trend"},
#'   or \code{"col"}. See the \strong{Plot types} section for guidance.
#' @param group_by The column name in the metadata to use for grouping cells
#'   (x-axis categories). Default: \code{"Sample"}. Only a single \code{group_by}
#'   column is supported.
#' @param groups The specific groups (levels of \code{group_by}) to include.
#'   Default: \code{NULL} (all groups included). If a named vector, names are used
#'   as display labels (e.g. \code{c(B = "P17B", L = "P17L")} renames "P17B" to "B").
#'   For \code{"chord"}/\code{"circos"}, exactly 2 groups are required. For
#'   \code{"box"}, \code{"violin"}, \code{"heatmap"}, \code{"pies"}, \code{"sankey"},
#'   and \code{"trend"}, at least 2 groups are required.
#' @param subgroup_by The column name in the metadata for subgrouping. Interpretation
#'   varies by plot type: for \code{"box"}/\code{"violin"}, it controls fill grouping;
#'   for \code{"heatmap"} with \code{"pies"}, it defines the pie chart composition;
#'   for \code{"heatmap"} without \code{"pies"}, it colors row labels. Not supported
#'   for \code{"bar"}, \code{"trend"}, or \code{"col"}. Default: \code{NULL}.
#' @param subgroups The specific subgroups (levels of \code{subgroup_by}) to include.
#'   Default: \code{NULL} (all subgroups included). If a vector, the same subgroups
#'   are applied to all groups. If a named list, different subgroups can be specified
#'   per group (names match \code{group_by} levels).
#' @param within_subgroup Whether clone selection (\code{clones}) should be performed
#'   within each subgroup separately. Default: \code{TRUE} for most plot types,
#'   \code{FALSE} for \code{"pies"}. When \code{TRUE}, clone selectors like
#'   \code{top(10)} select the top 10 clones within each subgroup rather than
#'   across all subgroups combined.
#' @param order A list specifying the order of levels for \code{group_by}. Default:
#'   \code{NULL} (uses the order present in the data). Lower priority than \code{groups}.
#' @param facet_by A column name to facet the plot into separate panels. Default:
#'   \code{NULL}. Not supported for \code{"col"}, \code{"heatmap"}, or \code{"pies"}
#'   plot types (use \code{split_by} instead).
#' @param split_by A column name to split the plot into separate subplots (via
#'   \pkg{patchwork}). Default: \code{NULL}. Unlike \code{facet_by}, splitting
#'   creates independent plots that can have different scales.
#' @param y The y-axis variable. Default: \code{NULL} (auto-determined from
#'   \code{values_by}). For \code{"bar"} plots, can be \code{"TotalSize"}
#'   (total cells in selected clones) or \code{"Count"} (number of selected clones).
#' @param xlab Custom x-axis label. Default: \code{NULL} (auto-generated).
#' @param ylab Custom y-axis label. Default: \code{NULL} (auto-generated based on
#'   \code{values_by}: "Clone Size", "Relative Abundance", or "Number of Clones").
#' @param ... Additional arguments passed to the underlying plot function from
#'   \pkg{plotthis}. For example:
#'   \itemize{
#'     \item For \code{"bar"}: see \code{\link[plotthis:BarPlot]{plotthis::BarPlot()}} (e.g. \code{position}, \code{palette})
#'     \item For \code{"box"}: see \code{\link[plotthis:BoxPlot]{plotthis::BoxPlot()}} (e.g. \code{add_box}, \code{comparison})
#'     \item For \code{"violin"}: see \code{\link[plotthis:ViolinPlot]{plotthis::ViolinPlot()}} (e.g. \code{add_box}, \code{comparison})
#'     \item For \code{"heatmap"} and \code{"pies"}: see \code{\link[plotthis:Heatmap]{plotthis::Heatmap()}} (e.g. \code{palette}, \code{show_row_names})
#'     \item For \code{"sankey"}: see \code{\link[plotthis:SankeyPlot]{plotthis::SankeyPlot()}} (e.g. \code{flow}, \code{node_palette})
#'     \item For \code{"trend"}: see \code{\link[plotthis:TrendPlot]{plotthis::TrendPlot()}} (e.g. \code{line_type}, \code{palette})
#'     \item For \code{"chord"}: see \code{\link[plotthis:ChordPlot]{plotthis::ChordPlot()}}
#'     \item For \code{"col"}: see \code{\link[plotthis:BarPlot]{plotthis::BarPlot()}} (used internally with faceting)
#'   }
#'   Additional arguments for \code{"col"} plot include \code{fill_by}, \code{fill_name},
#'   \code{facet_scale}, \code{facet_ncol}, \code{x_text_angle}, \code{aspect.ratio},
#'   \code{legend.position}, and \code{theme_args}.
#' @return A \code{ggplot} object (or a \code{patchwork} object if \code{split_by}
#'   is used) invisibly.
#' @importFrom rlang parse_expr as_label enexpr
#' @importFrom dplyr %>% summarise arrange desc mutate ungroup across slice_head n filter
#' @importFrom tidyr pivot_wider
#' @importFrom plotthis SankeyPlot TrendPlot
#' @note
#' \itemize{
#'   \item ClonalStatPlot requires at least 2 groups for \code{"box"}, \code{"violin"},
#'     \code{"heatmap"}, \code{"pies"}, \code{"sankey"}, and \code{"trend"} plot types.
#'     Only \code{"bar"} and \code{"col"} work with a single group.
#'   \item \code{"chord"}/\code{"circos"} plots are limited to exactly 2 groups.
#'     For more groups, use \code{"sankey"} instead.
#'   \item \code{facet_by} is not supported for \code{"col"}, \code{"heatmap"}, and
#'     \code{"pies"} plot types because these plots use internal faceting. Use
#'     \code{split_by} as an alternative for creating separate subplots.
#'   \item \code{subgroup_by} is not supported for \code{"bar"}, \code{"trend"}, and
#'     \code{"col"} plot types.
#'   \item When using clone groups (a named list for \code{clones}), the
#'     \code{relabel} argument has no effect since group names are used directly.
#'   \item Clone selection expressions are evaluated after the data is filtered to
#'     the specified \code{groups}. If you reference group names in your expression
#'     (e.g., \code{"sel(P17B > 10)"}), ensure those groups are included in
#'     \code{groups} if they differ from the display groups.
#'   \item For \code{"pies"} plots, \code{within_subgroup} defaults to \code{FALSE},
#'     meaning clone selection occurs across all subgroups combined. Set to
#'     \code{TRUE} to select clones within each subgroup independently.
#' }
#' @seealso
#' \itemize{
#'   \item \code{\link{CloneSelectors}} for the full clone selection expression system
#'   \item \code{\link{ClonalCompositionPlot}} for visualizing clonal space composition (homeostasis)
#'   \item \code{\link{ClonalDiversityPlot}} for clonal diversity metrics
#'   \item \code{\link{ClonalGeneUsagePlot}} for V(D)J gene segment usage
#'   \item \code{\link{ClonalPositionalPlot}} for CDR3 positional analysis
#'   \item \code{\link{ClonalKmerPlot}} for CDR3 k-mer motif analysis
#' }
#' @export
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
#' # add a fake variable (e.g. cell type from scRNA-seq)
#' data <- lapply(data, function(x) {
#'     x$CellType <- factor(
#'         sample(c("CD4", "CD8", "B", "NK"), nrow(x), replace = TRUE),
#'        levels = c("CD8", "CD4", "B", "NK")
#'     )
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
#' # using pies with subgroups for groups of clones
#' ClonalStatPlot(data, group_by = "Sample", groups = c("P17B", "P17L"),
#'     clones = list(
#'         ExpandedClonesInP17L = "sel(P17L > 20, group_by = 'Sample')",
#'         ExpandedClonesInP17B = "sel(P17B > 20, group_by = 'Sample')"
#'     ), subgroup_by = "CellType", pie_size = sqrt,
#'     plot_type = "pies", show_row_names = TRUE, show_column_names = TRUE,
#'     title = "Clones larger than 20 in P17L and P17B (pies with subgroups by CellType)")
#' # using heatmap with subgroups for groups of clones
#' ClonalStatPlot(data, group_by = "Sample", groups = c("P17B", "P17L"),
#'    clones = list(
#'        ExpandedClonesInP17L = "sel(P17L > 20, group_by = 'Sample')",
#'        ExpandedClonesInP17B = "sel(P17B > 20, group_by = 'Sample')"
#'    ), subgroup_by = "CellType", pie_size = sqrt, within_subgroup = FALSE,
#'    plot_type = "heatmap", show_row_names = TRUE, show_column_names = TRUE,
#'    title = "Clones larger than 20 in P17L and P17B (pies with subgroups by CellType)")
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
#' # Rename groups
#' ClonalStatPlot(data, plot_type = "col", groups = c(P17B = "B", P17L = "L"),
#'     facet_ncol = 1, legend.position = "right",
#'     relabel = TRUE, fill_by = ".Clones", fill_name = "Clones")
#' }
ClonalStatPlot <- function(
    data, clones = "top(10)", clone_call = "aa", chain = "both", values_by = c("count", "fraction", "n"),
    plot_type = c("bar", "box", "violin", "heatmap", "pies", "circos", "chord", "sankey", "alluvial", "trend", "col"),
    group_by = "Sample", groups = NULL, subgroup_by = NULL, subgroups = NULL, order = NULL,
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
    data <- clonal_size_data(data, clone_call, chain, all_groupings, order)
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
        groups <- unlist(groups)
        if (is.null(names(groups))) {
            names(groups) <- groups
        }
        nonexist_groups <- setdiff(names(groups), levels(data[[group_by]]))
        if (length(nonexist_groups) > 0) {
            stop(paste("The following groups do not exist in the data:", paste(nonexist_groups, collapse = ", ")))
        }
        data <- data %>% filter(!!sym(group_by) %in% names(groups))
        # Reverse the names and values for fct_recode to rename the groups
        groups <- stats::setNames(names(groups), groups)
        data[[group_by]] <- forcats::fct_recode(data[[group_by]], !!!groups)
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
        data[[subgroup_by]] <- factor(data[[subgroup_by]], levels = subgroups)
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
                        "See ?CloneSelectors for available selector functions.\n",
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
        args$columns_by <- xlab
        args$name <- ylab
        args$split_by <- split_by
        args$show_row_names <- args$show_row_names %||% TRUE
        args$show_column_names <- args$show_column_names %||% TRUE
        if (!is.null(subgroup_by)) {
            args$rows_split_by <- group_by
            args$rows_by <- subgroup_by
        } else {
            args$rows_by <- group_by
        }

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
