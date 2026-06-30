#' Visualize TCR/BCR gene segment usage
#'
#' Adaptive immune receptors (TCRs and BCRs) are assembled through V(D)J recombination,
#' where variable (V), diversity (D), and joining (J) gene segments are randomly selected
#' and rearranged. The frequency with which different gene segments are used — termed
#' *gene usage* — provides insight into immune repertoire composition, T/B cell
#' development, and antigen-driven selection. Skewed gene usage can indicate clonal
#' expansion, immune aging, or disease-associated repertoire bias.
#'
#' `ClonalGeneUsagePlot` visualizes the usage frequency of TCR or BCR gene segments
#' across samples and conditions. It supports two modes of analysis:
#'
#' * **Single gene analysis** (`genes = "TRBV"`): The second axis is the `group_by`
#'   variable (typically Sample), producing a gene-by-sample matrix suitable for
#'   comparing usage across conditions.
#' * **Dual gene analysis** (`genes = c("TRBV", "TRBJ")`): Both axes are gene segments,
#'   and the `group_by` variable is used for faceting. This reveals preferential
#'   V-J pairings, which can reflect structural constraints in the receptor or
#'   antigen-driven selection for specific chain combinations.
#'
#' Gene usage can be displayed as raw counts or as proportions scaled within each
#' group (when `scale = TRUE`), the latter being more appropriate for comparing
#' samples of different sequencing depths.
#'
#' The function uses \code{\link[scRepertoire:vizGenes]{scRepertoire::vizGenes()}}
#' internally to compute gene usage frequencies, then passes the result to the
#' appropriate \pkg{plotthis} visualization function.
#'
#' @section Single vs. dual gene analysis:
#' The `genes` parameter determines the dimensionality of the analysis:
#'
#' * **Single gene prefix** (e.g., `"TRBV"`, `"TRBJ"`, `"IGHV"`): The x-axis shows
#'   gene segments (e.g., TRBV5-1, TRBV7-2), and the y-axis (or columns in a heatmap)
#'   shows the `group_by` variable, typically Sample. This is the standard mode for
#'   comparing gene usage across samples or conditions. Use a bar plot to quickly
#'   identify over- or under-represented genes, or a heatmap for a compact overview
#'   of many genes across many samples.
#'
#' * **Two gene prefixes** (e.g., `c("TRBV", "TRBJ")`): The x-axis shows the first
#'   gene set and the y-axis (or columns) shows the second gene set. Each cell
#'   represents a specific V-J pair, and the `group_by` variable is used for
#'   faceting. This mode is essential for studying chain pairing — certain V segments
#'   preferentially pair with specific J segments, and deviations from expected
#'   pairing frequencies can indicate structural constraints or disease-associated
#'   repertoire features. The sankey/alluvial plot type is particularly effective
#'   for dual gene analysis, showing the flow from one gene set to the other.
#'
#' @section Plot types:
#' Five visualization types are available, each suited to different analytical goals:
#'
#' * **`"bar"`** (default): Horizontal bar chart showing gene usage per sample/group.
#'   Best for comparing a moderate number of genes across a few samples. Gene labels
#'   are rotated 90 degrees for readability. The `aspect.ratio` parameter controls
#'   bar height (defaults to `2 / top` to automatically scale with the number of genes).
#'
#' * **`"heatmap"`**: Matrix heatmap with genes on rows and samples/conditions on
#'   columns. Ideal for surveying many genes across many samples simultaneously.
#'   When a single gene prefix is used, a "Total Usage" row annotation is
#'   automatically added showing the aggregate usage per gene as a line plot.
#'   Supports custom row annotations, annotation types, and aggregation functions.
#'
#' * **`"circos"` / `"chord"`**: Chord diagram showing flow from gene segments to
#'   samples (single gene) or from one gene set to another (dual gene). Chord
#'   diagrams excel at revealing broad patterns of connectivity. In dual gene mode,
#'   the plot is automatically split by the `group_by` variable.
#'
#' * **`"alluvial"` / `"sankey"`**: Sankey diagram showing the flow from gene
#'   segments to samples (single gene) or from one gene set to another (dual gene).
#'   Unlike chord diagrams, sankey diagrams preserve the ordering of categories
#'   and are often easier to read when there are many connections. In dual gene
#'   mode, the plot is faceted by `group_by`. Note that `"alluvial"` is
#'   automatically mapped to `"sankey"`.
#'
#' @section Gene segment naming:
#' Gene segment prefixes follow standard IMGT nomenclature:
#'
#' * **Human TCR**: `TRBV`, `TRBD`, `TRBJ` (beta chain); `TRAV`, `TRAJ` (alpha chain);
#'   `TRGV`, `TRGJ` (gamma chain); `TRDV`, `TRDJ` (delta chain)
#' * **Human BCR**: `IGHV`, `IGHD`, `IGHJ` (heavy chain); `IGKV`, `IGKJ` (kappa chain);
#'   `IGLV`, `IGLJ` (lambda chain)
#' * **Mouse**: Prefixes are similar but use `Trbv`, `Ighv`, etc. (sentence case)
#'
#' The gene prefix is used to identify matching columns in the data. Only genes
#' matching the prefix will be included in the analysis. The `top` parameter
#' selects the most frequently used genes/genepairs for display.
#'
#' @param data The product of \code{\link[scRepertoire:combineTCR]{scRepertoire::combineTCR()}},
#'  \code{\link[scRepertoire:combineBCR]{scRepertoire::combineBCR()}}, or
#'  \code{\link[scRepertoire:combineExpression]{scRepertoire::combineExpression()}}.
#'  Must contain columns with gene segment information (e.g., TRBV_, TRBJ_ for TCR beta chain).
#' @param genes A character vector of gene segment prefixes to analyze. Default is `"TRBV"`.
#'  If a single prefix is provided (e.g., `"TRBV"`), the second dimension of the plot
#'  will be the `group_by` variable. If two prefixes are provided
#'  (e.g., `c("TRBV", "TRBJ")`), both axes represent gene segments, revealing
#'  V-J pairing frequencies, and `group_by` is used for faceting. Must be of
#'  length 1 or 2.
#' @param scale Logical; whether to normalize gene usage as proportions within each
#'  group (`TRUE`, default), or display raw counts (`FALSE`). Setting `scale = TRUE`
#'  is recommended when comparing samples of different sequencing depths, as it
#'  accounts for differences in total clone counts. Raw counts may be preferred
#'  when absolute clone numbers are biologically meaningful.
#' @param top Integer; the number of top genes (single gene mode) or top gene pairs
#'  (dual gene mode) to display. Genes are ranked by total usage across all groups.
#'  Default is `20`. Increase this value to show more genes, or decrease it to
#'  focus on the most dominant segments.
#' @param plot_type Character; the type of visualization. One of `"bar"` (default),
#'  `"heatmap"`, `"circos"`, `"chord"`, `"alluvial"`, or `"sankey"`. Note that
#'  `"alluvial"` and `"sankey"` are treated identically. See the **Plot types**
#'  section for guidance on selecting the appropriate visualization.
#' @param group_by Character vector; the column(s) in `data` to group by.
#'  Default is `"Sample"`. For single gene analysis, this becomes the second
#'  dimension of the plot. For dual gene analysis, it is used for faceting.
#'  Set to `NULL` to disable grouping.
#' @param order A named list specifying the order of factor levels for grouping
#'  variables. For example, `list(Sample = c("P17B", "P17L"))` or
#'  `list(Type = c("L", "B"))`. Default is `NULL`, which uses the order
#'  present in the data. Names in the list must match column names in `data`.
#' @param facet_by A character vector of column names to facet the plots by.
#'  Default is `NULL`. This parameter is typically set automatically by the
#'  function and should not be specified manually — providing a value will
#'  raise an error.
#' @param facet_ncol Integer; the number of columns in the facet grid when
#'  faceting is applied. Default is `1`. Only relevant for plot types that
#'  produce faceted output (e.g., sankey in dual gene mode).
#' @param split_by A character vector of column names to split the plots by.
#'  Default is `NULL`. When provided, separate plots are generated for each
#'  combination of values in the specified columns. For dual gene analysis,
#'  `split_by` is used in single gene mode; in dual gene mode, the plot is
#'  automatically split by `group_by` and specifying `split_by` will raise
#'  an error.
#' @param aspect.ratio Numeric; the aspect ratio (height / width) of bar plot
#'  panels. Default is `2 / top`, which automatically scales the ratio based
#'  on the number of genes displayed. Only applicable for `plot_type = "bar"`.
#' @param theme_args A named list of arguments passed to
#'  \code{\link[ggplot2:theme]{ggplot2::theme()}} for customizing the plot
#'  appearance. For bar plots, `panel.grid.major.y` defaults to
#'  `element_blank()` for a cleaner look.
#' @param ylab Character; the y-axis label. Default is `NULL`, which
#'  automatically uses `"Gene Usage Fraction"` when `scale = TRUE` or
#'  `"Gene Usage Count"` when `scale = FALSE`.
#' @param row_annotation A named list specifying row annotations for heatmap
#'  plots. Each element should be a column name in the data to use as
#'  annotation data. Default is `NULL`. When a single gene prefix is used,
#'  a `"Total Usage"` annotation is automatically added showing the aggregate
#'  usage per gene. Only applicable for `plot_type = "heatmap"`.
#' @param row_annotation_type A named list specifying the annotation type for
#'  each row annotation. For example, `list("Total Usage" = "lines")`. Default
#'  is an empty list. The `"Total Usage"` annotation defaults to `"lines"`.
#'  Only applicable for `plot_type = "heatmap"`.
#' @param row_annotation_side Character; the side of the heatmap where row
#'  annotations are placed. One of `"right"` (default) or `"left"`. Only
#'  applicable for `plot_type = "heatmap"`.
#' @param row_annotation_agg A named list of aggregation functions for row
#'  annotations when multiple values exist per row. For example,
#'  `list("Total Usage" = function(x) ifelse(length(x) > 1, x[1], 0))`.
#'  Default is an empty list. Only applicable for `plot_type = "heatmap"`.
#' @param show_row_names Logical; whether to display row names (gene segment
#'  names) in the heatmap. Default is `TRUE`. Set to `FALSE` to hide gene
#'  labels when there are too many to display legibly. Only applicable for
#'  `plot_type = "heatmap"`.
#' @param show_column_names Logical; whether to display column names (sample/
#'  group names or second gene set names) in the heatmap. Default is `TRUE`.
#'  Only applicable for `plot_type = "heatmap"`.
#' @param ... Additional arguments passed to the underlying \pkg{plotthis}
#'  visualization function, depending on `plot_type`:
#'  * For `"bar"`: \code{\link[plotthis:BarPlot]{plotthis::BarPlot()}}
#'  * For `"heatmap"`: \code{\link[plotthis:Heatmap]{plotthis::Heatmap()}}
#'  * For `"circos"` / `"chord"`: \code{\link[plotthis:ChordPlot]{plotthis::ChordPlot()}}
#'  * For `"alluvial"` / `"sankey"`: \code{\link[plotthis:SankeyPlot]{plotthis::SankeyPlot()}}
#'
#'  Common arguments include `title` (plot title), `legend.position`, and
#'  color palette parameters. See the respective \pkg{plotthis} documentation
#'  for available options.
#' @return A `ggplot` object (or a list of `ggplot` objects if `combine = FALSE`
#'  is passed via `...`)
#' @note
#' * Gene usage data is computed using
#'   \code{\link[scRepertoire:vizGenes]{scRepertoire::vizGenes()}}, which
#'   requires the input data to contain gene segment columns (e.g., `TRBV_gene`,
#'   `TRBJ_gene`). These columns are automatically created by
#'   \code{\link[scRepertoire:combineTCR]{scRepertoire::combineTCR()}} and
#'   \code{\link[scRepertoire:combineBCR]{scRepertoire::combineBCR()}}.
#' * For bar plots, `aspect.ratio` is automatically set to `2 / top` to
#'   accommodate varying numbers of genes. You can override this by passing
#'   a custom value.
#' * The `facet_by` parameter is set internally by the function and should not
#'   be provided by the user. Doing so will raise an error.
#' * Dual gene analysis with `split_by` in chord/circos mode is not supported —
#'   the plot is automatically split by `group_by` instead.
#' * When using `scale = TRUE`, proportions sum to 1 within each group for
#'   each gene set, making it possible to compare relative usage across
#'   groups of different sizes. However, this normalization can mask
#'   differences in absolute clone counts — use `scale = FALSE` when
#'   total clone numbers are of interest.
#' * Gene segment names follow IMGT nomenclature. Ensure your gene prefixes
#'   match the naming convention in your dataset (e.g., `"TRBV"` for human
#'   TRB V genes, `"Trbv"` for mouse).
#' @seealso
#' * \code{\link{ClonalCompositionPlot}} for analyzing clonal homeostasis and
#'   expansion/contraction composition
#' * \code{\link{ClonalDiversityPlot}} for analyzing clonal diversity metrics
#' * \code{\link{ClonalPositionalPlot}} for analyzing amino acid positional
#'   distributions within CDR3 sequences
#' * \code{\link{ClonalKmerPlot}} for analyzing k-mer motifs in CDR3 sequences
#' * \code{\link[scRepertoire:vizGenes]{scRepertoire::vizGenes()}} for the
#'   underlying gene usage computation
#' @importFrom rlang := syms
#' @importFrom dplyr mutate rename slice_max summarise pull filter all_of ungroup
#' @importFrom tidyr separate
#' @importFrom ggplot2 unit element_blank
#' @importFrom scRepertoire vizGenes
#' @importFrom plotthis BarPlot Heatmap ChordPlot SankeyPlot
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
#'
#' ClonalGeneUsagePlot(data)
#' ClonalGeneUsagePlot(data, genes = c("TRBJ", "TRBV"))
#' ClonalGeneUsagePlot(data, top = 40, plot_type = "heatmap")
#' ClonalGeneUsagePlot(data, genes = c("TRBV", "TRBJ"), plot_type = "heatmap")
#' ClonalGeneUsagePlot(data, genes = "TRBV", group_by = "Type", plot_type = "chord")
#' ClonalGeneUsagePlot(data, genes = c("TRBV", "TRBJ"), group_by = "Type", plot_type = "chord")
#' ClonalGeneUsagePlot(data, genes = c("TRBV", "TRBJ"), plot_type = "alluvial",
#'      facet_scales = "free_y")
#' ClonalGeneUsagePlot(data, genes = c("TRBV", "TRBJ"), plot_type = "alluvial",
#'      group_by = NULL)
#' }
ClonalGeneUsagePlot <- function(
    data, genes = "TRBV", scale = TRUE, top = 20, order = NULL,
    plot_type = c("bar", "heatmap", "circos", "chord", "alluvial", "sankey"),
    group_by = "Sample", facet_by = NULL, facet_ncol = 1, split_by = NULL,
    aspect.ratio = 2 / top, theme_args = list(), ylab = NULL,
    show_row_names = TRUE, show_column_names = TRUE,
    row_annotation = NULL, row_annotation_type = list(), row_annotation_side = "right",
    row_annotation_agg = list(), ...
) {
    plot_type <- match.arg(plot_type)
    if (plot_type == "alluvial") { plot_type <- "sankey" }
    if (!is.null(facet_by)) {
        stop("'facet_by' should not be specified in ClonalGeneUsagePlot.")
    }

    if (length(genes) > 1) {
        genes2 <- genes[2]
        genes <- genes[1]
    } else {
        genes2 <- NULL
    }

    all_groupings <- unique(c(group_by, split_by))
    grouping_levels <- get_clonal_grouping_levels(data, all_groupings, order)
    data <- merge_clonal_groupings(data, all_groupings)
    data <- vizGenes(data, x.axis = genes, y.axis = genes2, group.by = ".group", scale = scale, exportTable = TRUE)
    if (is.null(genes2)) {
        stopifnot("[ClonalGeneUsagePlot] 'genes' must be of length 2." = plot_type != "sankey")
        axis1 <- genes
        axis2 <- group_by

        data <- separate(data, "y.axis", into = all_groupings, sep = " // ") %>%
            rename(!!sym(axis1) := "x.axis")

        for (gl in names(grouping_levels)) {
            if (!is.null(data[[gl]])) {
                data[[gl]] <- factor(data[[gl]], levels = grouping_levels[[gl]])
            }
        }

        selected_genes <- data %>% dplyr::group_by(!!sym(genes)) %>%
            summarise(total = sum(!!sym("count")), .groups = "drop") %>%
            slice_max(order_by = !!sym("total"), n = top) %>%
            pull(axis1)%>%
            as.character()

        selected_genes_levels <- if (is.factor(data[[axis1]])) {
            levels(data[[axis1]])
        } else {
            unique(data[[axis1]])
        }
        selected_genes_levels <- selected_genes_levels[selected_genes_levels %in% selected_genes]

        data <- data %>% filter(!!sym(axis1) %in% selected_genes) %>%
            mutate(!!sym(axis1) := factor(!!sym(axis1), levels = selected_genes_levels))
    } else {
        axis1 <- genes
        axis2 <- genes2

        if (!is.null(all_groupings)) {
            data <- separate(data, "element.names", into = all_groupings, sep = " // ")
            for (gl in names(grouping_levels)) {
                if (!is.null(data[[gl]])) {
                    data[[gl]] <- factor(data[[gl]], levels = grouping_levels[[gl]])
                }
            }
        }
        data <- data %>%
            rename(!!sym(genes) := "x.axis", !!sym(genes2) := "y.axis") %>%
            unite("GenePairs", c(genes, genes2), sep = " // ", remove = FALSE)

        genepairs <- data %>% dplyr::group_by(!!!syms(c(genes, genes2))) %>%
            summarise(total = sum(!!sym("count")), .groups = "drop") %>%
            slice_max(order_by = !!sym("total"), n = top) %>%
            unite("GenePairs", c(genes, genes2), sep = " // ") %>%
            pull("GenePairs")

        data <- data %>% filter(!!sym("GenePairs") %in% genepairs)
        data$GenePairs <- NULL
    }

    if (plot_type == "bar") {
        # theme_args$panel.spacing <- theme_args$panel.spacing %||% unit(-0.1, "lines")
        theme_args$panel.grid.major.y <- theme_args$panel.grid.major.y %||% element_blank()
        ylab <- ylab %||% ifelse(scale, "Gene Usage Fraction", "Gene Usage Count")

        BarPlot(data, x = axis1, y = ifelse(scale, "proportion", "count"), facet_by = axis2, facet_ncol = facet_ncol,
            split_by = split_by, legend.position = "none", x_text_angle = 90, aspect.ratio = aspect.ratio,
            facet_args = list(strip.position = "right"), theme_args = theme_args, ylab = ylab, ...)
    } else if (plot_type == "heatmap") {
        row_annotation_type[["Total Usage"]] <- row_annotation_type[["Total Usage"]] %||% "lines"
        row_annotation_agg[["Total Usage"]] <- row_annotation_agg[["Total Usage"]] %||% function(x) ifelse(length(x) > 1, x[1], 0)
        if (is.null(genes2)) {
            data <- data %>% dplyr::group_by(!!!syms(c(axis1, split_by))) %>%
                mutate(.total = sum(!!sym(ifelse(scale, "proportion", "count")))) %>%
                ungroup()
            row_annotation <- row_annotation %||% list(`Total Usage` = ".total")
        }

        Heatmap(
            data,
            values_by = ifelse(scale, "proportion", "count"), values_fill = 0,
            rows_by = axis1, columns_by = axis2, split_by = split_by,
            rows_name = axis1, name = ifelse(scale, "Gene Usage Fraction", "Gene Usage Count"),
            row_annotation = row_annotation, row_annotation_type = row_annotation_type, row_annotation_agg = row_annotation_agg,
            row_annotation_side = row_annotation_side,
            show_row_names = show_row_names, show_column_names = show_column_names,
            ...)
    } else if (plot_type %in% c("circos", "chord")) {
        if (is.null(genes2)) {
            ChordPlot(data, from = axis1, to = axis2, y = ifelse(scale, "proportion", "count"),
                split_by = split_by, theme_args = theme_args, ...)
        } else {
            if (!is.null(split_by)) {
                stop("[ClonalGeneUsagePlot] 'split_by' should not be specified when 'genes' has length 2, since the plot will be split by 'group_by'.")
            }
            ChordPlot(data, from = axis1, to = axis2, y = ifelse(scale, "proportion", "count"),
                split_by = group_by, theme_args = theme_args, ...)
        }
    } else {  # alluvial / sankey
        SankeyPlot(data, x = c(axis1, axis2), y = ifelse(scale, "proportion", "count"),
            links_fill_by = axis1, facet_by = group_by,
            split_by = split_by, theme_args = theme_args, ...)
    }
}
