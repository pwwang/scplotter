#' Cell statistics plot
#'
#' @description
#' Visualizes cell-level statistics — counts, fractions, and composition —
#' across cell identities and metadata groupings. This is the primary function
#' for exploring the distribution of cell types, clusters, and categorical
#' metadata in single-cell transcriptomics datasets. It answers questions such
#' as: "What proportion of each cell type is in each condition?", "How do
#' cluster abundances change across samples?", and "What is the clonal
#' composition within each cell type?"
#'
#' `CellStatPlot` serves as a unified interface across 15+ visualization
#' types, all driven by a common data aggregation and fraction-calculation
#' pipeline. It supports four single-cell data containers:
#' \itemize{
#'   \item **Seurat objects** — Extracts `@meta.data`; uses `Idents()` as
#'         the default identity when `ident = NULL`.
#'   \item **Giotto objects** — Extracts cell metadata via
#'         `getCellMetadata()` using `spat_unit` and `feat_type`.
#'   \item **h5ad files** (.h5ad or opened `H5File`) — Reads from `obs`
#'         via `h5group_to_dataframe()`.
#'   \item **Data frames** — Internal method; all other methods ultimately
#'         delegate here after metadata extraction.
#' }
#'
#' @section Plot types:
#' The `plot_type` parameter selects the visualization. Plot types fall into
#' several conceptual categories:
#'
#' **Composition within groups** (requires `group_by`):
#' \itemize{
#'   \item `"bar"` — Grouped or stacked bar chart. The default and most
#'         common choice. Shows count or fraction per ident, colored by
#'         group.
#'   \item `"trend"` — Line chart connecting points across idents. Best
#'         for showing trends across ordered categories (e.g., pseudotime
#'         bins, dose levels).
#'   \item `"area"` — Stacked area chart. Similar to trend but emphasizes
#'         the cumulative composition.
#'   \item `"ring"` / `"donut"` — Ring (donut) chart. A radial alternative
#'         to stacked bars. `frac` is forced to `"group"`.
#'   \item `"radar"` — Radar chart showing each ident as an axis. Compact
#'         for comparing group profiles across multiple idents.
#'   \item `"spider"` — Spider chart with filled polygons. Similar to radar
#'         but emphasizes the area covered by each group.
#' }
#'
#' **Flow between categories** (requires `group_by`):
#' \itemize{
#'   \item `"circos"` — Circos plot showing directed edges from one
#'         category to another (e.g., cluster to condition).
#'   \item `"sankey"` / `"alluvial"` — Sankey (alluvial) diagram for
#'         multi-step categorical flows. Supports multiple `group_by`
#'         columns to show cascading relationships.
#' }
#'
#' **Single-group composition** (no `group_by`):
#' \itemize{
#'   \item `"pie"` — Pie chart of cell counts per ident. Use for a quick
#'         overview of cluster proportions.
#' }
#'
#' **Matrix views** (requires `group_by`):
#' \itemize{
#'   \item `"heatmap"` — Heatmap of counts or fractions with idents as
#'         rows and `group_by` categories as columns. Supports row
#'         splitting via `rows_split_by` and column splitting via
#'         `columns_split_by`.
#'   \item `"pies"` — Heatmap where each cell is a pie chart showing
#'         sub-composition (e.g., clone distribution within each
#'         cluster-by-condition intersection). Requires `rows_by`.
#' }
#'
#' **Distribution within idents** (requires `group_by`):
#' \itemize{
#'   \item `"violin"` — Violin plot showing the distribution of values per
#'         ident, grouped or split by metadata.
#'   \item `"box"` — Box plot with the same structure as violin.
#' }
#'
#' @section Fraction calculation:
#' The `frac` parameter controls how cell counts are normalized:
#' \itemize{
#'   \item `"none"` — Raw cell counts (default). The y-axis shows the
#'         absolute number of cells.
#'   \item `"group"` — Fraction within each group. The total across all
#'         idents within each group sums to 1. Answers: "Of the cells in
#'         condition X, what fraction are each cell type?"
#'   \item `"ident"` — Fraction within each ident. The total across all
#'         groups within each ident sums to 1. Answers: "Of the Beta cells,
#'         what fraction are in each condition?" Requires `group_by`.
#'   \item `"cluster"` — Alias for `"ident"`.
#'   \item `"all"` — Fraction of all cells in the data (or split/facet
#'         subset). Answers: "What fraction of all cells are Beta cells in
#'         condition X?"
#' }
#' Fractions are calculated independently within each `split_by` and/or
#' `facet_by` subset.
#'
#' @section Custom aggregation:
#' By default, `CellStatPlot` counts the number of cells in each
#' ident-by-group intersection (`agg = "n()"`). The `agg` parameter accepts
#' any expression that can be passed to `dplyr::summarise()`, enabling
#' custom metrics. Examples:
#' \itemize{
#'   \item `"sum(hasTCR) / n()"` — Fraction of cells with a TCR in each
#'         group
#'   \item `"mean(expression_score)"` — Mean of a numeric metadata column
#'   \item `"sum(clone_size > 1)"` — Number of expanded clones
#' }
#' Note: `agg` is ignored for `"circos"` and `"pies"` plot types, which
#' always count cells.
#'
#' @section The ident/group_by duality:
#' `CellStatPlot` is built around two categorical axes:
#' \itemize{
#'   \item **`ident`** — The primary cell identity (clusters, cell types).
#'         This typically forms the x-axis, pie slices, or heatmap rows.
#'   \item **`group_by`** — A secondary grouping (conditions, samples,
#'         time points). This typically forms the fill colors, bar stacks,
#'         or heatmap columns.
#' }
#' The `swap` parameter exchanges these two roles: when `swap = TRUE`, the
#' `group_by` variable is used as the x-axis and `ident` provides the fill
#' colors. The exact behavior varies by plot type — see the `swap`
#' parameter documentation for details.
#'
#' @param object A Seurat object, a Giotto object, a path to an `.h5ad`
#'   file, an opened `H5File` from the \pkg{hdf5r} package, or a data
#'   frame (internal method) containing cell metadata.
#' @param ident The metadata column containing cell identities (clusters,
#'   cell types). If `NULL` for a Seurat object, the active identity is
#'   used and stored under the column name `"Identity"`. Required for
#'   Giotto and h5ad objects. This column forms the primary categorical
#'   axis of the plot — typically the x-axis, pie slices, or heatmap rows.
#' @param group_by The metadata column(s) used for secondary grouping
#'   (conditions, samples, stimulation status). Default is `NULL`.
#'   Behavior varies by plot type:
#'   \itemize{
#'     \item For most plots, multiple columns are concatenated into a
#'           single grouping variable (separated by `group_by_sep`).
#'     \item For `"sankey"` and `"heatmap"` plots, multiple columns are
#'           NOT concatenated — each column forms a separate node or
#'           column grouping.
#'     \item For `"violin"` and `"box"` plots, at most 2 columns are
#'           allowed: the first determines the x-axis breakdown, the second
#'           is passed as `group_by` to the underlying plot function.
#'   }
#' @param group_by_sep Separator used when concatenating multiple
#'   `group_by` columns into a single variable. Default is `"_"`. Ignored
#'   for `"sankey"` and `"heatmap"` plots where columns are not combined.
#' @param spat_unit Spatial unit name for Giotto objects (e.g., `"cell"`).
#'   Ignored for non-Giotto inputs. If `NULL`, auto-detected via
#'   `GiottoClass::set_default_spat_unit()`.
#' @param feat_type Feature type name for Giotto objects (e.g., `"rna"`,
#'   `"dna"`, `"protein"`). Ignored for non-Giotto inputs. If `NULL`,
#'   auto-detected via `GiottoClass::set_default_feat_type()`.
#' @param split_by Metadata column(s) used to split the data into separate
#'   plots. Each unique value (or combination) produces an independent
#'   plot. Multiple columns are concatenated using `split_by_sep`. Default
#'   is `NULL` (no splitting).
#' @param split_by_sep Separator used when concatenating multiple
#'   `split_by` columns. Default is `"_"`.
#' @param columns_split_by Metadata column used to split the columns of
#'   `"heatmap"` or `"pies"` plots. This adds an additional level of column
#'   faceting beyond what `group_by` provides. Default is `NULL`.
#' @param rows_by Metadata column(s) used as the rows for `"heatmap"` or
#'   `"pies"` plots. For `"pies"`, this defines what each pie chart
#'   represents (e.g., clones). For `"heatmap"`, this overrides the default
#'   row grouping (which would otherwise be `ident`). When multiple columns
#'   are provided, they must be logical or numeric. Default is `NULL`.
#' @param facet_by Metadata column(s) used to facet the plots (separate
#'   panels within the same output). Not available for `"circos"`,
#'   `"sankey"`, and `"heatmap"` plot types. Default is `NULL`.
#' @param plot_type The visualization type. One of:
#'   `"bar"` (default), `"circos"`, `"pie"`, `"pies"`, `"ring"`/`"donut"`,
#'   `"trend"`, `"area"`, `"sankey"`/`"alluvial"`, `"heatmap"`,
#'   `"radar"`, `"spider"`, `"violin"`, or `"box"`.
#'   See the \emph{Plot types} section for guidance on choosing a type.
#'   Note: `"donut"` is an alias for `"ring"`; `"alluvial"` is an alias
#'   for `"sankey"`.
#' @param agg An expression string passed to `dplyr::summarise()` to
#'   compute the value for each ident-by-group intersection. Default is
#'   `"n()"` (count of cells). For example,
#'   `"sum(hasTCR) / n()"` computes the fraction of TCR-positive cells.
#'   Ignored for `"circos"` and `"pies"` plot types. See the
#'   \emph{Custom aggregation} section for more examples.
#' @param frac The fraction normalization mode. One of `"none"` (default),
#'   `"group"`, `"ident"`, `"cluster"` (alias for `"ident"`), or `"all"`.
#'   Fractions are calculated within each `split_by` and `facet_by` subset.
#'   See the \emph{Fraction calculation} section for detailed semantics.
#' @param swap Whether to exchange the roles of `ident` and `group_by`.
#'   Default is `FALSE`. Behavior by plot type:
#'   \itemize{
#'     \item **`"bar"`, `"trend"`, `"area"`, `"ring"`, `"radar"`,
#'           `"spider"`:** Puts `group_by` on the x-axis and uses `ident`
#'           for fill/color.
#'     \item **`"circos"`:** Reverses the direction of edges (from ident
#'           to group_by instead of group_by to ident).
#'     \item **`"pies"`:** Swaps `group_by` and `ident` for the column
#'           grouping and pie grouping respectively.
#'     \item **`"heatmap"`:** Swaps `group_by` and `columns_split_by`
#'           (requires `columns_split_by` to be set).
#'     \item **`"violin"`, `"box"`:** Swaps the x-axis and `group_by`
#'           assignments from the two `group_by` columns.
#'   }
#'   Only works when `group_by` is specified. Note that `swap` is different
#'   from `flip` (which transposes coordinates).
#' @param name Display name for the main legend / value scale in
#'   `"heatmap"` and `"pies"` plots. Default is `NULL`, which auto-generates
#'   a name based on the fraction mode.
#' @param ylab Y-axis label. Default is `NULL`, which auto-generates a
#'   label based on the fraction mode (e.g., `"Number of cells"` or
#'   `"Fraction of cells"`).
#' @param ... Additional arguments passed to the underlying
#'   \pkg{plotthis} plotting function. The target function depends on
#'   `plot_type`:
#'   \itemize{
#'     \item `"bar"` — \code{\link[plotthis:BarPlot]{plotthis::BarPlot()}}
#'           (`position`, `palette`, `alpha`, `x_text_angle`, ...)
#'     \item `"circos"` — \code{\link[plotthis:CircosPlot]{plotthis::CircosPlot()}}
#'           (`labels_rot`, `links_alpha`, ...)
#'     \item `"pie"` — \code{\link[plotthis:PieChart]{plotthis::PieChart()}}
#'     \item `"ring"` / `"donut"` — \code{\link[plotthis:RingPlot]{plotthis::RingPlot()}}
#'           (`palette`, `alpha`, ...)
#'     \item `"trend"` — \code{\link[plotthis:TrendPlot]{plotthis::TrendPlot()}}
#'           (`palette`, `point_size`, `line_width`, ...)
#'     \item `"area"` — \code{\link[plotthis:AreaPlot]{plotthis::AreaPlot()}}
#'           (`palette`, `alpha`, ...)
#'     \item `"sankey"` / `"alluvial"` — \code{\link[plotthis:SankeyPlot]{plotthis::SankeyPlot()}}
#'           (`links_alpha`, `links_fill_by`, ...)
#'     \item `"heatmap"` — \code{\link[plotthis:Heatmap]{plotthis::Heatmap()}}
#'           (`palette`, `cluster_rows`, `cluster_columns`, `show_row_names`,
#'           `show_column_names`, `cell_type`, ...)
#'     \item `"pies"` — \code{\link[plotthis:Heatmap]{plotthis::Heatmap()}}
#'           (with `cell_type = "pie"`; `pie_size`, `pie_values`,
#'           `row_names_side`, `column_names_side`, ...)
#'     \item `"radar"` — \code{\link[plotthis:RadarPlot]{plotthis::RadarPlot()}}
#'     \item `"spider"` — \code{\link[plotthis:SpiderPlot]{plotthis::SpiderPlot()}}
#'     \item `"violin"` — \code{\link[plotthis:ViolinPlot]{plotthis::ViolinPlot()}}
#'           (`add_box`, `comparisons`, `aspect.ratio`, ...)
#'     \item `"box"` — \code{\link[plotthis:BoxPlot]{plotthis::BoxPlot()}}
#'           (`comparisons`, `aspect.ratio`, ...)
#'   }
#'   Common layout parameters (`split_by`, `facet_by`, `combine`, `nrow`,
#'   `ncol`) are handled by `CellStatPlot` directly.
#' @return A `ggplot` object, or a list of `ggplot` objects if
#'   `combine = FALSE` is passed via `...`.
#' @note
#' **Default y-axis label:** When `ylab = NULL`, the label is
#' auto-generated as `"Number of cells"` (when `frac = "none"`) or
#' `"Fraction of cells"` (otherwise). Override with `ylab` for custom
#' labels.
#'
#' **Heatmap/pies metadata columns:** When multiple `rows_by` columns are
#' specified for `"heatmap"` or `"pies"` plots, each column must be logical
#' or numeric. The function filters cells where the column is `TRUE` or
#' greater than zero.
#'
#' **Facet restrictions:** `facet_by` is not supported for `"circos"`,
#' `"sankey"`, and `"heatmap"` plot types. For heatmaps, use
#' `split_by` to create separate plots.
#'
#' **Factor ordering:** The order of categories in plots follows the factor
#' levels of the metadata columns. Set factor levels before plotting to
#' control the display order of cell types and groups.
#'
#' **Giotto spatial units and feature types:** The `spat_unit` and
#' `feat_type` parameters are required to locate the correct metadata
#' within Giotto's hierarchical spatial data structure. When `NULL`,
#' Giotto's own default resolution is used.
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{CellDimPlot}} — Dimension reduction visualization of
#'         cell clusters
#'   \item \code{\link{FeatureStatPlot}} — Feature expression statistics
#'   \item \code{\link[plotthis:BarPlot]{plotthis::BarPlot()}},
#'         \code{\link[plotthis:Heatmap]{plotthis::Heatmap()}},
#'         \code{\link[plotthis:SankeyPlot]{plotthis::SankeyPlot()}} —
#'         Underlying \pkg{plotthis} functions for individual plot types
#' }
#'
#' @export
#' @importFrom rlang sym syms parse_expr
#' @importFrom SeuratObject Idents
#' @importFrom dplyr %>% summarise mutate ungroup n
#' @importFrom tidyr drop_na pivot_wider pivot_longer
#' @importFrom plotthis BarPlot CircosPlot PieChart RingPlot TrendPlot AreaPlot SankeyPlot Heatmap RadarPlot SpiderPlot ViolinPlot BoxPlot
#' @details See:
#' * <https://pwwang.github.io/scplotter/articles/Giotto_Xenium.html>
#'
#' for examples of using this function with a Giotto object.
#'
#' And see:
#' * <https://pwwang.github.io/scplotter/articles/Working_with_anndata_h5ad_files.html>
#'
#' for examples of using this function with .h5ad files.
#' @examples
#' \donttest{
#' # library(patchwork)
#' data(ifnb_sub)
#'
#' # Bar plot
#' p1 <- CellStatPlot(ifnb_sub)
#' p2 <- CellStatPlot(ifnb_sub, group_by = "stim")
#' p3 <- CellStatPlot(ifnb_sub, group_by = "stim", palette = "Set2")
#' p4 <- CellStatPlot(ifnb_sub, group_by = "stim", position = "stack")
#' (p1 | p2) / (p3 | p4)
#'
#' # Fraction of cells
#' p1 <- CellStatPlot(ifnb_sub, group_by = "stim", frac = "group",
#'                    ident = "seurat_annotations", x_text_angle = 60)
#' p2 <- CellStatPlot(ifnb_sub, group_by = "stim", frac = "group")
#' p3 <- CellStatPlot(ifnb_sub, group_by = "stim", frac = "ident",
#'                    position = "stack", alpha = .6)
#' p4 <- CellStatPlot(ifnb_sub, group_by = "stim", frac = "group",
#'                    swap = TRUE, position = "stack")
#' (p1 | p2) / (p3 | p4)
#'
#' # Splitting/Facetting the plot
#' CellStatPlot(ifnb_sub, split_by = "stim")
#' CellStatPlot(ifnb_sub, facet_by = "stim")
#' CellStatPlot(ifnb_sub, facet_by = "stim", facet_nrow = 2)
#'
#' # Circos plot
#' CellStatPlot(ifnb_sub, group_by = "stim", plot_type = "circos")
#' CellStatPlot(ifnb_sub, group_by = "stim", ident = "seurat_annotations",
#'              plot_type = "circos", labels_rot = TRUE)
#'
#' # Pie plot
#' CellStatPlot(ifnb_sub, plot_type = "pie")
#' CellStatPlot(ifnb_sub, plot_type = "pie", split_by = "stim")
#'
#' # Ring plot
#' CellStatPlot(ifnb_sub, plot_type = "ring", group_by = "stim",
#'              palette = "Spectral")
#'
#' # Trend plot
#' CellStatPlot(ifnb_sub, plot_type = "trend", frac = "group",
#'              x_text_angle = 90, group_by = c("stim", "seurat_annotations"))
#'
#' # Sankey plot
#' CellStatPlot(ifnb_sub, plot_type = "sankey", group_by = c("seurat_clusters", "stim"),
#'              links_alpha = .6)
#' CellStatPlot(ifnb_sub, plot_type = "sankey", links_alpha = .6,
#'              group_by = c("stim", "seurat_annotations", "orig.ident"))
#' CellStatPlot(ifnb_sub, plot_type = "sankey", links_alpha = .6,
#'              group_by = c("seurat_clusters", "stim", "seurat_annotations", "orig.ident"))
#'
#' # Area plot
#' CellStatPlot(ifnb_sub, plot_type = "area", frac = "group", x_text_angle = 90,
#'              group_by = "seurat_annotations", split_by = "stim")
#'
#' # Heatmap
#' CellStatPlot(ifnb_sub, plot_type = "heatmap", group_by = "stim", palette = "Blues")
#' CellStatPlot(ifnb_sub, plot_type = "heatmap", group_by = "stim",
#'    frac = "group", columns_split_by = "seurat_annotations", swap = TRUE)
#'
#' # Pies
#' # Simulate some sets of cells (e.g. clones)
#' ifnb_sub$r1 <- ifelse(ifnb_sub$seurat_clusters %in% c("0", "1", "2"), 1, 0)
#' ifnb_sub$r2 <- sample(c(1, 0), ncol(ifnb_sub), prob = c(0.5, 0.5), replace = TRUE)
#' ifnb_sub$r3 <- sample(c(1, 0), ncol(ifnb_sub), prob = c(0.7, 0.3), replace = TRUE)
#' CellStatPlot(ifnb_sub, plot_type = "pies", group_by = "stim", rows_name = "Clones",
#'    rows_by = c("r1", "r2", "r3"), column_names_side = "top", cluster_columns = FALSE,
#'    row_names_side = "right", pie_size = "identity", pie_values = "sum")
#'
#' # Expand pies into heatmap
#' CellStatPlot(ifnb_sub, plot_type = "heatmap", group_by = "stim", rows_name = "Clones",
#'    rows_by = c("r1", "r2", "r3"), column_names_side = "top", cluster_columns = FALSE,
#'    row_names_side = "right")
#'
#' # Rows split by clones instead of idents and show fraction of cells in each clone.
#' CellStatPlot(ifnb_sub, plot_type = "heatmap", group_by = "stim", frac = "group",
#'   rows_split_by = c("r1", "r2", "r3"), column_names_side = "top", cluster_columns = FALSE,
#'   row_names_side = "right", label = function(x) scales::number(x, accuracy = 0.01),
#'   cell_type = "label")
#'
#' # Radar plot/Spider plot
#' pr <- CellStatPlot(ifnb_sub, plot_type = "radar", group_by = "stim")
#' ps <- CellStatPlot(ifnb_sub, plot_type = "spider", group_by = "stim")
#' pr | ps
#'
#' # Box/Violin plot
#' ifnb_sub$group <- sample(paste0("g", 1:10), nrow(ifnb_sub), replace = TRUE)
#' CellStatPlot(ifnb_sub, group_by = c("group", "stim"), frac = "group",
#'    plot_type = "violin", add_box = TRUE, ident = "seurat_annotations",
#'    x_text_angle = 60, comparisons = TRUE, aspect.ratio = 0.8)
#'
#' # Use different agg other than counting the number of cells.
#' # Let's say we do the fraction of g1 in each stim group.
#' CellStatPlot(ifnb_sub, agg = "sum(group == 'g1') / n()",
#'    plot_type = "bar", ylab = "Fraction of g1 cells")
#' CellStatPlot(ifnb_sub, group_by = "stim", agg = "sum(group == 'g1') / n()",
#'    plot_type = "bar", ylab = "Fraction of g1 cells")
#' }
CellStatPlot <- function(
    object, ident = NULL, group_by = NULL, group_by_sep = "_", spat_unit = NULL, feat_type = NULL,
    split_by = NULL, split_by_sep = "_", facet_by = NULL, rows_by = NULL, columns_split_by = NULL,
    frac = c("none", "group", "ident", "cluster", "all"), name = NULL, agg = "n()",
    plot_type = c("bar", "circos", "pie", "pies", "ring", "donut", "trend", "area", "sankey", "alluvial", "heatmap", "radar", "spider", "violin", "box"),
    swap = FALSE, ylab = NULL, ...
) {
    UseMethod("CellStatPlot")
}

#' @export
CellStatPlot.giotto <- function(
    object, ident = NULL, group_by = NULL, group_by_sep = "_", spat_unit = NULL, feat_type = NULL,
    split_by = NULL, split_by_sep = "_", facet_by = NULL, rows_by = NULL, columns_split_by = NULL,
    frac = c("none", "group", "ident", "cluster", "all"), name = NULL, agg = "n()",
    plot_type = c("bar", "circos", "pie", "pies", "ring", "donut", "trend", "area", "sankey", "alluvial", "heatmap", "radar", "spider", "violin", "box"),
    swap = FALSE, ylab = NULL, ...
) {
    stopifnot("[CellStatPlot] 'ident' is required for Giotto object." = !is.null(ident))

    spat_unit <- GiottoClass::set_default_spat_unit(
        gobject = object,
        spat_unit = spat_unit
    )

    feat_type <- GiottoClass::set_default_feat_type(
        gobject = object,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    data <- GiottoClass::getCellMetadata(
        gobject = object,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "data.table"
    )

    CellStatPlot.data.frame(
        data, ident = ident, group_by = group_by, group_by_sep = group_by_sep,
        split_by = split_by, split_by_sep = split_by_sep, facet_by = facet_by,
        rows_by = rows_by, columns_split_by = columns_split_by, frac = frac,
        name = name, plot_type = plot_type, agg = agg,
        swap = swap, ylab = ylab, ...
    )
}

#' @export
CellStatPlot.Seurat <- function(
    object, ident = NULL, group_by = NULL, group_by_sep = "_", spat_unit = NULL, feat_type = NULL,
    split_by = NULL, split_by_sep = "_", facet_by = NULL, rows_by = NULL, columns_split_by = NULL,
    frac = c("none", "group", "ident", "cluster", "all"), name = NULL, agg = "n()",
    plot_type = c("bar", "circos", "pie", "pies", "ring", "donut", "trend", "area", "sankey", "alluvial", "heatmap", "radar", "spider", "violin", "box"),
    swap = FALSE, ylab = NULL, ...
) {
    data <- object@meta.data
    if (is.null(ident)) {
        ident <- "Identity"
        data[[ident]] <- Idents(object)
    }

    CellStatPlot.data.frame(
        data, ident = ident, group_by = group_by, group_by_sep = group_by_sep,
        split_by = split_by, split_by_sep = split_by_sep, facet_by = facet_by,
        rows_by = rows_by, columns_split_by = columns_split_by, frac = frac,
        name = name, plot_type = plot_type, agg = agg,
        swap = swap, ylab = ylab, ...
    )
}

#' @export
CellStatPlot.character <- function(
    object, ident = NULL, group_by = NULL, group_by_sep = "_", spat_unit = NULL, feat_type = NULL,
    split_by = NULL, split_by_sep = "_", facet_by = NULL, rows_by = NULL, columns_split_by = NULL,
    frac = c("none", "group", "ident", "cluster", "all"), name = NULL, agg = "n()",
    plot_type = c("bar", "circos", "pie", "pies", "ring", "donut", "trend", "area", "sankey", "alluvial", "heatmap", "radar", "spider", "violin", "box"),
    swap = FALSE, ylab = NULL, ...
) {
    if (!endsWith(object, ".h5ad")) {
        stop("[CellStatPlot] Currently only supports .h5ad files when called with a string/path.")
    }

    object <- hdf5r::H5File$new(object, mode = "r")
    on.exit(object$close_all())

    CellStatPlot.H5File(
        object, ident = ident, group_by = group_by, group_by_sep = group_by_sep,
        spat_unit = spat_unit, feat_type = feat_type, agg = agg,
        split_by = split_by, split_by_sep = split_by_sep, facet_by = facet_by,
        rows_by = rows_by, columns_split_by = columns_split_by, frac = frac,
        name = name, plot_type = plot_type,
        swap = swap, ylab = ylab, ...
    )
}

#' @export
CellStatPlot.H5File <- function(
    object, ident = NULL, group_by = NULL, group_by_sep = "_", spat_unit = NULL, feat_type = NULL,
    split_by = NULL, split_by_sep = "_", facet_by = NULL, rows_by = NULL, columns_split_by = NULL,
    frac = c("none", "group", "ident", "cluster", "all"), name = NULL, agg = "n()",
    plot_type = c("bar", "circos", "pie", "pies", "ring", "donut", "trend", "area", "sankey", "alluvial", "heatmap", "radar", "spider", "violin", "box"),
    swap = FALSE, ylab = NULL, ...
) {
    stopifnot("[CellStatPlot] 'ident' is required for anndata (h5ad) object." = !is.null(ident))

    object <- h5group_to_dataframe(object[["obs"]])

    CellStatPlot.data.frame(
        object, ident = ident, group_by = group_by, group_by_sep = group_by_sep,
        spat_unit = spat_unit, feat_type = feat_type, agg = agg,
        split_by = split_by, split_by_sep = split_by_sep, facet_by = facet_by,
        rows_by = rows_by, columns_split_by = columns_split_by, frac = frac,
        name = name, plot_type = plot_type,
        swap = swap, ylab = ylab, ...
    )
}

#' @export
CellStatPlot.data.frame <- function(
    object, ident = NULL, group_by = NULL, group_by_sep = "_", spat_unit = NULL, feat_type = NULL,
    split_by = NULL, split_by_sep = "_", facet_by = NULL, rows_by = NULL, columns_split_by = NULL,
    frac = c("none", "group", "ident", "cluster", "all"), name = NULL, agg = "n()",
    plot_type = c("bar", "circos", "pie", "pies", "ring", "donut", "trend", "area", "sankey", "alluvial", "heatmap", "radar", "spider", "violin", "box"),
    swap = FALSE, ylab = NULL, ...
) {
    plot_type <- match.arg(plot_type)
    if (plot_type == "donut") plot_type <- "ring"
    if (plot_type == "alluvial") plot_type <- "sankey"
    if (isFALSE(swap) && plot_type %in% c("sankey", "heatmap", "violin", "box")) {
        group_by <- check_columns(object, group_by, force_factor = TRUE,
            allow_multi = TRUE)
    } else {
        group_by <- check_columns(object, group_by, force_factor = TRUE,
            allow_multi = TRUE, concat_multi = TRUE, concat_sep = group_by_sep)
    }
    split_by <- check_columns(object, split_by, force_factor = TRUE,
        allow_multi = TRUE, concat_multi = TRUE, concat_sep = split_by_sep)
    facet_by <- check_columns(object, facet_by, force_factor = TRUE,
        allow_multi = TRUE)
    rows_by <- check_columns(object, rows_by, allow_multi = TRUE)
    rows_split_by <- list(...)[["rows_split_by"]]
    rows_split_by <- check_columns(object, rows_split_by, allow_multi = TRUE)
    stopifnot("[CellStatPlot] Can't have both 'rows_by' and 'rows_split_by' for plot_type = 'heatmap'" = !(plot_type == "heatmap" && !is.null(rows_by) && !is.null(rows_split_by)))
    stopifnot("[CellStatPlot] Can't have a single 'rows_by' column as the heatmap rows for plot_type will be idents. Do you mean 'rows_split_by'?" = !(plot_type == "heatmap" && length(rows_by) == 1))

    frac <- match.arg(frac)
    if (frac == "cluster") frac <- "ident"
    if (identical(plot_type, "ring") && !identical(frac, "group")) {
        message("'frac' is forced to 'group' for 'ring' plot.")
        frac <- "group"
    }
    if (is.null(group_by) && frac == "ident") {
        stop("Cannot calculate the fraction by 'ident' without specifying 'group_by'.")
    }

    if (isTRUE(swap) && is.null(group_by)) {
        stop("Cannot swap the 'ident' and 'group_by' without specifying 'group_by'.")
    }

    object <- object %>% drop_na(!!sym(ident))
    if (nrow(object) == 0) {
        stop("No cells found with ident:'", ident, "'")
    }

    if (plot_type != "pies") {
        # Heatmap itself will handle the data for cell_type = pie
        if (is.null(group_by)) {
            object <- object %>%
                dplyr::group_by(!!!syms(unique(c(split_by, facet_by, columns_split_by, ident)))) %>%
                summarise(.n = !!parse_expr(agg), .groups = "drop")

            if (frac != "none") {
                object <- object %>%
                    dplyr::group_by(!!!syms(unique(c(split_by, facet_by, columns_split_by)))) %>%
                    mutate(.frac = !!sym(".n") / sum(!!sym(".n")))
            } else {
                object <- object %>% mutate(.frac = 1)  # not used
            }
        } else if (length(group_by) > 1 && !plot_type %in% c("sankey", "violin", "box")) {
            # calculate the frac for each group. we don't want them to be concatenated.
            object <- do_call(rbind, lapply(group_by, function(g) {
                dat <- object %>%
                    dplyr::group_by(!!!syms(unique(c(split_by, facet_by, columns_split_by, ident, g)))) %>%
                    summarise(.n = !!parse_expr(agg), .groups = "drop")
                dat <- dat[!is.na(dat[[g]]), , drop = FALSE]

                if (frac == "group") {
                    dat <- dat %>%
                        dplyr::group_by(!!!syms(unique(c(split_by, facet_by, g, columns_split_by)))) %>%
                        mutate(.frac = !!sym(".n") / sum(!!sym(".n"))) %>%
                        ungroup()
                } else if (frac == "ident") {
                    dat <- dat %>%
                        dplyr::group_by(!!!syms(unique(c(split_by, facet_by, columns_split_by, ident)))) %>%
                        mutate(.frac = !!sym(".n") / sum(!!sym(".n"))) %>%
                        ungroup()
                } else if (frac == "all") {
                    dat <- dat %>%
                        dplyr::group_by(!!!syms(unique(c(split_by, facet_by, columns_split_by)))) %>%
                        mutate(.frac = !!sym(".n") / sum(!!sym(".n"))) %>%
                        ungroup()
                } else {
                    dat <- dat %>% mutate(.frac = 1)  # not used
                }

                dat[, setdiff(group_by, g)] <- NA
                dat
            }))
        } else if (plot_type == "heatmap" && (length(c(rows_by, rows_split_by)) > 1)) {
            is_rows_by <- length(rows_by) > 0
            row_cols <- if (is_rows_by) rows_by else rows_split_by
            # if multiple columns provided, they must be logical or numeric columns
            rc_obj <- NULL
            for (col in row_cols) {
                col_agg <- if (is.list(agg)) agg[[col]] else agg
                col_agg <- col_agg %||% "n()"
                if (!is.logical(object[[col]]) && !is.numeric(object[[col]])) {
                    stop("[CellStatPlot] For 'heatmap' plot, multiple columns specified in 'rows_by' and 'rows_split_by' must be logical or numeric. But column '", col, "' is not.")
                }
                tmp <- object %>% dplyr::filter(!is.na(!!sym(col)) & (!!sym(col) == TRUE | !!sym(col) > 0)) %>%
                    dplyr::group_by(!!!syms(unique(c(split_by, facet_by, group_by, columns_split_by, ident)))) %>%
                    summarise(.n = !!parse_expr(col_agg), .groups = "drop")

                tmp[[".rows_by_or_rows_split_by"]] <- as.character(col)
                if (frac == "group") {
                    tmp <- tmp %>%
                        dplyr::group_by(!!!syms(unique(c(split_by, facet_by, group_by, columns_split_by)))) %>%
                        mutate(.frac = !!sym(".n") / sum(!!sym(".n"))) %>%
                        ungroup()
                } else if (frac == "ident") {
                    tmp <- tmp %>%
                        dplyr::group_by(!!!syms(unique(c(split_by, facet_by, columns_split_by, ident)))) %>%
                        mutate(.frac = !!sym(".n") / sum(!!sym(".n"))) %>%
                        ungroup()
                } else if (frac == "all") {
                    tmp <- tmp %>%
                        dplyr::group_by(!!!syms(unique(c(split_by, facet_by, columns_split_by)))) %>%
                        mutate(.frac = !!sym(".n") / sum(!!sym(".n"))) %>%
                        ungroup()
                } else {
                    tmp <- tmp %>% mutate(.frac = 1)  # not used
                }
                rc_obj <- rbind(rc_obj, tmp)
            }
            rc_obj[[".rows_by_or_rows_split_by"]] <- factor(rc_obj[[".rows_by_or_rows_split_by"]], levels = row_cols)
            object <- rc_obj
            rm(rc_obj)
            if (is_rows_by) {
                rows_by <- ".rows_by_or_rows_split_by"
            } else {
                rows_split_by <- ".rows_by_or_rows_split_by"
            }
        } else {
            row_cols <- if (plot_type == "heatmap") c(rows_by, rows_split_by) else NULL
            object <- object %>%
                dplyr::group_by(!!!syms(unique(c(split_by, facet_by, group_by, row_cols, columns_split_by, ident)))) %>%
                summarise(.n = !!parse_expr(agg), .groups = "drop")

            if (frac == "group") {
                object <- object %>%
                    dplyr::group_by(!!!syms(unique(c(split_by, facet_by, group_by, row_cols, columns_split_by)))) %>%
                    mutate(.frac = !!sym(".n") / sum(!!sym(".n"))) %>%
                    ungroup()
            } else if (frac == "ident") {
                object <- object %>%
                    dplyr::group_by(!!!syms(unique(c(split_by, facet_by, row_cols, columns_split_by, ident)))) %>%
                    mutate(.frac = !!sym(".n") / sum(!!sym(".n"))) %>%
                    ungroup()
            } else if (frac == "all") {
                object <- object %>%
                    dplyr::group_by(!!!syms(unique(c(split_by, facet_by, row_cols, columns_split_by)))) %>%
                    mutate(.frac = !!sym(".n") / sum(!!sym(".n"))) %>%
                    ungroup()
            } else {
                object <- object %>% mutate(.frac = 1)  # not used
            }
        }
    }

    ylab <- ylab %||% paste0(ifelse(identical(frac, "none"), "Number", "Fraction"), " of cells")
    if (plot_type == "bar") {
        if (is.null(group_by)) {
            BarPlot(
                object, x = ident, y = ifelse(identical(frac, "none"), ".n", ".frac"),
                split_by = split_by, facet_by = facet_by, ylab = ylab, ...)
        } else {
            BarPlot(
                object,
                x = ifelse(swap, group_by, ident),
                y = ifelse(identical(frac, "none"), ".n", ".frac"),
                group_by = ifelse(swap, ident, group_by),
                split_by = split_by, facet_by = facet_by, ylab = ylab, ...)
        }
    } else if (plot_type == "circos") {
        if (is.null(group_by)) {
            stop("Cannot create a circos plot without specifying 'group_by'.")
        }
        CircosPlot(
            object,
            from = ifelse(swap, ident, group_by),
            to = ifelse(swap, group_by, ident),
            split_by = split_by, facet_by = facet_by, ...)
    } else if (plot_type == "pie") {
        if (!is.null(group_by)) {
            stop("Cannot create a pie plot with 'group_by'. You may want to use 'ring' plot instead.")
        }
        PieChart(
            object, x = ident, y = ifelse(identical(frac, "none"), ".n", ".frac"),
            split_by = split_by, facet_by = facet_by, ...)
    } else if (plot_type == "ring") {
        if (!identical(frac, "group")) {
            stop("'frac' must be 'group' for 'ring' plot.")
        }
        RingPlot(
            object,
            x = ifelse(swap, ident, group_by),
            y = ifelse(identical(frac, "none"), ".n", ".frac"),
            group_by = ifelse(swap, group_by, ident),
            split_by = split_by, facet_by = facet_by, ...)
    } else if (plot_type == "trend") {
        if (is.null(group_by)) {
            stop("Cannot create a trend plot without specifying 'group_by'.")
        }
        TrendPlot(
            object,
            x = ifelse(swap, ident, group_by),
            y = ifelse(identical(frac, "none"), ".n", ".frac"),
            group_by = ifelse(swap, group_by, ident),
            ylab = ylab, split_by = split_by, facet_by = facet_by, ...)
    } else if (plot_type == "area") {
        if (is.null(group_by)) {
            stop("Cannot create an area plot without specifying 'group_by'.")
        }
        AreaPlot(
            object,
            x = ifelse(swap, ident, group_by),
            y = ifelse(identical(frac, "none"), ".n", ".frac"),
            group_by = ifelse(swap, group_by, ident),
            ylab = ylab, split_by = split_by, facet_by = facet_by, ...)
    } else if (plot_type == "sankey") {
        if (is.null(group_by)) {
            stop("Cannot create a sankey plot without specifying 'group_by'.")
        }
        if (frac == "ident") {
            stop("Cannot calculate the fraction by 'ident' for 'sankey' plot.")
        }
        if (isTRUE(swap)) {
            stop("'swap = TRUE' is not supported for 'sankey' plot.")
        }
        SankeyPlot(
            object, x = group_by, links_fill_by = ident, xlab = "", ylab = "", flow = TRUE,
            y = ifelse(identical(frac, "none"), ".n", ".frac"),
            split_by = split_by, facet_by = facet_by, ...)
    } else if (plot_type == "pies") {
        if (is.null(group_by)) {
            stop("Cannot create a heatmap (cell_type = 'pie') plot without specifying 'group_by' (works as columns of the plot).")
        }
        if (is.null(rows_by)) {
            stop("Cannot create a heatmap (cell_type = 'pie') plot without specifying 'rows_by' (works as rows of the plot).")
        }
        if (!is.null(facet_by)) {
            stop("Cannot create a heatmap (cell_type = 'pie') plot with 'facet_by'.")
        }
        args <- rlang::dots_list(...)
        args$data <- object
        args$rows_by <- rows_by
        args$cell_type <- "pie"
        args$values_by <- args$values_by %||% name
        args$columns_by <- if (swap) ident else group_by
        args$values_fill <- args$values_fill %||% 0
        args$pie_group_by <- if (swap) group_by else ident
        args$columns_split_by <- columns_split_by
        args$split_by <- split_by
        args$show_row_names <- args$show_row_names %||% TRUE
        args$show_column_names <- args$show_column_names %||% TRUE
        args$add_reticle <- args$add_reticle %||% TRUE

        do_call(Heatmap, args)
    } else if (plot_type == "heatmap") {
        if (is.null(group_by)) {
            stop("Cannot create a heatmap plot without specifying 'group_by', which should work as the columns of the heatmap.")
        }
        if (!is.null(facet_by)) {
            stop("Cannot create a heatmap plot with 'facet_by'.")
        }
        if (swap && is.null(columns_split_by)) {
            stop("Cannot swap between 'group_by' and 'columns_split_by' without specifying 'columns_split_by'.")
        }

        args <- rlang::dots_list(...)
        args$data <- object
        args$name <- args$name %||% ifelse(identical(frac, "none"), "Number of cells", "Fraction of cells")
        args$values_by <- ifelse(identical(frac, "none"), ".n", ".frac")
        args$in_form <- "long"
        args$columns_by <- if (swap) columns_split_by else group_by
        args$values_fill <- args$values_fill %||% 0
        args$columns_split_by <- if (swap) group_by else columns_split_by
        args$split_by <- split_by
        args$show_row_names <- args$show_row_names %||% TRUE
        args$show_column_names <- args$show_column_names %||% TRUE
        if (!is.null(rows_by)) {
            args$rows_split_by <- ident
            args$rows_by <- rows_by
        } else if (!is.null(rows_split_by)) {
            args$rows_split_by <- rows_split_by
            args$rows_split_name <- args$rows_split_name %||% " "
            args$rows_by <- ident
        } else {
            args$rows_by <- ident
        }

        do_call(Heatmap, args)
    } else if (plot_type %in% c("violin", "box")) {
        if (is.null(group_by)) {
            stop("Cannot create a 'violin'/'box' plot without specifying 'group_by'.")
        }
        if (length(group_by) > 2) {
            stop("Cannot create a 'violin'/'box' plot with more than 2 'group_by'.")
        }
        fn <- ifelse(plot_type == "violin", ViolinPlot, BoxPlot)
        if (length(group_by) == 1) {
            fn(
                object,
                x = ifelse(swap, group_by, ident),
                y = ifelse(identical(frac, "none"), ".n", ".frac"),
                in_form = "long",
                ylab = ylab, split_by = split_by, facet_by = facet_by, ...)
        } else {
            fn(
                object,
                x = ifelse(swap, group_by[2], ident),
                y = ifelse(identical(frac, "none"), ".n", ".frac"),
                group_by = ifelse(swap, ident, group_by[2]),
                ylab = ylab, split_by = split_by, facet_by = facet_by, ...)
        }
    } else {  # radar/spider plot
        if (is.null(group_by)) {
            stop("Cannot create a '", plot_type, "' plot without specifying 'group_by'.")
        }
        if (plot_type == "radar") {
            RadarPlot(
                object,
                x = ifelse(swap, group_by, ident),
                y = ifelse(identical(frac, "none"), ".n", ".frac"),
                group_by = ifelse(swap, ident, group_by),
                split_by = split_by, facet_by = facet_by, ...)
        } else if (plot_type == "spider") {
            SpiderPlot(
                object,
                x = ifelse(swap, group_by, ident),
                y = ifelse(identical(frac, "none"), ".n", ".frac"),
                group_by = ifelse(swap, ident, group_by),
                split_by = split_by, facet_by = facet_by, ...)
        }
    }
}
