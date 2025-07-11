% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/featurestatplot.R
\name{FeatureStatPlot}
\alias{FeatureStatPlot}
\title{Feature statistic plot}
\usage{
FeatureStatPlot(
  object,
  features,
  plot_type = c("violin", "box", "bar", "ridge", "dim", "cor", "heatmap", "dot"),
  spat_unit = NULL,
  feat_type = NULL,
  reduction = NULL,
  graph = NULL,
  bg_cutoff = 0,
  dims = 1:2,
  rows_name = "Features",
  ident = NULL,
  assay = NULL,
  layer = NULL,
  agg = mean,
  group_by = NULL,
  split_by = NULL,
  facet_by = NULL,
  xlab = NULL,
  ylab = NULL,
  x_text_angle = NULL,
  ...
)
}
\arguments{
\item{object}{A seurat object, a giotto object, a path to an .h5ad file or an opened \code{H5File} by \code{hdf5r} package.}

\item{features}{A character vector of feature names}

\item{plot_type}{Type of the plot. It can be "violin", "box", "bar", "ridge", "dim", "cor", "heatmap" or "dot"}

\item{spat_unit}{The spatial unit to use for the plot. Only applied to Giotto objects.}

\item{feat_type}{feature type of the features (e.g. "rna", "dna", "protein"), only applied to Giotto objects.}

\item{reduction}{Name of the reduction to plot (for example, "umap"), only used when \code{plot_type} is "dim" or you can to use the reduction as feature.}

\item{graph}{Specify the graph name to add edges between cell neighbors to the plot, only used when \code{plot_type} is "dim".}

\item{bg_cutoff}{Background cutoff for the dim plot, only used when \code{plot_type} is "dim".}

\item{dims}{Dimensions to plot, only used when \code{plot_type} is "dim".}

\item{rows_name}{The name of the rows in the heatmap, only used when \code{plot_type} is "heatmap".}

\item{ident}{The column name in the meta data to identify the cells.}

\item{assay}{The assay to use for the feature data.}

\item{layer}{The layer to use for the feature data.}

\item{agg}{The aggregation function to use for the bar plot.}

\item{group_by}{The column name in the meta data to group the cells.}

\item{split_by}{Column name in the meta data to split the cells to different plots.
If TRUE, the cells are split by the features.}

\item{facet_by}{Column name in the meta data to facet the plots. Should be always NULL.}

\item{xlab}{The x-axis label.}

\item{ylab}{The y-axis label.}

\item{x_text_angle}{The angle of the x-axis text. Only used when \code{plot_type} is "violin", "bar", or "box".}

\item{...}{Other arguments passed to the plot functions.
\itemize{
\item For \code{plot_type} "violin", the arguments are passed to \code{\link[plotthis:boxviolinplot]{plotthis::ViolinPlot()}}.
\item For \code{plot_type} "box", the arguments are passed to \code{\link[plotthis:boxviolinplot]{plotthis::BoxPlot()}}.
\item For \code{plot_type} "bar", the arguments are passed to \code{\link[plotthis:barplot]{plotthis::BarPlot()}}.
\item For \code{plot_type} "ridge", the arguments are passed to \code{\link[plotthis:RidgePlot]{plotthis::RidgePlot()}}.
\item For \code{plot_type} "dim", the arguments are passed to \code{\link[plotthis:dimplot]{plotthis::FeatureDimPlot()}}.
\item For \code{plot_type} "heatmap", the arguments are passed to \code{\link[plotthis:Heatmap]{plotthis::Heatmap()}}.
\item For \code{plot_type} "cor" with 2 features, the arguments are passed to \code{\link[plotthis:CorPlot]{plotthis::CorPlot()}}.
\item For \code{plot_type} "cor" with more than 2 features, the arguments are passed to \code{\link[plotthis:CorPairsPlot]{plotthis::CorPairsPlot()}}.
\item For \code{plot_type} "dot", the arguments are passed to \code{\link[plotthis:Heatmap]{plotthis::Heatmap()}} with \code{cell_type} set to "dot".
}}
}
\value{
A ggplot object or a list if \code{combine} is FALSE
}
\description{
This function creates various types of feature statistic plots for a Seurat object, a Giotto object,
a path to an .h5ad file or an opened \code{H5File} by \code{hdf5r} package.
It allows for plotting features such as gene expression, scores, or other metadata across different groups or conditions.
The function supports multiple plot types including violin, box, bar, ridge, dimension reduction, correlation, heatmap, and dot plots.
It can also handle multiple features and supports faceting, splitting, and grouping by metadata columns.
}
\details{
See:
\itemize{
\item \url{https://pwwang.github.io/scplotter/articles/Giotto_Visium.html}
\item \url{https://pwwang.github.io/scplotter/articles/Giotto_Xenium.html}
}

for examples of using this function with Giotto objects.

And see:
\itemize{
\item \url{https://pwwang.github.io/scplotter/articles/Working_with_anndata_h5ad_files.html}
}

for examples of using this function with .h5ad files.
}
\examples{
\donttest{
data(pancreas_sub)

FeatureStatPlot(pancreas_sub, features = c("G2M_score", "nCount_RNA"),
   ident = "SubCellType", facet_scales = "free_y")
FeatureStatPlot(pancreas_sub, features = c("G2M_score", "nCount_RNA"),
   ident = "SubCellType", plot_type = "box", facet_scales = "free_y")
FeatureStatPlot(pancreas_sub, features = c("G2M_score", "nCount_RNA"),
   ident = "SubCellType", plot_type = "bar", facet_scales = "free_y")
FeatureStatPlot(pancreas_sub, features = c("G2M_score", "nCount_RNA"),
   ident = "SubCellType", plot_type = "ridge", flip = TRUE, facet_scales = "free_y")
FeatureStatPlot(pancreas_sub, features = c("G2M_score", "nCount_RNA"),
   ident = "SubCellType", facet_scales = "free_y", add_point = TRUE)
FeatureStatPlot(pancreas_sub, features = c("G2M_score", "nCount_RNA"),
   ident = "SubCellType", facet_scales = "free_y", add_trend = TRUE)
FeatureStatPlot(pancreas_sub, features = c("G2M_score", "nCount_RNA"),
   ident = "SubCellType", facet_scales = "free_y", add_stat = mean)
FeatureStatPlot(pancreas_sub, features = c("G2M_score", "nCount_RNA"),
   ident = "SubCellType", facet_scales = "free_y", group_by = "Phase")
FeatureStatPlot(
   subset(pancreas_sub,
       subset = SubCellType \%in\% c("Ductal", "Ngn3 low EP", "Ngn3 high EP")),
   features = c("G2M_score"),
   ident = "SubCellType", group_by = "Phase", comparisons = TRUE)
FeatureStatPlot(pancreas_sub, features = c("Rbp4", "Pyy"), ident = "SubCellType",
   add_bg = TRUE, add_box = TRUE, stack = TRUE)
FeatureStatPlot(pancreas_sub, features = c(
       "Sox9", "Anxa2", "Bicc1", # Ductal
       "Neurog3", "Hes6", # EPs
       "Fev", "Neurod1", # Pre-endocrine
       "Rbp4", "Pyy", # Endocrine
       "Ins1", "Gcg", "Sst", "Ghrl" # Beta, Alpha, Delta, Epsilon
   ), ident = "SubCellType", add_bg = TRUE, stack = TRUE,
   legend.position = "top", legend.direction = "horizontal")
FeatureStatPlot(pancreas_sub, plot_type = "box", features = c(
      "Sox9", "Anxa2", "Bicc1", # Ductal
      "Neurog3", "Hes6", # EPs
      "Fev", "Neurod1", # Pre-endocrine
      "Rbp4", "Pyy", # Endocrine
      "Ins1", "Gcg", "Sst", "Ghrl" # Beta, Alpha, Delta, Epsilon
   ), ident = "SubCellType", add_bg = TRUE, stack = TRUE, flip = TRUE,
   legend.position = "top", legend.direction = "horizontal")
# Use splitting instead of facetting
FeatureStatPlot(pancreas_sub, features = c("Neurog3", "Rbp4", "Ins1"),
   ident = "CellType", split_by = TRUE)

FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "G2M_score", reduction = "UMAP")
FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "G2M_score", reduction = "UMAP",
   bg_cutoff = -Inf)
FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "G2M_score", reduction = "UMAP",
   theme = "theme_blank")
FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "G2M_score", reduction = "UMAP",
   theme = ggplot2::theme_classic, theme_args = list(base_size = 16))

# Label and highlight cell points
FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "Rbp4", reduction = "UMAP",
   highlight = 'SubCellType == "Delta"')
FeatureStatPlot(pancreas_sub, plot_type = "dim",
   features = "Rbp4", split_by = "Phase", reduction = "UMAP",
   highlight = TRUE, theme = "theme_blank")

# Add a density layer
FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "Rbp4", reduction = "UMAP",
   add_density = TRUE)
FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "Rbp4", reduction = "UMAP",
   add_density = TRUE, density_filled = TRUE)

# Change the plot type from point to the hexagonal bin
FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "Rbp4", reduction = "UMAP",
   hex = TRUE)
FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "Rbp4", reduction = "UMAP",
   hex = TRUE, hex_bins = 20)

# Show lineages on the plot based on the pseudotime
FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "Lineage3", reduction = "UMAP",
   lineages = "Lineage3")
FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "Lineage3", reduction = "UMAP",
   lineages = "Lineage3", lineages_whiskers = TRUE)
FeatureStatPlot(pancreas_sub, plot_type = "dim", features = "Lineage3", reduction = "UMAP",
   lineages = "Lineage3", lineages_span = 0.1)

FeatureStatPlot(pancreas_sub, plot_type = "dim",
  features = c("Sox9", "Anxa2", "Bicc1"), reduction = "UMAP",
  theme = "theme_blank",
  theme_args = list(plot.subtitle = ggplot2::element_text(size = 10),
     strip.text = ggplot2::element_text(size = 8))
)

# Plot multiple features with different scales
endocrine_markers <- c("Ins1", "Gcg", "Sst", "Ghrl")
FeatureStatPlot(pancreas_sub, endocrine_markers, reduction = "UMAP", plot_type = "dim")
FeatureStatPlot(pancreas_sub, endocrine_markers, reduction = "UMAP", lower_quantile = 0,
   upper_quantile = 0.8, plot_type = "dim")
FeatureStatPlot(pancreas_sub, endocrine_markers, reduction = "UMAP",
   lower_cutoff = 1, upper_cutoff = 4, plot_type = "dim")
FeatureStatPlot(pancreas_sub, endocrine_markers, reduction = "UMAP", bg_cutoff = 2,
   lower_cutoff = 2, upper_cutoff = 4, plot_type = "dim")
FeatureStatPlot(pancreas_sub, c("Sst", "Ghrl"), split_by = "Phase", reduction = "UMAP",
   plot_type = "dim")
FeatureStatPlot(pancreas_sub, features = c("G2M_score", "nCount_RNA"),
   ident = "SubCellType", plot_type = "dim", facet_by = "Phase", split_by = TRUE, ncol = 1)

# Heatmap
features <- c(
   "Sox9", "Anxa2", "Bicc1", # Ductal
   "Neurog3", "Hes6", # EPs
   "Fev", "Neurod1", # Pre-endocrine
   "Rbp4", "Pyy", # Endocrine
   "Ins1", "Gcg", "Sst", "Ghrl" # Beta, Alpha, Delta, Epsilon
)
rows_data <- data.frame(
   Features = features,  # 'rows_name' default is "Features"
   group = c(
       "Ductal", "Ductal", "Ductal", "EPs", "EPs", "Pre-endocrine",
       "Pre-endocrine", "Endocrine", "Endocrine", "Beta", "Alpha", "Delta", "Epsilon"),
   TF = c(TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE,
       TRUE, TRUE, TRUE),
   CSPA = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE,
       FALSE, FALSE, FALSE)
)
FeatureStatPlot(pancreas_sub, features = features, ident = "SubCellType",
   plot_type = "heatmap", name = "Expression Level")
FeatureStatPlot(pancreas_sub, features = features, ident = "Phase",
   plot_type = "heatmap", name = "Expression Level", columns_split_by = "SubCellType")
FeatureStatPlot(pancreas_sub, features = features, ident = "SubCellType",
   plot_type = "heatmap", cell_type = "bars", name = "Expression Level")
FeatureStatPlot(pancreas_sub, features = features, ident = "SubCellType", cell_type = "dot",
   plot_type = "heatmap", name = "Expression Level", dot_size = function(x) sum(x > 0) / length(x),
   dot_size_name = "Percent Expressed", add_bg = TRUE, rows_data = rows_data,
   show_row_names = TRUE, rows_split_by = "group", cluster_rows = FALSE,
   column_annotation = c("Phase", "G2M_score"),
   column_annotation_type = list(Phase = "pie", G2M_score = "violin"),
   column_annotation_params = list(G2M_score = list(show_legend = FALSE)),
   row_annotation = c("TF", "CSPA"),
   row_annotation_side = "right",
   row_annotation_type = list(TF = "simple", CSPA = "simple"))
FeatureStatPlot(pancreas_sub, features = features, ident = "SubCellType", cell_type = "dot",
   plot_type = "heatmap", name = "Expression Level", dot_size = function(x) sum(x > 0) / length(x),
   dot_size_name = "Percent Expressed", add_bg = TRUE,
   rows_data = rows_data, show_column_names = TRUE, rows_split_by = "group",
   cluster_rows = FALSE, flip = TRUE, palette = "YlOrRd",
   column_annotation = c("Phase", "G2M_score"),
   column_annotation_type = list(Phase = "pie", G2M_score = "violin"),
   column_annotation_params = list(G2M_score = list(show_legend = FALSE)),
   row_annotation = c("TF", "CSPA"),
   row_annotation_side = "right",
   row_annotation_type = list(TF = "simple", CSPA = "simple"))
FeatureStatPlot(pancreas_sub, features = features, ident = "SubCellType", cell_type = "violin",
   plot_type = "heatmap", name = "Expression Level", show_row_names = TRUE,
   cluster_columns = FALSE, rows_split_by = "group", rows_data = rows_data)
FeatureStatPlot(pancreas_sub, features = features, ident = "SubCellType", cell_type = "dot",
   plot_type = "heatmap", dot_size = function(x) sum(x > 0) / length(x),
   dot_size_name = "Percent Expressed", palette = "viridis", add_reticle = TRUE,
   rows_data = rows_data, name = "Expression Level", show_row_names = TRUE,
   rows_split_by = "group")
# Visualize the markers for each sub-cell type (the markers can overlap)
# Say: markers <- Seurat::FindAllMarkers(pancreas_sub, ident = "SubCellType")
markers <- data.frame(
    avg_log2FC = c(
         3.44, 2.93, 2.72, 2.63, 2.13, 1.97, 2.96, 1.92, 5.22, 3.91, 3.64, 4.52,
         3.45, 2.45, 1.75, 2.08, 9.10, 4.45, 3.61, 6.30, 4.96, 3.49, 3.91, 3.90,
         10.58, 5.84, 4.73, 3.34, 7.22, 4.52, 10.10, 4.25),
    cluster = factor(rep(
         c("Ductal", "Ngn3 low EP", "Ngn3 high EP", "Pre-endocrine", "Beta",
           "Alpha", "Delta", "Epsilon"), each = 4),
         levels = levels(pancreas_sub$SubCellType)),
    gene = c(
         "Cyr61", "Adamts1", "Anxa2", "Bicc1", "1700011H14Rik", "Gsta3", "8430408G22Rik",
         "Anxa2", "Ppp1r14a", "Btbd17", "Neurog3", "Gadd45a", "Fev", "Runx1t1", "Hmgn3",
         "Cryba2", "Ins2", "Ppp1r1a", "Gng12", "Sytl4", "Irx1", "Tmem27", "Peg10", "Irx2",
         "Sst", "Ptprz1", "Arg1", "Frzb", "Irs4", "Mboat4", "Ghrl", "Arg1"
    )
)
FeatureStatPlot(pancreas_sub,
  features = unique(markers$gene), ident = "SubCellType", cell_type = "bars",
  plot_type = "heatmap", rows_data = markers, rows_name = "gene", rows_split_by = "cluster",
  show_row_names = TRUE, show_column_names = TRUE, name = "Expression Level",
  cluster_rows = FALSE, cluster_columns = FALSE, rows_split_palette = "Paired")

# Use plot_type = "dot" to as a shortcut for heatmap with cell_type = "dot"
FeatureStatPlot(pancreas_sub, features = features, ident = "SubCellType", plot_type = "dot")

named_features <- list(
   Ductal = c("Sox9", "Anxa2", "Bicc1"),
   EPs = c("Neurog3", "Hes6"),
   `Pre-endocrine` = c("Fev", "Neurod1"),
   Endocrine = c("Rbp4", "Pyy"),
   Beta = "Ins1", Alpha = "Gcg", Delta = "Sst", Epsilon = "Ghrl"
)
FeatureStatPlot(pancreas_sub, features = named_features, ident = "SubCellType",
   plot_type = "heatmap", name = "Expression Level", show_row_names = TRUE)

# Correlation plot
FeatureStatPlot(pancreas_sub, features = c("Pyy", "Rbp4"), plot_type = "cor",
   anno_items = c("eq", "r2", "spearman"))
FeatureStatPlot(pancreas_sub, features = c("Ins1", "Gcg", "Sst", "Ghrl"),
   plot_type = "cor")
}
}
