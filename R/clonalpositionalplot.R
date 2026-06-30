#' Visualize positional properties of CDR3 sequences
#'
#' The complementarity-determining region 3 (CDR3) is the most variable region of
#' T cell and B cell receptors, and is the primary determinant of antigen
#' specificity. Analyzing how amino acid composition, diversity, and
#' physicochemical properties vary across CDR3 positions provides insight into
#' repertoire structure, selection pressures, and the biophysical constraints
#' that shape antigen recognition.
#'
#' `ClonalPositionalPlot` computes and visualizes position-level properties of
#' CDR3 sequences by aligning them at their N-terminus (Cys residue) and analyzing
#' each position independently. Three fundamentally different categories of
#' analysis are available through the `method` parameter:
#'
#' * **Amino acid frequency** (`method = "AA"`): Raw frequency of each amino acid
#'   at each position, displayed as stacked bars or a pie-chart heatmap. This
#'   reveals position-specific amino acid preferences (e.g., hydrophobic residues
#'   in the center of the CDR3 loop, or glycine enrichment at turn positions).
#'
#' * **Positional entropy** (`method = "shannon"`, `"inv.simpson"`, or
#'   `"norm.entropy"`): Quantifies the diversity of amino acid usage at each
#'   position. Low entropy at a position indicates conservation (often structural
#'   constraint), while high entropy indicates variability (potential antigen
#'   contact). Supports bar, line, heatmap, box, and violin plot types, with
#'   box/violin requiring a `group_by` variable for comparisons.
#'
#' * **Physicochemical properties** (`method = "Atchley"`, `"Kidera"`,
#'   `"stScales"`, `"tScales"`, or `"VHSE"`): Each amino acid is represented by
#'   a set of numeric factor scores that capture its physicochemical properties
#'   (e.g., hydrophobicity, size, charge, secondary structure propensity). These
#'   are averaged at each position to reveal trends along the CDR3 sequence.
#'   Supports bar and line plot types.
#'
#' The function delegates computation to three \pkg{scRepertoire} functions:
#' \code{\link[scRepertoire:percentAA]{scRepertoire::percentAA()}} for amino acid
#' frequency, \code{\link[scRepertoire:positionalEntropy]{scRepertoire::positionalEntropy()}}
#' for entropy metrics, and \code{\link[scRepertoire:positionalProperty]{scRepertoire::positionalProperty()}}
#' for physicochemical property scores.
#'
#' @section Analysis methods:
#' The `method` parameter selects both the computational method and determines
#' which plot types are available:
#'
#' **Amino acid frequency (`"AA"`):**
#' Computes the frequency of each of the 20 standard amino acids at every CDR3
#' position up to `aa_length`. Plot types:
#' * `"bar"` — Stacked bar chart with one bar per position, colored by amino acid.
#'   Faceted by `group_by`. `facet_by` is not available (set internally).
#' * `"heatmap"` — Pie-chart heatmap where each cell is a position and the pie
#'   segments show amino acid proportions at that position. Groups are on rows.
#'
#' **Entropy methods (`"shannon"`, `"inv.simpson"`, `"norm.entropy"`):**
#' Each produces a single diversity value per position per sample. Plot types:
#' * `"bar"` — Bar chart of diversity per position, faceted by `group_by`.
#' * `"line"` — Line plot connecting positions, with each group as a separate
#'   line. Useful for comparing diversity trends between conditions.
#' * `"heatmap"` — Position-by-group matrix of diversity values.
#' * `"box"` — Box plot summarizing diversity across positions for each group.
#'   Requires `group_by` to be provided (not `"Sample"`).
#' * `"violin"` — Violin plot variant of the box plot. Requires `group_by`.
#'
#' **Property methods (`"Atchley"`, `"Kidera"`, `"stScales"`, `"tScales"`, `"VHSE"`):**
#' Each yields multiple numeric property scores per position. Plot types:
#' * `"bar"` — Mean property values per position, faceted by property and
#'   `group_by`. `facet_by` is not available (set internally).
#' * `"line"` — Line plot of mean property values by position, with `group_by`
#'   as line grouping and `"property"` as faceting.
#'
#' @section Physicochemical property sets:
#' Each property method represents amino acids using a different set of derived
#' factors. The choice depends on the biological question:
#'
#' * **Atchley factors** (`"Atchley"`): Five factors derived from 494 amino acid
#'   properties — Factor I (polarity / hydrophobicity), Factor II (secondary
#'   structure), Factor III (molecular size), Factor IV (codon diversity), and
#'   Factor V (electrostatic charge). Widely used and interpretable.
#'
#' * **Kidera factors** (`"Kidera"`): Ten factors derived from 188 physical
#'   properties of amino acids, capturing more granular property variation
#'   than Atchley factors.
#'
#' * **VHSE** (`"VHSE"`): Vectors of Hydrophobic, Steric, and Electronic
#'   properties — eight descriptors derived from 50 physicochemical properties.
#'
#' * **stScales** (`"stScales"`): A set of scales derived from the structural
#'   topology of proteins (ST-scale), emphasizing local structural context.
#'
#' * **tScales** (`"tScales"`): T-scales derived from the topological properties
#'   of amino acids (T-scale), focusing on structural and folding properties.
#'
#' @section CDR3 length and alignment:
#' The `aa_length` parameter controls how many positions from the CDR3 are
#' analyzed. CDR3 sequences are aligned at their N-terminus (the conserved
#' cysteine residue), and positions beyond `aa_length` are excluded. Sequences
#' shorter than a given position contribute no data to that position. This means:
#'
#' * A larger `aa_length` includes more positions but the C-terminal positions
#'   may have sparser data (fewer sequences reach those lengths).
#' * The default of 20 amino acids captures the majority of typical TCR beta
#'   chain CDR3 lengths, though BCR heavy chain CDR3s can be substantially
#'   longer.
#' * For entropy and property methods, sparse positions may show artificially
#'   low or noisy values — consider the distribution of CDR3 lengths in your
#'   dataset when setting `aa_length`.
#'
#' @param data The product of \code{\link[scRepertoire:combineTCR]{scRepertoire::combineTCR()}},
#'  \code{\link[scRepertoire:combineBCR]{scRepertoire::combineBCR()}}, or
#'  \code{\link[scRepertoire:combineExpression]{scRepertoire::combineExpression()}}.
#'  Must contain columns with CDR3 amino acid sequences (e.g., `TRA_cdr3_aa`,
#'  `TRB_cdr3_aa`).
#' @param chain Character; the TCR or BCR chain to analyze. Default is `"TRB"`.
#'  Common values include `"TRA"` (TCR alpha), `"TRB"` (TCR beta), `"TRG"`
#'  (TCR gamma), `"TRD"` (TCR delta), `"IGH"` (BCR heavy), `"IGL"` (BCR
#'  lambda), and `"IGK"` (BCR kappa). The chain determines which CDR3 column
#'  is used for analysis.
#' @param aa_length Integer; the number of CDR3 amino acid positions to analyze,
#'  counting from the N-terminus. Default is `20`. Sequences shorter than this
#'  length contribute data only to the positions they span. Increase for
#'  receptors with longer CDR3s (e.g., BCR heavy chain).
#' @param group_by Character vector; the column(s) in `data` to group samples by.
#'  Default is `"Sample"`. Each unique combination of grouping values becomes
#'  a distinct line/bar/facet in the plot. For box and violin plot types in
#'  entropy mode, `group_by` must be explicitly provided (cannot be the default
#'  `"Sample"`).
#' @param group_by_sep Character; the separator used when concatenating multiple
#'  `group_by` columns into a single identifier. Default is `"_"`.
#' @param split_by Character vector; column(s) in `data` to split the plot by,
#'  producing separate sub-plots for each unique combination. Default is `NULL`.
#' @param method Character; the analysis method. One of:
#'  * `"AA"` — Amino acid frequency at each position (default).
#'  * `"shannon"` — Shannon entropy (diversity).
#'  * `"inv.simpson"` — Inverse Simpson index (diversity, emphasizes dominant
#'    amino acids).
#'  * `"norm.entropy"` — Normalized entropy (diversity, scaled to the range 0–1).
#'  * `"Atchley"` — Atchley factor scores (5 properties).
#'  * `"Kidera"` — Kidera factor scores (10 properties).
#'  * `"stScales"` — ST-scale property scores.
#'  * `"tScales"` — T-scale property scores.
#'  * `"VHSE"` — VHSE property scores (8 descriptors).
#'
#'  See the **Analysis methods** and **Physicochemical property sets**
#'  sections for details on each method.
#' @param order A named list specifying the order of factor levels for grouping
#'  variables. For example, `list(Type = c("L", "B"))`. Default is `NULL`,
#'  which uses the order present in the data. Names in the list must match
#'  column names in `data`.
#' @param plot_type Character; the type of visualization. One of `"bar"`
#'  (default), `"line"`, `"heatmap"`, `"box"`, or `"violin"`. Not all plot
#'  types are available for all methods — see the **Analysis methods**
#'  section for the full mapping.
#' @param theme_args A named list of arguments passed to
#'  \code{\link[ggplot2:theme]{ggplot2::theme()}} for customizing the plot
#'  appearance. For bar plots, `panel.grid.major.y` defaults to
#'  `element_blank()` for a cleaner look.
#' @param xlab Character; the x-axis label. Default is `NULL`, which
#'  automatically uses `"Position"`.
#' @param ylab Character; the y-axis label. Default is `NULL`, which
#'  automatically uses a method-appropriate label (`"Amino Acid Frequency"`
#'  for AA, the method name for entropy, or `"Mean Values"` for properties).
#' @param facet_by A character vector of column names to facet the plots by.
#'  Default is `NULL`. Note that for AA and property bar plots, `facet_by`
#'  is set internally and should not be specified manually (doing so will
#'  raise an error). For entropy line plots, it can be used to add an
#'  additional faceting dimension.
#' @param facet_ncol Integer; the number of columns in the facet grid.
#'  Default is `NULL`, which uses `1` for bar plots and automatic
#'  determination for other plot types.
#' @param facet_nrow Integer; the number of rows in the facet grid.
#'  Default is `NULL`. For property bar plots, defaults to the number of
#'  properties in the selected factor set.
#' @param aspect.ratio Numeric; the aspect ratio (height / width) of plot
#'  panels. Default is `NULL`, which uses method- and plot_type-appropriate
#'  defaults (e.g., `2 / aa_length` for entropy bar plots,
#'  `6 / aa_length` for entropy line plots, `10 / aa_length` for box/violin
#'  plots).
#' @param ... Additional arguments passed to the underlying \pkg{plotthis}
#'  visualization function, depending on `plot_type`:
#'  * For `"bar"`: \code{\link[plotthis:BarPlot]{plotthis::BarPlot()}}
#'  * For `"line"`: \code{\link[plotthis:LinePlot]{plotthis::LinePlot()}}
#'  * For `"heatmap"`: \code{\link[plotthis:Heatmap]{plotthis::Heatmap()}}
#'  * For `"box"`: \code{\link[plotthis:BoxPlot]{plotthis::BoxPlot()}}
#'  * For `"violin"`: \code{\link[plotthis:ViolinPlot]{plotthis::ViolinPlot()}}
#'
#'  Common arguments include `title`, `legend.position`, and color palette
#'  parameters. See the respective \pkg{plotthis} documentation for
#'  available options.
#' @return A `ggplot` object (or a list of `ggplot` objects if `combine = FALSE`
#'  is passed via `...`)
#' @note
#' * Positional analysis requires CDR3 amino acid sequences. These are
#'   automatically extracted by \pkg{scRepertoire} during data processing.
#'   Ensure your input data has been processed with
#'   \code{\link[scRepertoire:combineTCR]{scRepertoire::combineTCR()}} or
#'   \code{\link[scRepertoire:combineBCR]{scRepertoire::combineBCR()}}.
#' * For the `"AA"` method, `facet_by` is internally set to `group_by` and
#'   must not be provided by the user. The same applies to property bar plots.
#' * Box and violin plot types (entropy methods only) require `group_by` to
#'   be explicitly set to a value other than the default `"Sample"`. This is
#'   because these plots need meaningful groups to compare distributions.
#' * Property methods (`"Atchley"`, `"Kidera"`, etc.) only support `"bar"`,
#'   `"line"`, and in some cases `"heatmap"` plot types. Box and violin plots
#'   are not available for property methods.
#' * CDR3 sequences shorter than `aa_length` will have missing data at
#'   C-terminal positions. This can introduce noise at those positions,
#'   especially for entropy and property analyses. Consider the length
#'   distribution of your CDR3s when choosing `aa_length`.
#' * The `"AA"` heatmap uses a pie-chart cell type (`cell_type = "pie"`),
#'   where each cell represents a position and pie segments show amino acid
#'   composition. This is a distinctive visualization unique to this method.
#' @seealso
#' * \code{\link{ClonalKmerPlot}} for k-mer (short motif) analysis of CDR3
#'   sequences
#' * \code{\link{ClonalLengthPlot}} for CDR3 length distribution analysis
#' * \code{\link{ClonalDiversityPlot}} for repertoire-level diversity metrics
#' * \code{\link{ClonalGeneUsagePlot}} for V(D)J gene segment usage analysis
#' * \code{\link[scRepertoire:percentAA]{scRepertoire::percentAA()}},
#'   \code{\link[scRepertoire:positionalEntropy]{scRepertoire::positionalEntropy()}},
#'   \code{\link[scRepertoire:positionalProperty]{scRepertoire::positionalProperty()}}
#'   for the underlying computational methods
#' @export
#' @importFrom ggplot2 element_blank
#' @importFrom scRepertoire percentAA positionalEntropy positionalProperty
#' @importFrom plotthis BarPlot LinePlot Heatmap BoxPlot ViolinPlot
#' @examples
#' \donttest{
#' set.seed(8525)
#' data(contig_list, package = "scRepertoire")
#' data <- scRepertoire::combineTCR(contig_list,
#'    samples = c("P17B", "P17L", "P18B", "P18L", "P19B","P19L", "P20B", "P20L"))
#' data <- scRepertoire::addVariable(data,
#'   variable.name = "Type",
#'   variables = factor(rep(c("B", "L"), 4), levels = c("L", "B"))
#' )
#'
#' ClonalPositionalPlot(data)
#' ClonalPositionalPlot(data, method = "shannon")
#' ClonalPositionalPlot(data, method = "norm.entropy", plot_type = "heatmap")
#' ClonalPositionalPlot(data, method = "Atchley", group_by = "Type", plot_type = "bar")
#' ClonalPositionalPlot(data, method = "Atchley", plot_type = "line")
#' }
ClonalPositionalPlot <- function (
    data, chain = "TRB", aa_length = 20, group_by = "Sample", group_by_sep = "_", split_by = NULL,
    method = c("AA", "shannon", "inv.simpson", "norm.entropy", "Atchley",
        "Kidera", "stScales", "tScales", "VHSE"), order = NULL,
    plot_type = c("bar", "line", "heatmap", "box", "violin"), theme_args = list(),
    xlab = NULL, ylab = NULL, facet_by = NULL, facet_ncol = NULL, facet_nrow = NULL,
    aspect.ratio = NULL,
    ...
) {
    method <- match.arg(method)
    plot_type <- match.arg(plot_type)
    if (plot_type %in% c("box", "violin")) {
        if (is.null(group_by) || identical(group_by, "Sample")) {
            stop("'group_by' must be provided for box/violin ClonalPositionalPlot")
        }
        all_groupings <- unique(c("Sample", group_by, facet_by, split_by))
    } else {
        all_groupings <- unique(c(group_by, facet_by, split_by))
    }
    grouping_levels <- get_clonal_grouping_levels(data, all_groupings, order)
    data <- merge_clonal_groupings(data, all_groupings)

    if (method == "AA") {
        if (!is.null(facet_by)) {
            stop("'facet_by' should not be specified for AA bar plot in ClonalPositionalPlot.")
        }
        data <- percentAA(data, chain = chain, aa.length = aa_length, group.by = ".group",
            exportTable = TRUE)
        data <- separate(data, "group", into = all_groupings, sep = " // ")
        for (gl in names(grouping_levels)) {
            if (!is.null(data[[gl]])) {
                data[[gl]] <- factor(data[[gl]], levels = grouping_levels[[gl]])
            }
        }

        if (plot_type == "bar") {
            theme_args$panel.grid.major.y <- theme_args$panel.grid.major.y %||% element_blank()

            BarPlot(data, x = "variable", y = "value", group_by = "AA", position = "stack",
                xlab = xlab %||% "Position", ylab = ylab %||% "Amino Acid Frequency",
                split_by = split_by, facet_by = group_by, facet_ncol = facet_ncol %||% 1, facet_nrow = facet_nrow,
                x_text_angle = 90, facet_args = list(strip.position = "right"),
                aspect.ratio = aspect.ratio %||% (2 / aa_length), theme_args = theme_args, ...
            )
        } else if (plot_type == "heatmap") {
            data <- data %>% unite(".group", !!!syms(group_by), sep = group_by_sep)
            allgroups <- unique(data$.group)
            data <- data %>%
                pivot_wider(names_from = ".group", values_from = "value") %>%
                rename(Position = "variable")

            Heatmap(data, columns_by = "Position", rows_by = allgroups, rows_name = paste(group_by, collapse = group_by_sep),
                cell_type = "pie", pie_group_by = "AA", cluster_rows = FALSE, cluster_columns = FALSE,
                pie_values = "sum", ...)
        } else {
            stop("Only 'bar' and 'heatmap' plot types are supported for AA in ClonalPositionalPlot.")
        }
    } else if (method %in% c("shannon", "inv.simpson", "norm.entropy")) {
        data <- positionalEntropy(data, chain = chain, aa.length = aa_length, group.by = ".group",
            method = method, exportTable = TRUE) %>%
            separate("Var1", into = all_groupings, sep = " // ") %>%
            rename(Position = "Var2") %>%
            unite(".group", !!!syms(group_by), sep = group_by_sep)

        for (gl in names(grouping_levels)) {
            if (!is.null(data[[gl]])) {
                data[[gl]] <- factor(data[[gl]], levels = grouping_levels[[gl]])
            }
        }

        group_by <- paste(group_by, sep = group_by_sep)
        data <- rename(data, !!sym(group_by) := ".group")

        if (plot_type == "bar") {
            if (!is.null(facet_by)) {
                stop("'facet_by' should not be specified for entropy bar plot in ClonalPositionalPlot.")
            }
            theme_args$panel.grid.major.y <- theme_args$panel.grid.major.y %||% element_blank()

            BarPlot(data, x = "Position", y = "value",
                xlab = xlab %||% "Position", ylab = ylab %||% method, split_by = split_by,
                facet_by = group_by, facet_ncol = facet_ncol %||% 1, facet_nrow = facet_nrow,
                x_text_angle = 90, legend.position = "none", facet_args = list(strip.position = "right"),
                aspect.ratio = aspect.ratio %||% (2 / aa_length), theme_args = theme_args, ...
            )
        } else if (plot_type == "line") {
            LinePlot(data, x = "Position", y = "value", group_by = group_by, pt_size = 2,
                xlab = xlab %||% "Position", ylab = ylab %||% method, split_by = split_by,
                facet_by = facet_by, facet_ncol = facet_ncol, facet_nrow = facet_nrow, x_text_angle = 90,
                facet_args = list(strip.position = "right"), aspect.ratio = aspect.ratio %||% (6 / aa_length),
                theme_args = theme_args, ...
            )
        } else if (plot_type == "heatmap") {
            allgroups <- unique(data[[group_by]])
            data <- data %>% pivot_wider(names_from = group_by, values_from = "value")

            Heatmap(data, columns_by = "Position", rows_by = allgroups, rows_name = group_by,
                name = method, cluster_columns = FALSE, show_column_names = TRUE, show_row_names = TRUE,
                ...)
        } else if (plot_type == "box") {
            BoxPlot(data, x = "Position", y = "value", xlab = xlab %||% "Position",
                ylab = ylab %||% method, split_by = split_by, group_by = group_by, facet_ncol = facet_ncol,
                facet_nrow = facet_nrow, x_text_angle = 90, theme_args = theme_args,
                aspect.ratio = aspect.ratio %||% (10 / aa_length), ...
            )
        } else if (plot_type == "violin") {
            ViolinPlot(data, x = "Position", y = "value", xlab = xlab %||% "Position",
                ylab = ylab %||% method, split_by = split_by, group_by = group_by, facet_ncol = facet_ncol,
                facet_nrow = facet_nrow, x_text_angle = 90, theme_args = theme_args,
                aspect.ratio = aspect.ratio %||% (10 / aa_length), ...
            )
        }
    } else {
        # https://github.com/ncborcherding/scRepertoire/issues/420
        data <- positionalProperty(data, chain = chain, aa.length = aa_length, group.by = ".group",
            method = method)$data %>%
            separate("group", into = all_groupings, sep = " // ") %>%
            rename(Position = "position") %>%
            unite(".group", !!!syms(group_by), sep = group_by_sep)

        for (gl in names(grouping_levels)) {
            if (!is.null(data[[gl]])) {
                data[[gl]] <- factor(data[[gl]], levels = grouping_levels[[gl]])
            }
        }

        group_by <- paste(group_by, sep = group_by_sep)
        data <- rename(data, !!sym(group_by) := ".group")

        n_properties <- length(unique(data$property))

        if (plot_type == "bar") {
            if (!is.null(facet_by)) {
                stop("'facet_by' should not be specified for property bar plot in ClonalPositionalPlot.")
            }
            theme_args$panel.grid.major.y <- theme_args$panel.grid.major.y %||% element_blank()

            BarPlot(data, x = "Position", y = "mean",
                xlab = xlab %||% "Position", ylab = ylab %||% "Mean Values", split_by = split_by,
                facet_by = c("property", group_by), facet_ncol = facet_ncol, facet_nrow = facet_nrow %||% n_properties,
                x_text_angle = 90, legend.position = "none",
                aspect.ratio = aspect.ratio %||% (4 / aa_length), theme_args = theme_args, ...
            )
        } else if (plot_type == "line") {
            theme_args$panel.grid.major.y <- theme_args$panel.grid.major.y %||% element_blank()
            LinePlot(data, x = "Position", y = "mean", group_by = group_by, pt_size = 2,
                xlab = xlab %||% "Position", ylab = ylab %||% "Mean Values", split_by = split_by,
                facet_by = "property", facet_ncol = facet_ncol %||% 1, facet_nrow = facet_nrow, x_text_angle = 90,
                facet_args = list(strip.position = "right"), aspect.ratio = aspect.ratio %||% (6 / aa_length),
                theme_args = theme_args, ...
            )
        } else {
            stop("Only 'bar' and 'line' plot types are supported for property in ClonalPositionalPlot.")
        }
    }
}

#' Visualize CDR3 k-mer (motif) frequency
#'
#' Short amino acid motifs within CDR3 sequences — termed *k-mers* — can reveal
#' shared binding specificities, common structural elements, and repertoire-level
#' sequence features that are not apparent from full-length sequence analysis
#' alone. Specific k-mers may be enriched in responses to particular antigens,
#' represent public TCR/BCR motifs shared across individuals, or reflect
#' convergent recombination events.
#'
#' `ClonalKmerPlot` extracts k-mers of length `k` from CDR3 amino acid (or
#' nucleotide) sequences and visualizes their frequency across samples and
#' conditions. The analysis uses \code{\link[scRepertoire:percentKmer]{scRepertoire::percentKmer()}}
#' to identify the `top` most frequent k-mers, then displays them using bar,
#' line, or heatmap visualizations.
#'
#' K-mer analysis is complementary to positional analysis
#' (\code{\link{ClonalPositionalPlot}}): while positional analysis asks "what
#' happens at position X?", k-mer analysis asks "what short sequence motifs
#' appear, regardless of their exact position?" This position-independent
#' perspective can capture sequence features that are distributed across
#' different CDR3 locations.
#'
#' @section K-mer length selection:
#' The choice of `k` represents a fundamental trade-off:
#'
#' * **k = 2** (dipeptides): Highly recurrent but often non-specific — many
#'   dipeptides appear in functionally unrelated receptors. Best for broad
#'   surveys of amino acid pairing preferences.
#' * **k = 3** (tripeptides, default): The most commonly used length. Tripeptides
#'   are long enough to capture meaningful motifs (e.g., the "CASS" motif at the
#'   start of TRB CDR3s) while being short enough to recur frequently across
#'   samples for robust statistical analysis.
#' * **k = 4 or 5**: More specific motifs that are more likely to reflect
#'   genuine functional constraints or antigen-driven selection. However,
#'   longer k-mers appear less frequently, requiring larger datasets for
#'   reliable frequency estimation.
#'
#' The `top` parameter controls how many of the most frequent k-mers are
#' displayed. The default of 25 provides a manageable view for bar and line
#' plots; increase for heatmap visualizations that can accommodate many more
#' motifs.
#'
#' @section Plot types:
#' Three visualization types are available:
#'
#' * **`"bar"`** (default): Bar chart of k-mer frequencies, faceted by
#'   `group_by`. Best for comparing motif usage across a moderate number of
#'   samples or conditions. `facet_by` is not available (set internally).
#'
#' * **`"line"`**: Line plot connecting k-mer frequencies, with each group
#'   as a separate line. Useful for visualizing trends across motifs or
#'   comparing the frequency profile shape between conditions. Supports
#'   `facet_by` for additional faceting dimensions.
#'
#' * **`"heatmap"`**: K-mer by group matrix showing frequency as color
#'   intensity. Ideal for surveying many motifs across many samples
#'   simultaneously, revealing clusters of co-enriched motifs.
#'
#' @param data The product of \code{\link[scRepertoire:combineTCR]{scRepertoire::combineTCR()}},
#'  \code{\link[scRepertoire:combineBCR]{scRepertoire::combineBCR()}}, or
#'  \code{\link[scRepertoire:combineExpression]{scRepertoire::combineExpression()}}.
#'  Must contain columns with CDR3 sequences (amino acid or nucleotide).
#' @param chain Character; the TCR or BCR chain to analyze. Default is `"TRB"`.
#'  Common values include `"TRA"`, `"TRB"`, `"TRG"`, `"TRD"`, `"IGH"`,
#'  `"IGL"`, and `"IGK"`.
#' @param clone_call Character; the column name specifying which CDR3 sequence
#'  to use for k-mer extraction. Default is `"aa"` (amino acid). Use `"nt"`
#'  for nucleotide sequences or `"strict"` for strict clonotype calls.
#' @param k Integer; the length of k-mers (motifs) to extract. Default is `3`.
#'  See the **K-mer length selection** section for guidance on choosing an
#'  appropriate value.
#' @param top Integer; the number of most frequent k-mers to display. Default
#'  is `25`. K-mers are ranked by total frequency across all groups. Increase
#'  for heatmap visualizations; decrease for focused bar charts.
#' @param group_by Character vector; the column(s) in `data` to group samples by.
#'  Default is `"Sample"`.
#' @param group_by_sep Character; the separator used when concatenating multiple
#'  `group_by` columns into a single identifier. Default is `"_"`.
#' @param order A named list specifying the order of factor levels for grouping
#'  variables. For example, `list(Type = c("L", "B"))`. Default is `NULL`,
#'  which uses the order present in the data.
#' @param facet_by A character vector of column names to facet the plots by.
#'  Default is `NULL`. For bar plots, `facet_by` is set internally to
#'  `group_by` and must not be specified manually (doing so will raise an
#'  error). For line plots, it can be used for additional faceting.
#' @param split_by Character vector; column(s) in `data` to split the plot by,
#'  producing separate sub-plots for each unique combination. Default is `NULL`.
#' @param plot_type Character; the type of visualization. One of `"bar"`
#'  (default), `"line"`, or `"heatmap"`. See the **Plot types** section
#'  for guidance.
#' @param theme_args A named list of arguments passed to
#'  \code{\link[ggplot2:theme]{ggplot2::theme()}} for customizing the plot
#'  appearance. For bar and line plots, `panel.grid.major.y` defaults to
#'  `element_blank()`.
#' @param aspect.ratio Numeric; the aspect ratio (height / width) of plot
#'  panels. Default is `NULL`, which uses `4 / length(motifs)` for bar
#'  plots and `8 / length(motifs)` for line plots, automatically scaling
#'  with the number of k-mers.
#' @param facet_ncol Integer; the number of columns in the facet grid.
#'  Default is `NULL`, which uses `1` for bar plots.
#' @param ... Additional arguments passed to the underlying \pkg{plotthis}
#'  visualization function, depending on `plot_type`:
#'  * For `"bar"`: \code{\link[plotthis:BarPlot]{plotthis::BarPlot()}}
#'  * For `"line"`: \code{\link[plotthis:LinePlot]{plotthis::LinePlot()}}
#'  * For `"heatmap"`: \code{\link[plotthis:Heatmap]{plotthis::Heatmap()}}
#'
#'  Common arguments include `title`, `legend.position`, `show_row_names`,
#'  `show_column_names`, and color palette parameters. See the respective
#'  \pkg{plotthis} documentation for available options.
#' @return A `ggplot` object (or a list of `ggplot` objects if `combine = FALSE`
#'  is passed via `...`)
#' @note
#' * K-mer extraction is performed by
#'   \code{\link[scRepertoire:percentKmer]{scRepertoire::percentKmer()}}, which
#'   uses a sliding window across CDR3 sequences. The resulting frequencies
#'   represent the percentage of sequences containing each k-mer at least once.
#' * For bar plots, `facet_by` is internally set to `group_by` and must not
#'   be provided by the user.
#' * The number of possible k-mers grows exponentially with `k` (20^k for
#'   amino acids), but only the `top` most frequent are displayed. Ensure
#'   `top` is set high enough to capture motifs of interest.
#' * K-mer frequencies can be noisy when sample sizes are small. Consider
#'   using `split_by` or `facet_by` to disaggregate data rather than relying
#'   on small within-group sample sizes.
#' * For nucleotide k-mers (`clone_call = "nt"`), the alphabet size is 4
#'   rather than 20, so shorter k values (2-3) are generally appropriate.
#' @seealso
#' * \code{\link{ClonalPositionalPlot}} for position-specific CDR3 analysis
#'   (amino acid frequency, entropy, and physicochemical properties)
#' * \code{\link{ClonalLengthPlot}} for CDR3 length distribution analysis
#' * \code{\link{ClonalGeneUsagePlot}} for V(D)J gene segment usage
#' * \code{\link{ClonalDiversityPlot}} for repertoire-level diversity metrics
#' * \code{\link[scRepertoire:percentKmer]{scRepertoire::percentKmer()}} for
#'   the underlying k-mer computation
#' @export
#' @importFrom tidyr pivot_longer separate unite
#' @importFrom dplyr %>% rename
#' @importFrom scRepertoire percentKmer
#' @importFrom plotthis BarPlot Heatmap
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
#' ClonalKmerPlot(data)
#' ClonalKmerPlot(data, group_by = "Type")
#' ClonalKmerPlot(data, group_by = "Type", plot_type = "line")
#' ClonalKmerPlot(data, group_by = "Type", plot_type = "heatmap")
#' }
ClonalKmerPlot <- function (
    data, chain = "TRB", clone_call = "aa", k = 3, top = 25, group_by = "Sample",
    group_by_sep = "_", facet_by = NULL, split_by = NULL, order = NULL,
    plot_type = c("bar", "line", "heatmap"), theme_args = list(), aspect.ratio = NULL,
    facet_ncol = NULL, ...
) {
    plot_type <- match.arg(plot_type)
    all_groupings <- unique(c(group_by, split_by))
    grouping_levels <- get_clonal_grouping_levels(data, all_groupings, order)
    data <- merge_clonal_groupings(data, all_groupings)
    data <- percentKmer(data, chain = chain, cloneCall = clone_call, motif.length = k,
        top.motifs = top, group.by = ".group", exportTable = TRUE)
    data <- as.data.frame(data)
    motifs <- colnames(data)
    data$.group <- rownames(data)
    data <- data %>%
        separate(".group", into = all_groupings, sep = " // ") %>%
        unite(".group", !!!syms(group_by), sep = group_by_sep) %>%
        rename(!!sym(paste(group_by, sep = group_by_sep)) := ".group")

    for (gl in names(grouping_levels)) {
        if (!is.null(data[[gl]])) {
            data[[gl]] <- factor(data[[gl]], levels = grouping_levels[[gl]])
        }
    }

    group_by <- paste(group_by, sep = group_by_sep)

    if (plot_type == "bar") {
        if (!is.null(facet_by)) {
            stop("'facet_by' should not be specified in bar ClonalKmerPlot.")
        }

        data <- data %>% pivot_longer(cols = motifs, names_to = "Motifs", values_to = "Frequency")
        theme_args$panel.grid.major.y <- theme_args$panel.grid.major.y %||% element_blank()

        BarPlot(data, x = "Motifs", y = "Frequency", facet_by = group_by,
            xlab = "Motifs", ylab = "Frequency", split_by = split_by,
            facet_ncol = facet_ncol %||% 1, x_text_angle = 90, facet_args = list(strip.position = "right"),
            aspect.ratio = aspect.ratio %||% (4 / length(motifs)), legend.position = "none",
            theme_args = theme_args, ...
        )
    } else if (plot_type == "line") {
        data <- data %>% pivot_longer(cols = motifs, names_to = "Motifs", values_to = "Frequency")
        theme_args$panel.grid.major.y <- theme_args$panel.grid.major.y %||% element_blank()

        LinePlot(data, x = "Motifs", y = "Frequency", group_by = group_by, pt_size = 2,
            xlab = "Motifs", ylab = "Frequency", split_by = split_by, facet_by = facet_by,
            facet_ncol = facet_ncol, x_text_angle = 90,
            aspect.ratio = aspect.ratio %||% (8 / length(motifs)), theme_args = theme_args, ...
        )
    } else if (plot_type == "heatmap") {
        args <- rlang::dots_list(...)
        args$data <- data
        args$columns_by <- group_by
        args$rows_by <- motifs
        args$rows_name <- "Motifs"
        args$name <- "Frequency"
        args$show_row_names <- args$show_row_names %||% TRUE
        args$show_column_names <- args$show_column_names %||% TRUE
        do_call(Heatmap, args)
    }
}
