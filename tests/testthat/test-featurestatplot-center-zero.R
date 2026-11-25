# Toy object with a PCA reduction and scale.data for every gene, so that
# FeatureStatPlot.Seurat() sees the same values the tests derive.
make_center_zero_obj <- function() {
    skip_if_not_installed("Seurat")
    suppressPackageStartupMessages(library(Seurat))

    set.seed(42)
    mat <- matrix(
        rpois(100 * 40, lambda = 5),
        nrow = 100,
        dimnames = list(paste0("gene", 1:100), paste0("cell", 1:40))
    )
    obj <- suppressWarnings(CreateSeuratObject(mat))
    obj$Sample <- rep(c("A", "B"), each = 20)
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, nfeatures = 50, verbose = FALSE)
    # Scale all genes; the default only scales VariableFeatures and the Seurat
    # method silently drops features missing from the layer.
    obj <- ScaleData(obj, features = rownames(obj), verbose = FALSE)
    obj <- RunPCA(obj, npcs = 5, verbose = FALSE)
    obj
}

# Feature values as FeatureStatPlot.Seurat() assembles them.
pooled_values <- function(obj, features, layer = "scale.data", reduction = "pca") {
    cells <- rownames(Embeddings(obj, reduction))
    assay_data <- GetAssayData(obj, layer = layer)
    as.numeric(t(as.matrix(assay_data[features, cells, drop = FALSE])))
}

# Run `expr` with the plotthis plotting functions mocked, returning the
# captured argument list of the plotting call.
capture_plot_args <- function(expr) {
    captured <- NULL
    local_mocked_bindings(
        Heatmap = function(...) {
            captured <<- list(...)
            "heatmap"
        },
        FeatureDimPlot = function(...) {
            captured <<- list(...)
            "dim"
        },
        .package = "scplotter"
    )
    force(expr)
    captured
}

test_that("center_zero uses all values when no quantiles are supplied", {
    obj <- make_center_zero_obj()
    features <- c("gene1", "gene2", "gene3")
    pooled <- pooled_values(obj, features)
    expected <- max(abs(pooled[is.finite(pooled)]))

    for (plot_type in c("dim", "heatmap", "dot")) {
        # layer = "scale.data" also makes center_zero default to TRUE for dim
        # and avoids the "intended for 0-centered layers" warning.
        args <- capture_plot_args(FeatureStatPlot(
            obj, features,
            plot_type = plot_type, layer = "scale.data", ident = "Sample"
        ))
        expect_equal(args$lower_cutoff, -expected, info = plot_type)
        expect_equal(args$upper_cutoff, expected, info = plot_type)
        expect_null(args$lower_quantile, info = plot_type)
        expect_null(args$upper_quantile, info = plot_type)
    }
})

test_that("center_zero bounds the data by the supplied quantiles", {
    obj <- make_center_zero_obj()
    features <- c("gene1", "gene2", "gene3")
    pooled <- pooled_values(obj, features)
    probs <- c(0.05, 0.95)
    expected <- max(abs(quantile(pooled[is.finite(pooled)], probs)))

    for (plot_type in c("dim", "heatmap", "dot")) {
        args <- capture_plot_args(FeatureStatPlot(
            obj, features,
            plot_type = plot_type, layer = "scale.data", ident = "Sample",
            lower_quantile = probs[1], upper_quantile = probs[2]
        ))
        expect_equal(args$lower_cutoff, -expected, info = plot_type)
        expect_equal(args$upper_cutoff, expected, info = plot_type)
        # quantiles are still forwarded to plotthis
        expect_equal(args$lower_quantile, probs[1], info = plot_type)
        expect_equal(args$upper_quantile, probs[2], info = plot_type)
    }
})

test_that("center_zero defaults an unspecified upper quantile to 0.99", {
    obj <- make_center_zero_obj()
    features <- c("gene1", "gene2", "gene3")
    pooled <- pooled_values(obj, features)
    expected <- max(abs(quantile(pooled[is.finite(pooled)], c(0.5, 0.99))))

    args <- capture_plot_args(FeatureStatPlot(
        obj, features, plot_type = "heatmap", ident = "Sample",
        lower_quantile = 0.5
    ))
    expect_equal(args$lower_cutoff, -expected)
    expect_equal(args$upper_cutoff, expected)
})

test_that("center_zero defaults an unspecified lower quantile to 0", {
    # Positive-skew data so the upper bound dominates (the scale.data fixture
    # above is negative-skew, where dropping the lower tail changes nothing).
    data <- data.frame(Sample = rep("A", 6), f1 = c(-1, 0, 1, 2, 3, 10))
    expected <- max(abs(quantile(data$f1, c(0, 0.5))))

    args <- capture_plot_args(.feature_stat_plot(
        data,
        features = "f1", plot_type = "heatmap",
        should_shrink = FALSE, should_pivot = FALSE, downsample = NULL,
        ident = "Sample", center_zero = TRUE, upper_quantile = 0.5
    ))
    expect_equal(args$lower_cutoff, -expected)  # 1.5; all-data behavior gives 10
    expect_equal(args$upper_cutoff, expected)
})

test_that("center_zero keeps explicit cutoffs untouched", {
    obj <- make_center_zero_obj()
    args <- capture_plot_args(FeatureStatPlot(
        obj, c("gene1", "gene2"), plot_type = "heatmap", ident = "Sample",
        lower_cutoff = -3, upper_cutoff = 3
    ))
    expect_equal(args$lower_cutoff, -3)
    expect_equal(args$upper_cutoff, 3)
})

test_that("center_zero falls back when there are no finite values", {
    data <- data.frame(Sample = c("A", "A"), f1 = c(NA_real_, NA_real_))
    expect_warning(
        args <- capture_plot_args(.feature_stat_plot(
            data,
            features = "f1", plot_type = "heatmap",
            should_shrink = FALSE, should_pivot = FALSE, downsample = NULL,
            ident = "Sample", center_zero = TRUE, lower_quantile = 0.1
        )),
        "all feature values are NA or infinite"
    )
    expect_equal(args$lower_cutoff, -1)
    expect_equal(args$upper_cutoff, 1)
})
