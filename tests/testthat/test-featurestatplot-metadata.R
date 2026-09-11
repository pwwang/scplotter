test_that("FeatureStatPlot() handles metadata features on objects without a reduction", {
    skip_if_not_installed("Seurat")
    suppressPackageStartupMessages(library(Seurat))

    set.seed(42)
    mat <- matrix(
        rpois(100 * 40, lambda = 5),
        nrow = 100,
        dimnames = list(paste0("gene", 1:100), paste0("cell", 1:40))
    )
    obj <- CreateSeuratObject(mat)
    obj$Sample <- rep(c("A", "B"), each = 20)

    # No reduction on purpose. Features that live in meta.data have no assay
    # data, and a NULL `reduction` selects the branch that cbind()s the assay
    # data onto the metadata. cbind() does not drop a NULL operand, so this used
    # to fail with "arguments imply differing number of rows: 40, 0".
    expect_no_error(
        p <- FeatureStatPlot(
            obj,
            features = "nCount_RNA",
            ident = "Sample",
            group_by = "Sample",
            plot_type = "violin"
        )
    )
    expect_false(is.null(p))

    # gene features still need an assay layer: without scale.data the message
    # should name the missing layer rather than fail inside cbind()
    expect_error(
        FeatureStatPlot(
            obj,
            features = "gene1",
            ident = "Sample",
            group_by = "Sample",
            plot_type = "violin"
        ),
        "does not have any data in layer 'scale.data'"
    )
})

test_that("FeatureStatPlot() still plots metadata features when a reduction exists", {
    skip_if_not_installed("Seurat")
    suppressPackageStartupMessages(library(Seurat))

    set.seed(42)
    mat <- matrix(
        rpois(100 * 40, lambda = 5),
        nrow = 100,
        dimnames = list(paste0("gene", 1:100), paste0("cell", 1:40))
    )
    obj <- CreateSeuratObject(mat)
    obj$Sample <- rep(c("A", "B"), each = 20)
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, nfeatures = 50, verbose = FALSE)
    obj <- ScaleData(obj, verbose = FALSE)
    obj <- RunPCA(obj, npcs = 5, verbose = FALSE)

    expect_no_error(
        p <- FeatureStatPlot(
            obj,
            features = "nCount_RNA",
            ident = "Sample",
            group_by = "Sample",
            plot_type = "violin"
        )
    )
    expect_false(is.null(p))
})
