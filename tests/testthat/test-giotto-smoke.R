# First automated coverage of the Giotto code paths.
# Uses the mini giotto objects bundled with GiottoData (offline, no downloads).

need <- function(pkg, min_version) {
    skip_if_not_installed(pkg)
    if (utils::packageVersion(pkg) < min_version) {
        skip(paste0(pkg, " >= ", min_version, " required, have ", utils::packageVersion(pkg)))
    }
}

first_categorical <- function(g) {
    md <- GiottoClass::getCellMetadata(g, output = "data.table")
    cand <- names(md)[vapply(md, function(x) is.character(x) || is.factor(x), logical(1))]
    setdiff(cand, c("cell_ID", "fov"))[1]
}

built_rows <- function(p) sum(vapply(ggplot2::ggplot_build(p)$data, nrow, integer(1)))

test_that("SpatDimPlot() renders a visium giotto object with and without its image", {
    need("GiottoClass", "0.5.1")
    need("GiottoData", "0.3.8")
    g <- GiottoData::loadGiottoMini("visium")
    grp <- first_categorical(g)
    skip_if(is.na(grp), "no categorical metadata column in the visium mini")
    p <- SpatDimPlot(g, group_by = grp)
    expect_s3_class(p, "ggplot")
    expect_gt(built_rows(p), 0)
    p_img <- SpatDimPlot(g, group_by = grp, image = TRUE)
    expect_s3_class(p_img, "ggplot")
    expect_gt(built_rows(p_img), built_rows(p))
})

test_that("SpatFeaturePlot() maps expression onto spatial coordinates", {
    need("GiottoClass", "0.5.1")
    need("GiottoData", "0.3.8")
    g <- GiottoData::loadGiottoMini("visium")
    feats <- head(GiottoClass::featIDs(g, feat_type = "rna"), 2)
    skip_if(length(feats) < 2, "mini object has fewer than 2 features")
    p <- SpatFeaturePlot(g, features = feats)
    expect_s3_class(p, "ggplot")
    expect_gt(built_rows(p), 0)
})

test_that("CellDimPlot()/CellStatPlot()/FeatureStatPlot() work on a giotto object", {
    need("GiottoClass", "0.5.1")
    need("GiottoData", "0.3.8")
    g <- GiottoData::loadGiottoMini("visium")
    grp <- first_categorical(g)
    skip_if(is.na(grp), "no categorical metadata column in the visium mini")
    expect_s3_class(CellDimPlot(g, group_by = grp), "ggplot")
    expect_s3_class(CellStatPlot(g, ident = grp, plot_type = "bar"), "ggplot")
    feats <- head(GiottoClass::featIDs(g, feat_type = "rna"), 2)
    skip_if(length(feats) < 1, "mini object has no features")
    expect_s3_class(
        FeatureStatPlot(g, features = feats, plot_type = "violin", ident = grp),
        "ggplot"
    )
})

test_that("molecule and shapes layers run on an FOV-style mini", {
    need("GiottoClass", "0.5.1")
    need("GiottoData", "0.3.8")
    g <- GiottoData::loadGiottoMini("cosmx")
    mols <- head(GiottoClass::featIDs(g, feat_type = "rna"), 2)
    skip_if(length(mols) < 1, "cosmx mini has no features")
    p <- SpatDimPlot(g, group_by = "molecules", features = mols)
    expect_s3_class(p, "ggplot")
    expect_gt(built_rows(p), 0)
    grp <- first_categorical(g)
    skip_if(is.na(grp), "no categorical metadata column in the cosmx mini")
    expect_s3_class(SpatDimPlot(g, group_by = grp, shapes = TRUE), "ggplot")
})

test_that("the spatial network overlay keeps the networkDT column contract", {
    need("GiottoClass", "0.5.1")
    need("GiottoData", "0.3.8")
    g <- GiottoData::loadGiottoMini("visium")
    net <- GiottoClass::getSpatialNetwork(g, output = "networkDT", verbose = FALSE)
    expect_true(all(
        c("sdimx_begin", "sdimy_begin", "sdimx_end", "sdimy_end", "weight") %in% names(net)
    ))
    grp <- first_categorical(g)
    skip_if(is.na(grp), "no categorical metadata column in the visium mini")
    expect_s3_class(SpatDimPlot(g, group_by = grp, graph = TRUE), "ggplot")
})
