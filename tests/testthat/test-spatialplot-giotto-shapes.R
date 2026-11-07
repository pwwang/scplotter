# Regression: shapes_fill_by must accept a cell metadata column for the polygon
# fill. featIDs(g, feat_type = "cell") is empty on CosMx objects, so the old
# gate (`shapes_fill_by %in% featIDs(...)`) never merged metadata columns and
# every polygon was filled black while the legend listed the levels.

need <- function(pkg, min_version) {
    skip_if_not_installed(pkg)
    if (utils::packageVersion(pkg) < min_version) {
        skip(paste0(pkg, " >= ", min_version, " required, have ", utils::packageVersion(pkg)))
    }
}

test_that("shapes_fill_by fills polygons from a cell metadata column", {
    need("GiottoClass", "0.5.1")
    need("GiottoData", "0.3.8")
    g <- GiottoData::loadGiottoMini("cosmx")
    poly_rows <- tryCatch(
        {
            poly <- GiottoClass::getPolygonInfo(gobject = g)
            if (is.list(poly)) poly <- poly[[1]]
            nrow(poly)
        },
        error = function(e) NA_integer_
    )
    skip_if(is.na(poly_rows) || poly_rows < 1, "cosmx mini has no usable polygon geometry")
    md <- GiottoClass::pDataDT(g, feat_type = "rna")
    cand <- names(md)[vapply(md, function(x) is.character(x) || is.factor(x), logical(1))]
    cand <- setdiff(cand, c("cell_ID", "fov"))
    cand <- cand[vapply(cand, function(nm) !anyNA(md[[nm]]), logical(1))]
    skip_if(!length(cand), "no categorical metadata column without NAs in the cosmx mini")
    grp <- cand[1]
    levels <- unique(as.character(md[[grp]]))
    skip_if(length(levels) < 2, "metadata column has only one level")
    skip_if(grp %in% GiottoClass::featIDs(g, feat_type = "rna"),
        "metadata column is also an rna feature ID; metadata path not exercised")

    p <- SpatDimPlot(g, shapes = TRUE, points = FALSE, shapes_fill_by = grp)
    expect_s3_class(p, "ggplot")

    # the fill actually varies: one fill per level, not a single uniform fill
    built <- ggplot2::ggplot_build(p)
    fills <- unique(stats::na.omit(unlist(lapply(built$data, function(d) d$fill))))
    expect_length(fills, length(levels))
    # and the fill scale has one entry per level
    scale <- built$plot$scales$get_scales("fill")
    expect_false(is.null(scale))
    expect_setequal(as.character(scale$get_limits()), levels)
})

# Regression: the same gate also blocked FEATURE fills. It tested membership against
# featIDs(g, feat_type = shapes_feat_type), and shapes_feat_type is "cell", where featIDs() is
# length 0 -- so a gene name never reached spatValues() and the polygons were filled black with
# no fill scale at all. This is observable on the mini, unlike the metadata path.
test_that("shapes_fill_by fills polygons from an rna feature", {
    need("GiottoClass", "0.5.1")
    need("GiottoData", "0.3.8")
    g <- GiottoData::loadGiottoMini("cosmx")
    poly_rows <- tryCatch(
        {
            poly <- GiottoClass::getPolygonInfo(gobject = g)
            if (is.list(poly)) poly <- poly[[1]]
            nrow(poly)
        },
        error = function(e) NA_integer_
    )
    skip_if(is.na(poly_rows) || poly_rows < 1, "cosmx mini has no usable polygon geometry")
    expr_ok <- tryCatch({
        GiottoClass::getExpression(g, values = "normalized", feat_type = "rna")
        TRUE
    }, error = function(e) FALSE)
    skip_if(!expr_ok, "cosmx mini has no normalized expression values")
    feat <- GiottoClass::featIDs(g, feat_type = "rna")[1]
    skip_if(is.na(feat), "no rna features")

    p <- SpatDimPlot(g, shapes = TRUE, points = FALSE, shapes_fill_by = feat)
    expect_s3_class(p, "ggplot")
    built <- ggplot2::ggplot_build(p)
    # a continuous fill scale must exist; before the fix there was none at all
    expect_false(is.null(built$plot$scales$get_scales("fill")))
    fills <- unique(stats::na.omit(unlist(lapply(built$data, function(d) d$fill))))
    expect_gt(length(fills), 1)
})
