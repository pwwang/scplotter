set.seed(8525)
data(enrich_example, package = "plotthis")
enrich_example$Group <- sample(LETTERS[1:3], nrow(enrich_example), replace = TRUE)

test_that("heatmap colors by one metric and marks by another", {
    p <- EnrichmentPlot(
        enrich_example, plot_type = "heatmap", group_by = "Group",
        metric = "zScore:p.adjust", cutoff = 0.05
    )
    expect_s3_class(p, "ggplot")
})

test_that("heatmap with a single metric keeps the old behavior", {
    p <- EnrichmentPlot(
        enrich_example, plot_type = "heatmap", group_by = "Group",
        metric = "p.adjust", cutoff = 0.05
    )
    expect_s3_class(p, "ggplot")
})

test_that("dot and lollipop grey out by the significance metric", {
    p <- EnrichmentPlot(
        enrich_example, plot_type = "dot",
        metric = "zScore:p.adjust", cutoff = 0.05
    )
    expect_s3_class(p, "ggplot")
    p <- EnrichmentPlot(
        enrich_example, plot_type = "lollipop",
        metric = "zScore:p.adjust", cutoff = 0.05
    )
    expect_s3_class(p, "ggplot")
    p <- EnrichmentPlot(
        enrich_example, plot_type = "comparison", group_by = "Group",
        metric = "zScore:p.adjust", cutoff = 0.05
    )
    expect_s3_class(p, "ggplot")
})

test_that("a term significant in one group is kept in the other groups", {
    d <- data.frame(
        Description = rep(paste0("Term", 1:2), each = 2),
        Group       = rep(c("A", "B"), times = 2),
        GeneRatio   = "10/100", BgRatio = "50/1000", Count = 10, geneID = "G1/G2",
        zScore      = c(5, 4, 3, 2),
        # Term1 passes only in A; Term2 fails everywhere
        p.adjust    = c(1e-6, 0.5, 0.9, 0.9),
        pvalue      = c(1e-6, 0.5, 0.9, 0.9),
        qvalue      = c(1e-6, 0.5, 0.9, 0.9)
    )
    captured <- NULL
    local_mocked_bindings(
        Heatmap = function(...) {
            captured <<- list(...)
            "heatmap"
        },
        .package = "scplotter"
    )
    EnrichmentPlot(d, plot_type = "heatmap", group_by = "Group",
        metric = "zScore:p.adjust", cutoff = 0.05)

    expect_setequal(
        as.character(captured$data$Description),
        c("Term1", "Term1")
    )
    expect_setequal(as.character(captured$data$Group), c("A", "B"))
})

test_that("dot and lollipop keep insignificant terms", {
    d <- data.frame(
        Description = paste0("Term", 1:3),
        GeneRatio   = "10/100", BgRatio = "50/1000", Count = 10, geneID = "G1/G2",
        zScore      = c(5, 4, 3),
        p.adjust    = c(1e-6, 0.5, 0.9),
        pvalue      = c(1e-6, 0.5, 0.9),
        qvalue      = c(1e-6, 0.5, 0.9)
    )
    captured <- NULL
    local_mocked_bindings(
        DotPlot = function(data, ...) {
            captured <<- data
            "dot"
        },
        LollipopPlot = function(data, ...) {
            captured <<- data
            "lollipop"
        },
        .package = "scplotter"
    )
    EnrichmentPlot(d, plot_type = "dot", metric = "zScore:p.adjust", cutoff = 0.05)
    expect_setequal(as.character(captured$Description), paste0("Term", 1:3))

    EnrichmentPlot(d, plot_type = "lollipop", metric = "zScore:p.adjust", cutoff = 0.05)
    expect_setequal(as.character(captured$Description), paste0("Term", 1:3))
})

test_that("bar skips the cutoff line when the metrics differ", {
    expect_warning(
        p <- EnrichmentPlot(enrich_example, plot_type = "bar", metric = "zScore:p.adjust", cutoff = 0.05),
        "No cutoff line is drawn"
    )
    expect_s3_class(p, "ggplot")
})

test_that("heatmap marks are looked up by cell name, not display position", {
    d <- data.frame(
        Description = rep(paste0("Term", 1:4), each = 2),
        Group       = rep(c("A", "B"), times = 4),
        GeneRatio   = "10/100", BgRatio = "50/1000", Count = 10, geneID = "G1/G2",
        zScore      = c(5, -5, 4, -4, 3, -3, 2, -2),
        pvalue      = c(1e-6, 0.5, 0.2, 1e-6, 0.3, 0.4, 1e-6, 0.9),
        # significant: Term1/A, Term2/B, Term4/A
        p.adjust    = c(1e-6, 0.5, 0.9, 1e-6, 0.8, 0.7, 1e-6, 0.9),
        qvalue      = c(1e-6, 0.5, 0.9, 1e-6, 0.8, 0.7, 1e-6, 0.9)
    )
    captured <- NULL
    local_mocked_bindings(
        Heatmap = function(...) {
            captured <<- list(...)
            "heatmap"
        },
        .package = "scplotter"
    )
    EnrichmentPlot(d, plot_type = "heatmap", group_by = "Group",
        metric = "zScore:p.adjust", cutoff = 0.05)

    label <- captured$label
    mark <- function(rname, cname) {
        # indices are deliberately arbitrary: only the names may drive the mark
        r <- label(1, 1, 1, rname, cname)
        if (is.list(r)) r[[1]] else r
    }
    expect_identical(mark("Term1", "A"), "*")
    expect_true(is.na(mark("Term1", "B")))
    expect_identical(mark("Term2", "B"), "*")
    expect_true(is.na(mark("Term2", "A")))
    expect_identical(mark("Term4", "A"), "*")
    expect_true(is.na(mark("Term4", "B")))
})

test_that("top_term selects by the fill metric, after filtering by the cutoff", {
    d <- data.frame(
        Description = paste0("Term", c("A", "B", "C", "D")),
        Group       = "G1",
        GeneRatio   = "10/100", BgRatio = "50/1000", Count = 10, geneID = "G1/G2",
        # largest zScore: A > B > C > D
        zScore      = c(10, 9, 8, 7),
        # only B and C pass the cutoff
        pvalue      = c(0.5, 1e-6, 1e-4, 0.2),
        qvalue      = c(0.5, 1e-6, 1e-4, 0.2),
        p.adjust    = c(0.5, 1e-6, 1e-4, 0.2)
    )
    captured <- list()
    local_mocked_bindings(
        Heatmap = function(...) {
            captured[[length(captured) + 1]] <<- list(...)
            "heatmap"
        },
        .package = "scplotter"
    )
    selected <- function() as.character(captured[[length(captured)]]$data$Description)

    # a separate significance metric filters first, then ranks the survivors by the fill metric
    EnrichmentPlot(d, plot_type = "heatmap", group_by = "Group",
        metric = "zScore:p.adjust", cutoff = 0.05, top_term = 1)
    expect_identical(selected(), "TermB")

    # without a modifier, the fill metric is also the significance metric
    EnrichmentPlot(d, plot_type = "heatmap", group_by = "Group",
        metric = "zScore", top_term = 1)
    expect_identical(selected(), "TermA")

    # p-value metrics select the smallest values
    EnrichmentPlot(d, plot_type = "heatmap", group_by = "Group",
        metric = "pvalue", top_term = 2)
    expect_setequal(selected(), c("TermB", "TermC"))

    # same rule without grouping
    local_mocked_bindings(
        DotPlot = function(data, ...) {
            captured[[length(captured) + 1]] <<- list(data = data, ...)
            "dot"
        },
        .package = "scplotter"
    )
    EnrichmentPlot(d, plot_type = "dot", metric = "zScore", top_term = 1)
    expect_identical(selected(), "TermA")
})

test_that("invalid metrics are rejected", {
    expect_error(
        EnrichmentPlot(enrich_example, plot_type = "heatmap", group_by = "Group", metric = "zScore:p.adjust:qvalue"),
        "must be 'column' or 'column:column'"
    )
    expect_error(
        EnrichmentPlot(enrich_example, plot_type = "heatmap", group_by = "Group", metric = "zScore:not_a_column"),
        "not found in data"
    )
})
