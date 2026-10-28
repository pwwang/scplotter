test_that("CellStatPlot() expands a single 'group_by' column for sankey plots", {
    df <- data.frame(
        cluster = c(rep("A", 12), rep("B", 8), rep("C", 5)),
        stim = c(rep("CTRL", 10), rep("STIM", 2), rep("STIM", 8), rep("CTRL", 5))
    )
    p <- CellStatPlot(df, ident = "cluster", group_by = "stim:expand", plot_type = "sankey")
    expect_s3_class(p, "ggplot")

    built <- ggplot2::ggplot_build(p)
    idx <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomStratum"), logical(1)))
    st <- built$data[[idx[1]]]
    # stages are ordered by factor levels: 1 = CTRL, 2 = STIM
    ctrl_strata <- as.character(st$stratum[st$x == 1])
    stim_strata <- as.character(st$stratum[st$x == 2])
    expect_true(all(c("A", "C") %in% ctrl_strata))
    expect_false("B" %in% ctrl_strata)
    expect_true(all(c("A", "B") %in% stim_strata))
    expect_false("C" %in% stim_strata)

    expect_no_error(
        CellStatPlot(df, ident = "cluster", group_by = "stim:expand",
                     plot_type = "sankey", frac = "group")
    )
})

test_that("CellStatPlot() validates the ':expand' modifier", {
    df <- data.frame(
        cluster = rep(c("A", "B"), each = 5),
        stim = rep(c("CTRL", "STIM"), each = 5)
    )
    expect_error(
        CellStatPlot(df, ident = "cluster", group_by = c("stim:expand", "cluster"),
                     plot_type = "sankey"),
        "single 'group_by' column"
    )
    expect_error(
        CellStatPlot(df, ident = "stim", group_by = "stim:expand", plot_type = "sankey"),
        "cannot be the same as 'ident'"
    )
})
