test_that("uniq() returns selected elements", {
    df <- data.frame(
        x = c("A", "B", "C", "D", "E"),
        g = c("X", "Y", "X", "Y", "Z"),
        g1 = c(10, 20, 30, 40, 50),
        g2 = c(5, 15, 25, 0, 0),
        g3 = c(1, 2, 3, 0, 5)
    )
    result <- uniq(g1, g2, data = df, in_form = "wide", return_ids = FALSE)
    expect_equal(result$x, c("D", "E"))

    result <- uniq(g1, g2, groups = g, data = df, in_form = "wide", return_ids = FALSE)
    expect_equal(result$x, c("D", "E"))

    result <- uniq(g1, g2, g3, data = df, in_form = "wide", return_ids = FALSE)
    expect_equal(result$x, c("D"))
})

test_that("uniq() respects return_ids", {
    df <- data.frame(
        x = c("A", "B", "C", "D"),
        g = c("X", "Y", "X", "Y"),
        g1 = c(10, 20, 30, 40),
        g2 = c(5, 15, 25, 0)
    )
    result <- uniq(g1, g2, id = x,
                   data = df, in_form = "wide", return_ids = TRUE)
    expect_equal(result, c(NA, NA, NA, "D"))
})

test_that("uniq() works with dplyr::mutate()", {
    df <- data.frame(
        CTaa = c(rep("A", 10), rep("B", 20), rep("C", 30)),
        group = c(rep("X", 10), rep("Y", 50)),
        group2 = c(rep("M", 20), rep("N", 20), rep("O", 20))
    )
    # CTaa      X     Y
    # <chr> <int> <int>
    # 1 A        10     0
    # 2 B         0    20
    # 3 C         0    30
    result <- dplyr::mutate(df, UniqueClones = uniq(X, Y, groups = "group"))
    expect_equal(result$UniqueClones, c(rep("A", 10), rep(NA, 50)))

    result <- dplyr::mutate(df, UniqueClones = uniq(X, Y, groups = c("group", "group2")))
    expect_equal(result$UniqueClones, c(rep("A", 10), rep(NA, 50)))
})

test_that("uniq() works with given environment", {
    data <- data.frame(
        group = c("A", "A", "B", "B"),
        x = c("A", "B", "C", "D"),
        y = c(10, 20, 30, 40),
        z = c(1, 2, 3, 0)
    )
    result <- uniq(y, z)
    expect_equal(result$x, c("D"))

    facet_by <- "group"
    result <- uniq(y, z)
    expect_equal(result$x, c("D"))
})
