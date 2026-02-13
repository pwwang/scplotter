test_that("shared() returns selected elements", {
    df <- data.frame(
        x = c("A", "B", "C", "D", "E"),
        g = c("X", "Y", "X", "Y", "Z"),
        g1 = c(10, 20, 30, 40, 50),
        g2 = c(5, 15, 25, 0, 30),
        g3 = c(1, 2, 3, 0, 0)
    )
    result <- shared(g1, g2, data = df, in_form = "wide", output = "data")
    expect_equal(result$x, c("A", "B", "C", "E"))

    result <- shared(g1, g2, groups = g, data = df, in_form = "wide", output = "data")
    expect_equal(result$x, c("A", "B", "C", "E"))

    result <- shared(g1, g2, g3, data = df, in_form = "wide", output = "data")
    expect_equal(result$x, c("A", "B", "C"))
})

test_that("shared() respects output", {
    df <- data.frame(
        x = c("A", "B", "C", "D"),
        g = c("X", "Y", "X", "Y"),
        g1 = c(10, 20, 30, 40),
        g2 = c(5, 15, 25, 0)
    )
    result <- shared(g1, g2, id = x,
                    data = df, in_form = "wide", output = "id")
    expect_equal(result, c("A", "B", "C", NA))
})

test_that("shared() works with dplyr::mutate()", {
    set.seed(8525)
    df <- data.frame(
        CTaa = c(rep("A", 10), rep("B", 20), rep("C", 30)),
        group = c(sample(c("X", "Y"), 30, replace = TRUE), rep("Y", 30)),
        group2 = c(rep("M", 20), rep("N", 20), rep("O", 20))
    )
    result <- dplyr::mutate(df, SharedClones = shared(X, Y, groups = "group"))
    expect_equal(result$SharedClones, c(rep("A", 10), rep("B", 20), rep(NA, 30)))

    result <- dplyr::mutate(df, SharedClones = shared(X, Y, groups = c("group", "group2")))
    expect_equal(result$SharedClones, c(rep("A", 10), rep("B", 20), rep(NA, 30)))
})

test_that("shared() works with given environment", {
    data <- data.frame(
        group = c("A", "A", "B", "B"),
        x = c("A", "B", "C", "D"),
        y = c(10, 20, 30, 40),
        z = c(1, 2, 3, 4)
    )
    result <- shared(y, z)
    expect_equal(result$x, c("A", "B", "C", "D"))

    facet_by <- "group"
    result <- shared(y, z)
    expect_equal(result$x, c("A", "B", "C", "D"))
})