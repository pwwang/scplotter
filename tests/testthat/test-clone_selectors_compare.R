test_that("gt()/ge()/lt()/le() returns selected elements", {
    df <- data.frame(
        x = c("A", "B", "C", "D", "E"),
        g = c("X", "Y", "X", "Y", "Z"),
        g1 = c(10, 20, 30, 40, 50),
        g2 = c(5, 15, 35, 40, 30)
    )
    result <- gt(g1, g2, data = df, in_form = "wide", return_ids = FALSE)
    expect_equal(result$x, c("A", "B", "E"))

    result <- gt(g1, g2, groups = g, data = df, in_form = "wide", return_ids = FALSE)
    expect_equal(result$x, c("A", "B", "E"))

    result <- ge(g1, g2, data = df, in_form = "wide", return_ids = FALSE)
    expect_equal(result$x, c("A", "B", "D", "E"))

    result <- ge(g1, g2, groups = g, data = df, in_form = "wide", return_ids = FALSE)
    expect_equal(result$x, c("A", "B", "D", "E"))

    result <- lt(g1, g2, data = df, in_form = "wide", return_ids = FALSE)
    expect_equal(result$x, c("C"))

    result <- lt(g1, g2, groups = g, data = df, in_form = "wide", return_ids = FALSE)
    expect_equal(result$x, c("C"))

    result <- le(g1, g2, data = df, in_form = "wide", return_ids = FALSE)
    expect_equal(result$x, c("C", "D"))

    result <- le(g1, g2, groups = g, data = df, in_form = "wide", return_ids = FALSE)
    expect_equal(result$x, c("C", "D"))

    result <- eq(g1, g2, data = df, in_form = "wide", return_ids = FALSE)
    expect_equal(result$x, c("D"))

    result <- eq(g1, g2, groups = g, data = df, in_form = "wide", return_ids = FALSE)
    expect_equal(result$x, c("D"))

    result <- ne(g1, g2, data = df, in_form = "wide", return_ids = FALSE)
    expect_equal(result$x, c("A", "B", "C", "E"))

    result <- ne(g1, g2, groups = g, data = df, in_form = "wide", return_ids = FALSE)
    expect_equal(result$x, c("A", "B", "C", "E"))
})

test_that("gt()/ge()/lt()/le() respects return_ids", {
    df <- data.frame(
        x = c("A", "B", "C", "D"),
        g = c("X", "Y", "X", "Y"),
        g1 = c(10, 20, 30, 40),
        g2 = c(5, 15, 30, 45)
    )
    result <- gt(g1, g2, id = x,
                data = df, in_form = "wide", return_ids = TRUE)
    expect_equal(result, c("A", "B", NA, NA))

    result <- ge(g1, g2, id = x,
                data = df, in_form = "wide", return_ids = TRUE)
    expect_equal(result, c("A", "B", "C", NA))

    result <- lt(g1, g2, id = x,
                data = df, in_form = "wide", return_ids = TRUE)
    expect_equal(result, c(NA, NA, NA, "D"))

    result <- le(g1, g2, id = x,
                data = df, in_form = "wide", return_ids = TRUE)
    expect_equal(result, c(NA, NA, "C", "D"))

    result <- eq(g1, g2, id = x,
                data = df, in_form = "wide", return_ids = TRUE)
    expect_equal(result, c(NA, NA, "C", NA))

    result <- ne(g1, g2, id = x,
                data = df, in_form = "wide", return_ids = TRUE)
    expect_equal(result, c("A", "B", NA, "D"))
})

test_that("gt()/ge()/lt()/le() works with dplyr::mutate()", {
    set.seed(8525)
    df <- data.frame(
        CTaa = c(rep("A", 10), rep("B", 20), rep("C", 30)),
        group = c(sample(c("X", "Y"), 30, replace = TRUE), c(rep("X", 15), rep("Y", 15))),
        group2 = c(rep("M", 20), rep("N", 20), rep("O", 20))
    )
    result <- dplyr::mutate(df, HighClones = gt(Y, X, groups = "group"))
    expect_equal(result$HighClones, c(rep(NA, 10), rep("B", 20), rep(NA, 30)))

    result <- dplyr::mutate(df, HighClones = ge(Y, X, groups = "group"))
    expect_equal(result$HighClones, c(rep(NA, 10), rep("B", 20), rep("C", 30)))

    result <- dplyr::mutate(df, HighClones = gt(Y, X, groups = c("group", "group2")))
    expect_equal(result$HighClones, c(rep(NA, 10), rep("B", 10), rep(NA, 20), rep("C", 20)))

    result <- dplyr::mutate(df, HighClones = ge(Y, X, groups = c("group", "group2")))
    expect_equal(result$HighClones, c(rep(NA, 10), rep("B", 20), rep(NA, 10), rep("C", 20)))

    result <- dplyr::mutate(df, LowClones = lt(Y, X, groups = "group"))
    expect_equal(result$LowClones, c(rep("A", 10), rep(NA, 50)))

    result <- dplyr::mutate(df, LowClones = le(Y, X, groups = "group"))
    expect_equal(result$LowClones, c(rep("A", 10), rep(NA, 20), rep("C", 30)))

    result <- dplyr::mutate(df, LowClones = lt(Y, X, groups = c("group", "group2")))
    expect_equal(result$LowClones, c(rep("A", 10), rep(NA, 20), rep("C", 10), rep(NA, 20)))

    result <- dplyr::mutate(df, LowClones = le(Y, X, groups = c("group", "group2")))
    expect_equal(result$LowClones, c(rep("A", 10), rep(NA, 10), rep("B", 10), rep("C", 10), rep(NA, 20)))

    result <- dplyr::mutate(df, EqualClones = eq(Y, X, groups = "group"))
    expect_equal(result$EqualClones, c(rep(NA, 30), rep("C", 30)))

    result <- dplyr::mutate(df, EqualClones = ne(Y, X, groups = "group"))
    expect_equal(result$EqualClones, c(rep("A", 10), rep("B", 20), rep(NA, 30)))

    result <- dplyr::mutate(df, EqualClones = eq(Y, X, groups = c("group", "group2")))
    expect_equal(result$EqualClones, c(rep(NA, 20), rep("B", 10), rep(NA, 30)))

    result <- dplyr::mutate(df, EqualClones = ne(Y, X, groups = c("group", "group2")))
    expect_equal(result$EqualClones, c(rep("A", 10), rep("B", 10), rep(NA, 10), rep("C", 30)))
})

test_that("gt()/ge() works with given environment", {
    data <- data.frame(
        group = c("A", "A", "B", "B"),
        x = c("A", "B", "C", "D"),
        y = c(10, 20, 30, 40),
        z = c(5, 15, 30, 45)
    )
    result <- gt(y, z)
    expect_equal(result$x, c("A", "B"))

    result <- lt(y, z)
    expect_equal(result$x, c("D"))

    result <- eq(y, z)
    expect_equal(result$x, c("C"))

    result <- ne(y, z)
    expect_equal(result$x, c("A", "B", "D"))

    facet_by <- "group"
    result <- ge(y, z)
    expect_equal(result$x, c("A", "B", "C"))

    result <- le(y, z)
    expect_equal(result$x, c("C", "D"))

    result <- eq(y, z)
    expect_equal(result$x, c("C"))

    result <- ne(y, z)
    expect_equal(result$x, c("A", "B", "D"))
})
