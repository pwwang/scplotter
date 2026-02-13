test_that("sel() returns selected elements without groups", {
    df <- data.frame(
        x = c("A", "B", "C", "D", "E"),
        y = c(10, 20, 30, 40, 50)
    )
    result <- sel(y < 20, data = df, in_form = "wide", output = "data")
    expect_equal(result$x, c("A"))

    # works with a character expression as well
    result <- sel("y == 10", data = df, in_form = "wide", output = "data")
    expect_equal(result$x, c("A"))
})

test_that("sel() returns selected elements with groups", {
    df <- data.frame(
        group = c("A", "A", "B", "B", "C"),
        x = c("A", "B", "C", "D", "E"),
        y = c(50, 20, 30, 40, 50)
    )
    result <- sel(y < 30, data = df, groups = "group", in_form = "wide", output = "data")
    expect_equal(result$x, c("B"))

    result <- sel(mean(y) == 35, data = df, groups = group, in_form = "wide", output = "data")
    expect_equal(result$x, c("A", "B", "C", "D"))
})

test_that("sel() respects output", {
    df <- data.frame(
        group = c("A", "A", "B", "B"),
        x = c("A", "B", "C", "D"),
        y = c(10, 20, 30, 40)
    )
    result <- sel(y < 30, data = df, groups = "group", id = x,
                  in_form = "wide", output = "id")
    expect_equal(result, c("A", "B", NA, NA))
})

test_that("sel() works with dplyr::mutate()", {
    df <- data.frame(
        CTaa = c(rep("A", 10), rep("B", 20), rep("C", 30)),
        group = c(rep("X", 10), rep("Y", 50)),
        y = c(rep(10, 10), rep(20, 20), rep(30, 30))
    )
    # CTaa
    result <- dplyr::mutate(df, SelectedClones = sel(.n > 20))
    expect_equal(result$SelectedClones, c(rep(NA, 30), rep("C", 30)))

    result <- dplyr::mutate(df, SelectedClones = sel(X == 0, groups = "group"))
    expect_equal(result$SelectedClones, c(rep(NA, 10), rep("B", 20), rep("C", 30)))
})

test_that("sel() works with given environment", {
    data <- data.frame(
        group = c("A", "A", "B", "B"),
        x = c("A", "B", "C", "D"),
        y = c(10, 20, 30, 40)
    )
    result <- sel(y < 30)
    expect_equal(result$x, c("A", "B"))

    facet_by <- "group"
    result <- sel(mean(y) >= 25)
    expect_equal(result$x, c("C", "D"))
})
