test_that("top() returns the top n elements without groups", {
    df <- data.frame(
        x = c("A", "B", "C", "D", "E"),
        y = c(10, 20, 30, 40, 50)
    )
    result <- top(3, data = df, in_form = "wide", output = "data")
    expect_equal(result$x, c("A", "B", "C"))
})

test_that("top() returns the top n elements with groups", {
    df <- data.frame(
        group = c("A", "A", "B", "B", "C"),
        x = c("A", "B", "C", "D", "E"),
        y = c(10, 20, 30, 40, 50)
    )
    result <- top(1, data = df, groups = "group", in_form = "wide", output = "data")
    expect_equal(result$x, c("A", "C", "E"))

    result <- top(1, data = df, groups = group, in_form = "wide", output = "data")
    expect_equal(result$x, c("A", "C", "E"))
})

test_that("top() works with ordering", {
    df <- data.frame(
        group = c("A", "A", "B", "B", "B"),
        x = c("A", "B", "C", "D", "E"),
        y = c(10, 20, 30, 40, 50)
    )
    result <- top(1, data = df, group = group, order = desc(y),
        in_form = "wide", output = "data")
    expect_equal(result$x, c("E", "B"))
})

test_that("top() respects output", {
    df <- data.frame(
        group = c("A", "A", "B", "B"),
        x = c("A", "B", "C", "D"),
        y = c(10, 20, 30, 40)
    )
    result <- top(1, data = df, groups = group, id = x,
        in_form = "wide", output = "id")
    expect_equal(result, c("A", NA, "C", NA))
})

test_that("top() works with dplyr::mutate()", {
    df <- data.frame(
        CTaa = c(rep("A", 10), rep("B", 20), rep("C", 30)),
        group = c(rep("A", 10), rep("B", 50))
    )
    result <- dplyr::mutate(df, TopClones = top(1))
    expect_equal(result$TopClone, c(rep(NA, 30), rep("C", 30)))

    result <- dplyr::mutate(df, TopClones = top(1, groups = "group"))
    expect_equal(result$TopClones, c(rep("A", 10), rep(NA, 20), rep("C", 30)))
})

test_that("top() works with given environment", {
    data <- data.frame(
        group = c("A", "A", "B", "B"),
        x = c("A", "B", "C", "D"),
        y = c(10, 20, 30, 40)
    )
    result <- top(1)
    expect_equal(result$x, c("A"))

    facet_by <- "group"
    result <- top(1, group = group)
    expect_equal(result$x, c("A", "C"))
})
