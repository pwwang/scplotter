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

test_that("top() respects within", {
    df <- data.frame(
        CTaa = c(rep("A", 10), rep("B", 20), rep("C", 30)),
        group = c(rep("A", 25), rep("B", 35))
    )
    result <- unique(top(1, data = df, id = "CTaa", in_form = "long",
        within = group == "A", output = "id"))
    expect_equal(result, c(NA, "B"))  # not C
})

test_that("top() respects output_within", {
    df <- data.frame(
        CTaa = c(rep("A", 10), rep("B", 20), rep("C", 30)),
        group = c(rep("A", 25), rep("B", 35))
    )
    result <- dplyr::mutate(df, TopClones = top(1, within = group == "A"))
    result <- dplyr::distinct(result, CTaa, group, TopClones)
    # without output_within, "B" from both groups A and B are returned
    expect_equal(result$TopClone, c(NA, "B", "B", NA))

    result <- dplyr::mutate(df, TopClones = top(1, within = group == "A", output_within = TRUE))
    result <- dplyr::distinct(result, CTaa, group, TopClones)
    # with output_within, only "B" from group A is returned, and "B" from group B is not returned
    expect_equal(result$TopClone, c(NA, "B", NA, NA))

    result <- dplyr::mutate(df, TopClones = top(1, within = group == "A", output_within = FALSE))
    result <- dplyr::distinct(result, CTaa, group, TopClones)
    # with output_within = FALSE, "B" from both groups A and B are returned, same as the case without output_within
    expect_equal(result$TopClone, c(NA, "B", "B", NA))

    result <- dplyr::mutate(df, TopClones = top(1, within = group == "A", output_within = group == "B"))
    result <- dplyr::distinct(result, CTaa, group, TopClones)
    # with output_within = group == "B", only "B" from group B is returned, and "B" from group A is not returned
    expect_equal(result$TopClone, c(NA, NA, "B", NA))
})
