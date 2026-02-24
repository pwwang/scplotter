test_that("sel() returns selected elements without group_by", {
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

test_that("sel() returns selected elements with group_by", {
    df <- data.frame(
        group = c("A", "A", "B", "B", "C"),
        x = c("A", "B", "C", "D", "E"),
        y = c(50, 20, 30, 40, 50)
    )
    result <- sel(y < 30, data = df, group_by = "group", in_form = "wide", output = "data")
    expect_equal(result$x, c("B"))

    result <- sel(mean(y) == 35, data = df, group_by = group, in_form = "wide", output = "data")
    expect_equal(result$x, c("A", "B", "C", "D"))
})

test_that("sel() respects output", {
    df <- data.frame(
        group = c("A", "A", "B", "B"),
        x = c("A", "B", "C", "D"),
        y = c(10, 20, 30, 40)
    )
    result <- sel(y < 30, data = df, group_by = "group", id = x,
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

    result <- dplyr::mutate(df, SelectedClones = sel(X == 0, group_by = "group"))
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

test_that("sel() respects within", {
    df <- data.frame(
        CTaa = c(rep("A", 10), rep("B", 20), rep("C", 30)),
        group = c(rep("A", 25), rep("B", 35))
    )

    result <- dplyr::mutate(df, SelectedClones = sel(.n > 10))
    result <- dplyr::distinct(result, CTaa, group, SelectedClones)
    expect_equal(result$SelectedClones, c(NA, "B", "B", "C"))

    result <- dplyr::mutate(df, SelectedClones = sel(.n > 10, within = group == "B"))
    result <- dplyr::distinct(result, CTaa, group, SelectedClones)
    expect_equal(result$SelectedClones, c(NA, NA, NA, "C"))

    result <- dplyr::mutate(df, SelectedClones = sel(.n > 10, within = group == "A"))
    result <- dplyr::distinct(result, CTaa, group, SelectedClones)
    expect_equal(result$SelectedClones, c(NA, "B", "B", NA))
})

test_that("sel() respects output_within", {
    df <- data.frame(
        CTaa = c(rep("A", 10), rep("B", 20), rep("C", 30)),
        group = c(rep("A", 25), rep("B", 35))
    )

    result <- dplyr::mutate(df, SelectedClones = sel(.n > 10, within = group == "A", output_within = TRUE))
    result <- dplyr::distinct(result, CTaa, group, SelectedClones)
    expect_equal(result$SelectedClones, c(NA, "B", NA, NA))

    result <- dplyr::mutate(df, SelectedClones = sel(.n > 10, within = group == "A", output_within = FALSE))
    result <- dplyr::distinct(result, CTaa, group, SelectedClones)
    expect_equal(result$SelectedClones, c(NA, "B", "B", NA))

    result <- dplyr::mutate(df, SelectedClones = sel(.n > 10, within = group == "A", output_within = group == "B"))
    result <- dplyr::distinct(result, CTaa, group, SelectedClones)
    expect_equal(result$SelectedClones, c(NA, NA, "B", NA))
})

test_that("sel() top argument limits number of selected clones in wide format", {
    # y > 15 matches B(40), C(20), D(50), E(30) -- 4 elements
    df <- data.frame(
        x = c("A", "B", "C", "D", "E"),
        y = c(10, 40, 20, 50, 30)
    )
    # top = 2 without order: first 2 matching rows in original order = B, C
    result <- sel("y > 15", data = df, id = "x", top = 2, in_form = "wide", output = "data")
    expect_equal(result$x, c("B", "C"))

    # top = 3: B, C, D
    result <- sel("y > 15", data = df, id = "x", top = 3, in_form = "wide", output = "data")
    expect_equal(result$x, c("B", "C", "D"))

    # top >= number of matches: all matching rows B, C, D, E
    result <- sel("y > 15", data = df, id = "x", top = 10, in_form = "wide", output = "data")
    expect_equal(result$x, c("B", "C", "D", "E"))
})

test_that("sel() order argument selects top n clones by specified order in wide format", {
    # y > 15 matches B(40), C(20), D(50), E(30)
    df <- data.frame(
        x = c("A", "B", "C", "D", "E"),
        y = c(10, 40, 20, 50, 30)
    )
    # top = 2, order = desc(y): largest y among matches = D(50), B(40)
    # result is in original row order: B then D
    result <- sel("y > 15", data = df, id = "x", top = 2, order = "desc(y)",
                  in_form = "wide", output = "data")
    expect_equal(result$x, c("B", "D"))

    # top = 2, order = y: smallest y among matches = C(20), E(30)
    # result is in original row order: C then E
    result <- sel("y > 15", data = df, id = "x", top = 2, order = "y",
                  in_form = "wide", output = "data")
    expect_equal(result$x, c("C", "E"))

    # top = 2, order = desc(y), output = "id": B and D selected, rest NA
    result <- sel("y > 15", data = df, id = "x", top = 2, order = "desc(y)",
                  in_form = "wide", output = "id")
    expect_equal(result[!is.na(result)], c("B", "D"))
    expect_equal(sum(is.na(result)), 3L)
})

test_that("sel() top argument works with group_by in wide format", {
    df <- data.frame(
        group = c("G1", "G1", "G1", "G2", "G2", "G2"),
        x     = c("A",  "B",  "C",  "D",  "E",  "F"),
        y     = c(50,    10,   30,   20,   40,   15)
    )
    # y > 12 matches all; top = 1 per group: first matching row in each group
    # Group G1 row order: A(50,T), B(10,T), C(30,T) -> top 1 = A
    # Group G2 row order: D(20,T), E(40,T), F(15,T) -> top 1 = D
    result <- sel("y > 12", data = df, id = "x", group_by = "group", top = 1,
                  in_form = "wide", output = "data")
    expect_equal(result$x, c("A", "D"))

    # top = 1, order = desc(y): largest y per group
    # Group G1: A(50) is largest -> top 1 = A
    # Group G2: E(40) is largest -> top 1 = E
    result <- sel("y > 12", data = df, id = "x", group_by = "group", top = 1,
                  order = "desc(y)", in_form = "wide", output = "data")
    expect_equal(result$x, c("A", "E"))
})
