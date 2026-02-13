test_that("and()/or() works with dplyr::mutate()", {
    set.seed(8525)
    df <- data.frame(
        CTaa = c(rep("A", 10), rep("B", 20), rep("C", 30)),
        group = c(sample(c("X", "Y"), 30, replace = TRUE), c(rep("X", 15), rep("Y", 15))),
        group2 = c(rep("M", 20), rep("N", 20), rep("O", 20))
    )
    # CTaa      X     Y
    # <chr> <int> <int>
    # 1 A         7     3
    # 2 B         9    11
    # 3 C        15    15
    # and()
    result <- dplyr::mutate(df, AndClones = and(eq(Y, X, groups = "group"), ge(Y, X, groups = "group")))
    expect_equal(result$AndClones, c(rep(NA, 10), rep(NA, 20), rep("C", 30)))

    # or()
    result <- dplyr::mutate(df, OrClones = or(eq(Y, X, groups = "group"), ge(Y, X, groups = "group")))
    expect_equal(result$OrClones, c(rep(NA, 10), rep("B", 20), rep("C", 30)))
})

test_that("and()/or() reports error when elements have different lengths", {
    expect_error(and(c(1, 2), c(1, 2, 3)))
    expect_error(or(c(1, 2), c(1, 2, 3)))
})

test_that("and()/or() respects the output of the input selectors", {
    expect_equal(and(c(TRUE, FALSE, TRUE), c(FALSE, TRUE, TRUE)), c(FALSE, FALSE, TRUE))
    expect_equal(or(c(TRUE, FALSE, TRUE), c(FALSE, TRUE, TRUE)), c(TRUE, TRUE, TRUE))
    expect_equal(and(c("A", "B", "C"), c("A", "B", NA)), c("A", "B", NA))
    expect_equal(or(c("A", NA, "C"), c("A", "B", NA)), c("A", "B", "C"))
})
