# Unit tests for the named-category interface to CMICE.
#
# Verify that mice.impute.cons_named():
#   1. Validates the exclude argument is character (not c(int, char) coercion)
#   2. Errors clearly on unknown labels with a helpful message
#   3. Excludes the named category from imputed values (no leakage)
#   4. Preserves the standard mice.impute.* interface (length == sum(wy))
#   5. Returns values that are all in levels(y) \ exclude

test_that("mice.impute.cons_named catches mixed integer/character exclude (R coerces to char)", {
  # In R, c(3L, "D") coerces to c("3", "D"). Our validator correctly
  # surfaces this as "label '3' not in levels(y)" -- which is the right
  # behaviour: if you pass integer positions, they will not match factor labels.
  set.seed(1)
  y  <- factor(sample(c("A", "B", "C", "D"), 100, replace = TRUE))
  ry <- !is.na(y); ry[1:20] <- FALSE
  x  <- matrix(rnorm(100 * 2), ncol = 2)

  expect_error(
    mice.impute.cons_named(y, ry, x, exclude = c(3L, "D")),
    regexp = "labels not in levels\\(y\\)"
  )
})

test_that("mice.impute.cons_named rejects truly non-character exclude (e.g., integer-only)", {
  set.seed(1)
  y  <- factor(sample(c("A", "B", "C", "D"), 100, replace = TRUE))
  ry <- !is.na(y); ry[1:20] <- FALSE
  x  <- matrix(rnorm(100 * 2), ncol = 2)

  expect_error(
    mice.impute.cons_named(y, ry, x, exclude = c(1L, 2L)),
    regexp = "must be a character vector"
  )
})

test_that("mice.impute.cons_named errors on unknown label", {
  set.seed(2)
  y  <- factor(sample(c("ADC", "SQC", "SMC"), 100, replace = TRUE))
  ry <- !is.na(y); ry[1:20] <- FALSE
  x  <- matrix(rnorm(100 * 2), ncol = 2)

  expect_error(
    mice.impute.cons_named(y, ry, x, exclude = c("SMc")),
    regexp = "labels not in levels\\(y\\).*SMc"
  )
})

test_that("mice.impute.cons_named excludes named category (no leakage)", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("mice")

  set.seed(3)
  n <- 200
  y  <- factor(sample(c("ADC", "SQC", "LAC", "SMC"), n, replace = TRUE),
               levels = c("ADC", "SQC", "LAC", "SMC"))
  ry <- !is.na(y); ry[1:50] <- FALSE
  x  <- matrix(rnorm(n * 3), ncol = 3)

  imp <- mice.impute.cons_named(y, ry, x, exclude = c("SMC"))

  expect_equal(length(imp), 50)
  expect_false(any(as.character(imp) == "SMC"))
  expect_true(all(as.character(imp) %in% c("ADC", "SQC", "LAC")))
})

test_that("mice.impute.cons_named handles empty exclude (= unconstrained)", {
  skip_if_not_installed("nnet")
  set.seed(4)
  n <- 150
  y  <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  ry <- !is.na(y); ry[1:30] <- FALSE
  x  <- matrix(rnorm(n * 2), ncol = 2)

  imp <- mice.impute.cons_named(y, ry, x, exclude = character(0))
  expect_equal(length(imp), 30)
  expect_true(all(as.character(imp) %in% levels(y)))
})
