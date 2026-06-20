# Unit tests for the core (positional/global constraint) engine
# mice.impute.cons(). The named-category interface (test-cons-named.R)
# exercises this engine only indirectly via delegation; these tests target
# its own branches directly:
#
#   1. Standard mice.impute.* contract: length == sum(wy); values in levels(y)
#   2. wy defaults to !ry, and an explicit wy is respected
#   3. Early return when there is nothing to impute (!any(wy))
#   4. Global constraint masking: constrained categories never imputed
#   5. Single observed class -> that class is imputed (unconstrained)
#   6. Single observed class that IS constrained -> NA_character_
#   7. All categories constrained (zero-row fallback) -> NA_character_
#   8. Binary outcome path (predict() may return a probability vector)
#   9. Single missing row (sum(wy) == 1) matrix-shape edge case
#  10. Reproducibility under a fixed seed

test_that("mice.impute.cons respects the standard interface (length, levels)", {
  skip_if_not_installed("nnet")
  set.seed(101)
  n  <- 200
  y  <- factor(sample(c("ADC", "SQC", "LAC", "SMC"), n, replace = TRUE),
               levels = c("ADC", "SQC", "LAC", "SMC"))
  ry <- rep(TRUE, n); ry[1:40] <- FALSE
  x  <- matrix(rnorm(n * 3), ncol = 3)

  imp <- mice.impute.cons(y, ry, x)

  expect_equal(length(imp), sum(!ry))
  expect_type(imp, "character")
  expect_true(all(imp %in% levels(y)))
})

test_that("mice.impute.cons defaults wy to !ry and honours explicit wy", {
  skip_if_not_installed("nnet")
  set.seed(102)
  n  <- 150
  y  <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  ry <- rep(TRUE, n); ry[1:30] <- FALSE
  x  <- matrix(rnorm(n * 2), ncol = 2)

  # Default: imputes the 30 missing
  imp_default <- mice.impute.cons(y, ry, x, wy = NULL)
  expect_equal(length(imp_default), 30)

  # Explicit wy targeting only the first 10 missing rows
  wy <- rep(FALSE, n); wy[1:10] <- TRUE
  imp_explicit <- mice.impute.cons(y, ry, x, wy = wy)
  expect_equal(length(imp_explicit), 10)
})

test_that("mice.impute.cons returns empty when there is nothing to impute", {
  set.seed(103)
  n  <- 50
  y  <- factor(sample(c("A", "B"), n, replace = TRUE))
  ry <- rep(TRUE, n)            # nothing missing
  x  <- matrix(rnorm(n * 2), ncol = 2)

  imp <- mice.impute.cons(y, ry, x)   # wy <- !ry is all FALSE
  expect_equal(length(imp), 0L)
})

test_that("mice.impute.cons masks globally constrained categories (no leakage)", {
  skip_if_not_installed("nnet")
  set.seed(104)
  n  <- 300
  y  <- factor(sample(c("ADC", "SQC", "LAC", "SMC"), n, replace = TRUE),
               levels = c("ADC", "SQC", "LAC", "SMC"))
  ry <- rep(TRUE, n); ry[1:80] <- FALSE
  x  <- matrix(rnorm(n * 3), ncol = 3)

  imp <- mice.impute.cons(y, ry, x, constrain = c("SMC", "LAC"))

  expect_equal(length(imp), 80)
  expect_false(any(imp %in% c("SMC", "LAC")))
  expect_true(all(imp %in% c("ADC", "SQC")))
})

test_that("mice.impute.cons imputes the only observed class when unconstrained", {
  set.seed(105)
  n  <- 60
  y  <- factor(rep("A", n), levels = c("A", "B", "C"))
  y[1:15] <- NA
  ry <- !is.na(y)
  x  <- matrix(rnorm(n * 2), ncol = 2)

  imp <- mice.impute.cons(y, ry, x)
  expect_equal(length(imp), 15)
  expect_true(all(imp == "A"))
})

test_that("mice.impute.cons returns NA when the only observed class is constrained", {
  set.seed(106)
  n  <- 60
  y  <- factor(rep("A", n), levels = c("A", "B", "C"))
  y[1:15] <- NA
  ry <- !is.na(y)
  x  <- matrix(rnorm(n * 2), ncol = 2)

  imp <- mice.impute.cons(y, ry, x, constrain = "A")
  expect_equal(length(imp), 15)
  expect_true(all(is.na(imp)))
})

test_that("mice.impute.cons returns NA when every category is constrained (zero-row fallback)", {
  skip_if_not_installed("nnet")
  set.seed(107)
  n  <- 120
  y  <- factor(sample(c("A", "B", "C"), n, replace = TRUE),
               levels = c("A", "B", "C"))
  ry <- rep(TRUE, n); ry[1:30] <- FALSE
  x  <- matrix(rnorm(n * 2), ncol = 2)

  imp <- mice.impute.cons(y, ry, x, constrain = c("A", "B", "C"))
  expect_equal(length(imp), 30)
  expect_true(all(is.na(imp)))
})

test_that("mice.impute.cons handles a binary outcome (predict() vector path)", {
  skip_if_not_installed("nnet")
  set.seed(108)
  n  <- 200
  y  <- factor(sample(c("yes", "no"), n, replace = TRUE), levels = c("no", "yes"))
  ry <- rep(TRUE, n); ry[1:50] <- FALSE
  x  <- matrix(rnorm(n * 2), ncol = 2)

  imp <- mice.impute.cons(y, ry, x)
  expect_equal(length(imp), 50)
  expect_true(all(imp %in% c("no", "yes")))
})

test_that("mice.impute.cons handles a single missing row (sum(wy) == 1)", {
  skip_if_not_installed("nnet")
  set.seed(109)
  n  <- 120
  y  <- factor(sample(c("A", "B", "C"), n, replace = TRUE),
               levels = c("A", "B", "C"))
  ry <- rep(TRUE, n); ry[1] <- FALSE      # exactly one missing
  x  <- matrix(rnorm(n * 2), ncol = 2)

  imp <- mice.impute.cons(y, ry, x)
  expect_equal(length(imp), 1L)
  expect_true(imp %in% levels(y))
})

test_that("mice.impute.cons is reproducible under a fixed seed", {
  skip_if_not_installed("nnet")
  n  <- 200
  y  <- factor(sample(c("A", "B", "C", "D"), n, replace = TRUE),
               levels = c("A", "B", "C", "D"))
  ry <- rep(TRUE, n); ry[1:50] <- FALSE
  x  <- matrix(rnorm(n * 3), ncol = 3)

  set.seed(2026); a <- mice.impute.cons(y, ry, x)
  set.seed(2026); b <- mice.impute.cons(y, ry, x)
  expect_identical(a, b)
})
