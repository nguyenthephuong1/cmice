# dev/project-specific-and-legacy.R
#
# NOT part of the installed package (this whole dev/ directory is .Rbuildignore'd).
# Preserved for provenance / reuse. Two kinds of content:
#
#   1. Project-specific convenience wrappers (mice.impute.cons_na1 / _na2) that
#      hardcode the constraint sets used in the CMICE Japanese lung-cancer
#      analysis. These are NOT general package API — the admissible-set
#      constraints are dataset-specific — so they live here rather than in R/.
#      To use them, source this file in the analysis project, or call
#      mice.impute.cons_named() / mice.impute.cons() directly with the
#      appropriate `exclude` / `constrain` argument.
#
#   2. The superseded dplyr-based implementation of mice.impute.cons(), kept
#      for historical reference. The current engine (R/mice-impute-cons.R)
#      replaced it to drop the tidyverse + mice:::augment dependencies.

## ── 1. Project-specific constraint wrappers ──────────────────────────────────

### Exclude categories {6, "NA2"} (CMICE Japan lung-cancer config A)
mice.impute.cons_na1 <- function (y, ry, x, wy = NULL, nnet.maxit = 100,
                                  nnet.trace = FALSE,
                                  nnet.MaxNWts = 1500, ...) {
  mice.impute.cons(y, ry, x, wy,
                   constrain = c(6, "NA2"),
                   nnet.maxit = nnet.maxit,
                   nnet.trace = nnet.trace,
                   nnet.MaxNWts = nnet.MaxNWts)
}

### Exclude categories {1, 6, "NA1"} (CMICE Japan lung-cancer config B)
mice.impute.cons_na2 <- function (y, ry, x, wy = NULL, nnet.maxit = 100,
                                  nnet.trace = FALSE,
                                  nnet.MaxNWts = 1500, ...) {
  mice.impute.cons(y, ry, x, wy,
                   constrain = c(1, 6, "NA1"),
                   nnet.maxit = nnet.maxit,
                   nnet.trace = nnet.trace,
                   nnet.MaxNWts = nnet.MaxNWts)
}


## ── 2. Legacy dplyr-based implementation (superseded) ────────────────────────
# Superseded by the dependency-free engine in R/mice-impute-cons.R.
# Depended on dplyr (%>%, select, mutate) and the non-exported mice:::augment().
#
# mice.impute.cons <- function (y, ry, x, wy = NULL, nnet.maxit = 100,
#                               nnet.trace = FALSE,
#                               constrain = NULL,
#                               nnet.MaxNWts = 1500, ...) {
#   install.on.demand("nnet", "tidyverse",...)
#   if (is.null(wy)) {
#     wy <- !ry
#   }
#   x <- as.matrix(x)
#   aug <- augment(y, ry, x, wy)
#   x <- aug$x
#   y <- aug$y
#   ry <- aug$ry
#   wy <- aug$wy
#   w <- aug$w
#   t <- y %ni% constrain
#   fy <- as.factor(as.character(y)[t])
#   nc <- length(levels(fy))
#   un <- rep(runif(sum(wy)), each = nc)
#   xy <- cbind.data.frame(y = y, x = x)
#   if (ncol(x) == 0L) {
#     xy <- data.frame(xy, int = 1)
#   }
#   cat.has.all.obs <- table(y[ry]) == sum(ry)
#   if (any(cat.has.all.obs)) {
#     return(rep(levels(fy)[cat.has.all.obs], sum(wy)))
#   }
#   fit <- nnet::multinom(formula(xy), data = xy[ry, , drop = FALSE],
#                         weights = w[ry],
#                         maxit = nnet.maxit,
#                         trace = nnet.trace,
#                         MaxNWts = nnet.MaxNWts,
#                         ...
#   )
#   post <- predict(fit, xy[wy, , drop = FALSE], type = "probs")
#   if (sum(wy) == 1) {
#     post <- matrix(post, nrow = 1, ncol = length(post))
#   }
#   if (is.vector(post)) {
#     post <- matrix(c(1 - post, post), ncol = 2)
#   }
#   post.df <- as.data.frame(post) %>%
#     dplyr::select(-constrain) %>%
#     mutate(sum = rowSums(.)) %>%
#     mutate_all(~./sum) %>%
#     dplyr::select(-sum)
#
#   draws <- un > apply(post.df, 1, cumsum)
#   idx <- 1 + apply(draws, 2, sum)
#   levels(fy)[idx]
# }
