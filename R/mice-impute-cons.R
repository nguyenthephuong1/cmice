####### Functions of CMICE
####### For data imputation
####### Updated: 2025/05/08




################################################################################
#########           METHOD CMICE -  PREDICTION  CONSTRAINS            ##########
################################################################################




### Internal data-augmentation helper (White, Daniel & Royston, CSDA 2010),
### bundled here so the package does not depend on the non-exported mice:::augment().
.cmice_augment <- function(y, ry, x, wy, maxcat = 50) {
  icod <- sort(unique(unclass(y)))
  k <- length(icod)
  if (k > maxcat) stop("Maximum number of categories (", maxcat, ") exceeded")
  p <- ncol(x)
  if (p == 0) return(list(y = y, ry = ry, x = x, wy = wy, w = rep(1, length(y))))
  if (sum(!ry) == 1) return(list(y = y, ry = ry, x = x, wy = wy, w = rep(1, length(y))))

  mean <- apply(x, 2, mean, na.rm = TRUE)
  sd   <- sqrt(apply(x, 2, var, na.rm = TRUE))
  minx <- apply(x, 2, min, na.rm = TRUE)
  maxx <- apply(x, 2, max, na.rm = TRUE)
  nr <- 2 * p * k
  a <- matrix(mean, nrow = nr, ncol = p, byrow = TRUE)
  b <- matrix(rep(c(rep.int(c(0.5, -0.5), k), rep.int(0, nr)), length = nr * p),
              nrow = nr, ncol = p, byrow = FALSE)
  cc <- matrix(sd, nrow = nr, ncol = p, byrow = TRUE)
  d <- a + b * cc
  d <- pmax(matrix(minx, nrow = nr, ncol = p, byrow = TRUE), d, na.rm = TRUE)
  d <- pmin(matrix(maxx, nrow = nr, ncol = p, byrow = TRUE), d, na.rm = TRUE)
  e <- rep(rep(icod, each = 2), p)

  dimnames(d) <- list(paste0("AUG", seq_len(nrow(d))), dimnames(x)[[2]])
  xa  <- rbind.data.frame(x, d)
  ya  <- if (is.factor(y)) as.factor(levels(y)[c(y, e)]) else c(y, e)
  rya <- c(ry, rep.int(TRUE, nr))
  wya <- c(wy, rep.int(FALSE, nr))
  wa  <- c(rep.int(1, length(y)), rep.int((p + 1) / nr, nr))
  list(y = ya, ry = rya, x = xa, w = wa, wy = wya)
}

### Core function (global constraint)
mice.impute.cons <- function(y, ry, x, wy = NULL,
                             nnet.maxit = 100,
                             nnet.trace = FALSE,
                             constrain = NULL,
                             nnet.MaxNWts = 1500,
                             ...) {
  # --- Safety checks ---
  if (is.null(wy)) wy <- !ry
  if (!any(wy)) return(y[wy]) # nothing to impute
  
  # Ensure y is character/factor-like
  y_chr <- as.character(y)
  
  # Ensure constrain is character vector (match y levels by label)
  if (!is.null(constrain)) constrain <- as.character(constrain)
  
  # Load required dependency (do not install)
  if (!requireNamespace("nnet", quietly = TRUE)) {
    stop("Package 'nnet' is required for mice.impute.cons().")
  }
  
  # Data augmentation via the bundled helper (avoids the non-exported mice:::augment()).
  aug <- .cmice_augment(y_chr, ry, as.matrix(x), wy)
  
  x   <- aug$x
  y_c <- aug$y
  ry  <- aug$ry
  wy  <- aug$wy
  w   <- aug$w
  
  # Define outcome levels from observed data (but keep stable ordering)
  y_obs <- y_c[ry]
  lev_all <- sort(unique(y_obs))
  
  # If observed data has only one class, impute that class (unless constrained)
  if (length(lev_all) == 1L) {
    if (!is.null(constrain) && lev_all %in% constrain) {
      # No valid class to impute; return NA to make failure explicit
      return(rep(NA_character_, sum(wy)))
    }
    return(rep(lev_all, sum(wy)))
  }
  
  # Build modeling data frame
  xy <- data.frame(y = factor(y_c, levels = lev_all), x, check.names = FALSE)
  
  # If no predictors, add intercept
  if (ncol(x) == 0L) xy$int <- 1
  
  # Fit multinomial on observed rows
  # Use y ~ . (all predictors in xy except y)
  fit <- nnet::multinom(
    y ~ .,
    data    = xy[ry, , drop = FALSE],
    weights = w[ry],
    maxit   = nnet.maxit,
    trace   = nnet.trace,
    MaxNWts = nnet.MaxNWts,
    ...
  )
  
  # Predict class probabilities for missing rows
  post <- predict(fit, newdata = xy[wy, , drop = FALSE], type = "probs")
  
  # Normalize shapes for edge cases
  if (sum(wy) == 1L) {
    post <- matrix(post, nrow = 1L)
  }
  if (is.vector(post)) {
    # Binary case may return vector of P(class2); reconstruct 2-column matrix
    # Determine column names from model levels if possible
    if (length(lev_all) != 2L) {
      stop("Unexpected post vector with non-binary outcome.")
    }
    post <- cbind(1 - post, post)
    colnames(post) <- lev_all
  } else {
    # Ensure columns aligned with lev_all
    # predict(multinom) often returns columns named by levels(fit)
    cn <- colnames(post)
    if (!is.null(cn)) {
      # Build full probability matrix in lev_all order
      P <- matrix(0, nrow = nrow(post), ncol = length(lev_all))
      colnames(P) <- lev_all
      P[, intersect(lev_all, cn)] <- post[, intersect(lev_all, cn), drop = FALSE]
      post <- P
    } else {
      # If no colnames, assume same order as lev_all
      if (ncol(post) != length(lev_all)) {
        stop("Predicted probability matrix has unexpected dimension.")
      }
      colnames(post) <- lev_all
    }
  }
  
  # Apply global constraints by zeroing disallowed categories
  if (!is.null(constrain) && length(constrain) > 0) {
    bad <- intersect(colnames(post), constrain)
    if (length(bad) > 0) post[, bad] <- 0
  }
  
  # Renormalize row-wise; if a row sums to 0, fallback to uniform over allowed
  rs <- rowSums(post)
  zero_rows <- which(rs <= 0)
  
  if (length(zero_rows) > 0) {
    allowed <- colnames(post)
    if (!is.null(constrain) && length(constrain) > 0) {
      allowed <- setdiff(allowed, intersect(allowed, constrain))
    }
    if (length(allowed) == 0L) {
      return(rep(NA_character_, sum(wy)))
    }
    post[zero_rows, ] <- 0
    post[zero_rows, allowed] <- 1 / length(allowed)
    rs <- rowSums(post)
  }
  
  post <- post / rs
  
  # Sample categories row-wise
  draws <- apply(post, 1L, function(p) {
    sample(colnames(post), size = 1L, prob = p)
  })
  
  as.character(draws)
}





### ADVANCE METHODS

### Edit the methods: specific functions
mice.impute.cons_na1 <- function (y, ry, x, wy = NULL, nnet.maxit = 100,
                                  nnet.trace = FALSE,
                                  nnet.MaxNWts = 1500, ...) {
  mice.impute.cons(y, ry, x, wy,
                   constrain = c(6, "NA2"),
                   nnet.maxit = nnet.maxit,
                   nnet.trace = nnet.trace,
                   nnet.MaxNWts = nnet.MaxNWts)
}


mice.impute.cons_na2 <- function (y, ry, x, wy = NULL, nnet.maxit = 100,
                                  nnet.trace = FALSE,
                                  nnet.MaxNWts = 1500, ...) {
  mice.impute.cons(y, ry, x, wy,
                   constrain = c(1, 6, "NA1"),
                   nnet.maxit = nnet.maxit,
                   nnet.trace = nnet.trace,
                   nnet.MaxNWts = nnet.MaxNWts)
}


### End







# ### Edit the methods: General functions
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
