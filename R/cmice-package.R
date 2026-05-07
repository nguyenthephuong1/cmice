#' Constrained Multiple Imputation: named-category interface
#'
#' Imputes missing categorical values under structural constraints. The
#' \code{exclude} argument names categories that are inadmissible and will
#' receive zero predicted probability before sampling. This is the recommended
#' interface; it validates labels against \code{levels(y)} at call time.
#'
#' @param y A factor vector to impute. Must have \code{levels(y)} defined.
#' @param ry Logical vector indicating observed (TRUE) vs missing (FALSE).
#' @param x Predictor matrix (numeric or coerced).
#' @param wy Optional logical vector indicating which observations to impute
#'   (defaults to \code{!ry}).
#' @param exclude Character vector of category labels in \code{levels(y)} that
#'   are structurally inadmissible for the missing observations being imputed.
#'   Every label must appear in \code{levels(y)}; mismatches raise an error.
#' @param nnet.maxit Maximum iterations for \code{nnet::multinom} (default 100).
#' @param nnet.MaxNWts Maximum number of weights for \code{nnet::multinom}
#'   (default 1500).
#' @param ... Additional arguments passed to \code{nnet::multinom}.
#'
#' @return A character vector of imputed values for the missing entries.
#'
#' @seealso \code{\link{mice.impute.cons}} (legacy positional interface).
#'
#' @examples
#' \dontrun{
#' # Define a custom imputation method that excludes "SMC" by name
#' mice.impute.cons_no_smc <- function(y, ry, x, wy = NULL, ...) {
#'   mice.impute.cons_named(y, ry, x, wy, exclude = c("SMC"), ...)
#' }
#' }
#' @export
mice.impute.cons_named <- function(y, ry, x, wy = NULL,
                                   exclude = character(0),
                                   nnet.maxit = 100,
                                   nnet.MaxNWts = 1500,
                                   ...) {
  if (is.null(wy)) wy <- !ry

  # --- Validate exclude argument ---
  if (!is.character(exclude)) {
    stop("`exclude` must be a character vector of category labels. ",
         "For positional/integer constraints use mice.impute.cons() instead.")
  }
  lev <- levels(as.factor(y))
  bad <- setdiff(exclude, lev)
  if (length(bad)) {
    stop("`exclude` contains labels not in levels(y): ",
         paste(shQuote(bad), collapse = ", "), ". ",
         "Available levels: ", paste(shQuote(lev), collapse = ", "), ".")
  }
  # Forward to the legacy implementation with character constraint
  mice.impute.cons(y = y, ry = ry, x = x, wy = wy,
                   constrain = exclude,
                   nnet.maxit = nnet.maxit,
                   nnet.MaxNWts = nnet.MaxNWts,
                   ...)
}
