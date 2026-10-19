# Internal helpers ---------------------------------------------------------

#' @importFrom stats pf lm.fit setNames
NULL

.lagmat <- function(x, p) {
  n <- length(x)
  out <- matrix(NA_real_, n - p, p)
  for (l in seq_len(p)) out[, l] <- x[(p - l + 1):(n - l)]
  out
}

.rss <- function(y, X) {
  X <- cbind(1, X)
  fit <- lm.fit(X, y)
  sum(fit$residuals^2)
}

.as_series_frame <- function(x, y = NULL, xname = "x", yname = "y") {
  to_df <- function(v, nm) {
    if (is.data.frame(v)) return(v)
    if (is.matrix(v)) {
      v <- as.data.frame(v)
      if (is.null(colnames(v)) || any(colnames(v) == "")) {
        colnames(v) <- paste0(nm, seq_len(ncol(v)))
      }
      return(v)
    }
    if (inherits(v, "ts") || is.numeric(v)) {
      d <- data.frame(as.numeric(v)); names(d) <- nm; return(d)
    }
    stop("Input must be a numeric vector, matrix or data.frame.", call. = FALSE)
  }
  X <- to_df(x, xname)
  Y <- if (is.null(y)) NULL else to_df(y, yname)
  list(X = X, Y = Y)
}

.prep_na <- function(v) {
  v <- as.numeric(v)
  bad <- !is.finite(v)
  if (!any(bad)) return(v)
  if (all(bad)) stop("A series contains only missing values.", call. = FALSE)
  idx <- seq_along(v)
  v[bad] <- stats::approx(idx[!bad], v[!bad], xout = idx[bad], rule = 2)$y
  v
}
