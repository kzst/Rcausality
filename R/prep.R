#' Prepare a time-series data frame
#'
#' Coerces input to a numeric data frame, linearly interpolates missing
#' values, and optionally differences and/or standardises the columns.
#'
#' @param x numeric vector, matrix or data frame.
#' @param diff integer; order of differencing (0 = none).
#' @param scale logical; standardise to zero mean, unit variance?
#' @param impute logical; linearly interpolate missing values?
#' @return A data frame of numeric columns.
#' @export
#' @examples
#' prep_series(c(1, NA, 3, 4, 6), scale = TRUE)
prep_series <- function(x, diff = 0L, scale = FALSE, impute = TRUE) {
  d <- .as_series_frame(x)$X
  d[] <- lapply(d, function(v) {
    v <- as.numeric(v)
    if (impute) v <- .prep_na(v)
    v
  })
  if (diff > 0L) {
    d <- as.data.frame(lapply(d, function(v) base::diff(v, differences = diff)))
  }
  if (scale) d[] <- lapply(d, function(v) (v - mean(v)) / stats::sd(v))
  d
}

#' Select the autoregressive lag order
#'
#' Chooses the lag order of a (multivariate) VAR fitted by least squares using
#' an information criterion.
#'
#' @param data numeric matrix / data frame (columns are series).
#' @param max_lag largest lag order considered.
#' @param criterion one of `"AIC"`, `"BIC"`, `"HQ"`.
#' @return Integer lag order.
#' @export
#' @examples
#' d <- sim_var(300, scenario = "x_to_y", seed = 1)
#' select_lag(d)
select_lag <- function(data, max_lag = 8L, criterion = c("BIC", "AIC", "HQ")) {
  criterion <- match.arg(criterion)
  M <- as.matrix(data)
  K <- ncol(M); n <- nrow(M)
  max_lag <- max(1L, min(as.integer(max_lag), floor((n - 1) / (2 * K + 2))))
  # common effective sample (n - max_lag) for comparability across orders
  Nn <- n - max_lag
  Y <- M[(max_lag + 1):n, , drop = FALSE]
  crit <- rep(Inf, max_lag)
  for (p in seq_len(max_lag)) {
    Z <- do.call(cbind, lapply(seq_len(p), function(l) M[(max_lag + 1 - l):(n - l), , drop = FALSE]))
    Z <- cbind(1, Z)
    R <- lm.fit(Z, Y)$residuals
    S <- crossprod(R) / Nn
    ld <- tryCatch(as.numeric(determinant(S, logarithm = TRUE)$modulus),
                   error = function(e) NA_real_)
    if (!is.finite(ld)) next
    k <- p * K * K + K
    pen <- switch(criterion,
                  AIC = 2 * k / Nn,
                  BIC = log(Nn) * k / Nn,
                  HQ  = 2 * log(log(Nn)) * k / Nn)
    crit[p] <- ld + pen
  }
  which.min(crit)
}

#' Time-delay embedding
#'
#' @param x numeric vector.
#' @param m embedding dimension.
#' @param tau delay.
#' @return Matrix with `m` columns; row `t` is `(x[t], x[t+tau], ...)`.
#' @export
#' @examples
#' embed_delay(1:10, m = 3, tau = 2)
embed_delay <- function(x, m = 3L, tau = 1L) {
  x <- as.numeric(x); n <- length(x)
  rows <- n - (m - 1L) * tau
  if (rows < 1L) stop("Series too short for this embedding.", call. = FALSE)
  out <- vapply(seq_len(m), function(j) x[(1 + (j - 1) * tau):(rows + (j - 1) * tau)],
                numeric(rows))
  matrix(out, nrow = rows, ncol = m)
}

#' Binarise a series at its mean (or other threshold)
#'
#' @param x numeric vector.
#' @param threshold cut-off; default is the mean of `x`.
#' @return Integer vector of 0/1.
#' @export
#' @examples
#' binarize(c(1, 5, 2, 8))
binarize <- function(x, threshold = mean(x, na.rm = TRUE)) {
  as.integer(x > threshold)
}

#' Ordinal (permutation) patterns
#'
#' Maps a series onto the index (1..m!) of its ordinal pattern of order `m`.
#'
#' @param x numeric vector.
#' @param m pattern length (2..6).
#' @param tau delay.
#' @return Integer vector of pattern indices.
#' @export
#' @examples
#' ordinal_pattern(c(1, 3, 2, 5, 4, 6), m = 3)
ordinal_pattern <- function(x, m = 3L, tau = 1L) {
  E <- embed_delay(x, m, tau)
  key <- apply(E, 1, function(r) paste(order(r), collapse = ""))
  as.integer(factor(key, levels = sort(unique(key))))
}
