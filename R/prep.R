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
#' # the missing value is interpolated, then the series is standardised
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
#' @references
#' Luetkepohl, H. (2005). *New Introduction to Multiple Time Series Analysis*. Springer, Berlin.
#'
#' Schwarz, G. (1978). Estimating the dimension of a model. *The Annals of Statistics*, 6(2), 461-464. \doi{10.1214/aos/1176344136}
#' @export
#' @examples
#' # the VAR order of a simulated VAR(2) selected by BIC (the default)
#' d <- sim_var(300, scenario = "x_to_y", seed = 1)
#' select_lag(d)
#' select_lag(d, criterion = "AIC")     # AIC tends to pick longer lags
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

# ---- the memory of the effect ----------------------------------------------------------------------------------------
# A test of "the cause does not drive the effect" conditions on a finite past of the effect, and its null hypothesis
# is the one with that past.  Without the link it is true only when the past holds what the effect remembers of
# itself: when the effect depends on earlier values of its own than the test conditions on, and the cause depends on
# them too (as in the reverse direction of a link), the cause is informative given the short past although it does
# not drive the effect.  Example: x[t] = 0.4 x[t-1] + 0.4 x[t-2] + noise drives y; the hypothesis "y does not help to
# predict x given one value of each" is false.  The order of the VAR by the BIC is no safeguard: the past of the driven
# series stands in for the missing lag of the driver, the first-order VAR fits nearly as well, and the criterion
# charges four coefficients for the one that matters (it chose order one in 82%, 27% and 0% of such pairs of 100, 300
# and 1000 time points, and the F test of the reverse direction then rejected in 56%, 31% and 8% of all pairs).
# .own_order(): the order of the autoregression of the effect alone, orders 0 to max_lag compared on the same sample
# by the corrected AIC (Hurvich and Tsai 1989).  It does not depend on the cause, so that it selects nothing in favour
# of the tested link; and the AIC rather than the BIC, because an order that is too long costs a little power and one
# that is too short makes the null hypothesis false (a series observed with noise has an autoregression of infinite
# order, for which the BIC stops early).
# Every test keeps its statistic for the lags of the cause and conditions on the past of the effect up to that order:
# the regression tests include that many lags of the effect itself (.granger_pair(): own lags), and the tests with
# conditional surrogates let the surrogates of the cause depend on the longer past (.cond_extra()).  The order of the
# test, the lags of the cause, is chosen as before.  A numeric `lag` of the user gives the regression of that order.
# options(Rcausality.own_memory = FALSE) switches the rule off (order 0: the tests of the versions up to 0.0.11); it
# is not meant for analyses, the simulation study of the accompanying paper uses it to show what the rule changes.
.own_order <- function(e, max_lag = 8L) {
  if (isFALSE(getOption("Rcausality.own_memory"))) return(0L)
  e <- as.numeric(e); n <- length(e)
  if (n < 12L || anyNA(e) || !is.finite(stats::sd(e)) || stats::sd(e) == 0) return(0L)
  pm <- max(1L, min(as.integer(max_lag), (n - 1L) %/% 4L)); m <- n - pm
  Lm <- cbind(1, .lagmat(e, pm)); y <- e[(pm + 1L):n]
  crit <- vapply(0:pm, function(p) {
    rss <- sum(stats::lm.fit(Lm[, seq_len(p + 1L), drop = FALSE], y)$residuals^2); k <- p + 2L
    if (rss <= 0 || m - k - 1L < 1L) Inf else m * log(rss / m) + 2 * k * m / (m - k - 1L)
  }, 0)
  as.integer(which.min(crit) - 1L)
}
# the lags `from`..`to` of v at the time points t_idx; time points before the start of the series take its first value
.hist_pad <- function(v, t_idx, from, to) {
  if (to < from) return(matrix(0, length(t_idx), 0L))
  do.call(cbind, lapply(from:to, function(l) v[pmax(t_idx - l, 1L)]))
}
# the conditioning vector of the conditional surrogates: `W`, the conditioning vector of the statistic, which holds
# `have` lags of the effect, with the further lags of the effect up to its own order `own` (standardized)
.cond_extra <- function(W, effect, t_idx, have, own) {
  if (is.null(own) || !is.finite(own) || own <= have) return(W)
  cbind(W, .hist_pad(.std(effect), t_idx, have + 1L, as.integer(own)))
}

#' Time-delay embedding
#'
#' @param x numeric vector.
#' @param m embedding dimension.
#' @param tau delay.
#' @return Matrix with `m` columns; row `t` is `(x[t], x[t+tau], ...)`.
#' @export
#' @examples
#' # rows are the delay vectors (x[t-4], x[t-2], x[t]) for tau = 2, m = 3
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
#' # 1 above the mean (4), 0 below
#' binarize(c(1, 5, 2, 8))
binarize <- function(x, threshold = mean(x, na.rm = TRUE)) {
  as.integer(x > threshold)
}

#' Ordinal (permutation) patterns
#'
#' Maps a series onto the index (1..m!) of its ordinal pattern of order `m`:
#' the Lehmer code of the rank permutation of every delay vector (ties are
#' broken by position). The code of a pattern does not depend on which
#' patterns occur in the series, so codes are comparable across series
#' (the increasing pattern is always 1, the decreasing one m!).
#'
#' @param x numeric vector.
#' @param m pattern length (2..10).
#' @param tau delay.
#' @return Integer vector of pattern indices.
#' @export
#' @examples
#' # each window of three values is coded by the rank order of its entries
#' # (1 = increasing ... 6 = decreasing), the symbols CCDH compresses
#' ordinal_pattern(c(1, 3, 2, 5, 4, 6), m = 3)
#' ordinal_pattern(c(3, 2, 1), m = 3)        # 6 = the decreasing pattern
ordinal_pattern <- function(x, m = 3L, tau = 1L) {
  m <- as.integer(m)
  if (m < 2L || m > 10L) stop("m must be between 2 and 10.", call. = FALSE)
  as.integer(.lehmer_codes(embed_delay(x, m, tau)))
}

# Lehmer code (1..m!) of the rank permutation of every row, ties broken by
# column position; vectorised over the rows, double precision (exact up to m = 18)
.lehmer_codes <- function(E) {
  E <- as.matrix(E); m <- ncol(E); n <- nrow(E)
  if (m == 1L) return(rep(1, n))
  code <- rep(0, n)
  for (j in seq_len(m - 1L)) {
    smaller_after <- rep(0, n)
    for (l in (j + 1L):m) smaller_after <- smaller_after + (E[, l] < E[, j])
    code <- code + smaller_after * factorial(m - j)
  }
  code + 1
}
