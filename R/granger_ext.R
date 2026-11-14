# Generalisations of the linear Granger test ----------------------------------
#  * Toda-Yamamoto (1995): lag-augmented VAR, valid for integrated and
#    cointegrated series without pre-testing (the d_max extra lags are not
#    tested);
#  * Diks-Panchenko (2006): nonparametric test of the conditional
#    independence  y[t] _||_ x-past | y-past  with kernel density
#    estimates (nonlinear Granger causality);
#  * Breitung-Candelon (2006) and Geweke (1982): Granger causality at a
#    frequency (spectral decomposition of the predictive gain).

# ---- Toda-Yamamoto -------------------------------------------------------------

.ty_pair <- function(effect, cause, controls = NULL, lag, d_max = 1L, prior = "bic", prior_par = NULL, vcov = "classical") {
  L <- lag + d_max; n <- length(effect); N <- n - L
  yd <- effect[(L + 1):n]
  own <- .lagmat(effect, L)
  cz <- if (is.null(controls)) NULL else do.call(cbind, lapply(seq_len(ncol(controls)), function(j) .lagmat(controls[, j], L)))
  xc <- .lagmat(cause, L)
  Xr <- cbind(own, cz, xc[, setdiff(seq_len(L), seq_len(lag)), drop = FALSE])   # the d_max augmenting lags stay
  Xu <- cbind(own, cz, xc)
  rss_r <- .rss(yd, Xr); rss_u <- .rss(yd, Xu)
  df2 <- N - ncol(Xu) - 1L
  if (df2 < 1 || rss_u <= 0)
    return(list(statistic = NA_real_, df1 = lag, df2 = df2, p_value = NA_real_, bf10 = NA_real_, effect_size = NA_real_,
                r2_gain = NA_real_, sign = NA_real_, n_eff = N))
  Fst <- ((rss_r - rss_u) / lag) / (rss_u / df2)
  hc <- .hc_wald(yd, Xu, ncol(Xu) - L + seq_len(lag), hac = identical(vcov, "hac"))
  het_p <- if (is.null(hc)) NA_real_ else .het_test(hc$residuals, Xu)
  if (vcov %in% c("hc3", "hac") && !is.null(hc) && is.finite(hc$statistic)) Fst <- hc$statistic
  log_bf <- if (prior == "bic") (N * log(rss_r / rss_u) - lag * log(N)) / 2 else {
    rss_0 <- sum((yd - mean(yd))^2)
    .log_bf_nested(N, ncol(Xr), ncol(Xu), 1 - rss_r / rss_0, 1 - rss_u / rss_0, prior, prior_par)
  }
  cf <- stats::lm.fit(cbind(1, Xu), yd)$coefficients
  cc <- cf[ncol(Xu) - L + 1L + seq_len(lag)]
  list(statistic = Fst, df1 = lag, df2 = df2, p_value = stats::pf(Fst, lag, df2, lower.tail = FALSE),
       bf10 = exp(min(log_bf, 700)), effect_size = log(rss_r / rss_u), r2_gain = 1 - rss_u / rss_r,
       sign = sign(sum(cc, na.rm = TRUE)), n_eff = N, het_p = het_p)
}

# ---- Diks-Panchenko ------------------------------------------------------------

# residuals of the VAR(p) of the involved series (linear part removed)
.var_filter <- function(M, p) {
  M <- as.matrix(M); n <- nrow(M)
  Z <- cbind(1, do.call(cbind, lapply(seq_len(ncol(M)), function(j) .lagmat(M[, j], p))))
  stats::lm.fit(Z, M[(p + 1):n, , drop = FALSE])$residuals
}

# the test statistic of Diks and Panchenko (2006, eq. 13) and its asymptotic
# standard error from the projection of the third-order U-statistic kernel;
# the long-run variance of the projection uses Bartlett weights
# Reference distribution.  The asymptotic normality of the statistic is proved for the bivariate first-order case
# (one lag of the cause, one of the effect).  With controls or more lags the bias of the kernel estimates does not
# vanish fast enough (Diks and Wolski 2016), and with the wide bandwidths of short series the conditioning is coarse,
# on a control and on the past of the effect alike.  The reference distribution therefore comes from conditional
# shifts of the cause (.cond_split()) in all cases; the normal reference is an option (see .dp_stat()).

# ---- conditional shifts ----------------------------------------------------------------------------------------------
# A surrogate of the cause for a conditional null hypothesis must depend on the conditioning vector w as the cause
# itself does: the statistics condition coarsely (kernel windows, nearest-neighbour balls, trees), so that a surrogate
# that is a worse proxy of a function of w than the cause is makes the cause look informative.  The lags of the cause
# are therefore written as
#     c = g(w) + s(w) e,
# location g, scale s and a standardized residual e; a surrogate keeps g and s and shifts e circularly (all lags by the
# same offset), which keeps the serial dependence of e.  The surrogate has the conditional distribution of the cause
# when e is independent of w (a location-scale model) and g and s are estimated well.  Both are estimated flexibly:
#   location: least squares on the polynomial basis in w (linear, quadratic with and without products, cubic) with the
#     smallest leave-one-out error, plus the mean of its residuals over the nearest neighbours in w (the point itself
#     excluded) when that mean predicts the residuals out of sample, with the number of neighbours that does it best;
#   scale: the same for the absolute residuals a = |c - g(w)|, starting from a constant (which gives the plain location
#     model), with two changes (.scale_bases(), .scale_fit()).  The candidate bases include the absolute values of w:
#     a scale often depends on the size of a conditioning variable and not on its sign, a V shape that a polynomial
#     of low order follows badly.  And the selected basis is refitted twice by weighted least squares with the
#     weights 1 / s^2: a = s |e| has a standard deviation proportional to s, so that an unweighted fit is imprecise
#     where the scale is small, and a scale that is too small there gives standardized residuals that are too large
#     everywhere they are shifted to.  The scale is bounded from below by a tenth of the mean absolute residual.
#     (With a cause x = z * noise and the control z, series of 100, 300 and 1000 time points, the kernel test
#     rejected a true conditional null hypothesis in 10% to 12% of the series at the 5% level with the unweighted
#     polynomial fit, in 7% to 9% of the same series with the two changes, and in 3% with the true scale.  A
#     multiplicative model, the fit of the logarithm of the absolute residuals, was much worse: 10% to 25%.)
# A dependence of the shape of the cause on w beyond location and scale is not reproduced; no test of conditional
# independence is valid for all null distributions (Shah and Peters 2020).
.poly_bases <- function(W) {
  N <- nrow(W); d <- ncol(W); one <- matrix(1, N, 1L)
  if (d == 0L) return(list(one))
  out <- list(one, cbind(one, W), cbind(one, W, W^2))
  if (d > 1L) { ij <- utils::combn(d, 2L); out[[4L]] <- cbind(out[[3L]], W[, ij[1L, ], drop = FALSE] * W[, ij[2L, ], drop = FALSE]) }
  out[[length(out) + 1L]] <- cbind(out[[length(out)]], W^3)
  out[vapply(out, ncol, 0L) <= max(2L, N %/% 8L)]
}
# the candidate bases of the scale: those of the location and, after the constant, the absolute values of w without
# and with w itself
.scale_bases <- function(W, bases) {
  if (!ncol(W)) return(bases)
  A <- abs(W); add <- list(cbind(1, A), cbind(1, A, W))
  c(bases[1L], add[vapply(add, ncol, 0L) <= max(2L, nrow(W) %/% 8L)], bases[-1L])
}
.scale_floor <- 0.1
# the least-squares fit with the smallest leave-one-out error (PRESS) among the bases from number `from` on
# (.press_best(): the fit, its error and the number of its basis; NULL when no basis has full rank)
.press_best <- function(bases, y, from = 1L) {
  best <- NULL
  for (b in seq_along(bases)) {
    if (b < from) next
    q <- qr(bases[[b]]); if (q$rank < ncol(bases[[b]])) next
    Q <- qr.Q(q); fit <- as.numeric(Q %*% crossprod(Q, y)); h <- pmin(rowSums(Q^2), 1 - 1e-8)
    press <- sum(((y - fit) / (1 - h))^2)
    if (is.null(best) || press < best$press) best <- list(fit = fit, press = press, basis = b)
  }
  best
}
.press_fit <- function(bases, y, from = 1L) {
  best <- .press_best(bases, y, from)
  if (is.null(best)) rep(mean(y), length(y)) else best$fit
}
# the fitted scale: the candidate basis with the smallest leave-one-out error for the absolute residuals `a`, refitted
# twice with the weights 1 / s^2 (s bounded from below by `lo`); a constant scale needs no weights
.scale_fit <- function(sbases, a, lo) {
  best <- .press_best(sbases, a)
  if (is.null(best)) return(rep(mean(a), length(a)))
  s <- best$fit; B <- sbases[[best$basis]]
  if (ncol(B) > 1L) for (i in 1:2) {
    s2 <- tryCatch(stats::lm.wfit(B, a, 1 / pmax(s, lo)^2)$fitted.values, error = function(e) NULL)
    if (is.null(s2) || any(!is.finite(s2))) break
    s <- as.numeric(s2)
  }
  s
}
# the nearest-neighbour mean of r (nn: neighbours without the point itself) with the number of neighbours that
# predicts r best; zero when no mean predicts r better than zero does
.knn_smooth <- function(nn, r, ks) {
  best <- numeric(length(r)); err <- sum(r^2)
  for (k in ks) {
    f <- rowMeans(matrix(r[nn[, seq_len(k), drop = FALSE]], nrow(nn))); e <- sum((r - f)^2)
    if (e < err) { best <- f; err <- e }
  }
  best
}
.cond_split <- function(W, X) {
  X <- as.matrix(X); N <- nrow(X); W <- .std_mat(W)
  bases <- .poly_bases(W); sbases <- .scale_bases(W, bases)
  k1 <- max(5L, as.integer(floor(sqrt(N)))); ks <- unique(pmax(1L, pmin(c(k1, 2L * k1, 4L * k1), (N - 1L) %/% 4L)))
  nn <- if (ncol(W) && N > 8L) knn_idx_cpp(.jit(W), max(ks)) else NULL
  fit <- X; sc <- X
  for (j in seq_len(ncol(X))) {
    f <- .press_fit(bases, X[, j], from = min(2L, length(bases)))
    if (!is.null(nn)) f <- f + .knn_smooth(nn, X[, j] - f, ks)
    a <- abs(X[, j] - f); lo <- .scale_floor * mean(a)
    s <- .scale_fit(sbases, a, lo)
    if (!is.null(nn)) s <- s + .knn_smooth(nn, a - s, ks)
    fit[, j] <- f; sc[, j] <- pmax(s, lo, 1e-12)
  }
  list(fit = fit, scale = sc, res = (X - fit) / sc, r2 = .cond_r2(X, fit))
}
# the share of the variance of the lags of the cause that the conditioning vector explains (the largest over the lags).
# A share near one means that the past of the cause is a function of the conditioning past: the cause is redundant,
# or the series behave as a deterministic system in which the driver is reconstructed from the driven series
# (Takens 1981).  The surrogates then equal the cause, and the conditional test is not informative: see
# .determined_check().
.cond_r2 <- function(X, fit) {
  v <- colSums((X - rep(colMeans(X), each = nrow(X)))^2)
  max(ifelse(v > 0, 1 - colSums((X - fit)^2) / v, NA_real_), na.rm = TRUE)
}
# A cause that the conditioning vector determines (99% of the variance of a lag or more) has no surrogates that differ
# from it by more than the error of the fitted location: what a test would measure is that error, not a property of
# the data.  No test is carried out then (p-value 1, Bayes factor 1) and a warning explains it.
.determined <- function(r2) isTRUE(is.finite(r2) && r2 > 0.99)
# the surrogate of the cause for the rows `rows` of the standardized residuals
.cond_surr <- function(sp, rows) sp$fit + sp$scale * sp$res[rows, , drop = FALSE]

# Reference distribution (`reference`): "conditional", conditional shifts of the cause (see .cond_split()); "normal",
# the asymptotic distribution of Diks and Panchenko (2006).  The normal reference presupposes that the bias of the
# kernel estimates is negligible beside the standard error of the statistic.  At the bandwidths of the rule (1.5
# standard deviations up to about 350 observations) it is not when the past of the cause depends on the conditioning
# past and the effect is persistent, as in the reverse direction of a link whose driver is autocorrelated: the window
# of the kernel does not resolve the past of the effect, and the cause, which is a function of that past, keeps
# information about the present of the effect.  A surrogate that keeps the dependence of the cause on the conditioning
# past has the same bias, which is why the conditional shifts are the default reference (the normal reference is used
# on request and when n_perm = 0).  Rejection rates of true null hypotheses for both references are in the simulation
# study of the accompanying paper.
.dp_stat <- function(effect, cause, controls = NULL, lag = 1L, eps = NULL, C = 8, n_perm = 0L, cores = 1L,
                     reference = "conditional") {
  e <- .std(effect); x <- .std(cause); n <- length(e); L <- max(1L, as.integer(lag))
  X <- .lagmat(x, L); Y <- .lagmat(e, L)
  if (!is.null(controls)) Y <- cbind(Y, do.call(cbind, lapply(seq_len(ncol(controls)), function(j) .lagmat(.std(controls[, j]), L))))
  Z <- e[(L + 1):n]; N <- length(Z)
  if (is.null(eps)) eps <- min(1.5, C * N^(-2 / 7))
  cc <- dp_counts_cpp(X, Y, Z, eps)
  term <- cc$cxyz * cc$cy - cc$cxy * cc$cyz
  Tn <- (N - 1) / (N * (N - 2)) * sum(term)
  r <- (term + cc$s) / 3 - Tn
  K <- max(1L, floor(N^(1 / 4)))
  lrv <- mean(r^2) + 2 * sum(vapply(seq_len(K), function(k) (1 - k / (K + 1)) * sum(r[seq_len(N - k)] * r[(k + 1):N]) / (N - k), 0))
  S <- sqrt(9 * max(lrv, .Machine$double.eps))
  z <- sqrt(N) * Tn / S
  if (n_perm >= 1L && reference != "normal") {
    sp <- .cond_split(Y, X); lo <- max(L + 1L, floor(0.1 * N))
    Ts <- unlist(.plapply(.shift_draws(lo, N - lo, n_perm), function(s) {
      c2 <- dp_counts_cpp(.cond_surr(sp, ((seq_len(N) + s - 1L) %% N) + 1L), Y, Z, eps)
      (N - 1) / (N * (N - 2)) * sum(c2$cxyz * c2$cy - c2$cxy * c2$cyz)
    }, cores = cores))
    sdT <- stats::sd(Ts)
    return(list(statistic = Tn, z = if (is.finite(sdT) && sdT > 0) (Tn - mean(Ts)) / sdT else NA_real_,
                p_value = if (.determined(sp$r2)) 1 else (1 + sum(Ts >= Tn)) / (length(Ts) + 1), eps = eps, n_eff = N, null = Ts,
                reference = "conditional shifts",
                cond_r2 = sp$r2))
  }
  list(statistic = Tn, z = z, p_value = stats::pnorm(z, lower.tail = FALSE), eps = eps, n_eff = N, null = NULL,
       reference = "normal", cond_r2 = NA_real_)
}

.dp_ref_arg <- function(v) {
  v <- tolower(if (is.null(v)) "conditional" else as.character(v)[1L])
  if (!v %in% c("conditional", "normal")) stop("`reference` must be \"conditional\" or \"normal\".", call. = FALSE)
  v
}

.dp_pair <- function(effect, cause, controls = NULL, lag = 1L, eps = NULL, C = 8, filter = "none", filter_lag = 1L,
                     n_perm = 0L, cores = 1L, reference = "conditional") {
  if (filter == "var") {
    R <- .var_filter(cbind(effect, cause, controls), filter_lag)
    effect <- R[, 1]; cause <- R[, 2]; controls <- if (ncol(R) > 2) R[, -(1:2), drop = FALSE] else NULL
  }
  .dp_stat(effect, cause, controls, lag, eps, C, n_perm, cores, reference)
}

#' Diks-Panchenko nonparametric Granger causality test
#'
#' Tests the conditional independence `y[t] _||_ (x[t-1], ..., x[t-lag]) |
#' (y[t-1], ..., y[t-lag])` without a model for the dependence, so that
#' nonlinear Granger causality is detected (Diks and Panchenko, 2006). The
#' statistic is the U-statistic
#' \deqn{T_n(\varepsilon) = \frac{n-1}{n(n-2)} \sum_i \left[\hat f_{XYZ}(W_i)\hat f_Y(W_i) - \hat f_{XY}(W_i)\hat f_{YZ}(W_i)\right]}{T_n = (n-1)/(n(n-2)) sum_i [f_XYZ f_Y - f_XY f_YZ]}
#' of box-kernel density estimates of bandwidth `eps` (`X` = the cause's
#' past, `Y` = the effect's past, `Z` = the effect's present, series
#' standardised); `sqrt(n) T_n / S_n` is asymptotically standard normal, with
#' `S_n` from the first-order projection of the kernel and a Bartlett
#' long-run variance. The test is one-sided. Unlike the Hiemstra-Jones
#' test, it is consistent because the kernel is weighted by the density of
#' the conditioning variable.
#'
#' @section Reference distribution:
#' By default the statistic is referred to conditional shifts of the cause
#' (`reference = "conditional"`): the lags of the cause are split into a
#' location and a scale, both estimated as functions of the conditioning
#' vector, and standardized residuals, which are shifted circularly (see
#' [rcausality()]). The surrogates keep the dependence of the cause's past
#' on the effect's past, which the null hypothesis allows. When the
#' conditioning vector determines the cause (99% of its variance or more,
#' as in noise-free deterministic systems), no test is carried out and the
#' p-value is 1.
#'
#' `reference = "normal"` is the asymptotic test of Diks and Panchenko
#' (2006), which needs no surrogates. Its normal limit is proved for one
#' lag of each series, and it presupposes that the bias of the kernel
#' estimates is negligible beside the standard error of the statistic. At
#' the bandwidths of the rule this fails when the cause's past depends on
#' the effect's past and the effect is autocorrelated, as in the reverse
#' direction of a link whose driver is persistent: a window of 1.5 standard
#' deviations does not resolve the effect's past, the cause is a function of
#' that past and keeps information about the effect's present, and the
#' normal reference rejects a true null hypothesis in a large share of the
#' cases (the simulation study of the accompanying paper gives the rates).
#' Use it for series whose pasts are at most weakly dependent, for example
#' after `filter = "var"` when the dependence is linear, and with one lag.
#'
#' @param x,y numeric vectors (both directions are tested).
#' @param lag length of the two pasts (`lX = lY`; the bandwidth rule is for
#'   `lag = 1`).
#' @param eps bandwidth in standard deviations; `NULL` = `min(1.5, C *
#'   n^(-2/7))`, the rule of Diks and Panchenko (2006, Sect. 4).
#' @param C constant of the bandwidth rule (about 8 for ARCH-type dynamics).
#' @param filter `"none"` or `"var"`: apply the test to the residuals of a
#'   VAR(`filter_lag`) of the two series, so that only the *nonlinear* part
#'   of the dependence is tested (the usual practice after a linear Granger
#'   test).
#' @param filter_lag lag order of the filtering VAR.
#' @param n_perm,seed number of conditional shifts and their seed
#'   (`n_perm = 0` gives the normal reference).
#' @param reference `"conditional"` (default) or `"normal"`: see the section
#'   on the reference distribution. Diks and Wolski (2016) discuss the
#'   failure of the normal reference with more than one lag or with further
#'   conditioning variables.
#' @return An object of class `diks_panchenko`: a data frame `tests` (one
#'   row per direction: statistic `T_n`, z-score (with conditional shifts:
#'   the statistic standardized by mean and standard deviation of its
#'   surrogate values), one-sided p-value, bandwidth, effective sample size,
#'   reference) and the settings. Method: `print`.
#' @references
#' Diks, C. and Panchenko, V. (2006). A new statistic and practical guidelines for nonparametric Granger causality testing. *Journal of Economic Dynamics and Control*, 30(9-10), 1647-1669. \doi{10.1016/j.jedc.2005.08.008}
#'
#' Hiemstra, C. and Jones, J. D. (1994). Testing for linear and nonlinear Granger causality in the stock price-volume relation. *The Journal of Finance*, 49(5), 1639-1664. \doi{10.1111/j.1540-6261.1994.tb04776.x}
#' @export
#' @examples
#' # x drives y through a nonlinear (logistic-map) coupling: the test rejects
#' # x -> y and keeps y -> x
#' d <- sim_logistic(400, "x_to_y", seed = 1)
#' diks_panchenko(d$x, d$y)
#' \donttest{
#' # only the nonlinear part: the VAR(2) residuals are tested
#' diks_panchenko(d$x, d$y, filter = "var", filter_lag = 2)
#' }
diks_panchenko <- function(x, y, lag = 1L, eps = NULL, C = 8, filter = c("none", "var"), filter_lag = 1L,
                           n_perm = 199L, seed = NULL, reference = c("conditional", "normal")) {
  filter <- match.arg(filter); reference <- match.arg(reference)
  .local_seed(seed)
  x <- .prep_na(x); y <- .prep_na(y)
  if (length(x) != length(y)) stop("x and y must have the same length.", call. = FALSE)
  if (length(x) < 30L) stop("At least 30 observations are needed.", call. = FALSE)
  a <- .dp_pair(y, x, NULL, lag, eps, C, filter, filter_lag, n_perm, reference = reference)
  b <- .dp_pair(x, y, NULL, lag, eps, C, filter, filter_lag, n_perm, reference = reference)
  tests <- data.frame(direction = c("x -> y", "y -> x"), statistic = c(a$statistic, b$statistic), z = c(a$z, b$z),
                      p_value = c(a$p_value, b$p_value), eps = c(a$eps, b$eps), n_eff = c(a$n_eff, b$n_eff),
                      reference = c(a$reference, b$reference), stringsAsFactors = FALSE)
  structure(list(tests = tests, lag = as.integer(lag), filter = filter, filter_lag = filter_lag, reference = a$reference),
            class = "diks_panchenko")
}

#' @export
print.diks_panchenko <- function(x, digits = 3, ...) {
  .out(sprintf("Diks-Panchenko nonparametric Granger causality test (lag %d, bandwidth %.3g sd%s)\n", x$lag, x$tests$eps[1],
              if (x$filter == "var") sprintf(", VAR(%d) residuals", x$filter_lag) else ""))
  if (!is.null(x$reference)) .out(sprintf("  reference distribution: %s\n", if (x$reference == "normal") "normal (asymptotic)" else x$reference))
  t <- x$tests
  print(data.frame(direction = t$direction, T.n = formatC(t$statistic, digits = digits, format = "g"),
                   z = formatC(t$z, digits = digits, format = "f"), p = .fmt_p(t$p_value), n = t$n_eff), row.names = FALSE,
        right = FALSE)
  .out("  one-sided test of H0: the cause's past adds no information about the effect's present\n")
  invisible(x)
}

# ---- frequency-domain Granger causality ----------------------------------------

#' Granger causality in the frequency domain
#'
#' Decomposes the predictive gain of a bivariate VAR(p) over the frequencies.
#' Two complementary quantities are returned for both directions:
#' * Geweke's (1982) spectral measure
#'   \deqn{f_{x \to y}(\omega) = \ln \frac{S_{yy}(\omega)}{S_{yy}(\omega) - (\Sigma_{xx} - \Sigma_{xy}^2/\Sigma_{yy})\,|H_{yx}(\omega)|^2},}{f_x->y(w) = log(S_yy(w) / (S_yy(w) - (Sigma_xx - Sigma_xy^2 / Sigma_yy) |H_yx(w)|^2)),}
#'   where `H` is the transfer function and `S` the spectral matrix of the
#'   VAR; its average over `(-pi, pi)` equals the time-domain measure
#'   `log(sigma^2_restricted / sigma^2_full)` for a stationary VAR;
#' * the Breitung-Candelon (2006) test of `H0: no causality at frequency w`,
#'   the linear restrictions `sum_l b_l cos(l w) = sum_l b_l sin(l w) = 0` on
#'   the cause coefficients of the effect equation (an F-test with two
#'   numerator degrees of freedom, one at `w = 0` and `w = pi`; informative
#'   for `p >= 3`).
#'
#' @param x,y numeric vectors.
#' @param lag VAR order or `"auto"` (selected by `lag_criterion` up to
#'   `max_lag`).
#' @param max_lag,lag_criterion lag selection, as in [rcausality()].
#' @param n_freq number of frequencies on `[0, pi]`.
#' @param alpha level of the pointwise tests (for the plot and the
#'   significant bands).
#' @return An object of class `granger_spectral`: `spectrum` (one row per
#'   direction and frequency: `freq` (radians), `cycles` (per time unit),
#'   `period`, `geweke`, `F`, `df1`, `p_value`), `time_domain` (Geweke's
#'   measure and the F-test of each direction), `bands` (frequency bands
#'   where the pointwise test rejects at `alpha`), `lag`, `n`. Methods:
#'   `print`, `plot`.
#' @references
#' Geweke, J. (1982). Measurement of linear dependence and feedback between multiple time series. *Journal of the American Statistical Association*, 77(378), 304-313. \doi{10.1080/01621459.1982.10477803}
#'
#' Breitung, J. and Candelon, B. (2006). Testing for short- and long-run causality: a frequency-domain approach. *Journal of Econometrics*, 132(2), 363-378. \doi{10.1016/j.jeconom.2005.02.004}
#'
#' Ding, M., Chen, Y. and Bressler, S. L. (2006). Granger causality: basic theory and application to neuroscience. In: *Handbook of Time Series Analysis*, Wiley-VCH, 437-460. \doi{10.1002/9783527609970.ch17}
#' @export
#' @examples
#' # y depends on x[t-1] + x[t-2]: the transfer gain 0.8 |cos(w/2)| is
#' # largest in the long run and vanishes at the period of 2 time units, so
#' # the pointwise tests reject at low but not at the highest frequencies;
#' # y does not drive x at any frequency
#' set.seed(2); n <- 600; x <- rnorm(n); y <- numeric(n)
#' for (t in 3:n) y[t] <- 0.3 * y[t - 1] + 0.4 * (x[t - 1] + x[t - 2]) + rnorm(1)
#' gs <- granger_spectral(x, y, lag = 3)
#' gs
#' \donttest{
#' plot(gs)
#' }
granger_spectral <- function(x, y, lag = "auto", max_lag = 8L, lag_criterion = c("BIC", "AIC", "HQ"), n_freq = 256L,
                             alpha = 0.05) {
  lag_criterion <- match.arg(lag_criterion)
  x <- .prep_na(x); y <- .prep_na(y); n <- length(x)
  if (length(y) != n) stop("x and y must have the same length.", call. = FALSE)
  M <- cbind(x = x, y = y)
  p <- if (identical(lag, "auto")) select_lag(M, max_lag, lag_criterion) else as.integer(lag)
  p <- max(1L, p); N <- n - p
  if (N < 4 * p + 5) stop("Series too short for a VAR of this order.", call. = FALSE)
  Z <- cbind(1, .lagmat(x, p), .lagmat(y, p))           # columns: x lags 1..p, then y lags 1..p
  Yd <- M[(p + 1):n, , drop = FALSE]
  fit <- stats::lm.fit(Z, Yd); B <- fit$coefficients; E <- fit$residuals
  k <- ncol(Z); Sigma <- crossprod(E) / (N - k)
  ZtZi <- chol2inv(chol(crossprod(Z)))
  A <- lapply(seq_len(p), function(l) rbind(c(B[1 + l, "x"], B[1 + p + l, "x"]), c(B[1 + l, "y"], B[1 + p + l, "y"])))
  w <- seq(0, pi, length.out = max(16L, as.integer(n_freq)))
  # transfer function H(w) = (I - sum_l A_l e^{-i l w})^{-1}
  spec <- lapply(w, function(om) {
    Aw <- diag(2) + 0i; for (l in seq_len(p)) Aw <- Aw - A[[l]] * exp(-1i * l * om)
    H <- solve(Aw); S <- H %*% Sigma %*% Conj(t(H))
    list(H = H, S = Re(diag(S)))
  })
  # cause / effect columns: index 1 = x, 2 = y
  one_dir <- function(ic, ie) {
    s_part <- Sigma[ic, ic] - Sigma[ic, ie]^2 / Sigma[ie, ie]
    ge <- vapply(spec, function(s) { see <- s$S[ie]; log(see / max(see - s_part * Mod(s$H[ie, ic])^2, 1e-300)) }, 0)
    cols <- if (ic == 1) 1 + seq_len(p) else 1 + p + seq_len(p)       # cause lags in the effect equation
    b <- B[cols, ie]; V <- Sigma[ie, ie] * ZtZi[cols, cols, drop = FALSE]
    bc <- t(vapply(w, function(om) {
      R <- rbind(cos(seq_len(p) * om), sin(seq_len(p) * om))
      if (max(abs(R[2, ])) < 1e-8) R <- R[1, , drop = FALSE]
      if (p == 1 && nrow(R) == 2) R <- R[1, , drop = FALSE]
      Rb <- R %*% b; RV <- R %*% V %*% t(R)
      Fs <- as.numeric(t(Rb) %*% solve(RV, Rb)) / nrow(R)
      c(Fs, nrow(R), stats::pf(Fs, nrow(R), N - k, lower.tail = FALSE))
    }, numeric(3)))
    data.frame(direction = if (ic == 1) "x -> y" else "y -> x", freq = w, cycles = w / (2 * pi),
               period = ifelse(w > 0, 2 * pi / w, Inf), geweke = ge, F = bc[, 1], df1 = bc[, 2], p_value = bc[, 3],
               stringsAsFactors = FALSE)
  }
  sp <- rbind(one_dir(1, 2), one_dir(2, 1))
  td <- do.call(rbind, lapply(list(c(1, 2), c(2, 1)), function(ij) {
    g <- .granger_pair(M[, ij[2]], M[, ij[1]], NULL, p)
    data.frame(direction = if (ij[1] == 1) "x -> y" else "y -> x", geweke = g$effect_size, F = g$statistic,
               df1 = g$df1, df2 = g$df2, p_value = g$p_value, stringsAsFactors = FALSE)
  }))
  bands <- do.call(rbind, lapply(split(sp, sp$direction), function(s) {
    sig <- s$p_value < alpha; if (!any(sig)) return(NULL)
    rl <- rle(sig); end <- cumsum(rl$lengths); start <- end - rl$lengths + 1
    data.frame(direction = s$direction[1], from = s$cycles[start[rl$values]], to = s$cycles[end[rl$values]],
               stringsAsFactors = FALSE)
  }))
  structure(list(spectrum = sp, time_domain = td, bands = bands, lag = p, n = n, alpha = alpha),
            class = "granger_spectral")
}

#' @export
print.granger_spectral <- function(x, digits = 3, ...) {
  .out(sprintf("Granger causality in the frequency domain (bivariate VAR(%d), n = %d)\n", x$lag, x$n))
  td <- x$time_domain
  sp <- x$spectrum
  mean_ge <- vapply(td$direction, function(d) mean(sp$geweke[sp$direction == d]), 0)
  print(data.frame(direction = td$direction, geweke.time = round(td$geweke, 4), geweke.mean.spectrum = round(mean_ge, 4),
                   F = round(td$F, 2), p = .fmt_p(td$p_value)), row.names = FALSE, right = FALSE)
  .out(sprintf("  frequency bands with causality (Breitung-Candelon, pointwise level %g), in cycles per time unit:\n", x$alpha))
  if (is.null(x$bands)) .out("    none\n") else
    for (i in seq_len(nrow(x$bands))) .out(sprintf("    %s: %.3f - %.3f (periods %s - %s)\n", x$bands$direction[i], x$bands$from[i],
                                                  x$bands$to[i], if (x$bands$to[i] > 0) formatC(1 / x$bands$to[i], digits = 3, format = "g") else "Inf",
                                                  if (x$bands$from[i] > 0) formatC(1 / x$bands$from[i], digits = 3, format = "g") else "Inf"))
  invisible(x)
}

#' @rdname granger_spectral
#' @param main optional title.
#' @param ... unused.
#' @export
plot.granger_spectral <- function(x, main = NULL, ...) {
  sp <- x$spectrum; dirs <- unique(sp$direction); cols <- c("#2B5C8A", "#B2403A")
  op <- graphics::par(mfrow = c(2, 1), mar = c(2.2, 4.2, 3, 1), oma = c(2, 0, 0, 0)); on.exit(graphics::par(op))
  yl <- range(0, sp$geweke, na.rm = TRUE)
  graphics::plot(NA, xlim = c(0, 0.5), ylim = yl, xlab = "", ylab = "Geweke measure f(w)",
                 main = if (is.null(main)) sprintf("Spectral Granger causality, VAR(%d)", x$lag) else main, font.main = 2,
                 cex.main = 0.95)
  graphics::abline(h = 0, col = .pal$grid)
  for (i in seq_along(dirs)) {
    s <- sp[sp$direction == dirs[i], ]
    graphics::lines(s$cycles, s$geweke, col = cols[i], lwd = 2)
    graphics::abline(h = x$time_domain$geweke[x$time_domain$direction == dirs[i]], col = cols[i], lty = 3)
  }
  graphics::legend("topright", legend = c(dirs, "time-domain measure (dotted)"), col = c(cols[seq_along(dirs)], "grey40"),
                   lty = c(rep(1, length(dirs)), 3), lwd = c(rep(2, length(dirs)), 1), bty = "n", cex = 0.75)
  lp <- -log10(pmax(sp$p_value, 1e-300))
  graphics::plot(NA, xlim = c(0, 0.5), ylim = c(0, max(lp, -log10(x$alpha) * 1.5, na.rm = TRUE)), xlab = "",
                 ylab = "-log10 p (Breitung-Candelon)", main = "Pointwise test of no causality at frequency w",
                 font.main = 2, cex.main = 0.9)
  if (!is.null(x$bands)) for (i in seq_len(nrow(x$bands))) {
    ci <- cols[match(x$bands$direction[i], dirs)]
    graphics::rect(x$bands$from[i], 0, x$bands$to[i], graphics::par("usr")[4], col = grDevices::adjustcolor(ci, 0.12), border = NA)
  }
  for (i in seq_along(dirs)) {
    s <- sp[sp$direction == dirs[i], ]
    graphics::lines(s$cycles, -log10(pmax(s$p_value, 1e-300)), col = cols[i], lwd = 2)
  }
  graphics::abline(h = -log10(x$alpha), lty = 2, col = "grey40")
  graphics::mtext("frequency (cycles per time unit; period = 1 / frequency)", side = 1, line = 0.6, outer = TRUE, cex = 0.8)
  invisible(x)
}

# ---- pair functions for rcausality() --------------------------------------------

# largest lag order that leaves degrees of freedom for the unrestricted
# regression of the effect on `n_vars` lagged series (plus `extra` lags)
.cap_lag <- function(n, n_vars, extra = 0L) max(1L, floor((n - 2 - extra) / (2 + 2 * n_vars)))

# the "first significant order" rule of the Granger family: the smallest order
# whose test is significant at level alpha / (number of orders scanned)
# (Bonferroni over the scan); the reported p-value is the Bonferroni-corrected
# minimum over the orders, a valid p-value of "no link at any order".  For
# Bayes-factor decisions: the first order with BF10 above the threshold
# (otherwise the order with the largest Bayes factor); the reported Bayes
# factor is the average over the orders, i.e. the Bayes factor of "a link at
# one of the scanned orders" with equal prior weights (the Bayesian
# counterpart of the Bonferroni correction).
.first_lag <- function(lags, p, bf, alpha, bf_threshold, by_bf = FALSE) {
  ok <- is.finite(p); nl <- max(1L, sum(ok))
  p_min <- if (any(ok)) min(1, min(p[ok]) * nl) else NA_real_
  if (by_bf && any(is.finite(bf))) {
    hit <- which(is.finite(bf) & bf > bf_threshold)
    l <- if (length(hit)) lags[hit[1L]] else lags[which.max(ifelse(is.finite(bf), bf, -Inf))]
  } else {
    hit <- which(ok & p < alpha / nl)
    l <- if (length(hit)) lags[hit[1L]] else if (any(ok)) lags[which.min(ifelse(ok, p, Inf))] else lags[1L]
  }
  list(lag = l, p_value = p_min, bf10 = if (any(is.finite(bf))) mean(bf[is.finite(bf)]) else NA_real_, n_orders = nl)
}

.ty_pair_fun <- function(fixed_lag, max_lag, lag_criterion, d_max, prior, first = FALSE, alpha = 0.05, prior_par = NULL,
                         vcov = "classical") {
  function(effect, cause, ctrl, cores = 1L) {
    inv <- cbind(effect, cause, ctrl); n <- nrow(inv)
    cap <- .cap_lag(n, ncol(inv), d_max)
    lags <- seq_len(max(1L, min(as.integer(max_lag), cap)))
    scan <- lapply(lags, function(l) .ty_pair(effect, cause, ctrl, l, d_max, prior, prior_par, vcov))
    bl <- data.frame(lag = lags, statistic = vapply(scan, `[[`, 0, "statistic"), p = vapply(scan, `[[`, 0, "p_value"))
    if (first) {
      fs <- .first_lag(lags, bl$p, vapply(scan, `[[`, 0, "bf10"), alpha, 3)
      p <- fs$lag; g <- scan[[match(p, lags)]]; g$p_value <- fs$p_value; g$bf10 <- fs$bf10
    } else {
      p <- if (is.null(fixed_lag)) select_lag(inv, max_lag, lag_criterion) else fixed_lag
      p <- max(1L, min(p, cap))
      g <- if (p %in% lags) scan[[match(p, lags)]] else .ty_pair(effect, cause, ctrl, p, d_max, prior, prior_par, vcov)
    }
    list(row = c(list(lag = p), g[c("statistic", "df1", "df2", "p_value", "bf10", "effect_size", "r2_gain", "sign", "n_eff")]),
         by_lag = bl, extra = data.frame(het_p = g$het_p, vcov = vcov, stringsAsFactors = FALSE))
  }
}

.dp_pair_fun <- function(fixed_lag, max_lag, lag_criterion, eps, C, filter, n_perm = 200L, reference = "conditional") {
  function(effect, cause, ctrl, cores = 1L) {
    L <- if (is.null(fixed_lag)) 1L else max(1L, fixed_lag)
    fl <- if (filter == "var") max(1L, select_lag(cbind(effect, cause, ctrl), max_lag, lag_criterion)) else 1L
    g <- .dp_pair(effect, cause, ctrl, L, eps, C, filter, fl, n_perm, cores, reference)
    # (No profile over further lags: every lag needs its own surrogate test.  A numeric `lag` tests longer pasts.)
    bl <- data.frame(lag = L, statistic = g$z, p = g$p_value)
    list(row = list(lag = L, statistic = g$statistic, p_value = g$p_value, effect_size = g$z, n_eff = g$n_eff),
         by_lag = bl,
         null = if (is.null(g$null)) list(statistic = g$z, null = list(mean = 0, sd = 1), alternative = "greater")
                else list(statistic = g$statistic, null = g$null, alternative = "greater", cond_r2 = g$cond_r2),
         extra = data.frame(z = g$z, eps = g$eps, filter = filter, reference = g$reference, cond_r2 = g$cond_r2))
  }
}
