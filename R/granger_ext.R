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

.ty_pair <- function(effect, cause, controls = NULL, lag, d_max = 1L, prior = "bic") {
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
  log_bf <- if (prior == "bic") (N * log(rss_r / rss_u) - lag * log(N)) / 2 else {
    rss_0 <- sum((yd - mean(yd))^2)
    .log_bf_nested(N, ncol(Xr), ncol(Xu), 1 - rss_r / rss_0, 1 - rss_u / rss_0, prior)
  }
  cf <- stats::lm.fit(cbind(1, Xu), yd)$coefficients
  cc <- cf[ncol(Xu) - L + 1L + seq_len(lag)]
  list(statistic = Fst, df1 = lag, df2 = df2, p_value = stats::pf(Fst, lag, df2, lower.tail = FALSE),
       bf10 = exp(min(log_bf, 700)), effect_size = log(rss_r / rss_u), r2_gain = 1 - rss_u / rss_r,
       sign = sign(sum(cc, na.rm = TRUE)), n_eff = N)
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
.dp_stat <- function(effect, cause, controls = NULL, lag = 1L, eps = NULL, C = 8) {
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
  list(statistic = Tn, z = z, p_value = stats::pnorm(z, lower.tail = FALSE), eps = eps, n_eff = N)
}

.dp_pair <- function(effect, cause, controls = NULL, lag = 1L, eps = NULL, C = 8, filter = "none", filter_lag = 1L) {
  if (filter == "var") {
    R <- .var_filter(cbind(effect, cause, controls), filter_lag)
    effect <- R[, 1]; cause <- R[, 2]; controls <- if (ncol(R) > 2) R[, -(1:2), drop = FALSE] else NULL
  }
  .dp_stat(effect, cause, controls, lag, eps, C)
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
#' @param x,y numeric vectors (both directions are tested).
#' @param lag length of the two pasts (`lX = lY`; the theory and the
#'   bandwidth rule are for `lag = 1`, larger lags give a conservative test).
#' @param eps bandwidth in standard deviations; `NULL` = `min(1.5, C *
#'   n^(-2/7))`, the rule of Diks and Panchenko (2006, Sect. 4).
#' @param C constant of the bandwidth rule (about 8 for ARCH-type dynamics).
#' @param filter `"none"` or `"var"`: apply the test to the residuals of a
#'   VAR(`filter_lag`) of the two series, so that only the *nonlinear* part
#'   of the dependence is tested (the usual practice after a linear Granger
#'   test).
#' @param filter_lag lag order of the filtering VAR.
#' @return An object of class `diks_panchenko`: a data frame `tests` (one
#'   row per direction: statistic `T_n`, z-score, one-sided p-value,
#'   bandwidth, effective sample size) and the settings. Method: `print`.
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
diks_panchenko <- function(x, y, lag = 1L, eps = NULL, C = 8, filter = c("none", "var"), filter_lag = 1L) {
  filter <- match.arg(filter)
  x <- .prep_na(x); y <- .prep_na(y)
  if (length(x) != length(y)) stop("x and y must have the same length.", call. = FALSE)
  if (length(x) < 30L) stop("At least 30 observations are needed.", call. = FALSE)
  a <- .dp_pair(y, x, NULL, lag, eps, C, filter, filter_lag)
  b <- .dp_pair(x, y, NULL, lag, eps, C, filter, filter_lag)
  tests <- data.frame(direction = c("x -> y", "y -> x"), statistic = c(a$statistic, b$statistic), z = c(a$z, b$z),
                      p_value = c(a$p_value, b$p_value), eps = c(a$eps, b$eps), n_eff = c(a$n_eff, b$n_eff),
                      stringsAsFactors = FALSE)
  structure(list(tests = tests, lag = as.integer(lag), filter = filter, filter_lag = filter_lag), class = "diks_panchenko")
}

#' @export
print.diks_panchenko <- function(x, digits = 3, ...) {
  cat(sprintf("Diks-Panchenko nonparametric Granger causality test (lag %d, bandwidth %.3g sd%s)\n", x$lag, x$tests$eps[1],
              if (x$filter == "var") sprintf(", VAR(%d) residuals", x$filter_lag) else ""))
  t <- x$tests
  print(data.frame(direction = t$direction, T.n = formatC(t$statistic, digits = digits, format = "g"),
                   z = formatC(t$z, digits = digits, format = "f"), p = .fmt_p(t$p_value), n = t$n_eff), row.names = FALSE,
        right = FALSE)
  cat("  one-sided test of H0: the cause's past adds no information about the effect's present\n")
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
  cat(sprintf("Granger causality in the frequency domain (bivariate VAR(%d), n = %d)\n", x$lag, x$n))
  td <- x$time_domain
  sp <- x$spectrum
  mean_ge <- vapply(td$direction, function(d) mean(sp$geweke[sp$direction == d]), 0)
  print(data.frame(direction = td$direction, geweke.time = round(td$geweke, 4), geweke.mean.spectrum = round(mean_ge, 4),
                   F = round(td$F, 2), p = .fmt_p(td$p_value)), row.names = FALSE, right = FALSE)
  cat(sprintf("  frequency bands with causality (Breitung-Candelon, pointwise level %g), in cycles per time unit:\n", x$alpha))
  if (is.null(x$bands)) cat("    none\n") else
    for (i in seq_len(nrow(x$bands))) cat(sprintf("    %s: %.3f - %.3f (periods %s - %s)\n", x$bands$direction[i], x$bands$from[i],
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

.ty_pair_fun <- function(fixed_lag, max_lag, lag_criterion, d_max, prior) {
  function(effect, cause, ctrl, cores = 1L) {
    inv <- cbind(effect, cause, ctrl); n <- nrow(inv)
    cap <- .cap_lag(n, ncol(inv), d_max)
    p <- if (is.null(fixed_lag)) select_lag(inv, max_lag, lag_criterion) else fixed_lag
    p <- max(1L, min(p, cap))
    g <- .ty_pair(effect, cause, ctrl, p, d_max, prior)
    bl <- do.call(rbind, lapply(seq_len(max(1L, min(as.integer(max_lag), cap))), function(l) {
      gl <- .ty_pair(effect, cause, ctrl, l, d_max, "bic")
      data.frame(lag = l, statistic = gl$statistic, p = gl$p_value)
    }))
    list(row = c(list(lag = p), g[c("statistic", "df1", "df2", "p_value", "bf10", "effect_size", "r2_gain", "sign", "n_eff")]),
         by_lag = bl)
  }
}

.dp_pair_fun <- function(fixed_lag, max_lag, lag_criterion, eps, C, filter) {
  function(effect, cause, ctrl, cores = 1L) {
    L <- if (is.null(fixed_lag)) 1L else max(1L, fixed_lag)
    fl <- if (filter == "var") max(1L, select_lag(cbind(effect, cause, ctrl), max_lag, lag_criterion)) else 1L
    g <- .dp_pair(effect, cause, ctrl, L, eps, C, filter, fl)
    lags <- if (is.null(fixed_lag)) seq_len(max(1L, min(3L, as.integer(max_lag)))) else L
    bl <- do.call(rbind, lapply(lags, function(l) {
      gl <- if (l == L) g else .dp_pair(effect, cause, ctrl, l, eps, C, filter, fl)
      data.frame(lag = l, statistic = gl$z, p = gl$p_value)
    }))
    list(row = list(lag = L, statistic = g$statistic, p_value = g$p_value, effect_size = g$z, n_eff = g$n_eff),
         by_lag = bl, null = list(statistic = g$z, null = list(mean = 0, sd = 1), alternative = "greater"),
         extra = data.frame(z = g$z, eps = g$eps, filter = filter))
  }
}
