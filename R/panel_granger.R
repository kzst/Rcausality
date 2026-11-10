# Granger causality in panels: the estimators and their null distributions -----
#
# Unit level: the F-test of the cause's lags in every unit.
# Pooled level: the test of Juodis, Karavias and Sarafidis (2021).  Under the
#   null hypothesis the coefficients of the cause's lags are zero in every
#   unit, hence equal across the units, and can be pooled although the
#   dynamics differ: the intercept, the effect's own lags and the controls'
#   lags get unit-specific coefficients (they are projected out unit by
#   unit), only the cause's lags are pooled.  The pooled estimator has a bias
#   of order 1/T (Nickell 1981), which makes the test of the reverse
#   direction of a true link over-reject in panels with many units; the
#   half-panel jackknife (Dhaene and Jochmans 2015) removes it:
#     beta = 2 * beta(all periods) - (beta(first half) + beta(second half)) / 2.
# Null distributions: in short series the F distribution is only an
#   approximation (the regressors are lagged values of persistent series), so
#   the tests are referred to a bootstrap under the null hypothesis, as
#   Dumitrescu and Hurlin (2012) recommend for their panel test: in every unit
#   the effect is regenerated from its own estimated autoregression (with the
#   controls) with resampled residuals, the cause is kept as it is.  The
#   bootstrap effect has the dynamics of the observed one (with the
#   small-sample bias of the estimated persistence corrected, Kilian 1998)
#   and no link from the cause.  (Shifting the cause in time is not an alternative in short
#   panels: within a short, demeaned series, values far apart are negatively
#   correlated, about -(1 + rho) / ((1 - rho) T), so that a shifted cause still
#   carries a true link and the surrogate distribution is contaminated.)
# One pass computes, for the observed panel and `B` bootstrap panels, the
# per-unit cross-products from which the unit F statistics and the pooled
# statistic follow.

# bootstrap effect series of one unit under H0: n x B matrix; the first p0
# values are the observed ones (conditioning on the initial values)
.boot_effect <- function(e, Zfix, p0, B) {
  n <- length(e); m <- n - p0; kf <- ncol(Zfix)
  f <- stats::lm.fit(cbind(Zfix, .lagmat(e, p0)), e[(p0 + 1L):n]); cf <- f$coefficients; cf[is.na(cf)] <- 0
  r <- f$residuals - mean(f$residuals); r <- r * sqrt(m / max(1, m - f$rank))
  fx <- as.numeric(Zfix %*% cf[seq_len(kf)]); phi <- cf[kf + seq_len(p0)]
  # Least squares underestimates the persistence of a short series, E(rho_hat) - rho ~ -(1 + 3 rho) / T (Kendall
  # 1954), and a bootstrap process that is less persistent than the data understates the null variability of the
  # tests (Kilian 1998): the sum of the autoregressive coefficients is corrected before the series are generated
  rho <- sum(phi)
  if (m > 8L && rho > 0.05) phi <- phi * min(0.98, rho + (1 + 3 * rho) / m) / rho
  sp <- sum(abs(phi)); if (sp >= 0.99) phi <- phi * 0.98 / sp          # keep the bootstrap process stable
  out <- matrix(0, n, B); out[seq_len(p0), ] <- e[seq_len(p0)]
  R <- matrix(r[sample.int(m, m * B, replace = TRUE)], m, B)
  for (t in seq_len(m)) {
    v <- fx[t] + R[t, ]
    for (l in seq_len(p0)) v <- v + phi[l] * out[p0 + t - l, ]
    out[p0 + t, ] <- v
  }
  out
}

# cross-products of one unit for the lag order l and every column of E (the
# observed effect and its bootstrap copies): after the unit-specific
# regressors Z = [1, controls' lags, augmenting lags of the cause, own lags]
# are projected out of the cause's lags X and of the response y,
#   A = X'X, b = X'y, yy = y'y, df = residual degrees of freedom,
# for all periods and for the two halves of the periods
.reg_unit_contrib <- function(E, x, Zc, l, aug = 0L) {
  n <- nrow(E); nb <- ncol(E); LL <- l + aug; m <- n - LL
  if (m < 4L) return(NULL)
  Xc <- .lagmat(x, LL); X <- Xc[, seq_len(l), drop = FALSE]
  Zf <- matrix(1, m, 1L)
  if (!is.null(Zc)) for (j in seq_len(ncol(Zc))) Zf <- cbind(Zf, .lagmat(Zc[, j], LL))
  if (aug > 0L) Zf <- cbind(Zf, Xc[, (l + 1L):LL, drop = FALSE])
  kz <- ncol(Zf) + LL
  if (m < kz + l + 2L) return(NULL)
  half <- floor(m / 2) >= kz + l + 1L
  sets <- list(full = seq_len(m)); if (half) { sets$h1 <- seq_len(floor(m / 2)); sets$h2 <- (floor(m / 2) + 1L):m }
  if (LL == 1L) {                                        # one own lag: all columns at once (Frisch-Waugh)
    G0 <- E[seq_len(m), , drop = FALSE]; Y0 <- E[2:n, , drop = FALSE]
    res <- lapply(sets, function(rows) {
      q <- qr(Zf[rows, , drop = FALSE]); X1 <- qr.resid(q, X[rows, , drop = FALSE])
      G <- qr.resid(q, G0[rows, , drop = FALSE]); Y1 <- qr.resid(q, Y0[rows, , drop = FALSE])
      gg <- pmax(colSums(G^2), 1e-12); gx <- as.numeric(crossprod(G, X1)); gy <- colSums(G * Y1)
      list(A = array(sum(X1^2) - gx^2 / gg, c(1L, 1L, nb)), b = matrix(as.numeric(crossprod(X1, Y1)) - gx * gy / gg, 1L, nb),
           yy = colSums(Y1^2) - gy^2 / gg, df = length(rows) - q$rank - 1L)
    })
  } else {
    res <- lapply(sets, function(rows) list(A = array(0, c(l, l, nb)), b = matrix(0, l, nb), yy = numeric(nb), df = 0L))
    for (b in seq_len(nb)) {
      Z <- cbind(Zf, .lagmat(E[, b], LL)); y <- E[(LL + 1L):n, b]
      for (nm in names(sets)) {
        rows <- sets[[nm]]; q <- qr(Z[rows, , drop = FALSE])
        Xt <- qr.resid(q, X[rows, , drop = FALSE]); yt <- qr.resid(q, y[rows])
        res[[nm]]$A[, , b] <- crossprod(Xt); res[[nm]]$b[, b] <- crossprod(Xt, yt); res[[nm]]$yy[b] <- sum(yt^2)
        res[[nm]]$df <- length(rows) - q$rank
      }
    }
  }
  y1 <- E[(LL + 1L):n, 1L]
  res$half <- half; res$m <- m; res$tss <- sum((y1 - mean(y1))^2); res$kr <- kz - 1L
  res
}

# The engine: for one link and the lag orders `lags`, the unit F statistics
# and the pooled test of the observed panel (column 1) and of B bootstrap
# panels (columns 2..B+1).
.reg_engine <- function(units, ef, cs, ctrl = character(0), lags = 1L, aug = 0L, B = 99L) {
  N <- length(units); p0 <- max(lags) + aug; nb <- B + 1L
  con <- vector("list", N)
  for (i in seq_len(N)) {
    u <- units[[i]]; n <- nrow(u); if (n < p0 + 5L) next
    e <- .std(u[[ef]]); x <- .std(u[[cs]])
    Zc <- if (length(ctrl)) apply(as.matrix(u[, ctrl, drop = FALSE]), 2, .std) else NULL
    Zf0 <- matrix(1, n - p0, 1L); if (!is.null(Zc)) for (j in seq_len(ncol(Zc))) Zf0 <- cbind(Zf0, .lagmat(Zc[, j], p0))
    E <- if (B > 0L) cbind(e, .boot_effect(e, Zf0, p0, B)) else matrix(e, ncol = 1L)
    con[[i]] <- lapply(lags, function(l) .reg_unit_contrib(E, x, Zc, l, aug))
  }
  quad <- function(A, b) {                               # b' A^-1 b for every column
    if (dim(A)[1L] == 1L) return(b[1L, ]^2 / pmax(A[1L, 1L, ], 1e-12))
    vapply(seq_len(ncol(b)), function(k) tryCatch(sum(b[, k] * solve(A[, , k], b[, k])), error = function(e) NA_real_), 0)
  }
  out <- lapply(seq_along(lags), function(k) {
    l <- lags[k]; cc <- lapply(con, function(ci) if (is.null(ci)) NULL else ci[[k]]); ok <- !vapply(cc, is.null, TRUE)
    # ---- units
    unitF <- matrix(NA_real_, N, nb); udf2 <- rep(NA_real_, N)
    uinfo <- data.frame(m = rep(NA_integer_, N), kr = NA_integer_, tss = NA_real_, rss_r = NA_real_, rss_u = NA_real_, sign = NA_real_)
    for (i in which(ok)) {
      r <- cc[[i]]$full; d2 <- r$df - l
      if (d2 >= 1L) {
        ssr <- quad(r$A, r$b); unitF[i, ] <- (ssr / l) / (pmax(r$yy - ssr, 1e-12) / d2); udf2[i] <- d2
        b1 <- tryCatch(as.numeric(solve(matrix(r$A[, , 1L], l, l), r$b[, 1L])), error = function(e) NA_real_)
        uinfo[i, ] <- list(cc[[i]]$m, cc[[i]]$kr, cc[[i]]$tss, r$yy[1L], max(r$yy[1L] - ssr[1L], 1e-12), sign(sum(b1)))
      }
    }
    # ---- pooled
    halfok <- vapply(cc, function(r) !is.null(r) && isTRUE(r$half), TRUE)
    if (sum(ok) < 2L) return(list(lag = l, obs = NULL, p = rep(NA_real_, nb), unitF = unitF, unit_df2 = udf2, unit = uinfo))
    hpj <- sum(halfok) >= 2L && sum(halfok) >= sum(ok) / 2; use <- which(if (hpj) halfok else ok)
    add <- function(set, what) Reduce(`+`, lapply(cc[use], function(r) r[[set]][[what]]))
    Af <- add("full", "A"); bf <- add("full", "b"); df2 <- max(1, sum(vapply(cc[use], function(r) r$full$df, 0)) - l)
    if (hpj) { A1 <- add("h1", "A"); b1 <- add("h1", "b"); A2 <- add("h2", "A"); b2 <- add("h2", "b") }
    sol <- function(A, b, j) tryCatch(as.numeric(solve(matrix(A[, , j], l, l), b[, j])), error = function(e) rep(NA_real_, l))
    p <- rep(NA_real_, nb); obs <- NULL
    for (j in seq_len(nb)) {
      b0 <- sol(Af, bf, j); beta <- if (hpj) 2 * b0 - (sol(A1, b1, j) + sol(A2, b2, j)) / 2 else b0
      if (anyNA(beta)) next
      V <- matrix(0, l, l); rss_u <- 0; rss_r <- 0
      for (r in cc[use]) {
        Ai <- matrix(r$full$A[, , j], l, l); bi <- r$full$b[, j]
        rss <- r$full$yy[j] - 2 * sum(beta * bi) + sum(beta * (Ai %*% beta))
        V <- V + max(rss, 0) / max(r$full$df, 1) * Ai
        if (j == 1L) { rss_r <- rss_r + r$full$yy[j]; rss_u <- rss_u + r$full$yy[j] - 2 * sum(b0 * bi) + sum(b0 * (Ai %*% b0)) }
      }
      Afj <- matrix(Af[, , j], l, l)
      C <- tryCatch({ Ai <- solve(Afj); Ai %*% V %*% Ai }, error = function(e) NULL); if (is.null(C)) next
      W <- tryCatch(sum(beta * solve(C, beta)), error = function(e) NA_real_); if (!is.finite(W)) next
      p[j] <- stats::pf(W / l, l, df2, lower.tail = FALSE)
      if (j == 1L) obs <- list(beta = beta, cov = C, W = W, statistic = W / l, df1 = l, df2 = df2, p_value = p[j],
                               effect_size = log(rss_r / max(rss_u, 1e-12)), r2_gain = 1 - rss_u / rss_r, sign = sign(sum(beta)),
                               n_eff = sum(vapply(cc[use], function(r) r$m, 0L)), n_units = length(use), hpj = hpj)
    }
    list(lag = l, obs = obs, p = p, unitF = unitF, unit_df2 = udf2, unit = uinfo)
  })
  names(out) <- as.character(lags)
  out
}

# ---- empirical null (Efron 2004) ----------------------------------------------------
# z-scores of p-values; location and scale of the null z-scores (median and
# upper half-spread, robust to the mass at p = 1 of Bonferroni-corrected
# p-values), never below 0 and 1: the calibration can only make a test more
# conservative than its nominal distribution
.zc <- function(p) pmin(pmax(stats::qnorm(p, lower.tail = FALSE), -6), 8)
.empirical_null <- function(p0, min_n = 40L) {
  p0 <- as.numeric(p0); p0 <- p0[is.finite(p0)]
  if (length(p0) < min_n) return(c(mu0 = 0, sigma0 = 1, n = length(p0)))
  z0 <- .zc(p0)
  md <- stats::median(z0)
  c(mu0 = max(0, md), sigma0 = max(1, stats::quantile(z0, stats::pnorm(1), names = FALSE) - md), n = length(z0))
}
.p_calibrate <- function(p, en) {
  if (en[["mu0"]] <= 0 && en[["sigma0"]] <= 1) return(p)
  pmax(p, stats::pnorm((.zc(p) - en[["mu0"]]) / en[["sigma0"]], lower.tail = FALSE))
}

# Bayes factor of "the pooled coefficients of the cause are zero" from the
# (bias-corrected) estimate and its covariance, i.e. with the normal
# approximation of the likelihood beta_hat ~ N(beta, C).  The series are
# standardised within the units, so the coefficients are in standard-deviation
# units.  "bic": unit-information approximation; "g" / "hyper_g": Zellner's
# g-prior and the hyper-g prior through the R-squared implied by the Wald
# statistic; "normal" / "spike_slab": beta ~ N(0, scale^2 I) under the
# alternative; "minnesota": beta_l ~ N(0, (lambda * cross / l^decay)^2), the
# tightness lambda integrated over its hyperprior.
.wald_bf <- function(h, ctl) {
  L <- h$df1; W <- h$W; prior <- ctl$prior
  lbf <- if (prior == "bic") (W - L * log(h$n_eff)) / 2
  else if (prior %in% c("g", "hyper_g")) .log_bf_g0(h$df2 + L + 1, L, W / (W + h$df2), prior, a = ctl$hyper_a)
  else {
    ln <- function(D) { S <- h$cov + D; -0.5 * (as.numeric(determinant(S, logarithm = TRUE)$modulus) + sum(h$beta * solve(S, h$beta))) }
    l0 <- ln(matrix(0, L, L))
    if (prior == "minnesota") {
      base <- (ctl$cross * seq_len(L)^(-ctl$decay))^2
      if (identical(ctl$lambda, "auto")) {
        g <- .gamma_from_mode_sd(ctl$lambda_hyper[["mode"]], ctl$lambda_hyper[["sd"]])
        lam <- exp(seq(log(0.005), log(10), length.out = 80L))
        lw <- stats::dgamma(lam, shape = g[["shape"]], scale = g[["scale"]], log = TRUE) + log(lam); lw <- lw - .logsumexp(lw)
        .logsumexp(lw + vapply(lam, function(l) ln(diag(l^2 * base, L)), 0)) - l0
      } else ln(diag(as.numeric(ctl$lambda)^2 * base, L)) - l0
    } else ln(diag(ctl$scale^2, L)) - l0
  }
  exp(max(min(lbf, 700), -700))
}

# The pooled result of one link from the engine output: the test at the chosen
# order (or the first significant order of the scan), its p-value calibrated
# with the empirical null of the bootstrap panels, the Bayes factor of the
# control's prior (lowered by the change of the likelihood ratio that the
# calibration implies)
.pooled_reg_result <- function(eng, A, lag_rule) {
  ctl <- A$control; lags <- vapply(eng, `[[`, 0, "lag"); have <- !vapply(eng, function(e) is.null(e$obs), TRUE)
  if (!any(have)) return(NULL)
  bonf <- function(p) { p <- p[is.finite(p)]; if (!length(p)) NA_real_ else min(1, min(p) * length(p)) }
  bfs <- vapply(eng, function(e) if (is.null(e$obs)) NA_real_ else .wald_bf(e$obs, ctl), 0)
  P <- do.call(cbind, lapply(eng, `[[`, "p"))            # (B + 1) x lags
  if (lag_rule == "first") {
    fs <- .first_lag(lags, P[1L, ], bfs, A$alpha, A$bf_threshold, by_bf = isTRUE(A$bayes))
    k <- match(fs$lag, lags); p_nom <- fs$p_value; bf <- fs$bf10
    p0 <- if (nrow(P) > 1L) apply(P[-1L, , drop = FALSE], 1L, bonf) else numeric(0)
    bl <- data.frame(lag = lags, statistic = vapply(eng, function(e) if (is.null(e$obs)) NA_real_ else e$obs$statistic, 0), p = P[1L, ])
  } else { k <- which(have)[1L]; p_nom <- P[1L, k]; bf <- bfs[k]; p0 <- P[-1L, k]; bl <- NULL }
  g <- eng[[k]]$obs; if (is.null(g)) return(NULL)
  en <- .empirical_null(p0); p_cal <- .p_calibrate(p_nom, en)
  dl <- (stats::qchisq(p_nom, g$df1, lower.tail = FALSE) - stats::qchisq(p_cal, g$df1, lower.tail = FALSE)) / 2
  if (is.finite(dl) && dl > 0) bf <- bf * exp(-dl)
  detail <- sprintf("pooled cause lags, unit-specific dynamics%s (%d units), F(%d, %d); bootstrap null: %s",
                    if (g$hpj) ", half-panel jackknife" else ", no bias correction (units too short)", g$n_units, as.integer(g$df1),
                    as.integer(g$df2), if (p_cal > p_nom) sprintf("location %.2f, scale %.2f (F-distribution p = %s)", en[["mu0"]], en[["sigma0"]],
                                                                  format.pval(p_nom, digits = 2, eps = 1e-4)) else "no correction")
  list(row = list(lag = lags[k], statistic = g$statistic, df1 = g$df1, df2 = g$df2, p_value = p_cal, bf10 = bf,
                  effect_size = g$effect_size, r2_gain = g$r2_gain, sign = g$sign, n_eff = g$n_eff),
       by_lag = bl, detail = detail, gate = TRUE)
}

# The bootstrap null of the unit statistic of one link: per unit a list
# (statistic, null, alternative, p0), with the F statistic at the order of
# the test, or -- for the scan over the orders -- minus the logarithm of the
# Bonferroni-corrected smallest p-value
.unit_reg_nulls <- function(eng, lag_rule, unit_lag = NULL) {
  lags <- vapply(eng, `[[`, 0, "lag"); N <- nrow(eng[[1L]]$unitF)
  lapply(seq_len(N), function(i) {
    if (lag_rule == "first") {
      Pm <- vapply(seq_along(lags), function(k) stats::pf(eng[[k]]$unitF[i, ], lags[k], eng[[k]]$unit_df2[i], lower.tail = FALSE), numeric(ncol(eng[[1L]]$unitF)))
      Pm <- matrix(Pm, ncol = length(lags))
      pb <- apply(Pm, 1L, function(p) { p <- p[is.finite(p)]; if (!length(p)) NA_real_ else min(1, min(p) * length(p)) })
      if (!is.finite(pb[1L])) return(NULL)
      list(statistic = -log(max(pb[1L], 1e-300)), null = -log(pmax(pb[-1L], 1e-300)), alternative = "greater", p0 = pb[-1L])
    } else {
      k <- if (is.null(unit_lag)) 1L else match(unit_lag[i], lags)
      if (is.na(k) || !is.finite(eng[[k]]$unitF[i, 1L])) return(NULL)
      Fv <- eng[[k]]$unitF[i, ]
      list(statistic = Fv[1L], null = Fv[-1L], alternative = "greater", p0 = stats::pf(Fv[-1L], lags[k], eng[[k]]$unit_df2[i], lower.tail = FALSE))
    }
  })
}

# The unit-level tests of one link from the engine output (Granger family with
# a closed-form Bayes factor): one row per unit, as rcausality() would give
# for the unit alone with the same lag order
.unit_reg_rows <- function(eng, rule, A, prior, bayes) {
  lags <- vapply(eng, `[[`, 0, "lag"); N <- nrow(eng[[1L]]$unitF)
  per <- lapply(seq_along(lags), function(k) {
    e <- eng[[k]]; l <- lags[k]; Fv <- e$unitF[, 1L]; d2 <- e$unit_df2; u <- e$unit
    lbf <- if (prior %in% c("g", "hyper_g")) vapply(seq_len(N), function(i) if (!is.finite(Fv[i])) NA_real_ else
             .log_bf_nested(u$m[i], u$kr[i], u$kr[i] + l, 1 - u$rss_r[i] / u$tss[i], 1 - u$rss_u[i] / u$tss[i], prior), 0)
           else (u$m * log(u$rss_r / u$rss_u) - l * log(u$m)) / 2
    data.frame(lag = l, statistic = as.numeric(Fv), df1 = l, df2 = as.numeric(d2), p_value = as.numeric(stats::pf(Fv, l, d2, lower.tail = FALSE)),
               bf10 = as.numeric(exp(pmin(lbf, 700))), effect_size = as.numeric(log(u$rss_r / u$rss_u)), sign = as.numeric(u$sign),
               row.names = NULL)
  })
  if (rule != "first") return(per[[1L]])
  do.call(rbind, lapply(seq_len(N), function(i) {
    p <- vapply(per, function(d) d$p_value[i], 0); bf <- vapply(per, function(d) d$bf10[i], 0)
    if (all(!is.finite(p))) return(per[[1L]][i, ])
    fs <- .first_lag(lags, p, bf, A$alpha, A$bf_threshold, by_bf = bayes)
    r <- per[[match(fs$lag, lags)]][i, ]; r$p_value <- fs$p_value; r$bf10 <- fs$bf10; r
  }))
}
