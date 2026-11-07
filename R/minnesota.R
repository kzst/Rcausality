# Bayesian Granger regression with a hierarchical Minnesota prior --------------
# For short series the unrestricted regression of the effect on the lags of
# all involved series has few residual degrees of freedom (or none), and the
# F-test loses its power or cannot be computed.  The Minnesota prior
# (Litterman 1986) shrinks the coefficients towards zero, more strongly for
# distant lags; its overall tightness is integrated out under the hyperprior
# of Giannone, Lenza and Primiceri (2015), so that the Bayes factor of
# "cause lags in" against "cause lags out" is a ratio of two marginal
# likelihoods that both average over the tightness:
#   p(y | M) = int p(y | M, lambda) p(lambda) d lambda,
# with p(y | M, lambda) the closed-form normal-inverse-gamma marginal
# likelihood.  Exact (independent) posterior draws of lambda, sigma^2 and the
# coefficients follow from the same grid.

# Gamma(shape, scale) with the given mode and standard deviation
.gamma_from_mode_sd <- function(mode, sd) {
  a <- 2 * sd^2 + mode^2
  shape <- (a + sqrt(a^2 - 4 * sd^4)) / (2 * sd^2)
  c(shape = shape, scale = mode / (shape - 1))
}

# residual sd of an AR(p) of a series (the Minnesota scale); the sd of the
# series when the autoregression has no residual degrees of freedom
.ar_scale <- function(v, p) {
  n <- length(v); N <- n - p
  s <- if (N - p - 1 >= 2) sqrt(.rss(v[(p + 1):n], .lagmat(v, p)) / (N - p - 1)) else stats::sd(v)
  if (!is.finite(s) || s <= 0) 1 else s
}

.logsumexp <- function(v) { m <- max(v); if (!is.finite(m)) m else m + log(sum(exp(v - m))) }

.bayes_granger_minnesota <- function(effect, cause, controls, lag, control, prior_prob, cores = 1L) {
  ys <- .std(effect); xs <- .std(cause); n <- length(ys)
  cs <- if (is.null(controls)) NULL else apply(as.matrix(controls), 2, .std)
  own <- .lagmat(ys, lag); Xc <- .lagmat(xs, lag); colnames(Xc) <- paste0("lag", seq_len(lag))
  cz <- if (is.null(cs)) NULL else do.call(cbind, lapply(seq_len(ncol(cs)), function(j) .lagmat(cs[, j], lag)))
  yd <- ys[(lag + 1):n] - control$own_mean * own[, 1]          # prior mean of the first own lag moved into y
  l_dec <- seq_len(lag)^(-2 * control$decay)
  ds <- c(l_dec / .ar_scale(ys, lag)^2,
          if (!is.null(cs)) unlist(lapply(seq_len(ncol(cs)), function(j) control$cross^2 * l_dec / .ar_scale(cs[, j], lag)^2)))
  dc <- control$cross^2 * l_dec / .ar_scale(xs, lag)^2
  out <- .minnesota_core(yd, matrix(1, length(yd), 1L), cbind(own, cz), ds, Xc, dc, control, prior_prob, cores)
  out$settings$lag <- lag
  out
}

# the Minnesota comparison of "cause lags in" (Xc) against "out" for the
# response yd: Xfix = regressors with a vague prior (intercept, unit
# effects), Xs = shrunk nuisance regressors with base prior variances ds,
# Xc = the cause's lags with base variances dc; the prior variance of a
# shrunk coefficient is lambda^2 times its base variance (times sigma^2)
.minnesota_core <- function(yd, Xfix, Xs, ds, Xc, dc, control, prior_prob, cores = 1L, draws = TRUE) {
  N <- length(yd); a0 <- 0.01; b0 <- 0.01
  if (identical(control$lambda, "auto")) {
    g <- .gamma_from_mode_sd(control$lambda_hyper[["mode"]], control$lambda_hyper[["sd"]])
    lam <- exp(seq(log(0.005), log(10), length.out = 80L))
    lw <- stats::dgamma(lam, shape = g[["shape"]], scale = g[["scale"]], log = TRUE) + log(lam)  # log grid: d lambda = lambda d log lambda
  } else { lam <- as.numeric(control$lambda); lw <- 0 }
  X0 <- cbind(Xfix, Xs); XU <- cbind(X0, Xc)
  if (is.null(colnames(Xc))) colnames(Xc) <- paste0("lag", seq_len(ncol(Xc)))
  dvec <- function(l2, with_cause) c(rep(1e4, ncol(Xfix)), l2 * ds, if (with_cause) l2 * dc)
  mlU <- vapply(lam, function(l) .log_ml_nig(yd, XU, dvec(l^2, TRUE), a0, b0), 0)
  mlR <- vapply(lam, function(l) .log_ml_nig(yd, X0, dvec(l^2, FALSE), a0, b0), 0)
  lbf <- .logsumexp(lw + mlU) - .logsumexp(lw + mlR)
  bf <- exp(max(min(lbf, 700), -700)); post <- bf * prior_prob / (bf * prior_prob + 1 - prior_prob)
  wl <- exp(lw + mlU - max(lw + mlU)); wl <- wl / sum(wl)            # posterior of lambda under H1
  lam_post <- c(posterior_mean = sum(wl * lam), posterior_mode = lam[which.max(wl)])
  if (!draws) return(list(bf10 = bf, log_bf10 = lbf, post_prob = post, lambda = lam_post))
  q <- ncol(Xc); k0 <- ncol(X0)
  n_keep <- ceiling(control$n_iter / control$thin)
  run_chain <- function(ch) {
    li <- sample.int(length(lam), n_keep, replace = TRUE, prob = wl)
    out <- matrix(0, n_keep, q + 3L); colnames(out) <- c(colnames(Xc), "sigma", "geweke", "lambda")
    for (u in unique(li)) {
      rows <- which(li == u); d <- dvec(lam[u]^2, TRUE)
      A <- crossprod(XU) + diag(1 / d, ncol(XU)); V <- chol2inv(chol(A)); Xy <- crossprod(XU, yd); m <- V %*% Xy
      b <- b0 + (sum(yd^2) - sum(m * Xy)) / 2
      LV <- t(chol(V))
      for (r in rows) {
        sig2 <- 1 / stats::rgamma(1, a0 + N / 2, rate = max(b, 1e-12))
        beta <- as.numeric(m + sqrt(sig2) * LV %*% stats::rnorm(ncol(XU)))
        bc <- beta[k0 + seq_len(q)]
        out[r, ] <- c(bc, sqrt(sig2), log(1 + stats::var(as.numeric(Xc %*% bc)) / sig2), lam[u])
      }
    }
    out
  }
  ch <- .plapply(seq_len(control$chains), run_chain, cores = cores)
  dr <- do.call(rbind, ch)
  ci <- function(v) c(mean = mean(v), lower = unname(stats::quantile(v, (1 - control$credible) / 2)),
                      upper = unname(stats::quantile(v, 1 - (1 - control$credible) / 2)))
  structure(list(bf10 = bf, post_prob = post, inclusion = NULL, geweke = ci(dr[, "geweke"]),
                 coefficients = t(apply(dr[, seq_len(q), drop = FALSE], 2, ci)), draws = dr, ess = nrow(dr),
                 rhat = .gelman_rhat(lapply(ch, function(d) d[, "geweke"])), n = N, lambda = lam_post,
                 settings = list(prior = "minnesota", lag = q, prior_prob = prior_prob, chains = control$chains,
                                 n_iter = control$n_iter, burnin = 0L, thin = control$thin, credible = control$credible,
                                 sampler = "exact (independent draws)", lambda = control$lambda, decay = control$decay,
                                 cross = control$cross)),
            class = "bayes_granger")
}
