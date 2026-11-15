# Granger-type directed test (one ordered pair, optional controls) ----------

# Returns F test (classical) and BIC-approximated Bayes factor (BF10, evidence
# for "cause helps predict effect").
# Heteroskedasticity.  The F test presupposes errors whose variance does not depend on the regressors.  When it does --
# for example when the effect depends on the square of the cause, so that the linear regression is misspecified -- the
# F distribution understates the variation of the estimated coefficients and the test rejects without a linear link.
#  * .het_test(): the test of White (1980) in the studentized form of Koenker (1981), without cross products: n R^2 of
#    the regression of the squared residuals on the regressors and their squares, referred to chi-square;
#  * .hc_wald(): the Wald test of the tested coefficients with the heteroskedasticity-consistent covariance HC3
#    (MacKinnon and White 1985), divided by the number of restrictions and referred to F with the residual degrees of
#    freedom, as the classical statistic is.  With hac = TRUE the autocovariances of the scores (regressors times
#    leverage-adjusted residuals) are added with Bartlett weights up to lag floor(4 (N / 100)^(2 / 9)) (Newey and West
#    1987, 1994): the errors of a misspecified regression are not martingale differences, the scores are then serially
#    correlated, and a covariance for independent scores understates the variation of the coefficients as well.
.het_test <- function(e, X) {
  Z <- cbind(1, X, X^2); u <- e^2; q <- qr(Z); k <- q$rank - 1L
  if (k < 1L || length(u) <= q$rank + 1L) return(NA_real_)
  r2 <- 1 - sum(qr.resid(q, u)^2) / sum((u - mean(u))^2)
  if (!is.finite(r2)) return(NA_real_)
  stats::pchisq(length(u) * r2, k, lower.tail = FALSE)
}
.hac_lag <- function(N) max(0L, min(N - 2L, as.integer(floor(4 * (N / 100)^(2 / 9)))))
.hc_wald <- function(y, X, idx, hac = FALSE) {
  X1 <- cbind(1, X); q <- qr(X1)
  if (q$rank < ncol(X1)) return(NULL)
  b <- qr.coef(q, y); e <- as.numeric(y - X1 %*% b); h <- pmin(rowSums(qr.Q(q)^2), 1 - 1e-8)
  S <- X1 * (e / (1 - h)); M <- crossprod(S); N <- length(y)
  if (hac) for (k in seq_len(.hac_lag(N))) {
    G <- crossprod(S[(k + 1L):N, , drop = FALSE], S[seq_len(N - k), , drop = FALSE])
    M <- M + (1 - k / (.hac_lag(N) + 1)) * (G + t(G))
  }
  Xi <- chol2inv(qr.R(q)); V <- Xi %*% M %*% Xi
  i <- idx + 1L
  W <- tryCatch(sum(b[i] * solve(V[i, i, drop = FALSE], b[i])), error = function(e) NA_real_)
  list(statistic = W / length(i), residuals = e)
}

# `own`: the number of lags of the effect itself (at least `lag`, the number of lags of the cause and of the controls).
# With the order chosen by the package it is the order of the autoregression of the effect (.own_order()), so that the
# regression holds what the effect remembers of itself; a numeric `lag` of the user gives the regression of that order.
.lags_upto <- function(v, lag, Q) .lagmat(v, Q)[, seq_len(lag), drop = FALSE]     # lags 1..lag on the time points Q + 1..n
.granger_pair <- function(effect, cause, controls = NULL, lag, prior = "bic", prior_par = NULL, vcov = "classical", own = lag) {
  n <- length(effect); Q <- max(as.integer(lag), as.integer(own), na.rm = TRUE)
  N <- n - Q
  yd <- effect[(Q + 1):n]
  own <- .lagmat(effect, Q)
  cz <- if (is.null(controls)) NULL else
    do.call(cbind, lapply(seq_len(ncol(controls)), function(j) .lags_upto(controls[, j], lag, Q)))
  Xr <- cbind(own, cz)
  Xu <- cbind(Xr, .lags_upto(cause, lag, Q))
  rss_r <- .rss(yd, Xr); rss_u <- .rss(yd, Xu)
  ku <- ncol(Xu) + 1L; kr <- ncol(Xr) + 1L
  q <- ku - kr; df2 <- N - ku
  if (df2 < 1 || rss_u <= 0) {
    return(list(statistic = NA_real_, df1 = q, df2 = df2, p_value = NA_real_,
                bf10 = NA_real_, effect_size = NA_real_, r2_gain = NA_real_,
                sign = NA_real_, n_eff = N, own = Q))
  }
  Fst <- ((rss_r - rss_u) / q) / (rss_u / df2)
  pv <- stats::pf(Fst, q, df2, lower.tail = FALSE)
  hc <- .hc_wald(yd, Xu, ncol(Xr) + seq_len(lag), hac = identical(vcov, "hac"))
  het_p <- if (is.null(hc)) NA_real_ else .het_test(hc$residuals, Xu)
  if (vcov %in% c("hc3", "hac") && !is.null(hc) && is.finite(hc$statistic)) {
    Fst <- hc$statistic; pv <- stats::pf(Fst, q, df2, lower.tail = FALSE)
  }
  log_bf <- if (prior == "bic") {
    bic_r <- N * log(rss_r / N) + kr * log(N)
    bic_u <- N * log(rss_u / N) + ku * log(N)
    (bic_r - bic_u) / 2
  } else {
    rss_0 <- sum((yd - mean(yd))^2)
    .log_bf_nested(N, k_r = ncol(Xr), k_u = ncol(Xu), R2_r = 1 - rss_r / rss_0, R2_u = 1 - rss_u / rss_0, prior = prior,
                   par = prior_par)
  }
  # sign of the summed cause coefficients (sign of the net influence)
  cf <- lm.fit(cbind(1, Xu), yd)$coefficients
  cc <- cf[(length(cf) - lag + 1):length(cf)]
  list(statistic = Fst, df1 = q, df2 = df2, p_value = pv,
       bf10 = exp(min(log_bf, 700)), log_bf10 = log_bf,
       effect_size = log(rss_r / rss_u),          # Geweke's F-measure
       r2_gain = 1 - rss_u / rss_r,
       sign = sign(sum(cc, na.rm = TRUE)), n_eff = N, het_p = het_p, own = Q)
}

#' Evidence category of a Bayes factor
#'
#' The categories of Jeffreys (1961) in the wording of Lee and Wagenmakers
#' (2013): "anecdotal" between 1 and 3, "moderate" between 3 and 10, "strong"
#' between 10 and 30, "very strong" between 30 and 100 and "decisive" above
#' 100, for `H1` when `BF10 > 1` and against it (with the reciprocal limits)
#' when `BF10 < 1`. On the scale of Kass and Raftery (1995) values between 3
#' and 20 are "positive", between 20 and 150 "strong" and above 150 "very
#' strong" evidence.
#'
#' @param bf numeric vector of Bayes factors (BF10).
#' @return Factor of evidence labels.
#' @references
#' Lee, M. D. and Wagenmakers, E.-J. (2013). *Bayesian Cognitive Modeling: A Practical Course*. Cambridge University Press, Cambridge. \doi{10.1017/CBO9781139087759}
#'
#' Kass, R. E. and Raftery, A. E. (1995). Bayes factors. *Journal of the American Statistical Association*, 90(430), 773-795. \doi{10.1080/01621459.1995.10476572}
#'
#' Jeffreys, H. (1961). *Theory of Probability* (3rd ed.). Oxford University Press, Oxford.
#'
#' Schwarz, G. (1978). Estimating the dimension of a model. *The Annals of Statistics*, 6(2), 461-464. \doi{10.1214/aos/1176344136}
#'
#' Wagenmakers, E.-J. (2007). A practical solution to the pervasive problems of p values. *Psychonomic Bulletin & Review*, 14, 779-804. \doi{10.3758/BF03194105}
#' @export
#' @examples
#' # BF10 < 1 favours H0, BF10 > 1 favours H1; the labels change at 3, 10, 30 and 100
#' bf_evidence(c(0.02, 0.5, 2, 15, 200))
bf_evidence <- function(bf) {
  cut(bf, c(-Inf, 1/100, 1/30, 1/10, 1/3, 1, 3, 10, 30, 100, Inf),
      labels = c("decisive against", "very strong against", "strong against", "moderate against",
                 "anecdotal against", "anecdotal for", "moderate for",
                 "strong for", "very strong for", "decisive for"),
      right = FALSE)
}


# Exact Bayes factors for nested Gaussian linear models ----------------------
# log BF of a model with k regressors (plus intercept) against the intercept-only
# model under Zellner's g-prior (g = N, unit-information) or the hyper-g prior
# (Liang et al. 2008, a = 3), from the coefficient of determination R2.
.log_bf_g0 <- function(N, k, R2, prior, g = N, a = 3) {
  R2 <- min(max(R2, 0), 1 - 1e-12)
  if (prior == "g") return((N - 1 - k) / 2 * log1p(g) - (N - 1) / 2 * log1p(g * (1 - R2)))
  # hyper-g: BF = (a-2)/2 * int_0^inf (1+g)^((N-1-k)/2 - a/2) (1 + g(1-R2))^(-(N-1)/2) dg,
  # with u = g/(1+g) on (0, 1): integrand (1-u)^((k+a)/2 - 2) (1 - u R2)^(-(N-1)/2)
  logf <- function(u) ((k + a) / 2 - 2) * log1p(-u) - (N - 1) / 2 * log1p(-u * R2)
  M <- max(logf(seq(0, 1 - 1e-9, length.out = 2001)))
  I <- stats::integrate(function(u) exp(logf(u) - M), 0, 1, rel.tol = 1e-8, subdivisions = 500L)$value
  log((a - 2) / 2) + M + log(I)
}
# Bayes factor of the unrestricted against the restricted regression under the g-prior family, with the restricted
# model as the base model (Liang et al. 2008; Bayarri et al. 2012): flat priors on the coefficients that the two
# models share (intercept, own lags, controls) and the g prior on the tested coefficients after the common regressors
# have been projected out, b | sigma^2 ~ N(0, g sigma^2 (Xc' M0 Xc)^{-1}).  Then
#   BF10 = (1 + g)^((nu0 - p) / 2) (1 + g (1 - R2p))^(-nu0 / 2),
# with nu0 = N - 1 - k_r the residual degrees of freedom of the restricted model, p = k_u - k_r tested coefficients
# and the partial coefficient of determination R2p = 1 - RSS_u / RSS_r, so that the Bayes factor depends on the data
# through the F statistic and its degrees of freedom only.  (The ratio of the Bayes factors of the two models against
# the intercept-only model, with the g prior on *all* coefficients, depends in addition on the fit of the restricted
# model: for a persistent effect series it tends to (1 + g)^(-p/2) whatever the evidence about the tested lags.)
# `par`: the settings of the prior, list(g_scale = , a = ): g = g_scale * N for the g-prior (g_scale = 1 is the
# unit-information prior) and the parameter a of the hyper-g prior; NULL = the defaults (g = N, a = 3)
.log_bf_nested <- function(N, k_r, k_u, R2_r, R2_u, prior, par = NULL) {
  g <- N * (if (is.null(par$g_scale)) 1 else par$g_scale); a <- if (is.null(par$a)) 3 else par$a
  R2p <- 1 - (1 - min(R2_u, 1)) / max(1 - R2_r, 1e-300)
  .log_bf_g0(N - k_r, k_u - k_r, R2p, prior, g = g, a = a)
}
# the prior settings of a bayes_control() list in that form
.prior_par <- function(control) list(g_scale = if (is.null(control$scale)) 1 else control$scale,
                                     a = if (is.null(control$hyper_a)) 3 else control$hyper_a)

#' Calibrated Bayes-factor bound from a p value
#'
#' The Sellke-Bayarri-Berger calibration: for a p value `p < 1/e` the Bayes
#' factor in favour of the alternative is at most `1 / (-e p log p)` over a
#' wide class of alternatives; the bound is `1` for larger p values. It
#' turns the permutation / surrogate p values of the model-free methods
#' into a Bayes factor (an upper bound, i.e. the most favourable reading of
#' the evidence) and a posterior probability of a link.
#'
#' @param p numeric vector of p values.
#' @param prior_prob prior probability of the alternative (a causal link).
#' @return A data frame with `bf10` (the bound) and `post_prob`.
#' @references
#' Sellke, T., Bayarri, M. J. and Berger, J. O. (2001). Calibration of p values for testing precise null hypotheses. *The American Statistician*, 55(1), 62-71. \doi{10.1198/000313001300339950}
#' @export
#' @examples
#' # the Sellke-Bayarri-Berger bound: p = 0.05 is at most a Bayes factor of 2.5
#' bf_from_p(c(0.001, 0.01, 0.05, 0.2, 0.5))
bf_from_p <- function(p, prior_prob = 0.5) {
  bf <- ifelse(is.na(p), NA_real_, ifelse(p < exp(-1), 1 / (-exp(1) * pmax(p, 1e-300) * log(pmax(p, 1e-300))), 1))
  data.frame(bf10 = bf, post_prob = bf * prior_prob / (bf * prior_prob + 1 - prior_prob))
}
