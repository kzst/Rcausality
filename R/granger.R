# Granger-type directed test (one ordered pair, optional controls) ----------

# Returns F-test (classical) and BIC-approximated Bayes factor (BF10, evidence
# for "cause helps predict effect").
.granger_pair <- function(effect, cause, controls = NULL, lag, prior = "bic") {
  n <- length(effect)
  N <- n - lag
  yd <- effect[(lag + 1):n]
  own <- .lagmat(effect, lag)
  cz <- if (is.null(controls)) NULL else
    do.call(cbind, lapply(seq_len(ncol(controls)), function(j) .lagmat(controls[, j], lag)))
  Xr <- cbind(own, cz)
  Xu <- cbind(Xr, .lagmat(cause, lag))
  rss_r <- .rss(yd, Xr); rss_u <- .rss(yd, Xu)
  ku <- ncol(Xu) + 1L; kr <- ncol(Xr) + 1L
  q <- ku - kr; df2 <- N - ku
  if (df2 < 1 || rss_u <= 0) {
    return(list(statistic = NA_real_, df1 = q, df2 = df2, p_value = NA_real_,
                bf10 = NA_real_, effect_size = NA_real_, r2_gain = NA_real_,
                sign = NA_real_, n_eff = N))
  }
  Fst <- ((rss_r - rss_u) / q) / (rss_u / df2)
  pv <- stats::pf(Fst, q, df2, lower.tail = FALSE)
  log_bf <- if (prior == "bic") {
    bic_r <- N * log(rss_r / N) + kr * log(N)
    bic_u <- N * log(rss_u / N) + ku * log(N)
    (bic_r - bic_u) / 2
  } else {
    rss_0 <- sum((yd - mean(yd))^2)
    .log_bf_nested(N, k_r = ncol(Xr), k_u = ncol(Xu), R2_r = 1 - rss_r / rss_0, R2_u = 1 - rss_u / rss_0, prior = prior)
  }
  # sign of the summed cause coefficients (sign of the net influence)
  cf <- lm.fit(cbind(1, Xu), yd)$coefficients
  cc <- cf[(length(cf) - lag + 1):length(cf)]
  list(statistic = Fst, df1 = q, df2 = df2, p_value = pv,
       bf10 = exp(min(log_bf, 700)), log_bf10 = log_bf,
       effect_size = log(rss_r / rss_u),          # Geweke's F-measure
       r2_gain = 1 - rss_u / rss_r,
       sign = sign(sum(cc, na.rm = TRUE)), n_eff = N)
}

#' Evidence category of a Bayes factor (Jeffreys/Kass-Raftery style)
#'
#' @param bf numeric vector of Bayes factors (BF10).
#' @return Character vector of evidence labels.
#' @references
#' Kass, R. E. and Raftery, A. E. (1995). Bayes factors. *Journal of the American Statistical Association*, 90(430), 773-795. \doi{10.1080/01621459.1995.10476572}
#'
#' Jeffreys, H. (1961). *Theory of Probability* (3rd ed.). Oxford University Press, Oxford.
#'
#' Schwarz, G. (1978). Estimating the dimension of a model. *The Annals of Statistics*, 6(2), 461-464. \doi{10.1214/aos/1176344136}
#'
#' Wagenmakers, E.-J. (2007). A practical solution to the pervasive problems of p values. *Psychonomic Bulletin & Review*, 14, 779-804. \doi{10.3758/BF03194105}
#' @export
#' @examples
#' # Kass-Raftery labels: BF10 < 1 favours H0, > 3 is positive evidence for H1
#' bf_evidence(c(0.02, 0.5, 2, 15, 200))
bf_evidence <- function(bf) {
  cut(bf, c(-Inf, 1/30, 1/10, 1/3, 1, 3, 10, 30, 100, Inf),
      labels = c("very strong against", "strong against", "moderate against",
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
.log_bf_nested <- function(N, k_r, k_u, R2_r, R2_u, prior) {
  .log_bf_g0(N, k_u, R2_u, prior) - .log_bf_g0(N, k_r, R2_r, prior)
}

#' Calibrated Bayes-factor bound from a p-value
#'
#' The Sellke-Bayarri-Berger calibration: for a p-value `p < 1/e` the Bayes
#' factor in favour of the alternative is at most `1 / (-e p log p)` over a
#' wide class of alternatives; the bound is `1` for larger p-values. It
#' turns the permutation / surrogate p-values of the model-free methods
#' into a Bayes factor (an upper bound, i.e. the most favourable reading of
#' the evidence) and a posterior probability of a link.
#'
#' @param p numeric vector of p-values.
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
