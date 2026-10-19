# Granger-type directed test (one ordered pair, optional controls) ----------

# Returns F-test (classical) and BIC-approximated Bayes factor (BF10, evidence
# for "cause helps predict effect").
.granger_pair <- function(effect, cause, controls = NULL, lag) {
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
  bic_r <- N * log(rss_r / N) + kr * log(N)
  bic_u <- N * log(rss_u / N) + ku * log(N)
  log_bf <- (bic_r - bic_u) / 2
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
#' @export
#' @examples
#' bf_evidence(c(0.02, 0.5, 2, 15, 200))
bf_evidence <- function(bf) {
  cut(bf, c(-Inf, 1/30, 1/10, 1/3, 1, 3, 10, 30, 100, Inf),
      labels = c("very strong against", "strong against", "moderate against",
                 "anecdotal against", "anecdotal for", "moderate for",
                 "strong for", "very strong for", "decisive for"),
      right = FALSE)
}
