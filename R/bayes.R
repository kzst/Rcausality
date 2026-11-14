# Bayesian variants of the causality methods ---------------------------------
#
# Three families:
#  (1) regression (Granger) methods: Bayesian linear regression of the effect
#      on its own past and the cause's past -- closed-form Bayes factors (BIC,
#      Zellner g-prior, hyper-g), a block comparison with a normal slab
#      ("normal") and stochastic-search variable selection over the cause's
#      lags ("spike_slab", George & McCulloch 1993, 1997; exact enumeration of
#      the lag subsets, Gibbs sampling for many lags), with posterior draws of
#      the coefficients and of Geweke's measure;
#  (2) statistic-versus-surrogates methods (CMI, transfer entropy, MAR-PID,
#      df, CCDH, recurrence plots, CCM, Markov tests): the observed statistic
#      is compared with the surrogate / permutation sample through a Bayes
#      factor with an explicit prior on the effect size (a Bayesian version of
#      the surrogate test: H0 = the statistic comes from the surrogate
#      distribution, H1 = it is shifted by an effect delta with a half-normal
#      or half-Cauchy prior scaled by the surrogate spread);
#  (3) posterior-based methods (dimensional causality, df-entropy): a prior
#      over the five causal cases replaces the uniform one.

#' Settings of the Bayesian variants of the causality methods
#'
#' Collects the options of `rcausality(bayes = TRUE)`: the prior of the
#' regression (Granger) family, the effect-size prior and null model of the
#' surrogate family, the prior over the causal cases of the posterior-based
#' methods, and the Monte Carlo / MCMC sizes.
#'
#' @param prior prior of the Granger regression comparison: `"spike_slab"`
#'   (stochastic-search variable selection over the cause's lags: every lag
#'   has a point-mass spike at zero with prior inclusion probability
#'   `inclusion_prob` and a normal slab of standard deviation `scale * sigma`;
#'   the link is the event that at least one lag is included), `"normal"`
#'   (block comparison: all cause lags with the normal slab against none),
#'   `"g"` (Zellner g-prior, `g = scale * n` when `scale` is given, else
#'   `g = n`), `"hyper_g"` (Liang et al. 2008, parameter `hyper_a`) or
#'   `"bic"` (Schwarz approximation).
#' @param effect_prior prior of the standardised effect size `delta / sd0`
#'   of the surrogate family (`sd0` = standard deviation of the surrogate
#'   statistics): `"half_cauchy"` (default; heavy tails make the Bayes
#'   factor and the effect posterior robust when the effect is far larger
#'   than the surrogate spread, as in the JZS t-test of Rouder et al. 2009),
#'   `"half_normal"` (one-sided; the alternative shifts the statistic in the
#'   direction of the effect) or `"normal"`, `"cauchy"` (two-sided).
#' @param scale prior scale: the slab standard deviation relative to `sigma`
#'   (`"spike_slab"`, `"normal"`), `g / n` for `"g"`, and the scale of the
#'   effect-size prior in units of the surrogate spread (1 = a shift of one
#'   surrogate standard deviation is the typical effect; 0.5 "narrow", 2
#'   "wide").
#' @param inclusion_prob prior probability that a given cause lag enters the
#'   regression (`"spike_slab"`).
#' @param null_model distribution of the statistic under the null of the
#'   surrogate family: `"kde"` (kernel density of the surrogate sample,
#'   default; robust to skewed nulls) or `"normal"` (normal fit; closed-form
#'   Bayes factor with the half-normal prior).
#' @param case_prior prior probabilities of the five causal cases
#'   (`x_to_y`, `y_to_x`, `bidirectional`, `common_cause`, `independent`) for
#'   the dimensional-causality and df-entropy methods and for the Bayesian
#'   Markov decision table. It replaces the prior of the method, which is
#'   used when `case_prior = NULL`: uniform over the five cases for the
#'   df-entropy classifier, uniform over the nine elementary events of Benko
#'   et al. (2024) for dimensional causality (weights 1 : 1 : 1 : 3 : 3 of
#'   the five cases), and the prior implied by equally probable outcomes of
#'   the six tests for the Markov table (2 : 2 : 2 : 1 : 1).
#' @param n_iter,burnin,thin,chains MCMC settings of the regression family:
#'   posterior draws of the coefficients, of `sigma` and of Geweke's measure
#'   (`n_iter` draws per chain after `burnin`, every `thin`-th kept; the
#'   inclusion indicators are enumerated exactly when the cause has at most
#'   12 lags and Gibbs-sampled otherwise).
#' @param n_draws Monte Carlo draws of the effect-size prior used to
#'   integrate the surrogate-family Bayes factor (and to summarise the
#'   posterior of the effect).
#' @param hyper_a parameter `a` of the hyper-g prior (2 < a <= 4).
#' @param credible level of the credible intervals.
#' @param lambda,decay,cross,own_mean,lambda_hyper settings of the
#'   Minnesota prior (`prior = "minnesota"`) for short series: the prior
#'   standard deviation of the coefficient of series `j` at lag `l` is
#'   `lambda * cross_j / l^decay * sigma / s_j` (Litterman 1986; `s_j` = the
#'   residual sd of an autoregression of series `j`, `cross_j = 1` for the
#'   effect's own lags and `cross` for the other series), the prior mean of
#'   the effect's first own lag is `own_mean` (0 for stationary,
#'   standardised series; 1 = random-walk prior), and the overall tightness
#'   `lambda` is either fixed (a number) or `"auto"`: integrated out under a
#'   Gamma hyperprior with the mode and standard deviation `lambda_hyper`
#'   (0.2 and 0.4, as in Giannone, Lenza and Primiceri 2015), on a grid, so
#'   that the Bayes factor and the posterior draws average over the
#'   tightness that the data support.
#' @param decision how a link is declared: `"bf"` (Bayes factor above
#'   `bf_threshold` of [rcausality()]) or `"posterior"` (posterior probability
#'   above `prob_threshold`).
#' @param effect_scale unit of the scale of the effect prior of the surrogate
#'   Bayes factor: `"local"` (default; `scale` null standard deviations of the
#'   statistic) or `"unit"` (`scale` null standard deviations of the
#'   statistic times the root of the number of observations, the scale of a
#'   unit-information prior for a statistic whose null standard deviation
#'   falls with the root of the length; see [surrogate_bf()]). With `"unit"`
#'   the evidence for the absence of a link grows with the length of the
#'   series.
#' @return A list of class `bayes_control`.
#' @references
#' George, E. I. and McCulloch, R. E. (1993). Variable selection via Gibbs sampling. *Journal of the American Statistical Association*, 88(423), 881-889. \doi{10.1080/01621459.1993.10476353}
#'
#' George, E. I. and McCulloch, R. E. (1997). Approaches for Bayesian variable selection. *Statistica Sinica*, 7, 339-373.
#'
#' Liang, F., Paulo, R., Molina, G., Clyde, M. A. and Berger, J. O. (2008). Mixtures of g priors for Bayesian variable selection. *Journal of the American Statistical Association*, 103(481), 410-423. \doi{10.1198/016214507000001337}
#'
#' Rouder, J. N., Speckman, P. L., Sun, D., Morey, R. D. and Iverson, G. (2009). Bayesian t tests for accepting and rejecting the null hypothesis. *Psychonomic Bulletin & Review*, 16, 225-237. \doi{10.3758/PBR.16.2.225}
#'
#' Kass, R. E. and Raftery, A. E. (1995). Bayes factors. *Journal of the American Statistical Association*, 90(430), 773-795. \doi{10.1080/01621459.1995.10476572}
#'
#' Litterman, R. B. (1986). Forecasting with Bayesian vector autoregressions: five years of experience. *Journal of Business & Economic Statistics*, 4(1), 25-38. \doi{10.1080/07350015.1986.10509491}
#'
#' Giannone, D., Lenza, M. and Primiceri, G. E. (2015). Prior selection for vector autoregressions. *The Review of Economics and Statistics*, 97(2), 436-451. \doi{10.1162/REST_a_00483}
#' @export
#' @examples
#' # the defaults: spike-and-slab regression prior, half-Cauchy effect prior,
#' # two chains of 2000 iterations, decision by the Bayes factor
#' bayes_control()
#' # a wider slab, longer chains, decision by the posterior probability
#' bayes_control(prior = "spike_slab", scale = 2, n_iter = 4000, chains = 4, decision = "posterior")
#' # the hierarchical Minnesota prior for short series (tightness integrated out)
#' bayes_control(prior = "minnesota")
#' # a prior over the causal cases (for dc, dfe and the Bayesian Markov table)
#' bayes_control(case_prior = c(x_to_y = 0.3, y_to_x = 0.3, bidirectional = 0.1,
#'                              common_cause = 0.1, independent = 0.2))
bayes_control <- function(prior = c("spike_slab", "normal", "minnesota", "g", "hyper_g", "bic"),
                          effect_prior = c("half_cauchy", "half_normal", "normal", "cauchy"),
                          scale = 1, inclusion_prob = 0.5, null_model = c("kde", "normal"),
                          case_prior = NULL, n_iter = 2000L, burnin = 500L, thin = 1L, chains = 2L,
                          n_draws = 4000L, hyper_a = 3, credible = 0.95,
                          lambda = "auto", decay = 1, cross = 1, own_mean = 0, lambda_hyper = c(mode = 0.2, sd = 0.4),
                          decision = c("bf", "posterior"), effect_scale = c("local", "unit")) {
  prior <- match.arg(prior); effect_prior <- match.arg(effect_prior); null_model <- match.arg(null_model)
  decision <- match.arg(decision); effect_scale <- match.arg(effect_scale)
  .check_number(scale, "scale", lower = 0, strict = TRUE)
  .check_number(inclusion_prob, "inclusion_prob", lower = 0, upper = 1, strict = TRUE)
  .check_number(credible, "credible", lower = 0, upper = 1, strict = TRUE)
  if (!identical(lambda, "auto")) .check_number(lambda, "lambda", lower = 0, strict = TRUE)
  .check_number(decay, "decay", lower = 0); .check_number(cross, "cross", lower = 0, strict = TRUE)
  if (length(lambda_hyper) != 2L || any(!is.finite(lambda_hyper)) || any(lambda_hyper <= 0))
    stop("`lambda_hyper` must hold the positive mode and standard deviation of the tightness hyperprior.", call. = FALSE)
  if (hyper_a <= 2 || hyper_a > 4) stop("`hyper_a` must lie in (2, 4].", call. = FALSE)
  if (!is.null(case_prior)) {
    if (length(case_prior) != 5L || any(case_prior < 0) || sum(case_prior) <= 0)
      stop("case_prior must hold five nonnegative probabilities.", call. = FALSE)
    case_prior <- case_prior / sum(case_prior)
    names(case_prior) <- c("x_to_y", "y_to_x", "bidirectional", "common_cause", "independent")
  }
  structure(list(prior = prior, effect_prior = effect_prior, scale = scale, effect_scale = effect_scale, inclusion_prob = inclusion_prob,
                 null_model = null_model, case_prior = case_prior, n_iter = as.integer(n_iter),
                 burnin = as.integer(burnin), thin = max(1L, as.integer(thin)), chains = max(1L, as.integer(chains)),
                 n_draws = as.integer(n_draws), hyper_a = hyper_a, credible = credible,
                 lambda = lambda, decay = decay, cross = cross, own_mean = own_mean,
                 lambda_hyper = stats::setNames(as.numeric(lambda_hyper), c("mode", "sd")), decision = decision),
            class = "bayes_control")
}

#' @export
print.bayes_control <- function(x, ...) {
  .out("Bayesian settings of rcausality():\n")
  .out(sprintf("  regression family: prior = %s, scale = %g, inclusion_prob = %g, hyper_a = %g\n", x$prior, x$scale, x$inclusion_prob, x$hyper_a))
  if (identical(x$prior, "minnesota"))
    .out(sprintf("  Minnesota prior: lambda = %s, decay = %g, cross = %g, own_mean = %g\n",
                if (identical(x$lambda, "auto")) sprintf("auto (Gamma hyperprior, mode %g, sd %g)", x$lambda_hyper[["mode"]], x$lambda_hyper[["sd"]])
                else format(x$lambda), x$decay, x$cross, x$own_mean))
  .out(sprintf("  surrogate family:  effect_prior = %s (scale %g x surrogate sd), null_model = %s, n_draws = %d\n",
              x$effect_prior, x$scale, x$null_model, x$n_draws))
  .out(sprintf("  case prior (dc, dfe, markov): %s\n", if (is.null(x$case_prior)) "uniform" else paste(round(x$case_prior, 3), collapse = " ")))
  .out(sprintf("  MCMC: %d chains x %d iterations, burn-in %d, thinning %d; credible level %g; decision by %s\n",
              x$chains, x$n_iter, x$burnin, x$thin, x$credible, x$decision))
  invisible(x)
}

# ---- (2) Bayes factor of a statistic against its surrogate sample ------------
# the unit of the effect scale of a control list and the number of observations behind a statistic
.eff_scale <- function(control) if (is.null(control$effect_scale)) "local" else control$effect_scale
.eff_n <- function(n_eff, n) { v <- suppressWarnings(as.numeric(n_eff)); if (length(v) == 1L && is.finite(v) && v >= 1) v else as.numeric(n) }

#' Bayes factor of an observed statistic against surrogate / permutation draws
#'
#' A Bayesian version of the surrogate test used by the model-free causality
#' methods. Under `H0` the observed statistic `T` is a draw from the
#' surrogate distribution `f0` (a kernel density or a normal fit of the
#' surrogate sample); under `H1` it is shifted by an effect `delta` with a
#' prior on the standardised effect `delta / sd0` (`sd0` = surrogate standard
#' deviation): half-Cauchy (default) or half-normal for a one-sided
#' alternative, Cauchy or normal for a two-sided one. The Bayes factor is `f1(T) / f0(T)` with
#' `f1(T) = E_delta[f0(T - delta)]`, integrated by Monte Carlo (closed form for
#' the normal null with the half-normal prior). The posterior of the effect
#' under `H1` is summarised from the importance-weighted prior draws.
#'
#' The construction is a Bayes factor based on a test statistic (Johnson,
#' 2005): it conditions on the statistic and not on the data, and it assumes
#' that the effect shifts the null distribution without changing its shape.
#' Outside the range of the surrogate statistics the kernel density is
#' continued by a normal tail (the element `extrapolated` of the result is
#' then `TRUE`): the replicates carry no information about the null density
#' there, so that such Bayes factors rest on the normal tail.
#'
#' @param statistic the observed statistic.
#' @param null numeric vector of surrogate statistics (at least 10), or a
#'   list `list(mean = , sd = )` describing a normal null.
#' @param alternative direction of the effect: `"greater"` (the statistic
#'   is larger under coupling), `"less"` or `"two.sided"`.
#' @param effect_prior,scale,null_model,n_draws,credible see
#'   [bayes_control()].
#' @param prior_prob prior probability of `H1`.
#' @param seed optional RNG seed for the Monte Carlo integration.
#' @param effect_scale,n_eff the unit of the prior scale of the effect:
#'   `"local"` (default), `scale` null standard deviations of the statistic,
#'   or `"unit"`, `scale * sd0 * sqrt(n_eff)`, where `n_eff` is the number
#'   of observations behind the statistic. For a statistic whose null
#'   standard deviation falls like `1 / sqrt(n_eff)` this is the scale of a
#'   unit-information prior (Kass and Wasserman 1995); the plug-in statistics
#'   of the nonlinear methods need not have that rate, and the option is
#'   then a wider prior scale that grows with the length, not a
#'   unit-information prior in the strict sense. With the
#'   local scale the test is most sensitive, but the evidence it can give for
#'   `H0` is bounded (about `BF10 = 0.5` at the centre of the null
#'   distribution, whatever the length of the series); with the unit scale
#'   the Bayes factor in favour of `H0` grows like `sqrt(n_eff)` when `H0`
#'   holds, as for the BIC, at the price of a later detection of weak effects.
#' @return An object of class `surrogate_bf`: `bf10`, `post_prob`,
#'   `effect` (posterior mean and credible interval of the effect under `H1`,
#'   in the units of the statistic), `z` (the statistic in surrogate standard
#'   deviations), `statistic`, `null_mean` and `null_sd` (mean and standard
#'   deviation of the surrogate statistics), `extrapolated` (`TRUE` when the
#'   statistic lies outside the range of the surrogate statistics, where the
#'   density of the null distribution is the normal tail and not an estimate
#'   from them), `bound` (the Bayes factor at the most extreme surrogate
#'   statistic on the side of the observed one: with `extrapolated = TRUE`
#'   the surrogates support `bf10 > bound` and no more; `NA` otherwise), and
#'   `settings`. Method: `print`, which shows an extrapolated Bayes factor as
#'   the bound.
#' @references
#' Rouder, J. N., Speckman, P. L., Sun, D., Morey, R. D. and Iverson, G. (2009). Bayesian t tests for accepting and rejecting the null hypothesis. *Psychonomic Bulletin & Review*, 16, 225-237. \doi{10.3758/PBR.16.2.225}
#'
#' Kass, R. E. and Raftery, A. E. (1995). Bayes factors. *Journal of the American Statistical Association*, 90(430), 773-795. \doi{10.1080/01621459.1995.10476572}
#'
#' Kass, R. E. and Wasserman, L. (1995). A reference Bayesian test for nested hypotheses and its relationship to the Schwarz criterion. *Journal of the American Statistical Association*, 90(431), 928-934. \doi{10.1080/01621459.1995.10476592}
#'
#' Johnson, V. E. (2005). Bayes factors based on test statistics. *Journal of the Royal Statistical Society B*, 67(5), 689-701. \doi{10.1111/j.1467-9868.2005.00521.x}
#' @export
#' @examples
#' # A statistic observed 3.2 surrogate standard deviations above a standard
#' # normal null: strong evidence for an effect (BF10 >> 1); the posterior of
#' # the effect size (in surrogate sd units) is reported with it.
#' set.seed(1)
#' surrogate_bf(3.2, rnorm(200))
#' # A statistic inside the null distribution: BF10 < 1 is evidence *for* H0,
#' # which a p-value cannot express.
#' surrogate_bf(0.3, rnorm(200), scale = 0.5)
surrogate_bf <- function(statistic, null, alternative = c("greater", "less", "two.sided"),
                         effect_prior = c("half_cauchy", "half_normal", "normal", "cauchy"), scale = 1,
                         null_model = c("kde", "normal"), prior_prob = 0.5, n_draws = 4000L, credible = 0.95,
                         seed = NULL, effect_scale = c("local", "unit"), n_eff = NULL) {
  alternative <- match.arg(alternative); effect_prior <- match.arg(effect_prior); null_model <- match.arg(null_model)
  effect_scale <- match.arg(effect_scale)
  if (effect_scale == "unit" && (is.null(n_eff) || !is.finite(n_eff) || n_eff < 1))
    stop("effect_scale = \"unit\" needs the number of observations behind the statistic (n_eff).", call. = FALSE)
  .local_seed(seed)
  if (is.list(null)) { m0 <- null$mean; s0 <- null$sd; null_model <- "normal"; draws <- NULL }
  else {
    draws <- as.numeric(null); draws <- draws[is.finite(draws)]
    if (length(draws) < 10L) stop("At least 10 surrogate statistics are needed.", call. = FALSE)
    m0 <- mean(draws); s0 <- stats::sd(draws)
  }
  s0 <- max(s0, 1e-8 * (abs(m0) + 1))
  T <- as.numeric(statistic)
  if (!is.finite(T)) return(structure(list(bf10 = NA_real_, post_prob = NA_real_, effect = c(mean = NA, lower = NA, upper = NA),
                                           z = NA_real_, settings = list()), class = "surrogate_bf"))
  sgn <- switch(alternative, greater = 1, less = -1, two.sided = 0)
  one_sided <- sgn != 0 && effect_prior %in% c("half_normal", "half_cauchy")
  if (sgn == 0 && effect_prior %in% c("half_normal", "half_cauchy")) effect_prior <- sub("half_", "", effect_prior)
  # scale of the effect prior: `scale` null standard deviations ("local": alternatives of the size of the sampling
  # error, the most sensitive choice, under which the evidence for H0 stays bounded), or `scale` null standard
  # deviations of ONE observation, s0 * sqrt(n_eff) ("unit": a unit-information prior, under which the Bayes factor
  # in favour of H0 grows like sqrt(n_eff) when H0 holds, as for the BIC)
  tau <- scale * s0 * (if (effect_scale == "unit") sqrt(n_eff) else 1)
  # null density.  Inside the range of the replicates: their kernel density.  Outside it the replicates carry no
  # information, and the tails of the kernel estimate are those of the kernel with the bandwidth as scale, which
  # decay far too fast (Bayes factors above the bound of the shift model); there the estimate is continued by the
  # tail of the normal distribution with the mean and the standard deviation of the replicates, matched at the
  # extreme replicate.
  extrapolated <- FALSE
  f0 <- if (null_model == "normal") function(v) stats::dnorm(v, m0, s0) else {
    bw <- stats::bw.nrd0(draws); lo_d <- min(draws); hi_d <- max(draws)
    kde <- function(u) mean(stats::dnorm(u, draws, bw))
    c_lo <- kde(lo_d) / stats::dnorm(lo_d, m0, s0); c_hi <- kde(hi_d) / stats::dnorm(hi_d, m0, s0)
    extrapolated <- T > hi_d || T < lo_d
    function(v) vapply(v, function(u) if (u > hi_d) c_hi * stats::dnorm(u, m0, s0) else if (u < lo_d) c_lo * stats::dnorm(u, m0, s0)
                                      else kde(u), 0)
  }
  f0T <- f0(T)
  # H1 density: closed form for the normal null with the half-normal / normal prior
  if (null_model == "normal" && effect_prior %in% c("half_normal", "normal")) {
    sd1 <- sqrt(s0^2 + tau^2)
    f1T <- if (effect_prior == "half_normal") 2 * stats::dnorm(T, m0, sd1) * stats::pnorm(sgn * (T - m0) * tau / (s0 * sd1))
           else stats::dnorm(T, m0, sd1)
    # posterior of delta under H1 (normal part): N(mu_d, v_d), truncated for the half-normal
    v_d <- 1 / (1 / s0^2 + 1 / tau^2); mu_d <- v_d * (T - m0) / s0^2
    dd <- stats::rnorm(n_draws, mu_d, sqrt(v_d)); if (one_sided) dd <- dd[sgn * dd >= 0]
    w <- rep(1, length(dd))
  } else {
    # importance sampling of f1(T) = E_prior[f0(T - delta)] with a defensive
    # mixture proposal (half prior, half a normal around the observed shift),
    # so that effects far outside the prior scale are integrated accurately
    dprior <- function(v) switch(effect_prior,
      half_normal = ifelse(sgn * v >= 0, 2 * stats::dnorm(v, 0, tau), 0), normal = stats::dnorm(v, 0, tau),
      half_cauchy = ifelse(sgn * v >= 0, 2 * stats::dcauchy(v, 0, tau), 0), cauchy = stats::dcauchy(v, 0, tau))
    rprior <- function(k) switch(effect_prior,
      half_normal = abs(stats::rnorm(k, 0, tau)) * sgn, normal = stats::rnorm(k, 0, tau),
      half_cauchy = abs(stats::rcauchy(k, 0, tau)) * sgn, cauchy = stats::rcauchy(k, 0, tau))
    k1 <- n_draws %/% 2; k2 <- n_draws - k1; shift <- T - m0
    if (is.finite(T) && !is.finite(f0T)) f0T <- 0
    dd <- c(rprior(k1), stats::rnorm(k2, shift, s0))
    qd <- 0.5 * dprior(dd) + 0.5 * stats::dnorm(dd, shift, s0)
    w <- f0(T - dd) * dprior(dd) / qd; w[!is.finite(w)] <- 0
    f1T <- mean(w)
  }
  bf <- if (f0T > 0) f1T / f0T else if (f1T > 0) Inf else 1
  bf <- min(bf, 1e12)
  post <- bf * prior_prob / (bf * prior_prob + 1 - prior_prob)
  eff <- if (length(dd) && sum(w) > 0) {
    ws <- w / sum(w); mu <- sum(ws * dd)
    o <- order(dd); cw <- cumsum(ws[o]); q <- function(p) dd[o][which(cw >= p)[1]]
    c(mean = mu, lower = q((1 - credible) / 2), upper = q(1 - (1 - credible) / 2))
  } else c(mean = NA_real_, lower = NA_real_, upper = NA_real_)
  # Outside the range of the replicates the Bayes factor rests on the continuation of the tail, a convention.  What
  # the replicates themselves support is the Bayes factor at the extreme replicate, a lower bound of the Bayes factor
  # of a statistic beyond it (the generator of the session is left as it was).
  bound <- NA_real_
  if (isTRUE(extrapolated)) {
    st <- .rng_get()
    bound <- tryCatch(surrogate_bf(if (T > max(draws)) max(draws) else min(draws), draws, alternative = alternative,
                                   effect_prior = effect_prior, scale = scale, null_model = null_model, prior_prob = prior_prob,
                                   n_draws = n_draws, credible = credible, effect_scale = effect_scale, n_eff = n_eff)$bf10,
                      error = function(e) NA_real_)
    .rng_put(st)
  }
  structure(list(bf10 = bf, post_prob = post, effect = eff, z = (T - m0) / s0, statistic = T, null_mean = m0, null_sd = s0,
                 extrapolated = extrapolated, bound = bound,
                 settings = list(alternative = alternative, effect_prior = effect_prior, scale = scale, null_model = null_model,
                                 effect_scale = effect_scale, n_eff = n_eff,
                                 n_surr = length(draws), n_draws = n_draws, prior_prob = prior_prob, credible = credible)),
            class = "surrogate_bf")
}

#' @export
print.surrogate_bf <- function(x, digits = 3, ...) {
  s <- x$settings
  .out("Bayes factor of a statistic against its surrogate distribution\n")
  .out(sprintf("  statistic %s = %.3f surrogate sd from the null mean (%s surrogates, %s null, %s prior of scale %g, %s)\n",
              formatC(x$statistic, digits = digits, format = "g"), x$z, if (is.null(s$n_surr)) "-" else s$n_surr,
              s$null_model, s$effect_prior, s$scale, s$alternative))
  if (isTRUE(x$extrapolated) && is.finite(x$bound)) {
    .out(sprintf("  BF10 > %s (evidence: at least \"%s\" an effect). The statistic lies outside the range of the %s surrogates,\n",
                 .fmt_bf(x$bound), as.character(bf_evidence(x$bound)), if (is.null(s$n_surr)) "" else s$n_surr))
    .out(sprintf("   which support the Bayes factor at their extreme value only; with the tail of a normal distribution beyond it\n   BF10 = %s and the posterior probability of an effect is %.3f (prior %.2f).\n",
                 .fmt_bf(x$bf10), x$post_prob, s$prior_prob))
  } else .out(sprintf("  BF10 = %s (%s); posterior probability of an effect = %.3f (prior %.2f)\n", .fmt_bf(x$bf10),
              as.character(bf_evidence(x$bf10)), x$post_prob, s$prior_prob))
  .out(sprintf("  effect under H1: %.3g [%.3g, %.3g] (%.0f%% credible)\n", x$effect[["mean"]], x$effect[["lower"]],
              x$effect[["upper"]], 100 * s$credible))
  invisible(x)
}

# ---- (1) Bayesian Granger regression -----------------------------------------

# log marginal likelihood of y ~ N(X_S b, sigma^2) with b ~ N(0, sigma^2 D),
# sigma^2 ~ IG(a0, b0), up to a constant common to all subsets of the same y
# (`n`: the number of observations that enters the exponent; smaller than length(y) when regressors with flat
# priors have been projected out of y and X)
.log_ml_nig <- function(y, X, d, a0 = 0.01, b0 = 0.01, n = length(y)) {
  if (!ncol(X)) { b <- b0 + sum(y^2) / 2; return(-(a0 + n / 2) * log(b)) }
  A <- crossprod(X) + diag(1 / d, ncol(X))
  R <- tryCatch(chol(A), error = function(e) NULL); if (is.null(R)) return(-Inf)
  Xy <- crossprod(X, y); m <- backsolve(R, forwardsolve(t(R), Xy))
  b <- b0 + (sum(y^2) - sum(m * Xy)) / 2
  if (b <= 0) return(-Inf)
  -sum(log(diag(R))) - 0.5 * sum(log(d)) - (a0 + n / 2) * log(b)
}

#' Bayesian Granger causality for one directed pair
#'
#' Bayesian linear regression of the effect on its own `lag` past values
#' (and the control series) and on the cause's `lag` past values, all
#' standardised. The link "cause -> effect" is assessed by
#' * `prior = "spike_slab"`: stochastic-search variable selection over the
#'   cause's lags (spike = exclusion, slab = `N(0, (scale * sigma)^2)`,
#'   prior inclusion probability `inclusion_prob` per lag); the 2^q lag
#'   subsets are enumerated exactly when `q <= 12` (Gibbs sampling otherwise)
#'   and the Bayes factor compares "at least one lag" with "no lag";
#' * `prior = "normal"`: the block of all cause lags with the normal slab
#'   against none. For these two priors the coefficients that both
#'   regressions share (intercept, own lags, lags of the controls) have flat
#'   priors: the restricted regression is the base model, and the marginal
#'   likelihoods have `n - k0` degrees of freedom. A regression whose
#'   always-included coefficients leave no degrees of freedom stops with an
#'   error (use `"minnesota"`);
#' * `prior = "minnesota"`: the hierarchical Minnesota prior for short
#'   series (lag-decaying shrinkage of all coefficients, the overall
#'   tightness integrated out under a Gamma hyperprior; see [bayes_control()]);
#'   closed-form Bayes factor and exact independent posterior draws. The
#'   posterior is proper even when the regression has more coefficients than
#'   observations;
#' * `prior = "g"`, `"hyper_g"`, `"bic"`: the closed-form Bayes factors of
#'   [rcausality()].
#' Posterior draws of the coefficients, of `sigma` and of Geweke's measure
#' `log(1 + var(X_c b_c) / sigma^2)` (the Bayesian analogue of
#' `log(RSS_r / RSS_u)`; `X_c`: the lags of the cause after the
#' always-included regressors have been projected out) are generated for
#' `prior = "spike_slab"` and
#' `"normal"` (`chains` x `n_iter`, burn-in, thinning).
#'
#' @param effect,cause numeric vectors.
#' @param controls optional matrix of control series (conditional Granger).
#' @param lag lag order.
#' @param control a [bayes_control()] list.
#' @param prior_prob prior probability of the link (for `"normal"`, `"g"`,
#'   `"hyper_g"` and `"bic"`; the spike-and-slab prior implies its own,
#'   `1 - (1 - inclusion_prob)^lag`).
#' @param seed optional RNG seed.
#' @param cores CPU cores: the chains run in parallel (see [Rcausality-parallel]).
#' @return An object of class `bayes_granger`: `bf10`, `post_prob` (computed
#'   from the Bayes factor and `prior_prob`), `any_included` (posterior and
#'   prior probability that at least one lag of the cause is included; under
#'   the spike-and-slab prior the prior value, `1 - (1 - inclusion_prob)^lag`,
#'   grows with the order, which is why this probability is not the one that
#'   decides), `inclusion` (posterior inclusion probability of every cause lag, spike
#'   and slab only), `geweke` (posterior mean and credible interval of the
#'   measure), `coefficients` (posterior mean and interval of the cause's
#'   coefficients), `draws` (matrix of posterior draws, the chains stacked),
#'   `ess` (effective sample size of the measure), `lambda` (posterior mean
#'   and mode of the tightness, Minnesota prior only), `rhat` (Gelman-Rubin
#'   potential scale reduction factor of the measure over the chains; values
#'   close to 1 indicate that the chains agree), `n`, `settings`. Method:
#'   `print`.
#' @references
#' George, E. I. and McCulloch, R. E. (1993). Variable selection via Gibbs sampling. *Journal of the American Statistical Association*, 88(423), 881-889. \doi{10.1080/01621459.1993.10476353}
#'
#' Geweke, J. (1982). Measurement of linear dependence and feedback between multiple time series. *Journal of the American Statistical Association*, 77(378), 304-313. \doi{10.1080/01621459.1982.10477803}
#' @export
#' @examples
#' # Does x Granger-cause y?  With the spike-and-slab prior every lag of x may
#' # enter the regression of y or not; BF10 compares "at least one lag" with
#' # "no lag", and the posterior inclusion probabilities say which lags carry
#' # the link.
#' d <- sim_var(400, "x_to_y", seed = 1)
#' b <- bayes_granger(d$y, d$x, lag = 2, control = bayes_control(n_iter = 300, chains = 2), seed = 1)
#' b
#' b$inclusion    # posterior inclusion probability of each lag of x
#' b$geweke       # Geweke's measure: posterior mean and credible interval
#' \donttest{
#' # the reverse direction: BF10 below 1 (evidence against y -> x)
#' bayes_granger(d$x, d$y, lag = 2, control = bayes_control(n_iter = 300, chains = 2), seed = 1)
#' # the closed-form g-prior Bayes factor needs no sampling
#' bayes_granger(d$y, d$x, lag = 2, control = bayes_control(prior = "g"))
#' # more chains in parallel (see ?Rcausality-parallel); R-hat checks their agreement
#' bayes_granger(d$y, d$x, lag = 2, control = bayes_control(n_iter = 1000, chains = 4),
#'               seed = 1, cores = 2)
#' }
bayes_granger <- function(effect, cause, controls = NULL, lag = 1L, control = bayes_control(), prior_prob = 0.5,
                          seed = NULL, cores = getOption("Rcausality.cores", 1L)) {
  .local_seed(seed)
  lag <- max(1L, as.integer(lag)); n <- length(effect)
  if (is.null(control$prior)) control <- bayes_control()
  if (identical(control$prior, "minnesota")) return(.bayes_granger_minnesota(effect, cause, controls, lag, control, prior_prob, cores))
  if (control$prior %in% c("g", "hyper_g", "bic")) {
    g <- .granger_pair(effect, cause, controls, lag, prior = control$prior, prior_par = .prior_par(control))
    bf <- g$bf10
    return(structure(list(bf10 = bf, post_prob = bf * prior_prob / (bf * prior_prob + 1 - prior_prob), inclusion = NULL,
                          geweke = c(mean = g$effect_size, lower = NA_real_, upper = NA_real_), coefficients = NULL,
                          draws = NULL, ess = NA_real_, n = g$n_eff,
                          settings = list(prior = control$prior, lag = lag, prior_prob = prior_prob)), class = "bayes_granger"))
  }
  yd <- .std(effect)[(lag + 1):n]
  own <- .lagmat(.std(effect), lag)
  cz <- if (is.null(controls)) NULL else do.call(cbind, lapply(seq_len(ncol(controls)), function(j) .lagmat(.std(controls[, j]), lag)))
  Xc <- .lagmat(.std(cause), lag); colnames(Xc) <- paste0("lag", seq_len(lag))
  X0 <- cbind(1, own, cz); d0 <- NULL                                    # intercept and always-in terms: flat priors
  .bayes_reg_core(yd, X0, d0, Xc, control, prior_prob, cores, lag)
}

# The spike-and-slab / normal-slab comparison of "cause columns Xc in" against
# "out" for the response yd, given the always-in regressors X0 with prior
# variances d0 (times sigma^2).  Used with the lag design of one series pair
# (bayes_granger()) and with the stacked, within-unit standardised design of a
# panel (the pooled Bayesian panel test).
.bayes_reg_core <- function(yd, X0, d0, Xc, control, prior_prob = 0.5, cores = 1L, lag = ncol(Xc)) {
  if (is.null(colnames(Xc))) colnames(Xc) <- paste0("lag", seq_len(ncol(Xc)))
  q <- ncol(Xc); tau2 <- control$scale^2; pin <- control$inclusion_prob
  N <- length(yd); a0 <- 0.01; b0 <- 0.01
  # The always-in regressors (intercept, own lags, lags of the controls) are common to all models and have flat
  # priors, with p(sigma^2) vague: the restricted regression is the base model, as for the g priors.  They are
  # projected out, which leaves nu0 = N - rank(X0) degrees of freedom for the lags of the cause and the error
  # variance.  (A proper normal prior with variance sigma^2 d0 on them, however vague, keeps N in the exponent of the
  # marginal likelihood although the common coefficients have used rank(X0) degrees of freedom: with 9 coefficients
  # and 13 observations the Bayes factor then exceeded 3 for 29% of the series without a link, and with more
  # coefficients than degrees of freedom for 80%.  `d0` is kept in the signature and no longer used.)
  q0 <- qr(X0); nu0 <- N - q0$rank
  if (nu0 < 2L)
    stop(sprintf(paste0("The regression has %d always-included coefficients for %d observations and leaves no degrees of freedom ",
                        "for the lags of the cause; prior = \"minnesota\" shrinks all coefficients and remains defined."),
                 q0$rank, N), call. = FALSE)
  yp <- qr.resid(q0, yd); Xp <- qr.resid(q0, Xc); colnames(Xp) <- colnames(Xc)
  ml <- function(gam) .log_ml_nig(yp, Xp[, gam, drop = FALSE], rep(tau2, sum(gam)), a0, b0, n = nu0)
  n_keep <- ceiling(control$n_iter / control$thin)                    # kept draws per chain
  exact <- control$prior == "normal" || q <= 12L
  if (control$prior == "normal") {
    lbf <- ml(rep(TRUE, q)) - ml(rep(FALSE, q)); bf <- exp(min(lbf, 700))
    post <- bf * prior_prob / (bf * prior_prob + 1 - prior_prob); incl <- NULL
    G <- matrix(TRUE, 1L, q); w <- 1
  } else if (exact) {                                                   # spike and slab, enumerated
    G <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), q)))
    lm_ <- apply(G, 1, ml); lp <- rowSums(G) * log(pin) + (q - rowSums(G)) * log(1 - pin)
    lw <- lm_ + lp; w <- exp(lw - max(lw)); w <- w / sum(w)
    p_none <- w[1]; incl <- colSums(G * w); names(incl) <- colnames(Xc)
    prior_none <- (1 - pin)^q
    bf <- ((1 - p_none) / max(p_none, 1e-300)) / ((1 - prior_none) / prior_none)
    post <- 1 - p_none
  }
  # one chain: the lag-subset indicators (exact draws, or a Gibbs sampler over
  # the indicators when the cause has more than 12 lags), then (b, sigma) given
  # the subset and Geweke's measure for every kept draw
  run_chain <- function(ch) {
    gam_keep <- if (exact) {
      G[sample.int(nrow(G), n_keep, replace = TRUE, prob = w), , drop = FALSE]
    } else {
      gam <- stats::runif(q) < pin; if (!any(gam)) gam[1L] <- TRUE
      cur <- ml(gam); keep <- matrix(FALSE, n_keep, q); kk <- 0L
      for (it in seq_len(control$burnin + control$n_iter)) {
        for (j in seq_len(q)) {
          g1 <- gam; g1[j] <- TRUE; g0 <- gam; g0[j] <- FALSE
          l1 <- if (gam[j]) cur else ml(g1); l0 <- if (gam[j]) ml(g0) else cur
          p1 <- 1 / (1 + exp(l0 - l1 + log(1 - pin) - log(pin)))
          nj <- stats::runif(1) < p1; if (nj != gam[j]) { gam[j] <- nj; cur <- if (nj) l1 else l0 }
        }
        if (it > control$burnin && (it - control$burnin) %% control$thin == 0L && kk < n_keep) { kk <- kk + 1L; keep[kk, ] <- gam }
      }
      keep[seq_len(kk), , drop = FALSE]
    }
    draws <- matrix(0, nrow(gam_keep), q + 2L); colnames(draws) <- c(colnames(Xc), "sigma", "geweke")
    for (s in seq_len(nrow(gam_keep))) {
      # the conjugate posterior of the cause coefficients and of sigma^2 given the subset, on the projected data
      gam <- gam_keep[s, ]; XS <- Xp[, gam, drop = FALSE]; bc <- numeric(q); a <- a0 + nu0 / 2
      if (any(gam)) {
        A <- crossprod(XS) + diag(1 / tau2, ncol(XS)); V <- solve(A); Xy <- crossprod(XS, yp); m <- V %*% Xy
        sig2 <- 1 / stats::rgamma(1, a, rate = max(b0 + (sum(yp^2) - sum(m * Xy)) / 2, 1e-12))
        bc[gam] <- as.numeric(m + t(chol(V * sig2)) %*% stats::rnorm(ncol(XS)))
        contrib <- mean(as.numeric(XS %*% bc[gam])^2)              # variance explained beyond the always-in terms
      } else { sig2 <- 1 / stats::rgamma(1, a, rate = max(b0 + sum(yp^2) / 2, 1e-12)); contrib <- 0 }
      draws[s, ] <- c(bc, sqrt(sig2), log(1 + contrib / sig2))
    }
    list(draws = draws, gamma = gam_keep)
  }
  chains <- .plapply(seq_len(control$chains), run_chain, cores = cores)
  draws <- do.call(rbind, lapply(chains, `[[`, "draws"))
  if (!exact) {                                                         # Gibbs: evidence from the pooled indicator draws
    keep <- do.call(rbind, lapply(chains, `[[`, "gamma"))
    incl <- colMeans(keep); names(incl) <- colnames(Xc)
    p_none <- (sum(!apply(keep, 1, any)) + 0.5) / (nrow(keep) + 1)
    prior_none <- (1 - pin)^q
    bf <- ((1 - p_none) / p_none) / ((1 - prior_none) / prior_none); post <- 1 - p_none
  }
  ci <- function(v) c(mean = mean(v), lower = unname(stats::quantile(v, (1 - control$credible) / 2)),
                      upper = unname(stats::quantile(v, 1 - (1 - control$credible) / 2)))
  ess <- {
    v <- draws[, "geweke"]; nd <- nrow(draws); if (stats::sd(v) == 0) nd else {
      ac <- stats::acf(v, plot = FALSE, lag.max = min(50, nd - 1))$acf[-1]
      first_neg <- which(ac < 0)[1]; if (!is.na(first_neg)) ac <- ac[seq_len(first_neg - 1)]
      nd / (1 + 2 * sum(ac)) } }
  rhat <- .gelman_rhat(lapply(chains, function(ch) ch$draws[, "geweke"]))
  # P(H1 | data) from the Bayes factor and the prior probability of the link, as for every other prior.  Under the
  # prior of the indicators "some lag is included" has the prior probability 1 - (1 - inclusion_prob)^lag, which
  # grows with the order (0.996 at lag 8); the posterior of that event is reported separately (`any_included`).
  any_included <- c(posterior = post, prior = if (control$prior == "normal") prior_prob else 1 - (1 - pin)^q)
  post <- bf * prior_prob / (bf * prior_prob + 1 - prior_prob)
  structure(list(bf10 = bf, post_prob = post, any_included = any_included, inclusion = incl, geweke = ci(draws[, "geweke"]),
                 coefficients = t(apply(draws[, seq_len(q), drop = FALSE], 2, ci)), draws = draws, ess = ess, rhat = rhat, n = N,
                 settings = list(prior = control$prior, lag = lag, scale = control$scale, inclusion_prob = pin,
                                 prior_prob = prior_prob,
                                 chains = control$chains, n_iter = control$n_iter, burnin = control$burnin, thin = control$thin,
                                 credible = control$credible, sampler = if (exact) "exact" else "Gibbs")), class = "bayes_granger")
}

# Gelman-Rubin potential scale reduction factor of one scalar quantity over
# the chains (NA with a single chain or constant draws)
.gelman_rhat <- function(chains) {
  m <- length(chains); n <- min(lengths(chains))
  if (m < 2L || n < 2L) return(NA_real_)
  M <- vapply(chains, function(v) v[seq_len(n)], numeric(n))
  W <- mean(apply(M, 2, stats::var)); B <- n * stats::var(colMeans(M))
  if (!is.finite(W) || W <= 0) return(NA_real_)
  sqrt(((n - 1) / n * W + B / n) / W)
}

#' @export
print.bayes_granger <- function(x, digits = 3, ...) {
  s <- x$settings
  .out(sprintf("Bayesian Granger causality: second argument (cause) -> first argument (effect); prior: %s, lag %d, n = %d\n",
               s$prior, s$lag, x$n))
  .out(sprintf("  BF10 = %s (%s); posterior probability of the link = %.3f (prior %.2f)\n", .fmt_bf(x$bf10),
              as.character(bf_evidence(x$bf10)), x$post_prob, s$prior_prob))
  if (!is.null(x$inclusion)) {
    .out("  posterior inclusion probability per cause lag:", paste(sprintf("%s %.2f", names(x$inclusion), x$inclusion), collapse = ", "), "\n")
    if (!is.null(x$any_included))
      .out(sprintf("  probability that some lag is included: %.3f (%.3f a priori, inclusion probability %g per lag)\n",
                   x$any_included[["posterior"]], x$any_included[["prior"]], s$inclusion_prob))
  }
  if (!is.na(x$geweke[["lower"]]))
    .out(sprintf("  Geweke measure log(1 + explained / sigma^2): %.3f [%.3f, %.3f] (%.0f%% credible; ESS %.0f of %d draws%s)\n",
                x$geweke[["mean"]], x$geweke[["lower"]], x$geweke[["upper"]], 100 * s$credible, x$ess, nrow(x$draws),
                if (is.null(x$rhat) || is.na(x$rhat)) "" else sprintf("; %d chains, R-hat %.3f", s$chains, x$rhat)))
  else .out(sprintf("  Geweke measure: %.3f\n", x$geweke[["mean"]]))
  if (!is.null(x$coefficients)) {
    .out("  cause coefficients (standardised series):\n")
    print(round(x$coefficients, digits))
  }
  invisible(x)
}

# ---- (3) prior over the causal cases ------------------------------------------
# The posterior / score vector `p` of the five cases was computed under the
# prior `base` of the method (NULL = uniform over the five cases): replace
# that prior by `case_prior`, i.e. posterior ~ p * case_prior / base.
# Dimensional causality: uniform over the nine elementary events of Benko et
# al. (2024), which gives the cases the weights (1, 1, 1, 3, 3) / 9.
.dc_case_base <- c(x_to_y = 1, y_to_x = 1, bidirectional = 1, common_cause = 3, independent = 3) / 9
.bayes_cases <- function(p, case_prior, base = NULL) {
  if (is.null(case_prior)) return(p)
  q <- p * case_prior[names(p)]
  if (!is.null(base)) q <- q / base[names(p)]
  q / sum(q)
}

# Bayesian Markov decision table: the six binary events (x Markov, y Markov,
# joint Markov, x->y, y->x, contemporaneous) get posterior probabilities from
# their surrogate Bayes factors; the 64 combinations are weighted by their
# (independent) probabilities and mapped through the decision table of
# Stippinger et al. (2023) to posterior probabilities of the two directions.
.bayes_markov <- function(pM, pD, iid, case_prior = NULL) {
  ev <- c("MA", "MB", "MJ", "sAB", "sBA", "sI0"); pr <- c(pM, pD)
  G <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), 6))); colnames(G) <- ev
  w <- apply(G, 1, function(g) prod(ifelse(g, pr, 1 - pr)))
  dec <- lapply(seq_len(nrow(G)), function(i) .markov_decide(G[i, "MA"], G[i, "MB"], G[i, "MJ"], G[i, "sAB"], G[i, "sBA"], G[i, "sI0"], iid))
  ab <- vapply(dec, function(d) d$ab != "n", TRUE); ba <- vapply(dec, function(d) d$ba != "n", TRUE)
  cls <- ifelse(ab & ba, "bidirectional", ifelse(ab, "x_to_y", ifelse(ba, "y_to_x", ifelse(G[, "sI0"], "common_cause", "independent"))))
  lev <- c("x_to_y", "y_to_x", "bidirectional", "common_cause", "independent")
  pc <- tapply(w, factor(cls, lev), sum); pc[is.na(pc)] <- 0
  if (!is.null(case_prior)) {
    # the prior of the cases that the table implies when every event has probability 1/2
    base <- tapply(rep(1 / nrow(G), nrow(G)), factor(cls, lev), sum); base[is.na(base)] <- 0
    pc <- .bayes_cases(pc, case_prior, pmax(base, 1e-12))
  }
  list(p_ab = sum(pc[c("x_to_y", "bidirectional")]), p_ba = sum(pc[c("y_to_x", "bidirectional")]), cases = pc)
}

# ---- post-processing of a causres under bayes = TRUE --------------------------
.bayes_postprocess <- function(r, control, prior_prob, bf_threshold, prob_threshold, seed = NULL) {
  s <- r$settings; P <- r$pairs; key <- paste(P$cause, "->", P$effect)
  rows <- list()
  bf_row <- function(k, bf, post, eff, how, bound = NA_real_)
    data.frame(pair = k, bf10 = bf, post_prob = post, effect_mean = eff[["mean"]], effect_lower = eff[["lower"]],
               effect_upper = eff[["upper"]], model = how, bf_bound = bound, stringsAsFactors = FALSE)
  na3 <- c(mean = NA_real_, lower = NA_real_, upper = NA_real_)
  if (s$method %in% c("dc", "dfe")) {
    # posterior over the five cases with the case prior; directed posterior = P(x->y) + P(x<->y)
    tab <- if (s$method == "dc") r$dc else r$dfe
    if (!is.null(tab) && !is.null(control$case_prior)) {
      for (i in seq_len(nrow(tab))) {
        p <- unlist(tab[i, c("x_to_y", "y_to_x", "bidirectional", "common_cause", "independent")])
        q <- .bayes_cases(p, control$case_prior, if (s$method == "dc") .dc_case_base else NULL); tab[i, names(q)] <- q
        a <- tab$var1[i]; b <- tab$var2[i]; prior_dir <- sum(control$case_prior[c("x_to_y", "bidirectional")])
        for (dir in c("ab", "ba")) {
          kk <- if (dir == "ab") paste(a, "->", b) else paste(b, "->", a)
          post <- if (dir == "ab") q[["x_to_y"]] + q[["bidirectional"]] else q[["y_to_x"]] + q[["bidirectional"]]
          pri <- if (dir == "ab") prior_dir else sum(control$case_prior[c("y_to_x", "bidirectional")])
          j <- match(kk, key); if (is.na(j)) next
          P$post_prob[j] <- post; P$bf10[j] <- (post / (1 - post)) / (pri / (1 - pri)); P$statistic[j] <- post
          rows[[length(rows) + 1L]] <- bf_row(kk, P$bf10[j], post, na3, "case prior")
        }
      }
      if (s$method == "dc") r$dc <- tab else r$dfe <- tab
    } else for (j in seq_along(key)) rows[[length(rows) + 1L]] <- bf_row(key[j], P$bf10[j], P$post_prob[j], na3, "posterior of the cases")
    P$evidence <- as.character(bf_evidence(P$bf10))
    P$significant <- if (control$decision == "bf") P$bf10 > bf_threshold else P$post_prob > prob_threshold
    P$significant[is.na(P$significant)] <- FALSE
  } else if (s$method == "markov" && !is.null(r$null_markov)) {
    M <- r$markov
    for (i in seq_len(nrow(M))) {
      nm <- r$null_markov[[paste(M$var1[i], "~", M$var2[i])]]
      bfs <- lapply(nm, function(e) surrogate_bf(e$statistic, e$null, alternative = e$alternative, effect_prior = control$effect_prior,
                                                 scale = control$scale, null_model = control$null_model, prior_prob = 0.5,
                                                 n_draws = control$n_draws, credible = control$credible,
                                                 effect_scale = .eff_scale(control), n_eff = .eff_n(NA, s$n)))
      # Markov tests: small statistic = not Markov; the event "Markov" has posterior 1 - P(effect)
      pM <- c(MA = 1 - bfs$x$post_prob, MB = 1 - bfs$y$post_prob, MJ = 1 - bfs$j$post_prob)
      pD <- c(sAB = bfs$xy$post_prob, sBA = bfs$yx$post_prob, sI0 = bfs$i0$post_prob)
      bm <- .bayes_markov(pM, pD, if (is.null(r$settings$iid)) "allowed" else r$settings$iid, control$case_prior)
      M$p_x_markov[i] <- pM[["MA"]]; M$p_y_markov[i] <- pM[["MB"]]; M$p_joint_markov[i] <- pM[["MJ"]]
      for (dir in c("ab", "ba")) {
        kk <- if (dir == "ab") paste(M$var1[i], "->", M$var2[i]) else paste(M$var2[i], "->", M$var1[i])
        post <- if (dir == "ab") bm$p_ab else bm$p_ba; j <- match(kk, key); if (is.na(j)) next
        pri <- if (is.null(control$case_prior)) 0.5 else
          sum(control$case_prior[c(if (dir == "ab") "x_to_y" else "y_to_x", "bidirectional")])
        pri <- min(max(pri, 1e-12), 1 - 1e-12)
        P$post_prob[j] <- post; P$bf10[j] <- min(1e12, (post / max(1 - post, 1e-12)) / (pri / (1 - pri)))
        e <- if (dir == "ab") bfs$xy$effect else bfs$yx$effect
        rows[[length(rows) + 1L]] <- bf_row(kk, P$bf10[j], post, e, "Bayesian decision table")
      }
      if (!is.null(r$instant)) {
        ii <- which(r$instant$var1 == M$var1[i] & r$instant$var2 == M$var2[i])
        if (length(ii)) { r$instant$post_prob[ii] <- bfs$i0$post_prob; r$instant$significant[ii] <- bfs$i0$bf10 > bf_threshold }
      }
    }
    r$markov <- M
    P$evidence <- as.character(bf_evidence(P$bf10))
    P$significant <- if (control$decision == "bf") P$bf10 > bf_threshold else P$post_prob > prob_threshold
    P$significant[is.na(P$significant)] <- FALSE
  } else if (!is.null(r$null)) {
    conv <- if (!is.null(r$ccm)) r$ccm$converging[match(key, paste(r$ccm$cause, "->", r$ccm$effect))] else rep(TRUE, nrow(P))
    for (j in seq_along(key)) {
      e <- r$null[[key[j]]]; if (is.null(e)) next
      if (.determined(e$cond_r2)) {                     # no test for a cause that the conditioning past determines
        P$bf10[j] <- 1; P$post_prob[j] <- prior_prob
        rows[[length(rows) + 1L]] <- bf_row(key[j], 1, prior_prob, na3, "not tested (cause determined by the conditioning past)")
        next
      }
      b <- surrogate_bf(e$statistic, e$null, alternative = e$alternative, effect_prior = control$effect_prior,
                        scale = control$scale, null_model = control$null_model, prior_prob = prior_prob,
                        n_draws = control$n_draws, credible = control$credible,
                        effect_scale = .eff_scale(control), n_eff = .eff_n(P$n_eff[j], s$n))
      P$bf10[j] <- b$bf10; P$post_prob[j] <- b$post_prob
      rows[[length(rows) + 1L]] <- bf_row(key[j], b$bf10, b$post_prob, b$effect, sprintf("surrogate BF (%s)", control$effect_prior),
                                          bound = if (isTRUE(b$extrapolated)) b$bound else NA_real_)
    }
    P$evidence <- as.character(bf_evidence(P$bf10))
    sig <- if (control$decision == "bf") P$bf10 > bf_threshold else P$post_prob > prob_threshold
    # the method-specific side conditions stay: CCM convergence, the coverage
    # index of the recurrence criterion, the reconstruction index of CCDH
    gate <- conv %in% TRUE
    if (s$method == "hirata" && !is.null(s$theta)) gate <- gate & P$statistic >= s$theta
    if (s$method == "ccdh" && !is.null(s$prob_threshold)) gate <- gate & (1 - P$statistic) <= s$prob_threshold
    if (s$method == "rf_granger") gate <- gate & P$statistic > 0          # an improvement of the prediction
    P$significant <- sig %in% TRUE & gate
    if (!is.null(r$instant) && !is.null(r$null_instant)) {
      for (i in seq_len(nrow(r$instant))) {
        e <- r$null_instant[[paste(r$instant$var1[i], "~", r$instant$var2[i])]]; if (is.null(e)) next
        b <- surrogate_bf(e$statistic, e$null, alternative = "greater", effect_prior = control$effect_prior, scale = control$scale,
                          null_model = control$null_model, prior_prob = prior_prob, n_draws = control$n_draws,
                          effect_scale = .eff_scale(control), n_eff = .eff_n(NA, s$n))
        r$instant$bf10[i] <- b$bf10; r$instant$post_prob[i] <- b$post_prob
        r$instant$significant[i] <- if (control$decision == "bf") b$bf10 > bf_threshold else b$post_prob > prob_threshold
      }
    }
  } else if (s$method == "toda_yamamoto") {
    # closed-form Bayes factor of the lag-augmented regressions (BIC, g or hyper-g)
    for (j in seq_along(key)) rows[[length(rows) + 1L]] <- bf_row(key[j], P$bf10[j], P$post_prob[j],
                                                                  c(mean = P$effect_size[j], lower = NA, upper = NA),
                                                                  "lag-augmented regression")
    P$significant <- if (control$decision == "bf") P$bf10 > bf_threshold else P$post_prob > prob_threshold
    P$significant[is.na(P$significant)] <- FALSE
  } else if (s$method == "pcmci") {
    # partial correlations: Fisher z against its sampling sd, Savage-Dickey with the effect prior
    for (j in seq_along(key)) {
      rho <- P$statistic[j]; ne <- P$n_eff[j]; if (is.na(rho) || is.na(ne) || ne < 10) next
      z <- atanh(max(min(rho, 0.999), -0.999)); sdz <- 1 / sqrt(ne - 3)
      b <- surrogate_bf(z, list(mean = 0, sd = sdz), alternative = "two.sided", effect_prior = "normal", scale = control$scale,
                        prior_prob = prior_prob, n_draws = control$n_draws, credible = control$credible,
                        effect_scale = .eff_scale(control), n_eff = ne)
      P$bf10[j] <- b$bf10; P$post_prob[j] <- b$post_prob
      rows[[length(rows) + 1L]] <- bf_row(key[j], b$bf10, b$post_prob, tanh(b$effect), "Fisher-z normal approximation")
    }
    P$evidence <- as.character(bf_evidence(P$bf10))
    P$significant <- if (control$decision == "bf") P$bf10 > bf_threshold else P$post_prob > prob_threshold
    P$significant[is.na(P$significant)] <- FALSE
  }
  r <- .rebuild_causres(r, P)
  r$bayes <- if (length(rows)) do.call(rbind, rows) else NULL
  r$settings$bayes <- TRUE; r$settings$bayes_control <- control; r$settings$decision <- control$decision
  r$settings$bf_threshold <- bf_threshold; r$settings$prob_threshold <- prob_threshold
  r
}

# Bayesian Granger regression for every directed pair of a Granger run: the
# F-test run supplies the selected lags; the sampler supplies the evidence
.bayes_granger_pairs <- function(r, control, prior_prob, bf_threshold, prob_threshold, conditional, cores = 1L,
                                 lag = "auto") {
  P <- r$pairs; D <- r$data; rows <- list()
  # the Minnesota posterior is proper for any lag order: a requested lag is
  # used as it is, even where the F-test run had to shorten it
  if (identical(control$prior, "minnesota") && !identical(lag, "auto")) P$lag <- rep(max(1L, as.integer(lag)), nrow(P))
  sc <- .split_cores(nrow(P), cores)                 # pairs in parallel, or the chains of each pair
  res <- .plapply(seq_len(nrow(P)), function(j) {
    cs <- P$cause[j]; ef <- P$effect[j]
    ctrl_names <- if (conditional) setdiff(names(D), c(cs, ef)) else character(0)
    bayes_granger(D[[ef]], D[[cs]], if (length(ctrl_names)) as.matrix(D[, ctrl_names, drop = FALSE]) else NULL,
                  lag = P$lag[j], control = control, prior_prob = prior_prob, cores = sc$inner)
  }, cores = sc$outer)
  for (j in seq_len(nrow(P))) {
    b <- res[[j]]
    if (identical(lag, "first") && !is.null(r$by_lag)) {
      # the order was searched: the Bayes factor of the selected order divided by the number of orders scanned is a
      # lower bound of the Bayes factor of "a link at one of the orders" (their average)
      K <- max(1L, sum(r$by_lag$cause == P$cause[j] & r$by_lag$effect == P$effect[j]))
      b$bf10 <- b$bf10 / K; b$post_prob <- b$bf10 * prior_prob / (b$bf10 * prior_prob + 1 - prior_prob)
    }
    P$bf10[j] <- b$bf10; P$post_prob[j] <- b$post_prob
    rows[[j]] <- data.frame(pair = paste(P$cause[j], "->", P$effect[j]), bf10 = b$bf10, post_prob = b$post_prob,
                            effect_mean = b$geweke[["mean"]], effect_lower = b$geweke[["lower"]], effect_upper = b$geweke[["upper"]],
                            ess = b$ess, rhat = if (is.null(b$rhat)) NA_real_ else b$rhat,
                            model = sprintf("Bayesian regression (%s)", control$prior), stringsAsFactors = FALSE)
  }
  P$evidence <- as.character(bf_evidence(P$bf10))
  P$significant <- if (control$decision == "bf") P$bf10 > bf_threshold else P$post_prob > prob_threshold
  P$significant[is.na(P$significant)] <- FALSE
  out <- .rebuild_causres(r, P)
  out$bayes <- do.call(rbind, rows); out$bayes_fits <- stats::setNames(res, out$bayes$pair)
  out$settings$method <- "bayes_granger"; out$settings$bf_threshold <- bf_threshold; out$settings$prob_threshold <- prob_threshold
  out$settings$bayes <- TRUE; out$settings$prior <- control$prior; out$settings$bayes_control <- control
  out
}
