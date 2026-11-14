#' Detect causal relations between time series or in panel data
#'
#' Single entry point of the package. The causality method, its Bayesian
#' variant and parallel execution are options of the same call,
#' `rcausality(x, method = , bayes = , cores = )`, for every kind of input:
#' * two series: `rcausality(x, y)`;
#' * one data frame / matrix: `rcausality(df)` -- all ordered pairs;
#' * two blocks (MIMO): `rcausality(X, Y)` with several columns each --
#'   every column of `X` is tested as a cause of every column of `Y`;
#' * one explained variable: `rcausality(df, target = "y")` -- which series
#'   drive `y`, at which lag, directly or through another series, and how the
#'   explanatory series drive each other (see [causal_drivers()]);
#' * panel data in long format: `rcausality(df, id = "unit", time = "year")`
#'   -- the same analysis in every unit, at the group level and, where the
#'   method allows, pooled over the units (see [panel_causality()]).
#'
#' @section Methods:
#' [causal_methods()] lists the methods; `causal_methods("<name>")` shows the
#' settings of one of them.
#' * Granger family: `"granger"` (F-test, default), `"bayes_granger"` (the
#'   Bayesian Granger regression, i.e. `method = "granger", bayes = TRUE`),
#'   `"toda_yamamoto"` (lag-augmented test for integrated series; setting
#'   `d_max`), `"diks_panchenko"` (nonparametric, nonlinear; settings `eps`,
#'   `C`, `filter`), `"rf_granger"` (random-forest Granger-like test;
#'   settings `num_trees`, `mtry`, `min_node_size`, `folds`, `shadow`,
#'   `importance`; needs \pkg{ranger}).
#' * Information-theoretic: `"cmi"` (kNN conditional mutual information
#'   `I(x[t-lag]; y[t] | history of y, controls)` with a surrogate test
#'   by conditional shifts of the cause, see below), `"cmi_lag"` (Zlatniczki et al. 2021: also conditions on the
#'   cause's own past and adds a contemporaneous test whose significance,
#'   without lagged links, indicates dependence within the sampling interval:
#'   a common shock or an instantaneous relation; class `common_cause`), `"te"` (binned
#'   transfer entropy, circular-shift surrogates), `"pid"` (MAR-PID: unique
#'   + synergistic information of the cause's past).
#' * State-space and dynamical-systems methods: `"ccm"` (convergent cross
#'   mapping), `"dc"` (Bayesian dimensional causality), `"ccdh"`
#'   (complexity-based causal discovery), `"hirata"` (recurrence plots),
#'   `"cd"` (correlation dimensions).
#' * Markov-chain and order-based methods: `"markov"` (Stippinger et al.
#'   2023), `"df"` (degrees of freedom, Telcs et al. 2025), `"dfe"`
#'   (df-entropy with a trained random forest).
#' * `"pcmci"` (PCMCI of Runge et al. 2019 through the Python package
#'   tigramite; always multivariate).
#'
#' @section Method settings:
#' Method-specific settings go in `method_args` (a named list; unknown names
#' are an error). The settings of the methods of version 0.0.1 may also be
#' given directly in `...`, as before.
#' * `"cmi"`, `"cmi_lag"`: `k` (neighbours of the estimator, default 5),
#'   `history` (order of the conditioning history of the effect; `"auto"` =
#'   the VAR order of the pair, at most 5; too short a history lets the
#'   cause's lags stand in for the effect's own deeper past and produces
#'   spurious reverse links), `reference` (the reference distribution:
#'   `"shift"`, the default, uses conditional shifts: the lags of the cause
#'   are written as location plus scale times a standardized residual, with
#'   location and scale estimated as functions of the conditioning vector
#'   (see the section on conditional shifts), and the residuals are shifted
#'   circularly, which keeps the serial dependence of the cause;
#'   `"local"` uses the local permutations of Runge (2018) among the
#'   `k_perm` nearest neighbours in the conditioning space, which do not
#'   keep it and reject too often for strongly autocorrelated series). With
#'   `lag = "auto"` the lags `1..max_lag` are tested at once (with
#'   `reference = "shift"` the scan is limited as for `"te"`): the reported
#'   lag is the one
#'   with the largest statistic and the p-value is that of the maximal
#'   statistic against the maximum of the replicated statistics over the
#'   lags (the same max-statistic test is used by `"te"`, `"pid"` and
#'   `"df"`). Series are standardised.
#' * `"te"` (see [transfer_entropy()]): `bins` (quantile bins), `history`
#'   (the effect's own history); the cause's lag is scanned. The circular
#'   shifts of a scan over `L` lags are spaced by `L`, which leaves about
#'   `0.8 n / L` replicates; the scan is limited to the lags for which
#'   enough of them remain (at least 40, and as many as the smallest
#'   p-value needs to pass a Holm or Bonferroni correction over the tests of
#'   the analysis; a message reports the limit). A numeric `lag` tests one
#'   lag with all shifts. `reference`: `"circular"` (default), circular
#'   shifts of the coded cause, which generate the hypothesis that the two
#'   series are independent; `"conditional"`, conditional shifts of the cause
#'   given the past of the effect, coded after the shift, for series on a
#'   continuous scale. In the reverse direction of a strong link whose
#'   driver is persistent the past of the cause is a function of the past
#'   of the effect, which a few symbols do not resolve: with circular
#'   shifts the test then rejects too often in long series (a warning names
#'   the situation), with conditional shifts it does not.
#' * `"ccm"` (see [ccm()]): `m`, `tau` (embedding); a link needs a
#'   significant and converging cross-map skill.
#' * `"pid"` (see [marpid()]): `k_target`, `k_source` (binarisation depths),
#'   `history`. `reference`: `"permutation"` (default), the cause is permuted
#'   within the levels of the coded own past of the effect, an exact test
#'   for the coded variables; `"conditional"`, conditional shifts of the
#'   cause given the past of the effect, coded after the shift, for series
#'   on a continuous scale. The code does not resolve a persistent past:
#'   in the reverse direction of a strong link whose driver is persistent
#'   the permutation test rejects too often, the more so the longer the
#'   series (a warning names the situation); the test with conditional
#'   shifts does not.
#' * `"dc"` (see [dimensional_causality()]): `m`, `tau`, `k_range`, `prep`;
#'   decisions by posterior probabilities (`prob_threshold`).
#' * `"hirata"` (see [hirata_rp()]): `m`, `tau`, `rate` (recurrence rate),
#'   `theta` (coverage-index threshold). `"cd"` (see [krakovska_cd()]): `m`,
#'   `tau`; bootstrap z-tests at level `alpha`. Both give verdicts per pair.
#' * `"ccdh"` (see [ccdh()]): `m`, `tau`, `K` (NSRPS steps), `n_surr`
#'   (circular-shift surrogates), `theta` (reconstruction-index threshold);
#'   a numeric `lag` aligns the pair.
#' * `"dfe"` (see [dfe_classify()]): `model`, a trained `dfe_model` from
#'   [dfe_train()] (`NULL` trains and caches a default one).
#' * `"df"`: `k` (default 10), `k_perm`, `history`, `max_constr` (largest
#'   number of constrained past lags). The orders (`$df`) give the Table-2
#'   reading of the paper (including the common-driver case), which decides
#'   a pair whenever it is determinate; the other pairs are decided by the
#'   directed variance-reduction tests (`$pairs`).
#' * `"markov"`: `memory` (`"auto"` from the PACF), `iid_injection`
#'   (`"allowed"` or `"forbidden"` in the decision table), `k_nn`, `n_bins`,
#'   `n_markov`, `k`, `k_perm`. Lags are fixed to 1.
#' * `"toda_yamamoto"`: `d_max` (default 1), the number of extra lags of the
#'   lag-augmented VAR (the maximal order of integration); the Wald (F) test
#'   concerns the first `lag` lags of the cause, so no unit-root pre-test is
#'   needed.
#' * `"diks_panchenko"` (see [diks_panchenko()]): `eps` (bandwidth in
#'   standard deviations, `NULL` = `min(1.5, C n^(-2/7))`), `C` (default 8),
#'   `filter` (`"none"` or `"var"`: only the nonlinear part is tested),
#'   `reference` (`"conditional"`, conditional shifts of the cause, the
#'   default; `"normal"`, the asymptotic test of Diks and Panchenko 2006,
#'   which rejects too often when the past of the cause depends on the past
#'   of an autocorrelated effect, as in the reverse direction of a link); a
#'   numeric `lag` is the length of the pasts, `"auto"` uses 1.
#' * `"rf_granger"` (see [rf_granger()]): `num_trees`, `mtry`,
#'   `min_node_size`, `folds`, `shadow`, `importance`; `n_perm` is the number
#'   of surrogates (conditional circular shifts).
#' * `"granger"`, `"toda_yamamoto"` (time series): `vcov`, `"classical"`
#'   (default), `"HC3"` or `"HAC"`. `"HC3"` is the Wald test with the
#'   heteroskedasticity-consistent covariance of MacKinnon and White (1985),
#'   divided by the number of restrictions and referred to the F
#'   distribution; `"HAC"` adds the autocovariances of the scores with
#'   Bartlett weights up to lag `floor(4 (n / 100)^(2/9))` (Newey and West
#'   1987), for regressions whose errors are not martingale differences, as
#'   those of a linear regression on a nonlinear process (it rejects a true
#'   null hypothesis somewhat too often in short series). The regressions
#'   are tested for heteroskedasticity
#'   (White's test in the form of Koenker 1981; `$extra$het_p`), and a
#'   warning names these options when the test rejects at the 0.1% level:
#'   the classical F test can then reject without a linear link, for
#'   example when the dependence is nonlinear.
#'
#' @section Conditional shifts:
#' The tests of `"cmi"`, `"cmi_lag"`, `"rf_granger"` and
#' `"diks_panchenko"` (and of `"te"` and `"pid"` with `reference =
#' "conditional"`) refer
#' their statistic to surrogates of the cause that keep its location and scale
#' given the conditioning vector. Every
#' lag of the cause is written as `g(w) + s(w) e`, with the conditioning
#' vector `w` (the past of the effect and of the controls), a location `g`,
#' a scale `s` and a standardized residual `e`; a surrogate keeps `g` and
#' `s` and shifts `e` circularly. Location and scale are estimated flexibly
#' (the polynomial basis in `w` with the smallest leave-one-out error, plus
#' a nearest-neighbour mean of its residuals where that predicts them out of
#' sample; an out-of-bag forest for `"rf_granger"`; for the scale the
#' absolute values of `w` are among the candidate bases, since a spread
#' often depends on the size of a variable and not on its sign, and the
#' fit is weighted by the inverse squared scale), so that a cause that is
#' a nonlinear function of a control, or whose spread a control drives, has
#' surrogates that depend on the control as it does. A dependence of the
#' shape of the cause on `w` beyond location and scale is not reproduced (no
#' test of conditional independence is valid for all null distributions;
#' Shah and Peters 2020).
#'
#' When `w` explains 99% or more of the variance of a lag of the cause, no
#' surrogate differs from the cause by more than the error of the fit: no
#' test is carried out (p-value 1, Bayes factor 1), and a warning says so.
#' This is the case of redundant series, of series with a common trend, and
#' of noise-free deterministic systems, in which the driver is a function of
#' the past of the driven series (Takens 1981). For deterministic systems
#' `"te"` and the state-space methods (`"ccm"`, `"dc"`, `"cd"`, `"ccdh"`,
#' `"hirata"`) are the tools.
#'
#' @section Bayesian variants:
#' With `bayes = TRUE` (settings in `control`, see [bayes_control()])
#' `"granger"` becomes a Bayesian regression: closed-form Bayes factors
#' (`prior = "bic"`, `"g"`, `"hyper_g"`), the hierarchical Minnesota prior
#' for short series (`"minnesota"`) or the spike-and-slab / normal-slab
#' sampler of [bayes_granger()]. The surrogate / permutation methods compare
#' the observed statistic with the surrogate sample through [surrogate_bf()];
#' `"markov"` combines the surrogate Bayes factors of its tests in a Bayesian
#' decision table; `"dc"` and `"dfe"` apply the prior over the causal cases;
#' `"pcmci"` uses a Savage-Dickey Bayes factor. Without `bayes = TRUE` the
#' regression tests are accompanied by the Bayes factor of the
#' unit-information prior (`prior = "g"` with its default scale; columns
#' `bf10`, `post_prob`), which does not enter the decision. `bf10`, `post_prob` and
#' `evidence` are filled for every directed test, the details are in
#' `$bayes`, and a link is declared when `bf10 > bf_threshold` (or, with
#' `control$decision = "posterior"`, when `post_prob > prob_threshold`).
#'
#' @param x numeric vector, matrix or data frame (the putative causes); for
#'   panel data a data frame in long format (see `id`).
#' @param y optional numeric vector, matrix or data frame (the putative
#'   effects). If `NULL` and `x` has several columns, all ordered pairs of the
#'   columns of `x` are tested.
#' @param method the causality method (see Methods and [causal_methods()]).
#' @param bayes logical; the Bayesian variant of the method (see Bayesian
#'   variants).
#' @param cores number of CPU cores (`1` = serial, an integer, `"auto"` = all
#'   physical cores but one, or a cluster from `parallel::makeCluster()`);
#'   default `getOption("Rcausality.cores", 1L)`. The work (pairs, panel
#'   units, permutations, surrogates, chains) is distributed over the cores;
#'   every unit of work has its own random stream, so the results are
#'   *identical* for any number of cores given the `seed` (see
#'   [Rcausality-parallel]).
#' @param lag integer lag order, `"auto"` (selected by `lag_criterion` on
#'   the VAR of the involved series; the lag with the largest statistic for
#'   the methods that scan the lags) or `"first"` (Granger family: the
#'   smallest order at which the test is significant at level
#'   `alpha / (number of orders)`, i.e. the delay after which the cause first
#'   helps; the reported p-value is the Bonferroni-corrected minimum over the
#'   orders `1..max_lag`; for Bayes-factor decisions the first order with
#'   `bf10 > bf_threshold`, and the reported Bayes factor is that of "a link
#'   at one of the orders", the average over the orders). In panels `"auto"`
#'   chooses the order of every link by the BIC of its pooled regression.
#'   `NULL` is `"auto"`, except in the time-series analysis of a `target`
#'   with a method of the Granger family, where it is `"first"`.
#' @param max_lag largest lag order considered when `lag` is not a number.
#' @param lag_criterion `"BIC"`, `"AIC"` or `"HQ"`.
#' @param conditional logical; if `TRUE` every test conditions on all other
#'   series in the data (full conditional / multivariate analysis, for the
#'   methods that support it), otherwise tests are pairwise.
#' @param directions for two series / two blocks: `"both"` (default) tests
#'   x -> y and y -> x; `"forward"` only the causes-in-`x`, effects-in-`y`
#'   direction. Ignored when `y` is `NULL` (all ordered pairs are tested).
#' @param alpha significance level (applied to adjusted p-values).
#' @param adjust multiple-testing correction passed to
#'   [stats::p.adjust()] (default `"holm"`; `"none"` disables it).
#' @param bf_threshold Bayes factor above which a link is declared in the
#'   Bayesian variants.
#' @param prob_threshold posterior probability above which a link is
#'   declared by the methods that decide by posterior probabilities (`"dc"`,
#'   `"dfe"`, and the Bayesian variants with `control$decision =
#'   "posterior"`).
#' @param prior_prob prior probability of a causal link (for the posterior).
#' @param prior shortcut for `control$prior` (the prior of the Granger
#'   regression comparison; see [bayes_control()]); `NULL` = the control's.
#' @param control a [bayes_control()] list with the priors, the Monte Carlo
#'   and MCMC sizes and the decision rule of the Bayesian variants.
#' @param n_perm number of permutations / surrogates of the methods with a
#'   resampling test (raised automatically when it is too small for the
#'   multiplicity correction to reach significance).
#' @param seed optional RNG seed (the caller's generator is left untouched).
#' @param impute logical (default `TRUE`): missing and non-finite values are
#'   replaced by linear interpolation (the ends are carried), and a message
#'   reports how many were replaced in which series; `FALSE` stops with an
#'   error when a value is missing.
#' @param target optional name of one explained variable: the analysis is
#'   organised around it (drivers, lags, importance, direct or mediated
#'   links, explanatory graph); see [causal_drivers()].
#' @param id,time panel data: the names of the unit and (optional) time
#'   columns of the long data frame `x`; see [panel_causality()].
#' @param panel a [panel_control()] list with the settings of the panel
#'   tests.
#' @param method_args named list of method-specific settings (see Method
#'   settings).
#' @param ... the settings of the 0.0.1 methods given directly (e.g.
#'   `bins = 6`), as an alternative to `method_args`.
#' @param m embedding dimension of the state-space methods, as in
#'   `method_args` (an argument of its own only because the name would
#'   otherwise be taken for an abbreviation of `method`).
#' @return An object of class [`causres`][causres-class]; with `target` also
#'   of class `causdrivers`, with `id` of class `causpanel`.
#' @references
#' Granger, C. W. J. (1969). Investigating causal relations by econometric models and cross-spectral methods. *Econometrica*, 37(3), 424-438. \doi{10.2307/1912791}
#'
#' Liang, F., Paulo, R., Molina, G., Clyde, M. A. and Berger, J. O. (2008). Mixtures of g priors for Bayesian variable selection. *Journal of the American Statistical Association*, 103(481), 410-423. \doi{10.1198/016214507000001337}
#'
#' Sellke, T., Bayarri, M. J. and Berger, J. O. (2001). Calibration of p values for testing precise null hypotheses. *The American Statistician*, 55(1), 62-71. \doi{10.1198/000313001300339950}
#'
#' Geweke, J. (1982). Measurement of linear dependence and feedback between multiple time series. *Journal of the American Statistical Association*, 77(378), 304-313. \doi{10.1080/01621459.1982.10477803}
#'
#' Kass, R. E. and Raftery, A. E. (1995). Bayes factors. *Journal of the American Statistical Association*, 90(430), 773-795. \doi{10.1080/01621459.1995.10476572}
#'
#' Schwarz, G. (1978). Estimating the dimension of a model. *The Annals of Statistics*, 6(2), 461-464. \doi{10.1214/aos/1176344136}
#'
#' Holm, S. (1979). A simple sequentially rejective multiple test procedure. *Scandinavian Journal of Statistics*, 6(2), 65-70.
#'
#' Frenzel, S. and Pompe, B. (2007). Partial mutual information for coupling analysis of multivariate time series. *Physical Review Letters*, 99, 204101. \doi{10.1103/PhysRevLett.99.204101}
#'
#' Runge, J. (2018). Conditional independence testing based on a nearest-neighbor estimator of conditional mutual information. *Proceedings of the 21st International Conference on Artificial Intelligence and Statistics (AISTATS)*, PMLR 84, 938-947.
#'
#' Shah, R. D. and Peters, J. (2020). The hardness of conditional independence testing and the generalised covariance measure. *The Annals of Statistics*, 48(3), 1514-1538. \doi{10.1214/19-AOS1857}
#'
#' Takens, F. (1981). Detecting strange attractors in turbulence. In *Dynamical Systems and Turbulence, Warwick 1980*, Lecture Notes in Mathematics 898, 366-381. Springer, Berlin. \doi{10.1007/BFb0091924}
#'
#' MacKinnon, J. G. and White, H. (1985). Some heteroskedasticity-consistent covariance matrix estimators with improved finite sample properties. *Journal of Econometrics*, 29(3), 305-325. \doi{10.1016/0304-4076(85)90158-7}
#'
#' Newey, W. K. and West, K. D. (1987). A simple, positive semi-definite, heteroskedasticity and autocorrelation consistent covariance matrix. *Econometrica*, 55(3), 703-708. \doi{10.2307/1913610}
#'
#' Koenker, R. (1981). A note on studentizing a test for heteroscedasticity. *Journal of Econometrics*, 17(1), 107-112. \doi{10.1016/0304-4076(81)90062-2}
#'
#' Zlatniczki, A., Stippinger, M., Benko, Z., Somogyvari, Z. and Telcs, A. (2021). Relaxation of some confusions about confounders. *Entropy*, 23(11), 1450. \doi{10.3390/e23111450}
#'
#' Stippinger, M., Bencze, A., Zlatniczki, A., Somogyvari, Z. and Telcs, A. (2023). Causal discovery of stochastic dynamical systems: a Markov chain approach. *Mathematics*, 11(4), 852. \doi{10.3390/math11040852}
#'
#' Toda, H. Y. and Yamamoto, T. (1995). Statistical inference in vector autoregressions with possibly integrated processes. *Journal of Econometrics*, 66(1-2), 225-250. \doi{10.1016/0304-4076(94)01616-8}
#'
#' Diks, C. and Panchenko, V. (2006). A new statistic and practical guidelines for nonparametric Granger causality testing. *Journal of Economic Dynamics and Control*, 30(9-10), 1647-1669. \doi{10.1016/j.jedc.2005.08.008}
#'
#' Kosztyan, Z. T., Matrai, R. and Kiss, D. (2026). Is the world getting smaller? Network-driven heterogeneity in gravity-model elasticities. *The World Economy*, 49(11), 2478-2525. \doi{10.1111/twec.70148}
#' @seealso [causal_methods()], [causal_drivers()], [panel_causality()],
#'   [compare_methods()], [placebo_test()], [robustness()],
#'   [predict.causres()]
#' @export
#' @examples
#' # --- two series, Granger causality --------------------------------------
#' # A VAR where x drives y.  The print shows both directed tests (lag chosen
#' # by BIC, F statistic, p-value, Bayes factor) and the class of the pair.
#' d <- sim_var(400, scenario = "x_to_y", seed = 1)
#' r <- rcausality(d$x, d$y)
#' r
#' r$pairs          # one row per directed test
#' r$classes        # the verdict per pair: x_to_y / y_to_x / bidirectional / ...
#'
#' # --- the three options: method, Bayesian variant, cores ---------------------
#' d3 <- sim_var(400, n_series = 3, seed = 2)
#' rcausality(d3, method = "granger", bayes = TRUE, prior = "bic")$classes
#' \donttest{
#' summary(r)
#' plot(r, "all")                    # network, p-value matrix and lag scan
#' rcausality(d3, method = "bayes_granger", control = bayes_control(n_iter = 300),
#'            seed = 1)$bayes                         # spike-and-slab: BF10, effect interval, R-hat
#'
#' # --- model-free methods on a non-linear coupling --------------------------
#' # (a second or two each; `cores = 2` halves the time, see ?Rcausality-parallel)
#' nl <- sim_nonlinear(300, "x_to_y", seed = 2)       # stochastic, no linear trace of the link
#' rcausality(nl, method = "cmi", max_lag = 2, n_perm = 50, seed = 1)$classes   # kNN conditional MI
#' lg <- sim_logistic(200, "x_to_y", seed = 2)        # deterministic maps
#' rcausality(lg, method = "te", n_perm = 100, seed = 1, cores = 2)$classes$relation
#' rcausality(lg, method = "te", method_args = list(bins = 6), n_perm = 100, seed = 1)$pairs
#' rcausality(lg, method = "ccm", n_perm = 50, seed = 1)$classes$relation  # convergent cross mapping
#' rcausality(sim_logistic_bdc(500, "x_to_y", seed = 1), method = "dc")$classes$relation
#'
#' # --- conditional (multivariate) Granger causality and blocks of series ---
#' rcausality(d3, conditional = TRUE)$classes$relation
#' rcausality(d3[, 1:2], d3[, 3, drop = FALSE], directions = "forward")
#'
#' # --- the Bayesian variant of a surrogate method ---------------------------
#' rcausality(lg, method = "te", bayes = TRUE, n_perm = 100, seed = 1)$bayes
#'
#' # --- generalisations of the Granger test ----------------------------------
#' rcausality(d, method = "toda_yamamoto", method_args = list(d_max = 1))$classes  # integrated series
#' rcausality(lg, method = "diks_panchenko")$classes                             # nonlinear
#' rcausality(d, lag = "first")$pairs[, c("cause", "effect", "lag", "p_value")]  # first order
#' # short series: the hierarchical Minnesota prior (tightness integrated out)
#' rcausality(sim_var(30, seed = 3), method = "bayes_granger", prior = "minnesota",
#'            control = bayes_control(n_iter = 200), seed = 1)$pairs[, c("cause", "effect", "bf10")]
#'
#' # --- one explained variable: which series drive v2? -----------------------
#' rcausality(d3, target = "v2")
#'
#' # --- panel data (long format): the same options ---------------------------
#' pd <- sim_panel(n_units = 20, n_time = 25, seed = 1)
#' rcausality(pd, id = "unit", time = "time", method = "granger", bayes = TRUE, prior = "bic")
#' }
rcausality <- function(x, y = NULL, method = "granger", bayes = FALSE, cores = getOption("Rcausality.cores", 1L),
                       lag = "auto", max_lag = 8L, lag_criterion = c("BIC", "AIC", "HQ"),
                       conditional = FALSE, directions = c("both", "forward"),
                       alpha = 0.05, adjust = "holm", bf_threshold = 3, prob_threshold = 0.5, prior_prob = 0.5,
                       prior = NULL, control = bayes_control(), n_perm = 200L, seed = NULL, impute = TRUE,
                       target = NULL, id = NULL, time = NULL, panel = panel_control(), method_args = list(), ..., m = NULL) {
  cl <- match.call()
  method <- match.arg(method, .method_names)
  lag_criterion <- match.arg(lag_criterion); directions <- match.arg(directions)
  lag_default <- missing(lag) || is.null(lag); if (is.null(lag)) lag <- "auto"
  if (is.numeric(lag) && (length(lag) != 1L || !is.finite(lag) || lag < 1 || lag != round(lag)))
    stop("`lag` must be a positive integer, \"auto\" or \"first\".", call. = FALSE)
  if (is.null(id) && !is.null(time))
    stop("`time` names the time index of a panel and needs `id`; for a time-series analysis remove the time column from the data.", call. = FALSE)
  if (is.null(id) && is.null(target) && is.data.frame(x)) {
    tl <- names(x)[tolower(names(x)) %in% c("time", "year", "date", "period", "month", "quarter")]
    tl <- tl[vapply(x[tl], function(v) is.numeric(v) && !anyNA(v) && length(v) > 2L && length(unique(diff(v))) == 1L, TRUE)]
    if (length(tl))
      message(sprintf("rcausality: the column(s) %s look like a time index and are analysed as variables; remove them from the data if that is not intended.",
                      paste(tl, collapse = ", ")))
  }
  if (identical(method, "bayes_granger") && !missing(bayes) && identical(bayes, FALSE))
    message("rcausality: method = \"bayes_granger\" is the Bayesian regression; bayes = FALSE is ignored.")
  st <- .resolve_settings(method, method_args, c(list(...), if (!is.null(m)) list(m = m)))
  A <- c(mget(setdiff(names(formals(rcausality)), c("x", "y", "method_args", "...", "m"))), st$tune,
         list(method_args = st$own))
  if (!is.null(target)) return(.drivers_run(x, y, A, cl, lag_default))
  if (!is.null(id)) return(.panel_run(x, y, A, cl))
  .rcausality_ts(x, y, A, cl)
}

# One time-series analysis from the complete argument list `A` (the formals of
# rcausality() other than the data, the tuning settings and the method's own
# settings).  Also the work horse of the panel and the target analyses.
.rcausality_ts <- function(x, y, A, cl = NULL) {
  method <- A$method; bayes <- isTRUE(A$bayes); control <- A$control
  if (!inherits(control, "bayes_control")) control <- do.call(bayes_control, as.list(control))
  if (!is.null(A[["prior"]])) control$prior <- match.arg(A[["prior"]], c("bic", "g", "hyper_g", "spike_slab", "normal", "minnesota"))
  prior <- control$prior
  if (method == "bayes_granger") bayes <- TRUE
  sampler <- bayes && method %in% c("granger", "bayes_granger") && prior %in% c("spike_slab", "normal", "minnesota")
  if (bayes && method == "granger" && !sampler) method <- "bayes_granger"
  if (sampler) method <- "granger"                      # the F-test run supplies the lags; the sampler adds the evidence
  old_rng <- .rng_setup(A$seed); on.exit(.rng_restore(old_rng), add = TRUE)   # the seed also covers the Bayesian stage
  run_args <- A[intersect(names(A), names(formals(.rcausality_run)))]
  # The Bayes factor that accompanies a frequentist regression test is the exact one of the unit-information prior
  # (the g prior with g = number of usable observations), which the BIC approximates: the approximation overstates
  # the evidence of regressions with few residual degrees of freedom.  An explicit `prior` is respected.
  run_args$method <- method
  info_default <- !bayes && !prior %in% c("bic", "g", "hyper_g")
  run_args$prior <- if (info_default) "g" else if (prior %in% c("bic", "g", "hyper_g")) prior else "bic"
  run_args$prior_par <- if (info_default) list(g_scale = 1, a = 3) else .prior_par(control)
  # (the call object must not travel through do.call(), which would evaluate it)
  r <- do.call(function(...) .rcausality_run(x = x, y = y, cl = cl, ...), run_args)
  r <- .restrict_to_grid(r)
  if (!bayes) .resolution_check(r)
  .persistence_check(r, bayes)
  if (!bayes) .het_check(r)
  .determined_check(r)
  if (!bayes) .coarse_check(r)
  if (sampler) r <- .bayes_granger_pairs(r, control, A$prior_prob, A$bf_threshold, A$prob_threshold, A$conditional, A$cores,
                                         lag = A$lag)
  else if (bayes && method != "bayes_granger") r <- .bayes_postprocess(r, control, A$prior_prob, A$bf_threshold, A$prob_threshold)
  # a requested order that the data cannot carry is reduced by the regression methods: say so
  if (is.numeric(A$lag) && method %in% c("granger", "bayes_granger", "toda_yamamoto", "rf_granger") &&
      any(r$pairs$lag != A$lag, na.rm = TRUE))
    warning(sprintf("%s: the requested order lag = %d was changed to %s (the series of length %d do not carry the requested order).",
                    method, as.integer(A$lag), paste(sort(unique(stats::na.omit(r$pairs$lag))), collapse = ", "), r$settings$n),
            call. = FALSE)
  r$settings$bayes <- bayes
  if (!is.null(r$pairs$bf10)) r$pairs$bf10[is.nan(r$pairs$bf10)] <- NA_real_     # (tests that involve a constant series)
  # the BIC approximation of a Bayes factor presupposes many more observations than coefficients
  ty_fallback <- bayes && method == "toda_yamamoto" && !prior %in% c("bic", "g", "hyper_g")
  if (ty_fallback)
    message(sprintf(paste0("toda_yamamoto: the lag-augmented regression has closed-form Bayes factors only (prior = \"bic\", \"g\", ",
                           "\"hyper_g\"); the BIC form is used instead of prior = \"%s\"."), prior))
  if (bayes && (prior == "bic" || ty_fallback) && method %in% c("granger", "bayes_granger", "toda_yamamoto") &&
      any(r$pairs$df2 < 10, na.rm = TRUE))
    warning(sprintf(paste0("prior = \"bic\": with %d residual degrees of freedom the BIC approximation overstates the evidence; ",
                           "prior = \"g\" gives the exact Bayes factor of the unit-information prior that it approximates."),
                    as.integer(min(r$pairs$df2, na.rm = TRUE))), call. = FALSE)
  if (bayes) {
    # the lag-augmented regressions have closed-form Bayes factors only: any other prior falls back to the BIC one
    r$settings$prior <- if (ty_fallback) "bic" else prior
    r$settings$bayes_control <- control; r$settings$decision <- control$decision
  }
  # the full argument list, so that the analysis can be repeated on other data
  # (predict(refit = TRUE), placebo_test(), robustness(), benchmark_methods())
  A$method <- if (sampler) "granger" else method; A[["prior"]] <- prior; A$bayes <- bayes; A$control <- control
  r$args <- .flat_args(A)
  r
}

# Notes on the data.  Missing values are interpolated by default; the user is told how many and where, and warned when
# a large share of a series is interpolated (interpolation manufactures serial dependence).
.note_missing <- function(D, impute, who, unit = NULL) {
  bad <- vapply(D, function(v) sum(!is.finite(suppressWarnings(as.numeric(v)))), 0L)
  if (!any(bad > 0L) || !isTRUE(impute)) return(invisible(bad))
  where <- paste(sprintf("%s: %d", names(bad)[bad > 0L], bad[bad > 0L]), collapse = ", ")
  message(sprintf("%s: %d missing or non-finite value(s) were replaced by linear interpolation (%s)%s; impute = FALSE stops instead.",
                  who, sum(bad), where, if (is.null(unit)) "" else paste0(" in ", unit)))
  share <- bad / nrow(D)
  share[share >= 1] <- 0                                # a series without any finite value stops with an error below
  if (any(share > 0.1))
    warning(sprintf("%s: more than 10%% of the values of %s were interpolated; interpolation creates serial dependence and can produce spurious links.",
                    who, paste(names(share)[share > 0.1], collapse = ", ")), call. = FALSE)
  invisible(bad)
}
# constant and exactly collinear series, and series shorter than the advisable minimal length of the method
.note_series <- function(D, method) {
  if (nrow(D) < 8L) stop(sprintf("The series have %d time points; no method of the package can be applied to fewer than 8.", nrow(D)), call. = FALSE)
  s <- vapply(D, function(v) stats::sd(v), 0)
  if (any(!is.finite(s) | s == 0))
    warning(sprintf("Constant series (%s): no test involving them is meaningful.", paste(names(D)[!is.finite(s) | s == 0], collapse = ", ")),
            call. = FALSE)
  ok <- which(is.finite(s) & s > 0)
  if (length(ok) >= 2L && nrow(D) >= 3L) {
    R <- abs(stats::cor(D[, ok, drop = FALSE])); diag(R) <- 0
    if (any(R > 1 - 1e-12)) {
      ij <- which(R > 1 - 1e-12 & upper.tri(R), arr.ind = TRUE)
      stop(sprintf(paste0("Exactly collinear series (%s): one is a linear function of the other, so that neither can add information ",
                          "about the other and the tests between them are not defined; remove one of them."),
                   paste(sprintf("%s ~ %s", names(D)[ok][ij[, 1]], names(D)[ok][ij[, 2]]), collapse = ", ")), call. = FALSE)
    }
  }
  # "dfe": the classifier knows binary Markov chains of short memory.  A continuous series that is cut at its mean
  # is a binary sequence with long memory, which the classifier reads as the trace of a hidden process: 93 of 100
  # pairs of independent first-order autoregressions (coefficient 0.5, 300 time points) were classified as common_cause.
  if (method == "dfe" && any(vapply(D, function(v) length(unique(v[is.finite(v)])) > 2L, TRUE)))
    warning(paste0("dfe: the series are not binary and are cut at their mean. The classifier is trained on binary Markov chains with a ",
                   "memory of one or two steps; a continuous autocorrelated series that is cut at its mean has longer memory, and ",
                   "unrelated series of this kind are usually classified as common_cause."), call. = FALSE)
  mn <- .method_info(method)$min_n
  if (!is.null(mn) && nrow(D) < mn)
    warning(sprintf("%s: the series have %d time points, fewer than the advisable minimum of %d for this method (see causal_methods()); the result is unreliable.",
                    method, nrow(D), as.integer(mn)), call. = FALSE)
  invisible(NULL)
}

# Resolution of the resampling tests.  A Monte Carlo p value cannot be smaller than 1 / (B + 1) with B replicates, and
# the circular-shift surrogates of a series of length n are a finite set of about 0.8 n distinct shifts (they are
# drawn without replacement, see .shift_draws()).  When this bound cannot pass a Holm or Bonferroni correction over the
# M tests of the analysis, no link can be declared whatever the data: the user is told so, with the remedy.
.resolution_check <- function(r) {
  s <- r$settings; adjust <- s$adjust; alpha <- s$alpha
  stepup <- c("BH", "BY", "hochberg", "hommel")
  if (is.null(r$null) || is.null(adjust) || !adjust %in% c("holm", "bonferroni", stepup) ||
      !s$method %in% c("cmi", "cmi_lag", "te", "ccm", "pid", "rf_granger", "markov", "df", "ccdh", "hirata", "diks_panchenko"))
    return(invisible(NULL))
  B <- vapply(r$null, function(z) { v <- if (is.list(z)) z$null else z; if (is.numeric(v)) length(v) else 0L }, 0L); B <- B[B > 0L]
  M <- nrow(r$pairs)
  if (!length(B) || M / (min(B) + 1) <= alpha) return(invisible(NULL))
  if (adjust %in% stepup) {
    # a step-up correction can declare links only when several tests attain the smallest p value together
    need <- as.integer(ceiling(M / ((min(B) + 1) * alpha)))
    warning(sprintf(paste0("%s: with the smallest attainable p-value, 1/%d, the %s correction over %d tests can declare a link at alpha = %g ",
                           "only when at least %d tests attain that value together%s. Use longer series, more replicates where the ",
                           "number of circular shifts allows, or fewer series at once."),
                    s$method, min(B) + 1L, adjust, M, alpha, need, if (need > M) " (which is impossible here)" else ""), call. = FALSE)
    return(invisible(NULL))
  }
  shift <- .uses_shifts(r)
  L <- if (shift && !is.null(r$by_lag) && nrow(r$by_lag)) max(1L, max(table(paste(r$by_lag$cause, r$by_lag$effect)))) else 1L   # lags of a scan
  warning(sprintf(paste0("%s: the smallest attainable p-value, 1/%d, cannot pass the %s correction over %d tests at alpha = %g, ",
                         "so that no link can be declared. %s"),
                  s$method, min(B) + 1L, adjust, M, alpha,
                  if (shift) sprintf(paste0("A series of length %d has only %d admissible circular shifts%s; about %d time points ",
                                            "would be needed. Use longer series, test fewer series at once%s, or use a method ",
                                            "with a permutation or analytic reference distribution%s."),
                                     s$n, min(B), if (L > 1L) sprintf(" for a scan of %d lags", L) else "",
                                     ceiling(1.25 * M * L / alpha), if (L > 1L) ", fix the lag" else "",
                                     if (s$method %in% c("cmi", "cmi_lag"))
                                       " (for weakly autocorrelated series: method_args = list(reference = \"local\"))"
                                     else if (s$method == "diks_panchenko")
                                       " (for series whose pasts are at most weakly dependent: method_args = list(reference = \"normal\"))"
                                     else "")
                  else "Increase n_perm."), call. = FALSE)
  invisible(NULL)
}

# The largest number of lags L of a scan that leaves at least `n_need` circular shifts spaced by L (the numbers of
# admissible shifts are those of .te_scan() and of .cmi_pair(); `history`: the longest conditioning history).
.scan_cap <- function(n, n_need, method, history = 1L) {
  count <- function(L) {
    if (method %in% c("te", "te_cond")) { N <- if (method == "te") n else n - max(L, history); lo <- max(2L * max(L, history) + 1L, N %/% 10L) }
    else { N <- n - if (method == "cmi_lag") L + history else max(L, history); lo <- max(2L, N %/% 10L) }
    if (N - 2L * lo < 0L) 0L else (N - 2L * lo) %/% L + 1L
  }
  L <- 1L
  while (L < n && count(L + 1L) >= n_need) L <- L + 1L
  L
}

# Does the reference distribution of the analysis come from circular shifts?
.uses_shifts <- function(r) {
  s <- r$settings
  s$method %in% c("te", "ccm", "rf_granger", "ccdh", "hirata", "markov") ||
    (s$method %in% c("cmi", "cmi_lag") && !identical(.method_args(s$method, s$method_args)$reference, "local")) ||
    (s$method %in% c("diks_panchenko", "pid") && any(r$extra$reference == "conditional shifts", na.rm = TRUE))
}

# Persistence of the series and circular shifts.  Neighbouring shifts of autocorrelated series give dependent
# replicates: for two series with first-order autocorrelations r1 and r2 a cross-dependence statistic decorrelates over
# about (1 + r1 r2) / (1 - r1 r2) shifts (the variance inflation of a cross-correlation, Bartlett 1946), so that the
# S admissible shifts are worth about S (1 - r1 r2) / (1 + r1 r2) independent ones.  A test at level alpha needs
# 1 / alpha of them at the very least.  In simulations (two independent first-order autoregressions with coefficient
# 0.9) the tests reject in 6 to 10% of the cases at the 5% level with about 8 independent shifts (100 time points), in
# 5 to 8% with about 25 (300 time points) and in 4 to 6% with about 84 (1000 time points): the warning is given below
# 2 / alpha.
.acf1 <- function(D) vapply(D, function(v) {
  a <- tryCatch(stats::acf(v, 1L, plot = FALSE, na.action = stats::na.pass)$acf[2L], error = function(e) NA_real_)
  if (is.finite(a)) min(max(a, 0), 0.999) else 0 }, 0)
.shift_eff <- function(D, cause, effect) {
  rho <- .acf1(D)
  pr <- rho[cause] * rho[effect]
  unname(0.8 * nrow(D) * (1 - pr) / (1 + pr))
}
.persistence_check <- function(r, bayes = FALSE) {
  s <- r$settings; D <- r$data; P <- r$pairs
  if (is.null(D) || is.null(P) || !nrow(P) || !all(c("cause", "effect") %in% names(P))) return(invisible(NULL))
  ok <- P$cause %in% names(D) & P$effect %in% names(D)
  if (!any(ok)) return(invisible(NULL))
  if (s$method %in% c("granger", "bayes_granger") && bayes) return(invisible(NULL))   # (the Bayes-factor rules stay near their rates)
  if (s$method == "granger") {
    # The F test presupposes stationary series.  For nearly integrated series it rejects too often in samples of any
    # practical length (two independent first-order autoregressions with coefficient 0.9: 10% at the 5% level for 50
    # observations, 7% for 200; coefficient 0.98 and 200 observations: 12%); the lag-augmented test does not.
    rho <- .acf1(D); hi <- pmin(rho[P$cause[ok]], rho[P$effect[ok]]) > 0.8
    if (any(hi)) {
      nm <- unique(c(P$cause[ok][hi], P$effect[ok][hi]))
      warning(sprintf(paste0("%s: the series %s are strongly autocorrelated (first-order autocorrelation %s). For nearly integrated ",
                             "or trending series the F test rejects too often; method = \"toda_yamamoto\" or differenced series are advisable."),
                      s$method, paste(nm, collapse = ", "), paste(sprintf("%.2f", rho[nm]), collapse = ", ")), call. = FALSE)
    }
    return(invisible(NULL))
  }
  if (s$method == "diks_panchenko" && all(c("cause", "effect", "reference", "filter") %in% names(r$extra))) {
    k <- match(paste(P$cause, P$effect), paste(r$extra$cause, r$extra$effect))
    P$reference <- r$extra$reference[k]; P$filter <- r$extra$filter[k]
    # The normal reference of the kernel statistic (on request: reference = "normal") estimates its long-run variance
    # with n^(1/4) autocovariances, which is too few for strongly autocorrelated series (two independent first-order
    # autoregressions with coefficient 0.9: 6 to 8% rejections at the 5% level for 100 to 1000 observations), and it
    # is biased towards rejection when the past of the cause depends on the past of a persistent effect (see
    # .dp_stat()).  The test with this reference is meant for filtered series (Diks and Panchenko 2006 apply it to the
    # residuals of a vector autoregression): 1.4 to 3.3% with filter = "var".
    nrm <- ok & P$reference %in% "normal" & P$filter %in% "none"
    if (any(nrm)) {
      rho <- .acf1(D); hi <- pmin(rho[P$cause[nrm]], rho[P$effect[nrm]]) > 0.8
      if (any(hi)) {
        nm <- unique(c(P$cause[nrm][hi], P$effect[nrm][hi]))
        warning(sprintf(paste0("%s: the series %s are strongly autocorrelated (first-order autocorrelation %s). The normal reference ",
                               "of the kernel statistic then rejects too often, in particular in the reverse direction of a link; ",
                               "the default reference (conditional shifts) does not, and method_args = list(filter = \"var\") applies ",
                               "the test to the residuals of a vector autoregression, as is usual for this test (it then tests for ",
                               "nonlinear dependence beyond the linear one, which the filter removes)."),
                        s$method, paste(nm, collapse = ", "), paste(sprintf("%.2f", rho[nm]), collapse = ", ")), call. = FALSE)
      }
    }
  }
  if (!.uses_shifts(r) && s$method != "pid") return(invisible(NULL))
  if (s$method == "diks_panchenko" && !is.null(P$filter)) ok <- ok & P$filter %in% "none"   # (filtered series are not persistent)
  if (!any(ok)) return(invisible(NULL))
  # (short series that are not autocorrelated are the business of the length and the resolution warnings)
  rho <- .acf1(D); pr <- unname(rho[P$cause[ok]] * rho[P$effect[ok]])
  se <- .shift_eff(D, P$cause[ok], P$effect[ok]); low <- se < 2 / s$alpha & pr > 0.2
  if (!any(low)) return(invisible(NULL))
  lab <- sprintf("%s -> %s", P$cause[ok][low], P$effect[ok][low])
  if (s$method == "pid" && !identical(s$method_args$reference, "conditional")) {
    # the permutation within the levels of the coded own past keeps the conditioning but not the serial dependence of
    # the cause (two independent first-order autoregressions with coefficient 0.9: 10% and 7% rejections at the 5%
    # level for 100 and 300 time points)
    warning(sprintf(paste0("pid: the series of %s are strongly autocorrelated for their length (the diagnostic of the shift-based tests ",
                           "gives about %d independent replicates, fewer than 2/alpha = %d). The permutation of this test does not keep ",
                           "the serial dependence of the cause, and its p-values can then be too small. Longer series, differenced or ",
                           "filtered series, or a regression-based method are advisable."),
                    if (length(lab) > 3L) sprintf("%s and %d further pairs", paste(lab[1:3], collapse = ", "), length(lab) - 3L)
                    else paste(lab, collapse = ", "), as.integer(floor(min(se[low]))), as.integer(ceiling(2 / s$alpha))), call. = FALSE)
    return(invisible(NULL))
  }
  warning(sprintf(paste0("%s: the series of %s are strongly autocorrelated for their length: their circular shifts are worth about %d ",
                         "independent replicates, fewer than the 2/alpha = %d from which the tests were close to their level in ",
                         "simulations, and the p-values of %s can be too small. Longer series, differenced or filtered series%s, or a regression-based method ",
                         "(\"granger\", \"toda_yamamoto\") are advisable."),
                  s$method, if (length(lab) > 3L) sprintf("%s and %d further pairs", paste(lab[1:3], collapse = ", "), length(lab) - 3L)
                            else paste(lab, collapse = ", "),
                  as.integer(floor(min(se[low]))), as.integer(ceiling(2 / s$alpha)),
                  if (length(lab) > 1L) "these pairs" else "this pair",
                  if (s$method == "diks_panchenko") " (method_args = list(filter = \"var\"))" else ""), call. = FALSE)
  invisible(NULL)
}

# the stored argument list in the form accepted by rcausality(): the method's
# own settings and the tuning settings that differ from their defaults are
# folded into `method_args`
.flat_args <- function(A) {
  tn <- names(.tuning_defaults)
  ch <- tn[!vapply(tn, function(nm) identical(A[[nm]], .tuning_defaults[[nm]]), TRUE)]
  rel <- intersect(ch, .method_info(A$method)$tuning)
  out <- A[setdiff(names(A), c(tn, "method_args", "m_family"))]
  out$method_args <- c(A$method_args, A[rel])
  c(out, A[setdiff(ch, rel)])                           # settings foreign to the method travel through `...`
}

.rcausality_run <- function(x, y = NULL, method, lag, max_lag, lag_criterion, alpha, bf_threshold, adjust,
                            conditional, prior_prob, impute, directions, k, k_perm, n_perm, history, seed, memory,
                            iid_injection, k_nn, n_bins, n_markov, cores, max_constr, m, tau, k_range, prep,
                            prob_threshold, K, n_surr, theta, model, bins, k_target, k_source, rate, prior,
                            method_args = list(), cl, prior_par = NULL, m_family = NULL) {
  lag_criterion <- match.arg(lag_criterion, c("BIC", "AIC", "HQ"))
  directions <- match.arg(directions, c("both", "forward"))
  .check_number(alpha, "alpha", lower = 0, upper = 1, strict = TRUE)
  .check_number(max_lag, "max_lag", lower = 1, integer = TRUE)
  .check_number(n_perm, "n_perm", lower = 1, integer = TRUE)
  .check_number(prior_prob, "prior_prob", lower = 0, upper = 1, strict = TRUE)
  .check_number(bf_threshold, "bf_threshold", lower = 0, strict = TRUE)
  first_lag <- identical(lag, "first")
  if (!identical(lag, "auto") && !first_lag) {
    if (is.character(lag)) stop("`lag` must be an integer, \"auto\" or \"first\".", call. = FALSE)
    .check_number(lag, "lag", lower = 0, integer = TRUE)
  }
  if (!identical(history, "auto")) .check_number(history, "history", lower = 1, integer = TRUE)
  if (!adjust %in% c(stats::p.adjust.methods)) stop("`adjust` must be one of ", paste(stats::p.adjust.methods, collapse = ", "), ".", call. = FALSE)
  iid_injection <- match.arg(iid_injection, c("allowed", "forbidden"))
  fr <- .as_series_frame(x, y, xname = if (is.null(y)) "x" else "x", yname = "y")
  if (is.null(fr$Y) && ncol(fr$X) < 2L) {
    stop("Supply two series (x, y) or a data frame with at least two columns.",
         call. = FALSE)
  }
  X <- fr$X; Y <- fr$Y
  if (!is.null(Y) && nrow(X) != nrow(Y)) stop("x and y must have the same length.", call. = FALSE)
  # make unique names across blocks
  all <- if (is.null(Y)) X else cbind(X, Y)
  if (anyDuplicated(names(all))) {
    if (!is.null(Y)) {
      names(X) <- paste0("x.", names(X)); names(Y) <- paste0("y.", names(Y))
      all <- cbind(X, Y)
    } else stop("Column names must be unique.", call. = FALSE)
  }
  .note_missing(all, impute, "rcausality")
  all[] <- lapply(all, function(v) { v <- as.numeric(v); if (impute) .prep_na(v) else v })
  if (any(!is.finite(as.matrix(all)))) stop("Missing/non-finite values present (use impute = TRUE).", call. = FALSE)
  .note_series(all, method)
  causes <- if (is.null(Y)) names(all) else names(X)
  effects <- if (is.null(Y)) names(all) else names(Y)
  grid <- expand.grid(cause = causes, effect = effects, stringsAsFactors = FALSE)
  if (!is.null(Y) && directions == "both") {
    grid <- rbind(grid, setNames(grid[, 2:1], c("cause", "effect")))
    causes <- effects <- names(all)
  }
  grid <- grid[grid$cause != grid$effect, , drop = FALSE]
  if (nrow(grid) == 0L) stop("No pairs to test.", call. = FALSE)
  if (isTRUE(conditional) && ncol(all) <= 2L) {          # nothing to condition on
    message("rcausality: conditional = TRUE has no effect with two series (there is no other series to condition on).")
    conditional <- FALSE
  }

  .n_cores(cores)                                        # validates the argument
  old_rng <- .rng_setup(seed); on.exit(.rng_restore(old_rng), add = TRUE)
  sc <- .split_cores(nrow(grid), cores); pair_cores <- sc$outer; inner_cores <- sc$inner
  if (method %in% c("cmi", "cmi_lag", "markov", "df", "ccdh", "te", "ccm", "pid", "rf_granger") ||
      (method == "diks_panchenko" && .dp_ref_arg(.method_args("diks_panchenko", method_args)$reference) != "normal")) {
    # enough replicates for the smallest p value, 1 / (B + 1), to pass the correction over the tests of the analysis
    # (circular-shift surrogates are limited by the number of distinct shifts: see .shift_draws(), .resolution_check())
    need <- as.integer(ceiling(nrow(grid) / alpha))        # (also with adjust = "none": a target analysis adjusts afterwards)
    if (method == "ccdh") {
      if (n_surr < need) {
        message("ccdh: n_surr raised from ", n_surr, " to ", need,
                " so that the surrogate p-values can reach significance after multiple-testing correction.")
        n_surr <- need
      }
    } else if (n_perm < need) {
      # (when the number of admissible circular shifts is the limit, more replicates do not help: .resolution_check())
      shifts_bind <- method %in% c("te", "ccm", "rf_granger", "cmi", "cmi_lag", "diks_panchenko") && 0.8 * nrow(all) < need
      if (!shifts_bind)
        message(method, ": n_perm raised from ", n_perm, " to ", need,
                " so that the smallest attainable p-value can pass the multiple-testing correction.")
      n_perm <- need
    }
  }
  fixed_lag <- if (identical(lag, "auto") || first_lag) NULL else as.integer(lag)
  cand_lags <- if (is.null(fixed_lag)) seq_len(max(1L, as.integer(max_lag))) else fixed_lag
  # conditioning history of the effect: the VAR order of the pair when "auto"
  hist_vec <- if (identical(history, "auto")) {
    vapply(seq_len(nrow(grid)), function(i)
      max(1L, min(5L, as.integer(max_lag), select_lag(all[, c(grid$effect[i], grid$cause[i])], max_lag, lag_criterion))), 0L)
  } else rep(max(1L, as.integer(history)), nrow(grid))
  if (is.null(k)) k <- if (method == "df") 10L else 5L
  if (is.null(theta)) theta <- if (method == "hirata") 0.4 else 0.45
  if (length(method_args) && !method %in% names(.method_arg_defaults))
    warning(sprintf("`method_args` is not used by method \"%s\".", method), call. = FALSE)
  if (method %in% c("toda_yamamoto", "diks_panchenko", "rf_granger")) {
    ma <- .method_args(method, method_args)
    if (!.method_info(method)$conditional && conditional)
      stop(sprintf("method = \"%s\" supports pairwise analysis only.", method), call. = FALSE)
    pair_fun <- switch(method,
      toda_yamamoto = .ty_pair_fun(fixed_lag, max_lag, lag_criterion, d_max = as.integer(ma$d_max), prior = prior,
                                   first = first_lag, alpha = alpha, prior_par = prior_par, vcov = .vcov_arg(ma$vcov)),
      diks_panchenko = .dp_pair_fun(fixed_lag, max_lag, lag_criterion, eps = ma$eps, C = ma$C,
                                    filter = match.arg(ma$filter, c("none", "var")), n_perm = n_perm,
                                    reference = .dp_ref_arg(ma$reference)),
      rf_granger = .rf_pair_fun(fixed_lag, max_lag, lag_criterion, n_perm = n_perm, settings = ma))
    out <- .rcausality_directed(all, grid, causes, effects, method, pair_fun, alpha, adjust, prior_prob, bf_threshold,
                                conditional, cl, cores, rule = "p")
    out$settings$method_args <- ma
    return(out)
  }
  if (method == "pcmci") {
    r <- pcmci(all, tau_max = if (is.null(fixed_lag)) as.integer(max_lag) else fixed_lag,
               tau_min = if (is.null(fixed_lag)) 1L else fixed_lag, alpha = alpha,
               adjust = if (adjust == "none") "none" else "fdr_bh")
    keep <- paste(r$pairs$cause, r$pairs$effect) %in% paste(grid$cause, grid$effect)
    if (!all(keep)) {
      bk <- paste(r$by_lag$cause, r$by_lag$effect) %in% paste(grid$cause, grid$effect)
      pm <- r$pcmci
      r <- .new_causres(r$pairs[keep, , drop = FALSE], series = names(all), causes = causes, effects = effects,
                        method = "pcmci", conditional = TRUE, alpha = alpha, bf_threshold = NA_real_,
                        adjust = r$settings$adjust, n = nrow(all), call = cl, data = all, instant = r$instant,
                        by_lag = r$by_lag[bk, , drop = FALSE])
      r$pcmci <- pm; r$settings$test <- pm$test; r$settings$plus <- pm$plus
    }
    r$call <- cl
    return(r)
  }
  if (method == "hirata") {
    if (conditional) stop("method = \"hirata\" supports pairwise analysis only.", call. = FALSE)
    rh <- .rcausality_verdict_pairs(all, grid, causes, effects, alpha, "hirata",
      fun = function(a, b, cores) hirata_rp(a, b, E = as.integer(m), tau = as.integer(tau), rate = rate, n_surr = n_perm,
                                            theta = theta, alpha = alpha, cores = cores),
      cl = cl, cores = cores,
      stat_fun = function(r) list(statistic = r$index, p = r$p_values[c("x_to_y", "y_to_x")], n = r$n,
                                  extra = data.frame(cov_x_to_y = r$coverage[["x_to_y"]], cov_y_to_x = r$coverage[["y_to_x"]],
                                                     excess = r$excess, p_excess = r$p_values[["excess"]])),
      null_fun = function(r, k) list(statistic = r$coverage[[k]], null = r$surr[, if (k == "x_to_y") "cov_xy" else "cov_yx"],
                                     alternative = "greater"))
    rh$settings$theta <- theta
    return(rh)
  }
  if (method == "cd") {
    if (conditional) stop("method = \"cd\" supports pairwise analysis only.", call. = FALSE)
    return(.rcausality_verdict_pairs(all, grid, causes, effects, alpha, "cd",
      fun = function(a, b, cores) krakovska_cd(a, b, E = as.integer(m), tau = as.integer(tau), alpha = alpha, cores = cores),
      cl = cl, cores = cores,
      stat_fun = function(r) list(statistic = c(x_to_y = r$d2[["y"]] - r$d2[["x"]], y_to_x = r$d2[["x"]] - r$d2[["y"]]),
                                  p = c(x_to_y = NA_real_, y_to_x = NA_real_), n = NA_integer_,
                                  extra = data.frame(D2_x = r$d2[["x"]], D2_y = r$d2[["y"]], D2_joint = r$d2[["joint"]],
                                                     D2_indep = r$d2[["indep"]], se_joint = r$se[["joint"]])),
      null_fun = function(r, k) list(statistic = if (k == "x_to_y") r$d2[["y"]] - r$d2[["x"]] else r$d2[["x"]] - r$d2[["y"]],
                                     null = list(mean = 0, sd = sqrt(r$se[["x"]]^2 + r$se[["y"]]^2)), alternative = "greater")))
  }
  pid_ref <- if (method == "pid") .pid_ref_arg(.method_args("pid", method_args)$reference) else "permutation"
  if (method == "pid" && conditional) stop("method = \"pid\" supports pairwise analysis only.", call. = FALSE)
  if (method == "pid" && pid_ref == "permutation")
    return(.rcausality_pid(all, grid, causes, effects, alpha, adjust, lags = cand_lags, history = hist_vec,
                           k_target = as.integer(k_target), k_source = as.integer(k_source), n_perm = n_perm,
                           cl = cl, cores = cores))
  cmi_ref <- if (method %in% c("cmi", "cmi_lag")) match.arg(.method_args(method, method_args)$reference, c("shift", "local")) else "local"
  te_ref <- if (method == "te") .te_ref_arg(.method_args("te", method_args)$reference) else "circular"
  pair_lags <- rep(list(cand_lags), nrow(grid))                    # the lags scanned for every pair
  if (method %in% c("te", "pid") || (method %in% c("cmi", "cmi_lag") && cmi_ref == "shift")) {
    # The shifts of a scan over L lags are spaced by L (see .te_scan()), so that a scan leaves S / L replicates of the S
    # admissible shifts.  The scan is limited so that enough of them remain: at least 40, and as many as the smallest
    # p-value needs to pass a Holm or Bonferroni correction over the tests of the analysis.
    # (m_family: a target analysis corrects within its families of tests -- the links to the target, from the target,
    # among the drivers -- and passes a function that gives the size of the family of every pair, so that the scan of
    # a small family is not limited by a large one; a single number applies to all pairs)
    m_pair <- if (is.function(m_family)) as.integer(m_family(grid$cause, grid$effect))
              else rep(if (!is.null(m_family)) as.integer(m_family) else if (adjust %in% c("holm", "bonferroni")) nrow(grid) else 0L,
                       nrow(grid))
    need_pair <- pmax(40L, as.integer(ceiling(m_pair / alpha)))
    if (is.null(fixed_lag)) {
      full <- max(cand_lags)
      for (nn in sort(unique(need_pair))) {
        te_cap <- .scan_cap(nrow(all), nn, if (method == "te" && te_ref == "conditional") "te_cond" else method, max(hist_vec))
        sel <- need_pair == nn
        if (full > te_cap)
          message(sprintf(paste0("%s: %s%s: a scan up to max_lag = %d would leave fewer than the %d circular ",
                                 "shifts with disjoint windows of lags that %s for series of length %d."),
                          method, if (te_cap > 1L) sprintf("the lags 1 to %d are scanned", te_cap) else "only lag 1 is tested",
                          if (length(unique(need_pair)) > 1L) sprintf(" in the families of %d tests", m_pair[sel][1]) else "",
                          full, nn, if (nn > 40L) sprintf("%d tests at alpha = %g need", m_pair[sel][1], alpha) else "a test needs",
                          nrow(all)))
        pair_lags[sel] <- list(cand_lags[cand_lags <= te_cap])
      }
      if (length(unique(need_pair)) == 1L) cand_lags <- pair_lags[[1L]]
    }
  }
  if (method == "pid")                                             # (with conditional shifts: the scan is limited as above)
    return(.rcausality_pid(all, grid, causes, effects, alpha, adjust, lags = pair_lags, history = hist_vec,
                           k_target = as.integer(k_target), k_source = as.integer(k_source), n_perm = n_perm,
                           cl = cl, cores = cores, reference = pid_ref))
  if (method == "te") {
    if (conditional) stop("method = \"te\" supports pairwise analysis only.", call. = FALSE)
    return(.rcausality_te(all, grid, causes, effects, alpha, adjust, conditional, lags = pair_lags,
                          history = hist_vec, bins = as.integer(bins), n_perm = n_perm, cl = cl, cores = cores, reference = te_ref))
  }
  if (method == "ccm") {
    if (conditional) stop("method = \"ccm\" supports pairwise analysis only.", call. = FALSE)
    return(.rcausality_ccm(all, grid, causes, effects, alpha, adjust, m = as.integer(m), tau = as.integer(tau),
                           n_perm = n_perm, cl = cl, cores = cores))
  }
  if (method == "dfe") {
    if (conditional) stop("method = \"dfe\" supports pairwise analysis only.", call. = FALSE)
    return(.rcausality_dfe(all, grid, causes, effects, alpha, model = model,
                           lag = if (is.null(fixed_lag)) 0L else fixed_lag,
                           prob_threshold = prob_threshold, cl = cl, cores = cores))
  }
  if (method == "ccdh") {
    if (conditional) stop("method = \"ccdh\" supports pairwise analysis only.", call. = FALSE)
    return(.rcausality_ccdh(all, grid, causes, effects, alpha, adjust, m = as.integer(m), tau = as.integer(tau),
                            K = as.integer(K), lag = if (is.null(fixed_lag)) 0L else fixed_lag, n_surr = n_surr,
                            theta = theta, cl = cl, cores = cores))
  }
  if (method == "dc") {
    if (conditional) stop("method = \"dc\" supports pairwise analysis only.", call. = FALSE)
    return(.rcausality_dc(all, grid, causes, effects, alpha, adjust, m = as.integer(m), tau = as.integer(tau),
                          k_range = k_range, prep = match.arg(prep, c("rank", "scale", "none")), eps = 0.05, downsample = 1L,
                          prob_threshold = prob_threshold, cl = cl, cores = cores))
  }
  if (method == "df") {
    return(.rcausality_df(all, grid, causes, effects, alpha, adjust, conditional, lags = cand_lags,
                          history = hist_vec, k = k, k_perm = k_perm, n_perm = n_perm,
                          max_constr = as.integer(max_constr), cl = cl, cores = cores))
  }
  if (method == "markov") {
    if (conditional) stop("method = \"markov\" supports pairwise analysis only.", call. = FALSE)
    return(.rcausality_markov(all, grid, causes, effects, alpha, adjust, memory, iid_injection,
                              k_nn, n_bins, n_markov, k, k_perm, n_perm, cl, directions,
                              cores = cores))
  }
  by_lag <- list()
  one_pair <- function(i) {
    cs <- grid$cause[i]; ef <- grid$effect[i]
    ctrl_names <- if (conditional) setdiff(names(all), c(cs, ef)) else character(0)
    if (method %in% c("cmi", "cmi_lag")) {
      g <- .cmi_pair(all[[ef]], all[[cs]],
                     if (length(ctrl_names)) as.matrix(all[, ctrl_names, drop = FALSE]) else NULL,
                     lags = pair_lags[[i]], history = hist_vec[i], k = k, k_perm = k_perm, n_perm = n_perm,
                     cause_past = (method == "cmi_lag"), cores = inner_cores, reference = cmi_ref)
      return(list(row = data.frame(cause = cs, effect = ef, lag = g$lag,
                                   statistic = g$statistic, df1 = NA_real_, df2 = NA_real_,
                                   p_value = g$p_value, bf10 = NA_real_,
                                   effect_size = g$statistic, r2_gain = NA_real_,
                                   sign = NA_real_, n_eff = g$n_eff, stringsAsFactors = FALSE),
                  by_lag = cbind(cause = cs, effect = ef, g$by_lag),
                  null = list(statistic = g$statistic, null = g$null, alternative = "greater", cond_r2 = g$cond_r2)))
    }
    inv <- all[, c(ef, cs, ctrl_names), drop = FALSE]
    cap <- .cap_lag(nrow(all), ncol(inv))
    ctrl <- if (length(ctrl_names)) as.matrix(all[, ctrl_names, drop = FALSE]) else NULL
    # the test at every candidate order (cheap): the lag plot and the "first" rule
    lags <- seq_len(max(1L, min(as.integer(max_lag), cap)))
    scan <- lapply(lags, function(l) .granger_pair(all[[ef]], all[[cs]], ctrl, l, prior = prior, prior_par = prior_par, vcov = vc))
    bl <- data.frame(cause = cs, effect = ef, lag = lags, statistic = vapply(scan, `[[`, 0, "statistic"),
                     p = vapply(scan, `[[`, 0, "p_value"), stringsAsFactors = FALSE)
    if (first_lag) {
      fs <- .first_lag(lags, bl$p, vapply(scan, `[[`, 0, "bf10"), alpha, bf_threshold, by_bf = method == "bayes_granger")
      p <- fs$lag; g <- scan[[match(p, lags)]]; g$p_value <- fs$p_value; g$bf10 <- fs$bf10
    } else {
      p <- if (is.null(fixed_lag)) select_lag(inv, max_lag, lag_criterion) else fixed_lag
      p <- max(1L, min(p, cap))
      g <- if (p %in% lags) scan[[match(p, lags)]] else .granger_pair(all[[ef]], all[[cs]], ctrl, p, prior = prior, prior_par = prior_par, vcov = vc)
    }
    list(row = data.frame(cause = cs, effect = ef, lag = p,
                          statistic = g$statistic, df1 = g$df1, df2 = g$df2,
                          p_value = g$p_value, bf10 = g$bf10,
                          effect_size = g$effect_size, r2_gain = g$r2_gain,
                          sign = g$sign, n_eff = g$n_eff,
                          stringsAsFactors = FALSE),
         by_lag = bl, het = g$het_p)
  }
  vc <- if (method == "granger") .vcov_arg(.method_args("granger", method_args)$vcov) else "classical"
  res <- .plapply(seq_len(nrow(grid)), one_pair, cores = pair_cores)
  rows <- lapply(res, `[[`, "row")
  by_lag <- do.call(rbind, lapply(res, `[[`, "by_lag"))
  names(by_lag)[names(by_lag) == "cmi"] <- "statistic"
  rownames(by_lag) <- NULL
  pairs <- do.call(rbind, rows)
  instant <- NULL; null_instant <- NULL
  if (method == "cmi_lag") {
    up <- .unordered_pairs(grid, all)
    ins <- .plapply(seq_len(nrow(up)), function(i) {
      ctrl_names <- if (conditional) setdiff(names(all), c(up$v1[i], up$v2[i])) else character(0)
      hp <- max(hist_vec[(grid$cause == up$v1[i] & grid$effect == up$v2[i]) | (grid$cause == up$v2[i] & grid$effect == up$v1[i])])
      g <- .cmi_instant(all[[up$v1[i]]], all[[up$v2[i]]],
                        if (length(ctrl_names)) as.matrix(all[, ctrl_names, drop = FALSE]) else NULL,
                        hp, k, k_perm, n_perm, cores = inner_cores, reference = cmi_ref)
      list(row = data.frame(var1 = up$v1[i], var2 = up$v2[i], cmi = g$cmi, p_value = g$p_value, stringsAsFactors = FALSE),
           null = list(statistic = g$cmi, null = g$null, alternative = "greater"))
    }, cores = pair_cores)
    null_instant <- stats::setNames(lapply(ins, `[[`, "null"), paste(up$v1, "~", up$v2))
    instant <- do.call(rbind, lapply(ins, `[[`, "row"))
    instant$p_adj <- if (adjust == "none") instant$p_value else
      stats::p.adjust(instant$p_value, method = adjust)
    instant$significant <- instant$p_adj < alpha
  }
  pairs <- .finish_pairs(pairs, adjust, alpha, prior_prob, bf_threshold,
                         rule = if (method %in% c("granger", "cmi", "cmi_lag")) "p" else "bf")
  out <- .new_causres(pairs, series = names(all), causes = causes, effects = effects,
                      method = method, conditional = conditional, alpha = alpha,
                      bf_threshold = bf_threshold, adjust = adjust, n = nrow(all),
                      call = cl, data = all, instant = instant, by_lag = by_lag)
  if (method %in% c("cmi", "cmi_lag")) {
    out$null <- stats::setNames(lapply(res, `[[`, "null"), paste(pairs$cause, "->", pairs$effect))
    out$null_instant <- null_instant
    out$settings$method_args <- list(reference = cmi_ref)
  }
  if (method == "granger") {
    out$extra <- data.frame(cause = pairs$cause, effect = pairs$effect, vcov = vc, stringsAsFactors = FALSE,
                            het_p = vapply(res, function(z) if (is.null(z$het)) NA_real_ else z$het, 0))
    out$settings$method_args <- list(vcov = vc)
  }
  out
}
.vcov_arg <- function(v) {
  v <- tolower(if (is.null(v)) "classical" else as.character(v)[1L])
  if (!v %in% c("classical", "hc3", "hac")) stop("`vcov` must be \"classical\", \"HC3\" or \"HAC\".", call. = FALSE)
  v
}
# The past of the cause is a function of the conditioning past (see .cond_r2()): the conditional tests with
# conditional shifts cannot be informative, and the user is told why and what to use instead.
.determined_check <- function(r) {
  s <- r$settings
  if (!s$method %in% c("cmi", "cmi_lag", "rf_granger", "diks_panchenko", "te", "pid")) return(invisible(NULL))
  v <- if (s$method == "diks_panchenko") {
    E <- r$extra; if (is.null(E) || is.null(E$cond_r2)) return(invisible(NULL))
    stats::setNames(E$cond_r2, paste(E$cause, "->", E$effect))
  } else {
    if (is.null(r$null)) return(invisible(NULL))
    vapply(r$null, function(z) if (is.list(z) && !is.null(z$cond_r2) && length(z$cond_r2) == 1L) as.numeric(z$cond_r2) else NA_real_, 0)
  }
  bad <- which(v > 0.99)
  if (!length(bad)) return(invisible(NULL))
  lab <- names(v)[bad]
  warning(sprintf(paste0("%s: for %s the past of the cause is a function of the conditioning past (the past of the effect%s explains ",
                         "%s of its variance), so that no surrogate of the cause differs from it. No test is carried out ",
                         "(p-value 1; Bayes factor 1): a conditional test of predictability is not informative here. Either the cause ",
                         "is redundant given the conditioning series, or the series share a trend (difference them), or they behave as ",
                         "a deterministic system, in which the driver is reconstructed from the driven series. For deterministic ",
                         "systems the state-space methods (\"ccm\", \"dc\", \"cd\", \"ccdh\", \"hirata\") or \"te\" with circular ",
                         "shifts are the tools; a shorter history%s leaves the reconstruction incomplete."),
                  s$method, if (length(lab) > 3L) sprintf("%s and %d further links", paste(lab[1:3], collapse = ", "), length(lab) - 3L)
                            else paste(lab, collapse = ", "),
                  if (isTRUE(s$conditional)) " and of the other series" else "",
                  if (min(v[bad]) >= 0.9995) "all" else sprintf("at least %.1f%%", floor(1000 * min(v[bad])) / 10),
                  if (s$method %in% c("cmi", "cmi_lag", "te", "pid")) " (history = 1)" else " (lag = 1)"), call. = FALSE)
  invisible(NULL)
}
# Coarse conditioning.  "te" conditions on a few symbols of the past of the effect, "pid" on its binary code, and the
# kernel statistic on a window of up to 1.5 standard deviations.  When the past of the cause is itself a function of
# the past of the effect, as in the reverse direction of a link, and the effect is persistent, such a code does not
# resolve that past: the cause keeps information about the present of the effect, and a reference distribution that
# does not keep the dependence of the cause on the past of the effect (circular shifts, the permutation of "pid", the
# normal reference of the kernel test) rejects a true null hypothesis, more often the longer the series.  (Reverse
# direction of y[t] = 0.5 y[t-1] + 1.5 (x[t-1]^2 - 1) + e[t] with a first-order autoregression x of coefficient 0.9:
# "te" with circular shifts rejects in 6% of the series of 300 and in 16 to 19% of the series of 1000 time points, with
# conditional shifts in 4% and 6%; the kernel test with its normal reference in 51% and 77%, with conditional shifts
# in 4%.)  A significant test of such a pair gets a warning: the lag-one autocorrelation of the effect is above
# .coarse_acf and the past of the effect explains more than .coarse_r2 of the variance of the past of the cause.
.coarse_acf <- 0.6; .coarse_r2 <- 0.2
.coarse_check <- function(r) {
  s <- r$settings; D <- r$data; P <- r$pairs
  if (!s$method %in% c("te", "pid", "diks_panchenko") || is.null(D) || is.null(P) || !nrow(P) ||
      !all(c("cause", "effect", "p_value") %in% names(P))) return(invisible(NULL))
  ok <- P$cause %in% names(D) & P$effect %in% names(D) & is.finite(P$p_value) & P$p_value < s$alpha
  {
    E <- r$extra; if (is.null(E) || !all(c("cause", "effect", "reference") %in% names(E))) return(invisible(NULL))
    k <- match(paste(P$cause, P$effect), paste(E$cause, E$effect))
    ok <- ok & E$reference[k] %in% c("circular shifts", "normal", "permutation")
    if (s$method == "diks_panchenko" && !is.null(E$filter)) ok <- ok & E$filter[k] %in% "none"
  }
  if (!any(ok)) return(invisible(NULL))
  rho <- .acf1(D); ok <- ok & unname(rho[P$effect]) > .coarse_acf
  bins <- if (s$method == "pid") 2L else 4L
  r2 <- vapply(which(ok), function(i) {
    e <- D[[P$effect[i]]]; x <- D[[P$cause[i]]]; n <- length(e)
    if (min(length(unique(e)), length(unique(x))) <= bins) return(0)              # (symbolic series are not coded)
    tryCatch(.cond_split(cbind(.std(e)[-n]), cbind(.std(x)[-n]))$r2, error = function(z) 0)
  }, 0)
  bad <- which(ok)[is.finite(r2) & r2 > .coarse_r2]
  if (!length(bad)) return(invisible(NULL))
  lab <- sprintf("%s -> %s", P$cause[bad], P$effect[bad])
  how <- switch(s$method,
    te = paste0("Four symbols do not resolve such a past, and the test with circular shifts can then reject without a link, more often ",
                "the longer the series. method_args = list(reference = \"conditional\") refers the statistic to surrogates that keep ",
                "the dependence of the cause on the past of the effect; method = \"cmi\" conditions on the values themselves."),
    pid = paste0("The binary code does not resolve such a past, and the permutation of this test can then reject without a link, more ",
                 "often the longer the series. method_args = list(reference = \"conditional\") refers the statistic to surrogates ",
                 "that keep the dependence of the cause on the past of the effect; method = \"cmi\" conditions on the values themselves."),
    paste0("The kernel window does not resolve such a past, and the normal reference then rejects without a link. The default ",
           "reference (conditional shifts) keeps the dependence of the cause on the past of the effect."))
  warning(sprintf(paste0("%s: for %s the past of the cause is largely a function of the past of the effect (which explains %s of its ",
                         "variance), and the effect is strongly autocorrelated, as in the reverse direction of a link whose driver is ",
                         "persistent. %s"),
                  s$method, if (length(lab) > 3L) sprintf("%s and %d further links", paste(lab[1:3], collapse = ", "), length(lab) - 3L)
                            else paste(lab, collapse = ", "),
                  sprintf("%.0f%%", 100 * min(r2[is.finite(r2) & r2 > .coarse_r2])), how), call. = FALSE)
  invisible(NULL)
}
# Heteroskedastic errors of the regression tests (see .het_test()): the user is told, with the remedies.
.het_check <- function(r) {
  s <- r$settings; E <- r$extra
  if (!s$method %in% c("granger", "toda_yamamoto") || is.null(E) || !all(c("het_p", "vcov") %in% names(E))) return(invisible(NULL))
  bad <- which(E$het_p < 0.001 & E$vcov == "classical")
  if (!length(bad)) return(invisible(NULL))
  lab <- if (all(c("cause", "effect") %in% names(E))) sprintf("%s -> %s", E$cause[bad], E$effect[bad]) else sprintf("pair %d", bad)
  warning(sprintf(paste0("%s: the errors of the regression%s for %s are heteroskedastic (White test, smallest p = %.1g). The F test ",
                         "presupposes errors whose variance does not depend on the regressors and can otherwise reject without a linear ",
                         "link, for example when the dependence is nonlinear. method_args = list(vcov = \"HC3\") uses a ",
                         "heteroskedasticity-consistent covariance, and vcov = \"HAC\" one that is also consistent when the ",
                         "regression is misspecified and its errors are serially dependent; the nonlinear methods ",
                         "(\"rf_granger\", \"cmi\", \"te\", \"diks_panchenko\") test for dependence of any form."),
                  s$method, if (length(bad) > 1L) "s" else "",
                  if (length(lab) > 3L) sprintf("%s and %d further pairs", paste(lab[1:3], collapse = ", "), length(lab) - 3L) else paste(lab, collapse = ", "),
                  min(E$het_p[bad])), call. = FALSE)
  invisible(NULL)
}

# Markov-chain method driver ------------------------------------------------
.rcausality_markov <- function(all, grid, causes, effects, alpha, adjust, memory, iid,
                               k_nn, n_bins, n_markov, k, k_perm, n_perm, cl, directions,
                               cores = 1L) {
  up <- .unordered_pairs(grid, all)
  # keep the (var1, var2) = (earlier name, later name) orientation used by causres
  sc <- .split_cores(nrow(up), cores); pair_cores <- sc$outer; inner_cores <- sc$inner
  res <- .plapply(seq_len(nrow(up)), function(i)
    .markov_pair(all[[up$v1[i]]], all[[up$v2[i]]], memory, k_nn, n_bins, n_markov, k, k_perm, n_perm,
                 cores = inner_cores), cores = pair_cores)
  rows <- list(); ins <- list(); mk <- list()
  for (i in seq_len(nrow(up))) {
    r <- res[[i]]; a <- up$v1[i]; b <- up$v2[i]; n_eff <- r$xy$n_eff
    for (dir in c("ab", "ba")) {
      g <- if (dir == "ab") r$xy else r$yx
      rows[[length(rows) + 1L]] <- data.frame(
        cause = if (dir == "ab") a else b, effect = if (dir == "ab") b else a,
        lag = 1L, statistic = g$statistic, df1 = NA_real_, df2 = NA_real_,
        p_value = g$p_value, bf10 = NA_real_, effect_size = g$statistic, r2_gain = NA_real_,
        sign = NA_real_, n_eff = n_eff, stringsAsFactors = FALSE)
    }
    ins[[i]] <- data.frame(var1 = a, var2 = b, cmi = r$i0$cmi, p_value = r$i0$p_value,
                           stringsAsFactors = FALSE)
  }
  pairs <- do.call(rbind, rows)
  pairs$p_adj <- if (adjust == "none") pairs$p_value else stats::p.adjust(pairs$p_value, method = adjust)
  pairs$post_prob <- NA_real_; pairs$bf10 <- NA_real_; pairs$evidence <- NA_character_
  instant <- do.call(rbind, ins)
  instant$p_adj <- if (adjust == "none") instant$p_value else stats::p.adjust(instant$p_value, method = adjust)
  instant$significant <- instant$p_adj < alpha
  pairs$significant <- FALSE
  for (i in seq_len(nrow(up))) {
    r <- res[[i]]; a <- up$v1[i]; b <- up$v2[i]
    ja <- which(pairs$cause == a & pairs$effect == b); jb <- which(pairs$cause == b & pairs$effect == a)
    dec <- .markov_decide(r$p_x >= alpha, r$p_y >= alpha, r$p_j >= alpha,
                          pairs$p_adj[ja] < alpha, pairs$p_adj[jb] < alpha,
                          instant$significant[i], iid)
    pairs$significant[ja] <- dec$ab != "n"; pairs$significant[jb] <- dec$ba != "n"
    mk[[i]] <- data.frame(var1 = a, var2 = b, memory = r$m, p_markov_x = r$p_x, p_markov_y = r$p_y,
                          p_markov_joint = r$p_j, p_joint_x_earlier = r$p_jm1, p_joint_x_later = r$p_jp1,
                          case = dec$case, inj_x_to_y = dec$ab, inj_y_to_x = dec$ba,
                          hidden = dec$hidden, consistent = dec$consistent, stringsAsFactors = FALSE)
  }
  pairs <- pairs[, .pairs_cols]
  rownames(pairs) <- NULL
  out <- .new_causres(pairs, series = names(all), causes = causes, effects = effects,
                      method = "markov", conditional = FALSE, alpha = alpha, bf_threshold = NA_real_,
                      adjust = adjust, n = nrow(all), call = cl, data = all, instant = instant,
                      markov = do.call(rbind, mk))
  out$null_markov <- stats::setNames(lapply(res, `[[`, "null"), paste(up$v1, "~", up$v2))
  out$settings$iid <- iid
  out
}

# df (degrees of freedom) method driver -------------------------------------
.rcausality_df <- function(all, grid, causes, effects, alpha, adjust, conditional, lags, history,
                           k, k_perm, n_perm, max_constr, cl, cores) {
  n <- nrow(all)
  sc <- .split_cores(nrow(grid), cores); pair_cores <- sc$outer; inner_cores <- sc$inner
  res <- .plapply(seq_len(nrow(grid)), function(i) {
    cs <- grid$cause[i]; ef <- grid$effect[i]
    ctrl_names <- if (conditional) setdiff(names(all), c(cs, ef)) else character(0)
    g <- .df_pair(all[[ef]], all[[cs]],
                  if (length(ctrl_names)) as.matrix(all[, ctrl_names, drop = FALSE]) else NULL,
                  lags, history[i], k, n_perm, k_perm, cores = inner_cores)
    list(row = data.frame(cause = cs, effect = ef, lag = g$lag, statistic = g$statistic,
                          df1 = NA_real_, df2 = NA_real_, p_value = g$p_value, bf10 = NA_real_,
                          effect_size = g$statistic, r2_gain = g$statistic, sign = NA_real_,
                          n_eff = g$n_eff, stringsAsFactors = FALSE),
         by_lag = cbind(cause = cs, effect = ef, g$by_lag),
         null = list(statistic = g$statistic, null = g$null, alternative = "greater", cond_r2 = g$cond_r2))
  }, cores = pair_cores)
  pairs <- do.call(rbind, lapply(res, `[[`, "row"))
  by_lag <- do.call(rbind, lapply(res, `[[`, "by_lag")); rownames(by_lag) <- NULL
  pairs$p_adj <- if (adjust == "none") pairs$p_value else stats::p.adjust(pairs$p_value, method = adjust)
  pairs$post_prob <- NA_real_; pairs$evidence <- NA_character_
  pairs$significant <- pairs$p_adj < alpha
  # orders of x, y and the joint (Table 2 of the paper), per unordered pair
  up <- .unordered_pairs(grid, all)
  t_idx <- (max_constr + 1):n
  n_perm_o <- min(n_perm, 49L)             # the order tests only need alpha-level resolution
  ord <- .plapply(seq_len(nrow(up)), function(i) {
    x <- all[[up$v1[i]]]; y <- all[[up$v2[i]]]
    xs <- x / stats::sd(x); ys <- y / stats::sd(y)
    cands <- stats::setNames(list(xs, ys), c(up$v1[i], up$v2[i]))
    ox <- .df_order(cbind(xs[t_idx]), cands, t_idx, max_constr, k, alpha, 0.02, n_perm_o, k_perm, inner_cores)
    oy <- .df_order(cbind(ys[t_idx]), cands, t_idx, max_constr, k, alpha, 0.02, n_perm_o, k_perm, inner_cores)
    oj <- .df_order(cbind(xs[t_idx], ys[t_idx]), cands, t_idx, max_constr, k, alpha, 0.02, n_perm_o, k_perm, inner_cores)
    uses <- function(o, nm) any(startsWith(o$chosen, paste0(nm, "[")))
    data.frame(var1 = up$v1[i], var2 = up$v2[i], df_x = ox$df, df_y = oy$df, df_joint = oj$df,
               conditions_x = paste(ox$chosen, collapse = " "), conditions_y = paste(oy$chosen, collapse = " "),
               conditions_joint = paste(oj$chosen, collapse = " "),
               table2 = .df_decide(ox$df, oy$df, oj$df, x_drives_y = uses(oy, up$v1[i]), y_drives_x = uses(ox, up$v2[i])),
               stringsAsFactors = FALSE)
  }, cores = pair_cores)
  dfo <- do.call(rbind, ord)
  instant <- data.frame(var1 = dfo$var1, var2 = dfo$var2, cmi = NA_real_, p_value = NA_real_,
                        p_adj = NA_real_, significant = dfo$table2 == "common_cause",
                        stringsAsFactors = FALSE)
  # Table 2 of the paper decides a pair whenever the orders give a definite case
  # (in a deterministic system the driven series reconstructs its driver, so
  # the directed tests alone read a one-way link as bidirectional).  The order
  # tests carry no multiple-testing control, so a direction claimed by the
  # table is kept only when its directed test is at least nominally
  # significant (raw p < alpha); undetermined pairs are decided by the
  # directed tests alone.
  dfo$decided_by <- ifelse(dfo$table2 == "undetermined", "directed tests", "Table 2")
  for (i in seq_len(nrow(dfo))) {
    if (dfo$table2[i] == "undetermined") next
    a <- dfo$var1[i]; b <- dfo$var2[i]
    ab <- pairs$cause == a & pairs$effect == b; ba <- pairs$cause == b & pairs$effect == a
    pairs$significant[ab] <- dfo$table2[i] %in% c("x_to_y", "bidirectional") & pairs$p_value[ab] < alpha
    pairs$significant[ba] <- dfo$table2[i] %in% c("y_to_x", "bidirectional") & pairs$p_value[ba] < alpha
  }
  pairs <- pairs[, .pairs_cols]
  rownames(pairs) <- NULL
  out <- .new_causres(pairs, series = names(all), causes = causes, effects = effects,
                      method = "df", conditional = conditional, alpha = alpha, bf_threshold = NA_real_,
                      adjust = adjust, n = n, call = cl, data = all, instant = instant,
                      by_lag = by_lag, df = dfo)
  out$null <- stats::setNames(lapply(res, `[[`, "null"), paste(pairs$cause, "->", pairs$effect))
  out
}
