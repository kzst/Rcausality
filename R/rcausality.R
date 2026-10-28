#' Detect causal relations between time series
#'
#' Single entry point of the package. Works with
#' * two series: `rcausality(x, y)`;
#' * one data frame / matrix: `rcausality(df)` -- all ordered pairs;
#' * two blocks (MIMO): `rcausality(X, Y)` with several columns each --
#'   every column of `X` is tested as a cause of every column of `Y`.
#'
#' @param x numeric vector, matrix or data frame (the putative causes).
#' @param y optional numeric vector, matrix or data frame (the putative
#'   effects). If `NULL` and `x` has several columns, all ordered pairs of the
#'   columns of `x` are tested.
#' @param method `"granger"` (F-test, default), `"bayes_granger"` (the
#'   Bayesian Granger regression, i.e. `bayes = TRUE`: decision by the Bayes
#'   factor of the spike-and-slab regression by default, or by the
#'   closed-form BIC / g / hyper-g factor with `prior = "bic"`, `"g"`,
#'   `"hyper_g"`; p-value *and* Bayes factor are both reported) or `"cmi"`
#'   (model-free: kNN conditional mutual information
#'   `I(x[t-lag]; y[t] | history of y, controls)` with a local-permutation
#'   p-value; detects non-linear dependence) or `"cmi_lag"` (the
#'   lag-shifted CMI scheme of Zlatniczki et al. (2021): like `"cmi"` but
#'   also conditioning on the cause's own past, plus a contemporaneous
#'   test `I(x[t]; y[t] | pasts)` whose significance, without lagged links,
#'   indicates a *common driver*; `common_cause` class available).
#' @param lag integer lag order or `"auto"` (selected by `lag_criterion` on
#'   the VAR of the involved series).
#' @param max_lag largest lag order considered when `lag = "auto"`.
#' @param lag_criterion `"BIC"`, `"AIC"` or `"HQ"`.
#' @param alpha significance level (applied to adjusted p-values).
#' @param bf_threshold Bayes factor above which a link is declared
#'   (used for `method = "bayes_granger"`).
#' @param adjust multiple-testing correction passed to
#'   [stats::p.adjust()] (default `"holm"`; `"none"` disables it).
#' @param conditional logical; if `TRUE` every test conditions on all other
#'   series in the data (full conditional / multivariate Granger), otherwise
#'   tests are pairwise.
#' @param directions for two series / two blocks: `"both"` (default) tests
#'   x -> y and y -> x; `"forward"` only the causes-in-`x`, effects-in-`y`
#'   direction. Ignored when `y` is `NULL` (all ordered pairs are tested).
#' @param k,k_perm,n_perm,history,seed settings of `method = "cmi"`: number
#'   of neighbours of the estimator (`NULL` = 5; 10 for `method = "df"`),
#'   size of the local-permutation neighbourhood, number of permutations,
#'   order of the conditioning history of the effect series (`"auto"` =
#'   the VAR order of the pair selected by `lag_criterion`, at most 5; too
#'   short a history lets the cause's lags stand in for the effect's own
#'   deeper past and produces spurious reverse links), and RNG seed. With
#'   `lag = "auto"` every lag `1..max_lag` is tested at once: the reported
#'   lag is the one with the largest statistic and the p-value is that of
#'   the maximal statistic against the maximum of the permuted statistics
#'   over the lags (the same max-statistic test is used by `"te"`, `"pid"`
#'   and `"df"`). Series are standardised before the neighbour searches.
#' @param m,tau,k_range,prep,prob_threshold settings of `method = "dc"`
#'   (Bayesian dimensional causality of Benko et al., 2024, see
#'   [dimensional_causality()]): embedding dimension and delay,
#'   neighbourhood sizes averaged over, preprocessing, and the posterior
#'   probability above which a directed link (direct + bidirectional) is
#'   reported as significant. Pairwise only; the method assumes
#'   deterministic dynamics (no p-values, decisions by posterior
#'   probabilities; `bf10` is the posterior odds against the uniform prior).
#' @param rate `method = "hirata"` (recurrence-plot criteria of Hirata and
#'   Aihara, 2010, see [hirata_rp()]): recurrence rate; `m`/`tau` give the
#'   embedding, `n_perm` the surrogates, `theta` the coverage-index
#'   threshold. `method = "cd"` (correlation-dimension relations of
#'   Krakovska, 2019, see [krakovska_cd()]) uses `m`/`tau` and bootstrap
#'   z-tests at level `alpha`. Both pairwise; verdicts are per pair.
#' @param k_target,k_source `method = "pcmci"` (PCMCI of Runge et al., 2019,
#'   via the Python package tigramite, see [pcmci()]; needs reticulate +
#'   tigramite): linear partial-correlation MCI tests over lags 1..`max_lag`
#'   (or the fixed `lag`), Benjamini-Hochberg adjusted unless
#'   `adjust = "none"`; always multivariate (conditional on the other
#'   series). `method = "pid"` (finite-resolution MAR-PID,
#'   see [marpid()]): binarisation depths of the effect's future (target)
#'   and of the two sources (the cause's past at the candidate lag and the
#'   effect's own past); the directed statistic is the unique + synergistic
#'   ("novel") information of the cause's past, tested with `n_perm`
#'   surrogates. Pairwise only.
#' @param bins `method = "te"` (binned transfer entropy, see
#'   [transfer_entropy()]): number of quantile bins; `history` is the effect's
#'   own history, the cause's history is scanned over the candidate lags,
#'   `n_perm` circular-shift surrogates give the p-values. `method = "ccm"`
#'   (convergent cross mapping, see [ccm()]) uses `m`/`tau` as embedding and
#'   `n_perm` surrogates; a link needs a significant and converging cross-map
#'   skill. Both are pairwise.
#' @param model `method = "dfe"` (df-entropy method with random-forest
#'   classification, Jakovac et al.; see [dfe_classify()]): a trained
#'   `dfe_model` from [dfe_train()]; `NULL` trains and caches a default one.
#'   Decisions by class scores as for `"dc"` (`prob_threshold`); a numeric
#'   `lag` shifts the pairs. Pairwise only.
#' @param K,n_surr,theta settings of `method = "ccdh"` (complexity-based
#'   causal discovery of Stippinger et al., see [ccdh()]): number of NSRPS
#'   steps, number of circular-shift surrogates and the reconstruction-index
#'   threshold; `m` and `tau` give the embedding, a numeric `lag` the
#'   alignment of the pair (`"auto"` = 0 for this method). Pairwise only.
#' @param max_constr `method = "df"` (degrees-of-freedom / conditional-variance
#'   method of Telcs, Kurbucz and Jakovac, 2025): largest number of past lags
#'   constrained when estimating the orders of x, y and of the joint process;
#'   `k` is the neighbourhood size of the adaptive windows (default 10 for
#'   this method), `history` the number of own-past lags conditioned on
#'   in the directed tests, `n_perm`/`k_perm` the local-permutation settings.
#'   The orders (`$df`) give the Table-2 reading of the paper (including the
#'   common-driver case), which decides a pair whenever it is determinate;
#'   pairs the table leaves undetermined are decided by the directed
#'   variance-reduction tests (`$pairs`), which are always reported.
#' @param memory,iid_injection,k_nn,n_bins,n_markov settings of
#'   `method = "markov"` (Markov-chain approach of Stippinger et al., 2023):
#'   memory (embedding dimension; `"auto"` from the PACF), whether i.i.d.
#'   injection of one series into the other is `"allowed"` (default, Occam's
#'   razor) or `"forbidden"` in the decision table, local-neighbourhood size,
#'   number of bins of the local contingency tables and number of local
#'   permutations of the Markov tests. Pairwise only (`conditional` is not
#'   supported). Lags are fixed to 1 in this version.
#' @param cores number of CPU cores (`1` = serial, an integer, `"auto"` = all
#'   physical cores but one, or a cluster from `parallel::makeCluster()`);
#'   default `getOption("Rcausality.cores", 1L)`. Pairs are distributed over
#'   the cores when there are at least as many pairs as cores, otherwise the
#'   permutations / surrogates / chains of each test are. Every unit of work
#'   has its own random stream, so the results are *identical* for any number
#'   of cores given the `seed` (see [Rcausality-parallel]).
#' @param prior_prob prior probability of a causal link (for the posterior).
#' @param bayes logical; the Bayesian variant of the method (settings in
#'   `control`, see [bayes_control()]). `"granger"` becomes a Bayesian
#'   regression: closed-form Bayes factors (`prior = "bic"`, `"g"`,
#'   `"hyper_g"`) or the spike-and-slab / normal-slab sampler of
#'   [bayes_granger()] with posterior inclusion probabilities and credible
#'   intervals of Geweke's measure. The surrogate / permutation methods
#'   (`cmi`, `cmi_lag`, `df`, `te`, `ccm`, `pid`, `ccdh`, `hirata`, `cd`)
#'   compare the observed statistic with the surrogate sample through
#'   [surrogate_bf()] (an explicit effect-size prior); `markov` combines the
#'   surrogate Bayes factors of its six tests in a Bayesian version of the
#'   decision table; `dc` and `dfe` apply the prior over the causal cases;
#'   `pcmci` uses a Savage-Dickey Bayes factor of the Fisher-transformed
#'   partial correlation. `bf10`, `post_prob` and `evidence` are filled for
#'   every directed test, the per-pair details are in `$bayes`, and a link
#'   is declared when `bf10 > bf_threshold` (or, with
#'   `control$decision = "posterior"`, when `post_prob > prob_threshold`).
#' @param prior shortcut for `control$prior` (the prior of the Granger
#'   regression comparison; see [bayes_control()]); `NULL` = the control's.
#' @param control a [bayes_control()] list with the priors, the Monte Carlo
#'   and MCMC sizes and the decision rule of the Bayesian variants.
#' @param impute logical; interpolate missing values?
#' @return An object of class [`causres`][causres-class].
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
#' Zlatniczki, A., Stippinger, M., Benko, Z., Somogyvari, Z. and Telcs, A. (2021). Relaxation of some confusions about confounders. *Entropy*, 23(11), 1450. \doi{10.3390/e23111450}
#'
#' Stippinger, M., Bencze, A., Zlatniczki, A., Somogyvari, Z. and Telcs, A. (2023). Causal discovery of stochastic dynamical systems: a Markov chain approach. *Mathematics*, 11(4), 852. \doi{10.3390/math11040852}
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
#' # --- all pairs of a data frame, decision by the Bayes factor ------------
#' # (closed-form BIC Bayes factor; without `prior` the spike-and-slab
#' # regression of bayes_granger() is sampled for every pair)
#' d3 <- sim_var(400, n_series = 3, seed = 2)
#' rcausality(d3, method = "bayes_granger", prior = "bic")$classes
#' \donttest{
#' summary(r)
#' plot(r, "all")                    # network, p-value matrix and lag scan
#' rcausality(d3, method = "bayes_granger", control = bayes_control(n_iter = 300),
#'            seed = 1)$bayes                         # spike-and-slab: BF10, effect interval, R-hat
#'
#' # --- model-free methods on a non-linear coupling --------------------------
#' # (a second or two each; `cores = 2` halves the time, see ?Rcausality-parallel)
#' lg <- sim_logistic(200, "x_to_y", seed = 2)
#' rcausality(lg, method = "cmi", max_lag = 2, n_perm = 50, seed = 1)$classes   # kNN conditional MI
#' rcausality(lg, method = "te", n_perm = 100, seed = 1)$classes$relation  # binned transfer entropy
#' rcausality(lg, method = "ccm", n_perm = 50, seed = 1)$classes$relation  # convergent cross mapping
#' rcausality(sim_logistic_bdc(500, "x_to_y", seed = 1), method = "dc")$classes$relation
#'
#' # --- conditional (multivariate) Granger causality and blocks of series ---
#' rcausality(d3, conditional = TRUE)$classes$relation
#' rcausality(d3[, 1:2], d3[, 3, drop = FALSE], directions = "forward")
#'
#' # --- the Bayesian variant of a surrogate method ---------------------------
#' rcausality(lg, method = "te", bayes = TRUE, n_perm = 100, seed = 1)$bayes
#' }
rcausality <- function(x, y = NULL, method = c("granger", "bayes_granger", "cmi", "cmi_lag", "markov", "df", "dc", "ccdh", "dfe",
                                  "te", "ccm", "pid", "hirata", "cd", "pcmci"),
                       lag = "auto", max_lag = 8L,
                       lag_criterion = c("BIC", "AIC", "HQ"),
                       alpha = 0.05, bf_threshold = 3,
                       adjust = "holm", conditional = FALSE,
                       prior_prob = 0.5, impute = TRUE,
                       directions = c("both", "forward"),
                       k = NULL, k_perm = 5L, n_perm = 200L, history = "auto",
                       seed = NULL, memory = "auto", iid_injection = c("allowed", "forbidden"),
                       k_nn = NULL, n_bins = 2L, n_markov = 99L, cores = getOption("Rcausality.cores", 1L),
                       max_constr = 4L, m = 3L, tau = 1L, k_range = seq(10L, 40L, by = 5L),
                       prep = c("rank", "scale", "none"), prob_threshold = 0.5,
                       K = 10L, n_surr = 50L, theta = NULL, model = NULL, bins = 4L,
                       k_target = 2L, k_source = 2L, rate = 0.05,
                       bayes = FALSE, prior = NULL, control = bayes_control()) {
  cl <- match.call()
  method <- match.arg(method)
  if (!inherits(control, "bayes_control")) control <- do.call(bayes_control, as.list(control))
  if (!is.null(prior)) control$prior <- match.arg(prior, c("bic", "g", "hyper_g", "spike_slab", "normal"))
  prior <- control$prior
  if (method == "bayes_granger") bayes <- TRUE
  sampler <- bayes && method %in% c("granger", "bayes_granger") && prior %in% c("spike_slab", "normal")
  if (bayes && method == "granger" && !sampler) method <- "bayes_granger"
  if (sampler) method <- "granger"                      # the F-test run supplies the lags; the sampler adds the evidence
  old_rng <- .rng_setup(seed); on.exit(.rng_restore(old_rng), add = TRUE)   # the seed also covers the Bayesian stage
  r <- .rcausality_run(x = x, y = y, method = method, lag = lag, max_lag = max_lag, lag_criterion = lag_criterion,
                       alpha = alpha, bf_threshold = bf_threshold, adjust = adjust, conditional = conditional,
                       prior_prob = prior_prob, impute = impute, directions = directions, k = k, k_perm = k_perm,
                       n_perm = n_perm, history = history, seed = seed, memory = memory, iid_injection = iid_injection,
                       k_nn = k_nn, n_bins = n_bins, n_markov = n_markov, cores = cores, max_constr = max_constr,
                       m = m, tau = tau, k_range = k_range, prep = prep, prob_threshold = prob_threshold, K = K,
                       n_surr = n_surr, theta = theta, model = model, bins = bins, k_target = k_target,
                       k_source = k_source, rate = rate, prior = if (prior %in% c("bic", "g", "hyper_g")) prior else "bic", cl = cl)
  r <- .restrict_to_grid(r)
  if (sampler) r <- .bayes_granger_pairs(r, control, prior_prob, bf_threshold, prob_threshold, conditional, cores)
  else if (bayes && method != "bayes_granger") r <- .bayes_postprocess(r, control, prior_prob, bf_threshold, prob_threshold)
  r$settings$bayes <- bayes
  if (bayes) { r$settings$prior <- prior; r$settings$bayes_control <- control; r$settings$decision <- control$decision }
  # the full argument list, so that the analysis can be repeated on a subset
  # (predict.causres(refit = TRUE), benchmark_methods())
  r$args <- mget(setdiff(names(formals(rcausality)), c("x", "y")))
  r$args$method <- if (sampler) "granger" else method; r$args$prior <- prior; r$args$bayes <- bayes; r$args$control <- control
  r
}

.rcausality_run <- function(x, y = NULL, method, lag, max_lag, lag_criterion, alpha, bf_threshold, adjust,
                            conditional, prior_prob, impute, directions, k, k_perm, n_perm, history, seed, memory,
                            iid_injection, k_nn, n_bins, n_markov, cores, max_constr, m, tau, k_range, prep,
                            prob_threshold, K, n_surr, theta, model, bins, k_target, k_source, rate, prior, cl) {
  lag_criterion <- match.arg(lag_criterion, c("BIC", "AIC", "HQ"))
  directions <- match.arg(directions, c("both", "forward"))
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
  all[] <- lapply(all, function(v) { v <- as.numeric(v); if (impute) .prep_na(v) else v })
  if (any(!is.finite(as.matrix(all)))) stop("Missing/non-finite values present (use impute = TRUE).", call. = FALSE)
  causes <- if (is.null(Y)) names(all) else names(X)
  effects <- if (is.null(Y)) names(all) else names(Y)
  grid <- expand.grid(cause = causes, effect = effects, stringsAsFactors = FALSE)
  if (!is.null(Y) && directions == "both") {
    grid <- rbind(grid, setNames(grid[, 2:1], c("cause", "effect")))
    causes <- effects <- names(all)
  }
  grid <- grid[grid$cause != grid$effect, , drop = FALSE]
  if (nrow(grid) == 0L) stop("No pairs to test.", call. = FALSE)

  .n_cores(cores)                                        # validates the argument
  old_rng <- .rng_setup(seed); on.exit(.rng_restore(old_rng), add = TRUE)
  sc <- .split_cores(nrow(grid), cores); pair_cores <- sc$outer; inner_cores <- sc$inner
  if (method %in% c("cmi", "cmi_lag", "markov", "df", "ccdh", "te", "ccm", "pid")) {
    need <- min(999L, ceiling(nrow(grid) / alpha))
    if (method == "ccdh") {
      if (n_surr < need) {
        message("ccdh: n_surr raised from ", n_surr, " to ", need,
                " so that the surrogate p-values can reach significance after multiple-testing correction.")
        n_surr <- need
      }
    } else if (n_perm < need) {
      message(method, ": n_perm raised from ", n_perm, " to ", need,
              " so that the permutation p-values can reach significance after ",
              "multiple-testing correction.")
      n_perm <- need
    }
  }
  fixed_lag <- if (identical(lag, "auto")) NULL else as.integer(lag)
  cand_lags <- if (is.null(fixed_lag)) seq_len(max(1L, as.integer(max_lag))) else fixed_lag
  # conditioning history of the effect: the VAR order of the pair when "auto"
  hist_vec <- if (identical(history, "auto")) {
    vapply(seq_len(nrow(grid)), function(i)
      max(1L, min(5L, as.integer(max_lag), select_lag(all[, c(grid$effect[i], grid$cause[i])], max_lag, lag_criterion))), 0L)
  } else rep(max(1L, as.integer(history)), nrow(grid))
  if (is.null(k)) k <- if (method == "df") 10L else 5L
  if (is.null(theta)) theta <- if (method == "hirata") 0.4 else 0.45
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
  if (method == "pid") {
    if (conditional) stop("method = \"pid\" supports pairwise analysis only.", call. = FALSE)
    return(.rcausality_pid(all, grid, causes, effects, alpha, adjust, lags = cand_lags, history = hist_vec,
                           k_target = as.integer(k_target), k_source = as.integer(k_source), n_perm = n_perm,
                           cl = cl, cores = cores))
  }
  if (method == "te") {
    if (conditional) stop("method = \"te\" supports pairwise analysis only.", call. = FALSE)
    return(.rcausality_te(all, grid, causes, effects, alpha, adjust, conditional, lags = cand_lags,
                          history = hist_vec, bins = as.integer(bins), n_perm = n_perm, cl = cl, cores = cores))
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
                     lags = cand_lags, history = hist_vec[i], k = k, k_perm = k_perm, n_perm = n_perm,
                     cause_past = (method == "cmi_lag"), cores = inner_cores)
      return(list(row = data.frame(cause = cs, effect = ef, lag = g$lag,
                                   statistic = g$statistic, df1 = NA_real_, df2 = NA_real_,
                                   p_value = g$p_value, bf10 = NA_real_,
                                   effect_size = g$statistic, r2_gain = NA_real_,
                                   sign = NA_real_, n_eff = g$n_eff, stringsAsFactors = FALSE),
                  by_lag = cbind(cause = cs, effect = ef, g$by_lag),
                  null = list(statistic = g$statistic, null = g$null, alternative = "greater")))
    }
    inv <- all[, c(ef, cs, ctrl_names), drop = FALSE]
    p <- if (is.null(fixed_lag)) select_lag(inv, max_lag, lag_criterion) else fixed_lag
    p <- max(1L, min(p, floor((nrow(all) - 2) / (2 + 2 * ncol(inv)))))
    ctrl <- if (length(ctrl_names)) as.matrix(all[, ctrl_names, drop = FALSE]) else NULL
    g <- .granger_pair(all[[ef]], all[[cs]], ctrl, p, prior = prior)
    # F statistic over all candidate lags (cheap) for the lag plot
    lags <- seq_len(max(1L, min(as.integer(max_lag), floor((nrow(all) - 2) / (2 + 2 * ncol(inv))))))
    bl <- do.call(rbind, lapply(lags, function(l) {
      gl <- .granger_pair(all[[ef]], all[[cs]], ctrl, l)
      data.frame(cause = cs, effect = ef, lag = l, statistic = gl$statistic, p = gl$p_value)
    }))
    list(row = data.frame(cause = cs, effect = ef, lag = p,
                          statistic = g$statistic, df1 = g$df1, df2 = g$df2,
                          p_value = g$p_value, bf10 = g$bf10,
                          effect_size = g$effect_size, r2_gain = g$r2_gain,
                          sign = g$sign, n_eff = g$n_eff,
                          stringsAsFactors = FALSE),
         by_lag = bl)
  }
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
                        hp, k, k_perm, n_perm, cores = inner_cores)
      list(row = data.frame(var1 = up$v1[i], var2 = up$v2[i], cmi = g$cmi, p_value = g$p_value, stringsAsFactors = FALSE),
           null = list(statistic = g$cmi, null = g$null, alternative = "greater"))
    }, cores = pair_cores)
    null_instant <- stats::setNames(lapply(ins, `[[`, "null"), paste(up$v1, "~", up$v2))
    instant <- do.call(rbind, lapply(ins, `[[`, "row"))
    instant$p_adj <- if (adjust == "none") instant$p_value else
      stats::p.adjust(instant$p_value, method = adjust)
    instant$significant <- instant$p_adj < alpha
  }
  pairs$p_adj <- if (adjust == "none") pairs$p_value else
    stats::p.adjust(pairs$p_value, method = adjust)
  pairs$post_prob <- (pairs$bf10 * prior_prob) /
    (pairs$bf10 * prior_prob + (1 - prior_prob))
  pairs$evidence <- as.character(bf_evidence(pairs$bf10))
  pairs$significant <- if (method %in% c("granger", "cmi", "cmi_lag")) {
    !is.na(pairs$p_adj) & pairs$p_adj < alpha
  } else {
    !is.na(pairs$bf10) & pairs$bf10 > bf_threshold
  }
  pairs <- pairs[, c("cause", "effect", "lag", "statistic", "df1", "df2", "p_value",
                     "p_adj", "bf10", "post_prob", "evidence", "effect_size",
                     "r2_gain", "sign", "n_eff", "significant")]
  rownames(pairs) <- NULL
  out <- .new_causres(pairs, series = names(all), causes = causes, effects = effects,
                      method = method, conditional = conditional, alpha = alpha,
                      bf_threshold = bf_threshold, adjust = adjust, n = nrow(all),
                      call = cl, data = all, instant = instant, by_lag = by_lag)
  if (method %in% c("cmi", "cmi_lag")) {
    out$null <- stats::setNames(lapply(res, `[[`, "null"), paste(pairs$cause, "->", pairs$effect))
    out$null_instant <- null_instant
  }
  out
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
  pairs <- pairs[, c("cause", "effect", "lag", "statistic", "df1", "df2", "p_value",
                     "p_adj", "bf10", "post_prob", "evidence", "effect_size",
                     "r2_gain", "sign", "n_eff", "significant")]
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
         null = list(statistic = g$statistic, null = g$null, alternative = "greater"))
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
  pairs <- pairs[, c("cause", "effect", "lag", "statistic", "df1", "df2", "p_value",
                     "p_adj", "bf10", "post_prob", "evidence", "effect_size",
                     "r2_gain", "sign", "n_eff", "significant")]
  rownames(pairs) <- NULL
  out <- .new_causres(pairs, series = names(all), causes = causes, effects = effects,
                      method = "df", conditional = conditional, alpha = alpha, bf_threshold = NA_real_,
                      adjust = adjust, n = n, call = cl, data = all, instant = instant,
                      by_lag = by_lag, df = dfo)
  out$null <- stats::setNames(lapply(res, `[[`, "null"), paste(pairs$cause, "->", pairs$effect))
  out
}
