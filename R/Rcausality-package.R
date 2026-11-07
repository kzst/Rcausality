#' Rcausality: causality analysis for time series
#'
#' One entry point, [rcausality()], for a wide range of causality tests on two
#' series, all pairs of a data frame or two blocks of series, returning a
#' [`causres`][causres-class] object with `print`, `summary` and `plot`
#' methods.
#'
#' @section Families of methods:
#' * Granger causality and its generalisations: the F-test (`"granger"`), the
#'   Bayesian regression comparison with spike-and-slab, Minnesota (short
#'   series), g, hyper-g or BIC priors (`"bayes_granger"`, [bayes_granger()],
#'   [bayes_control()]), the Toda-Yamamoto test for integrated series
#'   (`"toda_yamamoto"`), the frequency-domain decomposition
#'   ([granger_spectral()]), the nonparametric Diks-Panchenko test
#'   (`"diks_panchenko"`, [diks_panchenko()]) and the random-forest
#'   Granger-like test (`"rf_granger"`, [rf_granger()]); panels of short
#'   series: [panel_granger()].
#' * Information-theoretic and model-free methods: conditional mutual
#'   information (`"cmi"`, `"cmi_lag"`), transfer entropy (`"te"`), MAR-PID
#'   (`"pid"`), the Markov-chain approach (`"markov"`), degrees of freedom
#'   (`"df"`), dimensional causality (`"dc"`), CCDH (`"ccdh"`), df-entropy
#'   (`"dfe"`), convergent cross mapping (`"ccm"`), recurrence plots
#'   (`"hirata"`), correlation dimensions (`"cd"`) and PCMCI (`"pcmci"`,
#'   optional, through Python).
#' * Every method has a Bayesian variant (`bayes = TRUE`).
#'
#' @section Checking and using the links:
#' [placebo_test()] (placebo data and time reversal), [robustness()]
#' (settings, subsamples, methods), [compare_methods()] and
#' [benchmark_methods()] (triangulation and simulation benchmarks),
#' [predict.causres()], [forecast_causal()] and [forecast_compare()]
#' (forecasting with the detected causes).
#'
#' @section Parallel computation:
#' All repeated work runs in parallel through the `cores` argument with
#' results that do not depend on the number of cores; see
#' [Rcausality-parallel].
#'
#' @section Vignettes:
#' `vignette("rcausality")` (getting started), `vignette("methods")` (the
#' methods paper by paper) and `vignette("extensions")` (generalised,
#' nonlinear, Bayesian and panel tests, robustness and forecasting).
#'
#' @keywords internal
"_PACKAGE"
