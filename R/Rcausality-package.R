#' Rcausality: causality analysis for time series and panel data
#'
#' One entry point, [rcausality()], for a wide range of causality tests on two
#' series, all pairs of a data frame, two blocks of series, a panel of units
#' observed over time (`id`), or the causal drivers of one explained variable
#' (`target`), returning a [`causres`][causres-class] object with `print`,
#' `summary` and `plot` methods. The method (`method`), its Bayesian variant
#' (`bayes`) and the number of CPU cores (`cores`) are options of every call;
#' [causal_methods()] lists the methods and their settings.
#'
#' @section Families of methods:
#' * Granger causality and its generalisations: the F-test (`"granger"`), the
#'   Bayesian regression comparison with spike-and-slab, Minnesota (short
#'   series), g, hyper-g or BIC priors (`"bayes_granger"`, [bayes_granger()],
#'   [bayes_control()]), the Toda-Yamamoto test for integrated series
#'   (`"toda_yamamoto"`), the frequency-domain decomposition
#'   ([granger_spectral()]), the nonparametric Diks-Panchenko test
#'   (`"diks_panchenko"`, [diks_panchenko()]) and the random-forest
#'   Granger-like test (`"rf_granger"`, [rf_granger()]).
#' * Information-theoretic and model-free methods: conditional mutual
#'   information (`"cmi"`, `"cmi_lag"`), transfer entropy (`"te"`), MAR-PID
#'   (`"pid"`), the Markov-chain approach (`"markov"`), degrees of freedom
#'   (`"df"`), dimensional causality (`"dc"`), CCDH (`"ccdh"`), df-entropy
#'   (`"dfe"`), convergent cross mapping (`"ccm"`), recurrence plots
#'   (`"hirata"`), correlation dimensions (`"cd"`) and PCMCI (`"pcmci"`,
#'   optional, through Python).
#' * Every method has a Bayesian variant (`bayes = TRUE`).
#'
#' @section Panel data and one explained variable:
#' [panel_causality()] (every method at the level of the units, of the group
#' of units and of the pooled units, with Bayesian variants;
#' [panel_control()], [plot.causpanel()], [sim_panel()], [panel_long()]) and
#' [causal_drivers()] (drivers of a target with lags, importance, direct and
#' indirect links and the causal graph of the explanatory variables, for
#' time series and panels; [plot.causdrivers()], [sim_drivers()]), following
#' Kosztyan, Matrai and Kiss (2026).
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
#' methods paper by paper), `vignette("extensions")` (generalised,
#' nonlinear and Bayesian tests, robustness and forecasting) and
#' `vignette("panel")` (panel data and the causal drivers of one variable).
#'
#' @references
#' Kosztyan, Z. T., Matrai, R. and Kiss, D. (2026). Is the world getting smaller? Network-driven heterogeneity in gravity-model elasticities. *The World Economy*. \doi{10.1111/twec.70148}
#'
#' @keywords internal
"_PACKAGE"
