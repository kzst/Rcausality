# PCMCI wrapper (optional; needs reticulate + the Python package tigramite) ----

#' Is the PCMCI back-end (Python tigramite) available?
#'
#' @param python optional path to a Python executable or the name of a
#'   virtualenv/conda environment passed to `reticulate::use_python()` /
#'   `use_virtualenv()`.
#' @return `TRUE` if reticulate can import `tigramite`, otherwise `FALSE`.
#' @export
#' @examples
#' pcmci_available()
pcmci_available <- function(python = NULL) {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(FALSE)
  if (!is.null(python)) try(reticulate::use_python(python, required = FALSE), silent = TRUE)
  isTRUE(tryCatch(reticulate::py_module_available("tigramite"), error = function(e) FALSE))
}

#' PCMCI / PCMCI+ causal discovery (Runge et al., 2019) via tigramite
#'
#' A thin wrapper around the Python package *tigramite* (Runge et al., 2019,
#' 2020) that runs PCMCI (or PCMCI+) on a multivariate series and returns the
#' result as a `causres` object, so that it prints, plots and compares like
#' the native methods. The Python side is reached through **reticulate**; both
#' reticulate and `tigramite` (`pip install tigramite`) must be installed
#' (see [pcmci_available()]). Nothing is ported: tigramite does the work.
#'
#' PCMCI first runs the PC\eqn{_1} condition-selection stage for every
#' variable (parents are selected with significance `pc_alpha`), then the MCI
#' (momentary conditional independence) stage tests every lagged link
#' \eqn{X^i_{t-\tau} \to X^j_t} conditioned on the parents of both variables.
#' PCMCI+ (`plus = TRUE`) additionally orients contemporaneous links.
#'
#' @param x a data frame / matrix with the series in columns (at least two).
#' @param tau_max,tau_min maximal and minimal lag.
#' @param test conditional-independence test: `"parcorr"` (linear partial
#'   correlation, analytic p-values), `"robust_parcorr"` (rank-transformed),
#'   `"cmiknn"` (kNN conditional mutual information, permutation p-values;
#'   slow) or `"gpdc"` (Gaussian-process / distance correlation).
#' @param pc_alpha significance level(s) of the condition-selection stage;
#'   `NULL` lets tigramite pick it by model selection.
#' @param alpha significance level of the final (adjusted) MCI p-values.
#' @param adjust multiple-testing correction applied by tigramite to the MCI
#'   p-values: `"fdr_bh"` (Benjamini–Hochberg, tigramite's default) or `"none"`.
#' @param plus run PCMCI+ (`run_pcmciplus`) instead of PCMCI.
#' @param python see [pcmci_available()].
#' @param verbosity tigramite verbosity.
#' @return A `causres` object (method `"pcmci"`). `pairs` holds, for every
#'   ordered pair, the lag with the smallest MCI p-value; `by_lag` holds all
#'   lags (`statistic` = test statistic, `p` = raw MCI p-value, `q` =
#'   adjusted); `instant` holds the contemporaneous links (lag 0) when
#'   `plus = TRUE` or `tau_min = 0`; `pcmci` keeps tigramite's raw
#'   `val_matrix`, `p_matrix`, `q_matrix` and `graph`.
#' @references
#' Runge, J., Nowack, P., Kretschmer, M., Flaxman, S., & Sejdinovic, D. (2019). Detecting and quantifying causal associations in large nonlinear time series datasets. *Science Advances*, 5(11), eaau4996. \doi{10.1126/sciadv.aau4996}
#'
#' Runge, J. (2020). Discovering contemporaneous and lagged causal relations in autocorrelated nonlinear time series datasets. *Proceedings of the 36th Conference on Uncertainty in Artificial Intelligence (UAI)*, PMLR 124, 1388-1397.
#' @export
#' @examples
#' \donttest{
#' if (pcmci_available()) {           # needs reticulate + the Python package tigramite
#'   d <- sim_var(500, n_series = 3, seed = 1)
#'   r <- pcmci(d, tau_max = 3)
#'   r
#'   plot(r)
#' }
#' }
pcmci <- function(x, tau_max = 3L, tau_min = 1L, test = c("parcorr", "robust_parcorr", "cmiknn", "gpdc"),
                  pc_alpha = NULL, alpha = 0.05, adjust = c("fdr_bh", "none"), plus = FALSE,
                  python = NULL, verbosity = 0L) {
  test <- match.arg(test); adjust <- match.arg(adjust)
  if (!pcmci_available(python))
    stop("pcmci() needs the R package 'reticulate' and the Python package 'tigramite' ",
         "(pip install tigramite). See ?pcmci_available.", call. = FALSE)
  fr <- .as_series_frame(x, NULL, xname = "x", yname = "y")
  all <- fr$X
  if (ncol(all) < 2L) stop("Supply a data frame with at least two columns.", call. = FALSE)
  all[] <- lapply(all, function(v) .prep_na(as.numeric(v)))
  cl <- match.call()
  tg_pp <- reticulate::import("tigramite.data_processing", convert = FALSE)
  tg_pc <- reticulate::import("tigramite.pcmci", convert = FALSE)
  ci <- switch(test,
    parcorr = reticulate::import("tigramite.independence_tests.parcorr")$ParCorr(significance = "analytic"),
    robust_parcorr = reticulate::import("tigramite.independence_tests.robust_parcorr")$RobustParCorr(significance = "analytic"),
    cmiknn = reticulate::import("tigramite.independence_tests.cmiknn")$CMIknn(significance = "shuffle_test"),
    gpdc = reticulate::import("tigramite.independence_tests.gpdc")$GPDC(significance = "analytic"))
  np <- reticulate::import("numpy", convert = FALSE)
  dataframe <- tg_pp$DataFrame(np$array(as.matrix(all)), var_names = as.list(names(all)))
  pc <- tg_pc$PCMCI(dataframe = dataframe, cond_ind_test = ci, verbosity = as.integer(verbosity))
  pca <- if (is.null(pc_alpha)) NULL else pc_alpha
  res <- if (plus) pc$run_pcmciplus(tau_min = as.integer(tau_min), tau_max = as.integer(tau_max), pc_alpha = pca)
         else pc$run_pcmci(tau_min = as.integer(tau_min), tau_max = as.integer(tau_max), pc_alpha = pca)
  res <- reticulate::py_to_r(res)
  val <- res$val_matrix; pm <- res$p_matrix
  qm <- if (adjust == "none") pm else
    reticulate::py_to_r(pc$get_corrected_pvalues(p_matrix = pm, tau_min = as.integer(tau_min),
                                                 tau_max = as.integer(tau_max), fdr_method = adjust))
  vars <- names(all); N <- length(vars)
  # lagged links: val[i, j, tau] = X^i_{t-tau} -> X^j_t
  rows <- list(); bl <- list()
  for (i in seq_len(N)) for (j in seq_len(N)) {
    if (i == j) next
    taus <- seq.int(max(1L, tau_min), tau_max)
    b <- data.frame(cause = vars[i], effect = vars[j], lag = taus, statistic = val[i, j, taus + 1L],
                    p = pm[i, j, taus + 1L], q = qm[i, j, taus + 1L], stringsAsFactors = FALSE)
    bl[[length(bl) + 1L]] <- b
    k <- which.min(b$q + 1e-12 * b$p)
    rows[[length(rows) + 1L]] <- data.frame(cause = vars[i], effect = vars[j], lag = b$lag[k], statistic = b$statistic[k],
      df1 = NA_real_, df2 = NA_real_, p_value = b$p[k], p_adj = b$q[k], bf10 = NA_real_, post_prob = NA_real_,
      evidence = NA_character_, effect_size = b$statistic[k], r2_gain = if (test %in% c("parcorr", "robust_parcorr")) b$statistic[k]^2 else NA_real_,
      sign = sign(b$statistic[k]), n_eff = nrow(all) - tau_max, significant = b$q[k] < alpha, stringsAsFactors = FALSE)
  }
  pairs <- do.call(rbind, rows); by_lag <- do.call(rbind, bl); rownames(pairs) <- rownames(by_lag) <- NULL
  instant <- NULL
  if (tau_min == 0L || plus) {
    ir <- list()
    for (i in seq_len(N - 1L)) for (j in (i + 1L):N) {
      pij <- min(qm[i, j, 1L], qm[j, i, 1L])
      ir[[length(ir) + 1L]] <- data.frame(var1 = vars[i], var2 = vars[j], statistic = val[i, j, 1L],
                                          p_value = min(pm[i, j, 1L], pm[j, i, 1L]), p_adj = pij,
                                          significant = pij < alpha, stringsAsFactors = FALSE)
    }
    instant <- do.call(rbind, ir)
  }
  out <- .new_causres(pairs, series = vars, causes = vars, effects = vars, method = "pcmci", conditional = TRUE,
                      alpha = alpha, bf_threshold = NA_real_, adjust = adjust, n = nrow(all), call = cl, data = all,
                      instant = instant, by_lag = by_lag)
  out$pcmci <- list(val_matrix = val, p_matrix = pm, q_matrix = qm, graph = res$graph, test = test,
                    plus = plus, tau_max = tau_max, tau_min = tau_min, pc_alpha = pc_alpha)
  out$settings$test <- test; out$settings$plus <- plus
  out
}
