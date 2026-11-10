# Method registry and the shared assembly of the `pairs` table ---------------
#
# One entry per method of rcausality(): the labels used by print() and plot(),
# the name of the statistic column, whether the method can condition on the
# other series, the family that decides how its Bayesian variant works
# ("regression": closed-form or sampled Bayes factors of nested regressions;
# "surrogate": the statistic against its surrogate / asymptotic null through
# surrogate_bf(); "posterior": posterior over the causal cases; "verdict":
# one verdict per pair, no directed p-values; "external": PCMCI), the tuning
# settings the method reads, whether a pooled (stacked) panel estimator
# exists, and the reference shown by causal_methods().

.method_registry <- list(
  granger        = list(label = "Granger causality (F-test)", short = "Granger causality (F-test)",
                        stat = "F", family = "regression", conditional = TRUE, type = "linear", pooled = TRUE, min_n = 12L,
                        desc = "F-test of the cause's lags in the autoregression of the effect",
                        ref = "Granger (1969)", tuning = character(0)),
  bayes_granger  = list(label = "Bayesian Granger causality", short = "Bayesian Granger causality (Bayes factor)",
                        stat = "F", family = "regression", conditional = TRUE, type = "linear", pooled = TRUE, min_n = 8L,
                        desc = "Bayes factor of the Granger regressions (spike-and-slab, Minnesota, g, hyper-g, BIC priors)",
                        ref = "George & McCulloch (1993); Giannone et al. (2015)", tuning = character(0)),
  toda_yamamoto  = list(label = "Toda-Yamamoto lag-augmented Granger causality (integrated or cointegrated series)",
                        short = "Toda-Yamamoto Granger causality", stat = "F", family = "regression", conditional = TRUE,
                        type = "linear", pooled = TRUE, min_n = 15L,
                        desc = "lag-augmented Wald test, valid for integrated and cointegrated series",
                        ref = "Toda & Yamamoto (1995)", tuning = character(0)),
  diks_panchenko = list(label = "Diks-Panchenko nonparametric (nonlinear) Granger causality test",
                        short = "Diks-Panchenko nonparametric test", stat = "T.n", family = "surrogate", conditional = TRUE,
                        type = "nonlinear", pooled = TRUE, min_n = 30L,
                        desc = "nonparametric test of conditional independence with kernel densities",
                        ref = "Diks & Panchenko (2006)", tuning = character(0)),
  rf_granger     = list(label = "Random-forest Granger-like causality (out-of-bag gain, circular-shift surrogates)",
                        short = "random-forest Granger-like causality", stat = "OOB.gain", family = "surrogate",
                        conditional = TRUE, type = "nonlinear", pooled = TRUE, min_n = 25L, needs = "ranger",
                        desc = "predictive gain of a random forest from the cause's past, surrogate test",
                        ref = "Kosztyan, Matrai & Kiss (2026)", tuning = character(0)),
  cmi            = list(label = "Conditional mutual information (kNN, local permutation)",
                        short = "kNN conditional mutual information", stat = "CMI", family = "surrogate", conditional = TRUE,
                        type = "nonlinear", pooled = TRUE, min_n = 50L,
                        desc = "kNN conditional mutual information of the cause's past and the effect, local permutation test",
                        ref = "Frenzel & Pompe (2007); Runge (2018)", tuning = c("k", "k_perm", "history")),
  cmi_lag        = list(label = "Lag-shifted conditional mutual information (Zlatniczki et al., 2021)",
                        short = "lag-shifted CMI (Zlatniczki et al. 2021)", stat = "CMI", family = "surrogate", conditional = TRUE,
                        type = "nonlinear", pooled = TRUE, min_n = 50L,
                        desc = "CMI conditioned on both pasts, with a contemporaneous test for common drivers",
                        ref = "Zlatniczki et al. (2021)", tuning = c("k", "k_perm", "history")),
  markov         = list(label = "Markov-chain approach (Stippinger et al., 2023)",
                        short = "Markov-chain approach (Stippinger et al. 2023)", stat = "CMI", family = "surrogate",
                        conditional = FALSE, type = "nonlinear", pooled = FALSE, min_n = 200L,
                        desc = "Markov-property tests of the marginal and joint processes with a decision table",
                        ref = "Stippinger et al. (2023)",
                        tuning = c("memory", "iid_injection", "k_nn", "n_bins", "n_markov", "k", "k_perm")),
  df             = list(label = "Degrees-of-freedom / conditional-variance method (Telcs, Kurbucz & Jakovac, 2025)",
                        short = "degrees of freedom (Telcs et al. 2025)", stat = "sd.red", family = "surrogate", conditional = TRUE,
                        type = "nonlinear", pooled = FALSE, min_n = 100L,
                        desc = "reduction of the conditional variance and the orders (degrees of freedom) of the processes",
                        ref = "Telcs, Kurbucz & Jakovac (2025)", tuning = c("k", "k_perm", "history", "max_constr")),
  dc             = list(label = "Bayesian dimensional causality (Benko et al., 2024)",
                        short = "dimensional causality (Benko et al. 2024)", stat = "P.direct", family = "posterior",
                        conditional = FALSE, type = "nonlinear", pooled = FALSE, min_n = 300L,
                        desc = "posterior of the causal cases from the intrinsic dimensions of the embedded manifolds",
                        ref = "Benko et al. (2024)", tuning = c("m", "tau", "k_range", "prep")),
  ccdh           = list(label = "CCDH: complexity-based causal discovery (Stippinger et al., 2023)",
                        short = "CCDH (Stippinger et al. 2023)", stat = "recon", family = "surrogate", conditional = FALSE,
                        type = "nonlinear", pooled = FALSE, min_n = 300L,
                        desc = "compression-complexity based reconstruction of the driver from the driven series",
                        ref = "Stippinger et al. (2023)", tuning = c("m", "tau", "K", "n_surr", "theta")),
  dfe            = list(label = "df-entropy method with random-forest classification (Jakovac et al.)",
                        short = "df-entropy + random forest (Jakovac et al.)", stat = "P.direct", family = "posterior",
                        conditional = FALSE, type = "nonlinear", pooled = FALSE, min_n = 500L, needs = "ranger",
                        desc = "entropy tables of the binarised pair classified by a trained random forest",
                        ref = "Jakovac et al.", tuning = "model"),
  te             = list(label = "Transfer entropy, binned (Schreiber, 2000) with surrogate test",
                        short = "transfer entropy (Schreiber 2000)", stat = "TE.bits", family = "surrogate", conditional = FALSE,
                        type = "nonlinear", pooled = TRUE, min_n = 30L,
                        desc = "binned transfer entropy with circular-shift surrogates",
                        ref = "Schreiber (2000)", tuning = c("bins", "history")),
  ccm            = list(label = "Convergent cross mapping (Sugihara et al., 2012)",
                        short = "convergent cross mapping (Sugihara et al. 2012)", stat = "rho", family = "surrogate",
                        conditional = FALSE, type = "nonlinear", pooled = TRUE, min_n = 50L,
                        desc = "cross-map skill of the cause from the effect's delay embedding, convergence and surrogates",
                        ref = "Sugihara et al. (2012); Clark et al. (2015)", tuning = c("m", "tau")),
  pid            = list(label = "MAR-PID: unique + synergistic information of the cause's past about the effect's future (Telcs, 2026)",
                        short = "MAR-PID (Telcs 2026)", stat = "U+S.bits", family = "surrogate", conditional = FALSE,
                        type = "nonlinear", pooled = FALSE, min_n = 200L,
                        desc = "partial information decomposition: novel information of the cause's past",
                        ref = "Telcs (2026)", tuning = c("k_target", "k_source", "history")),
  hirata         = list(label = "Recurrence-plot criteria: coverage and simultaneous recurrences (Hirata & Aihara, 2010)",
                        short = "recurrence plots (Hirata & Aihara 2010)", stat = "cover.idx", family = "verdict",
                        conditional = FALSE, type = "nonlinear", pooled = FALSE, min_n = 200L,
                        desc = "coverage of recurrence plots and simultaneous recurrences",
                        ref = "Hirata & Aihara (2010)", tuning = c("m", "tau", "rate", "theta")),
  cd             = list(label = "Correlation-dimension relations (Krakovska, 2019)",
                        short = "correlation dimension (Krakovska 2019)", stat = "D2.diff", family = "verdict", conditional = FALSE,
                        type = "nonlinear", pooled = FALSE, min_n = 500L,
                        desc = "relations of the correlation dimensions of the marginal and joint embeddings",
                        ref = "Krakovska (2019)", tuning = c("m", "tau")),
  pcmci          = list(label = "PCMCI / PCMCI+ via tigramite (Runge et al., 2019)",
                        short = "PCMCI (Runge et al. 2019)", stat = "MCI.stat", family = "external", conditional = TRUE,
                        type = "linear", pooled = FALSE, min_n = 50L, needs = "reticulate + Python tigramite",
                        desc = "PC-stable condition selection and momentary conditional independence tests",
                        ref = "Runge et al. (2019)", tuning = character(0))
)

.method_names <- names(.method_registry)

.method_info <- function(method) {
  info <- .method_registry[[method]]
  if (is.null(info)) list(label = method, short = method, stat = "stat", family = "surrogate", conditional = FALSE,
                          type = "nonlinear", pooled = FALSE, min_n = 30L, tuning = character(0)) else info
}

# ---- method settings ------------------------------------------------------------
# Tuning settings shared by several of the methods of version 0.0.1 (they were
# top-level arguments of rcausality(); they are still accepted there through
# `...`), and the settings of the methods added later.  All of them can be
# given in `method_args`; unknown names are an error, so that a misspelt
# setting is not silently ignored.
.tuning_defaults <- list(k = NULL, k_perm = 5L, history = "auto", memory = "auto", iid_injection = "allowed",
                         k_nn = NULL, n_bins = 2L, n_markov = 99L, max_constr = 4L, m = 3L, tau = 1L,
                         k_range = seq(10L, 40L, by = 5L), prep = "rank", K = 10L, n_surr = 50L, theta = NULL,
                         model = NULL, bins = 4L, k_target = 2L, k_source = 2L, rate = 0.05)

.method_arg_defaults <- list(
  toda_yamamoto  = list(d_max = 1L),
  diks_panchenko = list(eps = NULL, C = 8, filter = "none"),
  rf_granger     = list(num_trees = 200L, mtry = NULL, min_node_size = 5L, folds = 5L, shadow = TRUE,
                        importance = TRUE)
)

# the settings of one method with their defaults (for causal_methods() and the checks)
.method_settings <- function(method) {
  c(.tuning_defaults[.method_info(method)$tuning], .method_arg_defaults[[method]])
}

.method_args <- function(method, args) {
  def <- .method_arg_defaults[[method]]
  if (is.null(args) || !length(args)) return(def)
  utils::modifyList(def, args[intersect(names(args), names(def))], keep.null = TRUE)
}

# merge the legacy top-level settings (`dots`) and `method_args` into the full
# tuning list and the method's own settings
.resolve_settings <- function(method, method_args = list(), dots = list()) {
  if (length(dots)) {
    if (is.null(names(dots)) || any(!nzchar(names(dots)))) stop("Unnamed arguments in `...` are not allowed.", call. = FALSE)
    bad <- setdiff(names(dots), names(.tuning_defaults))
    if (length(bad)) stop("Unused argument(s): ", paste(bad, collapse = ", "),
                          ". Method settings go in `method_args` (see causal_methods(\"", method, "\")).", call. = FALSE)
  }
  if (is.null(method_args)) method_args <- list()
  if (length(method_args)) {
    if (!is.list(method_args) || is.null(names(method_args)) || any(!nzchar(names(method_args))))
      stop("`method_args` must be a named list.", call. = FALSE)
    ok <- names(.method_settings(method))
    bad <- setdiff(names(method_args), ok)
    if (length(bad))
      stop(sprintf("Unknown setting(s) for method \"%s\" in `method_args`: %s (available: %s).", method,
                   paste(bad, collapse = ", "), if (length(ok)) paste(ok, collapse = ", ") else "none"), call. = FALSE)
  }
  tn <- names(.tuning_defaults)
  tune <- utils::modifyList(.tuning_defaults, dots, keep.null = TRUE)
  tune <- utils::modifyList(tune, method_args[intersect(names(method_args), tn)], keep.null = TRUE)
  list(tune = tune[tn], own = method_args[setdiff(names(method_args), tn)])
}

#' The causality methods of the package and their settings
#'
#' Without an argument, a table of the methods that can be chosen with the
#' `method` argument of [rcausality()] (and of the functions built on it):
#' what the method tests, whether it is linear, whether it can condition on
#' the other series (`conditional = TRUE`), how its Bayesian variant
#' (`bayes = TRUE`) is obtained, and which panel-data tests are available
#' (every method has unit-level and group-level panel tests; `pooled` marks
#' the methods with a pooled estimator on the stacked units). With a method
#' name, the settings of that method (to be given in `method_args`) and their
#' defaults.
#'
#' @param method optional method name.
#' @return A data frame of class `causal_methods` (one row per method), or,
#'   for one method, a list with its description, capabilities, `settings`
#'   (named list of defaults) and reference.
#' @seealso [rcausality()], [panel_causality()], [causal_drivers()]
#' @export
#' @examples
#' causal_methods()
#' causal_methods("te")          # settings of the transfer entropy and their defaults
causal_methods <- function(method = NULL) {
  bayes_how <- c(regression = "Bayes factor of nested regressions", surrogate = "surrogate Bayes factor",
                 posterior = "posterior of the causal cases", verdict = "surrogate Bayes factor",
                 external = "Savage-Dickey Bayes factor")
  if (!is.null(method)) {
    method <- match.arg(method, .method_names); info <- .method_info(method)
    return(structure(list(method = method, description = info$desc, label = info$label, type = info$type,
                          conditional = info$conditional, bayes = unname(bayes_how[info$family]),
                          panel = c("units", "group", if (isTRUE(info$pooled)) "pooled"), min_length = info$min_n,
                          needs = if (is.null(info$needs)) NA_character_ else info$needs,
                          settings = .method_settings(method), reference = info$ref), class = "causal_method"))
  }
  tab <- do.call(rbind, lapply(.method_names, function(m) {
    i <- .method_info(m)
    data.frame(method = m, description = i$desc, type = i$type, conditional = i$conditional,
               bayes = unname(bayes_how[i$family]), pooled = isTRUE(i$pooled), min_length = i$min_n,
               needs = if (is.null(i$needs)) "" else i$needs, reference = i$ref, stringsAsFactors = FALSE)
  }))
  structure(tab, class = c("causal_methods", "data.frame"))
}

#' @export
print.causal_methods <- function(x, ...) {
  if (!all(c("method", "type", "conditional", "pooled", "needs", "description") %in% names(x)))     # a subset of the columns
    return(invisible(print.data.frame(x, row.names = FALSE, right = FALSE)))
  cat("Causality methods of Rcausality: rcausality(x, method = <name>, bayes = TRUE / FALSE, cores = <n>)\n")
  tb <- data.frame(method = x$method, type = x$type, conditional = ifelse(x$conditional, "yes", "-"),
                   pooled.panel = ifelse(x$pooled, "yes", "-"), needs = x$needs, description = x$description,
                   stringsAsFactors = FALSE)
  print(tb, row.names = FALSE, right = FALSE)
  cat("  every method has a Bayesian variant (bayes = TRUE) and unit- and group-level panel tests (id = , time = );\n",
      " pooled.panel: a pooled estimator on the stacked units exists.  causal_methods(\"<name>\") lists the settings.\n", sep = "")
  invisible(x)
}

#' @export
print.causal_method <- function(x, ...) {
  cat(sprintf("%s -- %s\n", x$method, x$label))
  cat(sprintf("  %s; %s; %s\n", x$description, x$type, if (x$conditional) "pairwise or conditional" else "pairwise only"))
  cat(sprintf("  Bayesian variant: %s; panel tests: %s; advisable series length: at least %d\n", x$bayes,
              paste(x$panel, collapse = ", "), x$min_length))
  if (!is.na(x$needs)) cat("  needs:", x$needs, "\n")
  if (length(x$settings)) {
    cat("  settings (method_args = list(...)) and defaults:\n")
    for (nm in names(x$settings)) {
      v <- x$settings[[nm]]
      cat(sprintf("    %-14s %s\n", nm, if (is.null(v)) "NULL (automatic)" else paste(format(v), collapse = ", ")))
    }
  } else cat("  no method-specific settings\n")
  cat("  reference:", x$reference, "\n")
  invisible(x)
}

# ---- the `pairs` table --------------------------------------------------------

.pairs_cols <- c("cause", "effect", "lag", "statistic", "df1", "df2", "p_value", "p_adj", "bf10", "post_prob",
                 "evidence", "effect_size", "r2_gain", "sign", "n_eff", "significant")

# one directed test as a one-row data frame (missing fields are NA)
.pair_row <- function(cause, effect, lag = NA_integer_, statistic = NA_real_, df1 = NA_real_, df2 = NA_real_,
                      p_value = NA_real_, bf10 = NA_real_, effect_size = NA_real_, r2_gain = NA_real_,
                      sign = NA_real_, n_eff = NA_integer_) {
  data.frame(cause = cause, effect = effect, lag = as.integer(lag), statistic = as.numeric(statistic),
             df1 = as.numeric(df1), df2 = as.numeric(df2), p_value = as.numeric(p_value), bf10 = as.numeric(bf10),
             effect_size = as.numeric(effect_size), r2_gain = as.numeric(r2_gain), sign = as.numeric(sign),
             n_eff = as.integer(n_eff), stringsAsFactors = FALSE)
}

# multiplicity correction, posterior probability, evidence label and decision
# of a pairs table; rule = "p" (adjusted p < alpha) or "bf" (BF10 > threshold)
.finish_pairs <- function(pairs, adjust, alpha, prior_prob = 0.5, bf_threshold = 3, rule = c("p", "bf")) {
  rule <- match.arg(rule)
  pairs$p_adj <- if (adjust == "none") pairs$p_value else stats::p.adjust(pairs$p_value, method = adjust)
  pairs$post_prob <- ifelse(is.na(pairs$bf10), NA_real_, pairs$bf10 * prior_prob / (pairs$bf10 * prior_prob + 1 - prior_prob))
  pairs$evidence <- as.character(bf_evidence(pairs$bf10))
  pairs$significant <- if (rule == "p") !is.na(pairs$p_adj) & pairs$p_adj < alpha else !is.na(pairs$bf10) & pairs$bf10 > bf_threshold
  pairs <- pairs[, .pairs_cols]
  rownames(pairs) <- NULL
  pairs
}

# Driver of the directed pairwise methods added in 0.0.2: `pair_fun(effect,
# cause, controls, cores)` returns list(row = <named list of the .pair_row()
# fields>, by_lag = <data frame (lag, statistic, p) or NULL>, null = <list
# (statistic, null, alternative) or NULL>, extra = <one-row data frame or
# NULL>).  Pairs run in parallel (or the work inside each pair, see
# .split_cores()); the extra rows are collected in the result slot `extra`.
.rcausality_directed <- function(all, grid, causes, effects, method, pair_fun, alpha, adjust, prior_prob,
                                 bf_threshold, conditional, cl, cores, rule = "p", extra = NULL) {
  sc <- .split_cores(nrow(grid), cores)
  res <- .plapply(seq_len(nrow(grid)), function(i) {
    cs <- grid$cause[i]; ef <- grid$effect[i]
    ctrl_names <- if (conditional) setdiff(names(all), c(cs, ef)) else character(0)
    ctrl <- if (length(ctrl_names)) as.matrix(all[, ctrl_names, drop = FALSE]) else NULL
    pair_fun(all[[ef]], all[[cs]], ctrl, cores = sc$inner)
  }, cores = sc$outer)
  pairs <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i)
    do.call(.pair_row, c(list(cause = grid$cause[i], effect = grid$effect[i]), res[[i]]$row))))
  pairs <- .finish_pairs(pairs, adjust, alpha, prior_prob, bf_threshold, rule)
  bl <- lapply(seq_len(nrow(grid)), function(i) {
    b <- res[[i]]$by_lag
    if (is.null(b) || !nrow(b)) NULL else data.frame(cause = grid$cause[i], effect = grid$effect[i], b, stringsAsFactors = FALSE)
  })
  bl <- do.call(rbind, bl); if (!is.null(bl)) rownames(bl) <- NULL
  out <- .new_causres(pairs, series = names(all), causes = causes, effects = effects, method = method,
                      conditional = conditional, alpha = alpha, bf_threshold = bf_threshold, adjust = adjust,
                      n = nrow(all), call = cl, data = all, by_lag = bl)
  nulls <- lapply(res, `[[`, "null")
  if (!all(vapply(nulls, is.null, TRUE))) out$null <- stats::setNames(nulls, paste(pairs$cause, "->", pairs$effect))
  ex <- lapply(seq_len(nrow(grid)), function(i) {
    e <- res[[i]]$extra
    if (is.null(e)) NULL else data.frame(cause = grid$cause[i], effect = grid$effect[i], e, stringsAsFactors = FALSE,
                                         check.names = FALSE)
  })
  ex <- do.call(rbind, ex)
  if (!is.null(ex)) { rownames(ex) <- NULL; out$extra <- ex }
  out
}
