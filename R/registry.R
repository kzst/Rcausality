# Method registry and the shared assembly of the `pairs` table ---------------
#
# One entry per method of rcausality(): the labels used by print() and plot(),
# the name of the statistic column, whether the method can condition on the
# other series, and the family that decides how its Bayesian variant works
# ("regression": closed-form or sampled Bayes factors of nested regressions;
# "surrogate": the statistic against its surrogate / asymptotic null through
# surrogate_bf(); "posterior": posterior over the causal cases; "verdict":
# one verdict per pair, no directed p-values; "external": PCMCI).

.method_registry <- list(
  granger        = list(label = "Granger causality (F-test)", short = "Granger causality (F-test)",
                        stat = "F", family = "regression", conditional = TRUE),
  bayes_granger  = list(label = "Bayesian Granger causality", short = "Bayesian Granger causality (Bayes factor)",
                        stat = "F", family = "regression", conditional = TRUE),
  toda_yamamoto  = list(label = "Toda-Yamamoto lag-augmented Granger causality (integrated or cointegrated series)",
                        short = "Toda-Yamamoto Granger causality", stat = "F", family = "regression", conditional = TRUE),
  diks_panchenko = list(label = "Diks-Panchenko nonparametric (nonlinear) Granger causality test",
                        short = "Diks-Panchenko nonparametric test", stat = "T.n", family = "surrogate", conditional = TRUE),
  rf_granger     = list(label = "Random-forest Granger-like causality (out-of-bag gain, circular-shift surrogates)",
                        short = "random-forest Granger-like causality", stat = "OOB.gain", family = "surrogate",
                        conditional = TRUE),
  cmi            = list(label = "Conditional mutual information (kNN, local permutation)",
                        short = "kNN conditional mutual information", stat = "CMI", family = "surrogate", conditional = TRUE),
  cmi_lag        = list(label = "Lag-shifted conditional mutual information (Zlatniczki et al., 2021)",
                        short = "lag-shifted CMI (Zlatniczki et al. 2021)", stat = "CMI", family = "surrogate", conditional = TRUE),
  markov         = list(label = "Markov-chain approach (Stippinger et al., 2023)",
                        short = "Markov-chain approach (Stippinger et al. 2023)", stat = "CMI", family = "surrogate",
                        conditional = FALSE),
  df             = list(label = "Degrees-of-freedom / conditional-variance method (Telcs, Kurbucz & Jakovac, 2025)",
                        short = "degrees of freedom (Telcs et al. 2025)", stat = "sd.red", family = "surrogate", conditional = TRUE),
  dc             = list(label = "Bayesian dimensional causality (Benko et al., 2024)",
                        short = "dimensional causality (Benko et al. 2024)", stat = "P.direct", family = "posterior",
                        conditional = FALSE),
  ccdh           = list(label = "CCDH: complexity-based causal discovery (Stippinger et al., 2023)",
                        short = "CCDH (Stippinger et al. 2023)", stat = "recon", family = "surrogate", conditional = FALSE),
  dfe            = list(label = "df-entropy method with random-forest classification (Jakovac et al.)",
                        short = "df-entropy + random forest (Jakovac et al.)", stat = "P.direct", family = "posterior",
                        conditional = FALSE),
  te             = list(label = "Transfer entropy, binned (Schreiber, 2000) with surrogate test",
                        short = "transfer entropy (Schreiber 2000)", stat = "TE.bits", family = "surrogate", conditional = FALSE),
  ccm            = list(label = "Convergent cross mapping (Sugihara et al., 2012)",
                        short = "convergent cross mapping (Sugihara et al. 2012)", stat = "rho", family = "surrogate",
                        conditional = FALSE),
  pid            = list(label = "MAR-PID: unique + synergistic information of the cause's past about the effect's future (Telcs, 2026)",
                        short = "MAR-PID (Telcs 2026)", stat = "U+S.bits", family = "surrogate", conditional = FALSE),
  hirata         = list(label = "Recurrence-plot criteria: coverage and simultaneous recurrences (Hirata & Aihara, 2010)",
                        short = "recurrence plots (Hirata & Aihara 2010)", stat = "cover.idx", family = "verdict",
                        conditional = FALSE),
  cd             = list(label = "Correlation-dimension relations (Krakovska, 2019)",
                        short = "correlation dimension (Krakovska 2019)", stat = "D2.diff", family = "verdict", conditional = FALSE),
  pcmci          = list(label = "PCMCI / PCMCI+ via tigramite (Runge et al., 2019)",
                        short = "PCMCI (Runge et al. 2019)", stat = "MCI.stat", family = "external", conditional = TRUE)
)

.method_names <- names(.method_registry)

.method_info <- function(method) {
  info <- .method_registry[[method]]
  if (is.null(info)) list(label = method, short = method, stat = "stat", family = "surrogate", conditional = FALSE) else info
}

# settings of the methods added in version 0.0.2 (argument `method_args` of
# rcausality()): defaults, merged with the user's list; unknown names are an
# error, so that a misspelt setting is not silently ignored
.method_arg_defaults <- list(
  toda_yamamoto  = list(d_max = 1L),
  diks_panchenko = list(eps = NULL, C = 8, filter = "none"),
  rf_granger     = list(num_trees = 200L, mtry = NULL, min_node_size = 5L, folds = 5L, shadow = TRUE,
                        importance = TRUE)
)

.method_args <- function(method, args) {
  def <- .method_arg_defaults[[method]]
  if (is.null(args) || !length(args)) return(def)
  if (!is.list(args) || is.null(names(args)) || any(!nzchar(names(args))))
    stop("`method_args` must be a named list.", call. = FALSE)
  bad <- setdiff(names(args), names(def))
  if (length(bad))
    stop(sprintf("Unknown setting(s) for method \"%s\" in `method_args`: %s (available: %s).", method,
                 paste(bad, collapse = ", "), if (length(def)) paste(names(def), collapse = ", ") else "none"), call. = FALSE)
  utils::modifyList(def, args, keep.null = TRUE)
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
