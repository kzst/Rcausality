# Random-forest Granger-like causality ----------------------------------------
# Granger's predictive definition with a nonparametric regression: does a
# random forest predict the effect better when the cause's past is added to
# the effect's own past (and the controls)?  The evidence combines
#  * the out-of-bag gain  1 - MSE_oob(own + cause) / MSE_oob(own),
#    tested against conditional shifts of the cause
#    (the lags of the cause are split into their regression on the
#    conditioning features and residuals; the residuals are shifted
#    circularly: autocorrelation and dependence on the conditioning
#    features kept, the coupling destroyed);
#  * the honest rolling-origin (expanding-window) out-of-sample gain;
#  * the permutation importance of the cause's lags and a shadow-feature
#    (Boruta-type) rule: a lag of the cause must be more important than every
#    lag of a circularly shifted copy of the cause.
# This is the approach of Kosztyan, Matrai and Kiss (2026) for gravity-model
# panels, implemented here for time series (and for panels in panel_causality()).

.need_ranger <- function(what = "method = \"rf_granger\"") {
  if (!requireNamespace("ranger", quietly = TRUE))
    stop(what, " needs the package 'ranger' (install.packages(\"ranger\")).", call. = FALSE)
}

# a seed for ranger drawn from the current random stream (ranger has its own
# generator; seeding it from R's keeps the results reproducible and
# independent of the number of workers)
.draw_seed <- function() sample.int(.Machine$integer.max, 1L)

.lag_block <- function(v, L, prefix) { m <- .lagmat(v, L); colnames(m) <- paste0(prefix, ".l", seq_len(L)); m }

.circ_shift <- function(v, min_shift) {
  n <- length(v); s <- sample(min_shift:(n - min_shift), 1L)
  c(v[(s + 1):n], v[seq_len(s)])
}

.rf_fit <- function(X, y, settings, importance = "none") {
  # mtry = "all": every split considers all features (bagged trees).  With a random subset of a handful of features,
  # an additional feature improves a forest even when it is redundant or irrelevant, because it adds diversity to the
  # trees; bagged trees do not have this artefact, so that the gain of a cause without novel information is zero.
  mtry <- if (identical(settings$mtry, "all")) ncol(X) else settings$mtry
  ranger::ranger(x = as.data.frame(X), y = y, num.trees = as.integer(settings$num_trees), mtry = mtry,
                 min.node.size = .rf_node_size(X, y, settings), importance = importance, num.threads = 1L,
                 seed = .draw_seed(), verbose = FALSE)
}

# Minimal node size of the trees.  With few features -- a handful of lags --
# fully grown trees (ranger's default node size of 5) interpolate the noise:
# the forest on the own past of an AR(1) series then has an out-of-bag R^2
# near zero, and *any* added feature lowers its error, because the random
# choice of the split variable acts as a regulariser.  "auto" therefore
# selects the node size with the smallest out-of-bag error from a grid that
# grows with the sample (5 and 2%, 5%, 10%, 20% of the observations).
.rf_node_grid <- function(n) sort(unique(pmax(5L, pmin(as.integer(round(c(5, n * c(0.02, 0.05, 0.1, 0.2)))), max(5L, n %/% 3L)))))
.rf_node_size <- function(X, y, settings) {
  ms <- settings$min_node_size
  if (is.numeric(ms)) return(as.integer(ms))
  grid <- .rf_node_grid(length(y))
  if (length(grid) == 1L) return(grid)
  err <- vapply(grid, function(g) { s <- settings; s$min_node_size <- g; .rf_fit(X, y, s)$prediction.error }, 0)
  grid[which.min(err)]
}

# rolling-origin (expanding window) out-of-sample MSE of the own-past model
# and of the model with the cause's past: folds + 1 consecutive blocks, block
# f + 1 predicted from blocks 1..f.  With min_node_size = "auto" every forest
# selects its node size on its own training part: a size chosen on the whole
# sample is too coarse for the short early windows.
.rf_rolling <- function(R, U, y, settings, cores = 1L) {
  N <- length(y); folds <- as.integer(settings$folds)
  if (folds < 1L || N < 10L * (folds + 1L)) return(c(mse_own = NA_real_, mse_full = NA_real_))
  edges <- floor(seq(0, N, length.out = folds + 2L))
  err <- .plapply(seq_len(folds), function(f) {
    tr <- seq_len(edges[f + 1L]); te <- (edges[f + 1L] + 1L):edges[f + 2L]
    fr <- .rf_fit(R[tr, , drop = FALSE], y[tr], settings); fu <- .rf_fit(U[tr, , drop = FALSE], y[tr], settings)
    pr <- stats::predict(fr, data = as.data.frame(R[te, , drop = FALSE]), num.threads = 1L)$predictions
    pu <- stats::predict(fu, data = as.data.frame(U[te, , drop = FALSE]), num.threads = 1L)$predictions
    c(sum((y[te] - pr)^2), sum((y[te] - pu)^2), length(te))
  }, cores = cores)
  E <- do.call(rbind, err)
  c(mse_own = sum(E[, 1]) / sum(E[, 3]), mse_full = sum(E[, 2]) / sum(E[, 3]))
}

# The lags of the cause split into location, scale and standardized residuals given the conditioning features (see
# .cond_split() for the model); a surrogate block has the residuals shifted by s rows.  The location is the least-
# squares fit on the polynomial basis with the smallest leave-one-out error, plus the out-of-bag predictions of a
# forest fitted to its residuals (one per lag, node size selected), which adds what a polynomial misses and next to
# nothing otherwise; the scale is the fit of the absolute residuals of .scale_fit(): the basis of .scale_bases()
# (a constant, the absolute values of the features, the polynomial bases) with the smallest leave-one-out error,
# refitted with the weights 1 / s^2.
# (A forest alone approximates a smooth dependence by a step function, and the surrogate cause is then a worse copy
# of a correlated control than the cause itself: in simulations such surrogates rejected a true conditional null
# hypothesis far too often.)
.rf_split_cause <- function(R, C, settings) {
  C <- as.matrix(C); W <- .std_mat(R); bases <- .poly_bases(W); sbases <- .scale_bases(W, bases); fit <- C; sc <- C
  for (j in seq_len(ncol(C))) {
    f <- .press_fit(bases, C[, j], from = min(2L, length(bases)))
    f <- f + .rf_fit(R, C[, j] - f, settings)$predictions
    a <- abs(C[, j] - f); lo <- .scale_floor * mean(a)
    fit[, j] <- f; sc[, j] <- pmax(.scale_fit(sbases, a, lo), lo, 1e-12)
  }
  list(fit = fit, scale = sc, res = (C - fit) / sc, r2 = .cond_r2(C, fit))
}
.rf_surr_cause <- function(sp, s, rows = NULL) {
  if (is.null(rows)) { n <- nrow(sp$res); rows <- ((seq_len(n) + s - 1L) %% n) + 1L }
  out <- .cond_surr(sp, rows); dimnames(out) <- dimnames(sp$fit); out
}

.rf_pair <- function(effect, cause, controls = NULL, lag = 2L, n_surr = 99L, settings = .method_arg_defaults$rf_granger,
                     cores = 1L) {
  .need_ranger()
  L <- max(1L, as.integer(lag)); n <- length(effect)
  y <- effect[(L + 1):n]
  R <- .lag_block(effect, L, "own")
  if (!is.null(controls)) for (j in seq_len(ncol(controls)))
    R <- cbind(R, .lag_block(controls[, j], L, if (is.null(colnames(controls))) paste0("control", j) else colnames(controls)[j]))
  C <- .lag_block(cause, L, "cause")
  U <- cbind(R, C)
  min_shift <- max(L + 1L, floor(0.1 * n))
  # the node size is chosen on the own-past forest, i.e. without reference to the cause, and used for every forest of
  # the pair: the observed statistic stays exchangeable with its surrogates and the own-past forest is not handicapped
  settings_roll <- settings
  settings$min_node_size <- .rf_node_size(R, y, settings)
  f_own <- .rf_fit(R, y, settings)
  f_full <- .rf_fit(U, y, settings, if (isTRUE(settings$importance)) "permutation" else "none")
  mse_own <- f_own$prediction.error
  gain <- 1 - f_full$prediction.error / mse_own
  # Conditional shifts of the cause.  The null hypothesis is conditional: the past of the cause
  # does not help GIVEN the own past of the effect and the controls, on which it may depend -- on the controls when
  # they are correlated with the cause, on the own past of the effect in the reverse direction of a link.  A circular
  # shift of the cause would make it independent of all of them and test another hypothesis.  Every lag of the cause
  # is therefore split into its regression on the conditioning features (out-of-bag forest predictions) and a
  # residual; a surrogate keeps the first part and shifts the residuals circularly, all lags by the same offset, so
  # that their serial dependence is kept.  For a cause that is unrelated to the conditioning features this is the
  # circular shift of the cause itself.
  sp <- .rf_split_cause(R, C, settings_roll)
  n_row <- nrow(C); lo_s <- max(L + 1L, floor(0.1 * n_row))
  shifts <- .shift_draws(lo_s, n_row - lo_s, n_surr)
  # (a determined cause has no surrogates that differ from it: no test, see .determined(); the surrogate forests
  # are not grown and the replicates stand at the observed value)
  null <- if (.determined(sp$r2)) rep(gain, length(shifts)) else unlist(.plapply(shifts, function(s) {
    1 - .rf_fit(cbind(R, .rf_surr_cause(sp, s)), y, settings)$prediction.error / mse_own
  }, cores = cores))
  # Granger's definition asks for an improvement: a non-positive gain is no
  # evidence (in coupled deterministic systems the driven series carries
  # information about its driver, and the true past of the effect can harm
  # the prediction of the cause less than a surrogate does)
  p <- if (gain > 0 && !.determined(sp$r2)) (1 + sum(null >= gain)) / (1 + length(null)) else 1
  roll <- .rf_rolling(R, U, y, settings_roll, cores)
  vi <- if (isTRUE(settings$importance)) f_full$variable.importance else NULL
  vic <- if (is.null(vi)) NA_real_ else max(vi[colnames(C)])
  vlag <- if (is.null(vi)) NA_integer_ else which.max(vi[colnames(C)])
  shadow_max <- NA_real_; hit <- NA
  if (isTRUE(settings$shadow)) {
    S <- .rf_surr_cause(sp, .shift_draws(lo_s, n_row - lo_s, 1L)); colnames(S) <- sub("^cause", "shadow", colnames(C))
    fs <- .rf_fit(cbind(U, S), y, settings, "permutation")
    vs <- fs$variable.importance
    shadow_max <- max(vs[colnames(S)])
    hit <- max(vs[colnames(C)]) > 0 && max(vs[colnames(C)]) > shadow_max
  }
  oos <- if (is.na(roll[["mse_own"]]) || roll[["mse_own"]] <= 0) NA_real_ else 1 - roll[["mse_full"]] / roll[["mse_own"]]
  list(gain = gain, p_value = p, null = null, oos_gain = oos, mse_own = mse_own, mse_full = f_full$prediction.error,
       importance = vi, vi_cause = vic, vi_lag = unname(vlag), shadow_max = shadow_max, shadow_hit = hit, n_eff = length(y),
       lag = L, node_size = settings$min_node_size, cond_r2 = sp$r2)
}

.rf_pair_fun <- function(fixed_lag, max_lag, lag_criterion, n_perm, settings) {
  function(effect, cause, ctrl, cores = 1L) {
    inv <- cbind(effect, cause, ctrl)
    p <- if (is.null(fixed_lag)) select_lag(inv, max_lag, lag_criterion) else fixed_lag
    p <- max(1L, min(p, .cap_lag(nrow(inv), ncol(inv))))
    g <- .rf_pair(effect, cause, ctrl, p, n_perm, settings, cores)
    list(row = list(lag = p, statistic = g$gain, p_value = g$p_value,
                    effect_size = if (is.na(g$oos_gain)) g$gain else g$oos_gain, r2_gain = g$gain, n_eff = g$n_eff),
         null = list(statistic = g$gain, null = g$null, alternative = "greater", cond_r2 = g$cond_r2),
         extra = data.frame(oob_gain = g$gain, oos_gain = g$oos_gain, mse_own = g$mse_own, mse_full = g$mse_full,
                            vi_cause = g$vi_cause, vi_lag = g$vi_lag, shadow_max = g$shadow_max, shadow_hit = g$shadow_hit,
                            node_size = g$node_size))
  }
}

#' Random-forest Granger-like causality test
#'
#' A nonlinear, model-free reading of Granger's predictive definition: a
#' random forest (package \pkg{ranger}) predicts `y[t]` from the own past
#' `y[t-1..t-lag]` with and without the past of `x`. The out-of-bag gain
#' \deqn{G = 1 - \mathrm{MSE}_{oob}(\text{own} + x) / \mathrm{MSE}_{oob}(\text{own})}{G = 1 - MSE_oob(own + x) / MSE_oob(own)}
#' is tested against `n_surr` conditional shifts of `x`, surrogates that
#' keep what the null hypothesis allows: the lags of `x` are split into a location given
#' the own past of `y` (the polynomial with the smallest leave-one-out error,
#' plus the out-of-bag predictions of a forest for its residuals), a scale
#' (a weighted least-squares fit of the absolute residuals on a polynomial
#' of the own past or on its absolute values) and standardized residuals,
#' and the residuals are shifted circularly. The surrogates keep
#' the autocorrelation of both series and the dependence of `x` on the past
#' of `y` (which exists in the reverse direction of a link) and destroy the
#' coupling, so that the optimism of out-of-bag errors on dependent data
#' cancels; for a cause that is unrelated to the past of `y` they are
#' circular shifts of `x`. A non-positive
#' gain gives `p = 1`, since Granger's definition asks for an improvement of
#' the prediction (not merely a dependence); so does a cause that the past
#' of `y` determines (99% of the variance of a lag or more, as in noise-free
#' deterministic systems), for which no surrogate differs from the cause and
#' no test is carried out. Reported with
#' it: the rolling-origin (expanding-window) out-of-sample gain over `folds`
#' forecast blocks, the permutation importance of the lags of `x`, and the
#' shadow-feature rule of Kosztyan, Matrai and Kiss (2026): a link is
#' supported when some lag of `x` is more important than every lag of a
#' circularly shifted copy of `x` added as a shadow.
#'
#' @param x,y numeric vectors (both directions are tested).
#' @param lag number of past values of each series.
#' @param n_surr number of surrogates (conditional circular shifts, see
#'   above; at most the number of admissible shifts, about 0.8 times the
#'   length of the series).
#' @param num_trees,mtry forest settings (passed to [ranger::ranger()]).
#'   `mtry = "all"` (default): every split considers all features (bagged
#'   trees); `NULL` is the default of ranger, a random subset of about the
#'   square root of the number of features. With a random subset of a handful
#'   of lags, a redundant or even an irrelevant feature improves a forest
#'   because it adds diversity to the trees, which would count as a gain of
#'   the cause; bagged trees do not have this artefact.
#' @param min_node_size minimal node size of the trees: a number, or
#'   `"auto"` (default): the size with the smallest out-of-bag error of the
#'   forest on the own past of the effect, from the grid 5 and 2%, 5%, 10%,
#'   20% of the observations; it is then used for all forests of the
#'   out-of-bag statistic and its surrogates, while every forest of the
#'   rolling-origin evaluation selects its size on its own training part.
#'   With only a few lags as features, fully grown trees (node size 5)
#'   interpolate the noise, the own-past forest predicts poorly, and any added
#'   series appears to improve it.
#' @param folds number of rolling-origin forecast blocks (`0` = none).
#' @param shadow logical; compute the shadow-feature rule.
#' @param seed optional RNG seed (the forests are seeded from R's generator).
#' @param cores CPU cores: the surrogate forests and the forecast blocks run in
#'   parallel (see [Rcausality-parallel]); the results do not depend on it.
#' @return An object of class `rf_granger`: `tests` (one row per direction:
#'   out-of-bag gain, surrogate p-value, out-of-sample gain, importance of the
#'   best lag of the cause, shadow decision), `importance` (permutation
#'   importance of every feature, per direction), `null` (surrogate gains).
#'   Method: `print`.
#' @references
#' Kosztyan, Z. T., Matrai, R. and Kiss, D. (2026). Is the world getting smaller? Network-driven heterogeneity in gravity-model elasticities. *The World Economy*, 49(11), 2478-2525. \doi{10.1111/twec.70148}
#'
#' Breiman, L. (2001). Random forests. *Machine Learning*, 45(1), 5-32. \doi{10.1023/A:1010933404324}
#'
#' Wright, M. N. and Ziegler, A. (2017). ranger: a fast implementation of random forests for high dimensional data in C++ and R. *Journal of Statistical Software*, 77(1), 1-17. \doi{10.18637/jss.v077.i01}
#'
#' Kursa, M. B. and Rudnicki, W. R. (2010). Feature selection with the Boruta package. *Journal of Statistical Software*, 36(11), 1-13. \doi{10.18637/jss.v036.i11}
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("ranger", quietly = TRUE)) {
#'   # a nonlinear link x -> y without a linear trace: the forest gains from
#'   # the past of x (small surrogate p-value, positive out-of-sample gain),
#'   # not the other way round
#'   d <- sim_nonlinear(300, "x_to_y", seed = 1)
#'   rf_granger(d$x, d$y, lag = 1, n_surr = 49, num_trees = 100, seed = 1)
#' }
#' }
rf_granger <- function(x, y, lag = 2L, n_surr = 99L, num_trees = 200L, mtry = "all", min_node_size = "auto", folds = 5L,
                       shadow = TRUE, seed = NULL, cores = getOption("Rcausality.cores", 1L)) {
  .need_ranger("rf_granger()")
  .local_seed(seed)
  x <- .prep_na(x); y <- .prep_na(y)
  if (length(x) != length(y)) stop("x and y must have the same length.", call. = FALSE)
  if (!identical(min_node_size, "auto")) .check_number(min_node_size, "min_node_size", lower = 1, integer = TRUE)
  st <- list(num_trees = num_trees, mtry = mtry, min_node_size = min_node_size, folds = folds, shadow = shadow,
             importance = TRUE)
  res <- .plapply(list(c("x", "y"), c("y", "x")), function(d) {
    cs <- if (d[1] == "x") x else y; ef <- if (d[1] == "x") y else x
    .rf_pair(ef, cs, NULL, lag, n_surr, st, cores = 1L)
  }, cores = cores)
  tests <- data.frame(direction = c("x -> y", "y -> x"), oob_gain = vapply(res, `[[`, 0, "gain"),
                      p_value = vapply(res, `[[`, 0, "p_value"), oos_gain = vapply(res, `[[`, 0, "oos_gain"),
                      vi_cause = vapply(res, `[[`, 0, "vi_cause"), vi_lag = vapply(res, function(r) as.integer(r$vi_lag), 0L),
                      shadow_hit = vapply(res, function(r) as.logical(r$shadow_hit), NA),
                      node_size = vapply(res, function(r) as.integer(r$node_size), 0L), stringsAsFactors = FALSE)
  structure(list(tests = tests, importance = stats::setNames(lapply(res, `[[`, "importance"), tests$direction),
                 null = stats::setNames(lapply(res, `[[`, "null"), tests$direction), lag = as.integer(lag),
                 n_surr = n_surr, settings = st), class = "rf_granger")
}

#' @export
print.rf_granger <- function(x, digits = 3, ...) {
  .out(sprintf("Random-forest Granger-like causality (lag %d, %d trees, node size %s, %d conditional-shift surrogates)\n", x$lag,
              as.integer(x$settings$num_trees), paste(unique(x$tests$node_size), collapse = "/"), x$n_surr))
  t <- x$tests
  print(data.frame(direction = t$direction, OOB.gain = sprintf("%+.3f", t$oob_gain), p = .fmt_p(t$p_value),
                   OOS.gain = ifelse(is.na(t$oos_gain), "-", sprintf("%+.3f", t$oos_gain)),
                   VI.best = ifelse(is.na(t$vi_cause), "-", sprintf("%.3g (lag %d)", t$vi_cause, t$vi_lag)),
                   shadow = ifelse(is.na(t$shadow_hit), "-", ifelse(t$shadow_hit, "beats shadow", "not above shadow"))),
        row.names = FALSE, right = FALSE)
  .out("  OOB.gain = 1 - MSE(own + cause) / MSE(own) out of bag, p = against conditional shifts of the cause (location and scale given the own past kept, residuals shifted circularly);\n",
      " OOS.gain = the same ratio from rolling-origin forecasts; VI = permutation importance of the cause's best lag\n", sep = "")
  invisible(x)
}
