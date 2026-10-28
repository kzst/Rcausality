# Causality-informed forecasting ---------------------------------------------
# Uses the significant links of a causres as a restricted VAR: every series is
# regressed on its own lags and on the lags of its significant causes; the
# gain over the own-lags-only model is measured on a hold-out part of the
# sample, and h-step-ahead forecasts are produced recursively.

.lag_design <- function(D, p, target, preds) {
  n <- nrow(D)
  X <- do.call(cbind, lapply(preds, function(v) .lagmat(D[[v]], p)))
  colnames(X) <- as.vector(outer(paste0("lag", seq_len(p)), preds, function(a, b) paste(b, a, sep = ".")))
  list(y = D[[target]][(p + 1):n], X = X)
}

#' Causality-informed forecasts from a `causres` object
#'
#' Fits, for every series, a linear model on its own `p` lags (baseline) and
#' a model that adds the `p` lags of the series' significant causes (the
#' restricted VAR implied by the detected links), evaluates both by rolling
#' one-step forecasts on the last `holdout` fraction of the sample, and
#' produces `h`-step-ahead forecasts (recursive, all series jointly) with
#' normal prediction intervals from the residual standard deviation.
#'
#' @param object a `causres` object (any method; the links in
#'   `object$pairs$significant` define the causes).
#' @param h forecast horizon.
#' @param p lag order of the models (default: the largest selected lag of
#'   the significant links, at least 1).
#' @param holdout fraction of the sample used for the out-of-sample
#'   evaluation (`0` = no evaluation).
#' @param level level of the prediction intervals.
#' @param refit re-run the causality analysis (same method and settings) on
#'   the training part alone and use *those* links in the hold-out
#'   evaluation (default), so that the evaluation does not see the hold-out
#'   data when the links are chosen; `FALSE` evaluates the links of
#'   `object`, which were selected on the full sample and make the
#'   comparison optimistic. The final forecasts always use the links of
#'   `object`.
#' @param ... unused.
#' @return An object of class `causfc`: `forecasts` (series, horizon, point,
#'   lower, upper, se, for the causal model; the intervals follow the
#'   multi-step MSE recursion of the restricted VAR, Luetkepohl 2005, Sect.
#'   2.2.2), `evaluation` (per series: RMSE and MAE of the baseline and of
#'   the causal model on the hold-out, the relative RMSE improvement, the
#'   causes used), `holdout_fit` (observed and one-step predicted values of
#'   both models on the hold-out), `fitted_causes`, `coefficients` (the
#'   restricted VAR coefficient array, series x series x lag), `sigma`
#'   (residual covariance), `p`, `data`. Methods: `print`, `summary`, `plot`.
#' @references
#' Granger, C. W. J. (1969). Investigating causal relations by econometric models and cross-spectral methods. *Econometrica*, 37(3), 424-438. \doi{10.2307/1912791}
#'
#' Luetkepohl, H. (2005). *New Introduction to Multiple Time Series Analysis*. Springer, Berlin. \doi{10.1007/978-3-540-27752-1}
#' @export
#' @examples
#' # Granger analysis of a VAR where x drives y, then a 5-step forecast from
#' # the restricted VAR that keeps only the detected link (x -> y) and the own
#' # lags: y is predicted from its past and from x, x from its past alone.
#' d <- sim_var(400, "x_to_y", seed = 1)
#' r <- rcausality(d)
#' fc <- predict(r, h = 5)
#' fc              # point forecasts with interval, and the hold-out evaluation
#' \donttest{
#' plot(fc)        # the last observations, the hold-out fit and the forecast fan
#' # a longer hold-out, 80% intervals
#' predict(r, h = 10, holdout = 0.3, level = 0.8)
#' }
predict.causres <- function(object, h = 10L, p = NULL, holdout = 0.2, level = 0.95, refit = TRUE, ...) {
  D <- as.data.frame(object$data); vars <- names(D); n <- nrow(D)
  links_of <- function(r) {                        # significant causes of every series in a causres
    Pr <- r$pairs[r$pairs$significant %in% TRUE, , drop = FALSE]
    stats::setNames(lapply(vars, function(v) unique(Pr$cause[Pr$effect == v])), vars)
  }
  lag_of <- function(r, p) {
    Pr <- r$pairs[r$pairs$significant %in% TRUE, , drop = FALSE]
    if (is.null(p)) p <- suppressWarnings(max(Pr$lag, na.rm = TRUE))
    max(1L, as.integer(if (is.finite(p)) p else 1L))
  }
  p_full <- lag_of(object, p); causes_of <- links_of(object)
  holdout <- max(0, min(holdout, 0.5)); n_train <- floor(n * (1 - holdout))
  if (holdout > 0 && n_train <= 3 * p_full + 5) stop("Series too short for the hold-out evaluation.", call. = FALSE)
  fit <- function(target, preds, idx, p) {
    d <- .lag_design(D[idx, , drop = FALSE], p, target, preds)
    f <- stats::lm.fit(cbind(1, d$X), d$y)
    f$coefficients[is.na(f$coefficients)] <- 0        # aliased (e.g. constant) columns drop out
    f
  }
  pred_one <- function(coef, D_, t, preds, p) {       # one-step prediction for time t from rows < t
    x <- unlist(lapply(preds, function(v) D_[[v]][t - seq_len(p)]))
    sum(coef * c(1, x))
  }
  # ---- hold-out evaluation: links (and lag) chosen on the training part only
  ev <- NULL; hold <- NULL; eval_links <- NULL; eval_source <- "none"
  if (holdout > 0) {
    tr <- seq_len(n_train)
    if (refit && !is.null(object$args)) {
      blocks <- !setequal(object$causes, vars) || !setequal(object$effects, vars)
      xs <- if (blocks) D[tr, object$causes, drop = FALSE] else D[tr, , drop = FALSE]
      ys <- if (blocks) D[tr, setdiff(object$effects, object$causes), drop = FALSE] else NULL
      if (!is.null(ys) && !ncol(ys)) ys <- NULL
      r_tr <- suppressMessages(do.call(rcausality, c(list(x = xs, y = ys), object$args)))
      eval_links <- links_of(r_tr); p_ev <- lag_of(r_tr, p); eval_source <- "training part (refit)"
    } else { eval_links <- causes_of; p_ev <- p_full; eval_source <- "full sample (optimistic)" }
    if (n_train <= 3 * p_ev + 5) stop("Series too short for the hold-out evaluation.", call. = FALSE)
    ev <- list(); hold <- list(); tt <- (n_train + 1):n
    for (v in vars) {
      base <- fit(v, v, tr, p_ev); caus <- fit(v, c(v, eval_links[[v]]), tr, p_ev)
      pb <- vapply(tt, function(t) pred_one(base$coefficients, D, t, v, p_ev), 0)
      pc <- vapply(tt, function(t) pred_one(caus$coefficients, D, t, c(v, eval_links[[v]]), p_ev), 0)
      e_b <- D[[v]][tt] - pb; e_c <- D[[v]][tt] - pc
      hold[[v]] <- data.frame(series = v, time = tt, observed = D[[v]][tt], pred_own = pb, pred_causal = pc, stringsAsFactors = FALSE)
      ev[[v]] <- data.frame(series = v, causes = paste(eval_links[[v]], collapse = ","),
                            rmse_baseline = sqrt(mean(e_b^2)), rmse_causal = sqrt(mean(e_c^2)),
                            mae_baseline = mean(abs(e_b)), mae_causal = mean(abs(e_c)), stringsAsFactors = FALSE)
    }
    ev <- do.call(rbind, ev); rownames(ev) <- NULL
    ev$improvement <- ifelse(ev$rmse_baseline > 0, 1 - ev$rmse_causal / ev$rmse_baseline, 0)
    hold <- do.call(rbind, hold)
  }
  # ---- final model on the full sample with the links of `object`
  p <- p_full; coefs <- list(); resid <- list()
  for (v in vars) {
    full <- fit(v, c(v, causes_of[[v]]), seq_len(n), p)
    coefs[[v]] <- list(coef = full$coefficients, preds = c(v, causes_of[[v]]))
    resid[[v]] <- full$residuals
  }
  # restricted VAR in companion form: A_l[i, j] = coefficient of series j at lag l in the equation of series i
  K <- length(vars); A <- array(0, c(K, K, p), dimnames = list(vars, vars, NULL))
  for (v in vars) {
    cf <- coefs[[v]]$coef[-1L]; names(cf) <- as.vector(outer(seq_len(p), coefs[[v]]$preds, function(l, w) paste(w, l)))
    for (w in coefs[[v]]$preds) for (l in seq_len(p)) A[v, w, l] <- cf[[paste(w, l)]]
  }
  Fm <- matrix(0, K * p, K * p); for (l in seq_len(p)) Fm[seq_len(K), (l - 1) * K + seq_len(K)] <- A[, , l]
  if (p > 1) Fm[K + seq_len(K * (p - 1)), seq_len(K * (p - 1))] <- diag(K * (p - 1))
  R <- do.call(cbind, resid); Sigma <- stats::cov(R)          # residual covariance (same length for all series)
  J <- cbind(diag(K), matrix(0, K, K * (p - 1)))
  # recursive h-step point forecasts and the Luetkepohl MSE recursion
  Dx <- D; z <- stats::qnorm(1 - (1 - level) / 2)
  fc <- list(); MSE <- matrix(0, K, K); Fi <- diag(K * p)
  for (k in seq_len(h)) {
    new <- vapply(vars, function(v) pred_one(coefs[[v]]$coef, Dx, nrow(Dx) + 1L, coefs[[v]]$preds, p), 0)
    Dx <- rbind(Dx, as.list(new))
    Phi <- J %*% Fi %*% t(J); MSE <- MSE + Phi %*% Sigma %*% t(Phi); Fi <- Fi %*% Fm
    se <- sqrt(pmax(diag(MSE), 0))
    for (i in seq_len(K)) {
      v <- vars[i]
      fc[[length(fc) + 1L]] <- data.frame(series = v, horizon = k, point = new[[v]], lower = new[[v]] - z * se[i],
                                          upper = new[[v]] + z * se[i], se = se[i], stringsAsFactors = FALSE)
    }
  }
  structure(list(forecasts = do.call(rbind, fc), evaluation = ev, holdout_fit = hold, evaluation_links = eval_links,
                 evaluation_source = eval_source, fitted_causes = causes_of, coefficients = A, sigma = Sigma,
                 p = p, h = h, level = level, holdout = holdout, data = D, method = object$settings$method), class = "causfc")
}

#' @export
print.causfc <- function(x, digits = 3, ...) {
  cat(sprintf("Causality-informed forecasts (restricted VAR from the '%s' links; p = %d lags, h = %d)\n", x$method, x$p, x$h))
  if (is.null(x$evaluation)) cat("  no hold-out evaluation (holdout = 0)\n") else {
    cat(sprintf("  hold-out evaluation (last %.0f%% of the sample, rolling one-step forecasts; links selected on the %s):\n",
                100 * x$holdout, x$evaluation_source))
    E <- x$evaluation
    print(data.frame(series = E$series, causes = ifelse(E$causes == "", "-", E$causes),
                     RMSE.own = round(E$rmse_baseline, digits), RMSE.causal = round(E$rmse_causal, digits),
                     improvement = sprintf("%+.1f%%", 100 * E$improvement)), row.names = FALSE, right = FALSE)
  }
  cat("  forecast model: links of the full-sample analysis; causes per series:",
      paste(sprintf("%s <- {%s}", names(x$fitted_causes), vapply(x$fitted_causes, function(v) paste(v, collapse = ","), "")), collapse = "; "), "\n")
  cat(sprintf("\n  forecasts: point [%d%% interval], one column per series (first %d horizons; all in $forecasts):\n",
              round(100 * x$level), min(x$h, 10L)))
  Fc <- x$forecasts; vars <- names(x$fitted_causes)
  W <- do.call(cbind, lapply(vars, function(v) {
    f <- Fc[Fc$series == v, ]; f <- f[order(f$horizon), ]
    sprintf("%s [%s, %s]", formatC(f$point, digits = digits - 1, format = "f"),
            formatC(f$lower, digits = digits - 1, format = "f"), formatC(f$upper, digits = digits - 1, format = "f"))
  }))
  W <- data.frame(h = seq_len(x$h), W, stringsAsFactors = FALSE); names(W) <- c("h", vars)
  print(utils::head(W, 10L), row.names = FALSE, right = FALSE)
  invisible(x)
}

#' @export
summary.causfc <- function(object, ...) {
  E <- object$evaluation
  structure(list(evaluation = E, mean_improvement = if (is.null(E)) NA_real_ else mean(E$improvement), p = object$p, h = object$h,
                 best = if (is.null(E)) NA_character_ else E$series[which.max(E$improvement)]), class = "summary.causfc")
}

#' @export
print.summary.causfc <- function(x, ...) {
  cat(sprintf("Causal forecasting summary: p = %d, h = %d\n", x$p, x$h))
  if (is.null(x$evaluation)) { cat("  no hold-out evaluation\n"); return(invisible(x)) }
  cat(sprintf("  mean RMSE improvement of the causal model over own lags: %+.1f%% (largest for %s)\n",
              100 * x$mean_improvement, x$best))
  print(data.frame(series = x$evaluation$series, causes = ifelse(x$evaluation$causes == "", "-", x$evaluation$causes),
                   improvement = sprintf("%+.1f%%", 100 * x$evaluation$improvement)), row.names = FALSE, right = FALSE)
  invisible(x)
}

#' @rdname predict.causres
#' @param x a `causfc` object.
#' @param n_last number of observed values shown before the forecasts.
#' @param main title.
#' @export
plot.causfc <- function(x, n_last = 60L, main = NULL, ...) {
  D <- x$data; vars <- names(D); n <- nrow(D); Fc <- x$forecasts; H <- x$holdout_fit; E <- x$evaluation
  op <- graphics::par(mfrow = grDevices::n2mfrow(length(vars)), mar = c(4, 4, 3, 1)); on.exit(graphics::par(op))
  idx <- max(1, n - n_last + 1):n
  for (v in vars) {
    f <- Fc[Fc$series == v, ]; hv <- if (is.null(H)) H[0, ] else H[H$series == v & H$time >= min(idx), ]
    ylim <- range(c(D[[v]][idx], f$lower, f$upper, hv$pred_own, hv$pred_causal), na.rm = TRUE)
    impr <- if (is.null(E)) "" else sprintf("; hold-out RMSE %+.0f%%", 100 * E$improvement[E$series == v])
    graphics::plot(idx, D[[v]][idx], type = "l", xlim = c(min(idx), n + x$h), ylim = ylim, xlab = "time", ylab = v,
                   main = if (is.null(main)) sprintf("%s  (causes: %s%s)", v,
                                                     if (length(x$fitted_causes[[v]])) paste(x$fitted_causes[[v]], collapse = ", ") else "none", impr) else main,
                   font.main = 2, cex.main = 0.95, col = "grey30", lwd = 1.3)
    if (!is.null(hv) && nrow(hv)) {
      graphics::lines(hv$time, hv$pred_own, col = "#B2403A", lty = 2, lwd = 1.2)
      graphics::lines(hv$time, hv$pred_causal, col = "#2B5C8A", lwd = 1.5)
    }
    graphics::polygon(c(n + f$horizon, rev(n + f$horizon)), c(f$lower, rev(f$upper)),
                      col = grDevices::adjustcolor("#1B7A43", 0.2), border = NA)
    graphics::lines(c(n, n + f$horizon), c(D[[v]][n], f$point), col = "#1B7A43", lwd = 2)
    graphics::points(n + f$horizon, f$point, col = "#1B7A43", pch = 16, cex = 0.7)
    graphics::abline(v = n, lty = 3, col = "grey60")
    if (v == vars[1])
      graphics::legend("topleft", legend = c("observed", "own lags (one step)", "with causes (one step)",
                                             sprintf("forecast, %d%% interval", round(100 * x$level))),
                       col = c("grey30", "#B2403A", "#2B5C8A", "#1B7A43"), lty = c(1, 2, 1, 1), lwd = c(1.3, 1.2, 1.5, 2),
                       bty = "n", cex = 0.65)
  }
  invisible(x)
}

#' Causality-informed forecasts (convenience wrapper)
#'
#' Runs [rcausality()] with the given method and then [predict.causres()].
#' @param x,y data as for [rcausality()].
#' @param method causality method used to select the causes.
#' @param h forecast horizon.
#' @param ... further arguments to [rcausality()] (e.g. `seed`) and to
#'   [predict.causres()] (`p`, `holdout`, `level`).
#' @return A `causfc` object.
#' @export
#' @examples
#' # one call: causal analysis + restricted-VAR forecast of every series
#' d <- sim_var(400, "x_to_y", seed = 1)
#' forecast_causal(d, h = 5)
forecast_causal <- function(x, y = NULL, method = "granger", h = 10L, ...) {
  dots <- list(...)
  pa <- intersect(names(dots), c("p", "holdout", "level"))
  r <- do.call(rcausality, c(list(x = x, y = y, method = method), dots[setdiff(names(dots), pa)]))
  do.call(predict.causres, c(list(object = r, h = h), dots[pa]))
}
