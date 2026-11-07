# Causality-informed forecasting ---------------------------------------------
# The significant links of a causres define, for every series, its causes.
# Every series is then predicted from its own lags and the lags of its causes
# by one of three forecast models:
#   "var"  -- the restricted VAR (least squares, analytic MSE intervals);
#   "bvar" -- the same equations with the hierarchical Minnesota prior
#             (short series; posterior-predictive simulation intervals);
#   "rf"   -- random forests on the same lags (nonlinear; recursive
#             simulation with bootstrapped out-of-bag residuals).
# The gain over the own-lags-only model is measured on a hold-out part of the
# sample with links selected on the training part only (refit = TRUE).

# regressor rows for the times `t` (intercept, then the lags 1..p of every
# predictor series in turn)
.fc_X <- function(D, preds, t, p) {
  cbind(1, do.call(cbind, lapply(preds, function(w) matrix(vapply(seq_len(p), function(l) D[[w]][t - l], numeric(length(t))),
                                                            nrow = length(t)))))
}

# ---- the three forecast models: fit on rows `idx` (times), predict rows ---------

.fc_fit <- function(model, D, target, preds, idx, p, settings) {
  t <- idx[idx > p]; X <- .fc_X(D, preds, t, p); y <- D[[target]][t]
  switch(model,
    var = {
      f <- stats::lm.fit(X, y); cf <- f$coefficients; cf[is.na(cf)] <- 0      # aliased (e.g. constant) columns drop out
      list(coef = cf, resid = as.numeric(y - X %*% cf), preds = preds)
    },
    bvar = {
      # standardised series inside D (see predict.causres); Minnesota base
      # variances: own lags 1 / l^(2 decay), the causes' lags cross^2 times that
      ctl <- settings$control; l_dec <- seq_len(p)^(-2 * ctl$decay)
      dbase <- c(l_dec, rep(ctl$cross^2 * l_dec, length(preds) - 1L))
      b <- .minnesota_posterior(y, X[, 1, drop = FALSE], X[, -1, drop = FALSE], dbase, ctl)
      list(coef = b$mean, resid = as.numeric(y - X %*% b$mean), preds = preds, post = b)
    },
    rf = {
      fit <- .rf_fit(X[, -1, drop = FALSE], y, settings$rf)
      list(fit = fit, resid = y - fit$predictions, preds = preds)
    })
}

.fc_predict <- function(model, fit, X) {
  if (model == "rf") stats::predict(fit$fit, data = as.data.frame(X[, -1, drop = FALSE]), num.threads = 1L)$predictions
  else as.numeric(X %*% fit$coef)
}

# posterior of the coefficients under the hierarchical Minnesota prior: the
# mixture over the tightness grid of normal-inverse-gamma posteriors; mean
# and a sampler of (beta, sigma)
.minnesota_posterior <- function(y, Xfix, Xs, dbase, control) {
  N <- length(y); a0 <- 0.01; b0 <- 0.01; X <- cbind(Xfix, Xs)
  if (identical(control$lambda, "auto")) {
    g <- .gamma_from_mode_sd(control$lambda_hyper[["mode"]], control$lambda_hyper[["sd"]])
    lam <- exp(seq(log(0.005), log(10), length.out = 60L))
    lw <- stats::dgamma(lam, shape = g[["shape"]], scale = g[["scale"]], log = TRUE) + log(lam)
  } else { lam <- as.numeric(control$lambda); lw <- 0 }
  comp <- lapply(lam, function(l) {
    d <- c(rep(1e4, ncol(Xfix)), l^2 * dbase)
    A <- crossprod(X) + diag(1 / d, ncol(X)); R <- chol(A); Xy <- crossprod(X, y)
    m <- backsolve(R, forwardsolve(t(R), Xy)); b <- b0 + (sum(y^2) - sum(m * Xy)) / 2
    list(m = as.numeric(m), R = R, b = max(b, 1e-12),
         lml = -sum(log(diag(R))) - 0.5 * sum(log(d)) - (a0 + N / 2) * log(max(b, 1e-12)))
  })
  lml <- vapply(comp, `[[`, 0, "lml") + lw; w <- exp(lml - max(lml)); w <- w / sum(w)
  mean_b <- Reduce(`+`, Map(function(cmp, wi) wi * cmp$m, comp, w))
  draw <- function(k) {
    li <- sample.int(length(lam), k, replace = TRUE, prob = w)
    out <- matrix(0, k, ncol(X) + 1L)
    for (i in seq_len(k)) {
      cmp <- comp[[li[i]]]; s2 <- 1 / stats::rgamma(1, a0 + N / 2, rate = cmp$b)
      out[i, ] <- c(cmp$m + sqrt(s2) * backsolve(cmp$R, stats::rnorm(ncol(X))), sqrt(s2))
    }
    out
  }
  list(mean = mean_b, draw = draw, lambda = sum(w * lam))
}

#' Causality-informed forecasts from a `causres` object
#'
#' Predicts every series from its own `p` lags and the `p` lags of its
#' significant causes (the restricted model implied by the detected links),
#' compares this model with the own-lags-only model by rolling one-step
#' forecasts on the last `holdout` fraction of the sample, and produces
#' `h`-step-ahead forecasts with prediction intervals. The links may come
#' from any method of [rcausality()]; the forecast model is chosen by
#' `model`:
#' * `"var"`: the restricted VAR estimated by least squares; recursive point
#'   forecasts and normal intervals from the multi-step MSE recursion of the
#'   VAR (Luetkepohl 2005, Sect. 2.2.2);
#' * `"bvar"`: the same equations with the hierarchical Minnesota prior of
#'   [bayes_control()] (`prior = "minnesota"`; tightness integrated out),
#'   which stabilises the coefficients of short series; intervals from
#'   `n_sim` posterior-predictive paths (coefficients, error scale and
#'   correlated shocks simulated);
#' * `"rf"`: random forests (package \pkg{ranger}) on the same lags, for
#'   nonlinear dynamics; intervals from `n_sim` recursive paths with
#'   bootstrapped out-of-bag residuals (resampled jointly across the series).
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
#' @param model forecast model: `"var"`, `"bvar"` or `"rf"` (see Details).
#' @param control a [bayes_control()] list for `model = "bvar"` (its Minnesota
#'   settings are used).
#' @param rf_settings forest settings for `model = "rf"` (`num_trees`, `mtry`,
#'   `min_node_size`).
#' @param n_sim number of simulated paths for the intervals of `"bvar"` and
#'   `"rf"`.
#' @param ... unused.
#' @return An object of class `causfc`: `forecasts` (series, horizon, point,
#'   lower, upper, se), `evaluation` (per series: RMSE and MAE of the baseline
#'   and of the causal model on the hold-out, the relative RMSE improvement,
#'   the causes used), `holdout_fit` (observed and one-step predicted values
#'   of both models on the hold-out), `fitted_causes`, `coefficients` (the
#'   restricted VAR coefficient array, series x series x lag; `"var"` and
#'   `"bvar"`), `sigma` (residual covariance), `model`, `p`, `data`. Methods:
#'   `print`, `summary`, `plot`.
#' @references
#' Granger, C. W. J. (1969). Investigating causal relations by econometric models and cross-spectral methods. *Econometrica*, 37(3), 424-438. \doi{10.2307/1912791}
#'
#' Luetkepohl, H. (2005). *New Introduction to Multiple Time Series Analysis*. Springer, Berlin. \doi{10.1007/978-3-540-27752-1}
#'
#' Giannone, D., Lenza, M. and Primiceri, G. E. (2015). Prior selection for vector autoregressions. *The Review of Economics and Statistics*, 97(2), 436-451. \doi{10.1162/REST_a_00483}
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
#' # a short series: the Minnesota-prior model (simulated intervals)
#' ds <- sim_var(60, "x_to_y", seed = 2)
#' predict(rcausality(ds), h = 5, model = "bvar", n_sim = 200)
#' }
predict.causres <- function(object, h = 10L, p = NULL, holdout = 0.2, level = 0.95, refit = TRUE,
                            model = c("var", "bvar", "rf"), control = bayes_control(prior = "minnesota"),
                            rf_settings = list(num_trees = 300L, mtry = NULL, min_node_size = 5L), n_sim = 500L, ...) {
  model <- match.arg(model)
  if (model == "rf") .need_ranger("predict(model = \"rf\")")
  .check_number(h, "h", lower = 1, integer = TRUE); .check_number(level, "level", lower = 0, upper = 1, strict = TRUE)
  D0 <- as.data.frame(object$data); vars <- names(D0); n <- nrow(D0)
  # the Minnesota model works on standardised series (back-transformed at the end)
  mu <- vapply(D0, mean, 0); sdv <- vapply(D0, stats::sd, 0); sdv[!is.finite(sdv) | sdv <= 0] <- 1
  D <- if (model == "bvar") as.data.frame(Map(function(v, m, s) (v - m) / s, D0, mu, sdv)) else D0
  settings <- list(control = control, rf = utils::modifyList(list(num_trees = 300L, mtry = NULL, min_node_size = 5L),
                                                             rf_settings, keep.null = TRUE))
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
  back <- function(v, x) if (model == "bvar") x * sdv[[v]] + mu[[v]] else x
  # ---- hold-out evaluation: links (and lag) chosen on the training part only
  ev <- NULL; hold <- NULL; eval_links <- NULL; eval_source <- "none"
  if (holdout > 0) {
    tr <- seq_len(n_train)
    if (refit && !is.null(object$args)) {
      r_tr <- .rerun_object(object, D0[tr, , drop = FALSE])
      eval_links <- links_of(r_tr); p_ev <- lag_of(r_tr, p); eval_source <- "training part (refit)"
    } else { eval_links <- causes_of; p_ev <- p_full; eval_source <- "full sample (optimistic)" }
    if (n_train <= 3 * p_ev + 5) stop("Series too short for the hold-out evaluation.", call. = FALSE)
    tt <- (n_train + 1):n; ev <- list(); hold <- list()
    for (v in vars) {
      base <- .fc_fit(model, D, v, v, tr, p_ev, settings)
      caus <- if (length(eval_links[[v]])) .fc_fit(model, D, v, c(v, eval_links[[v]]), tr, p_ev, settings) else base
      pb <- back(v, .fc_predict(model, base, .fc_X(D, v, tt, p_ev)))
      pc <- back(v, .fc_predict(model, caus, .fc_X(D, c(v, eval_links[[v]]), tt, p_ev)))
      obs <- D0[[v]][tt]; e_b <- obs - pb; e_c <- obs - pc
      hold[[v]] <- data.frame(series = v, time = tt, observed = obs, pred_own = pb, pred_causal = pc, stringsAsFactors = FALSE)
      ev[[v]] <- data.frame(series = v, causes = paste(eval_links[[v]], collapse = ","),
                            rmse_baseline = sqrt(mean(e_b^2)), rmse_causal = sqrt(mean(e_c^2)),
                            mae_baseline = mean(abs(e_b)), mae_causal = mean(abs(e_c)), stringsAsFactors = FALSE)
    }
    ev <- do.call(rbind, ev); rownames(ev) <- NULL
    ev$improvement <- ifelse(ev$rmse_baseline > 0, 1 - ev$rmse_causal / ev$rmse_baseline, 0)
    hold <- do.call(rbind, hold); rownames(hold) <- NULL
  }
  # ---- final model on the full sample with the links of `object`
  p <- p_full; K <- length(vars)
  fits <- stats::setNames(lapply(vars, function(v) .fc_fit(model, D, v, c(v, causes_of[[v]]), seq_len(n), p, settings)), vars)
  R <- do.call(cbind, lapply(fits, `[[`, "resid")); Sigma <- stats::cov(R)
  A <- NULL
  if (model != "rf") {
    # restricted VAR coefficients: A_l[i, j] = coefficient of series j at lag l in the equation of series i
    A <- array(0, c(K, K, p), dimnames = list(vars, vars, NULL))
    for (v in vars) {
      cf <- fits[[v]]$coef[-1L]
      for (j in seq_along(fits[[v]]$preds)) A[v, fits[[v]]$preds[j], ] <- cf[(j - 1L) * p + seq_len(p)]
    }
  }
  z <- stats::qnorm(1 - (1 - level) / 2)
  if (model == "var") {
    Fm <- matrix(0, K * p, K * p); for (l in seq_len(p)) Fm[seq_len(K), (l - 1) * K + seq_len(K)] <- A[, , l]
    if (p > 1) Fm[K + seq_len(K * (p - 1)), seq_len(K * (p - 1))] <- diag(K * (p - 1))
    J <- cbind(diag(K), matrix(0, K, K * (p - 1)))
    Dx <- D; fc <- list(); MSE <- matrix(0, K, K); Fi <- diag(K * p)
    for (k in seq_len(h)) {
      new <- vapply(vars, function(v) .fc_predict(model, fits[[v]], .fc_X(Dx, fits[[v]]$preds, nrow(Dx) + 1L, p)), 0)
      Dx <- rbind(Dx, as.list(new))
      Phi <- J %*% Fi %*% t(J); MSE <- MSE + Phi %*% Sigma %*% t(Phi); Fi <- Fi %*% Fm
      se <- sqrt(pmax(diag(MSE), 0))
      fc[[k]] <- data.frame(series = vars, horizon = k, point = new, lower = new - z * se, upper = new + z * se, se = se,
                            stringsAsFactors = FALSE)
    }
    fc <- do.call(rbind, fc)
  } else {
    fc <- .fc_simulate(model, fits, D, vars, p, h, n_sim, R, level)
    for (v in vars) {                                # back to the original units
      k <- fc$series == v
      fc$point[k] <- back(v, fc$point[k]); fc$lower[k] <- back(v, fc$lower[k]); fc$upper[k] <- back(v, fc$upper[k])
      if (model == "bvar") fc$se[k] <- fc$se[k] * sdv[[v]]
    }
    if (model == "bvar") Sigma <- Sigma * outer(sdv, sdv)
  }
  fc <- fc[order(match(fc$series, vars), fc$horizon), ]; rownames(fc) <- NULL
  structure(list(forecasts = fc, evaluation = ev, holdout_fit = hold, evaluation_links = eval_links,
                 evaluation_source = eval_source, fitted_causes = causes_of, coefficients = A, sigma = Sigma,
                 model = model, p = p, h = h, level = level, holdout = holdout, data = D0, method = object$settings$method),
            class = "causfc")
}

# simulated forecast paths of the "bvar" and "rf" models: point forecast =
# mean of the paths, interval = quantiles
.fc_simulate <- function(model, fits, D, vars, p, h, n_sim, R, level) {
  K <- length(vars); n <- nrow(D)
  hist <- as.matrix(D[(n - p + 1):n, vars, drop = FALSE])                 # last p observations
  paths <- array(NA_real_, c(n_sim, p + h, K), dimnames = list(NULL, NULL, vars))
  for (s in seq_len(n_sim)) paths[s, seq_len(p), ] <- hist
  if (model == "bvar") {
    draws <- lapply(vars, function(v) fits[[v]]$post$draw(n_sim)); names(draws) <- vars
    Cr <- stats::cor(R); Cr[!is.finite(Cr)] <- 0; diag(Cr) <- 1
    Lc <- tryCatch(t(chol(Cr)), error = function(e) diag(K))
  } else Rc <- R[stats::complete.cases(R), , drop = FALSE]
  for (k in seq_len(h)) {
    tcur <- p + k
    shock <- if (model == "bvar") t(Lc %*% matrix(stats::rnorm(K * n_sim), K, n_sim)) else Rc[sample.int(nrow(Rc), n_sim, TRUE), , drop = FALSE]
    for (j in seq_len(K)) {
      v <- vars[j]; pr <- fits[[v]]$preds
      X <- cbind(1, do.call(cbind, lapply(pr, function(w) paths[, tcur - seq_len(p), w, drop = FALSE][, , 1, drop = TRUE])))
      X <- matrix(X, nrow = n_sim)
      if (model == "bvar") {
        B <- draws[[v]]; mean_k <- rowSums(X * B[, seq_len(ncol(X)), drop = FALSE])
        paths[, tcur, v] <- mean_k + B[, ncol(B)] * shock[, j]
      } else paths[, tcur, v] <- .fc_predict("rf", fits[[v]], X) + shock[, j]
    }
  }
  a <- (1 - level) / 2
  do.call(rbind, lapply(seq_len(h), function(k) do.call(rbind, lapply(vars, function(v) {
    x <- paths[, p + k, v]
    data.frame(series = v, horizon = k, point = mean(x), lower = unname(stats::quantile(x, a)),
               upper = unname(stats::quantile(x, 1 - a)), se = stats::sd(x), stringsAsFactors = FALSE)
  }))))
}

# the analysis of `object` repeated on new data with the same arguments
.rerun_object <- function(object, data) .rerun(object, data, seed = object$args$seed)

#' @export
print.causfc <- function(x, digits = 3, ...) {
  mdl <- c(var = "restricted VAR", bvar = "restricted Bayesian VAR, Minnesota prior", rf = "random forests")[[if (is.null(x$model)) "var" else x$model]]
  cat(sprintf("Causality-informed forecasts (%s from the '%s' links; p = %d lags, h = %d)\n", mdl, x$method, x$p, x$h))
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
#' Runs [rcausality()] with the given method and then [predict.causres()]
#' with the given forecast model: the causality method decides which series
#' enter the forecast equation of which, the forecast model how they are
#' combined.
#' @param x,y data as for [rcausality()].
#' @param method causality method used to select the causes.
#' @param h forecast horizon.
#' @param model forecast model: `"var"`, `"bvar"` or `"rf"` (see
#'   [predict.causres()]).
#' @param ... further arguments to [rcausality()] (e.g. `seed`) and to
#'   [predict.causres()] (`p`, `holdout`, `level`, `refit`, `control`,
#'   `rf_settings`, `n_sim`).
#' @return A `causfc` object.
#' @export
#' @examples
#' # one call: causal analysis + restricted-VAR forecast of every series
#' d <- sim_var(400, "x_to_y", seed = 1)
#' forecast_causal(d, h = 5)
#' \donttest{
#' # the links from the transfer entropy, a Minnesota-prior forecast model
#' forecast_causal(d, method = "te", model = "bvar", h = 5, n_perm = 100, seed = 1, n_sim = 200)
#' }
forecast_causal <- function(x, y = NULL, method = "granger", h = 10L, model = c("var", "bvar", "rf"), ...) {
  model <- match.arg(model)
  dots <- list(...)
  pa <- intersect(names(dots), c("p", "holdout", "level", "refit", "control", "rf_settings", "n_sim"))
  ra <- dots[setdiff(names(dots), pa)]
  if ("control" %in% names(dots) && isTRUE(dots$bayes)) ra$control <- dots$control
  r <- do.call(rcausality, c(list(x = x, y = y, method = method), ra))
  do.call(predict.causres, c(list(object = r, h = h, model = model), dots[pa]))
}

# Diebold-Mariano test of equal accuracy (squared-error loss) of two series
# of h-step forecast errors, with the Harvey-Leybourne-Newbold correction
.dm_test <- function(e1, e2, h = 1L) {
  d <- e1^2 - e2^2; n <- length(d)
  if (n < 5 || stats::sd(d) == 0) return(c(statistic = NA_real_, p_value = NA_real_))
  g <- function(k) sum((d[(k + 1):n] - mean(d)) * (d[seq_len(n - k)] - mean(d))) / n
  v <- (g(0) + 2 * sum(vapply(seq_len(h - 1L), g, 0))) / n
  if (!is.finite(v) || v <= 0) return(c(statistic = NA_real_, p_value = NA_real_))
  hln <- sqrt((n + 1 - 2 * h + h * (h - 1) / n) / n)
  s <- hln * mean(d) / sqrt(v)
  c(statistic = s, p_value = 2 * stats::pt(-abs(s), n - 1))
}

#' Compare causality methods by the forecasts they enable
#'
#' For every method, the causes of every series are selected on the training
#' part of the sample (the first `1 - holdout` share) and the forecast model
#' with those causes is evaluated by rolling one-step forecasts on the
#' hold-out part; the reference is the own-lags-only model (`"none"`). The
#' relative RMSE gain of every method and series is reported with the
#' Diebold-Mariano test of equal accuracy against the reference (squared-error
#' loss, Harvey-Leybourne-Newbold correction). Since the links are chosen
#' without the hold-out data, the comparison measures the out-of-sample
#' usefulness of the detected causal structure.
#'
#' @param x,y data as for [rcausality()]; for `plot()`, `x` is a `causfcomp`
#'   object.
#' @param methods causality methods (see [rcausality()]).
#' @param model forecast model (see [predict.causres()]).
#' @param holdout share of the sample in the hold-out.
#' @param p lag order of the forecast models (`NULL` = per method, the
#'   largest lag of its links).
#' @param args optional named list of per-method argument lists for
#'   [rcausality()].
#' @param cores CPU cores: the methods run in parallel (see
#'   [Rcausality-parallel]).
#' @param ... arguments to every [rcausality()] call (e.g. `seed`, `n_perm`)
#'   and to [predict.causres()] (`control`, `rf_settings`).
#' @return An object of class `causfcomp`: `table` (one row per method and
#'   series: causes, RMSE, gain over the own-lags model, Diebold-Mariano
#'   statistic and p-value), `summary` (mean gain per method), `holdout`,
#'   `model`. Methods: `print`, `plot`.
#' @references
#' Diebold, F. X. and Mariano, R. S. (1995). Comparing predictive accuracy. *Journal of Business & Economic Statistics*, 13(3), 253-263. \doi{10.1080/07350015.1995.10524599}
#'
#' Harvey, D., Leybourne, S. and Newbold, P. (1997). Testing the equality of prediction mean squared errors. *International Journal of Forecasting*, 13(2), 281-291. \doi{10.1016/S0169-2070(96)00719-4}
#' @export
#' @examples
#' d <- sim_var(300, n_series = 3, seed = 2)
#' fcmp <- forecast_compare(d, methods = c("granger", "toda_yamamoto"))
#' fcmp
#' \donttest{
#' plot(fcmp)
#' forecast_compare(d, methods = c("granger", "te", "bayes_granger"), n_perm = 100,
#'                  seed = 1, cores = 2)
#' }
forecast_compare <- function(x, y = NULL, methods = c("granger", "bayes_granger", "te"), model = c("var", "bvar", "rf"),
                             holdout = 0.25, p = NULL, args = list(), cores = getOption("Rcausality.cores", 1L), ...) {
  model <- match.arg(model)
  methods <- match.arg(methods, .method_names, several.ok = TRUE)
  .check_number(holdout, "holdout", lower = 0.05, upper = 0.5)
  dots <- list(...)
  pa <- intersect(names(dots), c("control", "rf_settings", "n_sim"))
  ra <- dots[setdiff(names(dots), pa)]
  fr <- .as_series_frame(x, y); D <- if (is.null(fr$Y)) fr$X else cbind(fr$X, fr$Y)
  D[] <- lapply(D, .prep_na); n <- nrow(D); n_train <- floor(n * (1 - holdout)); tr <- seq_len(n_train); tt <- (n_train + 1):n
  vars <- names(D)
  settings <- list(control = if (is.null(dots$control)) bayes_control(prior = "minnesota") else dots$control,
                   rf = utils::modifyList(list(num_trees = 300L, mtry = NULL, min_node_size = 5L),
                                          if (is.null(dots$rf_settings)) list() else dots$rf_settings, keep.null = TRUE))
  Ds <- if (model == "bvar") as.data.frame(lapply(D, function(v) (v - mean(v[tr])) / stats::sd(v[tr]))) else D
  back <- function(v, z) if (model == "bvar") z * stats::sd(D[[v]][tr]) + mean(D[[v]][tr]) else z
  sc <- .split_cores(length(methods), cores)
  runs <- .plapply(methods, function(m) {
    a <- c(list(x = D[tr, , drop = FALSE], method = m, cores = sc$inner), ra, if (!is.null(args[[m]])) args[[m]])
    a <- a[!duplicated(names(a), fromLast = TRUE)]
    r <- suppressMessages(do.call(rcausality, a))
    Pr <- r$pairs[r$pairs$significant %in% TRUE, , drop = FALSE]
    list(causes = stats::setNames(lapply(vars, function(v) unique(Pr$cause[Pr$effect == v])), vars),
         p = as.integer(max(1L, suppressWarnings(max(c(1L, Pr$lag), na.rm = TRUE)))))
  }, cores = sc$outer)
  names(runs) <- methods
  p_all <- if (is.null(p)) max(vapply(runs, `[[`, 0L, "p")) else as.integer(p)
  err_of <- function(v, preds) {
    f <- .fc_fit(model, Ds, v, preds, tr, p_all, settings)
    D[[v]][tt] - back(v, .fc_predict(model, f, .fc_X(Ds, preds, tt, p_all)))
  }
  base <- stats::setNames(lapply(vars, function(v) err_of(v, v)), vars)
  rows <- list()
  for (v in vars) rows[[length(rows) + 1L]] <- data.frame(method = "none", series = v, causes = "-",
                                                         rmse = sqrt(mean(base[[v]]^2)), gain = 0, dm = NA_real_, p_dm = NA_real_,
                                                         stringsAsFactors = FALSE)
  for (m in methods) for (v in vars) {
    cs <- runs[[m]]$causes[[v]]
    e <- if (length(cs)) err_of(v, c(v, cs)) else base[[v]]
    dm <- if (length(cs)) .dm_test(base[[v]], e) else c(statistic = NA_real_, p_value = NA_real_)
    rows[[length(rows) + 1L]] <- data.frame(method = m, series = v, causes = if (length(cs)) paste(cs, collapse = ",") else "-",
                                           rmse = sqrt(mean(e^2)), gain = 1 - sqrt(mean(e^2)) / sqrt(mean(base[[v]]^2)),
                                           dm = dm[["statistic"]], p_dm = dm[["p_value"]], stringsAsFactors = FALSE)
  }
  tab <- do.call(rbind, rows); rownames(tab) <- NULL
  sm <- stats::aggregate(gain ~ method, data = tab, FUN = mean)
  sm <- sm[match(c("none", methods), sm$method), ]; rownames(sm) <- NULL
  structure(list(table = tab, summary = sm, holdout = holdout, model = model, p = p_all, n = n, methods = methods),
            class = "causfcomp")
}

#' @export
print.causfcomp <- function(x, ...) {
  mdl <- c(var = "restricted VAR", bvar = "Minnesota-prior VAR", rf = "random forests")[[x$model]]
  cat(sprintf("Forecast comparison of causality methods (%s, p = %d; hold-out = last %.0f%% of n = %d, one-step forecasts)\n",
              mdl, x$p, 100 * x$holdout, x$n))
  T <- x$table
  print(data.frame(method = T$method, series = T$series, causes = T$causes, RMSE = signif(T$rmse, 4),
                   gain = sprintf("%+.1f%%", 100 * T$gain), DM = ifelse(is.na(T$dm), "-", sprintf("%.2f", T$dm)),
                   p.DM = ifelse(is.na(T$p_dm), "-", .fmt_p(T$p_dm))), row.names = FALSE, right = FALSE)
  cat("  gain = RMSE reduction against the own-lags model ('none'); DM = Diebold-Mariano statistic (positive: more accurate)\n")
  cat("  mean gain per method: ", paste(sprintf("%s %+.1f%%", x$summary$method, 100 * x$summary$gain), collapse = ", "), "\n", sep = "")
  invisible(x)
}

#' @rdname forecast_compare
#' @param main optional title.
#' @export
plot.causfcomp <- function(x, main = NULL, ...) {
  T <- x$table[x$table$method != "none", ]
  M <- 100 * tapply(T$gain, list(T$method, T$series), mean); M <- M[x$methods, , drop = FALSE]
  P <- tapply(T$p_dm, list(T$method, T$series), mean); P <- P[x$methods, , drop = FALSE]
  cols <- grDevices::hcl.colors(nrow(M), "Dark 3")
  op <- graphics::par(mar = c(4.2, 4.5, 3, 1)); on.exit(graphics::par(op))
  yl <- range(0, M, na.rm = TRUE); yl <- yl + c(-0.1, 0.15) * diff(yl)
  bp <- graphics::barplot(M, beside = TRUE, col = cols, border = .pal$border, ylim = yl, ylab = "RMSE gain over own lags (%)",
                          main = if (is.null(main)) sprintf("Out-of-sample value of the detected causes (%s)", x$model) else main,
                          font.main = 2, cex.main = 0.95, legend.text = rownames(M),
                          args.legend = list(bty = "n", cex = 0.75, x = "topright"))
  graphics::abline(h = 0, col = "grey40")
  star <- ifelse(is.na(P), "", ifelse(P < 0.01, "**", ifelse(P < 0.05, "*", "")))
  graphics::text(bp, M + ifelse(M >= 0, 1, -1) * 0.03 * diff(yl), star, cex = 1)
  graphics::mtext("* / ** : Diebold-Mariano p < 0.05 / 0.01 against the own-lags model", side = 1, line = 2.8, cex = 0.7, col = "grey30")
  invisible(x)
}
