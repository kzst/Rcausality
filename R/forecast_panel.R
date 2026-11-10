# Causality-informed forecasting of panel data -----------------------------------
# The links of the panel analysis define the causes of every series; the
# forecast equation of a series is pooled over the units: the series are
# centred and scaled within the units (unit fixed effects and unit scales),
# the coefficients of the own lags and of the causes' lags are common.  The
# same three forecast models as for time series ("var", "bvar", "rf").

# units standardised with the means and standard deviations of their first
# `n_train` rows (all rows when NULL)
.pfc_std <- function(U, vars, n_train = NULL) {
  out <- lapply(seq_along(U), function(i) {
    d <- U[[i]][, vars, drop = FALSE]; tr <- if (is.null(n_train)) seq_len(nrow(d)) else seq_len(n_train[i])
    mu <- vapply(d, function(v) mean(v[tr]), 0); s <- vapply(d, function(v) stats::sd(v[tr]), 0); s[!is.finite(s) | s <= 0] <- 1
    list(D = as.data.frame(Map(function(v, m, sv) (v - m) / sv, d, mu, s)), mu = mu, sd = s)
  })
  names(out) <- names(U); out
}

# stacked regressors and response of one forecast equation over the rows
# `rows[[i]]` (time indices) of every unit
.pfc_xy <- function(S, target, preds, rows, p) {
  parts <- lapply(seq_along(S), function(i) {
    t <- rows[[i]]; t <- t[t > p]; if (!length(t)) return(NULL)
    list(X = .fc_X(S[[i]]$D, preds, t, p), y = S[[i]]$D[[target]][t], unit = rep(i, length(t)), t = t)
  })
  parts <- parts[!vapply(parts, is.null, TRUE)]
  list(X = do.call(rbind, lapply(parts, `[[`, "X")), y = unlist(lapply(parts, `[[`, "y")),
       unit = unlist(lapply(parts, `[[`, "unit")), t = unlist(lapply(parts, `[[`, "t")))
}

.panel_links <- function(r, vars) {
  Pr <- r$pairs[r$pairs$significant %in% TRUE, , drop = FALSE]
  list(causes = stats::setNames(lapply(vars, function(v) unique(Pr$cause[Pr$effect == v])), vars),
       p = as.integer(max(1L, suppressWarnings(max(c(1L, Pr$lag), na.rm = TRUE)))))
}

# hold-out errors (original units) of the pooled model of `target` with the
# predictors `preds`: fitted on the first n_train[i] rows of every unit
.pfc_errors <- function(model, S, U, target, preds, n_train, p, settings) {
  tr <- lapply(n_train, seq_len); te <- lapply(seq_along(U), function(i) if (n_train[i] < nrow(U[[i]])) (n_train[i] + 1L):nrow(U[[i]]) else integer(0))
  a <- .pfc_xy(S, target, preds, tr, p); f <- .fc_fit_xy(model, a$X, a$y, preds, p, settings)
  b <- .pfc_xy(S, target, preds, te, p); pr <- .fc_predict(model, f, b$X)
  mu <- vapply(S, function(s) s$mu[[target]], 0)[b$unit]; sdv <- vapply(S, function(s) s$sd[[target]], 0)[b$unit]
  obs <- b$y * sdv + mu; pred <- pr * sdv + mu
  data.frame(unit = names(U)[b$unit], time = b$t, step = b$t - n_train[b$unit], observed = obs, pred = pred, stringsAsFactors = FALSE)
}

#' Causality-informed forecasts of panel data
#'
#' The panel version of [predict.causres()]. The significant links of the
#' panel analysis define the causes of every series. Every series is then
#' predicted, in all units with one pooled equation, from its own `p` lags
#' and the `p` lags of its causes; the series are centred and scaled within
#' the units (unit fixed effects and unit scales). The model is compared
#' with the own-lags-only model on the last `holdout` share of the periods
#' of every unit (rolling one-step forecasts; with `refit = TRUE` the links
#' are selected on the training periods only), and `h`-step forecasts with
#' prediction intervals are produced for every unit. Forecast models as in
#' [predict.causres()]: `"var"` (least squares), `"bvar"` (hierarchical
#' Minnesota prior) and `"rf"` (random forest).
#'
#' @param object a `causpanel` object from [panel_causality()] /
#'   [rcausality()] with `id`.
#' @param h forecast horizon.
#' @param p lag order (default: the largest lag of the significant links).
#' @param holdout share of the periods of every unit used for the
#'   out-of-sample evaluation (`0` = none).
#' @param level level of the prediction intervals.
#' @param refit select the links on the training periods for the evaluation
#'   (see [predict.causres()]).
#' @param model forecast model: `"var"`, `"bvar"` or `"rf"`.
#' @param control,rf_settings,n_sim as in [predict.causres()] (`n_sim` paths
#'   per unit).
#' @param ... unused.
#' @return An object of class `causfcpanel`: `forecasts` (unit, series,
#'   horizon, point, lower, upper, se), `evaluation` (per series, pooled over
#'   the units: RMSE of the own-lags and of the causal model, improvement,
#'   causes), `evaluation_units` (the same per unit), `holdout_fit`,
#'   `fitted_causes`, `coefficients` (`"var"`, `"bvar"`; standardised
#'   units), `model`, `p`. Methods: `print`, `summary`, `plot`.
#' @seealso [predict.causres()], [panel_causality()], [forecast_compare()]
#' @export
#' @examples
#' pd <- sim_panel(n_units = 20, n_time = 30, seed = 1)
#' pg <- panel_causality(pd, id = "unit", time = "time")
#' fc <- predict(pg, h = 3)
#' fc
#' \donttest{
#' plot(fc)
#' predict(pg, h = 3, model = "bvar", n_sim = 100)
#' }
predict.causpanel <- function(object, h = 5L, p = NULL, holdout = 0.2, level = 0.95, refit = TRUE,
                              model = c("var", "bvar", "rf"), control = bayes_control(prior = "minnesota"),
                              rf_settings = list(num_trees = 300L, mtry = NULL, min_node_size = 5L), n_sim = 200L, ...) {
  model <- match.arg(model)
  if (model == "rf") .need_ranger("predict(model = \"rf\")")
  .check_number(h, "h", lower = 1, integer = TRUE); .check_number(level, "level", lower = 0, upper = 1, strict = TRUE)
  vars <- object$series; D0 <- as.data.frame(object$data)
  PD <- .panel_units(D0, object$id, object$time, vars = vars); U <- PD$units; N <- length(U); n_i <- PD$n
  settings <- list(control = control, rf = utils::modifyList(list(num_trees = 300L, mtry = NULL, min_node_size = 5L),
                                                             rf_settings, keep.null = TRUE))
  full <- .panel_links(object, vars); causes_of <- full$causes
  p_full <- if (is.null(p)) full$p else max(1L, as.integer(p))
  holdout <- max(0, min(holdout, 0.5)); n_train <- pmax(floor(n_i * (1 - holdout)), 1L)
  ev <- NULL; evu <- NULL; hold <- NULL; eval_links <- NULL; eval_source <- "none"
  if (holdout > 0 && any(n_train < n_i)) {
    if (refit && !is.null(object$args)) {
      uid <- as.character(D0[[object$id]]); pos <- stats::ave(seq_along(uid), uid, FUN = seq_along)
      Dtr <- D0[pos <= n_train[match(uid, names(U))], , drop = FALSE]
      r_tr <- tryCatch(.rerun_object(object, Dtr), error = function(e) NULL)
      if (is.null(r_tr)) { el <- full; eval_source <- "full sample (the analysis could not be repeated on the training periods)" }
      else { el <- .panel_links(r_tr, vars); eval_source <- "training periods (refit)" }
    } else { el <- full; eval_source <- "full sample (optimistic)" }
    eval_links <- el$causes; p_ev <- if (is.null(p)) el$p else p_full
    if (sum(pmax(n_train - p_ev, 0L)) <= 3 * p_ev * length(vars) + 5) stop("Panel too short for the hold-out evaluation.", call. = FALSE)
    S <- .pfc_std(U, vars, n_train)
    res <- lapply(vars, function(v) {
      eb <- .pfc_errors(model, S, U, v, v, n_train, p_ev, settings)
      ec <- if (length(eval_links[[v]])) .pfc_errors(model, S, U, v, c(v, eval_links[[v]]), n_train, p_ev, settings) else eb
      data.frame(series = v, eb[, c("unit", "time", "observed")], pred_own = eb$pred, pred_causal = ec$pred, stringsAsFactors = FALSE)
    })
    hold <- do.call(rbind, res); rownames(hold) <- NULL
    rm <- function(e) sqrt(mean(e^2))
    ev <- do.call(rbind, lapply(vars, function(v) {
      k <- hold[hold$series == v, ]
      data.frame(series = v, causes = paste(eval_links[[v]], collapse = ","), rmse_baseline = rm(k$observed - k$pred_own),
                 rmse_causal = rm(k$observed - k$pred_causal), mae_baseline = mean(abs(k$observed - k$pred_own)),
                 mae_causal = mean(abs(k$observed - k$pred_causal)), stringsAsFactors = FALSE)
    }))
    ev$improvement <- ifelse(ev$rmse_baseline > 0, 1 - ev$rmse_causal / ev$rmse_baseline, 0)
    evu <- stats::aggregate(cbind(se_own = (observed - pred_own)^2, se_causal = (observed - pred_causal)^2) ~ unit + series, data = hold, FUN = mean)
    evu$rmse_baseline <- sqrt(evu$se_own); evu$rmse_causal <- sqrt(evu$se_causal)
    evu$improvement <- ifelse(evu$rmse_baseline > 0, 1 - evu$rmse_causal / evu$rmse_baseline, 0); evu$se_own <- evu$se_causal <- NULL
  }
  # ---- final model on all periods with the links of `object`
  p <- p_full; K <- length(vars); S <- .pfc_std(U, vars); all_rows <- lapply(n_i, seq_len)
  fits <- stats::setNames(lapply(vars, function(v) {
    a <- .pfc_xy(S, v, c(v, causes_of[[v]]), all_rows, p); f <- .fc_fit_xy(model, a$X, a$y, c(v, causes_of[[v]]), p, settings)
    f$key <- paste(a$unit, a$t); f
  }), vars)
  keys <- Reduce(intersect, lapply(fits, `[[`, "key"))
  R <- do.call(cbind, lapply(fits, function(f) f$resid[match(keys, f$key)])); colnames(R) <- vars
  Sigma <- stats::cov(R); A <- NULL
  if (model != "rf") {
    A <- array(0, c(K, K, p), dimnames = list(vars, vars, NULL))
    for (v in vars) { cf <- fits[[v]]$coef[-1L]; for (j in seq_along(fits[[v]]$preds)) A[v, fits[[v]]$preds[j], ] <- cf[(j - 1L) * p + seq_len(p)] }
  }
  ok <- n_i >= p; a <- (1 - level) / 2
  hist <- array(NA_real_, c(N, p, K), dimnames = list(NULL, NULL, vars))
  for (i in which(ok)) hist[i, , ] <- as.matrix(S[[i]]$D[(n_i[i] - p + 1L):n_i[i], vars, drop = FALSE])
  step <- function(paths, tcur, v) {                    # regressor rows of series v for every path at time tcur
    pr <- fits[[v]]$preds
    matrix(cbind(1, do.call(cbind, lapply(pr, function(w) matrix(paths[, tcur - seq_len(p), w], nrow = nrow(paths))))), nrow = nrow(paths))
  }
  if (model == "var") {
    paths <- array(NA_real_, c(N, p + h, K), dimnames = list(NULL, NULL, vars)); paths[, seq_len(p), ] <- hist
    Fm <- matrix(0, K * p, K * p); for (l in seq_len(p)) Fm[seq_len(K), (l - 1) * K + seq_len(K)] <- A[, , l]
    if (p > 1) Fm[K + seq_len(K * (p - 1)), seq_len(K * (p - 1))] <- diag(K * (p - 1))
    J <- cbind(diag(K), matrix(0, K, K * (p - 1))); MSE <- matrix(0, K, K); Fi <- diag(K * p); se <- matrix(0, h, K, dimnames = list(NULL, vars))
    for (k in seq_len(h)) {
      for (v in vars) paths[, p + k, v] <- .fc_predict(model, fits[[v]], step(paths, p + k, v))
      Phi <- J %*% Fi %*% t(J); MSE <- MSE + Phi %*% Sigma %*% t(Phi); Fi <- Fi %*% Fm; se[k, ] <- sqrt(pmax(diag(MSE), 0))
    }
    z <- stats::qnorm(1 - a)
    fc <- do.call(rbind, lapply(which(ok), function(i) do.call(rbind, lapply(vars, function(v) {
      pt <- paths[i, p + seq_len(h), v] * S[[i]]$sd[[v]] + S[[i]]$mu[[v]]; s <- se[, v] * S[[i]]$sd[[v]]
      data.frame(unit = names(U)[i], series = v, horizon = seq_len(h), point = pt, lower = pt - z * s, upper = pt + z * s, se = s, stringsAsFactors = FALSE)
    }))))
  } else {
    .check_number(n_sim, "n_sim", lower = 10, integer = TRUE)
    ui <- rep(which(ok), each = n_sim); M <- length(ui)
    paths <- array(NA_real_, c(M, p + h, K), dimnames = list(NULL, NULL, vars)); paths[, seq_len(p), ] <- hist[ui, , , drop = FALSE]
    if (model == "bvar") {
      draws <- lapply(vars, function(v) fits[[v]]$post$draw(n_sim)[rep(seq_len(n_sim), times = sum(ok)), , drop = FALSE]); names(draws) <- vars
      Cr <- stats::cor(R); Cr[!is.finite(Cr)] <- 0; diag(Cr) <- 1
      Lc <- tryCatch(t(chol(Cr)), error = function(e) diag(K))
    } else Rc <- R[stats::complete.cases(R), , drop = FALSE]
    for (k in seq_len(h)) {
      shock <- if (model == "bvar") t(Lc %*% matrix(stats::rnorm(K * M), K, M)) else Rc[sample.int(nrow(Rc), M, TRUE), , drop = FALSE]
      new <- vapply(seq_len(K), function(j) {
        v <- vars[j]; X <- step(paths, p + k, v)
        if (model == "bvar") { B <- draws[[v]]; rowSums(X * B[, seq_len(ncol(X)), drop = FALSE]) + B[, ncol(B)] * shock[, j] }
        else .fc_predict("rf", fits[[v]], X) + shock[, j]
      }, numeric(M))
      paths[, p + k, ] <- new
    }
    fc <- do.call(rbind, lapply(which(ok), function(i) do.call(rbind, lapply(vars, function(v) {
      x <- paths[ui == i, p + seq_len(h), v, drop = FALSE]; x <- matrix(x, ncol = h) * S[[i]]$sd[[v]] + S[[i]]$mu[[v]]
      data.frame(unit = names(U)[i], series = v, horizon = seq_len(h), point = colMeans(x),
                 lower = apply(x, 2, stats::quantile, a, names = FALSE), upper = apply(x, 2, stats::quantile, 1 - a, names = FALSE),
                 se = apply(x, 2, stats::sd), stringsAsFactors = FALSE)
    }))))
  }
  rownames(fc) <- NULL
  structure(list(forecasts = fc, evaluation = ev, evaluation_units = evu, holdout_fit = hold, evaluation_links = eval_links,
                 evaluation_source = eval_source, fitted_causes = causes_of, coefficients = A, sigma = Sigma, model = model, p = p,
                 h = h, level = level, holdout = holdout, data = D0, id = object$id, time = object$time, series = vars,
                 n_units = N, method = object$settings$method), class = "causfcpanel")
}

#' @export
print.causfcpanel <- function(x, digits = 3, ...) {
  mdl <- c(var = "pooled restricted VAR", bvar = "pooled Bayesian VAR, Minnesota prior", rf = "pooled random forests")[[x$model]]
  cat(sprintf("Causality-informed panel forecasts (%s from the '%s' links; %d units, p = %d lags, h = %d)\n", mdl, x$method, x$n_units, x$p, x$h))
  if (is.null(x$evaluation)) cat("  no hold-out evaluation\n") else {
    cat(sprintf("  hold-out evaluation (last %.0f%% of the periods of every unit, rolling one-step forecasts; links selected on the %s):\n",
                100 * x$holdout, x$evaluation_source))
    E <- x$evaluation; EU <- x$evaluation_units
    better <- vapply(E$series, function(v) mean(EU$improvement[EU$series == v] > 0), 0)
    print(data.frame(series = E$series, causes = ifelse(E$causes == "", "-", E$causes), RMSE.own = round(E$rmse_baseline, digits),
                     RMSE.causal = round(E$rmse_causal, digits), improvement = sprintf("%+.1f%%", 100 * E$improvement),
                     units.better = ifelse(E$causes == "", "-", sprintf("%.0f%%", 100 * better))), row.names = FALSE, right = FALSE)
  }
  cat("  forecast model: links of the full panel; causes per series:",
      paste(sprintf("%s <- {%s}", names(x$fitted_causes), vapply(x$fitted_causes, function(v) paste(v, collapse = ","), "")), collapse = "; "), "\n")
  Fc <- x$forecasts; u1 <- utils::head(unique(Fc$unit), 3L)
  cat(sprintf("\n  forecasts: point [%d%% interval]; first units (all in $forecasts):\n", round(100 * x$level)))
  f <- Fc[Fc$unit %in% u1 & Fc$horizon <= min(x$h, 5L), ]
  f$forecast <- sprintf("%s [%s, %s]", formatC(f$point, digits = digits - 1, format = "f"), formatC(f$lower, digits = digits - 1, format = "f"),
                        formatC(f$upper, digits = digits - 1, format = "f"))
  W <- stats::reshape(f[, c("unit", "series", "horizon", "forecast")], idvar = c("unit", "horizon"), timevar = "series", direction = "wide")
  names(W) <- sub("^forecast\\.", "", names(W)); names(W)[names(W) == "horizon"] <- "h"
  print(W, row.names = FALSE, right = FALSE)
  invisible(x)
}

#' @export
summary.causfcpanel <- function(object, ...) {
  E <- object$evaluation
  structure(list(evaluation = E, mean_improvement = if (is.null(E)) NA_real_ else mean(E$improvement), p = object$p, h = object$h,
                 best = if (is.null(E)) NA_character_ else E$series[which.max(E$improvement)]), class = "summary.causfc")
}

#' @rdname predict.causpanel
#' @param x a `causfcpanel` object.
#' @param series the series whose forecasts are drawn (default: the series
#'   with the largest hold-out improvement).
#' @param units units to draw (default: the first six).
#' @param n_last number of observed values shown before the forecasts.
#' @param main optional title.
#' @export
plot.causfcpanel <- function(x, series = NULL, units = NULL, n_last = 30L, main = NULL, ...) {
  E <- x$evaluation; EU <- x$evaluation_units; Fc <- x$forecasts
  if (is.null(series)) series <- if (is.null(E)) x$series[1L] else E$series[which.max(E$improvement)]
  if (!series %in% x$series) stop("Unknown series: ", series, call. = FALSE)
  if (is.null(units)) units <- utils::head(unique(Fc$unit), 6L)
  units <- intersect(units, unique(Fc$unit)); if (!length(units)) stop("No such unit.", call. = FALSE)
  np <- length(units) + !is.null(EU)
  op <- graphics::par(mfrow = grDevices::n2mfrow(np), mar = c(3.6, 4, 2.8, 1), mgp = c(2.1, 0.7, 0)); on.exit(graphics::par(op))
  if (!is.null(EU)) {
    imp <- split(100 * EU$improvement, factor(EU$series, levels = x$series))
    graphics::boxplot(imp, col = "#CFE8D5", border = "grey30", ylab = "RMSE improvement by unit (%)", las = 1, outline = TRUE,
                      main = if (is.null(main)) "Hold-out gain of the causal model" else main, font.main = 2, cex.main = 0.95)
    graphics::abline(h = 0, lty = 2, col = "grey40")
    graphics::points(seq_along(imp), 100 * E$improvement[match(names(imp), E$series)], pch = 23, bg = "#1B7837", cex = 1.3)
    graphics::mtext("diamond: pooled over the units", side = 3, line = 0.1, cex = 0.62, col = "grey30")
  }
  D <- x$data; uid <- as.character(D[[x$id]]); H <- x$holdout_fit
  for (u in units) {
    yv <- D[[series]][uid == u]; n <- length(yv); idx <- max(1, n - n_last + 1):n
    f <- Fc[Fc$unit == u & Fc$series == series, ]; hv <- if (is.null(H)) NULL else H[H$unit == u & H$series == series, ]
    graphics::plot(idx, yv[idx], type = "l", xlim = c(min(idx), n + x$h), ylim = range(c(yv[idx], f$lower, f$upper, hv$pred_causal, hv$pred_own), na.rm = TRUE),
                   xlab = "period", ylab = series, main = sprintf("%s: %s", u, series), font.main = 2, cex.main = 0.95, col = "grey30", lwd = 1.3, las = 1)
    if (!is.null(hv) && nrow(hv)) {
      graphics::lines(hv$time, hv$pred_own, col = "#B2403A", lty = 2, lwd = 1.2); graphics::lines(hv$time, hv$pred_causal, col = "#2B5C8A", lwd = 1.5)
    }
    graphics::polygon(c(n + f$horizon, rev(n + f$horizon)), c(f$lower, rev(f$upper)), col = grDevices::adjustcolor("#1B7A43", 0.2), border = NA)
    graphics::lines(c(n, n + f$horizon), c(yv[n], f$point), col = "#1B7A43", lwd = 2); graphics::abline(v = n, lty = 3, col = "grey60")
    if (u == units[1L]) graphics::legend("topleft", bty = "n", cex = 0.6, col = c("grey30", "#B2403A", "#2B5C8A", "#1B7A43"), lty = c(1, 2, 1, 1),
                                         lwd = c(1.3, 1.2, 1.5, 2), legend = c("observed", "own lags (one step)", "with causes (one step)", "forecast"))
  }
  invisible(x)
}

# forecast_compare() for panel data: the links of every method from the panel
# analysis of the training periods; pooled forecast equations; Diebold-Mariano
# test on the loss differentials averaged over the units per hold-out step
.forecast_compare_panel <- function(x, methods, model, holdout, p, args, cores, dots, ra) {
  id <- dots$id; time <- dots$time
  PD <- .panel_units(as.data.frame(x), id, time); U <- PD$units; vars <- PD$vars; D0 <- PD$data; n_i <- PD$n
  n_train <- pmax(floor(n_i * (1 - holdout)), 1L)
  uid <- as.character(D0[[id]]); pos <- stats::ave(seq_along(uid), uid, FUN = seq_along)
  Dtr <- D0[pos <= n_train[match(uid, names(U))], c(id, time, vars), drop = FALSE]
  settings <- list(control = if (is.null(dots$control)) bayes_control(prior = "minnesota") else dots$control,
                   rf = utils::modifyList(list(num_trees = 300L, mtry = NULL, min_node_size = 5L),
                                          if (is.null(dots$rf_settings)) list() else dots$rf_settings, keep.null = TRUE))
  sc <- .split_cores(length(methods), cores)
  runs <- .plapply(methods, function(m) {
    a <- c(list(x = Dtr, method = m, cores = sc$inner), ra, if (!is.null(args[[m]])) args[[m]]); a <- a[!duplicated(names(a), fromLast = TRUE)]
    r <- tryCatch(suppressWarnings(suppressMessages(do.call(rcausality, a))), error = function(e) NULL)
    if (is.null(r)) NULL else .panel_links(r, vars)
  }, cores = sc$outer)
  names(runs) <- methods
  bad <- vapply(runs, is.null, TRUE)
  if (all(bad)) stop("No method could be run on the training periods of the panel.", call. = FALSE)
  if (any(bad)) { warning("Method(s) left out: ", paste(methods[bad], collapse = ", "), call. = FALSE); runs <- runs[!bad]; methods <- methods[!bad] }
  p_all <- if (is.null(p)) max(vapply(runs, `[[`, 0L, "p")) else as.integer(p)
  S <- .pfc_std(U, vars, n_train)
  err_of <- function(v, preds) { e <- .pfc_errors(model, S, U, v, preds, n_train, p_all, settings); e$err <- e$observed - e$pred; e }
  # the test on the cross-sectional mean of the loss differential in every hold-out step
  dm_of <- function(eb, ec) {
    d <- tapply(eb$err^2 - ec$err^2, eb$step, mean); n <- length(d)
    if (n < 3 || stats::sd(d) == 0) return(c(statistic = NA_real_, p_value = NA_real_))
    s <- mean(d) / (stats::sd(d) / sqrt(n)); c(statistic = s, p_value = 2 * stats::pt(-abs(s), n - 1))
  }
  base <- stats::setNames(lapply(vars, function(v) err_of(v, v)), vars); rm <- function(e) sqrt(mean(e$err^2))
  rows <- lapply(vars, function(v) data.frame(method = "none", series = v, causes = "-", rmse = rm(base[[v]]), gain = 0, dm = NA_real_,
                                              p_dm = NA_real_, stringsAsFactors = FALSE))
  for (m in methods) for (v in vars) {
    cs <- runs[[m]]$causes[[v]]; e <- if (length(cs)) err_of(v, c(v, cs)) else base[[v]]
    dm <- if (length(cs)) dm_of(base[[v]], e) else c(statistic = NA_real_, p_value = NA_real_)
    rows[[length(rows) + 1L]] <- data.frame(method = m, series = v, causes = if (length(cs)) paste(cs, collapse = ",") else "-", rmse = rm(e),
                                           gain = 1 - rm(e) / rm(base[[v]]), dm = dm[["statistic"]], p_dm = dm[["p_value"]], stringsAsFactors = FALSE)
  }
  tab <- do.call(rbind, rows); rownames(tab) <- NULL
  sm <- stats::aggregate(gain ~ method, data = tab, FUN = mean); sm <- sm[match(c("none", methods), sm$method), ]; rownames(sm) <- NULL
  structure(list(table = tab, summary = sm, holdout = holdout, model = model, p = p_all, n = sum(n_i), n_units = length(U), methods = methods),
            class = "causfcomp")
}
