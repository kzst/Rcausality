# Granger causality in panels of (short) time series ---------------------------
# Many units observed over a short period (countries x years, firms x
# quarters): every unit alone has too few observations for a reliable test,
# but the panel carries the evidence.  panel_granger() triangulates
#  * unit-by-unit Granger tests with false-discovery-rate control across the
#    units (Benjamini-Hochberg);
#  * the Dumitrescu-Hurlin (2012) test for heterogeneous panels (average of
#    the unit Wald statistics, standardised with their exact finite-T
#    moments);
#  * the pooled fixed-effects Granger test (within transformation);
#  * the pooled Bayesian comparison (hierarchical Minnesota prior, or BIC /
#    g-prior), evidence for H0 as well as for H1;
#  * the pooled random-forest Granger-like test (out-of-bag gain against
#    within-unit circular-shift surrogates of the cause);
# and an optional placebo distribution of every pooled statistic from
# within-unit circular shifts of the cause (Kosztyan, Matrai and Kiss, 2026).

# split a long panel into a list of unit data frames ordered by time
.panel_units <- function(data, id, time, vars) {
  if (!is.data.frame(data)) stop("`data` must be a data frame in long format.", call. = FALSE)
  for (v in c(id, time, vars)) if (!v %in% names(data)) stop(sprintf("Column '%s' not found in `data`.", v), call. = FALSE)
  ord <- if (is.null(time)) order(data[[id]]) else order(data[[id]], data[[time]])
  D <- data[ord, , drop = FALSE]
  U <- split(D[, vars, drop = FALSE], as.character(D[[id]]))
  lapply(U, function(u) { u[] <- lapply(u, function(v) .prep_na(as.numeric(v))); u })
}

# stacked lag design of one directed pair: y, nuisance regressors (own and
# control lags), cause lags, unit index; within-unit demeaning removes the
# fixed effects.  `shift`: optional list of circular shifts of the cause per unit
.panel_design <- function(units, ef, cs, ctrl, L, shift = NULL, demean = TRUE) {
  parts <- lapply(seq_along(units), function(i) {
    u <- units[[i]]; n <- nrow(u)
    if (n < L + 3) return(NULL)
    x <- u[[cs]]; if (!is.null(shift)) x <- if (shift[i] > 0) c(x[(shift[i] + 1):n], x[seq_len(shift[i])]) else x
    Xs <- .lagmat(u[[ef]], L)
    for (c in ctrl) Xs <- cbind(Xs, .lagmat(u[[c]], L))
    list(y = u[[ef]][(L + 1):n], Xs = Xs, Xc = .lagmat(x, L), unit = rep(i, n - L))
  })
  parts <- parts[!vapply(parts, is.null, TRUE)]
  y <- unlist(lapply(parts, `[[`, "y")); Xs <- do.call(rbind, lapply(parts, `[[`, "Xs"))
  Xc <- do.call(rbind, lapply(parts, `[[`, "Xc")); unit <- unlist(lapply(parts, `[[`, "unit"))
  if (demean) {
    dm <- function(M) M - rowsum(M, unit)[match(unit, sort(unique(unit))), , drop = FALSE] / as.vector(table(unit)[as.character(unit)])
    y <- as.numeric(dm(matrix(y))); Xs <- dm(Xs); Xc <- dm(Xc)
  }
  list(y = y, Xs = Xs, Xc = Xc, unit = unit, n_units = length(unique(unit)))
}

.panel_shifts <- function(units, L) {
  vapply(units, function(u) { n <- nrow(u); m <- max(L + 1L, floor(0.1 * n)); if (n - m <= m) 0L else sample(m:(n - m), 1L) }, 0L)
}

# pooled statistics of one design: fixed-effects F-test, Bayes factor, RF gain
.panel_pooled_F <- function(ds) {
  N <- length(ds$y); q <- ncol(ds$Xc)
  rss_r <- sum(stats::lm.fit(ds$Xs, ds$y)$residuals^2)
  rss_u <- sum(stats::lm.fit(cbind(ds$Xs, ds$Xc), ds$y)$residuals^2)
  df2 <- N - ds$n_units - ncol(ds$Xs) - q
  Fs <- ((rss_r - rss_u) / q) / (rss_u / df2)
  list(statistic = Fs, df1 = q, df2 = df2, p_value = stats::pf(Fs, q, df2, lower.tail = FALSE),
       bf10_bic = exp(min((N * log(rss_r / rss_u) - q * log(N)) / 2, 700)), effect_size = log(rss_r / rss_u), n_eff = N)
}

.panel_bayes <- function(ds, L, control, prior_prob) {
  if (identical(control$prior, "minnesota")) {
    l_dec <- seq_len(L)^(-2 * control$decay); nb <- ncol(ds$Xs) / L
    ds_v <- c(l_dec, rep(control$cross^2 * l_dec, nb - 1L))
    cs <- function(M) { s <- apply(M, 2, stats::sd); s[!is.finite(s) | s <= 0] <- 1; sweep(M, 2, s, "/") }
    b <- .minnesota_core(ds$y / stats::sd(ds$y), matrix(0, length(ds$y), 0L), cs(ds$Xs), ds_v, cs(ds$Xc),
                         control$cross^2 * l_dec, control, prior_prob, draws = FALSE)
    return(b$log_bf10)
  }
  N <- length(ds$y)
  rss0 <- sum(ds$y^2); rss_r <- sum(stats::lm.fit(ds$Xs, ds$y)$residuals^2)
  rss_u <- sum(stats::lm.fit(cbind(ds$Xs, ds$Xc), ds$y)$residuals^2)
  if (control$prior %in% c("g", "hyper_g"))
    .log_bf_nested(N, ncol(ds$Xs), ncol(ds$Xs) + ncol(ds$Xc), 1 - rss_r / rss0, 1 - rss_u / rss0, control$prior)
  else (N * log(rss_r / rss_u) - ncol(ds$Xc) * log(N)) / 2
}

.panel_rf_gain <- function(ds, settings) {
  f_own <- .rf_fit(ds$Xs, ds$y, settings); f_full <- .rf_fit(cbind(ds$Xs, ds$Xc), ds$y, settings)
  1 - f_full$prediction.error / f_own$prediction.error
}

# Dumitrescu-Hurlin from the unit Wald statistics W_i = df1 * F_i; the exact
# moments of K * F(K, d2) under H0 give the finite-T standardisation Z-tilde
.dh_test <- function(W, K, d2) {
  ok <- is.finite(W) & d2 > 4
  W <- W[ok]; d2 <- d2[ok]; N <- length(W)
  if (N < 2L) return(c(W_bar = NA, Z_bar = NA, Z_tilde = NA, p_value = NA, n_units = N))
  EW <- K * d2 / (d2 - 2); VW <- 2 * K * d2^2 * (K + d2 - 2) / ((d2 - 2)^2 * (d2 - 4))
  Zb <- sqrt(N / (2 * K)) * (mean(W) - K)
  Zt <- (sum(W) - sum(EW)) / sqrt(sum(VW))
  c(W_bar = mean(W), Z_bar = Zb, Z_tilde = Zt, p_value = 2 * stats::pnorm(-abs(Zt)), n_units = N)
}

#' Granger causality in panels of short time series
#'
#' Tests every ordered pair of the variables of a long-format panel (units x
#' time) with several complementary procedures, so that the evidence can be
#' triangulated: the unit-by-unit tests are too weak when each unit has a
#' short series, the pooled tests assume a common effect, the
#' Dumitrescu-Hurlin test allows heterogeneous effects, and the Bayesian and
#' random-forest versions add evidence for the null and nonlinear
#' predictability.
#'
#' * `"units"`: Granger F-test in every unit (no pooling), Benjamini-Hochberg
#'   adjusted across the units of a pair; the statistic is the share of
#'   units with a significant link and the p-value the smallest adjusted
#'   unit p-value (Simes' test of the global null "no unit has a link");
#' * `"dh"`: Dumitrescu and Hurlin (2012): the mean of the unit Wald
#'   statistics `W_i = K F_i` standardised with the exact moments of `K F(K,
#'   d2)` under H0 (`Z-tilde`, valid for fixed T) and asymptotically
#'   (`Z-bar`); two-sided p-value of `Z-tilde`, as in \pkg{plm} and
#'   `xtgcause`;
#' * `"pooled"`: fixed-effects (within) regression of the effect on its own
#'   lags and the cause's lags pooled over the units; F-test with
#'   `N - n_units - k` residual degrees of freedom (the own-lag coefficients
#'   have the Nickell bias of order 1/T, the test of the cause's lags remains
#'   informative);
#' * `"bayes"`: the pooled comparison with a Bayes factor (`control$prior`:
#'   the hierarchical Minnesota prior by default, or `"g"`, `"hyper_g"`,
#'   `"bic"`);
#' * `"rf"`: pooled random-forest Granger-like test: out-of-bag gain of the
#'   forest with the cause's lags (within-unit demeaned data) against
#'   `n_surr` within-unit circular-shift surrogates of the cause.
#'
#' With `n_placebo > 0` every pooled statistic is recomputed on placebo
#' panels in which the cause is circularly shifted within every unit by an
#' independent random offset (the time-series structure of every unit kept,
#' the cross-variable timing destroyed), and the placebo p-value is reported.
#'
#' @param data a data frame in long format.
#' @param id name of the unit column.
#' @param time optional name of the time column (the rows are ordered by it
#'   within units; consecutive periods are assumed, missing values are
#'   interpolated).
#' @param vars the variables to analyse (default: every numeric column
#'   other than `id` and `time`).
#' @param lag lag order `K` (common to all units).
#' @param methods subset of `c("units", "dh", "pooled", "bayes", "rf")`.
#' @param conditional logical; condition every test on the lags of the
#'   other variables.
#' @param alpha significance level; `bf_threshold` Bayes-factor threshold.
#' @param bf_threshold Bayes factor above which `"bayes"` declares a link.
#' @param adjust correction across the directed pairs within a method
#'   (passed to [stats::p.adjust()], default `"BH"`).
#' @param control a [bayes_control()] list for `"bayes"` (and the prior of the
#'   unit-level Bayes factors).
#' @param n_surr surrogates of the `"rf"` test.
#' @param rf_settings forest settings of `"rf"` (`num_trees`, `mtry`,
#'   `min_node_size`), see [rf_granger()].
#' @param n_placebo number of placebo panels (`0` = none).
#' @param seed optional RNG seed.
#' @param cores CPU cores: pairs, surrogates and placebo panels run in
#'   parallel with results independent of the number of workers (see
#'   [Rcausality-parallel]).
#' @return An object of class `panel_granger`: `tests` (one row per directed
#'   pair and method: statistic, p-value, adjusted p-value, Bayes factor,
#'   placebo p-value, decision), `units` (the unit-level tests), `settings`.
#'   Methods: `print`, `plot`.
#' @references
#' Dumitrescu, E.-I. and Hurlin, C. (2012). Testing for Granger non-causality in heterogeneous panels. *Economic Modelling*, 29(4), 1450-1460. \doi{10.1016/j.econmod.2012.02.014}
#'
#' Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society B*, 57(1), 289-300. \doi{10.1111/j.2517-6161.1995.tb02031.x}
#'
#' Nickell, S. (1981). Biases in dynamic models with fixed effects. *Econometrica*, 49(6), 1417-1426. \doi{10.2307/1911408}
#'
#' Kosztyan, Z. T., Matrai, R. and Kiss, D. (2026). Is the world getting smaller? Network-driven heterogeneity in gravity-model elasticities. *The World Economy*. \doi{10.1111/twec.70148}
#' @export
#' @examples
#' # 30 units observed for 20 periods; x drives y in every unit with a
#' # unit-specific strength.  A single unit is too short for a reliable test,
#' # the panel procedures agree on x -> y.
#' pd <- sim_panel(n_units = 30, n_time = 20, scenario = "x_to_y", seed = 1)
#' pg <- panel_granger(pd, id = "unit", time = "time", methods = c("units", "dh", "pooled", "bayes"))
#' pg
#' \donttest{
#' plot(pg)
#' # with a placebo distribution of the pooled statistics
#' panel_granger(pd, id = "unit", time = "time", methods = c("dh", "pooled"), n_placebo = 99, seed = 1)
#' }
panel_granger <- function(data, id, time = NULL, vars = NULL, lag = 1L,
                          methods = c("units", "dh", "pooled", "bayes"), conditional = FALSE, alpha = 0.05,
                          bf_threshold = 3, adjust = "BH", control = bayes_control(prior = "minnesota"), n_surr = 99L,
                          rf_settings = list(num_trees = 200L, mtry = NULL, min_node_size = 5L), n_placebo = 0L,
                          seed = NULL, cores = getOption("Rcausality.cores", 1L)) {
  methods <- match.arg(methods, c("units", "dh", "pooled", "bayes", "rf"), several.ok = TRUE)
  .check_number(lag, "lag", lower = 1, integer = TRUE); L <- as.integer(lag)
  .check_number(n_placebo, "n_placebo", lower = 0, integer = TRUE)
  if ("rf" %in% methods) .need_ranger("panel_granger(methods = \"rf\")")
  if (!inherits(control, "bayes_control")) control <- do.call(bayes_control, as.list(control))
  .local_seed(seed)
  if (is.null(vars)) vars <- setdiff(names(data)[vapply(data, is.numeric, TRUE)], c(id, time))
  if (length(vars) < 2L) stop("At least two variables are needed.", call. = FALSE)
  units <- .panel_units(data, id, time, vars)
  T_i <- vapply(units, nrow, 0L)
  grid <- expand.grid(cause = vars, effect = vars, stringsAsFactors = FALSE)
  grid <- grid[grid$cause != grid$effect, , drop = FALSE]; rownames(grid) <- NULL
  rf_st <- utils::modifyList(list(num_trees = 200L, mtry = NULL, min_node_size = 5L), rf_settings, keep.null = TRUE)
  sc <- .split_cores(nrow(grid), cores)
  one_pair <- function(i) {
    cs <- grid$cause[i]; ef <- grid$effect[i]; ctrl <- if (conditional) setdiff(vars, c(cs, ef)) else character(0)
    out <- list(units = NULL, tests = list())
    if (any(c("units", "dh") %in% methods)) {
      ut <- do.call(rbind, lapply(names(units), function(nm) {
        u <- units[[nm]]
        g <- if (nrow(u) > (2 + length(ctrl)) * L + 3) .granger_pair(u[[ef]], u[[cs]], if (length(ctrl)) as.matrix(u[, ctrl, drop = FALSE]), L)
             else list(statistic = NA_real_, df1 = L, df2 = NA_real_, p_value = NA_real_, bf10 = NA_real_)
        data.frame(unit = nm, cause = cs, effect = ef, n = nrow(u), F = g$statistic, df1 = g$df1, df2 = g$df2,
                   W = g$df1 * g$statistic, p_value = g$p_value, bf10 = g$bf10, stringsAsFactors = FALSE)
      }))
      ut$p_adj <- stats::p.adjust(ut$p_value, method = "BH")
      ut$significant <- !is.na(ut$p_adj) & ut$p_adj < alpha
      out$units <- ut
      if ("units" %in% methods)
        out$tests$units <- list(statistic = mean(ut$significant), p_value = suppressWarnings(min(ut$p_adj, na.rm = TRUE)),
                                bf10 = NA_real_, detail = sprintf("%d of %d units significant (BH)", sum(ut$significant), sum(!is.na(ut$p_value))))
      if ("dh" %in% methods) {
        dh <- .dh_test(ut$W, L, ut$df2)
        out$tests$dh <- list(statistic = dh[["Z_tilde"]], p_value = dh[["p_value"]], bf10 = NA_real_,
                             detail = sprintf("W-bar = %.3f, Z-bar = %.2f, %d units", dh[["W_bar"]], dh[["Z_bar"]], as.integer(dh[["n_units"]])))
      }
    }
    ds <- .panel_design(units, ef, cs, ctrl, L)
    if ("pooled" %in% methods) {
      pf <- .panel_pooled_F(ds)
      out$tests$pooled <- list(statistic = pf$statistic, p_value = pf$p_value, bf10 = pf$bf10_bic,
                               detail = sprintf("F(%d, %d), %d obs. of %d units", pf$df1, pf$df2, pf$n_eff, ds$n_units))
    }
    if ("bayes" %in% methods) {
      lb <- .panel_bayes(ds, L, control, 0.5)
      out$tests$bayes <- list(statistic = lb, p_value = NA_real_, bf10 = exp(max(min(lb, 700), -700)),
                              detail = sprintf("log BF10 = %.2f (%s prior)", lb, control$prior))
    }
    if ("rf" %in% methods) {
      gain <- .panel_rf_gain(ds, rf_st)
      null <- unlist(.plapply(seq_len(n_surr), function(b)
        .panel_rf_gain(.panel_design(units, ef, cs, ctrl, L, shift = .panel_shifts(units, L)), rf_st), cores = sc$inner))
      out$tests$rf <- list(statistic = gain, p_value = if (gain > 0) (1 + sum(null >= gain)) / (1 + n_surr) else 1,
                           bf10 = NA_real_, detail = sprintf("OOB gain %+.3f vs %d within-unit shifts", gain, n_surr))
    }
    if (n_placebo > 0) {
      pm <- intersect(c("dh", "pooled", "bayes"), methods)
      if (length(pm)) {
        pl <- .plapply(seq_len(n_placebo), function(b) {
          sh <- .panel_shifts(units, L); res <- c()
          if ("dh" %in% pm) {
            W <- vapply(seq_along(units), function(k) {
              u <- units[[k]]; n <- nrow(u); s <- sh[k]; xs <- if (s > 0) c(u[[cs]][(s + 1):n], u[[cs]][seq_len(s)]) else u[[cs]]
              if (n <= (2 + length(ctrl)) * L + 3) return(NA_real_)
              g <- .granger_pair(u[[ef]], xs, if (length(ctrl)) as.matrix(u[, ctrl, drop = FALSE]), L); g$df1 * g$statistic
            }, 0)
            d2 <- vapply(units, function(u) nrow(u) - L - (2 + length(ctrl)) * L - 1, 0)
            res["dh"] <- .dh_test(W, L, d2)[["Z_tilde"]]
          }
          dsp <- .panel_design(units, ef, cs, ctrl, L, shift = sh)
          if ("pooled" %in% pm) res["pooled"] <- .panel_pooled_F(dsp)$statistic
          if ("bayes" %in% pm) res["bayes"] <- .panel_bayes(dsp, L, control, 0.5)
          res
        }, cores = sc$inner)
        PL <- do.call(rbind, pl)
        for (m in pm) {
          obs <- out$tests[[m]]$statistic
          out$tests[[m]]$p_placebo <- (1 + sum(PL[, m] >= obs, na.rm = TRUE)) / (1 + sum(!is.na(PL[, m])))
        }
      }
    }
    out
  }
  res <- .plapply(seq_len(nrow(grid)), one_pair, cores = sc$outer)
  tests <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) do.call(rbind, lapply(names(res[[i]]$tests), function(m) {
    t <- res[[i]]$tests[[m]]
    data.frame(cause = grid$cause[i], effect = grid$effect[i], method = m, statistic = t$statistic, p_value = t$p_value,
               bf10 = t$bf10, p_placebo = if (is.null(t$p_placebo)) NA_real_ else t$p_placebo, detail = t$detail,
               stringsAsFactors = FALSE)
  }))))
  tests$p_adj <- NA_real_
  for (m in unique(tests$method)) {
    k <- tests$method == m
    if (any(!is.na(tests$p_value[k]))) tests$p_adj[k] <- stats::p.adjust(tests$p_value[k], method = adjust)
  }
  tests$significant <- ifelse(tests$method == "bayes", tests$bf10 > bf_threshold, !is.na(tests$p_adj) & tests$p_adj < alpha)
  tests <- tests[, c("cause", "effect", "method", "statistic", "p_value", "p_adj", "bf10", "p_placebo", "significant", "detail")]
  units_tab <- do.call(rbind, lapply(res, `[[`, "units")); if (!is.null(units_tab)) rownames(units_tab) <- NULL
  structure(list(tests = tests, units = units_tab, vars = vars,
                 settings = list(lag = L, methods = methods, conditional = conditional, alpha = alpha, bf_threshold = bf_threshold,
                                 adjust = adjust, prior = control$prior, n_units = length(units), T_range = range(T_i),
                                 n_surr = n_surr, n_placebo = n_placebo)),
            class = "panel_granger")
}

#' @export
print.panel_granger <- function(x, digits = 3, ...) {
  s <- x$settings
  cat(sprintf("Panel Granger causality: %d units, T = %d..%d, lag %d, %s\n", s$n_units, s$T_range[1], s$T_range[2], s$lag,
              if (s$conditional) "conditional on the other variables" else "pairwise"))
  lab <- c(units = "unit tests (BH)", dh = "Dumitrescu-Hurlin Z~", pooled = "pooled FE F-test",
           bayes = "pooled Bayes factor", rf = "pooled random forest")
  T <- x$tests
  tb <- data.frame(link = paste(T$cause, "->", T$effect), method = lab[T$method],
                   statistic = formatC(T$statistic, digits = digits, format = "g"),
                   p.adj = ifelse(T$method == "bayes", "-", .fmt_p(T$p_adj)),
                   BF10 = ifelse(is.na(T$bf10), "-", .fmt_bf(T$bf10)), stringsAsFactors = FALSE)
  if (s$n_placebo > 0) tb$p.placebo <- ifelse(is.na(T$p_placebo), "-", .fmt_p(T$p_placebo))
  tb$sig <- ifelse(T$significant, "yes", "no"); tb$detail <- T$detail
  print(tb, row.names = FALSE, right = FALSE)
  cat(sprintf("  decisions: adjusted p < %g (%s across the directed pairs), Bayes factor > %g; 'units' = share of units\n",
              s$alpha, s$adjust, s$bf_threshold))
  cat("  with a significant unit-level test after Benjamini-Hochberg correction across the units\n")
  invisible(x)
}

#' @rdname panel_granger
#' @param x a `panel_granger` object.
#' @param main optional title.
#' @param ... unused.
#' @export
plot.panel_granger <- function(x, main = NULL, ...) {
  U <- x$units
  if (is.null(U)) stop("No unit-level tests (methods 'units' or 'dh') to plot.", call. = FALSE)
  links <- unique(paste(U$cause, "->", U$effect)); k <- length(links)
  lp <- -log10(pmax(U$p_value, 1e-12)); key <- paste(U$cause, "->", U$effect)
  op <- graphics::par(mar = c(4.2, 7, 3, 1)); on.exit(graphics::par(op))
  graphics::plot(NA, xlim = c(0, max(lp, 2, na.rm = TRUE) * 1.05), ylim = c(0.5, k + 0.5), yaxt = "n",
                 xlab = "-log10 p of the unit-level Granger test", ylab = "",
                 main = if (is.null(main)) sprintf("Unit-level evidence (%d units, lag %d)", x$settings$n_units, x$settings$lag) else main,
                 font.main = 2, cex.main = 0.95)
  graphics::axis(2, at = seq_len(k), labels = links, las = 1, cex.axis = 0.85)
  graphics::abline(v = -log10(x$settings$alpha), lty = 2, col = "grey50")
  for (j in seq_len(k)) {
    idx <- which(key == links[j]); yy <- j + stats::runif(length(idx), -0.25, 0.25)
    graphics::points(lp[idx], yy, pch = 21, bg = ifelse(U$significant[idx], .pal$causal, "grey85"), col = .pal$border, cex = 0.9)
    T <- x$tests[paste(x$tests$cause, "->", x$tests$effect) == links[j], ]
    txt <- paste(vapply(seq_len(nrow(T)), function(r) sprintf("%s:%s", T$method[r], if (T$significant[r]) "+" else "-"), ""), collapse = " ")
    graphics::text(graphics::par("usr")[2], j + 0.38, txt, adj = c(1, 0.5), cex = 0.65, col = "grey30")
  }
  graphics::legend("bottomright", legend = c("significant after BH (unit level)", "not significant", sprintf("p = %g", x$settings$alpha)),
                   pt.bg = c(.pal$causal, "grey85", NA), pch = c(21, 21, NA), lty = c(NA, NA, 2), col = c(.pal$border, .pal$border, "grey50"),
                   bty = "n", cex = 0.7)
  invisible(x)
}

#' Simulate a panel of short VAR(1) series
#'
#' Every unit follows `x[t] = a x[t-1] + e1`, `y[t] = a y[t-1] + b_i x[t-1] +
#' e2` (scenario `"x_to_y"`), with unit-specific coupling `b_i ~ N(coupling,
#' (heterogeneity * coupling)^2)` and unit fixed effects; `"independent"` sets
#' `b_i = 0`, `"bidirectional"` couples both ways.
#'
#' @param n_units number of units.
#' @param n_time number of periods per unit.
#' @param scenario `"x_to_y"`, `"bidirectional"` or `"independent"`.
#' @param coupling mean coupling strength; `heterogeneity` its relative
#'   standard deviation across the units.
#' @param heterogeneity relative standard deviation of the coupling.
#' @param ar autoregressive coefficient.
#' @param seed optional RNG seed.
#' @return A long data frame with columns `unit`, `time`, `x`, `y`.
#' @export
#' @examples
#' head(sim_panel(n_units = 3, n_time = 5, seed = 1))
sim_panel <- function(n_units = 30L, n_time = 20L, scenario = c("x_to_y", "bidirectional", "independent"),
                      coupling = 0.4, heterogeneity = 0.5, ar = 0.5, seed = NULL) {
  scenario <- match.arg(scenario); .local_seed(seed)
  out <- lapply(seq_len(n_units), function(i) {
    b <- stats::rnorm(1, coupling, heterogeneity * abs(coupling)); fe <- stats::rnorm(2)
    bxy <- if (scenario == "independent") 0 else b; byx <- if (scenario == "bidirectional") b else 0
    burn <- 50L; Tt <- n_time + burn; x <- y <- numeric(Tt)
    for (t in 2:Tt) {
      x[t] <- ar * x[t - 1] + byx * y[t - 1] + stats::rnorm(1)
      y[t] <- ar * y[t - 1] + bxy * x[t - 1] + stats::rnorm(1)
    }
    data.frame(unit = sprintf("u%02d", i), time = seq_len(n_time), x = x[-seq_len(burn)] + fe[1], y = y[-seq_len(burn)] + fe[2],
               stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, out); attr(d, "truth") <- scenario; d
}
