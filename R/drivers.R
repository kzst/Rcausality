# Explanatory causal analysis of one explained variable -------------------------
# Which series drive the target, at which lag, how strongly, directly or
# through another driver -- and how the explanatory variables drive each
# other.  Any method of rcausality(), its Bayesian variant, time series and
# panels.  The analysis generalises the Granger-based explanatory analysis
# of Kosztyan, Matrai and Kiss (2026) to all methods of the package.

.family_of <- function(cause, effect, target)
  ifelse(effect == target, "to_target", ifelse(cause == target, "feedback", "among"))

# family-wise multiplicity correction of a pairs-like table computed with
# adjust = "none": the p-values are corrected within the families (drivers ->
# target, target -> drivers, among the drivers); Bayes-factor, posterior and
# verdict decisions are kept
.family_adjust <- function(P, fam, adjust, alpha, use_bf) {
  P$p_adj <- NA_real_
  for (f in unique(fam)) {
    k <- fam == f
    if (any(!is.na(P$p_value[k]))) P$p_adj[k] <- if (adjust == "none") P$p_value[k] else stats::p.adjust(P$p_value[k], method = adjust)
  }
  if (!use_bf) P$significant <- P$significant %in% TRUE & (is.na(P$p_adj) | P$p_adj < alpha)
  P
}

# The delay at which a cause acts most strongly (Granger family).  The order
# of a Granger regression says how many lags are needed, not where the effect
# sits: an effect at lag 2 of an autocorrelated cause is already picked up by
# the order-1 test, because the first lag is a proxy of the second.  The peak
# lag is found by adding the cause's lags one at a time to the regression of
# the effect on its own (and the controls') past and taking the lag with the
# largest reduction of the residual sum of squares.  Panels: the reductions
# are added over the units, every unit with its own nuisance coefficients.
.peak_lag_unit <- function(e, x, Z = NULL, max_lag = 8L) {
  n <- length(e); nz <- if (is.null(Z)) 0L else ncol(Z)
  p <- max(1L, min(as.integer(max_lag), floor((n - 3) / (3 + nz))))
  Y <- e[(p + 1L):n]; X0 <- cbind(1, .lagmat(e, p))
  if (nz) for (j in seq_len(nz)) X0 <- cbind(X0, .lagmat(Z[, j], p))
  q <- qr(X0); yt <- qr.resid(q, Y); Xt <- qr.resid(q, .lagmat(x, p))
  list(num = as.numeric(crossprod(Xt, yt))^2, den = colSums(Xt^2), rss = sum(yt^2))
}
.peak_lag <- function(units, ef, cs, ctrl = character(0), max_lag = 8L) {
  n_min <- min(vapply(units, nrow, 0L)); nz <- length(ctrl)
  p <- max(1L, min(as.integer(max_lag), floor((n_min - 3) / (3 + nz))))
  gain <- numeric(p)
  for (u in units) {
    r <- .peak_lag_unit(.std(u[[ef]]), .std(u[[cs]]), if (nz) apply(as.matrix(u[, ctrl, drop = FALSE]), 2, .std) else NULL, p)
    g <- r$num / pmax(r$den, 1e-12); g[!is.finite(g)] <- 0
    gain <- gain + g
  }
  which.max(gain)
}

.drivers_run <- function(x, y, A, cl, lag_default = FALSE) {
  target <- A$target; panel <- !is.null(A$id)
  if (!is.character(target) || length(target) != 1L || is.na(target)) stop("`target` must be the name of one variable.", call. = FALSE)
  base <- if (A$method == "bayes_granger") "granger" else A$method; info <- .method_info(base)
  if (lag_default && !panel && base %in% c("granger", "toda_yamamoto")) A$lag <- "first"
  # ---- the data: all series in one table
  if (panel) {
    if (!is.null(y)) stop("Panel data: give one long data frame with the unit column `id`.", call. = FALSE)
    D <- as.data.frame(x)
    vars <- setdiff(names(D)[vapply(D, function(v) is.numeric(v) || is.logical(v), TRUE)], c(A$id, A$time))
  } else {
    X <- .as_series_frame(x, y, "x", "y")
    D <- if (is.null(X$Y)) X$X else {
      if (ncol(X$Y) != 1L) stop("With `target`, `y` must be one series (the target) or NULL.", call. = FALSE)
      if (!target %in% names(X$Y)) names(X$Y) <- target
      cbind(X$X[, setdiff(names(X$X), target), drop = FALSE], X$Y)
    }
    vars <- names(D)
  }
  if (!target %in% vars) stop(sprintf("Target '%s' is not a numeric column of the data.", target), call. = FALSE)
  drivers <- setdiff(vars, target)
  if (!length(drivers)) stop("No explanatory variable besides the target.", call. = FALSE)
  A1 <- A; A1["target"] <- list(NULL); A1$adjust <- "none"; A1$directions <- "both"
  fam_fun <- function(cause, effect) .family_of(cause, effect, target)
  old_rng <- .rng_setup(A$seed); on.exit(.rng_restore(old_rng), add = TRUE)
  A1["seed"] <- list(NULL)

  # ---- (1) all directed links among the series (the target included)
  if (panel) {
    r <- .panel_run(D, NULL, A1, cl, fam_fun = fam_fun, adjust = A$adjust)
    use_bf <- r$settings$panel$use_bf
  } else {
    r <- .rcausality_ts(D, NULL, A1, cl)
    use_bf <- isTRUE(r$settings$bayes) || r$settings$method %in% c("bayes_granger", "dc", "dfe")
    P <- .family_adjust(r$pairs, fam_fun(r$pairs$cause, r$pairs$effect), A$adjust, A$alpha, use_bf)
    r <- .rebuild_causres(r, P[, .pairs_cols])
  }
  r$settings$adjust <- A$adjust
  P <- r$pairs; to_t <- P[P$effect == target, , drop = FALSE]; to_t <- to_t[match(drivers, to_t$cause), , drop = FALSE]
  fb <- P[P$cause == target, , drop = FALSE]; fb <- fb[match(drivers, fb$effect), , drop = FALSE]
  among <- P[P$cause != target & P$effect != target, , drop = FALSE]

  # ---- (2) direct or mediated: the drivers of the target given all other series
  direct <- rep(NA, length(drivers)); p_dir <- bf_dir <- rep(NA_real_, length(drivers)); cond_note <- NULL
  light <- isTRUE(getOption("Rcausality.drivers_light"))        # re-runs of placebo_test() / robustness(): the links only
  if (isTRUE(A$conditional)) {
    direct <- to_t$significant; p_dir <- to_t$p_adj; bf_dir <- to_t$bf10
  } else if (light) {
  } else if (info$conditional && length(drivers) >= 2L) {
    A2 <- A1; A2$conditional <- TRUE; A2$directions <- "forward"
    rc <- tryCatch(suppressWarnings(suppressMessages({
      if (panel) .panel_run(D, NULL, A2, cl, xvars = drivers, yvars = target, fam_fun = fam_fun, adjust = A$adjust)
      else {
        z <- .rcausality_ts(D[, drivers, drop = FALSE], D[, target, drop = FALSE], A2, cl)
        Pz <- .family_adjust(z$pairs, fam_fun(z$pairs$cause, z$pairs$effect), A$adjust, A$alpha, use_bf); z$pairs <- Pz; z
      }})), error = function(e) conditionMessage(e))
    if (is.character(rc)) cond_note <- sprintf("the conditional tests could not be computed (%s)", rc)
    else {
      C <- rc$pairs[rc$pairs$effect == target, , drop = FALSE]; C <- C[match(drivers, C$cause), , drop = FALSE]
      direct <- C$significant %in% TRUE; p_dir <- C$p_adj; bf_dir <- C$bf10
      r$conditional <- rc$pairs
    }
  } else if (!info$conditional) cond_note <- sprintf("method \"%s\" cannot condition on the other series: direct and mediated links are not separated", base)

  # ---- (3) the drivers table
  # Granger family: the delay of the strongest effect, next to the order of the test
  peak <- function(cs, ef) NA_integer_
  if (base %in% c("granger", "toda_yamamoto") && !light) {
    UU <- if (panel) .panel_units(D, A$id, A$time, vars = vars, impute = A$impute)$units else list(D)
    peak <- function(cs, ef) tryCatch(.peak_lag(UU, ef, cs, if (isTRUE(A$conditional)) setdiff(vars, c(cs, ef)) else character(0), A$max_lag),
                                      error = function(e) NA_integer_)
  }
  sig <- to_t$significant %in% TRUE
  ev <- if (use_bf) pmax(0, log10(pmax(to_t$bf10, 1e-300))) else -log10(pmax(to_t$p_value, 1e-300))
  if (all(is.na(ev))) ev <- ifelse(sig, abs(to_t$statistic), 0)
  ev[!is.finite(ev)] <- 0
  imp <- if (sum(ev) > 0) ev / sum(ev) else rep(0, length(ev))
  edge <- function(a, b) any(among$cause == a & among$effect == b & among$significant %in% TRUE)
  via <- vapply(seq_along(drivers), function(i) {
    if (!sig[i]) return("")
    cand <- setdiff(drivers[sig & (is.na(direct) | direct %in% TRUE)], drivers[i])
    m <- cand[vapply(cand, function(mm) edge(drivers[i], mm), TRUE)]
    paste(m, collapse = ", ")
  }, "")
  role <- ifelse(!sig, "not a driver",
          ifelse(is.na(direct), ifelse(nzchar(via), "driver (possibly via another driver)", "driver"),
          ifelse(direct, "direct driver", ifelse(nzchar(via), "indirect driver", "driver, not direct"))))
  pk <- vapply(seq_along(drivers), function(i) if (sig[i]) as.integer(peak(drivers[i], target)) else NA_integer_, 0L)
  dt <- data.frame(driver = drivers, lag = ifelse(is.na(pk), to_t$lag, pk), order = to_t$lag, statistic = to_t$statistic, p_value = to_t$p_value, p_adj = to_t$p_adj,
                   bf10 = to_t$bf10, post_prob = to_t$post_prob, significant = sig, importance = imp, effect_size = to_t$effect_size,
                   sign = to_t$sign, feedback = fb$significant %in% TRUE, direct = direct, p_direct = p_dir, bf_direct = bf_dir,
                   via = ifelse(sig & !(direct %in% TRUE), via, ""), role = role, stringsAsFactors = FALSE)
  dt <- dt[order(!dt$significant, -dt$importance), , drop = FALSE]; rownames(dt) <- NULL
  r$drivers <- dt
  r$graph <- among[among$significant %in% TRUE, c("cause", "effect", "lag", "statistic", "p_adj", "bf10", "sign"), drop = FALSE]
  if (nrow(r$graph)) {
    pg <- vapply(seq_len(nrow(r$graph)), function(i) as.integer(peak(r$graph$cause[i], r$graph$effect[i])), 0L)
    r$graph$order <- r$graph$lag; r$graph$lag <- ifelse(is.na(pg), r$graph$lag, pg)
  }
  rownames(r$graph) <- NULL
  r$target <- target; r$settings$target <- target; r$settings$use_bf <- use_bf; r$settings$drivers_note <- cond_note
  A[["prior"]] <- r$args[["prior"]]; A$method <- r$args$method; A$bayes <- r$args$bayes
  if (!is.null(r$args$control)) A$control <- r$args$control
  if (panel) A$panel <- r$args$panel
  r$args <- .flat_args(A)
  if (!panel) r$args$x <- NULL
  class(r) <- c("causdrivers", class(r))
  r
}

#' Causal drivers of one explained variable
#'
#' An explanatory causal analysis organised around one target: which of the
#' other series drive it, at which lag, with how much evidence, directly or
#' through another driver, whether the target feeds back on them -- and which
#' causal links exist among the explanatory variables themselves. It works
#' with every method of [rcausality()], with the Bayesian variants, and with
#' time series as well as panels (`id`).
#'
#' Three families of directed links are tested and corrected for
#' multiplicity separately (`adjust`): drivers -> target, target -> drivers
#' (feedback) and the links among the drivers (the explanatory graph).
#' For the methods that can condition on other series, every driver is
#' tested again given all other series: a driver that stays significant is a
#' **direct** driver; one that does not is **indirect** when it drives a
#' direct driver (the mediator is named in `via`), and otherwise shares its
#' information with the other drivers. The **importance** of a driver is its
#' share in the evidence for the links into the target (`-log10 p`, or
#' `log10` Bayes factor for the Bayesian variants).
#'
#' For the Granger family the time-series test is by default that of the
#' first significant lag order (`lag = "first"`: the orders `1..max_lag` are
#' scanned, with a Bonferroni correction over the scan), as in the lag
#' ladder of Kosztyan, Matrai and Kiss (2026); in panels the order of every
#' link is chosen by the BIC of its pooled regression. The order of a Granger
#' regression says how many lags are needed, not where the effect sits (the
#' first lag of an autocorrelated cause is a proxy of its second lag), so
#' the table also gives the **lag** at which the driver acts most strongly:
#' the cause's lags are added one at a time to the regression of the target
#' on its own past, and the lag with the largest gain is reported. The
#' model-free methods report the lag with the largest statistic.
#'
#' @param data a data frame (or matrix) with the series in columns; panel
#'   data in long format with `id`.
#' @param target name of the explained variable.
#' @param drivers optional names of the explanatory variables (default: all
#'   other numeric columns).
#' @param method,bayes,cores the causality method, its Bayesian variant and
#'   the number of CPU cores, as in [rcausality()].
#' @param id,time panel data: names of the unit and time columns (see
#'   [panel_causality()]).
#' @param lag lag order; `NULL` = `"first"` for the Granger family on time
#'   series and `"auto"` otherwise (see [rcausality()]).
#' @param ... further arguments of [rcausality()] (`max_lag`, `alpha`,
#'   `adjust`, `conditional`, `control`, `prior`, `n_perm`, `seed`, `panel`,
#'   `method_args`, ...).
#' @return An object of class `causdrivers`, inheriting from
#'   [`causres`][causres-class] (and from `causpanel` for panel data), with
#'   * `drivers`: one row per explanatory variable -- `lag` (the delay of
#'     the strongest effect), `order` (the lag order of the test),
#'     `statistic`, `p_value`, `p_adj`, `bf10`, `significant`, `importance`,
#'     `effect_size`, `sign`, `feedback` (the target drives it), `direct`,
#'     `p_direct` / `bf_direct` (the conditional test), `via` (mediators)
#'     and `role`;
#'   * `graph`: the significant links among the explanatory variables;
#'   * `pairs`, `classes`, ...: all directed tests, as in a `causres`;
#'   * `conditional`: the conditional tests of the drivers.
#'
#'   Methods: `print`, `summary`, `plot` (see [plot.causdrivers()]) and the
#'   methods of `causres` ([placebo_test()], [robustness()], `predict`).
#' @references
#' Kosztyan, Z. T., Matrai, R. and Kiss, D. (2026). Is the world getting smaller? Network-driven heterogeneity in gravity-model elasticities. *The World Economy*. \doi{10.1111/twec.70148}
#'
#' Granger, C. W. J. (1969). Investigating causal relations by econometric models and cross-spectral methods. *Econometrica*, 37(3), 424-438. \doi{10.2307/1912791}
#'
#' Geweke, J. (1984). Measures of conditional linear dependence and feedback between time series. *Journal of the American Statistical Association*, 79(388), 907-915. \doi{10.1080/01621459.1984.10477110}
#' @seealso [rcausality()], [plot.causdrivers()], [sim_drivers()],
#'   [panel_causality()]
#' @export
#' @examples
#' # x1 and x3 drive y (at lags 1 and 2), x2 drives y only through x1, x4 is noise
#' d <- sim_drivers(300, seed = 1)
#' cd <- causal_drivers(d, target = "y")
#' cd
#' cd$drivers[, c("driver", "lag", "significant", "importance", "direct", "via")]
#' \donttest{
#' plot(cd)                  # the lag ladder with the explanatory graph
#' plot(cd, "importance")
#' # the same analysis with another method and its Bayesian variant
#' causal_drivers(d[, c("x1", "x3", "y")], target = "y", method = "te", bayes = TRUE,
#'                n_perm = 100, seed = 1)
#' # panel data: 15 units observed for 30 periods
#' pd <- sim_panel(n_units = 15, n_time = 30, scenario = "drivers", seed = 1)
#' cp <- causal_drivers(pd, target = "y", drivers = c("x1", "x2"), id = "unit", time = "time",
#'                      seed = 1)
#' cp
#' plot(cp, "heatmap")       # drivers by unit, with the panel-level tests
#' }
causal_drivers <- function(data, target, drivers = NULL, method = "granger", bayes = FALSE,
                           cores = getOption("Rcausality.cores", 1L), id = NULL, time = NULL, lag = NULL, ...) {
  cl <- match.call()
  if (!is.data.frame(data)) data <- as.data.frame(data)
  if (!is.null(drivers)) {
    bad <- setdiff(c(drivers, target), names(data))
    if (length(bad)) stop("Variable(s) not found: ", paste(bad, collapse = ", "), call. = FALSE)
    data <- data[, unique(c(id, time, drivers, target)), drop = FALSE]
  }
  out <- rcausality(data, method = method, bayes = bayes, cores = cores, lag = lag, target = target, id = id, time = time, ...)
  out$call <- cl
  out
}

#' Simulate a system with direct and indirect drivers of a target
#'
#' Five series: `x1 -> y` at lag 1 and `x3 -> y` at lag 2 are direct links,
#' `x2 -> x1` at lag 1 makes `x2` an indirect driver of `y` (through `x1`,
#' with a total delay of two steps), and `x4` is unrelated noise. Every
#' series has an autoregressive part.
#'
#' @param n number of observations.
#' @param coupling strength of the links.
#' @param nonlinear logical; couple through the centred square of the cause
#'   (invisible to linear tests) instead of linearly.
#' @param ar autoregressive coefficient.
#' @param seed optional RNG seed.
#' @return A data frame with the columns `x1`, `x2`, `x3`, `x4`, `y`.
#' @seealso [causal_drivers()], [sim_panel()] (`scenario = "drivers"`)
#' @export
#' @examples
#' head(sim_drivers(100, seed = 1))
sim_drivers <- function(n = 300L, coupling = 0.5, nonlinear = FALSE, ar = 0.5, seed = NULL) {
  .local_seed(seed)
  .check_number(n, "n", lower = 10, integer = TRUE)
  as.data.frame(.sim_drivers_unit(as.integer(n), coupling, nonlinear, ar))
}

.sim_drivers_unit <- function(n, b, nonlinear = FALSE, ar = 0.5, burn = 50L) {
  f <- if (nonlinear) function(v) (v^2 - 1) / sqrt(2) else identity
  Tt <- n + burn; E <- matrix(stats::rnorm(Tt * 5L), Tt, 5L)
  X <- matrix(0, Tt, 5L, dimnames = list(NULL, c("x1", "x2", "x3", "x4", "y")))
  for (t in 3:Tt) {
    X[t, "x2"] <- ar * X[t - 1, "x2"] + E[t, 2]
    X[t, "x1"] <- ar * X[t - 1, "x1"] + b * f(X[t - 1, "x2"]) + E[t, 1]
    X[t, "x3"] <- ar * X[t - 1, "x3"] + E[t, 3]
    X[t, "x4"] <- ar * X[t - 1, "x4"] + E[t, 4]
    X[t, "y"]  <- ar * X[t - 1, "y"] + b * f(X[t - 1, "x1"]) + b * f(X[t - 2, "x3"]) + E[t, 5]
  }
  X[-seq_len(burn), , drop = FALSE]
}

# ---- print / summary ------------------------------------------------------------

#' @export
print.causdrivers <- function(x, digits = 3, ...) {
  s <- x$settings; d <- x$drivers; panel <- inherits(x, "causpanel")
  base <- if (is.null(s$base_method)) s$method else s$base_method
  cat(sprintf("Causal drivers of '%s': %s%s\n", x$target, .method_info(base)$short, if (isTRUE(s$bayes)) ", Bayesian variant" else ""))
  if (panel) {
    ps <- s$panel
    cat(sprintf("  panel of %d units (%s periods), decision: %s\n", ps$n_units,
                if (ps$n_range[1] == ps$n_range[2]) ps$n_range[1] else paste(ps$n_range, collapse = "-"),
                c(group = "group level", pooled = "pooled level", both = "group or pooled level")[[ps$decision]]))
  } else cat(sprintf("  %d explanatory series, n = %d\n", nrow(d), s$n))
  cat(sprintf("  link if %s%s\n", if (isTRUE(s$use_bf)) sprintf("BF10 > %g", s$bf_threshold) else sprintf("adj. p < %g (%s)", s$alpha, s$adjust),
              if (isTRUE(s$use_bf)) "" else ";\n  correction within the families: drivers -> target, target -> drivers, among the drivers"))
  show <- data.frame(driver = d$driver, lag = ifelse(is.na(d$lag), "-", as.character(d$lag)), stringsAsFactors = FALSE)
  if (any(d$significant & !is.na(d$order) & d$order != d$lag)) show$order <- ifelse(is.na(d$order), "-", as.character(d$order))
  if (isTRUE(s$use_bf)) show$BF10 <- .fmt_bf(d$bf10) else show$p.adj <- .fmt_p(d$p_adj)
  show$importance <- sprintf("%4.0f%%", abs(round(100 * d$importance)))
  show$sign <- ifelse(is.na(d$sign), "", ifelse(d$sign >= 0, "+", "-"))
  if (any(!is.na(d$direct))) show$direct <- ifelse(!d$significant, "", ifelse(is.na(d$direct), "?", ifelse(d$direct, "yes", "no")))
  show$feedback <- ifelse(d$feedback, "yes", "")
  show$role <- ifelse(nzchar(d$via), sprintf("%s (via %s)", d$role, d$via), d$role)
  cat("\n"); print(show, row.names = FALSE, right = FALSE)
  g <- x$graph
  cat("\nLinks among the explanatory variables:")
  if (!nrow(g)) cat(" none\n") else {
    cat("\n"); for (i in seq_len(nrow(g))) cat(sprintf("  %s -> %s%s\n", g$cause[i], g$effect[i], if (is.na(g$lag[i])) "" else sprintf(" (lag %d)", g$lag[i])))
  }
  if (!is.null(s$drivers_note)) cat(sprintf("\nNote: %s.\n", s$drivers_note))
  cat("\n  lag: delay of the strongest effect", if (!is.null(show$order)) " (order: lag order of the test)" else "",
      ";\n  importance: share in the evidence for the links into the target;\n",
      "  direct: still a driver given all other series; feedback: the target drives the series.\n",
      "  plot(x) draws the lag ladder.\n", sep = "")
  invisible(x)
}

#' @export
summary.causdrivers <- function(object, ...) {
  structure(list(target = object$target, drivers = object$drivers, graph = object$graph, settings = object$settings,
                 panel = if (inherits(object, "causpanel")) .panel_table(object) else NULL), class = "summary.causdrivers")
}

#' @export
print.summary.causdrivers <- function(x, digits = 3, ...) {
  cat(sprintf("Causal drivers of '%s'\n\n", x$target))
  d <- x$drivers; num <- vapply(d, is.numeric, TRUE); d[num] <- lapply(d[num], function(v) signif(v, digits))
  print(d, row.names = FALSE)
  cat("\nLinks among the explanatory variables:\n")
  if (nrow(x$graph)) { g <- x$graph; num <- vapply(g, is.numeric, TRUE); g[num] <- lapply(g[num], function(v) signif(v, digits)); print(g, row.names = FALSE) }
  else cat("  none\n")
  if (!is.null(x$panel)) {
    cat("\nPanel levels of the links into the target:\n")
    p <- x$panel[grepl(paste0("-> ", x$target, "$"), x$panel$link), , drop = FALSE]
    num <- vapply(p, is.numeric, TRUE); p[num] <- lapply(p[num], function(v) signif(v, digits)); print(p, row.names = FALSE)
  }
  invisible(x)
}

# ---- plots -----------------------------------------------------------------------

#' Plot the causal drivers of a target
#'
#' @param x a `causdrivers` object.
#' @param type one of
#'   * `"ladder"`: the lag ladder. The target stands on the right; every
#'     driver is placed in the column of its lag (the further left, the
#'     longer the delay), from top to bottom by importance. Solid arrows are
#'     direct links, dashed arrows links that vanish given the other series
#'     (indirect or shared); the bar in a box is the importance of the
#'     driver; a returning arrow marks feedback from the target; the thin
#'     arrows between the boxes are the causal links among the explanatory
#'     variables. Series that do not drive the target are listed below;
#'   * `"importance"`: importance of every explanatory variable with its
#'     evidence, role and lag;
#'   * `"graph"`: the causal network of all series (target included), see
#'     [plot.causres()];
#'   * `"heatmap"`, `"units"`, `"prevalence"` (panel data): the drivers by
#'     unit, see [plot.causpanel()];
#'   * any other type of [plot.causres()].
#' @param graph logical; draw the links among the explanatory variables in
#'   the ladder?
#' @param main optional title.
#' @param ... passed on to [plot.causres()] / [plot.causpanel()].
#' @return Invisibly, `x`.
#' @references
#' Kosztyan, Z. T., Matrai, R. and Kiss, D. (2026). Is the world getting smaller? Network-driven heterogeneity in gravity-model elasticities. *The World Economy*. \doi{10.1111/twec.70148}
#' @export
#' @examples
#' cd <- causal_drivers(sim_drivers(300, seed = 1), target = "y")
#' plot(cd)
#' \donttest{
#' plot(cd, "importance")
#' plot(cd, "graph")
#' }
plot.causdrivers <- function(x, type = c("ladder", "importance", "graph", "heatmap", "units", "prevalence", "network", "matrix",
                                         "bars", "lags"), graph = TRUE, main = NULL, ...) {
  type <- match.arg(type)
  if (type == "ladder") return(invisible(.plot_ladder(x, graph, main)))
  if (type == "importance") return(invisible(.plot_importance(x, main)))
  if (type == "graph") type <- "network"
  if (type %in% c("heatmap", "units", "prevalence")) {
    if (!inherits(x, "causpanel")) stop(sprintf("type = \"%s\" needs panel data.", type), call. = FALSE)
    dots <- list(...)
    if (is.null(dots$links)) dots$links <- .link_key(x$drivers$driver, x$target)
    y <- x; class(y) <- setdiff(class(x), "causdrivers")
    do.call(graphics::plot, c(list(y, type = type, main = if (is.null(main) && type == "heatmap") sprintf("Drivers of %s by unit and in the panel", x$target) else main), dots))
    return(invisible(x))
  }
  y <- x; class(y) <- setdiff(class(x), "causdrivers")
  plot(y, type = type, main = main, ...)
  invisible(x)
}

.box <- function(x, y, w, h, col, border = "#555555", lwd = 1.2) {
  r <- min(w, h) * 0.22; a <- seq(0, pi / 2, length.out = 8)
  cx <- c(x + w / 2 - r + r * cos(a), x - w / 2 + r - r * sin(a), x - w / 2 + r - r * cos(a), x + w / 2 - r + r * sin(a))
  cy <- c(y + h / 2 - r + r * sin(a), y + h / 2 - r + r * cos(a), y - h / 2 + r - r * sin(a), y - h / 2 + r - r * cos(a))
  graphics::polygon(cx, cy, col = col, border = border, lwd = lwd)
}

.plot_ladder <- function(x, graph = TRUE, main = NULL) {
  d <- x$drivers; s <- x$settings; tg <- x$target
  sg <- d[d$significant, , drop = FALSE]; ns <- d[!d$significant, , drop = FALSE]
  lags <- sort(unique(ifelse(is.na(sg$lag), 0L, sg$lag)), decreasing = TRUE)
  if (!length(lags)) lags <- 1L
  nc <- length(lags); nr <- max(1L, nrow(sg))
  op <- graphics::par(mar = c(0.6, 0.6, 3.4, 0.6), xpd = NA); on.exit(graphics::par(op))
  bw <- 0.74; bh <- 0.56
  ybot <- if (nrow(ns)) -1.15 else -0.2
  graphics::plot.new(); graphics::plot.window(c(0.35, nc + 1.75), c(ybot - 0.75, nr + 0.75))
  base <- if (is.null(s$base_method)) s$method else s$base_method
  graphics::title(main = if (is.null(main)) sprintf("Causal drivers of %s", tg) else main, font.main = 2, line = 2.1)
  graphics::mtext(sprintf("%s%s;  link if %s", .method_info(base)$short, if (isTRUE(s$bayes)) " (Bayesian)" else "",
                          if (isTRUE(s$use_bf)) sprintf("BF10 > %g", s$bf_threshold) else sprintf("adj. p < %g (%s)", s$alpha, s$adjust)),
                  side = 3, line = 0.8, cex = 0.72, col = "grey30")
  # lag columns
  for (j in seq_len(nc)) {
    graphics::rect(j - 0.47, 0.42, j + 0.47, nr + 0.58, col = if (j %% 2) "#F4F6F8" else "#FAFBFC", border = NA)
    graphics::text(j, nr + 0.72, if (lags[j] == 0L) "lag -" else sprintf("lag %d", lags[j]), cex = 0.8, font = 2, col = "grey25")
  }
  xt <- nc + 1.15; yt <- (nr + 1) / 2
  pos <- list()
  if (nrow(sg)) {
    px <- match(ifelse(is.na(sg$lag), 0L, sg$lag), lags); py <- nr + 1 - seq_len(nrow(sg))
    for (i in seq_len(nrow(sg))) pos[[sg$driver[i]]] <- c(px[i], py[i])
    # arrows to the target
    for (i in seq_len(nrow(sg))) {
      dirct <- sg$direct[i]; lwd <- 1.2 + 5 * sg$importance[i] / max(sg$importance, 1e-9)
      col <- if (is.na(dirct)) "#4C8C6B" else if (dirct) "#1B7837" else "#8AA6B8"
      y1 <- yt + (py[i] - yt) * 0.5 * min(1, 1.6 / nr * 2)
      graphics::arrows(px[i] + bw / 2, py[i], xt - 0.42, y1, length = 0.11, lwd = lwd, col = col, lty = if (isFALSE(dirct)) 2 else 1)
      if (sg$feedback[i]) graphics::arrows(xt - 0.42, y1 - 0.09, px[i] + bw / 2 + 0.03, py[i] - 0.12, length = 0.07, lwd = 1.1, col = "#D98000", lty = 3)
    }
    # links among the explanatory variables (over the arrows, behind the boxes)
    g <- x$graph
    if (graph && nrow(g)) for (i in seq_len(nrow(g))) {
      a <- pos[[g$cause[i]]]; b <- pos[[g$effect[i]]]
      if (is.null(a) || is.null(b)) next
      .curve_arrow(a[1], a[2], b[1], b[2], bw, bh, col = "#7B5EA7", lab = if (is.na(g$lag[i])) "" else sprintf("lag %d", g$lag[i]))
    }
    for (i in seq_len(nrow(sg))) {
      dirct <- sg$direct[i]
      .box(px[i], py[i], bw, bh, col = if (isFALSE(dirct)) "#DDE8EF" else "#CFE8D5")
      graphics::text(px[i], py[i] + 0.1, sg$driver[i], cex = 0.85, font = 2, col = .pal$text)
      graphics::rect(px[i] - bw * 0.38, py[i] - 0.2, px[i] + bw * 0.38, py[i] - 0.11, col = "white", border = "grey60", lwd = 0.6)
      graphics::rect(px[i] - bw * 0.38, py[i] - 0.2, px[i] - bw * 0.38 + bw * 0.76 * sg$importance[i], py[i] - 0.11, col = "#2E7D32", border = NA)
      lab <- sprintf("%s%s  %.0f%%", if (is.na(sg$sign[i])) "" else if (sg$sign[i] >= 0) "(+) " else "(-) ",
                     if (isTRUE(s$use_bf)) paste0("BF ", .fmt_bf(sg$bf10[i])) else paste0("p ", .fmt_p(sg$p_adj[i])), 100 * sg$importance[i])
      graphics::text(px[i], py[i] - 0.4, lab, cex = 0.58, col = "grey25")
    }
  } else graphics::text(nc / 2 + 0.5, yt, "no significant driver", col = "grey40", cex = 0.95)
  .box(xt, yt, 0.82, max(0.8, min(nr * 0.6, 1.6)), col = "#BFD7EA", lwd = 1.6)
  graphics::text(xt, yt, tg, cex = 1, font = 2, col = .pal$text)
  if (nrow(ns)) {
    graphics::text(0.55, -0.42, "not drivers:", adj = 0, cex = 0.72, col = "grey35")
    w <- (nc + 1.2) / max(nrow(ns), 4)
    for (i in seq_len(nrow(ns))) {
      xx <- 0.55 + (i - 0.5) * min(w, 0.95)
      .box(xx, -0.95, min(w, 0.95) * 0.86, 0.42, col = "#F1F1F1", border = "grey65", lwd = 0.8)
      graphics::text(xx, -0.95, ns$driver[i], cex = 0.72, col = "grey35")
    }
  }
  has_dir <- any(!is.na(sg$direct))
  graphics::legend(0.4, ybot - 0.3, horiz = FALSE, ncol = 3, bty = "n", cex = 0.62, xjust = 0, yjust = 1,
                   legend = c(if (has_dir) c("direct link", "not direct (given the other series)") else "link into the target",
                              "feedback from the target", "link among the explanatory variables", "bar: importance"),
                   col = c(if (has_dir) c("#1B7837", "#8AA6B8") else "#4C8C6B", "#D98000", "#7B5EA7", "#2E7D32"),
                   lty = c(if (has_dir) c(1, 2) else 1, 3, 1, 1), lwd = c(if (has_dir) c(2.5, 2) else 2.5, 1.5, 1.3, 4))
  x
}

# curved arrow between two boxes of the ladder
.curve_arrow <- function(x0, y0, x1, y1, bw, bh, col, lab = "") {
  if (x0 == x1) {                                   # same column: bow on the left side
    sx <- x0 - bw / 2; t <- seq(0, 1, length.out = 30)
    cx <- sx - 0.22 * sin(pi * t) * (1 + 0.08 * abs(y1 - y0)); cy <- y0 + (y1 - y0) * t
  } else {
    sgn <- if (x1 > x0) 1 else -1; t <- seq(0, 1, length.out = 30)
    xa <- x0 + sgn * bw / 2; xb <- x1 - sgn * bw / 2
    ya <- y0 + if (y1 >= y0) 0.12 else -0.12; yb <- y1 + if (y1 >= y0) -0.12 else 0.12
    if (y0 == y1) { ya <- y0 + 0.16; yb <- y1 + 0.16 }
    cx <- xa + (xb - xa) * t; cy <- ya + (yb - ya) * (3 * t^2 - 2 * t^3) + if (y0 == y1) 0.18 * sin(pi * t) else 0
  }
  n <- length(cx)
  graphics::lines(cx[-n], cy[-n], col = col, lwd = 1.4)
  graphics::arrows(cx[n - 1], cy[n - 1], cx[n], cy[n], length = 0.08, col = col, lwd = 1.4)
  if (nzchar(lab)) {
    m <- max(2L, round(n * 0.22)); w <- graphics::strwidth(lab, cex = 0.55) / 2 + 0.02; h <- graphics::strheight(lab, cex = 0.55) / 2 + 0.02
    graphics::rect(cx[m] - w, cy[m] - h, cx[m] + w, cy[m] + h, col = grDevices::adjustcolor("white", 0.9), border = NA)
    graphics::text(cx[m], cy[m], lab, cex = 0.55, col = col)
  }
}

.plot_importance <- function(x, main = NULL) {
  d <- x$drivers; s <- x$settings; d <- d[order(d$importance), , drop = FALSE]; n <- nrow(d)
  lab_w <- max(graphics::strwidth(d$driver, units = "inches", cex = 0.9)) + 0.3
  op <- graphics::par(mai = c(0.85, lab_w, 1.1, 2.3), xpd = NA); on.exit(graphics::par(op))
  col <- ifelse(!d$significant, "grey82", ifelse(d$direct %in% FALSE, "#8AA6B8", "#2E7D32"))
  xm <- max(d$importance, 0.05) * 1.05
  graphics::plot.new(); graphics::plot.window(c(0, xm), c(0.4, n + 0.6), xaxs = "i")
  graphics::axis(1, at = pretty(c(0, xm)), labels = paste0(100 * pretty(c(0, xm)), "%"), cex.axis = 0.8)
  graphics::title(xlab = "importance (share in the evidence for the links into the target)", line = 2.2, cex.lab = 0.8)
  graphics::rect(0, seq_len(n) - 0.32, d$importance, seq_len(n) + 0.32, col = col, border = NA)
  graphics::text(-0.01 * xm, seq_len(n), d$driver, adj = 1, cex = 0.9, font = ifelse(d$significant, 2, 1))
  ev <- if (isTRUE(s$use_bf)) paste0("BF ", .fmt_bf(d$bf10)) else paste0("p ", .fmt_p(d$p_adj))
  lab <- sprintf("%s%s%s", ev, ifelse(is.na(d$lag) | !d$significant, "", sprintf(", lag %d", d$lag)),
                 ifelse(d$significant, sprintf(" - %s%s", d$role, ifelse(nzchar(d$via), sprintf(" (via %s)", d$via), "")), ""))
  graphics::text(xm * 1.02, seq_len(n), lab, adj = 0, cex = 0.66, col = "grey20")
  base <- if (is.null(s$base_method)) s$method else s$base_method
  graphics::title(main = if (is.null(main)) sprintf("Importance of the drivers of %s", x$target) else main, font.main = 2, line = 3)
  graphics::mtext(sprintf("%s%s", .method_info(base)$short, if (isTRUE(s$bayes)) " (Bayesian)" else ""), side = 3, line = 1.9, cex = 0.72, col = "grey30")
  graphics::legend(0, n + 0.62, horiz = TRUE, bty = "n", cex = 0.66, xjust = 0, yjust = 0, fill = c("#2E7D32", "#8AA6B8", "grey82"), border = NA,
                   legend = c("driver", "driver, not direct", "not a driver"), x.intersp = 0.5)
  x
}
