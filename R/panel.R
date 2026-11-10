# Causality in panel data ------------------------------------------------------
# Many units (countries, firms, sectors, subjects) observed over time, often
# over a short period.  Every method of rcausality() is carried to panels at
# three levels:
#   units  -- the method in every unit on its own, with false-discovery-rate
#             control over all unit x link cells (Benjamini-Hochberg; Bayesian
#             FDR for the Bayesian variants);
#   group  -- a test of "no unit has the link" built from the unit results:
#             the mean of the standardised unit statistics against the
#             convolution of their null distributions (bootstrap for the
#             Granger family, surrogates for the model-free methods), the
#             Dumitrescu-Hurlin (2012) statistic, Stouffer / Fisher / Simes
#             combinations of the unit p-values, or the count of units with
#             the link; Bayesian variant: the Bayes factor of the group
#             statistic against its null distribution, or a prevalence model
#             on the unit Bayes factors;
#   pooled -- one estimate from all units: the heterogeneous-panel Granger
#             test of Juodis, Karavias and Sarafidis (2021) with the half-panel
#             jackknife for the Granger family (R/panel_granger.R); for the
#             model-free methods the statistic of the stacked, within-unit
#             standardised delay vectors (random forest, Diks-Panchenko, kNN
#             conditional mutual information, transfer entropy, multispatial
#             convergent cross mapping) against within-unit circular shifts
#             of the cause; Bayesian variant: Bayes factor of the pooled
#             estimate or surrogate Bayes factor.
# The approach follows the triangulation of traditional, Bayesian and
# random-forest panel Granger tests with FDR control and placebo tests of
# Kosztyan, Matrai and Kiss (2026).

# ---- data handling ---------------------------------------------------------------

# split a long panel into a list of unit data frames ordered by time
.panel_units <- function(data, id, time = NULL, vars = NULL, impute = TRUE, time_effects = FALSE) {
  if (!is.data.frame(data)) stop("Panel data must be a data frame in long format (one row per unit and period).", call. = FALSE)
  data <- as.data.frame(data)
  for (v in c(id, time)) if (!is.character(v) || length(v) != 1L || !v %in% names(data))
    stop(sprintf("Column '%s' not found in the data.", paste(v, collapse = ", ")), call. = FALSE)
  if (is.null(vars)) vars <- setdiff(names(data)[vapply(data, function(v) is.numeric(v) || is.logical(v), TRUE)], c(id, time))
  bad <- setdiff(vars, names(data)); if (length(bad)) stop("Variable(s) not found: ", paste(bad, collapse = ", "), call. = FALSE)
  if (length(vars) < 2L) stop("At least two numeric variables are needed besides `id` and `time`.", call. = FALSE)
  uid <- as.character(data[[id]])
  if (anyNA(uid)) stop("Missing values in the unit column `id`.", call. = FALSE)
  if (!is.null(time)) {
    if (any(duplicated(data.frame(uid, data[[time]])))) stop("Duplicated (id, time) rows in the panel.", call. = FALSE)
    data <- data[order(uid, data[[time]]), , drop = FALSE]; uid <- as.character(data[[id]])
  }
  ids <- unique(uid)
  if (isTRUE(time_effects)) {
    # common shocks: every variable is centred across the units within each period
    if (is.null(time)) stop("Removing time effects needs the `time` column.", call. = FALSE)
    tt <- data[[time]]
    for (v in vars) { z <- as.numeric(data[[v]]); data[[v]] <- z - stats::ave(z, tt, FUN = function(a) mean(a, na.rm = TRUE)) }
  }
  U <- lapply(ids, function(u) {
    d <- data[uid == u, vars, drop = FALSE]; rownames(d) <- NULL
    d[] <- lapply(d, function(v) {
      v <- as.numeric(v)
      if (all(!is.finite(v))) stop(sprintf("A variable of unit '%s' has no finite value.", u), call. = FALSE)
      if (impute) .prep_na(v) else v
    })
    if (any(!is.finite(as.matrix(d)))) stop("Missing/non-finite values present (use impute = TRUE).", call. = FALSE)
    d
  })
  names(U) <- ids
  list(units = U, vars = vars, ids = ids, n = vapply(U, nrow, 0L), data = data)
}

#' Long panel data from matrices or arrays
#'
#' Converts the wide layouts in which panel series are often stored to the
#' long data frame expected by [rcausality()] / [panel_causality()]: every
#' named argument in `...` is a `time x unit` matrix of one variable, and
#' `array` is an optional `variable x unit x time` array (the layout of the
#' scripts of Kosztyan, Matrai and Kiss, 2026).
#'
#' @param ... named `time x unit` matrices (or data frames), one per variable.
#' @param array optional three-dimensional array, `variable x unit x time`.
#' @param units,times optional unit and time labels (default: the dimnames,
#'   or `u1, u2, ...` and `1, 2, ...`).
#' @return A data frame with the columns `unit`, `time` and one column per
#'   variable.
#' @export
#' @examples
#' Y <- matrix(rnorm(40), 10, 4)                    # 10 periods, 4 units
#' X <- array(rnorm(80), c(2, 4, 10), list(c("a", "b"), NULL, NULL))
#' str(panel_long(y = Y, array = X))
panel_long <- function(..., array = NULL, units = NULL, times = NULL) {
  mats <- list(...)
  if (length(mats) && (is.null(names(mats)) || any(!nzchar(names(mats))))) stop("The matrices in `...` must be named.", call. = FALSE)
  mats <- lapply(mats, as.matrix)
  if (!is.null(array)) {
    if (length(dim(array)) != 3L) stop("`array` must be a variable x unit x time array.", call. = FALSE)
    vn <- dimnames(array)[[1]]; if (is.null(vn)) vn <- paste0("x", seq_len(dim(array)[1]))
    for (i in seq_len(dim(array)[1])) {
      m <- t(matrix(array[i, , ], dim(array)[2], dim(array)[3]))
      if (is.null(units) && !is.null(dimnames(array)[[2]])) colnames(m) <- dimnames(array)[[2]]
      mats[[vn[i]]] <- m
    }
  }
  if (!length(mats)) stop("No data given.", call. = FALSE)
  nt <- nrow(mats[[1]]); nu <- ncol(mats[[1]])
  if (any(vapply(mats, nrow, 0L) != nt) || any(vapply(mats, ncol, 0L) != nu))
    stop("All variables must have the same number of periods (rows) and units (columns).", call. = FALSE)
  cn <- Filter(Negate(is.null), lapply(mats, colnames))
  if (is.null(units)) units <- if (length(cn)) cn[[1]] else paste0("u", formatC(seq_len(nu), width = nchar(nu), flag = "0"))
  rn <- Filter(Negate(is.null), lapply(mats, rownames))
  if (is.null(times)) times <- if (length(rn) && !anyNA(suppressWarnings(as.numeric(rn[[1]])))) as.numeric(rn[[1]]) else seq_len(nt)
  out <- data.frame(unit = rep(as.character(units), each = nt), time = rep(times, times = nu), stringsAsFactors = FALSE)
  for (nm in names(mats)) out[[nm]] <- as.numeric(mats[[nm]])
  out
}

#' Settings of the panel-data tests
#'
#' @param levels which levels to compute: `"units"` (the method in every
#'   unit), `"group"` (the panel test built from the unit results) and
#'   `"pooled"` (one estimate from all units, for the methods marked `pooled`
#'   in [causal_methods()]).
#' @param decision the level that decides the links of the panel (`$pairs`,
#'   `$classes`): `"group"`, `"pooled"`, `"both"` or `"auto"` (= `"both"`
#'   when both levels are available, otherwise the available one). With
#'   `"both"` a link is declared when either level shows it: the smaller of
#'   the two p-values, doubled (a Bonferroni correction over the two tests),
#'   or -- for Bayes-factor decisions -- the average of the two Bayes
#'   factors. The two levels are sensitive to different alternatives: the
#'   pooled test to an effect that has the same direction in most units, the
#'   group test to an effect in some of the units, of any sign.
#' @param combine the group-level test:
#'   `"auto"` (the natural one for the method: `"mean"` for the methods whose
#'   unit statistic has a bootstrap or surrogate null distribution,
#'   `"stouffer"` for the methods with asymptotic p-values only, `"count"`
#'   for the methods whose decision is a table or a verdict), or one of
#'   * `"mean"`: the sum of the unit statistics, each standardised with the
#'     mean and the standard deviation of its own null sample (a
#'     nonparametric combination test, Pesarin and Salmaso 2010), against the
#'     convolution of the unit null distributions (Monte Carlo with
#'     `n_draws` draws);
#'   * `"dh"`: Dumitrescu and Hurlin (2012), the mean of the unit Wald
#'     statistics standardised with the moments of the F distribution
#'     (Granger family). The normal approximation is anti-conservative for
#'     short and persistent series; the statistic is always reported next to
#'     the `"mean"` test;
#'   * `"stouffer"`, `"fisher"`, `"simes"`: combinations of the unit
#'     p-values;
#'   * `"count"`: the number of units with a link at level `alpha` against
#'     the Binomial(`n_units`, `alpha`) distribution.
#' @param unit_adjust correction of the unit-level p-values over all unit x
#'   link cells (passed to [stats::p.adjust()]; default Benjamini-Hochberg).
#' @param prevalence_prior the two parameters of the Beta prior of the share
#'   of units with the link in the prevalence model (default uniform).
#' @param min_units smallest number of units for a group-level test.
#' @param min_length units shorter than this are not analysed on their own
#'   (`NULL` = the advisable minimal length of the method, see
#'   [causal_methods()]); they still enter the pooled estimate.
#' @param time_effects logical; remove common shocks first, i.e. centre every
#'   variable across the units within each period (needs `time`)? The tests
#'   treat the units as independent; with strong common factors this
#'   transformation is advisable.
#' @param n_draws Monte Carlo size of the `"mean"` group test.
#' @return A list of class `panel_control`.
#' @references
#' Dumitrescu, E.-I. and Hurlin, C. (2012). Testing for Granger non-causality in heterogeneous panels. *Economic Modelling*, 29(4), 1450-1460. \doi{10.1016/j.econmod.2012.02.014}
#'
#' Pesarin, F. and Salmaso, L. (2010). *Permutation Tests for Complex Data: Theory, Applications and Software*. Wiley, Chichester. \doi{10.1002/9780470689516}
#' @seealso [panel_causality()]
#' @export
#' @examples
#' panel_control(decision = "pooled", combine = "fisher")
panel_control <- function(levels = c("units", "group", "pooled"), decision = c("auto", "both", "group", "pooled"),
                          combine = c("auto", "dh", "mean", "stouffer", "fisher", "simes", "count"),
                          unit_adjust = "BH", prevalence_prior = c(1, 1), min_units = 3L, min_length = NULL,
                          time_effects = FALSE, n_draws = 10000L) {
  levels <- match.arg(levels, several.ok = TRUE); decision <- match.arg(decision); combine <- match.arg(combine)
  if (!unit_adjust %in% stats::p.adjust.methods) stop("`unit_adjust` must be a method of stats::p.adjust().", call. = FALSE)
  if (length(prevalence_prior) != 2L || any(!is.finite(prevalence_prior)) || any(prevalence_prior <= 0))
    stop("`prevalence_prior` must hold two positive Beta parameters.", call. = FALSE)
  .check_number(min_units, "min_units", lower = 2, integer = TRUE)
  if (!is.null(min_length)) .check_number(min_length, "min_length", lower = 4, integer = TRUE)
  .check_number(n_draws, "n_draws", lower = 100, integer = TRUE)
  if ("group" %in% levels) levels <- union("units", levels)
  structure(list(levels = levels, decision = decision, combine = combine, unit_adjust = unit_adjust,
                 prevalence_prior = as.numeric(prevalence_prior), min_units = as.integer(min_units),
                 min_length = min_length, time_effects = isTRUE(time_effects), n_draws = as.integer(n_draws)), class = "panel_control")
}

#' @export
print.panel_control <- function(x, ...) {
  cat("Panel-test settings\n")
  cat(sprintf("  levels: %s; decision: %s; group test: %s; unit-level correction: %s\n", paste(x$levels, collapse = ", "),
              x$decision, x$combine, x$unit_adjust))
  cat(sprintf("  prevalence prior: Beta(%g, %g); min_units = %d; min_length = %s; time effects removed: %s; n_draws = %d\n", x$prevalence_prior[1],
              x$prevalence_prior[2], x$min_units, if (is.null(x$min_length)) "method default" else x$min_length,
              if (isTRUE(x$time_effects)) "yes" else "no", x$n_draws))
  invisible(x)
}

# ---- group-level tests -----------------------------------------------------------

# Dumitrescu-Hurlin from the unit Wald statistics W_i = df1 * F_i; the exact
# moments of K * F(K, d2) under H0 give the finite-T standardisation Z-tilde
.dh_test <- function(W, K, d2) {
  K <- rep_len(K, length(W)); ok <- is.finite(W) & is.finite(d2) & d2 > 4
  W <- W[ok]; d2 <- d2[ok]; K <- K[ok]; N <- length(W)
  if (N < 2L) return(c(W_bar = NA_real_, Z_bar = NA_real_, Z_tilde = NA_real_, p_value = NA_real_, n_units = N))
  EW <- K * d2 / (d2 - 2); VW <- 2 * K * d2^2 * (K + d2 - 2) / ((d2 - 2)^2 * (d2 - 4))
  Zb <- sqrt(N / 2) * mean((W - K) / sqrt(K))
  Zt <- (sum(W) - sum(EW)) / sqrt(sum(VW))
  c(W_bar = mean(W), Z_bar = Zb, Z_tilde = Zt, p_value = 2 * stats::pnorm(-abs(Zt)), n_units = N)
}

# Group mean of the standardised unit statistics (a nonparametric combination
# test, Pesarin and Salmaso 2010).  `nulls`: per unit a list (statistic,
# null, alternative) with `null` either surrogate draws or list(mean, sd).
# In every unit the observed statistic and its surrogates form one pool; the
# statistic is standardised with the mean and the standard deviation of the
# pool, and the null distribution of the sum over the units is the
# convolution of the standardised pools (Monte Carlo).  Under the null
# hypothesis the observed value is exchangeable with its surrogates, so the
# test keeps its level for any shape of the unit distributions, any unit
# length and any number of surrogates.
.group_mean_test <- function(nulls, n_draws = 10000L) {
  std <- lapply(nulls, function(e) {
    if (is.null(e) || !is.finite(e$statistic)) return(NULL)
    dev <- switch(if (is.null(e$alternative)) "greater" else e$alternative, less = function(v) -v, two.sided = abs, identity)
    if (is.list(e$null)) {
      if (!is.finite(e$null$sd) || e$null$sd <= 0) return(NULL)
      return(list(z = dev((e$statistic - e$null$mean) / e$null$sd), draws = NULL, two = identical(e$alternative, "two.sided")))
    }
    pool <- c(e$statistic, e$null[is.finite(e$null)]); s <- stats::sd(pool)
    if (length(pool) < 3L || !is.finite(s) || s <= 0) return(NULL)
    zz <- dev((pool - mean(pool)) / s)
    list(z = zz[1L], draws = zz, two = FALSE)
  })
  std <- std[!vapply(std, is.null, TRUE)]; N <- length(std)
  if (N < 2L) return(list(statistic = NA_real_, p_value = NA_real_, n_units = N, z = numeric(0)))
  z <- vapply(std, `[[`, 0, "z"); Z <- sum(z) / sqrt(N)
  sim <- vapply(std, function(e) if (is.null(e$draws)) { v <- stats::rnorm(n_draws); if (e$two) abs(v) else v }
                                 else e$draws[sample.int(length(e$draws), n_draws, replace = TRUE)], numeric(n_draws))
  Zs <- rowSums(sim) / sqrt(N)
  list(statistic = Z, p_value = (1 + sum(Zs >= Z - 1e-12)) / (1 + n_draws), n_units = N, z = z, null = Zs)
}

# the lag order of one link, the same in every unit: the BIC of its pooled
# fixed-effects regression, all candidate orders compared on the same sample
.panel_pair_lag <- function(units, ef, cs, ctrl, max_lag, aug = 0L) {
  cap <- .pooled_cap(units, max_lag, aug, 2L + length(ctrl)); lags <- seq_len(cap)
  bic <- vapply(lags, function(l) tryCatch(.pooled_F(.panel_design(units, ef, cs, ctrl, l, aug = aug, start = cap + aug))$bic,
                                           error = function(e) NA_real_), 0)
  if (all(!is.finite(bic))) 1L else lags[which.min(bic)]
}

.combine_p <- function(p, how = c("stouffer", "fisher", "simes")) {
  how <- match.arg(how); p <- p[is.finite(p)]; N <- length(p)
  if (N < 2L) return(list(statistic = NA_real_, p_value = NA_real_, n_units = N))
  p <- pmin(pmax(p, 1e-15), 1 - 1e-15)
  switch(how,
    stouffer = { Z <- sum(stats::qnorm(1 - p)) / sqrt(N); list(statistic = Z, p_value = stats::pnorm(Z, lower.tail = FALSE), n_units = N) },
    fisher = { X <- -2 * sum(log(p)); list(statistic = X, p_value = stats::pchisq(X, 2 * N, lower.tail = FALSE), n_units = N) },
    simes = { ps <- min(sort(p) * N / seq_len(N)); list(statistic = -log10(ps), p_value = min(1, ps), n_units = N) })
}

# Bayesian prevalence model on the unit Bayes factors: the share `pi` of the
# units with the link has a Beta prior; the marginal likelihood of the panel
# relative to "no unit has the link" is  E_pi[ prod_i (1 - pi + pi BF_i) ].
# Returns the Bayes factor of "pi > 0", the posterior of pi and the posterior
# link probability of every unit (partial pooling: the units borrow the
# panel's estimate of the prevalence as their prior).
.prevalence_bf <- function(bf, prior = c(1, 1), credible = 0.95) {
  ok <- is.finite(bf); N <- sum(ok)
  out <- list(bf10 = NA_real_, mean = NA_real_, lower = NA_real_, upper = NA_real_, unit_prob = rep(NA_real_, length(bf)), n_units = N)
  if (N < 1L) return(out)
  lb <- log(pmin(pmax(bf[ok], 1e-300), 1e300))
  g <- (seq_len(1000L) - 0.5) / 1000
  lw <- stats::dbeta(g, prior[1], prior[2], log = TRUE); lw <- lw - .logsumexp(lw)
  em1 <- expm1(pmin(lb, 700))                                           # BF_i - 1
  ll <- vapply(g, function(p) sum(log1p(pmax(p * em1, -1 + 1e-300))), 0)
  lpost <- ll + lw; lm1 <- .logsumexp(lpost)
  post <- exp(lpost - lm1); cdf <- cumsum(post)
  a <- (1 - credible) / 2
  up <- vapply(seq_len(N), function(i) sum(post * exp(log(g) + lb[i] - log1p(pmax(g * em1[i], -1 + 1e-300)))), 0)
  out$bf10 <- exp(max(min(lm1, 700), -700)); out$mean <- sum(post * g)
  out$lower <- g[which(cdf >= a)[1]]; out$upper <- g[which(cdf >= 1 - a)[1]]
  out$unit_prob[ok] <- pmin(pmax(up, 0), 1)
  out
}

# Bayesian false-discovery-rate control (Newton et al. 2004): the largest set
# of the most probable cells whose mean posterior null probability is at most alpha
.bayes_fdr <- function(post, alpha) {
  sel <- rep(FALSE, length(post)); ok <- which(is.finite(post))
  if (!length(ok)) return(sel)
  o <- ok[order(post[ok], decreasing = TRUE)]
  k <- max(c(0L, which(cumsum(1 - post[o]) / seq_along(o) <= alpha)))
  if (k > 0L) sel[o[seq_len(k)]] <- TRUE
  sel
}

# ---- pooled (stacked) estimators ---------------------------------------------------

.unit_shift <- function(v, s) { n <- length(v); if (s <= 0L || s >= n) v else c(v[(s + 1):n], v[seq_len(s)]) }

# surrogate panel of the model-free pooled tests: in every unit the cause is
# delayed circularly by a random number of steps between the lag window (at
# least 10% of the unit) and 60% of its length (0 for units too short to
# shift).  Delaying does not move the future of the cause next to the present
# of the effect, which would bring a link in the reverse direction into the
# surrogates.
.panel_shifts <- function(units, L) {
  vapply(units, function(u) {
    n <- nrow(u); lo <- as.integer(max(L + 1L, ceiling(0.1 * n))); hi <- as.integer(floor(0.6 * n))
    if (hi < lo) 0L else as.integer(n - (lo + sample.int(hi - lo + 1L, 1L) - 1L))
  }, 0L)
}

# stacked lag design of one directed pair: y, nuisance regressors (own and
# control lags 1..L + aug, and the cause's augmenting lags L+1..L+aug), the
# cause's lags 1..L, unit index and position within the unit.  scale =
# "demean": exact within transformation (unit means of every column removed);
# "standardize": the unit series are z-scored first and the columns centred
# by unit; `shift`: circular shifts of the cause per unit (surrogates).
.panel_design <- function(units, ef, cs, ctrl = character(0), L = 1L, shift = NULL, scale = c("demean", "standardize"),
                          aug = 0L, start = NULL) {
  scale <- match.arg(scale); LL <- L + aug; skip <- if (is.null(start)) 0L else max(0L, as.integer(start) - LL)
  pr <- if (scale == "standardize") .std else as.numeric
  parts <- lapply(seq_along(units), function(i) {
    u <- units[[i]]; n <- nrow(u)
    if (n < LL + skip + 3L) return(NULL)
    x <- pr(u[[cs]]); if (!is.null(shift)) x <- .unit_shift(x, shift[i])
    e <- pr(u[[ef]])
    Xs <- .lagmat(e, LL)
    for (c in ctrl) Xs <- cbind(Xs, .lagmat(pr(u[[c]]), LL))
    Xc <- .lagmat(x, LL)
    if (aug > 0L) { Xs <- cbind(Xs, Xc[, (L + 1L):LL, drop = FALSE]); Xc <- Xc[, seq_len(L), drop = FALSE] }
    keep <- (skip + 1L):(n - LL)
    list(y = e[(LL + 1L):n][keep], Xs = Xs[keep, , drop = FALSE], Xc = Xc[keep, , drop = FALSE], unit = rep(i, length(keep)),
         pos = keep)
  })
  parts <- parts[!vapply(parts, is.null, TRUE)]
  if (!length(parts)) stop("Every unit is too short for this lag order.", call. = FALSE)
  y <- unlist(lapply(parts, `[[`, "y")); Xs <- do.call(rbind, lapply(parts, `[[`, "Xs"))
  Xc <- do.call(rbind, lapply(parts, `[[`, "Xc")); unit <- unlist(lapply(parts, `[[`, "unit"))
  f <- factor(unit); cnt <- as.vector(table(f))[as.integer(f)]
  dm <- function(M) { M <- as.matrix(M); M - rowsum(M, f)[as.integer(f), , drop = FALSE] / cnt }
  list(y = as.numeric(dm(y)), Xs = dm(Xs), Xc = dm(Xc), unit = unit, pos = unlist(lapply(parts, `[[`, "pos")),
       n_units = nlevels(f))
}

# largest lag order usable in the pooled design
# (the rule of the time-series tests, .cap_lag(), for a unit of median length: the
# unit regressions and their half-panel versions are then estimable)
.pooled_cap <- function(units, max_lag, aug = 0L, nvars = 2L)
  max(1L, min(as.integer(max_lag), .cap_lag(stats::median(vapply(units, nrow, 0L)), nvars, aug)))
# largest lag scanned by the pooled model-free statistics
.delay_cap <- function(units, max_lag) max(1L, min(as.integer(max_lag), floor(stats::median(vapply(units, nrow, 0L)) / 4)))

.pooled_F <- function(ds, prior = "bic") {
  N <- length(ds$y); q <- ncol(ds$Xc); k0 <- ncol(ds$Xs)
  fr <- stats::lm.fit(ds$Xs, ds$y); fu <- stats::lm.fit(cbind(ds$Xs, ds$Xc), ds$y)
  rss_r <- sum(fr$residuals^2); rss_u <- sum(fu$residuals^2)
  df2 <- N - ds$n_units - k0 - q
  if (df2 < 1 || rss_u <= 0) return(list(statistic = NA_real_, df1 = q, df2 = df2, p_value = NA_real_, bf10 = NA_real_,
                                         effect_size = NA_real_, r2_gain = NA_real_, sign = NA_real_, n_eff = N, bic = NA_real_))
  Fs <- ((rss_r - rss_u) / q) / (rss_u / df2)
  lbf <- if (prior %in% c("g", "hyper_g")) {
    rss0 <- sum(ds$y^2)
    .log_bf_nested(N - ds$n_units + 1L, k0, k0 + q, 1 - rss_r / rss0, 1 - rss_u / rss0, prior)
  } else (N * log(rss_r / rss_u) - q * log(N)) / 2
  cc <- fu$coefficients[k0 + seq_len(q)]
  list(statistic = Fs, df1 = q, df2 = df2, p_value = stats::pf(Fs, q, df2, lower.tail = FALSE), bf10 = exp(max(min(lbf, 700), -700)),
       effect_size = log(rss_r / rss_u), r2_gain = 1 - rss_u / rss_r, sign = sign(sum(cc, na.rm = TRUE)), n_eff = N,
       bic = log(rss_u / N) + (k0 + q + ds$n_units) * log(N) / N)
}

# surrogate test of a pooled statistic: `stat(shift)` computes the statistic of
# the stacked design with the cause shifted circularly within the units
.pooled_surrogate <- function(stat, units, L, n_perm, cores) {
  obs <- stat(NULL)
  null <- unlist(.plapply(seq_len(n_perm), function(b) stat(.panel_shifts(units, L)), cores = cores))
  list(statistic = obs, null = null, p_value = (1 + sum(null >= obs, na.rm = TRUE)) / (1 + sum(is.finite(null))))
}

.pooled_rf <- function(units, ef, cs, ctrl, A, cores = 1L) {
  .need_ranger("the pooled random-forest test")
  st <- .method_args("rf_granger", A$method_args)
  L <- if (is.numeric(A$lag)) max(1L, as.integer(A$lag)) else .panel_pair_lag(units, ef, cs, ctrl, A$max_lag)
  d0 <- .panel_design(units, ef, cs, ctrl, L, scale = "standardize")
  f_own <- .rf_fit(d0$Xs, d0$y, st); mse_own <- f_own$prediction.error
  s <- .pooled_surrogate(function(shift) {
    d <- if (is.null(shift)) d0 else .panel_design(units, ef, cs, ctrl, L, shift = shift, scale = "standardize")
    1 - .rf_fit(cbind(d$Xs, d$Xc), d$y, st)$prediction.error / mse_own
  }, units, L, A$n_perm, cores)
  p <- if (s$statistic > 0) s$p_value else 1
  list(row = list(lag = L, statistic = s$statistic, p_value = p, effect_size = s$statistic, r2_gain = s$statistic, n_eff = length(d0$y)),
       null = list(statistic = s$statistic, null = s$null, alternative = "greater"),
       detail = sprintf("pooled forest, out-of-bag gain %+.3f", s$statistic), gate = s$statistic > 0)
}

.pooled_dp <- function(units, ef, cs, ctrl, A, cores = 1L) {
  st <- .method_args("diks_panchenko", A$method_args)
  L <- if (is.numeric(A$lag)) max(1L, as.integer(A$lag)) else 1L
  stat <- function(shift) {
    d <- .panel_design(units, ef, cs, ctrl, L, shift = shift, scale = "standardize"); N <- length(d$y)
    eps <- if (is.null(st$eps)) min(1.5, st$C * N^(-2 / 7)) else st$eps
    cc <- dp_counts_cpp(d$Xc, d$Xs, d$y, eps)
    (N - 1) / (N * (N - 2)) * sum(cc$cxyz * cc$cy - cc$cxy * cc$cyz)
  }
  s <- .pooled_surrogate(stat, units, L, A$n_perm, cores)
  z <- (s$statistic - mean(s$null)) / stats::sd(s$null)
  list(row = list(lag = L, statistic = s$statistic, p_value = s$p_value, effect_size = z, n_eff = sum(pmax(vapply(units, nrow, 0L) - L, 0L))),
       null = list(statistic = s$statistic, null = s$null, alternative = "greater"),
       detail = sprintf("stacked delay vectors, z = %.2f against within-unit shifts", z), gate = TRUE)
}

# pooled kNN conditional mutual information: the delay vectors of the
# standardised unit series are stacked; local permutation test on the pooled
# sample, max-statistic over the candidate lags as in the time-series test
.pooled_cmi <- function(units, ef, cs, ctrl, A, cause_past = FALSE, cores = 1L) {
  lags <- if (is.numeric(A$lag)) max(1L, as.integer(A$lag)) else seq_len(.delay_cap(units, min(A$max_lag, 5L)))
  history <- if (identical(A$history, "auto")) {
    .panel_pair_lag(units, ef, cs, ctrl, min(5L, A$max_lag))
  } else max(1L, as.integer(A$history))
  k <- if (is.null(A$k)) 5L else as.integer(A$k)
  Lmax <- if (cause_past) max(lags) + history else max(lags, history)
  parts <- lapply(units, function(u) {
    n <- nrow(u); if (n < Lmax + 3L) return(NULL)
    e <- .std(u[[ef]]); x <- .std(u[[cs]]); t_idx <- (Lmax + 1L):n
    Z <- .hist_cols(e, t_idx, history)
    for (c in ctrl) Z <- cbind(Z, .hist_cols(.std(u[[c]]), t_idx, history))
    list(y = e[t_idx], Z = Z, X = lapply(lags, function(tau) x[t_idx - tau]),
         ZX = if (cause_past) lapply(lags, function(tau) .hist_cols(x, t_idx, history, shift = tau)))
  })
  parts <- parts[!vapply(parts, is.null, TRUE)]
  if (!length(parts)) stop("Every unit is too short for the pooled conditional mutual information.", call. = FALSE)
  yt <- unlist(lapply(parts, `[[`, "y")); Z <- do.call(rbind, lapply(parts, `[[`, "Z"))
  res <- lapply(seq_along(lags), function(j) {
    xt <- unlist(lapply(parts, function(p) p$X[[j]]))
    Zt <- if (cause_past) cbind(Z, do.call(rbind, lapply(parts, function(p) p$ZX[[j]]))) else Z
    local_permutation_test(xt, yt, Zt, k = k, k_perm = A$k_perm, n_perm = A$n_perm, cores = cores, standardize = FALSE)
  })
  obs <- vapply(res, `[[`, 0, "cmi"); null_max <- do.call(pmax, lapply(res, `[[`, "null")); b <- which.max(obs)
  list(row = list(lag = lags[b], statistic = obs[b], p_value = (1 + sum(null_max >= obs[b])) / (length(null_max) + 1),
                  effect_size = obs[b], n_eff = length(yt)),
       by_lag = data.frame(lag = lags, statistic = obs, p = vapply(res, `[[`, 0, "p_value")),
       null = list(statistic = obs[b], null = null_max, alternative = "greater"),
       detail = sprintf("stacked delay vectors (history %d), local permutation", history), gate = TRUE)
}

# pooled binned transfer entropy: symbols coded within every unit, the
# conditional entropies from the stacked symbol histories; within-unit
# circular shifts of the cause, max-statistic over the lags
.pooled_te <- function(units, ef, cs, A, cores = 1L) {
  bins <- as.integer(A$bins)
  lags <- if (is.numeric(A$lag)) max(1L, as.integer(A$lag)) else seq_len(.delay_cap(units, min(A$max_lag, 5L)))
  l <- if (identical(A$history, "auto")) {
    .panel_pair_lag(units, ef, cs, character(0), min(5L, A$max_lag))
  } else max(1L, as.integer(A$history))
  Lmax <- max(lags, l)
  ok <- vapply(units, nrow, 0L) >= Lmax + 3L
  if (!any(ok)) stop("Every unit is too short for the pooled transfer entropy.", call. = FALSE)
  U <- units[ok]
  sy <- lapply(U, function(u) .sym_bins(u[[ef]], bins)); sx <- lapply(U, function(u) .sym_bins(u[[cs]], bins))
  nb <- max(unlist(sy), unlist(sx)) + 1L
  idx <- lapply(U, function(u) (Lmax + 1L):nrow(u))
  ty <- unlist(Map(function(s, t) s[t], sy, idx))
  Hy <- do.call(rbind, Map(function(s, t) do.call(cbind, lapply(seq_len(l), function(j) s[t - j])), sy, idx))
  storage.mode(Hy) <- "integer"; h0 <- cond_entropy_sym_cpp(ty, Hy, nb, nb)
  te_all <- function(shift) vapply(lags, function(tau) {
    hx <- unlist(Map(function(s, t, sh) .unit_shift(s, sh)[t - tau], sx, idx, if (is.null(shift)) rep(0L, length(U)) else shift))
    Hxy <- cbind(Hy, hx); storage.mode(Hxy) <- "integer"
    h0 - cond_entropy_sym_cpp(ty, Hxy, nb, nb)
  }, 0)
  obs <- te_all(NULL); b <- which.max(obs)
  null <- do.call(rbind, .plapply(seq_len(A$n_perm), function(r) te_all(.panel_shifts(U, Lmax)), cores = cores))
  null_max <- apply(null, 1, max)
  list(row = list(lag = lags[b], statistic = obs[b], p_value = (1 + sum(null_max >= obs[b])) / (length(null_max) + 1),
                  effect_size = obs[b] - mean(null[, b]), n_eff = length(ty)),
       by_lag = data.frame(lag = lags, statistic = obs, p = (1 + colSums(sweep(null, 2, obs, ">="))) / (nrow(null) + 1)),
       null = list(statistic = obs[b], null = null_max, alternative = "greater"),
       detail = sprintf("stacked symbol histories (%d bins, history %d)", bins, l), gate = TRUE)
}

# multispatial convergent cross mapping (Clark et al. 2015): the delay
# embeddings of the effect in all units form one library
.pooled_ccm <- function(units, ef, cs, A, cores = 1L) {
  E <- as.integer(A$m); tau <- as.integer(A$tau); span <- (E - 1L) * tau
  ok <- vapply(units, nrow, 0L) >= span + 4L
  if (!any(ok)) stop("Every unit is too short for the embedding.", call. = FALSE)
  U <- units[ok]
  My <- .num_mat(do.call(rbind, lapply(U, function(u) embed_delay(.std(u[[ef]]), E, tau)))); ne <- nrow(My)
  xs <- lapply(U, function(u) .std(u[[cs]])[(span + 1L):nrow(u)])
  lib_sizes <- sort(unique(as.integer(round(exp(seq(log(E + 2L), log(ne), length.out = 10L))))))
  lib_sizes <- lib_sizes[lib_sizes >= E + 2L & lib_sizes <= ne]
  if (!length(lib_sizes)) stop("Panel too small for the embedding.", call. = FALSE)
  skill <- function(xv, Ls, n_samples = 20L) mean(vapply(seq_len(if (Ls >= ne) 1L else n_samples), function(s) {
    lib <- if (Ls >= ne) seq_len(ne) else sort(sample.int(ne, Ls))
    pr <- simplex_predict_cpp(My, lib, lib, xv, E + 1L)
    r <- suppressWarnings(stats::cor(xv[lib], pr, use = "complete.obs")); if (is.na(r)) 0 else r
  }, 0))
  xt <- unlist(xs)
  rho <- vapply(lib_sizes, function(Ls) skill(xt, Ls), 0); rho_max <- rho[length(rho)]
  kt <- if (length(rho) >= 3) suppressWarnings(stats::cor.test(lib_sizes, rho, method = "kendall", alternative = "greater")$p.value) else NA_real_
  null <- unlist(.plapply(seq_len(A$n_perm), function(b) {
    xv <- unlist(lapply(xs, function(v) { n <- length(v); lo <- max(span + 2L, ceiling(0.1 * n)); hi <- floor(0.6 * n)
      if (hi < lo) sample(v) else .unit_shift(v, n - (lo + sample.int(hi - lo + 1L, 1L) - 1L)) }))
    skill(xv, lib_sizes[length(lib_sizes)])
  }, cores = cores))
  conv <- !is.na(kt) && kt < 0.05 && rho[length(rho)] > rho[1]
  list(row = list(lag = NA_integer_, statistic = rho_max, p_value = (1 + sum(null >= rho_max)) / (length(null) + 1),
                  effect_size = rho_max, r2_gain = rho_max^2, n_eff = ne),
       by_lag = data.frame(lag = lib_sizes, statistic = rho, p = NA_real_),
       null = list(statistic = rho_max, null = null, alternative = "greater"),
       detail = sprintf("multispatial library of %d delay vectors, %s", ne, if (conv) "converging" else "not converging"), gate = conv)
}

.pooled_test <- function(method, units, ef, cs, ctrl, A, cores = 1L) {
  switch(method,
    rf_granger = .pooled_rf(units, ef, cs, ctrl, A, cores),
    diks_panchenko = .pooled_dp(units, ef, cs, ctrl, A, cores),
    cmi = .pooled_cmi(units, ef, cs, ctrl, A, FALSE, cores),
    cmi_lag = .pooled_cmi(units, ef, cs, ctrl, A, TRUE, cores),
    te = .pooled_te(units, ef, cs, A, cores),
    ccm = .pooled_ccm(units, ef, cs, A, cores),
    NULL)
}

# ---- the panel analysis --------------------------------------------------------------

.group_default <- c(granger = "mean", bayes_granger = "mean", toda_yamamoto = "mean", diks_panchenko = "mean", rf_granger = "mean",
                    cmi = "mean", cmi_lag = "mean", te = "mean", pid = "mean", ccm = "count", ccdh = "count", hirata = "count",
                    cd = "count", df = "count", markov = "count", dc = "count", dfe = "count", pcmci = "stouffer")

# `fam_fun(cause, effect)` (optional) names the family of every link: the
# multiplicity correction `adjust` is then applied within the families
.panel_run <- function(x, y = NULL, A, cl = NULL, xvars = NULL, yvars = NULL, fam_fun = NULL, adjust = A$adjust) {
  if (!is.null(y)) stop("Panel data: give one long data frame with the unit column `id`; use `target` for one explained variable.",
                        call. = FALSE)
  pc <- A$panel; if (!inherits(pc, "panel_control")) pc <- do.call(panel_control, as.list(pc))
  method <- A$method; ctl <- A$control
  if (!inherits(ctl, "bayes_control")) ctl <- do.call(bayes_control, as.list(ctl))
  if (!is.null(A[["prior"]])) ctl$prior <- match.arg(A[["prior"]], c("bic", "g", "hyper_g", "spike_slab", "normal", "minnesota"))
  A$control <- ctl; A[["prior"]] <- NULL
  if (method == "bayes_granger") { method <- "granger"; A$bayes <- TRUE }
  bayes <- isTRUE(A$bayes); A$method <- method
  info <- .method_info(method)
  if (isTRUE(A$conditional) && !info$conditional) stop(sprintf("method = \"%s\" supports pairwise analysis only.", method), call. = FALSE)
  .check_number(A$alpha, "alpha", lower = 0, upper = 1, strict = TRUE)
  .check_number(A$max_lag, "max_lag", lower = 1, integer = TRUE)
  .check_number(A$n_perm, "n_perm", lower = 1, integer = TRUE)
  if (!adjust %in% stats::p.adjust.methods) stop("`adjust` must be one of ", paste(stats::p.adjust.methods, collapse = ", "), ".", call. = FALSE)
  if (is.character(A$lag) && !A$lag %in% c("auto", "first")) stop("`lag` must be an integer, \"auto\" or \"first\".", call. = FALSE)
  PD <- .panel_units(x, A$id, A$time, vars = if (is.null(xvars)) NULL else c(xvars, yvars), impute = A$impute, time_effects = pc$time_effects)
  U <- PD$units; vars <- PD$vars
  if (is.null(xvars)) xvars <- vars
  if (length(U) < 2L) stop("A panel needs at least two units.", call. = FALSE)
  causes <- xvars; effects <- if (is.null(yvars)) xvars else yvars
  grid <- expand.grid(cause = causes, effect = effects, stringsAsFactors = FALSE)
  if (!is.null(yvars) && A$directions == "both") {
    grid <- rbind(grid, stats::setNames(grid[, 2:1], c("cause", "effect"))); causes <- effects <- vars
  }
  grid <- grid[grid$cause != grid$effect, , drop = FALSE]; rownames(grid) <- NULL
  key <- paste(grid$cause, "->", grid$effect)
  ctrl_of <- function(j) if (isTRUE(A$conditional)) setdiff(vars, c(grid$cause[j], grid$effect[j])) else character(0)
  old_rng <- .rng_setup(A$seed); on.exit(.rng_restore(old_rng), add = TRUE)
  .n_cores(A$cores)
  # resampling tests of the pooled level: enough surrogates for the multiplicity correction
  if (!method %in% c("granger", "toda_yamamoto") && isTRUE(info$pooled) && "pooled" %in% pc$levels && adjust != "none") {
    both <- pc$decision %in% c("auto", "both") && "group" %in% pc$levels        # the smaller p-value is doubled
    need <- min(999L, ceiling((1 + both) * nrow(grid) / A$alpha))
    if (A$n_perm < need) {
      message(method, ": n_perm raised from ", A$n_perm, " to ", need, " so that the surrogate p-values of the pooled tests can reach ",
              "significance after the multiplicity correction.")
      A$n_perm <- need
    }
  }
  use_bf <- bayes || info$family == "posterior"
  min_len <- if (is.null(pc$min_length)) info$min_n else pc$min_length
  long <- PD$n >= min_len

  # ---- the lag order: one order per link, the same in every unit, when it is chosen from the data
  fam_reg <- method %in% c("granger", "toda_yamamoto")
  aug <- if (method == "toda_yamamoto") as.integer(.method_args("toda_yamamoto", A$method_args)$d_max) else 0L
  first <- identical(A$lag, "first") && fam_reg; lag_pair <- rep(NA_integer_, nrow(grid))
  if (identical(A$lag, "auto") && (fam_reg || method == "rf_granger"))
    lag_pair <- vapply(seq_len(nrow(grid)), function(j) .panel_pair_lag(U, grid$effect[j], grid$cause[j], ctrl_of(j), A$max_lag), 0L)
  lag_set <- if (anyNA(lag_pair)) NA_integer_ else sort(unique(lag_pair))
  # the scan over the orders (lag = "first") uses the same orders in all units and in the pooled test
  n_ref <- if (any(long)) min(PD$n[long]) else as.integer(stats::median(PD$n))
  scan_lags <- seq_len(max(1L, min(as.integer(A$max_lag), .cap_lag(n_ref, 2L + if (isTRUE(A$conditional)) length(vars) - 2L else 0L, aug))))

  # ---- Granger family: the engine (unit F statistics and pooled test of the observed and of the bootstrap panels)
  eng <- NULL; B <- max(19L, min(as.integer(A$n_perm), 99L))
  if (fam_reg) {
    eng <- .plapply(seq_along(key), function(j) {
      lags <- if (first) scan_lags else if (!is.na(lag_pair[j])) lag_pair[j] else max(1L, as.integer(A$lag))
      tryCatch(.reg_engine(U, grid$effect[j], grid$cause[j], ctrl_of(j), lags, aug, B), error = function(e) NULL)
    }, cores = A$cores)
  }

  # ---- units: the method in every unit that is long enough
  unit_res <- list(); unit_err <- character(0); units_tab <- NULL; nulls <- list(); calib <- NULL; unit_ids <- character(0)
  fast_units <- fam_reg && (!bayes || ctl$prior %in% c("bic", "g", "hyper_g"))     # unit tests straight from the engine
  if ("units" %in% pc$levels && any(long) && fast_units) {
    rows <- lapply(seq_along(key), function(j) {
      if (is.null(eng[[j]])) return(NULL)
      r <- .unit_reg_rows(eng[[j]], if (first) "first" else "fixed", A, ctl$prior, bayes)
      pp <- r$bf10 * A$prior_prob / (r$bf10 * A$prior_prob + 1 - A$prior_prob)
      nom <- if (bayes) { if (identical(ctl$decision, "posterior")) pp > A$prob_threshold else r$bf10 > A$bf_threshold } else r$p_value < A$alpha
      data.frame(unit = names(U), cause = grid$cause[j], effect = grid$effect[j], n = PD$n, lag = as.integer(r$lag), statistic = r$statistic,
                 df1 = r$df1, df2 = r$df2, p_value = r$p_value, bf10 = r$bf10, post_prob = pp, effect_size = r$effect_size, sign = r$sign,
                 nominal = nom %in% TRUE, stringsAsFactors = FALSE)
    })
    units_tab <- do.call(rbind, rows)
    if (!is.null(units_tab)) {
      done <- tapply(is.finite(units_tab$statistic), units_tab$unit, any)
      unit_ids <- names(U)[long & names(U) %in% names(done)[done]]
      units_tab <- units_tab[units_tab$unit %in% unit_ids, , drop = FALSE]
      units_tab <- units_tab[order(match(units_tab$unit, unit_ids), match(paste(units_tab$cause, "->", units_tab$effect), key)), ]
      rownames(units_tab) <- NULL; if (!nrow(units_tab)) units_tab <- NULL
    }
  } else if ("units" %in% pc$levels && any(long)) {
    Au <- A; Au$id <- Au$time <- Au$target <- NULL; Au$adjust <- "none"; Au["seed"] <- list(NULL)
    Au$method <- method; Au$bayes <- bayes
    if (first) Au$max_lag <- max(scan_lags)
    sc <- .split_cores(sum(long), A$cores); Au$cores <- sc$inner
    rr <- .plapply(names(U)[long], function(u) {
      d <- U[[u]]
      tryCatch(lapply(lag_set, function(L) {
        if (!is.na(L)) Au$lag <- L
        suppressWarnings(suppressMessages(
          .rcausality_ts(d[, xvars, drop = FALSE], if (is.null(yvars)) NULL else d[, yvars, drop = FALSE], Au)))
      }), error = function(e) structure(conditionMessage(e), class = "unit_error"))
    }, cores = sc$outer)
    names(rr) <- names(U)[long]
    failed <- vapply(rr, inherits, TRUE, "unit_error")
    unit_err <- vapply(rr[failed], as.character, ""); unit_res <- rr[!failed]; unit_ids <- names(unit_res)
    grp <- if (anyNA(lag_pair)) rep(1L, nrow(grid)) else match(lag_pair, lag_set)   # the run that holds each link
    if (length(unit_res)) {
      units_tab <- do.call(rbind, lapply(unit_ids, function(u) {
        do.call(rbind, lapply(seq_along(lag_set), function(g) {
          jj <- which(grp == g); if (!length(jj)) return(NULL)
          P <- unit_res[[u]][[g]]$pairs; k <- match(key[jj], paste(P$cause, "->", P$effect))
          data.frame(unit = u, cause = grid$cause[jj], effect = grid$effect[jj], n = PD$n[[u]], lag = P$lag[k], statistic = P$statistic[k],
                     df1 = P$df1[k], df2 = P$df2[k], p_value = P$p_value[k], bf10 = P$bf10[k], post_prob = P$post_prob[k],
                     effect_size = P$effect_size[k], sign = P$sign[k], nominal = P$significant[k] %in% TRUE, stringsAsFactors = FALSE)
        }))
      }))
      units_tab <- units_tab[order(match(units_tab$unit, unit_ids), match(paste(units_tab$cause, "->", units_tab$effect), key)), ]
      rownames(units_tab) <- NULL
      nulls <- lapply(unit_res, function(r) stats::setNames(lapply(seq_along(key), function(j) r[[grp[j]]]$null[[key[j]]]), key))
    }
  }

  # ---- group level and the unit-level error control
  panel_rows <- list(); prev_rows <- list()
  comb <- if (pc$combine == "auto") .group_default[[method]] else pc$combine
  if (comb == "dh" && !fam_reg) stop("combine = \"dh\" (Dumitrescu-Hurlin) is defined for the Granger family (\"granger\", \"toda_yamamoto\").", call. = FALSE)
  if (!is.null(units_tab)) {
    units_tab$post_pooled <- NA_real_
    do_group <- "group" %in% pc$levels && length(unit_ids) >= pc$min_units
    rows_of <- function(j) which(units_tab$cause == grid$cause[j] & units_tab$effect == grid$effect[j])
    # Granger family: bootstrap null of the unit statistics, and the empirical null of the unit tests.  In short
    # series the F-test of a unit is not exact (persistent regressors); the bootstrap panels show the null
    # distribution of the unit z-scores of a link, and the unit p-values are recomputed against the normal
    # distribution with the location and scale of that empirical null (Efron 2004), never below zero / one: the
    # correction can only make a unit test more conservative.  The unit Bayes factors are lowered by the
    # corresponding change of the likelihood ratio.
    reg_nulls <- NULL
    if (fam_reg) {
      ui <- match(unit_ids, names(U))
      reg_nulls <- lapply(seq_along(key), function(j) {
        if (is.null(eng[[j]])) return(vector("list", length(ui)))
        ul <- if (first) NULL else { v <- rep(NA_integer_, length(U)); v[ui] <- units_tab$lag[rows_of(j)]; v }
        .unit_reg_nulls(eng[[j]], if (first) "first" else "fixed", ul)[ui]
      })
      units_tab$p_ftest <- units_tab$p_value; units_tab$bf10_raw <- units_tab$bf10
      calib <- do.call(rbind, lapply(seq_along(key), function(j) {
        en <- .empirical_null(unlist(lapply(reg_nulls[[j]], `[[`, "p0")))
        data.frame(cause = grid$cause[j], effect = grid$effect[j], n_null = en[["n"]], mu0 = en[["mu0"]], sigma0 = en[["sigma0"]])
      }))
      for (j in seq_along(key)) if (calib$mu0[j] > 0 || calib$sigma0[j] > 1) {
        rows <- rows_of(j); ut <- units_tab[rows, ]
        pc_ <- .p_calibrate(ut$p_ftest, c(mu0 = calib$mu0[j], sigma0 = calib$sigma0[j])); pc_[is.na(ut$p_ftest)] <- NA_real_
        Fn <- stats::qf(ut$p_ftest, ut$df1, ut$df2, lower.tail = FALSE); Fc <- stats::qf(pc_, ut$df1, ut$df2, lower.tail = FALSE)
        dl <- (ut$df2 + ut$df1) / 2 * (log1p(ut$df1 * Fn / ut$df2) - log1p(ut$df1 * Fc / ut$df2)); dl[!is.finite(dl)] <- 0
        units_tab$p_value[rows] <- pc_
        units_tab$bf10[rows] <- ut$bf10 * exp(-pmax(dl, 0))
        units_tab$post_prob[rows] <- units_tab$bf10[rows] * A$prior_prob / (units_tab$bf10[rows] * A$prior_prob + 1 - A$prior_prob)
        units_tab$nominal[rows] <- if (use_bf) { if (identical(ctl$decision, "posterior")) units_tab$post_prob[rows] > A$prob_threshold
                                                 else units_tab$bf10[rows] > A$bf_threshold } else pc_ < A$alpha
        units_tab$nominal[rows][is.na(units_tab$nominal[rows])] <- FALSE
      }
    }
    for (j in seq_along(key)) {
      rows <- rows_of(j); ut <- units_tab[rows, ]
      # a unit whose side condition fails (no convergence, no predictive gain) gives no evidence for the link
      bf_use <- ifelse(!ut$nominal & is.finite(ut$bf10) & ut$bf10 > A$bf_threshold & ctl$decision == "bf", 1, ut$bf10)
      pv <- .prevalence_bf(bf_use, pc$prevalence_prior, ctl$credible)
      units_tab$post_pooled[rows] <- pv$unit_prob
      prev_rows[[j]] <- data.frame(cause = grid$cause[j], effect = grid$effect[j], n_units = nrow(ut), n_nominal = sum(ut$nominal),
                                   share_nominal = mean(ut$nominal), prevalence_bf10 = pv$bf10, prevalence = pv$mean,
                                   prevalence_lower = pv$lower, prevalence_upper = pv$upper, stringsAsFactors = FALSE)
      if (!do_group) next
      g <- list(statistic = NA_real_, p_value = NA_real_, n_units = nrow(ut)); how <- comb; detail <- ""
      if (fam_reg) {
        d <- .dh_test(ut$df1 * ut$statistic, ut$df1, ut$df2)
        if (how == "dh") {
          g <- list(statistic = d[["Z_tilde"]], p_value = d[["p_value"]], n_units = as.integer(d[["n_units"]]))
          detail <- sprintf("Dumitrescu-Hurlin, W-bar = %.2f, Z-bar = %.2f, asymptotic p-value", d[["W_bar"]], d[["Z_bar"]])
        }
      }
      if (how == "mean") {
        nl <- if (fam_reg) reg_nulls[[j]] else lapply(nulls[ut$unit], function(e) e[[key[j]]])
        g <- .group_mean_test(nl, pc$n_draws)
        detail <- if (fam_reg) sprintf("group mean of the unit statistics against their bootstrap null (W-bar = %.2f; Dumitrescu-Hurlin Z-tilde = %.2f, asymptotic p = %s)",
                                       d[["W_bar"]], d[["Z_tilde"]], format.pval(d[["p_value"]], digits = 2, eps = 1e-4))
                  else "group mean of the standardised unit statistics against the unit surrogates"
        if (is.na(g$p_value) && any(is.finite(ut$p_value))) how <- "stouffer"
      }
      if (how %in% c("stouffer", "fisher", "simes")) {
        g <- .combine_p(ut$p_value, how)
        detail <- c(stouffer = "Stouffer combination of the unit p-values", fisher = "Fisher combination of the unit p-values",
                    simes = "Simes test of the unit p-values")[[how]]
      } else if (how == "count") {
        k <- sum(ut$nominal); N <- nrow(ut)
        g <- list(statistic = k / N, p_value = stats::pbinom(k - 1, N, A$alpha, lower.tail = FALSE), n_units = N)
        detail <- sprintf("%d of %d units with the link, binomial test against the rate %g", k, N, A$alpha)
      }
      # Bayesian variant: the group statistic against its null distribution (as for the surrogate methods of
      # rcausality()) when the group test is the calibrated mean; otherwise the prevalence model
      bf_g <- pv$bf10; bf_how <- "prevalence model on the unit Bayes factors"
      if (bayes && how == "mean" && is.finite(g$statistic) && length(g$null) > 100L) {
        b <- tryCatch(surrogate_bf(g$statistic, g$null[seq_len(min(length(g$null), 1000L))], alternative = "greater",
                                   effect_prior = ctl$effect_prior, scale = ctl$scale, null_model = ctl$null_model,
                                   prior_prob = A$prior_prob, n_draws = min(ctl$n_draws, 2000L), credible = ctl$credible), error = function(e) NULL)
        if (!is.null(b) && is.finite(b$bf10)) { bf_g <- b$bf10; bf_how <- "Bayes factor of the group statistic against its null distribution" }
      }
      if (use_bf) detail <- sprintf("%s; %s", bf_how, detail)
      sig_units <- ut[ut$nominal, , drop = FALSE]
      lag_g <- suppressWarnings(as.integer(round(stats::median(if (nrow(sig_units)) sig_units$lag else ut$lag, na.rm = TRUE))))
      sg <- suppressWarnings(sign(sum(if (nrow(sig_units)) sig_units$sign else ut$sign, na.rm = TRUE)))
      panel_rows[[length(panel_rows) + 1L]] <- data.frame(
        cause = grid$cause[j], effect = grid$effect[j], level = "group", test = how, lag = lag_g,
        statistic = g$statistic, df1 = NA_real_, df2 = NA_real_, p_value = if (use_bf) NA_real_ else g$p_value, p_freq = g$p_value,
        bf10 = bf_g, effect_size = mean(ut$effect_size, na.rm = TRUE), r2_gain = NA_real_,
        sign = if (all(is.na(ut$sign))) NA_real_ else sg, n_eff = g$n_units, gate = TRUE, detail = detail, stringsAsFactors = FALSE)
    }
    # error control over all unit x link cells
    if (use_bf) {
      pp <- ifelse(is.finite(units_tab$post_pooled), units_tab$post_pooled, units_tab$post_prob)
      units_tab$q_value <- NA_real_
      units_tab$significant <- .bayes_fdr(pp, A$alpha) & units_tab$nominal
    } else {
      units_tab$q_value <- stats::p.adjust(units_tab$p_value, method = pc$unit_adjust)
      units_tab$significant <- units_tab$nominal & !is.na(units_tab$q_value) & units_tab$q_value < A$alpha
      if (all(is.na(units_tab$p_value))) units_tab$significant <- units_tab$nominal
    }
  }

  # ---- pooled level
  pooled_by_lag <- NULL
  if ("pooled" %in% pc$levels && isTRUE(info$pooled)) {
    if (fam_reg) {
      pr <- lapply(seq_along(key), function(j) if (is.null(eng[[j]])) NULL else .pooled_reg_result(eng[[j]], A, if (first) "first" else "fixed"))
      if (all(vapply(pr, is.null, TRUE))) unit_err <- c(unit_err, pooled = "the units are too short for the pooled regression")
    } else {
      Ap <- A; Ap$method <- method; Ap$bayes <- bayes
      sc <- .split_cores(nrow(grid), A$cores)
      pr <- .plapply(seq_len(nrow(grid)), function(j) {
        if (!is.na(lag_pair[j])) Ap$lag <- lag_pair[j]
        tryCatch(.pooled_test(method, U, grid$effect[j], grid$cause[j], ctrl_of(j), Ap, sc$inner),
                 error = function(e) structure(conditionMessage(e), class = "unit_error"))
      }, cores = sc$outer)
    }
    for (j in seq_along(key)) {
      r <- pr[[j]]; if (is.null(r)) next
      if (inherits(r, "unit_error")) { unit_err <- c(unit_err, pooled = as.character(r)); next }
      row <- do.call(.pair_row, c(list(cause = grid$cause[j], effect = grid$effect[j]), r$row))
      if (bayes && !is.null(r$null)) {
        b <- surrogate_bf(r$null$statistic, r$null$null, alternative = r$null$alternative, effect_prior = ctl$effect_prior,
                          scale = ctl$scale, null_model = ctl$null_model, prior_prob = A$prior_prob, n_draws = ctl$n_draws,
                          credible = ctl$credible)
        row$bf10 <- b$bf10
      }
      panel_rows[[length(panel_rows) + 1L]] <- data.frame(
        cause = row$cause, effect = row$effect, level = "pooled", test = method, lag = row$lag, statistic = row$statistic,
        df1 = row$df1, df2 = row$df2, p_value = if (bayes) NA_real_ else row$p_value, p_freq = row$p_value, bf10 = row$bf10,
        effect_size = row$effect_size, r2_gain = row$r2_gain, sign = row$sign, n_eff = row$n_eff, gate = isTRUE(r$gate),
        detail = r$detail, stringsAsFactors = FALSE)
      if (!is.null(r$by_lag)) pooled_by_lag <- rbind(pooled_by_lag, data.frame(cause = row$cause, effect = row$effect, r$by_lag))
    }
  }
  if (!length(panel_rows)) {
    msg <- if (length(unit_err)) paste0(" First error: ", unit_err[[1]]) else ""
    stop(sprintf(paste0("No panel-level test could be computed for method \"%s\": fewer than %d units are long enough for a ",
                        "unit-level analysis (advisable length %d; see panel_control(min_length = , min_units = )) and the ",
                        "method has no pooled estimator.%s"), method, pc$min_units, min_len, msg), call. = FALSE)
  }
  PT <- do.call(rbind, panel_rows); rownames(PT) <- NULL
  PT$p_adj <- NA_real_; PT$significant <- FALSE
  fam <- if (is.null(fam_fun)) rep("all", nrow(PT)) else fam_fun(PT$cause, PT$effect)
  for (lv in unique(paste(PT$level, fam))) {
    k <- paste(PT$level, fam) == lv
    if (any(!is.na(PT$p_value[k]))) PT$p_adj[k] <- if (adjust == "none") PT$p_value[k] else stats::p.adjust(PT$p_value[k], method = adjust)
    PT$significant[k] <- if (use_bf) (!is.na(PT$bf10[k]) & PT$bf10[k] > A$bf_threshold & PT$gate[k])
                         else (!is.na(PT$p_adj[k]) & PT$p_adj[k] < A$alpha & PT$gate[k])
  }
  # ---- the decision of the panel: one level, or both
  have <- unique(PT$level); dec <- pc$decision
  if (dec == "auto") dec <- if (all(c("group", "pooled") %in% have)) "both" else have[1L]
  if (dec == "both" && !all(c("group", "pooled") %in% have)) dec <- have[1L]
  if (!dec %in% c(have, "both")) {
    warning(sprintf("The %s level is not available for this analysis; the %s level decides.", dec, have[1L]), call. = FALSE)
    dec <- have[1L]
  }
  pick <- function(lv) { D <- PT[PT$level == lv, ]; D[match(key, paste(D$cause, "->", D$effect)), ] }
  if (dec == "both") {
    # either level may show the link: the smaller of the two p-values, doubled (Bonferroni over the two tests);
    # Bayesian variant: the average of the two Bayes factors, the Bayes factor of "a link at one of the levels"
    G <- pick("group"); Pp <- pick("pooled"); D <- Pp
    pm <- cbind(G$p_value, Pp$p_value); bm <- cbind(G$bf10, Pp$bf10)
    D$p_value <- apply(pm, 1L, function(p) { p <- p[is.finite(p)]; if (!length(p)) NA_real_ else min(1, min(p) * length(p)) })
    D$bf10 <- apply(bm, 1L, function(b) { b <- b[is.finite(b)]; if (!length(b)) NA_real_ else mean(b) })
    lead <- ifelse(is.na(Pp$p_freq) | (!is.na(G$p_freq) & G$p_freq < Pp$p_freq), 1L, 2L)      # the level with the stronger evidence
    if (use_bf) lead <- ifelse(is.na(Pp$bf10) | (!is.na(G$bf10) & G$bf10 > Pp$bf10), 1L, 2L)
    for (cn in c("lag", "statistic", "df1", "df2", "effect_size", "r2_gain", "sign", "n_eff")) D[[cn]] <- ifelse(lead == 1L, G[[cn]], Pp[[cn]])
    D$lag <- ifelse(is.na(Pp$lag), G$lag, Pp$lag); D$sign <- ifelse(is.na(Pp$sign), G$sign, Pp$sign)
    D$effect_size <- ifelse(is.na(Pp$effect_size), G$effect_size, Pp$effect_size); D$r2_gain <- Pp$r2_gain
    D$gate <- ifelse(lead == 1L, G$gate, Pp$gate) %in% TRUE
    D$p_adj <- NA_real_
    famD <- if (is.null(fam_fun)) rep("all", nrow(D)) else fam_fun(D$cause, D$effect)
    for (f in unique(famD)) { k <- famD == f
      if (any(!is.na(D$p_value[k]))) D$p_adj[k] <- if (adjust == "none") D$p_value[k] else stats::p.adjust(D$p_value[k], method = adjust) }
    D$significant <- if (use_bf) (!is.na(D$bf10) & D$bf10 > A$bf_threshold & D$gate) else (!is.na(D$p_adj) & D$p_adj < A$alpha & D$gate)
  } else D <- pick(dec)
  pairs <- data.frame(cause = grid$cause, effect = grid$effect, lag = as.integer(D$lag), statistic = D$statistic, df1 = D$df1,
                      df2 = D$df2, p_value = D$p_value, p_adj = D$p_adj, bf10 = D$bf10,
                      post_prob = ifelse(is.na(D$bf10), NA_real_, D$bf10 * A$prior_prob / (D$bf10 * A$prior_prob + 1 - A$prior_prob)),
                      evidence = as.character(bf_evidence(D$bf10)), effect_size = D$effect_size, r2_gain = D$r2_gain,
                      sign = D$sign, n_eff = as.integer(D$n_eff), significant = D$significant %in% TRUE, stringsAsFactors = FALSE)
  m_out <- if (bayes && method == "granger") "bayes_granger" else method
  out <- .new_causres(pairs, series = vars, causes = causes, effects = effects, method = m_out,
                      conditional = isTRUE(A$conditional), alpha = A$alpha, bf_threshold = A$bf_threshold, adjust = adjust,
                      n = sum(PD$n), call = cl, data = PD$data[, c(A$id, A$time, vars), drop = FALSE], by_lag = pooled_by_lag,
                      prob_threshold = A$prob_threshold)
  out$units <- units_tab; out$panel <- PT
  out$prevalence <- if (length(prev_rows)) do.call(rbind, prev_rows) else NULL
  out$id <- A$id; out$time <- A$time; out$unit_errors <- unit_err
  out$settings$bayes <- bayes; out$settings$base_method <- method
  if (bayes) { out$settings$prior <- ctl$prior; out$settings$bayes_control <- ctl; out$settings$decision <- ctl$decision }
  out$settings$panel <- list(n_units = length(U), n_range = range(PD$n), units_analysed = length(unit_ids),
                             units_short = sum(!long), units_failed = sum(names(unit_err) != "pooled"), decision = dec,
                             levels = have, group_test = if ("group" %in% have) unique(PT$test[PT$level == "group"]) else NA_character_,
                             unit_adjust = pc$unit_adjust, use_bf = use_bf, min_length = min_len, calibration = calib,
                             lag_rule = if (!anyNA(lag_pair)) "pooled BIC" else if (first) "first" else NA_character_, control = pc)
  A$panel <- pc; A$method <- if (bayes && method == "granger" && ctl$prior %in% c("spike_slab", "normal", "minnesota")) "granger" else m_out
  A$bayes <- bayes; A[["prior"]] <- ctl$prior
  out$args <- .flat_args(A)
  class(out) <- c("causpanel", "causres")
  out
}

#' Causality analysis of panel data
#'
#' The panel version of [rcausality()]: the same `method`, `bayes` and
#' `cores` options for a long data frame with a unit column. Every directed
#' pair of the variables is analysed at up to three levels (see
#' [panel_control()]):
#' * **units**: the method in every unit that is long enough for it, with
#'   false-discovery-rate control over all unit x link cells
#'   (Benjamini-Hochberg; for the Bayesian variants the Bayesian FDR of
#'   Newton et al. 2004 on the partially pooled posterior probabilities);
#' * **group**: a test of "no unit has the link" from the unit results, by
#'   default the mean of the unit statistics, each standardised with its own
#'   null sample, against the convolution of the unit null distributions.
#'   The Bayesian variant is the Bayes factor of that group statistic against
#'   its null distribution (as for the surrogate methods of [rcausality()]);
#'   in addition, a prevalence model on the unit Bayes factors `BF_i` (the
#'   share of the units with the link has a Beta prior; the Bayes factor of
#'   "share > 0" is the prior mean of `prod_i (1 - share + share * BF_i)`)
#'   gives the posterior of the share and partially pooled posterior
#'   probabilities of the link in every unit;
#' * **pooled**: one estimate from all units, for the methods with a pooled
#'   estimator.
#'
#' **Granger family** (`"granger"`, `"bayes_granger"`, `"toda_yamamoto"`).
#' The pooled test is that of Juodis, Karavias and Sarafidis (2021): the
#' cause's lags have common coefficients (under the null hypothesis they are
#' zero in every unit), while the intercept, the effect's own lags and the
#' controls' lags have unit-specific coefficients; the bias of order `1/T`
#' of the pooled estimator (Nickell 1981), which makes the test of the
#' reverse direction of a true link over-reject, is removed by the
#' half-panel jackknife of Dhaene and Jochmans (2015). In short series the F
#' distribution is only an approximation, so the unit statistics and the
#' pooled statistic are referred to a bootstrap under the null hypothesis
#' (in every unit the effect is regenerated from its own estimated
#' autoregression with resampled residuals, the cause is kept; the
#' small-sample bias of the estimated persistence is corrected first, Kendall
#' 1954, Kilian 1998): the
#' bootstrap sample is the null sample of the group test, and it gives the
#' empirical null (Efron 2004) with which the unit p-values and the pooled
#' p-value are calibrated -- the calibration can only make a test more
#' conservative. The Bayes factors use the prior of [bayes_control()]: the
#' unit Bayes factors as in [rcausality()], the pooled Bayes factor from the
#' pooled estimate and its covariance (BIC, g, hyper-g, normal /
#' spike-and-slab or Minnesota prior on the pooled coefficients); both are
#' lowered by the change of the likelihood ratio that the calibration
#' implies. With `lag = "auto"` the order of every link is chosen by the BIC
#' of its pooled regression and used in all units (choosing the order unit
#' by unit would distort the unit tests); `lag = "first"` scans the orders
#' with a Bonferroni correction.
#'
#' **Model-free methods.** The pooled estimators use the stacked delay
#' vectors of the units, standardised within the units: the random forest
#' (`"rf_granger"`), the Diks-Panchenko statistic, the kNN conditional
#' mutual information (`"cmi"`, `"cmi_lag"`), the binned transfer entropy
#' (`"te"`) and multispatial convergent cross mapping (`"ccm"`, Clark et al.
#' 2015). Their null distribution comes from surrogate panels in which the
#' cause is delayed circularly within every unit; with `bayes = TRUE` the
#' statistic is compared with its surrogates by [surrogate_bf()]. These
#' estimators use all `units x periods` delay vectors, so they work with
#' units that are far too short for a unit-by-unit analysis (see
#' `min_length` in [causal_methods()]); within short units the shifted
#' cause is not independent of the original one, which makes these tests
#' conservative.
#'
#' The result is a [`causres`][causres-class] object whose `pairs` and
#' `classes` hold the panel-level decision, so that `plot()`,
#' [placebo_test()], [robustness()], [compare_methods()] and `predict()`
#' work as for time series, with the unit and panel tables in addition.
#' The units are treated as independent; see `time_effects` in
#' [panel_control()] for common shocks.
#'
#' @param data a data frame in long format (one row per unit and period); see
#'   [panel_long()] for matrices and arrays.
#' @param id name of the unit column.
#' @param time optional name of the time column (the rows are ordered by it
#'   within the units; consecutive periods are assumed, missing values are
#'   interpolated within the unit).
#' @param method,bayes,cores the causality method, its Bayesian variant and
#'   the number of CPU cores, as in [rcausality()] (units, pairs, surrogates
#'   and chains are distributed over the cores; the results do not depend on
#'   the number of cores).
#' @param target optional name of one explained variable (see
#'   [causal_drivers()]).
#' @param panel a [panel_control()] list.
#' @param ... further arguments of [rcausality()] (`lag`, `max_lag`,
#'   `conditional`, `alpha`, `adjust`, `control`, `prior`, `n_perm`, `seed`,
#'   `method_args`, ...). `n_perm` is also the number of bootstrap panels of
#'   the Granger family (at most 99).
#' @return An object of class `causpanel` (inheriting from
#'   [`causres`][causres-class]) with, in addition to the slots of a
#'   `causres`: `units` (one row per unit and directed pair: lag, statistic,
#'   p-value, Bayes factor, `nominal` decision of the unit, `q_value`,
#'   `post_pooled` and the error-controlled decision `significant`; for the
#'   Granger family also the uncalibrated `p_ftest` and `bf10_raw`), `panel`
#'   (the group- and pooled-level tests, with a `detail` text), `prevalence`
#'   (share of units with the link: observed, and the posterior under the
#'   prevalence model), `id`, `time`. Methods: `print`, `summary`, `plot`
#'   (see [plot.causpanel()]), `as.matrix`, `predict`.
#' @references
#' Juodis, A., Karavias, Y. and Sarafidis, V. (2021). A homogeneous approach to testing for Granger non-causality in heterogeneous panels. *Empirical Economics*, 60(1), 93-112. \doi{10.1007/s00181-020-01970-9}
#'
#' Dhaene, G. and Jochmans, K. (2015). Split-panel jackknife estimation of fixed-effect models. *The Review of Economic Studies*, 82(3), 991-1030. \doi{10.1093/restud/rdv007}
#'
#' Dumitrescu, E.-I. and Hurlin, C. (2012). Testing for Granger non-causality in heterogeneous panels. *Economic Modelling*, 29(4), 1450-1460. \doi{10.1016/j.econmod.2012.02.014}
#'
#' Kosztyan, Z. T., Matrai, R. and Kiss, D. (2026). Is the world getting smaller? Network-driven heterogeneity in gravity-model elasticities. *The World Economy*. \doi{10.1111/twec.70148}
#'
#' Efron, B. (2004). Large-scale simultaneous hypothesis testing: the choice of a null hypothesis. *Journal of the American Statistical Association*, 99(465), 96-104. \doi{10.1198/016214504000000089}
#'
#' Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society B*, 57(1), 289-300. \doi{10.1111/j.2517-6161.1995.tb02031.x}
#'
#' Newton, M. A., Noueiry, A., Sarkar, D. and Ahlquist, P. (2004). Detecting differential gene expression with a semiparametric hierarchical mixture method. *Biostatistics*, 5(2), 155-176. \doi{10.1093/biostatistics/5.2.155}
#'
#' Ince, R. A. A., Paton, A. T., Kay, J. W. and Schyns, P. G. (2021). Bayesian inference of population prevalence. *eLife*, 10, e62461. \doi{10.7554/eLife.62461}
#'
#' Pesarin, F. and Salmaso, L. (2010). *Permutation Tests for Complex Data: Theory, Applications and Software*. Wiley, Chichester. \doi{10.1002/9780470689516}
#'
#' Clark, A. T., Ye, H., Isbell, F., Deyle, E. R., Cowles, J., Tilman, G. D. and Sugihara, G. (2015). Spatial convergent cross mapping to detect causal relationships from short time series. *Ecology*, 96(5), 1174-1181. \doi{10.1890/14-1479.1}
#'
#' Nickell, S. (1981). Biases in dynamic models with fixed effects. *Econometrica*, 49(6), 1417-1426. \doi{10.2307/1911408}
#'
#' Kendall, M. G. (1954). Note on bias in the estimation of autocorrelation. *Biometrika*, 41(3-4), 403-404. \doi{10.1093/biomet/41.3-4.403}
#'
#' Kilian, L. (1998). Small-sample confidence intervals for impulse response functions. *The Review of Economics and Statistics*, 80(2), 218-230. \doi{10.1162/003465398557465}
#' @seealso [rcausality()], [panel_control()], [plot.causpanel()],
#'   [causal_drivers()], [sim_panel()], [panel_long()]
#' @export
#' @examples
#' # 30 units observed for 20 periods; x drives y in every unit with a
#' # unit-specific strength.  A single unit is too short for a reliable test;
#' # the group-level and the pooled tests find x -> y.
#' pd <- sim_panel(n_units = 30, n_time = 20, scenario = "x_to_y", seed = 1)
#' pg <- panel_causality(pd, id = "unit", time = "time", seed = 1)
#' pg
#' \donttest{
#' plot(pg)                       # heat map: links x units, with the panel tests
#' plot(pg, "prevalence")
#' # the same call through rcausality(), Bayesian variant (closed-form Bayes factors)
#' rcausality(pd, id = "unit", time = "time", bayes = TRUE, prior = "bic", seed = 1)
#' # a nonlinear coupling in a short panel: the pooled conditional mutual information
#' pn <- sim_panel(n_units = 30, n_time = 20, scenario = "x_to_y", nonlinear = TRUE,
#'                 coupling = 0.8, seed = 2)
#' rcausality(pn, id = "unit", time = "time", method = "cmi", n_perm = 100, seed = 1, cores = 2)
#' }
panel_causality <- function(data, id, time = NULL, method = "granger", bayes = FALSE,
                            cores = getOption("Rcausality.cores", 1L), target = NULL, panel = panel_control(), ...) {
  rcausality(data, method = method, bayes = bayes, cores = cores, id = id, time = time, target = target, panel = panel, ...)
}

#' @rdname panel_causality
#' @details `panel_granger()` is `panel_causality()` with `method =
#'   "granger"` (kept for the scripts written for version 0.0.2).
#' @export
panel_granger <- function(data, id, time = NULL, ...) panel_causality(data, id = id, time = time, method = "granger", ...)

#' Simulate panels of short series
#'
#' Every unit follows a VAR(1) with unit-specific coupling strengths
#' `b_i ~ N(coupling, (heterogeneity * coupling)^2)` and unit fixed effects.
#' * two series (`x`, `y`): `"x_to_y"`, `"bidirectional"`, `"independent"`;
#' * `scenario = "drivers"`: five series for the analysis of one explained
#'   variable `y` -- `x1 -> y` at lag 1 and `x3 -> y` at lag 2 (direct
#'   drivers), `x2 -> x1` (so `x2` drives `y` only through `x1`), `x4`
#'   unrelated.
#' With `prevalence < 1` only that share of the units has the links (a
#' heterogeneous panel); `nonlinear = TRUE` couples through the centred
#' square of the cause, which is invisible to linear tests.
#'
#' @param n_units number of units.
#' @param n_time number of periods per unit (a vector gives an unbalanced
#'   panel).
#' @param scenario `"x_to_y"`, `"bidirectional"`, `"independent"` or
#'   `"drivers"`.
#' @param coupling mean coupling strength; `heterogeneity` its relative
#'   standard deviation across the units.
#' @param heterogeneity relative standard deviation of the coupling.
#' @param prevalence share of the units in which the links exist.
#' @param nonlinear logical; quadratic instead of linear coupling.
#' @param ar autoregressive coefficient.
#' @param seed optional RNG seed.
#' @return A long data frame with the columns `unit`, `time` and the series;
#'   attribute `"linked"` names the units with the links, `"truth"` the
#'   scenario.
#' @export
#' @examples
#' head(sim_panel(n_units = 3, n_time = 5, seed = 1))
#' str(sim_panel(n_units = 4, n_time = 6, scenario = "drivers", seed = 1))
sim_panel <- function(n_units = 30L, n_time = 20L, scenario = c("x_to_y", "bidirectional", "independent", "drivers"),
                      coupling = 0.4, heterogeneity = 0.5, prevalence = 1, nonlinear = FALSE, ar = 0.5, seed = NULL) {
  scenario <- match.arg(scenario); .local_seed(seed)
  .check_number(prevalence, "prevalence", lower = 0, upper = 1)
  n_time <- rep_len(as.integer(n_time), n_units)
  linked <- seq_len(n_units) %in% sample.int(n_units, round(prevalence * n_units))
  f <- if (nonlinear) function(v) (v^2 - 1) / sqrt(2) else identity
  ids <- sprintf(paste0("u%0", nchar(n_units), "d"), seq_len(n_units))
  out <- lapply(seq_len(n_units), function(i) {
    b <- if (linked[i] && scenario != "independent") stats::rnorm(1, coupling, heterogeneity * abs(coupling)) else 0
    burn <- 50L; Tt <- n_time[i] + burn
    if (scenario == "drivers") {
      X <- .sim_drivers_unit(n_time[i], b, nonlinear, ar, burn)
    } else {
      bxy <- b; byx <- if (scenario == "bidirectional") b else 0
      X <- matrix(0, Tt, 2L, dimnames = list(NULL, c("x", "y"))); E <- matrix(stats::rnorm(Tt * 2L), Tt, 2L)
      for (t in 2:Tt) {
        X[t, "x"] <- ar * X[t - 1, "x"] + byx * f(X[t - 1, "y"]) + E[t, 1]
        X[t, "y"] <- ar * X[t - 1, "y"] + bxy * f(X[t - 1, "x"]) + E[t, 2]
      }
      X <- X[-seq_len(burn), , drop = FALSE]
    }
    X <- sweep(X, 2, stats::rnorm(ncol(X)), "+")                         # unit fixed effects
    data.frame(unit = ids[i], time = seq_len(n_time[i]), X, stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, out); rownames(d) <- NULL
  attr(d, "truth") <- scenario; attr(d, "linked") <- ids[linked]
  d
}
