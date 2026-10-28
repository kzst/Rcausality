# Recurrence-plot based detection of hidden common causes and of directional
# coupling (Hirata & Aihara, Phys. Rev. E 81, 016203, 2010) and the
# correlation-dimension method of Krakovska (Entropy 21, 818, 2019).

#' Recurrence-plot causality and hidden common cause (Hirata and Aihara, 2010)
#'
#' Two criteria on the recurrence plots `R_X`, `R_Y` of the embedded series:
#' (1) *simultaneous recurrences*: related series recur together far more
#' often than independent ones, measured by the joint recurrence rate
#' relative to the product of the single rates (`excess`), tested against
#' time-shifted surrogates; (2) *coverage*: if `X` drives `Y`, the state of
#' `Y` determines the state of `X` (forced Takens theorem), so every
#' recurrence of `Y` is also a recurrence of `X` -- the fraction of `Y`'s
#' (tight, rate `rate / 5`) recurrences covered by `X`'s (loose, rate
#' `rate`) recurrences, `P(R_X = 1 | R_Y = 1)`, is close to 1 -- whereas
#' non-coverage denies the directed coupling (with equal rates the two
#' conditional probabilities would coincide, hence the tight/loose pair).
#' The normalised coverage index `(coverage - rate) / (1 - rate)` is 0 for
#' unrelated series and 1 for perfect coverage. Verdict: a direction whose index exceeds
#' `theta` and is significant against the surrogates; both directions =
#' bidirectional; significant excess of simultaneous recurrences without a
#' covered direction = hidden common cause; otherwise independent.
#'
#' @param x,y numeric vectors (for `plot`, `x` is a `hirata` object).
#' @param E,tau embedding dimension and delay.
#' @param rate target recurrence rate (the thresholds are the corresponding
#'   quantiles of the pairwise distances).
#' @param theiler Theiler window (pairs closer in time are ignored).
#' @param n_surr number of circular-shift surrogates.
#' @param theta coverage-index threshold for a directed link.
#' @param alpha significance level.
#' @param seed optional RNG seed.
#' @param cores CPU cores for the surrogates (see [Rcausality-parallel]).
#' @return An object of class `hirata`: `coverage` (`x_to_y` =
#'   `P(R_X | R_Y)`, `y_to_x`), `index` (normalised coverages), `excess`
#'   (joint rate over the product of the rates), `p_values` (surrogate
#'   p-values of the two coverages and of the excess), `verdict`,
#'   `embedding` (the embedded series and thresholds). Methods: `print`,
#'   `plot` (the two recurrence plots and the coverage picture).
#' @references
#' Hirata, Y. and Aihara, K. (2010). Identifying hidden common causes from bivariate time series by using the recurrence plot. *Physical Review E*, 81, 016203. \doi{10.1103/PhysRevE.81.016203}
#'
#' Marwan, N., Romano, M. C., Thiel, M. and Kurths, J. (2007). Recurrence plots for the analysis of complex systems. *Physics Reports*, 438(5-6), 237-329. \doi{10.1016/j.physrep.2006.11.001}
#' @export
#' @examples
#' # x drives y: the recurrences of y (the driven system) are covered by those
#' # of x, so P(R_X | R_Y) is high and its coverage index exceeds theta, while
#' # the reverse coverage stays near the recurrence rate.
#' d <- sim_logistic_bdc(500, "x_to_y", seed = 1)
#' h <- hirata_rp(d$x, d$y, E = 3, n_surr = 20, seed = 1)
#' h
#' \donttest{
#' plot(h)      # the two recurrence plots and the coverage picture
#' # the default 100 surrogates, two workers
#' hirata_rp(d$x, d$y, E = 3, seed = 1, cores = 2)$p_values
#' }
hirata_rp <- function(x, y, E = 3L, tau = 1L, rate = 0.05, theiler = NULL, n_surr = 100L, theta = 0.4,
                      alpha = 0.05, seed = NULL, cores = getOption("Rcausality.cores", 1L)) {
  .local_seed(seed)
  x <- as.numeric(x); y <- as.numeric(y); n <- length(x)
  if (length(y) != n) stop("x and y must have the same length.", call. = FALSE)
  E <- as.integer(E); tau <- as.integer(tau)
  X <- .num_mat(scale(embed_delay(x, E, tau))); Y <- .num_mat(scale(embed_delay(y, E, tau))); ne <- nrow(X)
  if (is.null(theiler)) theiler <- E * tau
  theiler <- max(1L, as.integer(theiler))
  # thresholds = quantiles of the pairwise distances; for long series they are
  # estimated from a random subset of the points (the full set is O(n^2) memory)
  sub_pts <- function(M) if (nrow(M) > 2000L) M[sort(sample.int(nrow(M), 2000L)), , drop = FALSE] else M
  dX <- pair_dist_cpp(sub_pts(X)); dY <- pair_dist_cpp(sub_pts(Y))
  ex <- stats::quantile(dX, c(rate, rate / 5), names = FALSE); ey <- stats::quantile(dY, c(rate, rate / 5), names = FALSE)
  stats_of <- function(s) {
    v <- recurrence_counts_cpp(X, Y, ex[1], ex[2], ey[1], ey[2], as.integer(theiler), as.integer(s))
    np <- v[1]; rx <- v[2] / np; ry <- v[3] / np; rj <- v[6] / np
    c(cov_xy = if (v[5] > 0) v[7] / v[5] else 0,       # P(R_X^loose | R_Y^tight): X covers Y's recurrences -> x drives y
      cov_yx = if (v[4] > 0) v[8] / v[4] else 0,       # P(R_Y^loose | R_X^tight)
      excess = if (rx * ry > 0) rj / (rx * ry) else 0, rate_x = rx, rate_y = ry)
  }
  obs <- stats_of(0L)
  lo <- max(2L * E * tau + 1L, floor(ne / 10)); hi <- ne - lo
  if (hi <= lo) stop("Series too short.", call. = FALSE)
  S <- do.call(rbind, .plapply(seq_len(n_surr), function(b) stats_of(sample(lo:hi, 1L)), cores = cores))
  pv <- c(x_to_y = (1 + sum(S[, "cov_xy"] >= obs[["cov_xy"]])) / (n_surr + 1),
          y_to_x = (1 + sum(S[, "cov_yx"] >= obs[["cov_yx"]])) / (n_surr + 1),
          excess = (1 + sum(S[, "excess"] >= obs[["excess"]])) / (n_surr + 1))
  index <- c(x_to_y = (obs[["cov_xy"]] - obs[["rate_x"]]) / (1 - obs[["rate_x"]]),
             y_to_x = (obs[["cov_yx"]] - obs[["rate_y"]]) / (1 - obs[["rate_y"]]))
  dx <- pv[["x_to_y"]] < alpha && index[["x_to_y"]] >= theta
  dy <- pv[["y_to_x"]] < alpha && index[["y_to_x"]] >= theta
  verdict <- if (dx && dy) "bidirectional" else if (dx) "x_to_y" else if (dy) "y_to_x" else
    if (pv[["excess"]] < alpha) "common_cause" else "independent"
  structure(list(coverage = c(x_to_y = obs[["cov_xy"]], y_to_x = obs[["cov_yx"]]), index = index,
                 excess = obs[["excess"]], rates = c(x = obs[["rate_x"]], y = obs[["rate_y"]]), p_values = pv,
                 surr = S, verdict = verdict, n = ne,
                 settings = list(E = E, tau = tau, rate = rate, theiler = theiler, n_surr = n_surr, theta = theta, alpha = alpha),
                 embedding = list(X = X, Y = Y, eps_x = ex, eps_y = ey)),
            class = "hirata")
}

#' @rdname hirata_rp
#' @param max_points at most this many embedded points are drawn (the first
#'   ones), to keep the recurrence plots legible.
#' @param main title.
#' @param ... unused.
#' @export
plot.hirata <- function(x, max_points = 500L, main = NULL, ...) {
  em <- x$embedding; n <- min(nrow(em$X), max_points); i <- seq_len(n)
  DX <- as.matrix(stats::dist(em$X[i, , drop = FALSE])); DY <- as.matrix(stats::dist(em$Y[i, , drop = FALSE]))
  RX <- DX <= em$eps_x[1]; RY <- DY <= em$eps_y[1]; RYt <- DY <= em$eps_y[2]
  band <- abs(outer(i, i, "-")) < x$settings$theiler          # Theiler window: excluded from the statistics
  RX[band] <- FALSE; RY[band] <- FALSE; RYt[band] <- FALSE
  op <- graphics::par(mfrow = c(1, 3), mar = c(4, 4, 3.5, 1), pty = "s", oma = c(0, 0, 2.5, 0)); on.exit(graphics::par(op))
  draw <- function(M, ttl, col) {
    graphics::image(i, i, M, col = c("white", col), axes = FALSE, xlab = "time i", ylab = "time j", main = ttl, font.main = 2, cex.main = 0.95)
    graphics::axis(1); graphics::axis(2); graphics::box(col = "grey50")
  }
  draw(RX, sprintf("R_X (rate %.2f)", x$settings$rate), "#2B5C8A")
  draw(RY, sprintf("R_Y (rate %.2f)", x$settings$rate), "#B2403A")
  # coverage picture: Y's tight recurrences (dark) that X's loose recurrences cover (light) or miss (red)
  C <- matrix(0L, n, n); C[RX] <- 1L; C[RYt & RX] <- 2L; C[RYt & !RX] <- 3L
  graphics::image(i, i, C, col = c("white", "#D6E4F0", "#1B7A43", "#C0392B"), axes = FALSE, xlab = "time i", ylab = "time j",
                  main = sprintf("tight R_Y on R_X: covered %.0f%%", 100 * x$coverage[["x_to_y"]]), font.main = 2, cex.main = 0.95, zlim = c(0, 3))
  graphics::axis(1); graphics::axis(2); graphics::box(col = "grey50")
  graphics::legend("topleft", legend = c("R_X", "R_Y covered by R_X", "R_Y not covered"), fill = c("#D6E4F0", "#1B7A43", "#C0392B"),
                   bty = "n", cex = 0.7, bg = "white")
  graphics::mtext(if (is.null(main)) sprintf("Recurrence plots (Hirata & Aihara, 2010): verdict %s; coverage x->y %.2f, y->x %.2f, excess %.2f",
                                             x$verdict, x$coverage[["x_to_y"]], x$coverage[["y_to_x"]], x$excess) else main,
                  outer = TRUE, font = 2, cex = 0.9)
  invisible(x)
}

#' @export
print.hirata <- function(x, digits = 3, ...) {
  s <- x$settings
  cat("Recurrence-plot causality (Hirata & Aihara, 2010)\n")
  cat(sprintf("  embedding E = %d, tau = %d, recurrence rate %.2f, Theiler window %d, n = %d, %d surrogates\n",
              s$E, s$tau, s$rate, s$theiler, x$n, s$n_surr))
  tb <- data.frame(criterion = c("x -> y: P(R_X | R_Y) (X covers Y's recurrences)", "y -> x: P(R_Y | R_X)",
                                 "simultaneous recurrences: joint rate / (rate_x rate_y)"),
                   value = round(c(x$coverage, x$excess), digits),
                   index = c(round(x$index, digits), NA),
                   p.surr = .fmt_p(x$p_values))
  print(tb, row.names = FALSE, right = FALSE)
  cat(sprintf("  verdict: %s   (index >= %g and p < %g => directed link; excess significant without a link => common cause)\n",
              x$verdict, s$theta, s$alpha))
  invisible(x)
}

# ---- correlation dimension (Krakovska 2019) ----------------------------------
#' Correlation dimension (Grassberger-Procaccia) with a bootstrap standard error
#'
#' Slope of `log C(r)` against `log r` over the scaling region between the
#' `q[1]` and `q[2]` quantiles of the pairwise distances of the embedded
#' points; the standard error comes from `n_boot` resamples of the points.
#'
#' @param M numeric matrix (rows = points), e.g. a delay embedding.
#' @param q quantiles delimiting the scaling region.
#' @param n_r number of radii in the region.
#' @param n_boot bootstrap resamples (0 = none).
#' @param max_points at most this many points are used (random subset).
#' @param cores CPU cores for the bootstrap (see [Rcausality-parallel]).
#' @return A list with `d2`, `se`, `radii`, `corr_sum`.
#' @references
#' Grassberger, P. and Procaccia, I. (1983). Measuring the strangeness of strange attractors. *Physica D*, 9(1-2), 189-208. \doi{10.1016/0167-2789(83)90298-1}
#' @export
#' @examples
#' # a one-dimensional map embedded in three dimensions: D2 is about 1
#' d <- sim_logistic(1000, "independent", seed = 1)
#' cd <- corr_dim(embed_delay(d$x, 3))
#' cd$d2        # the estimate
#' cd$se        # its bootstrap standard error
corr_dim <- function(M, q = c(0.02, 0.2), n_r = 12L, n_boot = 20L, max_points = 1500L,
                     cores = getOption("Rcausality.cores", 1L)) {
  M <- .num_mat(M); n <- nrow(M)
  if (n > max_points) M <- M[sort(sample.int(n, max_points)), , drop = FALSE]
  d <- pair_dist_cpp(M); d <- d[d > 0]
  radii <- exp(seq(log(stats::quantile(d, q[1], names = FALSE)), log(stats::quantile(d, q[2], names = FALSE)), length.out = n_r))
  cs <- vapply(radii, function(r) mean(d < r), 0)
  fit <- function(cs_) stats::coef(stats::lm(log(cs_) ~ log(radii)))[2]
  d2 <- fit(cs)
  se <- NA_real_
  if (n_boot > 0) {
    nn <- nrow(M)
    b <- unlist(.plapply(seq_len(n_boot), function(i) {
      Mi <- M[sample.int(nn, nn, replace = TRUE), , drop = FALSE]
      di <- pair_dist_cpp(Mi); di <- di[di > 1e-12]
      fit(vapply(radii, function(r) max(mean(di < r), 1 / length(di)), 0))
    }, cores = cores))
    se <- stats::sd(b)
  }
  list(d2 = unname(d2), se = se, radii = radii, corr_sum = cs)
}

#' Correlation-dimension causality (Krakovska, 2019)
#'
#' Compares the correlation dimensions of the reconstructed state spaces of
#' `x`, of `y`, of their joint and of an independent joint: if `x` drives
#' `y` the joint has the dimension of `y` and `D2(x) < D2(y)`; mutual
#' coupling gives equal dimensions; a hidden common driver makes the joint
#' dimension exceed both but stay below the independent joint; independence
#' makes the joint as large as the independent joint. The joint is the
#' embedding of the additive observation `a x + y` (`a = sqrt(29/31)`) and
#' the independent joint `a X + Y*` with time-shuffled `Y*`, as in Benko et
#' al. (2024), so that all four dimensions carry the same embedding bias;
#' relations are tested with the bootstrap standard errors of the estimates
#' (z-tests at level `alpha`).
#'
#' @param x,y numeric vectors.
#' @param E,tau embedding dimension and delay.
#' @param alpha significance level of the equality tests.
#' @param tol relative tolerance: two dimensions are also taken as equal when
#'   they differ by less than `tol` times their mean (the bootstrap standard
#'   errors do not cover the embedding-dependent bias of the estimates).
#' @param n_boot bootstrap resamples per dimension estimate.
#' @param ... passed to [corr_dim()].
#' @param seed optional RNG seed (shuffle of the independent joint).
#' @param cores CPU cores for the bootstrap resamples (see [Rcausality-parallel]).
#' @return An object of class `krakovska`: `d2` (named vector: `x`, `y`,
#'   `joint`, `indep`), `se`, `verdict`, `tests` (outcomes of the relation
#'   tests). Method: `print`.
#' @references
#' Krakovska, A. (2019). Correlation dimension detects causal links in coupled dynamical systems. *Entropy*, 21(9), 818. \doi{10.3390/e21090818}
#'
#' Krakovska, A., Jakubik, J., Chvostekova, M., Coufal, D., Jajcay, N. and Palus, M. (2018). Comparison of six methods for the detection of causality in a bivariate time series. *Physical Review E*, 97, 042207. \doi{10.1103/PhysRevE.97.042207}
#' @export
#' @examples
#' # x drives y: the joint attractor has the dimension of y (the driven system
#' # already contains its driver) and D2(x) < D2(y), so the verdict is x_to_y
#' d <- sim_logistic_bdc(600, "x_to_y", seed = 1)
#' krakovska_cd(d$x, d$y, E = 3, n_boot = 5, seed = 1)
#' \donttest{
#' krakovska_cd(d$x, d$y, E = 3, seed = 1)   # the default 20 bootstrap resamples
#' }
krakovska_cd <- function(x, y, E = 3L, tau = 1L, alpha = 0.05, tol = 0.1, n_boot = 20L, seed = NULL,
                         cores = getOption("Rcausality.cores", 1L), ...) {
  .local_seed(seed)
  x <- .rank_norm(as.numeric(x)); y <- .rank_norm(as.numeric(y)); a <- sqrt(29 / 31)
  X <- embed_delay(x, E, tau); Y <- embed_delay(y, E, tau); J <- embed_delay(a * x + y, E, tau)
  Z <- a * X + Y[sample.int(nrow(Y)), , drop = FALSE]
  cx <- corr_dim(X, n_boot = n_boot, cores = cores, ...); cy <- corr_dim(Y, n_boot = n_boot, cores = cores, ...)
  cj <- corr_dim(J, n_boot = n_boot, cores = cores, ...); cz <- corr_dim(Z, n_boot = n_boot, cores = cores, ...)
  d2 <- c(x = cx$d2, y = cy$d2, joint = cj$d2, indep = cz$d2); se <- c(x = cx$se, y = cy$se, joint = cj$se, indep = cz$se)
  zc <- stats::qnorm(1 - alpha / 2)
  eq <- function(a, b, sa, sb) abs(a - b) <= max(zc * sqrt(sa^2 + sb^2), tol * (a + b) / 2)
  eq_jx <- eq(d2[["joint"]], d2[["x"]], se[["joint"]], se[["x"]])
  eq_jy <- eq(d2[["joint"]], d2[["y"]], se[["joint"]], se[["y"]])
  eq_xy <- eq(d2[["x"]], d2[["y"]], se[["x"]], se[["y"]])
  eq_jz <- eq(d2[["joint"]], d2[["indep"]], se[["joint"]], se[["indep"]])
  verdict <- if (eq_jx && eq_jy && eq_xy) "bidirectional" else
    if (eq_jy && d2[["x"]] < d2[["y"]]) "x_to_y" else if (eq_jx && d2[["y"]] < d2[["x"]]) "y_to_x" else
    if (eq_jz) "independent" else if (d2[["joint"]] > max(d2[["x"]], d2[["y"]]) && d2[["joint"]] < d2[["indep"]]) "common_cause" else
    "undetermined"
  structure(list(d2 = d2, se = se, verdict = verdict,
                 tests = c(joint_eq_x = eq_jx, joint_eq_y = eq_jy, x_eq_y = eq_xy, joint_eq_indep = eq_jz),
                 settings = list(E = E, tau = tau, alpha = alpha, tol = tol)), class = "krakovska")
}

#' @export
print.krakovska <- function(x, digits = 2, ...) {
  cat("Correlation-dimension causality (Krakovska, 2019)\n")
  cat(sprintf("  embedding E = %d, tau = %d\n", x$settings$E, x$settings$tau))
  print(data.frame(manifold = c("X", "Y", "joint (a x + y)", "independent joint"), D2 = round(x$d2, digits), se = round(x$se, digits)),
        row.names = FALSE, right = FALSE)
  cat(sprintf("  relations (z-tests): joint = X: %s, joint = Y: %s, X = Y: %s, joint = independent joint: %s\n",
              x$tests[["joint_eq_x"]], x$tests[["joint_eq_y"]], x$tests[["x_eq_y"]], x$tests[["joint_eq_indep"]]))
  cat(sprintf("  verdict: %s\n", x$verdict))
  invisible(x)
}

# ---- rcausality drivers --------------------------------------------------------
.rcausality_verdict_pairs <- function(all, grid, causes, effects, alpha, method, fun, cl, cores, stat_fun, null_fun = NULL) {
  up <- .unordered_pairs(grid, all)
  sc <- .split_cores(nrow(up), cores)
  res <- .plapply(seq_len(nrow(up)), function(i) fun(all[[up$v1[i]]], all[[up$v2[i]]], cores = sc$inner), cores = sc$outer)
  rows <- list(); ins <- list(); extra <- list(); nulls <- list()
  for (i in seq_len(nrow(up))) {
    r <- res[[i]]; a <- up$v1[i]; b <- up$v2[i]; st <- stat_fun(r)
    for (dir in c("ab", "ba")) {
      k <- if (dir == "ab") "x_to_y" else "y_to_x"
      if (!is.null(null_fun)) nulls[[paste(if (dir == "ab") a else b, "->", if (dir == "ab") b else a)]] <- null_fun(r, k)
      rows[[length(rows) + 1L]] <- data.frame(
        cause = if (dir == "ab") a else b, effect = if (dir == "ab") b else a, lag = NA_integer_,
        statistic = st$statistic[[k]], df1 = NA_real_, df2 = NA_real_, p_value = st$p[[k]], bf10 = NA_real_,
        effect_size = st$statistic[[k]], r2_gain = NA_real_, sign = NA_real_, n_eff = st$n,
        significant = r$verdict %in% c(k, "bidirectional"), stringsAsFactors = FALSE)
    }
    ins[[i]] <- data.frame(var1 = a, var2 = b, cmi = NA_real_, p_value = NA_real_, p_adj = NA_real_,
                           significant = r$verdict == "common_cause", stringsAsFactors = FALSE)
    extra[[i]] <- cbind(data.frame(var1 = a, var2 = b, stringsAsFactors = FALSE), st$extra, verdict = r$verdict)
  }
  pairs <- do.call(rbind, rows)
  pairs$p_adj <- pairs$p_value; pairs$post_prob <- NA_real_; pairs$evidence <- NA_character_
  pairs <- pairs[, c("cause", "effect", "lag", "statistic", "df1", "df2", "p_value",
                     "p_adj", "bf10", "post_prob", "evidence", "effect_size",
                     "r2_gain", "sign", "n_eff", "significant")]
  rownames(pairs) <- NULL
  out <- .new_causres(pairs, series = names(all), causes = causes, effects = effects,
                      method = method, conditional = FALSE, alpha = alpha, bf_threshold = NA_real_,
                      adjust = "none", n = nrow(all), call = cl, data = all, instant = do.call(rbind, ins),
                      ccdh = do.call(rbind, extra))
  if (length(nulls)) out$null <- nulls
  out
}
