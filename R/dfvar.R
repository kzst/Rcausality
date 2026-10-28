# Degrees-of-freedom (df) causality: conditional-variance analysis
# (Telcs, Kurbucz & Jakovac, "Unified causality analysis based on the degrees
#  of freedom", Physical Review E, 2025).
#
# The number of degrees of freedom of a series is the number of past values
# that must be constrained for the distribution (in practice: the variance) of
# the present value to stop shrinking.  The orders of x, of y and of the joint
# process (x, y) identify the causal case (Table 2 of the paper).  Constraints
# are imposed as in Eqs. (14)-(15): every observed past configuration serves
# as a constraint centre and the samples whose past lies within a window of
# width sigma_W (max-norm) around it form the conditioned sample; the pooled
# variance of the present values over all centres is sigma(sigma_W, n_constr).
# The decision level uses adaptive (k-nearest-neighbour) windows.

# Conditional-variance kernels (C++, src/knn.cpp): `C` is the constraint
# matrix (rows = centres, columns = constrained past values; 0 columns = no
# constraint), `pres` the matrix of present values.  .cond_var_knn returns the
# mean over the columns of pres of the log pooled within-window variance with
# adaptive windows (k nearest neighbours of every centre in the max-norm);
# .cond_var_radius the same with fixed-width windows (radius r).
.cond_var_knn <- function(pres, C, k) cond_var_knn_cpp(C, pres, as.integer(k))
.cond_var_radius <- function(pres, C, r, min_n = 5L) {
  v <- cond_var_radius_cpp(C, pres, r, as.integer(min_n)); c(sigma2 = v[1], n = v[2])
}
.empty_C <- function(n) matrix(numeric(0), n, 0)
.add_cols <- function(C, cols) do.call(cbind, c(list(C), cols))

# Does adding the constraint `newcols` (list of aligned lag vectors, e.g. the
# components of one block at one lag) reduce the conditional variance of
# `pres` given the constraints already in C?  Statistic: log-variance
# reduction; null: local permutation of the new columns among the k_perm
# nearest neighbours in the current constraint space (Runge, 2018; a plain
# permutation when there is no constraint yet).
.cv_test <- function(pres, C, newcols, k, n_perm, k_perm = 5L, cores = 1L) {
  ne <- nrow(pres)
  base <- .cond_var_knn(pres, C, k)
  obs <- base - .cond_var_knn(pres, .add_cols(C, newcols), k)
  nn <- if (ncol(C) > 0L) knn_idx_cpp(C, as.integer(k_perm)) else NULL
  null <- unlist(.plapply(seq_len(n_perm), function(b) {
    perm <- if (is.null(nn)) sample.int(ne) else .local_perm(nn, k_perm)
    base - .cond_var_knn(pres, .add_cols(C, lapply(newcols, function(v) v[perm])), k)
  }, cores = cores))
  list(statistic = obs, p_value = (1 + sum(null >= obs)) / (n_perm + 1),
       sigma_ratio = exp(-obs / 2), C_new = .add_cols(C, newcols), null = null)
}

# lag columns of a block (matrix) aligned to present index t_idx
.block_lags <- function(M, t_idx, lags) {
  out <- list()
  for (l in lags) for (j in seq_len(ncol(M))) out[[length(out) + 1L]] <- M[t_idx - l, j]
  out
}

#' Degrees of freedom of a series by conditional-variance analysis
#'
#' Estimates how many past values of a series (or block of series) must be
#' constrained before the variance of its present value stops decreasing --
#' the number of degrees of freedom (df) of Telcs, Kurbucz and Jakovac (2025).
#' Constraints are windows of width `sigma_W` (max-norm) around every observed
#' past configuration; `sigma(sigma_W, n_constr)` is the pooled standard
#' deviation of the present values inside the windows (Eqs. 14-15 of the
#' paper). The df is read off the adaptive-window (k nearest neighbours)
#' profile: constraints are added one lag at a time while the conditional
#' standard deviation still drops by more than `tol` (relative).
#'
#' @param x numeric vector or matrix (columns = components of one subsystem;
#'   a constraint then fixes the whole lag vector).
#' @param max_constr largest number of past lags constrained.
#' @param k neighbourhood size of the adaptive windows (default 10; small
#'   values keep the windows narrow, which matters for deterministic systems).
#' @param widths window widths (in units of the standard deviation of `x`)
#'   for the `sigma(sigma_W)` profile; default a geometric grid.
#' @param alpha significance level of the step tests.
#' @param tol minimal relative decrease of the conditional standard deviation
#'   that still counts as a new degree of freedom (guards against significant
#'   but negligible reductions in long series).
#' @param n_perm,k_perm number of local permutations and neighbourhood size of
#'   the step tests (see Details).
#' @param min_n minimum number of samples in a fixed-width window.
#' @param seed optional RNG seed.
#' @param cores CPU cores for the permutations (see [Rcausality-parallel]).
#' @details Each additional lag is accepted as a degree of freedom when the
#'   reduction of the conditional variance it brings is significant against a
#'   local-permutation null (the lag values are shuffled among the `k_perm`
#'   nearest neighbours in the space of the constraints already imposed, cf.
#'   Runge, 2018) and larger than `tol`. In deterministic systems the
#'   reductions are by orders of magnitude and the procedure reproduces the
#'   "count the saturated lines" reading of the paper's figures; in stochastic
#'   systems the test makes the reading objective.
#' @return An object of class `dfvar`: `df` (estimated degrees of freedom, in
#'   lags, from the step tests), `df_profile` (the reading of the
#'   sigma-versus-width profile as in the figures of the paper: the first
#'   number of constraints whose line is proportional to the window width,
#'   otherwise the number of distinct saturated levels), `steps` (per number of constraints: conditional sd, relative
#'   reduction, p-value), `profile` (data frame: `n_constr`, `width`, `sigma`,
#'   `n`), `n`, `k`, `call`. Methods: `print`, `plot` (the sigma-versus-width
#'   figure of the paper).
#' @references
#' Telcs, A., Kurbucz, M. T. and Jakovac, A. (2025). Unified causality analysis based on the degrees of freedom. *Physical Review E*, 112(1), 014204. \doi{10.1103/PhysRevE.112.014204}, Sections 3.1, 6.1 and 6.4.
#'
#' Runge, J. (2018). Conditional independence testing based on a nearest-neighbor estimator of conditional mutual information. *Proceedings of the 21st International Conference on Artificial Intelligence and Statistics (AISTATS)*, PMLR 84, 938-947.
#' @export
#' @examples
#' # An AR(2) process has two degrees of freedom: constraining x[t-1] and then
#' # x[t-2] reduces the conditional variance of x[t] significantly, a third
#' # lag does not.
#' set.seed(1)
#' x <- as.numeric(stats::arima.sim(list(ar = c(0.5, 0.3)), 300))
#' d <- df_variance(x, max_constr = 3)
#' d              # df = 2 (two accepted lags)
#' \donttest{
#' plot(d)        # the variance-reduction path and the permutation p-values
#' }
df_variance <- function(x, max_constr = 6L, k = NULL, widths = NULL, alpha = 0.05, tol = 0.02,
                        n_perm = 49L, k_perm = 5L, min_n = 5L, seed = NULL, cores = getOption("Rcausality.cores", 1L)) {
  cl <- match.call()
  old <- .rng_setup(seed); on.exit(.rng_restore(old), add = TRUE)
  if (n_perm < ceiling(1 / alpha)) {
    message("df_variance: n_perm raised to ", ceiling(1 / alpha), " so that p < alpha is attainable.")
    n_perm <- ceiling(1 / alpha)
  }
  M <- as.matrix(x); storage.mode(M) <- "double"
  n <- nrow(M); if (n < 50) stop("Series too short.", call. = FALSE)
  max_constr <- as.integer(max_constr)
  t_idx <- (max_constr + 1):n; ne <- length(t_idx)
  if (is.null(k)) k <- 10L
  sds <- apply(M, 2, stats::sd); sds[sds == 0] <- 1
  Ms <- sweep(M, 2, sds, "/")                       # standardised for the max-norm
  pres <- M[t_idx, , drop = FALSE]
  if (is.null(widths)) widths <- 2^seq(1, -7, by = -0.5)
  C <- .empty_C(ne)
  steps <- data.frame(n_constr = 0:max_constr, sigma = NA_real_, reduction = NA_real_, p_value = NA_real_)
  steps$sigma[1] <- sqrt(exp(.cond_var_knn(pres, C, k)))
  prof <- list(); df <- 0L; open <- TRUE
  for (j in 0:max_constr) {
    if (j > 0) {
      tst <- .cv_test(pres, C, .block_lags(Ms, t_idx, j), k, n_perm, k_perm, cores)
      C <- tst$C_new
      steps$sigma[j + 1L] <- sqrt(exp(.cond_var_knn(pres, C, k)))
      steps$reduction[j + 1L] <- 1 - tst$sigma_ratio
      steps$p_value[j + 1L] <- tst$p_value
      if (open && tst$p_value < alpha && (1 - tst$sigma_ratio) > tol) df <- j else open <- FALSE
    }
    for (w in widths) {
      r <- .cond_var_radius(pres, C, w, min_n)
      prof[[length(prof) + 1L]] <- data.frame(n_constr = j, width = w, sigma = sqrt(r[["sigma2"]]),
                                              n = r[["n"]])
    }
  }
  profile <- do.call(rbind, prof)
  structure(list(df = df, df_profile = .df_from_profile(profile, tol = 0.1), steps = steps, sigma_knn = stats::setNames(steps$sigma, paste0("n_constr=", 0:max_constr)),
                 profile = profile, n = n, k = k, alpha = alpha, tol = tol,
                 series = if (is.null(colnames(M))) "x" else paste(colnames(M), collapse = ","), call = cl),
            class = "dfvar")
}

# Reading of the sigma(sigma_W) profile as in the figures of the paper:
# the df is the first number of constraints whose line is proportional to the
# window width at small widths (the present value is then fully determined);
# if no line is proportional (stochastic system) it is the number of distinct
# saturated levels beyond the unconstrained one.
.df_from_profile <- function(P, tol = 0.1, min_pts = 3L) {
  nc <- sort(unique(P$n_constr)); level <- rep(NA_real_, length(nc)); prop <- logical(length(nc))
  for (i in seq_along(nc)) {
    p <- P[P$n_constr == nc[i] & is.finite(P$sigma) & P$sigma > 0, ]; p <- p[order(p$width), ]
    if (nrow(p) < min_pts) next
    level[i] <- p$sigma[1]
    q <- p[seq_len(min(4L, nrow(p))), ]
    prop[i] <- nc[i] > 0 && stats::coef(stats::lm(log(q$sigma) ~ log(q$width)))[2] > 0.5
  }
  if (any(prop)) return(nc[which(prop)[1]])
  df <- 0L
  for (i in seq_along(nc)[-1]) {
    if (is.na(level[i]) || is.na(level[i - 1])) break
    if (level[i] < (1 - tol) * level[i - 1]) df <- nc[i] else break
  }
  df
}

#' @export
print.dfvar <- function(x, digits = 3, ...) {
  cat("Degrees of freedom by conditional-variance analysis (Telcs, Kurbucz & Jakovac, 2025)\n")
  cat(sprintf("  series: %s, n = %d, k = %d, alpha = %g, tol = %g\n", x$series, x$n, x$k, x$alpha, x$tol))
  cat("  constraints on the past vs conditional sd of the present value:\n")
  S <- x$steps
  print(data.frame(n_constr = S$n_constr, cond.sd = round(S$sigma, digits),
                   reduction = ifelse(is.na(S$reduction), "", sprintf("%.1f%%", 100 * S$reduction)),
                   p = ifelse(is.na(S$p_value), "", .fmt_p(S$p_value)),
                   accepted = ifelse(is.na(S$p_value), "", ifelse(S$n_constr <= x$df, "yes", "no"))),
        row.names = FALSE)
  cat(sprintf("  estimated df (past lags that significantly reduce the sd, adaptive windows): %d\n", x$df))
  cat(sprintf("  profile reading (first line proportional to the width, else distinct saturated levels): %d\n",
              x$df_profile))
  invisible(x)
}

#' @rdname df_variance
#' @param main title.
#' @param ... further arguments (unused).
#' @export
plot.dfvar <- function(x, main = NULL, ...) {
  P <- x$profile; P <- P[is.finite(P$sigma) & P$sigma > 0, ]
  nc <- sort(unique(P$n_constr)); cols <- grDevices::hcl.colors(length(nc), "Zissou 1")
  op <- graphics::par(mar = c(4.5, 4.5, 3, 1)); on.exit(graphics::par(op))
  graphics::plot(NA, xlim = range(P$width), ylim = range(P$sigma), log = "xy",
                 xlab = expression(paste("window width ", sigma[W], " (in sd units)")),
                 ylab = expression(paste("conditional sd ", sigma, " of the present value")),
                 main = if (is.null(main)) sprintf("df profile: %s (df = %d, profile reading %d)", x$series, x$df, x$df_profile) else main,
                 font.main = 2)
  for (i in seq_along(nc)) {
    p <- P[P$n_constr == nc[i], ]; p <- p[order(p$width), ]
    graphics::lines(p$width, p$sigma, col = cols[i], lwd = 2)
    graphics::points(p$width, p$sigma, col = cols[i], pch = 16, cex = 0.7)
  }
  graphics::legend("bottomright", legend = paste0("n_constr = ", nc), col = cols, lwd = 2, bty = "n",
                   cex = 0.7, ncol = 2)
  graphics::mtext("flat line = saturated constraint set; line proportional to the width = present value determined",
                  side = 3, line = 0.2, cex = 0.7, col = "grey40")
  invisible(x)
}

# ---- joint df, own df and the directed variance-reduction test ----------------

# Order (df) of the present `pres` (matrix of standardised present values)
# with constraints drawn greedily from the pasts of the observed series in
# `cands` (named list of standardised vectors): at every step the candidate
# (next unused lag of each series) with the larger observed reduction is
# tested first; stop when no candidate brings a significant reduction.
# This is the "number of conditions" of Table 1 of Telcs et al. (2025): for
# y alone it may include conditions on x's past (y's order is (ox, oy) when
# x drives y).
.df_order <- function(pres, cands, t_idx, max_constr, k, alpha, tol, n_perm, k_perm, cores) {
  ne <- length(t_idx); C <- .empty_C(ne)
  used <- stats::setNames(integer(length(cands)), names(cands)); chosen <- character(0); pv <- numeric(0)
  repeat {
    cand <- list()
    for (nm in names(cands)) if (used[[nm]] < max_constr)
      cand[[nm]] <- list(cands[[nm]][t_idx - used[[nm]] - 1L])
    if (!length(cand)) break
    base <- .cond_var_knn(pres, C, k)
    obs <- vapply(cand, function(cl) base - .cond_var_knn(pres, .add_cols(C, cl), k), 0)
    accepted <- FALSE
    for (b in names(sort(obs, decreasing = TRUE))) {
      tst <- .cv_test(pres, C, cand[[b]], k, n_perm, k_perm, cores)
      if (tst$p_value < alpha && (1 - tst$sigma_ratio) > tol) {
        C <- tst$C_new; accepted <- TRUE; used[[b]] <- used[[b]] + 1L
        chosen <- c(chosen, sprintf("%s[t-%d]", b, used[[b]])); pv <- c(pv, tst$p_value)
        break
      }
    }
    if (!accepted) break
  }
  list(df = sum(used), lags = used, chosen = chosen, p = pv)
}

# Directed test cause -> effect at lag tau: does constraining cause[t - tau]
# reduce the conditional variance of effect[t] beyond the effect's own
# `history` past (and the controls)?  Statistic = relative reduction of the
# conditional sd; null = local permutation (Runge, 2018).
.df_pair <- function(effect, cause, controls, lags, history, k, n_perm, k_perm = 5L, cores = 1L) {
  n <- length(effect)
  Lmax <- max(lags, history); t_idx <- (Lmax + 1):n; ne <- length(t_idx)
  ys <- effect / stats::sd(effect); xs <- cause / stats::sd(cause)
  pres <- cbind(ys[t_idx])
  Cz <- .empty_C(ne)
  for (l in seq_len(history)) Cz <- cbind(Cz, ys[t_idx - l])
  if (!is.null(controls)) for (j in seq_len(ncol(controls))) {
    cs <- controls[, j] / stats::sd(controls[, j])
    for (l in seq_len(history)) Cz <- cbind(Cz, cs[t_idx - l])
  }
  tests <- lapply(lags, function(tau) .cv_test(pres, Cz, list(xs[t_idx - tau]), k, n_perm, k_perm, cores))
  obs <- vapply(tests, `[[`, 0, "statistic"); pl <- vapply(tests, `[[`, 0, "p_value")
  null_max <- do.call(pmax, lapply(tests, `[[`, "null"))        # max-statistic test over the lags
  b <- which.max(obs)
  list(lag = lags[b], statistic = 1 - exp(-obs[b] / 2), p_value = (1 + sum(null_max >= obs[b])) / (length(null_max) + 1),
       n_eff = ne, by_lag = data.frame(lag = lags, statistic = 1 - exp(-obs / 2), p = pl), null = 1 - exp(-null_max / 2))
}

# Table 2 of Telcs, Kurbucz & Jakovac (2025): orders of x, y and the joint.
# A series of order 0 (white noise given every candidate past) has no
# dynamics and cannot be driven; whether it drives the other series is read
# from the conditions chosen for the other series (`x_drives_y`: a past
# value of x was accepted when determining the order of y).
.df_decide <- function(dx, dy, dj, x_drives_y = NA, y_drives_x = NA) {
  if (dx == 0 && dy == 0) return(if (dj == 0) "independent" else "undetermined")
  if (dx == 0) return(if (isTRUE(x_drives_y)) "x_to_y" else if (isFALSE(x_drives_y)) "independent" else "undetermined")
  if (dy == 0) return(if (isTRUE(y_drives_x)) "y_to_x" else if (isFALSE(y_drives_x)) "independent" else "undetermined")
  if (dj == dx + dy) "independent"
  else if (dx < dy && dj == dy) "x_to_y"
  else if (dy < dx && dj == dx) "y_to_x"
  else if (dx == dy && dj == dx) "bidirectional"
  else if (dj > max(dx, dy) && dj < dx + dy) "common_cause"
  else "undetermined"
}
