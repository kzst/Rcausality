# Bayesian dimensional causality (Benko, Zlatniczki, Stippinger, Fabo, Solyom,
# Eross, Telcs & Somogyvari, "Bayesian inference of causal relations between
# dynamical systems", Chaos, Solitons & Fractals 185 (2024) 115142).
#
# Independent implementation written from the paper (Sections 2.1-2.3,
# Appendices A.2 and B.1).  An AGPL-3 reference implementation by the authors
# exists (github.com/adam-zlatniczki/dimensional_causality); no code was
# taken from it.

# Moore-Penrose inverse (fallback for a singular dimension covariance)
.ginv <- function(S) { e <- eigen(S, symmetric = TRUE); v <- e$values; v <- ifelse(v > 1e-10 * max(v), 1 / v, 0); e$vectors %*% diag(v, length(v)) %*% t(e$vectors) }
.rank_norm <- function(v) rank(v, ties.method = "first") / (length(v) + 1)   # ties broken by position: no coincident points

# the nine mutually exclusive events A_ij of Appendix A.2 (rows: D_X < D_Y,
# D_X > D_Y, D_X = D_Y; columns: direct drive, common cause, independent).
# `free` = names of the free coordinates, `map` = which free coordinate gives
# (D_X, D_Y, D_J, D_Z), `order` = the chain of strict inequalities.
.dc_events <- list(
  A11 = list(free = c("x", "y", "z"),      map = c("x", "y", "y", "z"), order = c("x", "y", "z")),
  A12 = list(free = c("x", "y", "j", "z"), map = c("x", "y", "j", "z"), order = c("x", "y", "j", "z")),
  A13 = list(free = c("x", "y", "j"),      map = c("x", "y", "j", "j"), order = c("x", "y", "j")),
  A21 = list(free = c("x", "y", "z"),      map = c("x", "y", "x", "z"), order = c("y", "x", "z")),
  A22 = list(free = c("x", "y", "j", "z"), map = c("x", "y", "j", "z"), order = c("y", "x", "j", "z")),
  A23 = list(free = c("x", "y", "j"),      map = c("x", "y", "j", "j"), order = c("y", "x", "j")),
  A31 = list(free = c("x", "z"),           map = c("x", "x", "x", "z"), order = c("x", "z")),
  A32 = list(free = c("x", "j", "z"),      map = c("x", "x", "j", "z"), order = c("x", "j", "z")),
  A33 = list(free = c("x", "j"),           map = c("x", "x", "j", "j"), order = c("x", "j")))

# likelihood of the measured dimension vector d (D_X, D_Y, D_J, D_Z) with
# covariance S under every event: the 4-variate normal density integrated
# over the event's pyramid (grid integration, +/- c standard deviations,
# `bins` points per free coordinate), Appendix A.2.
.dc_likelihoods <- function(d, S, c = 3, bins = 25L) {
  S <- S + diag(1e-8 * max(1, mean(diag(S))), 4)          # ridge: discrete / tied data can make S singular
  Sinv <- tryCatch(solve(S), error = function(e) .ginv(S))
  sds <- pmax(sqrt(diag(S)), 1e-6); names(d) <- names(sds) <- c("x", "y", "j", "z")
  vapply(.dc_events, function(ev) {
    rng <- lapply(ev$free, function(v) {
      i <- which(ev$map == v)
      c(max(0, min(d[i] - c * sds[i])), max(d[i] + c * sds[i]))
    })
    grids <- lapply(rng, function(r) seq(r[1], r[2], length.out = bins))
    dw <- prod(vapply(rng, function(r) diff(r) / (bins - 1), 0))
    G <- as.matrix(expand.grid(grids)); colnames(G) <- ev$free
    ok <- rep(TRUE, nrow(G))
    for (l in seq_len(length(ev$order) - 1L)) ok <- ok & G[, ev$order[l]] < G[, ev$order[l + 1L]]
    if (!any(ok)) return(0)
    W <- G[ok, ev$map, drop = FALSE]
    R <- sweep(W, 2, d)
    q <- rowSums((R %*% Sinv) * R)
    sum(exp(-0.5 * q)) * dw
  }, 0)
}

#' Bayesian dimensional causality (Benko et al., 2024)
#'
#' Infers the causal relation between two deterministic dynamical systems
#' from the intrinsic dimensions of their time-delay embedded manifolds
#' `X`, `Y`, their joint `J` (embedding of `a x + y`, `a = sqrt(29/31)`) and
#' the independent joint `Z` (`a X + Y*`, `Y*` time-shuffled): `X -> Y` iff
#' `D_X < D_Y = D_J`, `X <-> Y` iff `D_X = D_Y = D_J`, hidden common cause iff
#' `max(D_X, D_Y) < D_J < D_Z`, independence iff `D_J = D_Z` (Eqs. 1-4).
#' Local dimensions are estimated with the k-nearest-neighbour estimator of
#' Farahmand, Szepesvari and Audibert (2007), `log(2) / log(r_k / r_{k/2})`,
#' trimmed (`eps` tails) and averaged; the sampling distribution of the four
#' global estimates is taken as multivariate normal with the effective
#' sample size `n / (2k)` (Appendix A.2), and the posterior probability of
#' the five causal cases is obtained by integrating this density over the
#' dimension pyramids of the nine elementary events with a uniform prior,
#' averaged over the neighbourhood sizes in `k_range`.
#'
#' @param x,y numeric vectors (the two time series).
#' @param m embedding dimension (the paper used 3-4 for logistic maps).
#' @param tau embedding delay.
#' @param k_range neighbourhood sizes over which the posterior is averaged
#'   (choose a range where the dimension estimates form a plateau, see
#'   `plot()`).
#' @param prep preprocessing: `"rank"` (rank normalisation, the paper's
#'   choice), `"scale"` (z-scores) or `"none"`.
#' @param eps trimming fraction of the local dimension estimates (each tail).
#' @param downsample keep every `downsample`-th embedded point (for flows
#'   sampled much faster than their dynamics).
#' @param c,bins half-width (in standard deviations) and resolution of the
#'   integration grid.
#' @param seed optional RNG seed (the shuffle defining `Z`).
#' @return An object of class `dimcaus`: `probs` (named vector: `x_to_y`,
#'   `y_to_x`, `bidirectional`, `common_cause`, `independent`), `best`,
#'   `by_k` (per neighbourhood size: dimension estimates, their standard
#'   errors and the five probabilities), `settings`, `call`. Methods:
#'   `print`, `plot` (dimension-versus-k and posterior panels, Fig. 1B-C).
#' @details Independent implementation written from the paper; the authors'
#'   AGPL-3 reference implementation (`dimensional_causality` on GitHub) was
#'   not used. The method assumes deterministic, stationary systems;
#'   observational noise must be filtered beforehand.
#' @references
#' Benko, Z., Zlatniczki, A., Stippinger, M., Fabo, D., Solyom, A., Eross, L., Telcs, A. and Somogyvari, Z. (2024). Bayesian inference of causal relations between dynamical systems. *Chaos, Solitons & Fractals*, 185, 115142. \doi{10.1016/j.chaos.2024.115142}
#'
#' Farahmand, A. M., Szepesvari, C. and Audibert, J.-Y. (2007). Manifold-adaptive dimension estimation. *Proceedings of the 24th International Conference on Machine Learning (ICML)*, 265-272. \doi{10.1145/1273496.1273530}
#'
#' Takens, F. (1981). Detecting strange attractors in turbulence. In: *Dynamical Systems and Turbulence, Warwick 1980*, Lecture Notes in Mathematics 898, 366-381. Springer. \doi{10.1007/BFb0091924}
#' @export
#' @examples
#' # Logistic maps where x drives y.  The driven series inherits the dimension
#' # of its driver, so D_X < D_Y = D_J and the posterior favours x_to_y.
#' d <- sim_logistic_bdc(800, "x_to_y", seed = 1)
#' r <- dimensional_causality(d$x, d$y, m = 3, k_range = c(10, 15), seed = 1)
#' r
#' r$dims         # the four manifold dimensions (X, Y, joint, independent joint)
#' \donttest{
#' # a hidden common cause: the joint is more complex than either series but
#' # simpler than the independent joint; longer series and a plateau of k
#' d <- sim_logistic_bdc(3000, "common_cause", seed = 1)
#' r <- dimensional_causality(d$x, d$y, m = 4, k_range = seq(12, 36, by = 8), seed = 1)
#' r
#' plot(r)        # dimension estimates against k and the posterior of the five cases
#' }
dimensional_causality <- function(x, y, m = 3L, tau = 1L, k_range = seq(10L, 40L, by = 5L),
                                  prep = c("rank", "scale", "none"), eps = 0.05, downsample = 1L,
                                  c = 3, bins = 25L, seed = NULL) {
  cl <- match.call(); prep <- match.arg(prep)
  .local_seed(seed)
  x <- as.numeric(x); y <- as.numeric(y)
  if (length(x) != length(y)) stop("x and y must have the same length.", call. = FALSE)
  f <- switch(prep, rank = .rank_norm, scale = function(v) as.numeric(scale(v)), none = identity)
  x <- f(x); y <- f(y)
  a <- sqrt(29 / 31)
  X <- embed_delay(x, m, tau); Y <- embed_delay(y, m, tau); J <- embed_delay(a * x + y, m, tau)
  Z <- a * X + Y[sample.int(nrow(Y)), , drop = FALSE]
  keep <- seq(1L, nrow(X), by = max(1L, as.integer(downsample)))
  X <- X[keep, , drop = FALSE]; Y <- Y[keep, , drop = FALSE]; J <- J[keep, , drop = FALSE]; Z <- Z[keep, , drop = FALSE]
  n <- nrow(X); k_range <- sort(unique(as.integer(k_range))); kmax <- max(k_range)
  if (kmax >= n - 1) stop("k_range too large for the number of embedded points.", call. = FALSE)
  jit <- function(M) M + matrix(stats::rnorm(length(M), 0, 1e-10), nrow(M))
  mans <- list(X = X, Y = Y, J = J, Z = Z)
  Dk <- lapply(mans, function(M) knn_dist_cpp(.num_mat(jit(M)), kmax))
  rows <- list()
  for (k in k_range) {
    k2 <- max(1L, k %/% 2L)        # Farahmand et al. (2007): d = log(k / k2) / log(r_k / r_k2)
    L <- vapply(Dk, function(D) log(k / k2) / log(pmax(D[, k], 1e-300) / pmax(D[, k2], 1e-300)), numeric(n))
    L[!is.finite(L)] <- NA
    keep_i <- rep(TRUE, n)
    for (j in 1:4) {
      q <- stats::quantile(L[, j], c(eps, 1 - eps), na.rm = TRUE)
      keep_i <- keep_i & !is.na(L[, j]) & L[, j] >= q[1] & L[, j] <= q[2]
    }
    Lk <- L[keep_i, , drop = FALSE]; nk <- nrow(Lk)
    d <- colMeans(Lk); S <- stats::cov(Lk) * (2 * k / nk)
    lik <- .dc_likelihoods(d, S, c, bins)
    pA <- lik / sum(lik)
    pr <- c(x_to_y = pA[["A11"]], y_to_x = pA[["A21"]], bidirectional = pA[["A31"]],
            common_cause = pA[["A12"]] + pA[["A22"]] + pA[["A32"]],
            independent = pA[["A13"]] + pA[["A23"]] + pA[["A33"]])
    rows[[length(rows) + 1L]] <- data.frame(k = k, D_X = d[1], D_Y = d[2], D_J = d[3], D_Z = d[4],
                                            se_X = sqrt(S[1, 1]), se_Y = sqrt(S[2, 2]), se_J = sqrt(S[3, 3]),
                                            se_Z = sqrt(S[4, 4]), t(pr), n_used = nk)
  }
  by_k <- do.call(rbind, rows); rownames(by_k) <- NULL
  probs <- colMeans(by_k[, c("x_to_y", "y_to_x", "bidirectional", "common_cause", "independent")])
  structure(list(probs = probs, best = names(probs)[which.max(probs)], by_k = by_k,
                 dims = colMeans(by_k[, c("D_X", "D_Y", "D_J", "D_Z")]),
                 settings = list(m = m, tau = tau, k_range = k_range, prep = prep, eps = eps,
                                 downsample = downsample, n = n), call = cl),
            class = "dimcaus")
}

#' @export
print.dimcaus <- function(x, digits = 3, ...) {
  s <- x$settings
  cat("Bayesian dimensional causality (Benko et al., 2024)\n")
  cat(sprintf("  embedding m = %d, tau = %d, %s normalisation, n = %d points, k = %s\n",
              s$m, s$tau, s$prep, s$n, paste(range(s$k_range), collapse = "-")))
  cat("  mean intrinsic dimensions:  "); cat(sprintf("D_%s = %.2f", c("X", "Y", "J", "Z"), x$dims), sep = "  "); cat("\n")
  cat("  posterior probabilities of the causal cases (averaged over k):\n")
  lab <- c(x_to_y = "X -> Y", y_to_x = "X <- Y", bidirectional = "X <-> Y",
           common_cause = "hidden common cause", independent = "independent")
  print(data.frame(case = lab[names(x$probs)], probability = round(x$probs, digits),
                   best = ifelse(names(x$probs) == x$best, "<==", "")), row.names = FALSE, right = FALSE)
  invisible(x)
}

#' @rdname dimensional_causality
#' @param main title.
#' @param ... unused.
#' @export
plot.dimcaus <- function(x, main = NULL, ...) {
  B <- x$by_k
  op <- graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 1)); on.exit(graphics::par(op))
  cols <- c(X = "#C62828", Y = "#1565C0", J = "#212121", Z = "#F9A825")
  ylim <- range(c(B$D_X - B$se_X, B$D_Y - B$se_Y, B$D_J - B$se_J, B$D_Z - B$se_Z,
                  B$D_X + B$se_X, B$D_Y + B$se_Y, B$D_J + B$se_J, B$D_Z + B$se_Z), na.rm = TRUE)
  ylim[2] <- ylim[2] + 0.25 * diff(ylim)                    # room for the legend
  graphics::plot(NA, xlim = range(B$k), ylim = ylim, xlab = "neighbourhood size k",
                 ylab = "intrinsic dimension", main = if (is.null(main)) "Manifold dimensions" else main,
                 font.main = 2)
  for (nm in names(cols)) {
    dcol <- B[[paste0("D_", nm)]]; se <- B[[paste0("se_", nm)]]
    graphics::polygon(c(B$k, rev(B$k)), c(dcol - se, rev(dcol + se)),
                      col = grDevices::adjustcolor(cols[nm], 0.15), border = NA)
    graphics::lines(B$k, dcol, col = cols[nm], lwd = 2); graphics::points(B$k, dcol, col = cols[nm], pch = 16)
  }
  graphics::legend("topright", legend = c("X", "Y", "J (joint)", "Z (independent joint)"),
                   col = cols, lwd = 2, bty = "n", cex = 0.8)
  lab <- c("X -> Y", "X <- Y", "X <-> Y", "common\ncause", "indep.")
  bc <- ifelse(names(x$probs) == x$best, "#1B7A43", "#BFD7EA")
  bp <- graphics::barplot(x$probs, names.arg = rep("", length(lab)), col = bc, border = "grey40", ylim = c(0, 1.08),
                          ylab = "posterior probability", main = "Causal cases", font.main = 2)
  graphics::mtext(lab, side = 1, at = bp, line = 1.6, cex = 0.72, padj = 0)
  graphics::text(bp, x$probs, sprintf("%.2f", x$probs), pos = 3, cex = 0.75, col = "grey20", xpd = NA)
  invisible(x)
}

# rcausality driver: one dimensional-causality run per unordered pair --------
.rcausality_dc <- function(all, grid, causes, effects, alpha, adjust, m, tau, k_range, prep, eps,
                           downsample, prob_threshold, cl, cores) {
  up <- .unordered_pairs(grid, all)
  res <- .plapply(seq_len(nrow(up)), function(i)
    dimensional_causality(all[[up$v1[i]]], all[[up$v2[i]]], m = m, tau = tau, k_range = k_range,
                          prep = prep, eps = eps, downsample = downsample), cores = cores)
  rows <- list(); ins <- list(); dc <- list()
  prior_drive <- 2 / 9                      # P(A11) + P(A31) under the uniform prior over the 9 events
  for (i in seq_len(nrow(up))) {
    r <- res[[i]]; p <- r$probs; a <- up$v1[i]; b <- up$v2[i]
    for (dir in c("ab", "ba")) {
      p_dir <- if (dir == "ab") p[["x_to_y"]] else p[["y_to_x"]]
      post <- p_dir + p[["bidirectional"]]
      bf <- (post / (1 - post)) / (prior_drive / (1 - prior_drive))
      rows[[length(rows) + 1L]] <- data.frame(
        cause = if (dir == "ab") a else b, effect = if (dir == "ab") b else a, lag = NA_integer_,
        statistic = p_dir, df1 = NA_real_, df2 = NA_real_, p_value = NA_real_, bf10 = bf,
        effect_size = r$dims[["D_J"]] - r$dims[[if (dir == "ab") "D_X" else "D_Y"]], r2_gain = NA_real_,
        sign = NA_real_, n_eff = r$settings$n, post_prob = post, stringsAsFactors = FALSE)
    }
    ins[[i]] <- data.frame(var1 = a, var2 = b, cmi = NA_real_, p_value = NA_real_, p_adj = NA_real_,
                           significant = r$best == "common_cause", stringsAsFactors = FALSE)
    dc[[i]] <- data.frame(var1 = a, var2 = b, D_X = r$dims[["D_X"]], D_Y = r$dims[["D_Y"]],
                          D_J = r$dims[["D_J"]], D_Z = r$dims[["D_Z"]], t(p), best = r$best,
                          stringsAsFactors = FALSE)
  }
  pairs <- do.call(rbind, rows)
  pairs$p_adj <- NA_real_; pairs$evidence <- bf_evidence(pairs$bf10)
  pairs$significant <- pairs$post_prob > prob_threshold
  pairs <- pairs[, c("cause", "effect", "lag", "statistic", "df1", "df2", "p_value",
                     "p_adj", "bf10", "post_prob", "evidence", "effect_size",
                     "r2_gain", "sign", "n_eff", "significant")]
  rownames(pairs) <- NULL
  .new_causres(pairs, series = names(all), causes = causes, effects = effects,
               method = "dc", conditional = FALSE, alpha = alpha, bf_threshold = NA_real_,
               adjust = "none", n = nrow(all), call = cl, data = all, instant = do.call(rbind, ins),
               dc = do.call(rbind, dc), prob_threshold = prob_threshold)
}
