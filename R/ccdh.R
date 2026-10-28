# CCDH: complexity-based causal discovery of hidden common cause
# (Stippinger, Varga, Benko, Fabo, Eross, Telcs & Somogyvari, Chaos, Solitons
# & Fractals 176 (2023) 114054, doi:10.1016/j.chaos.2023.114054).
#
# Series are rank-order discretised after time-delay embedding (symbols
# 1..m!), compressed by non-sequential recursive pair substitution (NSRPS),
# and the compression rates cr_k of a series are compared with those of its
# "skewed joint" with the other series (Eqs. 7-8): if the other series can be
# reconstructed from the embedding (it drives), the joint is not more complex
# than the series itself (Eq. 10); a hidden common cause makes both joints
# more complex (Proposition 1).  The paper leaves the significance testing
# open ("surrogate data and test of significance are necessary"); here the
# excess complexity of a joint is compared with surrogates in which the added
# series is circularly shifted (no shared information).

.ordinal_codes <- function(E) {
  E <- as.matrix(E)
  if (ncol(E) > 12L) stop("embedding dimension too large for rank-order coding (at most 12).", call. = FALSE)
  code <- .lehmer_codes(E)                 # canonical pattern codes (ties broken by position)
  match(code, unique(code))                # renumbered by first appearance for the compressor
}

# compression-rate curve (k = 1..K) of the rank-order symbols of matrix E
.cr_curve <- function(E, K) nsrps_cr_cpp(.ordinal_codes(E), as.integer(K))

#' CCDH: complexity-based causal discovery with hidden common cause
#'
#' Compares the NSRPS compression rates of the rank-order patterns of the
#' embedded series `X` and `Y` with those of their skewed joints
#' `J_xxy = (x[t], ..., x[t-m+2], y[t-m+1])` (the oldest coordinate replaced
#' by the lagged other series, Appendix of the paper) and `J_yyx`.
#' If `y` drives `x`, the lagged `y` is a function of the embedded `x` (Takens),
#' so the joint is no more complex than `X` (`cr_k(J_xxy) = cr_k(X)`); if
#' `x` drives `y` the same holds for `J_yyx` and `Y`; both equalities hold
#' for bidirectional coupling, neither for a hidden common cause.
#'
#' The excess complexity `Delta_x = mean_k (cr_k(J_xxy) - cr_k(X))` is
#' compared with the excess obtained with `n_surr` circularly shifted
#' copies of `y` (no shared information): `p_share` is the surrogate
#' p-value that `y` carries information about `X` at all, and the
#' reconstruction index `s_x = Delta_x / mean(Delta_x*)` is 0 when `y` is
#' perfectly reconstructible from `X` and 1 when `y` is as uninformative as a
#' surrogate. A direction is accepted when the shared information is
#' significant and the index is below `theta`.
#'
#' @param x,y numeric vectors.
#' @param m embedding dimension (paper: 3 for maps, 4-5 for flows/EEG).
#' @param tau embedding delay.
#' @param K number of NSRPS steps whose compression rates are averaged.
#' @param lag shift applied to `y` (`y[t + lag]` is paired with `x[t]`); the
#'   paper aligns the series at the lag of maximal shared information.
#' @param n_surr number of circular-shift surrogates.
#' @param theta reconstruction-index threshold below which a joint is
#'   considered as complex as the single series (0.45 by default: on the
#'   coupled logistic maps of the paper, n = 1000, m = 3, the index of a
#'   reconstructible driver is typically 0.15-0.4, that of a hidden common
#'   driver 0.5-0.85, that of an unrelated series about 1).
#' @param alpha significance level for the shared-information test.
#' @param seed optional RNG seed.
#' @param cores CPU cores for the surrogates (see [Rcausality-parallel]).
#' @return An object of class `ccdh` with `cr` (K x 4 matrix of compression
#'   rates: X, Y, J_xxy, J_yyx), `delta` (excess complexity of the two
#'   joints), `delta_surr` (surrogate means), `index` (`s_x`, `s_y`),
#'   `p_share`, `verdict` (one of `x_to_y`, `y_to_x`, `bidirectional`,
#'   `common_cause`, `independent`), `settings`. Methods: `print`, `plot`
#'   (compression rates against the number of steps, Fig. 1 of the paper).
#' @references
#' Stippinger, M., Varga, B., Benko, Z., Fabo, D., Eross, L., Telcs, A. and Somogyvari, Z. (2023). CCDH: Complexity based causal discovery of hidden common cause in time series. *Chaos, Solitons & Fractals*, 176, 114054. \doi{10.1016/j.chaos.2023.114054}
#'
#' Kathpalia, A. and Nagaraj, N. (2019). Data-based intervention approach for Complexity-Causality measure. *PeerJ Computer Science*, 5, e196. \doi{10.7717/peerj-cs.196}
#'
#' Nagaraj, N., Balasubramanian, K. and Dey, S. (2013). A new complexity measure for time series analysis and classification. *The European Physical Journal Special Topics*, 222, 847-860. \doi{10.1140/epjst/e2013-01888-9}
#' @export
#' @examples
#' # x0 drives x1 in the coupled logistic maps.  The driver is reconstructible
#' # from the driven series: the joint of x1 with lagged x0 compresses as well
#' # as x1 alone (small excess complexity, index close to 0), while the joint of
#' # x0 with lagged x1 does not.  The verdict is x_to_y.
#' d <- sim_logistic3(1000, "driver", b = 0.4, c = 0.4, seed = 1)
#' r <- ccdh(d$x0, d$x1, m = 3, n_surr = 20, seed = 1)
#' r
#' \donttest{
#' plot(r)      # compression-rate curves of the four symbol sequences (Fig. 1 of the paper)
#' # x1 and x2 are both driven by x0: with x0 left out they share a hidden
#' # common cause, and both joints are more complex than the single series
#' ccdh(d$x1, d$x2, m = 3, n_surr = 50, seed = 1)$verdict
#' }
ccdh <- function(x, y, m = 3L, tau = 1L, K = 10L, lag = 0L, n_surr = 50L, theta = 0.45, alpha = 0.05,
                 seed = NULL, cores = getOption("Rcausality.cores", 1L)) {
  cl <- match.call()
  .local_seed(seed)
  x <- as.numeric(x); y <- as.numeric(y); n <- length(x)
  if (length(y) != n) stop("x and y must have the same length.", call. = FALSE)
  lag <- as.integer(lag)
  if (lag != 0L) {                       # pair x[t] with y[t + lag]
    i <- seq_len(n); j <- i + lag; ok <- j >= 1 & j <= n
    x <- x[i[ok]]; y <- y[j[ok]]; n <- length(x)
  }
  nz <- function(v) (v - min(v)) / (max(v) - min(v) + 1e-300)
  x <- nz(x); y <- nz(y)
  m <- as.integer(m); tau <- as.integer(tau); K <- as.integer(K)
  # skewed joint of Appendix 7.1-7.2 ("low dimensional joint"): the oldest
  # coordinate of the embedding of `a` is replaced by the lagged `b`, so the
  # joint has the same dimension m and the same alphabet (m! patterns) as the
  # single series:  (a[t], a[t-tau], ..., a[t-(m-2)tau], b[t-(m-1)tau])
  # (columns oldest first, as in embed_delay(), so that with b = a the joint
  # is exactly the embedding of a, ties included)
  joint <- function(a, b) {
    t_idx <- ((m - 1) * tau + 1):n
    cbind(b[t_idx - (m - 1) * tau], do.call(cbind, lapply((m - 2):0, function(l) a[t_idx - l * tau])))
  }
  cr <- cbind(X = .cr_curve(embed_delay(x, m, tau), K), Y = .cr_curve(embed_delay(y, m, tau), K),
              J_xxy = .cr_curve(joint(x, y), K), J_yyx = .cr_curve(joint(y, x), K))
  delta <- c(x = mean(cr[, "J_xxy"] - cr[, "X"]), y = mean(cr[, "J_yyx"] - cr[, "Y"]))
  lo <- max(2L * m * tau + 1L, floor(n / 10)); hi <- n - lo
  if (hi <= lo) stop("Series too short for the surrogates.", call. = FALSE)
  shift <- function(v, s) v[((seq_len(n) + s - 1L) %% n) + 1L]
  surr <- do.call(rbind, .plapply(seq_len(n_surr), function(b) {
    s1 <- sample(lo:hi, 1L); s2 <- sample(lo:hi, 1L)
    c(x = mean(.cr_curve(joint(x, shift(y, s1)), K) - cr[, "X"]),
      y = mean(.cr_curve(joint(y, shift(x, s2)), K) - cr[, "Y"]))
  }, cores = cores))
  p_share <- c(x = (1 + sum(surr[, "x"] <= delta[["x"]])) / (n_surr + 1),
               y = (1 + sum(surr[, "y"] <= delta[["y"]])) / (n_surr + 1))
  ds <- colMeans(surr)
  index <- pmax(delta / pmax(ds, 1e-12), 0)
  rec_x <- p_share[["x"]] < alpha && index[["x"]] <= theta       # y reconstructible from X: y -> x
  rec_y <- p_share[["y"]] < alpha && index[["y"]] <= theta       # x reconstructible from Y: x -> y
  shared <- p_share[["x"]] < alpha || p_share[["y"]] < alpha
  verdict <- if (rec_x && rec_y) "bidirectional" else if (rec_x) "y_to_x" else if (rec_y) "x_to_y" else
    if (shared) "common_cause" else "independent"
  structure(list(cr = cr, delta = delta, delta_surr = ds, index = index, p_share = p_share, surr = surr,
                 verdict = verdict, n = n,
                 settings = list(m = m, tau = tau, K = K, lag = lag, n_surr = n_surr, theta = theta, alpha = alpha),
                 call = cl), class = "ccdh")
}

#' @export
print.ccdh <- function(x, digits = 3, ...) {
  s <- x$settings
  cat("CCDH: complexity-based causal discovery (Stippinger et al.)\n")
  cat(sprintf("  embedding m = %d, tau = %d, %d NSRPS steps, lag = %d, n = %d, %d surrogates\n",
              s$m, s$tau, s$K, s$lag, x$n, s$n_surr))
  cat("  mean compression rates:  "); cat(sprintf("%s = %.3f", colnames(x$cr), colMeans(x$cr)), sep = "  "); cat("\n")
  tab <- data.frame(joint = c("J_xxy vs X  (y reconstructible from X: y -> x)", "J_yyx vs Y  (x reconstructible from Y: x -> y)"),
                    excess = round(x$delta, digits), surrogate = round(x$delta_surr, digits),
                    index = round(x$index, digits), p.share = .fmt_p(x$p_share),
                    reconstructible = ifelse(x$p_share < s$alpha & x$index <= s$theta, "yes", "no"))
  print(tab, row.names = FALSE, right = FALSE)
  cat(sprintf("  verdict: %s   (index <= %g and p.share < %g => reconstructible)\n", x$verdict, s$theta, s$alpha))
  invisible(x)
}

#' @rdname ccdh
#' @param main title.
#' @param ... unused.
#' @export
plot.ccdh <- function(x, main = NULL, ...) {
  cr <- x$cr; K <- nrow(cr)
  cols <- c(X = "#C62828", Y = "#1565C0", J_xxy = "#E57373", J_yyx = "#64B5F6")
  op <- graphics::par(mar = c(4.5, 4.5, 3.5, 1)); on.exit(graphics::par(op))
  graphics::plot(NA, xlim = c(0, K), ylim = c(min(cr) * 0.9, 1), xlab = "number of NSRPS compression steps",
                 ylab = "compression rate", main = if (is.null(main)) sprintf("CCDH compression rates (verdict: %s)", x$verdict) else main,
                 font.main = 2)
  for (nm in colnames(cr)) {
    graphics::lines(0:K, c(1, cr[, nm]), col = cols[nm], lwd = 2, lty = if (grepl("^J", nm)) 2 else 1)
    graphics::points(0:K, c(1, cr[, nm]), col = cols[nm], pch = if (grepl("^J", nm)) 1 else 16)
  }
  graphics::legend("topright", legend = c("X", "Y", "joint J_xxy (X + y[t-m])", "joint J_yyx (Y + x[t-m])"),
                   col = cols, lwd = 2, lty = c(1, 1, 2, 2), pch = c(16, 16, 1, 1), bty = "n", cex = 0.8)
  graphics::mtext("joint on its single series = the added series is reconstructible (it drives); both above = hidden common cause",
                  side = 3, line = 0.2, cex = 0.7, col = "grey40")
  invisible(x)
}

# rcausality driver ------------------------------------------------------------
.rcausality_ccdh <- function(all, grid, causes, effects, alpha, adjust, m, tau, K, lag, n_surr, theta,
                             cl, cores) {
  up <- .unordered_pairs(grid, all)
  sc <- .split_cores(nrow(up), cores)
  res <- .plapply(seq_len(nrow(up)), function(i)
    ccdh(all[[up$v1[i]]], all[[up$v2[i]]], m = m, tau = tau, K = K, lag = lag, n_surr = n_surr,
         theta = theta, alpha = alpha, cores = sc$inner), cores = sc$outer)
  rows <- list(); ins <- list(); cc <- list(); nulls <- list()
  for (i in seq_len(nrow(up))) {
    r <- res[[i]]; a <- up$v1[i]; b <- up$v2[i]
    # a -> b: b's joint with a (J_yyx) is as simple as b  <=>  index y small
    for (dir in c("ab", "ba")) {
      which_ <- if (dir == "ab") "y" else "x"
      # evidence of reconstructibility: the excess complexity is SMALL compared with the surrogates
      nulls[[paste(if (dir == "ab") a else b, "->", if (dir == "ab") b else a)]] <-
        list(statistic = r$delta[[which_]], null = r$surr[, which_], alternative = "less")
      rows[[length(rows) + 1L]] <- data.frame(
        cause = if (dir == "ab") a else b, effect = if (dir == "ab") b else a, lag = r$settings$lag,
        statistic = 1 - r$index[[which_]], df1 = NA_real_, df2 = NA_real_, p_value = r$p_share[[which_]],
        bf10 = NA_real_, effect_size = r$delta[[which_]], r2_gain = NA_real_, sign = NA_real_,
        n_eff = r$n, index = r$index[[which_]], stringsAsFactors = FALSE)
    }
    ins[[i]] <- data.frame(var1 = a, var2 = b, cmi = NA_real_, p_value = NA_real_, p_adj = NA_real_,
                           significant = r$verdict == "common_cause", stringsAsFactors = FALSE)
    cc[[i]] <- data.frame(var1 = a, var2 = b, cr_x = mean(r$cr[, "X"]), cr_y = mean(r$cr[, "Y"]),
                          cr_jxxy = mean(r$cr[, "J_xxy"]), cr_jyyx = mean(r$cr[, "J_yyx"]),
                          index_x = r$index[["x"]], index_y = r$index[["y"]], verdict = r$verdict,
                          stringsAsFactors = FALSE)
  }
  pairs <- do.call(rbind, rows)
  pairs$p_adj <- if (adjust == "none") pairs$p_value else stats::p.adjust(pairs$p_value, method = adjust)
  pairs$post_prob <- NA_real_; pairs$evidence <- NA_character_
  pairs$significant <- pairs$p_adj < alpha & pairs$index <= theta
  pairs <- pairs[, c("cause", "effect", "lag", "statistic", "df1", "df2", "p_value",
                     "p_adj", "bf10", "post_prob", "evidence", "effect_size",
                     "r2_gain", "sign", "n_eff", "significant")]
  rownames(pairs) <- NULL
  out <- .new_causres(pairs, series = names(all), causes = causes, effects = effects,
               method = "ccdh", conditional = FALSE, alpha = alpha, bf_threshold = NA_real_,
               adjust = adjust, n = nrow(all), call = cl, data = all, instant = do.call(rbind, ins),
               ccdh = do.call(rbind, cc), prob_threshold = theta)
  out$null <- nulls
  out
}
