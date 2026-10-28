# Baseline methods: transfer entropy (Schreiber 2000) with symbolic binning
# and convergent cross mapping (Sugihara et al. 2012).

# ---- transfer entropy ----------------------------------------------------------
# symbolic coding: discrete series (at most `bins` distinct values) keep their
# values as symbols; otherwise equal-frequency bins of the (average) ranks, so
# that tied values share a bin and heavy ties cannot collapse the coding
.sym_bins <- function(v, bins) {
  u <- sort(unique(v))
  if (length(u) <= bins) return(match(v, u) - 1L)
  r <- rank(v, ties.method = "average")
  s <- findInterval(r, seq(0, length(v), length.out = bins + 1), rightmost.closed = TRUE) - 1L
  as.integer(pmin(pmax(s, 0L), bins - 1L))
}

.te_sym <- function(sx, sy, k, l, bins, t_idx) {
  # TE_{x -> y} = H(y[t] | y past l) - H(y[t] | y past l, x past k)
  ty <- sy[t_idx]
  Hy <- do.call(cbind, lapply(seq_len(l), function(j) sy[t_idx - j]))
  Hxy <- cbind(Hy, do.call(cbind, lapply(seq_len(k), function(j) sx[t_idx - j])))
  storage.mode(Hy) <- "integer"; storage.mode(Hxy) <- "integer"
  cond_entropy_sym_cpp(ty, Hy, bins, bins) - cond_entropy_sym_cpp(ty, Hxy, bins, bins)
}

#' Transfer entropy (binned) with a surrogate test
#'
#' Schreiber's transfer entropy `TE(x -> y) = H(y[t] | y[t-1..t-l]) -
#' H(y[t] | y[t-1..t-l], x[t-1..t-k])` in bits, estimated from the
#' quantile-binned (symbolic) series. The null distribution comes from
#' circular-shift surrogates of `x` (its own dynamics kept, the coupling
#' destroyed); the *effective* transfer entropy is `TE - mean(TE of the
#' surrogates)` (Marschinski and Kantz, 2002). A kNN (continuous) transfer
#' entropy is available as `rcausality(method = "cmi")`, since the transfer
#' entropy is the conditional mutual information `I(x past; y[t] | y past)`.
#'
#' @param x,y numeric vectors (`x` the candidate cause).
#' @param k,l history lengths of the cause and of the effect.
#' @param bins number of quantile bins.
#' @param n_surr number of surrogates.
#' @param seed optional RNG seed.
#' @param cores CPU cores for the surrogates (see [Rcausality-parallel]).
#' @return A list with `te`, `te_eff` (effective TE), `p_value`, `null`.
#' @references
#' Schreiber, T. (2000). Measuring information transfer. *Physical Review Letters*, 85(2), 461-464. \doi{10.1103/PhysRevLett.85.461}
#'
#' Marschinski, R. and Kantz, H. (2002). Analysing the information flow between financial time series: an improved estimator for transfer entropy. *The European Physical Journal B*, 30, 275-281. \doi{10.1140/epjb/e2002-00379-2}
#' @export
#' @examples
#' # x drives y in the simulated VAR.  The transfer entropy x -> y exceeds
#' # its circular-shift surrogates (small p-value); the reverse direction sits
#' # at the surrogate level (large p-value), so there is no feedback.
#' d <- sim_var(500, "x_to_y", seed = 1)
#' te_xy <- transfer_entropy(d$x, d$y, n_surr = 50, seed = 1)
#' te_xy$te         # transfer entropy in bits
#' te_xy$te_eff     # effective TE: the surrogate mean subtracted (bias correction)
#' te_xy$p_value    # small: x -> y
#' transfer_entropy(d$y, d$x, n_surr = 50, seed = 1)$p_value   # large: no y -> x
transfer_entropy <- function(x, y, k = 1L, l = 1L, bins = 4L, n_surr = 100L, seed = NULL,
                             cores = getOption("Rcausality.cores", 1L)) {
  .local_seed(seed)
  x <- as.numeric(x); y <- as.numeric(y); n <- length(x)
  if (length(y) != n) stop("x and y must have the same length.", call. = FALSE)
  k <- as.integer(k); l <- as.integer(l)
  if (k < 1L || l < 1L) stop("k and l (the cause's and the effect's history lengths) must be at least 1.", call. = FALSE)
  sx <- .sym_bins(x, bins); sy <- .sym_bins(y, bins)
  nb <- max(sx, sy) + 1L
  L <- max(k, l); t_idx <- (L + 1):n
  te <- .te_sym(sx, sy, k, l, nb, t_idx)
  lo <- max(2L * L + 1L, floor(n / 10)); hi <- n - lo
  if (hi <= lo) stop("Series too short for the surrogate test.", call. = FALSE)
  null <- unlist(.plapply(seq_len(n_surr), function(b) {
    s <- sample(lo:hi, 1L); .te_sym(sx[((seq_len(n) + s - 1L) %% n) + 1L], sy, k, l, nb, t_idx)
  }, cores = cores))
  list(te = te, te_eff = te - mean(null), p_value = (1 + sum(null >= te)) / (n_surr + 1), null = null,
       settings = list(k = k, l = l, bins = bins))
}

# ---- convergent cross mapping ------------------------------------------------
#' Convergent cross mapping (Sugihara et al., 2012)
#'
#' Tests whether `x` causes `y` by cross-mapping `x` from the time-delay
#' embedding of `y` ("`y` xmap `x`"): if `x` drives `y`, the states of `x` are
#' encoded in `y`'s manifold and the simplex-projection estimate of `x` from
#' `y`'s neighbours improves (converges) as the library size grows. The
#' cross-map skill `rho(L)` (Pearson correlation between `x` and its
#' estimate) is computed for a grid of library sizes with random libraries;
#' convergence is tested with Kendall's tau of `rho` against `L`, and the
#' skill at the largest library against circular-shift surrogates of `x`.
#'
#' @param x,y numeric vectors (`x` the candidate cause, `y` the effect whose
#'   manifold is used).
#' @param E embedding dimension; `tau` the delay.
#' @param tau embedding delay.
#' @param lib_sizes library sizes (default: 10 values from `E + 2` to the
#'   number of embedded points).
#' @param n_samples random libraries per size.
#' @param n_surr number of surrogates for the skill test.
#' @param seed optional RNG seed.
#' @param cores CPU cores for the surrogates (see [Rcausality-parallel]).
#' @return A list with `rho` (data frame: `lib_size`, `rho`), `rho_max`,
#'   `p_value` (surrogate test of `rho_max`), `p_convergence` (Kendall test of
#'   increasing skill), `converging`.
#' @references
#' Sugihara, G., May, R., Ye, H., Hsieh, C., Deyle, E., Fogarty, M. and Munch, S. (2012). Detecting causality in complex ecosystems. *Science*, 338(6106), 496-500. \doi{10.1126/science.1227079}
#'
#' Sugihara, G. and May, R. M. (1990). Nonlinear forecasting as a way of distinguishing chaos from measurement error in time series. *Nature*, 344, 734-741. \doi{10.1038/344734a0}
#' @export
#' @examples
#' # Coupled logistic maps where x drives y.  Because the driven series y
#' # encodes the state of its driver, "y xmap x" (estimating x from the
#' # neighbours on y's manifold) is skilful and improves with the library size.
#' d <- sim_logistic(300, "x_to_y", seed = 1)
#' ccm(d$x, d$y, E = 2, n_surr = 10, seed = 1)$rho_max   # skill at the largest library
#' \donttest{
#' r <- ccm(d$x, d$y, E = 2, n_surr = 50, seed = 1)
#' r$rho              # cross-map skill against the library size (should increase)
#' r$p_convergence    # Kendall test of increasing skill
#' r$p_value          # skill at the largest library against circular-shift surrogates
#' r$converging       # both conditions of a causal link
#' ccm(d$y, d$x, E = 2, n_surr = 50, seed = 1)$rho_max   # the reverse cross map is weaker
#' }
ccm <- function(x, y, E = 3L, tau = 1L, lib_sizes = NULL, n_samples = 20L, n_surr = 100L, seed = NULL,
                cores = getOption("Rcausality.cores", 1L)) {
  .local_seed(seed)
  x <- as.numeric(x); y <- as.numeric(y); n <- length(x)
  if (length(y) != n) stop("x and y must have the same length.", call. = FALSE)
  E <- as.integer(E); tau <- as.integer(tau)
  My <- .num_mat(embed_delay(y, E, tau)); ne <- nrow(My)
  xt <- x[((E - 1L) * tau + 1L):n]                      # x aligned with the embedded y
  if (is.null(lib_sizes)) lib_sizes <- unique(round(exp(seq(log(E + 2L), log(ne), length.out = 10L))))
  lib_sizes <- sort(unique(as.integer(lib_sizes[lib_sizes >= E + 2L & lib_sizes <= ne])))
  if (!length(lib_sizes)) stop("Series too short for the embedding.", call. = FALSE)
  skill <- function(xv, L) {
    mean(vapply(seq_len(if (L >= ne) 1L else n_samples), function(s) {   # the full library is unique
      lib <- if (L >= ne) seq_len(ne) else sort(sample.int(ne, L))
      pr <- simplex_predict_cpp(My, lib, lib, xv, E + 1L)
      r <- suppressWarnings(stats::cor(xv[lib], pr)); if (is.na(r)) 0 else r
    }, 0))
  }
  rho <- vapply(lib_sizes, function(L) skill(xt, L), 0)
  rho_max <- rho[length(rho)]
  kt <- if (length(rho) >= 3) suppressWarnings(stats::cor.test(lib_sizes, rho, method = "kendall", alternative = "greater")$p.value) else NA_real_
  lo <- max(2L * E * tau + 1L, floor(ne / 10)); hi <- ne - lo
  if (hi <= lo) stop("Series too short for the surrogate test.", call. = FALSE)
  null <- unlist(.plapply(seq_len(n_surr), function(b) {
    s <- sample(lo:hi, 1L); skill(xt[((seq_len(ne) + s - 1L) %% ne) + 1L], lib_sizes[length(lib_sizes)])
  }, cores = cores))
  list(rho = data.frame(lib_size = lib_sizes, rho = rho), rho_max = rho_max,
       p_value = (1 + sum(null >= rho_max)) / (n_surr + 1), p_convergence = kt,
       converging = !is.na(kt) && kt < 0.05 && rho[length(rho)] > rho[1], null = null,
       settings = list(E = E, tau = tau))
}

# ---- rcausality drivers ----------------------------------------------------------
.rcausality_te <- function(all, grid, causes, effects, alpha, adjust, conditional, lags, history, bins,
                           n_perm, cl, cores) {
  sc <- .split_cores(nrow(grid), cores); pair_cores <- sc$outer; inner_cores <- sc$inner
  res <- .plapply(seq_len(nrow(grid)), function(i) {
    cs <- grid$cause[i]; ef <- grid$effect[i]
    x <- all[[cs]]; y <- all[[ef]]
    # lag tau: the cause at t - tau only (x shifted by tau - 1, cause history 1)
    nn <- length(x)
    tests <- lapply(lags, function(tau) transfer_entropy(x[seq_len(nn - tau + 1L)], y[tau:nn], k = 1L, l = history[i],
                                                         bins = bins, n_surr = n_perm, cores = inner_cores))
    bl <- data.frame(cause = cs, effect = ef, lag = lags, statistic = vapply(tests, `[[`, 0, "te"),
                     te_eff = vapply(tests, `[[`, 0, "te_eff"), p = vapply(tests, `[[`, 0, "p_value"))
    null_max <- do.call(pmax, lapply(tests, `[[`, "null"))      # max-statistic test over the lags
    b <- which.max(bl$statistic)
    list(row = data.frame(cause = cs, effect = ef, lag = bl$lag[b], statistic = bl$statistic[b],
                          df1 = NA_real_, df2 = NA_real_, p_value = (1 + sum(null_max >= bl$statistic[b])) / (length(null_max) + 1),
                          bf10 = NA_real_, effect_size = bl$te_eff[b], r2_gain = NA_real_, sign = NA_real_,
                          n_eff = length(x) - max(lags, history[i]), stringsAsFactors = FALSE),
         null = list(statistic = bl$statistic[b], null = null_max, alternative = "greater"),
         by_lag = bl[, c("cause", "effect", "lag", "statistic", "p")])
  }, cores = pair_cores)
  pairs <- do.call(rbind, lapply(res, `[[`, "row"))
  by_lag <- do.call(rbind, lapply(res, `[[`, "by_lag")); rownames(by_lag) <- NULL
  pairs$p_adj <- if (adjust == "none") pairs$p_value else stats::p.adjust(pairs$p_value, method = adjust)
  pairs$post_prob <- NA_real_; pairs$evidence <- NA_character_
  pairs$significant <- pairs$p_adj < alpha
  pairs <- pairs[, c("cause", "effect", "lag", "statistic", "df1", "df2", "p_value",
                     "p_adj", "bf10", "post_prob", "evidence", "effect_size",
                     "r2_gain", "sign", "n_eff", "significant")]
  rownames(pairs) <- NULL
  out <- .new_causres(pairs, series = names(all), causes = causes, effects = effects,
                      method = "te", conditional = FALSE, alpha = alpha, bf_threshold = NA_real_,
                      adjust = adjust, n = nrow(all), call = cl, data = all, by_lag = by_lag)
  out$null <- stats::setNames(lapply(res, `[[`, "null"), paste(pairs$cause, "->", pairs$effect))
  out
}

.rcausality_ccm <- function(all, grid, causes, effects, alpha, adjust, m, tau, n_perm, cl, cores) {
  sc <- .split_cores(nrow(grid), cores); pair_cores <- sc$outer; inner_cores <- sc$inner
  res <- .plapply(seq_len(nrow(grid)), function(i) {
    cs <- grid$cause[i]; ef <- grid$effect[i]
    r <- ccm(all[[cs]], all[[ef]], E = m, tau = tau, n_surr = n_perm, cores = inner_cores)
    list(null = list(statistic = r$rho_max, null = r$null, alternative = "greater"),
         row = data.frame(cause = cs, effect = ef, lag = NA_integer_, statistic = r$rho_max,
                          df1 = NA_real_, df2 = NA_real_, p_value = r$p_value, bf10 = NA_real_,
                          effect_size = r$rho_max, r2_gain = r$rho_max^2, sign = NA_real_,
                          n_eff = nrow(r$rho) * 0 + length(all[[cs]]), converging = r$converging,
                          p_convergence = r$p_convergence, stringsAsFactors = FALSE),
         by_lag = data.frame(cause = cs, effect = ef, lag = r$rho$lib_size, statistic = r$rho$rho, p = NA_real_))
  }, cores = pair_cores)
  pairs <- do.call(rbind, lapply(res, `[[`, "row"))
  by_lag <- do.call(rbind, lapply(res, `[[`, "by_lag")); rownames(by_lag) <- NULL
  pairs$p_adj <- if (adjust == "none") pairs$p_value else stats::p.adjust(pairs$p_value, method = adjust)
  pairs$post_prob <- NA_real_; pairs$evidence <- NA_character_
  pairs$significant <- pairs$p_adj < alpha & pairs$converging
  conv <- pairs[, c("cause", "effect", "converging", "p_convergence")]
  pairs <- pairs[, c("cause", "effect", "lag", "statistic", "df1", "df2", "p_value",
                     "p_adj", "bf10", "post_prob", "evidence", "effect_size",
                     "r2_gain", "sign", "n_eff", "significant")]
  rownames(pairs) <- NULL
  out <- .new_causres(pairs, series = names(all), causes = causes, effects = effects,
                      method = "ccm", conditional = FALSE, alpha = alpha, bf_threshold = NA_real_,
                      adjust = adjust, n = nrow(all), call = cl, data = all, by_lag = by_lag, ccm = conv)
  out$null <- stats::setNames(lapply(res, `[[`, "null"), paste(pairs$cause, "->", pairs$effect))
  out
}
