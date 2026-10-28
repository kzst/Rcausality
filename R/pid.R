# Empirical nonnegative finite-resolution MAR-PID for continuous variables
# (Telcs, Entropy 2026, 28, 641) with the Blackwell-order decomposition of
# Mages, Anastasiadi & Rohner (Entropy 2024, 26, 424), two-source case.
#
# A binary target C and two finite-valued sources V1, V2.  Each source is a
# binary-input channel kappa_S(s | t) = P(S = s | C = t) whose zonogon is the
# Minkowski sum of the column segments; the Blackwell join of two channels is
# the convex hull of the union of their zonogons (its upper boundary = the
# channel obtainable from either source).  With i(kappa) the mutual
# information of a channel for the target bias p:
#   union   = i(kappa_1 v kappa_2)
#   unique_1 = union - i(kappa_2),  unique_2 = union - i(kappa_1)
#   redundant = i(kappa_1) + i(kappa_2) - union
#   synergy  = i(kappa_{12}) - union       (joint source (V1, V2))
# all nonnegative (zonogon inclusion-exclusion).  Continuous variables are
# represented by recursive quantile (median) binarisation; target bits are
# decomposed one by one and the atoms summed over the resolved target bits.

# recursive quantile binarisation (Sect. 3): depth-k code of a numeric vector
.rq_bits <- function(v, k) {
  n <- length(v); B <- matrix(0L, n, k); cell <- rep(1L, n)
  r <- rank(v, ties.method = "average")              # ties stay on the same side of every split
  for (j in seq_len(k)) {
    for (cl in unique(cell)) {
      i <- which(cell == cl); if (!length(i)) next
      u <- unique(r[i]); if (length(u) < 2L) next
      m <- stats::median(u); B[i, j] <- as.integer(r[i] > m)
    }
    cell <- 2L * cell + B[, j]
  }
  B
}
.bits_code <- function(B) { if (!ncol(B)) return(rep(0L, nrow(B))); as.integer(B %*% (2L^(seq_len(ncol(B)) - 1L))) }

# channel columns v_s = (P(s | T = 1), P(s | T = 0)) from counts
.channel <- function(s, t) {
  tb <- table(factor(s), factor(t, levels = 0:1))
  cbind(p1 = tb[, "1"] / max(1, sum(tb[, "1"])), p0 = tb[, "0"] / max(1, sum(tb[, "0"])))
}
# Boundary chain of the zonogon of a channel in the (cum P(s|T=1), cum P(s|T=0))
# plane: the columns sorted by decreasing likelihood ratio p1/p0 give the
# chain from (0, 0) to (1, 1) with increasing slope, i.e. the LOWER (convex)
# boundary; the zonogon is its reflection through (1/2, 1/2) glued to it.
.zono_upper <- function(cols) {
  cols <- cols[rowSums(cols) > 0, , drop = FALSE]
  o <- order(-ifelse(cols[, 2] > 0, cols[, 1] / cols[, 2], Inf), -cols[, 1])
  rbind(c(0, 0), apply(cols[o, , drop = FALSE], 2, cumsum))          # vertices (x = p1 cum, y = p0 cum)
}
# Blackwell join of two channels (Mages et al. 2024): the smallest zonogon
# containing both, i.e. the convex hull of the union of the two zonogons.
# Its lower boundary is the lower convex envelope of the union of the two
# lower chains (the join is as informative as the more informative channel
# along every direction, never less).
.zono_join <- function(V1, V2) {
  P <- unique(rbind(V1, V2)); P <- P[order(P[, 1], P[, 2]), , drop = FALSE]
  hull <- list()
  for (i in seq_len(nrow(P))) {
    while (length(hull) >= 2) {
      a <- hull[[length(hull) - 1]]; b <- hull[[length(hull)]]; c <- P[i, ]
      cross <- (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])
      if (cross <= 1e-12) hull[[length(hull)]] <- NULL else break    # keep the lower (convex) chain
    }
    hull[[length(hull) + 1]] <- P[i, ]
  }
  do.call(rbind, hull)
}
# mutual information (bits) of the channel given by the edges of an upper boundary
.mi_upper <- function(V, p) {
  E <- diff(V); E <- E[rowSums(E) > 0, , drop = FALSE]
  s1 <- E[, 1]; s0 <- E[, 2]; mix <- p * s1 + (1 - p) * s0
  t1 <- ifelse(s1 > 0, p * s1 * log2(s1 / mix), 0); t0 <- ifelse(s0 > 0, (1 - p) * s0 * log2(s0 / mix), 0)
  sum(t1 + t0)
}

# two-source MAR-PID of one binary target
.marpid_bit <- function(ct, s1, s2) {
  p <- mean(ct == 1)
  if (p <= 0 || p >= 1) return(c(redundant = 0, unique_1 = 0, unique_2 = 0, synergy = 0, total = 0))
  Z1 <- .zono_upper(.channel(s1, ct)); Z2 <- .zono_upper(.channel(s2, ct))
  Z12 <- .zono_upper(.channel(paste(s1, s2), ct))
  i1 <- .mi_upper(Z1, p); i2 <- .mi_upper(Z2, p); i12 <- .mi_upper(Z12, p)
  iu <- .mi_upper(.zono_join(Z1, Z2), p)
  iu <- min(max(iu, i1, i2), i12)
  c(redundant = max(0, i1 + i2 - iu), unique_1 = max(0, iu - i2), unique_2 = max(0, iu - i1),
    synergy = max(0, i12 - iu), total = i12)
}

#' Finite-resolution MAR-PID of a target given two sources (Telcs, 2026)
#'
#' Nonnegative partial information decomposition (redundant, unique,
#' synergistic atoms) of the information two sources carry about a target,
#' using the Blackwell-order construction of Mages, Anastasiadi and Rohner
#' (2024) on a finite-resolution binary representation: the target is
#' recursively median-binarised to `k_target` bits, each source to
#' `k_source` bits (treated as one finite-valued source), MAR-PID is computed
#' for every target bit and the atoms are summed (and averaged) over the
#' resolved target bits (Sects. 3-5 of the paper). Two sources only; the
#' general antichain lattice is not implemented.
#'
#' @param target,s1,s2 numeric vectors (continuous or discrete) of equal
#'   length; the rows are the joint observations.
#' @param k_target,k_source binarisation depths (1 bit = median split).
#' @param n_surr number of surrogates (rows of `s1` and `s2` permuted
#'   independently) for the p-values of the unique and synergistic atoms and
#'   of the "novel" information of each source (unique + synergistic, i.e.
#'   what the source adds beyond the other one); 0 = no test.
#' @param seed optional RNG seed.
#' @param cores CPU cores for the surrogates (see [Rcausality-parallel]).
#' @return An object of class `marpid`: `atoms` (named vector in bits,
#'   summed over the target bits: `redundant`, `unique_1`, `unique_2`,
#'   `synergy`, `total`), `per_bit` (average per target bit), `by_bit`
#'   (matrix: one row per target bit), `p_values` (surrogate p-values of
#'   `unique_1`, `unique_2`, `synergy`, `novel_1`, `novel_2`), `null` (the
#'   surrogate draws of the same quantities, per target bit), `settings`.
#'   Method: `print`.
#' @references
#' Telcs, A. (2026). Empirical nonnegative finite-resolution MAR-PID for continuous variables. *Entropy*, 28(6), 641. \doi{10.3390/e28060641}
#'
#' Mages, T., Anastasiadi, E. and Rohner, C. (2024). Non-negative decomposition of multivariate information: from minimum to Blackwell-specific information. *Entropy*, 26(5), 424. \doi{10.3390/e26050424}
#'
#' Blackwell, D. (1953). Equivalent comparisons of experiments. *The Annals of Mathematical Statistics*, 24(2), 265-272. \doi{10.1214/aoms/1177729032}
#' @export
#' @examples
#' # Two textbook cases.  XOR: neither source alone tells anything about the
#' # target, together they determine it -- all information is synergistic.
#' # Copies: both sources carry the same information -- all of it is redundant.
#' set.seed(1); a <- rbinom(2000, 1, 0.5); b <- rbinom(2000, 1, 0.5)
#' marpid(target = (a + b) %% 2, s1 = a, s2 = b, k_target = 1, k_source = 1)   # pure synergy
#' marpid(target = a, s1 = a, s2 = a, k_target = 1, k_source = 1)             # pure redundancy
#' \donttest{
#' # continuous variables are binarised to k bits; surrogates give p-values of
#' # the unique and synergistic atoms (here s1 is informative, s2 is noise)
#' t <- rnorm(1000); s1 <- t + rnorm(1000); s2 <- rnorm(1000)
#' marpid(t, s1, s2, k_target = 2, k_source = 2, n_surr = 50, seed = 1)
#' }
marpid <- function(target, s1, s2, k_target = 2L, k_source = 2L, n_surr = 0L, seed = NULL,
                   cores = getOption("Rcausality.cores", 1L)) {
  .local_seed(seed)
  target <- as.numeric(target); s1 <- as.numeric(s1); s2 <- as.numeric(s2); n <- length(target)
  if (length(s1) != n || length(s2) != n) stop("target, s1 and s2 must have the same length.", call. = FALSE)
  nb <- function(v, k) { u <- length(unique(v)); if (u <= 2L) 1L else min(k, floor(log2(u))) }
  kt <- nb(target, k_target); ks1 <- nb(s1, k_source); ks2 <- nb(s2, k_source)
  Ct <- .rq_bits(target, kt); c1 <- .bits_code(.rq_bits(s1, ks1)); c2 <- .bits_code(.rq_bits(s2, ks2))
  decompose <- function(a, b) {
    M <- t(vapply(seq_len(kt), function(l) .marpid_bit(Ct[, l], a, b), numeric(5)))
    rownames(M) <- paste0("bit", seq_len(kt)); M
  }
  by_bit <- decompose(c1, c2)
  atoms <- colSums(by_bit)
  pv <- NULL
  if (n_surr > 0) {
    S <- do.call(rbind, .plapply(seq_len(n_surr), function(b) colSums(decompose(c1[sample.int(n)], c2[sample.int(n)])),
                                 cores = cores))
    pv <- c(unique_1 = (1 + sum(S[, "unique_1"] >= atoms[["unique_1"]])) / (n_surr + 1),
            unique_2 = (1 + sum(S[, "unique_2"] >= atoms[["unique_2"]])) / (n_surr + 1),
            synergy = (1 + sum(S[, "synergy"] >= atoms[["synergy"]])) / (n_surr + 1),
            novel_1 = (1 + sum(S[, "unique_1"] + S[, "synergy"] >= atoms[["unique_1"]] + atoms[["synergy"]])) / (n_surr + 1),
            novel_2 = (1 + sum(S[, "unique_2"] + S[, "synergy"] >= atoms[["unique_2"]] + atoms[["synergy"]])) / (n_surr + 1))
  }
  null <- if (n_surr > 0) list(unique_1 = S[, "unique_1"] / kt, unique_2 = S[, "unique_2"] / kt, synergy = S[, "synergy"] / kt,
                               novel_1 = (S[, "unique_1"] + S[, "synergy"]) / kt, novel_2 = (S[, "unique_2"] + S[, "synergy"]) / kt) else NULL
  structure(list(atoms = atoms, per_bit = atoms / kt, by_bit = by_bit, p_values = pv, null = null, n = n,
                 settings = list(k_target = kt, k_source = c(ks1, ks2), n_surr = n_surr)), class = "marpid")
}

#' @export
print.marpid <- function(x, digits = 3, ...) {
  s <- x$settings
  cat("Finite-resolution MAR-PID, two sources (Telcs, 2026; Mages et al., 2024)\n")
  cat(sprintf("  n = %d, target depth %d bit(s), source depths %d/%d bits\n", x$n, s$k_target, s$k_source[1], s$k_source[2]))
  tb <- data.frame(atom = c("redundant {{1},{2}}", "unique {{1}}", "unique {{2}}", "synergistic {{1,2}}", "total I(T; S1, S2)"),
                   bits = round(x$atoms, digits), per.target.bit = round(x$per_bit, digits))
  if (!is.null(x$p_values)) tb$p.surr <- c("", .fmt_p(x$p_values[["unique_1"]]), .fmt_p(x$p_values[["unique_2"]]),
                                           .fmt_p(x$p_values[["synergy"]]), "")
  print(tb, row.names = FALSE, right = FALSE)
  invisible(x)
}

# rcausality driver: sources = the past of the cause (at the candidate lag)
# and the past of the effect, target = the effect's future.  The directed
# statistic is the "novel" information of the cause's past: unique +
# synergistic atoms, i.e. what the cause adds beyond the effect's own past
# (the Blackwell-specific unique atom alone is zero whenever the cause's
# channel is a garbling of the effect's own channel, even if the cause
# carries independent information).
# The effect's own past enters MAR-PID as one (scalar) source: its last value
# when history = 1, otherwise the linear predictor of y[t] from y[t-1..t-h]
# (the own-past information collapsed to one variable; MAR-PID takes two
# scalar sources).
.own_past <- function(y, t_idx, h) {
  if (h <= 1L) return(y[t_idx - 1L])
  H <- do.call(cbind, lapply(seq_len(h), function(l) y[t_idx - l]))
  as.numeric(cbind(1, H) %*% stats::lm.fit(cbind(1, H), y[t_idx])$coefficients)
}

.rcausality_pid <- function(all, grid, causes, effects, alpha, adjust, lags, history, k_target, k_source, n_perm,
                            cl, cores) {
  n <- nrow(all)
  sc <- .split_cores(nrow(grid), cores); pair_cores <- sc$outer; inner_cores <- sc$inner
  res <- .plapply(seq_len(nrow(grid)), function(i) {
    cs <- grid$cause[i]; ef <- grid$effect[i]; x <- all[[cs]]; y <- all[[ef]]; h <- history[i]
    L <- max(max(lags), h); t_idx <- (L + 1):n; own <- .own_past(y, t_idx, h)
    tests <- lapply(lags, function(tau) marpid(y[t_idx], x[t_idx - tau], own, k_target, k_source, n_surr = n_perm,
                                               cores = inner_cores))
    bl <- do.call(rbind, lapply(seq_along(lags), function(j) { r <- tests[[j]]
      data.frame(cause = cs, effect = ef, lag = lags[j], statistic = r$per_bit[["unique_1"]] + r$per_bit[["synergy"]],
                 p = r$p_values[["novel_1"]], unique_cause = r$per_bit[["unique_1"]], redundant = r$per_bit[["redundant"]],
                 unique_effect = r$per_bit[["unique_2"]], synergy = r$per_bit[["synergy"]], total = r$per_bit[["total"]]) }))
    null_max <- do.call(pmax, lapply(tests, function(r) r$null[["novel_1"]]))   # max-statistic test over the lags
    b <- which.max(bl$statistic)
    list(null = list(statistic = bl$statistic[b], null = null_max, alternative = "greater"),
         row = data.frame(cause = cs, effect = ef, lag = bl$lag[b], statistic = bl$statistic[b],
                          df1 = NA_real_, df2 = NA_real_,
                          p_value = if (length(null_max)) (1 + sum(null_max >= bl$statistic[b])) / (length(null_max) + 1) else NA_real_,
                          bf10 = NA_real_, effect_size = bl$statistic[b] / max(bl$total[b], 1e-12), r2_gain = NA_real_,
                          sign = NA_real_, n_eff = n - L, stringsAsFactors = FALSE),
         by_lag = bl[, c("cause", "effect", "lag", "statistic", "p")],
         atoms = data.frame(cause = cs, effect = ef, lag = bl$lag[b], unique_cause = bl$unique_cause[b],
                            unique_effect = bl$unique_effect[b], redundant = bl$redundant[b],
                            synergy = bl$synergy[b], total = bl$total[b], stringsAsFactors = FALSE))
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
                      method = "pid", conditional = FALSE, alpha = alpha, bf_threshold = NA_real_,
                      adjust = adjust, n = n, call = cl, data = all, by_lag = by_lag,
                      pid = do.call(rbind, lapply(res, `[[`, "atoms")))
  out$null <- stats::setNames(lapply(res, `[[`, "null"), paste(pairs$cause, "->", pairs$effect))
  out
}
