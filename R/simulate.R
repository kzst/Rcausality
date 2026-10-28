#' Simulators with known causal ground truth
#'
#' All simulators return a data frame with columns `x`, `y` (and `z` for the
#' common-cause scenario) and attach the true class as attribute `"truth"`.
#'
#' Scenarios: `"x_to_y"`, `"y_to_x"`, `"bidirectional"`, `"common_cause"`
#' (a hidden driver `z` drives both `x` and `y`), `"independent"`.
#'
#' @param n number of observations returned.
#' @param scenario causal scenario (see Details).
#' @param coupling coupling strength (maps: coupled as
#'   `(1-c) f(y) + c f(x)`; strong coupling synchronises the systems).
#' @param noise standard deviation of additive observational noise
#'   (relative to the signal sd).
#' @param seed optional RNG seed.
#' @param burn burn-in length.
#' @return A data frame with attribute `"truth"`.
#' @name simulators
#' @examples
#' # two coupled logistic maps, x driving y with coupling strength 0.3;
#' # the true relation is kept as an attribute for benchmarks
#' d <- sim_logistic(500, scenario = "x_to_y", coupling = 0.3, seed = 1)
#' attr(d, "truth")
#' plot(d$x[-500], d$x[-1], pch = 20, cex = 0.4)   # the parabola of the logistic map
NULL

.scenarios <- c("x_to_y", "y_to_x", "bidirectional", "common_cause", "independent")

.finish <- function(M, scenario, noise) {
  if (noise > 0) M <- M + matrix(stats::rnorm(length(M), 0, noise * stats::sd(M)), nrow(M))
  df <- as.data.frame(M)
  attr(df, "truth") <- scenario
  df
}

.coupling_mat <- function(scenario, c) {
  # rows = target, cols = source for variables (x, y, z)
  A <- matrix(0, 3, 3)
  switch(scenario,
         x_to_y = A[2, 1] <- c, y_to_x = A[1, 2] <- c,
         bidirectional = { A[2, 1] <- c; A[1, 2] <- c },
         common_cause = { A[1, 3] <- c; A[2, 3] <- c },
         independent = NULL)
  A
}

#' @references
#' Luetkepohl, H. (2005). *New Introduction to Multiple Time Series Analysis*. Springer, Berlin.
#'
#' Zlatniczki et al. (2021) and Stippinger et al. (2023) (benchmark designs); see Zlatniczki, A., Stippinger, M., Benko, Z., Somogyvari, Z. and Telcs, A. (2021). Relaxation of some confusions about confounders. *Entropy*, 23(11), 1450. \doi{10.3390/e23111450}; Stippinger, M., Bencze, A., Zlatniczki, A., Somogyvari, Z. and Telcs, A. (2023). Causal discovery of stochastic dynamical systems: a Markov chain approach. *Mathematics*, 11(4), 852. \doi{10.3390/math11040852}
#' @rdname simulators
#' @param n_series number of series for `sim_var` (>2 gives a random sparse
#'   network ground truth, returned in attribute `"adjacency"`).
#' @export
sim_var <- function(n = 500, scenario = "x_to_y", coupling = 0.4, noise = 0,
                    n_series = 2L, seed = NULL, burn = 200L) {
  .local_seed(seed)
  if (n_series > 2L) {
    k <- n_series
    A <- matrix(0, k, k)
    cand <- which(row(A) != col(A))
    pick <- sample(cand, max(1L, round(0.25 * length(cand))))
    A[pick] <- sample(c(-1, 1), length(pick), TRUE) * coupling
    B <- diag(0.5, k) + A
    ev <- max(abs(eigen(B, only.values = TRUE)$values))
    if (ev > 0.9) B <- B * 0.9 / ev
    E <- matrix(stats::rnorm((n + burn) * k), ncol = k); X <- matrix(0, n + burn, k)
    for (t in 2:(n + burn)) X[t, ] <- B %*% X[t - 1, ] + E[t, ]
    X <- X[-seq_len(burn), , drop = FALSE]
    colnames(X) <- paste0("v", seq_len(k))
    df <- .finish(X, "network", noise)
    A_true <- t((B - diag(diag(B), k)) != 0)      # B[j, i] != 0  <=>  v_i -> v_j
    dimnames(A_true) <- list(cause = colnames(X), effect = colnames(X))
    attr(df, "adjacency") <- A_true
    return(df)
  }
  stopifnot(scenario %in% .scenarios)
  A <- .coupling_mat(scenario, coupling)
  B <- diag(0.5, 3) + A
  ev <- max(abs(eigen(B, only.values = TRUE)$values))
  if (ev >= 0.95) {                      # keep the process stationary (e.g. strong bidirectional coupling)
    warning("coupling too strong for a stationary VAR; coefficients rescaled by ", signif(0.9 / ev, 3), call. = FALSE)
    B <- B * 0.9 / ev
  }
  T <- n + burn; X <- matrix(0, T, 3)
  E <- matrix(stats::rnorm(T * 3), ncol = 3)
  for (t in 2:T) X[t, ] <- as.numeric(B %*% X[t - 1, ]) + E[t, ]
  X <- X[-seq_len(burn), 1:2, drop = FALSE]; colnames(X) <- c("x", "y")
  .finish(X, scenario, noise)
}

.coupled_map <- function(f, n, scenario, coupling, noise, seed, burn) {
  .local_seed(seed)
  stopifnot(scenario %in% .scenarios)
  A <- .coupling_mat(scenario, coupling)
  T <- n + burn; X <- matrix(NA_real_, T, 3)
  X[1, ] <- stats::runif(3, 0.2, 0.8)
  for (t in 2:T) {
    prev <- X[t - 1, ]
    for (i in 1:3) {
      src <- which(A[i, ] > 0)
      c_i <- if (length(src)) A[i, src[1]] else 0
      drive <- if (length(src)) prev[src[1]] else prev[i]
      X[t, i] <- (1 - c_i) * f(prev[i]) + c_i * f(drive)
      if (i == 3) X[t, i] <- f(prev[3])
    }
    X[t, ] <- pmin(pmax(X[t, ], 1e-6), 1 - 1e-6)
  }
  X <- X[-seq_len(burn), , drop = FALSE]
  X <- X[, 1:2]; colnames(X) <- c("x", "y")
  if (max(abs(X[, 1] - X[, 2])) < 1e-6)
    warning("The two systems synchronised; lower 'coupling'.", call. = FALSE)
  .finish(X, scenario, noise)
}

#' @rdname simulators
#' @param r growth parameter of the logistic / tent map.
#' @export
sim_logistic <- function(n = 500, scenario = "x_to_y", coupling = 0.2, noise = 0,
                         r = 3.9, seed = NULL, burn = 200L) {
  .coupled_map(function(u) r * u * (1 - u), n, scenario, coupling, noise, seed, burn)
}

#' @rdname simulators
#' @export
sim_tent <- function(n = 500, scenario = "x_to_y", coupling = 0.2, noise = 0,
                     seed = NULL, burn = 200L) {
  .coupled_map(function(u) 1.99 * (0.5 - abs(u - 0.5)) * 1, n, scenario, coupling, noise, seed, burn)
}

#' @rdname simulators
#' @export
sim_henon <- function(n = 500, scenario = "x_to_y", coupling = 0.3, noise = 0,
                      seed = NULL, burn = 500L) {
  .local_seed(seed)
  stopifnot(scenario %in% .scenarios)
  A <- .coupling_mat(scenario, coupling)
  T <- n + burn; X <- matrix(0, T, 3); X[1:2, ] <- stats::runif(6, -0.1, 0.1)
  a <- 1.4; b <- 0.3
  for (t in 3:T) {
    for (i in 1:3) {
      src <- which(A[i, ] > 0)
      c_i <- if (length(src)) A[i, src[1]] else 0
      s <- if (length(src)) (1 - c_i) * X[t - 1, i] + c_i * X[t - 1, src[1]] else X[t - 1, i]
      X[t, i] <- a - s^2 + b * X[t - 2, i]
    }
  }
  X <- X[-seq_len(burn), 1:2, drop = FALSE]; colnames(X) <- c("x", "y")
  if (any(!is.finite(X))) stop("Henon trajectory diverged; reduce coupling.", call. = FALSE)
  .finish(X, scenario, noise)
}

#' @rdname simulators
#' @param p transition parameter of the binary Markov chains.
#' @export
sim_markov01 <- function(n = 1000, scenario = "x_to_y", coupling = 0.5, noise = 0,
                         p = 0.8, seed = NULL, burn = 100L) {
  .local_seed(seed)
  stopifnot(scenario %in% .scenarios)
  A <- .coupling_mat(scenario, coupling)
  T <- n + burn; X <- matrix(0L, T, 3); X[1, ] <- sample(0:1, 3, TRUE)
  for (t in 2:T) for (i in 1:3) {
    src <- which(A[i, ] > 0)
    own <- X[t - 1, i]
    stay <- stats::runif(1) < p
    base <- if (stay) own else 1L - own
    X[t, i] <- if (length(src) && stats::runif(1) < A[i, src[1]]) X[t - 1, src[1]] else base
  }
  X <- X[-seq_len(burn), 1:2, drop = FALSE]; colnames(X) <- c("x", "y")
  .finish(X, scenario, 0)
}

#' Simulators of the Markov-chain paper (Stippinger et al., 2023, Sect. 6.1)
#'
#' `sim_turnwheel()` is the "turn-wheel" example: Alice is a two-state Markov
#' chain (stay probabilities `p`, `q`); Bob and Barry are independent
#' three-state chains that, with one step delay, turn their wheel by +1 (right)
#' with probability `pB` if Alice is in state 0 and by -1 otherwise (and the
#' opposite way with probability `1 - pB`).
#' `sim_logistic3()` is the system of three logistic maps with system noise
#' and reflective boundary: an autonomous driver `x0`, and `x1`, `x2` that are
#' either both driven by `x0` (`scenario = "driver"`, `b > 0`, `c = 0`) or drive
#' each other (`scenario = "mutual"`, `b = 0`, `c > 0`), always with one step
#' delay.
#'
#' @param n series length.
#' @param p,q Alice's probabilities of repeating state 0 / state 1.
#' @param pB probability of following Alice's advice.
#' @param scenario `"driver"` or `"mutual"` (`sim_logistic3`).
#' @param b,c coupling strengths (driver -> driven, mutual).
#' @param sigma standard deviation of the system (dynamical) noise.
#' @param obs_noise standard deviation of additive observational noise added
#'   to the returned series (the paper's Figure 3 has observational noise
#'   only; its level is not stated).
#' @param burn warm-up steps (default 100, as in the paper).
#' @param seed optional RNG seed.
#' @return A data frame (`alice`, `bob`, `barry` or `x0`, `x1`, `x2`) with
#'   `attr(, "truth")`.
#' @references
#' Stippinger, M., Bencze, A., Zlatniczki, A., Somogyvari, Z. and Telcs, A.
#' (2023). Causal discovery of stochastic dynamical systems: a Markov chain
#' approach. *Mathematics*, 11(4), 852. \doi{10.3390/math11040852}
#' (Section 6.1, Figures 2-5).
#' @export
#' @examples
#' # three-state wheels turned by the previous state of the other wheel
#' w <- sim_turnwheel(500, seed = 1)
#' head(w)
#' table(w$alice, w$bob)      # bob turns right (+1) after alice = 0
sim_turnwheel <- function(n = 1000, p = 0.75, q = 0.75, pB = 0.7, seed = NULL) {
  .local_seed(seed)
  a <- integer(n); a[1] <- sample(0:1, 1)
  for (t in 2:n) a[t] <- if (a[t - 1] == 0) (if (stats::runif(1) < p) 0L else 1L) else (if (stats::runif(1) < q) 1L else 0L)
  step <- function() {
    s <- integer(n); s[1] <- sample(0:2, 1)
    for (t in 2:n) {
      right <- if (a[t - 1] == 0) stats::runif(1) < pB else stats::runif(1) >= pB
      s[t] <- (s[t - 1] + if (right) 1L else -1L) %% 3L
    }
    s
  }
  d <- data.frame(alice = a, bob = step(), barry = step())
  attr(d, "truth") <- "alice -> bob, alice -> barry"
  d
}

#' @rdname sim_turnwheel
#' @export
sim_logistic3 <- function(n = 1000, scenario = c("driver", "mutual"), b = 1.5, c = 1.5,
                          sigma = 0, obs_noise = 0, burn = 100L, seed = NULL) {
  scenario <- match.arg(scenario)
  .local_seed(seed)
  if (scenario == "driver") c <- 0 else b <- 0
  r <- stats::runif(3, 3.8, 3.9)
  fold <- function(z) abs((z - 1) - 2 * floor((z - 1) / 2) - 1)
  N <- n + burn; X <- matrix(0, N, 3); X[1, ] <- stats::runif(3)
  for (t in 1:(N - 1)) {
    e <- stats::rnorm(3, 0, sigma); x <- X[t, ]
    X[t + 1, 1] <- fold(r[1] * x[1] * (1 - x[1] + e[1]))
    X[t + 1, 2] <- fold(r[2] * x[2] * (1 - x[2] - b * x[1] - c * x[3] + e[2]))
    X[t + 1, 3] <- fold(r[3] * x[3] * (1 - x[3] - b * x[1] - c * x[2] + e[3]))
  }
  d <- as.data.frame(X[(burn + 1):N, ]); names(d) <- c("x0", "x1", "x2")
  if (obs_noise > 0) d[] <- lapply(d, function(v) v + stats::rnorm(n, 0, obs_noise))
  attr(d, "truth") <- if (scenario == "driver") "x0 -> x1, x0 -> x2" else "x1 <-> x2"
  d
}

#' Linear rotation system of the df-causality paper (Telcs, Kurbucz & Jakovac, 2025, Sect. 6.4)
#'
#' Six-dimensional linear stochastic recursion `s[n] = Q s[n-1] + sigma * xi[n]`
#' with three two-dimensional rotating subsystems `x`, `y`, `z` (Eq. 25): `z`
#' is autonomous and drives both `x` and `y` through the coupling `g` (its
#' first component enters both components of `x` and of `y`), so `x` and `y`
#' have a hidden common cause and no direct link. Two parameter sets are
#' used in the paper: `"deterministic"` (Eq. 29: `sigma = 5e-4`,
#' `alpha = 1 - 1e-6`, `g = 0.3`) and `"stochastic"` (Eq. 30: `sigma = 0.01`,
#' `alpha = 1 - 1e-4`, `g = 0.1`), with rotation angles `phi = (0.5, 0.3, 0.4)`.
#'
#' @param n series length.
#' @param set parameter set, `"deterministic"` or `"stochastic"`.
#' @param sigma,g,alpha,phi override the parameters of the chosen set.
#' @param burn warm-up steps discarded.
#' @param seed optional RNG seed.
#' @return A data frame with columns `x1, x2, y1, y2, z1, z2` and
#'   `attr(, "truth")`.
#' @references
#' Telcs, A., Kurbucz, M. T. and Jakovac, A. (2025). Unified causality analysis based on the degrees of freedom. *Physical Review E*, 112(1), 014204. \doi{10.1103/PhysRevE.112.014204}, Section 6.4, Eqs. (25)-(30).
#' @export
#' @examples
#' # the rotation system of the df paper: z1 drives z2
#' d <- sim_df_rotation(2000, "stochastic", seed = 1)
#' plot(d$z1[1:200], type = "l")
sim_df_rotation <- function(n = 5000, set = c("deterministic", "stochastic"), sigma = NULL,
                            g = NULL, alpha = NULL, phi = c(0.5, 0.3, 0.4), burn = 1000L, seed = NULL) {
  set <- match.arg(set)
  .local_seed(seed)
  if (is.null(sigma)) sigma <- if (set == "deterministic") 5e-4 else 0.01
  if (is.null(g)) g <- if (set == "deterministic") 0.3 else 0.1
  if (is.null(alpha)) alpha <- if (set == "deterministic") 1 - 1e-6 else 1 - 1e-4
  rot <- function(a) alpha * matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2)  # columns: (c, s), (-s, c)
  Q <- matrix(0, 6, 6)
  Q[1:2, 1:2] <- rot(phi[1]); Q[3:4, 3:4] <- rot(phi[2]); Q[5:6, 5:6] <- rot(phi[3])
  Q[1:4, 5] <- g
  N <- n + burn; S <- matrix(0, N, 6); S[1, ] <- stats::rnorm(6)
  for (t in 2:N) S[t, ] <- Q %*% S[t - 1, ] + sigma * stats::rnorm(6)
  d <- as.data.frame(S[(burn + 1):N, ]); names(d) <- c("x1", "x2", "y1", "y2", "z1", "z2")
  attr(d, "truth") <- "z -> x, z -> y; x and y share the hidden common cause z"
  d
}

#' Coupled logistic maps of the Bayesian dimensional-causality paper (Benko et al., 2024, Eq. B.1)
#'
#' Three logistic maps `x_j[t+1] = r x_j[t] (1 - x_j[t] - sum_l beta_jl x_l[t])`
#' with r = 3.99; the first two are observed (`x`, `y`), the third is hidden.
#' Scenarios (coupling matrices (a)-(d) of Appendix B.2): `"x_to_y"`
#' (`beta_21`), `"y_to_x"`, `"bidirectional"` (`beta_12 = beta_21`),
#' `"common_cause"` (the hidden third map drives both, `beta_13 = beta_23`)
#' and `"independent"`. Values are kept in `[0, 1]` by the reflective fold of
#' Stippinger et al. (2023); with `beta = 0.5` it is rarely needed.
#'
#' @param n series length.
#' @param scenario coupling scenario.
#' @param beta coupling strength (paper: 0.5).
#' @param r logistic parameter (paper: 3.99).
#' @param burn warm-up steps.
#' @param seed optional RNG seed.
#' @return A data frame with `x`, `y` and `attr(, "truth")`.
#' @references
#' Benko, Z., Zlatniczki, A., Stippinger, M., Fabo, D., Solyom, A., Eross, L., Telcs, A. and Somogyvari, Z. (2024). Bayesian inference of causal relations between dynamical systems. *Chaos, Solitons & Fractals*, 185, 115142. \doi{10.1016/j.chaos.2024.115142}
#' @export
#' @examples
#' # the logistic-map systems of the dimensional-causality paper; here x and y
#' # are both driven by a hidden third map
#' d <- sim_logistic_bdc(500, "common_cause", seed = 1)
#' head(d)
#' attr(d, "truth")
sim_logistic_bdc <- function(n = 10000, scenario = c("x_to_y", "y_to_x", "bidirectional",
                                                     "common_cause", "independent"),
                             beta = 0.5, r = 3.99, burn = 500L, seed = NULL) {
  scenario <- match.arg(scenario)
  .local_seed(seed)
  B <- matrix(0, 3, 3)
  switch(scenario,
         x_to_y = B[2, 1] <- beta, y_to_x = B[1, 2] <- beta,
         bidirectional = { B[1, 2] <- beta; B[2, 1] <- beta },
         common_cause = { B[1, 3] <- beta; B[2, 3] <- beta }, independent = NULL)
  fold <- function(z) abs((z - 1) - 2 * floor((z - 1) / 2) - 1)
  N <- n + burn; X <- matrix(0, N, 3); X[1, ] <- stats::runif(3)
  for (t in 1:(N - 1)) {
    x <- X[t, ]
    X[t + 1, ] <- fold(r * x * (1 - x - as.vector(B %*% x)))
  }
  d <- data.frame(x = X[(burn + 1):N, 1], y = X[(burn + 1):N, 2])
  attr(d, "truth") <- scenario
  d
}
