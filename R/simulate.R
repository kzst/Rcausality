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
#' d <- sim_logistic(500, scenario = "x_to_y", coupling = 0.3, seed = 1)
#' attr(d, "truth")
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

#' @rdname simulators
#' @param n_series number of series for `sim_var` (>2 gives a random sparse
#'   network ground truth, returned in attribute `"adjacency"`).
#' @export
sim_var <- function(n = 500, scenario = "x_to_y", coupling = 0.4, noise = 0,
                    n_series = 2L, seed = NULL, burn = 200L) {
  if (!is.null(seed)) set.seed(seed)
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
    A_true <- (B - diag(diag(B), k)) != 0
    dimnames(A_true) <- list(cause = colnames(X), effect = colnames(X))
    attr(df, "adjacency") <- t(A_true)
    return(df)
  }
  stopifnot(scenario %in% .scenarios)
  A <- .coupling_mat(scenario, coupling)
  T <- n + burn; X <- matrix(0, T, 3)
  E <- matrix(stats::rnorm(T * 3), ncol = 3)
  for (t in 2:T) X[t, ] <- 0.5 * X[t - 1, ] + as.numeric(A %*% X[t - 1, ]) + E[t, ]
  X <- X[-seq_len(burn), 1:2, drop = FALSE]; colnames(X) <- c("x", "y")
  .finish(X, scenario, noise)
}

.coupled_map <- function(f, n, scenario, coupling, noise, seed, burn) {
  if (!is.null(seed)) set.seed(seed)
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
  if (!is.null(seed)) set.seed(seed)
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
  if (!is.null(seed)) set.seed(seed)
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
