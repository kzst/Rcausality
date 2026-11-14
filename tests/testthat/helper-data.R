# A stochastic pair with a nonlinear link from x to y that leaves no linear trace (the mean of y depends on the square
# of the past of x).  The tests of conditional independence are tried on such data: on noise-free deterministic maps
# the cause is a function of the past of the effect, and no conditional test is informative (see .determined()).
nl_pair <- function(n, seed) {
  set.seed(seed); x <- as.numeric(stats::arima.sim(list(ar = 0.5), n)); y <- numeric(n)
  for (t in 2:n) y[t] <- 0.3 * y[t - 1] + 0.8 * (x[t - 1]^2 - 1) + stats::rnorm(1)
  data.frame(x = x, y = y)
}
