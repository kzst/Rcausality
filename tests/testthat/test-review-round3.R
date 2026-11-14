# Version 0.0.9: a HAC covariance for the regression tests; the kernel test with conditional shifts as its default
# reference; conditional shifts for the transfer entropy.

quad_pair <- function(n, seed) {                                  # x drives y through its square: no linear effect
  set.seed(seed); x <- as.numeric(arima.sim(list(ar = 0.5), n)); y <- numeric(n)
  for (t in 2:n) y[t] <- 0.5 * y[t - 1] + 0.5 * (x[t - 1]^2 - 1) + rnorm(1)
  data.frame(x = x, y = y)
}
rev_pair <- function(n, seed, f = identity) {                     # y is driven by a persistent x: the test of y -> x is a true null
  set.seed(seed); x <- as.numeric(arima.sim(list(ar = 0.9), n)) * sqrt(0.19); y <- numeric(n)
  for (t in 2:n) y[t] <- 0.5 * y[t - 1] + 1.5 * f(x[t - 1]) + rnorm(1)
  data.frame(x = x, y = y)
}

test_that("the regression tests have a Wald form with a HAC covariance", {
  d <- quad_pair(300, 1)
  h <- suppressWarnings(rcausality(d, method = "granger", lag = 2, method_args = list(vcov = "HC3")))
  a <- expect_no_warning(rcausality(d, method = "granger", lag = 2, method_args = list(vcov = "HAC")))
  expect_true(all(a$extra$vcov == "hac"))
  expect_false(isTRUE(all.equal(a$pairs$statistic, h$pairs$statistic)))
  expect_error(rcausality(d, method = "granger", method_args = list(vcov = "NW")), "HAC")
  ty <- suppressWarnings(rcausality(d, method = "toda_yamamoto", lag = 1, method_args = list(vcov = "HAC")))
  expect_true(all(ty$extra$vcov == "hac"))
  # the statistic: Bartlett weights up to lag floor(4 (N / 100)^(2 / 9)) on the scores with leverage-adjusted residuals
  n <- nrow(d); y <- d$y[3:n]; X <- cbind(1, d$y[2:(n - 1)], d$y[1:(n - 2)], d$x[2:(n - 1)], d$x[1:(n - 2)]); N <- length(y)
  XtXi <- solve(crossprod(X)); b <- XtXi %*% crossprod(X, y); e <- as.numeric(y - X %*% b); hat <- rowSums((X %*% XtXi) * X)
  S <- X * (e / (1 - hat)); L <- floor(4 * (N / 100)^(2 / 9)); M <- matrix(0, 5, 5)
  for (k in -L:L) for (t in seq_len(N)) if (t - k >= 1 && t - k <= N) M <- M + (1 - abs(k) / (L + 1)) * tcrossprod(S[t, ], S[t - k, ])
  V <- XtXi %*% M %*% XtXi; W <- as.numeric(t(b[4:5]) %*% solve(V[4:5, 4:5]) %*% b[4:5]) / 2
  xy <- a$pairs$cause == "x"
  expect_equal(a$pairs$statistic[xy], W, tolerance = 1e-8)
  expect_equal(a$pairs$p_value[xy], stats::pf(W, 2, N - 5, lower.tail = FALSE), tolerance = 1e-8)
  expect_equal(Rcausality:::.hac_lag(300), 5L); expect_equal(Rcausality:::.hac_lag(3), 1L)
  # ... which is close to the Newey-West covariance of the sandwich package (plain residuals there)
  skip_if_not_installed("sandwich"); skip_if_not_installed("lmtest")
  e2 <- data.frame(y = y, y1 = X[, 2], y2 = X[, 3], x1 = X[, 4], x2 = X[, 5])
  fu <- lm(y ~ y1 + y2 + x1 + x2, e2); fr <- lm(y ~ y1 + y2, e2)
  w <- lmtest::waldtest(fr, fu, vcov = sandwich::NeweyWest(fu, lag = L, prewhite = FALSE, adjust = FALSE), test = "F")
  expect_equal(a$pairs$statistic[xy], w$F[2], tolerance = 0.1)
})

test_that("the kernel test with conditional shifts keeps the reverse direction of a link with a persistent driver", {
  # (the normal reference rejects this true null hypothesis in most series of this kind: the bandwidth does not
  # resolve the past of the effect)
  rej <- vapply(1:6, function(s) {
    d <- rev_pair(300, s)
    c(normal = diks_panchenko(d$y, d$x, reference = "normal")$tests$p_value[1] < 0.05,
      conditional = diks_panchenko(d$y, d$x, n_perm = 99, seed = s)$tests$p_value[1] < 0.05)
  }, c(normal = TRUE, conditional = TRUE))
  expect_gte(sum(rej["normal", ]), 4L); expect_lte(sum(rej["conditional", ]), 1L)
  d <- rev_pair(300, 1)
  r <- suppressWarnings(suppressMessages(rcausality(d, method = "diks_panchenko", n_perm = 99, seed = 1)))
  expect_true(all(r$extra$reference == "conditional shifts"))
  expect_equal(r$classes$class, "x_to_y")
  expect_output(print(diks_panchenko(d$x, d$y, n_perm = 40, seed = 1)), "conditional shifts")
  expect_output(print(diks_panchenko(d$x, d$y, reference = "normal")), "normal")
})

test_that("the transfer entropy can be referred to conditional shifts", {
  d <- sim_var(400, "x_to_y", seed = 1)
  r0 <- suppressMessages(rcausality(d, method = "te", n_perm = 60, seed = 1))
  expect_true(all(r0$extra$reference == "circular shifts"))
  r <- suppressMessages(rcausality(d, method = "te", n_perm = 60, seed = 1, method_args = list(reference = "conditional")))
  expect_true(all(r$extra$reference == "conditional shifts")); expect_true(all(is.finite(r$extra$cond_r2)))
  expect_identical(r$settings$method_args$reference, "conditional")
  expect_lt(r$pairs$p_value[r$pairs$cause == "x"], 0.05); expect_gt(r$pairs$p_value[r$pairs$cause == "y"], 0.05)
  expect_error(suppressMessages(rcausality(d, method = "te", method_args = list(reference = "shift"))), "conditional")
  # symbolic series are shifted circularly whatever the setting
  m <- sim_markov_pair(600, "x_to_y", seed = 1)[, c("x", "y")]
  rm <- suppressMessages(rcausality(m, method = "te", n_perm = 60, seed = 1, method_args = list(reference = "conditional")))
  expect_true(all(rm$extra$reference == "circular shifts"))
  a <- transfer_entropy(d$x, d$y, n_surr = 60, seed = 1, reference = "conditional")
  expect_lt(a$p_value, 0.05); expect_identical(a$settings$reference, "conditional shifts"); expect_length(a$null, 60L)
  b <- transfer_entropy(d$x, d$y, k = 2, l = 2, n_surr = 40, seed = 1, reference = "conditional")
  expect_true(is.finite(b$te)); expect_length(b$null, 40L)
  expect_identical(transfer_entropy(m$x, m$y, n_surr = 40, seed = 1, reference = "conditional")$settings$reference, "circular shifts")
})

test_that("a significant test that conditions coarsely on a persistent past gets a warning", {
  d <- rev_pair(1000, 2, function(v) v^2 - 1)                      # x drives y; the test of y -> x is a true null
  w <- character(0)
  r <- withCallingHandlers(suppressMessages(rcausality(d, method = "diks_panchenko", method_args = list(reference = "normal"))),
                           warning = function(e) { w <<- c(w, conditionMessage(e)); invokeRestart("muffleWarning") })
  expect_true(r$pairs$p_value[r$pairs$cause == "y"] < 0.05)       # the normal reference rejects the reverse direction
  expect_true(any(grepl("largely a function of the past of the effect", w)))
  w <- character(0)
  r2 <- withCallingHandlers(suppressMessages(rcausality(d, method = "diks_panchenko", n_perm = 99, seed = 1)),
                            warning = function(e) { w <<- c(w, conditionMessage(e)); invokeRestart("muffleWarning") })
  expect_false(any(grepl("largely a function", w)))
  expect_gt(r2$pairs$p_value[r2$pairs$cause == "y"], 0.05)
  # independent series: no such warning
  set.seed(5); e <- data.frame(x = as.numeric(arima.sim(list(ar = 0.9), 400)), y = as.numeric(arima.sim(list(ar = 0.9), 400)))
  w <- character(0)
  withCallingHandlers(suppressMessages(rcausality(e, method = "te", seed = 1)),
                      warning = function(z) { w <<- c(w, conditionMessage(z)); invokeRestart("muffleWarning") })
  expect_false(any(grepl("largely a function", w)))
})

test_that("the decomposition can be referred to conditional shifts (version 0.0.10)", {
  d <- sim_var(400, "x_to_y", seed = 1)
  r0 <- suppressMessages(rcausality(d, method = "pid", n_perm = 60, seed = 1))
  expect_true(all(r0$extra$reference == "permutation"))
  r <- suppressMessages(rcausality(d, method = "pid", n_perm = 60, seed = 1, method_args = list(reference = "conditional")))
  expect_true(all(r$extra$reference == "conditional shifts")); expect_true(all(is.finite(r$extra$cond_r2)))
  expect_identical(r$settings$method_args$reference, "conditional")
  expect_lt(r$pairs$p_value[r$pairs$cause == "x"], 0.05); expect_gt(r$pairs$p_value[r$pairs$cause == "y"], 0.05)
  expect_gte(length(r$null[[1]]$null), 40L)                       # the scan leaves enough spaced shifts
  expect_error(suppressMessages(rcausality(d, method = "pid", method_args = list(reference = "shift"))), "conditional")
  # the statistic is that of the permutation test at the same lag
  f <- suppressMessages(rcausality(d, method = "pid", lag = 1, n_perm = 40, seed = 1, method_args = list(reference = "conditional")))
  f0 <- suppressMessages(rcausality(d, method = "pid", lag = 1, n_perm = 40, seed = 1))
  expect_equal(f$pairs$statistic, f0$pairs$statistic)
  # binary series are permuted whatever the setting
  m <- sim_markov_pair(600, "x_to_y", seed = 1)[, c("x", "y")]
  rm <- suppressMessages(rcausality(m, method = "pid", n_perm = 40, seed = 1, method_args = list(reference = "conditional")))
  expect_true(all(rm$extra$reference == "permutation"))
})

test_that("the panel versions of cmi and of the kernel test warn for strongly persistent units (version 0.0.11)", {
  set.seed(3)
  mk <- function(a) lapply(1:10, function(i) data.frame(x = as.numeric(arima.sim(list(ar = a), 60)), y = as.numeric(arima.sim(list(ar = a), 60))))
  grid <- data.frame(cause = c("x", "y"), effect = c("y", "x"), stringsAsFactors = FALSE)
  chk <- Rcausality:::.panel_persistence_check
  for (m in c("cmi", "cmi_lag", "diks_panchenko")) expect_warning(chk(mk(0.95), grid, m), "group test adds this excess")
  expect_warning(chk(mk(0.95), grid, "cmi"), "decision = .pooled")
  expect_no_warning(chk(mk(0.5), grid, "cmi"))
  expect_no_warning(chk(mk(0.95), grid, "te"))                      # (the transfer entropy kept its level in such panels)
  expect_no_warning(chk(mk(0.95), grid, "cmi", bayes = TRUE))
  # the limit of the Granger family is unchanged (0.9), and its message is the former one
  expect_warning(chk(mk(0.99), grid, "granger"), "persistent or trending")
})
