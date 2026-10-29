# Regression tests for the issues found in the 0.0.1 code review (October 2026).

test_that("CMI decisions do not depend on the units of one series", {
  lg <- sim_logistic(400, "x_to_y", seed = 2)
  r0 <- rcausality(lg$x, lg$y, method = "cmi", lag = 1, n_perm = 60, seed = 1)
  r1 <- rcausality(1000 * lg$x, lg$y, method = "cmi", lag = 1, n_perm = 60, seed = 1)
  r2 <- rcausality(lg$x, 1000 * lg$y, method = "cmi", lag = 1, n_perm = 60, seed = 1)
  expect_equal(r0$pairs$statistic, r1$pairs$statistic, tolerance = 1e-6)
  expect_equal(r0$pairs$statistic, r2$pairs$statistic, tolerance = 1e-6)
  expect_equal(r1$classes$class, "x_to_y"); expect_equal(r2$classes$class, "x_to_y")
  set.seed(1); x <- rnorm(300); y <- x + rnorm(300); z <- rnorm(300)
  a <- local_permutation_test(x, y, z, n_perm = 30)$cmi
  b <- local_permutation_test(1000 * x, y, z, n_perm = 30)$cmi
  expect_equal(a, b, tolerance = 1e-6)
})

test_that("the conditioning history follows the effect's memory (no spurious reverse link)", {
  skip_on_cran()      # long-running; run by devtools::test() / NOT_CRAN=true
  set.seed(2); n <- 600
  y <- as.numeric(stats::arima.sim(list(ar = c(0.5, 0.3)), n)); x <- numeric(n)
  for (t in 2:n) x[t] <- 0.3 * x[t - 1] + 0.6 * y[t - 1] + rnorm(1)
  d <- data.frame(x = x, y = y)
  for (m in c("cmi", "te", "pid")) {
    r <- rcausality(d, method = m, seed = 1)
    expect_equal(r$classes$class, "y_to_x", info = m)
    expect_true(all(r$pairs$n_eff < n - 1))
  }
  expect_equal(rcausality(d, method = "granger")$classes$class, "y_to_x")
})

test_that("the Table-2 reading handles white-noise series and needs nominal directed evidence", {
  expect_equal(Rcausality:::.df_decide(0, 1, 1, x_drives_y = FALSE), "independent")
  expect_equal(Rcausality:::.df_decide(0, 2, 2, x_drives_y = TRUE), "x_to_y")
  expect_equal(Rcausality:::.df_decide(1, 0, 1, y_drives_x = TRUE), "y_to_x")
  expect_equal(Rcausality:::.df_decide(0, 0, 0), "independent")
  expect_equal(Rcausality:::.df_decide(2, 2, 4), "independent")
  expect_equal(Rcausality:::.df_decide(2, 3, 3), "x_to_y")
  set.seed(2); d <- data.frame(x = rnorm(400), y = as.numeric(stats::arima.sim(list(ar = 0.6), 400)))
  r <- rcausality(d, method = "df", lag = 1, n_perm = 39, seed = 1)
  expect_equal(r$classes$class, "independent")
})

test_that("transfer entropy works on binary / discrete series", {
  set.seed(1); x <- rbinom(500, 1, 0.7); y <- c(0L, x[-500])
  r <- transfer_entropy(x, y, n_surr = 40, seed = 1)
  expect_gt(r$te, 0.5); expect_lt(r$p_value, 0.05)
  expect_equal(r$te, transfer_entropy(1 - x, 1 - y, n_surr = 40, seed = 1)$te)
  expect_equal(unname(table(Rcausality:::.sym_bins(c(rep(1, 300), rnorm(200)), 4))[3]), 300L)
  expect_error(transfer_entropy(rnorm(100), rnorm(100), l = 0), "at least 1")
})

test_that("the MAR-PID Blackwell join gives two unique atoms (Z-channel and its mirror)", {
  set.seed(1); n <- 20000; T <- rbinom(n, 1, 0.5)
  s1 <- ifelse(T == 1, rbinom(n, 1, 0.5), 0L); s2 <- ifelse(T == 0, rbinom(n, 1, 0.5), 1L)
  a <- marpid(T, s1, s2, 1, 1)$atoms
  expect_gt(a[["unique_1"]], 0.15); expect_gt(a[["unique_2"]], 0.15)
  expect_lt(abs(a[["unique_1"]] - a[["unique_2"]]), 0.03)
  expect_lt(a[["redundant"]], 0.15)
  x <- rbinom(n, 1, 0.5); y <- rbinom(n, 1, 0.5)
  expect_equal(unname(marpid((x + y) %% 2, x, y, 1, 1)$atoms[["synergy"]]), 1, tolerance = 1e-3)
  expect_equal(unname(marpid(x, x, y, 1, 1)$atoms[["unique_1"]]), 1, tolerance = 1e-3)
})

test_that("the df-entropy method rejects short or incomplete input instead of crashing", {
  skip_if_not_installed("ranger")
  expect_error(Rcausality:::.ce_curves(rbinom(10, 1, .5), rbinom(10, 1, .5)), "too short")
  expect_error(Rcausality:::.ce_curves(c(NA, rbinom(99, 1, .5)), rbinom(100, 1, .5)), "Missing")
  expect_error(cond_entropy_bits_cpp(c(0L, NA), matrix(0L, 2, 1), 2L))
  m <- dfe_train(n_per_class = 3, lengths = 300L, seed = 1)
  expect_equal(m$lengths, 300L)
})

test_that("directions = 'forward' is honoured by every method", {
  skip_if_not_installed("ranger")
  net <- sim_var(300, n_series = 3, seed = 1)
  m_small <- dfe_train(n_per_class = 3, lengths = c(200, 300), seed = 1)
  for (m in c("markov", "hirata", "cd", "dc", "ccdh", "dfe", "te", "pid")) {
    r <- suppressMessages(rcausality(net[, 1:2], net[, 3, drop = FALSE], method = m, directions = "forward",
                                     seed = 1, n_perm = 40, model = m_small))
    expect_equal(nrow(r$pairs), 2L, info = m)
    expect_true(all(r$pairs$effect == "v3"), info = m)
    expect_false(any(r$classes$both_directions_tested), info = m)
  }
})

test_that("max_lag is honoured by the model-free methods (lag-5 link)", {
  skip_on_cran()      # long-running; run by devtools::test() / NOT_CRAN=true
  set.seed(3); n <- 600; x <- rnorm(n); y <- numeric(n)
  for (t in 6:n) y[t] <- 0.4 * y[t - 1] + 0.8 * x[t - 5] + rnorm(1)
  r <- rcausality(data.frame(x, y), method = "cmi", max_lag = 6, seed = 1)
  P <- r$pairs
  expect_equal(P$lag[P$cause == "x" & P$effect == "y"], 5L)
  expect_true(P$significant[P$cause == "x" & P$effect == "y"])
  expect_equal(nrow(r$by_lag), 12L)
})

test_that("forecast evaluation refits the links on the training part and copes with constant series", {
  d <- sim_var(400, "x_to_y", seed = 1); r <- rcausality(d)
  fc <- predict(r, h = 3)
  expect_equal(fc$evaluation_source, "training part (refit)")
  expect_equal(predict(r, h = 3, refit = FALSE)$evaluation_source, "full sample (optimistic)")
  fc0 <- predict(r, h = 3, holdout = 0)
  expect_null(fc0$evaluation); expect_output(print(fc0), "no hold-out")
  dc <- cbind(d, z = 1); rz <- rcausality(dc); fz <- predict(rz, h = 3)
  expect_true(all(is.finite(fz$forecasts$lower)))
  expect_equal(fz$evaluation$improvement[fz$evaluation$series == "z"], 0)
})

test_that("benchmark_methods orients the truth by the simulator's columns and reports failures", {
  simr <- function(n, scenario, seed) { d <- sim_var(n, scenario, seed = seed); names(d) <- c("y1", "x1"); d }
  b <- benchmark_methods("granger", scenarios = c("x_to_y", "y_to_x"), simulator = simr, n = 300, reps = 2, seed = 1)
  expect_true(all(b$accuracy == 1))
  simbad <- function(n, scenario, seed) data.frame(x = rnorm(n), y = NA_real_)   # all-missing series: the analysis fails
  expect_warning(bb <- benchmark_methods("granger", scenarios = "independent", simulator = simbad, n = 50, reps = 1),
                 "failed")
  expect_true(all(!is.na(bb$runs$error)))
})

test_that("ccdh and the pairwise drivers do not depend on the column names", {
  d <- sim_logistic3(800, "driver", b = 0.4, seed = 3)
  r1 <- rcausality(data.frame(p = d$x0, q = d$x1), method = "ccdh", seed = 1, n_surr = 40)
  r2 <- rcausality(data.frame(q = d$x1, p = d$x0), method = "ccdh", seed = 1, n_surr = 40)
  expect_equal(r1$classes$relation, "p -> q"); expect_equal(r2$classes$relation, "p -> q")
  set.seed(1); x <- round(rnorm(500), 1)
  expect_equal(unname(ccdh(x, x, m = 3, n_surr = 10, seed = 1)$delta), c(0, 0))
})

test_that("dimensional causality: odd k and discrete data", {
  db <- sim_logistic_bdc(1200, "x_to_y", seed = 1)
  expect_equal(dimensional_causality(db$x, db$y, m = 3, k_range = c(11, 15, 21), seed = 1)$best, "x_to_y")
  expect_no_error(dimensional_causality(round(db$x, 1), round(db$y, 1), m = 3, k_range = c(10, 20), seed = 1))
})

test_that("sim_var stays stationary and the simulators restore the RNG state", {
  expect_warning(d <- sim_var(300, "bidirectional", coupling = 0.6, seed = 1), "stationary")
  expect_lt(max(abs(as.matrix(d))), 20)
  set.seed(5); a <- runif(1); set.seed(5); invisible(sim_logistic(50, seed = 9)); b <- runif(1)
  expect_equal(a, b)
})

test_that("contemporaneous tables are matched whatever the order of the pair", {
  pairs <- data.frame(cause = c("b", "a"), effect = c("a", "b"), lag = 1L, statistic = 0, df1 = NA, df2 = NA,
                      p_value = 0.5, p_adj = 0.5, bf10 = NA, post_prob = NA, evidence = NA, effect_size = 0,
                      r2_gain = 0, sign = NA, n_eff = 100, significant = FALSE, stringsAsFactors = FALSE)
  inst <- data.frame(var1 = "b", var2 = "a", cmi = 0.3, p_value = 0.001, p_adj = 0.001, significant = TRUE)
  r <- Rcausality:::.new_causres(pairs, series = c("b", "a"), causes = c("b", "a"), effects = c("b", "a"), method = "cmi_lag",
                                 conditional = FALSE, alpha = 0.05, bf_threshold = NA, adjust = "holm", n = 100,
                                 call = NULL, data = NULL, instant = inst)
  expect_equal(r$classes$class, "common_cause")
})

test_that("recurrence counts never include the diagonal", {
  set.seed(1); h <- hirata_rp(rnorm(400), rnorm(400), E = 2, theiler = 0, n_surr = 10)
  expect_lt(max(abs(h$index)), 0.1)
  expect_equal(h$settings$theiler, 1L)
})

test_that("non-numeric columns and function inputs give clear errors", {
  expect_error(rcausality(data.frame(date = Sys.Date() + 1:50, x = rnorm(50), y = rnorm(50))), "Non-numeric")
  expect_error(rcausality(rnorm(50), df), "function")
})

test_that("ordinal patterns use canonical codes", {
  expect_equal(ordinal_pattern(c(1, 2, 3), 3), 1L)
  expect_equal(ordinal_pattern(c(3, 2, 1), 3), 6L)
  expect_equal(ordinal_pattern(c(1, 3, 2, 5, 4, 6), 3), c(2L, 3L, 2L, 3L))
  expect_error(Rcausality:::.ordinal_codes(matrix(rnorm(26), 2)), "at most 12")
})

test_that("ccm uses the largest library for rho_max and plots work for posterior / skill based results", {
  lg <- sim_logistic(300, "x_to_y", seed = 1)
  r <- ccm(lg$x, lg$y, E = 2, lib_sizes = c(100, 20, 50), n_surr = 10, seed = 1)
  expect_equal(r$rho$lib_size, c(20, 50, 100)); expect_equal(r$rho_max, r$rho$rho[3])
  pdf(NULL); on.exit(dev.off())
  rc <- rcausality(lg, method = "ccm", seed = 1, n_perm = 40)
  expect_no_error(plot(rc, "bars")); expect_no_error(plot(rc, "embedding"))
  rd <- rcausality(sim_logistic_bdc(600, "x_to_y", seed = 1), method = "dc")
  expect_no_error(plot(rd, "bars")); expect_no_error(plot(rd, "embedding"))
})

test_that("explicit theta / k are not overridden", {
  d <- sim_logistic_bdc(500, "x_to_y", seed = 1)
  r <- rcausality(d, method = "hirata", theta = 0.45, n_perm = 20, seed = 1)
  expect_equal(r$args$theta, 0.45)
  r2 <- rcausality(d, method = "hirata", n_perm = 20, seed = 1)
  expect_null(r2$args$theta)
})
