test_that("knn_cmi recovers known mutual information", {
  set.seed(1); x <- rnorm(500); y <- x + rnorm(500)
  expect_equal(knn_cmi(x, y), -0.5 * log(1 - 0.5), tolerance = 0.25)
  z <- rnorm(500)
  expect_lt(abs(knn_cmi(x, z, y)), 0.06)
})

test_that("local permutation test separates dependence", {
  set.seed(3); z <- rnorm(250); x <- z + rnorm(250, 0, .5); y <- z + rnorm(250, 0, .5)
  expect_gt(local_permutation_test(x, y, z, n_perm = 49)$p_value, 0.05)
  expect_lt(local_permutation_test(x, y, NULL, n_perm = 49)$p_value, 0.05)
})

test_that("method = 'cmi' detects non-linear coupling direction", {
  d <- sim_logistic(300, scenario = "x_to_y", coupling = 0.25, seed = 1)
  r <- rcausality(d$x, d$y, method = "cmi", seed = 1)
  expect_equal(r$classes$class, "x_to_y")
  expect_output(print(r), "mutual information")
  pdf(NULL); on.exit(dev.off()); expect_no_error(plot(r, "both"))
  expect_no_error(summary(r))
})

test_that("cmi_lag finds a common (shared-noise) driver, not a directed link", {
  set.seed(5); n <- 300
  e1 <- rnorm(n); e2 <- 0.8 * e1 + 0.6 * rnorm(n)
  x <- stats::filter(e1, 0.5, "recursive"); y <- stats::filter(e2, 0.5, "recursive")
  r <- rcausality(as.numeric(x), as.numeric(y), method = "cmi_lag", seed = 1)
  expect_equal(r$classes$class, "common_cause")
  expect_true(r$classes$common_driver)
  expect_output(print(r), "Contemporaneous")
  pdf(NULL); on.exit(dev.off()); expect_no_error(plot(r))
})

test_that("cmi_lag keeps direction for a driven pair and helpers work", {
  d <- sim_logistic(300, scenario = "x_to_y", coupling = 0.25, seed = 1)
  r <- rcausality(d$x, d$y, method = "cmi_lag", seed = 1)
  expect_true(r$classes$class %in% c("x_to_y", "bidirectional"))
  expect_true(r$pairs$significant[r$pairs$cause == "x"])
  expect_equal(nrow(cmi_convergence(d$x, d$y, 1, 1:2)), 2)
  expect_match(lag_triangle_check(2, 5, 1), "direct")
  expect_match(lag_triangle_check(2, 2, 1), "common")
})

test_that("C++ kNN kernels agree with the pure-R reference implementation", {
  set.seed(11); n <- 150
  x <- rnorm(n); y <- 0.7 * x + rnorm(n); z <- cbind(rnorm(n), x + rnorm(n))
  expect_equal(knn_cmi(x, y), Rcausality:::.knn_cmi_r(x, y), tolerance = 1e-6)
  expect_equal(knn_cmi(x, y, z), Rcausality:::.knn_cmi_r(x, y, z), tolerance = 1e-6)
  nn <- Rcausality:::knn_idx_cpp(z, 4L)
  ref <- t(apply(Rcausality:::.dist_block(z) + diag(Inf, n), 1, order))[, 1:4]
  expect_equal(nn, ref, ignore_attr = TRUE)
  p <- Rcausality:::.local_perm(nn, 4L)
  expect_true(all(p %in% seq_len(n)) && length(unique(p)) > 0.8 * n)   # (almost) a permutation ...
  expect_true(all(mapply(function(i, j) j %in% nn[i, 1:4], seq_len(n), p)))  # ... within the neighbourhoods
})
