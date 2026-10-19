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
