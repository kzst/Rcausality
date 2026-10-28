test_that("binary conditional entropy and the coupled-chain simulator behave", {
  skip_if_not_installed("ranger")
  set.seed(1)
  x <- rbinom(2000, 1, 0.5)
  H0 <- Rcausality:::cond_entropy_bits_cpp(x, matrix(0L, 2000, 0), 2L)
  expect_equal(H0, 1, tolerance = 0.01)                               # fair coin: 1 bit
  d <- sim_markov_pair(3000, "x_to_y", df_x = 1, df_y = 1, seed = 2)
  E <- Rcausality:::.ce_curves(d$x, d$y, 3L, 2L)
  expect_true(E$Y[2, 1] >= E$Y[2, 2] - 1e-9)                         # x's past reduces H(y | y past)
  expect_lt(abs(E$X[2, 1] - E$X[2, 3]), 0.05)                         # but y's past does not help x
})

test_that("dfe_train / dfe_classify recover simulated Markov classes", {
  skip_if_not_installed("ranger")
  m <- dfe_train(n_per_class = 40, lengths = c(500, 1000), num_trees = 50, seed = 1)
  expect_s3_class(m, "dfe_model"); expect_output(print(m), "random-forest")
  d <- sim_markov_pair(2000, "x_to_y", df_x = 1, df_y = 2, seed = 5)
  r <- dfe_classify(d$x, d$y, model = m)
  expect_s3_class(r, "dfe"); expect_equal(sum(r$scores), 1, tolerance = 1e-8)
  expect_equal(r$best, "x_to_y")
  expect_output(print(r), "class scores")
  pdf(NULL); on.exit(dev.off()); expect_no_error(plot(r))
  d2 <- sim_var(1500, "x_to_y", seed = 1)
  rc <- rcausality(d2, method = "dfe", model = m)
  expect_s3_class(rc, "causres"); expect_true(all(is.na(rc$pairs$p_value)))
  expect_output(print(rc), "class scores")
})
