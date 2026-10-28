test_that("Markov test: AR(1) accepted at m = 1, AR(2) rejected at m = 1", {
  set.seed(4)
  x1 <- as.numeric(stats::arima.sim(list(ar = .6), 500))
  expect_gt(markov_test(x1, m = 1, n_perm = 49, seed = 3)$p_value, 0.05)
  set.seed(1)
  x2 <- as.numeric(stats::arima.sim(list(ar = c(.2, .6)), 800))
  expect_lt(markov_test(x2, m = 1, n_perm = 49, seed = 1)$p_value, 0.05)
  expect_gt(markov_test(x2, m = 2, n_perm = 49, seed = 1)$p_value, 0.01)
})

test_that("turn-wheel: driver Markov, driven not (discrete exact test)", {
  w <- sim_turnwheel(5000, seed = 2)
  expect_true(all(w$bob %in% 0:2))
  expect_gt(markov_test(w$alice, m = 1, n_perm = 99, seed = 1)$p_value, 0.05)
  expect_lt(markov_test(w$bob, m = 1, n_perm = 99, seed = 1)$p_value, 0.05)
  expect_equal(dim(sim_logistic3(200, "mutual", seed = 1)), c(200L, 3L))
})

test_that("method = 'markov' returns a decision table in the causres", {
  d <- sim_markov01(300, scenario = "x_to_y", seed = 1)
  r <- rcausality(d, method = "markov", n_perm = 49, n_markov = 29, seed = 1)
  expect_s3_class(r, "causres"); expect_false(is.null(r$markov))
  expect_output(print(r), "Markov property")
})
