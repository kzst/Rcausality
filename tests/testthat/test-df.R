test_that("df_variance recovers the order of AR processes", {
  set.seed(2)
  x2 <- as.numeric(stats::arima.sim(list(ar = c(.5, .3)), 400))
  d <- df_variance(x2, max_constr = 3, n_perm = 49, seed = 1)
  expect_s3_class(d, "dfvar"); expect_equal(d$df, 2L)
  expect_output(print(d), "Degrees of freedom")
  pdf(NULL); on.exit(dev.off()); expect_no_error(plot(d))
  x1 <- as.numeric(stats::arima.sim(list(ar = .6), 400))
  expect_equal(df_variance(x1, max_constr = 3, n_perm = 49, seed = 1)$df, 1L)
})

test_that("deterministic rotation: z needs one constrained pair (profile reading)", {
  d <- sim_df_rotation(600, "deterministic", seed = 1)
  dz <- df_variance(d[, c("z1", "z2")], max_constr = 2, n_perm = 49, seed = 1)
  expect_equal(dz$df, 1L); expect_equal(dz$df_profile, 1L)
})

test_that("method = 'df' finds x -> y and fills the orders table", {
  d <- sim_var(300, "x_to_y", seed = 1)
  r <- rcausality(d, method = "df", n_perm = 60, seed = 1)
  expect_s3_class(r, "causres")
  expect_true(r$pairs$significant[r$pairs$cause == "x"])
  expect_false(r$pairs$significant[r$pairs$cause == "y"])
  expect_equal(r$df$table2, "x_to_y")
  expect_output(print(r), "Table-2")
  pdf(NULL); on.exit(dev.off()); expect_no_error(plot(r, "all"))
})

test_that("parallel runs reproduce with a seed and restore the RNG kind", {
  d <- sim_logistic(200, "x_to_y", seed = 1)
  kind <- RNGkind()[1]
  r1 <- rcausality(d, method = "cmi", n_perm = 60, seed = 1, cores = 2)
  r2 <- rcausality(d, method = "cmi", n_perm = 60, seed = 1, cores = 2)
  expect_equal(r1$pairs$p_value, r2$pairs$p_value)
  expect_equal(RNGkind()[1], kind)
})

test_that("new plot types and plot_embedding work", {
  d <- sim_var(200, "x_to_y", seed = 1); r <- rcausality(d)
  pdf(NULL); on.exit(dev.off())
  for (tp in c("matrix", "lags", "embedding", "all")) expect_no_error(plot(r, tp))
  expect_equal(ncol(plot_embedding(d$x, m = 3)), 3L)
  expect_no_error(plot_embedding(d$x, d$y))
  expect_false(is.null(r$by_lag))
})
