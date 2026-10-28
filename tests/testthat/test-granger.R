test_that("Granger F-test agrees with lmtest::grangertest", {
  skip_if_not_installed("lmtest")
  d <- sim_var(300, scenario = "x_to_y", seed = 11)
  r <- rcausality(d$x, d$y, lag = 2, adjust = "none")
  ref <- lmtest::grangertest(y ~ x, order = 2, data = d)
  row <- r$pairs[r$pairs$cause == "x" & r$pairs$effect == "y", ]
  expect_equal(row$statistic, ref$F[2], tolerance = 1e-6)
  expect_equal(row$p_value, ref$`Pr(>F)`[2], tolerance = 1e-6)
})

test_that("direction is recovered in VAR scenarios", {
  d <- sim_var(600, scenario = "x_to_y", seed = 1)
  r <- rcausality(d$x, d$y)
  expect_s3_class(r, "causres")
  expect_equal(r$classes$class, "x_to_y")
  d <- sim_var(600, scenario = "bidirectional", seed = 2)
  expect_equal(rcausality(d$x, d$y)$classes$class, "bidirectional")
  d <- sim_var(600, scenario = "independent", seed = 3)
  expect_equal(rcausality(d$x, d$y)$classes$class, "independent")
})

test_that("bayes_granger reports Bayes factors and decides by threshold", {
  d <- sim_var(500, scenario = "y_to_x", seed = 4)
  r <- rcausality(d$x, d$y, method = "bayes_granger")
  expect_true(all(c("bf10", "post_prob", "evidence") %in% names(r$pairs)))
  expect_equal(r$classes$class, "y_to_x")
})

test_that("data.frame (all pairs) and MIMO block modes work", {
  d <- sim_var(500, n_series = 4, seed = 5)
  r <- rcausality(d)
  expect_equal(nrow(r$pairs), 12)
  expect_equal(dim(as.matrix(r, "lag")), c(4, 4))
  r2 <- rcausality(d[, 1:2], d[, 3:4])
  expect_equal(nrow(r2$pairs), 8)
  expect_equal(nrow(rcausality(d[, 1:2], d[, 3:4], directions = "forward")$pairs), 4)
  rc <- rcausality(d, conditional = TRUE)
  expect_true(all(rc$pairs$df1 == rc$pairs$lag))
})

test_that("print / summary / plot run", {
  d <- sim_var(400, n_series = 3, seed = 6)
  r <- rcausality(d)
  expect_output(print(r), "Granger")
  expect_output(print(summary(r)), "summary")
  pdf(NULL); on.exit(dev.off())
  for (tp in c("network", "heatmap", "bars", "both")) expect_no_error(plot(r, type = tp))
  expect_no_error(plot(r, type = "heatmap", evidence = "bf"))
})

test_that("input validation", {
  expect_error(rcausality(1:50), "two series")
  expect_error(rcausality(1:50, 1:40), "same length")
})

test_that("simulators and prep helpers", {
  for (f in list(sim_var, sim_logistic, sim_tent, sim_henon, sim_markov01))
    for (sc in c("x_to_y", "y_to_x", "bidirectional", "common_cause", "independent")) {
      d <- f(200, scenario = sc, seed = 1)
      expect_equal(dim(d), c(200, 2)); expect_equal(attr(d, "truth"), sc)
      expect_true(all(is.finite(as.matrix(d))))
    }
  expect_equal(dim(embed_delay(1:10, 3, 2)), c(6, 3))
  expect_equal(binarize(c(1, 5, 2, 8)), c(0L, 1L, 0L, 1L))
  expect_equal(length(ordinal_pattern(rnorm(50), 3)), 48)
  expect_false(anyNA(prep_series(c(1, NA, 3))))
})

test_that("README quick-start works and a function input gives a helpful error", {
  d <- sim_var(300, scenario = "x_to_y", seed = 1)
  r1 <- rcausality(d$x, d$y); r2 <- rcausality(d)
  expect_equal(r1$classes$class, "x_to_y"); expect_equal(r2$classes$class, "x_to_y")
  expect_error(rcausality(stats::df), "function, not data")
  net <- sim_var(300, n_series = 4, seed = 42)
  expect_equal(nrow(rcausality(net[, 1:2], net[, 3:4])$pairs), 8)
})
