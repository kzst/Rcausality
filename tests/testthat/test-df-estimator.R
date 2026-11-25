# The Gaussian conditional variance and the minimal order search of method "df"

test_that("df with the Gaussian conditional variance finds the hidden driver of the rotation system", {
  # Telcs et al. (2025, Section 6.4): one component of each oscillator observed, orders 4, 4 and 6
  d <- sim_df_rotation(1000, "deterministic", seed = 1)[, c("x1", "y1")]
  g <- suppressWarnings(suppressMessages(rcausality(d, method = "df", seed = 1,
                                                     method_args = list(estimator = "gaussian"))))
  expect_equal(c(g$df$df_x, g$df$df_y, g$df$df_joint), c(4L, 4L, 6L))
  expect_equal(g$df$table2, "common_cause")
  expect_true(g$instant$significant[1]); expect_false(any(g$pairs$significant))
  expect_error(rcausality(d, method = "df", method_args = list(estimator = "normal")), "estimator")
  expect_error(rcausality(d, method = "df", method_args = list(search = "full")), "search")
  # the step tests of df_variance() with the same estimator: an AR(2) process has two degrees of freedom
  set.seed(1); x <- as.numeric(stats::arima.sim(list(ar = c(0.5, 0.3)), 300))
  v <- df_variance(x, max_constr = 3, estimator = "gaussian", seed = 1)
  expect_equal(v$df, 2L); expect_identical(v$estimator, "gaussian")
  expect_output(print(v), "Gaussian conditional variance")
})

test_that("the minimal order search needs fewer conditions than the greedy one where the other series carries a hidden driver", {
  d <- sim_df_rotation(1500, "deterministic", seed = 1)
  xs <- d$x1 / stats::sd(d$x1); ys <- d$y1 / stats::sd(d$y1); t_idx <- 5:1500
  cands <- list(x1 = xs, y1 = ys); P <- cbind(xs[t_idx], ys[t_idx])
  ord <- function(search) { set.seed(1)
    Rcausality:::.df_order(P, cands, t_idx, 4L, 10L, 0.05, 0.02, 49L, 5L, 1L, "gaussian", search) }
  gr <- ord("greedy"); mi <- ord("minimal")
  expect_equal(gr$df, 7L); expect_equal(mi$df, 6L)
  expect_true(all(mi$lags >= 2L))                    # (both pasts: the split with the smallest variance, here 4 and 2)
  expect_setequal(mi$chosen, c(sprintf("x1[t-%d]", seq_len(mi$lags[["x1"]])), sprintf("y1[t-%d]", seq_len(mi$lags[["y1"]]))))
  # the default of the windows is the greedy search; the order of an independent white noise is 0 under both
  set.seed(2); w <- stats::rnorm(400); z <- stats::rnorm(400); ti <- 5:400
  for (s in c("greedy", "minimal"))
    expect_equal(Rcausality:::.df_order(cbind(w[ti]), list(w = w, z = z), ti, 4L, 10L, 0.05, 0.02, 49L, 5L, 1L, "knn", s)$df, 0L)
  expect_identical(Rcausality:::.df_search_arg(NULL, "knn"), "greedy")
  expect_identical(Rcausality:::.df_search_arg("auto", "gaussian"), "minimal")
})
