# Changes of version 0.0.13

test_that("lag = \"first\" works with the Minnesota prior, and a numeric order is kept by it", {
  d <- sim_var(60, "x_to_y", seed = 1)
  r0 <- suppressWarnings(rcausality(d, lag = "first", seed = 1))
  w <- character(0)
  r <- withCallingHandlers(rcausality(d, lag = "first", bayes = TRUE, prior = "minnesota", seed = 1,
                                      control = bayes_control(n_iter = 200)),
                           warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  expect_false(any(grepl("NAs introduced|non-missing", w)))
  expect_equal(r$pairs$lag[order(r$pairs$cause)], r0$pairs$lag[order(r0$pairs$cause)])   # the order of the scan
  expect_true(all(is.finite(r$pairs$bf10)))
  # a requested order that the F test cannot carry: reduced with a warning, kept under the Minnesota prior
  d13 <- sim_var(13, "x_to_y", seed = 2)
  w <- character(0)
  f <- withCallingHandlers(rcausality(d13, lag = 2, seed = 1),
                           warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("requested order lag = 2", w))); expect_true(all(f$pairs$lag == 1L))
  w <- character(0)
  m <- withCallingHandlers(rcausality(d13, lag = 2, bayes = TRUE, prior = "minnesota", seed = 1,
                                      control = bayes_control(n_iter = 200)),
                           warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  expect_false(any(grepl("requested order", w))); expect_true(all(m$pairs$lag == 2L))
})

test_that("the Bayesian decision of the df method keeps its side condition (a positive statistic)", {
  d <- sim_var(80, "x_to_y", seed = 3)
  r <- suppressWarnings(suppressMessages(rcausality(d, method = "df", seed = 1, n_perm = 49)))
  r$df$table2 <- "undetermined"                     # (the table of the orders decides nothing: the tests do)
  k <- paste(r$pairs$cause, "->", r$pairs$effect)
  set.seed(1); null <- stats::rnorm(49, -0.5, 0.05)  # replicates far below the statistic
  for (j in seq_along(k)) { r$null[[k[j]]]$null <- null; r$null[[k[j]]]$alternative <- "greater" }
  r$null[[k[1]]]$statistic <- r$pairs$statistic[1] <- -0.2    # no fall of the conditional variance
  r$null[[k[2]]]$statistic <- r$pairs$statistic[2] <- 0.05     # a fall
  out <- suppressWarnings(Rcausality:::.bayes_postprocess(r, bayes_control(), prior_prob = 0.5, bf_threshold = 3,
                                                          prob_threshold = 0.5))
  expect_gt(out$pairs$bf10[1], 3); expect_gt(out$pairs$bf10[2], 3)
  expect_false(out$pairs$significant[1]); expect_true(out$pairs$significant[2])
})
