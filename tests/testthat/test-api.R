# The unified interface of version 0.0.3: method, bayes, cores as options;
# method settings through `method_args`; the method catalogue

test_that("causal_methods() lists the methods and their settings", {
  cm <- causal_methods()
  expect_s3_class(cm, "causal_methods")
  expect_setequal(cm$method, Rcausality:::.method_names)
  expect_true(all(c("type", "conditional", "pooled", "min_length") %in% names(cm)))
  expect_true(all(cm$pooled[cm$method %in% c("granger", "toda_yamamoto", "cmi", "te", "ccm", "rf_granger", "diks_panchenko")]))
  expect_output(print(cm), "granger")
  one <- causal_methods("te")
  expect_s3_class(one, "causal_method")
  expect_true(all(c("bins", "history") %in% names(one$settings)))
  expect_output(print(one), "bins")
  expect_error(causal_methods("nope"))
})

test_that("method settings: `method_args` and the settings given directly agree", {
  d <- sim_logistic(150, "x_to_y", seed = 2)
  a <- rcausality(d, method = "te", bins = 6, n_perm = 40, seed = 1)
  b <- rcausality(d, method = "te", method_args = list(bins = 6), n_perm = 40, seed = 1)
  expect_equal(a$pairs, b$pairs)
  expect_equal(a$args$method_args$bins, 6)
  expect_error(rcausality(d, method = "te", method_args = list(binz = 6)), "binz")
  expect_error(rcausality(d, method = "te", foo = 1), "foo")
  # a setting of another method is not accepted in method_args
  expect_error(rcausality(d, method = "granger", method_args = list(bins = 6)), "bins")
})

test_that("the stored arguments reproduce the analysis", {
  d <- sim_var(200, "x_to_y", seed = 1)
  r <- rcausality(d, method = "toda_yamamoto", method_args = list(d_max = 2))
  r2 <- do.call(rcausality, c(list(x = d), r$args))
  expect_equal(r$pairs, r2$pairs)
  rt <- rcausality(d, method = "te", n_perm = 40, seed = 3, method_args = list(bins = 3))
  rt2 <- do.call(rcausality, c(list(x = d), rt$args))
  expect_equal(rt$pairs, rt2$pairs)
})

test_that("bayes = TRUE and cores are options of every call", {
  d <- sim_var(200, "x_to_y", seed = 1)
  f <- rcausality(d); b <- rcausality(d, bayes = TRUE, prior = "bic")
  expect_false(isTRUE(f$settings$bayes)); expect_true(b$settings$bayes)
  expect_equal(b$classes$class, "x_to_y")
  expect_equal(rcausality(d, cores = 1)$pairs, rcausality(d, cores = 2)$pairs)
})

test_that("lag = 'first': first significant order, Bonferroni over the scan", {
  d <- sim_var(300, "x_to_y", seed = 1)
  r <- rcausality(d, lag = "first", max_lag = 4, adjust = "none")
  bl <- r$by_lag[r$by_lag$cause == "x" & r$by_lag$effect == "y", ]
  nl <- nrow(bl)
  expect_equal(r$pairs$p_value[r$pairs$cause == "x"], min(1, min(bl$p) * nl))
  expect_equal(r$pairs$lag[r$pairs$cause == "x"], bl$lag[which(bl$p < 0.05 / nl)[1]])
  # the Bayes factor is that of 'a link at one of the orders': the average over the orders
  fs <- Rcausality:::.first_lag(1:3, c(0.5, 0.001, 0.2), c(0.2, 40, 1), 0.05, 3)
  expect_equal(fs$lag, 2L); expect_equal(fs$p_value, 0.003); expect_equal(fs$bf10, mean(c(0.2, 40, 1)))
  fb <- Rcausality:::.first_lag(1:3, c(0.5, 0.001, 0.2), c(5, 40, 1), 0.05, 3, by_bf = TRUE)
  expect_equal(fb$lag, 1L)
  expect_error(rcausality(d, lag = "last"), "lag")
})
