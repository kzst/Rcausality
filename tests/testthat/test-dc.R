test_that("dimensional causality recovers the five logistic-map scenarios", {
  skip_on_cran()      # long-running; run by devtools::test() / NOT_CRAN=true
  truth <- c("x_to_y", "y_to_x", "bidirectional", "common_cause", "independent")
  for (sc in truth) {
    d <- sim_logistic_bdc(2000, sc, seed = 1)
    r <- dimensional_causality(d$x, d$y, m = 4, k_range = c(12, 20, 28, 36, 44), seed = 1)
    expect_s3_class(r, "dimcaus")
    expect_equal(r$best, sc, info = sc)
    expect_equal(sum(r$probs), 1, tolerance = 1e-8)
  }
  expect_output(print(r), "posterior")
  pdf(NULL); on.exit(dev.off()); expect_no_error(plot(r))
})

test_that("method = 'dc' fills a causres with posterior decisions", {
  d <- sim_logistic_bdc(1500, "x_to_y", seed = 3)
  r <- rcausality(d, method = "dc", m = 4, k_range = c(12, 20, 28), seed = 1)
  expect_s3_class(r, "causres")
  expect_equal(r$classes$class, "x_to_y")
  expect_true(all(is.na(r$pairs$p_value)))
  expect_true(r$pairs$post_prob[r$pairs$cause == "x"] > 0.5)
  expect_output(print(r), "Intrinsic dimensions")
  pdf(NULL); on.exit(dev.off()); expect_no_error(plot(r, "all"))
})
