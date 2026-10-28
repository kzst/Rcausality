test_that("NSRPS compression reproduces the worked example of the paper", {
  x <- c(1, 2, 2, 1, 2, 1, 1, 2, 2, 1, 2, 1, 2, 2, 1, 1)      # 1221211221212211
  cr <- Rcausality:::nsrps_cr_cpp(x, 3L)
  expect_equal(cr[1], 11 / 16)                                 # (12) -> new symbol: 32313233211
  expect_true(all(diff(cr) < 0))
  expect_equal(Rcausality:::nsrps_cr_cpp(c(2, 1, 1, 1, 1, 1, 2, 2, 1), 1L), 7 / 9)  # 211111221 -> 2(11)(11)1221
  expect_equal(Rcausality:::nsrps_cr_cpp(c(1, 1, 1, 1), 2L), c(1, 1))            # zero entropy: nothing to do
})

test_that("ccdh separates driving, hidden common driver and mutual coupling (logistic maps)", {
  d <- sim_logistic3(1500, "driver", b = 0.4, seed = 3)
  r <- ccdh(d$x0, d$x1, m = 3, n_surr = 40, seed = 1)
  expect_s3_class(r, "ccdh"); expect_equal(r$verdict, "x_to_y")
  expect_lt(r$index[["y"]], r$index[["x"]])
  expect_output(print(r), "verdict")
  pdf(NULL); on.exit(dev.off()); expect_no_error(plot(r))
  r2 <- ccdh(d$x1, d$x2, m = 3, n_surr = 40, seed = 1)       # both driven by the hidden x0
  expect_equal(r2$verdict, "common_cause")
  dm <- sim_logistic3(1500, "mutual", c = 0.4, seed = 3)
  expect_equal(ccdh(dm$x1, dm$x2, m = 3, n_surr = 40, seed = 1)$verdict, "bidirectional")
})

test_that("method = 'ccdh' integrates with causres", {
  d <- sim_logistic3(1500, "driver", b = 0.4, seed = 5)
  r <- rcausality(d[, c("x0", "x1")], method = "ccdh", m = 3, n_surr = 40, seed = 1)
  expect_s3_class(r, "causres"); expect_equal(r$classes$class, "x_to_y")
  expect_output(print(r), "compression rates")
})
