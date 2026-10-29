skip_on_cran()      # starts Python (reticulate); run by devtools::test()
skip_if_not(pcmci_available(), "tigramite (Python) not available")

test_that("pcmci() recovers a simple VAR link and returns a causres", {
  d <- sim_var(500, "x_to_y", seed = 3)
  r <- pcmci(d, tau_max = 2)
  expect_s3_class(r, "causres")
  expect_equal(r$settings$method, "pcmci")
  P <- r$pairs
  expect_true(P$significant[P$cause == "x" & P$effect == "y"])
  expect_false(P$significant[P$cause == "y" & P$effect == "x"])
  expect_equal(r$classes$class, "x_to_y")
  expect_equal(nrow(r$by_lag), 4L)
  expect_true(all(c("val_matrix", "p_matrix", "q_matrix") %in% names(r$pcmci)))
  expect_output(print(r), "PCMCI")
})

test_that("pcmci() through rcausality() and PCMCI+ contemporaneous links", {
  d <- sim_var(500, n_series = 3, seed = 1)
  r <- rcausality(d, method = "pcmci", max_lag = 3)
  expect_s3_class(r, "causres")
  A <- attr(d, "adjacency")
  found <- r$pairs[r$pairs$significant, c("cause", "effect")]
  truth <- which(A, arr.ind = TRUE)
  expect_true(all(paste(rownames(A)[truth[, 1]], colnames(A)[truth[, 2]]) %in% paste(found$cause, found$effect)))
  r2 <- pcmci(d, tau_max = 2, plus = TRUE)
  expect_false(is.null(r2$instant))
  expect_equal(nrow(r2$instant), 3L)
  # forward-only block analysis restricts the tested pairs
  r3 <- rcausality(d[, 1:2], d[, 3, drop = FALSE], method = "pcmci", directions = "forward", max_lag = 2)
  expect_equal(nrow(r3$pairs), 2L)
  expect_true(all(r3$pairs$effect == "v3"))
})
