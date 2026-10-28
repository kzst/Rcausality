test_that("recurrence-plot (Hirata) and correlation-dimension (Krakovska) methods recover logistic-map scenarios", {
  for (sc in c("x_to_y", "common_cause", "independent")) {
    d <- sim_logistic_bdc(1200, sc, seed = 3)
    h <- hirata_rp(d$x, d$y, E = 3, n_surr = 40, seed = 1)
    expect_s3_class(h, "hirata"); expect_equal(h$verdict, sc, info = paste("hirata", sc))
    k <- krakovska_cd(d$x, d$y, E = 3, n_boot = 10, seed = 1)
    expect_s3_class(k, "krakovska"); expect_equal(k$verdict, sc, info = paste("cd", sc))
  }
  expect_output(print(h), "verdict"); expect_output(print(k), "D2")
  cd <- corr_dim(embed_delay(sim_logistic(800, "independent", seed = 1)$x, 3), n_boot = 5)
  expect_gt(cd$d2, 0.8); expect_lt(cd$d2, 1.3)
  d <- sim_logistic_bdc(1200, "x_to_y", seed = 3)
  r1 <- rcausality(d, method = "hirata", m = 3, n_perm = 40, seed = 1); expect_equal(r1$classes$class, "x_to_y")
  r2 <- rcausality(d, method = "cd", m = 3); expect_equal(r2$classes$class, "x_to_y")
  expect_output(print(r1), "Per-pair"); expect_output(print(r2), "D2")
})
