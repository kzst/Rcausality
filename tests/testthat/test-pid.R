test_that("MAR-PID atoms: XOR is synergy, copies are redundancy, one informative source is unique", {
  set.seed(1); a <- rbinom(3000, 1, .5); b <- rbinom(3000, 1, .5)
  r <- marpid((a + b) %% 2, a, b, 1, 1)
  expect_gt(r$atoms[["synergy"]], 0.95); expect_lt(r$atoms[["redundant"]] + r$atoms[["unique_1"]] + r$atoms[["unique_2"]], 0.05)
  r <- marpid(a, a, a, 1, 1); expect_gt(r$atoms[["redundant"]], 0.99)
  r <- marpid(a, a, b, 1, 1); expect_gt(r$atoms[["unique_1"]], 0.99); expect_lt(r$atoms[["unique_2"]], 0.01)
  expect_true(all(r$atoms >= 0)); expect_output(print(r), "MAR-PID")
  # four-valued target (Sect. 6.3): bit 1 unique to C1, bit 2 unique to C2
  C1 <- rbinom(3000, 1, .5); C2 <- rbinom(3000, 1, .5)
  r4 <- marpid(2 * C1 + C2, C1, C2, k_target = 2, k_source = 1)
  expect_gt(r4$by_bit["bit1", "unique_1"], 0.95); expect_gt(r4$by_bit["bit2", "unique_2"], 0.95)
})

test_that("method = 'pid' finds the directed link of a VAR pair", {
  d <- sim_var(800, "x_to_y", seed = 2)
  r <- rcausality(d, method = "pid", n_perm = 60, seed = 1)
  expect_s3_class(r, "causres"); expect_true(r$pairs$significant[r$pairs$cause == "x"])
  expect_output(print(r), "MAR-PID atoms")
  expect_equal(rcausality(sim_var(800, "independent", seed = 3), method = "pid", n_perm = 60, seed = 1)$classes$class, "independent")
})
