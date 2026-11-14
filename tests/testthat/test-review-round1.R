# Corrections of version 0.0.4 (first review round of the paper)

test_that("circular-shift surrogates are distinct and the p-value respects their number", {
  set.seed(1)
  s <- Rcausality:::.shift_draws(10L, 60L, 20L)
  expect_length(s, 20L); expect_equal(anyDuplicated(s), 0L); expect_true(all(s >= 10L & s <= 60L))
  expect_identical(Rcausality:::.shift_draws(10L, 60L, 500L), 10:60)          # more replicates than shifts: all shifts
  x <- rnorm(60); y <- rnorm(60)
  te <- transfer_entropy(x, y, n_surr = 1000, seed = 1)
  expect_lt(length(te$null), 60L)                                              # the admissible shifts of 60 points
  expect_gte(te$p_value, 1 / (length(te$null) + 1))
  expect_equal(te$p_value, (1 + sum(te$null >= te$te)) / (length(te$null) + 1))
  cc <- ccm(x, y, E = 2, n_surr = 1000, seed = 1)
  expect_lt(length(cc$null), 60L)
})

test_that("the user is told when no link can pass the correction for multiplicity", {
  set.seed(2); d <- as.data.frame(matrix(rnorm(500), 100, 5))
  expect_warning(r <- suppressMessages(rcausality(d, method = "te", lag = 1, seed = 1)), "cannot pass the holm correction over 20 tests")
  expect_false(any(r$pairs$significant))
  # enough replicates are drawn for many tests (no cap at 999)
  d8 <- as.data.frame(matrix(rnorm(8 * 60), 60, 8))
  r8 <- suppressWarnings(suppressMessages(rcausality(d8, method = "pid", lag = 1, seed = 1, method_args = list(k_target = 1, k_source = 1))))
  expect_equal(length(r8$null[[1]]$null), ceiling(56 / 0.05))
})

test_that("the lag scans use joint replicates (te, cmi, pid)", {
  d <- sim_var(200, "x_to_y", seed = 1)
  r <- suppressMessages(rcausality(d$x, d$y, method = "te", max_lag = 3, n_perm = 60, seed = 1, directions = "forward"))
  nl <- r$null[[1]]
  expect_equal(r$pairs$p_value[1], (1 + sum(nl$null >= nl$statistic)) / (length(nl$null) + 1))
  sc <- Rcausality:::.te_scan(d$x, d$y, 1:3, 1L, 4L, 30L)
  expect_equal(dim(sc$null), c(30L, 3L)); expect_length(sc$obs, 3L)
  # the shifts of a scan are spaced by the number of lags (disjoint windows of lags), so that the maximum over all
  # shifts is not counted once per lag
  expect_true(all(diff(sort(Rcausality:::.shift_draws(20L, 180L, 1000L, by = 3L))) == 3L))
  expect_equal(nrow(Rcausality:::.te_scan(d$x, d$y, 1:3, 1L, 4L, 1000L)$null), length(seq(20L, 180L, by = 3L)))
  # the default scan is limited by the length of the series
  expect_message(r8 <- rcausality(d$x[1:120], d$y[1:120], method = "te", seed = 1, directions = "forward"), "lags 1 to 2 are scanned")
  expect_true(all(r8$by_lag$lag <= 2))
  lg <- nl_pair(200, 1)
  rc <- suppressMessages(rcausality(lg$x, lg$y, method = "cmi", max_lag = 2, n_perm = 30, seed = 1, directions = "forward"))
  expect_true(rc$pairs$significant[1]); expect_length(rc$null[[1]]$null, 30L); expect_equal(rc$pairs$lag[1], 1)
  rp <- suppressMessages(rcausality(d$x, d$y, method = "pid", max_lag = 2, n_perm = 40, seed = 1, directions = "forward"))
  expect_length(rp$null[[1]]$null, 40L)
})

test_that("the PID test permutes the cause within the levels of the other source", {
  set.seed(3); code <- sample(0:3, 200, replace = TRUE)
  p <- Rcausality:::.perm_within(code)
  expect_equal(sort(p), seq_along(code)); expect_identical(code[p], code)
  # s1 depends on s2 and the target on s2 only: no novel information of s1, although s1 and the target are dependent
  pv <- vapply(1:20, function(s) { set.seed(s); s2 <- rnorm(300); s1 <- s2 + 0.5 * rnorm(300); tg <- s2 + rnorm(300)
    marpid(tg, s1, s2, n_surr = 99, .first_only = TRUE)$p_values[["novel_1"]] }, 0)
  expect_lt(mean(pv < 0.05), 0.2)
  pw <- vapply(1:10, function(s) { set.seed(s); s2 <- rnorm(300); s1 <- rnorm(300); tg <- s1 + s2 + rnorm(300)
    marpid(tg, s1, s2, n_surr = 99, .first_only = TRUE)$p_values[["novel_1"]] }, 0)
  expect_gt(mean(pw < 0.05), 0.8)
})

test_that("the posterior probability follows from the Bayes factor and prior_prob for every prior", {
  set.seed(4); x <- rnorm(200); y <- rnorm(200)
  for (pr in c("spike_slab", "normal", "bic", "g")) {
    r <- rcausality(x, y, lag = 8, bayes = TRUE, prior = pr, seed = 1, control = bayes_control(prior = pr, decision = "posterior"))
    expect_equal(r$pairs$post_prob, r$pairs$bf10 / (1 + r$pairs$bf10), tolerance = 1e-8, info = pr)
    expect_false(any(r$pairs$significant), info = pr)        # white noise at order 8: no link by the posterior rule
  }
  b <- bayes_granger(y, x, lag = 3, seed = 1)
  expect_equal(b$post_prob, b$bf10 / (1 + b$bf10)); expect_equal(b$any_included[["prior"]], 1 - 0.5^3)
  expect_output(print(b), "probability that some lag is included")
  b2 <- bayes_granger(y, x, lag = 3, seed = 1, prior_prob = 0.2)
  expect_equal(b2$post_prob, 0.2 * b2$bf10 / (0.2 * b2$bf10 + 0.8))
})

test_that("the surrogate Bayes factor is bounded by its model and has a unit-information scale", {
  set.seed(5); null <- rnorm(200)
  for (z in c(4, 5, 6)) {
    b <- surrogate_bf(z, null, seed = 1)
    expect_true(b$extrapolated)
    expect_lt(b$bf10, 5 * exp(z^2 / 2))                                        # no shift alternative can exceed exp(z^2 / 2)
  }
  expect_false(surrogate_bf(0.5, null, seed = 1)$extrapolated)
  expect_error(surrogate_bf(0, null, effect_scale = "unit"), "n_eff")
  loc <- surrogate_bf(0, null, seed = 1)$bf10
  u1 <- surrogate_bf(0, null, seed = 1, effect_scale = "unit", n_eff = 100)$bf10
  u2 <- surrogate_bf(0, null, seed = 1, effect_scale = "unit", n_eff = 10000)$bf10
  expect_gt(loc, 0.3); expect_lt(u1, loc / 3); expect_lt(u2, u1 / 5)           # evidence for H0 grows with n
  d <- sim_var(150, "independent", seed = 2)
  r <- suppressMessages(rcausality(d, method = "te", lag = 1, bayes = TRUE, seed = 1, control = bayes_control(effect_scale = "unit")))
  r0 <- suppressMessages(rcausality(d, method = "te", lag = 1, bayes = TRUE, seed = 1))
  expect_true(all(r$pairs$bf10 < r0$pairs$bf10))
})

test_that("defective input is reported", {
  d <- sim_var(120, "x_to_y", seed = 1); dn <- d; dn$x[c(5, 6, 40)] <- NA
  expect_message(rcausality(dn), "3 missing or non-finite value\\(s\\) were replaced by linear interpolation \\(x: 3\\)")
  expect_error(rcausality(dn, impute = FALSE), "Missing")
  dm <- d; dm$y[seq(1, 120, by = 4)] <- NA
  expect_warning(suppressMessages(rcausality(dm)), "more than 10% of the values of y")
  expect_warning(suppressMessages(rcausality(d[1:100, ], method = "ccdh", n_surr = 20, seed = 1)), "fewer than the advisable minimum of 300")
  expect_warning(rcausality(d[1:14, ], lag = 5), "requested order lag = 5 was changed")
  expect_error(rcausality(data.frame(a = d$x, b = 2 * d$x + 1)), "Exactly collinear")
  p <- sim_panel(n_units = 6, n_time = 12, seed = 1)
  expect_message(suppressWarnings(rcausality(p, id = "unit", lag = 1, seed = 1)), "column\\(s\\) time are analysed as variables")
  pg <- p[!(p$unit == p$unit[1] & p$time %in% 5:6), ]
  expect_warning(rcausality(pg, id = "unit", time = "time", lag = 1, seed = 1), "gaps in `time`")
})

test_that("the forecast stage is seeded", {
  skip_if_not_installed("ranger")
  d <- sim_var(150, "x_to_y", seed = 1); r <- rcausality(d, lag = 1)
  f1 <- predict(r, h = 2, model = "rf", rf_settings = list(num_trees = 30), n_sim = 20, seed = 3)
  f2 <- predict(r, h = 2, model = "rf", rf_settings = list(num_trees = 30), n_sim = 20, seed = 3)
  expect_identical(f1$forecasts, f2$forecasts); expect_identical(f1$evaluation, f2$evaluation)
  c1 <- forecast_compare(d, methods = "granger", model = "rf", lag = 1, seed = 1, rf_settings = list(num_trees = 30))
  c2 <- forecast_compare(d, methods = "granger", model = "rf", lag = 1, seed = 1, rf_settings = list(num_trees = 30))
  expect_identical(c1$table, c2$table)
})

test_that("a method of a benchmark does not depend on the other methods of the call", {
  a <- suppressMessages(benchmark_methods(c("granger", "te"), scenarios = c("x_to_y", "independent"), n = 150, reps = 2, seed = 1, lag = 1))
  b <- suppressMessages(benchmark_methods("te", scenarios = c("x_to_y", "independent"), n = 150, reps = 2, seed = 1, lag = 1))
  expect_identical(a$runs$predicted[a$runs$method == "te"], b$runs$predicted)
})

test_that("the forest test keeps a true conditional null hypothesis (cause correlated with a control)", {
  skip_on_cran(); skip_if_not_installed("ranger")
  st <- Rcausality:::.method_arg_defaults$rf_granger; st$num_trees <- 60; st$folds <- 0; st$shadow <- FALSE; st$importance <- FALSE
  pv <- vapply(1:12, function(s) { set.seed(s); n <- 200
    z <- as.numeric(arima.sim(list(ar = 0.5), n)); y <- c(0, 0.7 * (z[-n]^2 - 1)) + rnorm(n); x <- z + 0.5 * rnorm(n)
    Rcausality:::.rf_pair(y, x, cbind(z = z), 1L, 19L, st)$p_value }, 0)
  expect_lte(sum(pv <= 0.05), 3)                             # a circular shift of the cause rejected in more than half of such data
  sp <- Rcausality:::.rf_split_cause(cbind(a = rnorm(50)), cbind(c1 = rnorm(50), c2 = rnorm(50)), st)
  expect_equal(Rcausality:::.rf_surr_cause(sp, 0L), sp$fit + sp$scale * sp$res)   # shift 0 returns the cause itself
})

test_that("the Diks-Panchenko test uses conditional shifts, and the normal reference on request", {
  set.seed(1); n <- 150; z <- as.numeric(arima.sim(list(ar = 0.5), n)); x <- z + 0.5 * rnorm(n); y <- c(0, 0.7 * z[-n]) + rnorm(n)
  r <- suppressMessages(rcausality(data.frame(x = x, z = z), data.frame(y = y), method = "diks_panchenko", directions = "forward",
                                   conditional = TRUE, n_perm = 60, seed = 1))
  expect_true(all(r$extra$reference == "conditional shifts"))
  expect_equal(r$pairs$significant, c(FALSE, TRUE))                             # x is a noisy copy of z; z drives y
  expect_length(r$null[["x -> y"]]$null, 60L)
  r1 <- rcausality(x, y, method = "diks_panchenko", seed = 1)                   # bivariate, one lag: conditional shifts too
  expect_true(all(r1$extra$reference == "conditional shifts"))
  r0 <- rcausality(x, y, method = "diks_panchenko", method_args = list(reference = "normal"))
  expect_true(all(r0$extra$reference == "normal"))
  expect_equal(r0$pairs$p_value, stats::pnorm(r0$extra$z, lower.tail = FALSE))
  expect_equal(diks_panchenko(x, y, n_perm = 0)$tests$p_value, diks_panchenko(x, y, reference = "normal")$tests$p_value)
  expect_error(rcausality(x, y, method = "diks_panchenko", method_args = list(reference = "shift")), "conditional")
  d2 <- diks_panchenko(x, y, lag = 2, n_perm = 40, seed = 1)
  expect_true(all(d2$tests$p_value >= 1 / 41))
  sp <- Rcausality:::.cond_split(cbind(z), cbind(x))
  expect_equal(sp$fit + sp$scale * sp$res, cbind(x), ignore_attr = TRUE)
  expect_gt(stats::cor(sp$fit[, 1], z), 0.9)                                    # the dependence on the control stays in the fit
})

test_that("the kNN conditional mutual information uses conditional shifts and keeps the local permutations as an option", {
  d <- sim_var(200, seed = 5)
  r <- suppressMessages(rcausality(d, method = "cmi", seed = 1, n_perm = 60))
  expect_identical(r$settings$method_args$reference, "shift")
  expect_equal(r$pairs$significant[r$pairs$cause == "x"], TRUE)
  expect_true(all(r$pairs$p_value >= 1 / 41))                                   # the shifts of the scan are spaced by its length
  l <- suppressMessages(rcausality(d, method = "cmi", seed = 1, n_perm = 60, method_args = list(reference = "local")))
  expect_identical(l$settings$method_args$reference, "local")
  expect_error(rcausality(d, method = "cmi", method_args = list(reference = "none")), "arg")
  # one set of shifts for all lags of a scan, also when the conditioning vector depends on the lag ("cmi_lag")
  set.seed(3)
  g <- Rcausality:::.cmi_pair(d$y, d$x, NULL, 1:3, 1L, 5L, 5L, 40L, cause_past = TRUE, reference = "shift")
  expect_length(g$null, 40L)
  expect_equal(g$statistic, max(g$by_lag$cmi))
  # shift 0 returns the observed statistic
  X <- cbind(d$x[1:199]); y <- Rcausality:::.jit(cbind(d$y[2:200])); Z <- cbind(d$y[1:199])
  s0 <- Rcausality:::.cmi_shift_null(X, y, Z, 5L, 0L)
  expect_equal(s0$null[1, ], s0$obs)
  expect_output(print(r), "conditional shifts of the cause")
})

test_that("short and strongly autocorrelated series get a warning from the tests that rely on circular shifts", {
  set.seed(1); n <- 100
  p <- data.frame(x = as.numeric(arima.sim(list(ar = 0.9), n)), y = as.numeric(arima.sim(list(ar = 0.9), n)))
  expect_warning(suppressMessages(rcausality(p, method = "te", lag = 1, seed = 1)), "strongly autocorrelated")
  w <- data.frame(x = rnorm(n), y = rnorm(n))
  expect_no_warning(suppressMessages(rcausality(w, method = "te", lag = 1, seed = 1)))
  expect_warning(rcausality(p, method = "granger"), "toda_yamamoto")            # the F test: its own warning
  expect_no_warning(rcausality(p, method = "toda_yamamoto"))                    # neither circular shifts nor the F test
  expect_no_warning(rcausality(w, method = "granger"))
  expect_equal(Rcausality:::.shift_eff(data.frame(a = rep(c(1, -1), 50), b = rnorm(100)), "a", "b"), 80)   # negative autocorrelation: no loss
  se <- Rcausality:::.shift_eff(p, c("x", "y"), c("y", "x"))
  expect_true(all(se < 19) && se[1] == se[2])
  # the warning is given below 2 / alpha effective shifts: 300 time points with coefficient 0.9 are about 25
  set.seed(2); q <- data.frame(x = as.numeric(arima.sim(list(ar = 0.9), 300)), y = as.numeric(arima.sim(list(ar = 0.9), 300)))
  expect_warning(suppressMessages(rcausality(q, method = "te", lag = 1, seed = 1)), "2/alpha")
  set.seed(2); q <- data.frame(x = as.numeric(arima.sim(list(ar = 0.5), 300)), y = as.numeric(arima.sim(list(ar = 0.5), 300)))
  expect_no_warning(suppressMessages(rcausality(q, method = "te", lag = 1, seed = 1)))
})

test_that("the scan of the lags leaves the shifts that the multiplicity correction needs", {
  d <- sim_var(300, seed = 2); d$z <- rnorm(300)
  expect_message(r <- rcausality(d, method = "te", seed = 1, conditional = FALSE), "6 tests at alpha = 0.05 need")
  expect_true(max(r$pairs$lag) <= 2L)
  expect_true(min(r$pairs$p_value) >= 1 / 122)                                # 121 shifts spaced by the two lags
  expect_true(r$pairs$significant[r$pairs$cause == "x" & r$pairs$effect == "y"])
  expect_message(rcausality(d, method = "te", seed = 1, conditional = FALSE, adjust = "none"), "the lags 1 to 6 are scanned")
})

test_that("the normal and spike-and-slab Bayes factors use the restricted regression as base model", {
  set.seed(4); n <- 40; y <- rnorm(n); X0 <- cbind(1, rnorm(n), rnorm(n)); Xc <- cbind(lag1 = rnorm(n), lag2 = rnorm(n))
  ctl <- bayes_control(prior = "normal", n_iter = 40, burnin = 5, chains = 1)
  core <- Rcausality:::.bayes_reg_core
  b1 <- core(y, X0, c(1e4, 100, 100), Xc, ctl)$bf10
  # flat priors on the always-in coefficients: no dependence on their scale, parametrisation or prior variances
  b2 <- core(y, X0 %*% matrix(c(2, 0, 0, 1, 5, 0, -3, 1, 0.1), 3), c(1, 1, 1), Xc, ctl)$bf10
  expect_equal(b1, b2, tolerance = 1e-8)
  # closed form on the projected data with nu0 = n - 3 degrees of freedom
  q <- qr(X0); yp <- qr.resid(q, y); Xp <- qr.resid(q, Xc); nu0 <- n - 3
  Q0 <- sum(yp^2); A <- crossprod(Xp) + diag(2); Q1 <- Q0 - sum(crossprod(Xp, yp) * solve(A, crossprod(Xp, yp)))
  expect_equal(log(b1), -0.5 * as.numeric(determinant(A)$modulus) - (0.01 + nu0 / 2) * (log(0.01 + Q1 / 2) - log(0.01 + Q0 / 2)), tolerance = 1e-8)
  # no degrees of freedom left: an error that names the remedy
  expect_error(core(y[1:4], X0[1:4, ], c(1, 1, 1), Xc[1:4, ], ctl), "minnesota")
  # short series with many always-in coefficients: BF > 3 is rare without a link (it held for 30% of such series
  # with proper vague priors on the common coefficients)
  bf <- vapply(1:60, function(s) { set.seed(s); Z <- matrix(rnorm(15 * 4), 15)
    bayes_granger(Z[, 1], Z[, 2], Z[, 3:4], lag = 2, control = ctl)$bf10 }, 0)
  expect_lte(mean(bf > 3), 0.2)
})

test_that("dfe warns for series that are not binary", {
  set.seed(1)
  expect_warning(Rcausality:::.note_series(data.frame(a = rnorm(600), b = rnorm(600)), "dfe"), "not binary")
  d <- sim_markov_pair(600, "x_to_y", seed = 1)
  expect_no_warning(Rcausality:::.note_series(d[, c("x", "y")], "dfe"))
})

test_that("the kernel test warns for strongly autocorrelated series unless they are filtered", {
  set.seed(3); n <- 300
  d <- data.frame(x = as.numeric(arima.sim(list(ar = 0.95), n)), y = as.numeric(arima.sim(list(ar = 0.95), n)))
  expect_warning(rcausality(d, method = "diks_panchenko", lag = 1, n_perm = 40, seed = 1), "filter = \"var\"")
  expect_warning(rcausality(d, method = "diks_panchenko", lag = 1, method_args = list(reference = "normal")), "normal reference")
  expect_no_warning(rcausality(d, method = "diks_panchenko", lag = 1, n_perm = 40, seed = 1, method_args = list(filter = "var")))
  expect_no_warning(rcausality(d, method = "diks_panchenko", lag = 1, method_args = list(filter = "var", reference = "normal")))
  e <- data.frame(x = rnorm(n), y = rnorm(n))
  expect_no_warning(rcausality(e, method = "diks_panchenko", lag = 1, n_perm = 40, seed = 1))
})

test_that("a target analysis limits the lag scan of a pair by the size of its own family of tests", {
  d <- sim_drivers(300, seed = 1); msgs <- character(0)
  a <- withCallingHandlers(causal_drivers(d, target = "y", method = "te", seed = 1),
                           message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  nl <- table(paste(a$by_lag$cause, a$by_lag$effect))
  expect_equal(unname(nl[["x1 y"]]), 3L)        # links to the target: 4 tests, 80 shifts, lags 1 to 3
  expect_equal(unname(nl[["y x1"]]), 3L)        # links from the target: 4 tests
  expect_equal(unname(nl[["x1 x2"]]), 1L)       # links among the drivers: 12 tests, 240 shifts, lag 1
  expect_true(any(grepl("families of 4 tests", msgs)) && any(grepl("families of 12 tests", msgs)))
  # an analysis without a target is not affected
  b <- suppressMessages(rcausality(d[, c("x1", "y")], method = "te", seed = 1))
  expect_equal(length(unique(table(paste(b$by_lag$cause, b$by_lag$effect)))), 1L)
})

test_that("the default dfe model does not depend on the state of the session", {
  skip_if_not_installed("ranger")
  e <- Rcausality:::.dfe_env; old <- e$model; on.exit(e$model <- old, add = TRUE)
  d <- sim_markov_pair(400, "x_to_y", seed = 1)
  sc <- function(seed) { e$model <- NULL; set.seed(seed)
    m <- suppressMessages(Rcausality:::.dfe_model(NULL, n_per_class = 6, lengths = 300L, num_trees = 30L))
    u <- stats::runif(1)                                   # drawn before the forest predicts (which draws a seed)
    list(scores = dfe_classify(d$x, d$y, model = m)$scores, next_draw = u) }
  a <- sc(1); b <- sc(99)
  expect_equal(a$scores, b$scores)                         # the same model whatever the state of the generator
  set.seed(1); expect_equal(a$next_draw, stats::runif(1))  # and the stream of the session is untouched
  m1 <- dfe_train(n_per_class = 6, lengths = 300L, num_trees = 30L, seed = 5, cores = 1)
  m2 <- dfe_train(n_per_class = 6, lengths = 300L, num_trees = 30L, seed = 5, cores = 2)
  expect_equal(dfe_classify(d$x, d$y, model = m1)$scores, dfe_classify(d$x, d$y, model = m2)$scores)
})
