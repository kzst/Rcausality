# Regression tests for the corrections of version 0.0.5 (second round of the review of the manuscript).

test_that("conditional shifts keep the location and the scale of the cause given the conditioning vector", {
  ns <- asNamespace("Rcausality")
  set.seed(1); n <- 400; z <- as.numeric(arima.sim(list(ar = 0.5), n))
  x1 <- z^2 + 0.3 * rnorm(n)                                   # a nonlinear function of the control
  sp <- ns$.cond_split(cbind(z), cbind(x1))
  expect_equal(sp$fit + sp$scale * sp$res, cbind(x1), ignore_attr = TRUE)
  expect_equal(ns$.cond_surr(sp, seq_len(n)), cbind(x1), ignore_attr = TRUE)      # shift 0 returns the cause
  expect_gt(stats::cor(sp$fit[, 1], z^2), 0.97); expect_gt(sp$r2, 0.8)
  s <- ns$.cond_surr(sp, c(101:n, 1:100))                       # a surrogate is as good a proxy of z^2 as the cause is
  expect_gt(stats::cor(s[, 1], z^2), 0.9 * stats::cor(x1, z^2))
  x2 <- (0.3 + abs(z)) * rnorm(n)                               # the control drives the scale of the cause only
  sp2 <- ns$.cond_split(cbind(z), cbind(x2))
  expect_gt(stats::cor(sp2$scale[, 1], abs(z)), 0.6)
  expect_lt(abs(stats::cor(abs(sp2$res[, 1]), abs(z))), 0.2)    # the standardized residual does not depend on it
  # the absolute values of the conditioning variables are among the candidate bases of the scale (version 0.0.7): a
  # scale that is proportional to |z| is followed closely, down to a tenth of the mean absolute residual
  B <- ns$.scale_bases(cbind(z), ns$.poly_bases(cbind(z)))
  expect_identical(ncol(B[[1]]), 1L); expect_equal(B[[2]], cbind(1, abs(z)), ignore_attr = TRUE)
  expect_identical(ns$.scale_bases(matrix(0, n, 0L), list(matrix(1, n, 1L))), list(matrix(1, n, 1L)))
  x4 <- z * rnorm(n); sp4 <- ns$.cond_split(cbind(z), cbind(x4))
  expect_gt(stats::cor(sp4$scale[, 1], abs(z)), 0.9)
  expect_gte(min(sp4$scale[, 1]), 0.1 * mean(abs(x4 - sp4$fit[, 1])) - 1e-12)
  expect_lt(abs(stats::cor(abs(sp4$res[, 1]), abs(z))), 0.2)
  # the weighted refit leaves a constant scale at the mean absolute residual
  expect_equal(ns$.scale_fit(list(matrix(1, n, 1L)), abs(x4), 0.1), rep(mean(abs(x4)), n))
  x3 <- rnorm(n); sp3 <- ns$.cond_split(cbind(z), cbind(x3))    # unrelated cause: (almost) the circular shift itself
  expect_lt(sp3$r2, 0.1); expect_gt(stats::cor(ns$.cond_surr(sp3, c(101:n, 1:100))[, 1], x3[c(101:n, 1:100)]), 0.95)
  expect_true(ns$.determined(0.995)); expect_false(ns$.determined(0.9)); expect_false(ns$.determined(NA_real_))
})

test_that("the conditional tests keep a true null hypothesis when the cause is a nonlinear or a scale copy of a control", {
  skip_on_cran()      # long-running; run by devtools::test() / NOT_CRAN=true
  gen <- function(s, kind, n = 200) { set.seed(s); z <- as.numeric(arima.sim(list(ar = 0.5), n))
    x <- if (kind == "nl") z^2 + 0.3 * rnorm(n) else (0.3 + abs(z)) * rnorm(n)
    y <- c(0, if (kind == "nl") 0.7 * (z[-n]^2 - 1) else 0.8 * abs(z[-n])) + rnorm(n)
    list(x = data.frame(x = x, z = z), y = data.frame(y = y)) }
  for (kind in c("nl", "scale")) {
    pv <- vapply(1:16, function(s) { d <- gen(s, kind)
      r <- suppressWarnings(suppressMessages(rcausality(d$x, d$y, method = "diks_panchenko", directions = "forward", conditional = TRUE,
                                                        n_perm = 60, seed = 1)))
      r$pairs$p_value[r$pairs$cause == "x"] }, 0)
    expect_lte(sum(pv <= 0.05), 4, label = paste("kernel test,", kind))   # linear conditional means: 85 to 99% of such data
  }
  pv <- vapply(1:10, function(s) { d <- gen(s, "nl")
    r <- suppressWarnings(suppressMessages(rcausality(d$x, d$y, method = "cmi", directions = "forward", conditional = TRUE, lag = 1,
                                                      n_perm = 40, seed = 1)))
    r$pairs$p_value[r$pairs$cause == "x"] }, 0)
  expect_lte(sum(pv <= 0.05), 3)
})

test_that("no conditional test is carried out for a cause that the past of the effect determines", {
  skip_on_cran()
  lg <- sim_logistic(300, "x_to_y", seed = 1)                   # noise-free maps: x[t] is a polynomial of y[t], y[t - 1]
  expect_warning(r <- suppressMessages(rcausality(lg, method = "cmi", seed = 1, n_perm = 40)), "No test is carried out")
  expect_equal(r$pairs$p_value[r$pairs$cause == "x"], 1); expect_false(any(r$pairs$significant))
  expect_gt(r$null[["x -> y"]]$cond_r2, 0.99)
  b <- suppressWarnings(suppressMessages(rcausality(lg, method = "cmi", seed = 1, n_perm = 40, bayes = TRUE)))
  expect_equal(b$pairs$bf10[b$pairs$cause == "x"], 1); expect_match(b$bayes$model[b$bayes$pair == "x -> y"], "not tested")
  # the test with plain shifts and the state-space methods are the tools for such systems ...
  expect_equal(suppressMessages(rcausality(lg, method = "te", seed = 1, method_args = list(reference = "circular")))$classes$class, "x_to_y")
  # (the default reference of "te", conditional shifts from version 0.0.12, is a conditional test as well)
  expect_warning(suppressMessages(rcausality(lg, method = "te", seed = 1)), "No test is carried out")
  # ... and observational noise makes the conditional test informative again
  ln <- sim_logistic(300, "x_to_y", seed = 1, noise = 0.05)
  expect_no_warning(rn <- suppressMessages(rcausality(ln, method = "cmi", seed = 1)))
  expect_equal(rn$classes$class, "x_to_y")
})

test_that("the panel bootstrap regenerates the cause with the effect and can resample common periods", {
  expect_identical(panel_control()$bootstrap, "joint"); expect_identical(panel_control()$resample, "units")
  expect_identical(panel_control(resample = "periods")$bootstrap, "effect")       # the Dumitrescu-Hurlin bootstrap
  expect_identical(panel_control(resample = "periods", bootstrap = "joint")$bootstrap, "joint")
  expect_error(panel_control(bootstrap = "none"), "arg")
  pd <- sim_panel(n_units = 12, n_time = 14, seed = 4)
  run <- function(...) suppressMessages(rcausality(pd, id = "unit", time = "time", seed = 1, n_perm = 39, panel = panel_control(...)))
  a <- run(); b <- run(bootstrap = "effect"); p <- run(resample = "periods")
  expect_identical(a$panel$p_value, run()$panel$p_value)                          # reproducible
  expect_equal(a$units$statistic, b$units$statistic); expect_equal(a$units$statistic, p$units$statistic)   # the same statistics,
  expect_false(isTRUE(all.equal(a$panel$p_value, b$panel$p_value)))               # other reference distributions
  expect_true(all(is.finite(p$panel$p_value)))
  pu <- pd[!(pd$unit == pd$unit[1] & pd$time > 10), ]
  expect_warning(suppressMessages(rcausality(pu, id = "unit", time = "time", seed = 1, n_perm = 39, panel = panel_control(resample = "periods"))),
                 "balanced")
})

test_that("the reverse direction of a strong link is not declared in panels of many short units", {
  skip_on_cran()      # long-running; run by devtools::test() / NOT_CRAN=true
  gen <- function(s, N = 60, TT = 10) { set.seed(s)
    do.call(rbind, lapply(seq_len(N), function(i) { x <- numeric(TT + 30); y <- x
      for (t in 2:(TT + 30)) { x[t] <- 0.5 * x[t - 1] + rnorm(1); y[t] <- 0.5 * y[t - 1] + 0.5 * x[t - 1] + rnorm(1) }
      data.frame(unit = i, time = seq_len(TT), x = x[-(1:30)], y = y[-(1:30)]) })) }
  pv <- vapply(1:12, function(s) {
    r <- suppressWarnings(suppressMessages(rcausality(gen(s), id = "unit", time = "time", lag = 1, seed = s, n_perm = 99, adjust = "none")))
    P <- r$panel; c(P$p_value[P$cause == "y" & P$level == "pooled"], P$p_value[P$cause == "x" & P$level == "pooled"]) }, c(0, 0))
  expect_lte(sum(pv[1, ] <= 0.05), 3)                           # y -> x: no link
  expect_gte(sum(pv[2, ] <= 0.05), 11)                          # x -> y: the link
})

test_that("time effects are removed from the series standardized within the units", {
  ns <- asNamespace("Rcausality"); set.seed(2)
  f <- as.numeric(arima.sim(list(ar = 0.5), 20)); sc <- c(1, 1, 1, 50)          # one unit fifty times as large
  d <- do.call(rbind, lapply(1:4, function(i) data.frame(unit = i, time = 1:20, x = sc[i] * (f + rnorm(20, 0, 0.3)),
                                                         y = sc[i] * (f + rnorm(20, 0, 0.3)))))
  U <- ns$.panel_units(d, "unit", "time", time_effects = TRUE)$units
  X <- vapply(U, function(u) u$x, numeric(20))
  expect_equal(rowMeans(X), rep(0, 20), tolerance = 1e-10)                      # centred within every period
  expect_lt(max(abs(stats::cor(X, f))), 0.75)                                   # the factor is removed from all units,
  expect_lt(diff(range(apply(X, 2, stats::sd))), 0.5)                           # ... which stay on one scale
})

test_that("placebo offsets keep every pair of series apart", {
  ns <- asNamespace("Rcausality"); set.seed(1)
  for (K in c(2L, 3L, 6L)) for (i in 1:50) {
    off <- ns$.spaced_offsets(K, 200L, 20L); dd <- abs(outer(off, off, "-")); dd <- pmin(dd, 200L - dd)
    expect_gte(min(dd[upper.tri(dd)]), 20L)
  }
  expect_equal(ns$.rot(1:5, 2L), c(3L, 4L, 5L, 1L, 2L)); expect_equal(ns$.rot(1:5, 0L), 1:5)
  D <- data.frame(a = 1:100, b = 1:100, c = 1:100); P <- ns$.placebo_series(D, 0.1)
  sh <- c(P$a[1] - P$b[1], P$a[1] - P$c[1], P$b[1] - P$c[1]) %% 100L
  expect_true(all(pmin(sh, 100L - sh) >= 10L))
})

test_that("the default classifier of dfe does not depend on the kind or the state of the generator", {
  skip_on_cran(); skip_if_not_installed("ranger")
  ns <- asNamespace("Rcausality"); old_kind <- RNGkind(); on.exit(do.call(RNGkind, as.list(old_kind)), add = TRUE)
  env <- ns$.dfe_env; keep <- env$model; on.exit(env$model <- keep, add = TRUE)
  d <- sim_markov_pair(600, "x_to_y", seed = 1)
  score <- function(kind, seed) { env$model <- NULL; RNGkind(kind); set.seed(seed); before <- .Random.seed
    m <- suppressMessages(ns$.dfe_model()); expect_identical(.Random.seed, before)       # the caller's stream is untouched
    set.seed(1); suppressMessages(dfe_classify(d$x, d$y, model = m))$scores }
  a <- score("Mersenne-Twister", 1); b <- score("L'Ecuyer-CMRG", 99)
  expect_identical(a, b)
})

test_that("the data arguments are evaluated before the seed of the analysis is set", {
  set.seed(5); a <- rcausality(rnorm(120), rnorm(120), method = "granger", lag = 1, seed = 1)$pairs$statistic
  set.seed(6); b <- rcausality(rnorm(120), rnorm(120), method = "granger", lag = 1, seed = 1)$pairs$statistic
  set.seed(5); x <- rnorm(120); y <- rnorm(120); cc <- rcausality(x, y, method = "granger", lag = 1, seed = 1)$pairs$statistic
  expect_false(isTRUE(all.equal(a, b))); expect_equal(a, cc)
  set.seed(5); f1 <- sim_var(100, seed = NULL); set.seed(5); f2 <- sim_var(100, seed = NULL); expect_identical(f1, f2)
})

test_that("the F test reports heteroskedastic errors and has a heteroskedasticity-consistent form", {
  d <- nl_pair(300, 1)                                          # no linear link; the variance of the errors depends on x
  expect_warning(g <- rcausality(d, method = "granger", lag = 1), "heteroskedastic")
  expect_lt(min(g$extra$het_p), 0.001)
  expect_no_warning(h <- rcausality(d, method = "granger", lag = 1, method_args = list(vcov = "HC3")))
  expect_true(all(h$extra$vcov == "hc3")); expect_gt(h$pairs$p_value[h$pairs$cause == "x"], g$pairs$p_value[g$pairs$cause == "x"])
  expect_error(rcausality(d, method = "granger", method_args = list(vcov = "HC1")), "vcov")
  # the statistic is the Wald statistic with the HC3 covariance, divided by the number of restrictions
  skip_if_not_installed("sandwich"); skip_if_not_installed("lmtest")
  n <- nrow(d); e <- data.frame(y = d$y[3:n], y1 = d$y[2:(n - 1)], y2 = d$y[1:(n - 2)], x1 = d$x[2:(n - 1)], x2 = d$x[1:(n - 2)])
  fu <- lm(y ~ y1 + y2 + x1 + x2, e); fr <- lm(y ~ y1 + y2, e)
  w <- lmtest::waldtest(fr, fu, vcov = sandwich::vcovHC(fu, type = "HC3"), test = "F")
  h2 <- rcausality(d, method = "granger", lag = 2, method_args = list(vcov = "HC3")); xy <- h2$pairs$cause == "x"
  expect_equal(h2$pairs$statistic[xy], w$F[2], tolerance = 1e-6); expect_equal(h2$pairs$p_value[xy], w$`Pr(>F)`[2], tolerance = 1e-6)
  # ... and the classical statistic is that of lmtest::grangertest()
  g2 <- suppressWarnings(rcausality(d, method = "granger", lag = 2))
  expect_equal(g2$pairs$statistic[g2$pairs$cause == "x"], lmtest::grangertest(y ~ x, order = 2, data = d)$F[2], tolerance = 1e-6)
  expect_equal(g2$pairs$statistic[g2$pairs$cause == "y"], lmtest::grangertest(x ~ y, order = 2, data = d)$F[2], tolerance = 1e-6)
  # Toda-Yamamoto accepts the same covariance
  ty <- suppressWarnings(rcausality(d, method = "toda_yamamoto", lag = 1, method_args = list(vcov = "HC3")))
  expect_true(all(ty$extra$vcov == "hc3"))
})

test_that("the unit statistics of a panel are those of plm::pgrangertest()", {
  skip_on_cran()
  skip_if_not_installed("plm")
  pd <- sim_panel(n_units = 10, n_time = 20, seed = 3)
  r <- suppressMessages(rcausality(pd, id = "unit", time = "time", lag = 1, seed = 1, n_perm = 19))
  pg <- plm::pgrangertest(y ~ x, data = plm::pdata.frame(pd, index = c("unit", "time")), order = 1L, test = "Zbar")
  u <- r$units[r$units$cause == "x" & r$units$effect == "y", ]
  expect_equal(u$statistic, unname(pg$indgranger$Chisq), tolerance = 1e-6)        # one restriction: Wald = F
})

test_that("Bayes factors beyond the range of the surrogates are reported as bounds", {
  set.seed(1); nl <- rnorm(100)
  b <- surrogate_bf(8, nl, seed = 1)
  expect_true(b$extrapolated); expect_true(is.finite(b$bound)); expect_gt(b$bf10, b$bound)
  expect_equal(b$bound, surrogate_bf(max(nl), nl, seed = 1)$bf10, tolerance = 0.05)   # (Monte Carlo integrals)
  expect_output(print(b), "BF10 > ")
  b0 <- surrogate_bf(0.5, nl, seed = 1); expect_false(b0$extrapolated); expect_output(print(b0), "BF10 = ")
  d <- nl_pair(400, 2); r <- suppressMessages(rcausality(d, method = "cmi", bayes = TRUE, seed = 1, lag = 1, n_perm = 60))
  expect_true("bf_bound" %in% names(r$bayes)); expect_true(is.finite(r$bayes$bf_bound[r$bayes$pair == "x -> y"]))
  expect_output(print(r), "at least")
})

test_that("defective calls are refused or reported", {
  d <- sim_var(100, seed = 1)
  expect_error(rcausality(d, lag = 0), "lag"); expect_error(rcausality(d, lag = 1.5), "lag")
  expect_error(rcausality(d, cores = 0), "cores"); expect_error(rcausality(d[1:6, ], lag = 1), "8")
  expect_error(rcausality(cbind(d, time = 1:100), time = "time"), "id")
  expect_message(suppressWarnings(rcausality(cbind(d, year = 1901:2000), lag = 1)), "year")
  expect_message(rcausality(d, method = "bayes_granger", bayes = FALSE), "ignored")
  M <- causal_methods(); nm <- stats::setNames(M$min_length, M$method)
  expect_equal(unname(nm[c("te", "rf_granger", "ccm")]), c(50, 50, 55))
  # a step-up correction over 12 tests with about 80 admissible shifts: several tests must attain the smallest p value
  d4 <- sim_var(100, n_series = 4, seed = 1)
  w <- capture_warnings(suppressMessages(rcausality(d4, method = "te", adjust = "BH", lag = 1, seed = 1)))
  expect_true(any(grepl("attain that value together", w)))
})

test_that("persistent series get warnings from pid and from the panel tests", {
  set.seed(1); n <- 120; x <- as.numeric(arima.sim(list(ar = 0.97), n)); y <- as.numeric(arima.sim(list(ar = 0.97), n))
  w <- capture_warnings(suppressMessages(rcausality(x, y, method = "pid", seed = 1, lag = 1)))
  expect_true(any(grepl("autocorrelated", w)))
  set.seed(2); pd <- do.call(rbind, lapply(1:12, function(i) data.frame(unit = i, time = 1:25, x = cumsum(rnorm(25) + 0.3),
                                                                        y = cumsum(rnorm(25) + 0.3))))
  expect_warning(suppressMessages(rcausality(pd, id = "unit", time = "time", lag = 1, seed = 1, n_perm = 19)), "persistent or trending")
  expect_no_warning(suppressMessages(rcausality(sim_panel(12, 25, seed = 1), id = "unit", time = "time", lag = 1, seed = 1, n_perm = 19)))
})

test_that("the print method names dependence without direction as such", {
  set.seed(5); n <- 300; e1 <- rnorm(n); e2 <- 0.8 * e1 + 0.6 * rnorm(n)
  x <- as.numeric(stats::filter(e1, 0.5, "recursive")); y <- as.numeric(stats::filter(e2, 0.5, "recursive"))
  r <- suppressMessages(rcausality(x, y, method = "cmi_lag", seed = 1))
  out <- paste(capture.output(print(r)), collapse = "\n")
  expect_match(out, "<~>", fixed = TRUE); expect_false(grepl("<-L->", out, fixed = TRUE))
})

test_that("the Bayes factors next to the frequentist tests are exact, and short regressions with prior = \"bic\" are reported", {
  # (correction 8 of version 0.0.4) the unit Bayes factor of a panel is the exact one of the unit-information prior
  ns <- asNamespace("Rcausality")
  pd <- sim_panel(n_units = 8, n_time = 14, seed = 2)
  r <- suppressMessages(rcausality(pd, id = "unit", time = "time", lag = 1, seed = 1, n_perm = 19))
  u1 <- pd[pd$unit == pd$unit[1], ]; g <- ns$.granger_pair(u1$y, u1$x, NULL, 1L, prior = "g", prior_par = list(g_scale = 1, a = 3))
  row <- r$units[r$units$unit == as.character(pd$unit[1]) & r$units$cause == "x", ]
  expect_equal(row$bf10_raw, g$bf10, tolerance = 1e-6); expect_equal(row$p_ftest, g$p_value, tolerance = 1e-6)
  gb <- ns$.granger_pair(u1$y, u1$x, NULL, 1L, prior = "bic")
  expect_false(isTRUE(all.equal(g$bf10, gb$bf10, tolerance = 1e-3)))             # the BIC form differs in short regressions
  expect_warning(suppressMessages(rcausality(u1$x[1:11], u1$y[1:11], bayes = TRUE, prior = "bic", lag = 1)), "degrees of freedom")
})

test_that("the persistence diagnostic covers the conditional form of the kernel test", {
  # (correction 12 of version 0.0.4)
  set.seed(1); n <- 100; x <- as.numeric(arima.sim(list(ar = 0.95), n)); y <- as.numeric(arima.sim(list(ar = 0.95), n))
  expect_warning(suppressMessages(rcausality(x, y, method = "diks_panchenko", lag = 2, n_perm = 40, seed = 1)), "autocorrelated")
  expect_no_warning(suppressMessages(rcausality(rnorm(n), rnorm(n), method = "diks_panchenko", lag = 2, n_perm = 40, seed = 1)))
})

test_that("the placebo test of two short series uses every relative offset once", {
  set.seed(3); n <- 60; x <- as.numeric(arima.sim(list(ar = 0.4), n)); y <- as.numeric(arima.sim(list(ar = 0.4), n))
  r <- rcausality(x, y, lag = 1)
  expect_no_warning(expect_message(pl <- placebo_test(r, n_placebo = 99, seed = 1), "49 relative circular offsets"))  # reruns are silent
  expect_identical(pl$settings$n_placebo, 49L); expect_true(pl$settings$all_offsets)
  expect_true(all(pl$links$p_placebo >= 1 / 50 - 1e-12))                       # the resolution of 49 placebo data sets
  expect_match(gsub("\\s+", " ", paste(capture.output(print(pl)), collapse = " ")), "all relative circular offsets")
  # fewer placebo data sets than offsets: random spaced offsets as before, no message
  expect_no_message(pl2 <- placebo_test(r, n_placebo = 19, seed = 1))
  expect_identical(pl2$settings$n_placebo, 19L); expect_false(pl2$settings$all_offsets)
  ns <- asNamespace("Rcausality")
  expect_identical(ns$.placebo_offsets(2L, 60L, 0.1), 6:54); expect_null(ns$.placebo_offsets(3L, 60L, 0.1))
  # the relative offset of the two series is the requested one
  D <- data.frame(a = 1:60, b = 1:60); P <- ns$.placebo_series(D, 0.1, rel = 7L)
  expect_identical((P$b[1] - P$a[1]) %% 60L, 7L)
})

test_that("dependence between the units of a panel is diagnosed", {
  ns <- asNamespace("Rcausality")
  # the statistic: hand computation for three units
  set.seed(1); E <- matrix(rnorm(60), 20, 3); R <- cor(E)
  cd <- sum(sqrt(20 - 2) * R[upper.tri(R)]) / sqrt(3); d <- 18 * 3 * 2 / 1
  z <- ns$.cd_stat(E)
  expect_equal(z[["cd"]], cd, tolerance = 1e-12); expect_equal(z[["mean_r"]], mean(R[upper.tri(R)]), tolerance = 1e-12)
  q <- d + cd * sqrt(2 * d); expect_equal(z[["p"]], min(1, 2 * min(pchisq(q, d, lower.tail = FALSE), pchisq(q, d))), tolerance = 1e-12)
  expect_true(all(is.na(ns$.cd_stat(E[, 1:2]))))
  # independent units: no warning, the slot holds the statistics; units with a common factor: warning
  pd <- sim_panel(n_units = 12, n_time = 25, seed = 1)
  r <- expect_no_warning(suppressMessages(rcausality(pd, id = "unit", time = "time", lag = 1, seed = 1, n_perm = 19)))
  expect_identical(rownames(r$dependence), c("x", "y")); expect_true(all(r$dependence$p > 0.001))
  set.seed(2); f <- as.numeric(arima.sim(list(ar = 0.5), 25)); pf <- pd
  for (v in c("x", "y")) pf[[v]] <- pf[[v]] + 1.5 * f[match(pf$time, sort(unique(pf$time)))] * stats::sd(pf[[v]])
  expect_warning(rf <- suppressMessages(rcausality(pf, id = "unit", time = "time", lag = 1, seed = 1, n_perm = 19)), "units move together")
  expect_true(all(rf$dependence$p < 0.001)); expect_true(all(rf$dependence$mean_r > 0.3))
  # settings that account for the dependence: no diagnostic
  for (pc in list(panel_control(time_effects = TRUE), panel_control(resample = "periods"))) {
    r2 <- expect_no_warning(suppressMessages(rcausality(pf, id = "unit", time = "time", lag = 1, seed = 1, n_perm = 19, panel = pc)))
    expect_null(r2$dependence)
  }
  # a common linear trend is not dependence: independent random walks with the same drift
  set.seed(3); rw <- do.call(rbind, lapply(1:12, function(i) data.frame(unit = i, time = 1:25, x = cumsum(rnorm(25, 0.5)), y = cumsum(rnorm(25, 0.5)))))
  cdw <- ns$.panel_cd(rw, rw$unit, rw$time, c("x", "y")); expect_true(all(cdw$p > 0.001))
  # without `time` the periods of a balanced panel are matched by position; an unbalanced one is not checked
  expect_equal(ns$.panel_cd(rw, rw$unit, NULL, c("x", "y")), cdw)
  expect_null(ns$.panel_cd(rw[-1, ], rw$unit[-1], NULL, c("x", "y")))
})

test_that("both tables of the print method show a Bayes factor beyond the replicates as a bound", {
  d <- nl_pair(300, seed = 1)
  r <- suppressMessages(suppressWarnings(rcausality(d, method = "cmi", lag = 1, bayes = TRUE, seed = 1, n_perm = 49)))
  skip_if(!any(is.finite(r$bayes$bf_bound)), "no statistic beyond its replicates in this data set")
  out <- capture.output(print(r)); b <- r$bayes$bf_bound[is.finite(r$bayes$bf_bound)][1]
  shown <- paste0(">", formatC(b, digits = 3, format = "g"))
  expect_gte(sum(grepl(shown, out, fixed = TRUE)), 2L)
})
