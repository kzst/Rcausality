# Panel data: all methods at the unit, group and pooled level

test_that("panel simulators and the long format", {
  pd <- sim_panel(n_units = 6, n_time = 12, seed = 1)
  expect_equal(dim(pd), c(72L, 4L)); expect_named(pd, c("unit", "time", "x", "y"))
  pu <- sim_panel(n_units = 3, n_time = c(8, 10, 12), scenario = "drivers", prevalence = 2 / 3, seed = 1)
  expect_equal(as.vector(table(pu$unit)), c(8L, 10L, 12L)); expect_named(pu, c("unit", "time", "x1", "x2", "x3", "x4", "y"))
  expect_length(attr(pu, "linked"), 2L)
  expect_equal(sim_panel(4, 8, seed = 3), sim_panel(4, 8, seed = 3))
  Y <- matrix(rnorm(40), 10, 4, dimnames = list(NULL, paste0("c", 1:4))); X <- array(rnorm(80), c(2, 4, 10), list(c("a", "b"), NULL, NULL))
  pl <- panel_long(y = Y, array = X)
  expect_equal(dim(pl), c(40L, 5L)); expect_named(pl, c("unit", "time", "y", "a", "b"))
  expect_equal(pl$y[pl$unit == "c2"], Y[, 2]); expect_equal(pl$a[pl$unit == "c3"], X["a", 3, ])
  expect_error(panel_long(Y), "named")
})

test_that("panel_control() checks its arguments", {
  pc <- panel_control(decision = "pooled", combine = "fisher")
  expect_s3_class(pc, "panel_control"); expect_output(print(pc), "fisher")
  expect_true("units" %in% panel_control(levels = "group")$levels)
  expect_error(panel_control(unit_adjust = "nope")); expect_error(panel_control(prevalence_prior = c(1, -1)))
  expect_error(panel_control(decision = "all"))
})

test_that("building blocks of the panel tests", {
  ns <- asNamespace("Rcausality")
  # Dumitrescu-Hurlin: moments of K F(K, d2)
  dh <- ns$.dh_test(W = rep(1, 10), K = 1, d2 = rep(20, 10))
  expect_equal(unname(dh["W_bar"]), 1); expect_lt(unname(dh["Z_tilde"]), 0)
  # prevalence model: no evidence -> BF 1; strong evidence in all units -> large BF, prevalence near 1
  p1 <- ns$.prevalence_bf(rep(1, 12)); expect_equal(p1$bf10, 1, tolerance = 1e-8)
  p2 <- ns$.prevalence_bf(rep(50, 12)); expect_gt(p2$bf10, 1e6); expect_gt(p2$mean, 0.8)
  p3 <- ns$.prevalence_bf(c(rep(100, 4), rep(0.2, 8)))
  expect_true(p3$lower < p3$mean && p3$mean < p3$upper); expect_true(all(p3$unit_prob >= 0 & p3$unit_prob <= 1))
  expect_gt(min(p3$unit_prob[1:4]), max(p3$unit_prob[5:12]))
  # Bayesian FDR: the most probable cells with mean posterior null probability <= alpha
  expect_equal(ns$.bayes_fdr(c(0.99, 0.97, 0.6, 0.2), 0.05), c(TRUE, TRUE, FALSE, FALSE))
  expect_false(any(ns$.bayes_fdr(c(0.6, 0.5), 0.05)))
  # combinations of p-values
  expect_lt(ns$.combine_p(rep(0.04, 10), "stouffer")$p_value, 1e-4)
  expect_equal(ns$.combine_p(c(0.01, 0.5, 0.9), "simes")$p_value, 0.03)
  # the group test: exchangeable pools give no evidence, a shifted observation does
  set.seed(1)
  nl0 <- lapply(1:20, function(i) { v <- rnorm(11); list(statistic = v[1], null = v[-1], alternative = "greater") })
  nl1 <- lapply(nl0, function(e) { e$statistic <- e$statistic + 2; e })
  expect_gt(ns$.group_mean_test(nl0, 2000)$p_value, 0.01); expect_lt(ns$.group_mean_test(nl1, 2000)$p_value, 0.01)
  # empirical null: centred null z-scores leave the p-values, a shifted or wider null makes them larger
  set.seed(2); p0 <- pnorm(rnorm(400), lower.tail = FALSE)
  en <- ns$.empirical_null(p0); expect_lt(en[["mu0"]], 0.25); expect_lt(en[["sigma0"]], 1.25)
  en2 <- ns$.empirical_null(pnorm(rnorm(400, 0.8, 1.5), lower.tail = FALSE)); expect_gt(en2[["mu0"]], 0.5); expect_gt(en2[["sigma0"]], 1.2)
  expect_true(all(ns$.p_calibrate(c(0.001, 0.04, 0.5), en2) >= c(0.001, 0.04, 0.5)))
  expect_equal(ns$.p_calibrate(0.02, c(mu0 = 0, sigma0 = 1)), 0.02)
  # the bootstrap effect keeps the initial values and has no link from the cause
  set.seed(3); e <- as.numeric(arima.sim(list(ar = 0.6), 40)); Eb <- ns$.boot_effect(e, matrix(1, 39, 1), 1L, 50L)
  expect_equal(dim(Eb), c(40L, 50L)); expect_true(all(Eb[1, ] == e[1])); expect_true(all(is.finite(Eb)))
})

test_that("the engine reproduces the unit F-tests and the pooled half-panel-jackknife estimate", {
  ns <- asNamespace("Rcausality")
  pd <- sim_panel(n_units = 8, n_time = 24, seed = 4); U <- ns$.panel_units(pd, "unit", "time")$units
  set.seed(1); eng <- ns$.reg_engine(U, "y", "x", character(0), 1:2, 0L, 9L)
  for (l in 1:2) {
    Fu <- vapply(U, function(u) ns$.granger_pair(u$y, u$x, NULL, l)$statistic, 0)
    expect_equal(unname(eng[[l]]$unitF[, 1]), unname(Fu), tolerance = 1e-8)
  }
  # pooled estimate at lag 1 by hand: unit-specific intercept and own lag projected out, cause lag pooled
  est <- function(rows_of) {
    A <- 0; b <- 0
    for (u in U) {
      e <- ns$.std(u$y); x <- ns$.std(u$x); n <- length(e); rows <- rows_of(n - 1)
      Z <- cbind(1, e[1:(n - 1)])[rows, ]; X <- x[1:(n - 1)][rows]; y <- e[2:n][rows]
      xt <- lm.fit(Z, X)$residuals; yt <- lm.fit(Z, y)$residuals; A <- A + sum(xt^2); b <- b + sum(xt * yt)
    }
    b / A
  }
  hpj <- 2 * est(function(m) seq_len(m)) - (est(function(m) seq_len(floor(m / 2))) + est(function(m) (floor(m / 2) + 1):m)) / 2
  expect_equal(eng[[1]]$obs$beta, hpj, tolerance = 1e-8)
  expect_true(eng[[1]]$obs$hpj); expect_equal(eng[[1]]$obs$n_units, 8L)
  expect_true(all(eng[[1]]$p >= 0 & eng[[1]]$p <= 1)); expect_length(eng[[1]]$p, 10L)
  # Bayes factors of the pooled estimate under the priors of bayes_control()
  for (pr in c("bic", "g", "hyper_g", "normal", "spike_slab", "minnesota"))
    expect_gt(ns$.wald_bf(eng[[1]]$obs, bayes_control(prior = pr)), 10)
})

test_that("Granger causality in a panel: units, group and pooled level", {
  pd <- sim_panel(n_units = 25, n_time = 20, scenario = "x_to_y", seed = 1)
  pg <- panel_causality(pd, id = "unit", time = "time", seed = 1)
  expect_s3_class(pg, c("causpanel", "causres"))
  expect_equal(pg$classes$class, "x_to_y")
  expect_setequal(unique(pg$panel$level), c("group", "pooled")); expect_equal(pg$settings$panel$decision, "both")
  P <- pg$panel; expect_true(all(P$significant[P$cause == "x"])); expect_false(any(P$significant[P$cause == "y"]))
  # the decision of both levels: the smaller p-value, doubled
  expect_equal(pg$pairs$p_value, vapply(split(P$p_value, factor(paste(P$cause, P$effect), unique(paste(P$cause, P$effect)))), function(p) min(1, 2 * min(p)), 0),
               ignore_attr = TRUE)
  U <- pg$units; expect_equal(nrow(U), 50L)
  expect_true(all(c("p_ftest", "p_value", "q_value", "nominal", "significant", "post_pooled", "bf10") %in% names(U)))
  expect_true(all(U$p_value >= U$p_ftest - 1e-12))              # the calibration is conservative
  expect_true(all(!U$significant | U$nominal))
  # a unit row is the time-series analysis of that unit at the panel's lag order
  u1 <- pd[pd$unit == "u07", c("x", "y")]; r1 <- rcausality(u1, lag = pg$pairs$lag[pg$pairs$cause == "x"], adjust = "none")
  k <- U[U$unit == "u07" & U$cause == "x", ]
  expect_equal(k$statistic, r1$pairs$statistic[r1$pairs$cause == "x"], tolerance = 1e-8)
  expect_equal(k$p_ftest, r1$pairs$p_value[r1$pairs$cause == "x"], tolerance = 1e-8)
  expect_equal(k$bf10_raw, r1$pairs$bf10[r1$pairs$cause == "x"], tolerance = 1e-6)
  expect_equal(pg$prevalence$n_units, c(25L, 25L)); expect_gt(pg$prevalence$prevalence[pg$prevalence$cause == "x"], 0.3)
  expect_output(print(pg), "Panel causality analysis"); expect_output(print(summary(pg)), "Panel-level tests")
  M <- as.matrix(pg, "p_value", by = "units"); expect_equal(dim(M), c(2L, 25L)); expect_equal(dim(as.matrix(pg)), c(2L, 2L))
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  for (tp in c("heatmap", "units", "prevalence", "network", "matrix")) expect_silent(plot(pg, tp))
  expect_silent(plot(pg, fill = "evidence", cluster = FALSE)); expect_silent(plot(pg, links = "significant"))
  # same seed, same result; the wrappers
  expect_equal(panel_causality(pd, id = "unit", time = "time", seed = 1)$pairs, pg$pairs)
  expect_equal(panel_granger(pd, id = "unit", time = "time", seed = 1)$pairs, pg$pairs)
  expect_equal(rcausality(pd, id = "unit", time = "time", seed = 1)$pairs, pg$pairs)
})

test_that("decision levels, group tests and lag rules", {
  pd <- sim_panel(n_units = 20, n_time = 20, scenario = "x_to_y", seed = 2)
  f <- function(...) panel_causality(pd, id = "unit", time = "time", seed = 1, ...)
  g <- f(panel = panel_control(decision = "group")); p <- f(panel = panel_control(decision = "pooled"))
  expect_equal(g$settings$panel$decision, "group"); expect_equal(p$settings$panel$decision, "pooled")
  expect_equal(g$pairs$p_value, g$panel$p_value[g$panel$level == "group"])
  expect_equal(p$pairs$statistic, p$panel$statistic[p$panel$level == "pooled"])
  expect_equal(g$classes$class, "x_to_y"); expect_equal(p$classes$class, "x_to_y")
  po <- f(panel = panel_control(levels = "pooled")); expect_null(po$units); expect_equal(po$settings$panel$decision, "pooled")
  for (cb in c("dh", "stouffer", "fisher", "simes", "count")) {
    r <- f(panel = panel_control(decision = "group", combine = cb))
    expect_equal(unique(r$panel$test[r$panel$level == "group"]), cb)
    expect_true(r$pairs$significant[r$pairs$cause == "x"])
  }
  expect_error(panel_causality(pd, id = "unit", method = "te", n_perm = 20, panel = panel_control(combine = "dh")), "Granger family")
  r1 <- f(lag = 1); expect_true(all(r1$pairs$lag == 1L)); expect_true(all(r1$units$lag == 1L))
  rf <- f(lag = "first", max_lag = 3); expect_equal(rf$settings$panel$lag_rule, "first"); expect_true(rf$pairs$significant[rf$pairs$cause == "x"])
  ty <- f(method = "toda_yamamoto"); expect_equal(ty$classes$class, "x_to_y")
  cn <- panel_causality(sim_panel(15, 25, "drivers", seed = 1)[, c("unit", "time", "x1", "x2", "y")], id = "unit", time = "time",
                        conditional = TRUE, seed = 1)
  expect_true(cn$pairs$significant[cn$pairs$cause == "x1" & cn$pairs$effect == "y"])
  expect_false(cn$pairs$significant[cn$pairs$cause == "x2" & cn$pairs$effect == "y"])     # x2 acts through x1 only
})

test_that("Bayesian panel variants", {
  pd <- sim_panel(n_units = 20, n_time = 20, scenario = "x_to_y", seed = 2)
  for (pr in c("bic", "hyper_g")) {
    b <- rcausality(pd, id = "unit", time = "time", bayes = TRUE, prior = pr, seed = 1)
    expect_true(b$settings$bayes); expect_true(b$settings$panel$use_bf)
    expect_equal(b$classes$class, "x_to_y")
    expect_true(all(is.na(b$pairs$p_value))); expect_true(all(is.finite(b$pairs$bf10)))
    expect_gt(b$pairs$bf10[b$pairs$cause == "x"], 100); expect_lt(b$pairs$bf10[b$pairs$cause == "y"], 3)
    expect_equal(b$pairs$bf10, vapply(split(b$panel$bf10, factor(paste(b$panel$cause, b$panel$effect), unique(paste(b$panel$cause, b$panel$effect)))), mean, 0),
                 ignore_attr = TRUE)
  }
  expect_output(print(b), "Bayesian variant")
  skip_on_cran()
  # the sampled priors: unit Bayes factors from the sampler, pooled Bayes factor from the pooled estimate
  s <- rcausality(pd, id = "unit", time = "time", method = "bayes_granger", control = bayes_control(n_iter = 200, burnin = 50), seed = 1)
  expect_equal(s$classes$class, "x_to_y"); expect_equal(s$settings$method, "bayes_granger")
  m <- rcausality(pd, id = "unit", time = "time", bayes = TRUE, prior = "minnesota", control = bayes_control(n_iter = 200), seed = 1)
  expect_equal(m$classes$class, "x_to_y")
})

test_that("model-free methods in short panels: pooled estimators", {
  skip_on_cran()
  pn <- sim_panel(n_units = 30, n_time = 20, scenario = "x_to_y", nonlinear = TRUE, coupling = 0.8, heterogeneity = 0.2, seed = 2)
  lin <- panel_causality(pn, id = "unit", time = "time", seed = 1)
  cm <- panel_causality(pn, id = "unit", time = "time", method = "cmi", n_perm = 60, seed = 1)
  expect_equal(cm$settings$panel$decision, "pooled"); expect_null(cm$units)           # units of 20 periods are not analysed one by one
  expect_true(cm$pairs$significant[cm$pairs$cause == "x"]); expect_false(cm$pairs$significant[cm$pairs$cause == "y"])
  # the squared coupling leaves no linear trace: the linear effect size is negligible next to the dependence found by the CMI
  expect_lt(lin$panel$r2_gain[lin$panel$level == "pooled" & lin$panel$cause == "x"], 0.05)
  te <- panel_causality(pn, id = "unit", time = "time", method = "te", n_perm = 60, seed = 1, method_args = list(bins = 3))
  expect_true(all(is.finite(te$pairs$p_value)))
  tb <- panel_causality(pn, id = "unit", time = "time", method = "te", bayes = TRUE, n_perm = 60, seed = 1, method_args = list(bins = 3))
  expect_true(all(is.finite(tb$pairs$bf10)))
  dp <- panel_causality(pn, id = "unit", time = "time", method = "diks_panchenko", n_perm = 40, seed = 1)
  expect_true(all(is.finite(dp$pairs$p_value)))
  cc <- panel_causality(pn, id = "unit", time = "time", method = "ccm", n_perm = 30, seed = 1)
  expect_true(all(is.finite(cc$pairs$statistic)))
  skip_if_not_installed("ranger")
  rf <- panel_causality(pn, id = "unit", time = "time", method = "rf_granger", n_perm = 30, seed = 1, method_args = list(num_trees = 50))
  expect_true(rf$pairs$significant[rf$pairs$cause == "x"])
  # a method without a pooled estimator needs units that are long enough
  expect_error(panel_causality(pn, id = "unit", time = "time", method = "pid", n_perm = 20), "no pooled estimator")
})

test_that("model-free methods with units analysed one by one: surrogate group test", {
  skip_on_cran()
  pl <- sim_panel(n_units = 6, n_time = 80, scenario = "x_to_y", coupling = 0.6, seed = 3)
  te <- panel_causality(pl, id = "unit", time = "time", method = "te", n_perm = 100, seed = 1)
  expect_setequal(unique(te$panel$level), c("group", "pooled")); expect_equal(unique(te$panel$test[te$panel$level == "group"]), "mean")
  expect_equal(nrow(te$units), 12L); expect_true(te$pairs$significant[te$pairs$cause == "x"])
  tb <- panel_causality(pl, id = "unit", time = "time", method = "te", bayes = TRUE, n_perm = 40, seed = 1)
  expect_true(all(is.finite(tb$units$bf10))); expect_true(all(is.finite(tb$prevalence$prevalence_bf10)))
})

test_that("unbalanced panels, missing values, time effects, input checks", {
  pd <- sim_panel(n_units = 12, n_time = c(rep(18, 6), rep(26, 6)), scenario = "x_to_y", coupling = 0.6, seed = 5)
  pd$y[c(5, 40)] <- NA
  r <- panel_causality(pd[sample(nrow(pd)), ], id = "unit", time = "time", seed = 1)      # row order does not matter
  expect_equal(r$settings$panel$n_range, c(18L, 26L)); expect_true(r$pairs$significant[r$pairs$cause == "x"])
  expect_error(panel_causality(pd, id = "unit", time = "time", impute = FALSE), "impute")
  rt <- panel_causality(pd, id = "unit", time = "time", seed = 1, panel = panel_control(time_effects = TRUE))
  expect_s3_class(rt, "causpanel")
  expect_error(panel_causality(pd, id = "unit", panel = panel_control(time_effects = TRUE)), "time")
  expect_error(panel_causality(pd, id = "nope"), "not found")
  expect_error(panel_causality(rbind(pd, pd[1, ]), id = "unit", time = "time"), "Duplicated")
  expect_error(panel_causality(pd[pd$unit == "u01", ], id = "unit", time = "time"), "two units")
  expect_error(rcausality(pd, y = 1:3, id = "unit"), "long data frame")
})

test_that("panel results do not depend on the number of cores", {
  skip_on_cran(); skip_if_not(.Platform$OS.type != "windows" && parallel::detectCores() >= 2)
  pd <- sim_panel(n_units = 10, n_time = 20, scenario = "drivers", seed = 1)[, c("unit", "time", "x1", "x3", "y")]
  a <- panel_causality(pd, id = "unit", time = "time", seed = 3, cores = 1); b <- panel_causality(pd, id = "unit", time = "time", seed = 3, cores = 2)
  expect_equal(a$pairs, b$pairs); expect_equal(a$units, b$units); expect_equal(a$panel, b$panel)
})

test_that("placebo test, robustness, comparison and forecasts of a panel analysis", {
  skip_on_cran()
  pd <- sim_panel(n_units = 15, n_time = 24, scenario = "x_to_y", coupling = 0.6, seed = 1)
  pg <- panel_causality(pd, id = "unit", time = "time", seed = 1)
  pl <- placebo_test(pg, n_placebo = 5, reverse = TRUE, seed = 1)
  expect_s3_class(pl, "causplacebo"); expect_true(pl$settings$panel); expect_length(pl$counts, 5L)
  expect_lte(mean(pl$counts), pl$observed); expect_output(print(pl), "placebo panels")
  # a placebo panel keeps every series and gives each unit its variables from different units
  set.seed(1); Dp <- Rcausality:::.placebo_panel(pg, pg$data, 0.1)
  expect_equal(sort(Dp$x), sort(pg$data$x)); expect_equal(Dp$unit, pg$data$unit)
  sx <- vapply(split(Dp$x, Dp$unit), function(v) names(which(vapply(split(pg$data$x, pg$data$unit), function(w) isTRUE(all.equal(w, v)), TRUE)))[1], "")
  sy <- vapply(split(Dp$y, Dp$unit), function(v) names(which(vapply(split(pg$data$y, pg$data$unit), function(w) isTRUE(all.equal(w, v)), TRUE)))[1], "")
  expect_false(any(sx == sy))
  rb <- robustness(pg, lags = 1:2, n_sub = 3, seed = 1)
  expect_true(all(c("lag", "subsample", "level", "adjust") %in% rb$grid$type))
  expect_gt(rb$summary$share[rb$summary$link == "x -> y"], 0.7); expect_output(print(rb), "Robustness")
  cmp <- compare_methods(pd, methods = c("granger", "toda_yamamoto"), id = "unit", time = "time", seed = 1)
  expect_equal(cmp$table$consensus, "x_to_y")
  expect_warning(c2 <- compare_methods(pd, methods = c("granger", "pid"), id = "unit", time = "time", seed = 1, n_perm = 20), "left out")
  expect_equal(c2$methods, "granger")
  fc <- predict(pg, h = 3)
  expect_s3_class(fc, "causfcpanel"); expect_equal(nrow(fc$forecasts), 15L * 2L * 3L)
  expect_true(all(fc$forecasts$lower <= fc$forecasts$point & fc$forecasts$point <= fc$forecasts$upper))
  expect_gt(fc$evaluation$improvement[fc$evaluation$series == "y"], 0)
  expect_output(print(fc), "panel forecasts"); expect_output(print(summary(fc)), "improvement")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off()); expect_silent(plot(fc)); expect_silent(plot(pl)); expect_silent(plot(rb))
  fb <- predict(pg, h = 2, model = "bvar", n_sim = 30); expect_equal(nrow(fb$forecasts), 15L * 2L * 2L)
  f2 <- forecast_causal(pd, id = "unit", time = "time", h = 2, seed = 1); expect_s3_class(f2, "causfcpanel")
  fcm <- forecast_compare(pd, methods = c("granger", "toda_yamamoto"), id = "unit", time = "time", seed = 1)
  expect_s3_class(fcm, "causfcomp"); expect_equal(fcm$n_units, 15L)
  expect_gt(fcm$table$gain[fcm$table$method == "granger" & fcm$table$series == "y"], 0); expect_output(print(fcm), "units")
})

test_that("benchmark of methods on simulated panels", {
  skip_on_cran()
  b <- benchmark_methods(c("granger", "toda_yamamoto"), scenarios = c("x_to_y", "independent"),
                         simulator = function(n, scenario, seed) sim_panel(15, n, scenario, coupling = 0.6, seed = seed),
                         n = 20, reps = 2, seed = 1, id = "unit", time = "time")
  expect_s3_class(b, "causbench"); expect_equal(dim(b$accuracy), c(2L, 2L)); expect_true(all(is.na(b$runs$error)))
  expect_gte(b$accuracy["granger", "x_to_y"], 0.5)
})
