# Changes of version 0.0.12 (third round of the review of the accompanying paper)

ar2_pair <- function(n, seed) {                                    # x, a second-order autoregression, drives y
  set.seed(seed); x <- as.numeric(arima.sim(list(ar = c(0.4, 0.4)), n)); y <- numeric(n)
  for (t in 2:n) y[t] <- 0.5 * y[t - 1] + 0.8 * x[t - 1] + rnorm(1)
  data.frame(x = x, y = y)
}
flat <- function(x) gsub("\\s+", " ", paste(utils::capture.output(print(x)), collapse = " "))

test_that("the regressions hold the lags of the effect up to its own order", {
  d <- ar2_pair(400, 1)
  expect_gte(Rcausality:::.own_order(d$x), 2L)
  expect_equal(Rcausality:::.own_order(rep(1, 50)), 0L)
  set.seed(2); expect_lte(Rcausality:::.own_order(rnorm(300)), 3L)
  # own lags beyond the order of the test: the same numerator, one coefficient more
  g1 <- Rcausality:::.granger_pair(d$x, d$y, NULL, 1L)
  g2 <- Rcausality:::.granger_pair(d$x, d$y, NULL, 1L, own = 2L)
  expect_equal(g1$df1, 1); expect_equal(g2$df1, 1); expect_equal(g2$df2, g1$df2 - 2); expect_equal(g2$own, 2L)
  n <- nrow(d); y <- d$x[3:n]; o1 <- d$x[2:(n - 1)]; o2 <- d$x[1:(n - 2)]; c1 <- d$y[2:(n - 1)]
  expect_equal(g2$p_value, stats::anova(stats::lm(y ~ o1 + o2), stats::lm(y ~ o1 + o2 + c1))[["Pr(>F)"]][2], tolerance = 1e-9)
  # an order chosen by the package: own lags of the effect; a numeric order: the regression of that order
  r <- suppressWarnings(rcausality(d, seed = 1))
  expect_true(all(r$extra$own_lags >= r$pairs$lag)); expect_gte(r$extra$own_lags[r$extra$effect == "x"], 2L)
  r1 <- suppressWarnings(rcausality(d, lag = 1, seed = 1))
  expect_true(all(r1$extra$own_lags == 1L)); expect_equal(r1$pairs$p_value[r1$pairs$cause == "y"], g1$p_value)
  ty <- suppressWarnings(rcausality(d, method = "toda_yamamoto", seed = 1))
  expect_true(all(ty$extra$own_lags >= ty$pairs$lag + 1L))
  # the Bayesian regressions take the same own lags
  b1 <- bayes_granger(d$x, d$y, lag = 1L, control = bayes_control(prior = "g"))
  b2 <- bayes_granger(d$x, d$y, lag = 1L, control = bayes_control(prior = "g"), own = 2L)
  expect_gt(b1$bf10, 3); expect_lt(b2$bf10, 1)                    # the first-order regression "finds" the reverse link
})

test_that("the reverse direction of a link with a second-order driver is not declared", {
  skip_on_cran()
  rej <- function(lag) mean(vapply(1:60, function(s) { d <- ar2_pair(200, 100 + s)
    suppressWarnings(suppressMessages(rcausality(data.frame(c = d$y), data.frame(e = d$x), directions = "forward",
                                                 lag = lag)))$pairs$p_value[1] < 0.05 }, TRUE))
  expect_lte(rej("auto"), 0.15); expect_gte(rej(1L), 0.5)
  rej_dp <- function(own) mean(vapply(1:40, function(s) { d <- ar2_pair(200, 200 + s); set.seed(s)
    Rcausality:::.dp_stat(d$x, d$y, NULL, 1L, n_perm = 39L, own = own)$p_value <= 0.05 }, TRUE))
  expect_lte(rej_dp(NULL), 0.2); expect_gte(rej_dp(0L), 0.5)
})

test_that("the conditional surrogates depend on the past of the effect up to its own order", {
  d <- ar2_pair(300, 3)
  set.seed(1); a <- Rcausality:::.dp_stat(d$x, d$y, NULL, 1L, n_perm = 39L)
  set.seed(1); b <- Rcausality:::.dp_stat(d$x, d$y, NULL, 1L, n_perm = 39L, own = 0L)
  expect_gte(a$own, 2L); expect_equal(b$own, 1L); expect_equal(a$statistic, b$statistic)   # the statistic is the same
  expect_equal(ncol(Rcausality:::.cond_extra(matrix(0, 5, 1), 1:10, 6:10, 1L, 3L)), 3L)
  expect_equal(ncol(Rcausality:::.cond_extra(matrix(0, 5, 1), 1:10, 6:10, 1L, 1L)), 1L)
  expect_equal(unname(Rcausality:::.hist_pad(1:10, 2:4, 2, 3)[1, ]), c(1, 1))               # before the start: the first value
  # the simpler estimators of the earlier versions are kept for the comparisons of the paper
  op <- options(Rcausality.cond_split = "linear"); on.exit(options(op), add = TRUE)
  W <- cbind(d$x[-nrow(d)]); X <- cbind(d$y[-nrow(d)]); sp <- Rcausality:::.cond_split(W, X)
  expect_equal(length(unique(round(sp$scale[, 1], 10))), 1L)
  expect_equal(sp$fit[, 1], unname(stats::fitted(stats::lm(X[, 1] ~ W[, 1]))), tolerance = 1e-8)
})

test_that("conditional shifts are the default reference of all tests of conditional hypotheses, under one name", {
  for (m in c("te", "pid", "cmi", "cmi_lag", "df", "diks_panchenko"))
    expect_identical(Rcausality:::.method_arg_defaults[[m]]$reference, "conditional")
  d <- sim_var(300, "x_to_y", seed = 1)
  r <- suppressMessages(rcausality(d, method = "te", n_perm = 60, seed = 1))
  expect_true(all(r$extra$reference == "conditional shifts"))
  expect_identical(eval(formals(transfer_entropy)$reference)[1], "conditional")
  rc <- suppressMessages(rcausality(d, method = "cmi", n_perm = 40, seed = 1, method_args = list(reference = "conditional")))
  expect_true(all(rc$extra$reference == "conditional shifts"))
  expect_output(print(causal_methods()), "min.length")
})

test_that("the directed test of df uses conditional shifts, and its Bayesian decision keeps the table of the orders", {
  skip_on_cran()
  d <- sim_var(300, "x_to_y", seed = 1)
  r <- suppressWarnings(suppressMessages(rcausality(d, method = "df", n_perm = 40, seed = 1)))
  expect_true(all(is.finite(r$pairs$p_value))); expect_true(all(r$pairs$p_value >= 1 / 41))
  l <- suppressWarnings(suppressMessages(rcausality(d, method = "df", n_perm = 40, seed = 1, method_args = list(reference = "local"))))
  expect_true(all(is.finite(l$pairs$p_value)))
  # statistics far beyond their replicates in both directions: the table of the orders decides what is declared
  big <- r; for (k in names(big$null)) big$null[[k]]$statistic <- max(big$null[[k]]$null) + 10 * stats::sd(big$null[[k]]$null)
  pp <- function(t2) { z <- big; z$df$table2 <- t2
    suppressWarnings(Rcausality:::.bayes_postprocess(z, bayes_control(), 0.5, 3, 0.5))$pairs }
  a <- pp("independent"); expect_true(all(a$bf10 > 3)); expect_false(any(a$significant))
  u <- pp("undetermined"); expect_true(all(u$significant))
  v1 <- big$df$var1[1]; o <- pp("x_to_y"); expect_true(o$significant[o$cause == v1]); expect_false(o$significant[o$effect == v1])
})

test_that("with jointly resampled periods the group statistic is referred to its values in the bootstrap panels", {
  set.seed(1); pool <- stats::rexp(200); obs <- unname(stats::quantile(pool, 0.9))
  nulls <- rep(list(list(statistic = obs, null = pool, alternative = "greater")), 20)   # 20 units that move as one
  a <- Rcausality:::.group_mean_test(nulls, aligned = TRUE); u <- Rcausality:::.group_mean_test(nulls, 2000L)
  expect_true(isTRUE(a$aligned)); expect_gt(a$p_value, 0.05); expect_lt(u$p_value, 0.01)
  pd <- sim_panel(n_units = 12, n_time = 20, scenario = "independent", seed = 1)
  r <- suppressWarnings(suppressMessages(panel_causality(pd, id = "unit", time = "time", seed = 1, panel = panel_control(resample = "periods"))))
  g <- r$panel[r$panel$level == "group", ]
  expect_true(all(grepl("199 bootstrap panels with jointly resampled periods", g$detail)))
  expect_true(all(g$p_freq >= 1 / 200 - 1e-12))
  r0 <- suppressWarnings(suppressMessages(panel_causality(pd, id = "unit", time = "time", seed = 1)))
  expect_true(all(grepl("against their bootstrap null", r0$panel$detail[r0$panel$level == "group"])))
  # more links than 199 bootstrap panels can resolve after the correction: more panels
  pd$z <- stats::rnorm(nrow(pd))
  expect_message(suppressWarnings(panel_causality(pd, id = "unit", time = "time", seed = 1, panel = panel_control(resample = "periods"))),
                 "240 bootstrap panels instead of 199")
})

test_that("panel warnings: persistent units, the resolution of the unit level, an unused covariance", {
  chk <- Rcausality:::.panel_persistence_check; set.seed(3)
  rw <- lapply(1:10, function(i) data.frame(x = cumsum(1 + stats::rnorm(20)), y = cumsum(1 + stats::rnorm(20))))
  grid <- data.frame(cause = c("x", "y"), effect = c("y", "x"), stringsAsFactors = FALSE)
  expect_warning(chk(rw, grid, "granger"), "tests reject then too often")
  expect_warning(chk(rw, grid, "granger", bayes = TRUE), "analysis declares links then too often")
  expect_warning(chk(rw, grid, "toda_yamamoto"), "lag-augmented panel")
  set.seed(4); st <- lapply(1:10, function(i) data.frame(x = stats::rnorm(20), y = stats::rnorm(20)))
  expect_no_warning(chk(st, grid, "granger")); expect_no_warning(chk(st, grid, "toda_yamamoto"))
  # unit tests with shifted series: the correction over the unit x link cells asks for more than they can give
  pd <- sim_panel(n_units = 12, n_time = 60, seed = 1); w <- character(0)
  withCallingHandlers(suppressMessages(panel_causality(pd, id = "unit", time = "time", method = "te", seed = 1)),
                      warning = function(e) { w <<- c(w, conditionMessage(e)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("unit level \\(te\\).*attain that value together", w)))
  # the covariance options have no panel version
  ps <- sim_panel(n_units = 8, n_time = 20, seed = 1); w <- character(0)
  withCallingHandlers(suppressMessages(panel_causality(ps, id = "unit", time = "time", seed = 1, method_args = list(vcov = "HC3"))),
                      warning = function(e) { w <<- c(w, conditionMessage(e)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("vcov = .HC3.. is not used", w)))
})

test_that("rows of tests with a constant series carry no evidence", {
  d <- sim_var(100, seed = 1); d$y <- 3
  expect_warning(r <- rcausality(d, seed = 1), "p = 1 and BF10 = 1")
  expect_true(all(r$pairs$p_value == 1)); expect_true(all(r$pairs$bf10 == 1)); expect_false(any(r$pairs$significant))
})

test_that("gaps in a daily date index of a panel are reported", {
  days <- as.Date("2024-01-01") + 0:69; wd <- days[!format(days, "%u") %in% c("6", "7")]
  mk <- function(dates) { set.seed(1)
    do.call(rbind, lapply(1:6, function(i) data.frame(unit = i, time = dates, x = stats::rnorm(length(dates)), y = stats::rnorm(length(dates))))) }
  warns <- function(p) { w <- character(0)
    withCallingHandlers(suppressMessages(panel_causality(p, id = "unit", time = "time", seed = 1)),
                        warning = function(e) { w <<- c(w, conditionMessage(e)); invokeRestart("muffleWarning") })
    w }
  expect_false(any(grepl("gaps in", warns(mk(wd)))))             # weekends are no gaps
  expect_true(any(grepl("gaps in", warns(mk(wd[-(15:24)])))))    # two missing weeks are
})

test_that("the output names the decision rule, the prior probability and the scale of the effect prior", {
  d <- sim_var(300, "x_to_y", seed = 1)
  expect_warning(r <- suppressMessages(rcausality(d, method = "te", bayes = TRUE, seed = 1, prior_prob = 1 - 0.5^(1 / 2),
                                                  control = bayes_control(decision = "posterior"))),
                 "beyond the range that the replicates support")
  out <- flat(r)
  expect_match(out, "posterior probability > 0.5 with the prior probability 0.293 of a link (BF10 > 2.41)", fixed = TRUE)
  expect_match(out, "effect_scale = \"local\"", fixed = TRUE)
  s <- flat(summary(r)); expect_match(s, "decision: posterior probability > 0.5", fixed = TRUE); expect_false(grepl("alpha =", s, fixed = TRUE))
  u <- suppressWarnings(suppressMessages(rcausality(d, method = "te", bayes = TRUE, seed = 1, n_perm = 60, control = bayes_control(effect_scale = "unit"))))
  expect_match(flat(u), "effect_scale = \"unit\"", fixed = TRUE)
  g <- flat(suppressWarnings(rcausality(d, seed = 1))); expect_match(g, "White test of the regression errors", fixed = TRUE)
})

test_that("the placebo test reports the share of placebo data sets with a declared link", {
  d <- sim_var(200, n_series = 3, seed = 2); r <- suppressWarnings(rcausality(d))
  pl <- suppressWarnings(suppressMessages(placebo_test(r, n_placebo = 19, seed = 1)))
  expect_equal(pl$any_rate, mean(pl$counts > 0)); expect_true(pl$any_rate >= 0 && pl$any_rate <= 1)
  expect_output(print(pl), "placebo data sets with at least one link")
})

test_that("mutually coupled maps are not returned as identical series without notice", {
  for (s in 1:5) {
    h <- sim_henon(500, "bidirectional", seed = s); expect_gt(max(abs(h$x - h$y)), 0.1)
    l <- sim_logistic(500, "bidirectional", seed = s); expect_gt(max(abs(l$x - l$y)), 0.1)
  }
  expect_warning(sim_henon(500, "bidirectional", coupling = 0.3, seed = 1), "synchronised")
  expect_warning(sim_logistic(500, "bidirectional", coupling = 0.3, seed = 1), "synchronised")
  for (s in c(1, 5, 7, 15, 25)) { b <- sim_logistic_bdc(1000, "bidirectional", seed = s); expect_gt(max(abs(b$x - b$y)), 1e-3) }
})

test_that("a step-up correction limits the scan of the lags as a Holm correction does", {
  set.seed(1); d <- as.data.frame(matrix(stats::rnorm(5 * 300), 300, 5))
  expect_message(r <- suppressWarnings(rcausality(d, method = "te", adjust = "BH", seed = 1)), "only lag 1 is tested")
  expect_true(all(r$pairs$lag == 1))
})

test_that("the package never starts more processes than `cores`", {
  ns <- asNamespace("Rcausality")
  # a call made inside a worker runs serially
  old <- options(Rcausality.in_worker = TRUE)
  expect_identical(ns$.plapply(1:3, function(i) Sys.getpid(), cores = 2), as.list(rep(Sys.getpid(), 3)))
  options(old)
  skip_on_os("windows"); skip_if_not(parallel::detectCores() >= 2)
  pids <- unlist(ns$.plapply(1:4, function(i) unlist(ns$.plapply(1:2, function(j) Sys.getpid(), cores = 2)), cores = 2))
  expect_false(Sys.getpid() %in% pids); expect_lte(length(unique(pids)), 2L)
  # the results do not depend on the nesting
  set.seed(3); a <- ns$.plapply(1:4, function(i) unlist(ns$.plapply(1:2, function(j) stats::runif(1), cores = 2)), cores = 2)
  set.seed(3); b <- ns$.plapply(1:4, function(i) unlist(ns$.plapply(1:2, function(j) stats::runif(1), cores = 1)), cores = 1)
  expect_identical(a, b)
})

test_that("cores = \"auto\" respects the limits of the environment", {
  ns <- asNamespace("Rcausality"); keys <- c("_R_CHECK_LIMIT_CORES_", "OMP_THREAD_LIMIT")
  oe <- Sys.getenv(keys, unset = NA)
  on.exit(for (k in keys) if (is.na(oe[[k]])) Sys.unsetenv(k) else do.call(Sys.setenv, as.list(oe[k])), add = TRUE)
  Sys.unsetenv("OMP_THREAD_LIMIT"); Sys.setenv("_R_CHECK_LIMIT_CORES_" = "TRUE")
  expect_lte(ns$.n_cores("auto"), 2L)
  Sys.setenv("_R_CHECK_LIMIT_CORES_" = "false", OMP_THREAD_LIMIT = "1")
  expect_identical(ns$.n_cores("auto"), 1L)
  expect_identical(ns$.n_cores(2), 2L)                # an explicit number is used as it is
})

test_that("the forests of \"dfe\" are given their number of threads", {
  src <- paste(deparse(Rcausality::dfe_train), collapse = " ")
  expect_match(src, "num.threads")
  expect_false(grepl("predict\\(rf, data = Fm\\)", src))
})

test_that("the option for the contrasts of the paper switches the own-order rule off", {
  set.seed(5); x <- as.numeric(arima.sim(list(ar = c(0.4, 0.4)), 300))
  expect_gte(Rcausality:::.own_order(x), 2L)
  old <- options(Rcausality.own_memory = FALSE); on.exit(options(old), add = TRUE)
  expect_identical(Rcausality:::.own_order(x), 0L)
  y <- numeric(300); for (t in 2:300) y[t] <- 0.5 * y[t - 1] + 0.8 * x[t - 1] + rnorm(1)
  r <- suppressWarnings(rcausality(data.frame(x = x, y = y)))
  expect_true(all(r$extra$own_lags == r$pairs$lag))
})


test_that("df: deterministic pairs are tested with local permutations, and a negative reduction is no evidence", {
  skip_on_cran()
  d <- sim_df_rotation(1500, "deterministic", seed = 1)
  expect_message(r <- rcausality(d[, c("z1", "x1")], method = "df", lag = 1, max_constr = 4, n_perm = 39, seed = 1),
                 "local permutations")
  expect_identical(r$classes$relation, "z1 -> x1")
  expect_match(r$extra$reference[r$extra$cause == "z1"], "determine")
  i <- r$pairs$cause == "x1"
  if (r$pairs$statistic[i] <= 0) expect_identical(r$pairs$p_value[i], 1)
  set.seed(2); a <- as.numeric(arima.sim(list(ar = 0.9), 300)); b <- as.numeric(arima.sim(list(ar = 0.9), 300))
  r2 <- suppressWarnings(suppressMessages(rcausality(data.frame(a = a, b = b), method = "df", seed = 1)))
  expect_true(all(r2$extra$reference == "conditional shifts"))
  r3 <- suppressWarnings(suppressMessages(rcausality(data.frame(a = a, b = b), method = "df", seed = 1,
                                                     method_args = list(reference = "local"))))
  expect_true(all(r3$extra$reference == "local permutations"))
})

test_that("a seeded call gives the same result wherever it is made", {
  d <- sim_var(200, "x_to_y", seed = 3)
  direct <- suppressWarnings(rcausality(d, method = "cmi", seed = 1, lag = 1, n_perm = 40))
  nested <- suppressWarnings(compare_methods(d, methods = c("granger", "cmi"), seed = 1, lag = 1, n_perm = 40))
  expect_identical(nested$results$cmi$pairs$p_value, direct$pairs$p_value)
  # data that are generated in the call are the data of the simulator and the seed
  inline <- suppressWarnings(compare_methods(sim_var(200, "x_to_y", seed = 3), methods = c("granger", "cmi"), seed = 1, lag = 1,
                                             n_perm = 40))
  expect_identical(inline$results$granger$pairs$statistic, nested$results$granger$pairs$statistic)
  expect_identical(inline$results$cmi$pairs$p_value, direct$pairs$p_value)
  # a simulator that is called inside a unit of work
  ns <- asNamespace("Rcausality")
  inside <- ns$.plapply(1:2, function(i) sim_var(50, "x_to_y", seed = 7)$x, cores = 1)
  expect_identical(inside[[1]], sim_var(50, "x_to_y", seed = 7)$x); expect_identical(inside[[1]], inside[[2]])
  # the benchmark: the data and the analysis of a replication can be repeated by hand
  bm <- suppressWarnings(benchmark_methods(methods = "cmi", scenarios = "x_to_y", simulator = sim_var, n = 200, reps = 1, seed = 1,
                                           args = list(cmi = list(lag = 1, n_perm = 40))))
  d1 <- sim_var(200, "x_to_y", seed = 1 + 1000 + 1)[, 1:2]
  r1 <- suppressWarnings(rcausality(d1, method = "cmi", lag = 1, n_perm = 40,
                                    seed = 1 + 1000 + 1 + 100000 * match("cmi", ns$.method_names)))
  expect_identical(bm$runs$predicted, r1$classes$class)
  # and the generator of the session is left as it was
  set.seed(11); a <- runif(1); set.seed(11); invisible(compare_methods(d, methods = "granger", seed = 5)); b <- runif(1)
  expect_identical(a, b); expect_identical(RNGkind()[1], "Mersenne-Twister")
})
