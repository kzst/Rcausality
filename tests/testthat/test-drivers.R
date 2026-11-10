# Explanatory causal analysis of one explained variable

test_that("sim_drivers() and the peak lag", {
  d <- sim_drivers(300, seed = 1)
  expect_equal(dim(d), c(300L, 5L)); expect_named(d, c("x1", "x2", "x3", "x4", "y"))
  expect_equal(sim_drivers(50, seed = 2), sim_drivers(50, seed = 2))
  ns <- asNamespace("Rcausality")
  expect_equal(ns$.peak_lag(list(d), "y", "x3", character(0), 6), 2L)       # x3 acts at lag 2
  expect_equal(ns$.peak_lag(list(d), "y", "x1", character(0), 6), 1L)
  expect_equal(ns$.peak_lag(list(d), "x1", "x2", character(0), 6), 1L)
})

test_that("causal drivers of a target: direct, indirect, not a driver", {
  d <- sim_drivers(300, seed = 1)
  cd <- causal_drivers(d, target = "y")
  expect_s3_class(cd, c("causdrivers", "causres"))
  D <- cd$drivers; rownames(D) <- D$driver
  expect_equal(D[c("x1", "x2", "x3", "x4"), "significant"], c(TRUE, TRUE, TRUE, FALSE))
  expect_equal(D[c("x1", "x3"), "role"], c("direct driver", "direct driver"))
  expect_equal(D["x2", "role"], "indirect driver"); expect_equal(D["x2", "via"], "x1"); expect_false(D["x2", "direct"])
  expect_equal(D["x4", "role"], "not a driver")
  expect_equal(D[c("x1", "x3", "x2"), "lag"], c(1L, 2L, 2L))
  expect_equal(sum(D$importance), 1); expect_false(any(D$feedback))
  expect_true(any(cd$graph$cause == "x2" & cd$graph$effect == "x1"))
  # the correction is applied within the three families of links
  P <- cd$pairs; fam <- Rcausality:::.family_of(P$cause, P$effect, "y")
  for (f in unique(fam)) expect_equal(P$p_adj[fam == f], p.adjust(P$p_value[fam == f], "holm"))
  expect_equal(as.vector(table(fam)[c("to_target", "feedback", "among")]), c(4L, 4L, 12L))
  expect_output(print(cd), "Causal drivers of 'y'"); expect_output(print(summary(cd)), "Links among")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  for (tp in c("ladder", "importance", "graph", "matrix")) expect_silent(plot(cd, tp))
  expect_silent(plot(cd, graph = FALSE)); expect_error(plot(cd, "heatmap"), "panel")
  # the same through rcausality(), a subset of drivers, and y given separately
  r <- rcausality(d, target = "y"); expect_equal(r$drivers, cd$drivers)
  c2 <- causal_drivers(d, target = "y", drivers = c("x1", "x4")); expect_setequal(c2$drivers$driver, c("x1", "x4"))
  c3 <- rcausality(d[, 1:4], d$y, target = "y"); expect_equal(c3$drivers$significant, cd$drivers$significant)
  expect_error(causal_drivers(d, target = "z"), "not found|not a numeric"); expect_error(causal_drivers(d[, "y", drop = FALSE], target = "y"), "explanatory")
})

test_that("drivers with other methods and Bayesian variants", {
  d <- sim_drivers(300, seed = 1)
  b <- causal_drivers(d, target = "y", bayes = TRUE, prior = "bic")
  expect_true(b$settings$use_bf); expect_setequal(b$drivers$driver[b$drivers$significant], c("x1", "x2", "x3"))
  expect_equal(b$drivers$role[b$drivers$driver == "x2"], "indirect driver")
  cc <- causal_drivers(d, target = "y", conditional = TRUE)                  # all links conditional: direct by definition
  expect_setequal(cc$drivers$driver[cc$drivers$significant], c("x1", "x3"))
  skip_on_cran()
  te <- causal_drivers(d, target = "y", method = "te", n_perm = 40, seed = 1)
  expect_true(all(is.na(te$drivers$direct))); expect_match(te$settings$drivers_note, "cannot condition")
  expect_true(te$drivers$significant[te$drivers$driver == "x1"]); expect_output(print(te), "Note")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off()); expect_silent(plot(te))
  # robustness tools on a target analysis
  cd <- causal_drivers(d, target = "y")
  pl <- placebo_test(cd, n_placebo = 5, seed = 1); expect_equal(pl$observed, sum(cd$pairs$significant))
  rb <- robustness(cd, lags = 1:2, n_sub = 2, seed = 1); expect_gt(rb$summary$share[rb$summary$link == "x1 -> y"], 0.8)
})

test_that("drivers of a target in a panel", {
  skip_on_cran()
  pd <- sim_panel(n_units = 20, n_time = 30, scenario = "drivers", coupling = 0.5, seed = 1)
  cp <- causal_drivers(pd, target = "y", id = "unit", time = "time", seed = 1)
  expect_s3_class(cp, c("causdrivers", "causpanel", "causres"))
  D <- cp$drivers; rownames(D) <- D$driver
  expect_true(all(D[c("x1", "x3"), "significant"])); expect_false(D["x4", "significant"])
  expect_equal(D["x3", "lag"], 2L); expect_true(all(D[c("x1", "x3"), "direct"]))
  expect_true(any(cp$graph$cause == "x2" & cp$graph$effect == "x1"))
  expect_output(print(cp), "panel of 20 units")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  for (tp in c("ladder", "importance", "heatmap", "units", "prevalence")) expect_silent(plot(cp, tp))
})
