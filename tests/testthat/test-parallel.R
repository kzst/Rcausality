# Parallel execution: results must not depend on the number of workers, and
# the session's generator must be left as a serial call leaves it.

two_cores_ok <- function() {
  if (.Platform$OS.type == "windows") return(FALSE)        # forked workers only in these tests (PSOCK start-up is slow)
  parallel::detectCores() >= 2
}

test_that("serial and parallel runs give identical results (model-free tests)", {
  skip_if_not(two_cores_ok())
  d <- sim_logistic(200, "x_to_y", seed = 1)
  for (m in c("cmi", "te", "ccm", "pid", "ccdh", "hirata")) {
    a <- rcausality(d, method = m, n_perm = 40, n_surr = 40, lag = 1, seed = 1, cores = 1)
    b <- rcausality(d, method = m, n_perm = 40, n_surr = 40, lag = 1, seed = 1, cores = 2)
    expect_equal(a$pairs, b$pairs, info = m)
    if (!is.null(a$instant)) expect_equal(a$instant, b$instant, info = m)
  }
})

test_that("pair-level parallelism (three series) is reproducible", {
  skip_if_not(two_cores_ok())
  d <- sim_var(200, n_series = 3, seed = 2)
  a <- rcausality(d, method = "te", lag = 1, n_perm = 40, seed = 3, cores = 1)
  b <- rcausality(d, method = "te", lag = 1, n_perm = 40, seed = 3, cores = 2)
  expect_equal(a$pairs, b$pairs)
  a <- rcausality(d, bayes = TRUE, control = bayes_control(n_iter = 100, chains = 2), seed = 4, cores = 1)
  b <- rcausality(d, bayes = TRUE, control = bayes_control(n_iter = 100, chains = 2), seed = 4, cores = 2)
  expect_equal(a$bayes, b$bayes)
})

test_that("bayes_granger chains run in parallel with identical draws; R-hat reported", {
  skip_if_not(two_cores_ok())
  d <- sim_var(200, "x_to_y", seed = 1)
  ctl <- bayes_control(n_iter = 150, chains = 3)
  a <- bayes_granger(d$y, d$x, lag = 2, control = ctl, seed = 1, cores = 1)
  b <- bayes_granger(d$y, d$x, lag = 2, control = ctl, seed = 1, cores = 2)
  expect_equal(a$draws, b$draws)
  expect_equal(nrow(a$draws), 3 * 150)
  expect_true(is.finite(a$rhat) && a$rhat > 0.9 && a$rhat < 1.3)
  expect_equal(a$settings$sampler, "exact")
  # the Gibbs sampler (more than 12 lags), thinning respected
  g <- bayes_granger(d$y, d$x, lag = 13, control = bayes_control(n_iter = 40, burnin = 10, thin = 2, chains = 2), seed = 2)
  expect_equal(g$settings$sampler, "Gibbs")
  expect_equal(nrow(g$draws), 2 * 20)
  expect_true(is.finite(g$bf10))
  # single chain: no R-hat
  s1 <- bayes_granger(d$y, d$x, lag = 2, control = bayes_control(n_iter = 50, chains = 1), seed = 1)
  expect_true(is.na(s1$rhat))
})

test_that("standalone tests and the benchmark are worker-count invariant", {
  skip_if_not(two_cores_ok())
  d <- sim_logistic(200, "x_to_y", seed = 1)
  set.seed(1); a <- local_permutation_test(d$x[-1], d$y[-200], n_perm = 30, cores = 1)
  set.seed(1); b <- local_permutation_test(d$x[-1], d$y[-200], n_perm = 30, cores = 2)
  expect_equal(a, b)
  set.seed(2); a <- hirata_rp(d$x, d$y, n_surr = 30, cores = 1)
  set.seed(2); b <- hirata_rp(d$x, d$y, n_surr = 30, cores = 2)
  expect_equal(a$surr, b$surr)
  a <- benchmark_methods(c("granger", "te"), scenarios = c("x_to_y", "independent"), n = 150, reps = 2, seed = 1, cores = 1)
  b <- benchmark_methods(c("granger", "te"), scenarios = c("x_to_y", "independent"), n = 150, reps = 2, seed = 1, cores = 2)
  expect_equal(a$runs[, c("scenario", "rep", "method", "predicted")], b$runs[, c("scenario", "rep", "method", "predicted")])
})

test_that("the session's generator is left where a serial call leaves it; kind unchanged", {
  d <- sim_logistic(150, "x_to_y", seed = 1)
  kind <- RNGkind()[1]
  set.seed(7); invisible(rcausality(d, method = "te", lag = 1, n_perm = 20, cores = 1)); r1 <- runif(1)
  expect_equal(RNGkind()[1], kind)
  # a call without a seed advances the caller's generator by the same amount
  # serially and in parallel (one draw per dispatch of a batch of work)
  if (two_cores_ok()) {
    set.seed(7); invisible(rcausality(d, method = "te", lag = 1, n_perm = 20, cores = 2)); r2 <- runif(1)
    expect_equal(r1, r2)
    expect_equal(RNGkind()[1], kind)
  }
  # with a seed the caller's state is untouched
  set.seed(11); u0 <- runif(1); set.seed(11)
  invisible(rcausality(d, method = "te", lag = 1, n_perm = 20, seed = 3)); expect_equal(runif(1), u0)
  expect_equal(RNGkind()[1], kind)
})

test_that("a user-supplied cluster is accepted and not stopped", {
  skip_on_cran()
  skip_if_not(parallel::detectCores() >= 2)
  d <- sim_logistic(150, "x_to_y", seed = 1)
  cl <- parallel::makeCluster(2)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  for (m in c("te", "cmi", "ccdh", "hirata")) {
    a <- rcausality(d, method = m, lag = 1, n_perm = 30, n_surr = 30, seed = 1, cores = 1)
    b <- rcausality(d, method = m, lag = 1, n_perm = 30, n_surr = 30, seed = 1, cores = cl)
    expect_equal(a$pairs, b$pairs, info = m)
  }
  # lazy arguments of the caller are forced before the closures travel to the nodes
  a <- benchmark_methods("granger", scenarios = c("x_to_y", "independent"), simulator = sim_logistic, n = 120, reps = 2, seed = 1)
  b <- benchmark_methods("granger", scenarios = c("x_to_y", "independent"), simulator = sim_logistic, n = 120, reps = 2, seed = 1, cores = cl)
  expect_equal(a$runs$predicted, b$runs$predicted)
  expect_true(isTRUE(parallel::clusterEvalQ(cl, TRUE)[[1]]))     # still alive
})

test_that("the Rcausality.cores option sets the default and bad values are rejected", {
  d <- sim_logistic(150, "x_to_y", seed = 1)
  old <- options(Rcausality.cores = 1); on.exit(options(old), add = TRUE)
  expect_equal(formals(rcausality)$cores, quote(getOption("Rcausality.cores", 1L)))
  expect_error(rcausality(d, method = "te", lag = 1, n_perm = 20, cores = "many"), "cores")
  expect_equal(Rcausality:::.n_cores(1), 1L)
  expect_true(Rcausality:::.n_cores("auto") >= 1L)
  expect_equal(Rcausality:::.split_cores(5, 2), list(outer = 2, inner = 1L))
  expect_equal(Rcausality:::.split_cores(1, 2), list(outer = 1L, inner = 2))
  expect_equal(Rcausality:::.split_cores(1, 1), list(outer = 1L, inner = 1L))
})

test_that("the socket-cluster path (used on Windows) leaves no connection open", {
  skip_on_cran()
  skip_if_not(parallel::detectCores() >= 2)
  d <- sim_logistic(150, "x_to_y", seed = 1)
  old <- options(Rcausality.psock = TRUE); on.exit(options(old), add = TRUE)
  before <- nrow(showConnections())
  a <- rcausality(d, method = "te", lag = 1, n_perm = 30, seed = 1, cores = 1)
  b <- rcausality(d, method = "te", lag = 1, n_perm = 30, seed = 1, cores = 2)
  expect_equal(a$pairs, b$pairs)
  expect_equal(nrow(showConnections()), before)
})
