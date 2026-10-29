test_that("exact g-prior and hyper-g Bayes factors agree with the BIC one in direction", {
  d <- sim_var(400, "x_to_y", seed = 1)
  bf <- sapply(c("bic", "g", "hyper_g"), function(pr) rcausality(d, method = "bayes_granger", prior = pr)$pairs$bf10)
  # rows: y -> x (no link), x -> y (link)
  expect_true(all(bf[2, ] > 1e6))
  expect_true(all(bf[1, ] < 1))
  r <- rcausality(d, method = "bayes_granger", prior = "hyper_g")
  expect_equal(r$settings$prior, "hyper_g")
  expect_true(r$settings$bayes)
  expect_output(print(r), "hyper-g prior")
})

test_that(".log_bf_g0 reproduces the closed form of the g-prior and is finite for hyper-g", {
  N <- 200; k <- 3; R2 <- 0.3
  expect_equal(Rcausality:::.log_bf_g0(N, k, R2, "g"),
               (N - 1 - k) / 2 * log(1 + N) - (N - 1) / 2 * log(1 + N * (1 - R2)))
  h <- Rcausality:::.log_bf_g0(N, k, R2, "hyper_g")
  expect_true(is.finite(h))
  expect_gt(h, 0)
  expect_lt(Rcausality:::.log_bf_g0(N, k, 0.01, "hyper_g"), h)
})

test_that("bf_from_p is the Sellke-Bayarri-Berger bound", {
  b <- bf_from_p(c(0.05, 0.5, NA), prior_prob = 0.5)
  expect_equal(b$bf10[1], 1 / (-exp(1) * 0.05 * log(0.05)))
  expect_equal(b$bf10[2], 1)
  expect_true(is.na(b$bf10[3]))
  expect_equal(b$post_prob[2], 0.5)
  expect_gt(bf_from_p(0.001)$bf10, 50)
})

test_that("bayes = TRUE calibrates the model-free methods and decides by BF", {
  skip_on_cran()      # long-running; run by devtools::test() / NOT_CRAN=true
  lg <- sim_logistic(500, "x_to_y", seed = 2)
  r <- rcausality(lg, method = "cmi", bayes = TRUE, seed = 1, n_perm = 100)
  P <- r$pairs
  expect_true(all(!is.na(P$bf10)))
  expect_true(all(!is.na(P$significant)))
  expect_true(P$significant[P$cause == "x" & P$effect == "y"])
  expect_false(P$significant[P$cause == "y" & P$effect == "x"])
  expect_equal(r$classes$class, "x_to_y")
  expect_output(print(r), "Bayesian surrogate test")
  expect_output(print(r), "BF10")
  # bayes = TRUE with granger is bayes_granger
  r2 <- rcausality(sim_var(300, "x_to_y", seed = 1), bayes = TRUE)
  expect_equal(r2$settings$method, "bayes_granger")
  # posterior-based methods are left alone
  r3 <- rcausality(sim_logistic_bdc(600, "x_to_y", seed = 1), method = "dc", bayes = TRUE)
  expect_equal(r3$settings$method, "dc")
})

# ---- full Bayesian variants -------------------------------------------------

test_that("bayes_control validates and prints", {
  bc <- bayes_control(prior = "normal", scale = 0.5, case_prior = c(1, 1, 1, 1, 6), n_iter = 100, chains = 1)
  expect_s3_class(bc, "bayes_control"); expect_equal(sum(bc$case_prior), 1)
  expect_output(print(bc), "regression family")
  expect_error(bayes_control(case_prior = 1:3), "five")
})

test_that("surrogate_bf: closed form agrees with Monte Carlo, direction and scale behave", {
  set.seed(1); nl <- rnorm(300)
  a <- surrogate_bf(3, nl, effect_prior = "half_normal", null_model = "normal")$bf10
  b <- surrogate_bf(3, nl, effect_prior = "half_normal", null_model = "kde", n_draws = 20000, seed = 2)$bf10
  expect_equal(log(a), log(b), tolerance = 0.25)
  expect_gt(surrogate_bf(3, nl)$bf10, 10)
  expect_lt(surrogate_bf(0, nl)$bf10, 1)
  expect_lt(surrogate_bf(-2, nl)$bf10, surrogate_bf(2, nl)$bf10)           # one-sided
  expect_gt(surrogate_bf(-3, nl, alternative = "less")$bf10, 10)
  two <- surrogate_bf(2.5, nl, alternative = "two.sided", seed = 1); expect_gt(two$bf10, 1)
  # a known normal null works without surrogate draws; prior probability enters the posterior
  n0 <- surrogate_bf(2.5, list(mean = 0, sd = 1), prior_prob = 0.2)
  expect_equal(n0$post_prob, n0$bf10 * 0.2 / (n0$bf10 * 0.2 + 0.8))
  expect_output(print(n0), "BF10")
  # effect posterior follows the data when the effect is far outside the prior scale (half-Cauchy)
  big <- surrogate_bf(12, nl)
  expect_gt(big$effect[["mean"]], 10); expect_lt(big$effect[["lower"]], 12); expect_gt(big$effect[["upper"]], 12)
  expect_error(surrogate_bf(1, rnorm(5)), "At least 10")
})

test_that("bayes_granger: spike-and-slab and normal-slab evidence, inclusion probabilities, draws", {
  d <- sim_var(400, "x_to_y", seed = 1)
  bc <- bayes_control(n_iter = 400, chains = 1)
  b1 <- bayes_granger(d$y, d$x, lag = 2, control = bc)
  b0 <- bayes_granger(d$x, d$y, lag = 2, control = bc)
  expect_gt(b1$bf10, 100); expect_lt(b0$bf10, 1)
  expect_gt(b1$inclusion[["lag1"]], 0.95); expect_lt(b0$inclusion[["lag1"]], 0.3)
  expect_equal(dim(b1$draws), c(400L, 4L)); expect_true(all(c("sigma", "geweke") %in% colnames(b1$draws)))
  expect_gt(b1$geweke[["lower"]], 0); expect_true(is.finite(b1$ess))
  expect_output(print(b1), "inclusion probability")
  bn <- bayes_granger(d$y, d$x, lag = 2, control = bayes_control(prior = "normal", n_iter = 200, chains = 1))
  expect_gt(bn$bf10, 100); expect_null(bn$inclusion)
  # conditional version with a control series and the closed-form priors through the same interface
  z <- rnorm(400)
  expect_gt(bayes_granger(d$y, d$x, controls = cbind(z), lag = 1, control = bc)$bf10, 100)
  expect_gt(bayes_granger(d$y, d$x, lag = 1, control = bayes_control(prior = "g"))$bf10, 100)
  # many lags: the Gibbs path
  bg <- bayes_granger(d$y, d$x, lag = 13, control = bayes_control(n_iter = 60, burnin = 20, chains = 1))
  expect_length(bg$inclusion, 13L); expect_true(is.finite(bg$bf10))
})

test_that("rcausality(bayes = TRUE) runs every family and fills the Bayesian columns", {
  skip_on_cran()      # long-running; run by devtools::test() / NOT_CRAN=true
  d <- sim_var(400, "x_to_y", seed = 1)
  r <- rcausality(d, bayes = TRUE, control = bayes_control(n_iter = 200, chains = 1))
  expect_equal(r$settings$method, "bayes_granger"); expect_equal(r$settings$prior, "spike_slab")
  expect_equal(r$classes$class, "x_to_y"); expect_false(is.null(r$bayes)); expect_s3_class(r$bayes_fits[[1]], "bayes_granger")
  expect_output(print(r), "spike-and-slab")
  rp <- rcausality(d, bayes = TRUE, control = bayes_control(decision = "posterior", n_iter = 100, chains = 1))
  expect_equal(rp$classes$class, "x_to_y"); expect_output(print(rp), "posterior probability >")
  # closed-form priors through `prior`
  expect_equal(rcausality(d, bayes = TRUE, prior = "hyper_g")$settings$prior, "hyper_g")
  # surrogate family
  rt <- rcausality(d, method = "te", bayes = TRUE, seed = 1)
  expect_true(all(!is.na(rt$pairs$bf10))); expect_equal(rt$classes$class, "x_to_y"); expect_match(rt$bayes$model[1], "surrogate")
  expect_output(print(rt), "Bayesian surrogate test")
  rc <- rcausality(d, method = "cmi_lag", bayes = TRUE, lag = 1, seed = 1)
  expect_true(all(!is.na(rc$instant$bf10)))
  # posterior-based methods with a case prior
  lg <- sim_logistic_bdc(800, "x_to_y", seed = 1)
  r0 <- rcausality(lg, method = "dc", bayes = TRUE)
  r1 <- rcausality(lg, method = "dc", bayes = TRUE, control = bayes_control(case_prior = c(1, 1, 1, 1, 20)))
  expect_lt(r1$pairs$post_prob[r1$pairs$cause == "x"], r0$pairs$post_prob[r0$pairs$cause == "x"])
  # the recurrence criterion keeps its coverage-index condition
  rh <- rcausality(lg, method = "hirata", bayes = TRUE, seed = 1)
  expect_equal(rh$classes$class, "x_to_y")
  # Bayesian decision table of the Markov method
  w <- sim_turnwheel(800, seed = 1)[, c("alice", "bob")]
  rm <- rcausality(w, method = "markov", bayes = TRUE, memory = 1, n_perm = 49, n_markov = 49, seed = 1)
  expect_equal(rm$classes$relation, "alice -> bob"); expect_match(rm$bayes$model[1], "decision table")
  pdf(NULL); on.exit(dev.off()); expect_no_error(plot(r, "posterior")); expect_no_error(plot(rt, "posterior"))
})
