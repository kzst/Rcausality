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
  lg <- nl_pair(500, 2)
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
  skip_on_cran()
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

test_that("the g-prior scale and the hyper-g parameter of bayes_control() take effect", {
  skip_on_cran()
  d <- sim_var(200, "x_to_y", seed = 1)
  bf <- function(...) bayes_granger(d$y, d$x, lag = 2, control = bayes_control(...))$bf10
  # g = scale * N: the closed form against its definition
  n <- 200; L <- 2; yd <- d$y[(L + 1):n]; N <- length(yd)
  lags <- function(v) sapply(1:L, function(l) v[(L + 1 - l):(n - l)])
  # the restricted regression is the base model: flat priors on its coefficients, the g prior on the lags of the cause
  # after the common regressors are projected out, so that the Bayes factor is a function of F and its degrees of
  # freedom only: BF10 = (1 + g)^((nu0 - p) / 2) (1 + g (1 - R2p))^(-nu0 / 2), R2p = p F / (p F + nu1)
  Xr <- lags(d$y); Xu <- cbind(Xr, lags(d$x))
  rss <- function(X) sum(stats::resid(stats::lm(yd ~ X))^2)
  nu0 <- N - 1 - 2; nu1 <- nu0 - 2; Fst <- ((rss(Xr) - rss(Xu)) / 2) / (rss(Xu) / nu1); R2p <- 2 * Fst / (2 * Fst + nu1)
  expect_equal(R2p, 1 - rss(Xu) / rss(Xr))
  for (s in c(1, 10))
    expect_equal(log(bf(prior = "g", scale = s)), (nu0 - 2) / 2 * log1p(s * N) - nu0 / 2 * log1p(s * N * (1 - R2p)), tolerance = 1e-8)
  # ... and against its definition: the marginal likelihoods integrated numerically over sigma^2 (the coefficients
  # are integrated in closed form for a given sigma^2)
  g <- N; S0 <- rss(Xr); S1 <- S0 - g / (1 + g) * (S0 - rss(Xu))
  m <- function(S, const) { u0 <- log(S / nu0)          # integration over u = log(sigma^2) around the mode
    stats::integrate(function(u) exp(const - nu0 / 2 * u - S * exp(-u) / 2 + nu0 / 2 * log(S0 / nu0) + nu0 / 2),
                     u0 - 3, u0 + 3, rel.tol = 1e-10, subdivisions = 400L)$value }
  expect_equal(bf(prior = "g"), m(S1, -log1p(g)) / m(S0, 0), tolerance = 1e-6)
  # the evidence does not depend on how well the own past predicts: the Bayes factor is the same for two data sets
  # with the same F statistic and degrees of freedom
  expect_equal(bf(prior = "g"), exp(Rcausality:::.log_bf_nested(N, 2, 4, 0.999, 1 - 0.001 * rss(Xu) / rss(Xr), "g")), tolerance = 1e-8)
  expect_false(isTRUE(all.equal(bf(prior = "hyper_g", hyper_a = 3), bf(prior = "hyper_g", hyper_a = 4))))
  # the same through rcausality(): Granger and Toda-Yamamoto regressions
  r1 <- rcausality(d$x, d$y, method = "bayes_granger", lag = 2, control = bayes_control(prior = "g"))
  r2 <- rcausality(d$x, d$y, method = "bayes_granger", lag = 2, control = bayes_control(prior = "g", scale = 10))
  expect_equal(r1$pairs$bf10[1], bf(prior = "g")); expect_equal(r2$pairs$bf10[1], bf(prior = "g", scale = 10))
  expect_output(print(r2), "g = 10 n")
  t1 <- rcausality(d$x, d$y, method = "toda_yamamoto", bayes = TRUE, lag = 2, control = bayes_control(prior = "hyper_g"))
  t2 <- rcausality(d$x, d$y, method = "toda_yamamoto", bayes = TRUE, lag = 2,
                   control = bayes_control(prior = "hyper_g", hyper_a = 4))
  expect_false(isTRUE(all.equal(t1$pairs$bf10, t2$pairs$bf10)))
  # a prior without a closed form for the lag-augmented regressions falls back to (and records) the BIC one
  expect_identical(rcausality(d$x, d$y, method = "toda_yamamoto", bayes = TRUE, lag = 2)$settings$prior, "bic")
  # panels: unit and pooled Bayes factors
  pd <- sim_panel(n_units = 8, n_time = 30, scenario = "x_to_y", seed = 1)
  p1 <- panel_causality(pd, id = "unit", time = "time", bayes = TRUE, lag = 1, control = bayes_control(prior = "g"), seed = 1)
  p2 <- panel_causality(pd, id = "unit", time = "time", bayes = TRUE, lag = 1, control = bayes_control(prior = "g", scale = 10),
                        seed = 1)
  expect_false(isTRUE(all.equal(p1$units$bf10_raw, p2$units$bf10_raw)))
  expect_false(isTRUE(all.equal(p1$panel$bf10[p1$panel$level == "pooled"], p2$panel$bf10[p2$panel$level == "pooled"])))
})

test_that("case_prior replaces the prior of the method over the causal cases", {
  p <- c(x_to_y = 0.5, y_to_x = 0.1, bidirectional = 0.1, common_cause = 0.2, independent = 0.1)
  base9 <- Rcausality:::.dc_case_base
  # the method's own prior as case_prior leaves the posterior unchanged
  expect_equal(Rcausality:::.bayes_cases(p, base9, base9), p)
  expect_equal(Rcausality:::.bayes_cases(p, rep(0.2, 5) |> stats::setNames(names(p))), p)
  # a uniform prior over the five cases removes the 1:1:1:3:3 weights of the nine events
  q <- Rcausality:::.bayes_cases(p, stats::setNames(rep(0.2, 5), names(p)), base9)
  expect_equal(unname(q), unname((p / base9) / sum(p / base9)))
  # Markov table: implied prior 2:2:2:1:1; giving it as case_prior changes nothing
  pM <- c(MA = 0.8, MB = 0.3, MJ = 0.9); pD <- c(sAB = 0.9, sBA = 0.2, sI0 = 0.1)
  b0 <- Rcausality:::.bayes_markov(pM, pD, "allowed")
  imp <- c(x_to_y = 2, y_to_x = 2, bidirectional = 2, common_cause = 1, independent = 1) / 8
  b1 <- Rcausality:::.bayes_markov(pM, pD, "allowed", imp)
  expect_equal(b1$cases, b0$cases); expect_equal(b1$p_ab, b0$p_ab)
  b2 <- Rcausality:::.bayes_markov(pM, pD, "allowed", stats::setNames(rep(0.2, 5), names(imp)))
  expect_equal(sum(b2$cases), 1); expect_false(isTRUE(all.equal(b2$cases, b0$cases)))
})
