# Reproduction of
#   Stippinger, M., Bencze, A., Zlatniczki, A., Somogyvari, Z. and Telcs, A. (2023).
#   Causal discovery of stochastic dynamical systems: a Markov chain approach.
#   Mathematics 11(4), 852. doi:10.3390/math11040852
# Sect. 6.1 ("turn-wheel", Figure 2): Alice (driver) is a two-state Markov chain,
# Bob and Barry are driven by Alice with one step delay. Expected: Alice Markov
# (large p); Bob, Barry NOT Markov (small p) because Alice is a hidden memory;
# the joint (Alice_t, Bob_{t+1}) IS Markov; the joint (Bob, Barry) is not.
library(Rcausality)
reps <- 20L; n <- 5000L                       # paper: 100 repetitions, n = 2000...50000
p <- t(sapply(seq_len(reps), function(r) {
  w <- sim_turnwheel(n, seed = r); N <- n
  mp <- function(x) markov_test(x, m = 1, n_perm = 99, seed = r)$p_value
  c(alice = mp(w$alice), bob = mp(w$bob), barry = mp(w$barry),
    `alice0 x bob1` = mp(cbind(w$alice[-N], w$bob[-1])),
    `bob0 x barry0` = mp(cbind(w$bob, w$barry)))
}))
cat("Median p-values of the Markov test over", reps, "repetitions, n =", n, "\n")
print(round(apply(p, 2, median), 3))
cat("Rejection rate at 5%:\n"); print(colMeans(p < 0.05))

# Decision table on one shorter realisation (the kNN-CMI part is O(n^2))
w <- sim_turnwheel(1000, seed = 1)
print(rcausality(w[, c("alice", "bob")], method = "markov", memory = 1,
                 n_perm = 99, n_markov = 99, seed = 1))

## --- Sect. 6.1, coupled logistic maps (Figure 3: no dynamical noise) --------
## x0 autonomous driver, x1 and x2 driven by x0 (b = 1.5), embedding m = 3.
## Expected (Fig. 3, n = 2000 column): driver Markov (median p about 0.45),
## driven1/driven2 not Markov (about 0.01), joint (driver_t, driven1_{t+1})
## Markov (0.2-0.3), joint (driven1, driven2) not Markov (about 0.02).
## The paper adds observational noise of unstated size; here obs_noise = 0.
reps <- 5L
q <- t(sapply(seq_len(reps), function(s) {
  d <- sim_logistic3(2000, "driver", sigma = 0, obs_noise = 0, seed = s)
  E <- lapply(d, embed_delay, m = 3, tau = 1); N <- nrow(E[[1]])
  mp <- function(A) markov_test(A, n_perm = 49, seed = s)$p_value
  c(driver = mp(E$x0), driven1 = mp(E$x1), driven2 = mp(E$x2),
    `driver:0 x driven1:1` = mp(cbind(E$x0[-N, ], E$x1[-1, ])),
    `driven1:0 x driven2:0` = mp(cbind(E$x1, E$x2)))
}))
cat("\nCoupled logistic maps, median p-values over", reps, "runs (n = 2000):\n")
print(round(apply(q, 2, median), 2))
