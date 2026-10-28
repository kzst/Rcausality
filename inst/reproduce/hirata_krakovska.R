# Checks of the two literature-extension methods on the coupled logistic maps
# of Benko et al. (2024) (sim_logistic_bdc, r = 3.99, beta = 0.5), n = 1500:
#   Hirata, Y. and Aihara, K. (2010). Identifying hidden common causes from
#     bivariate time series by using the recurrence plot. Phys. Rev. E 81, 016203.
#   Krakovska, A. (2019). Correlation dimension detects causal links in coupled
#     dynamical systems. Entropy 21(9), 818.
# The original papers' own test systems (wind data; Rossler/Lorenz systems)
# are not reproduced here; the scripts shows that the two criteria of Hirata &
# Aihara (coverage, simultaneous recurrences) and the dimension relations of
# Krakovska separate the five canonical cases on this benchmark.
library(Rcausality)
classes <- c("x_to_y", "y_to_x", "bidirectional", "common_cause", "independent")
res <- do.call(rbind, lapply(1:5, function(s) do.call(rbind, lapply(classes, function(sc) {
  d <- sim_logistic_bdc(1500, sc, seed = s)
  data.frame(seed = s, truth = sc,
             hirata = hirata_rp(d$x, d$y, E = 3, n_surr = 50, seed = 1)$verdict,
             cd = krakovska_cd(d$x, d$y, E = 3, seed = 1)$verdict)
}))))
cat("Hirata & Aihara recurrence criteria:\n"); print(table(truth = res$truth, verdict = res$hirata))
cat("\nKrakovska correlation dimensions:\n"); print(table(truth = res$truth, verdict = res$cd))
