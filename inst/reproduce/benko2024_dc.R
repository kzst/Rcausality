# Reproduction of
#   Benko, Z., Zlatniczki, A., Stippinger, M., Fabo, D., Solyom, A., Eross, L.,
#   Telcs, A. and Somogyvari, Z. (2024). Bayesian inference of causal relations
#   between dynamical systems. Chaos, Solitons & Fractals 185, 115142.
#   doi:10.1016/j.chaos.2024.115142
# Figure 1 and Appendix B.2: three coupled logistic maps (Eq. B.1, r = 3.99,
# beta = 0.5, N = 10000, rank normalisation, tau = 1, M = 4, k = 12..44),
# the third map hidden.  Expected: the highest posterior on the true case in
# all four coupling patterns (unidirectional, bidirectional, common cause,
# independent).
library(Rcausality)
for (sc in c("x_to_y", "bidirectional", "common_cause", "independent")) {
  d <- sim_logistic_bdc(10000, sc, seed = 1)
  r <- dimensional_causality(d$x, d$y, m = 4, tau = 1, k_range = seq(12, 44, by = 4), seed = 1)
  cat("\n==== truth:", sc, "====\n"); print(r)
}
# Fig. 1B-C style panels for the common-cause case
d <- sim_logistic_bdc(10000, "common_cause", seed = 1)
plot(dimensional_causality(d$x, d$y, m = 4, k_range = seq(12, 44, by = 4), seed = 1))
