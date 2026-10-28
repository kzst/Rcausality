# Reproduction of
#   Stippinger, M., Varga, B., Benko, Z., Fabo, D., Eross, L., Telcs, A. and
#   Somogyvari, Z. CCDH: Complexity based causal discovery of hidden common
#   cause in time series. Chaos, Solitons & Fractals 176, 114054. doi:10.1016/j.chaos.2023.114054
# Section 4, Figure 1: three coupled logistic maps (Eqs. 14-17; r random in
# [3.5699, 3.98), couplings in [0.3, 0.5)), n = 1000, m = 3, tau = 1, rank-order
# coding, compression rates after 1..9 NSRPS steps, averaged over 100 systems.
# Expected: (a) driver0 -> driven1: the joints are as compressible as the
# driven series; (b) hidden common driver (driven1, driven2): the joints are
# less compressible than the singles; (c) mutual driving: all curves coincide.
library(Rcausality)
reps <- 25L   # paper: 100
curves <- function(sc, pair) {
  cr <- Reduce(`+`, lapply(seq_len(reps), function(s) {
    set.seed(s)
    d <- sim_logistic3(1000, sc, b = runif(1, .3, .5), c = runif(1, .3, .5), seed = s)
    # sim_logistic3 draws r in [3.8, 3.9); the paper uses [3.5699, 3.98)
    ccdh(d[[pair[1]]], d[[pair[2]]], m = 3, K = 9, n_surr = 10, seed = s)$cr
  })) / reps
  cr
}
a <- curves("driver", c("x0", "x1")); b <- curves("driver", c("x1", "x2")); c <- curves("mutual", c("x1", "x2"))
op <- par(mfrow = c(1, 3))
for (p in list(list(a, "(a) driver0 and driven1"), list(b, "(b) hidden common driver"), list(c, "(c) bidirectional coupling"))) {
  matplot(0:9, rbind(1, p[[1]]), type = "b", pch = c(16, 16, 1, 1), lty = c(1, 1, 2, 2),
          col = c("#C62828", "#1565C0", "#E57373", "#64B5F6"), xlab = "number of compression steps",
          ylab = "compression rate", main = p[[2]], ylim = c(0.2, 1))
  legend("topright", colnames(p[[1]]), col = c("#C62828", "#1565C0", "#E57373", "#64B5F6"), lty = c(1, 1, 2, 2), bty = "n")
}
par(op)
cat("mean compression rates (steps 1-9):\n")
print(round(rbind(`(a) x0,x1` = colMeans(a), `(b) x1,x2 common` = colMeans(b), `(c) mutual` = colMeans(c)), 3))
# verdict frequencies with the default decision rule
tab <- sapply(list(a = c("driver", "x0", "x1"), b = c("driver", "x1", "x2"), c = c("mutual", "x1", "x2")), function(z) {
  v <- sapply(seq_len(reps), function(s) { set.seed(s)
    d <- sim_logistic3(1000, z[1], b = runif(1, .3, .5), c = runif(1, .3, .5), seed = s)
    ccdh(d[[z[2]]], d[[z[3]]], m = 3, n_surr = 30, seed = s)$verdict })
  table(factor(v, c("x_to_y", "y_to_x", "bidirectional", "common_cause", "independent")))
})
print(tab)
