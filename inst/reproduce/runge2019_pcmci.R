## Runge, Nowack, Kretschmer, Flaxman & Sejdinovic (2019), Science Advances
## 5(11), eaau4996 -- PCMCI.  Optional wrapper around the Python package
## tigramite (needs reticulate + `pip install tigramite`).
##
## Reproduced: the four-variable linear toy model of the tigramite tutorial
## (the authors' reference implementation accompanying the paper):
##   X0_t = 0.7 X0_{t-1} - 0.8 X1_{t-1} + e0
##   X1_t = 0.8 X1_{t-1} + 0.8 X3_{t-1} + e1
##   X2_t = 0.5 X2_{t-1} + 0.5 X1_{t-2} + 0.6 X3_{t-3} + e2
##   X3_t = 0.7 X3_{t-1} + e3
## PCMCI with ParCorr and tau_max = 3 should recover exactly the four lagged
## cross-links X1->X0 (lag 1), X3->X1 (lag 1), X1->X2 (lag 2), X3->X2 (lag 3)
## with the signs of the coefficients, and no spurious link.
library(Rcausality)
if (!pcmci_available()) stop("tigramite is not available; set RETICULATE_PYTHON or pip install tigramite")

sim_runge_toy <- function(n = 1000, seed = 42) {
  set.seed(seed); T <- n + 100
  X <- matrix(0, T, 4); E <- matrix(stats::rnorm(T * 4), T, 4)
  for (t in 4:T) {
    X[t, 1] <- 0.7 * X[t - 1, 1] - 0.8 * X[t - 1, 2] + E[t, 1]
    X[t, 2] <- 0.8 * X[t - 1, 2] + 0.8 * X[t - 1, 4] + E[t, 2]
    X[t, 3] <- 0.5 * X[t - 1, 3] + 0.5 * X[t - 2, 2] + 0.6 * X[t - 3, 4] + E[t, 3]
    X[t, 4] <- 0.7 * X[t - 1, 4] + E[t, 4]
  }
  d <- as.data.frame(X[-(1:100), ]); names(d) <- paste0("X", 0:3); d
}
d <- sim_runge_toy()
r <- pcmci(d, tau_max = 3)
print(r)
cat("\nSignificant lagged links (expected: X1->X0 @1 (-), X3->X1 @1 (+), X1->X2 @2 (+), X3->X2 @3 (+)):\n")
print(r$pairs[r$pairs$significant, c("cause", "effect", "lag", "statistic", "p_adj")], row.names = FALSE)
# PCMCI+ on the same data (no contemporaneous links in the model)
rp <- pcmci(d, tau_max = 3, plus = TRUE)
cat("\nPCMCI+ contemporaneous links found:", sum(rp$instant$significant), "(expected 0)\n")
plot(r, type = "both")
