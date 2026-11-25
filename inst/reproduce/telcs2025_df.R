# Reproduction of
#   Telcs, A., Kurbucz, M. T. and Jakovac, A. (2025). Unified causality analysis
#   based on the degrees of freedom. Physical Review E (manuscript, round II).
# Sect. 6.4 (linear rotation system, Eqs. 25-30, Figures 1-2 and 4) and Sect. 7
# (chickens and eggs, Figures 6-7).
library(Rcausality)
cores <- 2L                                   # raise on a multi-core machine

## --- Sect. 6.4, "deterministic" parameter set (Eq. 29) ---------------------
## z: autonomous 2-D rotation -> one constrained pair determines z[n]:
##    the N = 0 line is flat, the N >= 1 lines are proportional to sigma_W (Fig. 1)
## x: driven by z -> N = 0 and N = 1 saturate, N >= 2 proportional (Fig. 2):
##    df = 2 pairs (= 4 in the paper's scalar count)
d <- sim_df_rotation(1500, "deterministic", seed = 1)
dz <- df_variance(d[, c("z1", "z2")], max_constr = 4, seed = 1, cores = cores)
dx <- df_variance(d[, c("x1", "x2")], max_constr = 4, seed = 1, cores = cores)
dz; dx
op <- par(mfrow = c(1, 2)); plot(dz, main = "z subsystem (Fig. 1, left)")
plot(dx, main = "x subsystem (Fig. 2, left)"); par(op)

## --- Sect. 6.4, "stochastic" parameter set (Eq. 30) -------------------------
ds <- sim_df_rotation(1500, "stochastic", seed = 1)
dzs <- df_variance(ds[, c("z1", "z2")], max_constr = 4, seed = 1, cores = cores)
dxs <- df_variance(ds[, c("x1", "x2")], max_constr = 4, seed = 1, cores = cores)
dzs; dxs
op <- par(mfrow = c(1, 2)); plot(dzs, main = "z subsystem (Fig. 1, right)")
plot(dxs, main = "x subsystem (Fig. 2, right)"); par(op)

## --- Sect. 6.4.5: the x-y subsystem, a hidden common driver -----------------
## One component of each oscillator observed (x1, y1); the orders, in lags of
## the observed series: 4 for x1, 4 for y1 and 6 for the joint process, so that
## max(4, 4) < 6 < 4 + 4: a hidden common driver (Table 2, row 5; Fig. 4 counts
## the same six conditions as three pairs).  The paper computes this linear
## example with the conditional normal distribution; the Gaussian conditional
## variance resolves the orders, which the windows of the nearest neighbours
## cannot at n = 1500 (in the nearly deterministic set the joint order stops at
## four).
rg <- rcausality(d[, c("x1", "y1")], method = "df", seed = 1, cores = cores,
                 method_args = list(estimator = "gaussian"))
rg$df                                         # orders 4, 4, 6: common_cause
rgs <- rcausality(ds[, c("x1", "y1")], method = "df", seed = 1, cores = cores,
                  method_args = list(estimator = "gaussian"))
rgs$df                                        # stochastic set: orders 3, 4, 6, common_cause
rk <- rcausality(d[, c("x1", "y1")], method = "df", seed = 1, cores = cores)
rk$df                                         # the windows of the neighbours (default)

## --- Sect. 7: chickens and eggs (Thurman & Fisher 1988 data, lmtest::ChickEgg)
## The paper differentiates and normalises both series and uses lag 2;
## constraining the egg past changes the chicken distribution (egg -> chicken),
## constraining the chicken past leaves the egg distribution unchanged.
if (requireNamespace("lmtest", quietly = TRUE)) {
  data("ChickEgg", package = "lmtest")
  ce <- prep_series(as.data.frame(ChickEgg), diff = 1, scale = TRUE)
  r <- rcausality(ce, method = "df", lag = 2, history = 1, n_perm = 200, seed = 1, cores = cores)
  print(r)
}
