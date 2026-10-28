# Reproduction / illustration of
#   Zlatniczki, A., Stippinger, M., Benko, Z., Somogyvari, Z. and Telcs, A. (2021).
#   Relaxation of some confusions about confounders. Entropy 23(11), 1450.
#   doi:10.3390/e23111450
# The paper is theoretical (no numerical experiments); this script checks its
# statements on simulated systems with known ground truth.
library(Rcausality)

## (1) Relation 1 (Sect. 2.1): I_A, I_B (lagged CMI) and I^(k)_{A,B} (k < tau)
##     decide between A->B, B->A, no coupling and a COMMON DRIVER.
set.seed(5); n <- 400
common <- local({                       # shared (i.i.d.) innovation = common driver
  e1 <- rnorm(n); e2 <- 0.8 * e1 + 0.6 * rnorm(n)
  data.frame(x = as.numeric(stats::filter(e1, .5, "recursive")),
             y = as.numeric(stats::filter(e2, .5, "recursive")))
})
cases <- list(
  "A -> B (direct driving)"   = sim_var(n, "x_to_y", seed = 1),
  "no coupling"               = sim_var(n, "independent", seed = 2),
  "common driver"             = common)
for (nm in names(cases)) {
  cat("\n====", nm, "====\n")
  print(rcausality(cases[[nm]], method = "cmi_lag", n_perm = 100, seed = 1))
}

## (2) Sect. 2.7 / Appendix A: with a common driver the k = 0 CMI does NOT
##     vanish as the conditioning history d grows; without it, it stays ~ 0.
cat("\nCMI I(x_n; y_n | histories of depth d)\n")
print(cbind(common_driver = cmi_convergence(common$x, common$y, 0, 1:4)$cmi,
            direct_only   = cmi_convergence(cases[[1]]$x, cases[[1]]$y, 0, 1:4)$cmi))

## (3) Sects. 2.4-2.6: geodesic triangle inequality on lags
lag_triangle_check(tau_LX = 2, tau_LY = 5, tau_XY = 1)  # reversed strict: direct driving
lag_triangle_check(tau_LX = 2, tau_LY = 2, tau_XY = 1)  # strict triangle: common driver
lag_triangle_check(tau_LX = 2, tau_LY = 3, tau_XY = 1)  # equality: confounder, ambiguous
