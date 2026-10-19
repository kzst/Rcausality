# Rcausality 0.0.0.9000

* Initial skeleton (M0): `causres` class with `print`/`summary`/`plot`,
  Granger and Bayesian (BIC Bayes factor) Granger tests, bivariate, all-pairs
  and MIMO block modes, simulators, preprocessing helpers.
* `method = "cmi_lag"`: lag-shifted CMI scheme of Zlatniczki et al. (2021)
  with contemporaneous common-driver test (`common_cause` class);
  `cmi_convergence()`, `lag_triangle_check()`.
