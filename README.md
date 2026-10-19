# Rcausality

Causality analysis for time series in R -- a unified interface to Granger
causality, its Bayesian (Bayes-factor) variant and, in coming releases, the
Markov-chain, conditional-mutual-information, dimensional (Bayesian),
entropy/random-forest and PID-based methods.

> Status: **0.0.0.9000 - development skeleton (M0)**. License: GPL-3.

```r
library(Rcausality)
d <- sim_var(500, scenario = "x_to_y", seed = 1)

r <- rcausality(d$x, d$y)         # two series (both directions)
r <- rcausality(df)               # data frame: all ordered pairs
r <- rcausality(X, Y)             # MIMO: blocks of series
r <- rcausality(df, conditional = TRUE)   # conditional on the rest

r; summary(r); plot(r, "both")    # 'causres' object with print/summary/plot
as.matrix(r, "p_adj")             # cause x effect matrices
```

Every test reports lag, F statistic, df, p-value, adjusted p-value, BIC Bayes
factor, posterior probability, evidence label, effect size (Geweke F-measure,
R2 gain), sign of the net influence and the decision; pairs are classified
as `x_to_y`, `y_to_x`, `bidirectional`, `independent` (and `common_cause`
for methods that can detect it).

## Roadmap
M1 (done: kNN-CMI + local permutation `method = "cmi"`; lag-shifted CMI + common-driver test `method = "cmi_lag"`, `cmi_convergence()`, `lag_triangle_check()`) Markov test, lag-shifted CMI (Zlatniczki), df-variance; M2 Bayesian dimensional
causality, CCDH; M3 df entropy + random forest, TE/CCM baselines, `compare_methods()`;
M4 vignettes + CRAN checks; M5 Hirata, Krakovska, MAR-PID, PCMCI wrapper;
`forecast()`/`predict()` methods; Bayesian variants of all methods.

## Ported code
Any code ported from the authors' Python reference implementations
(`dimensional_causality`, `df_causality`) will be explicitly marked as
**ported by the Rcausality authors**, in the source header, in `?` help pages
and in `inst/COPYRIGHTS`.
