# Rcausality 0.0.11

* **Panels: a warning for strongly persistent units with `"cmi"`,
  `"cmi_lag"` and `"diks_panchenko"`.** The unit tests of these methods
  refer the statistic of a unit to shifted series of that unit. For units
  that are short relative to their persistence such tests reject somewhat
  too often, and the group test, which adds the standardised unit
  statistics, adds this excess over the units. In simulated panels of 30
  independent units of 50 periods with a persistence of 0.9 in both series
  the unit tests rejected in 9% (`"cmi"`) and 7% (`"diks_panchenko"`) of the
  units, the group tests in 47% and 18% of the panels, the pooled tests in
  2% and 8%, and the default decisions in 35% and 13%; with a persistence
  of 0.5 the group test of `"cmi"` rejected in 6% of the panels and its
  decision in 3%. A panel analysis with one of these methods now warns
  when the median first-order autocorrelation of the units, corrected for
  their length, exceeds 0.8 (the limit of the Granger family stays 0.9).
  Differenced series are advisable; for the two information measures the
  pooled level alone, `panel_control(decision = "pooled")`, kept its level
  in these panels. No result changes.
* `benchmark_methods()`: the documentation of `seed` says that the data
  sets are simulated under the generator of the parallel streams
  (`RNGkind("L'Ecuyer-CMRG")`), which is needed to obtain one of them
  outside the function.

# Rcausality 0.0.10

* **Partial information decomposition: conditional shifts as an option.**
  `method = "pid"` tests the novel information of the cause by permuting
  the cause within the levels of the coded own past of the effect, an
  exact test for the coded variables. A code of one bit does not resolve
  the past of a persistent effect: when the past of the cause is a
  function of that past, as in the reverse direction of a link, the coded
  variables are conditionally dependent although the series are not. In
  the reverse direction of a strong link whose driver has an
  autocorrelation of 0.9 the test rejected in 7% to 12% of series of 100
  and 300 time points and in 37% to 55% of series of 1000. `method_args =
  list(reference = "conditional")` uses conditional shifts of the lags of
  the cause given the past of the effect, coded after the shift, with the
  lag scan limited as for `"te"`: 1% and 3% of the same series at 300 time
  points, 2% (linear link) and 10% (nonlinear link, 60 series) at 1000,
  with the power of the default for a weak link (38% against 42%). The
  default is unchanged, and so is every result of the default settings;
  the warning that version 0.0.9 gives in this situation names the option.
  `$extra` of a `"pid"` analysis records the reference of every test.

# Rcausality 0.0.9

* **Kernel test (`"diks_panchenko"`): conditional shifts are the default
  reference in all cases.** Up to version 0.0.8 the test with one lag and
  without controls used the asymptotic normal reference of Diks and
  Panchenko (2006). That reference presupposes that the bias of the kernel
  estimates is negligible beside the standard error of the statistic. At
  the bandwidths of the rule (1.5 standard deviations up to about 350
  observations) it is not when the past of the cause depends on the past
  of the effect and the effect is autocorrelated, which is the reverse
  direction of every link whose driver is persistent: the kernel window
  does not resolve the past of the effect, and the cause, a function of
  that past, keeps information about the present of the effect. In
  simulations the normal reference rejected the true null hypothesis of
  the reverse direction in 10% to 30% of series of 300 time points for a
  moderate link whose driver has an autocorrelation of 0.8 to 0.9, and in
  51% to 82% for a strong one. Surrogates that keep the dependence of the
  cause on the past of the effect have the same bias: with conditional
  shifts the test rejected in 3% to 4% of the same series, and it was
  more powerful in the true direction. `method_args = list(reference =
  "normal")` (and `diks_panchenko(reference = "normal")` or `n_perm = 0`)
  gives the asymptotic test, which needs no surrogates and has no smallest
  attainable p-value; it is appropriate when the pasts of the two series
  are at most weakly dependent, and a warning names the situation
  otherwise. The p-values of the default kernel test change in every
  analysis, and the test is slower (the statistic is recomputed for every
  surrogate) and subject to the resolution of circular shifts (`n_perm`,
  number of admissible shifts) as the other tests with surrogates are.
* **Transfer entropy: conditional shifts as an option.** `method = "te"`
  refers the statistic to circular shifts of the coded cause, which
  generate the hypothesis that the series are independent. With four
  symbols the past of a persistent effect is not resolved either, and in
  the reverse direction of a strong link whose driver has an
  autocorrelation of 0.9 the test rejected in 16% to 25% of series of 1000
  time points (4% to 7% at 300). `method_args = list(reference =
  "conditional")` (and `transfer_entropy(reference = "conditional")`) uses
  conditional shifts of the cause given the past of the effect, coded
  after the shift: 3% to 6% in the same series, with the same power. The
  default stays the circular shift, which is the reference that works for
  noise-free deterministic systems and for symbolic series (these are
  always shifted circularly).
* A warning for the tests that condition coarsely (`"te"` with circular
  shifts, `"pid"`, the kernel test with its normal reference): when a test
  is significant, the effect is strongly autocorrelated and its past
  explains a large share of the variance of the past of the cause, the
  user is told that the rejection can be due to the coding and what to
  use instead.
* **Regression tests: `vcov = "HAC"`.** The Wald test of `"granger"` and
  `"toda_yamamoto"` can use a covariance that adds the autocovariances of
  the scores with Bartlett weights up to lag `floor(4 (n / 100)^(2/9))`
  (Newey and West 1987). The errors of a linear regression on a nonlinear
  process are not martingale differences: with a quadratic link and no
  linear effect the classical F test rejected in 35% to 39% of the
  series, the Wald test with the covariance HC3 in 13% to 14%, and with
  the HAC covariance in 11%, 9% and 7% of series of 100, 300 and 1000
  time points (it rejects somewhat too often in short series: 7% at 100
  time points for independent autoregressions).

# Rcausality 0.0.8

Wording only; no result of an analysis changes.

* The print methods, the help pages, the README and the vignettes no
  longer say that the surrogates of the conditional tests "respect the
  conditioning", which claims more than they do. They are named for what
  they keep: conditional shifts keep the location and the scale of the
  cause given the conditioning past and shift the residuals circularly.

# Rcausality 0.0.7

* **Scale of the conditional surrogates**: the conditional shifts of
  `"cmi"`, `"cmi_lag"` and `"rf_granger"`, and of `"diks_panchenko"` with
  controls or several lags, write the lags of the cause as a location plus
  a scale times a standardized residual. The scale was the unweighted
  least-squares fit of the absolute residuals on a polynomial of the
  conditioning variables. It is now estimated more precisely where it
  depends on the size of a conditioning variable: the absolute values of
  the conditioning variables are among the candidate bases, and the
  selected basis is refitted twice with the weights `1 / s^2` (an absolute
  residual has a standard deviation proportional to the scale, so that an
  unweighted fit is imprecise where the scale is small). The lower bound
  of the scale is a tenth of the mean absolute residual (before: a fifth).
  For a cause whose scale is set by a control (`x = z * noise`) the kernel
  test rejected a true conditional null hypothesis in 10% to 12% of series
  of 100, 300 and 1000 time points at the 5% level, and rejects in 7% to
  9% of the same series now; with the true scale it rejects in 3%, so that
  the rest is the error of the estimated scale, the limit of these
  surrogates. For a volatility that the cause shares with the effect the
  rate is about 5% (before: 6% to 7%). The rates of the other conditional
  tests did not change in the same comparison. The p-values of the tests
  named above can change slightly in any data, since the estimated scale
  enters every surrogate.
* The print method of a Bayesian panel analysis labels the prevalence as
  what it is: the model-based share of the units with a detectable link.
* The repository is `kzst/Rcausality` in the URLs of the package.

# Rcausality 0.0.6

Corrections found while the illustrations of the accompanying paper were
rerun with version 0.0.5. No test statistic and no reference distribution
of an analysis changes; the placebo test of two short series does.

* **Dependence between the units of a panel**: panel analyses compute the
  CD statistic of Pesaran (2021) for every variable (the scaled sum of the
  correlations between the units of the residuals of a first-order
  autoregression with trend) and warn below p = 0.001 that the units move
  together, unless `panel_control(time_effects = TRUE)` or, for the
  Granger family in balanced panels, `resample = "periods"` accounts for
  it. All levels of a panel analysis treat the units as independent, and
  common shocks or trends then produce links between unrelated series: for
  the Grunfeld data the default analysis declared links that 30 to 90% of
  placebo panels show as well, without a warning. The statistic is
  referred to a chi-square approximation instead of the normal one, which
  is too thin in the upper tail of short panels (0.3 to 1.1% of
  independent panels of 10 to 20 periods beyond its 0.1% point; 0.03 to
  0.17% with the chi-square form). The new slot `dependence` holds the
  statistics. A factor with loadings of both signs is not detected.
* **`placebo_test()` with two series**: two series of length `n` have only
  `n - 2 * floor(min_shift * n) + 1` relative circular offsets. When
  `n_placebo` reaches that number every offset is now used once, with a
  message. Before, the offsets were drawn with repetition and the
  smallest placebo p-value, `1 / (n_placebo + 1)`, stated a resolution
  that the data do not have (0.005 instead of 0.022 for 53 observations
  and `n_placebo = 199`).
* **Reruns are silent**: the warnings of the analyses that
  `placebo_test()`, `robustness()` and `predict()` run on placebo data,
  under alternative settings or on a training part were passed on for
  time series (not for panels), where they read as warnings about the
  data of the user.
* **Printing**: a Bayes factor of the surrogate test that lies beyond the
  range of the surrogates is shown as "> bound" in the main table as well
  (it showed the continued value there and the bound in the table of the
  evidence).
* The help text of `effect_scale = "unit"` no longer calls the scale a
  unit-information prior for statistics whose null standard deviation
  need not fall with the root of the length.

# Rcausality 0.0.5

Corrections and additions from the second round of the review of the
accompanying paper. Several of them change results.

## Reference distributions

* **Conditional shifts** (`"diks_panchenko"` with controls or more lags,
  `"cmi"`, `"cmi_lag"`, `"rf_granger"`): the surrogates of the cause keep
  its *location and scale* given the conditioning vector, both estimated
  flexibly (the polynomial basis with the smallest leave-one-out error,
  plus a nearest-neighbour mean of its residuals where that predicts them
  out of sample; for the forest test an out-of-bag forest on the residuals
  of the polynomial). Version 0.0.4 kept a *linear* conditional mean only
  (a linear fit plus a forest for the forest test). A cause that is a
  nonlinear or a scale function of a control was then a better proxy of
  that function than its surrogates, and a true conditional null
  hypothesis was rejected: with a cause `z^2 + noise` and an effect that
  depends on `z^2` the kernel test rejected in 99%, `"cmi"` in 90% and the
  forest test in 17% of the simulated data sets of 300 time points at the
  5% level (now 2.5%, 4% and 3%); with a control that drives the scale of
  the cause the kernel test rejected in 85 to 89% (now 7%). The power for
  a linear and a quadratic link is unchanged within simulation error.
* **A cause that the conditioning past determines**: when the past of the
  effect (and of the controls) explains 99% or more of the variance of a
  lag of the cause, no surrogate differs from the cause by more than the
  error of the fitted location, and what a test would measure is that
  error. No test is carried out then (p-value 1, Bayes factor 1) and a
  warning explains it. This is the case of noise-free deterministic
  systems, in which the driver is a function of the past of the driven
  series (Takens 1981; for coupled logistic maps the driver is a quadratic
  polynomial of two values of the response): the conditional tests
  `"cmi"` and `"rf_granger"`, which 0.0.4 let "detect" such links through
  surrogates that were worse proxies than the cause, are not informative
  for them. `"te"` (plain shifts) and the state-space methods (`"ccm"`,
  `"dc"`, `"cd"`, `"ccdh"`, `"hirata"`) are the tools for deterministic
  systems; with observational or dynamic noise the conditional tests are
  informative again.
* **Panel bootstrap** (`"granger"`, `"toda_yamamoto"`): the cause is
  regenerated together with the effect under the null hypothesis, from its
  own lags and the lags of the bootstrap effect
  (`panel_control(bootstrap = "joint")`, the new default). Keeping the
  cause as observed (0.0.4; now `bootstrap = "effect"`) removes the
  feedback from the effect to the cause from the bootstrap panels and with
  it the bias that this feedback leaves in the pooled estimator: in panels
  of 100 units of 10 periods the reverse direction of a link was declared
  in 11% of the simulated panels at the 5% level (now 3%).
* **Dependence between the units**: `panel_control(resample = "periods")`
  uses the same resampled periods in all units of a balanced panel, which
  keeps common shocks and factors with any loadings in the reference
  distribution (with the effect regenerated and the cause kept, this is
  the bootstrap of Dumitrescu and Hurlin 2012, Section 6.2). With a common
  factor of loading 0.7 the default tests rejected in 45 to 52% of the
  panels without a link; with common periods in 2 to 3%.
* **Time effects** (`panel_control(time_effects = TRUE)`): every variable
  is standardized within its unit before the period means are removed. On
  the raw scale the period means of units that differ in size are those of
  the largest units, and their removal put one common series into the
  small units.
* **`placebo_test()`**: the circular offsets of the series are spaced so
  that every pair of series is displaced by at least `min_shift`.
  Independent offsets left small relative offsets between two series, at
  which a link survives in the placebo data at a neighbouring lag.
* **`"dfe"`**: the default classifier is trained under a fixed kind of
  random number generator as well as a fixed seed. Workers of the package
  run under L'Ecuyer-CMRG, so that a first call from inside
  `benchmark_methods()` or with `cores > 1` trained another forest than a
  first call from the session.

## New diagnostics and options

* **Heteroskedastic errors** (`"granger"`, `"toda_yamamoto"` on time
  series): the regressions are tested for heteroskedasticity (White's test
  in Koenker's form); below p = 0.001 a warning says that the F test can
  reject without a linear link (a nonlinear link whose linear projection
  is zero gave 27% rejections) and names
  `method_args = list(vcov = "HC3")`, the Wald test with the
  heteroskedasticity-consistent covariance HC3 (5%). `$extra` holds the
  p-value of the test and the covariance used.
* **Bayes factors beyond the range of the surrogates**: `surrogate_bf()`
  returns `bound`, the Bayes factor at the most extreme surrogate, which is
  all that the surrogates support; `$bayes` has the column `bf_bound`, and
  the print methods show such Bayes factors as "> bound". The decisions
  are unchanged.
* **Persistence**: `"pid"` and the panel tests of the Granger family warn
  for strongly persistent series (panels: median first-order
  autocorrelation of the units, corrected for their length, above 0.9).
  The warning for short and strongly autocorrelated series is not raised
  for series without autocorrelation.
* **Resolution**: with a step-up correction (`"BH"`, `"BY"`, `"hochberg"`,
  `"hommel"`) a warning states how many tests must attain the smallest
  p-value together for a link to be declared. The advisable minimal
  lengths of `"te"`, `"rf_granger"` (50) and `"ccm"` (55) are those at
  which the default test of two directions can reject at all.
* **Printing**: dependence without direction is shown as `<~>` (it was
  `<-L->`), and the column of the classes is called `dependent`.

## Input checks

* `lag` below 1 or not an integer, `cores` below 1 and series of fewer than
  8 time points are errors; so are exactly collinear series. `time`
  without `id` is an error (a time-series call takes the rows in their
  order), and a column that looks like a time index is reported. Gaps in a
  `Date` index of weekly or lower frequency are detected.
* Arguments that hold the data are evaluated before the seed of the
  analysis is set: `rcausality(rnorm(100), rnorm(100), seed = 1)` analysed
  the same data in every call.
* `bayes = TRUE` with `"toda_yamamoto"` and a sampled prior says that the
  BIC form is used; `bayes = FALSE` with `"bayes_granger"` says that it is
  ignored.

## Other

* `FNN` is no longer suggested (it was not used); `plm` and `sandwich` are
  suggested for tests against reference implementations.

# Rcausality 0.0.4

Corrections and additions that came out of the review of the accompanying
paper. Several of them change results.

## Reference distributions

* **Circular-shift surrogates** (`"te"`, `"ccm"`, `"rf_granger"`, `"ccdh"`,
  `"hirata"`) are drawn *without replacement*, and all admissible shifts
  are used when more replicates are requested than there are shifts. A
  series of length n has only about 0.8 n distinct shifts; drawing them
  with replacement returned ties among the replicates and p-values below
  1 / (number of shifts + 1), which are not valid. With many tests on short
  series this inflated the family-wise error rate (five white-noise series
  of length 100, `"te"`: 25% instead of 5%).
* The automatic number of replicates is no longer capped at 999 (with 50 or
  more tests no link could be declared at the 5% level, silently). When
  the smallest attainable p-value cannot pass a Holm or Bonferroni
  correction over the tests of an analysis, a warning says so and names the
  required series length.
* **Lag scans**: for `"te"`, `"cmi"` and `"pid"` every replicate is now
  computed for all candidate lags from one surrogate (one shift, one
  permutation), so that its maximum over the lags is a replicate of the
  observed maximum, as documented. (`"df"` keeps independent replicates per
  lag, which is conservative.)
* **`"cmi"`, `"cmi_lag"`, `"markov"`**: the reference distribution of the
  kNN conditional mutual information comes from surrogates that respect
  the conditioning *and* keep the serial dependence of the cause (the lags
  of the cause are split into their regression on the conditioning vector
  and residuals; the residuals are shifted circularly). The local
  permutations of Runge (2018), used so far, exchange time points and do
  not keep it: two independent first-order autoregressions with
  coefficient 0.9 were declared linked in 17% (100 time points) and 10%
  (300 time points) of the simulated data sets at the 5% level, against 9%
  and 5% now. `method_args = list(reference = "local")` restores the local
  permutations. With the new reference the replicates of a scan are joint
  for `"cmi_lag"` as well, and the default scan of the lags is limited as
  for `"te"` (see below).
* **Short and strongly autocorrelated series**: neighbouring circular
  shifts of such series give dependent replicates, and a test at level
  alpha needs 1/alpha independent ones at the very least. `rcausality()`
  estimates their number from the first-order autocorrelations of the pair
  and warns below 2/alpha (two independent autoregressions with coefficient
  0.9: 6 to 10% rejections at the 5% level with about 8 independent shifts,
  5 to 8% with about 25, 4 to 6% with about 84).
* **`"te"`**: the circular shifts of a scan over L lags are spaced by L. The
  windows of lags of neighbouring shifts overlap, so that the largest value
  over all shifts entered up to L replicates; with 8 lags and 80 shifts a
  true null hypothesis was rejected in 9% of the cases at the 5% level. The
  default scan is limited to the lags for which enough such shifts remain:
  at least 40, and as many as the smallest p-value needs to pass a Holm or
  Bonferroni correction over the tests of the analysis (a message reports
  the limit). A target analysis corrects within its families of tests (the
  links to the target, from the target, among the drivers), and the scan of
  a pair is limited by the size of its own family: with four drivers and
  300 time points the links to the target are scanned over the lags 1 to 3
  and the twelve links among the drivers at lag 1. The same holds for
  `"cmi"` and `"cmi_lag"` with their new reference distribution.
* **`"diks_panchenko"`**: with controls or more than one lag the p-value
  comes from surrogates that respect the conditioning (the lags of the
  cause are split into their regression on the conditioning vector and
  residuals, and the residuals are shifted circularly). The normal
  reference distribution is justified for the bivariate first-order test
  only; in the conditional test a cause that was correlated with a control
  was declared in 16 to 89% of the simulated data sets without a link.
  `diks_panchenko()` has the arguments `n_perm` and `seed` for `lag > 1`.
* **`"pid"` / `marpid()`**: the permutation test of the novel information of
  a source permutes it *within the levels of the other source*, the
  conditional null hypothesis; both sources were permuted freely, which
  tests their joint independence of the target.
* **`"rf_granger"`**: the surrogates respect the conditioning. A circular
  shift of the cause makes it independent of the controls and of the own
  past of the effect, on which it may depend under the null hypothesis; a
  cause that was correlated with a control was then declared in more than
  half of the simulated data sets without a link. The lags of the cause
  are now split into their regression on the conditioning features (linear
  plus a forest on the residuals) and residuals, and only the residuals are
  shifted. The forests of the test use all features at every split
  (`mtry = "all"`, the new default): with a random subset of a handful of
  lags, a redundant feature improves a forest and counted as a gain.

## Bayesian variants

* **`prior = "g"` and `"hyper_g"`**: the Bayes factor of the Granger and
  Toda-Yamamoto regressions takes the restricted regression as the base
  model (flat priors on the coefficients of the own lags and the controls,
  the g prior on the lags of the cause after the common regressors are
  projected out). It is then a function of the F statistic and its degrees
  of freedom. The former ratio of two Bayes factors against the
  intercept-only model depended on how well the own past predicts: for a
  persistent effect series it gave evidence against a link whatever the
  data (chickens -> eggs under the hyper-g prior: 0.0005, now 0.35).
* **Posterior probability under the spike-and-slab prior**: `post_prob` is
  computed from the Bayes factor and `prior_prob`, as for every other
  prior. It was the posterior probability that some lag is included, whose
  prior value 1 - (1 - inclusion_prob)^lag grows with the order (0.996 at
  order 8), so that `decision = "posterior"` declared links between
  independent series. `bayes_granger()` reports that probability
  separately (`any_included`).
* **`surrogate_bf()`**: outside the range of the surrogate statistics the
  kernel density of the null distribution is continued by a normal tail
  (element `extrapolated`). The tails of the kernel estimate decay with the
  bandwidth and gave Bayes factors above the bound of the model (3e7 for a
  statistic four standard deviations above its null mean; now about 1e3).
  New `effect_scale = "unit"` (also in `bayes_control()`): a
  unit-information scale of the effect prior, under which the evidence for
  the absence of a link grows with the length of the series.
* `bf_evidence()`: symmetric categories ("decisive against" below 1/100);
  the labels are those of Jeffreys in the wording of Lee and Wagenmakers.

* **Normal-slab and spike-and-slab priors** (`prior = "normal"`,
  `"spike_slab"`, the default of `bayes = TRUE`): the coefficients that the
  restricted and the unrestricted regression share (intercept, own lags,
  lags of the controls) have flat priors and are projected out, as for the
  g priors, so that the marginal likelihoods have `n - k0` degrees of
  freedom. Before, they had proper "vague" normal priors, which left all
  `n` observations in the exponent of the marginal likelihood: in a
  regression of 9 coefficients on 13 observations the Bayes factor exceeded
  3 for 29% of simulated series without a link (now 6%), and in regressions
  with nearly as many coefficients as observations for 80%. When the
  always-included coefficients leave no degrees of freedom the two priors
  stop with a message that names `prior = "minnesota"`, which shrinks all
  coefficients and remains defined. Bayes factors of long series change
  little (chickens and eggs: 2.3 to 1.8 for the default prior).
* **Bayes factor next to frequentist tests**: without `bayes = TRUE` the
  regression tests are accompanied by the exact Bayes factor of the
  unit-information prior (the g prior with g = number of usable
  observations) instead of its BIC approximation.
* **Panel units**: the Bayes factors of the unit regressions that are
  reported with a frequentist panel analysis, and those of
  `prior = "bic"`, are the exact Bayes factors of the unit-information
  prior (the g prior with g = number of usable periods) instead of its BIC
  approximation. With few residual degrees of freedom the approximation
  overstates the evidence (16 periods and 11 coefficients: Bayes factors in
  the thousands next to p = 0.25), which distorted the posterior share of
  units with a link (`prevalence`) and the group-level Bayes factor of the
  frequentist output. For time series `bayes = TRUE, prior = "bic"` warns
  when a regression has fewer than 10 residual degrees of freedom.

* **Strongly autocorrelated series and the F test**: `method = "granger"`
  warns when both series of a pair have a first-order autocorrelation
  above 0.8 and names `"toda_yamamoto"`. For nearly integrated series the
  F test rejects too often in samples of any practical length (two
  independent first-order autoregressions with coefficient 0.9: 10% at the
  5% level with 50 observations, 7% with 200; coefficient 0.98 and 200
  observations: 12%).

* **Strongly autocorrelated series and the kernel test**:
  `method = "diks_panchenko"` with its normal reference (two series, one
  lag) warns in the same situation and names
  `method_args = list(filter = "var")`. The long-run variance of the
  statistic is estimated from `n^(1/4)` autocovariances, too few for such
  series (coefficient 0.9: 6 to 8% rejections at the 5% level for 100 to
  1000 observations; 1.4 to 3.3% for the residuals of a vector
  autoregression, to which the test is usually applied). The filtered test
  detects nonlinear dependence only: the filter removes a linear link.

* **`"common_cause"` in stochastic systems**: the print method, the plots,
  the README and the vignettes now describe the class of `"cmi_lag"` and
  `"markov"` as what it is in data, dependence within the sampling
  interval (a common shock, an instantaneous or a fast link). A hidden
  driver *with memory* predicts both series from the past and is reported
  as a directed or bidirectional link by every method that uses the two
  series alone.

* **`"dfe"`** warns when the series are not binary. Its classifier is
  trained on binary Markov chains of short memory; continuous
  autocorrelated series that are cut at their mean have long memory, and
  92 of 100 pairs of independent first-order autoregressions were
  classified as `common_cause`.

## Input checking and reproducibility

* The default classifier of `"dfe"` is trained with a fixed seed and is
  therefore the same in every session. Before, it was trained from the
  state of the random number generator at the first call of a session, so
  that results with the default model were not reproducible across
  sessions (about one verdict in ten changed for the continuous series of
  a simulation that was recomputed). The random number stream of the
  session is not affected by the training.

* Missing values: a message reports how many values were interpolated in
  which series; a warning when more than 10% of a series is interpolated.
* Warnings for series shorter than the advisable minimal length of the
  method, for a requested order that had to be reduced, for constant and
  exactly collinear series; in panels for gaps in `time`, and a message when
  a column that looks like a time index is analysed as a variable.
* `predict()`, `forecast_causal()` and `forecast_compare()`: `seed` also
  seeds the forecast models (the forests of `model = "rf"` and the
  simulated paths); `predict()` has a `seed` argument.
* `benchmark_methods()`: every (scenario, replication, method) has its own
  seed, so that the result of a method does not depend on the other methods
  of the call.
* `print(bayes_granger())` states the direction that is tested.

# Rcausality 0.0.3

## One interface: method, Bayesian variant and cores as options

* `rcausality(x, y, method, bayes, cores, ...)`: the causality method, its
  Bayesian variant and the number of CPU cores are the first options of
  every analysis, for time series and for panel data alike.
* `method_args`: the settings of a method in one named list (the settings
  of the 0.0.1 methods are still accepted directly, e.g. `bins = 6`);
  unknown settings are an error. `causal_methods()` lists the methods (type,
  conditional and pooled estimators, advisable series length, reference)
  and `causal_methods("te")` the settings of one method.
* `lag = "first"` (Granger family): the first significant lag order, with a
  Bonferroni correction over the scanned orders; the Bayes factor is that of
  "a link at one of the orders" (the average over the orders).

## Panel data: every method, with Bayesian variants

* `rcausality(x, id = , time = )` / `panel_causality()`: every method for a
  long data frame of units observed over time, at three levels --
  **units** (the method in every unit that is long enough, with
  Benjamini-Hochberg or Bayesian false-discovery-rate control over the unit
  x link cells), **group** (a test of "no unit has the link" from the unit
  results: the mean of the standardised unit statistics against the
  convolution of their null distributions; Dumitrescu-Hurlin, Stouffer,
  Fisher, Simes and count tests as alternatives) and **pooled** (one
  estimate from all units). The panel decision uses the group level, the
  pooled level or both (`panel_control(decision = )`).
* Granger family: the pooled test is the heterogeneous-panel test of Juodis,
  Karavias and Sarafidis (2021) with the half-panel jackknife of Dhaene and
  Jochmans (2015), which removes the bias that makes fixed-effects tests of
  the reverse direction of a true link over-reject. Unit and pooled
  statistics are referred to a bootstrap under the null hypothesis (with the
  small-sample bias of the persistence corrected) through an empirical null
  (Efron 2004), because the F distribution is anti-conservative in short
  panels of persistent series. With `lag = "auto"` the order of every link is
  chosen by the BIC of its pooled regression.
* Model-free methods: pooled estimators from the stacked, within-unit
  standardised delay vectors -- random forest, Diks-Panchenko, kNN
  conditional mutual information (`"cmi"`, `"cmi_lag"`), transfer entropy and
  multispatial convergent cross mapping (Clark et al. 2015) -- with
  surrogate panels (the cause delayed circularly within the units). They
  work with units that are far too short for a unit-by-unit analysis.
* Bayesian variants of all of them but the random forest, which is tested
  with surrogates only (`bayes = TRUE`): unit Bayes factors, a
  prevalence model (posterior share of the units with the link, partially
  pooled unit probabilities), the Bayes factor of the group statistic, and
  pooled Bayes factors (BIC, g, hyper-g, normal / spike-and-slab and
  Minnesota priors for the Granger family; surrogate Bayes factors for the
  model-free methods).
* `plot()` of a panel analysis: heat map of the links by unit with the
  panel-level tests (after Kosztyan, Matrai and Kiss 2026), unit evidence,
  prevalence, network; `summary()`, `as.matrix(by = "units")`.
* `panel_control()` (levels, decision, group test, unit-level correction,
  prevalence prior, common time effects), `panel_long()` (matrices and
  arrays to long format), `sim_panel()` (heterogeneous prevalence, nonlinear
  coupling, unbalanced panels, a five-variable `"drivers"` scenario).
* `panel_granger()` is kept as `panel_causality(method = "granger")`; its
  result is now a `causpanel` object (the 0.0.2 arguments `methods` and
  `n_placebo` are replaced by `panel_control()` and `placebo_test()`).

## The causal drivers of one explained variable

* `rcausality(x, target = )` / `causal_drivers()`: which series drive the
  target, at which lag, with how much evidence (importance), directly or
  through another driver (conditional re-tests; the mediators are named),
  with feedback from the target or not -- and the causal links among the
  explanatory variables. Three families of links are corrected for
  multiplicity separately. Any method, Bayesian variants, time series and
  panels.
* `plot()`: the lag ladder with the explanatory graph, the importance
  chart, the causal network; for panels the drivers-by-unit heat map.
  `sim_drivers()` simulates a system with direct and indirect drivers.

## Robustness, comparison and forecasting for panels and targets

* `placebo_test()` and `robustness()` accept panel and target analyses
  (shifts and time reversal within the units; random subsets of the units;
  the decision of each panel level).
* `compare_methods(..., id = )` compares methods on panel data; a method that
  cannot be run on the data is left out with a warning.
* `predict()` of a panel analysis: forecast equations pooled over the units
  (`"var"`, `"bvar"`, `"rf"`) with a hold-out evaluation by unit;
  `forecast_causal()` and `forecast_compare()` accept `id` (Diebold-Mariano
  test on the loss differentials averaged over the units).

## Corrections

* `bayes_control(prior = "g", scale = )` and `bayes_control(prior =
  "hyper_g", hyper_a = )`: the scale of the g-prior (`g = scale * n`) and the
  parameter `a` of the hyper-g prior were documented but not used -- the
  closed-form Bayes factors of the Granger and Toda-Yamamoto regressions
  were always computed with `g = n` and `a = 3` (and the pooled panel Bayes
  factor with `g = n`). The settings now take effect in `rcausality()`,
  `bayes_granger()` and the panel analysis, and `print` shows the values in
  use. Results with the default settings are unchanged.
* `rcausality(method = "toda_yamamoto", bayes = TRUE)` with a prior that has
  no closed form for the lag-augmented regressions (the default
  spike-and-slab, the normal slab, Minnesota) uses the BIC Bayes factor, as
  documented; the result now records this prior (`$settings$prior` was the
  requested one).

* Random forests (`"rf_granger"`, `rf_granger()`, the pooled panel forest,
  `predict(model = "rf")`, `forecast_compare(model = "rf")`): the minimal
  node size of the trees is now selected by the out-of-bag error
  (`min_node_size = "auto"`, the new default; grid: 5 and 2%, 5%, 10%, 20% of
  the observations) instead of being fixed at 5. With only a few lags as
  features, fully grown trees interpolate the noise: the forest on the own
  past of an AR(1) series had an out-of-bag R-squared near zero, and *any*
  added series lowered its error by about 10% (the random choice of the
  split variable acts as a regulariser). The surrogate test was valid all
  the same, because its surrogates share the artefact, but the reported
  out-of-bag and out-of-sample gains were inflated, and the forecast
  comparison with random forests credited unrelated series with a gain.
  In the test the node size is chosen on the own-past forest (without
  reference to the cause) and used for all forests of the out-of-bag
  statistic and its surrogates; in the rolling-origin evaluation and in the
  forecast models every forest is tuned on its own training data. The
  selected size is
  reported (`$tests$node_size`, `$extra$node_size`); `min_node_size = 5`
  restores the previous behaviour.
* `bayes_control(case_prior = )` now *replaces* the prior of the method over
  the five causal cases. It was multiplied onto the posterior probabilities,
  which is correct only when the method's own prior is uniform (the
  df-entropy classifier); for dimensional causality (uniform prior over the
  nine elementary events of Benko et al., i.e. weights 1 : 1 : 1 : 3 : 3)
  and for the Bayesian Markov table (2 : 2 : 2 : 1 : 1) the requested prior
  was therefore combined with the method's prior, and the Bayes factor of
  the Markov directions used the prior probability 1/2 whatever the
  requested prior. Results with `case_prior = NULL` are unchanged.

## Output

* The explanatory lines of all `print` methods are wrapped to the width of
  the console (`getOption("width")`), so that the output fits narrow
  consoles and printed pages; `options(Rcausality.guide = FALSE)` hides the
  reading guide that `print` appends to a result.
* Heat map of a panel analysis: the unit labels are scaled to the width of
  the cells and no longer overlap when there are many units.

## Documentation

* New vignette `vignette("panel")`; the package cites Kosztyan, Matrai and
  Kiss (2026, *The World Economy*, 49(11), 2478-2525,
  doi:10.1111/twec.70148), whose panel and explanatory Granger analyses are
  generalised here.
* Corrections in the vignettes: the power of the Minnesota Bayes factor for
  short series is described as the stored simulation shows it (equal to the
  F-test for the shortest series, higher from about 30 observations), and
  the dimension relation of a directed link in the dimensional-causality
  section reads `D_X < D_Y = D_J`.

# Rcausality 0.0.2

## Generalisations of the Granger test

* `rcausality(method = "toda_yamamoto")`: lag-augmented Granger test of
  Toda and Yamamoto (1995) for integrated or cointegrated series (`d_max`
  extra lags, not tested; BIC / g / hyper-g Bayes factors with
  `bayes = TRUE`).
* `rcausality(method = "diks_panchenko")` and `diks_panchenko()`: the
  nonparametric (nonlinear) Granger test of Diks and Panchenko (2006), with
  the bandwidth rule of the paper, a C++ kernel (O(n^2) time, O(n) memory),
  optional VAR prefiltering (only the nonlinear part tested) and a Bayesian
  variant through `surrogate_bf()`.
* `granger_spectral()`: Granger causality in the frequency domain -- Geweke's
  (1982) spectral measure and the Breitung-Candelon (2006) pointwise tests,
  with the frequency bands of significant causality and a `plot()` method.

## Machine-learning Granger-like causality

* `rcausality(method = "rf_granger")` and `rf_granger()`: random-forest
  Granger-like test after Kosztyan, Matrai and Kiss (2026): out-of-bag gain
  against circular-shift surrogates of the cause, rolling-origin
  out-of-sample gain, permutation importance of the cause's lags and the
  shadow-feature rule; reproducible in parallel (forests seeded from the
  worker's random stream).

## Bayesian analysis of short series

* `bayes_control(prior = "minnesota")`: hierarchical Minnesota prior for the
  Bayesian Granger regression (Litterman 1986), the overall tightness
  integrated out under the hyperprior of Giannone, Lenza and Primiceri
  (2015); closed-form Bayes factor, exact independent posterior draws of the
  coefficients, the error scale, the tightness and Geweke's measure. Works
  when the F-test has no residual degrees of freedom; a requested lag order
  is used as it is.

## Panels

* `panel_granger()`: Granger causality in panels of short series --
  unit-level tests with Benjamini-Hochberg control, the Dumitrescu-Hurlin
  (2012) test with exact finite-T moments, the pooled fixed-effects test, the
  pooled Bayes factor (Minnesota prior), the pooled random-forest test, and
  placebo distributions from within-unit circular shifts; `sim_panel()`.

## Robustness, placebo and triangulation

* `placebo_test()`: placebo distribution of the number of links and of the
  evidence of every link (independent circular shifts of the series), and
  the time-reversed Granger score (Haufe et al. 2013; Winkler et al. 2016).
* `robustness()`: the links under other lag orders, multiplicity
  corrections (Holm, Benjamini-Hochberg, none), contiguous subsamples
  (stability) and other methods; tile plot.
* `compare_methods()` runs the methods in parallel (`cores`).

## Forecasting with any causality method

* `predict.causres(model = c("var", "bvar", "rf"))`: forecast models built
  from the detected links -- the restricted VAR, the Minnesota-prior VAR
  (posterior-predictive intervals) and random forests (recursive simulation
  with bootstrapped out-of-bag residuals); `forecast_causal(method, model)`.
* `forecast_compare()`: the out-of-sample value of the links of several
  causality methods against the own-lags model (links chosen on the training
  part only), with Diebold-Mariano tests (Harvey-Leybourne-Newbold
  correction) and a `plot()` method.

## Code quality

* A method registry (labels, statistic names, families) replaces the label
  tables repeated in `print()` and `plot()`; the new methods share one
  pairwise driver and one `pairs` assembly (`.finish_pairs()`); settings of
  the new methods go through `method_args`, whose names are checked; argument
  checks of the numeric settings of `bayes_control()` and of the new
  functions.
* A new vignette, "Beyond linear Granger causality", with precomputed
  long-running results; package-level help page `?Rcausality`.

# Rcausality 0.0.1

* Initial skeleton (M0): `causres` class with `print`/`summary`/`plot`,
  Granger and Bayesian (BIC Bayes factor) Granger tests, bivariate, all-pairs
  and MIMO block modes, simulators, preprocessing helpers.
* `method = "cmi_lag"`: lag-shifted CMI scheme of Zlatniczki et al. (2021)
  with contemporaneous common-driver test (`common_cause` class);
  `cmi_convergence()`, `lag_triangle_check()`.
* `method = "markov"`: Markov-chain approach of Stippinger et al. (2023)
  (`markov_test()`, `select_embedding()`, decision tables (10)-(12)); exact
  contingency-table test for discrete states; `sim_turnwheel()`,
  `sim_logistic3()`.
* `method = "df"`: degrees-of-freedom / conditional-variance method of
  Telcs, Kurbucz & Jakovac (2025): `df_variance()` with `print`/`plot`
  (sigma-versus-window-width profile), orders of x, y and the joint process
  with the Table-2 reading, directed variance-reduction tests;
  `sim_df_rotation()`.
* Multi-core execution: `cores` argument of `rcausality()`,
  `local_permutation_test()`, `markov_test()`, `df_variance()`
  (reproducible with `seed`; L'Ecuyer-CMRG streams).
* Plots: `plot(r, "matrix")` (causal matrix, cause rows -> effect columns),
  `"lags"` (statistic per candidate lag), `"embedding"`, `"all"`;
  `plot_embedding()` (2-D/3-D time-delay embedding); `causres$by_lag`.
* Reproduction scripts in `inst/reproduce/` and a README section per paper.
* References with DOIs in all help pages; `stats::` qualified calls.
* C++ (Rcpp) back-end for the max-norm neighbour searches: `knn_cmi()`,
  `local_permutation_test()`, the Markov test and the df kernels are 6-10x
  faster; the pure-R estimator is kept as the reference in the tests.
* Markov test for continuous states rewritten after Algorithms 1-2 of
  Stippinger et al. (2023): local-permutation surrogates of the futures,
  2-bin local tables, aggregation by a clipped geometric mean (`aggregate`,
  `clip`); the coupled logistic-map pattern of the paper's Figure 3 is
  reproduced; `sim_logistic3(obs_noise = )`.
* `method = "dc"`: Bayesian dimensional causality of Benko et al. (2024)
  (`dimensional_causality()` with `print`/`plot`, `sim_logistic_bdc()`);
  independent implementation from the paper, posterior probabilities of the
  five causal cases averaged over the neighbourhood sizes.
* Author/maintainer: Zsolt T. Kosztyan.
* `method = "ccdh"`: complexity-based causal discovery of a hidden common
  cause (Stippinger et al.): NSRPS compression in C++, rank-order coding,
  skewed joints, circular-shift surrogates and a reconstruction index
  (`ccdh()` with `print`/`plot`).
* `method = "dfe"`: degrees-of-freedom method with conditional-entropy
  features and random-forest classification (Jakovac et al.):
  `sim_markov_pair()`, `dfe_train()`, `dfe_classify()` with `print`/`plot`;
  binary conditional entropies and the coupled-chain simulator in C++;
  `ranger` in Suggests.
* README: status table per milestone; PRE reference of the df method updated
  (112(1), 014204).
* Baselines: `method = "te"` (binned transfer entropy with circular-shift
  surrogates, `transfer_entropy()`) and `method = "ccm"` (convergent cross
  mapping with simplex projection in C++, convergence and surrogate tests,
  `ccm()`).
* `compare_methods()` (verdicts of several methods per pair, majority
  consensus, tile plot) and `benchmark_methods()` (accuracy of methods on
  simulated scenarios with known truth, bar plot).
* `method = "pid"` and `marpid()`: finite-resolution MAR-PID (Telcs 2026;
  Blackwell-order decomposition of Mages et al. 2024) for two sources --
  recursive quantile binarisation, zonogon join, nonnegative redundant /
  unique / synergistic atoms, surrogate tests.
* `method = "hirata"` / `hirata_rp()`: recurrence-plot criteria of Hirata &
  Aihara (2010) (coverage, simultaneous recurrences, surrogates);
  `method = "cd"` / `krakovska_cd()` and `corr_dim()`: correlation-dimension
  relations of Krakovska (2019) with bootstrap tests.
* `predict.causres()` / `forecast_causal()`: causality-informed forecasts
  (restricted VAR on the detected links, hold-out comparison with the
  own-lags model, h-step forecasts with intervals); class `causfc` with
  `print`, `summary`, `plot`.
* Optional PCMCI / PCMCI+ wrapper `pcmci()` / `method = "pcmci"` around the
  Python package tigramite (Runge et al. 2019) through reticulate (Suggests);
  `pcmci_available()`; results are mapped onto `causres`.
* `sim_var(n_series > 2)`: the `"adjacency"` attribute now carries the
  correct `cause` x `effect` dimnames.
* README: development notes for checking inside a cloud-synced folder.
* Bayesian variants: `prior = c("bic", "g", "hyper_g")` for
  `method = "bayes_granger"` (exact g-prior and hyper-g Bayes factors);
  `bayes = TRUE` calibrates the p-values of every model-free method into
  Bayes factors / posterior link probabilities (`bf_from_p()`, Sellke,
  Bayarri & Berger 2001) and decides by `bf_threshold`.
* `method = "df"`: a determinate Table-2 reading of the orders now decides
  the pair (the directed variance-reduction tests decide undetermined
  pairs and are always reported); `$df$decided_by`.
* `ccdh()`: rank-order coding vectorised (about 100x faster).
* Version 0.0.1; vignettes "Getting started with Rcausality" and "The
  methods, paper by paper"; `inst/CITATION`; URL/BugReports fields.
* Plots: `plot.hirata()` (recurrence plots and the coverage picture),
  `plot_series()` (stacked panels of the series), `plot.causfc()` shows
  the hold-out one-step predictions of both models and intervals from the
  multi-step MSE recursion of the restricted VAR; benchmark and
  dimensional-causality plots tidied.
* `ccm()`: the full-library skill is computed once per surrogate (20x
  faster surrogate tests).
* Vignettes rewritten: full references with DOIs per method, the idea and
  the location of the reproduced result in every paper, figures for every
  method, comparative analyses (`compare_methods()`, `benchmark_methods()`
  on linear and chaotic benchmarks); sidebar table of contents, vector
  (SVG) figures.
* Examples: the quick part runs, the longer illustrations are wrapped in
  `\donttest{}`.

## Fixes after the 0.0.1 code review (October 2026)

* kNN-CMI tests (`cmi`, `cmi_lag`, `local_permutation_test()`) standardise
  the series first: the max-norm neighbour search was not invariant to
  rescaling one variable alone (a `1000 * x` could flip the verdict).
* `history = "auto"` (new default): the effect's conditioning history is
  the VAR order of the pair (BIC); with the old default `history = 1` an
  effect with memory > 1 produced spurious reverse links in `cmi`,
  `cmi_lag`, `te`, `df` and `pid`. `pid` collapses a deeper own past to its
  linear predictor (MAR-PID takes two scalar sources).
* Lag scans of `cmi`, `cmi_lag`, `te`, `df` and `pid` honour `max_lag`
  (the silent cap at 3 is gone) and use a max-statistic permutation test
  over the lags instead of a Bonferroni factor (no extra permutations are
  needed for more lags); `te` scans single lags of the cause.
* `df`: the Table-2 reading handles series of order 0 (white noise is not
  driven; whether it drives is read from the other series' accepted
  conditions), and a direction claimed by the table is kept only when its
  directed test is nominally significant, so the uncontrolled order tests
  can no longer override non-significant directed tests.
* `transfer_entropy()`: discrete series keep their values as symbols and
  continuous ones are binned by average ranks (quantile bins collapsed on
  0/1 data); `l`, `k` must be >= 1; surrogate-shift range checked.
* `marpid()`: the Blackwell join is the lower convex envelope of the
  union of the two zonogon chains (the convex hull of the union of the
  zonogons); the previous upper envelope collapsed the join to the diagonal
  and made one unique atom always zero (min-MI PID). The Z-channel and its
  mirror now give unique 0.19 / 0.19, redundancy 0.12.
* df-entropy method: short series or NA input give an error instead of a
  segfault (bounds checks in `cond_entropy_bits_cpp()` and
  `cond_entropy_sym_cpp()`); `dfe_train(lengths = <one value>)` no longer
  samples from `1:lengths`.
* `directions = "forward"` is honoured by all methods (`markov`, `hirata`,
  `cd`, `dc`, `ccdh`, `dfe` used to report both directions and adjust over
  unrequested tests).
* `predict.causres()`: the hold-out evaluation re-selects the links on the
  training part (`refit = TRUE`, new default) instead of using links chosen
  on the full sample; constant series (aliased coefficients) no longer
  produce NA forecasts; `holdout = 0` skips the evaluation explicitly.
* `benchmark_methods()`: the truth is oriented by the simulator's columns
  (non-alphabetical names gave accuracy 0); failed runs count as wrong and
  raise a warning with the first error (`$runs$error`).
* Pairwise drivers (`ccdh`, `dfe`, `dc`, `hirata`, `cd`, `markov`) order the
  unordered pairs by column position, so a `lag` applied to the second
  series no longer depends on the names; `ccdh`'s skewed joint uses the
  column orientation of `embed_delay()` (so `ccdh(x, x)` gives exactly 0
  with ties).
* `dimensional_causality()`: the Farahmand estimator uses `log(k / k2)`
  (odd `k` gave a saw-tooth bias); ties are broken by position and a ridge /
  pseudo-inverse guards the dimension covariance (discrete data no longer
  stop with "computationally singular").
* `sim_var()` (two series) rescales the coefficients to keep the process
  stationary and warns; the simulators and the other seeded functions
  restore the caller's RNG state (`.local_seed()`), `rcausality(seed = )`
  as well.
* `pcmci()`: contemporaneous rows are matched whatever their order (the
  `common_cause` flag could be lost with non-alphabetical names).
* `hirata_rp()`: the diagonal is never a recurrence (`theiler = 0` gave
  index 0.28 on independent data); thresholds of long series are estimated
  from a subset of the points (O(n^2) memory).
* Non-numeric (date, factor, character) columns give an error; the
  function-input message now fires for `y` as well.
* `ordinal_pattern()` returns canonical (Lehmer) codes comparable across
  series; rank-order coding limited to m <= 12 (no integer overflow).
* `ccm()`: library sizes sorted, `rho_max` at the largest library,
  surrogate-shift range checked. `plot(r, "bars")` / `"embedding"` work for
  posterior- and skill-based results (`dc`, `ccm`). `compare_methods()`
  tolerates all-NA verdicts.
* Explicit `theta` / `k` are no longer overridden (`NULL` defaults with
  method-specific values). Documentation of `markov_test()` (aggregation
  default) and `dfe_train()` (in-sample confusion) corrected. CCDH is cited
  as published: Stippinger et al. (2023), *Chaos, Solitons & Fractals* 176,
  114054. Vignettes fall back to PNG figures without cairo.
* New regression tests `tests/testthat/test-review-fixes.R`; `FNN` in
  Suggests for the KSG reference test.
* Plots made more self-explanatory: the network names the method, the
  decision rule and n, labels every significant edge with its lag and
  evidence (p-value, posterior or Bayes factor), draws the tested but
  non-significant pairs as faint dotted arrows (up to five series), colours
  the nodes by role (driver / driven / mediator / isolated) and lays two
  series out side by side; the causal matrix prints the p-value / posterior
  and the lag of every directed test; the lag scan colours only the
  significant pairs, labels every line at its end and falls back to a bar
  chart for a single lag; embeddings carry a time colour bar and a
  one-line interpretation; the comparison tiles name the direction with the
  series names, show the agreement and carry a legend. `print(<causres>)`
  ends with a short reading guide.
* `print(<causfc>)` shows the forecasts in a compact table (one column per
  series, point [interval]); vignettes use `svglite` figures when available
  (installed documentation under 1 MB).

## Bayesian variants of every method

* `rcausality(bayes = TRUE, control = bayes_control(...))` gives every
  method a Bayesian variant with explicit priors, Monte Carlo / MCMC sizes
  and a decision rule: `bayes_control()` collects `prior` (`spike_slab`,
  `normal`, `g`, `hyper_g`, `bic`), `effect_prior` (half-Cauchy, half-normal,
  Cauchy, normal), `scale`, `inclusion_prob`, `null_model` (kde / normal),
  `case_prior`, `n_iter`, `burnin`, `thin`, `chains`, `n_draws`, `hyper_a`,
  `credible`, `decision` (`bf` / `posterior`).
* `bayes_granger()`: Bayesian Granger regression with stochastic-search
  variable selection over the cause's lags (exact enumeration up to 12
  lags, Gibbs beyond), a normal-slab block comparison, posterior draws of
  the coefficients and of Geweke's measure with credible intervals and an
  effective sample size.
* `surrogate_bf()`: Bayes factor of an observed statistic against its
  surrogate / permutation draws with an effect-size prior (importance
  sampling; closed form for the normal null and half-normal prior); used by
  `cmi`, `cmi_lag`, `te`, `pid`, `df`, `ccdh`, `hirata`, `cd`, `ccm` (the
  surrogate draws are now kept in `$null`). The Markov-chain method combines
  the Bayes factors of its six tests in a Bayesian decision table; `dc` and
  `dfe` take a prior over the causal cases; `pcmci` uses a Savage-Dickey
  Bayes factor of the partial correlation. Results carry `$bayes` (per-link
  Bayes factor, posterior, effect interval, model) and `$bayes_fits`;
  `plot(r, "posterior")` shows them. The Sellke-Bayarri-Berger calibration
  (`bf_from_p()`) is no longer used by `bayes = TRUE` but stays available.
* README: the status section was removed (the package is complete; open
  research items are listed under "Open items").

## Parallel computation with identical results, examples and README

* Parallel execution was reworked so that **the results never depend on the
  number of cores**: every unit of work (a permutation or surrogate, a
  series pair, an MCMC chain, a bootstrap resample, a benchmark
  replication) draws from its own L'Ecuyer-CMRG stream derived from the
  caller's generator, and the session's generator (kind and state) is
  restored afterwards. `cores = 1` and `cores = 8` give the same numbers
  for the same `seed`, on every platform, and a script mixing Rcausality
  calls with other random draws stays reproducible when `cores` changes
  (previously serial and parallel runs differed). `?Rcausality-parallel`
  documents the design.
* `cores` was added to `transfer_entropy()`, `ccm()`, `marpid()`,
  `ccdh()`, `hirata_rp()`, `corr_dim()`, `krakovska_cd()`,
  `bayes_granger()` (chains in parallel) and `benchmark_methods()`
  (replications in parallel); inside `rcausality()` every method now
  parallelises the surrogates of a single pair when there are fewer pairs
  than cores. `cores` accepts an integer, `"auto"` or a cluster object from
  `parallel::makeCluster()`; its default is `getOption("Rcausality.cores",
  1L)`. On Windows a socket cluster is started for the call and stopped when
  it returns (no connection is left open).
* `bayes_granger()` runs its chains independently and reports the
  Gelman-Rubin `rhat` of Geweke's measure (also in `$bayes` of a Bayesian
  Granger run, with the effective sample size); `thin` is honoured by the
  Gibbs sampler; `settings$sampler` tells whether the lag subsets were
  enumerated or Gibbs-sampled.
* Examples: every help page explains what its example computes and what to
  look for; the parts that take more than a moment are wrapped in
  `\donttest{}` (still shown on the help page and run by
  `example(..., run.donttest = TRUE)`), so that `R CMD check` runs the
  examples in a few seconds.
* README: the long script was replaced by a sectioned walkthrough that
  says, for every method, what is examined, what the code does and what to
  look for in the output ("How to read a result", one section per method,
  Bayesian variants, comparison and benchmark, forecasting, plots, speed,
  own data).
* Regression tests of the parallel guarantees (`test-parallel.R`).
