# Rcausality

Causality analysis for time series in R -- a unified interface to Granger
causality, kNN conditional mutual information (plain and lag-shifted with
common-driver detection), the Markov-chain approach, the degrees-of-freedom
(conditional-variance) method, Bayesian dimensional causality, the
complexity-based CCDH method, the df-entropy / random-forest method,
finite-resolution MAR-PID, recurrence-plot and correlation-dimension
criteria, transfer entropy and convergent cross mapping, an optional PCMCI
wrapper, causality-informed forecasting, tools to compare and benchmark the
methods -- and a **Bayesian variant of every method** (`bayes = TRUE`, with
the priors, iterations and decision rules set in `bayes_control()`).
Every method from the group's papers has a reproduction of a published
result in `inst/reproduce/` and a section below ("Reproductions of
published results") stating what was and what was not reproduced.

> Author and maintainer: Zsolt T. Kosztyán (<kosztyan.zsolt@gtk.uni-pannon.hu>),
> University of Pannonia. Bug reports and questions: please open an issue on
> GitHub or write to the maintainer. License: GPL-3.

## Install

```r
# from GitHub (development version; needs a C++ compiler: Rtools on Windows,
# Xcode command-line tools on macOS)
# install.packages("remotes")
remotes::install_github("kzst/rcausality", build_vignettes = TRUE)
```

The neighbour searches of the kNN-based methods (`cmi`, `cmi_lag`, `markov`,
`df`, `dc`) and the NSRPS compression of `ccdh` run in C++ (Rcpp); `cores`
distributes pairs, permutations, chains or replications over several CPU
cores without changing the results (see "Speed" below).

## Quick start: two series

The simulator `sim_var()` produces data with a *known* truth, so the first
steps can be checked against it. Here `x` drives `y` with lag 1.

```r
library(Rcausality)
d <- sim_var(500, scenario = "x_to_y", seed = 1)
attr(d, "truth")                      # "x_to_y"

r <- rcausality(d$x, d$y)             # Granger causality, lag order chosen by BIC
r
```

**What you see.** The header names the method, the number of series and of
directed tests and the decision rule (adjusted p < 0.05, Holm). Then one row
per *directed test* (`x -> y` and `y -> x`): the lag order used, the F
statistic with its degrees of freedom, the raw and the adjusted p-value, the
BIC-approximated Bayes factor `BF10`, the posterior probability of the link
and the decision (`sig`, starred by strength). Below the tests comes the
*classification of the pair* -- one of `x_to_y`, `y_to_x`, `bidirectional`,
`independent` (and `common_cause` for the methods that can detect it) -- and
a short reading guide. For the simulated pair the `x -> y` row is significant
(F = 118, p < 1e-4, `BF10 > 1e4`), the `y -> x` row is not (p = 0.77,
`BF10 = 0.05`, i.e. evidence *against* a link), and the class is `x_to_y`.
`summary(r)` adds the sign of the net influence, the effect size (the gain
in R^2), the Kass-Raftery evidence label and the *node roles* (driver,
target, hub/mediator, isolated).

```r
summary(r)                            # verdict per pair with the evidence behind it
r$pairs                               # the directed tests as a data frame
r$classes                             # the pair classes as a data frame
plot(r)                               # network: an arrow per significant link, labelled with lag and p
```

> **Note.** The data frame is called `d`, not `df`: `df` is a base R function
> (the F density), and `rcausality(df)` fails with "Input is a function, not
> data".

## How to read a result

Every method fills the same `causres` object, so the output always has the
same shape:

* `$pairs` -- one row per directed test: `cause`, `effect`, `lag`,
  `statistic`, `p_value`, `p_adj` (after the multiple-testing correction
  `adjust`, Holm by default), `bf10`, `post_prob`, `evidence`,
  `effect_size`, `r2_gain`, `sign`, `n_eff` and `significant`.
* `$classes` -- one row per *unordered* pair with its class (the five classes
  above) and a `relation` string such as `x -> y`.
* `$instant` -- for the methods that test a contemporaneous / common-driver
  component (`cmi_lag`, `markov`, `df`, `ccdh`, `hirata`, `cd`, `dc`).
* `$by_lag` -- the statistic of every directed test at every candidate lag
  (the lag-scan plot), `$settings` -- every setting that was used, and
  method-specific tables (`$markov`, `$df`, `$dc`, `$ccdh`, `$pid`, `$ccm`,
  `$bayes`).

`as.matrix(r, "p_adj")` gives the cause x effect matrix of any column of
`$pairs`, `as.data.frame(r)` the tests, and `plot(r, type)` draws the
network, the p-value matrix, the lag scan, the embeddings or the posterior
panel (see "Plots" below).

## A walkthrough on a known network

### The data: a sparse VAR with a known adjacency matrix

`sim_var(n_series > 2)` draws a random sparse VAR network and attaches the
true adjacency matrix (rows = cause, columns = effect), so that every verdict
below can be checked.

```r
net <- sim_var(600, n_series = 4, seed = 42)
attr(net, "adjacency")
plot_series(net, n_last = 150)        # the last 150 observations, one panel per series
```

### Granger causality on all ordered pairs

`rcausality(net)` tests every ordered pair of the four series; `lag = "auto"`
picks the VAR order of each pair by BIC and the p-values are Holm-adjusted
over all tests.

```r
g <- rcausality(net, method = "granger")
summary(g)                            # five-class verdict per unordered pair
as.matrix(g, "p_adj")                 # cause (rows) x effect (columns): adjusted p-values
plot(g, "matrix")                     # the same matrix as a picture, starred = significant
```

**What to look for.** Compare the starred cells of the matrix with the
`TRUE` entries of `attr(net, "adjacency")`: the detected links should be
the true ones (here `v1 -> v2`, `v2 -> v3`, `v4 -> v3`, all with
p < 1e-4). A pair with links in both directions is `bidirectional`; a pair
without any is `independent`. `summary(g)` lists the significant links
strongest first and the role of every series (`v3` is the target of three
links, `v2` a hub/mediator). The lag printed for a link is the VAR order
the test used, not necessarily the lag at which the coupling acts;
`plot(g, "lags")` shows the F statistic of every pair over the candidate
lags.

### Decision by Bayes factors instead of p-values

`method = "bayes_granger"` decides by the Bayes factor `BF10` of the link
(Kass & Raftery labels: 1-3 anecdotal, 3-20 positive, 20-150 strong,
> 150 very strong); a link needs `BF10 > bf_threshold` (3 by default).
With `prior = "bic"` the factor is the closed-form BIC approximation of the
same regressions as the F-test; without `prior` the spike-and-slab
regression of `bayes_granger()` is sampled for every pair (see "Bayesian
variants" below). Unlike a p-value, `BF10 < 1` is evidence *for* the
absence of a link.

```r
b <- rcausality(net, method = "bayes_granger", prior = "bic", bf_threshold = 3)
b
plot(b, "heatmap", evidence = "bf")   # log10 BF10 instead of -log10 p
```

### Conditional (multivariate) Granger causality

Pairwise tests see an *indirect* link `a -> b -> c` as a link `a -> c`.
With `conditional = TRUE` every test conditions on all the other series, so
only direct links remain.

```r
rcausality(net, conditional = TRUE)$classes
```

**What to look for.** Links that disappear when conditioning are mediated by
a third series; the pairwise network and the conditional network together
tell which paths are direct. In the network above the true links are
`v1 -> v2`, `v2 -> v3` and `v4 -> v3`; the pairwise test also reports
`v1 -> v3` (weak, "anecdotal" evidence), which is the indirect path
`v1 -> v2 -> v3`, and `conditional = TRUE` removes it.

### Blocks of causes and effects (MIMO)

Two blocks test every series of the first block as a cause of every series of
the second; `directions = "forward"` restricts the tests to that direction.

```r
rcausality(net[, 1:2], net[, 3:4])                          # both directions, all cross pairs
rcausality(net[, 1:2], net[, 3:4], directions = "forward")  # only block 1 -> block 2
```

## Non-linear coupling and model-free methods

Granger causality is linear: a coupling that acts through a non-linear map
can be invisible to it. The model-free methods estimate information-theoretic
or geometric quantities instead and test them against permutations or
surrogates. They are slower (the kNN searches are O(n^2) per permutation, in
C++), and `cores` distributes the work (see "Speed" below); `seed` makes a
run reproducible.

### kNN conditional mutual information (`"cmi"`)

The directed statistic is `I(x[t-lag]; y[t] | y's own past)`, the
Frenzel-Pompe kNN estimator, with the local-permutation null of Runge (2018)
(x's past is permuted within neighbourhoods of the conditioning set). With
`lag = "auto"` every lag `1..max_lag` is tested at once and the p-value is a
max-statistic test over the lags, so no further correction is needed.

```r
lg <- sim_logistic(300, scenario = "x_to_y", seed = 1)    # coupled logistic maps, x drives y
rc <- rcausality(lg, method = "cmi", max_lag = 3, n_perm = 100, seed = 1)
rc
plot(rc, "lags")                                           # the CMI over the candidate lags
```

**What to look for.** The `x -> y` row has a clearly positive CMI (in
nats) with a small p-value and the reported lag is the one with the largest
CMI; the `y -> x` row is at the permutation level (a CMI around zero, p near
1). The linear F-test may also react to a logistic coupling, but its F
statistic measures the fit of a linear model, not the amount of dependence;
the CMI quantifies the information flow directly and its null preserves the
effect's own dynamics. `history = "auto"` conditions on the VAR order of the
pair (at most 5 lags of the effect's own past); too short a history lets the
cause's lags stand in for the effect's deeper past and produces spurious
reverse links.

### Hidden common driver: the lag-shifted CMI (`"cmi_lag"`)

Zlatniczki et al. (2021): if neither lagged direction is significant but the
*contemporaneous* dependence `I(x[t]; y[t] | both pasts)` is, the two series
share a common driver. `cmi_lag` adds this test (the `$instant` table) and
the class `common_cause`.

```r
# two AR(1) series with correlated innovations: no directed link, a shared driver
set.seed(5); n <- 300
e1 <- rnorm(n); e2 <- 0.8 * e1 + 0.6 * rnorm(n)
x <- as.numeric(stats::filter(e1, 0.5, "recursive"))
y <- as.numeric(stats::filter(e2, 0.5, "recursive"))
cl <- rcausality(x, y, method = "cmi_lag", seed = 1)
cl$instant                                                 # the contemporaneous test
cl$classes                                                 # common_cause
```

### The Markov-chain approach (`"markov"`)

Stippinger et al. (2023): a driver is a Markov chain on its own, a driven
series is not (its driver acts as hidden memory), and the lagged joint
`(x[t], y[t+1])` is Markov again. `rcausality(method = "markov")` runs the
Markov-property tests of `x`, `y` and the joint, the two directed CMI tests
and the contemporaneous test, and maps their outcomes through the decision
table of the paper (`$markov` shows the six p-values and the case).

```r
mk <- sim_markov01(600, scenario = "x_to_y", seed = 1)     # coupled binary chains
rm <- rcausality(mk, method = "markov", n_perm = 100, n_markov = 49, seed = 1)
rm$markov                                                  # Markov p-values of x, y, joint; the case
rm$classes
markov_test(mk$x, m = 1, n_perm = 49)$p_value              # the driver alone is Markov: large p
```

### Degrees of freedom / conditional variance (`"df"`)

Telcs, Kurbucz & Jakovac (2025): the *order* of a series is the number of
past lags that still reduce the conditional variance of its present
(`df_variance()`), and Table 2 of the paper reads the causal case from the
orders of `x`, of `y` and of the joint process (`$df`). The directed
"constrain the cause's past, watch the effect's spread" tests (`$pairs`)
decide the pairs the table leaves open.

```r
rd <- rcausality(d, method = "df", n_perm = 100, seed = 1)
rd$df                                                      # orders of x, y, joint; the Table-2 case
rd$classes
dv <- df_variance(d$y); dv                                 # how many past lags determine y[t]?
plot(dv)                                                   # the conditional sd against the window width
```

### Bayesian dimensional causality (`"dc"`)

Benko et al. (2024), for deterministic systems: the intrinsic dimensions of
the reconstructed manifolds of `X`, `Y`, their joint `J` and an independent
joint `Z` decide the case (`X -> Y` iff `D_X < D_Y = D_J`; a hidden common
cause iff `max(D_X, D_Y) < D_J < D_Z`; independence iff `D_J = D_Z`), with a
posterior probability for each of the five cases.

```r
lb <- sim_logistic_bdc(3000, "common_cause", seed = 1)     # a hidden third map drives both
rcausality(lb, method = "dc", m = 4, k_range = seq(12, 44, by = 8), seed = 1)$dc
dc <- dimensional_causality(lb$x, lb$y, m = 4, k_range = seq(12, 44, by = 4), seed = 1)
dc
plot(dc)                                                   # dimension-vs-k plateaus + posterior bars
```

**What to look for.** `plot(dc)` shows the four dimension estimates against
the neighbourhood size `k`; choose `k_range` where they form a plateau. For
the common-cause system `D_X = D_Y < D_J < D_Z` and the posterior of
`common_cause` is close to 1.

### CCDH: compressibility of rank-order patterns (`"ccdh"`)

Stippinger et al. (2023, CSF 176): the ordinal patterns of a series are
compressed by NSRPS; a driver is reconstructible from the driven series, so
the joint "driven + lagged driver" compresses as well as the driven series
alone (reconstruction index near 0), whereas with a hidden common cause both
joints are harder to compress than the single series.

```r
l3 <- sim_logistic3(1000, "driver", b = 0.4, c = 0.4, seed = 1)   # x0 -> x1, x0 -> x2
rcausality(l3[, c("x0", "x1")], method = "ccdh", m = 3, seed = 1)$classes   # x0 -> x1
cc <- ccdh(l3$x1, l3$x2, m = 3, seed = 1)                # x1, x2: their common driver x0 is hidden
cc                                                       # verdict: common_cause
plot(cc)                                                 # compression-rate curves of the four sequences
```

### df-entropy with a random forest (`"dfe"`)

Jakovac et al.: the series are binarised, conditional-entropy curves are
computed for growing own and cross histories, and a random forest trained on
coupled binary Markov chains of the five classes scores the pair (needs the
`ranger` package). The model is trained once per session (or once, with
`dfe_train()`, and reused).

```r
mdl <- dfe_train(n_per_class = 200, seed = 1, cores = 2)   # about 20 s; done once
mk2 <- sim_markov_pair(2000, "common_cause", seed = 1)     # hidden z drives x and y
dfe_classify(mk2$x, mk2$y, model = mdl)                    # class scores of the five cases
rcausality(d, method = "dfe", model = mdl)$classes         # the VAR pair of the quick start
```

### Information atoms: MAR-PID (`"pid"`)

Telcs (2026): the information that the cause's past and the effect's own
past carry about the effect's future is decomposed into redundant, unique and
synergistic atoms. The directed statistic is the "novel" information of the
cause's past (unique + synergistic), tested with surrogates; `$pid` lists the
atoms of every pair.

```r
rp <- rcausality(d, method = "pid", seed = 1)
rp$pid                                                     # the atoms per directed pair
marpid(target = d$y[-1], s1 = d$x[-nrow(d)], s2 = d$y[-nrow(d)], n_surr = 100, seed = 1)
```

### Baselines: transfer entropy and convergent cross mapping

`"te"` is Schreiber's transfer entropy on quantile-binned series with
circular-shift surrogates (and the effective TE, the surrogate mean
subtracted); `"ccm"` is Sugihara's convergent cross mapping, which declares
`x -> y` when the cross-map skill of `y xmap x` is significant *and*
converges with the library size.

```r
rcausality(d, method = "te", seed = 1)$classes
transfer_entropy(d$x, d$y, seed = 1)[c("te", "te_eff", "p_value")]
rcausality(lg, method = "ccm", m = 2, seed = 1)$ccm        # skill, convergence test per direction
ccm(lg$x, lg$y, E = 2, seed = 1)$rho                       # cross-map skill against the library size
```

### Recurrence plots and correlation dimensions (`"hirata"`, `"cd"`)

Two further criteria for deterministic systems: Hirata & Aihara (2010) --
the recurrences of the driven series are covered by those of the driver,
and excess simultaneous recurrences without a direction point to a common
cause; Krakovska (2019) -- the correlation dimensions of `X`, `Y`, the joint
and an independent joint obey the same relations as in the dimensional
causality method.

```r
rcausality(lb, method = "hirata", m = 3, seed = 1)$classes
rcausality(lb, method = "cd", m = 3, seed = 1)$classes
h <- hirata_rp(lb$x, lb$y, E = 3, seed = 1); h; plot(h)   # recurrence plots and coverage picture
krakovska_cd(lb$x, lb$y, E = 3, seed = 1)
```

### PCMCI through tigramite (`"pcmci"`, optional)

`pcmci()` wraps the authors' Python package tigramite (reached through
reticulate; nothing is ported). It needs `install.packages("reticulate")`
and a Python with `pip install tigramite`; `pcmci_available()` tells whether
they are found. The result is an ordinary `causres`.

```r
if (pcmci_available()) {
  pcmci(net, tau_max = 3)                     # PCMCI, ParCorr MCI tests, BH-adjusted
  pcmci(net, tau_max = 3, plus = TRUE)        # PCMCI+ (also contemporaneous links)
  rcausality(net, method = "pcmci")           # the same through the dispatcher
}
```

## Bayesian variants of every method

`bayes = TRUE` switches any method to its Bayesian variant; the priors, the
Monte Carlo and MCMC sizes and the decision rule are collected in
`bayes_control()`. Two families:

* **Regression family** (`granger`, conditional Granger): Bayesian linear
  regression with a spike-and-slab prior over the cause's lags (the default;
  posterior inclusion probability of every lag, Bayes factor of "at least
  one lag" against "none"), a normal-slab block prior, or the closed-form
  `g`, `hyper_g` and `bic` factors. The MCMC chains give credible intervals
  of the coefficients and of Geweke's measure, an effective sample size and
  the Gelman-Rubin `R-hat` across chains.
* **Surrogate family** (`cmi`, `cmi_lag`, `te`, `pid`, `df`, `ccdh`,
  `hirata`, `cd`, `ccm`): the observed statistic is compared with its
  surrogate sample through `surrogate_bf()` -- a kernel-density or normal
  null against an alternative shifted by an effect with a half-Cauchy
  (default), half-normal, normal or Cauchy prior -- which yields a Bayes
  factor and a posterior of the effect size. `markov` combines its six
  surrogate Bayes factors in a Bayesian decision table; `dc` and `dfe` take
  a prior over the five causal cases; `pcmci` uses a Savage-Dickey factor.

```r
rb <- rcausality(net, bayes = TRUE, seed = 1)               # spike-and-slab Granger regression
rb$bayes                                                    # BF10, posterior probability, effect interval, ESS, R-hat per link
plot(rb, "posterior")                                       # posterior probabilities and effect intervals
rcausality(net, bayes = TRUE, control = bayes_control(prior = "normal", scale = 0.5, n_iter = 4000, chains = 4), seed = 1)
rcausality(net, bayes = TRUE, prior = "hyper_g")            # closed-form factors: "bic", "g", "hyper_g"
bayes_granger(net$v2, net$v1, lag = 2, seed = 1)            # one pair in detail

rcausality(lg, method = "cmi", bayes = TRUE, max_lag = 3, seed = 1)$bayes   # Bayesian surrogate test
rcausality(lg, method = "te", bayes = TRUE, seed = 1,
           control = bayes_control(effect_prior = "half_normal", scale = 2, null_model = "normal"))
rcausality(lb, method = "dc", bayes = TRUE, m = 4,
           control = bayes_control(case_prior = c(1, 1, 1, 4, 4)))          # prior over the five cases
rcausality(net, bayes = TRUE, control = bayes_control(decision = "posterior"), prob_threshold = 0.9, seed = 1)
surrogate_bf(3.1, rnorm(200))                               # the generic Bayes factor of a statistic vs its surrogates
```

**What to look for.** In `$bayes`, `BF10` far above 1 with a posterior
probability near 1 is a link; `BF10` below 1 is evidence for its absence (a
p-value cannot say that). The `effect` column is the posterior mean and
credible interval of the effect (Geweke's measure for the regression family,
the shift in surrogate standard deviations for the surrogate family). For
the sampled variants `R-hat` close to 1 means that the chains agree; raise
`n_iter` or `chains` if it is not.

## Comparing and benchmarking methods

`compare_methods()` runs several methods on the same data and tabulates the
verdict of every pair under every method with a majority-vote consensus;
`benchmark_methods()` simulates data with a known truth and scores the
methods against it.

```r
cmp <- compare_methods(net, methods = c("granger", "cmi", "te", "dc"),
                       args = list(cmi = list(max_lag = 3), dc = list(m = 3)), seed = 1, cores = 2)
cmp
plot(cmp)                                     # tiles: verdict per pair and method, majority vote
bm <- benchmark_methods(c("granger", "te", "cmi"), simulator = sim_var, n = 400, reps = 10,
                        seed = 1, cores = 2)   # replications in parallel
bm
plot(bm)                                      # accuracy per method and scenario
```

**What to look for.** In the comparison, pairs on which the methods
disagree deserve a look at the individual results (`cmp$results$cmi`, say).
On the VAR network above `granger`, `cmi` and `te` agree on the three true
links, except that the pairwise F-test also reports the indirect
`v1 -> v3`; `dc` calls every pair independent, because it assumes
deterministic dynamics and a stochastic VAR is outside its domain -- a
reminder that each method answers its own question. In the benchmark,
`accuracy` is the fraction of replications in which the method recovered
the true class, per scenario; `$confusion` shows what the errors were.

## Causality-informed forecasting

The significant links define a restricted VAR: every series is regressed on
its own lags and on the lags of its detected causes. `predict()` measures the
gain over the own-lags-only model on a hold-out part (rolling one-step
forecasts) and then forecasts `h` steps ahead with intervals;
`forecast_causal()` does both steps in one call.

```r
fc <- predict(rcausality(net), h = 8)        # class "causfc": print / summary / plot
fc
plot(fc)                                     # last observations, hold-out fit and the forecast fan
forecast_causal(net, method = "granger", h = 8)
```

**What to look for.** The hold-out table compares the RMSE of the causal
model with that of the own-lags model for every series: a driven series
gains, an undriven one does not (the gain is about 0).

## Plots

```r
r <- rcausality(net)
plot(r, "all")                        # network, causal matrix, lag scan, joint embedding
plot(r, "network")                    # arrows = significant links, labelled with lag and p; node roles
plot(r, "matrix")                     # cause (rows) -> effect (columns), starred = significant
plot(r, "lags")                       # statistic of every directed pair against the candidate lag
plot(r, "embedding")                  # time-delay embedding of every series
plot_embedding(lg$x, m = 3)           # 3-D delay reconstruction of one series
plot_embedding(lg$x, lg$y)            # joint state (cause[t-1], effect[t], cause[t])
```

## Speed: several cores

Every function with repeated independent work -- permutations and surrogates,
series pairs, MCMC chains, bootstrap resamples, benchmark replications --
takes `cores`. `rcausality()` distributes the pairs when there are at least
as many pairs as cores, otherwise the permutations of every test, so a
two-series analysis on eight cores also uses all eight. Every unit of work
has its own random stream, so **the results are identical for any number of
cores** given the seed: parallelisation changes the running time, never a
p-value, a Bayes factor or a verdict. See `?Rcausality-parallel`.

```r
rcausality(net, method = "cmi", cores = 4, seed = 1)       # four forked workers (socket workers on Windows)
rcausality(net, method = "cmi", cores = "auto", seed = 1)  # all physical cores but one
options(Rcausality.cores = 4)                              # the default for every call of the session
cl <- parallel::makeCluster(4)                             # or a cluster you manage yourself
rcausality(net, method = "cmi", cores = cl, seed = 1); parallel::stopCluster(cl)
```

## Your own data

```r
my <- read.csv("series.csv")          # numeric columns = time series (drop date / text columns)
my <- prep_series(my, diff = 1)       # optional: difference, impute missing values
res <- rcausality(my)                 # Granger on all pairs; add method = ... for the others
```

Series must be numeric and of equal length; missing values are interpolated
when `impute = TRUE` (the default). Non-stationary series should be
differenced or detrended first (`prep_series()`), and the kNN methods
standardise the series internally.

## Reproductions of published results

Each subsection names the paper, the exact location of the reproduced result
in it, the code (in `inst/reproduce/`, also installed with the package under
`system.file("reproduce", package = "Rcausality")`) and what we obtained.
We only claim what was actually run; open items are listed explicitly.

### Zlatniczki et al. (2021) -- lag-shifted conditional mutual information

* **Reference:** Zlatniczki, A., Stippinger, M., Benko, Z., Somogyvari, Z. and
  Telcs, A. (2021). Relaxation of some confusions about confounders.
  *Entropy*, 23(11), 1450. <https://doi.org/10.3390/e23111450>
* **Where in the paper:** Relation 1 (Sect. 2.1, "The decomposable case");
  Sect. 2.2 (confounder case); Sects. 2.4-2.6 (strict reversed / strict
  triangle inequality, equality); Sect. 2.7 and Appendix A (convergence of the
  conditional mutual information in the history depth).
* **Note:** the paper is theoretical and has no numerical experiments, so
  these are *illustrations of its propositions* on simulated systems with known
  truth, not reproductions of published figures.
* **Code:** `inst/reproduce/zlatniczki2021.R`
  (`rcausality(method = "cmi_lag")`, `cmi_convergence()`, `lag_triangle_check()`).

| Simulated system (truth) | Paper statement | `cmi_lag` result |
|---|---|---|
| x drives y (lag 1) | `I_B = 0`, `I_A > 0`, `I^(k) = 0` => A -> B | `x_to_y`, k = 0 CMI = 0.008 (n.s.) |
| independent AR(1) series | all quantities zero | `independent` |
| shared innovation (hidden common driver) | `I^(k) > 0`, no lagged link => common driver | `common_cause`, k = 0 CMI = 0.42 (p = 0.01) |

Appendix A: for the common-driver system the k = 0 CMI stays clearly positive
as the history depth grows (0.42, 0.33, 0.27, 0.23 for d = 1..4), while for the
directly driven system it stays at about 0 (0.003, 0.004, 0.007, -0.005).
The geodesic checks return "direct driving" for tau_LY > tau_LX + tau_XY,
"common driver" for <, and "ambiguous (confounder / chain L->X->Y)" for =.

### Stippinger et al. (2023) -- Markov-chain approach

* **Reference:** Stippinger, M., Bencze, A., Zlatniczki, A., Somogyvari, Z.
  and Telcs, A. (2023). Causal discovery of stochastic dynamical systems: a
  Markov chain approach. *Mathematics*, 11(4), 852.
  <https://doi.org/10.3390/math11040852>
* **Where in the paper:** Sect. 6.1 "Simulated Data", the "turn-wheel"
  example, Figure 2; decision tables (10)-(12) in Sect. 4.
* **Code:** `inst/reproduce/stippinger2023.R`
  (`sim_turnwheel()`, `markov_test()`, `rcausality(method = "markov")`).
* **What is reproduced:** Alice (driver) is a Markov chain (large p), Bob and
  Barry are not Markov (small p) because Alice acts as hidden memory, the
  lagged joint (Alice_t, Bob_{t+1}) is Markov again, the joint (Bob, Barry) is
  not. For discrete (symbolic) states the test uses exact contingency tables
  conditional on the present state.
* **Result** (20 repetitions, n = 5000, 99 local permutations; the paper uses
  100 repetitions and n = 2000-50000):

  | series | median p | rejection rate at 5% | paper (Fig. 2) |
  |---|---|---|---|
  | Alice (driver) | 0.58 | 0.10 | Markov (median p about 0.4) |
  | Bob | 0.010 | 1.00 | not Markov (median p about 0.01) |
  | Barry | 0.010 | 1.00 | not Markov (median p about 0.01) |
  | joint (Alice_t, Bob_{t+1}) | 0.63 | 0.00 | Markov (median p about 0.2-0.3) |
  | joint (Bob_t, Barry_t) | 0.020 | 0.90 | not Markov (median p about 0.01-0.03) |

  The paper's medians are read off the printed values above the box plots of
  Figure 2.
* **Coupled logistic maps (Figure 3, no dynamical noise; `sim_logistic3()`,
  m = 3, 2 bins, n = 2000, 5 runs, 49 surrogates):** the pattern of the
  figure is reproduced -- driver Markov (median p about 0.4; paper 0.45),
  driven1 and driven2 not Markov (0.05-0.09; paper 0.01), the lagged joint
  (driver_t, driven1_{t+1}) Markov (0.4-0.5; paper 0.22), the joint
  (driven1, driven2) not Markov (0.02; paper 0.02). The paper's series carry
  observational noise of unstated size; with `obs_noise = 0.02` our test
  rejects the driver as well, so this sensitivity is an open item. With
  dynamical noise sigma = 0.02 (Figure 4) the driver and the driven1 results
  hold, driven2 and the lagged joint are borderline.
* **Implementation notes:** continuous states use Algorithm 1-2 of the paper
  (kNN conditions of size sqrt(n), local random projections, 2 bins,
  Yates-corrected modified log-likelihood test, local-permutation surrogates
  of the futures). The local p-values are aggregated by a clipped geometric
  mean instead of the paper's harmonic mean: the harmonic mean is dominated
  by the few perfectly dependent neighbourhoods that smooth deterministic
  dynamics produce and then rejects the Markov property of the driver
  (`aggregate = "harmonic"` restores the paper's choice).

### Telcs, Kurbucz and Jakovac (2025) -- degrees of freedom (conditional variance)

* **Reference:** Telcs, A., Kurbucz, M. T. and Jakovac, A. (2025). Unified
  causality analysis based on the degrees of freedom. *Physical Review E*,
  112(1), 014204. <https://doi.org/10.1103/PhysRevE.112.014204>
* **Where in the paper:** Definition 1 and Sect. 6.1 (determining the df by
  constraining the past within windows of width sigma_W), Table 2 (orders of
  x, y and the joint process decide the causal case), Sect. 6.4 with
  Eqs. (25)-(30) and Figures 1-2 (linear rotation system, z drives x and y),
  Sect. 7 and Figures 6-7 (chickens and eggs).
* **Code:** `inst/reproduce/telcs2025_df.R` (`sim_df_rotation()`,
  `df_variance()`, `rcausality(method = "df")`).
* **Implementation note:** constraints are imposed as in Eqs. (14)-(15) with
  every observed past configuration as a centre (so no data are wasted on a
  single target value); the `sigma(sigma_W, n_constr)` profile is the figure
  of the paper, and the decision uses adaptive (k-nearest-neighbour) windows
  with a local-permutation test of every added constraint. In
  `rcausality()` the Table-2 reading of the orders decides a pair whenever
  it is determinate (on the rotation system, n = 1500, `max_constr = 4`:
  orders 4 for x, 2 for z and 4 for the joint -> z -> x, as in the paper);
  the directed "constrain the cause's past, watch the effect's spread" test
  of Sect. 7 decides the undetermined pairs and is always reported. The
  directed test alone reads the deterministic rotation as bidirectional,
  because the driven series reconstructs its driver.
* **Result, deterministic set (n = 1500):** z subsystem: the N = 0 line is
  flat and all N >= 1 lines are proportional to the window width -> one
  constrained pair determines z (df = 1 pair, the paper's "df = 2") -- both
  the step tests and the profile reading give 1. x subsystem: the N = 0 and
  N = 1 lines saturate and the N >= 2 lines are proportional to the width ->
  profile reading 2 pairs (the paper's "df = 4, two pairs"); the step tests
  with adaptive windows stop at 1 pair, because at n = 1500 the second pair
  only lowers the conditional sd by about 1 % (the adaptive windows are
  wider than the fixed widths at which the lines separate).
* **Result, chickens and eggs** (`lmtest::ChickEgg`, 54 yearly values,
  differenced and standardised, lag 2 as in the paper): constraining the egg
  past reduces the conditional sd of the chicken series by 7.2 %, the
  reverse by 4.6 % -- the same ordering as the paper's Figs. 6-7
  (egg -> chicken), but with 53 observations neither reduction is
  significant (p = 0.29 and 0.30), so we do not claim the paper's
  conclusion is reproduced, only its direction.
* **Not yet reproduced:** the joint x-y constraint scheme of Fig. 4
  (three pairs, hidden common cause); the Table-2 reading in `rcausality()`
  implements the common-cause row but has not been checked on this system.

### Benko et al. (2024) -- Bayesian dimensional causality

* **Reference:** Benko, Z., Zlatniczki, A., Stippinger, M., Fabo, D.,
  Solyom, A., Eross, L., Telcs, A. and Somogyvari, Z. (2024). Bayesian
  inference of causal relations between dynamical systems. *Chaos, Solitons
  & Fractals*, 185, 115142. <https://doi.org/10.1016/j.chaos.2024.115142>
* **Where in the paper:** Sects. 2.1-2.3 (dimensions of X, Y, the joint J
  and the independent joint Z; Eqs. 1-4), Appendix A.2 (the nine events and
  the Bayesian integrals), Appendix B.1 (workflow), Appendix B.2 and Figure 1
  (three coupled logistic maps, Eq. B.1, r = 3.99, beta = 0.5, N = 10000,
  rank normalisation, M = 4, k = 12..44).
* **Code:** `inst/reproduce/benko2024_dc.R` (`sim_logistic_bdc()`,
  `dimensional_causality()`, `rcausality(method = "dc")`).
* **Result** (N = 10000, M = 4, k = 12..44 step 4; posterior of the true
  case): X -> Y 0.83, X <-> Y 0.75, hidden common cause 1.00, independent
  0.89 -- the true case gets the highest posterior in all four patterns, as
  in Fig. 1C (the paper shows bars, not numbers). The dimension-versus-k
  plateaus (Fig. 1B) are reproduced by `plot()`; e.g. for the common-cause
  case D_X = D_Y = 2.3 < D_J = 2.9 < D_Z = 3.6.
* **Implementation note:** written from the paper; the authors' AGPL-3
  reference implementation (`adam-zlatniczki/dimensional_causality`) was
  not used, so no code is ported. The likelihood integrals of Appendix A.2
  are evaluated on a grid (+/- 3 standard deviations, 25 points per free
  dimension); the effective sample size is n/(2k).

### Stippinger et al. (2023) -- CCDH, complexity-based discovery of a hidden common cause

* **Reference:** Stippinger, M., Varga, B., Benko, Z., Fabo, D., Eross, L.,
  Telcs, A. and Somogyvari, Z. (2023). CCDH: Complexity based causal
  discovery of hidden common cause in time series. *Chaos, Solitons &
  Fractals*, 176, 114054. <https://doi.org/10.1016/j.chaos.2023.114054>
* **Where in the paper:** Sect. 2 (NSRPS compression, rank-order coding),
  Sect. 3 (skewed joints, Eqs. 7-10, the decision table and Propositions
  1-3), Appendix 7.1 (worked NSRPS example, low-dimensional joint), Sect. 4
  and Figure 1 (coupled logistic maps, Eqs. 14-17).
* **Code:** `inst/reproduce/stippinger_ccdh.R` (`ccdh()`,
  `rcausality(method = "ccdh")`, `sim_logistic3()`).
* **Result:** the worked example of Appendix 7.1 is reproduced exactly
  (`1221211221212211` -> 11 symbols after the first step). On the logistic
  maps (n = 1000, m = 3, 25 systems per scenario) the three panels of Fig. 1
  come out as in the paper: (a) the joint of the driven series with the
  driver is as compressible as the driven series while the joint of the
  driver with the driven is not (reconstruction index 0.24 vs 0.86 in the
  median); (b) with a hidden common driver both joints are less
  compressible than the singles (indices 0.64-0.67); (c) with mutual driving
  the joints stay close to the singles (indices about 0.3; mean compression
  rates 0.54 vs 0.48, the paper's curves coincide). With the default
  decision rule the verdicts over the 25 systems were: (a) 24 x x -> y;
  (b) 19 x common cause, 4 x independent, 2 others; (c) 23 x bidirectional.
* **Implementation note:** the paper leaves the significance test open
  ("surrogate data and test of significance are necessary"); `ccdh()` adds
  circular-shift surrogates of the added series: `p_share` tests that the
  two series share information at all, and the reconstruction index
  `Delta / Delta_surrogate` (0 = the added series is a function of the
  embedding, 1 = as uninformative as a surrogate) is thresholded at
  `theta = 0.45`, calibrated on the logistic-map scenarios (independent
  series give indices about 1 and no significant sharing).

### Jakovac et al. -- degrees-of-freedom method with conditional entropies and a random forest

* **Reference:** Jakovac, A., Fabo, D., Eross, L., Stippinger, M., Benko, Z.,
  Somogyvari, Z. and Telcs, A. Causality beyond the Reichenbach principle:
  the degrees-of-freedom method (manuscript; DOI to be added on
  publication). Builds on Telcs, Kurbucz and Jakovac (2025), *Physical
  Review E*, 112(1), 014204.
* **Where in the paper:** Sect. 2.1 (df relations, Table 1), Sect. 2.2-2.3
  (conditional-entropy curves, curvature features, surrogate correction,
  random-forest classification on simulated binary Markov chains), Sect. 3
  and Figures 2-3 (confusion matrix at length 2000, accuracy versus
  length), Appendix F (pipeline: mean-threshold binarisation, Lmax = 6,
  Kmax = 5, surrogate shift 10 Lmax + 1, 1000 training series per class of
  lengths 200-5000, 100 trees).
* **Code:** `inst/reproduce/jakovac_dfe.R` (`sim_markov_pair()`,
  `dfe_train()`, `dfe_classify()`, `rcausality(method = "dfe")`).
* **Result:** with 600 training series per class (lengths 200-5000) and
  the default forest, fresh chains of length 2000 are classified with
  91 % accuracy and chains of length 5000 with 93 % (40 per class; paper:
  92.4 % at length 2000, "no substantial improvement at 5000"). The
  remaining confusions are the same as in the paper's Fig. 2: independent
  vs common cause, and a few bidirectional pairs read as unidirectional.
  The default model trained at first use (200 per class) reaches about 85 %.
* **Beyond the paper's training domain:** on linear VAR pairs
  (`sim_var()`, n = 2000) the Markov-trained model gets the directed,
  bidirectional and common-cause cases right (10/10, 10/10, 10/10, 9/10)
  but calls 6 of 10 independent pairs "common cause"; on coupled logistic
  maps it is unreliable without the tuning the paper applies to its Lorenz
  and logistic benchmarks (not implemented).
* **Implementation notes:** the feature set follows Appendix F.0.5 (two
  largest curvatures of every entropy curve, the same for the
  surrogate-corrected curves, asymmetry measures) and adds the
  surrogate-corrected entropy surfaces themselves, which raised the test
  accuracy from 84 % to 91 %. The exact asymmetry measures of the paper are
  not specified, so ours are differences of coupling gains and entropy rates
  between X and Y. Needs the `ranger` package (Suggests).

### Telcs (2026) -- finite-resolution MAR-PID

* **Reference:** Telcs, A. (2026). Empirical nonnegative finite-resolution
  MAR-PID for continuous variables. *Entropy*, 28(6), 641.
  <https://doi.org/10.3390/e28060641>; the decomposition itself is the
  Blackwell-order PID of Mages, Anastasiadi and Rohner (2024), *Entropy*
  26(5), 424.
* **Where in the paper:** Sect. 3 (recursive quantile binarisation), Sect. 4-5
  (bit-level MAR-PID, aggregation over target bits), Sect. 6.1 and Table 1
  (XOR under the doubling map), Sect. 6.2 and Table 2 (mixed
  redundancy-synergy example), Sect. 6.3 (four-valued target).
* **Code:** `inst/reproduce/telcs2026_marpid.R` (`marpid()`,
  `rcausality(method = "pid")`).
* **Result** (n = 4000): Table 2 is reproduced within sampling error (e.g.
  q = 0.5: I(A;V1) = 0.173, I(A;V2) = 0.173, I(A;V1,V2) = 0.651 vs 0.189,
  0.189, 0.656), with the atoms moving from pure redundancy (q = 0, R = 1)
  to pure synergy (q = 1, S = 1); the four-valued target gives one unique bit
  per source as in Sect. 6.3; the doubling-map XOR at resolution 3 gives a
  synergistic atom of 2.7 bits (paper: 3; the 64-state joint channel is
  estimated from 4000 samples, the value grows with n) and all other atoms
  below 0.01.
* **Implementation note:** two sources only (the general antichain lattice
  and Moebius inversion over Dedekind-many atoms are not implemented); the
  Blackwell join of two binary-input channels is computed as the upper
  concave envelope of the union of their zonogon boundaries. In
  `rcausality(method = "pid")` the sources are the cause's past (at the
  candidate lag) and the effect's own past, the target is the effect's
  future, and the directed statistic is the unique + synergistic ("novel")
  information of the cause's past with a surrogate test -- the unique atom
  alone is zero whenever the cause's channel is a garbling of the effect's
  own channel.

### Hirata and Aihara (2010), Krakovska (2019) -- literature extensions

* **References:** Hirata, Y. and Aihara, K. (2010). Identifying hidden
  common causes from bivariate time series by using the recurrence plot.
  *Physical Review E*, 81, 016203. <https://doi.org/10.1103/PhysRevE.81.016203>
  -- Krakovska, A. (2019). Correlation dimension detects causal links in
  coupled dynamical systems. *Entropy*, 21(9), 818.
  <https://doi.org/10.3390/e21090818>
* **What is implemented:** `hirata_rp()` / `method = "hirata"` -- the two
  criteria of the abstract of Hirata and Aihara: excessive simultaneous
  recurrences (joint recurrence rate over the product of the rates, surrogate
  tested) and the coverage of one series' recurrences by the other's
  (`P(R_X | R_Y)` close to 1 when X drives Y, non-coverage denies the
  direction); tight (rate/5) recurrences are conditioned on and loose (rate)
  ones cover, because with equal rates the two conditional probabilities
  coincide. `krakovska_cd()` / `method = "cd"` -- Grassberger-Procaccia
  correlation dimensions of X, Y, the joint and an independent joint, with
  the dimension relations of Krakovska (the same table as the dimensional
  causality method) tested by bootstrap z-tests plus a 10 % relative
  tolerance.
* **Code:** `inst/reproduce/hirata_krakovska.R`.
* **Result:** on the five logistic-map scenarios of Benko et al. (n = 1500,
  5 seeds each) both methods return the true case in 24 of 25 runs (Hirata:
  one common-cause run read as bidirectional; Krakovska: one bidirectional
  run undetermined). The papers' own data (wind measurements; Rossler and
  Lorenz systems) were not reproduced, and the full texts were not available
  to us, so the implementations follow the abstracts and the standard
  recurrence-plot / correlation-dimension definitions.

### Runge et al. (2019) -- PCMCI (optional wrapper, nothing ported)

* **Reference:** Runge, J., Nowack, P., Kretschmer, M., Flaxman, S. and
  Sejdinovic, D. (2019). Detecting and quantifying causal associations in
  large nonlinear time series datasets. *Science Advances*, 5(11), eaau4996.
  <https://doi.org/10.1126/sciadv.aau4996>; PCMCI+: Runge, J. (2020),
  *UAI 2020*, PMLR 124, 1388-1397.
* **What is implemented:** `pcmci()` / `method = "pcmci"` is a wrapper
  around the authors' Python package **tigramite** (reached through
  reticulate; both are optional, `pcmci_available()` tells whether they are
  installed). The MCI test statistics, raw and Benjamini-Hochberg adjusted
  p-values and (for PCMCI+) the contemporaneous links are mapped onto the
  `causres` object, so the result prints, plots and compares like the native
  methods.
* **Code:** `inst/reproduce/runge2019_pcmci.R`.
* **Result:** the four-variable linear toy model of the tigramite tutorial
  (the reference implementation accompanying the paper; n = 1000, ParCorr,
  `tau_max = 3`) is recovered exactly: the four lagged cross-links
  X1 -> X0 (lag 1, negative), X3 -> X1 (lag 1), X1 -> X2 (lag 2) and
  X3 -> X2 (lag 3) are the only significant ones (adjusted p < 1e-40), and
  PCMCI+ finds no contemporaneous link, as in the model.

### Causality-informed forecasting (`predict()`)

* **Idea:** the detected links of any method define a restricted VAR
  (Luetkepohl, 2005, *New Introduction to Multiple Time Series Analysis*,
  Springer): each series is regressed on its own `p` lags and on the `p`
  lags of its significant causes. `predict(<causres>, h)` compares this with
  the own-lags-only model by rolling one-step forecasts on the last 20 % of
  the sample and then produces recursive `h`-step forecasts with normal
  intervals; `forecast_causal()` does both steps in one call.
* **Result:** on `sim_var(600, n_series = 3, seed = 42)` the causal model
  lowers the hold-out RMSE by 10-12 % for the two driven series and leaves
  the undriven one unchanged (0 %), as it should.

### Bayesian variants of every method

* **References:** George, E. I. and McCulloch, R. E. (1993). Variable
  selection via Gibbs sampling. *JASA*, 88(423), 881-889.
  <https://doi.org/10.1080/01621459.1993.10476353> -- Liang, F., Paulo, R.,
  Molina, G., Clyde, M. A. and Berger, J. O. (2008). Mixtures of g priors
  for Bayesian variable selection. *JASA*, 103(481), 410-423.
  <https://doi.org/10.1198/016214507000001337> -- Rouder, J. N., Speckman,
  P. L., Sun, D., Morey, R. D. and Iverson, G. (2009). Bayesian t tests for
  accepting and rejecting the null hypothesis. *Psychonomic Bulletin &
  Review*, 16, 225-237. <https://doi.org/10.3758/PBR.16.2.225> -- Kass,
  R. E. and Raftery, A. E. (1995). Bayes factors. *JASA*, 90(430), 773-795.
* **What is implemented** (`bayes = TRUE`, settings in `bayes_control()`):
  - *Regression family* (`granger`, conditional Granger): Bayesian linear
    regression of the effect on its own past, the controls and the cause's
    past. `prior = "spike_slab"` (default) is stochastic-search variable
    selection over the cause's lags (point-mass spike, normal slab of
    scale `scale * sigma`, prior inclusion probability `inclusion_prob`);
    the 2^q lag subsets are enumerated exactly for q <= 12 lags (Gibbs
    sampling beyond), giving the posterior inclusion probability of every
    lag and the Bayes factor of "at least one lag" against "none".
    `prior = "normal"` compares the whole block with the slab against none;
    `"g"`, `"hyper_g"` (parameter `hyper_a`) and `"bic"` are the closed-form
    factors. Posterior draws (`chains` x `n_iter`, `burnin`, `thin`) give
    credible intervals (`credible`) of the coefficients and of Geweke's
    measure `log(1 + explained / sigma^2)`, with an effective sample size
    (`bayes_granger()`).
  - *Surrogate family* (`cmi`, `cmi_lag`, `te`, `pid`, `df`, `ccdh`,
    `hirata`, `cd`, `ccm`): a Bayesian surrogate test (`surrogate_bf()`):
    under H0 the observed statistic comes from the surrogate distribution
    (kernel density or normal fit, `null_model`), under H1 it is shifted by
    an effect with a half-Cauchy (default), half-normal, Cauchy or normal
    prior (`effect_prior`) of scale `scale` surrogate standard deviations;
    the Bayes factor is integrated by importance sampling (`n_draws`) and
    the posterior of the effect is summarised with a credible interval. The
    method-specific side conditions (CCM convergence, the coverage index
    of the recurrence criterion, the reconstruction index of CCDH) are
    kept.
  - *Markov-chain approach*: the six tests of a pair (x, y and the joint
    Markov properties, the two directed and the contemporaneous CMI tests)
    get surrogate Bayes factors; the 64 combinations of their outcomes,
    weighted by their posterior probabilities, are mapped through the
    decision table of Stippinger et al. (2023) to posterior probabilities
    of the five causal cases and of the two directions.
  - *Posterior-based methods* (`dc`, `dfe`): a prior over the five causal
    cases (`case_prior`) replaces the uniform one; `pcmci`: Savage-Dickey
    Bayes factor of the Fisher-transformed partial correlation.
  - Every variant fills `bf10`, `post_prob` (with `prior_prob`) and
    `evidence` (Jeffreys / Kass-Raftery labels), stores the per-link
    details in `$bayes` (and the fitted regressions in `$bayes_fits`),
    decides by `bf10 > bf_threshold` or, with `decision = "posterior"`, by
    `post_prob > prob_threshold`, and plots with `plot(r, "posterior")`.
* **Result:** on `sim_var(400, "x_to_y")` the spike-and-slab regression
  gives BF10 > 1e4 with inclusion probability 0.99 for lag 1 of x -> y and
  BF10 = 0.09 (inclusion 0.07) for y -> x; the three closed-form factors
  agree (BF10 > 1e14 and < 1). The Z-channel examples and the surrogate
  tests behave as expected on independent pairs (Bayes factors below 1)
  and on coupled ones (decisive evidence), see `tests/testthat/test-bayes.R`.

## Methods

| `method` | Reference |
|---|---|
| `"granger"` | Granger (1969), *Econometrica* 37(3), 424-438; Geweke (1982), *JASA* 77(378), 304-313 |
| `"bayes_granger"` | spike-and-slab Bayesian regression: George & McCulloch (1993), *JASA* 88(423), 881-889; closed-form BIC / g / hyper-g factors (`prior =`): Schwarz (1978), Liang et al. (2008); evidence labels: Kass & Raftery (1995), *JASA* 90(430), 773-795 |
| `"cmi"`, `"cmi_lag"` | Frenzel & Pompe (2007), *PRL* 99, 204101; Runge (2018), AISTATS PMLR 84; Zlatniczki et al. (2021), *Entropy* 23(11), 1450, doi:10.3390/e23111450 |
| `"markov"` | Stippinger et al. (2023), *Mathematics* 11(4), 852, doi:10.3390/math11040852 |
| `"df"` | Telcs, Kurbucz & Jakovac (2025), *Physical Review E* 112(1), 014204, doi:10.1103/PhysRevE.112.014204 |
| `"dc"` | Benko et al. (2024), *Chaos, Solitons & Fractals* 185, 115142, doi:10.1016/j.chaos.2024.115142 |
| `"ccdh"` | Stippinger et al. (2023), *Chaos, Solitons & Fractals* 176, 114054, doi:10.1016/j.chaos.2023.114054 |
| `"dfe"` | Jakovac et al., Causality beyond the Reichenbach principle: the degrees-of-freedom method (manuscript); needs `ranger` |
| `"te"` | Schreiber (2000), *PRL* 85, 461; effective TE: Marschinski & Kantz (2002), *Eur. Phys. J. B* 30, 275 |
| `"ccm"` | Sugihara et al. (2012), *Science* 338, 496, doi:10.1126/science.1227079 |
| `"pid"` | Telcs (2026), *Entropy* 28(6), 641, doi:10.3390/e28060641; Mages, Anastasiadi & Rohner (2024), *Entropy* 26(5), 424 |
| `"hirata"` | Hirata & Aihara (2010), *Phys. Rev. E* 81, 016203, doi:10.1103/PhysRevE.81.016203 |
| `"cd"` | Krakovska (2019), *Entropy* 21(9), 818, doi:10.3390/e21090818 |
| `"pcmci"` (optional) | Runge et al. (2019), *Science Advances* 5(11), eaau4996, doi:10.1126/sciadv.aau4996; Runge (2020), UAI, PMLR 124; wrapper around the Python package tigramite |
| `bayes = TRUE`, `bayes_control()` | George & McCulloch (1993), *JASA* 88(423), 881-889; Liang et al. (2008), *JASA* 103(481), 410-423; Rouder et al. (2009), *Psychon. Bull. Rev.* 16, 225-237; Kass & Raftery (1995), *JASA* 90(430), 773-795 |
| `predict()` / `forecast_causal()` | restricted VAR forecasts from the detected links; Luetkepohl (2005), Springer, doi:10.1007/978-3-540-27752-1 |

## Open items
The Markov test for near-deterministic series with observational noise
(Stippinger et al. 2023, Figs. 3-5); a kd-tree for series longer than about
10^4 points; the joint x-y constraint scheme of the df method (Telcs et al.
2025, Fig. 4). Further papers can be added on request; the `causres`
object, `compare_methods()` and `benchmark_methods()` take any new method.

## Development notes (checking the package inside a cloud-synced folder)
`devtools::check()` writes `Rcausality.Rcheck/` and lock directories next to
the package; a syncing client (Dropbox, OneDrive, ...) can re-create files
while R deletes them, which shows up as "failed to lock directory",
"conflicted copy" files or "there is no package called 'Rcausality'". These
are not package errors. Run the check in a temporary directory and reinstall
afterwards:

```r
unlink(file.path(getwd(), "Rcausality.Rcheck"), recursive = TRUE)   # stale check output, if any
devtools::check(check_dir = tempdir())
devtools::install()        # the examples and scripts use the installed package
```

Also delete any `00LOCK-*` directory in the library path (`.libPaths()`) if
an installation was interrupted. The optional PCMCI wrapper needs
`install.packages("reticulate")` and a Python with `pip install tigramite`;
point reticulate at it with `Sys.setenv(RETICULATE_PYTHON = "<path to python>")`
or `pcmci_available(python = "<path>")`.

## Ported code
Any code ported from the authors' Python reference implementations
(`dimensional_causality`, `df_causality`) will be explicitly marked as
**ported by the Rcausality authors**, in the source header, in `?` help pages
and in `inst/COPYRIGHTS`.
