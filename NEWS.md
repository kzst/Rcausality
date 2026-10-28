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
  1L)`. On Windows the socket cluster is started once and kept for the
  session.
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
