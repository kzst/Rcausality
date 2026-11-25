## Rcausality 0.0.13 -- an update with two bug fixes and an option

We are sorry for an update so soon after 0.0.12 (submitted on 2026-11-21).
While checking the simulations of the paper that accompanies the package we
found two errors in documented options, and we would rather not leave them
to the users:

* `rcausality(..., bayes = TRUE, prior = "minnesota", lag = "first")`
  stopped with an error (`invalid 'nrow' value`).
* The Bayesian decision of the method `"df"` did not apply a condition of
  the frequentist test (a statistic that is not positive is no evidence of
  a link), so that it could declare links that the test rates with p = 1.

Both are corrected, with tests (`tests/testthat/test-fixes-0.0.13.R`). The
same check showed that the method `"df"` rarely finds the hidden common
driver of the example system of its own source paper with its estimator of
the conditional variance (windows of nearest neighbours, which cannot
resolve small conditional variances). The update adds the estimator with
which that paper computes the example (the Gaussian conditional variance)
and a search for the smallest number of conditions, as options; the
default results do not change (tests in
`tests/testthat/test-df-estimator.R`). Apart from these, only the
documentation of `lag` (the largest lag order of the regression methods),
the wording of a warning and the README change; see NEWS.md.

Local check (Ubuntu 24.04, R 4.3.3, two cores), `R CMD check --as-cran`:
0 errors. As for 0.0.12, the two warnings and four of the six notes come
from the check machine (locale `en_US.UTF-8` not available, missing
`inconsolata.sty` for the PDF manual and the `.tex` file it leaves behind,
no `tidy`, system clock not verifiable, a compiler flag); the fifth is the
installed size (5.1 Mb, of which the compiled code 2.0 Mb and the four
vignettes 1.5 Mb), the sixth the "New submission" note of the incoming
feasibility check. Run times: examples 14 s (with `--run-donttest` 44 s),
tests 44 s, re-building of the vignettes 44 s.

## Rcausality 0.0.12 -- resubmission after the review of 2026-11-14

Thank you for the review. The three points:

* *Single quotes in Title and Description.* Only the names of software are
  quoted now ('Python', 'tigramite'); the quotes around the class name and
  the method names are removed.
* *Acronyms in the Description.* All are written out or explained where
  they first occur (the complexity-based causal discovery of hidden common
  causes (CCDH); the combination of the Peter-Clark algorithm with momentary
  conditional independence tests (PCMCI); Bayesian information criterion,
  processor cores, partial information decomposition, degrees of freedom).
* *No more than 2 cores.* Examples, vignettes and tests ask for at most two
  workers. The review made us look for other sources: the forests of one
  method were predicted without an explicit number of threads ('ranger'
  then uses all cores), which is corrected, and a test asserts it. A call
  made inside a worker of the package now always runs serially, so that the
  number of processes cannot exceed `cores`, and `cores = "auto"` respects
  `_R_CHECK_LIMIT_CORES_` and `OMP_THREAD_LIMIT`.

Run times (local: Ubuntu 24.04, R 4.3.3, two cores; `R CMD check --as-cran`
about 3.5 minutes in total): examples 11 s, with `--run-donttest` 38 s, the
slowest example 2.3 s; tests 38 s (the long-running tests are skipped on
CRAN); re-building of the four vignettes 44 s (their long-running results
are precomputed in `vignettes/precomputed/`, the code shown is the code
that produced them).

`R CMD check --as-cran`: 0 errors. One warning and two notes come from the
check machine (locale, system clock not verifiable, a compiler flag). One
note on the installed size (5.1 Mb on this machine, of which the compiled
code 2.0 Mb and the four vignettes 1.5 Mb).

The other changes of this version are described in NEWS.md.

## Rcausality 0.0.3

A unified interface (the method, its Bayesian variant and the number of cores
as options of every call), every method for panel data (unit, group and
pooled level, with Bayesian variants), the analysis of the causal drivers of
one explained variable, panel-aware placebo / robustness / comparison /
forecasting tools, and a fourth vignette (see NEWS.md).

* local: Ubuntu 24.04, R 4.3.3, `R CMD check --as-cran`: 0 errors; the only
  warnings / notes come from the check machine (locale, missing
  `inconsolata.sty`, compiler flags, no `tidy`).
* The long-running tests are skipped on CRAN; the simulation table and the
  model-free panel examples of the new vignette are precomputed
  (`vignettes/precomputed/`, the code shown is the code that produced them).
* `panel_granger()` (0.0.2) is kept as a wrapper of the new
  `panel_causality()`; its result class changed (`causpanel`).

## Rcausality 0.0.2

New methods and tools (see NEWS.md): Toda-Yamamoto, Diks-Panchenko and
random-forest Granger tests, frequency-domain decomposition, the
hierarchical Minnesota prior for short series, panel Granger tests,
placebo / time-reversal / robustness analyses, forecast models built from
the links of any method, and a third vignette.

* local: Ubuntu 24.04, R 4.3.3, `R CMD check --as-cran`: 0 errors; the only
  warnings / notes come from the check machine (locale, missing
  `inconsolata.sty`, compiler flags, no `tidy`).
* Tests about 45 s, examples (with --run-donttest) about 40 s, the three
  vignettes about 40 s together (long-running results precomputed in
  `vignettes/precomputed/`).
* `ranger` stays in Suggests: the random-forest functions check for it and
  the tests and vignette chunks skip without it.

## Rcausality 0.0.1 -- first CRAN submission (resubmission)

### Changes after the incoming pre-tests of 2026-10-29

* Windows ERROR "connections left open" in the examples of
  `?Rcausality-parallel`: the socket cluster used for `cores > 1` on Windows
  was kept for the session; it is now started for the call and stopped when
  the call returns, so no connection outlives a call.
* Windows NOTE "Overall checktime 13 min > 10 min": the long-running tests
  are skipped on CRAN (`skip_on_cran()`; they run with `devtools::test()`),
  and the benchmark / method-comparison results of the vignettes are
  precomputed (shipped in `vignettes/precomputed/`, the code shown is the
  code that produced them). Locally the tests now take about 45 s and the
  two vignettes about 30 s.
* Debian NOTE on the example of `pcmci_available()` (CPU time > 2.5 x
  elapsed, from starting Python): the examples that start Python are now in
  `\dontrun{}`, and the vignette does not start Python during `R CMD check`.
* The words flagged as possibly misspelled in DESCRIPTION are author names
  (Aihara, Benko, Granger, Hirata, Jakovac, Krakovska, Kurbucz, Runge,
  Stippinger, Telcs, Zlatniczki), method acronyms (CCDH, PCMCI, PID, df) and
  "et al."; they are spelled correctly.

### Test environments

* local: Ubuntu 24.04, R 4.3.3 (R CMD check --as-cran)
* macOS (maintainer's machine), R 4.5 (devtools::check())

### R CMD check results

0 errors | 0 warnings | 0 notes (apart from the "new submission" NOTE).

### Comments

* The package compiles a small Rcpp back-end (max-norm nearest-neighbour
  searches); no system requirements beyond a C++ compiler.
* The optional PCMCI wrapper (`pcmci()`) uses the Python package tigramite
  through reticulate (Suggests). All tests and vignette chunks that need it
  are skipped when it is not available, so no Python is needed for the
  checks.
* `ranger` (Suggests) is needed only by the df-entropy method; the tests
  and the vignette chunk skip without it.
* Examples run in well under 5 s each (the longer illustrations are in
  \donttest); the test suite takes about 45 s on CRAN, the two vignettes
  about 30 s together.
* Parallel execution (`cores`) is off by default; the examples and tests
  that exercise it use at most two workers (forked on Unix-alikes; the
  socket-cluster tests are skipped on CRAN).
* The DOIs in the Description refer to the published methods implemented
  by the package.
* No code is ported from the authors' reference implementations; the
  methods are implemented from the papers (see README, "Ported code").
