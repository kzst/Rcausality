#' Parallel computation in Rcausality
#'
#' Every function that repeats independent units of work -- the permutations
#' and surrogates of the model-free tests, the series pairs of a multivariate
#' analysis, the MCMC chains of [bayes_granger()], the bootstrap resamples of
#' [corr_dim()], the replications of [benchmark_methods()] and the training
#' systems of [dfe_train()] -- takes a `cores` argument. Its default is
#' `getOption("Rcausality.cores", 1L)`, so that
#' `options(Rcausality.cores = 4)` switches the whole package to four
#' workers for the session.
#'
#' @section What `cores` accepts:
#' * `1` (the default): serial execution.
#' * an integer `> 1`: that many workers. On Unix-alikes (Linux, macOS)
#'   these are forked copies of the session (`parallel::mclapply()`), which
#'   need no set-up. On Windows a socket cluster
#'   (`parallel::makePSOCKcluster()`) is started on first use, kept for the
#'   session and closed when R exits; the package must be *installed* on the
#'   machine (the workers load it by name).
#' * `"auto"`: all physical cores but one (at most 2 under `R CMD check`).
#' * a cluster object from `parallel::makeCluster()`: used as it is and never
#'   stopped by the package. This is the way to parallelise over several
#'   machines, or to use socket workers on a platform where forking is
#'   discouraged (e.g. inside some IDEs).
#'
#' @section Where the work is split:
#' [rcausality()] distributes the *series pairs* over the workers when there
#' are at least as many pairs as workers, and otherwise the *permutations or
#' surrogates* of every test (so that a two-series analysis on eight cores
#' still uses all eight). The two levels are never nested. The standalone
#' tests (e.g. [local_permutation_test()], [transfer_entropy()], [ccm()],
#' [hirata_rp()], [ccdh()], [marpid()]) parallelise their surrogate loop;
#' [bayes_granger()] runs its chains in parallel; [benchmark_methods()]
#' distributes the simulated data sets. Parallel execution pays off when a
#' single unit of work takes more than a few milliseconds (kNN estimates,
#' recurrence counts, cross-map predictions, MCMC chains); for the cheapest
#' statistics the overhead of dispatching the work can exceed the gain.
#'
#' @section Reproducibility -- parallel results are *identical* to serial ones:
#' Every unit of work draws its random numbers from its own L'Ecuyer-CMRG
#' stream (`parallel::nextRNGStream()`), and the streams are derived from the
#' state of the session's generator at the time of the call. Consequently
#' * the result of a call does not depend on `cores`: `cores = 1`,
#'   `cores = 2` and `cores = 8` give the same numbers (to the last digit)
#'   for the same `seed`, on every platform, and
#' * the generator of the session is left exactly where a serial call would
#'   have left it, so a script that mixes Rcausality calls with other random
#'   draws stays reproducible when `cores` changes.
#' Parallelisation therefore never changes a p-value, a Bayes factor or a
#' verdict; it only changes the wall-clock time. Use `seed` (or `set.seed()`
#' before the call) for reproducible runs as usual.
#'
#' @section CRAN and continuous integration:
#' Package checks limit the number of processes to two; `cores = "auto"`
#' respects that limit automatically, and the examples and tests of the
#' package use at most two workers.
#'
#' @name Rcausality-parallel
#' @aliases parallel-Rcausality Rcausality.cores
#' @examples
#' d <- sim_logistic(300, "x_to_y", seed = 1)
#' # the same analysis serially and with two workers: identical results
#' r1 <- rcausality(d, method = "te", seed = 1, cores = 1)
#' r2 <- rcausality(d, method = "te", seed = 1, cores = 2)
#' all.equal(r1$pairs, r2$pairs)
#' \donttest{
#' # a kNN-based test, where the workers pay off
#' system.time(r1 <- rcausality(d, method = "cmi", max_lag = 2, n_perm = 60, seed = 1, cores = 1))
#' system.time(r2 <- rcausality(d, method = "cmi", max_lag = 2, n_perm = 60, seed = 1, cores = 2))
#' identical(r1$pairs, r2$pairs)
#'
#' # switch the whole session to two workers
#' options(Rcausality.cores = 2)
#' rcausality(d, method = "ccm", seed = 1)$classes
#' options(Rcausality.cores = 1)
#'
#' # a user-managed socket cluster (works on every platform; not stopped by the package)
#' cl <- parallel::makeCluster(2)
#' rcausality(d, method = "te", seed = 1, cores = cl)$classes
#' parallel::stopCluster(cl)
#' }
NULL
