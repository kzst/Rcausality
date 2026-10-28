# Internal helpers ---------------------------------------------------------

#' @importFrom stats pf lm.fit setNames
#' @importFrom Rcpp sourceCpp
#' @useDynLib Rcausality, .registration = TRUE
NULL

# ---- parallel execution -------------------------------------------------------
# Every unit of work (a permutation, a surrogate, a series pair, an MCMC
# chain, a benchmark run) gets its own L'Ecuyer-CMRG random stream, derived
# from the caller's RNG state (one draw).  The result therefore does not
# depend on how the units are distributed over the workers: cores = 1 and
# cores = 8 give identical numbers for the same seed, and the caller's RNG is
# left exactly where a serial run would leave it.
#
# cores = 1 -> plain lapply; an integer > 1 -> forked workers on Unix-alikes
# (parallel::mclapply) or a cached PSOCK cluster on Windows (the package must
# be installed there); a cluster object from parallel::makeCluster() is used
# as it is (and not stopped).

.par_env <- new.env(parent = emptyenv())

.plapply <- function(X, FUN, cores = 1L, ...) {
  n <- length(X)
  if (n == 0L) return(list())
  streams <- .rng_streams(n)
  f <- function(i) .with_stream(streams[[i]], function() FUN(X[[i]], ...))
  nc <- .n_cores(cores)
  if (nc == 1L || n < 2L) return(lapply(seq_len(n), f))
  nw <- min(nc, n)
  out <- if (inherits(cores, "cluster")) {
    .cluster_run(cores, nw, seq_len(n), f)
  } else if (.Platform$OS.type == "windows") {
    .cluster_run(.psock_cluster(nw), nw, seq_len(n), f)
  } else {
    parallel::mclapply(seq_len(n), f, mc.cores = nw, mc.set.seed = FALSE, mc.preschedule = TRUE)
  }
  bad <- vapply(out, inherits, logical(1), "try-error")
  if (any(bad)) stop("Parallel worker failed: ", conditionMessage(attr(out[[which(bad)[1]]], "condition")),
                     call. = FALSE)
  out
}

# the number of workers implied by `cores`: an integer, "auto" (all physical
# cores but one; at most 2 under R CMD check) or a cluster object
.n_cores <- function(cores) {
  if (inherits(cores, "cluster")) return(length(cores))
  if (is.null(cores) || identical(cores, "auto")) {
    cores <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
    chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    if (nzchar(chk) && !identical(chk, "false")) cores <- min(cores, 2L)
  }
  cores <- suppressWarnings(as.integer(cores)[1])
  if (is.na(cores)) stop("`cores` must be a positive integer, \"auto\" or a cluster from parallel::makeCluster().", call. = FALSE)
  max(1L, cores)
}

# snapshot / restore of the session's generator: the kind (as RNGkind() reports
# it) and the state vector.  The kind is restored explicitly: assigning
# .Random.seed alone leaves R's internal kind as it was until the next draw,
# and if the vector is then removed the session would silently continue with
# the kind last set (L'Ecuyer-CMRG after a stream), which is what must not
# happen to a caller that uses the default Mersenne-Twister
.rng_get <- function() {
  g <- globalenv()
  list(kind = RNGkind(), seed = if (exists(".Random.seed", envir = g, inherits = FALSE)) get(".Random.seed", envir = g, inherits = FALSE) else NULL)
}
.rng_put <- function(st) {
  g <- globalenv()
  suppressWarnings(try(do.call(RNGkind, as.list(st$kind)), silent = TRUE))
  if (!is.null(st$seed)) assign(".Random.seed", st$seed, envir = g)
  else if (exists(".Random.seed", envir = g, inherits = FALSE)) rm(".Random.seed", envir = g)
  invisible(NULL)
}

# n independent L'Ecuyer-CMRG streams seeded from one draw of the caller's RNG;
# the caller's generator (kind and state) is put back afterwards
.rng_streams <- function(n) {
  g <- globalenv()
  if (!exists(".Random.seed", envir = g, inherits = FALSE)) stats::runif(1L)
  base <- sample.int(.Machine$integer.max, 1L)
  saved <- .rng_get()
  on.exit(.rng_put(saved), add = TRUE)
  suppressWarnings(RNGkind("L'Ecuyer-CMRG"))
  set.seed(base)
  s <- get(".Random.seed", envir = g, inherits = FALSE)
  out <- vector("list", n)
  for (i in seq_len(n)) { out[[i]] <- s; s <- parallel::nextRNGStream(s) }
  out
}

# evaluate fun() with the RNG set to the given stream; the previous generator
# of the process (kind and state) is restored afterwards, so that in a serial
# run the units of work leave no trace in the caller's generator
.with_stream <- function(stream, fun) {
  saved <- .rng_get()
  on.exit(.rng_put(saved), add = TRUE)
  assign(".Random.seed", stream, envir = globalenv())
  fun()
}

# run f over idx on (at most nw nodes of) a cluster; the package namespace is
# loaded on the nodes (a no-op after the first call).  The closures are
# serialised to the nodes together with their enclosing frames, so every
# lazy argument in those frames is forced first: an unevaluated promise whose
# expression refers to the caller's workspace cannot be evaluated on a node
.cluster_run <- function(cl, nw, idx, f) {
  if (length(cl) > nw) cl <- structure(cl[seq_len(nw)], class = class(cl))
  .force_frames(environment(f))
  parallel::clusterEvalQ(cl, suppressWarnings(requireNamespace("Rcausality", quietly = TRUE)))
  parallel::parLapply(cl, idx, f)
}

.force_frames <- function(e, seen = list()) {
  while (!identical(e, globalenv()) && !identical(e, emptyenv()) && !isNamespace(e) &&
         !any(vapply(seen, identical, logical(1), e))) {
    seen[[length(seen) + 1L]] <- e
    for (nm in ls(e, all.names = TRUE)) {
      v <- try(get(nm, envir = e, inherits = FALSE), silent = TRUE)      # get() forces a promise
      if (is.function(v) && !is.primitive(v) && !is.null(environment(v))) seen <- .force_frames(environment(v), seen)
    }
    e <- parent.env(e)
  }
  seen
}

# PSOCK cluster kept for the session (starting one costs seconds on Windows);
# re-created when more workers are asked for or when it has died
.psock_cluster <- function(nw) {
  cl <- .par_env$cluster
  alive <- !is.null(cl) && length(cl) >= nw &&
    isTRUE(tryCatch({ parallel::clusterEvalQ(cl, TRUE); TRUE }, error = function(e) FALSE))
  if (alive) return(cl)
  .stop_cluster()
  cl <- parallel::makePSOCKcluster(nw)
  .par_env$cluster <- cl
  cl
}

.stop_cluster <- function() {
  cl <- .par_env$cluster
  if (!is.null(cl)) { try(parallel::stopCluster(cl), silent = TRUE); .par_env$cluster <- NULL }
  invisible(NULL)
}

# where to parallelise: over the n_units outer units (pairs) when there are
# at least as many of them as workers, otherwise inside each unit (over the
# permutations / surrogates / chains); never both at once
.split_cores <- function(n_units, cores) {
  nc <- .n_cores(cores)
  if (nc == 1L) return(list(outer = 1L, inner = 1L))
  if (n_units >= nc) list(outer = cores, inner = 1L) else list(outer = 1L, inner = cores)
}

# seed handling of the top-level functions: the seed (when given) applies to
# the call only; the caller's generator is restored by .rng_restore()
.rng_setup <- function(seed) {
  if (is.null(seed)) return(NULL)
  old <- .rng_get()
  set.seed(seed)
  old
}
.rng_restore <- function(old) {
  if (!is.null(old)) .rng_put(old)
  invisible(NULL)
}

.lagmat <- function(x, p) {
  n <- length(x)
  out <- matrix(NA_real_, n - p, p)
  for (l in seq_len(p)) out[, l] <- x[(p - l + 1):(n - l)]
  out
}

.rss <- function(y, X) {
  X <- cbind(1, X)
  fit <- lm.fit(X, y)
  sum(fit$residuals^2)
}

.as_series_frame <- function(x, y = NULL, xname = "x", yname = "y") {
  to_df <- function(v, nm) {
    if (is.data.frame(v)) return(v)
    if (is.matrix(v)) {
      v <- as.data.frame(v)
      if (is.null(colnames(v)) || any(colnames(v) == "")) {
        colnames(v) <- paste0(nm, seq_len(ncol(v)))
      }
      return(v)
    }
    if (inherits(v, "ts") || is.numeric(v) || is.logical(v)) {
      d <- data.frame(as.numeric(v)); names(d) <- nm; return(d)
    }
    if (is.function(v))
      stop("Input is a function, not data. Did you pass an undefined object called `df`/`data`/`t`? ",
           "(these names are base R functions) -- use e.g. `d <- data.frame(...)`.", call. = FALSE)
    stop("Input must be a numeric vector, matrix or data.frame.", call. = FALSE)
  }
  check_cols <- function(d) {
    if (is.null(d)) return(d)
    if (inherits(d, c("mts", "ts"))) d <- as.data.frame(unclass(d))
    bad <- names(d)[!vapply(d, function(v) is.numeric(v) || is.logical(v), TRUE)]
    if (length(bad))
      stop("Non-numeric column(s): ", paste(bad, collapse = ", "),
           ". Only numeric time series can be analysed; drop date, factor and character columns first ",
           "(e.g. `d[, sapply(d, is.numeric)]`).", call. = FALSE)
    d[] <- lapply(d, as.numeric)
    d
  }
  X <- check_cols(to_df(x, xname))
  Y <- if (is.null(y)) NULL else check_cols(to_df(y, yname))
  list(X = X, Y = Y)
}

.prep_na <- function(v) {
  v <- as.numeric(v)
  bad <- !is.finite(v)
  if (!any(bad)) return(v)
  if (all(bad)) stop("A series contains only missing values.", call. = FALSE)
  idx <- seq_along(v)
  v[bad] <- stats::approx(idx[!bad], v[!bad], xout = idx[bad], rule = 2)$y
  v
}


# unordered pairs of a grid of directed tests, (v1, v2) in the column order of
# the data (not alphabetical), so that a lag applied to "the second series"
# does not depend on the names
.unordered_pairs <- function(grid, all) {
  i1 <- match(grid$cause, names(all)); i2 <- match(grid$effect, names(all))
  v1 <- ifelse(i1 < i2, grid$cause, grid$effect); v2 <- ifelse(i1 < i2, grid$effect, grid$cause)
  up <- unique(data.frame(v1 = v1, v2 = v2, stringsAsFactors = FALSE))
  rownames(up) <- NULL; up
}


# seed the RNG for the calling function only: the caller's generator (kind
# and state) is restored when that function exits (registered as an on.exit
# of the caller)
.local_seed <- function(seed) {
  if (is.null(seed)) return(invisible(NULL))
  old <- .rng_get()
  set.seed(seed)
  do.call(on.exit, list(bquote(.rng_put(.(old))), add = TRUE), envir = parent.frame())
  invisible(NULL)
}
