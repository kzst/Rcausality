# Method comparison on one data set and simulation benchmarks ------------------

.class_levels <- c("x_to_y", "y_to_x", "bidirectional", "common_cause", "independent")
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

#' Run several causality methods on the same data and compare their verdicts
#'
#' Calls [rcausality()] once per method and tabulates the class of every
#' series pair under every method, with a majority-vote consensus.
#'
#' @param x,y data as for [rcausality()].
#' @param methods character vector of methods (see [rcausality()]).
#' @param args optional named list of per-method argument lists, e.g.
#'   `list(dc = list(m = 4), ccm = list(m = 3))`.
#' @param ... arguments passed to every method (e.g. `seed`, `cores`,
#'   `n_perm`, `alpha`).
#' @return An object of class `causcomp`: `table` (one row per pair, one
#'   column per method, plus `consensus` and `agreement`), `results` (the
#'   `causres` objects), `methods`. Methods: `print`, `plot` (tile plot of
#'   the verdicts).
#' @export
#' @examples
#' # The same pair under three methods.  A linear VAR coupling is found by the
#' # F-test, by the Bayes factor and by the transfer entropy; 'consensus' is
#' # the majority vote and 'agreement' the share of methods behind it.
#' d <- sim_var(300, "x_to_y", seed = 1)
#' cmp <- compare_methods(d, methods = c("granger", "bayes_granger", "te"), seed = 1)
#' cmp
#' \donttest{
#' plot(cmp)            # tile plot of the verdicts per pair and method
#' cmp$results$te       # the complete causres object of one of the methods
#' # method-specific arguments go in `args`; the dots reach every method
#' compare_methods(d, methods = c("granger", "cmi"), args = list(cmi = list(max_lag = 2)),
#'                 n_perm = 100, seed = 1)
#' }
compare_methods <- function(x, y = NULL, methods = c("granger", "bayes_granger", "cmi", "te"),
                            args = list(), ...) {
  dots <- list(...)
  results <- lapply(methods, function(m) {
    a <- c(list(x = x, y = y, method = m), dots, if (!is.null(args[[m]])) args[[m]])
    a <- a[!duplicated(names(a), fromLast = TRUE)]
    do.call(rcausality, a)
  })
  names(results) <- methods
  cls <- results[[1]]$classes[, c("var1", "var2")]
  for (m in methods) {
    C <- results[[m]]$classes
    cls[[m]] <- C$class[match(paste(cls$var1, cls$var2), paste(C$var1, C$var2))]
  }
  votes <- as.matrix(cls[, methods, drop = FALSE]); storage.mode(votes) <- "character"
  cls$consensus <- vapply(seq_len(nrow(votes)), function(i) {
    v <- votes[i, ]; v <- v[!is.na(v)]
    if (!length(v)) return(NA_character_)
    tb <- sort(table(v), decreasing = TRUE)
    if (length(tb) > 1 && tb[1] == tb[2]) "undetermined" else names(tb)[1]
  }, "")
  cls$agreement <- vapply(seq_len(nrow(votes)), function(i) {
    v <- votes[i, ]; v <- v[!is.na(v)]; if (!length(v)) NA_real_ else max(table(v)) / length(v)
  }, 0)
  structure(list(table = cls, results = results, methods = methods, n = results[[1]]$settings$n),
            class = "causcomp")
}

#' @export
print.causcomp <- function(x, ...) {
  cat("Rcausality method comparison (n = ", x$n, ", ", length(x$methods), " methods)\n", sep = "")
  tb <- x$table
  tb$pair <- paste(tb$var1, "~", tb$var2); tb$var1 <- tb$var2 <- NULL
  tb$agreement <- sprintf("%.0f%%", 100 * tb$agreement)
  print(tb[, c("pair", x$methods, "consensus", "agreement")], row.names = FALSE, right = FALSE)
  cat("  consensus = majority vote over the methods ('undetermined' on ties)\n")
  invisible(x)
}

#' @rdname compare_methods
#' @param main title.
#' @export
plot.causcomp <- function(x, main = NULL, ...) {
  tb <- x$table; M <- as.matrix(tb[, c(x$methods, "consensus"), drop = FALSE])
  pal <- c(x_to_y = "#1B7A43", y_to_x = "#C62828", bidirectional = "#6A1B9A",
           common_cause = "#E69F00", independent = "#BFD7EA", undetermined = "grey70")
  nr <- nrow(M); nc <- ncol(M)
  op <- graphics::par(mar = c(3.2, 7, 4.6, 4.5), xpd = NA); on.exit(graphics::par(op))
  graphics::plot(NA, xlim = c(0.5, nc + 0.5), ylim = c(0.5, nr + 0.5), axes = FALSE, xlab = "", ylab = "")
  graphics::title(main = if (is.null(main)) "Verdicts of the methods per series pair" else main, font.main = 2, line = 3)
  graphics::mtext(sprintf("n = %d; the last column is the majority vote", x$n), side = 3, line = 1.9, cex = 0.72, col = "grey30")
  for (i in seq_len(nr)) for (j in seq_len(nc)) {
    v <- M[i, j]; if (is.na(v)) v <- "undetermined"
    graphics::rect(j - .5, nr - i + .5, j + .5, nr - i + 1.5, col = pal[v], border = "white", lwd = 2)
    # labels are relative to the pair (var1 = x, var2 = y)
    lab <- c(x_to_y = sprintf("%s -> %s", tb$var1[i], tb$var2[i]), y_to_x = sprintf("%s -> %s", tb$var2[i], tb$var1[i]),
             bidirectional = sprintf("%s <-> %s", tb$var1[i], tb$var2[i]), common_cause = "common\ncause",
             independent = "indep.", undetermined = "?")[v]
    graphics::text(j, nr - i + 1, lab, cex = 0.72, col = if (v %in% c("x_to_y", "y_to_x", "bidirectional")) "white" else "black")
  }
  graphics::axis(3, seq_len(nc), colnames(M), tick = FALSE, line = -0.8, cex.axis = 0.8, font = c(rep(1, nc - 1), 2))
  graphics::axis(2, nr:1, paste(tb$var1, "~", tb$var2), tick = FALSE, las = 1, line = -0.8, cex.axis = 0.8)
  graphics::axis(4, nr:1, sprintf("%.0f%%", 100 * tb$agreement), tick = FALSE, las = 1, line = -0.8, cex.axis = 0.7, col.axis = "grey35")
  graphics::mtext("agreement", side = 4, line = 2.6, cex = 0.7, col = "grey35")
  graphics::legend(x = mean(c(0.5, nc + 0.5)), y = 0.35, xjust = 0.5, yjust = 1, horiz = TRUE, bty = "n", cex = 0.62,
                   fill = pal, border = "white", x.intersp = 0.6,
                   legend = c("x -> y", "y -> x", "x <-> y", "common cause", "indep.", "undetermined"))
  invisible(x)
}

#' Benchmark causality methods on simulated scenarios with known truth
#'
#' Simulates `reps` data sets per scenario with a simulator that takes
#' `(n, scenario, seed)` (e.g. [sim_var()], [sim_logistic()],
#' [sim_logistic_bdc()]), runs every method on the pair and scores the class
#' against the truth.
#'
#' @param methods methods to compare.
#' @param scenarios scenarios of the simulator (the truth labels).
#' @param simulator simulator function `function(n, scenario, seed, ...)`
#'   returning a two-column data frame (further columns are ignored).
#' @param n series length.
#' @param reps replications per scenario.
#' @param args per-method argument lists as in [compare_methods()].
#' @param seed base RNG seed (the data set of replication `r` of a scenario
#'   is simulated with `seed + 1000 * scenario index + r`).
#' @param cores CPU cores: the simulated data sets are analysed in parallel
#'   (every method runs serially inside a replication; see
#'   [Rcausality-parallel]).
#' @param ... passed to every method.
#' @return An object of class `causbench`: `accuracy` (method x scenario),
#'   `confusion` (per method), `runs` (one row per run and method), `time`
#'   (mean seconds per run and method). Methods: `print`, `plot`.
#' @export
#' @examples
#' # Two methods on two scenarios of a linear VAR, three replications each.
#' # The accuracy table shows how often each method recovered the true class;
#' # $runs has every run (including the error message of a failed run).
#' b <- benchmark_methods(c("granger", "bayes_granger"), scenarios = c("x_to_y", "independent"),
#'                        simulator = sim_var, n = 200, reps = 2, seed = 1)
#' b
#' \donttest{
#' # a benchmark with a model-free method, replications in parallel
#' b <- benchmark_methods(c("granger", "te"), scenarios = c("x_to_y", "independent"),
#'                        simulator = sim_var, n = 300, reps = 3, seed = 1, cores = 2, n_perm = 100)
#' plot(b)
#' }
benchmark_methods <- function(methods = c("granger", "cmi", "te"), scenarios = c("x_to_y", "y_to_x",
                              "bidirectional", "common_cause", "independent"), simulator = sim_var,
                              n = 500L, reps = 10L, args = list(), seed = 1L,
                              cores = getOption("Rcausality.cores", 1L), ...) {
  dots <- list(...); dots$cores <- NULL
  design <- expand.grid(rep = seq_len(reps), scenario = scenarios, stringsAsFactors = FALSE)[, c("scenario", "rep")]
  old <- .rng_setup(seed); on.exit(.rng_restore(old), add = TRUE)
  run_one <- function(i) {
    sc <- design$scenario[i]; r <- design$rep[i]
    d <- simulator(n, sc, seed = seed + 1000L * match(sc, scenarios) + r)
    d <- as.data.frame(d)[, 1:2]
    rows <- list()
    for (m in methods) {
      a <- c(list(x = d, method = m), dots, if (!is.null(args[[m]])) args[[m]], list(cores = 1L))
      a <- a[!duplicated(names(a), fromLast = TRUE)]
      t0 <- proc.time()[["elapsed"]]
      err <- NA_character_
      res <- tryCatch(do.call(rcausality, a), error = function(e) { err <<- conditionMessage(e); NULL })
      el <- proc.time()[["elapsed"]] - t0
      # the class of the pair oriented as (first column -> second column) of the simulator
      pred <- if (is.null(res)) NA_character_ else {
        cl <- res$classes$class[1]
        if (res$classes$var1[1] == names(d)[1]) cl else c(x_to_y = "y_to_x", y_to_x = "x_to_y")[cl] %||% cl
      }
      rows[[m]] <- data.frame(scenario = sc, rep = r, method = m, predicted = unname(pred),
                              seconds = el, error = err, stringsAsFactors = FALSE)
    }
    do.call(rbind, rows)
  }
  runs <- do.call(rbind, .plapply(seq_len(nrow(design)), run_one, cores = cores))
  rownames(runs) <- NULL
  runs$correct <- !is.na(runs$predicted) & runs$predicted == runs$scenario   # a failed run counts as wrong
  if (any(!is.na(runs$error)))
    warning(sum(!is.na(runs$error)), " run(s) failed and count as wrong; see $runs$error. First error: ",
            runs$error[which(!is.na(runs$error))[1]], call. = FALSE)
  acc <- tapply(runs$correct, list(runs$method, runs$scenario), mean)
  acc <- acc[methods, scenarios, drop = FALSE]
  conf <- lapply(methods, function(m) {
    R <- runs[runs$method == m, ]
    table(truth = factor(R$scenario, scenarios), predicted = factor(R$predicted, .class_levels))
  }); names(conf) <- methods
  structure(list(accuracy = acc, overall = rowMeans(acc, na.rm = TRUE), confusion = conf, runs = runs,
                 time = tapply(runs$seconds, runs$method, mean)[methods], n = n, reps = reps,
                 simulator = deparse(substitute(simulator))), class = "causbench")
}

#' @export
print.causbench <- function(x, digits = 2, ...) {
  cat(sprintf("Rcausality benchmark: simulator %s, n = %d, %d replications per scenario\n", x$simulator, x$n, x$reps))
  cat("  accuracy (fraction of runs with the true class):\n")
  A <- cbind(round(x$accuracy, digits), overall = round(x$overall, digits), `sec/run` = round(x$time, 2))
  print(A)
  invisible(x)
}

#' @rdname benchmark_methods
#' @param x a `causbench` object.
#' @param main title.
#' @export
plot.causbench <- function(x, main = NULL, ...) {
  A <- x$accuracy
  A <- cbind(A, overall = x$overall)
  op <- graphics::par(mar = c(4.5, 4.5, 4.5, 1)); on.exit(graphics::par(op))
  cols <- grDevices::hcl.colors(nrow(A), "Dark 3")
  bp <- graphics::barplot(A, beside = TRUE, col = cols, ylim = c(0, 1.05), ylab = "accuracy", border = "grey40",
                          cex.names = 0.8, xpd = NA)
  graphics::title(main = if (is.null(main)) sprintf("Benchmark on %s (n = %d, %d reps)", x$simulator, x$n, x$reps) else main,
                  font.main = 2, line = 3)
  graphics::legend(x = mean(range(bp)), y = 1.2, legend = rownames(A), fill = cols, border = "grey40",
                   horiz = TRUE, bty = "n", cex = 0.8, xjust = 0.5, xpd = NA)
  graphics::abline(h = 1 / 5, lty = 3, col = "grey50")
  graphics::mtext("dotted line = chance level for five classes; 'overall' = mean over the scenarios", side = 1, line = 3.2, cex = 0.7, col = "grey40")
  invisible(x)
}
