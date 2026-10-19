#' Detect causal relations between time series
#'
#' Single entry point of the package. Works with
#' * two series: `rcausality(x, y)`;
#' * one data frame / matrix: `rcausality(df)` -- all ordered pairs;
#' * two blocks (MIMO): `rcausality(X, Y)` with several columns each --
#'   every column of `X` is tested as a cause of every column of `Y`.
#'
#' @param x numeric vector, matrix or data frame (the putative causes).
#' @param y optional numeric vector, matrix or data frame (the putative
#'   effects). If `NULL` and `x` has several columns, all ordered pairs of the
#'   columns of `x` are tested.
#' @param method `"granger"` (F-test, default) or `"bayes_granger"`
#'   (decision by BIC-approximated Bayes factor). Both report p-value *and*
#'   Bayes factor.
#' @param lag integer lag order or `"auto"` (selected by `lag_criterion` on
#'   the VAR of the involved series).
#' @param max_lag largest lag order considered when `lag = "auto"`.
#' @param lag_criterion `"BIC"`, `"AIC"` or `"HQ"`.
#' @param alpha significance level (applied to adjusted p-values).
#' @param bf_threshold Bayes factor above which a link is declared
#'   (used for `method = "bayes_granger"`).
#' @param adjust multiple-testing correction passed to
#'   [stats::p.adjust()] (default `"holm"`; `"none"` disables it).
#' @param conditional logical; if `TRUE` every test conditions on all other
#'   series in the data (full conditional / multivariate Granger), otherwise
#'   tests are pairwise.
#' @param directions for two series / two blocks: `"both"` (default) tests
#'   x -> y and y -> x; `"forward"` only the causes-in-`x`, effects-in-`y`
#'   direction. Ignored when `y` is `NULL` (all ordered pairs are tested).
#' @param prior_prob prior probability of a causal link (for the posterior).
#' @param impute logical; interpolate missing values?
#' @return An object of class [`causres`][causres-class].
#' @export
#' @examples
#' d <- sim_var(400, scenario = "x_to_y", seed = 1)
#' r <- rcausality(d$x, d$y)
#' r
#' summary(r)
#'
#' # all pairs of a data frame
#' d3 <- sim_var(400, n_series = 4, seed = 2)
#' r3 <- rcausality(d3)
#' plot(r3)
rcausality <- function(x, y = NULL, method = c("granger", "bayes_granger"),
                       lag = "auto", max_lag = 8L,
                       lag_criterion = c("BIC", "AIC", "HQ"),
                       alpha = 0.05, bf_threshold = 3,
                       adjust = "holm", conditional = FALSE,
                       prior_prob = 0.5, impute = TRUE,
                       directions = c("both", "forward")) {
  cl <- match.call()
  method <- match.arg(method)
  lag_criterion <- match.arg(lag_criterion)
  directions <- match.arg(directions)
  fr <- .as_series_frame(x, y, xname = if (is.null(y)) "x" else "x", yname = "y")
  if (is.null(fr$Y) && ncol(fr$X) < 2L) {
    stop("Supply two series (x, y) or a data frame with at least two columns.",
         call. = FALSE)
  }
  X <- fr$X; Y <- fr$Y
  if (!is.null(Y) && nrow(X) != nrow(Y)) stop("x and y must have the same length.", call. = FALSE)
  # make unique names across blocks
  all <- if (is.null(Y)) X else cbind(X, Y)
  if (anyDuplicated(names(all))) {
    if (!is.null(Y)) {
      names(X) <- paste0("x.", names(X)); names(Y) <- paste0("y.", names(Y))
      all <- cbind(X, Y)
    } else stop("Column names must be unique.", call. = FALSE)
  }
  all[] <- lapply(all, function(v) { v <- as.numeric(v); if (impute) .prep_na(v) else v })
  if (any(!is.finite(as.matrix(all)))) stop("Missing/non-finite values present (use impute = TRUE).", call. = FALSE)
  causes <- if (is.null(Y)) names(all) else names(X)
  effects <- if (is.null(Y)) names(all) else names(Y)
  grid <- expand.grid(cause = causes, effect = effects, stringsAsFactors = FALSE)
  if (!is.null(Y) && directions == "both") {
    grid <- rbind(grid, setNames(grid[, 2:1], c("cause", "effect")))
    causes <- effects <- names(all)
  }
  grid <- grid[grid$cause != grid$effect, , drop = FALSE]
  if (nrow(grid) == 0L) stop("No pairs to test.", call. = FALSE)

  fixed_lag <- if (identical(lag, "auto")) NULL else as.integer(lag)
  rows <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    cs <- grid$cause[i]; ef <- grid$effect[i]
    ctrl_names <- if (conditional) setdiff(names(all), c(cs, ef)) else character(0)
    inv <- all[, c(ef, cs, ctrl_names), drop = FALSE]
    p <- if (is.null(fixed_lag)) select_lag(inv, max_lag, lag_criterion) else fixed_lag
    p <- max(1L, min(p, floor((nrow(all) - 2) / (2 + 2 * ncol(inv)))))
    ctrl <- if (length(ctrl_names)) as.matrix(all[, ctrl_names, drop = FALSE]) else NULL
    g <- .granger_pair(all[[ef]], all[[cs]], ctrl, p)
    rows[[i]] <- data.frame(cause = cs, effect = ef, lag = p,
                            statistic = g$statistic, df1 = g$df1, df2 = g$df2,
                            p_value = g$p_value, bf10 = g$bf10,
                            effect_size = g$effect_size, r2_gain = g$r2_gain,
                            sign = g$sign, n_eff = g$n_eff,
                            stringsAsFactors = FALSE)
  }
  pairs <- do.call(rbind, rows)
  pairs$p_adj <- if (adjust == "none") pairs$p_value else
    stats::p.adjust(pairs$p_value, method = adjust)
  pairs$post_prob <- (pairs$bf10 * prior_prob) /
    (pairs$bf10 * prior_prob + (1 - prior_prob))
  pairs$evidence <- as.character(bf_evidence(pairs$bf10))
  pairs$significant <- if (method == "granger") {
    !is.na(pairs$p_adj) & pairs$p_adj < alpha
  } else {
    !is.na(pairs$bf10) & pairs$bf10 > bf_threshold
  }
  pairs <- pairs[, c("cause", "effect", "lag", "statistic", "df1", "df2", "p_value",
                     "p_adj", "bf10", "post_prob", "evidence", "effect_size",
                     "r2_gain", "sign", "n_eff", "significant")]
  rownames(pairs) <- NULL
  .new_causres(pairs, series = names(all), causes = causes, effects = effects,
               method = method, conditional = conditional, alpha = alpha,
               bf_threshold = bf_threshold, adjust = adjust, n = nrow(all),
               call = cl, data = all)
}
