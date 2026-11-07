# Robustness and placebo analyses of a causres --------------------------------
# A detected link is credible when it survives
#  * placebo data in which the timing between the series is destroyed while
#    every series keeps its own dynamics (independent circular shifts):
#    the number of links found in placebo data estimates the false-positive
#    level of the whole analysis, and the placebo distribution of a link's
#    evidence gives a placebo p-value (Theiler et al. 1992; Kosztyan, Matrai
#    and Kiss 2026);
#  * time reversal: for a genuine directed link the reversed series show the
#    opposite direction, whereas links produced by (non-directional)
#    dependence or by the noise level stay (Haufe et al. 2013; Winkler et al.
#    2016);
#  * other reasonable settings: lag orders, multiplicity corrections,
#    subsamples (contiguous blocks, in the spirit of stability selection),
#    and other methods (triangulation).

# repeat the analysis of `object` on new data (same method and settings);
# cores = 1 and seed = NULL inside, the replicate's random stream is used
.rerun <- function(object, data, ...) {
  a <- object$args
  if (is.null(a)) stop("The causres object has no stored arguments (create it with rcausality()).", call. = FALSE)
  a$cores <- 1L; a$seed <- NULL
  extra <- list(...); for (nm in names(extra)) a[[nm]] <- extra[[nm]]
  blocks <- !setequal(object$causes, object$series) || !setequal(object$effects, object$series)
  xs <- if (blocks) data[, object$causes, drop = FALSE] else data
  ys <- if (blocks) data[, setdiff(object$effects, object$causes), drop = FALSE] else NULL
  if (!is.null(ys) && !ncol(ys)) ys <- NULL
  suppressMessages(do.call(rcausality, c(list(x = xs, y = ys), a)))
}

# evidence of every directed test on a common scale: log10 BF10 for the
# Bayesian variants, -log10 adjusted p when the method has p-values,
# otherwise log10 BF10, otherwise the statistic
.evidence_of <- function(P, bayes = FALSE) {
  if (bayes && any(!is.na(P$bf10))) log10(pmax(P$bf10, 1e-300))
  else if (any(!is.na(P$p_adj))) -log10(pmax(P$p_adj, 1e-300))
  else if (any(!is.na(P$bf10))) log10(pmax(P$bf10, 1e-300))
  else P$statistic
}

#' Placebo test of a causality analysis
#'
#' Repeats the analysis of a [`causres`][causres-class] object on `n_placebo`
#' placebo data sets in which every series is circularly shifted by an
#' independent random offset (at least `min_shift` of the length from both
#' ends). The shift keeps the autocorrelation and the marginal distribution
#' of every series and destroys the lead-lag relations between them, so the
#' placebo runs show what the method finds when there is nothing to find:
#' * the global placebo p-value compares the number of links of the
#'   original analysis with the placebo distribution of the number of links;
#' * the link-level placebo p-value compares the evidence of a link
#'   (`-log10` adjusted p-value, or `log10` Bayes factor) with its placebo
#'   distribution; the placebo rate is the share of placebo runs in which the
#'   link is declared.
#' With `reverse = TRUE` the analysis is also repeated on the time-reversed
#' series and the time-reversed Granger score of Haufe et al. (2013) and
#' Winkler et al. (2016) is reported: the net directed score `s(x -> y) -
#' s(y -> x)` (with `s` the effect size of the method, e.g. Geweke's measure)
#' of the data minus the same net score of the reversed data. Reversal turns
#' a genuine lead-lag relation around, so a positive difference supports the
#' direction `x -> y`, whereas a link produced by symmetric dependence (for
#' instance a common driver observed with noise) does not change sign.
#'
#' @param object a `causres` object from [rcausality()].
#' @param n_placebo number of placebo data sets.
#' @param min_shift smallest shift as a fraction of the series length.
#' @param reverse logical; add the time-reversal check.
#' @param seed optional RNG seed.
#' @param cores CPU cores: the placebo runs are distributed over the cores
#'   (each run is serial); the result does not depend on the number of
#'   cores (see [Rcausality-parallel]).
#' @return An object of class `causplacebo`: `links` (one row per directed
#'   test: observed decision and evidence, placebo rate, placebo p-value and,
#'   with `reverse`, the net scores and their difference `trgc`), `counts` (the
#'   number of links in every placebo run), `observed` (number of links of the
#'   original analysis), `p_global`, `settings`. Methods: `print`, `plot`.
#' @references
#' Theiler, J., Eubank, S., Longtin, A., Galdrikian, B. and Farmer, J. D. (1992). Testing for nonlinearity in time series: the method of surrogate data. *Physica D*, 58(1-4), 77-94. \doi{10.1016/0167-2789(92)90102-S}
#'
#' Haufe, S., Nikulin, V. V., Mueller, K.-R. and Nolte, G. (2013). A critical assessment of connectivity measures for EEG data: a simulation study. *NeuroImage*, 64, 120-133. \doi{10.1016/j.neuroimage.2012.09.036}
#'
#' Winkler, I., Panknin, D., Bartz, D., Mueller, K.-R. and Haufe, S. (2016). Validity of time reversal for testing Granger causality. *IEEE Transactions on Signal Processing*, 64(11), 2746-2760. \doi{10.1109/TSP.2016.2531628}
#'
#' Kosztyan, Z. T., Matrai, R. and Kiss, D. (2026). Is the world getting smaller? Network-driven heterogeneity in gravity-model elasticities. *The World Economy*. \doi{10.1111/twec.70148}
#' @export
#' @examples
#' # Granger analysis of three series and its placebo distribution: in data
#' # with the lead-lag structure destroyed the analysis finds (almost) no
#' # link, so the detected links are not an artefact of autocorrelation
#' d <- sim_var(300, n_series = 3, seed = 2)
#' r <- rcausality(d)
#' pl <- placebo_test(r, n_placebo = 19, seed = 1)
#' pl
#' \donttest{
#' plot(pl)
#' placebo_test(r, n_placebo = 99, reverse = TRUE, seed = 1, cores = 2)
#' }
placebo_test <- function(object, n_placebo = 99L, min_shift = 0.1, reverse = FALSE, seed = NULL,
                         cores = getOption("Rcausality.cores", 1L)) {
  if (!inherits(object, "causres")) stop("`object` must be a causres object.", call. = FALSE)
  .check_number(n_placebo, "n_placebo", lower = 1, integer = TRUE)
  .check_number(min_shift, "min_shift", lower = 0, upper = 0.5, strict = TRUE)
  .local_seed(seed)
  D <- as.data.frame(object$data); n <- nrow(D); ms <- max(1L, floor(min_shift * n))
  P <- object$pairs; key <- paste(P$cause, "->", P$effect)
  is_bayes <- function(r) isTRUE(r$settings$bayes) || identical(r$settings$method, "bayes_granger")
  obs_ev <- .evidence_of(P, is_bayes(object))
  runs <- .plapply(seq_len(n_placebo), function(b) {
    Ds <- D; Ds[] <- lapply(D, function(v) .circ_shift(v, ms))
    rb <- .rerun(object, Ds)
    Pb <- rb$pairs; k <- match(key, paste(Pb$cause, "->", Pb$effect))
    list(sig = Pb$significant[k] %in% TRUE, ev = .evidence_of(Pb, is_bayes(rb))[k])
  }, cores = cores)
  S <- vapply(runs, `[[`, logical(length(key)), "sig"); E <- vapply(runs, `[[`, numeric(length(key)), "ev")
  if (is.null(dim(S))) { S <- matrix(S, nrow = 1L); E <- matrix(E, nrow = 1L) }
  counts <- colSums(S); n_obs <- sum(P$significant %in% TRUE)
  links <- data.frame(link = key, cause = P$cause, effect = P$effect, significant = P$significant %in% TRUE,
                      evidence = obs_ev, placebo_rate = rowMeans(S),
                      p_placebo = (1 + rowSums(E >= obs_ev, na.rm = TRUE)) / (1 + rowSums(!is.na(E))),
                      stringsAsFactors = FALSE)
  if (reverse) {
    # time-reversed Granger causality (Haufe et al. 2013; Winkler et al. 2016):
    # the net directed score s(x->y) - s(y->x) of the data minus that of the
    # time-reversed data; a genuine x -> y has a positive difference
    Dr <- D[rev(seq_len(n)), , drop = FALSE]; rownames(Dr) <- NULL
    rr <- .rerun(object, Dr)
    Pr <- rr$pairs; kr <- paste(Pr$cause, "->", Pr$effect)
    score <- function(Q, k) { v <- Q$effect_size; if (all(is.na(v))) v <- Q$statistic; v[match(k, paste(Q$cause, "->", Q$effect))] }
    rev_key <- paste(P$effect, "->", P$cause)
    net_f <- score(P, key) - score(P, rev_key); net_r <- score(Pr, key) - score(Pr, rev_key)
    links$net_forward <- net_f; links$net_reversed <- net_r; links$trgc <- net_f - net_r
    links$reversed_significant <- Pr$significant[match(key, kr)] %in% TRUE
  }
  structure(list(links = links, counts = counts, observed = n_obs,
                 p_global = (1 + sum(counts >= n_obs)) / (1 + n_placebo),
                 settings = list(method = object$settings$method, n_placebo = n_placebo, min_shift = min_shift,
                                 reverse = reverse, n = n)), class = "causplacebo")
}

#' @export
print.causplacebo <- function(x, ...) {
  s <- x$settings
  cat(sprintf("Placebo test of the '%s' analysis: %d placebo data sets (independent circular shifts, at least %.0f%% of n)\n",
              s$method, s$n_placebo, 100 * s$min_shift))
  cat(sprintf("  links found: %d in the data; in the placebo data mean %.2f, 95%% quantile %g; global placebo p = %s\n",
              x$observed, mean(x$counts), stats::quantile(x$counts, 0.95), .fmt_p(x$p_global)))
  L <- x$links
  tb <- data.frame(link = L$link, significant = ifelse(L$significant, "yes", "no"), evidence = round(L$evidence, 2),
                   placebo.rate = sprintf("%.2f", L$placebo_rate), p.placebo = .fmt_p(L$p_placebo), stringsAsFactors = FALSE)
  if (isTRUE(s$reverse)) {
    tb$TRGC <- ifelse(is.na(L$trgc), "-", sprintf("%+.3g", L$trgc))
    tb$time.reversal <- ifelse(is.na(L$trgc), "-", ifelse(L$trgc > 0, "supports the direction", "does not support"))
  }
  print(tb, row.names = FALSE, right = FALSE)
  cat("  evidence = -log10 adjusted p (or log10 BF10); placebo.rate = share of placebo runs declaring the link;\n",
      " p.placebo = share of placebo runs with at least the observed evidence", if (isTRUE(s$reverse))
        paste0(";\n TRGC = net score s(cause->effect) - s(effect->cause) in the data minus the same in the time-reversed data",
               "\n (s = effect size of the method; positive: the temporal order supports the direction)") else "", "\n", sep = "")
  invisible(x)
}

#' @rdname placebo_test
#' @param x a `causplacebo` object.
#' @param main optional title.
#' @param ... unused.
#' @export
plot.causplacebo <- function(x, main = NULL, ...) {
  op <- graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3, 1)); on.exit(graphics::par(op))
  tb <- table(factor(x$counts, levels = 0:max(c(x$counts, x$observed))))
  graphics::barplot(tb, col = "grey80", border = .pal$border, xlab = "number of links in a placebo run", ylab = "placebo runs",
                    main = if (is.null(main)) sprintf("Placebo distribution (global p = %s)", .fmt_p(x$p_global)) else main,
                    font.main = 2, cex.main = 0.9)
  bx <- graphics::barplot(tb, plot = FALSE)
  graphics::abline(v = bx[x$observed + 1], col = .pal$neg, lwd = 2)
  graphics::mtext("red line: links in the data", side = 3, line = -1, adj = 0.98, cex = 0.7, col = .pal$neg)
  L <- x$links[order(x$links$evidence), ]
  graphics::par(mar = c(4.2, 6, 3, 1))
  xr <- range(0, L$evidence, na.rm = TRUE); xr[2] <- xr[2] + 0.3 * diff(xr)
  graphics::plot(L$evidence, seq_len(nrow(L)), xlim = xr, yaxt = "n", ylab = "",
                 xlab = "evidence (-log10 adj. p or log10 BF)", pch = 21, bg = ifelse(L$significant, .pal$causal, "grey85"),
                 col = .pal$border, cex = 1.3, main = "Link evidence and placebo p-value", font.main = 2, cex.main = 0.9)
  graphics::axis(2, at = seq_len(nrow(L)), labels = L$link, las = 1, cex.axis = 0.8)
  graphics::text(L$evidence, seq_len(nrow(L)), sprintf("  p.plac %s", .fmt_p(L$p_placebo)), adj = 0, cex = 0.65, col = "grey30", xpd = NA)
  invisible(x)
}

#' Robustness of a causality analysis to its settings
#'
#' Repeats the analysis of a [`causres`][causres-class] object under
#' alternative settings and reports, for every directed link, in how many of
#' them it is declared:
#' * `lags`: fixed lag orders instead of the selected one;
#' * `adjust`: other multiplicity corrections of the same p-values (e.g.
#'   Benjamini-Hochberg FDR control instead of Holm's family-wise control;
#'   no re-run needed);
#' * `subsamples`: `n_sub` contiguous blocks of `fraction` of the
#'   observations at random positions (a block version of stability
#'   selection, Meinshausen and Buehlmann 2010; contiguous blocks keep the
#'   time-series structure);
#' * `methods`: other causality methods on the same data (triangulation).
#'
#' @param object a `causres` object from [rcausality()].
#' @param lags integer vector of lag orders (`NULL` = none).
#' @param adjust character vector of corrections for [stats::p.adjust()].
#' @param n_sub number of subsamples (`0` = none); `fraction` their length.
#' @param fraction length of a subsample as a fraction of `n`.
#' @param methods character vector of other methods of [rcausality()] (run
#'   with their default settings and the stored `alpha`, `seed`, `n_perm`).
#' @param seed optional RNG seed.
#' @param cores CPU cores: the alternative analyses run in parallel (see
#'   [Rcausality-parallel]).
#' @return An object of class `causrob`: `grid` (one row per directed link
#'   and setting: `type`, `setting`, `significant`, `evidence`), `summary`
#'   (per link: decision of the original analysis, share of the alternative
#'   settings that declare it, by type), `settings`. Methods: `print`,
#'   `plot`.
#' @references
#' Meinshausen, N. and Buehlmann, P. (2010). Stability selection. *Journal of the Royal Statistical Society B*, 72(4), 417-473. \doi{10.1111/j.1467-9868.2010.00740.x}
#'
#' Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society B*, 57(1), 289-300. \doi{10.1111/j.2517-6161.1995.tb02031.x}
#' @export
#' @examples
#' d <- sim_var(300, n_series = 3, seed = 2)
#' r <- rcausality(d)
#' rb <- robustness(r, lags = 1:3, adjust = c("holm", "BH", "none"), n_sub = 10, seed = 1)
#' rb
#' \donttest{
#' plot(rb)
#' robustness(r, lags = NULL, n_sub = 0, methods = c("toda_yamamoto", "te"), seed = 1)
#' }
robustness <- function(object, lags = 1:3, adjust = c("holm", "BH", "none"), n_sub = 20L, fraction = 0.7,
                       methods = NULL, seed = NULL, cores = getOption("Rcausality.cores", 1L)) {
  if (!inherits(object, "causres")) stop("`object` must be a causres object.", call. = FALSE)
  .check_number(n_sub, "n_sub", lower = 0, integer = TRUE)
  .check_number(fraction, "fraction", lower = 0.2, upper = 1)
  if (!is.null(methods)) methods <- match.arg(methods, .method_names, several.ok = TRUE)
  .local_seed(seed)
  D <- as.data.frame(object$data); n <- nrow(D); P <- object$pairs; key <- paste(P$cause, "->", P$effect)
  s <- object$settings
  is_bayes <- function(r) isTRUE(r$settings$bayes) || identical(r$settings$method, "bayes_granger")
  jobs <- c(lapply(lags, function(l) list(type = "lag", setting = sprintf("lag %d", l), args = list(lag = as.integer(l)))),
            if (n_sub > 0) lapply(seq_len(n_sub), function(b) list(type = "subsample", setting = sprintf("block %d", b))),
            lapply(methods, function(m) list(type = "method", setting = m, args = list(method = m, bayes = FALSE, method_args = list()))))
  runs <- .plapply(jobs, function(j) {
    if (j$type == "subsample") {
      len <- max(30L, floor(fraction * n)); st <- sample.int(n - len + 1L, 1L)
      rb <- .rerun(object, D[st:(st + len - 1L), , drop = FALSE])
    } else rb <- do.call(.rerun, c(list(object = object, data = D), j$args))
    Pb <- rb$pairs; k <- match(key, paste(Pb$cause, "->", Pb$effect))
    data.frame(link = key, type = j$type, setting = j$setting, significant = Pb$significant[k] %in% TRUE,
               evidence = .evidence_of(Pb, is_bayes(rb))[k], stringsAsFactors = FALSE)
  }, cores = cores)
  adj <- if (any(!is.na(P$p_value)) && !s$method %in% c("dc", "dfe", "hirata", "cd", "df", "markov") && !is_bayes(object))
    lapply(adjust, function(a) {
      pa <- if (a == "none") P$p_value else stats::p.adjust(P$p_value, method = a)
      data.frame(link = key, type = "adjust", setting = a, significant = !is.na(pa) & pa < s$alpha,
                 evidence = -log10(pmax(pa, 1e-300)), stringsAsFactors = FALSE)
    }) else list()
  G <- do.call(rbind, c(list(data.frame(link = key, type = "original", setting = "original",
                                        significant = P$significant %in% TRUE, evidence = .evidence_of(P, is_bayes(object)),
                                        stringsAsFactors = FALSE)), adj, runs))
  rownames(G) <- NULL
  alt <- G[G$type != "original", ]
  sm <- data.frame(link = key, original = P$significant %in% TRUE,
                   share = vapply(key, function(k) mean(alt$significant[alt$link == k]), 0), stringsAsFactors = FALSE)
  for (t in unique(alt$type)) sm[[paste0("share_", t)]] <- vapply(key, function(k) mean(alt$significant[alt$link == k & alt$type == t]), 0)
  rownames(sm) <- NULL
  structure(list(grid = G, summary = sm, settings = list(method = s$method, lags = lags, adjust = adjust, n_sub = n_sub,
                                                          fraction = fraction, methods = methods)), class = "causrob")
}

#' @export
print.causrob <- function(x, ...) {
  s <- x$settings
  cat(sprintf("Robustness of the '%s' analysis: %d alternative settings\n", s$method, sum(x$grid$type != "original") / nrow(x$summary)))
  S <- x$summary
  tb <- data.frame(link = S$link, original = ifelse(S$original, "yes", "no"), all = sprintf("%.0f%%", 100 * S$share),
                   stringsAsFactors = FALSE)
  for (nm in grep("^share_", names(S), value = TRUE)) tb[[sub("^share_", "", nm)]] <- sprintf("%.0f%%", 100 * S[[nm]])
  print(tb, row.names = FALSE, right = FALSE)
  cat("  share of the alternative settings (all / by type) in which the link is declared\n")
  invisible(x)
}

#' @rdname robustness
#' @param x a `causrob` object.
#' @param main optional title.
#' @param ... unused.
#' @export
plot.causrob <- function(x, main = NULL, ...) {
  G <- x$grid; links <- unique(G$link); st <- unique(paste(G$type, G$setting, sep = ": "))
  M <- matrix(NA, length(links), length(st), dimnames = list(links, st))
  for (i in seq_len(nrow(G))) M[G$link[i], paste(G$type[i], G$setting[i], sep = ": ")] <- G$significant[i]
  op <- graphics::par(mar = c(9, 7, 3, 1)); on.exit(graphics::par(op))
  graphics::image(seq_along(st), seq_along(links), t(M) * 1, col = c("grey90", .pal$causal), axes = FALSE, xlab = "", ylab = "",
                  main = if (is.null(main)) sprintf("Robustness of the '%s' links", x$settings$method) else main,
                  font.main = 2, cex.main = 0.95)
  graphics::axis(1, at = seq_along(st), labels = st, las = 2, cex.axis = 0.65)
  graphics::axis(2, at = seq_along(links), labels = links, las = 1, cex.axis = 0.8)
  graphics::abline(v = seq_along(st) + 0.5, col = "white"); graphics::abline(h = seq_along(links) + 0.5, col = "white")
  graphics::box(col = .pal$border)
  invisible(x)
}
