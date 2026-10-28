#' The `causres` class
#'
#' Result object returned by [rcausality()]. A list with
#' * `pairs`: one row per tested ordered pair (cause -> effect) with lag,
#'   test statistic, degrees of freedom, p-value, adjusted p-value, Bayes
#'   factor, posterior probability, evidence label, effect size (Geweke
#'   F-measure and R-squared gain), sign of the net influence and the
#'   `significant` decision;
#' * `classes`: one row per unordered pair with the 5-class verdict
#'   (`x_to_y`, `y_to_x`, `bidirectional`, `independent`, and -- for methods
#'   able to detect it -- `common_cause`), with evidence for each direction;
#' * `instant` (contemporaneous-dependence table, `cmi_lag`/`markov`),
#'   `markov` (Markov-property table, `markov`), `by_lag` (statistic and
#'   p-value of every directed pair at every candidate lag);
#' * `series`, `causes`, `effects`, `settings`, `call`, `data`.
#'
#' Methods: `print`, `summary`, `plot`, `as.data.frame`, `as.matrix`,
#' `coef`-style accessors via [as.matrix.causres()].
#'
#' @name causres-class
#' @aliases causres
NULL

.new_causres <- function(pairs, series, causes, effects, method, conditional,
                         alpha, bf_threshold, adjust, n, call, data, instant = NULL, markov = NULL, by_lag = NULL, df = NULL, dc = NULL,
                         ccdh = NULL, ccm = NULL, pid = NULL, prob_threshold = NULL) {
  # unordered-pair classification
  key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "\r")
  up <- unique(data.frame(v1 = pmin(pairs$cause, pairs$effect),
                          v2 = pmax(pairs$cause, pairs$effect),
                          stringsAsFactors = FALSE))
  cls <- lapply(seq_len(nrow(up)), function(i) {
    a <- up$v1[i]; b <- up$v2[i]
    ab <- pairs[pairs$cause == a & pairs$effect == b, , drop = FALSE]
    ba <- pairs[pairs$cause == b & pairs$effect == a, , drop = FALSE]
    s_ab <- if (nrow(ab)) ab$significant else NA
    s_ba <- if (nrow(ba)) ba$significant else NA
    cl <- if (isTRUE(s_ab) && isTRUE(s_ba)) "bidirectional" else
      if (isTRUE(s_ab)) "x_to_y" else if (isTRUE(s_ba)) "y_to_x" else
        if (is.na(s_ab) && is.na(s_ba)) NA_character_ else "independent"
    # one-sided blocks (x -> y only tested): the reverse is not assessed
    tested_both <- nrow(ab) > 0 && nrow(ba) > 0
    cd <- FALSE
    if (!is.null(instant)) {
      ii <- instant[(instant$var1 == a & instant$var2 == b) | (instant$var1 == b & instant$var2 == a), , drop = FALSE]
      cd <- nrow(ii) > 0 && isTRUE(ii$significant)
      if (cd && cl == "independent") cl <- "common_cause"
    }
    data.frame(var1 = a, var2 = b, class = cl, common_driver = cd,
               both_directions_tested = tested_both,
               p_adj_12 = if (nrow(ab)) ab$p_adj else NA_real_,
               p_adj_21 = if (nrow(ba)) ba$p_adj else NA_real_,
               bf10_12 = if (nrow(ab)) ab$bf10 else NA_real_,
               bf10_21 = if (nrow(ba)) ba$bf10 else NA_real_,
               stringsAsFactors = FALSE)
  })
  classes <- do.call(rbind, cls)
  classes$relation <- ifelse(is.na(classes$class), NA_character_,
    ifelse(classes$class == "x_to_y", paste(classes$var1, "->", classes$var2),
    ifelse(classes$class == "y_to_x", paste(classes$var2, "->", classes$var1),
    ifelse(classes$class == "bidirectional", paste(classes$var1, "<->", classes$var2),
    ifelse(classes$class == "common_cause", paste(classes$var1, "<-L->", classes$var2),
           paste(classes$var1, "|", classes$var2))))))
  classes$relation <- ifelse(!is.na(classes$class) & classes$common_driver &
                               classes$class != "common_cause",
                             paste0(classes$relation, " (+ common driver?)"), classes$relation)
  # direction labels are relative to (var1, var2): x = var1, y = var2
  structure(list(pairs = pairs, classes = classes, series = series,
                 causes = causes, effects = effects,
                 settings = list(method = method, conditional = conditional,
                                 alpha = alpha, bf_threshold = bf_threshold,
                                 adjust = adjust, n = n, prob_threshold = prob_threshold),
                 call = call, data = data, instant = instant, markov = markov, by_lag = by_lag, df = df, dc = dc, ccdh = ccdh, ccm = ccm, pid = pid),
            class = "causres")
}

# restrict a result to the directed tests that were asked for (blocks with
# directions = "forward"): the pairwise methods compute both directions of
# every unordered pair; the unrequested direction is dropped and the
# multiplicity correction is recomputed on the requested tests
.restrict_to_grid <- function(r) {
  asked <- paste(rep(r$causes, times = length(r$effects)), rep(r$effects, each = length(r$causes)))
  keep <- paste(r$pairs$cause, r$pairs$effect) %in% asked
  if (all(keep)) return(r)
  P <- r$pairs[keep, , drop = FALSE]; s <- r$settings
  if (!s$method %in% c("dc", "dfe", "hirata", "cd", "df") && any(!is.na(P$p_value))) {
    P$p_adj <- if (s$adjust == "none") P$p_value else stats::p.adjust(P$p_value, method = s$adjust)
    sig <- !is.na(P$p_adj) & P$p_adj < s$alpha
    if (s$method == "bayes_granger") sig <- !is.na(P$bf10) & P$bf10 > s$bf_threshold
    if (s$method == "ccm" && !is.null(r$ccm)) {
      cv <- r$ccm[match(paste(P$cause, P$effect), paste(r$ccm$cause, r$ccm$effect)), "converging"]
      sig <- sig & cv %in% TRUE
    }
    P$significant <- sig
  }
  sel <- function(tab) if (is.null(tab) || !all(c("cause", "effect") %in% names(tab))) tab else
    tab[paste(tab$cause, tab$effect) %in% asked, , drop = FALSE]
  r$by_lag <- sel(r$by_lag); r$ccm <- sel(r$ccm); r$pid <- sel(r$pid)
  rownames(P) <- NULL
  .rebuild_causres(r, P)
}

# rebuild a causres from edited pairs, keeping every other slot
.rebuild_causres <- function(r, pairs) {
  s <- r$settings
  out <- .new_causres(pairs, series = r$series, causes = r$causes, effects = r$effects, method = s$method,
                      conditional = s$conditional, alpha = s$alpha, bf_threshold = s$bf_threshold, adjust = s$adjust,
                      n = s$n, call = r$call, data = r$data, instant = r$instant, markov = r$markov, by_lag = r$by_lag,
                      df = r$df, dc = r$dc, ccdh = r$ccdh, ccm = r$ccm, pid = r$pid, prob_threshold = s$prob_threshold)
  extra <- setdiff(names(r), names(out)); for (e in extra) out[[e]] <- r[[e]]
  more <- setdiff(names(s), names(out$settings)); for (e in more) out$settings[[e]] <- s[[e]]
  out
}

#' @export
as.data.frame.causres <- function(x, ...) x$pairs

#' Matrix view of a `causres` object
#'
#' @param x a `causres` object.
#' @param what one of `"p_adj"`, `"p_value"`, `"bf10"`, `"post_prob"`,
#'   `"lag"`, `"statistic"`, `"effect_size"`, `"significant"`.
#' @param ... unused.
#' @return Matrix with causes in rows and effects in columns
#'   (`NA` where not tested).
#' @export
as.matrix.causres <- function(x, what = "p_adj", ...) {
  P <- x$pairs
  if (!what %in% names(P)) stop("Unknown column: ", what, call. = FALSE)
  M <- matrix(NA, length(x$series), length(x$series),
              dimnames = list(cause = x$series, effect = x$series))
  M[cbind(match(P$cause, x$series), match(P$effect, x$series))] <- P[[what]]
  M
}

.fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 1e-4, "<1e-04", formatC(p, digits = 3, format = "g")))
.fmt_bf <- function(b) ifelse(is.na(b), "NA", ifelse(b > 1e4, ">1e+04", ifelse(b < 1e-4, "<1e-04", formatC(b, digits = 3, format = "g"))))
.stars <- function(p) ifelse(is.na(p), "", ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", ifelse(p < .1, ".", "")))))

#' @export
print.causres <- function(x, digits = 3, ...) {
  s <- x$settings
  mname <- c(granger = "Granger causality (F-test)",
             bayes_granger = "Bayesian Granger causality",
             cmi = "Conditional mutual information (kNN, local permutation)",
             cmi_lag = "Lag-shifted conditional mutual information (Zlatniczki et al., 2021)",
             markov = "Markov-chain approach (Stippinger et al., 2023)",
             df = "Degrees-of-freedom / conditional-variance method (Telcs, Kurbucz & Jakovac, 2025)",
             dc = "Bayesian dimensional causality (Benko et al., 2024)",
             ccdh = "CCDH: complexity-based causal discovery (Stippinger et al., 2023)",
             dfe = "df-entropy method with random-forest classification (Jakovac et al.)",
             te = "Transfer entropy, binned (Schreiber, 2000) with surrogate test",
             ccm = "Convergent cross mapping (Sugihara et al., 2012)",
             pid = "MAR-PID: unique + synergistic information of the cause's past about the effect's future (Telcs, 2026)",
             hirata = "Recurrence-plot criteria: coverage and simultaneous recurrences (Hirata & Aihara, 2010)",
             cd = "Correlation-dimension relations (Krakovska, 2019)",
             pcmci = "PCMCI / PCMCI+ via tigramite (Runge et al., 2019)")[[s$method]]
  cat("\n  Rcausality: ", mname, "\n", sep = "")
  cat("  ", length(x$series), " series, n = ", s$n, ", ",
      nrow(x$pairs), " directed tests",
      if (s$conditional) " (conditional on the other series)" else " (pairwise)",
      "\n", sep = "")
  bc <- s$bayes_control
  bayes_rule <- function() {
    dec <- if (!is.null(bc) && identical(bc$decision, "posterior")) sprintf("posterior probability > %g", s$prob_threshold)
           else sprintf("BF10 > %g", s$bf_threshold)
    how <- if (s$method %in% c("granger", "bayes_granger"))
      switch(if (is.null(s$prior)) "bic" else s$prior, bic = "BIC approximation", g = "Zellner g-prior", hyper_g = "hyper-g prior",
             spike_slab = sprintf("spike-and-slab regression, slab scale %g, inclusion prior %g", bc$scale, bc$inclusion_prob),
             normal = sprintf("Bayesian regression, normal slab of scale %g", bc$scale))
      else if (s$method %in% c("dc", "dfe")) "posterior of the causal cases"
      else if (s$method == "markov") "Bayesian decision table from the surrogate Bayes factors"
      else if (s$method == "pcmci") "Savage-Dickey Bayes factor of the partial correlation"
      else sprintf("Bayesian surrogate test, %s effect prior of scale %g, %s null", bc$effect_prior, bc$scale, bc$null_model)
    sprintf("%s (%s)", dec, how)
  }
  cat("  decision: ",
      if (isTRUE(s$bayes) && !is.null(bc)) bayes_rule()
      else if (s$method %in% c("dc", "dfe")) sprintf("class score P(cause drives effect) > %g", s$prob_threshold)
      else if (s$method == "ccdh") sprintf("shared information: adjusted p < %g (%s); reconstruction index <= %g", s$alpha, s$adjust, s$prob_threshold)
      else if (s$method == "hirata") sprintf("coverage index >= threshold and surrogate p < %g; verdict per pair", s$alpha)
      else if (s$method == "cd") sprintf("bootstrap z-tests of the dimension relations at level %g; verdict per pair", s$alpha)
      else if (s$method == "pcmci") sprintf("MCI test '%s'%s, adjusted p < %g (%s)", s$test, if (isTRUE(s$plus)) " (PCMCI+)" else "", s$alpha, s$adjust)
      else if (s$method == "df") sprintf("Table-2 reading of the orders where determinate, else directed variance-reduction tests with adjusted p < %g (%s)", s$alpha, s$adjust)
      else if (s$method != "bayes_granger") sprintf("adjusted p < %g (%s)", s$alpha, s$adjust)
      else sprintf("BF10 (%s) > %g", switch(if (is.null(s$prior)) "bic" else s$prior, bic = "BIC approximation",
                                            g = "Zellner g-prior, g = n", hyper_g = "hyper-g prior, a = 3"), s$bf_threshold),
      "\n\n", sep = "")
  P <- x$pairs
  tab <- data.frame(
    link = paste(P$cause, "->", P$effect),
    lag = P$lag,
    F = formatC(P$statistic, digits = digits, format = "f"),
    df = ifelse(is.na(P$df1), "-", paste0(P$df1, ",", P$df2)),
    p = .fmt_p(P$p_value),
    p.adj = .fmt_p(P$p_adj),
    BF10 = .fmt_bf(P$bf10),
    P.post = formatC(P$post_prob, digits = digits, format = "f"),
    sig = ifelse(P$significant, paste0("yes ", .stars(P$p_adj)), "no"),
    check.names = FALSE, stringsAsFactors = FALSE)
  cmi_like <- s$method %in% c("cmi", "cmi_lag", "markov", "df", "ccdh", "te", "ccm", "pid", "hirata", "cd", "pcmci")
  names(tab)[names(tab) == "F"] <- switch(s$method, df = "sd.red", dc = , dfe = "P.direct", ccdh = "recon",
                                          te = "TE.bits", ccm = "rho", pid = "U+S.bits", hirata = "cover.idx",
                                          cd = "D2.diff", pcmci = "MCI.stat", if (cmi_like) "CMI" else "F")
  if (cmi_like) tab <- tab[, !names(tab) %in% c("df", if (!isTRUE(s$bayes)) c("BF10", "P.post"))]
  if (s$method %in% c("dc", "dfe")) {
    tab <- tab[, !names(tab) %in% c("df", "p", "p.adj", if (s$method == "dc") "lag")]
    tab$sig <- ifelse(P$significant, "yes", "no")
  }
  print(tab, row.names = FALSE, right = FALSE)
  if (s$method %in% c("dc", "dfe")) {
    cat("  P.direct = score of the unidirectional link, P.post = P(direct) + P(bidirectional), BF10 against the uniform prior\n")
    cat(if (s$method == "dc") "\n  Intrinsic dimensions and posterior probabilities of the five causal cases:\n"
        else "\n  Random-forest class scores of the five causal cases:\n")
    Dc <- x$dc
    tb <- data.frame(pair = paste(Dc$var1, "~", Dc$var2))
    if (s$method == "dc") tb <- cbind(tb, D.X = round(Dc$D_X, 2), D.Y = round(Dc$D_Y, 2), D.J = round(Dc$D_J, 2), D.Z = round(Dc$D_Z, 2))
    tb <- cbind(tb, `x->y` = round(Dc$x_to_y, 3), `y->x` = round(Dc$y_to_x, 3), `<->` = round(Dc$bidirectional, 3),
                common = round(Dc$common_cause, 3), indep = round(Dc$independent, 3), best = Dc$best)
    print(tb, row.names = FALSE, right = FALSE)
  }
  if (s$method == "df") cat("  sd.red = relative reduction of the conditional sd of the effect when the cause's past is constrained\n")
  if (s$method == "ccm") {
    tab$lag <- NULL
    cat("  rho = cross-map skill of 'effect xmap cause' at the largest library; p = surrogate test;\n")
    cat("  a link also requires convergence (Kendall test of rho against the library size):\n")
    Cv <- x$ccm
    print(data.frame(link = paste(Cv$cause, "->", Cv$effect), converging = ifelse(Cv$converging, "yes", "no"),
                     p.convergence = .fmt_p(Cv$p_convergence)), row.names = FALSE, right = FALSE)
  }
  if (s$method == "pid") {
    cat("  U+S.bits = unique + synergistic information (bits per target bit) of the cause's past about the effect's future beyond the effect's own past\n")
    cat("\n  MAR-PID atoms at the selected lag (bits per target bit):\n")
    Pd <- x$pid
    print(data.frame(link = paste(Pd$cause, "->", Pd$effect), lag = Pd$lag, unique.cause = round(Pd$unique_cause, 3),
                     unique.effect = round(Pd$unique_effect, 3), redundant = round(Pd$redundant, 3),
                     synergy = round(Pd$synergy, 3), total = round(Pd$total, 3)), row.names = FALSE, right = FALSE)
  }
  if (s$method == "te") cat("  TE.bits = transfer entropy at the selected lag (bits); effect_size column of $pairs = effective TE\n")
  if (s$method %in% c("hirata", "cd")) {
    if (s$method == "hirata") cat("  cover.idx = normalised coverage P(R_cause | R_effect); p = surrogate test of the coverage\n")
    else cat("  D2.diff = D2(effect) - D2(cause); no p-values, the verdict comes from the z-tests of the dimension relations\n")
    cat("\n  Per-pair quantities and verdict:\n")
    Cc <- x$ccdh; Cc$pair <- paste(Cc$var1, "~", Cc$var2); Cc$var1 <- Cc$var2 <- NULL
    num <- vapply(Cc, is.numeric, TRUE); Cc[num] <- lapply(Cc[num], round, 3)
    print(Cc[, c("pair", setdiff(names(Cc), "pair"))], row.names = FALSE, right = FALSE)
  }
  if (s$method == "ccdh") {
    cat("  recon = 1 - reconstruction index (1: the cause is perfectly reconstructible from the effect's embedding); p = shared-information surrogate test\n")
    cat("\n  Mean compression rates of the series and of their skewed joints:\n")
    Cc <- x$ccdh
    print(data.frame(pair = paste(Cc$var1, "~", Cc$var2), cr.x = round(Cc$cr_x, 3), cr.y = round(Cc$cr_y, 3),
                     cr.Jxxy = round(Cc$cr_jxxy, 3), cr.Jyyx = round(Cc$cr_jyyx, 3),
                     index.x = round(Cc$index_x, 3), index.y = round(Cc$index_y, 3), verdict = Cc$verdict,
                     check.names = FALSE), row.names = FALSE, right = FALSE)
  }
  if (!is.null(x$instant) && s$method == "pcmci") {
    cat("\n  Contemporaneous links (lag 0):\n")
    I <- x$instant
    print(data.frame(pair = paste(I$var1, "--", I$var2),
                     MCI.stat = formatC(I$statistic, digits = digits, format = "f"),
                     p = .fmt_p(I$p_value), p.adj = .fmt_p(I$p_adj),
                     significant = ifelse(I$significant, "yes", "no"),
                     check.names = FALSE), row.names = FALSE, right = FALSE)
  }
  if (!is.null(x$instant) && !s$method %in% c("df", "dc", "ccdh", "dfe", "hirata", "cd", "pcmci")) {
    cat("\n  Contemporaneous dependence given both pasts (k = 0; common-driver evidence):\n")
    I <- x$instant
    print(data.frame(pair = paste(I$var1, "~", I$var2),
                     CMI = formatC(I$cmi, digits = digits, format = "f"),
                     p = .fmt_p(I$p_value), p.adj = .fmt_p(I$p_adj),
                     common.driver = ifelse(I$significant, "yes", "no"),
                     check.names = FALSE), row.names = FALSE, right = FALSE)
  }
  if (!is.null(x$df)) {
    cat("\n  Orders (degrees of freedom = conditions on the observed pasts) and the Table-2 reading:\n")
    Dd <- x$df
    print(data.frame(pair = paste(Dd$var1, "~", Dd$var2), df.x = Dd$df_x, df.y = Dd$df_y,
                     df.joint = Dd$df_joint, conditions.joint = Dd$conditions_joint,
                     `Table 2` = Dd$table2, decided.by = Dd$decided_by, check.names = FALSE), row.names = FALSE, right = FALSE)
    cat("  a determinate Table-2 reading decides the pair; undetermined pairs are decided by the directed tests above\n")
  }
  if (!is.null(x$markov)) {
    cat("\n  Markov property (p small = not Markov) and decision-table case:\n")
    M <- x$markov
    print(data.frame(pair = paste(M$var1, "~", M$var2), m = M$memory,
                     p.x = .fmt_p(M$p_markov_x), p.y = .fmt_p(M$p_markov_y),
                     p.joint = .fmt_p(M$p_markov_joint),
                     case = M$case, `x->y` = M$inj_x_to_y, `y->x` = M$inj_y_to_x,
                     hidden = M$hidden, consistent = M$consistent, check.names = FALSE),
          row.names = FALSE, right = FALSE)
    cat("  injection: i.i. = i.i.d. injection, n.i. = non-i.i.d., d = driving; U_* = hidden own-source\n")
  }
  if (!is.null(x$bayes)) {
    B <- x$bayes
    cat("\n  Bayesian evidence per directed test (effect = posterior mean [credible interval] under H1):\n")
    tb <- data.frame(link = B$pair, BF10 = .fmt_bf(B$bf10), P.post = formatC(B$post_prob, digits = 3, format = "f"),
                     evidence = as.character(bf_evidence(B$bf10)),
                     effect = ifelse(is.na(B$effect_mean), "-", sprintf("%.3g [%.3g, %.3g]", B$effect_mean, B$effect_lower, B$effect_upper)),
                     model = B$model, check.names = FALSE, stringsAsFactors = FALSE)
    if (!is.null(B$rhat) && any(!is.na(B$rhat)))
      tb$`R-hat` <- ifelse(is.na(B$rhat), "-", formatC(B$rhat, digits = 3, format = "f"))
    print(tb, row.names = FALSE, right = FALSE)
    if (!is.null(B$rhat) && any(!is.na(B$rhat)))
      cat("  R-hat = Gelman-Rubin factor of the effect over the chains (close to 1: the chains agree)\n")
  }
  cat("\n  Classification of series pairs:\n")
  C <- x$classes
  ct <- data.frame(pair = paste(C$var1, "~", C$var2), class = C$class,
                   relation = C$relation,
                   stringsAsFactors = FALSE)
  print(ct, row.names = FALSE, right = FALSE)
  if (any(!C$both_directions_tested))
    cat("  (pairs with only one tested direction are labelled by that direction alone)\n")
  cat("  How to read: one row per directed test cause -> effect; 'sig' applies the decision rule above;\n",
      " the classification combines the two directions of a pair (x -> y, y -> x, x <-> y, independent;\n",
      " 'common_cause' = no directed link but evidence of a hidden common driver). plot(x, \"all\") shows the same.\n", sep = "")
  invisible(x)
}

#' Summary of a `causres` object
#'
#' @param object a `causres` object.
#' @param ... unused.
#' @return An object of class `summary.causres` (print it).
#' @references
#' Granger, C. W. J. (1969). Investigating causal relations by econometric models and cross-spectral methods. *Econometrica*, 37(3), 424-438. \doi{10.2307/1912791}
#'
#' Geweke, J. (1982). Measurement of linear dependence and feedback between multiple time series. *Journal of the American Statistical Association*, 77(378), 304-313. \doi{10.1080/01621459.1982.10477803}
#'
#' Holm, S. (1979). A simple sequentially rejective multiple test procedure. *Scandinavian Journal of Statistics*, 6(2), 65-70.
#' @export
summary.causres <- function(object, ...) {
  P <- object$pairs; s <- object$settings
  sig <- P[P$significant, , drop = FALSE]
  sig <- sig[order(-sig$effect_size), , drop = FALSE]
  vars <- object$series
  deg <- data.frame(series = vars,
                    out_degree = vapply(vars, function(v) sum(sig$cause == v), 0L),
                    in_degree = vapply(vars, function(v) sum(sig$effect == v), 0L),
                    stringsAsFactors = FALSE)
  deg$role <- ifelse(deg$out_degree > 0 & deg$in_degree == 0, "driver",
               ifelse(deg$in_degree > 0 & deg$out_degree == 0, "target",
               ifelse(deg$in_degree > 0 & deg$out_degree > 0, "hub/mediator", "isolated")))
  lagdist <- if (nrow(sig)) table(sig$lag) else table(integer(0))
  structure(list(settings = s, n_tests = nrow(P), n_sig = nrow(sig),
                 density = nrow(sig) / nrow(P), links = sig, degree = deg,
                 classes = table(factor(object$classes$class,
                                        levels = c("x_to_y", "y_to_x", "bidirectional",
                                                   "common_cause", "independent"))),
                 lags = lagdist, series = vars, call = object$call),
            class = "summary.causres")
}

#' @export
print.summary.causres <- function(x, digits = 3, ...) {
  s <- x$settings
  cat("\nRcausality summary\n", strrep("=", 50), "\n", sep = "")
  cat("Call:   ", deparse(x$call), "\n", sep = "")
  cat("Method: ", s$method, if (s$conditional) " (conditional)" else " (pairwise)",
      "   n = ", s$n, "   alpha = ", s$alpha, " (", s$adjust, ")\n", sep = "")
  cat("Tests:  ", x$n_tests, "   significant links: ", x$n_sig,
      sprintf("   (density %.0f%%)\n", 100 * x$density), sep = "")
  cat("\nPair classes (x = first, y = second variable of the pair):\n")
  print(x$classes)
  if (x$n_sig) {
    cat("\nSignificant links, strongest first:\n")
    L <- x$links
    tab <- data.frame(link = paste(L$cause, "->", L$effect), lag = L$lag,
                      sign = ifelse(L$sign > 0, "+", ifelse(L$sign < 0, "-", "0")),
                      F = formatC(L$statistic, digits = digits, format = "f"),
                      p.adj = .fmt_p(L$p_adj), BF10 = .fmt_bf(L$bf10),
                      P.post = formatC(L$post_prob, digits = digits, format = "f"),
                      R2.gain = formatC(L$r2_gain, digits = digits, format = "f"),
                      evidence = L$evidence, check.names = FALSE,
                      stringsAsFactors = FALSE)
    print(tab, row.names = FALSE, right = FALSE)
    cat("\nNode roles:\n")
    print(x$degree, row.names = FALSE, right = FALSE)
    cat("\nLag distribution of significant links:\n")
    print(x$lags)
  } else cat("\nNo significant causal links at the chosen threshold.\n")
  invisible(x)
}
