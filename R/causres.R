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
#' The explanatory lines of the `print` methods are wrapped to the width of
#' the console (`getOption("width")`); `options(Rcausality.guide = FALSE)`
#' hides the reading guide that `print` appends to a result.
#'
#' Two subclasses add slots and methods: `causpanel` (panel data, see
#' [panel_causality()]: `units`, `panel`, `prevalence`) and `causdrivers`
#' (the drivers of one explained variable, see [causal_drivers()]: `drivers`,
#' `graph`, `target`). Both keep the slots above, so that every function
#' that takes a `causres` works with them.
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
    ifelse(classes$class == "common_cause", paste(classes$var1, "<~>", classes$var2),
           paste(classes$var1, "|", classes$var2))))))
  classes$relation <- ifelse(!is.na(classes$class) & classes$common_driver &
                               classes$class != "common_cause",
                             paste0(classes$relation, " (+ <~>)"), classes$relation)
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
  info <- .method_info(s$method); mname <- info$label
  .out("\n  Rcausality: ", mname, "\n", sep = "")
  .out("  ", length(x$series), " series, n = ", s$n, ", ",
      nrow(x$pairs), " directed tests",
      if (s$conditional) " (conditional on the other series)" else " (pairwise)",
      "\n", sep = "")
  bc <- s$bayes_control
  g_txt <- function() if (is.null(bc) || isTRUE(all.equal(bc$scale, 1))) "n" else sprintf("%g n", bc$scale)
  bayes_rule <- function() {
    dec <- if (!is.null(bc) && identical(bc$decision, "posterior")) sprintf("posterior probability > %g", s$prob_threshold)
           else sprintf("BF10 > %g", s$bf_threshold)
    how <- if (s$method %in% c("granger", "bayes_granger"))
      switch(if (is.null(s$prior)) "bic" else s$prior, bic = "BIC approximation", g = sprintf("Zellner g-prior, g = %s", g_txt()),
             hyper_g = sprintf("hyper-g prior, a = %g", bc$hyper_a),
             spike_slab = sprintf("spike-and-slab regression, slab scale %g, inclusion prior %g", bc$scale, bc$inclusion_prob),
             normal = sprintf("Bayesian regression, normal slab of scale %g", bc$scale),
             minnesota = sprintf("hierarchical Minnesota prior, tightness %s, lag decay %g",
                                 if (identical(bc$lambda, "auto")) "integrated out" else format(bc$lambda), bc$decay))
      else if (s$method == "toda_yamamoto")
        sprintf("lag-augmented regressions, %s", switch(if (is.null(s$prior)) "bic" else s$prior, g = sprintf("Zellner g-prior, g = %s", g_txt()),
                                                       hyper_g = sprintf("hyper-g prior, a = %g", bc$hyper_a), "BIC approximation"))
      else if (s$method %in% c("dc", "dfe")) "posterior of the causal cases"
      else if (s$method == "markov") "Bayesian decision table from the surrogate Bayes factors"
      else if (s$method == "pcmci") "Savage-Dickey Bayes factor of the partial correlation"
      else sprintf("Bayesian surrogate test, %s effect prior of scale %g, %s null", bc$effect_prior, bc$scale, bc$null_model)
    sprintf("%s (%s)", dec, how)
  }
  .out("  decision: ",
      if (isTRUE(s$bayes) && !is.null(bc)) bayes_rule()
      else if (s$method %in% c("dc", "dfe")) sprintf("class score P(cause drives effect) > %g", s$prob_threshold)
      else if (s$method == "ccdh") sprintf("shared information: adjusted p < %g (%s); reconstruction index <= %g", s$alpha, s$adjust, s$prob_threshold)
      else if (s$method == "hirata") sprintf("coverage index >= threshold and surrogate p < %g; verdict per pair", s$alpha)
      else if (s$method == "cd") sprintf("bootstrap z-tests of the dimension relations at level %g; verdict per pair", s$alpha)
      else if (s$method == "pcmci") sprintf("MCI test '%s'%s, adjusted p < %g (%s)", s$test, if (isTRUE(s$plus)) " (PCMCI+)" else "", s$alpha, s$adjust)
      else if (s$method == "df") sprintf("Table-2 reading of the orders where determinate, else directed variance-reduction tests with adjusted p < %g (%s)", s$alpha, s$adjust)
      else if (s$method != "bayes_granger") sprintf("adjusted p < %g (%s)", s$alpha, s$adjust)
      else sprintf("BF10 (%s) > %g", switch(if (is.null(s$prior)) "bic" else s$prior, bic = "BIC approximation",
                                            g = sprintf("Zellner g-prior, g = %s", g_txt()),
                                            hyper_g = sprintf("hyper-g prior, a = %g", if (is.null(bc)) 3 else bc$hyper_a)), s$bf_threshold),
      "\n\n", sep = "")
  P <- x$pairs
  # Bayes factors beyond the range of the replicates are shown as "> bound", as in the table of the Bayesian evidence
  bfs <- .fmt_bf(P$bf10)
  if (!is.null(x$bayes) && !is.null(x$bayes$bf_bound)) {
    bd <- x$bayes$bf_bound[match(paste(P$cause, "->", P$effect), x$bayes$pair)]
    bfs <- ifelse(is.finite(bd), paste0(">", .fmt_bf(bd)), bfs)
  }
  tab <- data.frame(
    link = paste(P$cause, "->", P$effect),
    lag = P$lag,
    F = formatC(P$statistic, digits = digits, format = "f"),
    df = ifelse(is.na(P$df1), "-", paste0(P$df1, ",", P$df2)),
    p = .fmt_p(P$p_value),
    p.adj = .fmt_p(P$p_adj),
    BF10 = bfs,
    P.post = formatC(P$post_prob, digits = digits, format = "f"),
    sig = ifelse(P$significant, paste0("yes ", .stars(P$p_adj)), "no"),
    check.names = FALSE, stringsAsFactors = FALSE)
  cmi_like <- info$family != "regression" && !s$method %in% c("dc", "dfe")
  names(tab)[names(tab) == "F"] <- info$stat
  if (cmi_like) tab <- tab[, !names(tab) %in% c("df", if (!isTRUE(s$bayes)) c("BF10", "P.post"))]
  if (s$method %in% c("dc", "dfe")) {
    tab <- tab[, !names(tab) %in% c("df", "p", "p.adj", if (s$method == "dc") "lag")]
    tab$sig <- ifelse(P$significant, "yes", "no")
  }
  print(tab, row.names = FALSE, right = FALSE)
  if (s$method %in% c("dc", "dfe")) {
    .out("  P.direct = score of the unidirectional link, P.post = P(direct) + P(bidirectional), BF10 against the uniform prior\n")
    .out(if (s$method == "dc") "\n  Intrinsic dimensions and posterior probabilities of the five causal cases:\n"
        else "\n  Random-forest class scores of the five causal cases:\n")
    Dc <- x$dc
    tb <- data.frame(pair = paste(Dc$var1, "~", Dc$var2))
    if (s$method == "dc") tb <- cbind(tb, D.X = round(Dc$D_X, 2), D.Y = round(Dc$D_Y, 2), D.J = round(Dc$D_J, 2), D.Z = round(Dc$D_Z, 2))
    tb <- cbind(tb, `x->y` = round(Dc$x_to_y, 3), `y->x` = round(Dc$y_to_x, 3), `<->` = round(Dc$bidirectional, 3),
                common = round(Dc$common_cause, 3), indep = round(Dc$independent, 3), best = Dc$best)
    print(tb, row.names = FALSE, right = FALSE)
  }
  if (s$method == "df") .out("  sd.red = relative reduction of the conditional sd of the effect when the cause's past is constrained\n")
  if (s$method == "ccm") {
    tab$lag <- NULL
    .out("  rho = cross-map skill of 'effect xmap cause' at the largest library; p = surrogate test;\n")
    .out("  a link also requires convergence (Kendall test of rho against the library size):\n")
    Cv <- x$ccm
    print(data.frame(link = paste(Cv$cause, "->", Cv$effect), converging = ifelse(Cv$converging, "yes", "no"),
                     p.convergence = .fmt_p(Cv$p_convergence)), row.names = FALSE, right = FALSE)
  }
  if (s$method == "pid") {
    .out("  U+S.bits = unique + synergistic information (bits per target bit) of the cause's past about the effect's future beyond the effect's own past\n")
    .out("\n  MAR-PID atoms at the selected lag (bits per target bit):\n")
    Pd <- x$pid
    print(data.frame(link = paste(Pd$cause, "->", Pd$effect), lag = Pd$lag, unique.cause = round(Pd$unique_cause, 3),
                     unique.effect = round(Pd$unique_effect, 3), redundant = round(Pd$redundant, 3),
                     synergy = round(Pd$synergy, 3), total = round(Pd$total, 3)), row.names = FALSE, right = FALSE)
  }
  if (s$method == "te") .out("  TE.bits = transfer entropy at the selected lag (bits); effect_size column of $pairs = effective TE\n")
  if (s$method %in% c("cmi", "cmi_lag"))
    .out(sprintf("  CMI = kNN conditional mutual information at the selected lag (nats); p = against %s\n",
                 if (identical(s$method_args$reference, "local")) "local permutations of the cause (Runge 2018)"
                 else "conditional shifts of the cause (location and scale given the conditioning past kept, residuals shifted circularly)"))
  if (s$method == "toda_yamamoto")
    .out(sprintf("  F = Wald test of the first 'lag' cause lags in a VAR augmented by d_max = %d extra lags (not tested)\n",
                as.integer(s$method_args$d_max)))
  if (s$method == "diks_panchenko") {
    E <- x$extra
    .out(sprintf("  T.n = Diks-Panchenko statistic (bandwidth %.3g sd%s); p = %s\n",
                E$eps[1], if (identical(E$filter[1], "var")) ", VAR residuals" else "",
                if (any(E$reference == "conditional shifts", na.rm = TRUE))
                  "against conditional shifts of the cause (location and scale given the conditioning past kept, residuals shifted circularly)"
                else "one-sided asymptotic normal test of z = sqrt(n) T.n / S.n"))
  }
  if (s$method == "rf_granger") {
    E <- x$extra
    .out("  OOB.gain = 1 - MSE(own + cause) / MSE(own) out of bag, p = against conditional shifts of the cause (location and scale given the conditioning past kept, residuals shifted circularly)\n")
    .out("\n  Out-of-sample gain (rolling origin), permutation importance and shadow-feature rule:\n")
    print(data.frame(link = paste(E$cause, "->", E$effect), OOB.gain = sprintf("%+.3f", E$oob_gain),
                     OOS.gain = ifelse(is.na(E$oos_gain), "-", sprintf("%+.3f", E$oos_gain)),
                     VI.best = ifelse(is.na(E$vi_cause), "-", sprintf("%.3g (lag %d)", E$vi_cause, as.integer(E$vi_lag))),
                     shadow = ifelse(is.na(E$shadow_hit), "-", ifelse(E$shadow_hit, "beats shadow", "not above shadow"))),
          row.names = FALSE, right = FALSE)
  }
  if (s$method %in% c("hirata", "cd")) {
    if (s$method == "hirata") .out("  cover.idx = normalised coverage P(R_cause | R_effect); p = surrogate test of the coverage\n")
    else .out("  D2.diff = D2(effect) - D2(cause); no p-values, the verdict comes from the z-tests of the dimension relations\n")
    .out("\n  Per-pair quantities and verdict:\n")
    Cc <- x$ccdh; Cc$pair <- paste(Cc$var1, "~", Cc$var2); Cc$var1 <- Cc$var2 <- NULL
    num <- vapply(Cc, is.numeric, TRUE); Cc[num] <- lapply(Cc[num], round, 3)
    print(Cc[, c("pair", setdiff(names(Cc), "pair"))], row.names = FALSE, right = FALSE)
  }
  if (s$method == "ccdh") {
    .out("  recon = 1 - reconstruction index (1: the cause is perfectly reconstructible from the effect's embedding); p = shared-information surrogate test\n")
    .out("\n  Mean compression rates of the series and of their skewed joints:\n")
    Cc <- x$ccdh
    print(data.frame(pair = paste(Cc$var1, "~", Cc$var2), cr.x = round(Cc$cr_x, 3), cr.y = round(Cc$cr_y, 3),
                     cr.Jxxy = round(Cc$cr_jxxy, 3), cr.Jyyx = round(Cc$cr_jyyx, 3),
                     index.x = round(Cc$index_x, 3), index.y = round(Cc$index_y, 3), verdict = Cc$verdict,
                     check.names = FALSE), row.names = FALSE, right = FALSE)
  }
  if (!is.null(x$instant) && s$method == "pcmci") {
    .out("\n  Contemporaneous links (lag 0):\n")
    I <- x$instant
    print(data.frame(pair = paste(I$var1, "--", I$var2),
                     MCI.stat = formatC(I$statistic, digits = digits, format = "f"),
                     p = .fmt_p(I$p_value), p.adj = .fmt_p(I$p_adj),
                     significant = ifelse(I$significant, "yes", "no"),
                     check.names = FALSE), row.names = FALSE, right = FALSE)
  }
  if (!is.null(x$instant) && !s$method %in% c("df", "dc", "ccdh", "dfe", "hirata", "cd", "pcmci")) {
    .out("\n  Contemporaneous dependence given both pasts (k = 0; a common driver or an instantaneous relation):\n")
    I <- x$instant
    print(data.frame(pair = paste(I$var1, "~", I$var2),
                     CMI = formatC(I$cmi, digits = digits, format = "f"),
                     p = .fmt_p(I$p_value), p.adj = .fmt_p(I$p_adj),
                     dependent = ifelse(I$significant, "yes", "no"),
                     check.names = FALSE), row.names = FALSE, right = FALSE)
  }
  if (!is.null(x$df)) {
    .out("\n  Orders (degrees of freedom = conditions on the observed pasts) and the Table-2 reading:\n")
    Dd <- x$df
    print(data.frame(pair = paste(Dd$var1, "~", Dd$var2), df.x = Dd$df_x, df.y = Dd$df_y,
                     df.joint = Dd$df_joint, conditions.joint = Dd$conditions_joint,
                     `Table 2` = Dd$table2, decided.by = Dd$decided_by, check.names = FALSE), row.names = FALSE, right = FALSE)
    .out("  a determinate Table-2 reading decides the pair; undetermined pairs are decided by the directed tests above\n")
  }
  if (!is.null(x$markov)) {
    .out("\n  Markov property (p small = not Markov) and decision-table case:\n")
    M <- x$markov
    print(data.frame(pair = paste(M$var1, "~", M$var2), m = M$memory,
                     p.x = .fmt_p(M$p_markov_x), p.y = .fmt_p(M$p_markov_y),
                     p.joint = .fmt_p(M$p_markov_joint),
                     case = M$case, `x->y` = M$inj_x_to_y, `y->x` = M$inj_y_to_x,
                     hidden = M$hidden, consistent = M$consistent, check.names = FALSE),
          row.names = FALSE, right = FALSE)
    .out("  injection: i.i. = i.i.d. injection, n.i. = non-i.i.d., d = driving; U_* = hidden own-source\n")
  }
  if (!is.null(x$bayes)) {
    B <- x$bayes
    .out("\n  Bayesian evidence per directed test (effect = posterior mean [credible interval] under H1):\n")
    ext <- if (is.null(B$bf_bound)) rep(FALSE, nrow(B)) else is.finite(B$bf_bound)   # statistics outside the range of their replicates
    tb <- data.frame(link = B$pair, BF10 = ifelse(ext, paste0(">", .fmt_bf(B$bf_bound)), .fmt_bf(B$bf10)),
                     P.post = formatC(B$post_prob, digits = 3, format = "f"),
                     evidence = ifelse(ext, paste("at least", as.character(bf_evidence(B$bf_bound))), as.character(bf_evidence(B$bf10))),
                     effect = ifelse(is.na(B$effect_mean), "-", sprintf("%.3g [%.3g, %.3g]", B$effect_mean, B$effect_lower, B$effect_upper)),
                     model = B$model, check.names = FALSE, stringsAsFactors = FALSE)
    if (!is.null(B$rhat) && any(!is.na(B$rhat)))
      tb$`R-hat` <- ifelse(is.na(B$rhat), "-", formatC(B$rhat, digits = 3, format = "f"))
    print(tb, row.names = FALSE, right = FALSE)
    if (!is.null(B$rhat) && any(!is.na(B$rhat)))
      .out("  R-hat = Gelman-Rubin factor of the effect over the chains (close to 1: the chains agree)\n")
    if (any(ext))
      .out(paste0("  '>': the statistic lies outside the range of its replicates; the value is the Bayes factor at the extreme replicate,\n",
                  "   a lower bound (P.post and the decision use a normal continuation of the tail; see ?surrogate_bf)\n"))
  }
  .out("\n  Classification of series pairs:\n")
  C <- x$classes
  ct <- data.frame(pair = paste(C$var1, "~", C$var2), class = C$class,
                   relation = C$relation,
                   stringsAsFactors = FALSE)
  print(ct, row.names = FALSE, right = FALSE)
  if (any(!C$both_directions_tested))
    .out("  (pairs with only one tested direction are labelled by that direction alone)\n")
  if (!isFALSE(getOption("Rcausality.guide")))
    .out("  How to read: one row per directed test cause -> effect; 'sig' applies the decision rule above; ",
         if (!isTRUE(s$bayes) && s$method %in% c("granger", "toda_yamamoto"))
           "BF10 and P.post (Bayes factor of the unit-information prior, posterior probability for a prior of 1/2) are given for information; " else "",
         "the classification combines the two directions of a pair (x -> y, y -> x, x <-> y, independent; ",
         "'independent' = no relation detected; 'common_cause', printed as <~> = the criterion of the method for shared information ",
         "without a directed link holds: for \"cmi_lag\" and \"markov\" dependence within the sampling interval, which a hidden driver ",
         "with memory does not produce; see the help page of the method for the systems on which the class is a hidden common cause). ",
         "plot(x, \"all\") shows the same. ",
         "(options(Rcausality.guide = FALSE) hides this note.)\n", sep = "")
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
  .out("\nRcausality summary\n", strrep("=", 50), "\n", sep = "")
  .out("Call:   ", deparse(x$call), "\n", sep = "")
  .out("Method: ", s$method, if (s$conditional) " (conditional)" else " (pairwise)",
      "   n = ", s$n, "   alpha = ", s$alpha, " (", s$adjust, ")\n", sep = "")
  .out("Tests:  ", x$n_tests, "   significant links: ", x$n_sig,
      sprintf("   (density %.0f%%)\n", 100 * x$density), sep = "")
  .out("\nPair classes (x = first, y = second variable of the pair):\n")
  print(x$classes)
  if (x$n_sig) {
    .out("\nSignificant links, strongest first:\n")
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
    .out("\nNode roles:\n")
    print(x$degree, row.names = FALSE, right = FALSE)
    .out("\nLag distribution of significant links:\n")
    print(x$lags)
  } else .out("\nNo significant causal links at the chosen threshold.\n")
  invisible(x)
}
