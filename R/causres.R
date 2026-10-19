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
#' * `series`, `causes`, `effects`, `settings`, `call`, `data`.
#'
#' Methods: `print`, `summary`, `plot`, `as.data.frame`, `as.matrix`,
#' `coef`-style accessors via [as.matrix.causres()].
#'
#' @name causres-class
#' @aliases causres
NULL

.new_causres <- function(pairs, series, causes, effects, method, conditional,
                         alpha, bf_threshold, adjust, n, call, data, instant = NULL, markov = NULL) {
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
      ii <- instant[instant$var1 == a & instant$var2 == b, , drop = FALSE]
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
                                 adjust = adjust, n = n),
                 call = call, data = data, instant = instant, markov = markov),
            class = "causres")
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
             bayes_granger = "Bayesian Granger causality (BIC Bayes factor)",
             cmi = "Conditional mutual information (kNN, local permutation)",
             cmi_lag = "Lag-shifted conditional mutual information (Zlatniczki et al., 2021)",
             markov = "Markov-chain approach (Stippinger et al., 2023)")[[s$method]]
  cat("\n  Rcausality: ", mname, "\n", sep = "")
  cat("  ", length(x$series), " series, n = ", s$n, ", ",
      nrow(x$pairs), " directed tests",
      if (s$conditional) " (conditional on the other series)" else " (pairwise)",
      "\n", sep = "")
  cat("  decision: ",
      if (s$method != "bayes_granger") sprintf("adjusted p < %g (%s)", s$alpha, s$adjust)
      else sprintf("BF10 > %g", s$bf_threshold), "\n\n", sep = "")
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
  cmi_like <- s$method %in% c("cmi", "cmi_lag", "markov")
  names(tab)[names(tab) == "F"] <- if (cmi_like) "CMI" else "F"
  if (cmi_like) tab <- tab[, !names(tab) %in% c("df", "BF10", "P.post")]
  print(tab, row.names = FALSE, right = FALSE)
  if (!is.null(x$instant)) {
    cat("\n  Contemporaneous dependence given both pasts (k = 0; common-driver evidence):\n")
    I <- x$instant
    print(data.frame(pair = paste(I$var1, "~", I$var2),
                     CMI = formatC(I$cmi, digits = digits, format = "f"),
                     p = .fmt_p(I$p_value), p.adj = .fmt_p(I$p_adj),
                     common.driver = ifelse(I$significant, "yes", "no"),
                     check.names = FALSE), row.names = FALSE, right = FALSE)
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
  cat("\n  Classification of series pairs:\n")
  C <- x$classes
  ct <- data.frame(pair = paste(C$var1, "~", C$var2), class = C$class,
                   relation = C$relation,
                   stringsAsFactors = FALSE)
  print(ct, row.names = FALSE, right = FALSE)
  if (any(!C$both_directions_tested))
    cat("  (pairs with only one tested direction are labelled by that direction alone)\n")
  invisible(x)
}

#' Summary of a `causres` object
#'
#' @param object a `causres` object.
#' @param ... unused.
#' @return An object of class `summary.causres` (print it).
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
