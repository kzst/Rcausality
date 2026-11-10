# Methods of the `causpanel` class: print, summary, matrix view and plots ----------

.link_key <- function(cause, effect) paste(cause, "->", effect)

# one row per directed link with the three levels side by side
.panel_table <- function(x) {
  P <- x$pairs; PT <- x$panel; s <- x$settings; key <- .link_key(P$cause, P$effect)
  pick <- function(level, col) { k <- PT[PT$level == level, , drop = FALSE]; k[[col]][match(key, .link_key(k$cause, k$effect))] }
  U <- x$units
  cnt <- function(col) if (is.null(U)) rep(NA_integer_, length(key)) else
    vapply(key, function(k) sum(U[[col]][.link_key(U$cause, U$effect) == k] %in% TRUE), 0L)
  nun <- if (is.null(U)) rep(NA_integer_, length(key)) else vapply(key, function(k) sum(.link_key(U$cause, U$effect) == k), 0L)
  pv <- x$prevalence; pk <- if (is.null(pv)) NULL else match(key, .link_key(pv$cause, pv$effect))
  out <- data.frame(link = key, lag = P$lag, units = nun, units_nominal = cnt("nominal"), units_fdr = cnt("significant"),
                    group_stat = pick("group", "statistic"), group_p = pick("group", "p_freq"), group_bf = pick("group", "bf10"),
                    pooled_stat = pick("pooled", "statistic"), pooled_p = pick("pooled", "p_freq"), pooled_bf = pick("pooled", "bf10"),
                    prevalence = if (is.null(pk)) NA_real_ else pv$prevalence[pk],
                    prevalence_lower = if (is.null(pk)) NA_real_ else pv$prevalence_lower[pk],
                    prevalence_upper = if (is.null(pk)) NA_real_ else pv$prevalence_upper[pk],
                    group = pick("group", "significant"), pooled = pick("pooled", "significant"),
                    significant = P$significant, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}

.panel_header <- function(x) {
  s <- x$settings; ps <- s$panel
  bal <- if (ps$n_range[1] == ps$n_range[2]) sprintf("%d periods", ps$n_range[1]) else sprintf("%d-%d periods", ps$n_range[1], ps$n_range[2])
  cat(sprintf("Panel causality analysis: %s%s\n", .method_info(s$base_method)$short, if (isTRUE(s$bayes)) ", Bayesian variant" else ""))
  cat(sprintf("  %d units x %s, %d series%s%s\n", ps$n_units, bal, length(x$series),
              if (isTRUE(s$conditional)) ", conditional on the other series" else "",
              if (identical(ps$lag_rule, "pooled BIC")) "; lag order of every link by the BIC of its pooled regression"
              else if (identical(ps$lag_rule, "first")) "; first significant lag order" else ""))
  lv <- c(group = "group level (the unit results combined)", pooled = "pooled level (stacked units)",
          both = if (ps$use_bf) "group or pooled level (average of the two Bayes factors)" else "group or pooled level (smaller p-value, doubled)")
  cat(sprintf("  decision: %s; link if %s\n", lv[[ps$decision]],
              if (ps$use_bf) sprintf("BF10 > %g", s$bf_threshold) else sprintf("adj. p < %g (%s)", s$alpha, s$adjust)))
  if (ps$units_analysed > 0)
    cat(sprintf("  unit level: %d of %d units analysed%s; error control over the unit x link cells: %s\n", ps$units_analysed, ps$n_units,
                if (ps$units_short > 0) sprintf(" (%d shorter than %d periods)", ps$units_short, ps$min_length) else "",
                if (ps$use_bf) "Bayesian FDR" else sprintf("%s", ps$unit_adjust)))
  else cat(sprintf("  unit level: not available (units shorter than %d periods, see panel_control(min_length = ))\n", ps$min_length))
  if (ps$units_failed > 0) cat(sprintf("  note: the unit-level analysis failed in %d unit(s): %s\n", ps$units_failed, x$unit_errors[[1]]))
}

#' @export
print.causpanel <- function(x, digits = 3, ...) {
  .panel_header(x)
  tb <- .panel_table(x); ps <- x$settings$panel
  show <- data.frame(link = tb$link, lag = ifelse(is.na(tb$lag), "-", as.character(tb$lag)), stringsAsFactors = FALSE)
  if (ps$units_analysed > 0) show$units <- sprintf("%d/%d (%d)", tb$units_nominal, tb$units, tb$units_fdr)
  ev <- function(p, bf) if (ps$use_bf) ifelse(is.na(bf), "-", paste("BF", .fmt_bf(bf))) else ifelse(is.na(p), "-", .fmt_p(p))
  if ("group" %in% ps$levels) show$group <- paste0(ev(tb$group_p, tb$group_bf), ifelse(tb$group %in% TRUE, " *", ""))
  if ("pooled" %in% ps$levels) show$pooled <- paste0(ev(tb$pooled_p, tb$pooled_bf), ifelse(tb$pooled %in% TRUE, " *", ""))
  show$panel <- paste0(if (ps$use_bf) ifelse(is.na(x$pairs$bf10), "-", paste("BF", .fmt_bf(x$pairs$bf10))) else
    ifelse(is.na(x$pairs$p_adj), "-", .fmt_p(x$pairs$p_adj)), ifelse(tb$significant, " *", ""))
  if (any(is.finite(tb$prevalence)))
    show$prevalence <- ifelse(is.finite(tb$prevalence), sprintf("%.2f [%.2f, %.2f]", tb$prevalence, tb$prevalence_lower, tb$prevalence_upper), "-")
  cat("\n"); print(show, row.names = FALSE, right = FALSE)
  cat("\n")
  if (ps$units_analysed > 0) cat("  units: units with the link / units analysed (after the error control over all cells)\n")
  cat(sprintf("  %s\n", if (ps$use_bf) "group: Bayes factor of 'some units have the link'; pooled: Bayes factor of the stacked units"
                         else "group, pooled: p-values of the two panel-level tests"))
  cat(sprintf("  panel: the decision -- %s; * = link\n",
              if (ps$decision == "both") { if (ps$use_bf) "average of the two Bayes factors" else sprintf("the smaller p-value doubled, %s-adjusted over the links", x$settings$adjust) }
              else if (ps$use_bf) sprintf("the Bayes factor of the %s level", ps$decision) else sprintf("the %s-adjusted p-value of the %s level", x$settings$adjust, ps$decision)))
  if (any(is.finite(tb$prevalence))) cat("  prevalence: posterior share of the units with the link [credible interval]\n")
  C <- x$classes
  cat("\nRelations in the panel:\n")
  for (i in seq_len(nrow(C))) cat(sprintf("  %s\n", C$relation[i]))
  invisible(x)
}

#' Summary of a panel causality analysis
#'
#' @param object a `causpanel` object.
#' @param ... unused.
#' @return A list of class `summary.causpanel` with the link table (`links`:
#'   one row per directed link with the unit counts, the group- and
#'   pooled-level statistics, p-values and Bayes factors, the posterior
#'   prevalence and the decisions), the unit table restricted to the units
#'   with a link after the error control (`units`), and the settings.
#' @export
#' @examples
#' pd <- sim_panel(n_units = 15, n_time = 20, seed = 1)
#' summary(panel_causality(pd, id = "unit", time = "time"))
summary.causpanel <- function(object, ...) {
  U <- object$units
  structure(list(links = .panel_table(object), units = if (is.null(U)) NULL else U[U$significant %in% TRUE, , drop = FALSE],
                 panel = object$panel, settings = object$settings, object = object), class = "summary.causpanel")
}

#' @export
print.summary.causpanel <- function(x, digits = 3, ...) {
  .panel_header(x$object)
  cat("\nLinks (three levels):\n")
  tb <- x$links; num <- vapply(tb, is.numeric, TRUE)
  tb[num] <- lapply(tb[num], function(v) signif(v, digits))
  print(tb, row.names = FALSE)
  cat("\nPanel-level tests:\n")
  PT <- x$panel[, c("cause", "effect", "level", "test", "lag", "statistic", "p_freq", "bf10", "significant", "detail")]
  PT$statistic <- signif(PT$statistic, digits); PT$p_freq <- signif(PT$p_freq, digits); PT$bf10 <- signif(PT$bf10, digits)
  print(PT, row.names = FALSE, right = FALSE)
  if (!is.null(x$units)) {
    cat(sprintf("\nUnits with a link after the error control: %d cell(s)\n", nrow(x$units)))
    if (nrow(x$units)) {
      tu <- table(.link_key(x$units$cause, x$units$effect))
      for (k in names(tu)) cat(sprintf("  %s: %s\n", k, paste(utils::head(x$units$unit[.link_key(x$units$cause, x$units$effect) == k], 12), collapse = ", ")))
    }
  }
  invisible(x)
}

#' Matrix views of a panel causality analysis
#'
#' @param x a `causpanel` object.
#' @param what a column of the pairs table (`by = "pairs"`, see
#'   [as.matrix.causres()]) or of the unit table (`by = "units"`: `"lag"`,
#'   `"statistic"`, `"p_value"`, `"q_value"`, `"bf10"`, `"post_pooled"`,
#'   `"nominal"`, `"significant"`; `"first_lag"` gives the lag of the links
#'   that are significant in the unit and 0 otherwise).
#' @param by `"pairs"`: causes x effects matrix of the panel-level result;
#'   `"units"`: links x units matrix of the unit-level results.
#' @param ... unused.
#' @return A matrix.
#' @export
#' @examples
#' pd <- sim_panel(n_units = 8, n_time = 25, seed = 1)
#' pg <- panel_causality(pd, id = "unit", time = "time")
#' as.matrix(pg, "p_value", by = "units")
as.matrix.causpanel <- function(x, what = if (by == "pairs") "p_adj" else "first_lag", by = c("pairs", "units"), ...) {
  by <- match.arg(by)
  if (by == "pairs") return(NextMethod("as.matrix", x, what = what))
  U <- x$units; if (is.null(U)) stop("No unit-level results in this object.", call. = FALSE)
  if (what == "first_lag") { U$first_lag <- ifelse(U$nominal %in% TRUE, ifelse(is.na(U$lag), 1L, U$lag), 0L) }
  if (!what %in% names(U)) stop("Unknown column: ", what, call. = FALSE)
  key <- .link_key(x$pairs$cause, x$pairs$effect); ids <- unique(U$unit)
  M <- matrix(NA, length(key), length(ids), dimnames = list(key, ids))
  M[cbind(match(.link_key(U$cause, U$effect), key), match(U$unit, ids))] <- U[[what]]
  M
}

# ---- plots ---------------------------------------------------------------------

.lag_pal <- function(n) grDevices::colorRampPalette(c("#FFF2A8", "#FDC66B", "#F58B3C", "#D4502A", "#8F2D23", "#5A1A1A"))(max(n, 5L))[seq_len(max(n, 1L))]
# black or white text on a fill colour
.txt_on <- function(fill) { v <- grDevices::col2rgb(fill) / 255; if (sum(c(0.299, 0.587, 0.114) * v[, 1L]) < 0.5) "white" else "grey10" }
.ev_pal <- function(n) grDevices::colorRampPalette(c("#E3EEF9", "#9EC5E8", "#4F94CD", "#1F5FA8", "#0B3C78"))(n)

#' Plot a panel causality analysis
#'
#' @param x a `causpanel` object.
#' @param type one of
#'   * `"heatmap"`: links (rows) by units and the panel-level tests
#'     (columns). A coloured cell is a link at the 5% level (or with a Bayes
#'     factor above the threshold) in that unit; the colour is the lag of the
#'     link or the strength of the evidence (`fill`); a star marks the cells
#'     that remain after the error control over all cells (unit columns) and
#'     the links of the panel (columns `Group`, `Pooled`). The units are
#'     ordered by a hierarchical clustering of their link patterns;
#'   * `"units"`: the evidence of every unit for a link, sorted, with the
#'     decision line -- shows whether a panel-level link is carried by many
#'     units or by a few;
#'   * `"prevalence"`: share of the units with each link (observed, after the
#'     error control, and the posterior of the prevalence model with its
#'     credible interval);
#'   * `"network"`, `"matrix"`, `"bars"`, `"lags"`: the plots of
#'     [plot.causres()] for the panel-level result.
#' @param fill colour of the heat map cells: `"lag"`, `"evidence"`
#'   (`-log10 p` or `log10` Bayes factor) or `"auto"` (the lag when the lags
#'   of the links differ, otherwise the evidence).
#' @param links optional character vector of the links to show (as
#'   `"x -> y"`), or `"significant"` for the links of the panel.
#' @param cluster logical; order the units by clustering (`"heatmap"`)?
#' @param evidence,main as in [plot.causres()].
#' @param ... passed on to [plot.causres()].
#' @return Invisibly, `x`.
#' @references
#' Kosztyan, Z. T., Matrai, R. and Kiss, D. (2026). Is the world getting smaller? Network-driven heterogeneity in gravity-model elasticities. *The World Economy*. \doi{10.1111/twec.70148}
#' @export
#' @examples
#' pd <- sim_panel(n_units = 20, n_time = 25, prevalence = 0.5, coupling = 0.6, seed = 1)
#' pg <- panel_causality(pd, id = "unit", time = "time")
#' plot(pg)                    # in which units does x drive y?
#' \donttest{
#' plot(pg, "units")
#' plot(pg, "prevalence")
#' plot(pg, "network")
#' }
plot.causpanel <- function(x, type = c("heatmap", "units", "prevalence", "network", "matrix", "bars", "lags"),
                           fill = c("auto", "lag", "evidence"), links = NULL, cluster = TRUE, evidence = c("p", "bf"),
                           main = NULL, ...) {
  type <- match.arg(type); fill <- match.arg(fill); evidence <- match.arg(evidence)
  if (type %in% c("network", "matrix", "bars", "lags")) {
    if (type == "lags" && is.null(x$by_lag)) stop("No lag scan in this object (pooled level of a method that scans the lags).", call. = FALSE)
    y <- x; class(y) <- "causres"
    if (isTRUE(x$settings$panel$use_bf)) evidence <- "bf"
    plot(y, type = type, evidence = evidence, main = if (is.null(main)) sprintf("Panel-level causal %s", if (type == "network") "network" else type) else main, ...)
    return(invisible(x))
  }
  key <- .link_key(x$pairs$cause, x$pairs$effect)
  if (!is.null(links)) {
    sel <- if (identical(links, "significant")) key[x$pairs$significant] else intersect(key, links)
    if (!length(sel)) stop("No link to show.", call. = FALSE)
  } else sel <- key
  switch(type, heatmap = .plot_panel_heat(x, sel, fill, cluster, main), units = .plot_panel_units(x, sel, main),
         prevalence = .plot_panel_prev(x, sel, main))
  invisible(x)
}

# evidence of unit cells / panel rows on a common scale
.cell_evidence <- function(p, bf, use_bf) if (use_bf) log10(pmax(bf, 1e-300)) else -log10(pmax(p, 1e-300))

.plot_panel_heat <- function(x, sel, fill, cluster, main) {
  s <- x$settings; ps <- s$panel; U <- x$units; PT <- x$panel
  ids <- if (is.null(U)) character(0) else unique(U$unit)
  nl <- length(sel); nu <- length(ids)
  lagM <- sigM <- nomM <- evM <- matrix(NA_real_, nl, nu, dimnames = list(sel, ids))
  if (nu) {
    r <- match(.link_key(U$cause, U$effect), sel); cc <- match(U$unit, ids); ok <- !is.na(r)
    idx <- cbind(r[ok], cc[ok])
    lagM[idx] <- ifelse(is.na(U$lag[ok]), 1, U$lag[ok]); nomM[idx] <- U$nominal[ok]; sigM[idx] <- U$significant[ok]
    evM[idx] <- .cell_evidence(U$p_value[ok], U$bf10[ok], ps$use_bf)
  }
  lv <- intersect(c("group", "pooled"), ps$levels)
  pan <- lapply(lv, function(l) { k <- PT[PT$level == l, ]; k[match(sel, .link_key(k$cause, k$effect)), ] })
  names(pan) <- lv
  # order of the units: clustering of the link patterns, then the number of links
  ord <- seq_len(nu)
  if (nu > 2L) {
    B <- t(ifelse(is.na(nomM), 0, nomM) + ifelse(is.na(sigM), 0, sigM)); nlk <- rowSums(B > 0)
    ord <- order(-nlk, apply(B > 0, 1L, function(r) paste(1L - r, collapse = "")))
    if (cluster && nl > 2L && sum(nlk > 0) > 2L && nrow(unique(B[nlk > 0, , drop = FALSE])) > 1L) {
      w <- which(nlk > 0); ord <- c(w[stats::hclust(stats::dist(B[w, , drop = FALSE]), "average")$order], which(nlk == 0))
    }
  }
  all_lag <- c(lagM[nomM %in% 1], unlist(lapply(pan, function(k) k$lag[k$significant %in% TRUE])))
  all_lag <- all_lag[is.finite(all_lag)]
  if (fill == "auto") fill <- if (length(unique(all_lag)) > 1L) "lag" else "evidence"
  maxlag <- if (length(all_lag)) max(all_lag) else 1
  ev_all <- c(evM[nomM %in% 1], unlist(lapply(pan, function(k) .cell_evidence(k$p_freq, k$bf10, ps$use_bf)[k$significant %in% TRUE])))
  ev_cap <- max(4, min(12, stats::quantile(ev_all[is.finite(ev_all)], 0.95, na.rm = TRUE, names = FALSE)), na.rm = TRUE)
  ev0 <- if (ps$use_bf) log10(s$bf_threshold) else -log10(s$alpha)
  colf <- function(lag, ev) {
    if (fill == "lag") .lag_pal(maxlag)[pmax(1, pmin(maxlag, ifelse(is.na(lag), 1, lag)))]
    else .ev_pal(100)[pmax(1, pmin(100, ceiling(100 * (pmin(ev, ev_cap) - ev0) / max(ev_cap - ev0, 1e-9))))]
  }
  np <- length(lv)
  plab <- lapply(pan, function(k) paste0(if (ps$use_bf) ifelse(is.na(k$bf10), "", paste0("BF ", .fmt_bf(k$bf10))) else ifelse(is.na(k$p_adj), "", .fmt_p(k$p_adj)),
                                        ifelse(k$significant %in% TRUE, " *", "")))
  lab_w <- max(graphics::strwidth(sel, units = "inches", cex = 0.8)) + 0.25
  mai <- c(if (nu && nu <= 60) 0.9 else 0.45, lab_w, 1.0, 1.6)
  op <- graphics::par(mai = mai, xpd = NA); on.exit(graphics::par(op))
  graphics::plot.new()
  # width of a panel column in units of a unit cell, so that its label fits
  pw_in <- max(0.55, max(graphics::strwidth(unlist(plab), units = "inches", cex = 0.62, font = 2), 0) + 0.12)
  pin <- max(1, graphics::par("fin")[1] - mai[2] - mai[4])
  cell_in <- if (nu) max(0.02, (pin - np * pw_in * 1.06 - 0.12) / nu) else pw_in
  pw <- pw_in / cell_in; pgap <- 0.06 * pw; gap <- if (nu && np) 0.12 / cell_in else 0
  W <- nu + gap + np * (pw + pgap)
  graphics::plot.window(c(0, W), c(nl, 0), xaxs = "i", yaxs = "i")
  cellcex <- if (nu > 40) 0.7 else 0.95
  for (i in seq_len(nl)) {
    for (jj in seq_len(nu)) {
      j <- ord[jj]; an <- !is.na(nomM[i, j])
      col <- if (!an) "grey93" else if (nomM[i, j] == 1) colf(lagM[i, j], evM[i, j]) else "white"
      graphics::rect(jj - 1, i - 1, jj, i, col = col, border = "grey82", lwd = 0.5)
      if (an && sigM[i, j] %in% 1) graphics::text(jj - 0.5, i - 0.42, "*", cex = cellcex * 1.25, font = 2, col = .txt_on(col))
    }
    for (q in seq_len(np)) {
      k <- pan[[q]][i, ]; x0 <- nu + gap + (q - 1) * (pw + pgap)
      sig <- k$significant %in% TRUE
      col <- if (is.na(k$statistic) && is.na(k$bf10)) "grey93" else if (sig) colf(k$lag, .cell_evidence(k$p_freq, k$bf10, ps$use_bf)) else "white"
      graphics::rect(x0, i - 1, x0 + pw, i, col = col, border = "grey40", lwd = 1)
      graphics::text(x0 + pw / 2, i - 0.5, plab[[q]][i], cex = 0.62, col = .txt_on(col), font = if (sig) 2 else 1)
    }
  }
  graphics::text(-0.3, seq_len(nl) - 0.5, sel, adj = 1, cex = 0.8)
  if (nu && nu <= 60) graphics::text(seq_len(nu) - 0.5, nl + 0.15, ids[ord], srt = 90, adj = 1, cex = min(0.7, 28 / nu + 0.25))
  if (nu) graphics::text(nu / 2, -0.25, sprintf("units (%d)", nu), cex = 0.75, col = "grey30", adj = c(0.5, 0))
  for (q in seq_len(np)) graphics::text(nu + gap + (q - 1) * (pw + pgap) + pw / 2, -0.25, c(group = "Group", pooled = "Pooled")[[lv[q]]], cex = 0.75,
                                        font = if (ps$decision %in% c(lv[q], "both")) 2 else 1, adj = c(0.5, 0))
  graphics::title(main = if (is.null(main)) "Causal links by unit and in the panel" else main, line = 3, font.main = 2)
  graphics::mtext(sprintf("%s%s;  %d units", .method_info(s$base_method)$short, if (isTRUE(s$bayes)) " (Bayesian)" else "", ps$n_units),
                  side = 3, line = 1.9, cex = 0.72, col = "grey30")
  # legend
  xl <- W + 0.04 * W; cw <- 0.035 * W; yt <- 0
  dy <- 0.17 * nl / max(graphics::par("pin")[2], 0.1)                # legend line spacing: 0.17 inch, in row units
  if (fill == "lag") {
    graphics::text(xl, yt, "lag of the link", adj = c(0, 0.5), cex = 0.68, font = 2)
    for (l in seq_len(maxlag)) {
      graphics::rect(xl, yt + l * dy - dy * 0.4, xl + cw, yt + l * dy + dy * 0.4, col = .lag_pal(maxlag)[l], border = "grey60")
      graphics::text(xl + cw * 1.3, yt + l * dy, l, adj = c(0, 0.5), cex = 0.65)
    }
    yb <- yt + (maxlag + 1.2) * dy
  } else {
    graphics::text(xl, yt, if (ps$use_bf) "log10 BF10" else "-log10 p", adj = c(0, 0.5), cex = 0.68, font = 2)
    sq <- seq(ev0, ev_cap, length.out = 5)
    for (l in seq_along(sq)) {
      graphics::rect(xl, yt + l * dy - dy * 0.4, xl + cw, yt + l * dy + dy * 0.4, col = .ev_pal(100)[max(1, round(100 * (l - 1) / 4))], border = "grey60")
      graphics::text(xl + cw * 1.3, yt + l * dy, sprintf("%s%.1f", if (l == 5) ">= " else "", sq[l]), adj = c(0, 0.5), cex = 0.65)
    }
    yb <- yt + 6.2 * dy
  }
  graphics::rect(xl, yb - dy * 0.4, xl + cw, yb + dy * 0.4, col = "white", border = "grey60")
  graphics::text(xl + cw * 1.3, yb, "no link", adj = c(0, 0.5), cex = 0.65)
  if (any(is.na(nomM)) && nu) {
    graphics::rect(xl, yb + dy * 0.6, xl + cw, yb + dy * 1.4, col = "grey93", border = "grey60")
    graphics::text(xl + cw * 1.3, yb + dy, "not analysed", adj = c(0, 0.5), cex = 0.65); yb <- yb + dy
  }
  graphics::text(xl, yb + dy * 1.4, sprintf("* units: %s", if (ps$use_bf) "Bayesian FDR" else paste0(ps$unit_adjust, if (ps$unit_adjust %in% c("BH", "BY", "fdr")) "-FDR" else "")),
                 adj = c(0, 0.5), cex = 0.65)
  graphics::text(xl, yb + dy * 2.4, sprintf("* panel: %s", if (ps$use_bf) sprintf("BF10 > %g", s$bf_threshold) else sprintf("adj. p < %g", s$alpha)),
                 adj = c(0, 0.5), cex = 0.65)
}

.plot_panel_units <- function(x, sel, main) {
  s <- x$settings; ps <- s$panel; U <- x$units
  if (is.null(U)) stop("No unit-level results in this object.", call. = FALSE)
  if (length(sel) > 9L) { sel <- sel[order(-x$pairs$significant[match(sel, .link_key(x$pairs$cause, x$pairs$effect))])][1:9] }
  mf <- if (length(sel) <= 3L) c(1L, length(sel)) else grDevices::n2mfrow(length(sel))
  op <- graphics::par(mfrow = mf, mar = c(3.4, 4, 3.2, 1), mgp = c(2.1, 0.7, 0)); on.exit(graphics::par(op))
  thr <- if (ps$use_bf) log10(s$bf_threshold) else -log10(s$alpha)
  for (k in sel) {
    u <- U[.link_key(U$cause, U$effect) == k, , drop = FALSE]
    ev <- .cell_evidence(u$p_value, u$bf10, ps$use_bf)
    if (all(!is.finite(ev))) ev <- u$statistic
    o <- order(ev, decreasing = TRUE); ev <- ev[o]; u <- u[o, , drop = FALSE]
    col <- ifelse(u$significant %in% TRUE, "#1B7837", ifelse(u$nominal %in% TRUE, "#A6DBA0", "grey70"))
    graphics::plot(seq_along(ev), ev, type = "n", xlab = "units (sorted)", ylab = if (ps$use_bf) "log10 BF10" else "-log10 p",
                   ylim = range(c(0, ev, thr), finite = TRUE) * c(1, 1.18), xaxt = if (nrow(u) <= 25) "n" else "s", las = 1, bty = "l")
    if (nrow(u) <= 25) graphics::axis(1, at = seq_along(ev), labels = u$unit, las = 2, cex.axis = 0.6)
    graphics::abline(h = thr, lty = 2, col = "grey40"); graphics::abline(h = 0, col = "grey85")
    graphics::segments(seq_along(ev), 0, seq_along(ev), ev, col = col, lwd = 2)
    graphics::points(seq_along(ev), ev, pch = 21, bg = col, col = "grey25", cex = 1.05)
    pr <- x$pairs[.link_key(x$pairs$cause, x$pairs$effect) == k, , drop = FALSE]
    sub <- if (!nrow(pr)) "" else if (ps$use_bf) sprintf("panel: BF10 = %s", .fmt_bf(pr$bf10)) else sprintf("panel: adj. p = %s", .fmt_p(pr$p_adj))
    graphics::title(main = if (is.null(main)) k else paste(main, k, sep = ": "), line = 1.7, font.main = 2, cex.main = 1)
    graphics::mtext(sprintf("%d of %d units (%d after error control); %s", sum(u$nominal), nrow(u), sum(u$significant), sub),
                    side = 3, line = 0.4, cex = 0.62, col = "grey30")
  }
  graphics::legend("topright", bty = "n", cex = 0.62, pch = 21, pt.bg = c("#1B7837", "#A6DBA0", "grey70"), col = "grey25",
                   legend = c("link after error control", "link in the unit only", "no link"), inset = c(0, -0.02))
}

.plot_panel_prev <- function(x, sel, main) {
  s <- x$settings; ps <- s$panel; tb <- .panel_table(x); tb <- tb[match(sel, tb$link), , drop = FALSE]
  if (all(is.na(tb$units))) stop("No unit-level results in this object.", call. = FALSE)
  n <- nrow(tb); o <- order(tb$units_nominal / tb$units); tb <- tb[o, , drop = FALSE]
  lab_w <- max(graphics::strwidth(tb$link, units = "inches", cex = 0.85)) + 0.3
  op <- graphics::par(mai = c(0.9, lab_w, 0.95, 1.5), xpd = NA); on.exit(graphics::par(op))
  graphics::plot.new(); graphics::plot.window(c(0, 1), c(0.4, n + 0.6), xaxs = "i")
  graphics::axis(1, at = seq(0, 1, 0.2), labels = paste0(seq(0, 100, 20), "%"), cex.axis = 0.8)
  graphics::title(xlab = "share of the units with the link", line = 2.3, cex.lab = 0.85)
  y <- seq_len(n)
  graphics::rect(0, y - 0.3, tb$units_nominal / tb$units, y + 0.3, col = "#A6DBA0", border = NA)
  graphics::rect(0, y - 0.3, tb$units_fdr / tb$units, y + 0.3, col = "#1B7837", border = NA)
  if (!ps$use_bf) graphics::segments(s$alpha, 0.4, s$alpha, n + 0.6, lty = 3, col = "grey40", xpd = FALSE)
  hp <- is.finite(tb$prevalence)
  if (any(hp)) {
    graphics::segments(tb$prevalence_lower[hp], y[hp], tb$prevalence_upper[hp], y[hp], lwd = 2, col = "grey15")
    graphics::points(tb$prevalence[hp], y[hp], pch = 23, bg = "white", col = "grey15", cex = 1.3, lwd = 1.5)
  }
  graphics::text(-0.02, y, tb$link, adj = 1, cex = 0.85, font = ifelse(tb$significant, 2, 1))
  ev <- if (ps$use_bf) ifelse(is.na(tb$group_bf), "", paste0("BF ", vapply(tb$group_bf, .fmt_bf, ""))) else
    ifelse(is.na(tb$group_p), "", paste0("p ", vapply(tb$group_p, .fmt_p, "")))
  graphics::text(1.02, y, paste0(ev, ifelse(tb$group %in% TRUE, " *", "")), adj = 0, cex = 0.72)
  graphics::text(1.02, n + 0.75, "group test", adj = 0, cex = 0.68, col = "grey30")
  graphics::title(main = if (is.null(main)) "Prevalence of the links across the units" else main, line = 2.6, font.main = 2)
  graphics::mtext(sprintf("%s%s;  %d units", .method_info(s$base_method)$short, if (isTRUE(s$bayes)) " (Bayesian)" else "", ps$units_analysed),
                  side = 3, line = 1.5, cex = 0.72, col = "grey30")
  graphics::legend(0, 0.35, horiz = TRUE, bty = "n", cex = 0.68, xjust = 0, yjust = 1,
                   legend = c("units with the link", "after error control", if (any(hp)) "posterior prevalence [CI]"),
                   pch = c(22, 22, if (any(hp)) 23), pt.bg = c("#A6DBA0", "#1B7837", if (any(hp)) "white"), col = c(NA, NA, if (any(hp)) "grey15"),
                   pt.cex = 1.4)
}
