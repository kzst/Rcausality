.pal <- list(causal = "#A8D5BA", nocausal = "#F4A6A6", node = "#BFD7EA",
             border = "#555555", pos = "#2E7D32", neg = "#C62828",
             grid = "grey85", text = "#212121")

#' Plot a `causres` object
#'
#' @param x a `causres` object.
#' @param type one of
#'   * `"network"`: directed graph of the significant links (edge width ~
#'     evidence, label = lag, colour = sign of the net influence, dashed
#'     orange = common-driver evidence);
#'   * `"matrix"` (alias `"heatmap"`): causal matrix, rows = cause, columns =
#'     effect; every cell is a directed test "row -> column", coloured by
#'     evidence, significant cells framed and starred;
#'   * `"bars"`: evidence per directed pair, sorted, with the decision line;
#'   * `"lags"`: test statistic of every directed pair as a function of the
#'     candidate lag (the selected lag is marked);
#'   * `"embedding"`: time-delay (Takens) embedding of every series, see
#'     [plot_embedding()];
#'   * `"posterior"`: Bayesian variants (`bayes = TRUE`): posterior
#'     probability of every directed link with its Bayes factor, and the
#'     credible interval of the effect under the alternative;
#'   * `"both"`: network + matrix side by side; `"all"`: the four panels
#'     network, matrix, lags, embedding.
#' @param evidence evidence shown for widths/colours: `"p"`
#'   (-log10 adjusted p) or `"bf"` (log10 Bayes factor).
#' @param main optional title.
#' @param ... passed on to [plot_embedding()] for `type = "embedding"`.
#' @return Invisibly, `x`.
#' @references
#' Granger, C. W. J. (1969). Investigating causal relations by econometric models and cross-spectral methods. *Econometrica*, 37(3), 424-438. \doi{10.2307/1912791}
#'
#' Takens, F. (1981). Detecting strange attractors in turbulence. In: *Dynamical Systems and Turbulence, Warwick 1980*, Lecture Notes in Mathematics 898, 366-381. Springer. \doi{10.1007/BFb0091924}
#'
#' Zlatniczki, A., Stippinger, M., Benko, Z., Somogyvari, Z. and Telcs, A. (2021). Relaxation of some confusions about confounders. *Entropy*, 23(11), 1450. \doi{10.3390/e23111450}
#' @export
#' @examples
#' # a VAR where x drives y: the network shows one arrow x -> y with its lag
#' # and p-value, the matrix the p-value of every directed test, the lag scan
#' # the F statistic over the candidate lags
#' d <- sim_var(300, "x_to_y", seed = 1)
#' r <- rcausality(d)
#' plot(r)                 # the network
#' \donttest{
#' plot(r, "all")          # network, matrix and lag scan side by side
#' plot(r, "matrix")
#' plot(r, "lags")
#' # the posterior panel of a Bayesian run
#' rb <- rcausality(d, bayes = TRUE, control = bayes_control(n_iter = 300))
#' plot(rb, "posterior")
#' }
plot.causres <- function(x, type = c("network", "matrix", "heatmap", "bars", "lags",
                                     "embedding", "posterior", "both", "all"),
                         evidence = c("p", "bf"), main = NULL, ...) {
  type <- match.arg(type); evidence <- match.arg(evidence)
  if (type == "heatmap") type <- "matrix"
  if (type %in% c("both", "all")) {
    op <- graphics::par(mfrow = if (type == "both") c(1, 2) else c(2, 2)); on.exit(graphics::par(op))
    .plot_network(x, evidence, main); .plot_heat(x, evidence, main)
    if (type == "all") { .plot_lags(x, main); .plot_embed_all(x, main, single = TRUE, ...) }
  } else switch(type, network = .plot_network(x, evidence, main),
                matrix = .plot_heat(x, evidence, main),
                bars = .plot_bars(x, evidence, main),
                lags = .plot_lags(x, main),
                posterior = .plot_posterior(x, main),
                embedding = .plot_embed_all(x, main, ...))
  invisible(x)
}

# posterior probabilities of the directed links (Bayesian variants)
.plot_posterior <- function(x, main) {
  P <- x$pairs; s <- x$settings
  if (all(is.na(P$post_prob))) { graphics::plot.new(); graphics::title("No posterior probabilities: run rcausality(..., bayes = TRUE)"); return(invisible()) }
  o <- order(P$post_prob); P <- P[o, , drop = FALSE]; B <- x$bayes
  lab <- paste(P$cause, "->", P$effect)
  has_eff <- !is.null(B) && any(!is.na(B$effect_mean))
  op <- graphics::par(mfrow = if (has_eff) c(1, 2) else c(1, 1), mar = c(5.2, max(6, 0.6 * max(nchar(lab))), 4, 1.5)); on.exit(graphics::par(op))
  thr <- if (!is.null(s$bayes_control) && identical(s$bayes_control$decision, "posterior")) s$prob_threshold else NA
  bp <- graphics::barplot(P$post_prob, horiz = TRUE, names.arg = lab, las = 1, cex.names = 0.75, xlim = c(0, 1.18),
                          col = ifelse(P$significant %in% TRUE, .pal$causal, .pal$nocausal), border = .pal$border,
                          xlab = "posterior probability of the link", main = "", xpd = NA)
  graphics::title(main = if (is.null(main)) "Posterior probability per directed link" else main, font.main = 2, line = 2.2)
  graphics::mtext(sprintf("%s; %s", .method_label(s), .rule_label(s)), side = 3, line = 0.6, cex = 0.72, col = "grey30")
  graphics::text(pmin(P$post_prob, 1) + 0.02, bp, paste0("BF ", .fmt_bf(P$bf10)), adj = 0, cex = 0.68, col = "grey25", xpd = NA)
  if (is.finite(thr)) graphics::abline(v = thr, lty = 2, col = "grey40") else graphics::abline(v = 0.5, lty = 3, col = "grey75")
  graphics::mtext("green = link declared; label = Bayes factor BF10", side = 1, line = 4, cex = 0.65, col = "grey40")
  if (has_eff) {
    Bo <- B[match(paste(P$cause, "->", P$effect), B$pair), ]
    xr <- range(c(0, Bo$effect_lower, Bo$effect_upper), na.rm = TRUE)
    graphics::plot(NA, xlim = xr, ylim = c(0.5, nrow(Bo) + 0.5), yaxt = "n", ylab = "", xlab = "effect under H1 (posterior mean, interval)",
                   main = "", las = 1)
    graphics::title(main = "Effect size", font.main = 2, line = 2.2)
    graphics::mtext(sprintf("%.0f%% credible intervals", 100 * (if (is.null(s$bayes_control)) 0.95 else s$bayes_control$credible)),
                    side = 3, line = 0.6, cex = 0.72, col = "grey30")
    graphics::abline(v = 0, lty = 3, col = "grey60")
    graphics::segments(Bo$effect_lower, seq_len(nrow(Bo)), Bo$effect_upper, seq_len(nrow(Bo)),
                       col = ifelse(P$significant %in% TRUE, .pal$pos, "grey55"), lwd = 2.5)
    graphics::points(Bo$effect_mean, seq_len(nrow(Bo)), pch = 16, col = ifelse(P$significant %in% TRUE, .pal$pos, "grey40"), cex = 1.1)
    graphics::axis(2, seq_len(nrow(Bo)), lab, las = 1, cex.axis = 0.75, tick = FALSE)
  }
  invisible(x)
}

# evidence per directed pair: -log10 adjusted p, log10 BF10, or -- for the
# posterior / verdict based methods without p-values -- the posterior
# probability or the statistic itself
.evid <- function(P, evidence) {
  e <- if (evidence == "p") -log10(pmax(P$p_adj, 1e-300)) else log10(pmax(P$bf10, 1e-300))
  if (all(is.na(e))) e <- if (any(!is.na(P$post_prob))) P$post_prob else P$statistic
  e[is.na(e)] <- 0
  e
}
.evid_lab <- function(evidence, P = NULL) {
  if (!is.null(P)) {
    e <- if (evidence == "p") P$p_adj else P$bf10
    if (all(is.na(e))) return(if (any(!is.na(P$post_prob))) "posterior probability" else "statistic")
  }
  if (evidence == "p") "-log10 adj. p" else "log10 BF10"
}

.plot_network <- function(x, evidence, main) {
  P <- x$pairs; vars <- x$series; k <- length(vars); s <- x$settings
  # layout: two series side by side, otherwise on a circle starting at the top
  xy <- if (k == 2L) rbind(c(-0.75, 0), c(0.75, 0)) else {
    ang <- pi / 2 - 2 * pi * (seq_len(k) - 1) / k; cbind(cos(ang), sin(ang)) }
  rownames(xy) <- vars
  op <- graphics::par(mar = c(1.5, 1, 3.6, 1), xpd = NA); on.exit(graphics::par(op))
  graphics::plot.new()
  graphics::plot.window(c(-1.5, 1.5), c(-1.5, 1.5), asp = 1)
  graphics::title(main = if (is.null(main)) "Causal network" else main, font.main = 2, line = 2)
  graphics::mtext(sprintf("%s;  link if %s;  n = %d", .method_label(s), .rule_label(s), s$n),
                  side = 3, line = 0.5, cex = 0.72, col = "grey30")
  sig <- P$significant %in% TRUE
  ev <- .evid(P, evidence); rng <- max(c(ev[sig], 1), na.rm = TRUE)
  rad <- if (k <= 3) 0.2 else 0.17
  both <- function(a, b) any(P$cause == b & P$effect == a & sig) && any(P$cause == a & P$effect == b & sig)
  edge <- function(i, strong) {
    a <- xy[P$cause[i], ]; b <- xy[P$effect[i], ]
    d <- b - a; L <- sqrt(sum(d^2)); u <- d / L; nrm <- c(-u[2], u[1])
    off <- if (both(P$cause[i], P$effect[i]) || !strong) 0.09 else 0
    a2 <- a + u * rad + nrm * off; b2 <- b - u * rad * 1.15 + nrm * off
    if (strong) {
      lwd <- 1.2 + 4 * ev[i] / rng
      col <- if (is.na(P$sign[i]) || P$sign[i] >= 0) .pal$pos else .pal$neg
      graphics::arrows(a2[1], a2[2], b2[1], b2[2], length = 0.13, lwd = lwd, col = col)
      m <- (a2 + b2) / 2 + nrm * 0.1
      lab <- sprintf("%s%s", if (is.na(P$lag[i])) "" else sprintf("lag %d, ", P$lag[i]), .evid_text(P, i))
      w <- graphics::strwidth(lab, cex = 0.62) / 2 + 0.03; h <- graphics::strheight(lab, cex = 0.62) / 2 + 0.025
      graphics::rect(m[1] - w, m[2] - h, m[1] + w, m[2] + h, col = grDevices::adjustcolor("white", 0.85), border = NA)
      graphics::text(m[1], m[2], lab, cex = 0.62, col = "#333333")
    } else {
      graphics::arrows(a2[1], a2[2], b2[1], b2[2], length = 0.08, lwd = 0.8, col = "grey78", lty = 3)
    }
  }
  if (k <= 5) for (i in which(!sig)) edge(i, FALSE)      # tested, not significant (light, dotted)
  for (i in which(sig)) edge(i, TRUE)
  Cd <- x$classes[x$classes$common_driver %in% TRUE, , drop = FALSE]
  for (i in seq_len(nrow(Cd))) {
    a <- xy[Cd$var1[i], ]; b <- xy[Cd$var2[i], ]
    graphics::segments(a[1], a[2], b[1], b[2], lty = 2, lwd = 2.2, col = "#E69F00")
    m <- (a + b) / 2; graphics::text(m[1], m[2] + 0.1, "hidden\ncommon driver?", cex = 0.6, col = "#B36B00")
  }
  deg_out <- vapply(vars, function(v) sum(P$cause == v & sig), 0)
  deg_in <- vapply(vars, function(v) sum(P$effect == v & sig), 0)
  role <- ifelse(deg_in + deg_out == 0, "isolated", ifelse(deg_out > 0 & deg_in == 0, "driver",
                 ifelse(deg_in > 0 & deg_out == 0, "driven", "mediator")))
  fill <- c(isolated = .pal$nocausal, driver = .pal$causal, driven = .pal$node, mediator = "#D7C6E8")[role]
  graphics::symbols(xy[, 1], xy[, 2], circles = rep(rad, k), inches = FALSE, add = TRUE, bg = fill, fg = .pal$border, lwd = 1.3)
  graphics::text(xy[, 1], xy[, 2], vars, cex = 0.9, font = 2, col = .pal$text)
  graphics::legend("bottomleft", bty = "n", cex = 0.62, inset = c(-0.02, -0.04),
                   legend = c("significant link (+ net effect)", "significant link (- net effect)", "tested, not significant",
                              "common-driver evidence"),
                   col = c(.pal$pos, .pal$neg, "grey70", "#E69F00"), lty = c(1, 1, 3, 2), lwd = c(2.5, 2.5, 1, 2))
  graphics::legend("bottomright", bty = "n", cex = 0.62, inset = c(-0.02, -0.04), legend = c("driver", "driven", "mediator", "isolated"),
                   pt.bg = c(.pal$causal, .pal$node, "#D7C6E8", .pal$nocausal), pch = 22, pt.cex = 1.5, col = .pal$border)
  graphics::mtext(sprintf("arrow = cause -> effect; width ~ %s; label = lag, evidence", .evid_lab(evidence, P)),
                  side = 1, line = 0.2, cex = 0.62, col = "grey40")
}

.method_label <- function(s) {
  switch(s$method, granger = "Granger causality (F-test)", bayes_granger = "Bayesian Granger causality (Bayes factor)",
         cmi = "kNN conditional mutual information", cmi_lag = "lag-shifted CMI (Zlatniczki et al. 2021)",
         markov = "Markov-chain approach (Stippinger et al. 2023)", df = "degrees of freedom (Telcs et al. 2025)",
         dc = "dimensional causality (Benko et al. 2024)", ccdh = "CCDH (Stippinger et al. 2023)",
         dfe = "df-entropy + random forest (Jakovac et al.)", te = "transfer entropy (Schreiber 2000)",
         ccm = "convergent cross mapping (Sugihara et al. 2012)", pid = "MAR-PID (Telcs 2026)",
         hirata = "recurrence plots (Hirata & Aihara 2010)", cd = "correlation dimension (Krakovska 2019)",
         pcmci = "PCMCI (Runge et al. 2019)", s$method)
}
.rule_label <- function(s) {
  if (s$method == "bayes_granger") sprintf("BF10 > %g", s$bf_threshold)
  else if (isTRUE(s$bayes)) {
    if (!is.null(s$bayes_control) && identical(s$bayes_control$decision, "posterior")) sprintf("posterior > %g", s$prob_threshold)
    else sprintf("BF10 > %g (Bayesian)", s$bf_threshold) }
  else if (s$method %in% c("dc", "dfe")) sprintf("P(case) > %g", s$prob_threshold)
  else if (s$method %in% c("hirata", "cd")) "verdict per pair"
  else if (s$method == "df") sprintf("Table 2 / adj. p < %g", s$alpha)
  else sprintf("adj. p < %g (%s)", s$alpha, s$adjust)
}
# short evidence text of a directed test for labels
.evid_text <- function(P, i) {
  if (!is.na(P$p_adj[i])) { f <- .fmt_p(P$p_adj[i]); if (startsWith(f, "<")) paste("p", f) else sprintf("p = %s", f) }
  else if (!is.na(P$post_prob[i])) sprintf("P = %.2f", P$post_prob[i])
  else if (!is.na(P$bf10[i])) sprintf("BF = %s", .fmt_bf(P$bf10[i]))
  else sprintf("%.2f", P$statistic[i])
}

.plot_heat <- function(x, evidence, main) {
  M <- as.matrix(x, if (evidence == "p") "p_adj" else "bf10")
  V <- if (evidence == "p") pmax(-log10(pmax(M, 1e-300)), 0) else log10(pmax(M, 1e-300))
  Sg <- as.matrix(x, "significant"); Lg <- as.matrix(x, "lag")
  k <- nrow(V); vars <- rownames(V)
  op <- graphics::par(mar = c(4.5, 5.5, 4, 5), xpd = NA); on.exit(graphics::par(op))
  zmax <- max(abs(V), 1, na.rm = TRUE)
  zlim <- if (evidence == "p") c(0, zmax) else c(-zmax, zmax)
  cols <- if (evidence == "p") grDevices::colorRampPalette(c("#FFFFFF", "#A8D5BA", "#1B7A43"))(64) else
    grDevices::colorRampPalette(c("#C62828", "#FFFFFF", "#1B7A43"))(64)
  z <- t(V)[, k:1, drop = FALSE]                     # x = effect, y = cause (top-down)
  graphics::image(seq_len(k), seq_len(k), z, col = cols, zlim = zlim, axes = FALSE,
                  xlab = "", ylab = "", main = "")
  graphics::title(main = if (is.null(main)) "Causal matrix: row series -> column series" else main,
                  font.main = 2, line = 2.2)
  graphics::mtext(sprintf("%s;  decision: %s", .method_label(x$settings), .rule_label(x$settings)),
                  side = 3, line = 0.6, cex = 0.75, col = "grey30")
  graphics::axis(1, seq_len(k), paste0("-> ", vars), las = 1, cex.axis = 0.85, tick = FALSE, line = -0.5)
  graphics::axis(2, seq_len(k), paste0(rev(vars), " ->"), las = 1, cex.axis = 0.85, tick = FALSE, line = -0.5)
  graphics::mtext("effect  (to)", side = 1, line = 2.6, cex = 0.85)
  graphics::mtext("cause  (from)", side = 2, line = 4, cex = 0.85)
  for (i in seq_len(k)) for (j in seq_len(k)) {      # i = cause (row), j = effect (col)
    xc <- j; yc <- k - i + 1
    if (i == j || is.na(V[i, j])) {
      graphics::rect(xc - .5, yc - .5, xc + .5, yc + .5, col = "grey90", border = "white")
      graphics::text(xc, yc, if (i == j) "-" else "n/a", col = "grey55", cex = 0.9)
    } else {
      sig <- isTRUE(Sg[i, j])
      dark <- abs(V[i, j]) / zmax > 0.55
      tcol <- if (dark) "white" else if (sig) "#1B4332" else "grey45"
      graphics::rect(xc - .5, yc - .5, xc + .5, yc + .5, border = if (sig) "#1B4332" else "white",
                     lwd = if (sig) 2.5 else 1)
      ii <- which(x$pairs$cause == vars[i] & x$pairs$effect == vars[j])[1]
      graphics::text(xc, yc + 0.2, sprintf("%s -> %s%s", vars[i], vars[j], if (sig) " *" else ""), cex = 0.72,
                     col = tcol, font = if (sig) 2 else 1)
      graphics::text(xc, yc - 0.02, if (is.na(ii)) "" else .evid_text(x$pairs, ii), cex = 0.66, col = tcol)
      graphics::text(xc, yc - 0.22, if (is.na(Lg[i, j])) "" else sprintf("lag %s", Lg[i, j]), cex = 0.62, col = tcol)
    }
  }
  # colour bar (positions as fractions of the plot width so that it fits any k)
  u <- graphics::par("usr"); w <- diff(u[1:2])
  nb <- length(cols); ys <- seq(0.5, k + 0.5, length.out = nb + 1)
  xb <- u[2] + 0.025 * w; bw <- 0.035 * w
  graphics::rect(xb, ys[-length(ys)], xb + bw, ys[-1], col = cols, border = NA)
  graphics::rect(xb, 0.5, xb + bw, k + 0.5, border = "grey40")
  ticks <- pretty(zlim, 4); ticks <- ticks[ticks >= zlim[1] & ticks <= zlim[2]]
  ty <- 0.5 + (ticks - zlim[1]) / diff(zlim) * k
  graphics::segments(xb + bw, ty, xb + bw + 0.01 * w, ty, col = "grey40")
  graphics::text(xb + bw + 0.02 * w, ty, ticks, adj = 0, cex = 0.65)
  graphics::text(xb + bw / 2, k + 0.55, .evid_lab(evidence), adj = c(0.5, 0), cex = 0.65)
  graphics::mtext("cell = directed test row -> column; * = significant; colour = evidence",
                  side = 1, line = 3.6, cex = 0.62, col = "grey40")
}

.plot_bars <- function(x, evidence, main) {
  P <- x$pairs; ev <- .evid(P, evidence)
  o <- order(ev); P <- P[o, , drop = FALSE]; ev <- ev[o]
  lab <- paste(P$cause, "->", P$effect, paste0("(lag ", P$lag, ")"))
  op <- graphics::par(mar = c(4, max(8, 0.55 * max(nchar(lab))), 3, 1)); on.exit(graphics::par(op))
  thr <- if (x$settings$method != "bayes_granger") -log10(x$settings$alpha) else log10(x$settings$bf_threshold)
  if (evidence == "p" && x$settings$method != "granger") thr <- NA
  if (evidence == "bf" && x$settings$method != "bayes_granger") thr <- NA
  lab <- sub(" \\(lag NA\\)", "", lab)
  graphics::barplot(ev, horiz = TRUE, names.arg = lab, las = 1, cex.names = 0.75,
                    col = ifelse(P$significant %in% TRUE, .pal$causal, .pal$nocausal), border = .pal$border,
                    xlab = .evid_lab(evidence, P), main = if (is.null(main)) "Evidence per directed pair" else main)
  if (is.finite(thr)) graphics::abline(v = thr, lty = 2, col = "grey40")
}

# statistic as a function of the candidate lag, one panel per directed pair
.plot_lags <- function(x, main) {
  B <- x$by_lag
  if (is.null(B) && !is.null(x$dc)) {           # dimensional causality / dfe: class scores per pair
    Dc <- x$dc
    M <- t(as.matrix(Dc[, c("x_to_y", "y_to_x", "bidirectional", "common_cause", "independent")]))
    colnames(M) <- paste(Dc$var1, "~", Dc$var2)
    op <- graphics::par(mar = c(4, 4.5, 3, 1)); on.exit(graphics::par(op))
    graphics::barplot(M, beside = TRUE, col = c("#1B7A43", "#C62828", "#6A1B9A", "#E69F00", "#BFD7EA"),
                      ylim = c(0, 1), ylab = "posterior probability", las = 1, cex.names = 0.8,
                      main = if (is.null(main)) "Posterior of the causal cases" else main, font.main = 2,
                      legend.text = c("x -> y", "y -> x", "x <-> y", "common cause", "independent"),
                      args.legend = list(x = "top", bty = "n", cex = 0.7, ncol = 3))
    return(invisible())
  }
  if (is.null(B)) { graphics::plot.new(); graphics::title("No lag scan stored for this method"); return(invisible()) }
  if (length(unique(B$lag)) == 1L && x$settings$method != "ccm") {     # a single lag: bars of the statistic per pair
    P <- x$pairs; o <- order(P$statistic); lab <- paste(P$cause[o], "->", P$effect[o])
    op <- graphics::par(mar = c(4.5, max(7, 0.6 * max(nchar(lab))), 3.5, 1)); on.exit(graphics::par(op))
    stat_lab <- switch(x$settings$method, granger = , bayes_granger = "F statistic", te = "transfer entropy (bits)",
                       df = "sd reduction", pid = "unique + synergistic information (bits)", "CMI (nats)")
    graphics::barplot(P$statistic[o], horiz = TRUE, names.arg = lab, las = 1, cex.names = 0.75, border = .pal$border,
                      col = ifelse(P$significant[o] %in% TRUE, .pal$causal, .pal$nocausal), xlab = stat_lab,
                      main = if (is.null(main)) sprintf("Statistic per directed pair (lag %d)", B$lag[1]) else main, font.main = 2)
    graphics::mtext(sprintf("%s; green = significant link", .method_label(x$settings)), side = 3, line = 0.3, cex = 0.72, col = "grey30")
    return(invisible())
  }
  P <- x$pairs; key <- paste(B$cause, "->", B$effect); pk <- paste(P$cause, "->", P$effect)
  pairs <- unique(key); np <- length(pairs)
  stat_lab <- switch(x$settings$method, granger = , bayes_granger = "F statistic", te = "transfer entropy (bits)",
                     ccm = "cross-map skill rho", df = "sd reduction", pid = "unique information (bits)", "CMI (nats)")
  xlab_ <- if (x$settings$method == "ccm") "library size" else "candidate lag"
  ylim <- range(c(0, B$statistic), na.rm = TRUE)
  sigp <- P$significant[match(pairs, pk)] %in% TRUE
  ns <- sum(sigp); cols <- rep("grey72", np); cols[sigp] <- grDevices::hcl.colors(max(ns, 2), "Dark 3")[seq_len(ns)]
  op <- graphics::par(mar = c(5, 4.5, 3.5, max(5, 1 + 0.55 * max(nchar(pairs))))); on.exit(graphics::par(op))
  xr <- range(B$lag); is_ccm <- x$settings$method == "ccm"
  graphics::plot(NA, xlim = xr, ylim = ylim, xlab = xlab_, ylab = stat_lab, xaxt = "n", log = if (is_ccm) "x" else "",
                 main = if (is.null(main)) (if (is_ccm) "Convergence: cross-map skill vs library size" else "Lag scan: statistic per directed pair") else main,
                 font.main = 2, las = 1)
  graphics::mtext(.method_label(x$settings), side = 3, line = 0.3, cex = 0.72, col = "grey30")
  graphics::axis(1, at = if (is_ccm) grDevices::axisTicks(log10(xr), log = TRUE) else sort(unique(B$lag)))
  graphics::grid(nx = NA, ny = NULL, col = "grey92"); graphics::abline(h = 0, col = "grey70")
  ends <- list()
  for (i in order(sigp)) {                         # grey (non-significant) lines first, coloured on top
    b <- B[key == pairs[i], , drop = FALSE]; b <- b[order(b$lag), ]
    r <- P[pk == pairs[i], ][1, ]
    graphics::lines(b$lag, b$statistic, col = cols[i], lwd = if (sigp[i]) 2.6 else 1.1, lty = if (sigp[i]) 1 else 3)
    graphics::points(b$lag, b$statistic, col = cols[i], pch = 16, cex = if (sigp[i]) 0.9 else 0.6)
    if (sigp[i] && !is.na(r$lag) && any(b$lag == r$lag))
      graphics::points(r$lag, b$statistic[b$lag == r$lag][1], col = cols[i], pch = 21, bg = "white", cex = 2, lwd = 2.2)
    ends[[i]] <- c(b$lag[nrow(b)], b$statistic[nrow(b)])
  }
  # labels at the right end of every line, spread out vertically to avoid overlaps
  ye <- vapply(ends, `[`, 0, 2); oy <- order(ye); yl <- ye
  gap <- 0.045 * diff(graphics::par("usr")[3:4])
  for (j in seq_along(oy)[-1]) if (yl[oy[j]] - yl[oy[j - 1]] < gap) yl[oy[j]] <- yl[oy[j - 1]] + gap
  xe <- if (is_ccm) 10^(log10(xr[2]) + 0.03 * diff(log10(xr))) else xr[2] + 0.03 * diff(xr)
  graphics::text(xe, yl, sprintf("%s%s", pairs, ifelse(sigp, " *", "")), adj = 0, cex = 0.68, xpd = NA,
                 col = ifelse(sigp, cols, "grey45"), font = ifelse(sigp, 2, 1))
  graphics::mtext("coloured = significant (*), grey dotted = not;  circle = selected lag",
                  side = 1, line = 3.8, cex = 0.62, col = "grey40")
}

#' Time-delay (Takens) embedding plot
#'
#' Draws the delay-coordinate reconstruction of a series: `x[t]` against
#' `x[t - tau]` (`m = 2`) or the first three delay coordinates in an
#' orthographic 3-D projection (`m >= 3`), the trajectory drawn as a thin line
#' coloured by time. For two series, `y` is used as the second coordinate
#' (joint state `(x[t], y[t])` of Stippinger et al., 2023).
#'
#' @param x numeric vector.
#' @param y optional second series (then the joint embedding `(x[t], y[t],
#'   x[t - tau])` is drawn).
#' @param m embedding dimension (2 or 3 coordinates are drawn).
#' @param tau delay.
#' @param theta,phi viewing angles (degrees) of the 3-D projection.
#' @param main title.
#' @param col point colours (default: a time gradient).
#' @param trajectory draw the trajectory as a thin line (useful for flows,
#'   noisy for maps).
#' @param note optional one-line note printed under the plot.
#' @param ... passed to [graphics::plot()].
#' @return Invisibly, the embedded matrix.
#' @references
#' Takens, F. (1981). Detecting strange attractors in turbulence. In: *Dynamical Systems and Turbulence, Warwick 1980*, Lecture Notes in Mathematics 898, 366-381. Springer. \doi{10.1007/BFb0091924}
#'
#' Stippinger, M., Bencze, A., Zlatniczki, A., Somogyvari, Z. and Telcs, A. (2023). Causal discovery of stochastic dynamical systems: a Markov chain approach. *Mathematics*, 11(4), 852. \doi{10.3390/math11040852}
#' @export
#' @examples
#' # the reconstructed attractors of two coupled logistic maps: the driver x
#' # draws the parabola of the map, the driven y a thicker, deformed one
#' d <- sim_logistic(500, "x_to_y", seed = 1)
#' plot_embedding(d$x, m = 2)
#' \donttest{
#' plot_embedding(d$y, m = 3)
#' plot_embedding(d$x, d$y)      # joint state of the driver and the driven series
#' }
plot_embedding <- function(x, y = NULL, m = 3L, tau = 1L, theta = 35, phi = 25, main = NULL,
                           col = NULL, trajectory = FALSE, note = NULL, ...) {
  x <- as.numeric(x)
  E <- if (is.null(y)) embed_delay(x, max(2L, as.integer(m)), tau) else {
    y <- as.numeric(y); n <- length(x); i <- (tau + 1):n
    cbind(x[i - tau], y[i], x[i])          # (cause at t - tau, effect at t, cause at t)
  }
  n <- nrow(E)
  if (is.null(col)) col <- grDevices::hcl.colors(n, "Viridis", alpha = 0.8)
  labs <- if (is.null(y)) c("x[t]", sprintf("x[t-%d]", tau), sprintf("x[t-%d]", 2 * tau)) else
    c(sprintf("x[t-%d]", tau), "y[t]", "x[t]")
  if (is.null(main)) main <- if (is.null(y)) sprintf("Delay embedding (m = %d, tau = %d)", ncol(E), tau) else
    "Joint state embedding (x, y)"
  time_bar <- function() {                       # colour = time, from the first to the last point
    u <- graphics::par("usr"); w <- diff(u[1:2]); h <- diff(u[3:4])
    cb <- grDevices::hcl.colors(40, "Viridis"); xs <- seq(u[2] - 0.32 * w, u[2] - 0.04 * w, length.out = 41)
    graphics::rect(xs[-41], u[4] - 0.06 * h, xs[-1], u[4] - 0.035 * h, col = cb, border = NA, xpd = NA)
    graphics::text(c(xs[1], xs[41]), u[4] - 0.075 * h, c("start", "end"), cex = 0.58, col = "grey35", adj = c(0.5, 1))
    graphics::text(mean(xs), u[4] - 0.03 * h, "colour = time", cex = 0.58, col = "grey35", adj = c(0.5, 0))
  }
  if (ncol(E) == 2L || (is.null(y) && m == 2L)) {
    graphics::plot(E[, 2], E[, 1], type = "n", xlab = labs[2], ylab = labs[1], main = main, font.main = 2, las = 1, ...)
    if (trajectory) graphics::lines(E[, 2], E[, 1], col = "grey75", lwd = 0.5)
    graphics::points(E[, 2], E[, 1], col = col, pch = 16, cex = 0.6)
    if (is.null(match.call()$col)) time_bar()
    return(invisible(E))
  }
  # orthographic projection of the first three coordinates
  Z <- scale(E[, 1:3])
  t1 <- theta * pi / 180; p1 <- phi * pi / 180
  R <- matrix(c(cos(t1), -sin(t1), 0, sin(t1), cos(t1), 0, 0, 0, 1), 3, byrow = TRUE) %*%
       matrix(c(1, 0, 0, 0, cos(p1), -sin(p1), 0, sin(p1), cos(p1)), 3, byrow = TRUE)
  Pz <- Z %*% t(R)
  ax <- diag(3) * 2.2; Ax <- ax %*% t(R); o <- c(0, 0, 0) %*% t(R)
  xr <- range(c(Pz[, 1], Ax[, 1] * 1.25)); yr <- range(c(Pz[, 3], Ax[, 3] * 1.25))
  op <- graphics::par(mar = c(if (is.null(note)) 1 else 2, 1, 3, 1), xpd = FALSE); on.exit(graphics::par(op))
  graphics::plot(Pz[, 1], Pz[, 3], type = "n", axes = FALSE, xlab = "", ylab = "", main = main,
                 font.main = 2, asp = 1, xlim = xr, ylim = yr, ...)
  if (!is.null(note)) graphics::mtext(note, side = 1, line = 0.6, cex = 0.62, col = "grey40")
  for (a in 1:3) {
    graphics::arrows(o[1], o[3], Ax[a, 1], Ax[a, 3], length = 0.08, col = "grey50")
    graphics::text(Ax[a, 1] * 1.14, Ax[a, 3] * 1.14, labs[a], cex = 0.8, col = "grey30")
  }
  if (trajectory) graphics::lines(Pz[, 1], Pz[, 3], col = "grey75", lwd = 0.5)
  graphics::points(Pz[, 1], Pz[, 3], col = col, pch = 16, cex = 0.55)
  if (is.null(match.call()$col)) time_bar()
  invisible(E)
}

.plot_embed_all <- function(x, main, single = FALSE, ...) {
  s <- x$settings; D <- x$data; vars <- x$series; k <- length(vars)
  if (single) {                      # one panel: joint state of the strongest directed pair
    P <- x$pairs
    i <- order(ifelse(is.na(P$p_adj), 1, P$p_adj), -ifelse(is.na(P$statistic), 0, P$statistic))[1]
    tau <- if (is.na(P$lag[i])) max(1L, as.integer(.call_arg(x$call, "tau", 1L))) else max(1L, as.integer(P$lag[i]))
    plot_embedding(D[[P$cause[i]]], D[[P$effect[i]]], tau = tau,
                   main = sprintf("Joint state: %s (cause) x %s (effect), lag %d%s", P$cause[i], P$effect[i], tau,
                                  if (P$significant[i]) " *" else ""),
                   note = if (s$method %in% c("dc", "ccm", "ccdh", "hirata", "cd"))
                     "deterministic coupling: the driven state determines the driver (thin surface)" else
                     "stochastic coupling: dependence of the effect on the lagged cause", ...)
    return(invisible(x))
  }
  lag_max <- suppressWarnings(max(x$pairs$lag, na.rm = TRUE)); if (!is.finite(lag_max)) lag_max <- 2L
  m <- if (!is.null(x$markov)) max(x$markov$memory) else
       if (s$method %in% c("dc", "ccm", "ccdh", "hirata", "cd")) max(2L, as.integer(.call_arg(x$call, "m", 3L))) else
       max(2L, min(3L, lag_max + 1L))
  if (k > 1) { op <- graphics::par(mfrow = grDevices::n2mfrow(k), mar = c(4, 4, 3, 1)); on.exit(graphics::par(op)) }
  for (v in vars) plot_embedding(D[[v]], m = m, main = sprintf("%s: delay embedding (m = %d)", v, m), ...)
  invisible(x)
}
.call_arg <- function(cl, name, default) {
  v <- tryCatch(eval(cl[[name]], parent.frame(2)), error = function(e) NULL)
  if (is.null(v) || !is.numeric(v)) default else v
}


#' Plot several time series in stacked panels
#'
#' A quick look at the data before a causality analysis: one panel per
#' series, shared time axis, optional highlighting of a window.
#'
#' @param x a data frame or matrix with the series in columns, or a numeric
#'   vector.
#' @param n_last draw only the last `n_last` observations (`NULL` = all).
#' @param col line colours (recycled).
#' @param main overall title.
#' @param ... further arguments to [graphics::plot()].
#' @return `x`, invisibly.
#' @export
#' @examples
#' # the last 150 observations of three simulated series, one panel each
#' plot_series(sim_var(300, n_series = 3, seed = 1), n_last = 150)
plot_series <- function(x, n_last = NULL, col = NULL, main = NULL, ...) {
  D <- if (is.null(dim(x))) data.frame(x = as.numeric(x)) else as.data.frame(x)
  n <- nrow(D); idx <- if (is.null(n_last)) seq_len(n) else max(1, n - n_last + 1):n
  k <- ncol(D); if (is.null(col)) col <- grDevices::hcl.colors(max(k, 2), "Dark 3")[seq_len(k)]
  op <- graphics::par(mfrow = c(k, 1), mar = c(0.6, 4.2, 0.6, 1), oma = c(4, 0, if (is.null(main)) 0.5 else 2.5, 0))
  on.exit(graphics::par(op))
  for (j in seq_len(k)) {
    graphics::plot(idx, D[[j]][idx], type = "l", col = col[(j - 1) %% length(col) + 1], lwd = 1.3, xaxt = "n",
                   xlab = "", ylab = names(D)[j], las = 1, ...)
    graphics::grid(nx = NA, ny = NULL, col = "grey90")
  }
  graphics::axis(1); graphics::mtext("time", side = 1, line = 2.5)
  if (!is.null(main)) graphics::mtext(main, outer = TRUE, font = 2, cex = 1.1)
  invisible(x)
}
