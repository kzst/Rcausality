.pal <- list(causal = "#A8D5BA", nocausal = "#F4A6A6", node = "#BFD7EA",
             border = "#555555", pos = "#2E7D32", neg = "#C62828",
             grid = "grey85", text = "#212121")

#' Plot a `causres` object
#'
#' @param x a `causres` object.
#' @param type `"network"` (directed graph of significant links; width ~
#'   evidence, label = lag, colour = sign of the net influence), `"heatmap"`
#'   (cause x effect matrix of evidence, stars mark significant cells),
#'   `"bars"` (evidence per directed pair, sorted) or `"both"`
#'   (network and heatmap side by side).
#' @param evidence evidence shown for widths/colours: `"p"`
#'   (-log10 adjusted p) or `"bf"` (log10 Bayes factor).
#' @param main optional title.
#' @param ... unused.
#' @return Invisibly, `x`.
#' @export
plot.causres <- function(x, type = c("network", "heatmap", "bars", "both"),
                         evidence = c("p", "bf"), main = NULL, ...) {
  type <- match.arg(type); evidence <- match.arg(evidence)
  if (type == "both") {
    op <- graphics::par(mfrow = c(1, 2)); on.exit(graphics::par(op))
    .plot_network(x, evidence, main); .plot_heat(x, evidence, main)
  } else switch(type, network = .plot_network(x, evidence, main),
                heatmap = .plot_heat(x, evidence, main),
                bars = .plot_bars(x, evidence, main))
  invisible(x)
}

.evid <- function(P, evidence) {
  if (evidence == "p") -log10(pmax(P$p_adj, 1e-300)) else log10(pmax(P$bf10, 1e-300))
}
.evid_lab <- function(evidence) if (evidence == "p") "-log10 adj. p" else "log10 BF10"

.plot_network <- function(x, evidence, main) {
  P <- x$pairs; vars <- x$series; k <- length(vars)
  ang <- pi / 2 - 2 * pi * (seq_len(k) - 1) / k
  xy <- cbind(cos(ang), sin(ang)); rownames(xy) <- vars
  op <- graphics::par(mar = c(1, 1, 3, 1)); on.exit(graphics::par(op))
  graphics::plot.new()
  graphics::plot.window(c(-1.45, 1.45), c(-1.45, 1.45), asp = 1)
  graphics::title(main = if (is.null(main)) "Causal network" else main, font.main = 2)
  S <- P[P$significant, , drop = FALSE]
  ev <- .evid(P, evidence); evs <- ev[P$significant]
  rng <- max(c(evs, 1), na.rm = TRUE)
  rad <- 0.17
  both <- function(a, b) any(S$cause == b & S$effect == a)
  for (i in seq_len(nrow(S))) {
    a <- xy[S$cause[i], ]; b <- xy[S$effect[i], ]
    d <- b - a; L <- sqrt(sum(d^2)); u <- d / L; nrm <- c(-u[2], u[1])
    off <- if (both(S$cause[i], S$effect[i])) 0.08 else 0
    a2 <- a + u * rad + nrm * off; b2 <- b - u * rad * 1.15 + nrm * off
    lwd <- 1 + 4 * evs[i] / rng
    col <- if (is.na(S$sign[i]) || S$sign[i] >= 0) .pal$pos else .pal$neg
    graphics::arrows(a2[1], a2[2], b2[1], b2[2], length = 0.12, lwd = lwd, col = col)
    m <- (a2 + b2) / 2 + nrm * 0.07
    graphics::text(m[1], m[2], S$lag[i], cex = 0.75, col = "#444444")
  }
  deg_out <- vapply(vars, function(v) sum(S$cause == v), 0)
  deg_in <- vapply(vars, function(v) sum(S$effect == v), 0)
  fill <- ifelse(deg_in + deg_out == 0, .pal$nocausal,
          ifelse(deg_out > 0 & deg_in == 0, .pal$causal, .pal$node))
  graphics::symbols(xy[, 1], xy[, 2], circles = rep(rad, k), inches = FALSE,
                    add = TRUE, bg = fill, fg = .pal$border)
  graphics::text(xy[, 1], xy[, 2], vars, cex = 0.85, col = .pal$text)
  graphics::legend("bottomleft", bty = "n", cex = 0.7,
                   legend = c("positive net effect", "negative net effect",
                              "driver", "target/mediator", "isolated"),
                   col = c(.pal$pos, .pal$neg, NA, NA, NA), lty = c(1, 1, NA, NA, NA),
                   pt.bg = c(NA, NA, .pal$causal, .pal$node, .pal$nocausal),
                   pch = c(NA, NA, 22, 22, 22), pt.cex = 1.4)
  graphics::mtext(sprintf("edge width ~ %s; label = lag; %s", .evid_lab(evidence), x$settings$method),
                  side = 1, line = -0.2, cex = 0.7, col = "grey40")
}

.plot_heat <- function(x, evidence, main) {
  M <- as.matrix(x, if (evidence == "p") "p_adj" else "bf10")
  V <- if (evidence == "p") pmax(-log10(pmax(M, 1e-300)), 0) else log10(pmax(M, 1e-300))
  Sg <- as.matrix(x, "significant")
  k <- nrow(V)
  op <- graphics::par(mar = c(5, 5, 3.5, 5)); on.exit(graphics::par(op))
  zmax <- max(abs(V), 1, na.rm = TRUE)
  zlim <- if (evidence == "p") c(0, zmax) else c(-zmax, zmax)
  cols <- if (evidence == "p") grDevices::colorRampPalette(c("#FFFFFF", "#A8D5BA", "#1B7A43"))(64) else
    grDevices::colorRampPalette(c("#C62828", "#FFFFFF", "#1B7A43"))(64)
  graphics::image(seq_len(k), seq_len(k), V[, k:1, drop = FALSE], col = cols, zlim = zlim,
                  axes = FALSE, xlab = "cause", ylab = "effect",
                  main = if (is.null(main)) "Evidence matrix" else main)
  graphics::axis(1, seq_len(k), rownames(V), las = 2, cex.axis = 0.8)
  graphics::axis(2, seq_len(k), rev(colnames(V)), las = 1, cex.axis = 0.8)
  graphics::box()
  for (i in seq_len(k)) for (j in seq_len(k)) {
    if (is.na(V[i, j])) {
      graphics::rect(i - .5, k - j + 0.5, i + .5, k - j + 1.5, col = "grey92", border = NA)
    } else {
      lab <- sprintf("%.1f%s", V[i, j] + 0, if (isTRUE(Sg[i, j])) "*" else "")
      graphics::text(i, k - j + 1, lab, cex = 0.75, font = if (isTRUE(Sg[i, j])) 2 else 1)
    }
  }
  graphics::mtext(paste(.evid_lab(evidence), "(* = significant)"), side = 1, line = 4, cex = 0.7)
}

.plot_bars <- function(x, evidence, main) {
  P <- x$pairs; ev <- .evid(P, evidence)
  o <- order(ev); P <- P[o, , drop = FALSE]; ev <- ev[o]
  lab <- paste(P$cause, "->", P$effect, paste0("(lag ", P$lag, ")"))
  op <- graphics::par(mar = c(4, max(8, 0.55 * max(nchar(lab))), 3, 1)); on.exit(graphics::par(op))
  thr <- if (x$settings$method != "bayes_granger") -log10(x$settings$alpha) else log10(x$settings$bf_threshold)
  if (evidence == "p" && x$settings$method != "granger") thr <- NA
  if (evidence == "bf" && x$settings$method != "bayes_granger") thr <- NA
  graphics::barplot(ev, horiz = TRUE, names.arg = lab, las = 1, cex.names = 0.75,
                    col = ifelse(P$significant, .pal$causal, .pal$nocausal), border = .pal$border,
                    xlab = .evid_lab(evidence), main = if (is.null(main)) "Evidence per directed pair" else main)
  if (is.finite(thr)) graphics::abline(v = thr, lty = 2, col = "grey40")
}
