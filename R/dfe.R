# Degrees-of-freedom / conditional-entropy method with random-forest
# classification (Jakovac, Fabo, Eross, Stippinger, Benko, Somogyvari & Telcs,
# "Causality beyond the Reichenbach principle: the degrees-of-freedom method",
# manuscript).  Series are binarised at their mean; conditional entropies
# H(target[t+1] | own history of length l, other series' history of length k)
# are computed for l = 0..Lmax, k = 0..Kmax, for the targets X, Y and the joint
# (X, Y), on the original pair and on an independent (time-shifted) surrogate;
# curvature features of the entropy curves feed a random forest trained on
# simulated coupled binary Markov chains of the five causal classes.

.dfe_classes <- c("x_to_y", "y_to_x", "common_cause", "bidirectional", "independent")

.binarize_mean <- function(v) as.integer(v > mean(v))

# history matrix: columns = lags 0..l-1 of a (0/1), then lags 0..k-1 of b, aligned to t_idx
.hist_bits <- function(a, b, l, k, t_idx) {
  cols <- c(lapply(seq_len(l) - 1L, function(j) a[t_idx - j]), lapply(seq_len(k) - 1L, function(j) b[t_idx - j]))
  if (!length(cols)) return(matrix(0L, length(t_idx), 0))
  M <- do.call(cbind, cols); storage.mode(M) <- "integer"; M
}

# conditional-entropy surfaces H[l + 1, k + 1] for the targets X, Y and the joint
.ce_curves <- function(bx, by, Lmax = 6L, Kmax = 5L) {
  n <- length(bx); L <- max(Lmax, Kmax)
  if (length(by) != n) stop("x and y must have the same length.", call. = FALSE)
  if (anyNA(bx) || anyNA(by)) stop("Missing values in the series (impute first, e.g. with prep_series()).", call. = FALSE)
  if (n < 4L * L + 20L) stop(sprintf("Series too short for the df-entropy method: at least %d points are needed for Lmax = %d, Kmax = %d.", 4L * L + 20L, Lmax, Kmax), call. = FALSE)
  t_idx <- (L + 1):(n - 1)
  tx <- bx[t_idx + 1L]; ty <- by[t_idx + 1L]; tj <- 2L * tx + ty
  HX <- HY <- HJ <- matrix(NA_real_, Lmax + 1L, Kmax + 1L)
  for (l in 0:Lmax) for (k in 0:Kmax) {
    HX[l + 1L, k + 1L] <- cond_entropy_bits_cpp(tx, .hist_bits(bx, by, l, k, t_idx), 2L)
    HY[l + 1L, k + 1L] <- cond_entropy_bits_cpp(ty, .hist_bits(by, bx, l, k, t_idx), 2L)
    HJ[l + 1L, k + 1L] <- cond_entropy_bits_cpp(tj, .hist_bits(bx, by, l, k, t_idx), 4L)
  }
  list(X = HX, Y = HY, J = HJ)
}

.top2_curv <- function(h) {            # two largest |second differences| along l
  if (length(h) < 3) return(c(0, 0))
  d2 <- abs(diff(h, differences = 2)); d2 <- sort(d2, decreasing = TRUE)
  c(d2[1], if (length(d2) > 1) d2[2] else 0)
}

# feature vector of a binary pair (Appendix F.0.5): curvature of every curve,
# the same for the surrogate-corrected curves, and asymmetry / coupling-gain
# measures
.dfe_features <- function(bx, by, Lmax = 6L, Kmax = 5L) {
  n <- length(bx); s <- 10L * Lmax + 1L
  H <- .ce_curves(bx, by, Lmax, Kmax)
  # independent surrogate: x shifted forward, y backward by s (circularly)
  sh <- function(v, d) v[((seq_len(n) + d - 1L) %% n) + 1L]
  Hs <- .ce_curves(sh(bx, s), sh(by, -s), Lmax, Kmax)
  f <- c()
  for (nm in c("X", "Y", "J")) {
    D <- H[[nm]] - Hs[[nm]]
    for (k in 0:Kmax) {
      f[paste0("c_", nm, "_k", k, "_", 1:2)] <- .top2_curv(H[[nm]][, k + 1L])
      f[paste0("dc_", nm, "_k", k, "_", 1:2)] <- .top2_curv(D[, k + 1L])
    }
    f[paste0("gain_", nm)] <- H[[nm]][Lmax + 1L, 1L] - H[[nm]][Lmax + 1L, Kmax + 1L]   # other series' past helps?
    f[paste0("dgain_", nm)] <- D[Lmax + 1L, 1L] - D[Lmax + 1L, Kmax + 1L]
    f[paste0("rate_", nm)] <- H[[nm]][Lmax + 1L, Kmax + 1L]
    f[paste0("drate_", nm)] <- D[Lmax + 1L, Kmax + 1L]
    f[paste0("dmean_", nm)] <- mean(D)
    f[paste0("elbow_", nm)] <- which.max(abs(diff(H[[nm]][, 1L], differences = 2)))
    # the surrogate-corrected surface itself (entropy reduction attributable to coupling)
    f[paste0("d_", nm, "_l", rep(0:Lmax, Kmax + 1L), "_k", rep(0:Kmax, each = Lmax + 1L))] <- as.vector(D)
    # own-history profile (shape of the self-entropy decay)
    f[paste0("h_", nm, "_l", 0:Lmax)] <- H[[nm]][, 1L]
  }
  f["asym_gain"] <- f[["gain_X"]] - f[["gain_Y"]]
  f["asym_dgain"] <- f[["dgain_X"]] - f[["dgain_Y"]]
  f["asym_rate"] <- f[["rate_X"]] - f[["rate_Y"]]
  f["joint_excess"] <- f[["rate_J"]] - f[["rate_X"]] - f[["rate_Y"]]
  f["log_n"] <- log(n)
  f
}

.rand_p <- function(k) { p <- stats::runif(k, 0.001, 0.98); ifelse(p < 0.49, p, p + 0.02) }

#' Coupled binary Markov chains (training systems of the df-entropy method)
#'
#' Simulates the reference systems of Jakovac et al.: binary chains `x`, `y`
#' (and a hidden `z`) whose transition probability to state 1 is a random
#' table over the own history (memory `df_x`, `df_y`, `df_z`) and, when
#' coupled, the current state of the driver (the driver's state selects
#' between two transition tables, as in the paper). Transition probabilities
#' are drawn uniformly from `[0.001, 0.49] U [0.51, 0.99]`.
#'
#' @param n length.
#' @param class causal class: `"x_to_y"`, `"y_to_x"`, `"bidirectional"`,
#'   `"common_cause"` (hidden `z` drives both) or `"independent"`.
#' @param df_x,df_y,df_z memories (number of own past states) of the chains.
#' @param seed optional RNG seed.
#' @return A data frame with integer columns `x`, `y` (and `z` for the
#'   common-cause class) and `attr(, "truth")`.
#' @references
#' Jakovac, A., Fabo, D., Eross, L., Stippinger, M., Benko, Z., Somogyvari, Z. and Telcs, A. Causality beyond the Reichenbach principle: the degrees-of-freedom method (manuscript), Sect. 3 and Appendix F.0.6.
#' @export
#' @examples
#' # binary chains where x drives y; the class is recorded as an attribute
#' d <- sim_markov_pair(500, "x_to_y", seed = 1)
#' head(d)
#' attr(d, "truth")
#' table(x = d$x, y = d$y)      # the coupling shows in the joint frequencies
sim_markov_pair <- function(n = 2000, class = c("x_to_y", "y_to_x", "bidirectional", "common_cause",
                                                 "independent"),
                            df_x = sample(1:2, 1), df_y = sample(1:2, 1), df_z = 2L, seed = NULL) {
  class <- match.arg(class)
  .local_seed(seed)
  y2x <- class %in% c("y_to_x", "bidirectional"); x2y <- class %in% c("x_to_y", "bidirectional")
  cc <- class == "common_cause"
  px <- .rand_p(2^(df_x + y2x + cc)); py <- .rand_p(2^(df_y + x2y + cc)); pz <- .rand_p(2^df_z)
  M <- sim_coupled_binary_cpp(as.integer(n), as.integer(df_x), as.integer(df_y), as.integer(df_z),
                              y2x, x2y, cc, cc, px, py, pz, stats::runif(3L * n))
  d <- data.frame(x = M[, 1], y = M[, 2])
  if (cc) d$z <- M[, 3]
  attr(d, "truth") <- class
  d
}

.dfe_env <- new.env(parent = emptyenv())

#' Train the random-forest classifier of the df-entropy method
#'
#' Generates coupled binary Markov chains of the five causal classes
#' ([sim_markov_pair()]), extracts the conditional-entropy features and trains
#' a random forest (`ranger`). The model is cached for the session and used by
#' [dfe_classify()] and `rcausality(method = "dfe")` when no model is given.
#' The paper uses 1000 realisations per class with lengths 200-5000 and 100
#' trees; the default here is smaller for speed -- raise `n_per_class` for a
#' more accurate classifier.
#'
#' @param n_per_class number of training series per class.
#' @param lengths series lengths sampled for the training systems.
#' @param Lmax,Kmax largest self- and cross-history lengths of the entropy
#'   curves (paper: 6 and 5).
#' @param num_trees number of trees.
#' @param seed optional RNG seed.
#' @param cores CPU cores for the feature extraction and the forest (see
#'   [Rcausality-parallel]).
#' @return An object of class `dfe_model` (the `ranger` forest, the feature
#'   names and settings, the out-of-bag Brier error and the in-sample
#'   confusion matrix).
#' @references
#' Jakovac, A., Fabo, D., Eross, L., Stippinger, M., Benko, Z., Somogyvari, Z. and Telcs, A. Causality beyond the Reichenbach principle: the degrees-of-freedom method (manuscript), Sect. 2.2-2.3 and Appendix F.
#'
#' Wright, M. N. and Ziegler, A. (2017). ranger: A fast implementation of random forests for high dimensional data in C++ and R. *Journal of Statistical Software*, 77(1), 1-17. \doi{10.18637/jss.v077.i01}
#' @export
#' @examples
#' \donttest{
#' # a small model (25 training systems per class) trains in a few seconds;
#' # the print shows the out-of-bag error and the in-sample confusion matrix
#' m <- dfe_train(n_per_class = 25, lengths = c(300, 600), seed = 1)
#' m
#' # reuse it in the classifier and in rcausality()
#' d <- sim_markov_pair(1000, "x_to_y", seed = 2)
#' dfe_classify(d$x, d$y, model = m)$best
#' rcausality(d, method = "dfe", model = m)$classes
#' }
dfe_train <- function(n_per_class = 200L, lengths = c(200L, 500L, 1000L, 2000L, 5000L), Lmax = 6L,
                      Kmax = 5L, num_trees = 100L, seed = NULL, cores = getOption("Rcausality.cores", 1L)) {
  if (!requireNamespace("ranger", quietly = TRUE))
    stop("Package 'ranger' is required for the df-entropy method (install.packages(\"ranger\")).", call. = FALSE)
  old <- .rng_setup(seed); on.exit(.rng_restore(old), add = TRUE)
  design <- expand.grid(class = .dfe_classes, i = seq_len(n_per_class), stringsAsFactors = FALSE)
  lengths <- as.integer(lengths)
  design$n <- lengths[sample.int(length(lengths), nrow(design), replace = TRUE)]
  rows <- .plapply(seq_len(nrow(design)), function(i) {
    d <- sim_markov_pair(design$n[i], design$class[i])
    .dfe_features(d$x, d$y, Lmax, Kmax)
  }, cores = cores)
  Fm <- as.data.frame(do.call(rbind, rows)); Fm$class <- factor(design$class, levels = .dfe_classes)
  rf <- ranger::ranger(class ~ ., data = Fm, num.trees = num_trees, probability = TRUE,
                       num.threads = max(1L, .n_cores(cores)))
  oob <- stats::predict(rf, data = Fm)$predictions   # in-sample; the OOB error is in rf$prediction.error
  structure(list(rf = rf, features = setdiff(names(Fm), "class"), Lmax = Lmax, Kmax = Kmax,
                 n_per_class = n_per_class, lengths = lengths, oob_error = rf$prediction.error,
                 confusion = table(truth = Fm$class, predicted = factor(.dfe_classes[apply(oob, 1, which.max)], .dfe_classes))),
            class = "dfe_model")
}

#' @export
print.dfe_model <- function(x, ...) {
  cat("df-entropy random-forest model (Jakovac et al.)\n")
  cat(sprintf("  %d training series per class, lengths %s, Lmax = %d, Kmax = %d, %d trees, %d features\n",
              x$n_per_class, paste(x$lengths, collapse = "/"), x$Lmax, x$Kmax, x$rf$num.trees, length(x$features)))
  cat(sprintf("  out-of-bag prediction error (Brier) %.3f; in-sample confusion matrix:\n", x$oob_error))
  print(x$confusion)
  invisible(x)
}

.dfe_model <- function(model = NULL, ...) {
  if (!is.null(model)) return(model)
  if (is.null(.dfe_env$model)) {
    message("dfe: training the default random forest (200 series per class); pass `model = dfe_train(...)` to reuse or enlarge it.")
    .dfe_env$model <- dfe_train(...)
  }
  .dfe_env$model
}

#' Classify the causal relation of a pair with the df-entropy method
#'
#' Binarises `x` and `y` at their means, computes the conditional-entropy
#' features (see [dfe_train()]) and returns the class scores of the random
#' forest for the five causal classes.
#'
#' @param x,y numeric vectors.
#' @param model a `dfe_model` from [dfe_train()]; `NULL` uses (and trains if
#'   necessary) the session's default model.
#' @param lag shift applied to `y` (`y[t + lag]` paired with `x[t]`).
#' @param ... passed to [dfe_train()] when the default model is trained.
#' @return An object of class `dfe`: `scores` (named vector over
#'   `x_to_y`, `y_to_x`, `bidirectional`, `common_cause`, `independent`),
#'   `best`, `features`, `entropy` (the conditional-entropy surfaces of X, Y
#'   and the joint). Methods: `print`, `plot` (entropy curves, Fig. 1 of the
#'   paper).
#' @references
#' Jakovac, A., Fabo, D., Eross, L., Stippinger, M., Benko, Z., Somogyvari, Z. and Telcs, A. Causality beyond the Reichenbach principle: the degrees-of-freedom method (manuscript).
#' @export
#' @examples
#' \donttest{
#' # train a small forest, then classify a pair where x drives y: the class
#' # scores of the five causal cases should favour x_to_y
#' m <- dfe_train(n_per_class = 25, lengths = c(300, 600), seed = 1)
#' d <- sim_markov_pair(1000, "x_to_y", seed = 2)
#' r <- dfe_classify(d$x, d$y, model = m)
#' r
#' r$features    # the conditional-entropy features the forest sees
#' plot(r)       # the entropy curves and the class scores
#' }
dfe_classify <- function(x, y, model = NULL, lag = 0L, ...) {
  if (!requireNamespace("ranger", quietly = TRUE))
    stop("Package 'ranger' is required for the df-entropy method (install.packages(\"ranger\")).", call. = FALSE)
  model <- .dfe_model(model, ...)
  x <- as.numeric(x); y <- as.numeric(y); n <- length(x); lag <- as.integer(lag)
  if (lag != 0L) { i <- seq_len(n); j <- i + lag; ok <- j >= 1 & j <= n; x <- x[i[ok]]; y <- y[j[ok]] }
  bx <- .binarize_mean(x); by <- .binarize_mean(y)
  f <- .dfe_features(bx, by, model$Lmax, model$Kmax)
  Fm <- as.data.frame(as.list(f))[, model$features, drop = FALSE]
  sc <- stats::predict(model$rf, data = Fm)$predictions[1, ]
  sc <- sc[.dfe_classes]
  structure(list(scores = sc, best = names(sc)[which.max(sc)], features = f,
                 entropy = .ce_curves(bx, by, model$Lmax, model$Kmax), n = length(bx), lag = lag),
            class = "dfe")
}

#' @export
print.dfe <- function(x, digits = 3, ...) {
  cat("df-entropy method (Jakovac et al.): random-forest class scores\n")
  cat(sprintf("  n = %d (binarised at the mean), lag = %d\n", x$n, x$lag))
  lab <- c(x_to_y = "X -> Y", y_to_x = "X <- Y", bidirectional = "X <-> Y",
           common_cause = "hidden common cause", independent = "independent")
  print(data.frame(case = lab[names(x$scores)], score = round(x$scores, digits),
                   best = ifelse(names(x$scores) == x$best, "<==", "")), row.names = FALSE, right = FALSE)
  invisible(x)
}

#' @rdname dfe_classify
#' @param main title.
#' @export
plot.dfe <- function(x, main = NULL, ...) {
  E <- x$entropy; Kmax <- ncol(E$X) - 1L; Lmax <- nrow(E$X) - 1L
  op <- graphics::par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3.5, 1)); on.exit(graphics::par(op))
  cols <- grDevices::hcl.colors(Kmax + 1L, "Zissou 1")
  for (nm in c("X", "Y", "J")) {
    H <- E[[nm]]
    graphics::plot(NA, xlim = c(0, Lmax), ylim = range(H), xlab = "self-history length l",
                   ylab = "conditional entropy (bits)",
                   main = sprintf("H(%s[t+1] | own past l, other's past k)", if (nm == "J") "X,Y" else nm), font.main = 2)
    for (k in 0:Kmax) { graphics::lines(0:Lmax, H[, k + 1L], col = cols[k + 1L], lwd = 2); graphics::points(0:Lmax, H[, k + 1L], col = cols[k + 1L], pch = 16) }
    if (nm == "X") graphics::legend("topright", legend = paste0("k = ", 0:Kmax), col = cols, lwd = 2, bty = "n", cex = 0.8)
  }
  if (!is.null(main)) graphics::mtext(main, side = 3, outer = TRUE, line = -1.5, font = 2)
  invisible(x)
}

# rcausality driver ---------------------------------------------------------
.rcausality_dfe <- function(all, grid, causes, effects, alpha, model, lag, prob_threshold, cl, cores) {
  model <- .dfe_model(model, cores = cores)
  up <- .unordered_pairs(grid, all)
  res <- .plapply(seq_len(nrow(up)), function(i)
    dfe_classify(all[[up$v1[i]]], all[[up$v2[i]]], model = model, lag = lag), cores = cores)
  rows <- list(); ins <- list(); dd <- list()
  prior <- 2 / 5
  for (i in seq_len(nrow(up))) {
    r <- res[[i]]; p <- r$scores; a <- up$v1[i]; b <- up$v2[i]
    for (dir in c("ab", "ba")) {
      p_dir <- if (dir == "ab") p[["x_to_y"]] else p[["y_to_x"]]
      post <- p_dir + p[["bidirectional"]]
      rows[[length(rows) + 1L]] <- data.frame(
        cause = if (dir == "ab") a else b, effect = if (dir == "ab") b else a, lag = lag,
        statistic = p_dir, df1 = NA_real_, df2 = NA_real_, p_value = NA_real_,
        bf10 = (post / (1 - post)) / (prior / (1 - prior)), effect_size = NA_real_, r2_gain = NA_real_,
        sign = NA_real_, n_eff = r$n, post_prob = post, stringsAsFactors = FALSE)
    }
    ins[[i]] <- data.frame(var1 = a, var2 = b, cmi = NA_real_, p_value = NA_real_, p_adj = NA_real_,
                           significant = r$best == "common_cause", stringsAsFactors = FALSE)
    dd[[i]] <- data.frame(var1 = a, var2 = b, t(p), best = r$best, stringsAsFactors = FALSE)
  }
  pairs <- do.call(rbind, rows)
  pairs$p_adj <- NA_real_; pairs$evidence <- bf_evidence(pairs$bf10)
  pairs$significant <- pairs$post_prob > prob_threshold
  pairs <- pairs[, c("cause", "effect", "lag", "statistic", "df1", "df2", "p_value",
                     "p_adj", "bf10", "post_prob", "evidence", "effect_size",
                     "r2_gain", "sign", "n_eff", "significant")]
  rownames(pairs) <- NULL
  .new_causres(pairs, series = names(all), causes = causes, effects = effects,
               method = "dfe", conditional = FALSE, alpha = alpha, bf_threshold = NA_real_,
               adjust = "none", n = nrow(all), call = cl, data = all, instant = do.call(rbind, ins),
               dc = do.call(rbind, dd), prob_threshold = prob_threshold)
}
