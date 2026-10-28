# Reproduction of
#   Jakovac, A., Fabo, D., Eross, L., Stippinger, M., Benko, Z., Somogyvari, Z.
#   and Telcs, A. Causality beyond the Reichenbach principle: the
#   degrees-of-freedom method (manuscript).
# Sect. 3, Figures 2-3: random forest trained on coupled binary Markov chains
# (five classes, lengths 200-5000, transition probabilities uniform on
# [0.001, 0.49] U [0.51, 0.99], driver state selects the transition table),
# tested on fresh chains; the paper reports 92.4 % accuracy at length 2000
# (1000 training series per class, default scikit-learn random forest).
library(Rcausality)
cores <- 2L
m <- dfe_train(n_per_class = 600, lengths = c(200, 500, 1000, 2000, 5000), seed = 1, cores = cores)
m
classes <- c("x_to_y", "y_to_x", "bidirectional", "common_cause", "independent")
for (L in c(500, 1000, 2000, 5000)) {
  truth <- rep(classes, each = 40)
  pred <- sapply(seq_along(truth), function(i) {
    d <- sim_markov_pair(L, truth[i], seed = 1000 + i); dfe_classify(d$x, d$y, model = m)$best })
  cat(sprintf("\nlength %d: accuracy %.3f\n", L, mean(truth == pred)))
  print(table(truth = factor(truth, classes), predicted = factor(pred, classes)))
}
# entropy curves of one x -> y chain (Fig. 1 style)
d <- sim_markov_pair(2000, "y_to_x", df_x = 2, df_y = 1, seed = 3)
plot(dfe_classify(d$x, d$y, model = m))
