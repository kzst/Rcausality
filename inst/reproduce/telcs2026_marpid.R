# Reproduction of
#   Telcs, A. (2026). Empirical nonnegative finite-resolution MAR-PID for
#   continuous variables. Entropy 28(6), 641. doi:10.3390/e28060641
# Sect. 6.1 (XOR under the doubling map, Table 1), Sect. 6.2 (mixed
# redundancy-synergy example, Table 2) and Sect. 6.3 (four-valued target).
library(Rcausality)
set.seed(1); n <- 4000
## 6.1: X, Y uniform, Z's binary digits = XOR of the digits; resolution k = 3.
## Table 1: synergistic atom = k, all other atoms 0.
x0 <- runif(n); y0 <- runif(n)
bits <- function(u, k) sapply(1:k, function(j) floor(u * 2^j) %% 2)
z <- as.vector(((bits(x0, 3) + bits(y0, 3)) %% 2) %*% 2^(-(1:3)))
marpid(z, x0, y0, k_target = 3, k_source = 3)      # synergy about 2.7 of 3 bits at n = 4000 (finite-sample loss)
## 6.2: A ~ Bern(1/2), latent switch B ~ Bern(q): B = 0 -> V1 = V2 = A (redundancy),
## B = 1 -> V1 = U, V2 = U xor A (XOR synergy).  Table 2 lists I(A;V1), I(A;V2),
## I(A;V1,V2) for q = 0, .25, .5, .75, 1: (1,1,1), (.456,.456,.741), (.189,.189,.656),
## (.046,.046,.697), (0,0,1).
for (q in c(0, .25, .5, .75, 1)) {
  A <- rbinom(n, 1, .5); B <- rbinom(n, 1, q); U <- rbinom(n, 1, .5)
  V1 <- ifelse(B == 0, A, U); V2 <- ifelse(B == 0, A, (U + A) %% 2)
  r <- marpid(A, V1, V2, 1, 1)$atoms
  cat(sprintf("q = %.2f: I(A;V1) = %.3f  I(A;V2) = %.3f  I(A;V1,V2) = %.3f | R = %.3f  U1 = %.3f  U2 = %.3f  S = %.3f\n",
              q, r[["redundant"]] + r[["unique_1"]], r[["redundant"]] + r[["unique_2"]], r[["total"]],
              r[["redundant"]], r[["unique_1"]], r[["unique_2"]], r[["synergy"]]))
}
## 6.3: Y = 2 C1 + C2: bit 1 of the target is unique to C1, bit 2 unique to C2
C1 <- rbinom(n, 1, .5); C2 <- rbinom(n, 1, .5)
marpid(2 * C1 + C2, C1, C2, k_target = 2, k_source = 1)$by_bit
