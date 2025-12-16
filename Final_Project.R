set.seed(123)

############################
## 1. Observed data tables
############################

tec_77 <- rbind(
  c(66, 87, 25, 22,  4),
  c(13, 14, 15,  9,  4),
  c( 0,  4,  4,  9,  1),
  c( 0,  0,  4,  3,  1),
  c( 0,  0,  0,  1,  1),
  c( 0,  0,  0,  0,  0)
)

tec_80 <- rbind(
  c(44, 62, 47, 38,  9),
  c(10, 13,  8, 11,  5),
  c( 0,  9,  2,  7,  3),
  c( 0,  0,  3,  5,  1),
  c( 0,  0,  0,  1,  0),
  c( 0,  0,  0,  0,  1)
)

sea_B <- rbind(
  c( 9, 12, 18,  9,  4),
  c( 1,  6,  6,  4,  3),
  c( 0,  2,  3,  4,  0),
  c( 0,  0,  1,  3,  2),
  c( 0,  0,  0,  0,  0),
  c( 0,  0,  0,  0,  0)
)

sea_H1N1 <- rbind(
  c(15, 12,  4),
  c(11, 17,  4),
  c( 0, 21,  4),
  c( 0,  0,  5)
)

############################
## 2. Longini–Koopman w_{j,s}
############################

wjs_matrix <- function(qh, qc, max_s) {
  max_j <- max_s
  W <- matrix(0, nrow = max_j + 1, ncol = max_s)
  alpha <- numeric(max_j + 1)
  alpha[1] <- 1
  
  for (s in 1:max_s) {
    # j = 0
    W[1, s] <- qc^s
    
    # 0 < j < s
    if (s >= 2) {
      for (j in 1:(s - 1)) {
        W[j + 1, s] <- choose(s, j) * alpha[j + 1] *
          (qc * qh^j)^(s - j)
      }
    }
    
    # j = s
    W[s + 1, s] <- 1 - sum(W[1:s, s])
    alpha[s + 1] <- W[s + 1, s]
  }
  
  # tidy columns
  W[W < 0] <- 0
  cs <- colSums(W)
  for (s in 1:max_s) if (cs[s] > 0) W[, s] <- W[, s] / cs[s]
  W
}

############################
## 3. Simulate outbreak table
############################

simulate_counts <- function(qh, qc, obs_counts) {
  max_s <- ncol(obs_counts)
  max_j <- nrow(obs_counts) - 1
  W <- wjs_matrix(qh, qc, max_s)
  stopifnot(nrow(W) == max_j + 1)
  
  N_s <- colSums(obs_counts)
  D_star <- matrix(0, nrow = max_j + 1, ncol = max_s)
  
  for (s in 1:max_s) {
    if (N_s[s] == 0) next
    probs <- W[1:(s + 1), s]
    D_star[1:(s + 1), s] <- rmultinom(1, N_s[s], prob = probs)
  }
  D_star
}

############################
## 4. Distance between data & sim
############################

distance_two_outbreaks <- function(theta, obs1, obs2) {
  qh1 <- theta[1]; qc1 <- theta[2]
  qh2 <- theta[3]; qc2 <- theta[4]
  
  D1 <- simulate_counts(qh1, qc1, obs1)
  D2 <- simulate_counts(qh2, qc2, obs2)
  
  d1 <- sqrt(sum((obs1 - D1)^2))
  d2 <- sqrt(sum((obs2 - D2)^2))
  0.5 * (d1 + d2)
}

############################
## 5. ABC rejection (4 params)
############################

abc_rejection_4param <- function(obs1, obs2,
                                 N = 50000,
                                 accept_frac = 0.01) {
  # prior: U(0,1)^4
  theta <- matrix(runif(N * 4), ncol = 4)
  colnames(theta) <- c("qh1", "qc1", "qh2", "qc2")
  
  # simulate & compute distances
  dists <- apply(theta, 1, distance_two_outbreaks,
                 obs1 = obs1, obs2 = obs2)
  
  # choose epsilon as quantile
  eps <- quantile(dists, accept_frac)
  keep <- dists <= eps
  
  list(theta = theta[keep, , drop = FALSE],
       dist  = dists[keep],
       eps   = eps)
}

############################
## 6. Run ABC + plots
############################

set.seed(2025)

res_tec <- abc_rejection_4param(tec_77, tec_80,
                                N = 900000,
                                accept_frac = 0.0005)

res_sea <- abc_rejection_4param(sea_B, sea_H1N1,
                                N = 900000,
                                accept_frac = 0.0005)

par(mfrow = c(1, 2))

plot(res_tec$theta[, "qh1"], res_tec$theta[, "qc1"],
     xlab = expression(q[h]), ylab = expression(q[c]),
     xlim = c(0, 1), ylim = c(0, 1),
     pch = 16, col = "red",
     main = "(a) Tecumseh: H3N2")
points(res_tec$theta[, "qh2"], res_tec$theta[, "qc2"],
       pch = 1, col = "blue")
legend("bottomleft",
       legend = c("1977–78 (qh1,qc1)", "1980–81 (qh2,qc2)"),
       col = c("red", "blue"), pch = c(16, 1), bty = "n")

plot(res_sea$theta[, "qh1"], res_sea$theta[, "qc1"],
     xlab = expression(q[h]), ylab = expression(q[c]),
     xlim = c(0, 1), ylim = c(0, 1),
     pch = 16, col = "red",
     main = "(c) Seattle: B vs H1N1")
points(res_sea$theta[, "qh2"], res_sea$theta[, "qc2"],
       pch = 1, col = "blue")
legend("bottomleft",
       legend = c("Influenza B (qh1,qc1)", "H1N1 (qh2,qc2)"),
       col = c("red", "blue"), pch = c(16, 1), bty = "n")

cat("Tecumseh epsilon =", res_tec$eps, "\n")
cat("Seattle epsilon   =", res_sea$eps, "\n")



