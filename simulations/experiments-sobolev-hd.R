
library(sphunif)
library(progressr)

# Define a progress bar
handlers(handler_progress(
  format = paste("(:spin) [:bar] :percent Iter: :current/:total Rate:",
                 ":tick_rate iter/sec ETA: :eta Elapsed: :elapsedfull"),
  clear = FALSE))

# Simulation setup
cores <- 10
M <- 1e4
n <- c(5, 50, 100, 500)
p <- c(5, 50, 100, 500) + 1

# Asymptotic standardization of statistics
hd_std_asymp <- function(stat, n, p, vk2) {

  # Only accept scalar n and p
  stopifnot(length(n) == 1, length(p) == 1)

  # Remove biases to make it a U-statistic
  nonzero_vk2 <- which(vk2 != 0)
  if (length(nonzero_vk2) == 0) stop("vk2 must contain at least one non-zero weight")
  if (any(vk2 < 0)) stop("vk2 must be non-negative")

  bk <- vk2_to_bk(vk2, p = p)
  bias <- drop(bk[nonzero_vk2] %*%
                 Gegen_polyn(theta = 0, k = nonzero_vk2, p = p))

  # # \sigma_n^{-1} for a single k0 s.t. v_{k0} != 0
  # # \sqrt{(2k_0!) / d_n^{k_0}} / 2 = \sqrt{k_0! / (2 * d_n^{k_0})}
  # # The division by 2 is done because the computed statistic is 2/n \sum_{i<j}
  # k_0 <- which(vk2 != 0)
  # inv_sigma <- sqrt(factorial(k_0) / (2 * (p - 1)^k_0))

  # General \sigma_n^{-1}
  k <- seq_along(vk2)
  inv_sigma <- 1 / sqrt(2 * sum(exp(2 * log(vk2) +
                                      d_p_k(p = p, k = k, log = TRUE))))

  # Standardize
  return(inv_sigma * (stat - bias))

}

## Run simulations under H0
{

# Weights
vk2 <- c(1, 0, 0)

# # Loop on vk2's
# vk2s <- rbind(diag(rep(1, 3)), c(1, 1, 0), c(1, 1, 1))
# for (i in seq_len(nrow(vk2s))) {
# n <- c(5, 50, 100, 500)
# p <- c(5, 50, 100, 500) + 1
# vk2 <- vk2s[i, ]

# Loop on (n, p)'s
stats <- array(dim = c(length(p), length(n), M))
for (ni in seq_along(n)) {
  for (pi in seq_along(p)) {

    # Monte Carlo
    cat("n =", n[ni], "p =", p[pi], "\n")
    with_progress(
      stats[pi, ni, ] <- unif_stat_MC(n = n[ni], p = p[pi], type = "Sobolev",
                                      M = M, r_H1 = NULL,
                                      Sobolev_vk2 = vk2,
                                      chunks = M / 10,
                                      cores = cores)$stats_MC[, 1]
    )

  }
}

# Save results
file_name_0 <- paste0("check-Sobolev-", paste0(vk2, collapse = ""))
save(stats, file = paste0(file_name_0, ".RData"))

# }

}

## Run simulations under H1
{

# Weights
vk2 <- c(1, 1, 1)
vk2_orig <- vk2

# Use dn modification according to Remark 5?
dn_modif <- TRUE

# # Loop on vk2's
# vk2s <- rbind(diag(rep(1, 3)), c(1, 1, 0), c(1, 1, 1))
# for (i in seq_len(nrow(vk2s))) { for (dn_modif in c(TRUE, FALSE)) {
#
# n <- c(5, 50, 100, 500)
# p <- c(5, 50, 100, 500) + 1
# vk2 <- vk2s[i, ]
# vk2_orig <- vk2

# Loop on (n, p)'s
k <- seq_along(vk2_orig)
stats <- array(dim = c(length(p), length(n), M))
for (ni in seq_along(n)) {
  for (pi in seq_along(p)) {

    # Use dn modification according to Remark 5?
    if (dn_modif) {

      c_k <- ifelse(k == 1, 1, 1 / (p[pi] - 1))
      vk2 <- vk2_orig * sqrt(c_k * factorial(k) * (p[pi] - 1)^(-(k - 1)))

    } else {

      vk2 <- vk2_orig

    }

    # Deviation (in asymptotic mean: Gamma * tau^2, with Gamma = 1 / sqrt(2)
    # if not blind)
    tau2 <- sqrt(2)
    # tau2 <- 1
    kappa <- sqrt(tau2) * (p[pi] - 1)^(3/4) / sqrt(n[ni])

    # Monte Carlo
    cat("n =", n[ni], "p =", p[pi], "\n")
    with_progress(
      stats[pi, ni, ] <- unif_stat_MC(n = n[ni], p = p[pi], type = "Sobolev",
                                      crit_val = NULL, M = M,
                                      r_H1 = r_alt, alt = "vMF", kappa = kappa,
                                      Sobolev_vk2 = vk2,
                                      chunks = M / 10,
                                      cores = cores)$stats_MC[, 1]
    )

  }
}

# Save results
file_name_1 <- paste0("check-vMF-Sobolev-", paste0(vk2_orig, collapse = ""),
                      ifelse(dn_modif, "-dn", ""))
save(stats, file = paste0(file_name_1, ".RData"))

# }}

}

## Plot results H0
{

# Simulation setup
n <- c(5, 50, 100, 500)
p <- c(5, 50, 100, 500) + 1

# Weights
vk2 <- c(1, 0, 0)

# # Loop on vk2's
# vk2s <- rbind(diag(rep(1, 3)), c(1, 1, 0), c(1, 1, 1))
# for (i in seq_len(nrow(vk2s))) {
# n <- c(5, 50, 100, 500)
# p <- c(5, 50, 100, 500) + 1
# vk2 <- vk2s[i, ]

# Load results
file_name_0 <- paste0("check-Sobolev-", paste0(vk2, collapse = ""))
load(file = paste0(file_name_0, ".RData"))

# Plots
# png(paste0(file_name_0, ".png"), width = 10, height = 10, res = 300,
#     units = "in", bg = "transparent")
pdf(paste0(file_name_0, ".pdf"), width = 10, height = 10)
par(mfrow = c(4, 4), mar = c(3.5, 3, 3, 0) + 0.1)
for (ni in n) {
  for (pi in p) {

    # Statistics times standardizing factor
    stats_n_p <- hd_std_asymp(stat = stats[which(p == pi), which(n == ni), ],
                              n = ni, p = pi, vk2 = vk2)
    stats_n_p <- na.omit(stats_n_p)
    if (all(is.na(stats_n_p))) stats_n_p <- rep(1, 2)

    # Histogram
    stats_n_p_trunc <- stats_n_p[abs(stats_n_p) < 10]
    hist(stats_n_p_trunc, xlab = "", ylab = "", prob = TRUE,
         breaks = seq(-10, 10, l = 75),
         main = substitute(list(n == ni, d == di),
                           list(ni = ni, di = pi - 1)),
         xlim = c(-5, 5), ylim = c(0, 0.6))
    rug(stats_n_p_trunc)
    lines(density(stats_n_p), lwd = 2)

    # Add normal fit
    mu <- 0
    sigma <- 1
    ks <- ks.test(stats_n_p, "pnorm", mean = mu, sd = sigma)
    lillie <- tryCatch(nortest::lillie.test(stats_n_p),
                       error = function(e) list(statistic = NA, p.value = NA))
    curve(dnorm(x, mean = mu, sd = sigma), add = TRUE, col = "blue", lwd = 2)
    title(sub = substitute(list("KS p-value" == p1, "L p-value" == p2),
                           list(p1 = signif(ks$p.value, 2),
                                p2 = signif(lillie$p.value, 2))),
          line = 2)

  }
}
dev.off()

# }

}

## Plot results H1
{

# Simulation setup
n <- c(5, 50, 100, 500)
p <- c(5, 50, 100, 500) + 1

# Weights
vk2 <- c(1, 0, 0)

# Use dn modification according to Remark 5?
dn_modif <- TRUE

# # Loop on vk2's
# vk2s <- rbind(diag(rep(1, 3)), c(1, 1, 0), c(1, 1, 1))
# for (i in seq_len(nrow(vk2s))) { for (dn_modif in c(TRUE, FALSE)) {
# n <- c(5, 50, 100, 500)
# p <- c(5, 50, 100, 500) + 1
# vk2 <- vk2s[i, ]
# vk2_orig <- vk2

# Load results
vk2_orig <- vk2
k <- seq_along(vk2_orig)
file_name_1 <- paste0("check-vMF-Sobolev-", paste0(vk2_orig, collapse = ""),
                      ifelse(dn_modif, "-dn", ""))
load(file = paste0(file_name_1, ".RData"))

# Plots
# png(paste0(file_name_1, ".png"), width = 10, height = 10, res = 300,
#     units = "in", bg = "transparent")
pdf(paste0(file_name_1, ".pdf"), width = 10, height = 10)
par(mfrow = c(4, 4), mar = c(3, 2.5, 1.5, 0) + 0.1)
for (ni in n) {
  for (pi in p) {

    # Use dn modification according to Remark 5?
    if (dn_modif) {

      c_k <- ifelse(k == 1, 1, 1 / (pi - 1))
      vk2 <- vk2_orig * sqrt(c_k * factorial(k) * (pi - 1)^(-(k - 1)))

    } else {

      vk2 <- vk2_orig

    }

    # Statistics times standardizing factor
    stats_n_p <- hd_std_asymp(stat = stats[which(p == pi), which(n == ni), ],
                              n = ni, p = pi, vk2 = vk2)
    stats_n_p <- na.omit(stats_n_p)
    if (all(is.na(stats_n_p))) stats_n_p <- rep(1, 2)

    # Histogram
    stats_n_p_trunc <- stats_n_p[abs(stats_n_p) < 10]
    hist(stats_n_p_trunc, xlab = "", ylab = "", prob = TRUE,
         breaks = seq(-10, 10, l = 75),
         main = "",
         xlim = c(-5, 5), ylim = c(0, 0.6))
    rug(stats_n_p_trunc)
    lines(density(stats_n_p), lwd = 2)

    # Add normal fit
    sigma <- 1
    mu_0 <- 0
    tau2 <- sqrt(2)
    # tau2 <- 1
    Gamma <- 1 / sqrt(2) * (dn_modif || all(vk2 == c(1, 0, 0)))
    mu_1 <- tau2 * Gamma
    ks <- ks.test(stats_n_p, "pnorm", mean = mu_1, sd = sigma)
    lillie <- tryCatch(nortest::lillie.test(stats_n_p),
                       error = function(e) list(statistic = NA, p.value = NA))
    curve(dnorm(x, mean = mu_1, sd = sigma), add = TRUE,
          col = "red", lwd = 2)
    curve(dnorm(x, mean = mu_0, sd = sigma), add = TRUE,
          col = "blue", lwd = 2)
    title(main = substitute(list(n == ni, d == di),
                            list(ni = ni, di = pi - 1)), line = 0)
    title(sub = substitute(list("KS" == p1, "L" == p2),
                           list(p1 = signif(ks$p.value, 2),
                                p2 = signif(lillie$p.value, 2))),
          line = 2)

  }
}
dev.off()

# }}

}
