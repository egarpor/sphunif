library(sphunif)
library(progress)
library(progressr)

# Progress bar
handlers(handler_progress(
  format = ":spin [:bar] :percent Total: :elapsedfull End \u2248 :eta",
  clear = FALSE
))

p <- 3# Ambient space dimension
n <- 100# Sample size
M <- 1e4# Monte Carlo samples

alpha <- 0.05# Significance level
stat_list <- c("Rayleigh", "Bingham", "Softmax", "Poisson", "Stereo")
Stereo_a <- c(-1,0,1)
Softmax_kappa <- c(1, 5, 10)
Poisson_rho <- c(0.25, 0.5, 0.75)

# Monte Carlo M=10^5 exact-n alpha quantiles
q_H0 <- unif_stat_MC(n = n, p = p, type = c(stat_list),
                     Softmax_kappa = Softmax_kappa,
                     Poisson_rho = Poisson_rho,
                     Stereo_a = Stereo_a,
                     M = 1e5, alpha = alpha)$crit_val_MC

# UAD simulation parameters
radius_list <- c(1, 10, 20, 45, 90, 135, 180)

stat_names <- colnames(q_H0)
emp_rej_antipodal <- array(dim = c(length(radius_list), length(stat_names)),
                           dimnames = list(radius_list, stat_names))

with_progress({
  for (radius_deg in radius_list) {

    print(paste0("r = ", radius_deg, "º"))

    # Concentration parameter
    kappa <- pi - radius_deg * pi / 180

    # Statistics under UAD simulation
    stat <- unif_stat_MC(n = n , p = p, M = M, type = stat_list,
                         Softmax_kappa = Softmax_kappa,
                         Poisson_rho = Poisson_rho,
                         Stereo_a = Stereo_a,
                         r_H1 = r_alt, alt = "UAD", kappa = kappa)$stats_MC

    # Rejection?
    for (s in stat_names) {
      emp_rej_antipodal[as.character(radius_deg),
                        s] <- mean(stat[, s] > q_H0[as.character(alpha), s])
    }

  }
})

print(emp_rej_antipodal)
