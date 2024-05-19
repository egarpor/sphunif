library(sphunif)
library(progress)
library(progressr)

# Progress bar
handlers(handler_progress(
  format = ":spin [:bar] :percent Total: :elapsedfull End \u2248 :eta",
  clear = FALSE
))

seed <- 12345# Random seed

p <- 4# Ambient space dimension # TODO: Adapt to 3
n <- 100# Sample size
M <- 1e4# Monte Carlo samples

alpha <- 0.05# Significance level
stat_list <- c("Rayleigh", "Bingham", "Softmax", "Poisson", "Stereo")
Poisson_rho <- c(0.25, 0.5, 0.75)
Softmax_kappa <- c(1, 5, 10)
Stereo_a <- c(-1, 0, 1)

## CV tests
K <- 10# Number of folds
# Grids of parameters
Poisson_grid <- seq(0.02, 0.98, 0.02)
Softmax_grid <- c(1e-3, 1e-2, seq(0.1, 1, 0.1),
                  seq(1.2, 10, 0.2), seq(11, 30, 1),
                  seq(35, 60, 5))
Stereo_grid <- seq(-1, 1, 0.1)

# Monte Carlo M=10^5 exact-n alpha quantiles
with_progress({
  q_H0 <- unif_stat_MC(n = n, p = p, type = c(stat_list),
                       Poisson_rho = Poisson_rho,
                       Softmax_kappa = Softmax_kappa,
                       Stereo_a = Stereo_a,
                       M = 1e5, alpha = alpha)$crit_val_MC
})
# Monte Carlo M=10^5 exact-n samples to compute critical values in CV tests
with_progress({
  stats_MC <- unif_stat_MC(n = n, type = avail_sph_cv_tests,
                           p = p, M = 1e5, return_stats = TRUE,
                           stats_sorted = TRUE, Poisson_rho = Poisson_grid,
                           Softmax_kappa = Softmax_grid,
                           Stereo_a = Stereo_grid
                           )$stats_MC
})

# UAD simulation parameters
radius_list <- c(1, 10, 20, 45, 90, 135, 180)

stat_names <- colnames(q_H0)
emp_rej_antipodal <- array(dim = c(length(radius_list), length(stat_names)),
                           dimnames = list(radius_list, stat_names))
emp_rej_antipodal_cv <- array(dim = c(length(radius_list),
                                      length(avail_sph_cv_tests)),
                              dimnames = list(radius_list, avail_sph_cv_tests))

## Simulation
for (radius_deg in radius_list) {

  print(paste0("r = ", radius_deg, "º"))

  # Concentration parameter
  kappa <- pi - radius_deg * pi / 180

  # Simulation of UAD sample
  set.seed(seed)
  samp <- r_alt(n = n, p = p, M = M, alt = "UAD", kappa = kappa)

  # Statistics computation
  stat <- unif_stat(data = samp, type = stat_list,
                    Poisson_rho = Poisson_rho,
                    Softmax_kappa = Softmax_kappa,
                    Stereo_a = Stereo_a)

  # Rejection?
  for (s in stat_names) {
    emp_rej_antipodal[as.character(radius_deg),
                      s] <- mean(stat[, s] > q_H0[as.character(alpha), s])

  }

  # K-fold tests computation
  cv_reject <- array(dim = c(M, length(avail_sph_cv_tests)),
                     dimnames = list(1:M, avail_sph_cv_tests))
  for (i in 1:M) {

    # Progress
    if ((100 * i / M) %% 10 == 0) print(100 * i / M)

    cv_test <- unif_test_cv(data = samp[, , i], type = stat_list, K = K,
                            p_value = "MC", alpha = alpha,
                            stats_MC = stats_MC,
                            Poisson_rho = Poisson_grid,
                            Softmax_kappa = Softmax_grid,
                            Stereo_a = Stereo_grid,
                            seed_fold = seed)
    for (s in avail_sph_cv_tests){

      cv_reject[i, s] <- cv_test[[s]][["reject"]]

    }

  }

  for (s in avail_sph_cv_tests){

    emp_rej_antipodal_cv[as.character(radius_deg), s] <- mean(cv_reject[, s])

  }

}

## Show results
print(emp_rej_antipodal)

colnames(emp_rej_antipodal_cv) <- paste0(avail_sph_cv_tests," (", K, "-CV)")
print(emp_rej_antipodal_cv)
