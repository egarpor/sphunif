n <- 100
K_list <- c(1, 2, 5, 10, 23, 50)
rho_list <- c(0.1, 0.25, 0.75)
kappa_list <- c(0.1, 1, 10)
a_list <- c(-1, 0, 1)

test_that("K-fold partition selects every and just once each observation", {
  
  for (K in K_list){
    
    folds <- k_fold_split(n, K)
    expect_equal(n, length(unique(unlist(folds))))
    
  }
  
})

test_that("K-fold partitions (min, max, and mean subsizes) are well balanced", {
  
  for (K in K_list){
    
    folds <- k_fold_split(n, K)
    theo_length <- n / K
    fold_lengths <- unlist(lapply(folds, length))
    
    expect_equal(mean(fold_lengths), theo_length, tolerance = 0)
    expect_lte(max(fold_lengths), as.integer(theo_length) + 1)
    expect_equal(min(fold_lengths), as.integer(theo_length), tolerance = 0)
    
  }
  
})

test_that("K-fold partition when K > n/2 must throw an error", {
  
  expect_error(k_fold_split(n = n, K = n + 2))
  expect_error(k_fold_split(n = n, K = n / 2 + 1))
  
})

test_that("Variance under null hypothesis coincides with MC variance (p > 3)", {
  
  for (p in c(4, 10)){
    
    # Monte Carlo variance
    MC_var <- apply(unif_stat_MC(n = n, p = p,
                                 type = c("Poisson", "Softmax", "Stereo"),
                                 M = 1e4,
                                 Poisson_rho = rho_list,
                                 Softmax_kappa = kappa_list,
                                 Stereo_a = a_list)$stats_MC, 2, var)
    
    # Theoretical variance
    Poisson_var <- null_var(n = n, p = p, type = "Poisson",
                            lambda_grid = rho_list)
    Softmax_var <- null_var(n = n, p = p, type = "Softmax",
                            lambda_grid = kappa_list)
    Stereo_var <- null_var(n = n, p = p, type = "Stereo",
                           lambda_grid = a_list)
    
    expect_equal(Poisson_var, MC_var[seq_along(rho_list)],
                 tolerance = 4e-2, ignore_attr = TRUE)
    expect_equal(Softmax_var,
                 MC_var[(length(rho_list) + 1):(
                   length(rho_list) + length(kappa_list))],
                 tolerance = 4e-2, ignore_attr = TRUE)
    expect_equal(Stereo_var,
                 MC_var[(length(rho_list) + length(kappa_list) + 1):(
                   length(rho_list) + length(kappa_list) + length(a_list))],
                 tolerance = 4e-2, ignore_attr = TRUE)
    
  }
  
})

test_that("Variance under null hypothesis equals MC variance (p <= 3)", {
  
  for (p in c(2, 3)){
    
    # Monte Carlo variance
    MC_var <- apply(unif_stat_MC(n = n, p = p,
                                 type = c("Poisson", "Softmax", "Stereo"),
                                 M = 1e4,
                                 Poisson_rho = rho_list,
                                 Softmax_kappa = kappa_list,
                                 Stereo_a = a_list)$stats_MC, 2, var)
    
    # Theoretical variance
    Poisson_var <- null_var(n = n, p = p, type = "Poisson",
                            lambda_grid = rho_list)
    Softmax_var <- null_var(n = n, p = p, type = "Softmax",
                            lambda_grid = kappa_list)
    
    expect_equal(Poisson_var,
                 MC_var[seq_along(rho_list)],
                 tolerance = 4e-2, ignore_attr = TRUE)
    expect_equal(Softmax_var,
                 MC_var[(length(rho_list) + 1):(
                   length(rho_list) + length(kappa_list))],
                 tolerance = 4e-2, ignore_attr = TRUE)
    expect_error(null_var(n = n, p = p, type = "Stereo",
                          lambda_grid = a_list))
    
  }
  
})

# Hyperspherical
for (p in c(4, 10)){
  
  set.seed(12345)
  samp <- r_unif_sph(n = n, p = p)
  
  # Hide warnings because of precision losses in besselI for high K_max
  t_asymp <- suppressWarnings(unif_test_cv(data = samp, type = "all", K = 4,
                                           p_value = "asymp", seed_fold = 123, K_max = 1000))
  t_MC <- unif_test_cv(data = samp, type = "all", K = 4,
                       p_value = "MC", M = 1e3, seed_fold = 123)
  stats_MC <- unif_stat_MC(n = n - round(n / 4),
                           type = c("Softmax", "Poisson", "Stereo"),
                           p = p, M = 1e3, r_H1 = NULL, crit_val = NULL,
                           return_stats = TRUE, stats_sorted = TRUE,
                           Poisson_rho = seq(0.1, 0.9, 0.1),
                           Softmax_kappa = seq(0.1, 20, 1),
                           Stereo_a = seq(-1, 1, 0.25))$stats_MC
  t_stats_MC <- unif_test_cv(data = samp, type = "all", K = 4,
                             p_value = "MC", stats_MC = stats_MC,
                             seed_fold = 123)
  
  
  test_that("Asymptotic matches MC p-value", {
    
    expect_equal(t_asymp$Softmax$p.value, t_MC$Softmax$p.value,
                 tolerance = 0.2)
    expect_equal(t_asymp$Poisson$p.value, t_MC$Poisson$p.value,
                 tolerance = 0.2)
    expect_equal(t_asymp$Stereo$p.value, t_MC$Stereo$p.value,
                 tolerance = 0.2)
    
  })
  
  test_that("MC p-value equals p-value with stats_MC", {
    
    expect_equal(t_MC$Softmax$p.value, t_stats_MC$Softmax$p.value,
                 tolerance = 0.1)
    expect_equal(t_MC$Poisson$p.value, t_stats_MC$Poisson$p.value,
                 tolerance = 0.1)
    expect_equal(t_MC$Stereo$p.value, t_stats_MC$Stereo$p.value,
                 tolerance = 0.15)
    
  })
  
}

# Circular and spherical
for (p in c(2, 3)){
  
  set.seed(12345)
  samp <- r_unif_sph(n = n, p = p)
  
  t_asymp <- unif_test_cv(data = samp, type = c("Softmax", "Poisson"), K = 4,
                          p_value = "asymp", seed_fold = 123, K_max = 100)
  t_MC <- unif_test_cv(data = samp, type = c("Softmax", "Poisson"), K = 4,
                       p_value = "MC", M = 1e3, seed_fold = 123)
  stats_MC <- unif_stat_MC(n = n - round(n / 4), type = c("Softmax", "Poisson"),
                           p = p, M = 1e3, r_H1 = NULL, crit_val = NULL,
                           return_stats = TRUE, stats_sorted = TRUE,
                           Poisson_rho = seq(0.1, 0.9, 0.1),
                           Softmax_kappa = seq(0.1, 20, 1),
                           Stereo_a = seq(-1, 1, 0.25))$stats_MC
  t_stats_MC <- unif_test_cv(data = samp, type = c("Softmax", "Poisson"), K = 4,
                             p_value = "MC", stats_MC = stats_MC,
                             M = 1e3, seed_fold = 123)
  
  test_that("Asymptotic matches MC p-value", {
    
    expect_equal(t_asymp$Softmax$p.value, t_MC$Softmax$p.value,
                 tolerance = 0.2)
    expect_equal(t_asymp$Poisson$p.value, t_MC$Poisson$p.value,
                 tolerance = 0.2)
    
  })
  
  test_that("MC p-value equals p-value with stats_MC", {
    
    expect_equal(t_MC$Softmax$p.value, t_stats_MC$Softmax$p.value,
                 tolerance = 0.1)
    expect_equal(t_MC$Poisson$p.value, t_stats_MC$Poisson$p.value,
                 tolerance = 0.1)
    
  })
  
}

test_that("If null_variance is given, unif_test_cv checks every statistic
          and every parameter is given", {
            
            Poisson_grid <- c(0.1, 0.5, 0.8)
            n_v <- list("Poisson" = null_var(n = 10, p = 3, type = "Poisson",
                                             lambda_grid = Poisson_grid))
            expect_error(unif_test_cv(data = r_unif_sph(n = 10, p = 3), K = 3,
                                      type = c("Poisson", "Softmax"),
                                      null_variance = n_v, Poisson_rho = Poisson_grid))
            expect_error(unif_test_cv(data = r_unif_sph(n = 10, p = 3), K = 3,
                                      type = c("Poisson"), null_variance = n_v,
                                      Poisson_rho = Poisson_grid[1:2]))
            expect_no_error(unif_test_cv(data = r_unif_sph(n = 10, p = 3), K = 3,
                                         type = c("Poisson"), null_variance = n_v,
                                         Poisson_rho = Poisson_grid))
            
          })