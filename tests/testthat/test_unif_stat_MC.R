
# To prevent hanging with parallel computations
Sys.unsetenv("R_TESTS")

n <- 10
set.seed(1242333)
cir_0 <- unif_stat_MC(n = n, M = 5e2, type = "all", p = 2, seeds = 5)
sph_0 <- unif_stat_MC(n = n, M = 5e2, type = "all", p = 3, seeds = 5)
crit_val_bad <- sph_0$crit_val_MC
colnames(crit_val_bad) <- paste0(colnames(crit_val_bad), "_bad")

# Circular
dirs_2 <- r_unif_sph(n = 20, p = 2, M = 1)[, , 1]
cir_pow <- unif_stat_MC(n = n, M = 5e2, p = 2, r_H1 = r_alt, alt = "MvMF",
                        kappa = 0.5, crit_val = cir_0$crit_val_MC,
                        CCF09_dirs = dirs_2)

# Spherical
dirs_3 <- r_unif_sph(n = 20, p = 3, M = 1)[, , 1]
sph_pow <- unif_stat_MC(n = n, M = 5e2, p = 3, r_H1 = r_alt, alt = "MvMF",
                        kappa = 0.5, crit_val = sph_0$crit_val_MC,
                        CCF09_dirs = dirs_3)

# Statistics with vectorised parameters
cir_stats_vectorised <- c("Cressie", "Max_uncover", "Num_uncover", "Vacancy",
                          "Rayleigh", "Riesz", "Rothman", "PRt", "Poisson",
                          "Pycke_q", "Softmax", "Sobolev")
sph_stats_vectorised <- c("Riesz", "PRt", "Poisson", "Softmax", "Stereo",
                          "Sobolev")
t <- c(0.2, 0.3, 0.8)
m <- 1:3
s <- c(0, 1, 2)
kappa <- 1:3
rho <- seq(0.1, 0.9, l = 3)
vk2 <- rbind(1:3, 3:1, 3:5)

test_that("Parameter-vectorized statistics work for p = 2", {

  stats_1 <- as.matrix(
    unif_stat_MC(n = 5, p = 2, M = 2, type = cir_stats_vectorised,
                 Cressie_t = t, cov_a = t, Rayleigh_m = m, Riesz_s = s,
                 Rothman_t = t, Softmax_kappa = kappa, Poisson_rho = rho,
                 Pycke_q = t, Sobolev_vk2 = vk2, seeds = 1)$stats_MC)
  stats_2 <- lapply(1:3, function(i) as.matrix(
    unif_stat_MC(n = 5, p = 2, M = 2, type = cir_stats_vectorised,
                 Cressie_t = t[i], cov_a = t[i], Rayleigh_m = m[i],
                 Riesz_s = s[i], Rothman_t = t[i], Softmax_kappa = kappa[i],
                 Poisson_rho = rho[i], Pycke_q = t[i], Sobolev_vk2 = vk2[i, ],
                 seeds = 1)$stats_MC))
  stats_2 <- unname(do.call(cbind, stats_2))
  stats_1 <- unname(stats_1)
  stats_1 <- stats_1[, order(stats_1[1, ])]
  stats_2 <- stats_2[, order(stats_2[1, ])]
  expect_equal(stats_1, stats_2)

})

test_that("Parameter-vectorized statistics work for p = 4", {

  stats_1 <- as.matrix(
    unif_stat_MC(n = 5, p = 4, M = 2, type = sph_stats_vectorised,
                 Riesz_s = s, Rothman_t = t, Softmax_kappa = kappa,
                 Poisson_rho = rho, Stereo_a = t, Sobolev_vk2 = vk2,
                 seeds = 1)$stats_MC)
  stats_2 <- lapply(1:3, function(i) as.matrix(
    unif_stat_MC(n = 5, p = 4, M = 2, type = sph_stats_vectorised,
                 Riesz_s = s[i], Rothman_t = t[i], Softmax_kappa = kappa[i],
                 Poisson_rho = rho[i], Stereo_a = t[i], Sobolev_vk2 = vk2[i, ],
                 seeds = 1)$stats_MC))
  stats_2 <- unname(do.call(cbind, stats_2))
  stats_1 <- unname(stats_1)
  stats_1 <- stats_1[, order(stats_1[1, ])]
  stats_2 <- stats_2[, order(stats_2[1, ])]
  expect_equal(stats_1, stats_2)

})

test_that("Rejections for MvMF", {

  expect_true(all(sph_pow$power_MC[1, ] > 0.01))
  expect_true(all(cir_pow$power_MC[1, ] > 0.01))

})

test_that("Edge cases", {

  expect_error(unif_stat_MC(n = n, M = 1e2, type = "all", p = 1))
  expect_error(unif_stat_MC(n = n, M = 1e2, type = "all", p = 3,
                            crit_val = sph_0$crit_val_MC[, 1:3]))
  expect_error(unif_stat_MC(n = n, M = 1e2, type = "all", p = 3,
                            crit_val = crit_val_bad[, 1:3]))
  suppressWarnings(expect_warning(unif_stat_MC(n = n, M = 1e2, type = "all",
                                               p = 5, seeds = 1:3,
                                               chunks = 2)))
  expect_equal(unif_stat_MC(n = n, M = 1e2, type = c("PAD", "Ajne"), p = 3,
                            seeds = 5)$stats$PAD,
               unif_stat_MC(n = n, M = 1e2, type = "PAD", p = 3, seeds = 5,
                            crit_val = sph_0$crit_val_MC)$stats$PAD)

})

test_that("Several options", {

  expect_equal(unif_stat_MC(n = n, M = 10, type = "all", p = 2,
                            return_stats = FALSE)$stats, NA)
  suppressWarnings(expect_equal(unif_stat_MC(n = n, M = 10, type = "all", p = 3,
                                             return_stats = FALSE)$stats, NA))
  suppressWarnings(
    expect_equal(unname(as.matrix(unif_stat_MC(n = n, M = 10, type = "all",
                                               p = 5, stats_sorted = TRUE,
                                               seeds = 1)$stats)),
                 sort_each_col(as.matrix(unif_stat_MC(n = n, M = 10,
                                                      type = "all", p = 5,
                                                      seeds = 1)$stats)))
  )
  set.seed(1)
  CCF09_dirs <- r_unif_sph(n = 50, p = 3, M = 1)[, , 1]
  expect_equal(unif_stat_MC(n = n, M = 5, type = "CCF09", p = 3,
                            stats_sorted = TRUE, seeds = 1, chunks = 1)$CCF09,
               unif_stat_MC(n = n, M = 5, type = "CCF09", p = 3,
                            stats_sorted = TRUE, seeds = 1, chunks = 1,
                            CCF09_dirs = CCF09_dirs)$CCF09)

})

test_that("Progress bars", {

  skip_on_cran()
  o1_silent <- capture.output(s <- unif_stat_MC(n = n, M = 10, type = "all",
                                                p = 2, cores = 1, chunks = 10))
  o2_silent <- capture.output(s <- unif_stat_MC(n = n, M = 10, type = "all",
                                                p = 2, cores = 2, chunks = 10))
  expect_equal(length(o1_silent), 0)
  expect_equal(length(o2_silent), 0)

})

test_that("Different results with seeds = NULL (default)", {

  skip_on_cran()
  expect_false(all(unif_stat_MC(n = n, M = 10, type = "Rayleigh",
                                p = 2, chunks = 1, cores = 1,
                                seeds = NULL)$stats_MC$Rayleigh ==
                     unif_stat_MC(n = n, M = 10, type = "Rayleigh",
                                  p = 2, chunks = 1, cores = 1,
                                  seeds = NULL)$stats_MC$Rayleigh))
  expect_false(all(unif_stat_MC(n = n, M = 10, type = "Rayleigh",
                                p = 2, chunks = 2, cores = 1,
                                seeds = NULL)$stats_MC$Rayleigh ==
                     unif_stat_MC(n = n, M = 10, type = "Rayleigh",
                                  p = 2, chunks = 2, cores = 1,
                                  seeds = NULL)$stats_MC$Rayleigh))
  expect_false(all(unif_stat_MC(n = n, M = 10, type = "Rayleigh",
                                p = 2, chunks = 1, cores = 2,
                                seeds = NULL)$stats_MC$Rayleigh ==
                     unif_stat_MC(n = n, M = 10, type = "Rayleigh",
                                  p = 2, chunks = 1, cores = 2,
                                  seeds = NULL)$stats_MC$Rayleigh))
  expect_false(all(unif_stat_MC(n = n, M = 10, type = "Rayleigh",
                                p = 2, chunks = 2, cores = 2,
                                seeds = NULL)$stats_MC$Rayleigh ==
                     unif_stat_MC(n = n, M = 10, type = "Rayleigh",
                                  p = 2, chunks = 2, cores = 2,
                                  seeds = NULL)$stats_MC$Rayleigh))

})

test_that("Same results with cores = 1 and cores = 2 and fixed seeds", {

  skip_on_cran()
  expect_true(max(abs(
    unif_stat_MC(n = n, M = 10, type = "all", p = 2, chunks = 10,
                 cores = 1, seeds = 1:10)$stats_MC -
      unif_stat_MC(n = n, M = 10, type = "all", p = 2, chunks = 10,
                   cores = 2, seeds = 1:10)$stats_MC)) < 1e-10)

})

test_that("Parallelization is faster", {

  skip_on_ci()
  skip_on_cran()
  t1 <- system.time(unif_stat_MC(n = 100, M = 1e4, type = "all", p = 2,
                                 chunks = 10, cores = 1))[3]
  t2 <- system.time(unif_stat_MC(n = 100, M = 1e4, type = "all", p = 2,
                                 chunks = 10, cores = 2))[3]
  expect_gt(t1, t2)

})
