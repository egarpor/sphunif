
n <- 50
set.seed(1242333)
cir_0 <- unif_stat_MC(n = n, M = 1e3, type = "all", p = 2)
sph_0 <- unif_stat_MC(n = n, M = 1e3, type = "all", p = 3)

# Circular
rand_dirs_2 <- r_unif_sph(n = 20, p = 2, M = 1)[, , 1]
cir_pow <- unif_stat_MC(n = n, M = 1e3, p = 2, r_H1 = r_alt,
                        scenario = "MvMF", kappa = 0.5,
                        crit_val = cir_0$crit_val_MC,
                        Cuesta_Albertos_rand_dirs = rand_dirs_2)

# Spherical
rand_dirs_3 <- r_unif_sph(n = 20, p = 3, M = 1)[, , 1]
sph_pow <- unif_stat_MC(n = n, M = 1e3, p = 3, r_H1 = r_alt,
                        scenario = "MvMF", kappa = 0.5,
                        crit_val = sph_0$crit_val_MC,
                        Cuesta_Albertos_rand_dirs = rand_dirs_3)

test_that("Rejections for vMF", {

  expect_true(all(sph_pow$power_MC[1, ] > 0.01))
  expect_true(all(cir_pow$power_MC[1, ] > 0.01))

})
