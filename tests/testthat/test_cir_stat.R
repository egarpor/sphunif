
n <- 10
set.seed(123456789)
Theta_1 <- r_unif_cir(n = n, M = 1)
Theta_2 <- r_unif_cir(n = n, M = 2)

test_that("PCvM vs. Watson", {

  expect_equal(cir_stat_PCvM(Theta_2), 2 * cir_stat_Watson(Theta_2))
  expect_equal(sph_stat_PCvM(Theta_to_X(Theta_1)), 2 * cir_stat_Watson(Theta_1))

})

test_that("Hodges_Ajne use_Cressie and asymp_std", {

  expect_equal(cir_stat_Hodges_Ajne(Theta_1, use_Cressie = TRUE),
               cir_stat_Hodges_Ajne(Theta_1, use_Cressie = FALSE))
  expect_equal(cir_stat_Hodges_Ajne(Theta_2, use_Cressie = TRUE),
               cir_stat_Hodges_Ajne(Theta_2, use_Cressie = FALSE))
  expect_equal(cir_stat_Hodges_Ajne(Theta_2, asymp_std = TRUE),
               (2 * cir_stat_Hodges_Ajne(Theta_2) - n) / sqrt(n))

})

test_that("cir_stat_Range max_gap", {

  expect_equal(2 * pi - cir_stat_Range(Theta_1, max_gap = FALSE),
               cir_stat_Range(Theta_1, max_gap = TRUE))

})

