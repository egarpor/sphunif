
n <- 10
set.seed(123456789)
Theta_1 <- r_unif_cir(n = n, M = 1)
Theta_2 <- r_unif_cir(n = n, M = 2)
Psi_2 <- Psi_mat(Theta_to_X(Theta_2))

test_that("Kolmogorov-Smirnov", {

  expect_equal(drop(cir_stat_Kuiper(Theta_1, KS = TRUE)),
               as.numeric(sqrt(n) * ks.test(Theta_1, y = "punif",
                                            min = 0, max = 2 * pi)$statistic))
  expect_equal(drop(cir_stat_Kuiper(Theta_2, KS = TRUE)),
               sqrt(n) * as.numeric(c(
                   ks.test(Theta_2[, 1], y = "punif",
                           min = 0, max = 2 * pi)$statistic,
                   ks.test(Theta_2[, 2], y = "punif",
                           min = 0, max = 2 * pi)$statistic)))

})

test_that("Cramér-von Mises", {

  expect_equal(drop(cir_stat_Watson(Theta_1, CvM = TRUE)),
               as.numeric(goftest::cvm.test(Theta_1, null = "punif",
                                            min = 0, max = 2 * pi)$statistic))
  expect_equal(drop(cir_stat_Watson(Theta_2, CvM = TRUE)),
               as.numeric(c(
                 goftest::cvm.test(Theta_2[, 1], null = "punif",
                                   min = 0, max = 2 * pi)$statistic,
                 goftest::cvm.test(Theta_2[, 2], null = "punif",
                                   min = 0, max = 2 * pi)$statistic)))

})

test_that("Anderson-Darling", {

  expect_equal(drop(cir_stat_PAD(Theta_1, AD = TRUE)),
               as.numeric(goftest::ad.test(Theta_1, null = "punif",
                                            min = 0, max = 2 * pi)$statistic))
  expect_equal(drop(cir_stat_PAD(Theta_2, AD = TRUE)),
               as.numeric(c(
                 goftest::ad.test(Theta_2[, 1], null = "punif",
                                  min = 0, max = 2 * pi)$statistic,
                 goftest::ad.test(Theta_2[, 2], null = "punif",
                                  min = 0, max = 2 * pi)$statistic)))

})

test_that("Stephens modifications", {
  
  expect_equal(drop(cir_stat_Kuiper(Theta_1, KS = TRUE, Stephens = TRUE)),
               drop(cir_stat_Kuiper(Theta_1, KS = TRUE, Stephens = FALSE)) *
                 (1 + 0.12 / sqrt(n) + 0.21 / n))
  expect_equal(drop(cir_stat_Kuiper(Theta_1, Stephens = TRUE)),
               drop(cir_stat_Kuiper(Theta_1, Stephens = FALSE)) *
                 (1 + 0.155 / sqrt(n) + 0.24 / n))
  expect_equal(drop(cir_stat_Watson(Theta_1, CvM = TRUE, Stephens = TRUE)),
               (drop(cir_stat_Watson(Theta_1, CvM = TRUE, Stephens = FALSE)) -
                  0.4 / n + 0.6 / n^2) * (1 + 1 / n))
  expect_equal(drop(cir_stat_Watson(Theta_1, Stephens = TRUE)),
               (drop(cir_stat_Watson(Theta_1, Stephens = FALSE)) -
                  0.1 / n + 0.1 / n^2) * (1 + 0.8 / n))
  
})

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

test_that("Watson_1976 minus", {

  expect_equal(drop(cir_stat_Watson_1976(Theta = Theta_1, minus = TRUE)),
               as.numeric(sqrt(n) * ks.test(Theta_1, y = "punif",
                                            min = 0, max = 2 * pi,
                                            alternative = "less")$statistic -
                            sqrt(n) * (mean(Theta_1 / (2 * pi)) - 1 / 2)))
  expect_equal(drop(cir_stat_Watson_1976(Theta = Theta_2, minus = TRUE)),
               apply(Theta_2, 2, function(x) {
                 as.numeric(sqrt(n) * ks.test(x, y = "punif",
                                              min = 0, max = 2 * pi,
                                              alternative = "less")$statistic -
                              sqrt(n) * (mean(x / (2 * pi)) - 1 / 2))
                 }))

})

test_that("Gine_Fn constructed from Gine_Gn and Ajne", {

  expect_equal(cir_stat_Gine_Fn(Theta_1),
               cir_stat_Gine_Gn(Theta_1) + 4 * cir_stat_Ajne(Theta_1))
  expect_equal(cir_stat_Gine_Fn(Psi_2, Psi_in_Theta = TRUE),
               cir_stat_Gine_Gn(Psi_2, Psi_in_Theta = TRUE) +
                 4 * cir_stat_Ajne(Psi_2, Psi_in_Theta = TRUE))

})

test_that("Bakshaev is a particular case of Riesz", {

  expect_equal(cir_stat_Riesz(Theta_2, s = 1), cir_stat_Bakshaev(Theta_2))
  expect_equal(cir_stat_Riesz(Psi_2, s = 1, Psi_in_Theta = TRUE),
               cir_stat_Bakshaev(Psi_2, Psi_in_Theta = TRUE))
  expect_equal(cir_stat_Bakshaev(Theta_2),
               cir_stat_Bakshaev(Psi_2, Psi_in_Theta = TRUE))

})

