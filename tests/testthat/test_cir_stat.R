
n <- 10
set.seed(123456789)
Theta_1 <- r_unif_cir(n = n, M = 1)
Theta_2 <- r_unif_cir(n = n, M = 2)

test_that("Kolmogorov-Smirnov", {

  expect_equal(drop(cir_stat_Kuiper(Theta_1, KS = TRUE)),
               as.numeric(sqrt(n) * ks.test(Theta_1, y = "punif",
                                            0, 2 * pi)$statistic))
  expect_equal(drop(cir_stat_Kuiper(Theta_2, KS = TRUE)),
               sqrt(n) * as.numeric(c(
                   ks.test(Theta_2[, 1], y = "punif",
                           0, 2 * pi)$statistic,
                   ks.test(Theta_2[, 2], y = "punif",
                           0, 2 * pi)$statistic)))

})

test_that("Cramér-von Mises", {

  expect_equal(drop(cir_stat_Watson(Theta_1, CvM = TRUE)),
               as.numeric(goftest::cvm.test(Theta_1, null = "punif",
                                            0, 2 * pi)$statistic))
  expect_equal(drop(cir_stat_Watson(Theta_2, CvM = TRUE)),
               as.numeric(c(
                 goftest::cvm.test(Theta_2[, 1], null = "punif",
                                   0, 2 * pi)$statistic,
                 goftest::cvm.test(Theta_2[, 2], null = "punif",
                                   0, 2 * pi)$statistic)))

})

test_that("Anderson-Darling", {

  expect_equal(drop(cir_stat_PAD(Theta_1, AD = TRUE)),
               as.numeric(goftest::ad.test(Theta_1, null = "punif",
                                            0, 2 * pi)$statistic))
  expect_equal(drop(cir_stat_PAD(Theta_2, AD = TRUE)),
               as.numeric(c(
                 goftest::ad.test(Theta_2[, 1], null = "punif",
                                  0, 2 * pi)$statistic,
                 goftest::ad.test(Theta_2[, 2], null = "punif",
                                  0, 2 * pi)$statistic)))

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

