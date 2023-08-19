
n <- 10
set.seed(123456789)
Th1 <- r_unif_cir(n = n, M = 1)
Th2 <- r_unif_cir(n = n, M = 2)
Psi1 <- Psi_mat(Theta_to_X(Th1))
Psi2 <- Psi_mat(Theta_to_X(Th2))
Th1_rep <- rbind(Th1, Th1)
Th2_rep <- rbind(Th2, Th2)
Psi1_rep <- Psi_mat(Theta_to_X(Th1_rep))
Psi2_rep <- Psi_mat(Theta_to_X(Th2_rep))
h <- function(th) 0.5 * (th^2 / (4 * pi^2) - th / (2 * pi) + 1 / 6)

test_that("Kolmogorov-Smirnov", {

  expect_equal(drop(cir_stat_Kuiper(Th1, KS = TRUE)),
               as.numeric(sqrt(n) * ks.test(Th1, y = "punif",
                                            min = 0, max = 2 * pi)$statistic))
  expect_equal(drop(cir_stat_Kuiper(Th2, KS = TRUE)),
               sqrt(n) * as.numeric(c(
                   ks.test(Th2[, 1], y = "punif",
                           min = 0, max = 2 * pi)$statistic,
                   ks.test(Th2[, 2], y = "punif",
                           min = 0, max = 2 * pi)$statistic)))

})

test_that("Cramér-von Mises", {

  expect_equal(drop(cir_stat_Watson(Th1, CvM = TRUE)),
               as.numeric(goftest::cvm.test(Th1, null = "punif",
                                            min = 0, max = 2 * pi)$statistic))
  expect_equal(drop(cir_stat_Watson(Th2, CvM = TRUE)),
               as.numeric(c(
                 goftest::cvm.test(Th2[, 1], null = "punif",
                                   min = 0, max = 2 * pi)$statistic,
                 goftest::cvm.test(Th2[, 2], null = "punif",
                                   min = 0, max = 2 * pi)$statistic)))

})

test_that("Anderson-Darling", {

  expect_equal(drop(cir_stat_PAD(Th1, AD = TRUE)),
               as.numeric(goftest::ad.test(Th1, null = "punif",
                                            min = 0, max = 2 * pi)$statistic))
  expect_equal(drop(cir_stat_PAD(Th2, AD = TRUE)),
               as.numeric(c(
                 goftest::ad.test(Th2[, 1], null = "punif",
                                  min = 0, max = 2 * pi)$statistic,
                 goftest::ad.test(Th2[, 2], null = "punif",
                                  min = 0, max = 2 * pi)$statistic)))

})

test_that("Stephens modifications", {

  expect_equal(drop(cir_stat_Kuiper(Th1, KS = TRUE, Stephens = TRUE)),
               drop(cir_stat_Kuiper(Th1, KS = TRUE, Stephens = FALSE)) *
                 (1 + 0.12 / sqrt(n) + 0.21 / n))
  expect_equal(drop(cir_stat_Kuiper(Th1, Stephens = TRUE)),
               drop(cir_stat_Kuiper(Th1, Stephens = FALSE)) *
                 (1 + 0.155 / sqrt(n) + 0.24 / n))
  expect_equal(drop(cir_stat_Watson(Th1, CvM = TRUE, Stephens = TRUE)),
               (drop(cir_stat_Watson(Th1, CvM = TRUE, Stephens = FALSE)) -
                  0.4 / n + 0.6 / n^2) * (1 + 1 / n))
  expect_equal(drop(cir_stat_Watson(Th1, Stephens = TRUE)),
               (drop(cir_stat_Watson(Th1, Stephens = FALSE)) -
                  0.1 / n + 0.1 / n^2) * (1 + 0.8 / n))

})

test_that("PCvM vs. Watson", {

  expect_equal(cir_stat_PCvM(Th2), 2 * cir_stat_Watson(Th2))
  expect_equal(sph_stat_PCvM(Theta_to_X(Th1)), 2 * cir_stat_Watson(Th1))

})


test_that("Watson form in MJ (2000, page 111)", {

  expect_equal(drop(cir_stat_Watson(Th1)), 2 * sum(h(Psi1)) / n + 1 / 12)

})

test_that("Hodges_Ajne use_Cressie and asymp_std", {

  expect_equal(cir_stat_Hodges_Ajne(Th1, use_Cressie = TRUE),
               cir_stat_Hodges_Ajne(Th1, use_Cressie = FALSE))
  expect_equal(cir_stat_Hodges_Ajne(Th2, use_Cressie = TRUE),
               cir_stat_Hodges_Ajne(Th2, use_Cressie = FALSE))
  expect_equal(cir_stat_Hodges_Ajne(Th2, asymp_std = TRUE),
               (2 * cir_stat_Hodges_Ajne(Th2) - n) / sqrt(n))

})

test_that("cir_stat_Range max_gap", {

  expect_equal(2 * pi - cir_stat_Range(Th1, max_gap = FALSE),
               cir_stat_Range(Th1, max_gap = TRUE))

})

test_that("Watson_1976 minus", {

  expect_equal(drop(cir_stat_Watson_1976(Theta = Th1, minus = TRUE)),
               as.numeric(sqrt(n) * ks.test(Th1, y = "punif",
                                            min = 0, max = 2 * pi,
                                            alternative = "less")$statistic -
                            sqrt(n) * (mean(Th1 / (2 * pi)) - 1 / 2)))
  expect_equal(drop(cir_stat_Watson_1976(Theta = Th2, minus = TRUE)),
               apply(Th2, 2, function(x) {
                 as.numeric(sqrt(n) * ks.test(x, y = "punif",
                                              min = 0, max = 2 * pi,
                                              alternative = "less")$statistic -
                              sqrt(n) * (mean(x / (2 * pi)) - 1 / 2))
                 }))

})

test_that("Gine_Fn constructed from Gine_Gn and Ajne", {

  expect_equal(cir_stat_Gine_Fn(Th1),
               cir_stat_Gine_Gn(Th1) + 4 * cir_stat_Ajne(Th1))
  expect_equal(cir_stat_Gine_Fn(Psi2, Psi_in_Theta = TRUE),
               cir_stat_Gine_Gn(Psi2, Psi_in_Theta = TRUE) +
                 4 * cir_stat_Ajne(Psi2, Psi_in_Theta = TRUE))

})

test_that("Bakshaev is a particular case of Riesz", {

  expect_equal(cir_stat_Riesz(Th2, s = 1), cir_stat_Bakshaev(Th2))
  expect_equal(cir_stat_Riesz(Psi2, s = 1, Psi_in_Theta = TRUE),
               cir_stat_Bakshaev(Psi2, Psi_in_Theta = TRUE))
  expect_equal(cir_stat_Bakshaev(Th2),
               cir_stat_Bakshaev(Psi2, Psi_in_Theta = TRUE))

})

test_that("PAD with data repetitions is computable", {

  expect_warning(expect_true(is.finite(cir_stat_PAD(Th1_rep))))
  expect_warning(expect_warning(
    expect_true(all(is.finite(cir_stat_PAD(Th2_rep))))))

})
