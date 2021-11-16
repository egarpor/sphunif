
k <- 0:5
k1 <- 1:5
psi_PCvM_2 <- function(th) psi_Pn(theta = th, q = 1, type = "PCvM")
psi_PCvM_3 <- function(th) psi_Pn(theta = th, q = 2, type = "PCvM")
psi_PCvM_4 <- function(th) psi_Pn(theta = th, q = 3, type = "PCvM")
psi_PCvM_11 <- function(th) psi_Pn(theta = th, q = 10, type = "PCvM")
psi_PAD_2 <- function(th) psi_Pn(theta = th, q = 1, type = "PAD")
psi_PAD_3 <- function(th) psi_Pn(theta = th, q = 2, type = "PAD")
psi_PAD_4 <- function(th) psi_Pn(theta = th, q = 3, type = "PAD")
psi_PAD_11 <- function(th) psi_Pn(theta = th, q = 10, type = "PAD")
psi_PRt_2 <- function(th) psi_Pn(theta = th, q = 1, type = "PRt")
psi_PRt_3 <- function(th) psi_Pn(theta = th, q = 2, type = "PRt")
psi_PRt_4 <- function(th) psi_Pn(theta = th, q = 3, type = "PRt")
psi_PRt_11 <- function(th) psi_Pn(theta = th, q = 10, type = "PRt")

test_that("Gegen_coefs (Gauss = TRUE) vs. *_Pn (Gauss = FALSE) for PCvM", {

  expect_equal(Gegen_coefs(psi = psi_PCvM_2, k = k, p = 2, Gauss = TRUE),
               Gegen_coefs_Pn(k = k, p = 2, type = "PCvM", Gauss = FALSE),
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_PCvM_3, k = k, p = 3, Gauss = TRUE),
               Gegen_coefs_Pn(k = k, p = 3, type = "PCvM", Gauss = FALSE),
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_PCvM_4, k = k, p = 4, Gauss = TRUE),
               Gegen_coefs_Pn(k = k, p = 4, type = "PCvM", Gauss = FALSE),
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_PCvM_11, k = k, p = 11, Gauss = TRUE),
               Gegen_coefs_Pn(k = k, p = 11, type = "PCvM", Gauss = FALSE),
               tolerance = 1e-6)

})

test_that("Gegen_coefs (Gauss = FALSE) vs. *_Pn (Gauss = TRUE) for PCvM", {

  expect_equal(Gegen_coefs(psi = psi_PCvM_2, k = k, p = 2, Gauss = FALSE),
               Gegen_coefs_Pn(k = k, p = 2, type = "PCvM", Gauss = TRUE),
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_PCvM_3, k = k, p = 3, Gauss = FALSE),
               Gegen_coefs_Pn(k = k, p = 3, type = "PCvM", Gauss = TRUE),
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_PCvM_4, k = k, p = 4, Gauss = FALSE),
               Gegen_coefs_Pn(k = k, p = 4, type = "PCvM", Gauss = TRUE),
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_PCvM_11, k = k, p = 11, Gauss = FALSE),
               Gegen_coefs_Pn(k = k, p = 11, type = "PCvM", Gauss = TRUE),
               tolerance = 1e-6)

})

test_that("Gegen_coefs (Gauss = TRUE) vs. *_Pn (Gauss = FALSE) for PAD", {

  expect_equal(Gegen_coefs(psi = psi_PAD_2, k = k, p = 2, Gauss = TRUE),
               Gegen_coefs_Pn(k = k, p = 2, type = "PAD", Gauss = FALSE),
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_PAD_3, k = k, p = 3, Gauss = TRUE),
               Gegen_coefs_Pn(k = k, p = 3, type = "PAD", Gauss = FALSE),
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_PAD_4, k = k, p = 4, Gauss = TRUE),
               Gegen_coefs_Pn(k = k, p = 4, type = "PAD", Gauss = FALSE),
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_PAD_11, k = k, p = 11, Gauss = TRUE),
               Gegen_coefs_Pn(k = k, p = 11, type = "PAD", Gauss = FALSE),
               tolerance = 1e-6)

})

test_that("Gegen_coefs (Gauss = FALSE) vs. *_Pn (Gauss = TRUE) for PAD", {

  expect_equal(Gegen_coefs(psi = psi_PAD_2, k = k, p = 2, Gauss = FALSE),
               Gegen_coefs_Pn(k = k, p = 2, type = "PAD", Gauss = TRUE),
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_PAD_3, k = k, p = 3, Gauss = FALSE),
               Gegen_coefs_Pn(k = k, p = 3, type = "PAD", Gauss = TRUE),
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_PAD_4, k = k, p = 4, Gauss = FALSE),
               Gegen_coefs_Pn(k = k, p = 4, type = "PAD", Gauss = TRUE),
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_PAD_11, k = k, p = 11, Gauss = FALSE),
               Gegen_coefs_Pn(k = k, p = 11, type = "PAD", Gauss = TRUE),
               tolerance = 1e-6)

})

test_that("Gegen_coefs (Gauss = TRUE) vs. *_Pn (Gauss = FALSE) for PRt", {

  expect_equal(Gegen_coefs(psi = psi_PRt_2, k = k, p = 2, Gauss = TRUE),
               Gegen_coefs_Pn(k = k, p = 2, type = "PRt", Gauss = FALSE),
               tolerance = 1e-5)
  expect_equal(Gegen_coefs(psi = psi_PRt_3, k = k, p = 3, Gauss = TRUE),
               Gegen_coefs_Pn(k = k, p = 3, type = "PRt", Gauss = FALSE),
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_PRt_4, k = k, p = 4, Gauss = TRUE),
               Gegen_coefs_Pn(k = k, p = 4, type = "PRt", Gauss = FALSE),
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_PRt_11, k = k, p = 11, Gauss = TRUE),
               Gegen_coefs_Pn(k = k, p = 11, type = "PRt", Gauss = FALSE),
               tolerance = 1e-6)

})

test_that("Gegen_coefs (Gauss = FALSE) vs. *_Pn (Gauss = TRUE) for PRt", {

  expect_equal(Gegen_coefs(psi = psi_PRt_2, k = k, p = 2, Gauss = FALSE),
               Gegen_coefs_Pn(k = k, p = 2, type = "PRt", Gauss = TRUE),
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_PRt_3, k = k, p = 3, Gauss = FALSE),
               Gegen_coefs_Pn(k = k, p = 3, type = "PRt", Gauss = TRUE),
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_PRt_4, k = k, p = 4, Gauss = FALSE),
               Gegen_coefs_Pn(k = k, p = 4, type = "PRt", Gauss = TRUE),
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_PRt_11, k = k, p = 11, Gauss = FALSE),
               Gegen_coefs_Pn(k = k, p = 11, type = "PRt", Gauss = TRUE),
               tolerance = 1e-6)

})

test_that("Gegen_coefs_Pn edge cases", {

  expect_error(Gegen_coefs_Pn(k = k, p = 11, type = "PRt2", Gauss = TRUE))

})

test_that("Gegen_coefs_Pn verbose", {

  expect_message(Gegen_coefs_Pn(k = k, p = 2, type = "PCvM",
                                verbose = TRUE), NA)
  expect_message(Gegen_coefs_Pn(k = k, p = 3, type = "PCvM",
                                verbose = TRUE), NA)
  expect_message(Gegen_coefs_Pn(k = k, p = 4, type = "PCvM",
                                verbose = TRUE), NA)
  expect_message(Gegen_coefs_Pn(k = k, p = 5, type = "PCvM",
                                verbose = TRUE))
  expect_message(Gegen_coefs_Pn(k = k, p = 2, type = "PAD",
                                verbose = TRUE), NA)
  expect_message(Gegen_coefs_Pn(k = k, p = 4, type = "PAD",
                                verbose = TRUE))
  expect_message(Gegen_coefs_Pn(k = k, p = 2, type = "PRt",
                                verbose = TRUE), NA)
  expect_message(Gegen_coefs_Pn(k = k, p = 3, type = "PRt",
                                verbose = TRUE), NA)

})

n <- 5
set.seed(123456789)
X2 <- r_unif_sph(n = n, p = 2)
X3 <- r_unif_sph(n = n, p = 3)
X4 <- r_unif_sph(n = n, p = 4)
X5 <- r_unif_sph(n = n, p = 5)
X11 <- r_unif_sph(n = n, p = 11)
Psi2 <- Psi_mat(X2)
Psi3 <- Psi_mat(X3)
Psi4 <- Psi_mat(X4)
Psi5 <- Psi_mat(X5)
Psi11 <- Psi_mat(X11)
dim(Psi2) <- c(dim(Psi2), 1)
dim(Psi3) <- c(dim(Psi3), 1)
dim(Psi4) <- c(dim(Psi4), 1)
dim(Psi5) <- c(dim(Psi5), 1)
dim(Psi11) <- c(dim(Psi11), 1)
Th2 <- X_to_Theta(X2)
Psi_to_mat <- function(Psi) {

  Psi_mat <- matrix(0, nrow = n, ncol = n)
  Psi_mat[upper.tri(Psi_mat)] <- Psi
  return(Psi_mat + t(Psi_mat))

}

test_that("PCvM vs. psi_Pn", {

  expect_equal(2 / n * sum(psi_Pn(theta = Psi2, q = 1, type = "PCvM")) +
                 (3 - 2 * n) / 6,
               drop(sph_stat_PCvM(Psi2, Psi_in_X = TRUE, p = 2)),
               tolerance = 1e-5)
  expect_equal(2 / n * sum(psi_Pn(theta = Psi3, q = 2, type = "PCvM")) +
                 (3 - 2 * n) / 6,
               drop(sph_stat_PCvM(Psi3, Psi_in_X = TRUE, p = 3)),
               tolerance = 1e-5)
  expect_equal(2 / n * sum(psi_Pn(theta = Psi4, q = 3, type = "PCvM")) +
                 (3 - 2 * n) / 6,
               drop(sph_stat_PCvM(Psi4, Psi_in_X = TRUE, p = 4)),
               tolerance = 1e-5)
  expect_equal(2 / n * sum(psi_Pn(theta = Psi5, q = 4, type = "PCvM")) +
                 (3 - 2 * n) / 6,
               drop(sph_stat_PCvM(Psi5, Psi_in_X = TRUE, p = 5)),
               tolerance = 1e-5)
  expect_equal(2 / n * sum(psi_Pn(theta = Psi11, q = 10, type = "PCvM")) +
                 (3 - 2 * n) / 6,
               drop(sph_stat_PCvM(Psi11, Psi_in_X = TRUE, p = 11)),
               tolerance = 5e-5)

})

test_that("PCvM vs. psi_Pn (psi_Gauss = FALSE)", {

  expect_equal(2 / n * sum(psi_Pn(theta = Psi5, q = 4, type = "PCvM",
                                  psi_Gauss = FALSE)) +
                 (3 - 2 * n) / 6,
               drop(sph_stat_PCvM(Psi5, Psi_in_X = TRUE, p = 5)),
               tolerance = 1e-5)
  expect_equal(2 / n * sum(psi_Pn(theta = Psi11, q = 10, type = "PCvM",
                                  psi_Gauss = FALSE)) +
                 (3 - 2 * n) / 6,
               drop(sph_stat_PCvM(Psi11, Psi_in_X = TRUE, p = 11)),
               tolerance = 5e-5)

})

test_that("PAD vs. psi_Pn", {

  expect_equal(2 / n * sum(psi_Pn(theta = Psi2, q = 1, type = "PAD")) + n,
               drop(sph_stat_PAD(Psi2, Psi_in_X = TRUE, p = 2)),
               tolerance = 1e-5)
  expect_equal(2 / n * sum(psi_Pn(theta = Psi3, q = 2, type = "PAD")) + n,
               drop(sph_stat_PAD(Psi3, Psi_in_X = TRUE, p = 3)),
               tolerance = 1e-5)
  expect_equal(2 / n * sum(psi_Pn(theta = Psi4, q = 3, type = "PAD")) + n,
               drop(sph_stat_PAD(Psi4, Psi_in_X = TRUE, p = 4)),
               tolerance = 1e-5)
  expect_equal(2 / n * sum(psi_Pn(theta = Psi5, q = 4, type = "PAD")) + n,
               drop(sph_stat_PAD(Psi5, Psi_in_X = TRUE, p = 5)),
               tolerance = 1e-5)
  expect_equal(2 / n * sum(psi_Pn(theta = Psi11, q = 10, type = "PAD")) + n,
               drop(sph_stat_PAD(Psi11, Psi_in_X = TRUE, p = 11)),
               tolerance = 1e-5)

})

test_that("PAD vs. psi_Pn (psi_Gauss = FALSE)", {

  expect_equal(2 / n * sum(psi_Pn(theta = Psi3, q = 2, type = "PAD",
                                  psi_Gauss = FALSE)) + n,
               drop(sph_stat_PAD(Psi3, Psi_in_X = TRUE, p = 3)),
               tolerance = 1e-5)
  expect_equal(2 / n * sum(psi_Pn(theta = Psi4, q = 3, type = "PAD",
                                  psi_Gauss = FALSE)) + n,
               drop(sph_stat_PAD(Psi4, Psi_in_X = TRUE, p = 4)),
               tolerance = 1e-5)
  expect_equal(2 / n * sum(psi_Pn(theta = Psi5, q = 4, type = "PAD",
                                  psi_Gauss = FALSE)) + n,
               drop(sph_stat_PAD(Psi5, Psi_in_X = TRUE, p = 5)),
               tolerance = 1e-5)
  expect_equal(2 / n * sum(psi_Pn(theta = Psi11, q = 10, type = "PAD",
                                  psi_Gauss = FALSE)) + n,
               drop(sph_stat_PAD(Psi11, Psi_in_X = TRUE, p = 11)),
               tolerance = 1e-5)

})

test_that("PRt vs. psi_Pn", {

  expect_equal(2 / n * sum(psi_Pn(theta = Psi2, q = 1, type = "PRt")) +
                 (1 - n) / 2 + n * (1 / 3 * (1 - 1 / 3)),
               drop(sph_stat_PRt(Psi2, Psi_in_X = TRUE, p = 2)),
               tolerance = 1e-5)
  expect_equal(2 / n * sum(psi_Pn(theta = Psi3, q = 2, type = "PRt")) +
                 (1 - n) / 2 + n * (1 / 3 * (1 - 1 / 3)),
               drop(sph_stat_PRt(Psi3, Psi_in_X = TRUE, p = 3)),
               tolerance = 1e-5)
  expect_equal(2 / n * sum(psi_Pn(theta = Psi4, q = 3, type = "PRt")) +
                 (1 - n) / 2 + n * (1 / 3 * (1 - 1 / 3)),
               drop(sph_stat_PRt(Psi4, Psi_in_X = TRUE, p = 4)),
               tolerance = 1e-5)
  expect_equal(2 / n * sum(psi_Pn(theta = Psi5, q = 4, type = "PRt")) +
                 (1 - n) / 2 + n * (1 / 3 * (1 - 1 / 3)),
               drop(sph_stat_PRt(Psi5, Psi_in_X = TRUE, p = 5)),
               tolerance = 1e-5)
  expect_equal(2 / n * sum(psi_Pn(theta = Psi11, q = 10, type = "PRt")) +
                 (1 - n) / 2 + n * (1 / 3 * (1 - 1 / 3)),
               drop(sph_stat_PRt(Psi11, Psi_in_X = TRUE, p = 11)),
               tolerance = 1e-5)

})

test_that("PRt vs. psi_Pn (psi_Gauss = FALSE)", {

  expect_equal(2 / n * sum(psi_Pn(theta = Psi5, q = 4, type = "PRt",
                                  psi_Gauss = FALSE)) +
                 (1 - n) / 2 + n * (1 / 3 * (1 - 1 / 3)),
               drop(sph_stat_PRt(Psi5, Psi_in_X = TRUE, p = 5)),
               tolerance = 1e-6)
  expect_equal(2 / n * sum(psi_Pn(theta = Psi11, q = 10, type = "PRt",
                                  psi_Gauss = FALSE)) +
                 (1 - n) / 2 + n * (1 / 3 * (1 - 1 / 3)),
               drop(sph_stat_PRt(Psi11, Psi_in_X = TRUE, p = 11)),
               tolerance = 1e-6)

})

test_that("Gegen_coefs_Pn PRt for p = 2)", {

  expect_equal(Gegen_coefs_Pn(k = k1, p = 2, type = "PRt"),
               exp(log(2) - 2 * log(k1 * pi) + log(sin(k1 * pi * 1 / 3)^2)))
  expect_equal(Gegen_coefs_Pn(k = k1, p = 2, type = "PRt", Rothman_t = 1 / 2),
               exp(log(2) - 2 * log(k1 * pi) + log(sin(k1 * pi * 1 / 2)^2)))

})

test_that("PCvM vs. psi_Pn tilde", {

  expect_equal(sum(psi_Pn(theta = c(Psi_to_mat(Psi2)), q = 1, type = "PCvM",
                          tilde = TRUE)) / n,
               drop(sph_stat_PCvM(Psi2, Psi_in_X = TRUE, p = 2)),
               tolerance = 1e-5)
  expect_equal(sum(psi_Pn(theta = c(Psi_to_mat(Psi3)), q = 2, type = "PCvM",
                          tilde = TRUE)) / n,
               drop(sph_stat_PCvM(Psi3, Psi_in_X = TRUE, p = 3)),
               tolerance = 1e-5)
  expect_equal(sum(psi_Pn(theta = c(Psi_to_mat(Psi4)), q = 3, type = "PCvM",
                          tilde = TRUE)) / n,
               drop(sph_stat_PCvM(Psi4, Psi_in_X = TRUE, p = 4)),
               tolerance = 1e-5)
  expect_equal(sum(psi_Pn(theta = c(Psi_to_mat(Psi5)), q = 4, type = "PCvM",
                          tilde = TRUE)) / n,
               drop(sph_stat_PCvM(Psi5, Psi_in_X = TRUE, p = 5)),
               tolerance = 1e-5)
  expect_equal(sum(psi_Pn(theta = c(Psi_to_mat(Psi11)), q = 10, type = "PCvM",
                          tilde = TRUE)) / n,
               drop(sph_stat_PCvM(Psi11, Psi_in_X = TRUE, p = 11)),
               tolerance = 5e-5)

})

test_that("PAD vs. psi_Pn tilde", {

  expect_equal(sum(psi_Pn(theta = c(Psi_to_mat(Psi2)), q = 1, type = "PAD",
                          tilde = TRUE)) / n,
               drop(sph_stat_PAD(Psi2, Psi_in_X = TRUE, p = 2)),
               tolerance = 1e-5)
  expect_equal(sum(psi_Pn(theta = c(Psi_to_mat(Psi3)), q = 2, type = "PAD",
                          tilde = TRUE)) / n,
               drop(sph_stat_PAD(Psi3, Psi_in_X = TRUE, p = 3)),
               tolerance = 5e-5)
  expect_equal(sum(psi_Pn(theta = c(Psi_to_mat(Psi4)), q = 3, type = "PAD",
                          tilde = TRUE)) / n,
               drop(sph_stat_PAD(Psi4, Psi_in_X = TRUE, p = 4)),
               tolerance = 1e-5)
  expect_equal(sum(psi_Pn(theta = c(Psi_to_mat(Psi5)), q = 4, type = "PAD",
                          tilde = TRUE)) / n,
               drop(sph_stat_PAD(Psi5, Psi_in_X = TRUE, p = 5)),
               tolerance = 1e-5)
  expect_equal(sum(psi_Pn(theta = c(Psi_to_mat(Psi11)), q = 10, type = "PAD",
                          tilde = TRUE)) / n,
               drop(sph_stat_PAD(Psi11, Psi_in_X = TRUE, p = 11)),
               tolerance = 1e-5)

})

test_that("PRt vs. psi_Pn tilde", {

  expect_equal(sum(psi_Pn(theta = c(Psi_to_mat(Psi2)), q = 1, type = "PRt",
                          tilde = TRUE)) / n,
               drop(sph_stat_PRt(Psi2, Psi_in_X = TRUE, p = 2)),
               tolerance = 1e-5)
  expect_equal(sum(psi_Pn(theta = c(Psi_to_mat(Psi3)), q = 2, type = "PRt",
                          tilde = TRUE)) / n,
               drop(sph_stat_PRt(Psi3, Psi_in_X = TRUE, p = 3)),
               tolerance = 1e-5)
  expect_equal(sum(psi_Pn(theta = c(Psi_to_mat(Psi4)), q = 3, type = "PRt",
                          tilde = TRUE)) / n,
               drop(sph_stat_PRt(Psi4, Psi_in_X = TRUE, p = 4)),
               tolerance = 1e-5)
  expect_equal(sum(psi_Pn(theta = c(Psi_to_mat(Psi5)), q = 4, type = "PRt",
                          tilde = TRUE)) / n,
               drop(sph_stat_PRt(Psi5, Psi_in_X = TRUE, p = 5)),
               tolerance = 1e-5)
  expect_equal(sum(psi_Pn(theta = c(Psi_to_mat(Psi11)), q = 10, type = "PRt",
                          tilde = TRUE)) / n,
               drop(sph_stat_PRt(Psi11, Psi_in_X = TRUE, p = 11)),
               tolerance = 1e-5)

})

x <- seq(-0.9, 0.9, l = 5)
k0 <- 0

test_that("Gegen_coefs_Pn edge cases", {

  expect_equal(Gegen_coefs_Pn(k = k0, p = 2, type = "PCvM", Gauss = TRUE),
               1 / 3)
  expect_equal(Gegen_coefs_Pn(k = k0, p = 2, type = "PCvM", Gauss = FALSE),
               1 / 3)
  expect_equal(Gegen_coefs_Pn(k = k0, p = 3, type = "PCvM", Gauss = TRUE),
               1 / 3)
  expect_equal(Gegen_coefs_Pn(k = k0, p = 3, type = "PCvM", Gauss = FALSE),
               1 / 3)
  expect_equal(Gegen_coefs_Pn(k = k0, p = 5, type = "PCvM", Gauss = TRUE),
               1 / 3)
  expect_equal(Gegen_coefs_Pn(k = k0, p = 5, type = "PCvM", Gauss = FALSE),
               1 / 3)
  expect_equal(Gegen_coefs_Pn(k = k0, p = 2, type = "PAD", Gauss = TRUE), -1)
  expect_equal(Gegen_coefs_Pn(k = k0, p = 2, type = "PAD", Gauss = FALSE), -1)
  expect_equal(Gegen_coefs_Pn(k = k0, p = 3, type = "PAD", Gauss = TRUE), -1)
  expect_equal(Gegen_coefs_Pn(k = k0, p = 3, type = "PAD", Gauss = FALSE), -1)
  expect_equal(Gegen_coefs_Pn(k = k0, p = 5, type = "PAD", Gauss = TRUE), -1)
  expect_equal(Gegen_coefs_Pn(k = k0, p = 5, type = "PAD", Gauss = FALSE), -1)
  expect_equal(Gegen_coefs_Pn(k = k0, p = 2, type = "PRt", Gauss = TRUE),
               5 / 18)
  expect_equal(Gegen_coefs_Pn(k = k0, p = 2, type = "PRt", Gauss = FALSE),
               5 / 18)
  expect_equal(Gegen_coefs_Pn(k = k0, p = 3, type = "PRt", Gauss = TRUE),
               5 / 18)
  expect_equal(Gegen_coefs_Pn(k = k0, p = 3, type = "PRt", Gauss = FALSE),
               5 / 18)
  expect_equal(Gegen_coefs_Pn(k = k0, p = 5, type = "PRt", Gauss = TRUE),
               5 / 18)
  expect_equal(Gegen_coefs_Pn(k = k0, p = 5, type = "PRt", Gauss = FALSE),
               5 / 18)

})

test_that("psi_Pn edge cases", {

  expect_error(psi_Pn(theta = Psi2, q = 1, type = "PRt2"))

})

test_that("akx definition for p = 2", {

  expect_equal(drop(t(akx(x = x, p = 2, k = k1))),
               drop(1 - Gegen_polyn(theta = acos(x), k = 2 * k1, p = 2)) /
                 (pi * k1)^2, tolerance = 1e-3)

})

test_that("akx definition for p >= 2", {

  for (p in 3:5) {
    for (i in 1:5) {
      expect_equal(
        drop(akx(x = x[i], p = p, k = k1)),
        drop((1 - x[i]^2)^(p - 1) *
               Gegen_polyn(theta = acos(x[i]), k = k1 - 1, p = p + 2)^2 *
               1 / (1 + 2 * k1 / (p - 2)) *
               ((p - 2) / (Gegen_coefs(k = k1, p = p, only_const = TRUE) *
                             k1 * (k1 + p - 2)))^2)
        )
    }
  }
  skip_on_cran()
  for (p in 3:5) {
    for (i in 1:5) {
      expect_equal(
        drop(akx(x = x[i], p = p, k = k1, sqr = TRUE)) /
          sqrt(1 + 2 * k1 / (p - 2)),
        drop(2^(p - 2) * gamma(p / 2)^2 / pi *
               factorial(k1 - 1) / gamma(k1 + p - 1) *
               (1 - x[i]^2)^((p - 1) / 2) *
               Gegen_polyn(theta = acos(x[i]), k = k1 - 1, p = p + 2))
      )
    }
  }

})

test_that("akx edge cases", {

  expect_equal(akx(x = x, k = k, p = 2, sqr = TRUE)^2,
               akx(x = x, k = k, p = 2, sqr = FALSE))
  expect_equal(akx(x = x, k = k, p = 5, sqr = TRUE)^2,
               akx(x = x, k = k, p = 5, sqr = FALSE))
  expect_equal(akx(x = x[1], k = k, p = 2),
               akx(x = x, k = k, p = 2)[1, , drop = FALSE])
  expect_equal(akx(x = x[1], k = k, p = 5),
               akx(x = x, k = k, p = 5)[1, , drop = FALSE])
  expect_equal(akx(x = x, k = k[1], p = 2),
               akx(x = x, k = k, p = 2)[, 1, drop = FALSE])
  expect_equal(akx(x = x, k = k[1], p = 5),
               akx(x = x, k = k, p = 5)[, 1, drop = FALSE])
  expect_equal(akx(x = x[1], k = k[1], p = 2),
               akx(x = x, k = k, p = 2)[1, 1, drop = FALSE])
  expect_equal(akx(x = x[1], k = k[1], p = 5),
               akx(x = x, k = k, p = 5)[1, 1, drop = FALSE])

})

f_PRt_2 <- f_locdev_Pn(p = 2, type = "PRt")
f_PRt_3 <- f_locdev_Pn(p = 3, type = "PRt")
f_PRt_4 <- f_locdev_Pn(p = 4, type = "PRt")
f_PRt_11 <- f_locdev_Pn(p = 11, type = "PRt")
f_PCvM_2 <- f_locdev_Pn(p = 2, type = "PCvM")
f_PCvM_3 <- f_locdev_Pn(p = 3, type = "PCvM")
f_PCvM_4 <- f_locdev_Pn(p = 4, type = "PCvM")
f_PCvM_11 <- f_locdev_Pn(p = 11, type = "PCvM")
f_PAD_2 <- f_locdev_Pn(p = 2, type = "PAD")
f_PAD_3 <- f_locdev_Pn(p = 3, type = "PAD")
f_PAD_4 <- f_locdev_Pn(p = 4, type = "PAD")
f_PAD_11 <- f_locdev_Pn(p = 11, type = "PAD")

test_that("Integral one for f_locdev_Pn with PRt", {

  expect_equal(int_sph_MC(f = function(x) f_PRt_2(x[, 1]), p = 2, M = 1e3,
                          chunks = 1, seeds = 1), 1,
               tolerance = 3e-2)
  expect_equal(int_sph_MC(f = function(x) f_PRt_3(x[, 1]), p = 3, M = 1e3,
                          chunks = 1, seeds = 1), 1,
               tolerance = 3e-2)
  expect_equal(int_sph_MC(f = function(x) f_PRt_4(x[, 1]), p = 4, M = 1e3,
                          chunks = 1, seeds = 1), 1,
               tolerance = 3e-2)
  expect_equal(int_sph_MC(f = function(x) f_PRt_11(x[, 1]), p = 11, M = 1e3,
                          chunks = 1, seeds = 1), 1,
               tolerance = 3e-2)

})

test_that("Integral one for f_locdev_Pn with PCvM", {

  expect_equal(int_sph_MC(f = function(x) f_PCvM_2(x[, 1]), p = 2, M = 1e3,
                          chunks = 1, seeds = 1), 1,
               tolerance = 3e-2)
  expect_equal(int_sph_MC(f = function(x) f_PCvM_3(x[, 1]), p = 3, M = 1e3,
                          chunks = 1, seeds = 1), 1,
               tolerance = 3e-2)
  expect_equal(int_sph_MC(f = function(x) f_PCvM_4(x[, 1]), p = 4, M = 1e3,
                          chunks = 1, seeds = 1), 1,
               tolerance = 3e-2)
  expect_equal(int_sph_MC(f = function(x) f_PCvM_11(x[, 1]), p = 11, M = 1e3,
                          chunks = 1, seeds = 1), 1,
               tolerance = 3e-2)

})

test_that("Integral one for f_locdev_Pn with PAD", {

  expect_equal(int_sph_MC(f = function(x) f_PAD_2(x[, 1]), p = 2, M = 1e3,
                          chunks = 1, seeds = 1), 1,
               tolerance = 7e-2)
  expect_equal(int_sph_MC(f = function(x) f_PAD_3(x[, 1]), p = 3, M = 1e3,
                          chunks = 1, seeds = 1), 1,
               tolerance = 7e-2)
  expect_equal(int_sph_MC(f = function(x) f_PAD_4(x[, 1]), p = 4, M = 1e3,
                          chunks = 1, seeds = 1), 1,
               tolerance = 7e-2)
  expect_equal(int_sph_MC(f = function(x) f_PAD_11(x[, 1]), p = 11, M = 1e3,
                          chunks = 1, seeds = 1), 1,
               tolerance = 7e-2)

})
