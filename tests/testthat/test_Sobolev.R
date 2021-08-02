
K <- 7
k <- 1:K
psi_Ajne <- function(th) 1 / 2 - th / (2 * pi)
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

test_that("d_p_k", {

  expect_equal(d_p_k(p = 2, k = 1:3), rep(2, 3))
  expect_equal(d_p_k(p = 2, k = 1:3, log = TRUE), rep(log(2), 3))
  expect_equal(d_p_k(p = 30, k = 1), 30)

})

test_that("Gegen_coefs (Gauss = TRUE) vs weights_dfs_Sobolev for Ajne", {

  for (p in c(2, 3, 4, 11)) {
    expect_equal(Gegen_coefs(psi = psi_Ajne, k = k[-K], p = p, Gauss = TRUE),
                 switch((p == 2) + 1, (1 + 2 * k[-K] / (p - 2)), 2) *
                   weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                       type = "Ajne", verbose = FALSE)$weights,
                 tolerance = 1e-6)
  }

})

test_that("Gegen_coefs (Gauss = FALSE) vs weights_dfs_Sobolev for Ajne", {

  for (p in c(2, 3, 4, 11)) {
    expect_equal(Gegen_coefs(psi = psi_Ajne, k = k[-K], p = p, Gauss = FALSE),
                 switch((p == 2) + 1, (1 + 2 * k[-K] / (p - 2)), 2) *
                   weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                       type = "Ajne", verbose = FALSE)$weights,
                 tolerance = 1e-6)
  }

})

test_that("Gegen_coefs (Gauss = TRUE) for PCvM", {

  expect_equal(Gegen_coefs(psi = psi_PCvM_2, k = k, p = 2, Gauss = TRUE),
               2 * weights_dfs_Sobolev(p = 2, K_max = K, thre = 0,
                                       type = "PCvM", verbose = FALSE)$weights,
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_PCvM_3, k = k, p = 3, Gauss = TRUE),
               (1 + 2 * k / (3 - 2)) *
                 weights_dfs_Sobolev(p = 3, K_max = K, thre = 0,
                                     type = "PCvM", verbose = FALSE)$weights,
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_PCvM_4, k = k, p = 4, Gauss = TRUE),
               (1 + 2 * k / (4 - 2)) *
                 weights_dfs_Sobolev(p = 4, K_max = K, thre = 0,
                                     type = "PCvM", verbose = FALSE)$weights,
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_PCvM_11, k = k, p = 11, Gauss = TRUE),
               (1 + 2 * k / (11 - 2)) *
                 weights_dfs_Sobolev(p = 11, K_max = K, thre = 0,
                                     type = "PCvM", verbose = FALSE)$weights,
               tolerance = 1e-6)

})

test_that("Gegen_coefs (Gauss = FALSE) for PCvM", {

  expect_equal(Gegen_coefs(psi = psi_PCvM_2, k = k, p = 2, Gauss = FALSE),
               2 * weights_dfs_Sobolev(p = 2, K_max = K, thre = 0,
                                       type = "PCvM", verbose = FALSE)$weights,
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_PCvM_3, k = k, p = 3, Gauss = FALSE),
               (1 + 2 * k / (3 - 2)) *
                 weights_dfs_Sobolev(p = 3, K_max = K, thre = 0,
                                     type = "PCvM", verbose = FALSE)$weights,
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_PCvM_4, k = k, p = 4, Gauss = FALSE),
               (1 + 2 * k / (4 - 2)) *
                 weights_dfs_Sobolev(p = 4, K_max = K, thre = 0,
                                     type = "PCvM", verbose = FALSE)$weights,
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_PCvM_11, k = k, p = 11, Gauss = FALSE),
               (1 + 2 * k / (11 - 2)) *
                 weights_dfs_Sobolev(p = 11, K_max = K, thre = 0,
                                     type = "PCvM", verbose = FALSE)$weights,
               tolerance = 1e-6)

})

test_that("Gegen_coefs (Gauss = TRUE) for PAD", {

  expect_equal(Gegen_coefs(psi = psi_PAD_2, k = k, p = 2, Gauss = TRUE),
               2 * weights_dfs_Sobolev(p = 2, K_max = K, thre = 0,
                                       type = "PAD", verbose = FALSE)$weights,
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_PAD_3, k = k, p = 3, Gauss = TRUE),
               (1 + 2 * k / (3 - 2)) *
                 weights_dfs_Sobolev(p = 3, K_max = K, thre = 0,
                                     type = "PAD", verbose = FALSE)$weights,
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_PAD_4, k = k, p = 4, Gauss = TRUE),
               (1 + 2 * k / (4 - 2)) *
                 weights_dfs_Sobolev(p = 4, K_max = K, thre = 0,
                                     type = "PAD", verbose = FALSE)$weights,
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_PAD_11, k = k, p = 11, Gauss = TRUE),
               (1 + 2 * k / (11 - 2)) *
                 weights_dfs_Sobolev(p = 11, K_max = K, thre = 0,
                                     type = "PAD", verbose = FALSE)$weights,
               tolerance = 1e-6)

})

test_that("Gegen_coefs (Gauss = FALSE) for PAD", {

  expect_equal(Gegen_coefs(psi = psi_PAD_2, k = k, p = 2, Gauss = FALSE),
               2 * weights_dfs_Sobolev(p = 2, K_max = K, thre = 0,
                                       type = "PAD", verbose = FALSE)$weights,
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_PAD_3, k = k, p = 3, Gauss = FALSE),
               (1 + 2 * k / (3 - 2)) *
                 weights_dfs_Sobolev(p = 3, K_max = K, thre = 0,
                                     type = "PAD", verbose = FALSE)$weights,
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_PAD_4, k = k, p = 4, Gauss = FALSE),
               (1 + 2 * k / (4 - 2)) *
                 weights_dfs_Sobolev(p = 4, K_max = K, thre = 0,
                                     type = "PAD", verbose = FALSE)$weights,
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_PAD_11, k = k, p = 11, Gauss = FALSE),
               (1 + 2 * k / (11 - 2)) *
                 weights_dfs_Sobolev(p = 11, K_max = K, thre = 0,
                                     type = "PAD", verbose = FALSE)$weights,
               tolerance = 1e-6)

})

test_that("Gegen_coefs (Gauss = TRUE) for PRt", {

  expect_equal(Gegen_coefs(psi = psi_PRt_2, k = k, p = 2, Gauss = TRUE),
               2 * weights_dfs_Sobolev(p = 2, K_max = K, thre = 0,
                                       type = "PRt", verbose = FALSE)$weights,
               tolerance = 1e-5)
  expect_equal(Gegen_coefs(psi = psi_PRt_3, k = k, p = 3, Gauss = TRUE),
               (1 + 2 * k / (3 - 2)) *
                 weights_dfs_Sobolev(p = 3, K_max = K, thre = 0,
                                     type = "PRt", verbose = FALSE)$weights,
               tolerance = 1e-5)
  expect_equal(Gegen_coefs(psi = psi_PRt_4, k = k, p = 4, Gauss = TRUE),
               (1 + 2 * k / (4 - 2)) *
                 weights_dfs_Sobolev(p = 4, K_max = K, thre = 0,
                                     type = "PRt", verbose = FALSE)$weights,
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_PRt_11, k = k, p = 11, Gauss = TRUE),
               (1 + 2 * k / (11 - 2)) *
                 weights_dfs_Sobolev(p = 11, K_max = K, thre = 0,
                                     type = "PRt", verbose = FALSE)$weights,
               tolerance = 1e-6)

})

test_that("Gegen_coefs (Gauss = FALSE) for PRt", {

  expect_equal(Gegen_coefs(psi = psi_PRt_2, k = k, p = 2, Gauss = FALSE),
               2 * weights_dfs_Sobolev(p = 2, K_max = K, thre = 0,
                                       type = "PRt", verbose = FALSE)$weights,
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_PRt_3, k = k, p = 3, Gauss = FALSE),
               (1 + 2 * k / (3 - 2)) *
                 weights_dfs_Sobolev(p = 3, K_max = K, thre = 0,
                                     type = "PRt", verbose = FALSE)$weights,
               tolerance = 5e-5)
  expect_equal(Gegen_coefs(psi = psi_PRt_4, k = k, p = 4, Gauss = FALSE),
               (1 + 2 * k / (4 - 2)) *
                 weights_dfs_Sobolev(p = 4, K_max = K, thre = 0,
                                     type = "PRt", verbose = FALSE)$weights,
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_PRt_11, k = k, p = 11, Gauss = FALSE),
               (1 + 2 * k / (11 - 2)) *
                 weights_dfs_Sobolev(p = 11, K_max = K, thre = 0,
                                     type = "PRt", verbose = FALSE)$weights,
               tolerance = 1e-6)

})

test_that("Gegen_coefs_Pn vs weights_dfs_Sobolev for PCvM", {

  for (p in c(2, 3, 4, 11)) {
    expect_equal(Gegen_coefs_Pn(k = k, p = p, type = "PCvM"),
                 switch((p == 2) + 1, (1 + 2 * k / (p - 2)), 2) *
                   weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                       type = "PCvM", verbose = FALSE)$weights)
  }

})

test_that("Gegen_coefs_Pn vs weights_dfs_Sobolev for PAD", {

  for (p in c(2, 3, 4, 11)) {
    expect_equal(Gegen_coefs_Pn(k = k, p = p, type = "PAD"),
                 switch((p == 2) + 1, (1 + 2 * k / (p - 2)), 2) *
                   weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                       type = "PAD", verbose = FALSE)$weights)
  }

})

test_that("Gegen_coefs_Pn vs weights_dfs_Sobolev for PRt", {

  for (p in c(2, 3, 4, 11)) {
    expect_equal(Gegen_coefs_Pn(k = k, p = p, type = "PRt"),
                 switch((p == 2) + 1, (1 + 2 * k / (p - 2)), 2) *
                   weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                       type = "PRt", verbose = FALSE)$weights)
  }

})

test_that("weights_dfs_Sobolev for Ajne vs PRt with t = 1 / 2", {

  for (p in c(2, 3, 4, 11)) {
    expect_equal(weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                     type = "Ajne", verbose = FALSE)$weights,
                 weights_dfs_Sobolev(p = p, K_max = K - 1, thre = 0,
                                     type = "PRt", Rothman_t = 1 / 2,
                                     verbose = FALSE)$weights,
                 tolerance = 1e-6)
  }

})
