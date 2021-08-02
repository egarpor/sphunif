
Sys.unsetenv("R_TESTS")

set.seed(21332)
u <- seq(0, 1, l = 13)
x <- seq(-1, 1, l = 13)
v <- runif(1e3)
f0 <- function(x) rep(1, length(x))
f1 <- function(x, kappa) exp(kappa * x)
f2 <- function(x, kappa) exp(kappa * x^2)
f3 <- function(x, kappa, nu) exp(-kappa * (x - nu)^2)
f4 <- function(x, kappa, q) {
  rho <- ((2 * kappa + 1) - sqrt(4 * kappa + 1)) / (2 * kappa)
  (1 - rho^2) / (1 + rho^2 - 2 * rho * x)^((q + 1) / 2) / 
    rotasym::w_p(p = q + 1)
}

test_that("Correct integration of con_f", {

  for (p in c(2:4, 11)) {
    expect_equal(con_f(f = function(x)
      rotasym::g_vMF(t = x, p = p, kappa = 3, scaled = TRUE),
      p = p, N = 320), 1)
    expect_equal(con_f(f = function(x) f4(x, kappa = 1, q = p - 1),
                       p = p, N = 320), 1)
  }

})

test_that("d_locdev", {

  for (p in 2:4) {
    xp <- r_unif_sph(n = 5, p = p)[, , 1]
    mu <- c(rep(0, p - 1), 1)
    expect_equal(d_locdev(x = xp[1, , drop = FALSE], mu = mu, kappa = 0.25,
                          f = function(z)
                            rotasym::g_vMF(t = z, kappa = 3, p = p)),
                 unname(d_locdev(x = xp[1, ], mu = mu, kappa = 0.25,
                                 f = function(z) 
                                   rotasym::g_vMF(t = z, kappa = 3, p = p))))
    expect_equal(d_locdev(x = xp, mu = mu, kappa = 0, f = NULL),
                 rep(1 / rotasym::w_p(p = p), nrow(xp)))
    expect_equal(d_locdev(x = xp, mu = mu, kappa = 0.25,
                          f = function(z)
                            rotasym::g_vMF(t = z, kappa = 3, p = p)),
                 0.25 * rotasym::g_vMF(t = xp[, p], kappa = 3, p = p) +
                 0.75 / rotasym::w_p(p = p))
  }

})

test_that("r_locdev coherence with d_locdev", {

  for (p in 2:4) {
    mu <- c(rep(0, p - 1), 1)
    samp_1 <- r_locdev(n = 100, mu = mu, kappa = 0.25, 
                       f = function(z) f4(x = z, kappa = 3, q = p - 1))[, p]
    samp_2 <- F_inv_from_f(f = function(z)
      0.25 * f4(x = z, kappa = 3, q = p - 1) + 0.75 / rotasym::w_p(p = p),
      p = p)(runif(n = 100))
    expect_gt(ks.test(samp_1, samp_2)$p.value, 0.01)
    samp_1 <- r_locdev(n = 100, mu = mu, kappa = 0, f = NULL)[, 1]
    samp_2 <- r_unif_sph(n = 100, p = p)[, 1, 1]
    expect_gt(ks.test(samp_1, samp_2)$p.value, 0.01)
  }

})

test_that("Edge cases d_locdev and r_locdev", {

  expect_error(d_locdev(x = 1, mu = 1, kappa = -1, f = NULL))
  expect_error(d_locdev(x = 1:2, mu = 1:3, kappa = -1, f = NULL))
  expect_error(r_locdev(n = 1, mu = 1, kappa = -1))

})

test_that("F_from_f via Gauss--Legendre", {

  for (p in c(2:4, 11)) {
    expect_equal(F_from_f(f = f0, p = p, Gauss = TRUE, K = 1e2)(x),
                 drop(p_proj_unif(x = x, p = p)), tolerance = 1e-3)
  }

})

test_that("F_from_f via integrate()", {

  for (p in c(2:4, 11)) {
    expect_equal(F_from_f(f = f0, p = p, Gauss = FALSE, K = 1e2)(x),
                 drop(p_proj_unif(x = x, p = p)), tolerance = 1e-3)
  }

})

test_that("F_from_f for vMF", {

  for (p in c(2:4, 11)) {
    samp_g <- drop(rotasym::r_g_vMF(n = 100, p = p, kappa = 3))
    expect_gt(ks.test(x = F_from_f(f = f1, p = p, Gauss = TRUE,
                                   K = 1e2, kappa = 3)(samp_g),
                      y = "punif")$p.value, 0.01)
  }

})

test_that("F_inv_from_f via Gauss--Legendre", {

  for (p in c(2:4, 11)) {
    expect_equal(F_inv_from_f(f = f0, p = p, Gauss = TRUE, K = 1e2)(u),
                 drop(q_proj_unif(u = u, p = p)), tolerance = 5e-3)
  }

})

test_that("F_inv_from_f via integrate()", {

  for (p in c(2:4, 11)) {
    expect_equal(F_inv_from_f(f = f0, p = p, Gauss = FALSE, K = 1e2)(u),
                 drop(q_proj_unif(u = u, p = p)), tolerance = 5e-3)
  }

})

test_that("F_inv_from_f for vMF", {

  expect_gt(ks.test(x = F_inv_from_f(f = f1, p = 2, Gauss = TRUE,
                                     K = 1e2, kappa = 3)(v),
                    y = rotasym::r_g_vMF(n = 100, p = 2,
                                         kappa = 3))$p.value, 0.01)
  expect_gt(ks.test(x = F_inv_from_f(f = f1, p = 3, Gauss = TRUE,
                                      K = 1e2, kappa = 5)(v),
                    y = rotasym::r_g_vMF(n = 100, p = 3,
                                          kappa = 5))$p.value, 0.01)
  expect_gt(ks.test(x = F_inv_from_f(f = f1, p = 4, Gauss = TRUE,
                                     K = 1e2, kappa = 5)(v),
                    y = rotasym::r_g_vMF(n = 100, p = 4,
                                         kappa = 5))$p.value, 0.01)
  expect_gt(ks.test(x = F_inv_from_f(f = f1, p = 5, Gauss = TRUE,
                                      K = 1e2, kappa = 10)(v),
                    y = rotasym::r_g_vMF(n = 100, p = 5,
                                         kappa = 10))$p.value, 0.01)
  expect_gt(ks.test(x = F_inv_from_f(f = f1, p = 11, Gauss = TRUE,
                                     K = 1e2, kappa = 20)(v),
                    y = rotasym::r_g_vMF(n = 100, p = 11,
                                         kappa = 20))$p.value, 0.01)

})

test_that("r_alt rotationally symmetric", {

  for (p in 2:4) {
    samp_g <- r_alt(n = 100, p = p, M = 1, kappa = 2, scenario = "vMF")[, p, 1]
    expect_gt(ks.test(x = F_from_f(f = f1, p = p, kappa = 2)(samp_g),
                      y = "punif")$p.value, 0.01)
    samp_g <- r_alt(n = 100, p = p, M = 1, kappa = 2, scenario = "W")[, p, 1]
    expect_gt(ks.test(x = F_from_f(f = f2, p = p, kappa = 2)(samp_g),
                      y = "punif")$p.value, 0.01)
    samp_g <- r_alt(n = 100, p = p, M = 1, kappa = 2, nu = 0.5,
                    scenario = "SC")[, p, 1]
    expect_gt(ks.test(x = F_from_f(f = f3, p = p, kappa = 2, nu = 0.5)(samp_g),
                      y = "punif")$p.value, 0.01)
    samp_g <- r_alt(n = 100, p = p, M = 1, kappa = 2, scenario = "C")[, p, 1]
    expect_gt(ks.test(x = F_from_f(f = f4, p = p, kappa = 2, q = p - 1)(samp_g),
                      y = "punif")$p.value, 0.01)
  }

})

test_that("r_alt non-rotationally symmetric", {

  for (p in 2:4) {
    samp_1 <- r_alt(n = 100, p = p, M = 1, kappa = 1, scenario = "MvMF")[, p, 1]
    samp_2 <- apply(diag(rep(1, p)), 1, function(mu) 
      rotasym::r_vMF(n = round(100 / p), mu = mu, kappa = 1))[, p]
    expect_gt(ks.test(x = samp_1, y = samp_2)$p.value, 0.01)
    samp_1 <- r_alt(n = 100, p = p, M = 1, kappa = 1, scenario = "ACG")[, p, 1]
    samp_2 <- mvtnorm::rmvnorm(n = 100, mean = rep(0, p),
                               sigma = diag(c(rep(1, p - 1), 1 + 1)))
    samp_2 <- samp_2 / sqrt(rowSums(samp_2^2))
    samp_2 <- samp_2[, p]
    expect_gt(ks.test(x = samp_1, y = samp_2)$p.value, 0.01)
  }
  expect_error(r_alt(n = 100, p = p, M = 1, kappa = 1, scenario = "WC"))

})

set.seed(12311)
n <- 20
X_2 <- r_unif_sph(n = n, p = 2, M = 1)
X_3 <- r_unif_sph(n = n, p = 3, M = 1)
X_4 <- r_unif_sph(n = n, p = 4, M = 1)
X_11 <- r_unif_sph(n = n, p = 11, M = 1)
th_k <- Gauss_Legen_nodes(a = 0, b = 2 * pi, N = 320)
w_k <- Gauss_Legen_weights(a = 0, b = 2 * pi, N = 320)
z <- seq(-1, 1, l = 1e3)
z8 <- seq(-0.8, 1, l = 1e3)

# PCvM
k_PCvM <- 1:1e3
uk_PCvM_2 <- bk_to_uk(Gegen_coefs_Pn(k = k_PCvM, p = 2, type = "PCvM",
                                     N = 0), p = 2)
f_locdev_PCvM_2 <- function(z) f_locdev(z = z, p = 2, uk = uk_PCvM_2)
integrand_vec_PCvM_2 <- function(x) {
  f_gamma <- matrix(f_locdev_PCvM_2(c(X_2[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 2) * f_gamma - 1)^2 / n / rotasym::w_p(p = 2)
}
f_locdev_PCvM_2_anal <- f_locdev_Pn(p = 2, type = "PCvM")
integrand_vec_PCvM_2_anal <- function(x) {
  f_gamma <- matrix(f_locdev_PCvM_2_anal(c(X_2[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 2) * f_gamma - 1)^2 / n / rotasym::w_p(p = 2)
}
uk_PCvM_3 <- bk_to_uk(Gegen_coefs_Pn(k = k_PCvM, p = 3, type = "PCvM",
                                     N = 0), p = 3)
f_locdev_PCvM_3 <- function(z) f_locdev(z = z, p = 3, uk = uk_PCvM_3)
integrand_vec_PCvM_3 <- function(x) {
  f_gamma <- matrix(f_locdev_PCvM_3(c(X_3[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 3) * f_gamma - 1)^2 / n / rotasym::w_p(p = 3)
}
uk_PCvM_4 <- bk_to_uk(Gegen_coefs_Pn(k = k_PCvM, p = 4, type = "PCvM",
                                     N = 0), p = 4)
f_locdev_PCvM_4 <- function(z) f_locdev(z = z, p = 4, uk = uk_PCvM_4)
integrand_vec_PCvM_4 <- function(x) {
  f_gamma <- matrix(f_locdev_PCvM_4(c(X_4[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 4) * f_gamma - 1)^2 / n / rotasym::w_p(p = 4)
}
uk_PCvM_11 <- bk_to_uk(pmax(Gegen_coefs_Pn(k = k_PCvM, p = 11, type = "PCvM",
                                           N = 5120), 0), p = 11)
f_locdev_PCvM_11 <- function(z) f_locdev(z = z, p = 11, uk = uk_PCvM_11)
integrand_vec_PCvM_11 <- function(x) {
  f_gamma <- matrix(f_locdev_PCvM_11(c(X_11[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 11) * f_gamma - 1)^2 / n / rotasym::w_p(p = 11)
}

# PAD
k_PAD <- 1:1e3
uk_PAD_2 <- bk_to_uk(pmax(Gegen_coefs_Pn(k = k_PAD, p = 2, type = "PAD",
                                         N = 5120), 0), p = 2)
f_locdev_PAD_2 <- function(z) f_locdev(z = z, p = 2, uk = uk_PAD_2)
integrand_vec_PAD_2 <- function(x) {
  f_gamma <- matrix(f_locdev_PAD_2(c(X_2[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 2) * f_gamma - 1)^2 / n / rotasym::w_p(p = 2)
}
uk_PAD_3 <- bk_to_uk(pmax(Gegen_coefs_Pn(k = k_PAD, p = 3, type = "PAD",
                                         N = 5120), 0), p = 3)
f_locdev_PAD_3 <- function(z) f_locdev(z = z, p = 3, uk = uk_PAD_3)
integrand_vec_PAD_3 <- function(x) {
  f_gamma <- matrix(f_locdev_PAD_3(c(X_3[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 3) * f_gamma - 1)^2 / n / rotasym::w_p(p = 3)
}
uk_PAD_4 <- bk_to_uk(pmax(Gegen_coefs_Pn(k = k_PAD, p = 4, type = "PAD",
                                         N = 5120), 0), p = 4)
f_locdev_PAD_4 <- function(z) f_locdev(z = z, p = 4, uk = uk_PAD_4)
integrand_vec_PAD_4 <- function(x) {
  f_gamma <- matrix(f_locdev_PAD_4(c(X_4[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 4) * f_gamma - 1)^2 / n / rotasym::w_p(p = 4)
}
uk_PAD_11 <- bk_to_uk(pmax(Gegen_coefs_Pn(k = k_PAD, p = 11, type = "PAD",
                                          N = 5120), 0), p = 11)
f_locdev_PAD_11 <- function(z) f_locdev(z = z, p = 11, uk = uk_PAD_11)
integrand_vec_PAD_11 <- function(x) {
  f_gamma <- matrix(f_locdev_PAD_11(c(X_11[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 11) * f_gamma - 1)^2 / n / rotasym::w_p(p = 11)
}

# PRt 1 / 3
k_PRt <- 1:1e3
uk_PRt_2 <- bk_to_uk(Gegen_coefs_Pn(k = k_PRt, p = 2, type = "PRt",
                                    N = 0, Rothman_t = 1 / 3),
                     p = 2)
f_locdev_PRt_2 <- function(z) f_locdev(z = z, p = 2, uk = uk_PRt_2)
integrand_vec_PRt_2 <- function(x) {
  f_gamma <- matrix(f_locdev_PRt_2(c(X_2[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 2) * f_gamma - 1)^2 / n / rotasym::w_p(p = 2)
}
f_locdev_PRt_2_anal <- f_locdev_Pn(p = 2, type = "PRt", Rothman_t = 1 / 3)
integrand_vec_PRt_2_anal <- function(x) {
  f_gamma <- matrix(f_locdev_PRt_2_anal(c(X_2[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 2) * f_gamma - 1)^2 / n / rotasym::w_p(p = 2)
}
uk_PRt_3 <- bk_to_uk(Gegen_coefs_Pn(k = k_PRt, p = 3, type = "PRt",
                                    N = 0, Rothman_t = 1 / 3),
                     p = 3)
f_locdev_PRt_3 <- function(z) f_locdev(z = z, p = 3, uk = uk_PRt_3)
integrand_vec_PRt_3 <- function(x) {
  f_gamma <- matrix(f_locdev_PRt_3(c(X_3[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 3) * f_gamma - 1)^2 / n / rotasym::w_p(p = 3)
}
f_locdev_PRt_3_anal <- f_locdev_Pn(p = 3, type = "PRt", Rothman_t = 1 / 3)
integrand_vec_PRt_3_anal <- function(x) {
  f_gamma <- matrix(f_locdev_PRt_3_anal(c(X_3[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 3) * f_gamma - 1)^2 / n / rotasym::w_p(p = 3)
}
uk_PRt_4 <- bk_to_uk(Gegen_coefs_Pn(k = k_PRt, p = 4, type = "PRt",
                                    N = 0, Rothman_t = 1 / 3),
                     p = 4)
f_locdev_PRt_4 <- function(z) f_locdev(z = z, p = 4, uk = uk_PRt_4)
integrand_vec_PRt_4 <- function(x) {
  f_gamma <- matrix(f_locdev_PRt_4(c(X_4[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 4) * f_gamma - 1)^2 / n / rotasym::w_p(p = 4)
}
f_locdev_PRt_4_anal <- f_locdev_Pn(p = 4, type = "PRt", Rothman_t = 1 / 3)
integrand_vec_PRt_4_anal <- function(x) {
  f_gamma <- matrix(f_locdev_PRt_4_anal(c(X_4[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 4) * f_gamma - 1)^2 / n / rotasym::w_p(p = 4)
}
uk_PRt_11 <- bk_to_uk(Gegen_coefs_Pn(k = k_PRt, p = 11, type = "PRt",
                                     N = 0, Rothman_t = 1 / 3),
                      p = 11)
f_locdev_PRt_11 <- function(z) f_locdev(z = z, p = 11, uk = uk_PRt_11)
integrand_vec_PRt_11 <- function(x) {
  f_gamma <- matrix(f_locdev_PRt_11(c(X_11[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 11) * f_gamma - 1)^2 / n / rotasym::w_p(p = 11)
}
f_locdev_PRt_11_anal <- f_locdev_Pn(p = 11, type = "PRt", Rothman_t = 1 / 3)
integrand_vec_PRt_11_anal <- function(x) {
  f_gamma <- matrix(f_locdev_PRt_11_anal(c(X_11[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 11) * f_gamma - 1)^2 / n / rotasym::w_p(p = 11)
}

# Ajne
k_Ajne <- 1:1e3
vk2_Ajne_2 <- weights_dfs_Sobolev(p = 2, K_max = max(k_Ajne),
                                  type = "Ajne", thre = 0,
                                  verbose = FALSE)$weights
uk_Ajne_2 <- vk2_to_uk(vk2 = vk2_Ajne_2, p = 2,
                       signs = (-1)^((seq_along(vk2_Ajne_2) - 1) %/% 2))
f_locdev_Ajne_2 <- function(z) f_locdev(z = z, p = 2, uk = uk_Ajne_2)
integrand_vec_Ajne_2 <- function(x) {
  f_gamma <- matrix(f_locdev_Ajne_2(c(X_2[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 2) * f_gamma - 1)^2 / n / rotasym::w_p(p = 2)
}
vk2_Ajne_3 <- weights_dfs_Sobolev(p = 3, K_max = max(k_Ajne),
                                  type = "Ajne", thre = 0,
                                  verbose = FALSE)$weights
uk_Ajne_3 <- vk2_to_uk(vk2 = vk2_Ajne_3 , p = 3)
f_locdev_Ajne_3 <- function(z) f_locdev(z = z, p = 3, uk = uk_Ajne_3)
integrand_vec_Ajne_3 <- function(x) {
  f_gamma <- matrix(f_locdev_Ajne_3(c(X_3[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 3) * f_gamma - 1)^2 / n / rotasym::w_p(p = 3)
}
vk2_Ajne_4 <- weights_dfs_Sobolev(p = 4, K_max = max(k_Ajne),
                                  type = "Ajne", thre = 0,
                                  verbose = FALSE)$weights
uk_Ajne_4 <- vk2_to_uk(vk2 = vk2_Ajne_4 , p = 4)
f_locdev_Ajne_4 <- function(z) f_locdev(z = z, p = 4, uk = uk_Ajne_4)
integrand_vec_Ajne_4 <- function(x) {
  f_gamma <- matrix(f_locdev_Ajne_4(c(X_4[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 4) * f_gamma - 1)^2 / n / rotasym::w_p(p = 4)
}
vk2_Ajne_11 <- weights_dfs_Sobolev(p = 11, K_max = max(k_Ajne),
                                   type = "Ajne", thre = 0,
                                   verbose = FALSE)$weights
uk_Ajne_11 <- vk2_to_uk(vk2 = vk2_Ajne_11 , p = 11)
f_locdev_Ajne_11 <- function(z) f_locdev(z = z, p = 11, uk = uk_Ajne_11)
integrand_vec_Ajne_11 <- function(x) {
  f_gamma <- matrix(f_locdev_Ajne_11(c(X_11[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 11) * f_gamma - 1)^2 / n / rotasym::w_p(p = 11)
}

# Ajne analytical
k_Ajne <- 1:1e3
f_Ajne_2 <- function(z) (1 * (z >= 0) + 0.5) / rotasym::w_p(p = 2)
integrand_vec_f_Ajne_2  <- function(x) {
  f_gamma <- matrix(f_Ajne_2(c(X_2[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 2) * f_gamma - 1)^2 / n / rotasym::w_p(p = 2)
}
f_Ajne_3 <- function(z) (1 * (z >= 0) + 0.5) / rotasym::w_p(p = 3)
integrand_vec_f_Ajne_3 <- function(x) {
  f_gamma <- matrix(f_Ajne_3(c(X_3[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 3) * f_gamma - 1)^2 / n / rotasym::w_p(p = 3)
}
f_Ajne_4 <- function(z) (1 * (z >= 0) + 0.5) / rotasym::w_p(p = 4)
integrand_vec_f_Ajne_4 <- function(x) {
  f_gamma <- matrix(f_Ajne_4(c(X_4[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 4) * f_gamma - 1)^2 / n / rotasym::w_p(p = 4)
}
f_Ajne_11 <- function(z) (1 * (z >= 0) + 0.5) / rotasym::w_p(p = 11)
integrand_vec_f_Ajne_11 <- function(x) {
  f_gamma <- matrix(f_Ajne_11(c(X_11[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 11) * f_gamma - 1)^2 / n / rotasym::w_p(p = 11)
}

test_that("PCvM as the integral of its local alternative", {

  expect_equal(drop(sph_stat_PCvM(X = X_2)),
               sum(w_k * integrand_vec_PCvM_2(cbind(cos(th_k), sin(th_k)))),
               tolerance = 5e-2)
  expect_equal(drop(sph_stat_PCvM(X = X_2)),
               integrate(function(th)
                 integrand_vec_PCvM_2_anal(x = cbind(cos(th), sin(th))),
                 lower = 0, upper = 2 * pi, abs.tol = 1e-5,
                 subdivisions = 1e3)$value,
               tolerance = 5e-2)
  expect_equal(drop(sph_stat_PCvM(X = X_3)), {
    int_sph_MC(f = integrand_vec_PCvM_3, p = 3, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16, M = 1e4, verbose = FALSE)
  }, tolerance = 5e-2)
  expect_equal(drop(sph_stat_PCvM(X = X_4)), {
    int_sph_MC(f = integrand_vec_PCvM_4, p = 4, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16, M = 1e4, verbose = FALSE)
  }, tolerance = 5e-2)
  expect_equal(drop(sph_stat_PCvM(X = X_11)), {
    int_sph_MC(f = integrand_vec_PCvM_11, p = 11, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16 + 10, M = 1e4, verbose = FALSE)
  }, tolerance = 1e-1)

})

test_that("PAD as the integral of its local alternative", {

  expect_equal(drop(sph_stat_PAD(X = X_2)),
               sum(w_k * integrand_vec_PAD_2(cbind(cos(th_k), sin(th_k)))),
               tolerance = 1e-1)
  expect_equal(drop(sph_stat_PAD(X = X_3)), {
    int_sph_MC(f = integrand_vec_PAD_3, p = 3, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16, M = 1e4, verbose = FALSE)
  }, tolerance = 5e-2)
  expect_equal(drop(sph_stat_PAD(X = X_4)), {
    int_sph_MC(f = integrand_vec_PAD_4, p = 4, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16, M = 1e4, verbose = FALSE)
  }, tolerance = 2e-2)
  expect_equal(drop(sph_stat_PAD(X = X_11)), {
    int_sph_MC(f = integrand_vec_PAD_11, p = 11, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16 + 20, M = 1e4, verbose = FALSE)
  }, tolerance = 3e-2)

})

test_that("PRt t = 1 / 3 as the integral of its local alternative", {

  expect_equal(drop(sph_stat_PRt(X = X_2, t = 1 / 3)),
               sum(w_k * integrand_vec_PRt_2(cbind(cos(th_k), sin(th_k)))),
               tolerance = 5e-2)
  expect_equal(drop(sph_stat_PRt(X = X_2, t = 1 / 3)),
               integrate(function(th)
                 integrand_vec_PRt_2_anal(x = cbind(cos(th), sin(th))),
                 lower = 0, upper = 2 * pi, abs.tol = 1e-5,
                 subdivisions = 1e3)$value,
               tolerance = 1e-4)
  expect_equal(drop(sph_stat_PRt(X = X_3, t = 1 / 3)), {
    int_sph_MC(f = integrand_vec_PRt_3, p = 3, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16, M = 1e4, verbose = FALSE)
  }, tolerance = 5e-2)
  expect_equal(drop(sph_stat_PRt(X = X_3, t = 1 / 3)), {
    int_sph_MC(f = integrand_vec_PRt_3_anal, p = 3, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16, M = 1e4, verbose = FALSE)
  }, tolerance = 3e-2)
  expect_equal(drop(sph_stat_PRt(X = X_4, t = 1 / 3)), {
    int_sph_MC(f = integrand_vec_PRt_4, p = 4, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16, M = 1e4, verbose = FALSE)
  }, tolerance = 5e-2)
  expect_equal(drop(sph_stat_PRt(X = X_4, t = 1 / 3)), {
    int_sph_MC(f = integrand_vec_PRt_4_anal, p = 4, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16, M = 1e4, verbose = FALSE)
  }, tolerance = 3e-2)
  expect_equal(drop(sph_stat_PRt(X = X_11, t = 1 / 3)), {
    int_sph_MC(f = integrand_vec_PRt_11, p = 11, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16 + 10, M = 1e4, verbose = FALSE)
  }, tolerance = 5e-2)
  expect_equal(drop(sph_stat_PRt(X = X_11, t = 1 / 3)), {
    int_sph_MC(f = integrand_vec_PRt_11_anal, p = 11, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16, M = 1e4, verbose = FALSE)
  }, tolerance = 3e-2)

})

test_that("Ajne as the integral of its local alternative (series expansion)", {

  expect_equal(drop(sph_stat_Ajne(X = X_2)),
               sum(w_k * integrand_vec_Ajne_2(cbind(cos(th_k), sin(th_k)))),
               tolerance = 5e-2)
  expect_equal(drop(sph_stat_Ajne(X = X_3)), {
    int_sph_MC(f = integrand_vec_Ajne_3, p = 3, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16, M = 1e4, verbose = FALSE)
  }, tolerance = 5e-2)
  expect_equal(drop(sph_stat_Ajne(X = X_4)), {
    int_sph_MC(f = integrand_vec_Ajne_4, p = 4, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16, M = 1e4, verbose = FALSE)
  }, tolerance = 5e-2)
  expect_equal(drop(sph_stat_Ajne(X = X_11)), {
    int_sph_MC(f = integrand_vec_Ajne_11, p = 11, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16 + 10, M = 1e4, verbose = FALSE)
  }, tolerance = 1e-1)

})

test_that("Ajne as the integral of its local alternative (analytical)", {

  expect_equal(drop(sph_stat_Ajne(X = X_2)),
               sum(w_k * integrand_vec_f_Ajne_2(cbind(cos(th_k), sin(th_k)))),
               tolerance = 5e-2)
  expect_equal(drop(sph_stat_Ajne(X = X_3)), {
    int_sph_MC(f = integrand_vec_f_Ajne_3, p = 3, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16, M = 1e4, verbose = FALSE)
  }, tolerance = 5e-2)
  expect_equal(drop(sph_stat_Ajne(X = X_4)), {
    int_sph_MC(f = integrand_vec_f_Ajne_4, p = 4, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16, M = 1e4, verbose = FALSE)
  }, tolerance = 5e-2)
  expect_equal(drop(sph_stat_Ajne(X = X_11)), {
    int_sph_MC(f = integrand_vec_f_Ajne_11, p = 11, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16 + 10, M = 1e4, verbose = FALSE)
  }, tolerance = 3e-2)

})

test_that("Check positivity of f_PCvM and f_PAD", {

  expect_true(all(f_locdev_PCvM_2(z = z) > 0))
  expect_true(all(f_locdev_PCvM_3(z = z) > 0))
  expect_true(all(f_locdev_PCvM_4(z = z) > 0))
  expect_true(all(f_locdev_PCvM_11(z = z8) > 0))
  expect_true(all(f_locdev_PAD_2(z = z) > 0))
  expect_true(all(f_locdev_PAD_3(z = z) > 0))
  expect_true(all(f_locdev_PAD_4(z = z) > 0))
  expect_true(all(f_locdev_PAD_11(z = z8) > 0))

})

vk2 <- 1:5
bk <- 5:1
uk <- seq(-1, 1, l = 5)

test_that("Conversion between coefficients", {

  expect_equal(bk_to_vk2(vk2_to_bk(vk2, p = 2), p = 2), vk2)
  expect_equal(bk_to_vk2(vk2_to_bk(vk2, p = 3), p = 3), vk2)
  expect_equal(bk_to_uk(uk_to_bk(uk, p = 2), p = 2, signs = sign(uk)), uk)
  expect_equal(bk_to_uk(uk_to_bk(uk, p = 3), p = 3, signs = sign(uk)), uk)
  expect_equal(vk2_to_bk(bk_to_vk2(bk, p = 2), p = 2), bk)
  expect_equal(vk2_to_bk(bk_to_vk2(bk, p = 3), p = 3), bk)
  expect_equal(vk2_to_uk(uk_to_vk2(uk, p = 2), p = 2, signs = sign(uk)), uk)
  expect_equal(vk2_to_uk(uk_to_vk2(uk, p = 3), p = 3, signs = sign(uk)), uk)
  expect_equal(uk_to_vk2(vk2_to_uk(vk2, p = 2), p = 2), vk2)
  expect_equal(uk_to_vk2(vk2_to_uk(vk2, p = 3), p = 3), vk2)
  expect_equal(uk_to_bk(bk_to_uk(bk, p = 2), p = 2), bk)
  expect_equal(uk_to_bk(bk_to_uk(bk, p = 3), p = 3), bk)

})

test_that("Conversion of bk to vk2 in projected-ecdf statistics", {

  for (type in c("PCvM", "PAD", "PRt")) {
    for (p in c(2, 3, 4, 9)) {
      expect_equal(bk_to_vk2(Gegen_coefs_Pn(k = 1:9, p = p, type = type),
                             p = p),
                   weights_dfs_Sobolev(K_max = 9, thre = 0, p = p, type = type,
                                       verbose = FALSE)$weights)
    }
  }

})

test_that("Conversion of bk to uk in projected-ecdf statistics", {

  for (type in c("PCvM", "PAD", "PRt")) {
    for (p in c(2, 3, 4, 9)) {
      expect_equal(bk_to_uk(Gegen_coefs_Pn(k = 1:9, p = p, type = type, 
                                           Rothman_t = 0.1), p = p),
                   abs(cutoff_locdev(K_max = 9, thre = 0, p = p, type = type,
                                     Rothman_t = 0.1)))
    }
  }

})

test_that("cutoff_locdev verbose", {

  expect_message(cutoff_locdev(K_max = 9, thre = 0, p = 9, type = "PAD",
                               verbose = 1))
  suppressMessages(
    expect_message(cutoff_locdev(K_max = 1e1, thre = 1e-4, p = 4, type = "PCvM",
                                 verbose = 2)))

})

# MJ (2000) page 114 and applying modulus
f1_orig <- function(theta) {
  return(theta^2 / (2 * pi^2))
}
f1_mod <- function(theta) {
  theta <- (theta + pi) %% (2 * pi) - pi
  return(theta^2 / (2 * pi^2))
}

# f^CvM local deviation
f2 <- function(theta) {
  res <- (1 - log(2 * (1 - cos(theta))) * sqrt(2) / (2 * pi)) / (2 * pi)
  res[is.infinite(res)] <- 0
  return(res)
}

# Simulation
M <- 50
n <- 20
int1_orig <- int1_mod <- int2 <- wat <- numeric(M)
set.seed(1323131)
for (i in 1:M) {

  # Sample
  theta_i <- rnorm(n = n) %% (2 * pi)

  # MJ (2000) (6.3.70) + (6.3.60)
  int1_orig[i] <- integrate(f = function(th) {
    sapply(th, function(theta) (sum(f1_mod(theta + theta_i)) - n / (2 * pi))^2)
  }, lower = 0, upper = 2 * pi, abs.tol = 1e-6, subdivisions = 1e4,
  stop.on.error = TRUE)$value / (2 * pi) / (4 * n)

  # MJ (2000) (6.3.70) + (6.3.60) applying modulus
  int1_mod[i] <- integrate(f = function(th) {
    sapply(th, function(theta) (sum(f1_orig(theta + theta_i)) - n / (2 * pi))^2)
  }, lower = 0, upper = 2 * pi, abs.tol = 1e-6, subdivisions = 1e4,
  stop.on.error = TRUE)$value / (2 * pi) / (4 * n)

  # f^CvM local deviation
  int2[i] <- 0.5 * integrate(f = function(th) {
    sapply(th, function(theta) sum((2 * pi) * f2(theta + theta_i) - 1)^2)
  }, lower = 0, upper = 2 * pi, abs.tol = 1e-6, subdivisions = 1e4,
  stop.on.error = TRUE)$value / (2 * pi * n)

  # Watson statistic
  wat[i] <- sphunif::cir_stat_Watson(Theta = cbind(theta_i))

}

test_that("Watson vs. Sobolev statistic using f^PCvM", {

  expect_equal(wat, int2, tolerance = 1e-4)

})

test_that("Watson vs. Sobolev statistic using f in MJ (2000) page 114", {

  expect_false(isTRUE(all.equal(wat, int1_orig, tolerance = 1e-3)))
  expect_false(isTRUE(all.equal(diff(wat / int1_orig), rep(0, M - 1),
                                tolerance = 1e-3)))
  expect_false(isTRUE(all.equal(wat, int1_mod, tolerance = 1e-3)))
  expect_false(isTRUE(all.equal(diff(wat / int1_mod), rep(0, M - 1),
                                tolerance = 1e-3)))

})

# par(mfrow = c(2, 2))
# plot(wat, int1_orig, main = "MJ (2000) page 114 + (6.3.70) + (6.3.60)")
# plot(wat, int1_mod, main = "With modulus")
# plot(wat, int2, main = "f^CvM local alternative")
# curve(f1_orig, from = 0, to = 2 * pi, n = 1e4, ylim = c(0, 1), main = "f functions")
# curve(f1_mod, add = TRUE, col = 2, n = 1e4)
# curve(f2, add = TRUE, col = 3, n = 1e4)
