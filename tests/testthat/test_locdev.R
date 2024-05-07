
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

  skip_on_cran()
  for (p in c(2:4, 11)) {
    expect_equal(con_f(f = function(x)
      rotasym::g_vMF(t = x, p = p, kappa = 3, scaled = TRUE),
      p = p, N = 320), 1)
    expect_equal(con_f(f = function(x) f4(x, kappa = 1, q = p - 1),
                       p = p, N = 320), 1)
  }

})

test_that("d_locdev", {

  skip_on_cran()
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

  skip_on_cran()
  for (p in 2:4) {
    mu <- c(rep(0, p - 1), 1)
    samp_1 <- r_locdev(n = 1e3, mu = mu, kappa = 0.25,
                       f = function(z) f4(x = z, kappa = 3, q = p - 1))[, p]
    samp_2 <- F_inv_from_f(f = function(z)
      0.25 * f4(x = z, kappa = 3, q = p - 1) + 0.75 / rotasym::w_p(p = p),
      p = p)(runif(n = 1e3))
    expect_gt(ks.test(samp_1, samp_2)$p.value, 0.01)
    samp_1 <- r_locdev(n = 1e3, mu = mu, kappa = 0, f = NULL)[, 1]
    samp_2 <- r_unif_sph(n = 1e3, p = p)[, 1, 1]
    expect_gt(ks.test(samp_1, samp_2)$p.value, 0.01)
  }

})

test_that("Edge cases d_locdev and r_locdev", {

  skip_on_cran()
  expect_error(d_locdev(x = 1, mu = 1, kappa = -1, f = NULL))
  expect_error(d_locdev(x = 1:2, mu = 1:3, kappa = -1, f = NULL))
  expect_error(r_locdev(n = 1, mu = 1, kappa = -1))

})

test_that("F_from_f via Gauss--Legendre", {

  skip_on_cran()
  for (p in c(2:4, 11)) {
    expect_equal(F_from_f(f = f0, p = p, Gauss = TRUE, K = 1e2)(x),
                 drop(p_proj_unif(x = x, p = p)), tolerance = 1e-3)
  }

})

test_that("F_from_f via integrate()", {

  skip_on_cran()
  for (p in c(2:4, 11)) {
    expect_equal(F_from_f(f = f0, p = p, Gauss = FALSE, K = 1e2)(x),
                 drop(p_proj_unif(x = x, p = p)), tolerance = 1e-3)
  }

})

test_that("F_from_f for vMF", {

  skip_on_cran()
  for (p in c(2:4, 11)) {
    samp_g <- drop(rotasym::r_g_vMF(n = 100, p = p, kappa = 3))
    expect_gt(ks.test(x = F_from_f(f = f1, p = p, Gauss = TRUE,
                                   K = 1e2, kappa = 3)(samp_g),
                      y = "punif")$p.value, 0.01)
  }
  expect_error(F_from_f(f = f1, p = 2, kappa = 1e5))

})

test_that("F_inv_from_f via Gauss--Legendre", {

  skip_on_cran()
  for (p in c(2:4, 11)) {
    expect_equal(F_inv_from_f(f = f0, p = p, Gauss = TRUE, K = 1e2)(u),
                 drop(q_proj_unif(u = u, p = p)), tolerance = 5e-3)
  }

})

test_that("F_inv_from_f via integrate()", {

  skip_on_cran()
  for (p in c(2:4, 11)) {
    expect_equal(F_inv_from_f(f = f0, p = p, Gauss = FALSE, K = 1e2)(u),
                 drop(q_proj_unif(u = u, p = p)), tolerance = 5e-3)
  }

})

test_that("F_inv_from_f for vMF", {

  skip_on_cran()
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
  expect_error(F_inv_from_f(f = f1, p = 2, kappa = 1e5))

})
