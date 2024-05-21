
theta <- seq(0, pi, l = 5)
k <- 0:2
coefs_1 <- 4:2
f_1_2 <- function(th) Gegen_series(theta = th, coefs = coefs_1, k = k, p = 2)
f_1_3 <- function(th) Gegen_series(theta = th, coefs = coefs_1, k = k, p = 3)
f_1_4 <- function(th) Gegen_series(theta = th, coefs = coefs_1, k = k, p = 4)
f_1_11 <- function(th) Gegen_series(theta = th, coefs = coefs_1, k = k, p = 11)
I <- diag(1, nrow = length(k), ncol = length(k))

the_1 <- seq(0, pi, l = 4)
the_2 <- seq(0, pi, l = 7)
m <- 0:1
coefs_2 <- 3:1 %*% t(-1:0)
f_2_2 <- function(th_1, th_2) Gegen_series_2d(theta_1 = th_1, theta_2 = th_2,
                                              coefs = coefs_2, k = k, m = m,
                                              p = 2)
f_2_3 <- function(th_1, th_2) Gegen_series_2d(theta_1 = th_1, theta_2 = th_2,
                                              coefs = coefs_2, k = k, m = m,
                                              p = 3)
f_2_4 <- function(th_1, th_2) Gegen_series_2d(theta_1 = th_1, theta_2 = th_2,
                                              coefs = coefs_2, k = k, m = m,
                                              p = 4)
f_2_11 <- function(th_1, th_2) Gegen_series_2d(theta_1 = th_1, theta_2 = th_2,
                                               coefs = coefs_2, k = k, m = m,
                                               p = 11)

th_k <- drop(Gauss_Legen_nodes(a = 0, b = pi, N = 80))
w_k <- drop(Gauss_Legen_weights(a = 0, b = pi, N = 80))
w_k <- tcrossprod(w_k)
sin_k <- sin(th_k) %o% sin(th_k)

test_that("Gegen_coefs orthonormality", {

  for (p in c(2, 3, 4, 11)) {
    expect_equal(sapply(k, function(l) {
      Gegen_coefs(k = k, p = p, psi = function(th)
        drop(Gegen_polyn(theta = th, k = l, p = p)))
    }), I)
  }

})

test_that("Gegen_coefs normalizing constants", {

  for (p in c(2, 3, 4, 11)) {
    expect_equal(Gegen_coefs(k = k, p = p, only_const = TRUE),
                 sapply(k, function(m) {
                   f <- function(th)
                     drop(Gegen_polyn(theta = th, k = m, p = p))^2 *
                     sin(th)^(p - 2)
                   integrate(f = f, lower = 0, upper = pi)$value
                 }))
  }

})

test_that("Gegen_coefs_2d orthonormality", {

  for (p in c(2, 3, 4, 11)) {
    expect_equal(sapply(m, function(l) {
      sum(Gegen_coefs_2d(k = k, m = m, p = p, psi = function(th_1, th_2)
        Gegen_polyn_2d(theta_1 = th_1, theta_2 = th_2, k = l, m = l, p = p),
        N = 40))
    }), rep(1, length(m)), tolerance = 1e-4)
  }

})

test_that("Gegen_coefs_2d normalizing constants", {

  for (p in c(2, 3, 4, 11)) {
    expect_equal(Gegen_coefs_2d(k = k, m = m, p = p, only_const = TRUE),
                 Gegen_coefs(k = k, p = p, only_const = TRUE) %*%
                   t(Gegen_coefs(k = m, p = p, only_const = TRUE)))
  }

})

test_that("Gegen_series (Gauss = TRUE)", {

  expect_equal(Gegen_coefs(k = k, p = 2, psi = f_1_2), coefs_1)
  expect_equal(Gegen_coefs(k = k, p = 3, psi = f_1_3), coefs_1)
  expect_equal(Gegen_coefs(k = k, p = 4, psi = f_1_4), coefs_1)
  expect_equal(Gegen_coefs(k = k, p = 11, psi = f_1_11), coefs_1)

})

test_that("Gegen_series (Gauss = TRUE) with unnormalized coefficients", {

  expect_equal(Gegen_series(theta = theta, p = 2, k = k, normalize = FALSE,
                            coefs = Gegen_coefs(k = k, p = 2, psi = f_1_2,
                                                normalize = FALSE)),
               f_1_2(theta))
  expect_equal(Gegen_series(theta = theta, p = 3, k = k, normalize = FALSE,
                            coefs = Gegen_coefs(k = k, p = 3, psi = f_1_3,
                                                normalize = FALSE)),
               f_1_3(theta))
  expect_equal(Gegen_series(theta = theta, p = 4, k = k, normalize = FALSE,
                            coefs = Gegen_coefs(k = k, p = 4, psi = f_1_4,
                                                normalize = FALSE)),
               f_1_4(theta))
  expect_equal(Gegen_series(theta = theta, p = 11, k = k, normalize = FALSE,
                            coefs = Gegen_coefs(k = k, p = 11, psi = f_1_11,
                                                normalize = FALSE)),
               f_1_11(theta))

})

test_that("Gegen_series (Gauss = FALSE)", {

  expect_equal(Gegen_coefs(k = k, p = 2, psi = f_1_2, Gauss = FALSE), coefs_1)
  expect_equal(Gegen_coefs(k = k, p = 3, psi = f_1_3, Gauss = FALSE), coefs_1)
  expect_equal(Gegen_coefs(k = k, p = 4, psi = f_1_4, Gauss = FALSE), coefs_1)
  expect_equal(Gegen_coefs(k = k, p = 11, psi = f_1_11, Gauss = FALSE), coefs_1)

})

test_that("Gegen_series_2d (Gauss = TRUE)", {

  expect_equal(Gegen_coefs_2d(k = k, m = m, p = 2, psi = f_2_2, N = 40),
               coefs_2)
  expect_equal(Gegen_coefs_2d(k = k, m = m, p = 3, psi = f_2_3, N = 40),
               coefs_2)
  expect_equal(Gegen_coefs_2d(k = k, m = m, p = 4, psi = f_2_4, N = 40),
               coefs_2)
  expect_equal(Gegen_coefs_2d(k = k, m = m, p = 11, psi = f_2_11, N = 40),
               coefs_2)

})

test_that("Gegen_series_2d (Gauss = FALSE)", {

  expect_equal(Gegen_coefs_2d(k = k, m = m, p = 2, psi = f_2_2,
                              Gauss = FALSE, tol = 1e-2), coefs_2)
  expect_equal(Gegen_coefs_2d(k = k, m = m, p = 3, psi = f_2_3,
                              Gauss = FALSE, tol = 1e-2), coefs_2)
  expect_equal(Gegen_coefs_2d(k = k, m = m, p = 4, psi = f_2_4,
                              Gauss = FALSE, tol = 1e-2), coefs_2)
  expect_equal(Gegen_coefs_2d(k = k, m = m, p = 11, psi = f_2_11,
                              Gauss = FALSE, tol = 1e-2), coefs_2)

})

test_that("Gegen_series_2d (Gauss = TRUE) with unnormalized coefficients", {

  expect_equal(Gegen_series_2d(theta_1 = the_1, theta_2 = the_2, p = 2,
                               k = k, m = m, normalize = FALSE,
                               coefs = Gegen_coefs_2d(k = k, m = m, p = 2,
                                                      psi = f_2_2,
                                                      normalize = FALSE)),
               f_2_2(the_1, the_2))
  expect_equal(Gegen_series_2d(theta_1 = the_1, theta_2 = the_2, p = 3,
                               k = k, m = m, normalize = FALSE,
                               coefs = Gegen_coefs_2d(k = k, m = m, p = 3,
                                                      psi = f_2_3,
                                                      normalize = FALSE)),
               f_2_3(the_1, the_2))
  expect_equal(Gegen_series_2d(theta_1 = the_1, theta_2 = the_2, p = 4,
                               k = k, m = m, normalize = FALSE,
                               coefs = Gegen_coefs_2d(k = k, m = m, p = 4,
                                                      psi = f_2_4,
                                                      normalize = FALSE)),
               f_2_4(the_1, the_2))
  expect_equal(Gegen_series_2d(theta_1 = the_1, theta_2 = the_2, p = 11,
                               k = k, m = m, normalize = FALSE,
                               coefs = Gegen_coefs_2d(k = k, m = m, p = 11,
                                                      psi = f_2_11,
                                                      normalize = FALSE)),
               f_2_11(the_1, the_2))

})

test_that("Gegen_norm", {

  expect_equal(Gegen_norm(coefs = coefs_1, k = k, p = 2)^2,
               integrate(f = function(th) f_1_2(th)^2, lower = 0,
                         upper = pi)$value)
  expect_equal(Gegen_norm(coefs = coefs_1, k = k, p = 3)^2,
               integrate(f = function(th) f_1_3(th)^2 * sin(th), lower = 0,
                         upper = pi)$value)
  expect_equal(Gegen_norm(coefs = coefs_1, k = k, p = 4)^2,
               integrate(f = function(th) f_1_4(th)^2 * sin(th)^2, lower = 0,
                         upper = pi)$value)
  expect_equal(Gegen_norm(coefs = coefs_1, k = k, p = 11)^2,
               integrate(f = function(th) f_1_11(th)^2 * sin(th)^9, lower = 0,
                         upper = pi)$value)

})

test_that("Gegen_norm with unnormalized coefficients", {

  expect_equal(Gegen_norm(coefs = Gegen_coefs(k = k, p = 2, psi = f_1_2,
                                              normalize = FALSE),
                          k = k, p = 2, normalize = FALSE)^2,
               integrate(f = function(th) f_1_2(th)^2, lower = 0,
                         upper = pi)$value)
  expect_equal(Gegen_norm(coefs = Gegen_coefs(k = k, p = 3, psi = f_1_3,
                                              normalize = FALSE),
                          k = k, p = 3, normalize = FALSE)^2,
               integrate(f = function(th) f_1_3(th)^2 * sin(th), lower = 0,
                         upper = pi)$value)
  expect_equal(Gegen_norm(coefs = Gegen_coefs(k = k, p = 4, psi = f_1_4,
                                              normalize = FALSE),
                          k = k, p = 4, normalize = FALSE)^2,
               integrate(f = function(th) f_1_4(th)^2 * sin(th)^2, lower = 0,
                         upper = pi)$value)
  expect_equal(Gegen_norm(coefs = Gegen_coefs(k = k, p = 11, psi = f_1_11,
                                              normalize = FALSE),
                          k = k, p = 11, normalize = FALSE)^2,
               integrate(f = function(th) f_1_11(th)^2 * sin(th)^9, lower = 0,
                         upper = pi)$value)

})

test_that("Gegen_norm_2d", {

  expect_equal(Gegen_norm_2d(coefs = coefs_2, k = k, m = m, p = 2)^2,
               sum(w_k * f_2_2(th_k, th_k)^2))
  expect_equal(Gegen_norm_2d(coefs = coefs_2, k = k, m = m, p = 3)^2,
               sum(w_k * f_2_3(th_k, th_k)^2 * sin_k))
  expect_equal(Gegen_norm_2d(coefs = coefs_2, k = k, m = m, p = 4)^2,
               sum(w_k * f_2_4(th_k, th_k)^2 * sin_k^2))
  expect_equal(Gegen_norm_2d(coefs = coefs_2, k = k, m = m, p = 11)^2,
               sum(w_k * f_2_11(th_k, th_k)^2 * sin_k^9))

})

test_that("Gegen_norm_2d with unnormalized coefficients", {

  expect_equal(Gegen_norm_2d(coefs = Gegen_coefs_2d(k = k, m = m, p = 2,
                                                    psi = f_2_2,
                                                    normalize = FALSE,
                                                    N = 80),
                             k = k, m = m, p = 2, normalize = FALSE)^2,
               sum(w_k * f_2_2(th_k, th_k)^2))
  expect_equal(Gegen_norm_2d(coefs = Gegen_coefs_2d(k = k, m = m, p = 3,
                                                    psi = f_2_3,
                                                    normalize = FALSE,
                                                    N = 80),
                             k = k, m = m, p = 3, normalize = FALSE)^2,
               sum(w_k * f_2_3(th_k, th_k)^2 * sin_k))
  expect_equal(Gegen_norm_2d(coefs = Gegen_coefs_2d(k = k, m = m, p = 4,
                                                    psi = f_2_4,
                                                    normalize = FALSE,
                                                    N = 80),
                             k = k, m = m, p = 4, normalize = FALSE)^2,
               sum(w_k * f_2_4(th_k, th_k)^2 * sin_k^2))
  expect_equal(Gegen_norm_2d(coefs = Gegen_coefs_2d(k = k, m = m, p = 11,
                                                    psi = f_2_11,
                                                    normalize = FALSE,
                                                    N = 80),
                             k = k, m = m, p = 11, normalize = FALSE)^2,
               sum(w_k * f_2_11(th_k, th_k)^2 * sin_k^9))
})

test_that("Bound for Gegenbauer polynomials", {

  x <- seq(-1, 1, l = 200)
  for (p in c(3, 4, 5, 100)) {
    for (k in c(1:3, 10, 50, 100)) {
      expect_true(all(
        abs(Gegen_polyn(theta = acos(x), k = k, p = p)) <
          (2 * (1 - x^2)^(-(p - 2) / 4) * gamma(k + (p - 2) / 2) /
             (gamma((p - 2) / 2) * gamma(k + 1)))
      ))
      expect_true(all(
        abs(Gegen_polyn(theta = acos(x), k = k, p = p)) <
          (2 * (1 - x^2)^(-(p - 2) / 4) * gamma(k + (p - 2) / 2) /
             (gamma((p - 2) / 2) * gamma(k + 1)))
      ))
    }
  }

})

test_that("Gegenbauer polynomials are bounded by their value at 1", {

  x <- seq(-1, 1, l = 200)
  for (p in c(3, 4, 5, 100)) {
    for (k in c(1:3, 10, 50, 100)) {
      expect_true(all(
        abs(Gegen_polyn(theta = acos(x), k = k, p = p)) <=
          drop(Gegen_polyn(theta = acos(1), k = k, p = p))
      ))
    }
  }

})
