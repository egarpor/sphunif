
set.seed(12314)
int_prod <- function(x, i1, i2, k1, k2) {
  g_i_k(x = x, k = k1, i = i1) * g_i_k(x = x, k = k2, i = i2)
}

test_that("d_p_k for p = 2 and k = 0,1", {

  expect_equal(d_p_k(p = 2, k = 0:3), c(1, rep(2, 3)))
  expect_equal(d_p_k(p = 2, k = 0:3, log = TRUE), c(0, rep(log(2), 3)))
  expect_equal(d_p_k(p = 30, k = 0:1), c(1, 30))

})

test_that("d_p_k = (1 + (2 * k) / (p - 2)) * C_k^{(p - 2) / 2}(1)", {

  k <- 0:5
  for (p in 2:5) {
    expect_equal(d_p_k(p = p, k = k),
                 drop(Gegen_polyn(theta = 0, k = k, p = p)) *
                   switch((p == 2) + 1, 1 + (2 * k) / (p - 2), 2 - (k == 0)))
  }

})

test_that("omega{p - 2} / C_k^{(p - 2) / 2}(1) =
           1 / c_p_k * (1 + (2 * k) / (p - 2))^(-1) * omega_{p - 1}", {

  k <- 0:5
  for (p in 2:5) {
   expect_equal(rotasym::w_p(p = p - 1) /
                  drop(Gegen_polyn(theta = 0, k = k, p = p)),
                1 / drop(Gegen_coefs(only_const = TRUE, p = p, k = k)) *
                switch((p == 2) + 1, (1 + (2 * k) / (p - 2)),
                       2 - (k == 0))^(-1) * rotasym::w_p(p = p))
  }

})

test_that("c_p_k = (1 + (2 * k) / (p - 2))^(-2) * d_p_k *
           omega_{p - 1} / omega{p - 2}", {

  k <- 0:5
  for (p in 2:5) {
    expect_equal(drop(Gegen_coefs(only_const = TRUE, p = p, k = k)),
                 switch((p == 2) + 1, (1 + (2 * k) / (p - 2)),
                        2 - (k == 0))^(-2) * d_p_k(p = p, k = k) *
                   rotasym::w_p(p = p) / rotasym::w_p(p = p - 1))
  }

})

test_that("angles_to_sphere vs. sphere_to_angles", {

  n <- 10
  for (p in c(2:4, 11)) {
    X <- r_unif_sph(n = n, p = p, M = 1)[, , 1]
    th <- t(matrix(runif(n * (p - 1)), nrow = p - 1, ncol = n) *
              c(rep(pi, p - 2), 2 * pi))
    expect_equal(angles_to_sphere(sphere_to_angles(x = X)), X)
    expect_equal(sphere_to_angles(angles_to_sphere(theta = th)), th)
    expect_equal(angles_to_sphere(sphere_to_angles(x = X[1, ])),
                 X[1, , drop = FALSE])
    expect_equal(sphere_to_angles(angles_to_sphere(theta = th[1, ])),
                 th[1, , drop = FALSE])

  }

})

test_that("Harmonics with (i, k) vs. m", {

  skip_on_cran()
  for (p in 2:5) {
    repeat {
      m <- c(rpois(n = p - 1, lambda = 2), sample(c(0, 1), size = 1))
      if (sum(m) > 0) break
    }
    k <- sum(m)
    ms <- as.matrix(do.call(expand.grid, c(rep(list(0:k), p - 1), list(0:1))))
    ms <- ms[rowSums(ms) == k, ]
    i <- which(apply(ms, 1, function(x) all(x == m)))
    X <- r_unif_sph(n = 1, p = p, M = 1)[, , 1]
    expect_message(expect_equal(
      g_i_k(x = X, i = i, k = k, show_m = TRUE), g_i_k(x = X, m = m)))
  }

})

test_that("Orthonormality of harmonics among i's with common k", {

  skip_on_cran()
  set.seed(323231)
  M <- 1e3
  for (p in 2:5) {
    for (k in 1:2) {

      dpk <- d_p_k(p = p, k = k)
      intprods <- diag(rep(1, dpk))
      for (i1 in 1:dpk) {
        for (i2 in i1:dpk) {
          intprods[i1, i2] <- int_sph_MC(f = int_prod, i1 = i1, i2 = i2,
                                         k1 = k, k2 = k, p = p, M = M,
                                         cores = 1) /
            rotasym::w_p(p = p)
          image(abs(intprods), zlim = c(0, 1.5))
          expect_lt(max(abs(intprods - diag(rep(1, dpk)))), 0.1)
        }
      }

    }
  }

})

test_that("Orthonormality of harmonics among i's with different k's = 1,2", {

  skip_on_cran()
  set.seed(323231)
  M <- 1e3
  for (p in 2:5) {

    dpk1 <- d_p_k(p = p, k = 1)
    dpk2 <- d_p_k(p = p, k = 2)
    intprods <- matrix(0, nrow = dpk1, ncol = dpk2)
    for (i1 in 1:dpk1) {
      for (i2 in 1:dpk2) {
        intprods[i1, i2] <- int_sph_MC(f = int_prod, i1 = i1, i2 = i2,
                                       k1 = 1, k2 = 2, p = p, M = M,
                                       cores = 1) /
          rotasym::w_p(p = p)
        image(abs(intprods), zlim = c(0, 1.5))
        expect_lt(max(abs(intprods)), 0.1)
      }
    }

  }

})

test_that("Orthonormality of harmonics among i's with different k's = 0,1", {

  skip_on_cran()
  set.seed(323231)
  M <- 1e3
  for (p in 2:5) {

    dpk0 <- d_p_k(p = p, k = 0)
    dpk1 <- d_p_k(p = p, k = 1)
    intprods <- matrix(0, nrow = dpk0, ncol = dpk1)
    for (i0 in 1:dpk0) {
      for (i1 in 1:dpk1) {
        intprods[i0, i1] <- int_sph_MC(f = int_prod, i1 = i0, i2 = i1,
                                       k1 = 0, k2 = 1, p = p, M = M,
                                       cores = 1) /
          rotasym::w_p(p = p)
        image(abs(intprods), zlim = c(0, 1.5))
        expect_lt(max(abs(intprods)), 0.1)
      }
    }

  }

})

test_that("Specific cases of harmonics for x = e_1", {

  for (p in 2:5) {
    for (k in 1:3) {
      for (i in 1:d_p_k(p = p, k = k)) {
        expect_equal(g_i_k(x = c(1, rep(0, p - 1)), i = i, k = k),
                     (i == 1) * sqrt((gamma(k + p - 2) * (2 * k + p - 2)) /
                                       (gamma(k + 1) * gamma(p - 1))))
      }
    }
  }

})

test_that("Funk-Hecke formula in Dai and Xu (2013)", {

  skip_on_cran()
  set.seed(323231)
  M <- 1e3
  kappa <- 2
  i0 <- 1
  for (p in 2:5) {
    for (k0 in 1:3) {
      mu <- c(1, rep(0, p - 1))
      psi <- function(x) rotasym::g_vMF(t = x, p = p, kappa = kappa)
      expect_equal(
        mean(g_i_k(x = rotasym::r_vMF(n = M, mu = mu, kappa = kappa),
                     i = i0, k = k0)) / rotasym::w_p(p = p),
        rotasym::w_p(p = p - 1) / rotasym::w_p(p = p) *
          integrate(function(t) psi(x = t) *
                      Gegen_polyn(theta = acos(t), k = k0, p = p) *
                      (1 - t^2)^((p - 3) / 2), lower = -1, upper = 1,
                    stop.on.error = FALSE)$value /
          drop(Gegen_polyn(theta = 0, k = k0, p = p)) *
          drop(g_i_k(x = mu, i = i0, k = k0)),
        tolerance = 0.05
      )
    }
  }

})

test_that("Alternative version Funk-Hecke formula", {

  skip_on_cran()
  set.seed(323231)
  M <- 1e3
  kappa <- 2
  i0 <- 1
  for (p in 2:5) {
    for (k0 in 1:3) {
    mu <- c(1, rep(0, p - 1))
    phi <- function(x) rotasym::g_vMF(t = x, p = p, kappa = kappa)
    expect_equal(
      mean(g_i_k(x = rotasym::r_vMF(n = M, mu = mu, kappa = kappa),
                 i = i0, k = k0)) / rotasym::w_p(p = p),
      switch((p == 2) + 1, (1 + (2 * k0) / (p - 2)),
             2 - (k0 == 0))^(-1) *
        drop(Gegen_coefs(k = k0, p = p, psi = function(th) phi(x = cos(th)))) *
        drop(g_i_k(x = mu, i = i0, k = k0)),
      tolerance = 0.05)
    }
  }

})
