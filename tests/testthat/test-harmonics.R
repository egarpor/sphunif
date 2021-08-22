
set.seed(12314)
int_prod <- function(x, i1, i2, k1, k2) {
  g_i_k(x = x, k = k1, i = i1) * g_i_k(x = x, k = k2, i = i2)
}

test_that("angles_to_sphere vs. sphere_to_angles", {

  n <- 25
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

test_that("Harmonics with (i,k) vs. m", {

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

  M <- 1e4
  for (p in 2:5) {
    for (k in 1:2) {

      dpk <- d_p_k(p = p, k = k)
      intprods <- diag(rep(1, dpk))
      for (i1 in 1:dpk) {
        for (i2 in i1:dpk) {
          intprods[i1, i2] <- int_sph_MC(f = int_prod, i1 = i1, i2 = i2,
                                         k1 = k, k2 = k, p = p, M = M,
                                        cores = 1, verbose = FALSE) /
            rotasym::w_p(p = p)
          image(intprods)
          expect_lt(max(abs(intprods - diag(rep(1, dpk)))), 0.1)
        }
      }

    }
  }

})

test_that("Orthonormality of harmonics among i's with different k's = 1,2", {

  M <- 1e4
  for (p in 2:5) {

    dpk1 <- d_p_k(p = p, k = 1)
    dpk2 <- d_p_k(p = p, k = 2)
    intprods <- matrix(0, nrow = dpk1, ncol = dpk2)
    for (i1 in 1:dpk1) {
      for (i2 in 1:dpk2) {
        intprods[i1, i2] <- int_sph_MC(f = int_prod, i1 = i1, i2 = i2,
                                       k1 = 1, k2 = 2, p = p, M = M,
                                       cores = 1, verbose = FALSE) /
          rotasym::w_p(p = p)
        image(intprods)
        expect_lt(max(abs(intprods)), 0.1)
      }
    }

  }

})

test_that("Orthonormality of harmonics among i's with different k's = 0,1", {

  M <- 1e4
  for (p in 2:5) {

    dpk0 <- d_p_k(p = p, k = 0)
    dpk1 <- d_p_k(p = p, k = 1)
    intprods <- matrix(0, nrow = dpk0, ncol = dpk1)
    for (i0 in 1:dpk0) {
      for (i1 in 1:dpk1) {
        intprods[i0, i1] <- int_sph_MC(f = int_prod, i1 = i0, i2 = i1,
                                       k1 = 0, k2 = 1, p = p, M = M,
                                       cores = 1, verbose = FALSE) /
          rotasym::w_p(p = p)
        image(intprods)
        expect_lt(max(abs(intprods)), 0.1)
      }
    }

  }

})
