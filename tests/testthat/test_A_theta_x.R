
set.seed(124132)
th <- seq(0, pi, l = 100)
x <- runif(5, -1, 1)
k <- 0:100

test_that("Reconstruction of A_theta_x for p = 2", {

  for (i in 1:5) {
    expect_equal(drop(A_theta_x(theta = th, x = x[i], p = 2)),
                 drop(Gegen_series(theta = th, coefs =
                                     akx(x = x[i], p = 2, k = k),
                                   k = k, p = 2)),
                 tolerance = 1e-3)
  }

})

test_that("Reconstruction of A_theta_x for p >= 3", {

  for (p in 3:5) {
    for (i in 1:5) {
      expect_equal(drop(A_theta_x(theta = th, x = x[i], p = p)),
                   drop(Gegen_series(theta = th, coefs =
                                       akx(x = x[i], p = p, k = k),
                                     k = k, p = p)),
                   tolerance = 1e-3)
    }
  }

})

test_that("Gegenbauer coefficients of A_theta_x for p >= 2", {

  for (p in 2:5) {
    for (i in 1:5) {
      expect_equal(drop(akx(x = x[i], p = p, k = k)),
                   drop(Gegen_coefs(psi = function(th)
                     drop(A_theta_x(theta = th, x = x[i], p = p)),
                     k = k, p = p)),
                   tolerance = 1e-3)
    }
  }

})
