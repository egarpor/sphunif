
set.seed(1231323)
Theta <- r_unif_cir(n = 10, M = 2)
X <- r_unif_sph(n = 10, p = 2, M = 2)
Theta_s <- sphunif:::sort_each_col(Theta)
n <- 200
samp <- rnorm(n)
samp_s <- sort(samp)
x <- seq(-1, 1, l = 10)
n_dist <- n * (n - 1) / 2
y <- seq(0, 1, l = 20)

test_that("Theta_to_X and X_to_Theta", {

  expect_equal(Theta, X_to_Theta(Theta_to_X(Theta)))
  expect_equal(X, Theta_to_X(X_to_Theta(X)))
  expect_error(Theta_to_X(Theta[, 1]))
  expect_error(X_to_Theta(X[, , 1]))
  expect_error(X_to_Theta(r_unif_sph(n = 10, p = 3, M = 2)))
  expect_error(Theta_to_X(X))
  expect_error(X_to_Theta(Theta))

})

test_that("cir_gaps", {

  expect_equal(cir_gaps(Theta), cbind(cir_gaps(Theta[, 1, drop = FALSE]),
                                      cir_gaps(Theta[, 2, drop = FALSE])))
  expect_equal(cir_gaps(Theta), cir_gaps(Theta_s, sorted = TRUE))
  expect_equal(cir_gaps(Theta), rbind(diff(Theta_s),
                                      2 * pi - (Theta_s[10, ] - Theta_s[1, ])))

})

test_that("ecdf_bin", {

  expect_equal(sphunif:::ecdf_bin(samp, x),
               sphunif:::ecdf_bin(samp_s, x, data_sorted = TRUE))
  expect_equal(sphunif:::ecdf_bin(samp, x),
               sphunif:::ecdf_bin(samp_s, x, data_sorted = TRUE))
  expect_equal(ecdf(samp)(x), drop(sphunif:::ecdf_bin(samp, x, efic = TRUE)))
  expect_equal(n * ecdf(samp)(x),
               drop(sphunif:::ecdf_bin(samp, x, divide_n = FALSE)))
  expect_equal(sphunif:::ecdf_bin(samp, x, efic = FALSE),
               sphunif:::ecdf_bin(samp, x, efic = TRUE))

})

test_that("beta_inc and beta_inc_inv", {

  expect_equal(y,
               drop(sphunif:::beta_inc_inv(sphunif:::beta_inc(y, a = 1, b = 2),
                                           a = 1, b = 2)))
  expect_equal(pmax(x, 0),
               drop(sphunif:::beta_inc_inv(sphunif:::beta_inc(x, a = 1, b = 2),
                                           a = 1, b = 2)))
  expect_equal(y,
               drop(sphunif:::beta_inc(sphunif:::beta_inc_inv(y, a = 1, b = 2),
                                       a = 1, b = 2)))
  expect_equal(y,
               drop(sphunif:::beta_inc(sphunif:::beta_inc_inv(y, a = 1, b = 2,
                                                              lower_tail = TRUE),
                                       a = 1, b = 2, lower_tail = TRUE)))
  expect_equal(exp(sphunif:::beta_inc(y, a = 1, b = 2, log = TRUE)),
               sphunif:::beta_inc(y, a = 1, b = 2))
  expect_equal(sphunif:::beta_inc_inv(log(y), a = 1, b = 2, log = TRUE),
               sphunif:::beta_inc_inv(y, a = 1, b = 2))
  expect_equal(drop(sphunif:::beta_inc(x, 0.75, 2)), pbeta(x, 0.75, 2))
  expect_equal(drop(sphunif:::beta_inc_inv(y, 0.75, 2)), qbeta(y, 0.75, 2))

})

test_that("t_inv_sqrt_one and n_from_dist_vector", {

  expect_equal(drop(sphunif:::t_inv_sqrt_one(x)), x / sqrt(1 - x^2))
  expect_equal(drop(sphunif:::n_from_dist_vector(n_dist)),
               0.5 * (sqrt(8 * n_dist + 1) + 1))

})
