
set.seed(1231323)
Theta <- r_unif_cir(n = 10, M = 2)
X <- r_unif_sph(n = 10, p = 2, M = 2)
Theta_s <- sphunif:::sort_each_col(Theta)
n <- 200
samp <- rnorm(n)
samp_s <- sort(samp)
x <- seq(-1, 1, l = 10)
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

  expect_equal(y, drop(sphunif:::beta_inc_inv(
    sphunif:::beta_inc(y, a = 1, b = 2), a = 1, b = 2)))
  expect_equal(pmax(x, 0), drop(sphunif:::beta_inc_inv(
    sphunif:::beta_inc(x, a = 1, b = 2), a = 1, b = 2)))
  expect_equal(y, drop(sphunif:::beta_inc(
    sphunif:::beta_inc_inv(y, a = 1, b = 2), a = 1, b = 2)))
  expect_equal(y, drop(sphunif:::beta_inc(
    sphunif:::beta_inc_inv(y, a = 1, b = 2, lower_tail = TRUE),
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
  expect_equal(drop(sphunif:::n_from_dist_vector(
    sum(lower.tri(tcrossprod(1:n), diag = FALSE)))), n)

})

## log_besselI_scaled_asymp()

test_that("Correct vectorizations on nu and x", {
  xs <- 1:10
  nus <- c(1:9, 101)
  nus_a <- c(101:102)
  expect_equal(log_besselI_scaled_asymp(nu = nus, x = xs),
               sapply(seq_along(nus), function(i)
                 log_besselI_scaled_asymp(nu = nus[i], x = xs[i])))
  expect_equal(log_besselI_scaled_asymp(nu = nus, x = xs[4]),
               sapply(seq_along(nus), function(i)
                 log_besselI_scaled_asymp(nu = nus[i], x = xs)[4]))
  expect_equal(log_besselI_scaled_asymp(nu = nus[4], x = xs),
               sapply(seq_along(nus), function(i)
                 log_besselI_scaled_asymp(nu = nus, x = xs[i])[4]))
  expect_equal(log_besselI_scaled_asymp(nu = nus_a, x = xs[4]),
               sapply(seq_along(nus_a), function(i)
                 log_besselI_scaled_asymp(nu = nus_a[i], x = xs)[4]))
})

test_that("Correct vectorizations on nu and x for NA's", {
  xs <- c(1, NA, 3:4, NA)
  nus <- 1:5
  for (j in 1:2) {
    expect_equal(log_besselI_scaled_asymp(nu = nus, x = xs),
                 sapply(seq_along(nus), function(i)
                   log_besselI_scaled_asymp(nu = nus[i], x = xs[i])))
    expect_equal(log_besselI_scaled_asymp(nu = nus, x = xs[j]),
                 sapply(seq_along(nus), function(i)
                   log_besselI_scaled_asymp(nu = nus[i], x = xs)[j]))
    expect_equal(log_besselI_scaled_asymp(nu = nus[j], x = xs),
                 sapply(seq_along(nus), function(i)
                   log_besselI_scaled_asymp(nu = nus, x = xs[i])[j]))
  }
})

test_that("Accuracy of log_besselI_scaled_asymp(nu = seq(0, 6, by = 0.5)) with
          asymptotic approximations", {
            x <- seq(1e4, 1e5, l = 100)
            nus <- seq(0, 10, by = 1)
            for (nu in nus) {
              expect_equal(
                log_besselI_scaled_asymp(nu = nu, x = x),
                log_besselI_scaled_asymp(nu = nu, x = x),
                tolerance = 1e-9)
            }
          })

test_that("Asymptotic-kappa Bessel approximation", {
  paper_asymp <- function(x, d) {
    log1p(-d * (d - 2) / (8 * x)) - log(2 * pi * x) / 2
  }
  Bessel_asymp <- function(x, d) {
    Bessel::besselIasym(x = x, nu = (d - 1) / 2, expon.scaled = TRUE,
                        log = TRUE, k.max = 1)
  }
  for (d in 1:10) {
    expect_equal(paper_asymp(x = c(50:100, 1e4, 1e5), d = d),
                 Bessel_asymp(x = c(50:100, 1e4, 1e5), d = d))
  }
})

test_that("Asymptotic-d Bessel approximation", {
  for (x in c(0.1, 1, 10, 100)) {
    expect_lt(
      max(abs(diff(log_besselI_scaled_asymp(nu = 95:105, x = x),
                   differences = 2))),
      max(abs(diff(log_besselI_scaled_asymp(nu = 85:95, x = x),
                   differences = 2)))
    )
  }
})
