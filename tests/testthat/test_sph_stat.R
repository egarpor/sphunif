
n <- 10
set.seed(123456789)
X2 <- r_unif_sph(n = n, p = 2)
X3 <- r_unif_sph(n = n, p = 3)
X4 <- r_unif_sph(n = n, p = 4)
X5 <- r_unif_sph(n = n, p = 5)
X9 <- r_unif_sph(n = n, p = 9)
X200 <- r_unif_sph(n = n, p = 200)
X2_rep <- array(rep(X2, each = 2), dim = c(2 * n, 2, 1))
X3_rep <- array(rep(X3, each = 2), dim = c(2 * n, 3, 1))
X4_rep <- array(rep(X4, each = 2), dim = c(2 * n, 4, 1))
X5_rep <- array(rep(X5, each = 2), dim = c(2 * n, 5, 1))
X9_rep <- array(rep(X9, each = 2), dim = c(2 * n, 9, 1))
X200_rep <- array(rep(X200, each = 2), dim = c(2 * n, 200, 1))
Psi2 <- Psi_mat(X2)
Psi3 <- Psi_mat(X3)
Psi4 <- Psi_mat(X4)
Psi5 <- Psi_mat(X5)
Psi9 <- Psi_mat(X9)
Psi200 <- Psi_mat(X200)
Psi2_rep <- Psi_mat(X2_rep)
Psi3_rep <- Psi_mat(X3_rep)
Psi4_rep <- Psi_mat(X4_rep)
Psi5_rep <- Psi_mat(X5_rep)
Psi9_rep <- Psi_mat(X9_rep)
Psi200_rep <- Psi_mat(X200_rep)
dim(Psi2) <- c(dim(Psi2), 1)
dim(Psi3) <- c(dim(Psi3), 1)
dim(Psi4) <- c(dim(Psi4), 1)
dim(Psi5) <- c(dim(Psi5), 1)
dim(Psi9) <- c(dim(Psi9), 1)
dim(Psi200) <- c(dim(Psi200), 1)
dim(Psi2_rep) <- c(dim(Psi2_rep), 1)
dim(Psi3_rep) <- c(dim(Psi3_rep), 1)
dim(Psi4_rep) <- c(dim(Psi4_rep), 1)
dim(Psi5_rep) <- c(dim(Psi5_rep), 1)
dim(Psi9_rep) <- c(dim(Psi9_rep), 1)
dim(Psi200_rep) <- c(dim(Psi200_rep), 1)
Th2 <- X_to_Theta(X2)
rd3 <- X3[1:5, , 1]

test_that("PAD with X", {

  expect_equal(drop(sph_stat_PAD(X2)), 0.5752646, tolerance = 1e-6)
  expect_equal(drop(sph_stat_PAD(X3)), 0.7093266, tolerance = 1e-6)
  expect_equal(drop(sph_stat_PAD(X4)), 1.4088836, tolerance = 1e-6)
  expect_equal(drop(sph_stat_PAD(X5)), 0.7106791, tolerance = 1e-6)
  expect_equal(drop(sph_stat_PAD(X9)), 1.5825283, tolerance = 1e-6)
  expect_equal(drop(sph_stat_PAD(X200)), 0.9412612, tolerance = 1e-6)

})

test_that("PAD with Psi", {

  expect_equal(drop(sph_stat_PAD(Psi2, Psi_in_X = TRUE, p = 2)), 0.5752646,
               tolerance = 1e-6)
  expect_equal(drop(sph_stat_PAD(Psi3, Psi_in_X = TRUE, p = 3)), 0.7093266,
               tolerance = 1e-6)
  expect_equal(drop(sph_stat_PAD(Psi4, Psi_in_X = TRUE, p = 4)), 1.4088836,
               tolerance = 1e-6)
  expect_equal(drop(sph_stat_PAD(Psi5, Psi_in_X = TRUE, p = 5)), 0.7106791,
               tolerance = 1e-6)
  expect_equal(drop(sph_stat_PAD(Psi9, Psi_in_X = TRUE, p = 9)), 1.5825283,
               tolerance = 1e-6)
  expect_equal(drop(sph_stat_PAD(Psi200, Psi_in_X = TRUE, p = 200)), 0.9412612,
               tolerance = 1e-6)

})

test_that("PCvM with X", {

  expect_equal(drop(sph_stat_PCvM(X2)), 0.0751836, tolerance = 1e-6)
  expect_equal(drop(sph_stat_PCvM(X3)), 0.1027005, tolerance = 1e-6)
  expect_equal(drop(sph_stat_PCvM(X4)), 0.2406067, tolerance = 1e-6)
  expect_equal(drop(sph_stat_PCvM(X5)), 0.1124227, tolerance = 1e-6)
  expect_equal(drop(sph_stat_PCvM(X9)), 0.2766211, tolerance = 1e-6)
  expect_equal(drop(sph_stat_PCvM(X200)), 0.1553891, tolerance = 1e-6)

})

test_that("PCvM with Psi", {

  expect_equal(drop(sph_stat_PCvM(Psi2, Psi_in_X = TRUE, p = 2)), 0.0751836,
               tolerance = 1e-6)
  expect_equal(drop(sph_stat_PCvM(Psi3, Psi_in_X = TRUE, p = 3)), 0.1027005,
               tolerance = 1e-6)
  expect_equal(drop(sph_stat_PCvM(Psi4, Psi_in_X = TRUE, p = 4)), 0.2406067,
               tolerance = 1e-6)
  expect_equal(drop(sph_stat_PCvM(Psi5, Psi_in_X = TRUE, p = 5)), 0.1124227,
               tolerance = 1e-6)
  expect_equal(drop(sph_stat_PCvM(Psi9, Psi_in_X = TRUE, p = 9)), 0.2766211,
               tolerance = 1e-6)
  expect_equal(drop(sph_stat_PCvM(Psi200, Psi_in_X = TRUE, p = 200)), 0.1553891,
               tolerance = 1e-6)

})

test_that("PCvM vs. Bakshaev", {

  expect_equal(sph_stat_PCvM(X3), 0.125 * sph_stat_Bakshaev(X3))
  expect_equal(sph_stat_PCvM(Psi3, Psi_in_X = TRUE, p = 3),
               0.125 * sph_stat_Bakshaev(Psi3, Psi_in_X = TRUE, p = 3))

})

test_that("PRt with X", {

  expect_equal(drop(sph_stat_PRt(X2)), 0.0994200, tolerance = 1e-6)
  expect_equal(drop(sph_stat_PRt(X3)), 0.1188291, tolerance = 1e-6)
  expect_equal(drop(sph_stat_PRt(X4)), 0.3213707, tolerance = 1e-6)
  expect_equal(drop(sph_stat_PRt(X5)), 0.1477722, tolerance = 1e-6)
  expect_equal(drop(sph_stat_PRt(X9)), 0.3794311, tolerance = 1e-6)
  expect_equal(drop(sph_stat_PRt(X200)), 0.2059544, tolerance = 1e-6)

})

test_that("PRt with Psi", {

  expect_equal(drop(sph_stat_PRt(Psi2, Psi_in_X = TRUE, p = 2)), 0.0994200,
               tolerance = 1e-6)
  expect_equal(drop(sph_stat_PRt(Psi3, Psi_in_X = TRUE, p = 3)), 0.1188291,
               tolerance = 1e-6)
  expect_equal(drop(sph_stat_PRt(Psi4, Psi_in_X = TRUE, p = 4)), 0.3213707,
               tolerance = 1e-6)
  expect_equal(drop(sph_stat_PRt(Psi5, Psi_in_X = TRUE, p = 5)), 0.1477722,
               tolerance = 1e-6)
  expect_equal(drop(sph_stat_PRt(Psi9, Psi_in_X = TRUE, p = 9)), 0.3794311,
               tolerance = 1e-6)
  expect_equal(drop(sph_stat_PRt(Psi200, Psi_in_X = TRUE, p = 200)), 0.2059544,
               tolerance = 1e-6)

})

test_that("PRt vs. Ajne", {

  expect_equal(sph_stat_PRt(X2, t = 0.5), sph_stat_Ajne(X2))
  expect_equal(sph_stat_PRt(X3, t = 0.5), sph_stat_Ajne(X3))
  expect_equal(sph_stat_PRt(X4, t = 0.5), sph_stat_Ajne(X4))
  expect_equal(sph_stat_PRt(X5, t = 0.5), sph_stat_Ajne(X5))
  expect_equal(sph_stat_PRt(X9, t = 0.5), sph_stat_Ajne(X9))
  expect_equal(sph_stat_PRt(X200, t = 0.5), sph_stat_Ajne(X200))

})

test_that("Riesz vs. Bakshaev", {

  expect_equal(sph_stat_Bakshaev(X2), sph_stat_Riesz(X2, s = 1),
               tolerance = 1e-5)
  expect_equal(sph_stat_Bakshaev(X3), sph_stat_Riesz(X3, s = 1),
               tolerance = 1e-5)
  expect_equal(sph_stat_Bakshaev(X4), sph_stat_Riesz(X4, s = 1),
               tolerance = 1e-5)
  expect_equal(sph_stat_Bakshaev(X5), sph_stat_Riesz(X5, s = 1),
               tolerance = 1e-5)
  expect_equal(sph_stat_Bakshaev(X9), sph_stat_Riesz(X9, s = 1),
               tolerance = 1e-5)
  expect_equal(sph_stat_Bakshaev(X200), sph_stat_Riesz(X200, s = 1),
               tolerance = 1e-5)

})

test_that("Riesz vs. Pycke", {

  expect_equal((n - 1) / (2 * n) * sph_stat_Pycke(X2),
               sph_stat_Riesz(X2, s = 0))
  expect_equal((log(4) - 1) / 2 + (2 * pi * (n - 1) / n) * sph_stat_Pycke(X3),
               sph_stat_Riesz(X3, s = 0))
  expect_warning(expect_equal(sph_stat_Pycke(X4), sph_stat_Riesz(X4, s = 0)))
  expect_warning(expect_equal(sph_stat_Pycke(X5), sph_stat_Riesz(X5, s = 0)))
  expect_warning(expect_equal(sph_stat_Pycke(X9), sph_stat_Riesz(X9, s = 0)))
  expect_warning(expect_equal(sph_stat_Pycke(X200),
                              sph_stat_Riesz(X200, s = 0)))
  expect_warning(sph_stat_Pycke(X4))
  expect_warning(sph_stat_Pycke(X5))
  expect_warning(sph_stat_Pycke(X9))
  expect_warning(sph_stat_Pycke(X200))

})

test_that("Pycke with data repetitions is computable", {

  expect_warning(expect_true(is.finite(sph_stat_Pycke(X2_rep))))
  expect_warning(expect_true(is.finite(sph_stat_Pycke(X3_rep))))

})

test_that("Riesz with data repetitions is computable", {

  for (s in c(-0.5, 0)) {
    expect_warning(expect_true(is.finite(sph_stat_Riesz(X2_rep, s = s))))
    expect_warning(expect_true(is.finite(sph_stat_Riesz(X3_rep, s = s))))
    expect_warning(expect_true(is.finite(sph_stat_Riesz(X4_rep, s = s))))
    expect_warning(expect_true(is.finite(sph_stat_Riesz(X5_rep, s = s))))
    expect_warning(expect_true(is.finite(sph_stat_Riesz(X9_rep, s = s))))
    expect_warning(expect_true(is.finite(sph_stat_Riesz(X200_rep, s = s))))
  }

})


test_that("sph_stat with psi_in_X = TRUE", {

  expect_equal(sph_stat_Ajne(Psi2, Psi_in_X = TRUE), sph_stat_Ajne(X2))
  expect_equal(sph_stat_Ajne(Psi3, Psi_in_X = TRUE), sph_stat_Ajne(X3))
  expect_equal(sph_stat_Gine_Gn(Psi2, Psi_in_X = TRUE, p = 2),
               sph_stat_Gine_Gn(X2))
  expect_equal(sph_stat_Gine_Gn(Psi3, Psi_in_X = TRUE, p = 3),
               sph_stat_Gine_Gn(X3))
  expect_equal(sph_stat_Gine_Fn(Psi2, Psi_in_X = TRUE, p = 2),
               sph_stat_Gine_Fn(X2))
  expect_equal(sph_stat_Gine_Fn(Psi3, Psi_in_X = TRUE, p = 3),
               sph_stat_Gine_Fn(X3))
  expect_equal(sph_stat_Pycke(cos(Psi2), Psi_in_X = TRUE, p = 2),
               sph_stat_Pycke(X2))
  expect_equal(sph_stat_Pycke(cos(Psi3), Psi_in_X = TRUE, p = 3),
               sph_stat_Pycke(X3))

})

test_that("CJ12", {

  expect_equal(sph_stat_CJ12(X3, regime = 1), sph_stat_CJ12(X3, regime = 2))
  pn <- 3
  expect_equal(drop(sph_stat_CJ12(X3, regime = 1) -
                      sph_stat_CJ12(X3, regime = 3)),
               4 * log(n) - log(log(n)) -
                 (4 * pn * (pn - 2) * log(n) - log(pn)))

})

test_that("CCF09", {

  expect_equal(1 - p_Kolmogorov(sph_stat_CCF09(X3, dirs = rd3)),
               sph_stat_CCF09(X3, dirs = rd3, original = TRUE))
  expect_error(sph_stat_CCF09(X4, dirs = rd3))

})

test_that("Edge cases", {

  expect_error(sph_stat_Gine_Gn(Psi2, Psi_in_X = TRUE))
  expect_error(sph_stat_Gine_Fn(Psi2, Psi_in_X = TRUE))
  expect_error(sph_stat_Pycke(Psi2, Psi_in_X = TRUE))
  expect_error(sph_stat_Bakshaev(Psi2, Psi_in_X = TRUE))
  expect_error(sph_stat_PCvM(Psi2, Psi_in_X = TRUE))
  expect_error(sph_stat_PAD(Psi2, Psi_in_X = TRUE))
  expect_error(sph_stat_PRt(Psi2, Psi_in_X = TRUE))
  expect_error(sph_stat_CJ12(Psi2, Psi_in_X = TRUE))

})

test_that("Edge cases psi", {

  expect_error(sphunif:::sph_stat_PCvM_Psi(cbind(drop(Psi2)), n = n, p = 1,
                                           th_grid = 1:10, int_grid = 1:10))
  expect_error(sphunif:::sph_stat_PAD_Psi(cbind(drop(Psi2)), n = n, p = 1,
                                          th_grid = 1:10, int_grid = 1:10))
  expect_error(sphunif:::sph_stat_PRt_Psi(cbind(drop(Psi2)), n = n, p = 1,
                                          t_m = 0.5, theta_t_m = 1,
                                          th_grid = 1:10, int_grid = 1:10))
  expect_error(expect_warning(
  sphunif:::sph_stat_Pycke_Psi(cbind(drop(Psi2)), n = n, p = 1)))
  expect_error(expect_warning(
    sphunif:::sph_stat_Pycke_Psi(cbind(drop(Psi2)), n = n, p = 5)))

})

# Circular projected statistic with weight w
cir_stat_Pn <- function(samp, w, N = 2560, ...) {

  # Integration nodes and weighs
  x_k <- Gauss_Legen_nodes(a = -1, b = 1, N = N)
  wx_k <- Gauss_Legen_weights(a = -1, b = 1, N = N)
  s <- sort(x_k, index.return = TRUE)
  x_k <- x_k[s$ix]
  wx_k <- wx_k[s$ix]
  th_k <- Gauss_Legen_nodes(a = 0, b = 2 * pi, N = N)
  wth_k <- Gauss_Legen_weights(a = 0, b = 2 * pi, N = N)

  # Integral on x
  int <- function(gamma) {

    Fn <- sphunif:::ecdf_bin(data = drop(samp[, , 1] %*% gamma), sorted_x = x_k)
    Fp <- p_proj_unif(x = x_k, p = 2)
    fp <- d_proj_unif(x = x_k, p = 2)
    sum(wx_k * (Fn - Fp)^2 * w(Fp, ...) * fp)

  }

  # integral on gamma
  n * sum(wth_k * apply(Theta_to_X(th_k)[, , 1], 1, int)) / (2 * pi)

}

# Parameters
n <- 50
set.seed(1223213)
samp <- r_unif_sph(n = n, p = 2)
w_PCvM <- function(x) rep(1, length(x))
w_k <- function(x, k = 1) -2 * k^2 * pi^2 * cos(2 * pi * k * x)
w_k_plus <- function(x, k = 1) pmax(w_k(x, k = k), 0)
w_k_minus <- function(x, k = 1) pmax(-w_k(x, k = k), 0)

test_that("cir_stat_Pn equals PCvM", {

  skip_on_cran()
  expect_equal(cir_stat_Pn(samp, w = w_PCvM), drop(sph_stat_PCvM(samp)),
               tolerance = 1e-4)

})

test_that("cir_stat_Pn equals Rayleigh", {

  skip_on_cran()
  expect_equal(cir_stat_Pn(samp, w = w_k, k = 1),
               drop(sph_stat_Rayleigh(samp) / 2),
               tolerance = 1e-4)
  expect_equal(cir_stat_Pn(samp, w = w_k_plus, k = 1) -
                 cir_stat_Pn(samp, w = w_k_minus, k = 1),
               drop(sph_stat_Rayleigh(samp) / 2),
               tolerance = 1e-4)

})

test_that("cir_stat_Pn equals Bingham", {

  skip_on_cran()
  expect_equal(cir_stat_Pn(samp, w = w_k, k = 2),
               drop(sph_stat_Bingham(samp) / 2),
               tolerance = 1e-4)
  expect_equal(cir_stat_Pn(samp, w = w_k_plus, k = 2) -
                 cir_stat_Pn(samp, w = w_k_minus, k = 2),
               drop(sph_stat_Bingham(samp) / 2), tolerance = 1e-4)

})
