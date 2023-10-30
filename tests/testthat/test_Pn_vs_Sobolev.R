
# To prevent hanging with parallel computations
Sys.unsetenv("R_TESTS")

# Projected statistic with weight w
stat_Pn_w <- function(samp, w, circ = TRUE, N = 5120, M = 1e4, seeds = 1:4,
                      ...) {

  # Sample dimensions
  n <- nrow(samp)
  p <- ncol(samp)

  # Integration nodes and weights for x
  x_k <- Gauss_Legen_nodes(a = -1, b = 1, N = N)
  wx_k <- Gauss_Legen_weights(a = -1, b = 1, N = N)
  s <- sort(x_k, index.return = TRUE)
  x_k <- x_k[s$ix]
  wx_k <- wx_k[s$ix]

  # Sample-independent terms of the integrand
  Fp <- drop(p_proj_unif(x = x_k, p = p))
  wFp <- wx_k * w(Fp) * drop(d_proj_unif(x = x_k, p = p))

  # Integral on x
  int <- function(gamma) {

    apply(rbind(gamma), 1, function(gam) {
      Fn <- sphunif:::ecdf_bin(data = drop(samp[, , 1] %*% gam),
                               sorted_x = x_k)
      sum((Fn - Fp)^2 * wFp, na.rm = TRUE)
    })

  }

  # Integration on circle or sphere
  if (ncol(samp) == 2 && circ) {

    th_k <- Gauss_Legen_nodes(a = 0, b = 2 * pi, N = N)
    wth_k <- Gauss_Legen_weights(a = 0, b = 2 * pi, N = N)
    stat <- n * sum(wth_k * apply(Theta_to_X(th_k)[, , 1], 1, int)) / (2 * pi)

  } else {

    # Integral on gamma
    stat <- n * int_sph_MC(f = int, p = p, M = M, cores = 2, chunks = 4,
                           .export = c("samp", "Fp", "wFp"), seeds = seeds,
                           ...) / rotasym::w_p(p = p)

  }
  return(stat)

}

# Sobolev statistic with weights vk2
stat_Sobolev <- function(samp, vk2) {

  n <- nrow(samp)
  p <- ncol(samp)
  k <- seq_along(vk2)
  coefs <- switch((p == 2) + 1, (1 + (2 * k) / (p - 2)), 2) * vk2
  stat <- 2 * sum(Gegen_series(theta = drop(Psi_mat(data = samp)),
                               coefs = coefs, k = k, p = p)) / n
  stat <- stat +
    drop(Gegen_series(theta = 0, coefs = coefs, k = k, p = p))
  return(stat)

}

# w_k's for circular case
w_k <- function(x, k) -4 * k^2 * pi^2 * cos(2 * pi * k * x)
w_k_plus <- function(x, k) pmax(w_k(x, k = k), 0)
w_k_minus <- function(x, k) pmax(-w_k(x, k = k), 0)

# Log-Pochhammer
lpoch <- function(a, n) {

  lgamma(a + n) - lgamma(a)

}

# e_k_ell's
e_k_ell_q <- function(k, ell, q) {

  c2l2q1 <- log(Gegen_coefs(k = 2 * ell, p = 2 * (q + 1) + 1,
                            only_const = TRUE))
  cte <- lpoch(a = q + 1, n = k) - 2 * lfactorial(k) + lchoose(n = k, k = ell)
  num <- lfactorial(ell) + lpoch(a = ell + q + 1, n = k) +
    lpoch(a = (q + 1) / 2, n = ell) + lpoch(a = 1 / 2, n = ell) +
    lpoch(a = 1 / 2, n = k - ell)
  den <- lpoch(a = q + 1 / 2, n = 2 * ell) + lpoch(a = q / 2 + 1, n = ell) +
    lpoch(a = 2 * ell + q + 3 / 2, n = k - ell)
  return(exp(c2l2q1 + cte + num - den))

}

# b^vk's from vk2's
vk2_to_bvk <- function(vk2, p) {

  # c_m_q
  m <- seq_along(vk2)
  ell <- 0:length(m)
  c_m_q <- Gegen_coefs(k = m, p = p, only_const = TRUE)

  # L_m^2
  q <- p - 1
  L_m <- -(q - 1)^2 / (c_m_q * m * (m + q - 1) * (2 * m + q - 1))
  L_m2 <- L_m^2

  # e_k_ell
  e_k_ell <- t(sapply(m, function(k1) {
    e_k_ell_q(k = k1 - 1, ell = ell, q = q)
  }))
  e_k_ell[!is.finite(e_k_ell)] <- 0

  # b^vk's
  return(forwardsolve(l = L_m2 * e_k_ell, x = vk2, k = nrow(e_k_ell)))

}

# wk's from b^vk's
bvk_to_wk <- function(bvk, p, K = length(bvk)) {

  # c_k_q
  M <- length(bvk) - 1
  k <- 1:K
  m <- 0:M
  c_k_q <- Gegen_coefs(k = k, p = p, only_const = TRUE)

  # L_k^2
  q <- p - 1
  L_k <- -(q - 1)^2 / (c_k_q * k * (k + q - 1) * (2 * k + q - 1))
  L_k2 <- L_k^2

  # e_k_m
  e_k_m <- t(sapply(k, function(k1) {
    e_k_ell_q(k = k1 - 1, ell = m, q = q)
  }))
  e_k_m[!is.finite(e_k_m)] <- 0
  return(drop(L_k2 * e_k_m %*% bvk))

}

# Construct w(F_q(x)) for spherical case
w_M <- function(x, p, bvk) {

  q <- p - 1
  x <- drop(q_proj_unif(u = x, p = p))
  Gegen_series(theta = acos(x), coefs = bvk, k = 2 * (seq_along(bvk) - 1),
               p = 2 * (q + 1) + 1) / drop(d_proj_unif(x = x, p = p))

}

## Validation e_k_ell_q

test_that("e_k_ell_q definition", {

  skip_on_cran()
  for (k in c(5, 10, 100)) {
    for (p in 3:10) {
      ell <- 0:k
      expect_true(
        max(abs(e_k_ell_q(k = k, ell = ell, q = p - 1) /
                  Gegen_coefs(k = 2 * ell, p = 2 * p + 1, psi = function(th)
                    drop(Gegen_polyn(theta = th, k = k, p = p + 2))^2,
                    normalize = FALSE, N = 640)) - 1) < 1e-6)
    }
  }

})

test_that("e_k_ell_q / e_k_k_q increasing", {

  skip_on_cran()
  for (k in c(5, 10, 100, 200, 1000)) {
    for (p in 3:20) {
      ell <- 0:k
      expect_true(all(diff(e_k_ell_q(k = k, ell = ell, q = p - 1) /
                             e_k_ell_q(k = k, ell = k, q = p - 1)) > 0))
    }
  }

})

test_that("e_k_ell_q / e_k_k_q bounded by 1/2", {

  skip_on_cran()
  for (k in c(5, 10, 100, 200, 1000)) {
    for (p in 3:20) {
      ell <- 0:k
      expect_true(all(e_k_ell_q(k = k, ell = ell[-(k + 1)], q = p - 1) /
                        e_k_ell_q(k = k, ell = k, q = p - 1) < 0.5))
    }
  }

})

## Circular case

n <- 50
set.seed(42)
samp <- r_unif_sph(n = n, p = 2)

test_that("stat_Sobolev equals Rayleigh and Bingham for p = 2", {

  skip_on_cran()
  expect_equal(stat_Sobolev(samp, vk2 = c(1, rep(0, 99))),
               drop(sph_stat_Rayleigh(samp)),
               tolerance = 1e-3)
  expect_equal(stat_Sobolev(samp, vk2 = c(0, 1, rep(0, 98))),
               drop(sph_stat_Bingham(samp)),
               tolerance = 1e-3)

})

test_that("stat_Pn_w equals Rayleigh for p = 2", {

  skip_on_cran()
  expect_equal(stat_Pn_w(samp, w = function(x) w_k(x, k = 1)),
               drop(sph_stat_Rayleigh(samp)),
               tolerance = 1e-3)
  expect_equal(stat_Pn_w(samp, w = function(x) w_k(x, k = 1), circ = FALSE),
               drop(sph_stat_Rayleigh(samp)),
               tolerance = 1e-1)
  expect_equal(stat_Pn_w(samp, w = function(x) w_k_plus(x, k = 1)) -
                 stat_Pn_w(samp, w = function(x) w_k_minus(x, k = 1), k = 1),
               drop(sph_stat_Rayleigh(samp)),
               tolerance = 1e-3)

})

test_that("stat_Pn_w equals Bingham for p = 2", {

  skip_on_cran()
  expect_equal(stat_Pn_w(samp, w = function(x) w_k(x, k = 2)),
               drop(sph_stat_Bingham(samp)),
               tolerance = 1e-3)
  expect_equal(stat_Pn_w(samp, w = function(x) w_k(x, k = 2), circ = FALSE),
               drop(sph_stat_Bingham(samp)),
               tolerance = 1e-1)
  expect_equal(stat_Pn_w(samp, w = function(x) w_k_plus(x, k = 2)) -
                 stat_Pn_w(samp, w = function(x) w_k_minus(x, k = 2)),
               drop(sph_stat_Bingham(samp)), tolerance = 1e-3)

})

test_that("stat_Pn_w equals v_3 = 1 statistic for p = 2", {

  skip_on_cran()
  expect_equal(stat_Pn_w(samp, w = function(x) w_k(x, k = 3)),
               stat_Sobolev(samp = samp, vk2 = c(0, 0, 1, rep(0, 97))),
               tolerance = 1e-3)

})

## Spherical case

n <- 25
p <- 4
set.seed(42)
samp <- r_unif_sph(n = n, p = p)
K <- 200
vk2_Ray <- c(1, rep(0, K - 1))
bvk_Ray <- vk2_to_bvk(vk2 = vk2_Ray, p = p)
vk2_Bing <- c(0, 1, rep(0, K - 2))
bvk_Bing <- vk2_to_bvk(vk2 = vk2_Bing, p = p)
vk2_Ajne <- weights_dfs_Sobolev(p = p, K_max = K, type = "Ajne", thre = 0,
                                N = 5120, verbose = FALSE)$weights
bvk_Ajne <- vk2_to_bvk(vk2 = vk2_Ajne, p = p)
vk2_Gine_Gn <- weights_dfs_Sobolev(p = p, K_max = K, type = "Gine_Gn", thre = 0,
                                   N = 5120, verbose = FALSE)$weights
bvk_Gine_Gn <- vk2_to_bvk(vk2 = vk2_Gine_Gn, p = p)
vk2_Gine_Fn <- weights_dfs_Sobolev(p = p, K_max = K, type = "Gine_Fn", thre = 0,
                                   N = 5120, verbose = FALSE)$weights
bvk_Gine_Fn <- vk2_to_bvk(vk2 = vk2_Gine_Fn, p = p)
vk2_PCvM <- weights_dfs_Sobolev(p = p, K_max = K, type = "PCvM", thre = 0,
                                N = 5120, verbose = FALSE)$weights
bvk_PCvM <- vk2_to_bvk(vk2 = vk2_PCvM, p = p)
vk2_PAD <- weights_dfs_Sobolev(p = p, K_max = K, type = "PAD", thre = 0,
                               N = 5120, verbose = FALSE)$weights
bvk_PAD <- vk2_to_bvk(vk2 = vk2_PAD, p = p)

test_that("vk2_to_bvk vs. bvk_to_wk vs", {

  skip_on_cran()
  for (i in 1:10) {
    vk2r <- runif(10)
    bvkr <- runif(10)
    pr <- rpois(1, lambda = 5) + 3
    expect_equal(bvk_to_wk(bvk = vk2_to_bvk(vk2 = vk2r, p = pr), p = pr), vk2r)
    expect_equal(vk2_to_bvk(vk2 = bvk_to_wk(bvk = bvkr, p = pr), p = pr), bvkr)
  }

})

test_that("stat_Pn_w equals Rayleigh for p > 2", {

  skip_on_cran()
  expect_equal(stat_Sobolev(samp, vk2 = vk2_Ray),
               drop(sph_stat_Rayleigh(samp)), tolerance = 1e-3)
  expect_equal(stat_Pn_w(samp, w = function(x) w_M(x, p = p, bvk = bvk_Ray)),
               drop(sph_stat_Rayleigh(samp)),
               tolerance = 1e-1)

})

test_that("stat_Pn_w equals Bingham for p > 2", {

  skip_on_cran()
  expect_equal(stat_Sobolev(samp, vk2 = vk2_Bing),
               drop(sph_stat_Bingham(samp)), tolerance = 1e-3)
  expect_equal(stat_Pn_w(samp, w = function(x) w_M(x, p = p, bvk = bvk_Bing)),
               drop(sph_stat_Bingham(samp)),
               tolerance = 1e-1)

})

test_that("stat_Pn_w equals PCvM for p > 2", {

  skip_on_cran()
  expect_equal(stat_Sobolev(samp, vk2 = vk2_PCvM),
               drop(sph_stat_PCvM(samp)), tolerance = 1e-2)
  expect_equal(stat_Pn_w(samp, w = function(x) w_M(x, p = p, bvk = bvk_PCvM)),
               drop(sph_stat_PCvM(samp)),
               tolerance = 1e-1)

})

test_that("stat_Pn_w equals PAD for p > 2", {

  skip_on_cran()
  expect_equal(stat_Sobolev(samp, vk2 = vk2_PAD),
               drop(sph_stat_PAD(samp)), tolerance = 1e-2)
  expect_equal(stat_Pn_w(samp, w = function(x) w_M(x, p = p, bvk = bvk_PAD)),
               drop(sph_stat_PAD(samp)),
               tolerance = 1e-1)

})

test_that("stat_Pn_w equals Ajne for p > 2", {

  skip_on_cran()
  expect_equal(stat_Sobolev(samp, vk2 = vk2_Ajne),
               drop(sph_stat_Ajne(samp)), tolerance = 1e-2)
  expect_equal(stat_Pn_w(samp, w = function(x) w_M(x, p = p, bvk = bvk_Ajne)),
               drop(sph_stat_Ajne(samp)),
               tolerance = 1e-1)

})

test_that("stat_Pn_w equals Gine_Gn for p > 2", {

  skip_on_cran()
  expect_equal(stat_Sobolev(samp, vk2 = vk2_Gine_Gn),
               drop(sph_stat_Gine_Gn(samp)), tolerance = 1e-2)
  expect_equal(stat_Pn_w(samp, w = function(x) w_M(x, p = p,
                                                   bvk = bvk_Gine_Gn)),
               drop(sph_stat_Gine_Gn(samp)),
               tolerance = 1e-1)

})

test_that("stat_Pn_w equals Gine_Fn for p > 2", {

  skip_on_cran()
  expect_equal(stat_Sobolev(samp, vk2 = vk2_Gine_Fn),
               drop(sph_stat_Gine_Fn(samp)), tolerance = 1e-2)
  expect_equal(stat_Pn_w(samp, w = function(x) w_M(x, p = p,
                                                   bvk = bvk_Gine_Fn)),
               drop(sph_stat_Gine_Fn(samp)),
               tolerance = 1e-1)

})
