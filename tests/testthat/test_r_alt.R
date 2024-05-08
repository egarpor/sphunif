
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


test_that("r_alt for rotationally symmetric alternatives", {

  skip_on_cran()
  for (p in 2:4) {

    samp_g <- r_alt(n = 100, p = p, M = 1, kappa = 2, alt = "vMF")[, p, 1]
    expect_gt(ks.test(x = F_from_f(f = f1, p = p, kappa = 2)(samp_g),
                      y = "punif")$p.value, 0.01)

    samp_g <- r_alt(n = 100, p = p, M = 1, kappa = 2, alt = "W")[, p, 1]
    expect_gt(ks.test(x = F_from_f(f = f2, p = p, kappa = 2)(samp_g),
                      y = "punif")$p.value, 0.01)

    samp_g <- r_alt(n = 100, p = p, M = 1, kappa = 2, nu = 0.5,
                    alt = "SC")[, p, 1]
    expect_gt(ks.test(x = F_from_f(f = f3, p = p, kappa = 2, nu = 0.5)(samp_g),
                      y = "punif")$p.value, 0.01)

    samp_g <- r_alt(n = 100, p = p, M = 1, kappa = 2, alt = "C")[, p, 1]
    expect_gt(ks.test(x = F_from_f(f = f4, p = p, kappa = 2, q = p - 1)(samp_g),
                      y = "punif")$p.value, 0.01)

    samp_1 <- r_alt(n = 1e3, p = p, M = 1, kappa = 0, alt = "MvMF")[, p, 1]
    samp_2 <- r_unif_sph(n = 1e3, p = p, M = 1)[, p, 1]
    expect_gt(ks.test(x = samp_1, y = samp_2)$p.value, 0.01)

    samp_1 <- r_alt(n = 1e3, p = p, M = 1, kappa = 0, alt = "MC")[, p, 1]
    samp_2 <- r_unif_sph(n = 1e3, p = p, M = 1)[, p, 1]
    expect_gt(ks.test(x = samp_1, y = samp_2)$p.value, 0.01)

  }

})

test_that("r_alt for non-rotationally symmetric alternatives", {

  skip_on_cran()
  for (p in c(2:4, 11)) {

    samp_1a <- r_alt(n = 1e3, p = p, M = 1, kappa = 2, alt = "MvMF",
                     axial_mix = TRUE)[, p, 1]
    samp_1b <- r_alt(n = 1e3, p = p, M = 1, kappa = 2, alt = "MvMF",
                     axial_mix = FALSE)[, p, 1]
    samp_2b <- c(apply(diag(rep(1, p)), 1, function(mu)
      t(r_alt(n = round(1e3 / p), p = p, alt = "vMF",
              mu = mu, kappa = 2)[, , 1])))
    samp_2b <- matrix(samp_2b, ncol = p, byrow = TRUE)[, p]
    samp_2a <- samp_2b * sample(c(-1, 1), size = length(samp_2b),
                                replace = TRUE)
    expect_gt(ks.test(x = samp_1a, y = samp_2a)$p.value, 0.01)
    expect_gt(ks.test(x = samp_1b, y = samp_2b)$p.value, 0.01)

    samp_1a <- r_alt(n = 1e3, p = p, M = 1, kappa = 2, alt = "MC",
                     axial_mix = TRUE)[, p, 1]
    samp_1b <- r_alt(n = 1e3, p = p, M = 1, kappa = 2, alt = "MC",
                     axial_mix = FALSE)[, p, 1]
    samp_2b <- c(apply(diag(rep(1, p)), 1, function(mu)
      t(r_alt(n = round(1e3 / p), p = p, alt = "C",
              mu = mu, kappa = 2)[, , 1])))
    samp_2b <- matrix(samp_2b, ncol = p, byrow = TRUE)[, p]
    samp_2a <- samp_2b * sample(c(-1, 1), size = length(samp_2b),
                                replace = TRUE)
    expect_gt(ks.test(x = samp_1a, y = samp_2a)$p.value, 0.01)
    expect_gt(ks.test(x = samp_1b, y = samp_2b)$p.value, 0.01)

    samp_1 <- r_alt(n = 1e3, p = p, M = 1, kappa = 1, alt = "ACG")[, p, 1]
    samp_2 <- mvtnorm::rmvnorm(n = 1e3, mean = rep(0, p),
                               sigma = diag(c(rep(1, p - 1), 1 + 1)))
    samp_2 <- samp_2 / sqrt(rowSums(samp_2^2))
    samp_2 <- samp_2[, p]
    expect_gt(ks.test(x = samp_1, y = samp_2)$p.value, 0.01)

  }

})

test_that("Edge cases in r_alt", {

  skip_on_cran()
  for (p in 2:3) {

    expect_length(r_alt(n = 5, p = p, M = 1, alt = "MvMF"), 5 * p)
    expect_equal(dim(r_alt(n = 1, p = p, M = 1, alt = "vMF")), c(1, p, 1))
    expect_equal(dim(r_alt(n = 1, p = p, M = 1, alt = "MvMF")), c(1, p, 1))
    expect_equal(dim(r_alt(n = 1, p = p, M = 1, alt = "SC")), c(1, p, 1))
    expect_equal(dim(r_alt(n = 1, p = p, M = 1, alt = "C")), c(1, p, 1))
    expect_equal(dim(r_alt(n = 1, p = p, M = 1, alt = "W")), c(1, p, 1))
    expect_equal(dim(r_alt(n = 1, p = p, M = 1, alt = "ACG")), c(1, p, 1))
    expect_equal(dim(r_alt(n = 1, p = p, M = 1, alt = "MC")), c(1, p, 1))
    expect_equal(dim(r_alt(n = 1, p = p, M = 1, alt = "UAD")), c(1, p, 1))
    expect_error(r_alt(n = 100, p = p, M = 1, alt = "WC"))
    expect_error(r_alt(n = 100, p = p, M = 1, kappa = -1, alt = "C"))
    expect_error(r_alt(n = 0, p = p, M = 1, alt = "C"))
    expect_no_warning(r_alt(n = 1, p = p, M = 1, alt = "UAD"))
    expect_no_warning(r_alt(n = 2, p = p, M = 1, alt = "UAD"))
    expect_no_warning(r_alt(n = 3, p = p, M = 2, alt = "UAD"))
    expect_no_warning(r_alt(n = 4, p = p, M = 2, alt = "UAD"))

  }

})

