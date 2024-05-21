
K <- 8
k <- 1:K
psi_Ajne <- function(th) 1 / 2 - th / (2 * pi)
psi_PCvM <- function(th, q) psi_Pn(theta = th, q = q, type = "PCvM")
psi_PAD <- function(th, q) psi_Pn(theta = th, q = q, type = "PAD")
psi_PRt <- function(th, q) psi_Pn(theta = th, q = q, type = "PRt")
psi_Gine_Gn <- function(th, q) {
  -q / 4 * (gamma(q / 2) / gamma((q + 1) / 2))^2 * sin(th)
}
psi_Gine_Fn <- function(th, q) {
  psi_Gine_Gn(th = th, q = q) + 4 * psi_Ajne(th = th)
}
psi_Bakshaev <- function(th) -2 * sin(th / 2)
psi_Riesz <- function(th, s) {
  if (s == 0) {
    return(-0.5 * log(2 * (1 - cos(th))))
  } else {
    return(-sign(s) * 2^s * sin(th / 2)^s)
  }
}
psi_Watson <- function(th) 2 * (th^2 / (4 * pi^2) - th / (2 * pi))
psi_Rothman <- function(th, t = 1 / 3) {
  tm <- min(t, 1 - t)
  pmax(0, tm - th / (2 * pi)) - tm^2
}
psi_Hermans_Rasson <- function(th) {
  beta2 <- (pi^2 / 36) / (0.5 - 4 / pi^2)
  -(th + beta2 * sin(th))
}
psi_Pycke_q <- function(th, q = 0.5) {
  2 * (cos(th) - q) / (1 + q^2 - 2 * q * cos(th))
}
psi_Poisson <- function(th, rho, q) {
  ((1 - rho) / sqrt(1 - 2 * rho * cos(th) + rho^2))^(q + 1)
}
psi_Softmax <- function(th, kappa) {
  exp(kappa * (cos(th) - 1))
}
psi_Stereo <- function(th, a) {
  ta <- tan(th / 2)
  1 / ta + a * ta
}
alpha <- c(0.10, 0.05, 0.01)
x <- c(0.1, 0.15, 0.2)
eps <- 1e-9
x_eps1 <- x + eps
x_eps2 <- x - eps

## weights_dfs_Sobolev()

test_that("weights_dfs_Sobolev returns the same number of coefficients", {

  expect_true(all(lapply(sapply(avail_cir_tests, function(x)
    tryCatch(as.data.frame(weights_dfs_Sobolev(p = 2, K = 4, type = x,
                                               thre = 0, verbose = FALSE,
                                               Sobolev_vk2 = 1:4)),
             error = function(e) rbind(0))), nrow) %in% c(1, 4)))
  expect_true(all(lapply(sapply(avail_sph_tests, function(x)
    tryCatch(as.data.frame(weights_dfs_Sobolev(p = 3, K = 4, type = x,
                                               thre = 0, verbose = FALSE,
                                               Sobolev_vk2 = 1:4)),
             error = function(e) rbind(0))), nrow) %in% c(1, 4)))

})

test_that("Gegen_coefs vs. weights_dfs_Sobolev for Ajne", {

  for (p in c(2, 3, 4, 11)) {
    expect_equal(bk_to_vk2(Gegen_coefs(psi = psi_Ajne, k = k, p = p), p = p),
                 weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                     type = "Ajne", verbose = FALSE)$weights,
                 tolerance = 1e-6)
  }

})

test_that("Gegen_coefs vs. weights_dfs_Sobolev for PCvM", {

  for (p in c(2:4, 11)) {
    expect_equal(bk_to_vk2(Gegen_coefs(psi = psi_PCvM, k = k, p = p, q = p - 1),
                           p = p),
                 weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                     type = "PCvM", verbose = FALSE)$weights,
                 tolerance = 1e-6)
  }

})

test_that("Gegen_coefs vs. weights_dfs_Sobolev for PAD", {

  for (p in c(2:4, 11)) {
    expect_equal(bk_to_vk2(Gegen_coefs(psi = psi_PAD, k = k, p = p, q = p - 1),
                           p = p),
                 weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                     type = "PAD", verbose = FALSE)$weights,
                 tolerance = 1e-6)
  }

})

test_that("Gegen_coefs vs. weights_dfs_Sobolev for PRt", {

  for (p in c(2:4, 11)) {
    expect_equal(bk_to_vk2(Gegen_coefs(psi = psi_PRt, k = k, p = p, q = p - 1),
                           p = p),
                 weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                     type = "PRt", verbose = FALSE)$weights,
                 tolerance = 1e-5)
  }

})

test_that("Gegen_coefs vs. weights_dfs_Sobolev for Gine_Gn", {

  for (p in c(2:4, 11)) {
    expect_equal(bk_to_vk2(Gegen_coefs(psi = psi_Gine_Gn, k = k, p = p,
                                       q = p - 1), p = p),
                 weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                     type = "Gine_Gn", verbose = FALSE)$weights,
                 tolerance = 1e-5)
  }

})

test_that("Gegen_coefs vs. weights_dfs_Sobolev for Gine_Fn", {

  for (p in c(2:4, 11)) {
    expect_equal(bk_to_vk2(Gegen_coefs(psi = psi_Gine_Fn, k = k, p = p,
                                       q = p - 1), p = p),
                 weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                     type = "Gine_Fn", verbose = FALSE)$weights,
                 tolerance = 1e-5)
  }

})

test_that("Gegen_coefs vs. weights_dfs_Sobolev for Bakshaev", {

  for (p in c(2:4, 11)) {
    expect_equal(bk_to_vk2(Gegen_coefs(psi = psi_Bakshaev, k = k, p = p),
                           p = p),
                 weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                     type = "Bakshaev",
                                     verbose = FALSE)$weights,
                 tolerance = 1e-5)
  }

})

test_that("Gegen_coefs vs. weights_dfs_Sobolev for Riesz", {

  for (s in c(0, 1, 2, 1.5, 0.1, -0.5)) {
    for (p in c(2:4, 11)) {
      expect_equal(bk_to_vk2(Gegen_coefs(psi = psi_Riesz, k = k, p = p, s = s,
                                         N = 5120), p = p),
                   weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                       type = "Riesz", Riesz_s = s,
                                       verbose = FALSE)$weights,
                   tolerance = 1e-3)
    }
  }

})

test_that("Gegen_coefs vs. weights_dfs_Sobolev for Poisson", {

  for (p in c(2:4, 11)) {
    for (rho in seq(0.1, 0.9, by = 0.1)) {
      expect_equal(bk_to_vk2(Gegen_coefs(psi = psi_Poisson, k = k, p = p,
                                         rho = rho, q = p - 1), p = p),
                   weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                       type = "Poisson", Poisson_rho = rho,
                                       verbose = FALSE)$weights,
                   tolerance = 1e-5)
    }
  }

})

test_that("Gegen_coefs vs. weights_dfs_Sobolev for Softmax", {

  for (p in c(2:4, 11)) {
    for (kappa in 0:3) {
      expect_equal(bk_to_vk2(Gegen_coefs(psi = psi_Softmax, k = k, p = p,
                                         kappa = kappa), p = p),
                   weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                       type = "Softmax", Softmax_kappa = kappa,
                                       verbose = FALSE)$weights,
                   tolerance = 1e-5)
    }
  }

})

test_that("Gegen_coefs vs. weights_dfs_Sobolev for Stereo", {

  for (p in c(3:4, 11)) {
    for (a in seq(-1, 1, by = 0.5)) {
      expect_equal(bk_to_vk2(Gegen_coefs(psi = psi_Stereo, k = k, p = p, a = a),
                             p = p),
                   weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                       type = "Stereo", Stereo_a = a,
                                       verbose = FALSE)$weights,
                   tolerance = 1e-5)
    }
  }

})

test_that(paste("Gegen_coefs vs. weights_dfs_Sobolev for Watson, Rothman,",
                "Hermans_Rasson, Pycke_q"), {

  expect_equal(Gegen_coefs(psi = psi_Watson, k = k, p = 2),
               2 * weights_dfs_Sobolev(p = 2, K_max = K, thre = 0,
                                       type = "Watson",
                                       verbose = FALSE)$weights,
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_Rothman, k = k, p = 2, N = 640),
               2 * weights_dfs_Sobolev(p = 2, K_max = K, thre = 0,
                                       type = "Rothman",
                                       verbose = FALSE)$weights,
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_Hermans_Rasson, k = k, p = 2),
               2 * weights_dfs_Sobolev(p = 2, K_max = K, thre = 0,
                                       type = "Hermans_Rasson",
                                       verbose = FALSE)$weights,
               tolerance = 1e-6)
  expect_equal(Gegen_coefs(psi = psi_Pycke_q, k = k, p = 2),
               2 * weights_dfs_Sobolev(p = 2, K_max = K, thre = 0,
                                       type = "Pycke_q",
                                       verbose = FALSE)$weights,
               tolerance = 1e-6)
  expect_error(weights_dfs_Sobolev(p = 3, K_max = K, type = "Pycke_q"))

})

test_that("weights_dfs_Sobolev for Watson vs. PCvM with p = 2", {

  expect_equal(0.5 * weights_dfs_Sobolev(p = 2, K_max = K, thre = 0,
                                         type = "Watson",
                                         verbose = FALSE)$weights,
               weights_dfs_Sobolev(p = 2, K_max = K, thre = 0, type = "PCvM",
                                   verbose = FALSE)$weights,
               tolerance = 1e-6)

})

test_that("weights_dfs_Sobolev for Ajne vs. PRt with t = 1 / 2", {

  for (p in c(2, 3, 4, 11)) {
    expect_equal(weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                     type = "Ajne", verbose = FALSE)$weights,
                 weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                     type = "PRt", Rothman_t = 1 / 2,
                                     verbose = FALSE)$weights,
                 tolerance = 1e-6)
  }

})

## Gegen_coefs_Pn()

test_that("Gegen_coefs_Pn vs. weights_dfs_Sobolev for PCvM", {

  for (p in c(2, 3, 4, 11)) {
    expect_equal(bk_to_vk2(Gegen_coefs_Pn(k = k, p = p, type = "PCvM"), p = p),
                 weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                     type = "PCvM", verbose = FALSE)$weights)
  }

})

test_that("Gegen_coefs_Pn vs. weights_dfs_Sobolev for PAD", {

  for (p in c(2, 3, 4, 11)) {
    expect_equal(bk_to_vk2(Gegen_coefs_Pn(k = k, p = p, type = "PAD"), p = p),
                 weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                     type = "PAD", verbose = FALSE)$weights)
  }

})

test_that("Gegen_coefs_Pn vs. weights_dfs_Sobolev for PRt", {

  for (p in c(2, 3, 4, 11)) {
    expect_equal(bk_to_vk2(Gegen_coefs_Pn(k = k, p = p, type = "PRt"), p = p),
                 weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                     type = "PRt", verbose = FALSE)$weights)
  }

})

## [pq]_Sobolev

test_that("p_Sobolev vs. q_Sobolev", {

  for (p in c(2:4, 11)) {
    expect_message(
      expect_equal(p_Sobolev(x = q_Sobolev(u = alpha, p = p, type = "Gine_Fn",
                                           thre = 0, method = "I"),
                             p = p, type = "Gine_Fn", verbose = FALSE,
                             thre = 0, method = "I"),
                   alpha, tolerance = 1e-4)
    )
    expect_message(
      expect_equal(q_Sobolev(u = p_Sobolev(x = x + 1, p = p, type = "Gine_Fn",
                                           thre = 0, method = "I"),
                             p = p, type = "Gine_Fn", verbose = FALSE,
                             thre = 0, method = "I"),
                   x + 1, tolerance = 1e-6)
    )
    expect_message(
      expect_equal(p_Sobolev(x = q_Sobolev(u = alpha, p = p, type = "PCvM",
                                           thre = 0, method = "SW"),
                             p = p, type = "PCvM", verbose = FALSE,
                             thre = 0, method = "SW"),
                   alpha, tolerance = 1e-6)
    )
    expect_message(
      expect_equal(q_Sobolev(u = p_Sobolev(x = x, p = p, type = "PCvM",
                                           thre = 0, method = "SW"),
                             p = p, type = "PCvM", verbose = FALSE,
                             thre = 0, method = "SW"),
                   x, tolerance = 1e-6)
    )
    expect_message(
      expect_equal(p_Sobolev(x = q_Sobolev(u = alpha, p = p, type = "PRt",
                                           thre = 0, method = "HBE"),
                             p = p, type = "PRt", verbose = FALSE,
                             thre = 0, method = "HBE"),
                   alpha, tolerance = 1e-6)
    )
    expect_message(
      expect_equal(q_Sobolev(u = p_Sobolev(x = x, p = p, type = "PRt",
                                           thre = 0, method = "HBE"),
                             p = p, type = "PRt", verbose = FALSE,
                             thre = 0, method = "HBE"),
                   x, tolerance = 1e-6)
    )
  }

})

test_that("p_Sobolev vs. d_Sobolev", {

  for (p in c(2:4, 11)) {
    expect_equal((p_Sobolev(x = x_eps1 + 1, p = p, type = "Gine_Fn",
                            verbose = FALSE, thre = 0, method = "I") -
                    p_Sobolev(x = x_eps2 + 1, p = p, type = "Gine_Fn",
                              verbose = FALSE, thre = 0, method = "I")) /
                   (2 * eps),
                   d_Sobolev(x = x + 1, p = p, type = "Gine_Fn",
                             verbose = FALSE, thre = 0, method = "I"),
                   tolerance = 1e-4)
    expect_equal((p_Sobolev(x = x_eps1, p = p, type = "PCvM", verbose = FALSE,
                            thre = 0, method = "SW") -
                    p_Sobolev(x = x_eps2, p = p, type = "PCvM", verbose = FALSE,
                              thre = 0, method = "SW")) / (2 * eps),
                 d_Sobolev(x = x, p = p, type = "PCvM", verbose = FALSE,
                           thre = 0, method = "SW"),
                 tolerance = 1e-4)
    expect_equal((p_Sobolev(x = x_eps1, p = p, type = "PRt", verbose = FALSE,
                            thre = 0, method = "HBE") -
                    p_Sobolev(x = x_eps2, p = p, type = "PRt", verbose = FALSE,
                              thre = 0, method = "HBE")) / (2 * eps),
                 d_Sobolev(x = x, p = p, type = "PRt", verbose = FALSE,
                           thre = 0, method = "HBE"),
                 tolerance = 1e-4)
  }

})

## sph_stat_Sobolev() and cir_stat_Sobolev()

n <- 10
vk2 <- runif(4)
set.seed(46868)
X2 <- r_unif_sph(n = n, p = 2)
Theta2 <- X_to_Theta(X2)
X3 <- r_unif_sph(n = n, p = 3)
X4 <- r_unif_sph(n = n, p = 4)
X5 <- r_unif_sph(n = n, p = 5)
X9 <- r_unif_sph(n = n, p = 9)
X200 <- r_unif_sph(n = n, p = 200)
Psi2 <- Psi_mat(X2)
Psi3 <- Psi_mat(X3)
Psi4 <- Psi_mat(X4)
Psi5 <- Psi_mat(X5)
Psi9 <- Psi_mat(X9)
Psi200 <- Psi_mat(X200)
dim(Psi2) <- c(dim(Psi2), 1)
dim(Psi3) <- c(dim(Psi3), 1)
dim(Psi4) <- c(dim(Psi4), 1)
dim(Psi5) <- c(dim(Psi5), 1)
dim(Psi9) <- c(dim(Psi9), 1)
dim(Psi200) <- c(dim(Psi200), 1)

test_that("sph_stat_Sobolev for a single and several vk2's", {

  expect_equal(cir_stat_Sobolev(X_to_Theta(X2), vk2 = rbind(vk2, vk2 - 1))[, 1],
               drop(cir_stat_Sobolev(X_to_Theta(X2), vk2 = vk2)))
  expect_equal(sph_stat_Sobolev(X2, vk2 = rbind(vk2, vk2 + 1))[, 1],
               drop(sph_stat_Sobolev(X2, vk2 = vk2)))
  expect_equal(sph_stat_Sobolev(X3, vk2 = rbind(vk2, vk2 + 1))[, 1],
               drop(sph_stat_Sobolev(X3, vk2 = vk2)))
  expect_equal(sph_stat_Sobolev(X4, vk2 = rbind(vk2, vk2 + 1))[, 1],
               drop(sph_stat_Sobolev(X4, vk2 = vk2)))

})

test_that("sph_stat_Sobolev vs. cir_stat_Sobolev", {

  expect_equal(sph_stat_Sobolev(X2, vk2 = vk2),
               cir_stat_Sobolev(X_to_Theta(X2), vk2 = vk2),
               tolerance = 1e-6)
  expect_equal(sph_stat_Sobolev(Psi2, Psi_in_X = TRUE, p = 2, vk2 = vk2),
               cir_stat_Sobolev(Psi2, Psi_in_Theta = TRUE, vk2 = vk2),
               tolerance = 1e-6)

})

test_that("sph_stat_Sobolev with X and Psi", {

  expect_equal(sph_stat_Sobolev(Psi2, Psi_in_X = TRUE, p = 2, vk2 = vk2),
               sph_stat_Sobolev(X2, vk2 = vk2), tolerance = 1e-6)
  expect_equal(sph_stat_Sobolev(Psi3, Psi_in_X = TRUE, p = 3, vk2 = vk2),
               sph_stat_Sobolev(X3, vk2 = vk2), tolerance = 1e-6)
  expect_equal(sph_stat_Sobolev(Psi4, Psi_in_X = TRUE, p = 4, vk2 = vk2),
               sph_stat_Sobolev(X4, vk2 = vk2), tolerance = 1e-6)
  expect_equal(sph_stat_Sobolev(Psi5, Psi_in_X = TRUE, p = 5, vk2 = vk2),
               sph_stat_Sobolev(X5, vk2 = vk2), tolerance = 1e-6)
  expect_equal(sph_stat_Sobolev(Psi9, Psi_in_X = TRUE, p = 9, vk2 = vk2),
               sph_stat_Sobolev(X9, vk2 = vk2), tolerance = 1e-6)
  expect_equal(sph_stat_Sobolev(Psi200, Psi_in_X = TRUE, p = 200, vk2 = vk2),
               sph_stat_Sobolev(X200, vk2 = vk2), tolerance = 1e-6)

})

test_that("sph_stat_Sobolev(vk2 = 1) is a linear form of Rayleigh statistic", {

  for (p in 2:9) {
    stats <- unif_stat_MC(n = 5, type = c("Rayleigh", "Sobolev"), p = p, M = 5,
                          return_stats = TRUE, Sobolev_vk2 = 1)
    expect_equal(drop(cor(stats$stats_MC$Sobolev, stats$stats_MC$Rayleigh)), 1)
  }

})

test_that("sph_stat_Sobolev(vk2 = c(0, 1)) is a linear form of Bingham
          statistic", {

  for (p in 2:9) {
    stats <- unif_stat_MC(n = 5, type = c("Bingham", "Sobolev"), p = p, M = 5,
                          return_stats = TRUE, Sobolev_vk2 = c(0, 1))
    expect_equal(drop(cor(stats$stats_MC$Sobolev, stats$stats_MC$Bingham)), 1)
  }

})

test_that("Edge cases sph_stat_Sobolev()", {

  expect_error(sph_stat_Sobolev(Psi2, vk2 = vk2, Psi_in_X = TRUE))

})

## Coefficients conversion and generating functions

set.seed(12311)
n <- 20
X_2 <- r_unif_sph(n = n, p = 2, M = 1)
X_3 <- r_unif_sph(n = n, p = 3, M = 1)
X_4 <- r_unif_sph(n = n, p = 4, M = 1)
X_11 <- r_unif_sph(n = n, p = 11, M = 1)
th_k <- Gauss_Legen_nodes(a = 0, b = 2 * pi, N = 320)
w_k <- Gauss_Legen_weights(a = 0, b = 2 * pi, N = 320)
z <- seq(-1, 1, l = 1e3)
z8 <- seq(-0.8, 1, l = 1e3)

# PCvM
k_PCvM <- 1:1e3
uk_PCvM_2 <- bk_to_uk(Gegen_coefs_Pn(k = k_PCvM, p = 2, type = "PCvM",
                                     N = 0), p = 2)
f_locdev_PCvM_2 <- function(z) f_locdev(z = z, p = 2, uk = uk_PCvM_2)
integrand_vec_PCvM_2 <- function(x) {
  f_gamma <- matrix(f_locdev_PCvM_2(c(X_2[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 2) * f_gamma - 1)^2 / n / rotasym::w_p(p = 2)
}
f_locdev_PCvM_2_anal <- f_locdev_Pn(p = 2, type = "PCvM")
integrand_vec_PCvM_2_anal <- function(x) {
  f_gamma <- matrix(f_locdev_PCvM_2_anal(c(X_2[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 2) * f_gamma - 1)^2 / n / rotasym::w_p(p = 2)
}
uk_PCvM_3 <- bk_to_uk(Gegen_coefs_Pn(k = k_PCvM, p = 3, type = "PCvM",
                                     N = 0), p = 3)
f_locdev_PCvM_3 <- function(z) f_locdev(z = z, p = 3, uk = uk_PCvM_3)
integrand_vec_PCvM_3 <- function(x) {
  f_gamma <- matrix(f_locdev_PCvM_3(c(X_3[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 3) * f_gamma - 1)^2 / n / rotasym::w_p(p = 3)
}
uk_PCvM_4 <- bk_to_uk(Gegen_coefs_Pn(k = k_PCvM, p = 4, type = "PCvM",
                                     N = 0), p = 4)
f_locdev_PCvM_4 <- function(z) f_locdev(z = z, p = 4, uk = uk_PCvM_4)
integrand_vec_PCvM_4 <- function(x) {
  f_gamma <- matrix(f_locdev_PCvM_4(c(X_4[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 4) * f_gamma - 1)^2 / n / rotasym::w_p(p = 4)
}
uk_PCvM_11 <- bk_to_uk(pmax(Gegen_coefs_Pn(k = k_PCvM, p = 11, type = "PCvM",
                                           N = 5120), 0), p = 11)
f_locdev_PCvM_11 <- function(z) f_locdev(z = z, p = 11, uk = uk_PCvM_11)
integrand_vec_PCvM_11 <- function(x) {
  f_gamma <- matrix(f_locdev_PCvM_11(c(X_11[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 11) * f_gamma - 1)^2 / n / rotasym::w_p(p = 11)
}

# PAD
k_PAD <- 1:1e3
uk_PAD_2 <- bk_to_uk(pmax(Gegen_coefs_Pn(k = k_PAD, p = 2, type = "PAD",
                                         N = 5120), 0), p = 2)
f_locdev_PAD_2 <- function(z) f_locdev(z = z, p = 2, uk = uk_PAD_2)
integrand_vec_PAD_2 <- function(x) {
  f_gamma <- matrix(f_locdev_PAD_2(c(X_2[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 2) * f_gamma - 1)^2 / n / rotasym::w_p(p = 2)
}
uk_PAD_3 <- bk_to_uk(pmax(Gegen_coefs_Pn(k = k_PAD, p = 3, type = "PAD",
                                         N = 5120), 0), p = 3)
f_locdev_PAD_3 <- function(z) f_locdev(z = z, p = 3, uk = uk_PAD_3)
integrand_vec_PAD_3 <- function(x) {
  f_gamma <- matrix(f_locdev_PAD_3(c(X_3[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 3) * f_gamma - 1)^2 / n / rotasym::w_p(p = 3)
}
uk_PAD_4 <- bk_to_uk(pmax(Gegen_coefs_Pn(k = k_PAD, p = 4, type = "PAD",
                                         N = 5120), 0), p = 4)
f_locdev_PAD_4 <- function(z) f_locdev(z = z, p = 4, uk = uk_PAD_4)
integrand_vec_PAD_4 <- function(x) {
  f_gamma <- matrix(f_locdev_PAD_4(c(X_4[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 4) * f_gamma - 1)^2 / n / rotasym::w_p(p = 4)
}
uk_PAD_11 <- bk_to_uk(pmax(Gegen_coefs_Pn(k = k_PAD, p = 11, type = "PAD",
                                          N = 5120), 0), p = 11)
f_locdev_PAD_11 <- function(z) f_locdev(z = z, p = 11, uk = uk_PAD_11)
integrand_vec_PAD_11 <- function(x) {
  f_gamma <- matrix(f_locdev_PAD_11(c(X_11[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 11) * f_gamma - 1)^2 / n / rotasym::w_p(p = 11)
}

# PRt 1 / 3
k_PRt <- 1:1e3
uk_PRt_2 <- bk_to_uk(Gegen_coefs_Pn(k = k_PRt, p = 2, type = "PRt",
                                    N = 0, Rothman_t = 1 / 3),
                     p = 2)
f_locdev_PRt_2 <- function(z) f_locdev(z = z, p = 2, uk = uk_PRt_2)
integrand_vec_PRt_2 <- function(x) {
  f_gamma <- matrix(f_locdev_PRt_2(c(X_2[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 2) * f_gamma - 1)^2 / n / rotasym::w_p(p = 2)
}
f_locdev_PRt_2_anal <- f_locdev_Pn(p = 2, type = "PRt", Rothman_t = 1 / 3)
integrand_vec_PRt_2_anal <- function(x) {
  f_gamma <- matrix(f_locdev_PRt_2_anal(c(X_2[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 2) * f_gamma - 1)^2 / n / rotasym::w_p(p = 2)
}
uk_PRt_3 <- bk_to_uk(Gegen_coefs_Pn(k = k_PRt, p = 3, type = "PRt",
                                    N = 0, Rothman_t = 1 / 3),
                     p = 3)
f_locdev_PRt_3 <- function(z) f_locdev(z = z, p = 3, uk = uk_PRt_3)
integrand_vec_PRt_3 <- function(x) {
  f_gamma <- matrix(f_locdev_PRt_3(c(X_3[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 3) * f_gamma - 1)^2 / n / rotasym::w_p(p = 3)
}
f_locdev_PRt_3_anal <- f_locdev_Pn(p = 3, type = "PRt", Rothman_t = 1 / 3)
integrand_vec_PRt_3_anal <- function(x) {
  f_gamma <- matrix(f_locdev_PRt_3_anal(c(X_3[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 3) * f_gamma - 1)^2 / n / rotasym::w_p(p = 3)
}
uk_PRt_4 <- bk_to_uk(Gegen_coefs_Pn(k = k_PRt, p = 4, type = "PRt",
                                    N = 0, Rothman_t = 1 / 3),
                     p = 4)
f_locdev_PRt_4 <- function(z) f_locdev(z = z, p = 4, uk = uk_PRt_4)
integrand_vec_PRt_4 <- function(x) {
  f_gamma <- matrix(f_locdev_PRt_4(c(X_4[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 4) * f_gamma - 1)^2 / n / rotasym::w_p(p = 4)
}
f_locdev_PRt_4_anal <- f_locdev_Pn(p = 4, type = "PRt", Rothman_t = 1 / 3)
integrand_vec_PRt_4_anal <- function(x) {
  f_gamma <- matrix(f_locdev_PRt_4_anal(c(X_4[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 4) * f_gamma - 1)^2 / n / rotasym::w_p(p = 4)
}
uk_PRt_11 <- bk_to_uk(Gegen_coefs_Pn(k = k_PRt, p = 11, type = "PRt",
                                     N = 0, Rothman_t = 1 / 3),
                      p = 11)
f_locdev_PRt_11 <- function(z) f_locdev(z = z, p = 11, uk = uk_PRt_11)
integrand_vec_PRt_11 <- function(x) {
  f_gamma <- matrix(f_locdev_PRt_11(c(X_11[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 11) * f_gamma - 1)^2 / n / rotasym::w_p(p = 11)
}
f_locdev_PRt_11_anal <- f_locdev_Pn(p = 11, type = "PRt", Rothman_t = 1 / 3)
integrand_vec_PRt_11_anal <- function(x) {
  f_gamma <- matrix(f_locdev_PRt_11_anal(c(X_11[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 11) * f_gamma - 1)^2 / n / rotasym::w_p(p = 11)
}

# Ajne
k_Ajne <- 1:1e3
vk2_Ajne_2 <- weights_dfs_Sobolev(p = 2, K_max = max(k_Ajne),
                                  type = "Ajne", thre = 0,
                                  verbose = FALSE)$weights
uk_Ajne_2 <- vk2_to_uk(vk2 = vk2_Ajne_2, p = 2,
                       signs = (-1)^((seq_along(vk2_Ajne_2) - 1) %/% 2))
f_locdev_Ajne_2 <- function(z) f_locdev(z = z, p = 2, uk = uk_Ajne_2)
integrand_vec_Ajne_2 <- function(x) {
  f_gamma <- matrix(f_locdev_Ajne_2(c(X_2[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 2) * f_gamma - 1)^2 / n / rotasym::w_p(p = 2)
}
vk2_Ajne_3 <- weights_dfs_Sobolev(p = 3, K_max = max(k_Ajne),
                                  type = "Ajne", thre = 0,
                                  verbose = FALSE)$weights
uk_Ajne_3 <- vk2_to_uk(vk2 = vk2_Ajne_3, p = 3)
f_locdev_Ajne_3 <- function(z) f_locdev(z = z, p = 3, uk = uk_Ajne_3)
integrand_vec_Ajne_3 <- function(x) {
  f_gamma <- matrix(f_locdev_Ajne_3(c(X_3[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 3) * f_gamma - 1)^2 / n / rotasym::w_p(p = 3)
}
vk2_Ajne_4 <- weights_dfs_Sobolev(p = 4, K_max = max(k_Ajne),
                                  type = "Ajne", thre = 0,
                                  verbose = FALSE)$weights
uk_Ajne_4 <- vk2_to_uk(vk2 = vk2_Ajne_4, p = 4)
f_locdev_Ajne_4 <- function(z) f_locdev(z = z, p = 4, uk = uk_Ajne_4)
integrand_vec_Ajne_4 <- function(x) {
  f_gamma <- matrix(f_locdev_Ajne_4(c(X_4[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 4) * f_gamma - 1)^2 / n / rotasym::w_p(p = 4)
}
vk2_Ajne_11 <- weights_dfs_Sobolev(p = 11, K_max = max(k_Ajne),
                                   type = "Ajne", thre = 0,
                                   verbose = FALSE)$weights
uk_Ajne_11 <- vk2_to_uk(vk2 = vk2_Ajne_11, p = 11)
f_locdev_Ajne_11 <- function(z) f_locdev(z = z, p = 11, uk = uk_Ajne_11)
integrand_vec_Ajne_11 <- function(x) {
  f_gamma <- matrix(f_locdev_Ajne_11(c(X_11[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 11) * f_gamma - 1)^2 / n / rotasym::w_p(p = 11)
}

# Ajne analytical
k_Ajne <- 1:1e3
f_Ajne_2 <- function(z) (1 * (z >= 0) + 0.5) / rotasym::w_p(p = 2)
integrand_vec_f_Ajne_2  <- function(x) {
  f_gamma <- matrix(f_Ajne_2(c(X_2[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 2) * f_gamma - 1)^2 / n / rotasym::w_p(p = 2)
}
f_Ajne_3 <- function(z) (1 * (z >= 0) + 0.5) / rotasym::w_p(p = 3)
integrand_vec_f_Ajne_3 <- function(x) {
  f_gamma <- matrix(f_Ajne_3(c(X_3[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 3) * f_gamma - 1)^2 / n / rotasym::w_p(p = 3)
}
f_Ajne_4 <- function(z) (1 * (z >= 0) + 0.5) / rotasym::w_p(p = 4)
integrand_vec_f_Ajne_4 <- function(x) {
  f_gamma <- matrix(f_Ajne_4(c(X_4[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 4) * f_gamma - 1)^2 / n / rotasym::w_p(p = 4)
}
f_Ajne_11 <- function(z) (1 * (z >= 0) + 0.5) / rotasym::w_p(p = 11)
integrand_vec_f_Ajne_11 <- function(x) {
  f_gamma <- matrix(f_Ajne_11(c(X_11[, , 1] %*% t(x))),
                    nrow = n, ncol = nrow(x))
  colSums(rotasym::w_p(p = 11) * f_gamma - 1)^2 / n / rotasym::w_p(p = 11)
}

test_that("PCvM as the integral of its local alternative", {

  skip_on_cran()
  expect_equal(drop(sph_stat_PCvM(X = X_2)),
               sum(w_k * integrand_vec_PCvM_2(cbind(cos(th_k), sin(th_k)))),
               tolerance = 5e-2)
  expect_equal(drop(sph_stat_PCvM(X = X_2)),
               integrate(function(th)
                 integrand_vec_PCvM_2_anal(x = cbind(cos(th), sin(th))),
                 lower = 0, upper = 2 * pi, abs.tol = 1e-5,
                 subdivisions = 1e3)$value,
               tolerance = 5e-2)
  expect_equal(drop(sph_stat_PCvM(X = X_3)), {
    int_sph_MC(f = integrand_vec_PCvM_3, p = 3, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16, M = 1e4)
  }, tolerance = 5e-2)
  expect_equal(drop(sph_stat_PCvM(X = X_4)), {
    int_sph_MC(f = integrand_vec_PCvM_4, p = 4, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16, M = 1e4)
  }, tolerance = 5e-2)
  expect_equal(drop(sph_stat_PCvM(X = X_11)), {
    int_sph_MC(f = integrand_vec_PCvM_11, p = 11, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16 + 10, M = 1e4)
  }, tolerance = 1e-1)

})

test_that("PAD as the integral of its local alternative", {

  skip_on_cran()
  expect_equal(drop(sph_stat_PAD(X = X_2)),
               sum(w_k * integrand_vec_PAD_2(cbind(cos(th_k), sin(th_k)))),
               tolerance = 1e-1)
  expect_equal(drop(sph_stat_PAD(X = X_3)), {
    int_sph_MC(f = integrand_vec_PAD_3, p = 3, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16, M = 1e4)
  }, tolerance = 5e-2)
  expect_equal(drop(sph_stat_PAD(X = X_4)), {
    int_sph_MC(f = integrand_vec_PAD_4, p = 4, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16, M = 1e4)
  }, tolerance = 2e-2)
  expect_equal(drop(sph_stat_PAD(X = X_11)), {
    int_sph_MC(f = integrand_vec_PAD_11, p = 11, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16 + 20, M = 1e4)
  }, tolerance = 3e-2)

})

test_that("PRt t = 1 / 3 as the integral of its local alternative", {

  skip_on_cran()
  expect_equal(drop(sph_stat_PRt(X = X_2, t = 1 / 3)),
               sum(w_k * integrand_vec_PRt_2(cbind(cos(th_k), sin(th_k)))),
               tolerance = 5e-2)
  expect_equal(drop(sph_stat_PRt(X = X_2, t = 1 / 3)),
               integrate(function(th)
                 integrand_vec_PRt_2_anal(x = cbind(cos(th), sin(th))),
                 lower = 0, upper = 2 * pi, abs.tol = 1e-5,
                 subdivisions = 1e3)$value,
               tolerance = 1e-4)
  expect_equal(drop(sph_stat_PRt(X = X_3, t = 1 / 3)), {
    int_sph_MC(f = integrand_vec_PRt_3, p = 3, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16, M = 1e4)
  }, tolerance = 5e-2)
  expect_equal(drop(sph_stat_PRt(X = X_3, t = 1 / 3)), {
    int_sph_MC(f = integrand_vec_PRt_3_anal, p = 3, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16, M = 1e4)
  }, tolerance = 3e-2)
  expect_equal(drop(sph_stat_PRt(X = X_4, t = 1 / 3)), {
    int_sph_MC(f = integrand_vec_PRt_4, p = 4, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16, M = 1e4)
  }, tolerance = 5e-2)
  expect_equal(drop(sph_stat_PRt(X = X_4, t = 1 / 3)), {
    int_sph_MC(f = integrand_vec_PRt_4_anal, p = 4, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16, M = 1e4)
  }, tolerance = 3e-2)
  expect_equal(drop(sph_stat_PRt(X = X_11, t = 1 / 3)), {
    int_sph_MC(f = integrand_vec_PRt_11, p = 11, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16 + 10, M = 1e4)
  }, tolerance = 5e-2)
  expect_equal(drop(sph_stat_PRt(X = X_11, t = 1 / 3)), {
    int_sph_MC(f = integrand_vec_PRt_11_anal, p = 11, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16, M = 1e4)
  }, tolerance = 3e-2)

})

test_that("Ajne as the integral of its local alternative (series expansion)", {

  expect_equal(drop(sph_stat_Ajne(X = X_2)),
               sum(w_k * integrand_vec_Ajne_2(cbind(cos(th_k), sin(th_k)))),
               tolerance = 5e-2)
  expect_equal(drop(sph_stat_Ajne(X = X_3)), {
    int_sph_MC(f = integrand_vec_Ajne_3, p = 3, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16, M = 1e4)
  }, tolerance = 5e-2)
  expect_equal(drop(sph_stat_Ajne(X = X_4)), {
    int_sph_MC(f = integrand_vec_Ajne_4, p = 4, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16, M = 1e4)
  }, tolerance = 5e-2)
  expect_equal(drop(sph_stat_Ajne(X = X_11)), {
    int_sph_MC(f = integrand_vec_Ajne_11, p = 11, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16 + 10, M = 1e4)
  }, tolerance = 1e-1)

})

test_that("Ajne as the integral of its local alternative (analytical)", {

  skip_on_cran()
  expect_equal(drop(sph_stat_Ajne(X = X_2)),
               sum(w_k * integrand_vec_f_Ajne_2(cbind(cos(th_k), sin(th_k)))),
               tolerance = 5e-2)
  expect_equal(drop(sph_stat_Ajne(X = X_3)), {
    int_sph_MC(f = integrand_vec_f_Ajne_3, p = 3, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16, M = 1e4)
  }, tolerance = 5e-2)
  expect_equal(drop(sph_stat_Ajne(X = X_4)), {
    int_sph_MC(f = integrand_vec_f_Ajne_4, p = 4, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16, M = 1e4)
  }, tolerance = 5e-2)
  expect_equal(drop(sph_stat_Ajne(X = X_11)), {
    int_sph_MC(f = integrand_vec_f_Ajne_11, p = 11, cores = 2, chunks = 16,
               .export = ls(), seeds = 1:16 + 10, M = 1e4)
  }, tolerance = 3e-2)

})

test_that("Check positivity of f_PCvM and f_PAD", {

  skip_on_cran()
  expect_true(all(f_locdev_PCvM_2(z = z) > 0))
  expect_true(all(f_locdev_PCvM_3(z = z) > 0))
  expect_true(all(f_locdev_PCvM_4(z = z) > 0))
  expect_true(all(f_locdev_PCvM_11(z = z8) > 0))
  expect_true(all(f_locdev_PAD_2(z = z) > 0))
  expect_true(all(f_locdev_PAD_3(z = z) > 0))
  expect_true(all(f_locdev_PAD_4(z = z) > 0))
  expect_true(all(f_locdev_PAD_11(z = z8) > 0))

})

vk2 <- 1:5
bk <- 5:1
uk <- seq(-1, 1, l = 5)

test_that("Conversion between coefficients", {

  skip_on_cran()
  expect_equal(bk_to_vk2(vk2_to_bk(vk2, p = 2), p = 2), vk2)
  expect_equal(bk_to_vk2(vk2_to_bk(vk2, p = 3), p = 3), vk2)
  expect_equal(bk_to_vk2(vk2_to_bk(vk2, p = 2, log = TRUE), p = 2, log = TRUE),
               vk2)
  expect_equal(bk_to_vk2(vk2_to_bk(vk2, p = 3, log = TRUE), p = 3, log = TRUE),
               vk2)
  expect_equal(bk_to_uk(uk_to_bk(uk, p = 2), p = 2, signs = sign(uk)), uk)
  expect_equal(bk_to_uk(uk_to_bk(uk, p = 3), p = 3, signs = sign(uk)), uk)
  expect_equal(vk2_to_bk(bk_to_vk2(bk, p = 2), p = 2), bk)
  expect_equal(vk2_to_bk(bk_to_vk2(bk, p = 3), p = 3), bk)
  expect_equal(vk2_to_bk(bk_to_vk2(bk, p = 2, log = TRUE), p = 2, log = TRUE),
               bk)
  expect_equal(vk2_to_bk(bk_to_vk2(bk, p = 3, log = TRUE), p = 3, log = TRUE),
               bk)
  expect_equal(vk2_to_uk(uk_to_vk2(uk, p = 2), p = 2, signs = sign(uk)), uk)
  expect_equal(vk2_to_uk(uk_to_vk2(uk, p = 3), p = 3, signs = sign(uk)), uk)
  expect_equal(uk_to_vk2(vk2_to_uk(vk2, p = 2), p = 2), vk2)
  expect_equal(uk_to_vk2(vk2_to_uk(vk2, p = 3), p = 3), vk2)
  expect_equal(uk_to_bk(bk_to_uk(bk, p = 2), p = 2), bk)
  expect_equal(uk_to_bk(bk_to_uk(bk, p = 3), p = 3), bk)

})

vk2_mat <- rbind(vk2, vk2 + 1)
log_vk2_mat <- log(vk2_mat)
bk_mat <- rbind(bk, bk + 1)
log_bk_mat <- log(bk_mat)
uk_mat <- rbind(uk + 1, uk + 2)
log_uk_mat <- log(uk_mat)

test_that("Conversion between coefficients, matrix form and log", {

  skip_on_cran()
  expect_equal(unname(vk2_to_bk(log_vk2_mat, p = 2, log = TRUE)),
               log(rbind(vk2_to_bk(vk2_mat[1, ], p = 2),
                         vk2_to_bk(vk2_mat[2, ], p = 2))))
  expect_equal(unname(vk2_to_bk(log_vk2_mat, p = 3, log = TRUE)),
               log(rbind(vk2_to_bk(vk2_mat[1, ], p = 3),
                         vk2_to_bk(vk2_mat[2, ], p = 3))))
  expect_equal(unname(vk2_to_uk(vk2_mat, p = 2)),
               rbind(vk2_to_uk(vk2_mat[1, ], p = 2),
                     vk2_to_uk(vk2_mat[2, ], p = 2)))
  expect_equal(unname(vk2_to_uk(vk2_mat, p = 3)),
               rbind(vk2_to_uk(vk2_mat[1, ], p = 3),
                     vk2_to_uk(vk2_mat[2, ], p = 3)))

  expect_equal(unname(bk_to_vk2(log_bk_mat, p = 2, log = TRUE)),
               log(rbind(bk_to_vk2(bk_mat[1, ], p = 2),
                         bk_to_vk2(bk_mat[2, ], p = 2))))
  expect_equal(unname(bk_to_vk2(log_bk_mat, p = 3, log = TRUE)),
               log(rbind(bk_to_vk2(bk_mat[1, ], p = 3),
                         bk_to_vk2(bk_mat[2, ], p = 3))))
  expect_equal(unname(bk_to_uk(bk_mat, p = 2)),
               rbind(bk_to_uk(bk_mat[1, ], p = 2),
                     bk_to_uk(bk_mat[2, ], p = 2)))
  expect_equal(unname(bk_to_uk(bk_mat, p = 3)),
               rbind(bk_to_uk(bk_mat[1, ], p = 3),
                     bk_to_uk(bk_mat[2, ], p = 3)))

  expect_equal(unname(uk_to_vk2(uk_mat, p = 2)),
               rbind(uk_to_vk2(uk_mat[1, ], p = 2),
                     uk_to_vk2(uk_mat[2, ], p = 2)))
  expect_equal(unname(uk_to_vk2(uk_mat, p = 3)),
               rbind(uk_to_vk2(uk_mat[1, ], p = 3),
                     uk_to_vk2(uk_mat[2, ], p = 3)))
  expect_equal(unname(uk_to_bk(uk_mat, p = 2)),
               rbind(uk_to_bk(uk_mat[1, ], p = 2),
                     uk_to_bk(uk_mat[2, ], p = 2)))
  expect_equal(unname(uk_to_bk(uk_mat, p = 3)),
               rbind(uk_to_bk(uk_mat[1, ], p = 3),
                     uk_to_bk(uk_mat[2, ], p = 3)))

})

test_that("Conversion of bk to vk2 in projected-ecdf statistics", {

  skip_on_cran()
  for (type in c("PCvM", "PAD", "PRt")) {
    for (p in c(2, 3, 4, 9)) {
      expect_equal(bk_to_vk2(Gegen_coefs_Pn(k = 1:9, p = p, type = type),
                             p = p),
                   weights_dfs_Sobolev(K_max = 9, thre = 0, p = p, type = type,
                                       verbose = FALSE)$weights)
    }
  }

})

test_that("Conversion of bk to uk in projected-ecdf statistics", {

  skip_on_cran()
  for (type in c("PCvM", "PAD", "PRt")) {
    for (p in c(2, 3, 4, 9)) {
      expect_equal(bk_to_uk(Gegen_coefs_Pn(k = 1:9, p = p, type = type,
                                           Rothman_t = 0.1), p = p),
                   abs(cutoff_locdev(K_max = 9, thre = 0, p = p, type = type,
                                     Rothman_t = 0.1)))
    }
  }

})

test_that("cutoff_locdev verbose", {

  skip_on_cran()
  expect_message(cutoff_locdev(K_max = 9, thre = 0, p = 9, type = "PAD",
                               verbose = 1))
  suppressMessages(
    expect_message(cutoff_locdev(K_max = 1e1, thre = 1e-4, p = 4, type = "PCvM",
                                 verbose = 2)))

})

# MJ (2000) page 114 and applying modulus
f1_orig <- function(theta) {
  return(theta^2 / (2 * pi^2))
}
f1_mod <- function(theta) {
  theta <- (theta + pi) %% (2 * pi) - pi
  return(theta^2 / (2 * pi^2))
}

# f^CvM local deviation
f2 <- function(theta) {
  res <- (1 - log(2 * (1 - cos(theta))) * sqrt(2) / (2 * pi)) / (2 * pi)
  res[is.infinite(res)] <- 0
  return(res)
}

# Simulation
M <- 50
n <- 20
int1_orig <- int1_mod <- int2 <- wat <- numeric(M)
set.seed(1323131)
for (i in 1:M) {

  # Sample
  theta_i <- rnorm(n = n) %% (2 * pi)

  # MJ (2000) (6.3.70) + (6.3.60)
  int1_orig[i] <- integrate(f = function(th) {
    sapply(th, function(theta) (sum(f1_mod(theta + theta_i)) - n / (2 * pi))^2)
  }, lower = 0, upper = 2 * pi, abs.tol = 1e-6, subdivisions = 1e4,
  stop.on.error = TRUE)$value / (2 * pi) / (4 * n)

  # MJ (2000) (6.3.70) + (6.3.60) applying modulus
  int1_mod[i] <- integrate(f = function(th) {
    sapply(th, function(theta) (sum(f1_orig(theta + theta_i)) - n / (2 * pi))^2)
  }, lower = 0, upper = 2 * pi, abs.tol = 1e-6, subdivisions = 1e4,
  stop.on.error = TRUE)$value / (2 * pi) / (4 * n)

  # f^CvM local deviation
  int2[i] <- 0.5 * integrate(f = function(th) {
    sapply(th, function(theta) sum((2 * pi) * f2(theta + theta_i) - 1)^2)
  }, lower = 0, upper = 2 * pi, abs.tol = 1e-6, subdivisions = 1e4,
  stop.on.error = TRUE)$value / (2 * pi * n)

  # Watson statistic
  wat[i] <- sphunif::cir_stat_Watson(Theta = cbind(theta_i))

}

test_that("Watson vs. Sobolev statistic using f^PCvM", {

  expect_equal(wat, int2, tolerance = 1e-4)

})

test_that("Watson vs. Sobolev statistic using f in MJ (2000) page 114", {

  skip_on_cran()
  expect_false(isTRUE(all.equal(wat, int1_orig, tolerance = 1e-3)))
  expect_false(isTRUE(all.equal(diff(wat / int1_orig), rep(0, M - 1),
                                tolerance = 1e-3)))
  expect_false(isTRUE(all.equal(wat, int1_mod, tolerance = 1e-3)))
  expect_false(isTRUE(all.equal(diff(wat / int1_mod), rep(0, M - 1),
                                tolerance = 1e-3)))

})

# par(mfrow = c(2, 2))
# plot(wat, int1_orig, main = "MJ (2000) page 114 + (6.3.70) + (6.3.60)")
# plot(wat, int1_mod, main = "With modulus")
# plot(wat, int2, main = "f^CvM local alternative")
# curve(f1_orig, from = 0, to = 2 * pi, n = 1e4, ylim = c(0, 1),
#        main = "f functions")
# curve(f1_mod, add = TRUE, col = 2, n = 1e4)
# curve(f2, add = TRUE, col = 3, n = 1e4)
