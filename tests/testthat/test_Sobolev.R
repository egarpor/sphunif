
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
alpha <- c(0.10, 0.05, 0.01)
x <- c(0.1, 0.15, 0.2)
eps <- 1e-9
x_eps1 <- x + eps
x_eps2 <- x - eps

test_that("d_p_k", {

  expect_equal(d_p_k(p = 2, k = 0:3), c(1, rep(2, 3)))
  expect_equal(d_p_k(p = 2, k = 0:3, log = TRUE), c(0, rep(log(2), 3)))
  expect_equal(d_p_k(p = 30, k = 0:1), c(1, 30))

})

## Gegen_coefs() and weights_dfs_Sobolev()

test_that("Gegen_coefs vs. weights_dfs_Sobolev for Ajne", {

  for (p in c(2, 3, 4, 11)) {
    expect_equal(Gegen_coefs(psi = psi_Ajne, k = k, p = p),
                 switch((p == 2) + 1, (1 + 2 * k / (p - 2)), 2) *
                   weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                       type = "Ajne", verbose = FALSE)$weights,
                 tolerance = 1e-6)
  }

})

test_that("Gegen_coefs vs. weights_dfs_Sobolev for PCvM", {

  for (p in c(2:4, 11)) {
    expect_equal(Gegen_coefs(psi = psi_PCvM, k = k, p = p, q = p - 1),
                 switch((p == 2) + 1, (1 + 2 * k / (p - 2)), 2) *
                   weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                       type = "PCvM", verbose = FALSE)$weights,
                 tolerance = 1e-6)
  }

})

test_that("Gegen_coefs vs. weights_dfs_Sobolev for PAD", {

  for (p in c(2:4, 11)) {
    expect_equal(Gegen_coefs(psi = psi_PAD, k = k, p = p, q = p - 1),
                 switch((p == 2) + 1, (1 + 2 * k / (p - 2)), 2) *
                   weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                       type = "PAD", verbose = FALSE)$weights,
                 tolerance = 1e-6)
  }

})

test_that("Gegen_coefs vs. weights_dfs_Sobolev for PRt", {

  for (p in c(2:4, 11)) {
    expect_equal(Gegen_coefs(psi = psi_PRt, k = k, p = p, q = p - 1),
                 switch((p == 2) + 1, (1 + 2 * k / (p - 2)), 2) *
                   weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                       type = "PRt", verbose = FALSE)$weights,
                 tolerance = 1e-5)
  }

})

test_that("Gegen_coefs vs. weights_dfs_Sobolev for Gine_Gn", {

  for (p in c(2:4, 11)) {
    expect_equal(Gegen_coefs(psi = psi_Gine_Gn, k = k, p = p, q = p - 1),
                 switch((p == 2) + 1, (1 + 2 * k / (p - 2)), 2) *
                   weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                       type = "Gine_Gn",
                                       verbose = FALSE)$weights,
                 tolerance = 1e-5)
  }

})

test_that("Gegen_coefs vs. weights_dfs_Sobolev for Gine_Fn", {

  for (p in c(2:4, 11)) {
    expect_equal(Gegen_coefs(psi = psi_Gine_Fn, k = k, p = p, q = p - 1),
                 switch((p == 2) + 1, (1 + 2 * k / (p - 2)), 2) *
                   weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                       type = "Gine_Fn",
                                       verbose = FALSE)$weights,
                 tolerance = 1e-5)
  }

})

test_that("Gegen_coefs vs. weights_dfs_Sobolev for Bakshaev", {

  for (p in c(2:4, 11)) {
    expect_equal(Gegen_coefs(psi = psi_Bakshaev, k = k, p = p),
                 switch((p == 2) + 1, (1 + 2 * k / (p - 2)), 2) *
                   weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                       type = "Bakshaev",
                                       verbose = FALSE)$weights,
                 tolerance = 1e-5)
  }

})

test_that("Gegen_coefs vs. weights_dfs_Sobolev for Riesz", {

  for (p in c(2:4, 11)) {
    expect_equal(Gegen_coefs(psi = psi_Riesz, k = k, p = p, s = 0,
                             N = 5120),
                 switch((p == 2) + 1, (1 + 2 * k / (p - 2)), 2) *
                   weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                       type = "Riesz", Riesz_s = 0,
                                       verbose = FALSE)$weights,
                 tolerance = 1e-5)
    expect_equal(Gegen_coefs(psi = psi_Riesz, k = k, p = p, s = 1),
                 switch((p == 2) + 1, (1 + 2 * k / (p - 2)), 2) *
                   weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                       type = "Riesz", Riesz_s = 1,
                                       verbose = FALSE)$weights,
                 tolerance = 1e-5)
    expect_equal(Gegen_coefs(psi = psi_Riesz, k = k, p = p, s = 2),
                 switch((p == 2) + 1, (1 + 2 * k / (p - 2)), 2) *
                   weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                       type = "Riesz", Riesz_s = 2,
                                       verbose = FALSE)$weights,
                 tolerance = 1e-5)
    expect_equal(Gegen_coefs(psi = psi_Riesz, k = k, p = p, s = 1.5),
                 switch((p == 2) + 1, (1 + 2 * k / (p - 2)), 2) *
                   weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                       type = "Riesz", Riesz_s = 1.5,
                                       verbose = FALSE)$weights,
                 tolerance = 1e-5)
    expect_equal(Gegen_coefs(psi = psi_Riesz, k = k, p = p, s = 0.5),
                 switch((p == 2) + 1, (1 + 2 * k / (p - 2)), 2) *
                   weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                       type = "Riesz", Riesz_s = 0.5,
                                       verbose = FALSE)$weights,
                 tolerance = 1e-5)
    expect_equal(Gegen_coefs(psi = psi_Riesz, k = k, p = p, s = 0.1,
                             N = 5120),
                 switch((p == 2) + 1, (1 + 2 * k / (p - 2)), 2) *
                   weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                       type = "Riesz", Riesz_s = 0.1,
                                       verbose = FALSE)$weights,
                 tolerance = 1e-5)
    expect_equal(Gegen_coefs(psi = psi_Riesz, k = k, p = p, s = -0.5,
                             N = 5120),
                 switch((p == 2) + 1, (1 + 2 * k / (p - 2)), 2) *
                   weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                       type = "Riesz", Riesz_s = -0.5,
                                       verbose = FALSE)$weights,
                 tolerance = 1e-3)
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
               2 * weights_dfs_Sobolev(p = 2, K_max = K %/% 2, thre = 0,
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
    expect_equal(Gegen_coefs_Pn(k = k, p = p, type = "PCvM"),
                 switch((p == 2) + 1, (1 + 2 * k / (p - 2)), 2) *
                   weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                       type = "PCvM", verbose = FALSE)$weights)
  }

})

test_that("Gegen_coefs_Pn vs. weights_dfs_Sobolev for PAD", {

  for (p in c(2, 3, 4, 11)) {
    expect_equal(Gegen_coefs_Pn(k = k, p = p, type = "PAD"),
                 switch((p == 2) + 1, (1 + 2 * k / (p - 2)), 2) *
                   weights_dfs_Sobolev(p = p, K_max = K, thre = 0,
                                       type = "PAD", verbose = FALSE)$weights)
  }

})

test_that("Gegen_coefs_Pn vs. weights_dfs_Sobolev for PRt", {

  for (p in c(2, 3, 4, 11)) {
    expect_equal(Gegen_coefs_Pn(k = k, p = p, type = "PRt"),
                 switch((p == 2) + 1, (1 + 2 * k / (p - 2)), 2) *
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
