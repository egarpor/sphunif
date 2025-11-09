
# To prevent hanging with parallel computations
Sys.unsetenv("R_TESTS")

## Correct asymptotic distributions

n <- 1e2
M <- 1e4
x <- c(seq(-5, 10, by = 0.2), n %/% 2 + 1:3)

test_that("Asymptotic distributions available in unif_stat_distr() give a proper
          match to the empirical distributions, for p = 2", {

  set.seed(45673)
  suppressMessages(suppressWarnings({
    distr_as_2 <- unif_stat_distr(x = x, p = 2, n = n, type = "all",
                                  K_max = 1e3, approx = "asymp", method = "HBE")
    distr_mc_2 <- unif_stat_distr(x = x, p = 2, n = n, type = names(distr_as_2),
                                  approx = "MC", cores = 2, chunks = 100)
  }))
  expect_equal(distr_as_2, distr_mc_2, tolerance = 0.05)
  expect_true(all(apply(abs(distr_as_2 - distr_mc_2), 2, max) < 0.1))

})

test_that("Asymptotic distributions available in unif_stat_distr() give a proper
          match to the empirical distributions, for p = 3", {

  set.seed(45673)
  suppressMessages(suppressWarnings({
    distr_as_3 <- unif_stat_distr(x = x, p = 3, n = n, type = "all",
                                  K_max = 1e3, approx = "asymp", method = "HBE")
    distr_mc_3 <- unif_stat_distr(x = x, p = 3, n = n, type = names(distr_as_3),
                                  approx = "MC", cores = 2, chunks = 100)
  }))
  distr_as_3 <- subset(distr_as_3, select = -c(CJ12, Rayleigh_HD))
  distr_mc_3 <- subset(distr_mc_3, select = -c(CJ12, Rayleigh_HD))
  expect_equal(distr_as_3, distr_mc_3, tolerance = 0.05)
  expect_true(all(apply(abs(distr_as_3 - distr_mc_3), 2, max) < 0.1))

})

test_that("Asymptotic distributions available in unif_stat_distr() give a proper
          match to the empirical distributions, for p = 5", {

  set.seed(45673)
  suppressMessages(suppressWarnings({
    distr_as_5 <- unif_stat_distr(x = x, p = 5, n = n, type = "all",
                                  K_max = 1e3, approx = "asymp", method = "HBE")
    distr_mc_5 <- unif_stat_distr(x = x, p = 5, n = n, type = names(distr_as_5),
                                  approx = "MC", cores = 2, chunks = 100)
  }))
  distr_as_5 <- subset(distr_as_5, select = -c(CJ12, Rayleigh_HD))
  distr_mc_5 <- subset(distr_mc_5, select = -c(CJ12, Rayleigh_HD))
  expect_equal(distr_as_5, distr_mc_5, tolerance = 0.05)
  expect_true(all(apply(abs(distr_as_5 - distr_mc_5), 2, max) < 0.1))

})

test_that("Statistic arguments of unif_stat_distr are contained in those of
          unif_test", {

  expect_true(all(
    names(formals(unif_stat_distr))[-c(1:7, 10, 16, 18:20)] %in%
      names(formals(unif_test))[-c(1:8)]
  ))

})

## Errors/warnings in edge cases

x <- 1:5

test_that("Errors in edge cases", {

  expect_error(unif_stat_distr(x = c(0, NA)))
  expect_error(unif_stat_distr(x = 1))
  expect_error(unif_stat_distr(x = cbind(1, 3), n = 1, p = 2,
                               type = "Rayleigh"))
  expect_error(unif_stat_distr(x = 1, p = 2, n = 1, type = "Invent"))
  expect_error(unif_stat_distr(x = 1, p = 2, n = 1, type = 1e3))
  expect_error(unif_stat_distr(x = 1, p = 2, type = function(x) x))
  expect_error(unif_stat_distr(x = 0.5, p = 2, type = "Range"))
  expect_error(unif_stat_distr(x = 0.5, p = 2, type = "Kuiper"))

})

test_that("Warnings in edge cases", {

  expect_warning(unif_stat_distr(x = 0.5, p = 2, type = "Watson"))

})

## Coherency of unif_stat_distr() vs. p_Sobolev()

test_that("Coherency of unif_stat_distr() vs. p_Sobolev()", {

  # Sobolev tests
  type_Sobolev_2 <- c("Watson", "Rothman", "Hermans_Rasson", "Pycke_q")
  type_Sobolev_p <- c("Ajne", "Gine_Gn", "Gine_Fn", "Bakshaev", "Riesz",
                      "PCvM", "PAD", "PRt", "Poisson", "Softmax",
                      "Stein", "Stereo", "Sobolev")

  x <- c(0, 0.1, 0.5, 2, 5)
  tail_asymp_p2 <- unif_stat_distr(x = x, p = 2, n = 0,
                                   type = type_Sobolev_2,
                                   approx = "asymp",
                                   K_max = 1e2, method = "HBE")
  tail_Sobolev_2 <- sapply(type_Sobolev_2, function(t) {
    p_Sobolev(x = x, p = 2, K_max = 1e2, method = "HBE", thre = 0, type = t,
              verbose = FALSE)
  })
  tail_Sobolev_2 <- as.data.frame(tail_Sobolev_2)
  expect_equal(tail_asymp_p2, tail_Sobolev_2, tolerance = 1e-2)

  x <- c(0, 0.1, 0.5, 2, 5)
  tail_asymp_p4 <- unif_stat_distr(x = x, p = 4,
                                   type = type_Sobolev_p,
                                   approx = "asymp",
                                   K_max = 1e2, method = "HBE")
  tail_Sobolev_p4 <- sapply(type_Sobolev_p, function(t) {

    if (t %in% c("Poisson", "Softmax", "Stereo")) {

      v_k2 <- weights_dfs_Sobolev(p = 4, K_max = 1e2, thre = 0, type = t,
                                  verbose = FALSE)$weights
      psi_tilde_0 <- Gegen_series(theta = 0,
                                  coefs = vk2_to_bk(vk2 = v_k2, p = 4),
                                  k = seq_along(v_k2), p = 4)

    } else {

      psi_tilde_0 <- 0

    }

    p_Sobolev(x = x + psi_tilde_0, p = 4, K_max = 1e2, method = "HBE",
              thre = 0, type = t, verbose = FALSE)
  })
  tail_Sobolev_p4 <- as.data.frame(tail_Sobolev_p4)
  expect_equal(tail_asymp_p4, tail_Sobolev_p4, tolerance = 1e-2)

})
