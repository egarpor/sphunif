
# To prevent hanging with parallel computations
Sys.unsetenv("R_TESTS")

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
    names(formals(unif_stat_distr))[-c(1:7, 20:30)] %in%
      names(formals(unif_test))[-c(1:8)]
  ))

})

x <- 1:5

test_that("Errors in edge cases", {

  expect_error(unif_stat_distr(x = c(0, NA)))
  expect_error(unif_stat_distr(x = 1))
  expect_error(unif_stat_distr(x = cbind(1, 3), n = 1, p = 2,
                               type = "Rayleigh"))
  expect_error(unif_stat_distr(x = 1, p = 2, n = 1, type = "Invent"))
  expect_error(unif_stat_distr(x = 1, p = 2, n = 1, type = 1e3))
  expect_error(expect_warning(unif_stat_distr(x = 1, p = 2,
                                              type = function(x) x)))

})

