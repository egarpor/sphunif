
test_that("Statistic arguments of unif_stat_distr are contained in those of
          unif_test", {

  expect_true(all(
    names(formals(unif_stat_distr))[-c(1:7, 20:30)] %in% names(formals(unif_test))[-c(1:8)]
  ))

})

x <- 1:5

test_that("Errors in edge cases", {

  expect_error(unif_stat_distr(x = c(0, NA)))
  expect_error(unif_stat_distr(x = 1))
  expect_error(unif_stat_distr(x = cbind(1, 3), p = 2, type = "Rayleigh"))
  expect_error(unif_stat_distr(x = 1, p = 2, type = "Invent"))
  expect_error(unif_stat_distr(x = 1, p = 2, type = 1e3))
  expect_error(unif_stat_distr(x = 1, p = 2, type = function(x) x))

})