
set.seed(14323)
K <- rpois(n = 1, lambda = 9) + 1
x <- seq(0, 3, l = 10)
weights <- abs(rnorm(n = K))
weights0 <- c(K:2, 0)
weights1 <- c(K:2, 1e-15)
dfs <- rpois(n = K, lambda = 4) + 1
ncps <- runif(n = K, min = 0, max = 2)^2
M <- 1e4

test_that("r_wschisq vs. r_wschisq_Cpp and definition", {

  samp1 <- r_wschisq(n = M, weights = weights, dfs = dfs, ncps = ncps)
  samp2 <- r_wschisq_Cpp(n = M, weights = weights, dfs = dfs, ncps = ncps)
  expect_gt(ks.test(samp1, samp2)$p.value, 0.10)
  samp3 <- colSums(weights * t(sapply(1:K, function(i) {
    lambda <- sqrt(ncps[i] / dfs[i])
    rowSums(matrix((rnorm(n = M * dfs[i]) + lambda)^2,
                   nrow = M, ncol = dfs[i]))
  })))
  expect_gt(ks.test(samp1, samp3)$p.value, 0.10)

})

test_that("Wrong method", {

  expect_error(d_wschisq(x = 1, weights = weights, dfs = dfs, method = "A"))
  expect_error(p_wschisq(x = 1, weights = weights, dfs = dfs, method = "A"))
  expect_error(q_wschisq(u = 0.5, weights = weights, dfs = dfs, method = "A"))

})

test_that("Single chi square", {

  expect_equal(d_wschisq(x = 1, weights = 2, dfs = 3, ncps = 1),
               drop(d_chisq(x = 1 / 2, df = 3, ncp = 1) / 2))
  expect_equal(p_wschisq(x = 1, weights = 2, dfs = 3, ncps = 1),
               drop(p_chisq(x = 1 / 2, df = 3, ncp = 1)))
  expect_equal(q_wschisq(u = 0.5, weights = 2, dfs = 3, ncps = 1),
               drop(2 * qchisq(p = 0.5, df = 3, ncp = 1)))

})

test_that("Zero weights method", {

  expect_equal(p_wschisq(x = 2:1, weights = weights0, dfs = dfs, ncps = ncps,
                         method = "HBE"),
               p_wschisq(x = 2:1, weights = weights1, dfs = dfs, ncps = ncps,
                         method = "HBE"),
               tolerance = 1e-6)
  expect_equal(q_wschisq(u = 0.5, weights = weights0, dfs = dfs, ncps = ncps,
                         method = "HBE"),
               q_wschisq(u = 0.5, weights = weights1, dfs = dfs, ncps = ncps,
                         method = "HBE"),
               tolerance = 1e-6)
  expect_equal(d_wschisq(x = 2:1, weights = weights0, dfs = dfs, ncps = ncps,
                         method = "HBE"),
               d_wschisq(x = 2:1, weights = weights1, dfs = dfs, ncps = ncps,
                         method = "HBE"))
  expect_equal({set.seed(1); r_wschisq(n = 5, weights = c(5, 0),
                                       dfs = 2:1, ncps = 2:1)},
               {set.seed(1); r_wschisq(n = 5, weights = c(5, 1e-10),
                                       dfs = 2:1, ncps = 2:1)})

})

test_that("p_wschisq methods", {

  expect_lt(max(abs(p_wschisq(x = x, weights = weights, dfs = dfs,
                              ncps = ncps, method = "HBE") -
                      p_wschisq(x = x, weights = weights, dfs = dfs,
                                ncps = ncps, method = "MC")
  )), 1e-5)
  expect_lt(max(abs(p_wschisq(x = x, weights = weights, dfs = dfs,
                              ncps = ncps, method = "MC") -
                      p_wschisq(x = x, weights = weights, dfs = dfs,
                                ncps = ncps, method = "MC", MC_sample = 
                                  r_wschisq(n = 1e4, weights = weights,
                                            dfs = dfs, ncps = ncps))
  )), 1e-5)
  expect_lt(max(abs(p_wschisq(x = x, weights = weights, dfs = dfs,
                              ncps = ncps, method = "HBE") -
                      p_wschisq(x = x, weights = weights, dfs = dfs,
                                ncps = ncps, method = "SW")
  )), 1e-5)

})
