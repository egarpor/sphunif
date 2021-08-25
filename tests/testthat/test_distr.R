
set.seed(123456)
x <- cbind(seq(-1, 1, l = 20))
u <- cbind(seq(0, 1, l = 20))
K <- rpois(n = 1, lambda = 9) + 1
weights <- rnorm(n = K)
dfs <- rpois(n = K, lambda = 4) + 1
ncps <- runif(n = K, min = 0, max = 2)^2

test_that("r_proj_unif vs. p_proj_unif", {

  for (p in c(2:4, 11)) {
    expect_gt(ks.test(r_proj_unif(n = 1e3, p = p),
                      "p_proj_unif", p = p)$p.value, 0.05)
  }

})

test_that("r_proj_unif vs. r_unif_sph", {

  for (p in c(2:4, 11)) {
    expect_gt(ks.test(r_proj_unif(n = 1e3, p = p),
                      r_unif_sph(n = 1e3, p = p)[, 1, 1])$p.value, 0.05)
  }

})

test_that("sorting in r_unif_cir", {

  expect_equal({set.seed(1); drop(r_unif_cir(n = 10, sorted = TRUE))},
               {set.seed(1); sort(r_unif_cir(n = 10, sorted = FALSE))})

})

test_that("p_proj_unif vs. q_proj_unif", {

  for (p in c(2:4, 11)) {
    expect_equal(q_proj_unif(u = p_proj_unif(x = x, p = p), p = p), x)
    expect_equal(p_proj_unif(x = q_proj_unif(u = u, p = p), p = p), u)
  }

})

test_that("Wrong dimensions", {

  expect_error(d_proj_unif(x = 0:1, p = 1))
  expect_error(p_proj_unif(x = 0:1, p = 1))
  expect_error(q_proj_unif(u = 0:1, p = 1))
  expect_error(r_proj_unif(n = 1, p = 1))

})

test_that("d_chisq and p_chisq", {

  expect_equal(d_chisq(x = x, df = 3), dchisq(x = x, df = 3))
  expect_equal(p_chisq(x = x, df = 3), pchisq(q = x, df = 3))

})

test_that("r_wschisq_Cpp vs. p_wschisq_MC", {

  samp1 <- r_wschisq_Cpp(n = 1e3, weights = weights, dfs = dfs, ncps = ncps)
  samp2 <- r_wschisq_Cpp(n = 1e4, weights = weights, dfs = dfs, ncps = ncps)
  expect_gt(ks.test(samp1, "p_wschisq_MC", weights = weights, dfs = dfs,
                    ncps = ncps)$p.value, 0.10)
  expect_gt(ks.test(samp1, "p_wschisq_MC", sample = samp2,
                    use_sample = TRUE)$p.value, 0.10)
  expect_gt(ks.test(samp1, "p_wschisq_MC", sample = samp2,
                    use_sample = TRUE, x_sorted = TRUE)$p.value, 0.10)

})
