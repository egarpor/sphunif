
p <- 4
f <- function(th) drop(-sin(th / 2) * sin(th)^(p - 2))
I <- integrate(f, lower = 0, upper = pi, rel.tol = 1e-12)$value

test_that("Gauss_Legen for N = 5, 10, 20, ..., 2560, 5120", {

  expect_equal(sum(Gauss_Legen_weights(a = 0, b = pi, N = 5) *
                     f(Gauss_Legen_nodes(a = 0, b = pi, N = 5))), I,
               tolerance = 1e-3)
  expect_equal(sum(Gauss_Legen_weights(a = 0, b = pi, N = 10) *
                     f(Gauss_Legen_nodes(a = 0, b = pi, N = 10))), I,
               tolerance = 1e-12)
  for (N in c(10 * 2^(1:9))) {
    expect_equal(sum(Gauss_Legen_weights(a = 0, b = pi, N = N) *
                       f(Gauss_Legen_nodes(a = 0, b = pi, N = N))), I,
                 tolerance = 3e-14)
  }

})

test_that("Gauss_Legen edge cases", {

  expect_error(Gauss_Legen_weights(a = 0, b = pi, N = 19))
  expect_error(Gauss_Legen_nodes(a = 0, b = pi, N = 19))

})