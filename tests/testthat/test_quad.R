
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
  expect_equal(sum(Gauss_Legen_weights(a = 0, b = pi, N = 20) *
                     f(Gauss_Legen_nodes(a = 0, b = pi, N = 20))), I,
               tolerance = 1e-14)
  expect_equal(sum(Gauss_Legen_weights(a = 0, b = pi, N = 40) *
                     f(Gauss_Legen_nodes(a = 0, b = pi, N = 40))), I,
               tolerance = 1e-14)
  expect_equal(sum(Gauss_Legen_weights(a = 0, b = pi, N = 80) *
                     f(Gauss_Legen_nodes(a = 0, b = pi, N = 80))), I,
               tolerance = 1e-14)
  expect_equal(sum(Gauss_Legen_weights(a = 0, b = pi, N = 160) *
                     f(Gauss_Legen_nodes(a = 0, b = pi, N = 160))), I,
               tolerance = 1e-14)
  expect_equal(sum(Gauss_Legen_weights(a = 0, b = pi, N = 320) *
                     f(Gauss_Legen_nodes(a = 0, b = pi, N = 320))), I,
               tolerance = 1e-14)
  expect_equal(sum(Gauss_Legen_weights(a = 0, b = pi, N = 640) *
                     f(Gauss_Legen_nodes(a = 0, b = pi, N = 640))), I,
               tolerance = 1e-14)
  expect_equal(sum(Gauss_Legen_weights(a = 0, b = pi, N = 640) *
                     f(Gauss_Legen_nodes(a = 0, b = pi, N = 640))), I,
               tolerance = 1e-14)
  expect_equal(sum(Gauss_Legen_weights(a = 0, b = pi, N = 1280) *
                     f(Gauss_Legen_nodes(a = 0, b = pi, N = 1280))), I,
               tolerance = 1e-13)
  expect_equal(sum(Gauss_Legen_weights(a = 0, b = pi, N = 2560) *
                     f(Gauss_Legen_nodes(a = 0, b = pi, N = 2560))), I,
               tolerance = 1e-13)
  expect_equal(sum(Gauss_Legen_weights(a = 0, b = pi, N = 5120) *
                     f(Gauss_Legen_nodes(a = 0, b = pi, N = 5120))), I,
               tolerance = 1e-14)

})

test_that("Gauss_Legen edge cases", {

  expect_error(Gauss_Legen_weights(a = 0, b = pi, N = 19))
  expect_error(Gauss_Legen_nodes(a = 0, b = pi, N = 19))

})