
Sys.unsetenv("R_TESTS")

one <- function(x) rep(1, nrow(x))
id <- function(x) x[, 1]
quad <- function(x) rowSums(x^2)
cubic <- function(x, a = 0) a + x[, 1] * x[, 2] * x[, 3]
seeds <- 1:5

test_that("Integration of one", {

  expect_equal(int_sph_MC(f = one, p = 2, M = 10, verbose = FALSE),
               rotasym::w_p(p = 2))
  expect_equal(int_sph_MC(f = one, p = 3, M = 10, verbose = FALSE),
               rotasym::w_p(p = 3))
  expect_equal(int_sph_MC(f = one, p = 4, M = 10, verbose = FALSE),
               rotasym::w_p(p = 4))
  expect_equal(int_sph_MC(f = one, p = 11, M = 10, verbose = FALSE),
               rotasym::w_p(p = 11))

})

test_that("Integration of id", {

  expect_equal(int_sph_MC(f = id, p = 2, M = 1e3, verbose = FALSE,
                          seeds = 1:5, chunks = 5), 0,
               tolerance = 3e-1)
  expect_equal(int_sph_MC(f = id, p = 3, M = 1e3, verbose = FALSE,
                          seeds = 1:5, chunks = 5), 0,
               tolerance = 3e-1)
  expect_equal(int_sph_MC(f = id, p = 4, M = 1e3, verbose = FALSE,
                          seeds = 1:5, chunks = 5), 0,
               tolerance = 3e-1)
  expect_equal(int_sph_MC(f = id, p = 11, M = 1e3, verbose = FALSE,
                          seeds = 1:5, chunks = 5), 0,
               tolerance = 3e-1)

})

test_that("Integration of quad", {

  expect_equal(int_sph_MC(f = quad, p = 2, M = 10, verbose = FALSE,
                          seeds = 1:5, chunks = 5),
               rotasym::w_p(p = 2))
  expect_equal(int_sph_MC(f = quad, p = 3, M = 10, verbose = FALSE,
                          seeds = 1:5, chunks = 5),
               rotasym::w_p(p = 3))
  expect_equal(int_sph_MC(f = quad, p = 4, M = 10, verbose = FALSE,
                          seeds = 1:5, chunks = 5),
               rotasym::w_p(p = 4))
  expect_equal(int_sph_MC(f = quad, p = 11, M = 10, verbose = FALSE,
                          seeds = 1:5, chunks = 5),
               rotasym::w_p(p = 11))

})

test_that("Integration of cubic", {

  expect_equal(int_sph_MC(f = cubic, p = 3, M = 1e3, verbose = FALSE,
                          seeds = 1:5, chunks = 5), 0,
               tolerance = 1e-1)
  expect_equal(int_sph_MC(f = cubic, p = 4, M = 1e3, verbose = FALSE,
                          seeds = 1:5, chunks = 5), 0,
               tolerance = 1e-1)
  expect_equal(int_sph_MC(f = cubic, p = 11, M = 1e3, verbose = FALSE,
                          seeds = 1:5, chunks = 5), 0,
               tolerance = 1e-1)

})

test_that("cores = 1 vs. cores = 2", {

  expect_equal(int_sph_MC(f = cubic, p = 3, M = 10, cores = 1, seeds = seeds,
                          chunks = 5, verbose = FALSE),
               int_sph_MC(f = cubic, p = 3, M = 10, cores = 2, seeds = seeds,
                          chunks = 5, verbose = FALSE))
  expect_equal(int_sph_MC(f = id, p = 3, M = 10, cores = 1, seeds = seeds,
                          chunks = 5, verbose = FALSE),
               int_sph_MC(f = id, p = 3, M = 10, cores = 2, seeds = seeds,
                          chunks = 5, verbose = FALSE))

})

test_that("Optional argument to f", {

  expect_equal(int_sph_MC(f = cubic, p = 5, M = 10, cores = 1,
                          seeds = seeds, chunks = 5, verbose = FALSE),
               int_sph_MC(f = cubic, p = 5, M = 10, cores = 1, a = 0,
                          seeds = seeds, chunks = 5, verbose = FALSE))
  expect_equal(int_sph_MC(f = cubic, p = 4, M = 10, cores = 1, a = -5,
                          seeds = seeds, chunks = 5, verbose = FALSE),
               -5 * rotasym::w_p(p = 4) +
                 int_sph_MC(f = cubic, p = 4, M = 10, cores = 1, a = 0,
                            seeds = seeds, chunks = 5, verbose = FALSE))

})

test_that("int_sph_MC edge cases", {

  expect_warning(capture.output(int_sph_MC(f = cubic, p = 3, M = 10, cores = 1,
                                           seeds = seeds[1:2], chunks = 5,
                                           verbose = TRUE)))

})

