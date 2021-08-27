
# To prevent hanging with parallel computations
Sys.unsetenv("R_TESTS")

one <- function(x) rep(1, nrow(x))
id <- function(x) x[, 1]
quad <- function(x) rowSums(x^2)
cubic <- function(x, a = 0) a + x[, 1] * x[, 2] * x[, 3]
time_consuming <- function(x) {Sys.sleep(0.5); 1}
seeds <- 1:5

test_that("Integration of one", {

  skip_on_cran()
  for (p in c(2, 3, 4, 11)) {
    expect_equal(int_sph_MC(f = one, p = p, M = 10, verbose = FALSE),
                 rotasym::w_p(p = p))
  }

})

test_that("Integration of id", {

  skip_on_cran()
  for (p in c(2, 3, 4, 11)) {
    expect_equal(int_sph_MC(f = id, p = p, M = 1e3, verbose = FALSE,
                            seeds = 1:5, chunks = 5), 0,
                 tolerance = 3e-1)
  }

})

test_that("Integration of quad", {

  skip_on_cran()
  for (p in c(2, 3, 4, 11)) {
    expect_equal(int_sph_MC(f = quad, p = p, M = 10, verbose = FALSE,
                            seeds = 1:5, chunks = 5),
                 rotasym::w_p(p = p))
  }

})

test_that("Integration of cubic", {

  skip_on_cran()
  for (p in c(3, 4, 11)) {
    expect_equal(int_sph_MC(f = cubic, p = p, M = 1e3, verbose = FALSE,
                            seeds = 1:5, chunks = 5), 0,
                 tolerance = 1e-1)
  }

})

test_that("Optional argument to f", {

  skip_on_cran()
  expect_equal(int_sph_MC(f = cubic, p = 5, M = 10, seeds = seeds,
                          chunks = 5, verbose = FALSE),
               int_sph_MC(f = cubic, p = 5, M = 10, a = 0, seeds = seeds,
                          chunks = 5, verbose = FALSE))
  expect_equal(int_sph_MC(f = cubic, p = 4, M = 10, a = -5, seeds = seeds,
                          chunks = 5, verbose = FALSE),
               -5 * rotasym::w_p(p = 4) +
                 int_sph_MC(f = cubic, p = 4, M = 10, a = 0, seeds = seeds,
                            chunks = 5, verbose = FALSE))

})

test_that("Edge cases", {

  skip_on_cran()
  expect_warning(int_sph_MC(f = cubic, p = 3, M = 10, seeds = seeds[1:2],
                            chunks = 5))

})

test_that("Progress bars", {

  skip_on_cran()
  o1_verbose <- capture.output(i <- int_sph_MC(f = cubic, p = 3, M = 10,
                                               cores = 1, verbose = TRUE))
  o1_silent <- capture.output(i <- int_sph_MC(f = cubic, p = 3, M = 10,
                                              cores = 1, verbose = FALSE))
  o2_verbose <- capture.output(i <- int_sph_MC(f = cubic, p = 3, M = 10,
                                               cores = 2, verbose = TRUE))
  o2_silent <- capture.output(i <- int_sph_MC(f = cubic, p = 3, M = 10,
                                              cores = 2, verbose = FALSE))
  expect_gt(length(o1_verbose), 0)
  expect_gt(length(o2_verbose), 0)
  expect_equal(length(o1_silent), 0)
  expect_equal(length(o2_silent), 0)

})

test_that("Same results with cores = 1 and cores = 2", {

  skip_on_cran()
  expect_equal(int_sph_MC(f = cubic, p = 3, M = 10, cores = 1, seeds = seeds,
                          chunks = 5, verbose = FALSE),
               int_sph_MC(f = cubic, p = 3, M = 10, cores = 2, seeds = seeds,
                          chunks = 5, verbose = FALSE))

})

test_that("Parallelization is faster", {

  skip_on_cran()
  t1 <- system.time(int_sph_MC(f = time_consuming, p = 3, M = 1e5, cores = 1,
                               chunks = 10, verbose = TRUE))[3]
  t2 <- system.time(int_sph_MC(f = time_consuming, p = 3, M = 1e5, cores = 2,
                               chunks = 10, verbose = TRUE))[3]
  expect_gt(t1, t2)

})
