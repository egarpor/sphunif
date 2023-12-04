
n <- 10
set.seed(1234567)
X2 <- r_unif_sph(n = n, p = 2)
X3 <- r_unif_sph(n = n, p = 3)
X4 <- r_unif_sph(n = n, p = 4)
X5 <- r_unif_sph(n = n, p = 5)
Sobolev_vk2 <- diag(1, nrow = 2)

test_that("All test work for p = 2 and return valid class information", {

  for (type in avail_cir_tests) {
    expect_no_error(test <- unif_test(X2, type = type, p_value = "MC", M = 2))
    expect_true(is.numeric(test$statistic))
    expect_true(is.numeric(test$p.value))
    expect_false(is.null(test$method))
    expect_false(is.null(test$alternative))
  }

})

test_that("All test work for p = 3 and return valid class information", {

  for (type in avail_sph_tests) {
    expect_no_error(test <- unif_test(X3, type = type, p_value = "MC", M = 2))
    expect_true(is.numeric(test$statistic))
    expect_true(is.numeric(test$p.value))
    expect_false(is.null(test$method))
    expect_false(is.null(test$alternative))
  }

})

test_that("All test work for p = 4 and return valid class information", {

  for (type in avail_sph_tests) {
    if (type == "Pycke") {
      suppressWarnings(expect_no_error(expect_warning(
        test <- unif_test(X4, type = type, p_value = "MC", M = 2))))
    } else {
      expect_no_error(test <- unif_test(X4, type = type, p_value = "MC", M = 2))
    }
    expect_true(is.numeric(test$statistic))
    expect_true(is.numeric(test$p.value))
    expect_false(is.null(test$method))
    expect_false(is.null(test$alternative))
  }

})

test_that("All test work for p = 5 and return valid class information", {

  for (type in avail_sph_tests) {
    if (type == "Pycke") {
      suppressWarnings(expect_no_error(expect_warning(
        test <- unif_test(X5, type = type, p_value = "MC", M = 2))))
    } else {
      expect_no_error(test <- unif_test(X5, type = type, p_value = "MC", M = 2))
    }
    expect_true(is.numeric(test$statistic))
    expect_true(is.numeric(test$p.value))
    expect_false(is.null(test$method))
    expect_false(is.null(test$alternative))
  }

})

test_that("Vectorization works for Sobolev test", {

    expect_no_error(test_mc <- unif_test(X2, type = "Sobolev",
                                         p_value = "MC", M = 100,
                                         Sobolev_vk2 = Sobolev_vk2))
    expect_no_error(test_as <- unif_test(X2, type = "Sobolev",
                                         p_value = "asymp",
                                         Sobolev_vk2 = Sobolev_vk2))
    expect_equal(test_mc$p.value, test_as$p.value, tolerance = 0.1)
    expect_no_error(test_mc <- unif_test(X3, type = "Sobolev",
                                         p_value = "MC", M = 100,
                                         Sobolev_vk2 = Sobolev_vk2))
    expect_no_error(test_as <- unif_test(X3, type = "Sobolev",
                                         p_value = "asymp",
                                         Sobolev_vk2 = Sobolev_vk2))
    expect_equal(test_mc$p.value, test_as$p.value, tolerance = 0.1)

})

test_that("Vectorization works for Sobolev test with other tests", {

  expect_no_error(test_mc <- unif_test(X2, type = c("Sobolev", "Ajne"),
                                       p_value = "MC", M = 100,
                                       Sobolev_vk2 = Sobolev_vk2))
  expect_no_error(test_as <- unif_test(X2, type = c("Sobolev", "Ajne"),
                                       p_value = "asymp",
                                       Sobolev_vk2 = Sobolev_vk2))
  expect_equal(test_mc$p.value, test_as$p.value, tolerance = 0.1)
  expect_no_error(test_mc <- unif_test(X3, type = c("Sobolev", "Ajne"),
                                       p_value = "MC", M = 100,
                                       Sobolev_vk2 = Sobolev_vk2))
  expect_no_error(test_as <- unif_test(X3, type = c("Sobolev", "Ajne"),
                                       p_value = "asymp",
                                       Sobolev_vk2 = Sobolev_vk2))
  expect_equal(test_mc$p.value, test_as$p.value, tolerance = 0.1)

})

test_that("Errors in edge cases", {

  expect_error(unif_test(data = rbind(X2, c(NA, NA))))
  expect_error(unif_test(data = rbind(X3, c(0, NA, 1))))
  expect_error(unif_test(data = X3, type = "Invent"))
  expect_error(unif_test(data = X3, type = 1e3))
  expect_error(unif_test(data = X3, type = function(x) x))

})
