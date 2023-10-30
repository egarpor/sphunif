
n <- 10
set.seed(1234567)
X2 <- r_unif_sph(n = n, p = 2)
X3 <- r_unif_sph(n = n, p = 3)
X4 <- r_unif_sph(n = n, p = 4)
X5 <- r_unif_sph(n = n, p = 5)

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
