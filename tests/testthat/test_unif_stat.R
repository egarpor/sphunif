
set.seed(131121)
Theta_1 <- r_unif_cir(n = 33, M = 1)
Theta_2 <- r_unif_cir(n = 10, M = 2)
X_2 <- r_unif_sph(n = 33, p = 2, M = 1)
X_3 <- r_unif_sph(n = 10, p = 3, M = 2)
X_4 <- r_unif_sph(n = 10, p = 4, M = 2)
r_d_2 <- r_unif_sph(n = 5, p = 2, M = 1)[, , 1]
r_d_3 <- r_unif_sph(n = 5, p = 3, M = 1)[, , 1]
rand_dirs_4 <- r_unif_sph(n = 5, p = 4, M = 1)[, , 1]
avail_cir_tests_s1 <- c("Gine_Fn", "Pycke", "Watson", "PCvM", "Cressie")
avail_cir_tests_s2 <- c("Gine_Gn", "Rothman", "Watson_1976", "PRt")
avail_sph_tests_s1 <- c("Gine_Fn", "Pycke", "PAD", "PCvM", "Cai")
avail_sph_tests_s2 <- c("Gine_Gn", "Cuesta_Albertos", "Ajne", "PRt")

test_that("unif_stat doing all or separate statistics", {

  expect_equal(as.numeric(unif_stat(data = Theta_1, type = "all",
                                    Cuesta_Albertos_rand_dirs = r_d_2)),
               as.numeric(sapply(avail_cir_tests, function(test)
                 unif_stat(data = Theta_1, type = test,
                           Cuesta_Albertos_rand_dirs = r_d_2))))
  expect_equal(as.matrix(unif_stat(data = Theta_2, type = "all",
                                   Cuesta_Albertos_rand_dirs = r_d_2)),
               sapply(avail_cir_tests, function(test)
                 as.matrix(unif_stat(data = Theta_2, type = test,
                                     Cuesta_Albertos_rand_dirs = r_d_2))))
  expect_equal(as.numeric(unif_stat(data = X_2, type = "all",
                                   Cuesta_Albertos_rand_dirs = r_d_2)),
               as.numeric(sapply(avail_cir_tests, function(test)
                 unif_stat(data = X_2, type = test,
                           Cuesta_Albertos_rand_dirs = r_d_2))))
  expect_equal(as.matrix(unif_stat(data = X_3, type = "all",
                                   Cuesta_Albertos_rand_dirs = r_d_3)),
               sapply(avail_sph_tests, function(test)
                 as.matrix(unif_stat(data = X_3, type = test,
                                     Cuesta_Albertos_rand_dirs = r_d_3))))
  expect_equal(as.matrix(unif_stat(data = X_4, type = "all",
                                   Cuesta_Albertos_rand_dirs = rand_dirs_4)),
               sapply(avail_sph_tests, function(test)
                 as.matrix(unif_stat(data = X_4, type = test,
                                     Cuesta_Albertos_rand_dirs = rand_dirs_4))))

})

test_that("unif_stat with randomly-chosen and separate statistics", {

  expect_equal(as.numeric(unif_stat(data = Theta_1, type = avail_cir_tests_s1,
                                    Cuesta_Albertos_rand_dirs = r_d_2)),
               as.numeric(sapply(avail_cir_tests_s1, function(test)
                 unif_stat(data = Theta_1, type = test,
                           Cuesta_Albertos_rand_dirs = r_d_2))))
  expect_equal(as.numeric(unif_stat(data = Theta_1, type = avail_cir_tests_s2,
                                    Cuesta_Albertos_rand_dirs = r_d_2)),
               as.numeric(sapply(avail_cir_tests_s2, function(test)
                 unif_stat(data = Theta_1, type = test,
                           Cuesta_Albertos_rand_dirs = r_d_2))))
  expect_equal(as.numeric(unif_stat(data = X_2, type = avail_sph_tests_s2,
                                   Cuesta_Albertos_rand_dirs = r_d_2)),
               as.numeric(sapply(avail_sph_tests_s2, function(test)
                 unif_stat(data = X_2, type = test,
                           Cuesta_Albertos_rand_dirs = r_d_2))))
  expect_equal(as.matrix(unif_stat(data = X_3, type = avail_sph_tests_s2,
                                   Cuesta_Albertos_rand_dirs = r_d_3)),
               sapply(avail_sph_tests_s2, function(test)
                 as.matrix(unif_stat(data = X_3, type = test,
                                     Cuesta_Albertos_rand_dirs = r_d_3))))
  expect_equal(as.matrix(unif_stat(data = X_4, type = avail_sph_tests_s2,
                                   Cuesta_Albertos_rand_dirs = rand_dirs_4)),
               sapply(avail_sph_tests_s2, function(test)
                 as.matrix(unif_stat(data = X_4, type = test,
                                     Cuesta_Albertos_rand_dirs = rand_dirs_4))))

})

test_that("unif_stat with data_sorted = TRUE and Cuesta_Albertos = NULL", {

  expect_equal(as.numeric(unif_stat(data = sort_each_col(Theta_1),
                                    type = avail_cir_tests_s1,
                                    Cuesta_Albertos_rand_dirs = r_d_2,
                                    data_sorted = TRUE)),
               as.numeric(sapply(avail_cir_tests_s1, function(test)
                 unif_stat(data = sort_each_col(Theta_1), type = test,
                           Cuesta_Albertos_rand_dirs = r_d_2,
                           data_sorted = TRUE))))
  expect_equal({
    set.seed(131312)
    rd <-  Cuesta_Albertos_rand_dirs <- r_unif_sph(n = 50, p = 2, M = 1)[, , 1]
    as.matrix(unif_stat(data = Theta_2, type = c("Range", "Rao",
                                                 "Cuesta_Albertos"),
                         Cuesta_Albertos_rand_dirs = rd))
    }, {
    set.seed(131312)
    as.matrix(unif_stat(data = Theta_2, type = c("Range", "Rao",
                                                 "Cuesta_Albertos")))
    })
  expect_equal({
    set.seed(131312)
    r_d <-  Cuesta_Albertos_rand_dirs <- r_unif_sph(n = 50, p = 3, M = 1)[, , 1]
    as.matrix(unif_stat(data = X_3, type = c("PCvM", "Bakshaev",
                                             "Cuesta_Albertos"),
                        Cuesta_Albertos_rand_dirs = r_d))
    }, {
    set.seed(131312)
    as.matrix(unif_stat(data = X_3, type = c("PCvM", "Bakshaev",
                                             "Cuesta_Albertos")))
    })

})

test_that("unif_stat equality for data = Theta and data = Theta_to_X(Theta)", {

  expect_equal(unif_stat(data = Theta_1, type = avail_cir_tests,
                         Cuesta_Albertos_rand_dirs = r_d_2),
               unif_stat(data = Theta_to_X(Theta_1), type = avail_cir_tests,
                         Cuesta_Albertos_rand_dirs = r_d_2))
  expect_equal(unif_stat(data = X_2, type = avail_sph_tests,
                         Cuesta_Albertos_rand_dirs = r_d_2),
               unif_stat(data = X_to_Theta(X_2), type = avail_sph_tests,
                         Cuesta_Albertos_rand_dirs = r_d_2))

})

test_that("Errors in edge cases", {

  expect_error(unif_stat(data = rbind(Theta_2, c(0, NA))))
  expect_error(unif_stat(data = rbind(X_3, c(0, NA, 1))))
  expect_error(unif_stat(data = rbind(c(1, 0))))
  expect_error(unif_stat(data = 1))
  expect_error(unif_stat(data = 1:3, type = function(x) x^2))

})

test_that("Passing edge cases", {

  expect_equal(unif_stat(data = 1:6, type = c("Ajne", "Cressie")),
               unif_stat(data = cbind(1:6), type = c("Ajne", "Cressie")))
  expect_equal(unif_stat(data = cbind(1:5, 2:6), type = c("Ajne", "Cressie")),
               {A <- array(dim = c(5, 1, 2)); A[, 1, ] <- cbind(1:5, 2:6);
               unif_stat(data = A, type = c("Ajne", "Cressie"))})
  expect_equal(c(as.matrix(unif_stat(data = Theta_2,
                                     type = avail_cir_tests[2:5],
                                     Cuesta_Albertos_rand_dirs = r_d_2))),
               c(sapply(2:5, function(test)
                 as.matrix(unif_stat(data = Theta_2, type = test,
                                     Cuesta_Albertos_rand_dirs = r_d_2)))))
  expect_equal(c(as.matrix(unif_stat(data = X_3, type = avail_sph_tests[2:5],
                                     Cuesta_Albertos_rand_dirs = r_d_3))),
               c(sapply(2:5, function(test)
                 as.matrix(unif_stat(data = X_3, type = test,
                                     Cuesta_Albertos_rand_dirs = r_d_3)))))
  expect_equal(unif_stat(data = 1:6, type = c("Ajne", "Ajne", "Watson")),
               unif_stat(data = 1:6, type = c("Ajne", "Watson")))

})

test_that("KS and CvM", {

  expect_equal(as.numeric(unif_stat(data = Theta_1, type = "KS")),
               as.numeric(cir_stat_Kuiper(Theta = Theta_1, KS = TRUE)))
  expect_equal(as.numeric(unif_stat(data = Theta_1, type = "CvM")),
               as.numeric(cir_stat_Watson(Theta = Theta_1, CvM = TRUE)))
  expect_equal(unif_stat(data = Theta_2, type = "KS")$KS,
               cir_stat_Kuiper(Theta = Theta_2, KS = TRUE))
  expect_equal(unif_stat(data = Theta_2, type = "CvM")$CvM,
               cir_stat_Watson(Theta = Theta_2, CvM = TRUE))

})
