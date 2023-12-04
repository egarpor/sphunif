
set.seed(131121)
Theta_1 <- r_unif_cir(n = 33, M = 1)
Theta_2 <- r_unif_cir(n = 10, M = 2)
X_2 <- r_unif_sph(n = 33, p = 2, M = 1)
X_3 <- r_unif_sph(n = 10, p = 3, M = 2)
X_4 <- r_unif_sph(n = 10, p = 4, M = 2)
X_9 <- r_unif_sph(n = 10, p = 9, M = 2)
r_d_2 <- r_unif_sph(n = 5, p = 2, M = 1)[, , 1]
r_d_3 <- r_unif_sph(n = 5, p = 3, M = 1)[, , 1]
r_d_4 <- r_unif_sph(n = 5, p = 4, M = 1)[, , 1]
r_d_9 <- r_unif_sph(n = 5, p = 9, M = 1)[, , 1]
avail_cir_tests_s1 <- c("Gine_Fn", "Pycke", "Watson", "PCvM", "Cressie")
avail_cir_tests_s2 <- c("Gine_Gn", "Rothman", "Watson_1976", "PRt")
avail_sph_tests_s1 <- c("Gine_Fn", "Pycke", "PAD", "PCvM", "Cai")
avail_sph_tests_s2 <- c("Gine_Gn", "CCF09", "Ajne", "PRt")

test_that("Statistic arguments of unif_stat are the same of unif_stat_MC", {

  expect_true(all(
    names(formals(unif_stat))[-c(1:3)] ==
      names(formals(unif_stat_MC))[-c(1:12, length(formals(unif_stat_MC)))]
  ))

})

test_that("Statistic arguments of unif_stat are contained in those of
          unif_test", {

  expect_true(all(
    names(formals(unif_stat))[-c(1:3)] %in% names(formals(unif_test))[-c(1:8)]
  ))

})

test_that("unif_stat doing all or separate statistics", {

  expect_equal(as.numeric(unif_stat(data = Theta_1, type = "all",
                                    CCF09_dirs = r_d_2)),
               as.numeric(sapply(avail_cir_tests, function(test)
                 unif_stat(data = Theta_1, type = test,
                           CCF09_dirs = r_d_2))))
  expect_equal(as.matrix(unif_stat(data = Theta_2, type = "all",
                                   CCF09_dirs = r_d_2)),
               sapply(avail_cir_tests, function(test)
                 as.matrix(unif_stat(data = Theta_2, type = test,
                                     CCF09_dirs = r_d_2))))
  expect_equal(as.numeric(unif_stat(data = X_2, type = "all",
                                   CCF09_dirs = r_d_2)),
               as.numeric(sapply(avail_cir_tests, function(test)
                 unif_stat(data = X_2, type = test,
                           CCF09_dirs = r_d_2))))
  expect_equal(as.matrix(unif_stat(data = X_3, type = "all",
                                   CCF09_dirs = r_d_3)),
               sapply(avail_sph_tests, function(test)
                 as.matrix(unif_stat(data = X_3, type = test,
                                     CCF09_dirs = r_d_3))))
  suppressWarnings(expect_warning(
    expect_equal(as.matrix(unif_stat(data = X_4, type = "all",
                                     CCF09_dirs = r_d_4)),
                 sapply(avail_sph_tests, function(test)
                   as.matrix(unif_stat(data = X_4, type = test,
                                       CCF09_dirs = r_d_4))))
  ))
  suppressWarnings(expect_warning(
    expect_equal(as.matrix(unif_stat(data = X_9, type = "all",
                                     CCF09_dirs = r_d_9)),
                 sapply(avail_sph_tests, function(test)
                   as.matrix(unif_stat(data = X_9, type = test,
                                       CCF09_dirs = r_d_9))))
  ))

})

test_that("unif_stat with randomly-chosen and separate statistics", {

  expect_equal(as.numeric(unif_stat(data = Theta_1, type = avail_cir_tests_s1,
                                    CCF09_dirs = r_d_2)),
               as.numeric(sapply(avail_cir_tests_s1, function(test)
                 unif_stat(data = Theta_1, type = test,
                           CCF09_dirs = r_d_2))))
  expect_equal(as.numeric(unif_stat(data = Theta_1, type = avail_cir_tests_s2,
                                    CCF09_dirs = r_d_2)),
               as.numeric(sapply(avail_cir_tests_s2, function(test)
                 unif_stat(data = Theta_1, type = test,
                           CCF09_dirs = r_d_2))))
  expect_equal(as.numeric(unif_stat(data = X_2, type = avail_sph_tests_s2,
                                   CCF09_dirs = r_d_2)),
               as.numeric(sapply(avail_sph_tests_s2, function(test)
                 unif_stat(data = X_2, type = test,
                           CCF09_dirs = r_d_2))))
  expect_equal(as.matrix(unif_stat(data = X_3, type = avail_sph_tests_s2,
                                   CCF09_dirs = r_d_3)),
               sapply(avail_sph_tests_s2, function(test)
                 as.matrix(unif_stat(data = X_3, type = test,
                                     CCF09_dirs = r_d_3))))
  expect_equal(as.matrix(unif_stat(data = X_4, type = avail_sph_tests_s2,
                                   CCF09_dirs = r_d_4)),
               sapply(avail_sph_tests_s2, function(test)
                 as.matrix(unif_stat(data = X_4, type = test,
                                     CCF09_dirs = r_d_4))))

})

test_that("unif_stat with data_sorted = TRUE and CCF09 = NULL", {

  expect_equal(as.numeric(unif_stat(data = sort_each_col(Theta_1),
                                    type = avail_cir_tests_s1,
                                    CCF09_dirs = r_d_2,
                                    data_sorted = TRUE)),
               as.numeric(sapply(avail_cir_tests_s1, function(test)
                 unif_stat(data = sort_each_col(Theta_1), type = test,
                           CCF09_dirs = r_d_2,
                           data_sorted = TRUE))))
  expect_equal({
    set.seed(131312)
    rd <-  CCF09_dirs <- r_unif_sph(n = 50, p = 2, M = 1)[, , 1]
    as.matrix(unif_stat(data = Theta_2, type = c("Range", "Rao",
                                                 "CCF09"),
                         CCF09_dirs = rd))
    }, {
    set.seed(131312)
    as.matrix(unif_stat(data = Theta_2, type = c("Range", "Rao",
                                                 "CCF09")))
    })
  expect_equal({
    set.seed(131312)
    r_d <-  CCF09_dirs <- r_unif_sph(n = 50, p = 3, M = 1)[, , 1]
    as.matrix(unif_stat(data = X_3, type = c("PCvM", "Bakshaev",
                                             "CCF09"),
                        CCF09_dirs = r_d))
    }, {
    set.seed(131312)
    as.matrix(unif_stat(data = X_3, type = c("PCvM", "Bakshaev",
                                             "CCF09")))
    })

})

test_that("unif_stat equality for data = Theta and data = Theta_to_X(Theta)", {

  expect_equal(unif_stat(data = Theta_1, type = avail_cir_tests,
                         CCF09_dirs = r_d_2),
               unif_stat(data = Theta_to_X(Theta_1), type = avail_cir_tests,
                         CCF09_dirs = r_d_2))
  expect_equal(unif_stat(data = X_2, type = avail_sph_tests,
                         CCF09_dirs = r_d_2),
               unif_stat(data = X_to_Theta(X_2), type = avail_sph_tests,
                         CCF09_dirs = r_d_2))

})

test_that("KS, CvM, and AD", {

  expect_equal(as.numeric(unif_stat(data = Theta_1, type = "KS")),
               as.numeric(cir_stat_Kuiper(Theta = Theta_1, KS = TRUE)))
  expect_equal(as.numeric(unif_stat(data = Theta_1, type = "CvM")),
               as.numeric(cir_stat_Watson(Theta = Theta_1, CvM = TRUE)))
  expect_equal(as.numeric(unif_stat(data = Theta_1, type = "AD")),
               as.numeric(goftest::ad.test(Theta_1, null = "punif",
                                           min = 0, max = 2 * pi)$statistic))
  expect_equal(unif_stat(data = Theta_2, type = "KS")$KS,
               cir_stat_Kuiper(Theta = Theta_2, KS = TRUE))
  expect_equal(unif_stat(data = Theta_2, type = "CvM")$CvM,
               cir_stat_Watson(Theta = Theta_2, CvM = TRUE))
  expect_equal(as.numeric(unif_stat(data = Theta_2, type = "AD")$AD),
               apply(Theta_2, 2, function(x)
                 as.numeric(goftest::ad.test(x, null = "punif", min = 0,
                                             max = 2 * pi)$statistic)))

})

test_that("Riesz vs. Rayleigh", {

  expect_equal(
    as.numeric(unif_stat(data = Theta_2, type = c("Riesz", "Rayleigh"),
                         Riesz_s = 2)$Riesz),
    as.numeric(unif_stat(data = Theta_2, type = "Riesz", Riesz_s = 2)$Riesz)
  )
  expect_equal(
    as.numeric(unif_stat(data = Theta_2, type = c("Riesz", avail_sph_tests_s2),
                         Riesz_s = 2)$Riesz),
    as.numeric(unif_stat(data = Theta_2, type = "Riesz", Riesz_s = 2)$Riesz)
  )
  expect_equal(
    as.numeric(unif_stat(data = X_2, type = c("Riesz", "Rayleigh"),
                         Riesz_s = 2)$Riesz),
    as.numeric(unif_stat(data = X_2, type = "Riesz", Riesz_s = 2)$Riesz)
  )
  expect_equal(
    as.numeric(unif_stat(data = X_2, type = c("Riesz", avail_sph_tests_s2),
                         Riesz_s = 2)$Riesz),
    as.numeric(unif_stat(data = X_2, type = "Riesz", Riesz_s = 2)$Riesz)
  )
  expect_equal(
    as.numeric(unif_stat(data = X_3, type = c("Riesz", "Rayleigh"),
                         Riesz_s = 2)$Riesz),
    as.numeric(unif_stat(data = X_3, type = "Riesz", Riesz_s = 2)$Riesz)
  )
  expect_equal(
    as.numeric(unif_stat(data = X_3, type = c("Riesz", avail_sph_tests_s2),
                         Riesz_s = 2)$Riesz),
    as.numeric(unif_stat(data = X_3, type = "Riesz", Riesz_s = 2)$Riesz)
  )
  expect_equal(
    as.numeric(unif_stat(data = X_4, type = c("Riesz", "Rayleigh"),
                         Riesz_s = 2)$Riesz),
    as.numeric(unif_stat(data = X_4, type = "Riesz", Riesz_s = 2)$Riesz)
  )
  expect_equal(
    as.numeric(unif_stat(data = X_4, type = c("Riesz", avail_sph_tests_s2),
                         Riesz_s = 2)$Riesz),
    as.numeric(unif_stat(data = X_4, type = "Riesz", Riesz_s = 2)$Riesz)
  )
  expect_equal(
    as.numeric(unif_stat(data = X_9, type = c("Riesz", "Rayleigh"),
                         Riesz_s = 2)$Riesz),
    as.numeric(unif_stat(data = X_9, type = "Riesz", Riesz_s = 2)$Riesz)
  )
  expect_equal(
    as.numeric(unif_stat(data = X_9, type = c("Riesz", avail_sph_tests_s2),
                         Riesz_s = 2)$Riesz),
    as.numeric(unif_stat(data = X_9, type = "Riesz", Riesz_s = 2)$Riesz)
  )

})

test_that("Riesz vs. Pycke", {

  expect_equal(
    as.numeric(unif_stat(data = Theta_2, type = c("Riesz", "Pycke"),
                         Riesz_s = 0)$Pycke),
    as.numeric(unif_stat(data = Theta_2, type = "Pycke")$Pycke)
  )
  expect_equal(
    as.numeric(unif_stat(data = X_2, type = c("Riesz", "Pycke"),
                         Riesz_s = 0)$Pycke),
    as.numeric(unif_stat(data = X_2, type = "Pycke")$Pycke)
  )
  expect_equal(
    as.numeric(unif_stat(data = X_3, type = c("Riesz", "Pycke"),
                         Riesz_s = 0)$Pycke),
    as.numeric(unif_stat(data = X_3, type = "Pycke")$Pycke)
  )
  suppressWarnings(expect_warning(
    expect_equal(
      as.numeric(unif_stat(data = X_4, type = c("Riesz", "Pycke"),
                           Riesz_s = 0)$Riesz),
      as.numeric(unif_stat(data = X_4, type = "Pycke")$Pycke)
    )))
  suppressWarnings(expect_warning(
    expect_equal(
      as.numeric(unif_stat(data = X_9, type = c("Riesz", "Pycke"),
                           Riesz_s = 0)$Riesz),
      as.numeric(unif_stat(data = X_9, type = "Pycke")$Pycke)
    )))

})

test_that("Rayleigh vs. Softmax and Poisson", {
  for (m in 1:2) {
    expect_equal(
      unif_stat(data = Theta_2, type = c("Softmax", "Poisson"),
                Softmax_kappa = 0, Poisson_rho = 0, Rayleigh_m = m)[, 1:2],
      unif_stat(data = Theta_2, type = c("Softmax", "Poisson", "Rayleigh"),
                Softmax_kappa = 0, Poisson_rho = 0, Rayleigh_m = m)[, 1:2]
    )
  }
  expect_equal(
    unif_stat(data = X_2, type = c("Softmax", "Poisson"),
              Softmax_kappa = 0, Poisson_rho = 0)[, 1:2],
    unif_stat(data = X_2, type = c("Softmax", "Poisson", "Rayleigh"),
              Softmax_kappa = 0, Poisson_rho = 0)[, 1:2]
  )
  expect_equal(
    unif_stat(data = X_3, type = c("Softmax", "Poisson"),
              Softmax_kappa = 0, Poisson_rho = 0)[, 1:2],
    unif_stat(data = X_3, type = c("Softmax", "Poisson", "Rayleigh"),
              Softmax_kappa = 0, Poisson_rho = 0)[, 1:2]
  )
  expect_equal(
    unif_stat(data = X_4, type = c("Softmax", "Poisson"),
              Softmax_kappa = 0, Poisson_rho = 0)[, 1:2],
    unif_stat(data = X_4, type = c("Softmax", "Poisson", "Rayleigh"),
              Softmax_kappa = 0, Poisson_rho = 0)[, 1:2]
  )
})

test_that("Vectorization works in Stereo test", {
  expect_equal(unname(as.matrix(unif_stat(data = X_3, type = "Stereo",
                                Stereo_a = c(-0.5, 0, 0.5)))),
               unname(cbind(
                 unif_stat(data = X_3, type = "Stereo", Stereo_a = -0.5)$Stereo,
                 unif_stat(data = X_3, type = "Stereo", Stereo_a = 0)$Stereo,
                 unif_stat(data = X_3, type = "Stereo", Stereo_a = 0.5)$Stereo)
                 ))
})

# Statistics with vectorised parameters
cir_stats_vectorised <- c("Cressie", "Max_uncover", "Num_uncover", "Vacancy",
                          "Rayleigh", "Riesz", "Rothman", "PRt", "Poisson",
                          "Pycke_q", "Softmax", "Sobolev")
sph_stats_vectorised <- c("Riesz", "PRt", "Poisson", "Softmax", "Stereo",
                          "Sobolev")
t <- c(0.2, 0.3, 0.8)
m <- 1:3
s <- c(0, 1, 2)
kappa <- 1:3
rho <- seq(0.1, 0.9, l = 3)
vk2 <- rbind(1:3, 3:1, 3:5)

test_that("Parameter-vectorized statistics work for p = 2", {
  stats_1 <- as.matrix(
    unif_stat(data = X_2, type = cir_stats_vectorised,
              Cressie_t = t, cov_a = t, Rayleigh_m = m, Riesz_s = s,
              Rothman_t = t, Softmax_kappa = kappa, Poisson_rho = rho,
              Pycke_q = t, Sobolev_vk2 = vk2))
  stats_2 <- lapply(1:3, function(i) as.matrix(
    unif_stat(data = X_2, type = cir_stats_vectorised,
              Cressie_t = t[i], cov_a = t[i], Rayleigh_m = m[i], Riesz_s = s[i],
              Rothman_t = t[i], Softmax_kappa = kappa[i], Poisson_rho = rho[i],
              Pycke_q = t[i], Sobolev_vk2 = vk2[i, ])))
  stats_2 <- unname(do.call(cbind, stats_2))
  stats_1 <- unname(stats_1)
  stats_1 <- stats_1[, order(stats_1[1, ])]
  stats_2 <- stats_2[, order(stats_2[1, ])]
  expect_equal(stats_1, stats_2)
})

test_that("Parameter-vectorized statistics work for p = 4", {
  stats_1 <- as.matrix(
    unif_stat(data = X_4, type = sph_stats_vectorised, Riesz_s = s,
              Rothman_t = t, Softmax_kappa = kappa, Poisson_rho = rho,
              Stereo_a = t, Sobolev_vk2 = vk2))
  stats_2 <- lapply(1:3, function(i) as.matrix(
    unif_stat(data = X_4, type = sph_stats_vectorised, Riesz_s = s[i],
              Rothman_t = t[i], Softmax_kappa = kappa[i],
              Poisson_rho = rho[i], Stereo_a = t[i],
              Sobolev_vk2 = vk2[i, ])))
  stats_2 <- unname(do.call(cbind, stats_2))
  stats_1 <- unname(stats_1)
  stats_1 <- stats_1[, order(stats_1[1, ])]
  stats_2 <- stats_2[, order(stats_2[1, ])]
  expect_equal(stats_1, stats_2)
})

test_that("Errors in edge cases", {

  expect_error(unif_stat(data = rbind(Theta_2, c(0, NA))))
  expect_error(unif_stat(data = rbind(X_3, c(0, NA, 1))))
  expect_error(unif_stat(data = rbind(c(1, 0))))
  expect_error(unif_stat(data = 1))
  expect_error(unif_stat(data = X_3, type = "Invent"))
  expect_error(unif_test(data = X_3, type = 1e3))
  expect_error(unif_test(data = X_3, type = function(x) x))

})

test_that("Passing edge cases", {

  expect_equal(unif_stat(data = 1:6, type = c("Ajne", "Cressie")),
               unif_stat(data = cbind(1:6), type = c("Ajne", "Cressie")))
  expect_equal(unif_stat(data = cbind(1:5, 2:6), type = c("Ajne", "Cressie")),
               {A <- array(dim = c(5, 1, 2)); A[, 1, ] <- cbind(1:5, 2:6);
               unif_stat(data = A, type = c("Ajne", "Cressie"))})
  expect_equal(c(as.matrix(unif_stat(data = Theta_2,
                                     type = avail_cir_tests[2:5],
                                     CCF09_dirs = r_d_2))),
               c(sapply(2:5, function(test)
                 as.matrix(unif_stat(data = Theta_2, type = test,
                                     CCF09_dirs = r_d_2)))))
  expect_equal(c(as.matrix(unif_stat(data = X_3, type = avail_sph_tests[2:5],
                                     CCF09_dirs = r_d_3))),
               c(sapply(2:5, function(test)
                 as.matrix(unif_stat(data = X_3, type = test,
                                     CCF09_dirs = r_d_3)))))
  expect_equal(unif_stat(data = 1:6, type = c("Ajne", "Ajne", "Watson")),
               unif_stat(data = 1:6, type = c("Ajne", "Watson")))

})
