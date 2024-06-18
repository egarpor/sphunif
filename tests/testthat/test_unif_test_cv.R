n <- 100
K_list <- c(1, 2, 5, 10, 23, 50)
rho_list <- c(0.1, 0.25, 0.75)
kappa_list <- c(0.1, 1, 10)
a_list <- c(-1, 0, 1)

test_that("K-fold partition selects every and just once each observation", {

  for (K in K_list){

    folds <- k_fold_split(n, K)
    expect_equal(n, length(unique(unlist(folds))))

  }

})

test_that("K-fold partitions (min, max, and mean subsizes) are well balanced", {

  for (K in K_list){

    folds <- k_fold_split(n, K)
    theo_length <- n / K
    fold_lengths <- unlist(lapply(folds, length))

    expect_equal(mean(fold_lengths), theo_length, tolerance = 0)
    expect_lte(max(fold_lengths), as.integer(theo_length) + 1)
    expect_equal(min(fold_lengths), as.integer(theo_length), tolerance = 0)

  }

})

test_that("K-fold partition when K > n/2 must throw an error", {

  expect_error(k_fold_split(n = n, K = n + 2))
  expect_error(k_fold_split(n = n, K = n / 2 + 1))

})

test_that("Prepare test data returns correct n, p, and correct sizes", {

  n <- 10
  # Spherical structure
  for (p in c(2:4, 11)) {

    data <- r_unif_sph(n = n, p = p, M = 1)
    d <- prepare_test_data(data)

    expect_equal(d[["n"]], n)
    expect_equal(d[["p"]], p)
    expect_equal(nrow(d[["data"]]), n)
    expect_equal(ncol(d[["data"]]), ifelse(p == 2, 1, p))
    if (p == 2) {
      expect_true(is.matrix(d[["data"]]))
      expect_equal(d[["avail_stats"]], avail_cir_cv_tests)
    } else {
      expect_true(is.array(d[["data"]]))
      expect_equal(d[["avail_stats"]], avail_sph_cv_tests)
      expect_equal(length(dim(d[["data"]])), 3)
    }

    data_multiple <- r_unif_sph(n = n, p = p, M = 10)
    expect_message(prepare_test_data(data_multiple))

  }

  # Circular structure
  data_cir <- r_unif_cir(n = n)
  d <- prepare_test_data(data_cir)
  expect_equal(d[["n"]], n)
  expect_equal(d[["p"]], 2)
  expect_equal(nrow(d[["data"]]), n)
  expect_equal(ncol(d[["data"]]), 1)
  expect_equal(d[["avail_stats"]], avail_cir_cv_tests)
  expect_true(is.matrix(d[["data"]]))

})

# K-fold statistic and optimal lambda computation
n <- 10
K <- 3
lambda_grid <- list("Poisson" = seq(0.1, 0.9, 0.1),
                    "Softmax" = c(0.1, 1, 3, 10, 20, 60),
                    "Stereo" = seq(-1, 1, 0.25))
for (K in c(2, 3)){
  for (p in c(2:4, 11)) {

    avail_cv_tests <- if (p <= 3) avail_cir_cv_tests else avail_sph_cv_tests

    data <- r_unif_sph(n = n, p = p)
    folds <- k_fold_split(n = n, K = K, seed = 1234)

    # Multiple statistics computation
    # lambda_hat
    suppressWarnings(
      opt_lambda_all <- lambda_hat(data = data,
                                   type = avail_cv_tests,
                                   lambda_grid = lambda_grid, folds = folds,
                                   verbose = FALSE)
    )
    # k_fold
    ## given precomputed folds
    suppressWarnings(
      k_given_fold_all <- k_fold_stat(data = data, type = avail_cv_tests,
                                lambda_grid = lambda_grid, folds = folds,
                                verbose = FALSE)
    )
    ## computing the folds
    suppressWarnings(
      k_fold_all <- k_fold_stat(data = data, type = avail_cv_tests,
                                      lambda_grid = lambda_grid, K = K,
                                      seed = 1234, verbose = FALSE)
    )

    test_that("lambda_hat and k_fold_stat multiple vs. single statistics", {

      # Single statistic computation
      for (type in avail_cv_tests){
        suppressWarnings(
          opt_lambda_single <- lambda_hat(data = data, type = type,
                                          lambda_grid = lambda_grid,
                                          folds = folds,
                                          verbose = FALSE)
        )
        suppressWarnings(
          k_given_fold_stat_single <- k_fold_stat(data = data, type = type,
                                                  lambda_grid = lambda_grid,
                                                  folds = folds,
                                                  verbose = FALSE)
        )
        expect_equal(opt_lambda_single[, type], opt_lambda_all[, type])
        expect_equal(k_given_fold_stat_single$lambda_hat[, type],
                     k_given_fold_all$lambda_hat[, type])
        expect_equal(k_given_fold_stat_single$statistic[, type],
                     k_given_fold_all$statistic[, type])
      }

    })

    test_that("lambda_hat returns correct structure", {

      expect_true(is.matrix(opt_lambda_all))
      expect_equal(nrow(opt_lambda_all), K)
      expect_equal(colnames(opt_lambda_all), avail_cv_tests)

    })

    test_that("k_fold_stat returns correct structure", {

      expect_true(is.data.frame(k_fold_all$statistic))
      expect_equal(nrow(k_fold_all$statistic), K)
      expect_equal(colnames(k_fold_all$statistic), avail_cv_tests)

    })

    test_that("optimal lambda in k_fold_stat vs. lambda_hat", {
      expect_equal(opt_lambda_all, k_given_fold_all$lambda_hat)
    })


    test_that("k_fold_stat with folds = NULL vs. given folds", {
      expect_equal(k_fold_all$statistic, k_given_fold_all$statistic)
    })

  }
}

# HMP p_value
p_val_multiple <- cbind(runif(n = 4), runif(n = 4))
colnames(p_val_multiple) <- c("Poisson", "Softmax")
hmp_multiple <- p_val_hmp(p_val_multiple, M = 1e3)
test_that("p_val_hmp returns correct structure for multiple tests", {

  expect_true(is.matrix(hmp_multiple))
  expect_equal(nrow(hmp_multiple), 1)
  expect_equal(colnames(hmp_multiple), colnames(p_val_multiple))

})

p_val_single <- cbind(runif(n = 6))
colnames(p_val_single) <- c("Poisson")
hmp_single <- p_val_hmp(p_val_single, M = 1e3)
test_that("p_val_hmp returns correct structure for single test", {

  expect_true(is.matrix(hmp_single))
  expect_equal(nrow(hmp_single), 1)
  expect_equal(colnames(hmp_single), colnames(p_val_single))
  expect_error(p_val_hmp(runif(4), M = 1e3))

})

# unif_test_cv
n <- 100
K <- 4
# Hyperspherical
for (p in c(4, 10)) {

  set.seed(1234)
  samp <- r_unif_sph(n = n, p = p)

  # Hide warnings because of precision losses in besselI for high K_max
  suppressWarnings(
    t_asymp <- unif_test_cv(data = samp, type = "all", K = K,
                            p_value = "asymp", seed_fold = 123, K_max = 1000,
                            verbose = FALSE)
  )
  suppressWarnings(
    t_MC <- unif_test_cv(data = samp, type = "all", K = K,
                         p_value = "MC", M = 1e3, seed_fold = 123, K_max = 1000,
                         verbose = FALSE)
  )
  stats_MC <- unif_stat_MC(n = n - round(n / 4),
                           type = c("Softmax", "Poisson", "Stereo"),
                           p = p, M = 1e3, r_H1 = NULL, crit_val = NULL,
                           return_stats = TRUE, stats_sorted = TRUE,
                           Poisson_rho = seq(0.1, 0.9, 0.1),
                           Softmax_kappa = c(0.1, seq(1, 5, 1), seq(10, 30, 5)),
                           Stereo_a = seq(-1, 1, 0.25))$stats_MC
  suppressWarnings(
    t_stats_MC <- unif_test_cv(data = samp, type = "all", K = K,
                               p_value = "MC", stats_MC = stats_MC,
                               K_max = 1000, seed_fold = 123, verbose = FALSE)
  )

  test_that("(hypsph) unif_test_cv computations MC are equal to asymp", {
    for (t in c("Softmax", "Poisson", "Stereo")){
      expect_equal(t_MC[[t]][["fold_statistics"]],
                   t_asymp[[t]][["fold_statistics"]])
      expect_equal(t_MC[[t]][["fold_params"]],
                   t_asymp[[t]][["fold_params"]])
      expect_equal(t_MC[[t]][["fold_p.values"]],
                   t_asymp[[t]][["fold_p.values"]], tolerance = 0.2)
      expect_equal(t_MC[[t]][["reject"]],
                   t_asymp[[t]][["reject"]])
    }
  })

  test_that("(hypsph) unif_test_cv computations are equal to stats_MC = NULL", {
    for (t in c("Softmax", "Poisson", "Stereo")){
      expect_equal(t_MC[[t]][["fold_statistics"]],
                   t_stats_MC[[t]][["fold_statistics"]])
      expect_equal(t_MC[[t]][["fold_params"]],
                   t_stats_MC[[t]][["fold_params"]])
      expect_equal(t_MC[[t]][["fold_p.values"]],
                   t_stats_MC[[t]][["fold_p.values"]], tolerance = 5e-2)
      expect_equal(t_MC[[t]][["reject"]],
                   t_stats_MC[[t]][["reject"]])
    }
  })

  test_that("(hypsph) unif_test_cv with multiple vs. single statistics", {
    for (t in c("Softmax", "Poisson", "Stereo")){
      suppressWarnings(
        t_asymp_single <- unif_test_cv(data = samp, type = t, K = K,
                                p_value = "asymp", seed_fold = 123,
                                K_max = 1000, verbose = FALSE)
      )
      suppressWarnings(
        t_stats_MC_single <- unif_test_cv(data = samp, type = t, K = K,
                                   p_value = "MC", stats_MC = stats_MC,
                                   K_max = 1000, seed_fold = 123,
                                   verbose = FALSE)
      )
    }
    for (output in c("fold_statistics", "fold_params",
                     "fold_p.values", "p.value", "reject")){
      expect_equal(t_asymp[[t]][[output]],
                   t_asymp_single[[output]])
      expect_equal(t_stats_MC[[t]][[output]],
                   t_stats_MC_single[[output]])
    }
  })

  test_that("(hypsph) unif_test_cv returns adequate structure", {

    for (t in c("Softmax", "Poisson", "Stereo")){
      for (output in c("fold_statistics", "fold_params",
                       "fold_p.values")){
        expect_true(is.vector(t_asymp[[t]][[output]]))
        expect_true(is.vector(t_MC[[t]][[output]]))
        expect_true(is.vector(t_stats_MC[[t]][[output]]))

        expect_equal(length(t_asymp[[t]][[output]]), K)
        expect_equal(length(t_MC[[t]][[output]]), K)
        expect_equal(length(t_stats_MC[[t]][[output]]), K)
      }
    }

  })

  test_that("(hypsph) Asymptotic matches MC p-value", {

    expect_equal(t_asymp$Softmax$p.value, t_MC$Softmax$p.value,
                 tolerance = 0.2)
    expect_equal(t_asymp$Poisson$p.value, t_MC$Poisson$p.value,
                 tolerance = 0.2)
    expect_equal(t_asymp$Stereo$p.value, t_MC$Stereo$p.value,
                 tolerance = 0.2)

  })

  test_that("(hypsph) MC p-value equals p-value with stats_MC", {

    expect_equal(t_MC$Softmax$p.value, t_stats_MC$Softmax$p.value,
                 tolerance = 0.1)
    expect_equal(t_MC$Poisson$p.value, t_stats_MC$Poisson$p.value,
                 tolerance = 0.1)
    expect_equal(t_MC$Stereo$p.value, t_stats_MC$Stereo$p.value,
                 tolerance = 0.15)

  })

}

# Circular and spherical
for (p in c(2, 3)) {

  set.seed(1234)
  samp <- r_unif_sph(n = n, p = p)

  suppressWarnings(
    t_asymp <- unif_test_cv(data = samp, type = c("Softmax", "Poisson"), K = K,
                            p_value = "asymp", seed_fold = 123, K_max = 100,
                            verbose = FALSE)
  )
  suppressWarnings(
    t_MC <- unif_test_cv(data = samp, type = c("Softmax", "Poisson"), K = K,
                         p_value = "MC", M = 1e3, seed_fold = 123, K_max = 100,
                         verbose = FALSE)
  )
  stats_MC <- unif_stat_MC(n = n - round(n / 4), type = c("Softmax", "Poisson"),
                           p = p, M = 1e3, r_H1 = NULL, crit_val = NULL,
                           return_stats = TRUE, stats_sorted = TRUE,
                           Poisson_rho = seq(0.1, 0.9, 0.1),
                           Softmax_kappa = c(0.1, seq(1, 5, 1), seq(10, 30, 5)),
                           Stereo_a = seq(-1, 1, 0.25))$stats_MC
  suppressWarnings(
    t_stats_MC <- unif_test_cv(data = samp, type = c("Softmax", "Poisson"),
                               K = K, p_value = "MC", stats_MC = stats_MC,
                               K_max = 100, M = 1e3, seed_fold = 123,
                               verbose = FALSE)
  )

  test_that("(cir-sph) unif_test_cv computations MC are equal to asymp", {
    for (t in c("Softmax", "Poisson")){
      expect_equal(t_MC[[t]][["fold_statistics"]],
                   t_asymp[[t]][["fold_statistics"]])
      expect_equal(t_MC[[t]][["fold_params"]],
                   t_asymp[[t]][["fold_params"]])
      expect_equal(t_MC[[t]][["fold_p.values"]],
                   t_asymp[[t]][["fold_p.values"]], tolerance = 5e-2)
      expect_equal(t_MC[[t]][["reject"]],
                   t_asymp[[t]][["reject"]])
    }
  })

  test_that("(cir-sph)unif_test_cv computations are equal to stats_MC = NULL", {
    for (t in c("Softmax", "Poisson")){
      expect_equal(t_MC[[t]][["fold_statistics"]],
                   t_stats_MC[[t]][["fold_statistics"]])
      expect_equal(t_MC[[t]][["fold_params"]],
                   t_stats_MC[[t]][["fold_params"]])
      expect_equal(t_MC[[t]][["fold_p.values"]],
                   t_stats_MC[[t]][["fold_p.values"]], tolerance = 5e-2)
      expect_equal(t_MC[[t]][["reject"]],
                   t_stats_MC[[t]][["reject"]])
    }
  })

  test_that("(cir-sph) unif_test_cv with multiple vs. single statistics", {
    for (t in c("Softmax", "Poisson")){
      suppressWarnings(
        t_asymp_single <- unif_test_cv(data = samp, type = t, K = K,
                                       p_value = "asymp", seed_fold = 123,
                                       K_max = 100, verbose = FALSE)
      )
      suppressWarnings(
        t_stats_MC_single <- unif_test_cv(data = samp, type = t, K = K,
                                          p_value = "MC", stats_MC = stats_MC,
                                          K_max = 100, seed_fold = 123,
                                          verbose = FALSE)
      )
    }
    for (output in c("fold_statistics", "fold_params",
                     "fold_p.values", "p.value", "reject")){
      expect_equal(t_asymp[[t]][[output]],
                   t_asymp_single[[output]])
      expect_equal(t_stats_MC[[t]][[output]],
                   t_stats_MC_single[[output]])
    }
  })

  test_that("(cir-sph) unif_test_cv returns adequate structure", {

    for (t in c("Softmax", "Poisson")){
      for (output in c("fold_statistics", "fold_params",
                       "fold_p.values")){
        expect_true(is.vector(t_asymp[[t]][[output]]))
        expect_true(is.vector(t_MC[[t]][[output]]))
        expect_true(is.vector(t_stats_MC[[t]][[output]]))

        expect_equal(length(t_asymp[[t]][[output]]), K)
        expect_equal(length(t_MC[[t]][[output]]), K)
        expect_equal(length(t_stats_MC[[t]][[output]]), K)
      }
    }

  })

  test_that("Stereo not available for p in {2,3}", {
    expect_error(unif_test_cv(data = samp, type = c("Stereo"),
                              K = K, p_value = "MC", M = 1e3, seed_fold = 123,
                              verbose = FALSE))
  })

  test_that("(cir-sph) Asymptotic matches MC p-value", {

    expect_equal(t_asymp$Softmax$p.value, t_MC$Softmax$p.value,
                 tolerance = 0.2)
    expect_equal(t_asymp$Poisson$p.value, t_MC$Poisson$p.value,
                 tolerance = 0.2)

  })

  test_that("(cir-sph) MC p-value equals p-value with stats_MC", {

    expect_equal(t_MC$Softmax$p.value, t_stats_MC$Softmax$p.value,
                 tolerance = 0.1)
    expect_equal(t_MC$Poisson$p.value, t_stats_MC$Poisson$p.value,
                 tolerance = 0.1)

  })

}

test_that("If null_variance is given, unif_test_cv checks every statistic
          and every parameter is given", {

    Poisson_grid <- c(0.1, 0.5, 0.8)
    K <- 3
    n_v <- list("Poisson" = null_var_Sobolev(n = 10, p = 3,
                                             type = "Poisson",
                                             lambda_grid = Poisson_grid,
                                             verbose = FALSE)
    )
    expect_error(unif_test_cv(data = r_unif_sph(n = 10, p = 3), K = K,
                              type = c("Poisson", "Softmax"),
                              null_variance = n_v, K_max = 10,
                              Poisson_rho = Poisson_grid))
    expect_error(unif_test_cv(data = r_unif_sph(n = 10, p = 3), K = K,
                              type = c("Poisson"), null_variance = n_v,
                              K_max = 10, Poisson_rho = Poisson_grid[1:2]))
    expect_no_error(unif_test_cv(data = r_unif_sph(n = 10, p = 3),
                                 K = K, type = c("Poisson"),
                                 null_variance = n_v, K_max = 10,
                                 Poisson_rho = Poisson_grid))

})

test_that("null_variance = NULL equals given null_variance", {
  #TODO
})
