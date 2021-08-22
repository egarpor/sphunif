
n <- 50
set.seed(1242333)
cir_0 <- unif_stat_MC(n = n, M = 5e2, type = "all", p = 2, seeds = 5,
                      verbose = FALSE)
sph_0 <- unif_stat_MC(n = n, M = 5e2, type = "all", p = 3, seeds = 5,
                      verbose = FALSE)
crit_val_bad <- sph_0$crit_val_MC
colnames(crit_val_bad) <- paste0(colnames(crit_val_bad), "_bad")

# Circular
dirs_2 <- r_unif_sph(n = 20, p = 2, M = 1)[, , 1]
cir_pow <- unif_stat_MC(n = n, M = 5e2, p = 2, r_H1 = r_alt, alt = "MvMF",
                        kappa = 0.5, crit_val = cir_0$crit_val_MC,
                        CCF09_dirs = dirs_2, verbose = FALSE)

# Spherical
dirs_3 <- r_unif_sph(n = 20, p = 3, M = 1)[, , 1]
sph_pow <- unif_stat_MC(n = n, M = 5e2, p = 3, r_H1 = r_alt, alt = "MvMF",
                        kappa = 0.5, crit_val = sph_0$crit_val_MC,
                        CCF09_dirs = dirs_3, verbose = FALSE)

test_that("Rejections for MvMF", {

  expect_true(all(sph_pow$power_MC[1, ] > 0.01))
  expect_true(all(cir_pow$power_MC[1, ] > 0.01))

})

test_that("Edge cases", {

  expect_error(unif_stat_MC(n = n, M = 1e2, type = "all", p = 1))
  expect_error(unif_stat_MC(n = n, M = 1e2, type = "all", p = 3,
                            crit_val = sph_0$crit_val_MC[, 1:3],
                            verbose = FALSE))
  expect_error(unif_stat_MC(n = n, M = 1e2, type = "all", p = 3,
                            crit_val = crit_val_bad[, 1:3],
                            verbose = FALSE))
  expect_warning(unif_stat_MC(n = n, M = 1e2, type = "all", p = 5, seeds = 1:3,
                              chunks = 2, verbose = FALSE))
  expect_equal(unif_stat_MC(n = n, M = 1e2, type = c("PAD", "Ajne"), p = 3,
                            seeds = 5, verbose = FALSE)$stats$PAD,
               unif_stat_MC(n = n, M = 1e2, type = "PAD", p = 3,
                            seeds = 5, verbose = FALSE,
                            crit_val = sph_0$crit_val_MC)$stats$PAD)

})

test_that("Several options", {

  expect_equal(unif_stat_MC(n = n, M = 10, type = "all", p = 2,
                            return_stats = FALSE, verbose = FALSE)$stats, NA)
  expect_equal(unif_stat_MC(n = n, M = 10, type = "all", p = 3,
                            return_stats = FALSE, verbose = FALSE)$stats, NA)
  expect_equal(unname(as.matrix(unif_stat_MC(n = n, M = 10, type = "all", p = 5,
                                             verbose = FALSE,
                                             stats_sorted = TRUE,
                                             seeds = 1)$stats)),
               sort_each_col(as.matrix(unif_stat_MC(n = n, M = 10, type = "all",
                                                    p = 5, verbose = FALSE,
                                                    seeds = 1)$stats)))
  set.seed(1)
  CCF09_dirs <- r_unif_sph(n = 50, p = 3, M = 1)[, , 1]
  expect_equal(unif_stat_MC(n = n, M = 5, type = "CCF09", p = 3,
                            verbose = FALSE, stats_sorted = TRUE, seeds = 1,
                            chunks = 1)$CCF09,
               unif_stat_MC(n = n, M = 5, type = "CCF09", p = 3,
                            verbose = FALSE, stats_sorted = TRUE, seeds = 1,
                            chunks = 1, CCF09_dirs = CCF09_dirs)$CCF09)
  expect_true(!is.null(capture.output(
    a <- unif_stat_MC(n = n, M = 10, type = "all", p = 2, verbose = TRUE)
  )))

})
