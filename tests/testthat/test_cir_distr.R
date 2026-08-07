
test_that("p_cir_stat_Watson is coherent with d_cir_stat_Watson", {

  # The pdf must be the derivative of the cdf, also under the Stephens (1970)
  # modification. Regression: the Stephens branch of the cdf used to apply the
  # Kuiper transform, breaking the coherence.
  n <- 25
  xs <- c(0.05, 0.1, 0.15, 0.2, 0.3)
  h <- 1e-5
  for (stephens in c(FALSE, TRUE)) {
    num <- (p_cir_stat_Watson(xs + h, n = n, Stephens = stephens) -
              p_cir_stat_Watson(xs - h, n = n, Stephens = stephens)) / (2 * h)
    ana <- d_cir_stat_Watson(xs, n = n, Stephens = stephens)
    expect_equal(num, ana, tolerance = 1e-4)
  }

})
