
# Clean workspace
rm(list = ls())

# Load packages
library(sphunif)
library(rotasym)
stopifnot(packageVersion("sphunif") >= "1.3.0")

# Maximum likelihood estimation of rotationally symmetric models like
# f(x) = c(kappa) * exp(kappa * (x'mu)^b)
mle_f_b <- function(data, b = 1, N = 5120) {

  # Data dimensions
  n <- nrow(data)
  p <- ncol(data)

  # Use Gauss--Legendre quadrature to integrate the normalizing constant,
  # saving the common terms
  z_k <- drop(Gauss_Legen_nodes(a = -1, b = 1, N = N))
  w_k <- drop(Gauss_Legen_weights(a = -1, b = 1, N = N))
  omega_p <- w_p(p = p - 1)
  z_k_b <- z_k^b
  common_term <- omega_p * w_k * (1 - z_k^2)^((p - 3) / 2)
  log_const <- function(kappa) {

    const <- sum(common_term * exp(kappa * z_k_b), na.rm = TRUE)
    return(-log(const))

  }

  # Minus log-likelihood
  minus_loglik <- function(param) {

    kappa <- abs(param[1])
    mu <- param[-1]
    mu <- mu / sqrt(sum(mu^2))
    return(-(n * log_const(kappa) + kappa * sum((data %*% mu)^b)))

  }

  # Specific approaches
  if (b == 1) {

    kappa <- DirStats::kappa_ml(data)
    mu <- DirStats::mu_ml(data)
    log_const <- function(kappa) c_vMF(p = p, kappa = kappa, log = TRUE)
    return(list("kappa" = kappa, "mu" = mu,
                "mll" = minus_loglik(c(kappa, mu))))

  }

  # Starting values
  start_vmf <- c(DirStats::kappa_ml(data), DirStats::mu_ml(data))
  start_vmf <- rbind(start_vmf, 0.1 * start_vmf, 10 * start_vmf)
  eig_vecs <- eigen(crossprod(data))$vectors
  start_wat <- rbind(c(3, eig_vecs[, 1]),
                     c(3, eig_vecs[, p]),
                     c(3, spherical_loc_PCA(data)))

  # Optimization
  opt <- sdetorus::mleOptimWrapper(minusLogLik = minus_loglik,
                                   start = rbind(start_vmf, start_wat),
                                   optMethod = "nlm", penalty = 0,
                                   selectSolution = "lowest")
  ind <- unlist(lapply(opt$solutionsOutput, function(x) x$minimum))
  ind[ind == 0] <- Inf
  ind_min <- which.min(ind)
  opt$solutionsOutput[[ind_min]]$estimate[1] <-
    abs(opt$solutionsOutput[[ind_min]]$estimate[1])
  opt$solutionsOutput[[ind_min]]$estimate[-1] <-
    opt$solutionsOutput[[ind_min]]$estimate[-1] /
    sqrt(sum(opt$solutionsOutput[[ind_min]]$estimate[-1]^2))
  result <- list("kappa" = opt$solutionsOutput[[ind_min]]$estimate[1],
                 "mu" = opt$solutionsOutput[[ind_min]]$estimate[-1],
                 "mll" = opt$solutionsOutput[[ind_min]]$minimum)
  # opt$par[1] <- abs(opt$par[1])
  # opt$par[-1] <- opt$par[-1] / sqrt(sum(opt$par[-1]^2))
  # result <- list("kappa" = opt$par[1], "mu" = opt$par[-1], "mll" = opt$value)
  message(paste0("kappa = ", result$kappa, ", mu = ",
                 paste(result$mu, collapse = ", "),
                 ", mll = ", result$mll, "."))
  return(result)

}

# Simulation of rotationally symmetric models like
# f(x) = c(kappa) * exp(kappa * (x'mu)^b)
r_f_b <- function(n, mu = c(0, 0, 1), kappa = 1, b = 1) {

  samp <- sphunif::r_locdev(n = n, mu = mu, f = function(t) exp(kappa * t^b),
                            kappa = 1)
  return(samp)

}

# Computes the local power for angular function f(s) = exp(s^b) of a finite
# Sobolev test with a single non-null weight at position k_v
power_exp_b <- function(tau, k_v, b, p, alpha = 0.05) {

  # Critical value
  df_k_v <- d_p_k(p = p, k = k_v)
  c_alpha <- qchisq(1 - alpha, df = df_k_v)

  # Derivatives
  if (b == 1) {

    f_k_v_0 <- 1

  } else if (b %in% 2:10) {

    # f^(k)(0) for f(s) = exp(s^b), for k = 1, ..., 10 (rows) and
    # b = 2, ..., 10 (columns)
    f_k_b <- matrix(c(
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      2, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 6, 0, 0, 0, 0, 0, 0, 0,
      12, 0, 24, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 120, 0, 0, 0, 0, 0,
      120, 360, 0, 0, 720, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 5040, 0, 0, 0,
      1680, 0, 20160, 0, 0, 0, 40320, 0, 0,
      0, 60480, 0, 0, 0, 0, 0, 362880, 0,
      30240, 0, 0, 1814400, 0, 0, 0, 0, 3628800
    ), nrow = 10, ncol = 9, byrow = TRUE)

    # Access for (k_v, b)
    f_k_v_0 <- f_k_b[k_v, b - 1]

  } else {

    stop("b must be 1, 2, 3, 4, 5, or 6!")

  }

  # Non-centrality parameter
  ncp <- df_k_v * f_k_v_0^2 * tau^(2 * k_v) / prod(p + 2 * (0:(k_v - 1)))^2

  # Power
  return(1 - pchisq(q = c_alpha, df = df_k_v, ncp = ncp))

}

# Computes the power of a finite Sobolev test with a single non-null weight at
# position k_v against a model f(x) = c(kappa) * exp(kappa * (x'mu)^b) where
# (kappa, mu) are estimated from the data
mle_power_exp_b <- function(data, k_v = 1:3, b = 1:3, N = 5120, alpha = 0.05,
                            sort_lik = TRUE) {

  # Data dimensions
  n <- nrow(data)
  p <- ncol(data)

  # kappa_n = n^{-1 / (2 * k_v)} * tau
  kappa_n <- sapply(b, function(bi) {
    mle <- mle_f_b(data = data, b = bi, N = N)
    return(c(mle$kappa, -mle$mll)) # Return loglik instead of -loglik
  })
  tau_n <- n^(1 / (2 * k_v)) %o% kappa_n[1, ]

  # Power
  power <- matrix(nrow = length(k_v), ncol = length(b))
  for (ki in seq_along(k_v)) {
    for (bi in seq_along(b)) {

      power[ki, bi] <- power_exp_b(tau = tau_n[ki, bi], k_v = k_v[ki],
                                   b = b[bi], p = p, alpha = alpha)

    }
  }
  power <- rbind(power, kappa_n[2, ])
  rownames(power) <- c(paste0("k=", k_v), "loglik")
  colnames(power) <- paste0("b=", b)

  # Sort according to likelihood?
  if (sort_lik) {

    power <- power[, order(kappa_n[2, ], decreasing = TRUE)]

  }
  return(power)

}

# Examples
run <- FALSE
if (run) {

  library(testthat)
  library(rgl)

  # Check MLE
  mu <- c(0, 0, 1)
  kappa <- 3
  for (b in 1:6) {
    samp <- r_f_b(n = 1000, mu = mu, kappa = kappa, b = b)
    est <- mle_f_b(data = samp, b = b)
    test_that("Estimation working properly", {
      expect_lt(abs(est$kappa - kappa), 5e-1)
      expect_lt(sqrt(sum((abs(est$mu) - abs(mu))^2)), 5e-1)
    })
  }
  plot3d(samp, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1))

  # Gives blue curve in Figure 2
  tau <- seq(0, 4, by = 0.1)
  plot(tau, power_exp_b(tau = tau, k_v = 2, b = 2, p = 3), ylim = c(0, 1))

  # Check powers
  p <- 3
  data <- r_alt(n = 500, p = p, alt = "vMF", kappa = 1)[, , 1]
  mle_power_exp_b(data = data, k_v = 1:3, b = 1:3)

  data <- r_alt(n = 500, p = p, alt = "W", kappa = 2)[, , 1]
  mle_power_exp_b(data = data, k_v = 1:3, b = 1:3)

  data <- r_alt(n = 500, p = p, alt = "MvMF", kappa = 2)[, , 1]
  mle_power_exp_b(data = data, k_v = 1:3, b = 1:3)

}

### Uniformity of comets

# Load data
data("comets", package = "sphunif")

# Add normal
comets$normal <- cbind(sin(comets$i) * sin(comets$om),
                       -sin(comets$i) * cos(comets$om),
                       cos(comets$i))

# Oort cloud (long-period comets)
comets_oort <- subset(x = comets,
                      subset = !(class %in% c("HYP", "PAR")) & per_y >= 200)
nrow(comets_oort)

# Comets without fragments
comets_oort_clean <- comets_oort[!comets_oort$frag, ]
nrow(comets_oort_clean)

# Kuiper belt (short-period comets)
comets_kuiper <- subset(x = comets, subset = !(class %in% c("HYP", "PAR")) &
                          per_y < 200)
nrow(comets_kuiper)

# Cleaned data
comets_kuiper_clean <- comets_kuiper[!comets_kuiper$frag, ]
nrow(comets_kuiper_clean)

## Tests Oort

# Test uniformity with all records
test_oort <-
  unif_test(data = comets_oort$normal, type = c("Sobolev", "PAD"),
            Sobolev_vk2 = diag(1, nrow = 6), p_value = "asymp")

# Test uniformity with clean records
test_oort_clean <-
  unif_test(data = comets_oort_clean$normal, type = c("Sobolev", "PAD"),
            Sobolev_vk2 = diag(1, nrow = 6), p_value = "asymp")

# p-values
round(unlist(lapply(test_oort, function(test) test$p.value)), 4)
round(unlist(lapply(test_oort_clean, function(test) test$p.value)), 4)

## Tests Kuiper

# Test uniformity with all records
test_kuiper <-
  unif_test(data = comets_kuiper$normal, type = c("Sobolev", "PAD"),
            Sobolev_vk2 = diag(1, nrow = 6), p_value = "asymp")

# Test uniformity with clean records
test_kuiper_clean <-
  unif_test(data = comets_kuiper_clean$normal, type = c("Sobolev", "PAD"),
            Sobolev_vk2 = diag(1, nrow = 6), p_value = "asymp")

# p-values
round(unlist(lapply(test_kuiper, function(test) test$p.value)), 4)
round(unlist(lapply(test_kuiper_clean, function(test) test$p.value)), 4)

## Power computation

pow_oort <- round(mle_power_exp_b(data = comets_oort$normal,
                                  k_v = 1:6, b = 1:6), 4)
pow_oort_clean <- round(mle_power_exp_b(data = comets_oort_clean$normal,
                                        k_v = 1:6, b = 1:6), 4)
pow_kuiper <- round(mle_power_exp_b(data = comets_kuiper$normal,
                                    k_v = 1:6, b = 1:6), 4)
pow_kuiper_clean <- round(mle_power_exp_b(data = comets_kuiper_clean$normal,
                                          k_v = 1:6, b = 1:6), 4)
pow_oort_clean
apply(pow_oort, 1, max)
apply(pow_oort_clean, 1, max)
apply(pow_kuiper, 1, max)
apply(pow_kuiper_clean, 1, max)

### Rotational symmetry of comets

# Obtain the signs with respect to the north pole c(0, 0, 1) -- these are the
# positions in S^1 of the longitudes associated to the data. The latitudes are
# disregarded to investigate if the longitudes are uniformly distributed, which
# is implied if rotational symmetry about (0, 0, 1) holds.
signs_oort <- signs(comets_oort$normal, theta = c(0, 0, 1))
signs_oort_clean <- signs(comets_oort_clean$normal, theta = c(0, 0, 1))
signs_kuiper <- signs(comets_kuiper$normal, theta = c(0, 0, 1))
signs_kuiper_clean <- signs(comets_kuiper_clean$normal, theta = c(0, 0, 1))

## Tests Oort

# Test uniformity with all records
test_rot_oort <-
  unif_test(data = signs_oort, type = c("Sobolev", "PAD"),
            Sobolev_vk2 = diag(1, nrow = 6), p_value = "asymp")

# Test uniformity with clean records
test_rot_oort_clean <-
  unif_test(data = signs_oort_clean, type = c("Sobolev", "PAD"),
            Sobolev_vk2 = diag(1, nrow = 6), p_value = "asymp")

# p-values
round(unlist(lapply(test_rot_oort, function(test) test$p.value)), 4)
round(unlist(lapply(test_rot_oort_clean, function(test) test$p.value)), 4)

# Unspecified theta
test_rotasym(comets_oort$normal, type = "loc_vMF")
test_rotasym(comets_oort$normal, type = "sc")
test_rotasym(comets_oort_clean$normal, type = "loc_vMF")
test_rotasym(comets_oort_clean$normal, type = "sc")

## Tests Kuiper

# Test uniformity with all records
test_rot_kuiper <-
  unif_test(data = signs_kuiper, type = c("Sobolev", "PAD"),
            Sobolev_vk2 = diag(1, nrow = 6), p_value = "asymp")

# Test uniformity with clean records
test_rot_kuiper_clean <-
  unif_test(data = signs_kuiper_clean, type = c("Sobolev", "PAD"),
            Sobolev_vk2 = diag(1, nrow = 6), p_value = "asymp")

# p-values
round(unlist(lapply(test_rot_kuiper, function(test) test$p.value)), 4)
round(unlist(lapply(test_rot_kuiper_clean, function(test) test$p.value)), 4)

# Unspecified theta
test_rotasym(comets_kuiper$normal, type = "loc_vMF")
test_rotasym(comets_kuiper$normal, type = "sc")
test_rotasym(comets_kuiper_clean$normal, type = "loc_vMF")
test_rotasym(comets_kuiper_clean$normal, type = "sc")

## Power computation

pow_rot_oort <- round(mle_power_exp_b(data = signs_oort,
                                      k_v = 1:6, b = 1:6), 4)
pow_rot_oort_clean <- round(mle_power_exp_b(data = signs_oort_clean,
                                            k_v = 1:6, b = 1:6), 4)
pow_rot_kuiper <- round(mle_power_exp_b(data = signs_kuiper,
                                        k_v = 1:6, b = 1:6), 4)
pow_rot_kuiper_clean <- round(mle_power_exp_b(data = signs_kuiper_clean,
                                              k_v = 1:6, b = 1:6), 4)
pow_rot_oort_clean
apply(pow_rot_oort, 1, max)
apply(pow_rot_oort_clean, 1, max)
apply(pow_rot_kuiper, 1, max)
apply(pow_rot_kuiper_clean, 1, max)
