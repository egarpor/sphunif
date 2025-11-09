

#' @title Asymptotic distributions for spherical uniformity statistics
#'
#' @description Computation of the asymptotic null distributions of
#' spherical uniformity statistics.
#'
#' @inheritParams cir_stat_distr
#' @inheritParams Sobolev
#' @inheritParams r_unif
#' @inheritParams wschisq
#' @inheritParams unif_stat_distr
#' @param vk2 weights for the finite Sobolev test. A non-negative vector or
#' matrix. Defaults to \code{c(0, 0, 1)}.
#' @param regime type of asymptotic regime for the CJ12 test, either \code{1}
#' (sub-exponential regime), \code{2} (exponential), or \code{3}
#' (super-exponential; default).
#' @param beta \eqn{\beta} parameter in the exponential regime of the CJ12
#' test, a non-negative real. Defaults to \code{0}.
#' @inheritParams cir_stat
#' @inheritParams unif_stat
#' @param ... further parameters passed to \code{\link{p_Sobolev}} or
#' \code{\link{d_Sobolev}} (such as \code{x_tail}).
#' @return
#' \itemize{
#'   \item \code{r_sph_stat_*}: a matrix of size \code{c(n, 1)} containing
#'   the sample.
#'   \item \code{p_sph_stat_*}, \code{d_sph_stat_*}: a matrix of size
#'   \code{c(nx, 1)} with the evaluation of the distribution or density
#'   functions at \code{x}.
#' }
#' @details
#' Descriptions and references on most of the asymptotic distributions
#' are available in García-Portugués and Verdebout (2018).
#' @examples
#' # Ajne
#' curve(d_sph_stat_Ajne(x, p = 3, method = "HBE"), n = 2e2, ylim = c(0, 4))
#' curve(p_sph_stat_Ajne(x, p = 3, method = "HBE"), n = 2e2, col = 2,
#'       add = TRUE)
#'
#' # Bakshaev
#' curve(d_sph_stat_Bakshaev(x, p = 3, method = "HBE"), to = 5, n = 2e2,
#'       ylim = c(0, 2))
#' curve(p_sph_stat_Bakshaev(x, p = 3, method = "HBE"), n = 2e2, col = 2,
#'       add = TRUE)
#'
#' # Bingham
#' curve(d_sph_stat_Bingham(x, p = 3), to = 20, n = 2e2, ylim = c(0, 1))
#' curve(p_sph_stat_Bingham(x, p = 3), n = 2e2, col = 2, add = TRUE)
#'
#' # CJ12
#' curve(d_sph_stat_CJ12(x, regime = 1), from = -10, to = 10, n = 2e2,
#'       ylim = c(0, 1))
#' curve(d_sph_stat_CJ12(x, regime = 2, beta = 0.1), n = 2e2, col = 2,
#'       add = TRUE)
#' curve(d_sph_stat_CJ12(x, regime = 3), n = 2e2, col = 3, add = TRUE)
#' curve(p_sph_stat_CJ12(x, regime = 1), n = 2e2, col = 1, add = TRUE)
#' curve(p_sph_stat_CJ12(x, regime = 2, beta = 0.1), n = 2e2, col = 2,
#'       add = TRUE)
#' curve(p_sph_stat_CJ12(x, regime = 3), col = 3, add = TRUE)
#'
#' # Gine Fn
#' curve(d_sph_stat_Gine_Fn(x, p = 3, method = "HBE"), to = 2, n = 2e2,
#'       ylim = c(0, 2))
#' curve(p_sph_stat_Gine_Fn(x, p = 3, method = "HBE"), n = 2e2, col = 2,
#'       add = TRUE)
#'
#' # Gine Gn
#' curve(d_sph_stat_Gine_Gn(x, p = 3, method = "HBE"), to = 1.5, n = 2e2,
#'       ylim = c(0, 2.5))
#' curve(p_sph_stat_Gine_Gn(x, p = 3, method = "HBE"), n = 2e2, col = 2,
#'       add = TRUE)
#'
#' # PAD
#' curve(d_sph_stat_PAD(x, p = 3, method = "HBE"), to = 3, n = 2e2,
#'       ylim = c(0, 1.5))
#' curve(p_sph_stat_PAD(x, p = 3, method = "HBE"), n = 2e2, col = 2,
#'       add = TRUE)
#'
#' # PCvM
#' curve(d_sph_stat_PCvM(x, p = 3, method = "HBE"), to = 0.6, n = 2e2,
#'       ylim = c(0, 7))
#' curve(p_sph_stat_PCvM(x, p = 3, method = "HBE"), n = 2e2, col = 2,
#'       add = TRUE)
#'
#' # Poisson
#' curve(d_sph_stat_Poisson(x, p = 3, method = "HBE"), to = 2, n = 2e2,
#'       ylim = c(0, 2))
#' curve(p_sph_stat_Poisson(x, p = 3, method = "HBE"), n = 2e2, col = 2,
#'       add = TRUE)
#'
#' # PRt
#' curve(d_sph_stat_PRt(x, p = 3, method = "HBE"), n = 2e2, ylim = c(0, 5))
#' curve(p_sph_stat_PRt(x, p = 3, method = "HBE"), n = 2e2, col = 2, add = TRUE)
#'
#' # Rayleigh
#' curve(d_sph_stat_Rayleigh(x, p = 3), to = 15, n = 2e2, ylim = c(0, 1))
#' curve(p_sph_stat_Rayleigh(x, p = 3), n = 2e2, col = 2, add = TRUE)
#'
#' # HD-standardized Rayleigh
#' curve(d_sph_stat_Rayleigh_HD(x, p = 3), from = -4, to = 4, n = 2e2,
#'       ylim = c(0, 1))
#' curve(p_sph_stat_Rayleigh_HD(x, p = 3), n = 2e2, col = 2, add = TRUE)
#'
#' # Riesz
#' curve(d_sph_stat_Riesz(x, p = 3, method = "HBE"), n = 2e2, from = 0, to = 5,
#'       ylim = c(0, 2))
#' curve(p_sph_stat_Riesz(x, p = 3, method = "HBE"), n = 2e2, col = 2,
#'       add = TRUE)
#'
#' # Sobolev
#' x <- seq(-1, 5, by = 0.05)
#' vk2 <- diag(rep(0.3, 2))
#' matplot(x, d_sph_stat_Sobolev(x = x, vk2 = vk2, p = 3), type = "l",
#'         ylim = c(0, 1), lty = 1)
#' matlines(x, p_sph_stat_Sobolev(x = x, vk2 = vk2, p = 3), lty = 1)
#' matlines(x, d_sph_stat_Sobolev(x = x, vk2 = vk2 + 0.01, p = 3), lty = 2)
#' matlines(x, p_sph_stat_Sobolev(x = x, vk2 = vk2 + 0.01, p = 3), lty = 2)
#'
#' # Softmax
#' curve(d_sph_stat_Softmax(x, p = 3, method = "HBE"), to = 2, n = 2e2,
#'       ylim = c(0, 2))
#' curve(p_sph_stat_Softmax(x, p = 3, method = "HBE"), n = 2e2, col = 2,
#'       add = TRUE)
#'
#' # Stereo
#' curve(d_sph_stat_Stereo(x, p = 4, method = "HBE"), from=-5,to = 10, n = 2e2,
#'       ylim = c(0, 2))
#' curve(p_sph_stat_Stereo(x, p = 4, method = "HBE"), n = 2e2, col = 2,
#'       add = TRUE)
#' @name sph_stat_distr
NULL


#' @rdname sph_stat_distr
#' @export
p_sph_stat_Ajne <- function(x, p, K_max = 1e3, thre = 0, method = "I", ...) {

  cbind(p_Sobolev(x = x, p = p, type = "Ajne", K_max = K_max, thre = thre,
                  method = method, ...))

}


#' @rdname sph_stat_distr
#' @export
d_sph_stat_Ajne <- function(x, p, K_max = 1e3, thre = 0, method = "I", ...) {

  cbind(d_Sobolev(x = x, p = p, type = "Ajne", K_max = K_max, thre = thre,
                  method = method, ...))

}


#' @rdname sph_stat_distr
#' @export
p_sph_stat_Bakshaev <- function(x, p, K_max = 1e3, thre = 0, method = "I",
                                ...) {

  cbind(p_Sobolev(x = x, p = p, type = "Bakshaev", K_max = K_max, thre = thre,
                  method = method, ...))

}


#' @rdname sph_stat_distr
#' @export
d_sph_stat_Bakshaev <- function(x, p, K_max = 1e3, thre = 0, method = "I",
                                ...) {

  cbind(d_Sobolev(x = x, p = p, type = "Bakshaev", K_max = K_max, thre = thre,
                  method = method, ...))

}


#' @rdname sph_stat_distr
#' @export
p_sph_stat_Gine_Fn <- function(x, p, K_max = 1e3, thre = 0, method = "I", ...) {

  cbind(p_Sobolev(x = x, p = p, type = "Gine_Fn", K_max = K_max, thre = thre,
                  method = method, ...))

}


#' @rdname sph_stat_distr
#' @export
d_sph_stat_Gine_Fn <- function(x, p, K_max = 1e3, thre = 0, method = "I", ...) {

  cbind(d_Sobolev(x = x, p = p, type = "Gine_Fn", K_max = K_max, thre = thre,
                  method = method, ...))

}


#' @rdname sph_stat_distr
#' @export
p_sph_stat_Gine_Gn <- function(x, p, K_max = 1e3, thre = 0, method = "I", ...) {

  cbind(p_Sobolev(x = x, p = p, type = "Gine_Gn", K_max = K_max, thre = thre,
                  method = method, ...))

}


#' @rdname sph_stat_distr
#' @export
d_sph_stat_Gine_Gn <- function(x, p, K_max = 1e3, thre = 0, method = "I", ...) {

  cbind(d_Sobolev(x = x, p = p, type = "Gine_Gn", K_max = K_max, thre = thre,
                  method = method, ...))

}


#' @rdname sph_stat_distr
#' @export
p_sph_stat_PAD <- function(x, p, K_max = 1e3, thre = 0, method = "I", ...) {

  cbind(p_Sobolev(x = x, p = p, type = "PAD", K_max = K_max, thre = thre,
                  method = method, ...))

}


#' @rdname sph_stat_distr
#' @export
d_sph_stat_PAD <- function(x, p, K_max = 1e3, thre = 0, method = "I", ...) {

  cbind(d_Sobolev(x = x, p = p, type = "PAD", K_max = K_max, thre = thre,
                  method = method, ...))

}


#' @rdname sph_stat_distr
#' @export
p_sph_stat_PCvM <- function(x, p, K_max = 1e3, thre = 0, method = "I", ...) {

  cbind(p_Sobolev(x = x, p = p, type = "PCvM", K_max = K_max, thre = thre,
                  method = method, ...))

}


#' @rdname sph_stat_distr
#' @export
d_sph_stat_PCvM <- function(x, p, K_max = 1e3, thre = 0, method = "I", ...) {

  cbind(d_Sobolev(x = x, p = p, type = "PCvM", K_max = K_max, thre = thre,
                  method = method, ...))

}


#' @rdname sph_stat_distr
#' @export
p_sph_stat_Poisson <- function(x, p, rho = 0.5, K_max = 1e3, thre = 0,
                               method = "I", ...) {

  # psi_tilde_0 is the shift from U-stat to V-stat
  b_0 <- (1 - rho)^(p - 1) / (1 + rho)
  psi_tilde_0 <- ((1 - rho) / sqrt(1 - 2 * rho + rho^2))^p - b_0
  cbind(p_Sobolev(x = x + psi_tilde_0, p = p, type = "Poisson",
                  Poisson_rho = rho, K_max = K_max, thre = thre,
                  method = method, ...))

}


#' @rdname sph_stat_distr
#' @export
d_sph_stat_Poisson <- function(x, p, rho = 0.5, K_max = 1e3, thre = 0,
                               method = "I", ...) {

  # psi_tilde_0 is the shift from U-stat to V-stat
  b_0 <- (1 - rho)^(p - 1) / (1 + rho)
  psi_tilde_0 <- ((1 - rho) / sqrt(1 - 2 * rho + rho^2))^p - b_0
  cbind(d_Sobolev(x = x + psi_tilde_0, p = p, type = "Poisson",
                  Poisson_rho = rho, K_max = K_max, thre = thre,
                  method = method, ...))

}


#' @rdname sph_stat_distr
#' @export
p_sph_stat_PRt <- function(x, p, t = 1 / 3, K_max = 1e3, thre = 0,
                           method = "I", ...) {

  cbind(p_Sobolev(x = x, p = p, type = "PRt", Rothman_t = t,
                  K_max = K_max, thre = thre, method = method, ...))

}


#' @rdname sph_stat_distr
#' @export
d_sph_stat_PRt <- function(x, p, t = 1 / 3, K_max = 1e3, thre = 0,
                           method = "I", ...) {

  cbind(d_Sobolev(x = x, p = p, type = "PRt", Rothman_t = t,
                  K_max = K_max, thre = thre, method = method, ...))

}


#' @rdname sph_stat_distr
#' @export
p_sph_stat_Riesz <- function(x, p, s = 1, K_max = 1e3, thre = 0,  method = "I",
                             ...) {

  cbind(p_Sobolev(x = x, p = p, type = "Riesz", Riesz_s = s, K_max = K_max,
                  thre = thre, method = method, ...))

}


#' @rdname sph_stat_distr
#' @export
d_sph_stat_Riesz <- function(x, p, s = 1, K_max = 1e3, thre = 0, method = "I",
                             ...) {

  cbind(d_Sobolev(x = x, p = p, type = "Riesz", Riesz_s = s, K_max = K_max,
                  thre = thre, method = method, ...))

}


#' @rdname sph_stat_distr
#' @export
p_sph_stat_Sobolev <- function(x, p, vk2 = c(0, 0, 1), method = "I", ...) {

  # As a matrix
  if (is.vector(x)) {

    x <- matrix(x, ncol = 1)

  }

  # Check x and vk2 are compatible
  vk2 <- rbind(vk2)
  n_vk2 <- nrow(vk2)
  n_col_x <- ncol(x)
  if (ncol(x) != n_vk2) {

    if (n_col_x == 1) {

      x <- matrix(x, nrow = nrow(x), ncol = n_vk2, byrow = FALSE)

    } else {

      stop("The number of columns of x must match the number of rows of vk2.")

    }

  }

  # Loop for different vk2's
  cdf <- matrix(nrow = nrow(x), ncol = n_vk2)
  for (j in seq_len(n_vk2)) {

    # Check if the asymptotic distribution is a single chi-squared
    nonzero_vk2_j <- which(vk2[j, ] != 0)
    if (length(nonzero_vk2_j) == 1) {

      cdf[, j] <- pchisq(q = x[, j] / vk2[j, nonzero_vk2_j],
                         df = d_p_k(p = p, k = nonzero_vk2_j))

    } else {

      cdf[, j] <- p_Sobolev(x = x[, j], p = p, type = "Sobolev",
                            Sobolev_vk2 = vk2[j, ], K_max = max(nonzero_vk2_j),
                            thre = 0, method = method, ...)

    }

  }
  return(cdf)

}


#' @rdname sph_stat_distr
#' @export
d_sph_stat_Sobolev <- function(x, p, vk2 = c(0, 0, 1), method = "I", ...) {

  # As a matrix
  if (is.vector(x)) {

    x <- matrix(x, ncol = 1)

  }

  # Check x and vk2 are compatible
  vk2 <- rbind(vk2)
  n_vk2 <- nrow(vk2)
  n_col_x <- ncol(x)
  if (ncol(x) != n_vk2) {

    if (n_col_x == 1) {

      x <- matrix(x, nrow = nrow(x), ncol = n_vk2, byrow = FALSE)

    } else {

      stop("The number of columns of x must match the number of rows of vk2.")

    }

  }

  # Loop for different vk2's
  pdf <- matrix(nrow = nrow(x), ncol = n_vk2)
  for (j in seq_len(n_vk2)) {

    # Check if the asymptotic distribution is a single chi-squared
    nonzero_vk2_j <- which(vk2[j, ] != 0)
    if (length(nonzero_vk2_j) == 1) {

      pdf[, j] <- dchisq(x = x[, j] / vk2[j, nonzero_vk2_j],
                         df = d_p_k(p = p, k = nonzero_vk2_j)) /
        vk2[j, nonzero_vk2_j]

    } else {

      pdf[, j] <- d_Sobolev(x = x[, j], p = p, type = "Sobolev",
                            Sobolev_vk2 = vk2[j, ], K_max = max(nonzero_vk2_j),
                            thre = 0, method = method, ...)

    }

  }
  return(pdf)

}


#' @rdname sph_stat_distr
#' @export
p_sph_stat_Softmax <- function(x, p, kappa = 1, K_max = 1e3, thre = 0,
                               method = "I", ...) {

  # psi_tilde_0 is the shift from U-stat to V-stat
  b_0 <- exp(switch((p > 2) + 1,
                    log(besselI(x = kappa, nu = 0, expon.scaled = TRUE)),
                    (0.5 * p - 1) * (log(2) - log(kappa)) +
                      lgamma(0.5 * p - 1) + log(0.5 * p - 1) +
                      log(besselI(x = kappa, nu = 0.5 * p - 1,
                                  expon.scaled = TRUE))))
  psi_tilde_0 <- 1 - b_0
  cbind(p_Sobolev(x = x + psi_tilde_0, p = p, type = "Softmax",
                  Softmax_kappa = kappa, K_max = K_max, thre = thre,
                  method = method, ...))

}


#' @rdname sph_stat_distr
#' @export
d_sph_stat_Softmax <- function(x, p, kappa = 1, K_max = 1e3, thre = 0,
                               method = "I", ...) {

  # psi_tilde_0 is the shift from U-stat to V-stat
  b_0 <- exp(switch((p > 2) + 1,
                    log(besselI(x = kappa, nu = 0, expon.scaled = TRUE)),
                    (0.5 * p - 1) * (log(2) - log(kappa)) +
                      lgamma(0.5 * p - 1) + log(0.5 * p - 1) +
                      log(besselI(x = kappa, nu = 0.5 * p - 1,
                                  expon.scaled = TRUE))))
  psi_tilde_0 <- 1 - b_0
  cbind(d_Sobolev(x = x + psi_tilde_0, p = p, type = "Softmax",
                  Softmax_kappa = kappa, K_max = K_max, thre = thre,
                  method = method, ...))

}


#' @rdname sph_stat_distr
#' @export
p_sph_stat_Stein <- function(x, p, Stein_K = 10, Stein_cf = FALSE, method = "I",
                             ...) {

  Stein_vk2 <- weights_dfs_Sobolev(p = p, K_max = Stein_K, thre = 0,
                                   type = "Stein", Stein_cf = Stein_cf,
                                   verbose = FALSE)$weights
  p_sph_stat_Sobolev(x = x, p = p, vk2 = Stein_vk2, method = method, ...)

}


#' @rdname sph_stat_distr
#' @export
d_sph_stat_Stein <- function(x, p, Stein_K = 10, Stein_cf = FALSE,
                             method = "I", ...) {

  Stein_vk2 <- weights_dfs_Sobolev(p = p, K_max = Stein_K, thre = 0,
                                   type = "Stein", Stein_cf = Stein_cf,
                                   verbose = FALSE)$weights
  d_sph_stat_Sobolev(x = x, p = p, vk2 = Stein_vk2, method = method, ...)

}


#' @rdname sph_stat_distr
#' @export
p_sph_stat_Stereo <- function(x, p, a = 0, K_max = 1e3, method = "I", ...) {

  if (p <= 3) {

    stop("The asymptotic distribution is only available for p > 3.")

  }

  # psi_tilde_0 is the shift from U-stat to V-stat
  v_k2 <- weights_dfs_Sobolev(p = p, K_max = K_max, thre = 0, type = "Stereo",
                              Stereo_a = a, ...)$weights
  psi_tilde_0 <- Gegen_series(theta = 0, coefs = vk2_to_bk(vk2 = v_k2, p = p),
                              k = 1:K_max, p = p)
  cbind(p_Sobolev(x = x + psi_tilde_0, p = p, type = "Stereo", Stereo_a = a,
                  K_max = K_max, thre = 0, method = method, ...))

}


#' @rdname sph_stat_distr
#' @export
d_sph_stat_Stereo <- function(x, p, a = 0, K_max = 1e3, method = "I", ...) {

  if (p <= 3) {

    stop("The asymptotic distribution is only available for p > 3.")

  }

  # psi_tilde_0 is the shift from U-stat to V-stat
  v_k2 <- weights_dfs_Sobolev(p = p, K_max = K_max, thre = 0, type = "Stereo",
                              Stereo_a = a, ...)$weights
  psi_tilde_0 <- Gegen_series(theta = 0, coefs = vk2_to_bk(vk2 = v_k2, p = p),
                              k = 1:K_max, p = p)
  cbind(d_Sobolev(x = x + psi_tilde_0, p = p, type = "Stereo", Stereo_a = a,
                  K_max = K_max, thre = 0, method = method, ...))

}
