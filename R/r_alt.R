

#' @title Sample non-uniformly distributed spherical data
#'
#' @description Simple simulation of prespecified non-uniform spherical
#' distributions: von Mises--Fisher (vMF), Mixture of vMF (MvMF),
#' Angular Central Gaussian (ACG), Small Circle (SC), Watson (W),
#' Cauchy-like (C), Mixture of Cauchy-like (MC), or Uniform distribution with
#' Antipodal-Dependent observations (UAD).
#'
#' @inheritParams r_unif
#' @param alt alternative, must be \code{"vMF"}, \code{"MvMF"}, \code{"ACG"},
#' \code{"SC"}, \code{"W"}, \code{"C"}, \code{"MC"}, or \code{"UAD"}. See
#' details below.
#' @param mu location parameter for \code{"vMF"}, \code{"SC"}, \code{"W"}, and
#' \code{"C"}. Defaults to \code{c(rep(0, p - 1), 1)}.
#' @param kappa non-negative parameter measuring the strength of the deviation
#' with respect to uniformity (obtained with \eqn{\kappa = 0}).
#' @param nu projection along \eqn{{\bf e}_p}{e_p} controlling the modal
#' strip of the small circle distribution. Must be in (-1, 1). Defaults to
#' \code{0.5}.
#' @inheritParams unif_cap
#' @param F_inv quantile function returned by \code{\link{F_inv_from_f}}. Used
#' for \code{"SC"}, \code{"W"}, and \code{"C"}. Computed by internally if
#' \code{NULL} (default).
#' @inheritParams F_inv_from_f
#' @param axial_mix use a mixture of von Mises--Fisher or Cauchy-like that is
#' axial (i.e., symmetrically distributed about the origin)? Defaults to
#' \code{TRUE}.
#' @details
#' The parameter \code{kappa} is used as \eqn{\kappa} in the following
#' distributions:
#' \itemize{
#'   \item \code{"vMF"}: von Mises--Fisher distribution with concentration
#'   \eqn{\kappa} and directional mean \eqn{\boldsymbol{\mu}}.
#'   \item \code{"MvMF"}: equally-weighted mixture of \eqn{p} von Mises--Fisher
#'   distributions with common concentration \eqn{\kappa} and directional means
#'   \eqn{\pm{\bf e}_1, \ldots, \pm{\bf e}_p}{±e_1, \ldots, ±e_p} if
#'   \code{axial_mix = TRUE}. If \code{axial_mix = FALSE}, then only means
#'   with positive signs are considered.
#'   \item \code{"ACG"}: Angular Central Gaussian distribution with diagonal
#'   shape matrix with diagonal given by
#'   \deqn{(1, \ldots, 1, 1 + \kappa) / (p + \kappa).}
#'   \item \code{"SC"}: Small Circle distribution with axis mean
#'   \eqn{\boldsymbol{\mu}} and concentration \eqn{\kappa} about the projection
#'   along the mean, \eqn{\nu}.
#'   \item \code{"W"}: Watson distribution with axis mean \eqn{\boldsymbol{\mu}}
#'   and concentration \eqn{\kappa}. The Watson distribution is a particular
#'   case of the Bingham distribution.
#'   \item \code{"C"}: Cauchy-like distribution with directional mode
#'   \eqn{\boldsymbol{\mu}} and concentration
#'   \eqn{\kappa = \rho / (1 - \rho^2)}. The circular Wrapped Cauchy
#'   distribution is a particular case of this Cauchy-like distribution.
#'   \item \code{"MC"}: equally-weighted mixture of \eqn{p} Cauchy-like
#'   distributions with common concentration \eqn{\kappa} and directional means
#'   \eqn{\pm{\bf e}_1, \ldots, \pm{\bf e}_p}{±e_1, \ldots, ±e_p} if
#'   \code{axial_mix = TRUE}. If \code{axial_mix = FALSE}, then only means
#'   with positive signs are considered.
#' }
#' The alternative \code{"UAD"} generates a sample formed by
#' \eqn{\lceil n/2\rceil} observations drawn uniformly on \eqn{S^{p-1}}
#' and the remaining observations drawn from a uniform spherical cap
#' distribution of angle \eqn{\pi-\kappa} about each of the
#' \eqn{\lceil n/2\rceil} observations (see \code{\link{unif_cap}}). Hence,
#' \code{kappa = 0} corresponds to a spherical cap covering the whole sphere and
#' \code{kappa = pi} is a one-point degenerate spherical cap.
#' @return An \bold{array} of size \code{c(n, p, M)} with \code{M} random
#' samples of size \code{n} of non-uniformly-generated directions on
#' \eqn{S^{p-1}}.
#' @details
#' Much faster sampling for \code{"SC"}, \code{"W"}, \code{"C"}, and \code{"MC"}
#' is achieved providing \code{F_inv}; see examples.
#' @examples
#' ## Simulation with p = 2
#'
#' p <- 2
#' n <- 50
#' kappa <- 20
#' nu <- 0.5
#' angle <- pi / 10
#' rho <- ((2 * kappa + 1) - sqrt(4 * kappa + 1)) / (2 * kappa)
#' F_inv_SC_2 <- F_inv_from_f(f = function(z) exp(-kappa * (z - nu)^2), p = 2)
#' F_inv_W_2 <- F_inv_from_f(f = function(z) exp(kappa * z^2), p = 2)
#' F_inv_C_2 <- F_inv_from_f(f = function(z) (1 - rho^2) /
#'                             (1 + rho^2 - 2 * rho * z)^(p / 2), p = 2)
#' x1 <- r_alt(n = n, p = p, alt = "vMF", kappa = kappa)[, , 1]
#' x2 <- r_alt(n = n, p = p, alt = "C", F_inv = F_inv_C_2)[, , 1]
#' x3 <- r_alt(n = n, p = p, alt = "SC", F_inv = F_inv_SC_2)[, , 1]
#' x4 <- r_alt(n = n, p = p, alt = "ACG", kappa = kappa)[, , 1]
#' x5 <- r_alt(n = n, p = p, alt = "W", F_inv = F_inv_W_2)[, , 1]
#' x6 <- r_alt(n = n, p = p, alt = "MvMF", kappa = kappa)[, , 1]
#' x7 <- r_alt(n = n, p = p, alt = "MC", kappa = kappa)[, , 1]
#' x8 <- r_alt(n = n, p = p, alt = "UAD", kappa = 1 - angle)[, , 1]
#' r <- runif(n, 0.95, 1.05) # Radius perturbation to improve visualization
#' plot(r * x1, pch = 16, xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1), col = 1)
#' points(r * x2, pch = 16, col = 2)
#' points(r * x3, pch = 16, col = 3)
#' plot(r * x4, pch = 16, xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1), col = 1)
#' points(r * x5, pch = 16, col = 2)
#' points(r * x6, pch = 16, col = 3)
#' points(r * x7, pch = 16, col = 4)
#' col <- rep(rainbow(n / 2), 2)
#' plot(r * x8, pch = 16, xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1), col = col)
#' for (i in seq(1, n, by = 2)) lines((r * x8)[i + 0:1, ], col = col[i])
#'
#' ## Simulation with p = 3
#'
#' n <- 50
#' p <- 3
#' kappa <- 20
#' angle <- pi / 10
#' nu <- 0.5
#' rho <- ((2 * kappa + 1) - sqrt(4 * kappa + 1)) / (2 * kappa)
#' F_inv_SC_3 <- F_inv_from_f(f = function(z) exp(-kappa * (z - nu)^2), p = 3)
#' F_inv_W_3 <- F_inv_from_f(f = function(z) exp(kappa * z^2), p = 3)
#' F_inv_C_3 <- F_inv_from_f(f = function(z) (1 - rho^2) /
#'                             (1 + rho^2 - 2 * rho * z)^(p / 2), p = 3)
#' x1 <- r_alt(n = n, p = p, alt = "vMF", kappa = kappa)[, , 1]
#' x2 <- r_alt(n = n, p = p, alt = "C", F_inv = F_inv_C_3)[, , 1]
#' x3 <- r_alt(n = n, p = p, alt = "SC", F_inv = F_inv_SC_3)[, , 1]
#' x4 <- r_alt(n = n, p = p, alt = "ACG", kappa = kappa)[, , 1]
#' x5 <- r_alt(n = n, p = p, alt = "W", F_inv = F_inv_W_3)[, , 1]
#' x6 <- r_alt(n = n, p = p, alt = "MvMF", kappa = kappa)[, , 1]
#' x7 <- r_alt(n = n, p = p, alt = "MC", kappa = kappa)[, , 1]
#' x8 <- r_alt(n = n, p = p, alt = "UAD", kappa = 1 - angle)[, , 1]
#' s3d <- scatterplot3d::scatterplot3d(x1, pch = 16, xlim = c(-1.1, 1.1),
#'                                     ylim = c(-1.1, 1.1), zlim = c(-1.1, 1.1))
#' s3d$points3d(x2, pch = 16, col = 2)
#' s3d$points3d(x3, pch = 16, col = 3)
#' s3d <- scatterplot3d::scatterplot3d(x4, pch = 16, xlim = c(-1.1, 1.1),
#'                                     ylim = c(-1.1, 1.1), zlim = c(-1.1, 1.1))
#' s3d$points3d(x5, pch = 16, col = 2)
#' s3d$points3d(x6, pch = 16, col = 3)
#' s3d$points3d(x7, pch = 16, col = 4)
#' col <- rep(rainbow(n / 2), 2)
#' s3d <- scatterplot3d::scatterplot3d(x8, pch = 16, xlim = c(-1.1, 1.1),
#'                                     ylim = c(-1.1, 1.1), zlim = c(-1.1, 1.1),
#'                                     color = col)
#' for (i in seq(1, n, by = 2)) s3d$points3d(x8[i + 0:1, ], col = col[i],
#'                                           type = "l")
#' @export
r_alt <- function(n, p, M = 1, alt = "vMF", mu = c(rep(0, p - 1), 1),
                  kappa = 1, nu = 0.5, F_inv = NULL, K = 1e3,
                  axial_mix = TRUE) {

  # Check location and concentration parameters, and sample size
  mu <- rotasym::check_unit_norm(x = mu, warnings = TRUE)
  stopifnot(kappa >= 0)
  stopifnot(n >= 1)

  # Sampling from uniform
  if (kappa == 0) {

    return(r_unif_sph(n = n, p = p, M = M))

  }

  # Choose alternative
  if (alt == "vMF") {

    long_samp <- rotasym::r_vMF(n = n * M, mu = mu, kappa = kappa)

  } else if (alt == "MvMF") {

    # Mixture components
    j <- sample(x = 1:p, size = n * M, replace = TRUE)
    nM_j <- tabulate(bin = j, nbins = p)
    mu_j <- diag(1, nrow = p, ncol = p)

    # Sample components
    long_samp <- matrix(nrow = n * M, ncol = p)
    for (k in which(nM_j > 0)) {

      long_samp[j == k, ] <- rotasym::r_vMF(n = nM_j[k], mu = mu_j[k, ],
                                            kappa = kappa)

    }

    # Add plus and minus means
    if (axial_mix) {

      long_samp <- sample(x = c(-1, 1), size = n * M, replace = TRUE) *
        long_samp

    }

    # Shuffle data
    long_samp <- long_samp[sample(x = n * M), , drop = FALSE]

  } else if (alt == "ACG") {

    Lambda <- diag(c(rep(1 / (p + kappa), p - 1), (1 + kappa) / (p + kappa)),
                   nrow = p, ncol = p)
    long_samp <- rotasym::r_ACG(n = n * M, Lambda = Lambda)

  } else if (alt == "SC") {

    # Compute the inverse of the distribution function F?
    if (is.null(F_inv)) {

      stopifnot(-1 < nu & nu < 1)
      f <- function(z) exp(-kappa * (z - nu)^2)
      F_inv <- F_inv_from_f(f = f, p = p, K = K)

    }

    # Sample the small circle distribution
    r_U <- function(n) r_unif_sph(n = n, p = p - 1, M = 1)[, , 1]
    r_V <- function(n) F_inv(runif(n = n))
    long_samp <- rotasym::r_tang_norm(n = n * M, theta = mu,
                                      r_U = r_U, r_V = r_V)

  } else if (alt == "W") {

    # Compute the inverse of the distribution function F?
    if (is.null(F_inv)) {

      f <- function(z) exp(kappa * (z^2 - 1))
      F_inv <- F_inv_from_f(f = f, p = p, K = K)

    }

    # Sample the small circle distribution
    r_U <- function(n) r_unif_sph(n = n, p = p - 1, M = 1)[, , 1]
    r_V <- function(n) F_inv(runif(n = n))
    long_samp <- rotasym::r_tang_norm(n = n * M, theta = mu,
                                      r_U = r_U, r_V = r_V)

  } else if (alt == "C") {

    # Compute the inverse of the distribution function F?
    if (is.null(F_inv)) {

      rho <- ifelse(kappa == 0, 0,
                    ((2 * kappa + 1) - sqrt(4 * kappa + 1)) / (2 * kappa))
      f <- function(z) ((1 - rho) / sqrt(1 + rho^2 - 2 * rho * z))^p
      F_inv <- F_inv_from_f(f = f, p = p, K = K)

    }

    # Sample the small circle distribution
    r_U <- function(n) r_unif_sph(n = n, p = p - 1, M = 1)[, , 1]
    r_V <- function(n) F_inv(runif(n = n))
    long_samp <- rotasym::r_tang_norm(n = n * M, theta = mu,
                                      r_U = r_U, r_V = r_V)

  } else if (alt == "MC") {

    # Mixture components
    j <- sample(x = 1:p, size = n * M, replace = TRUE)
    nM_j <- tabulate(bin = j, nbins = p)
    mu_j <- diag(1, nrow = p, ncol = p)

    # Compute the inverse of the distribution function F?
    if (is.null(F_inv)) {

      rho <- ifelse(kappa == 0, 0,
                    ((2 * kappa + 1) - sqrt(4 * kappa + 1)) / (2 * kappa))
      f <- function(z) ((1 - rho) / sqrt(1 + rho^2 - 2 * rho * z))^p
      F_inv <- F_inv_from_f(f = f, p = p, K = K)

    }

    # Tangent-normal component distributions
    r_U <- function(n) r_unif_sph(n = n, p = p - 1, M = 1)[, , 1]
    r_V <- function(n) F_inv(runif(n = n))

    # Sample components
    long_samp <- matrix(nrow = n * M, ncol = p)
    for (k in which(nM_j > 0)) {

      long_samp[j == k, ] <- rotasym::r_tang_norm(n = nM_j[k],
                                                  theta = mu_j[k, ],
                                                  r_U = r_U, r_V = r_V)

    }

    # Add plus and minus means
    if (axial_mix) {

      long_samp <- sample(x = c(-1, 1), size = n * M, replace = TRUE) *
        long_samp

    }

    # Shuffle data
    long_samp <- long_samp[sample(x = n * M), , drop = FALSE]

  } else if (alt == "UAD") {

    # Check kappa
    if (kappa < 0) {

      stop("kappa must be non-negative.")

    }
    if (kappa > pi) {

      stop("kappa must be lower than pi.")

    }

    # Sample uniform
    n_2 <- ceiling(n / 2)
    n_ant <- n - n_2
    u <- rbind(r_unif_sph(n = n_2 * M, p = p, M = 1)[, , 1])

    # Sample antipodal
    ant <- array(dim = c(n_ant * M, p))
    for (i in seq_len(n_ant * M)) {

      ant[i, ] <- r_unif_cap(n = 1, mu = -u[i, ], angle = pi - kappa)

    }

    # Intercalate rows
    long_samp <- rbind(u, ant)
    if (n > 1) {

      ind_u <- seq_len(n_2 * M)
      ind_ant <- seq(n_2 * M + 1, n * M, by = 1)
      if (n %% 2 != 0) ind_ant <- c(ind_ant, rep(NA, M))
      long_samp <- na.omit(long_samp[c(rbind(ind_u, ind_ant)), ])

    }

  } else {

    stop(paste("Wrong alt; must be \"vMF\", \"MvMF\", \"Bing\"",
               "\"ACG\", \"SC\", \"W\", \"C\", \"MC\", or \"UAD\"."))

  }

  # As an array
  samp <- array(dim = c(n, p, M))
  for (j in 1:M) {

    samp[, , j] <- long_samp[(1 + (j - 1) * n):(j * n), , drop = FALSE]

  }
  return(samp)

}
