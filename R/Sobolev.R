

#' @title Asymptotic distributions of Sobolev statistics of spherical uniformity
#'
#' @description Approximated density, distribution, and quantile functions, and
#' variance for the asymptotic null distributions of Sobolev statistics of
#' uniformity on \eqn{S^{p-1}:=\{{\bf x}\in R^p:||{\bf x}||=1\}}{S^{p-1}:=
#' \{x\in R^p:||x||=1\}}. These asymptotic distributions are infinite
#' weighted sums of (central) chi squared random variables:
#' \deqn{\sum_{k = 1}^\infty v_k^2 \chi^2_{d_{p, k}},}
#' where
#' \deqn{d_{p, k} := {{p + k - 3}\choose{p - 2}} + {{p + k - 2}\choose{p - 2}}}
#' is the dimension of the space of eigenfunctions of the Laplacian on
#' \eqn{S^{p-1}}, \eqn{p\ge 2}, associated to the \eqn{k}-th
#' eigenvalue, \eqn{k\ge 1}.
#'
#' @inheritParams r_unif
#' @param k sequence of integer indexes.
#' @param method method for approximating the density, distribution, or
#' quantile function of the weighted sum of chi squared random variables. Must
#' be \code{"I"} (Imhof), \code{"SW"} (Satterthwaite--Welch), \code{"HBE"}
#' (Hall--Buckley--Eagleson), or \code{"MC"} (Monte Carlo; only for distribution
#' or quantile functions). Defaults to \code{"I"}.
#' @param K_max integer giving the truncation of the series that compute the
#' asymptotic p-value of a Sobolev test. Defaults to \code{1e3}.
#' @param K_start either \code{0} or \code{1}, giving the order of the first
#' coefficient of the series that compute the asymptotic p-value of a Sobolev
#' test. If equals \code{0}, the bias coefficient is included.
#' Defaults to \code{1}.
#' @param thre error threshold for the tail probability given by the
#' the first terms of the truncated series of a Sobolev test. Defaults to
#' \code{1e-3}.
#' @param type name of the Sobolev statistic, using the naming from
#' \code{\link{avail_cir_tests}} and \code{\link{avail_sph_tests}}.
#' @param log compute the logarithm of \eqn{d_{p,k}}? Defaults to
#' \code{FALSE}.
#' @param lambda_grid a single parameter or vector of parameters for the
#' statistic. (See \code{Poisson_rho}, \code{Softmax_kappa}, or \code{Stereo_a}
#' for additional requirements).
#' @param verbose output information about the truncation? Defaults to
#' \code{TRUE}.
#' @inheritParams unif_stat
#' @inheritParams wschisq
#' @inheritParams Gegenbauer
#' @inheritParams Pn
#' @param force_positive set negative weights to zero? Defaults to \code{TRUE}.
#' @param ... further parameters passed to \code{*_\link{wschisq}}.
#' @return
#' \itemize{
#'   \item \code{d_p_k}: a vector of size \code{length(k)} with the
#'   evaluation of \eqn{d_{p,k}}.
#'   \item \code{weights_dfs_Sobolev}: a list with entries \code{weights} and
#'   \code{dfs}, automatically truncated according to \code{K_max} and
#'   \code{thre} (see details).
#'   \item \code{d_Sobolev}: density function evaluated at \code{x}, a vector.
#'   \item \code{p_Sobolev}: distribution function evaluated at \code{x},
#'   a vector.
#'   \item \code{q_Sobolev}: quantile function evaluated at \code{u}, a vector.
#'   \item \code{null_var_Sobolev}: variance of Sobolev statistic of sample size
#'   \code{n} and parameter(s) \code{lambda_grid}.
#' }
#' @author Eduardo García-Portugués and Paula Navarro-Esteban.
#' @details
#' The truncation of \eqn{\sum_{k = 1}^\infty v_k^2 \chi^2_{d_{p, k}}} is
#' done to the first \code{K_max} terms and then up to the index such that
#' the first terms explain the tail probability at the \code{x_tail} with
#' an absolute error smaller than \code{thre} (see details in
#' \code{\link{cutoff_wschisq}}). This automatic truncation takes place when
#' calling \code{*_Sobolev}. Setting \code{thre = 0} truncates to \code{K_max}
#' terms exactly. If the series only contains odd or even non-zero terms, then
#' only \code{K_max / 2} addends are \emph{effectively} taken into account
#' in the first truncation.
#' @examples
#' # Circular-specific statistics
#' curve(p_Sobolev(x = x, p = 2, type = "Watson", method = "HBE"),
#'       n = 2e2, ylab = "Distribution", main = "Watson")
#' curve(p_Sobolev(x = x, p = 2, type = "Rothman", method = "HBE"),
#'       n = 2e2, ylab = "Distribution", main = "Rothman")
#' curve(p_Sobolev(x = x, p = 2, type = "Pycke_q", method = "HBE"), to = 10,
#'       n = 2e2, ylab = "Distribution", main = "Pycke_q")
#' curve(p_Sobolev(x = x, p = 2, type = "Hermans_Rasson", method = "HBE"),
#'       to = 10, n = 2e2, ylab = "Distribution", main = "Hermans_Rasson")
#'
#' # Statistics for arbitrary dimensions
#' test_statistic <- function(type, to = 1, pmax = 5, M = 1e3, ...) {
#'
#'   col <- viridisLite::viridis(pmax - 1)
#'   curve(p_Sobolev(x = x, p = 2, type = type, method = "MC", M = M,
#'                   ...), to = to, n = 2e2, col = col[pmax - 1],
#'                   ylab = "Distribution", main = type, ylim = c(0, 1))
#'   for (p in 3:pmax) {
#'     curve(p_Sobolev(x = x, p = p, type = type, method = "MC", M = M,
#'                     ...), add = TRUE, n = 2e2, col = col[pmax - p + 1])
#'   }
#'   legend("bottomright", legend = paste("p =", 2:pmax), col = rev(col),
#'          lwd = 2)
#'
#' }
#'
#' # Ajne
#' test_statistic(type = "Ajne")
#' \donttest{
#' # Gine_Gn
#' test_statistic(type = "Gine_Gn", to = 1.5)
#'
#' # Gine_Fn
#' test_statistic(type = "Gine_Fn", to = 2)
#'
#' # Bakshaev
#' test_statistic(type = "Bakshaev", to = 3)
#'
#' # Riesz
#' test_statistic(type = "Riesz", Riesz_s = 0.5, to = 3)
#'
#' # PCvM
#' test_statistic(type = "PCvM", to = 0.6)
#'
#' # PAD
#' test_statistic(type = "PAD", to = 3)
#'
#' # PRt
#' test_statistic(type = "PRt", Rothman_t = 0.5)
#'
#' # Quantiles
#' p <- c(2, 3, 4, 11)
#' t(sapply(p, function(p) q_Sobolev(u = c(0.10, 0.05, 0.01), p = p,
#'                                   type = "PCvM")))
#' t(sapply(p, function(p) q_Sobolev(u = c(0.10, 0.05, 0.01), p = p,
#'                                   type = "PAD")))
#' t(sapply(p, function(p) q_Sobolev(u = c(0.10, 0.05, 0.01), p = p,
#'                                   type = "PRt")))
#'
#' # Series truncation for thre = 1e-5
#' sapply(p, function(p) length(weights_dfs_Sobolev(p = p, type = "PCvM")$dfs))
#' sapply(p, function(p) length(weights_dfs_Sobolev(p = p, type = "PRt")$dfs))
#' sapply(p, function(p) length(weights_dfs_Sobolev(p = p, type = "PAD")$dfs))
#'
#' # Variance
#' n <- 100
#' p <- 4
#' (var_Poisson <- null_var_Sobolev(n = n, p = p, type = "Poisson", K_max = 200,
#'                                  lambda_grid = seq(0.1, 0.9, 0.1)))
#' (var_Softmax <- null_var_Sobolev(n = n, p = p, type = "Softmax", K_max = 20,
#'                                  lambda_grid = c(0.1, 1, 5)))
#' (var_Stereo <- null_var_Sobolev(n = n, p = p, type = "Stereo",
#'                                 lambda_grid = seq(-1, 1, 0.25)))
#'
#' }
#' @name Sobolev


#' @rdname Sobolev
#' @export
d_p_k <- function(p, k, log = FALSE) {

  # log(nu_{p, k})
  p <- p - 2
  log_dfs <- lchoose(n = p + k - 1, k = p) + log(2 + p / k)
  log_dfs[k == 0] <- 0

  # nu_{p, k}
  if (!log) {

    log_dfs <- round(exp(log_dfs))

  }
  return(log_dfs)

}


#' @rdname Sobolev
#' @export
weights_dfs_Sobolev <- function(p, K_max = 1e3, thre = 1e-3, type,
                                Rothman_t = 1 / 3, Pycke_q = 0.5, Riesz_s = 1,
                                Poisson_rho = 0.5, Softmax_kappa = 1,
                                Stereo_a = 0, Sobolev_vk2 = c(0, 0, 1),
                                K_start = 1, log = FALSE, verbose = TRUE,
                                Gauss = TRUE, N = 320, tol = 1e-6,
                                force_positive = TRUE, x_tail = NULL) {

  # Check K_start lower than K_max
  if (K_start > 1) {

    stop(paste("K_start =", K_start, "must be either 0 or 1"))

  } else if (K_start >= K_max) {

    stop(paste("K_max =", K_max, "must be larger than K_start = ", K_start))

  }

  # alpha
  alpha <- 0.5 * p - 1

  # Sobolev weights and dfs
  if (p == 2 && type %in% c("Watson", "Rothman", "Hermans_Rasson", "Pycke_q")) {

    if (type == "Watson") {

      # Sequence of indexes
      k <- 1:K_max

      # log(v_k^2)
      log_vk2 <- -2 * log(k * pi)

      # log(d_{2, k})
      log_dk <- d_p_k(p = 2, k = k, log = TRUE)

      # Log weights and dfs
      log_weights <- log_vk2
      log_dfs <- log_dk

    } else if (type == "Rothman") {

      # Sequence of indexes
      k <- 1:K_max

      # log(v_k^2)
      log_vk2 <- log((sin(k * pi * Rothman_t) / (pi * k))^2)

      # log(d_{2, k})
      log_dk <- d_p_k(p = 2, k = k, log = TRUE)

      # Log weights and dfs
      log_weights <- log_vk2
      log_dfs <- log_dk

    } else if (type == "Hermans_Rasson") {

      # Halve K_max since we compute the odd/even coefficients separately
      K_max <- K_max %/% 2

      # Sequence of indexes
      k <- 1:K_max

      # log(v_{2 * k - 1}^2)
      log_v2km12 <- log(2 / pi) - log((2 * k - 1)^2)

      # log(v_{2 * k}^2)
      beta2 <- (pi^2 / 36) / (0.5 - 4 / pi^2)
      log_v2k2 <- log(2 * beta2 / pi) - log((2 * k)^2 - 1)

      # log(d_{p, 2 * k - 1})
      log_d2km1 <- d_p_k(p = 2, k = 2 * k - 1, log = TRUE)

      # log(d_{p, 2 * k})
      log_dk <- d_p_k(p = 2, k = 2 * k, log = TRUE)

      # Log weights and dfs
      log_weights <- c(rbind(log_v2km12, log_v2k2))
      log_dfs <- c(rbind(log_d2km1, log_dk))

    } else if (type == "Pycke_q") {

      # Sequence of indexes
      k <- 1:K_max

      # log(v_k^2)
      log_vk2 <- (k - 1) * log(Pycke_q)

      # log(d_{2, k})
      log_dk <- d_p_k(p = 2, k = k, log = TRUE)

      # Log weights and dfs
      log_weights <- log_vk2
      log_dfs <- log_dk

    }

  } else {

    if (type == "Ajne") {

      # Halve K_max since the even coefficients are zero
      K_max <- K_max %/% 2

      # Sequence of indexes
      k <- 1:K_max

      # log(v_{2 * k - 1}^2)
      log_v2km12 <- (p - 2) * log(2) + lgamma(alpha + 1) + lgamma(k + alpha) +
        lgamma(2 * k - 1) - (log(pi) + lgamma(k) + lgamma(2 * k + p - 2))
      log_v2km12 <- 2 * log_v2km12

      # Add the zero coefficients for the even terms
      log_v2km12 <- c(rbind(log_v2km12, -Inf))

      # log(d_{p, k})
      log_dk <- d_p_k(p = p, k = 1:(2 * K_max), log = TRUE)

      # Log weights and dfs
      log_weights <- log_v2km12
      log_dfs <- log_dk

    } else if (type == "Gine_Gn") {

      # Halve K_max since the odd coefficients are zero
      K_max <- K_max %/% 2

      # Sequence of indexes
      k <- 1:K_max

      # log(v_{2 * k}^2)
      log_v2k2 <- log((p - 1) * (2 * k - 1) / (8 * pi * (2 * k + p - 1))) +
        2 * (lgamma(alpha + 0.5) + lgamma(k - 0.5) - lgamma(k + alpha + 0.5))

      # Add the zero coefficients for the odd terms
      log_v2k2 <- c(rbind(-Inf, log_v2k2))

      # log(d_{p, k})
      log_dk <- d_p_k(p = p, k = 1:(2 * K_max), log = TRUE)

      # Log weights and dfs
      log_weights <- log_v2k2
      log_dfs <- log_dk

    } else if (type == "Gine_Fn") {

      # Halve K_max since we sum 4 * Ajne and Gine_Gn
      K_max <- K_max %/% 2

      # Sequence of indexes
      k <- 1:K_max

      # log(v_{2 * k - 1}^2)
      log_v2km12 <- (p - 2) * log(2) + lgamma(alpha + 1) + lgamma(k + alpha) +
        lgamma(2 * k - 1) - (log(pi) + lgamma(k) + lgamma(2 * k + p - 2))
      log_v2km12 <- 2 * log_v2km12

      # log(v_{2 * k}^2)
      log_v2k2 <- log((p - 1) * (2 * k - 1) / (8 * pi * (2 * k + p - 1))) +
        2 * (lgamma(alpha + 0.5) + lgamma(k - 0.5) - lgamma(k + alpha + 0.5))

      # log(d_{p, k})
      log_dk <- d_p_k(p = p, k = 1:(2 * K_max), log = TRUE)

      # Log weights and dfs
      log_weights <- c(rbind(log(4) + log_v2km12, log_v2k2))
      log_dfs <- log_dk

    } else if (type == "Bakshaev") {

      # Sequence of indexes
      k <- 1:K_max

      # log(b_k)
      log_vk2 <- switch((p > 2) + 1, log(8) - log(pi) - log(4 * k^2 - 1),
                        (p - 3) * log(2) + log(2 * k + p - 2) +
                          lgamma(k - 1 / 2) + lgamma(alpha) +
                          lgamma(p / 2) - log(pi) - lgamma(k + p - 1 / 2))

      # Switch from bk to vk2
      log_vk2 <- bk_to_vk2(bk = log_vk2, p = p, log = TRUE)

      # log(d_{p, k})
      log_dk <- d_p_k(p = p, k = k, log = TRUE)

      # Log weights and dfs
      log_weights <- log_vk2
      log_dfs <- log_dk

    } else if (type == "Riesz") {

      # Sequence of indexes
      k <- 1:K_max

      if (Riesz_s == 0) {

        if (p == 2) {

          log_vk2 <- -log(k)

        } else if (p == 3) {

          log_vk2 <- log(1 + 2 * k) - log(2 * k) - log(k + 1)

        } else {

          log_vk2 <- log(2 * k + p - 2) + lgamma(k) + lgamma(p - 2) -
            log(2) - lgamma(k + p - 1)

        }

      } else if (Riesz_s == 2) {

        if (p == 2) {

          log_vk2 <- ifelse(k == 1, log(2), -Inf)

        } else {

          log_vk2 <- ifelse(k == 1, log(2) - log(p - 2), -Inf)

        }

      } else {

        # tau_{s, k, p}
        if (p == 2) {

          tau <- 2^(Riesz_s + 1) / (1 + (k == 0))

        } else {

          tau <- 2^(p + Riesz_s - 3) * (p + 2 * k - 2) * gamma((p - 2) / 2)

        }

        # log(b_k)
        log_vk2 <- log(tau / sqrt(pi)) + lgamma((p - 1 + Riesz_s) / 2) +
          lgamma(-Riesz_s / 2 + k) - lgamma(p - 1 + Riesz_s / 2 + k) -
          lgamma(-Riesz_s / 2)

      }

      # Switch from bk to vk2
      log_vk2 <- bk_to_vk2(bk = log_vk2, p = p, log = TRUE)

      # log(d_{p, k})
      log_dk <- d_p_k(p = p, k = k, log = TRUE)

      # Log weights and dfs
      log_weights <- log_vk2
      log_dfs <- log_dk

    } else if (type == "PCvM") {

      # Sequence of indexes
      k <- 1:K_max

      # log(b_k)
      log_vk2 <- log(Gegen_coefs_Pn(k = k, p = p, type = "PCvM",
                                    Gauss = Gauss, N = N, tol = tol))

      # Switch from bk to vk2
      log_vk2 <- bk_to_vk2(bk = log_vk2, p = p, log = TRUE)

      # log(d_{p, k})
      log_dk <- d_p_k(p = p, k = k, log = TRUE)

      # Log weights and dfs
      log_weights <- log_vk2
      log_dfs <- log_dk

    } else if (type == "PAD") {

      # Sequence of indexes
      k <- 1:K_max

      # log(b_k)
      log_vk2 <- log(Gegen_coefs_Pn(k = k, p = p, type = "PAD",
                                    Gauss = Gauss, N = N, tol = tol))

      # Switch from bk to vk2
      log_vk2 <- bk_to_vk2(bk = log_vk2, p = p, log = TRUE)

      # log(d_{p, k})
      log_dk <- d_p_k(p = p, k = k, log = TRUE)

      # Log weights and dfs
      log_weights <- log_vk2
      log_dfs <- log_dk

    } else if (type == "PRt") {

      # Sequence of indexes
      k <- 1:K_max

      # log(b_k)
      log_vk2 <- log(Gegen_coefs_Pn(k = k, p = p, type = "PRt",
                                    Rothman_t = Rothman_t, Gauss = Gauss,
                                    N = N, tol = tol))

      # Switch from bk to vk2
      log_vk2 <- bk_to_vk2(bk = log_vk2, p = p, log = TRUE)

      # log(d_{p, k})
      log_dk <- d_p_k(p = p, k = k, log = TRUE)

      # Log weights and dfs
      log_weights <- log_vk2
      log_dfs <- log_dk

    } else if (type == "Poisson") {

      # Sequence of indexes
      k <- seq(K_start, K_max, 1)

      # log(b_k)
      if (p == 2) {

        log_vk2 <- log(2 - (k == 0)) + k * log(Poisson_rho)

      } else {

        log_vk2 <- log(2 * k + p - 2) - log(p - 2) + k * log(Poisson_rho)

      }
      log_vk2 <- log_vk2 + (p - 1) * log(1 - Poisson_rho) - log(1 + Poisson_rho)

      # Switch from bk to vk2
      log_vk2 <- bk_to_vk2(bk = log_vk2, p = p, log = TRUE, K_start = K_start)

      # log(d_{p, k})
      log_dk <- d_p_k(p = p, k = k, log = TRUE)

      # Log weights and dfs
      log_weights <- log_vk2
      log_dfs <- log_dk

    } else if (type == "Softmax") {

      # Sequence of indexes
      k <- seq(K_start, K_max, 1)

      # log(b_k)
      if (p == 2) {

        log_vk2 <- log(2 - (k == 0)) +
          log(besselI(x = Softmax_kappa, nu = k, expon.scaled = TRUE))

      } else {

        log_vk2 <- alpha * log(2 / Softmax_kappa) + lgamma(alpha) +
          log(k + alpha) + log(besselI(x = Softmax_kappa, nu = k + alpha,
                                       expon.scaled = TRUE))

      }

      # Switch from bk to vk2
      log_vk2 <- bk_to_vk2(bk = log_vk2, p = p, log = TRUE, K_start = K_start)

      # log(d_{p, k})
      log_dk <- d_p_k(p = p, k = k, log = TRUE)

      # Log weights and dfs
      log_weights <- log_vk2
      log_dfs <- log_dk

    } else if (type == "Stereo") {

      # Halve K_max since we compute the odd/even coefficients separately
      K_max <- K_max %/% 2

      # Sequence of indexes
      k_odd <- seq_len(K_max)
      k_even <- seq(K_start, K_max, 1)

      # log(v_{2 * k - 1}^2)
      log_v2km12 <- log(2 * k_odd - 1) + log(4 * (k_odd - 1) + p) -
        log(2 * k_odd + p - 3) + log(1 - Stereo_a)
      log_alpha_2km1 <- 2 * (lgamma(k_odd - 0.5) + lgamma(alpha) -
                               lgamma(k_odd + alpha - 0.5)) - log(2 * pi)

      # log(v_{2 * k}^2)
      log_v2k2 <- log(4 * k_even + p - 2) + log(1 + Stereo_a)
      log_alpha_2k <- 2 * (lgamma(k_even + 0.5) + lgamma(alpha) -
                             lgamma(k_even + 0.5 * (p - 1))) - log(2 * pi)

      # log(b_k)
      if (K_start == 0) {
        log_vk2 <- c(log_alpha_2k[1] + log_v2k2[1],
                     rbind(log_alpha_2km1 + log_v2km12,
                           log_alpha_2k[seq(2, length(log_alpha_2k), 1)] +
                             log_v2k2[seq(2, length(log_alpha_2k), 1)]))
      } else {
        log_vk2 <- c(rbind(log_alpha_2km1 + log_v2km12,
                           log_alpha_2k + log_v2k2))
      }

      # Switch from bk to vk2
      log_vk2 <- bk_to_vk2(bk = log_vk2, p = p, log = TRUE, K_start = K_start)

      # log(d_{p, k})
      log_dk <- d_p_k(p = p, k = seq(K_start, 2 * K_max, 1), log = TRUE)

      # Log weights and dfs
      log_weights <- log_vk2
      log_dfs <- log_dk

    } else if (type == "Sobolev") {

      # log(d_{p, k})
      log_dk <- d_p_k(p = p, k = seq_along(Sobolev_vk2), log = TRUE)

      # Log weights and dfs
      log_weights <- log(Sobolev_vk2)
      log_dfs <- log_dk

    } else {

      stop("Incompatible choice of p and type.")

    }

  }

  # Set NaNs (coming from negative weights) to 0?
  if (force_positive) {

    log_weights[is.nan(log_weights)] <- -Inf

  }

  # Automatic truncation based on the tail probability mean of the sum of
  # weighted chi squared
  K_max <- length(log_weights)
  cutoff <- cutoff_wschisq(thre = thre, weights = log_weights,
                           dfs = log_dfs, log = TRUE, x_tail = x_tail)$prob[1]

  # Truncate displaying optional information
  log_weights <- log_weights[1:cutoff]
  log_dfs <- log_dfs[1:cutoff]
  if (verbose) {

    message("Series truncated from ", K_max, " to ", cutoff,
            " terms (difference <= ", thre, " with the HBE tail probability;",
            " last weight = ", sprintf("%.3e", exp(log_weights[cutoff])), ").")

  }

  # Return weights and dfs?
  if (!log) {

    log_weights <- exp(log_weights)
    log_dfs <- exp(log_dfs)

  }

  return(list("weights" = log_weights, "dfs" = log_dfs))

}

#' @rdname Sobolev
#' @export
d_Sobolev <- function(x, p, type, method = c("I", "SW", "HBE")[1], K_max = 1e3,
                      thre = 1e-3, Rothman_t = 1 / 3, Pycke_q = 0.5,
                      Riesz_s = 1, Poisson_rho = 0.5, Softmax_kappa = 1,
                      Stereo_a = 0, Sobolev_vk2 = c(0, 0, 1), ncps = 0,
                      verbose = TRUE, N = 320, x_tail = NULL, ...) {

  weights_dfs <- weights_dfs_Sobolev(p = p, K_max = K_max, thre = thre,
                                     type = type, Rothman_t = Rothman_t,
                                     Pycke_q = Pycke_q, Riesz_s = Riesz_s,
                                     Poisson_rho = Poisson_rho,
                                     Softmax_kappa = Softmax_kappa,
                                     Stereo_a = Stereo_a,
                                     Sobolev_vk2 = Sobolev_vk2,
                                     verbose = verbose, Gauss = TRUE, N = N,
                                     x_tail = x_tail)
  d_wschisq(x = x, weights = weights_dfs$weights, dfs = weights_dfs$dfs,
            ncps = ncps, method = method, ...)

}


#' @rdname Sobolev
#' @export
p_Sobolev <- function(x, p, type, method = c("I", "SW", "HBE", "MC")[1],
                      K_max = 1e3, thre = 1e-3, Rothman_t = 1 / 3,
                      Pycke_q = 0.5, Riesz_s = 1, Poisson_rho = 0.5,
                      Softmax_kappa = 1, Stereo_a = 0, Sobolev_vk2 = c(0, 0, 1),
                      ncps = 0, verbose = TRUE, N = 320, x_tail = NULL, ...) {

  weights_dfs <- weights_dfs_Sobolev(p = p, K_max = K_max, thre = thre,
                                     type = type, Rothman_t = Rothman_t,
                                     Pycke_q = Pycke_q, Riesz_s = Riesz_s,
                                     Poisson_rho = Poisson_rho,
                                     Softmax_kappa = Softmax_kappa,
                                     Stereo_a = Stereo_a,
                                     Sobolev_vk2 = Sobolev_vk2,
                                     verbose = verbose, Gauss = TRUE, N = N,
                                     x_tail = x_tail)
  p_wschisq(x = x, weights = weights_dfs$weights, dfs = weights_dfs$dfs,
            ncps = ncps, method = method, ...)

}


#' @rdname Sobolev
#' @export
q_Sobolev <- function(u, p, type, method = c("I", "SW", "HBE", "MC")[1],
                      K_max = 1e3, thre = 1e-3, Rothman_t = 1 / 3,
                      Pycke_q = 0.5, Riesz_s = 1, Poisson_rho = 0.5,
                      Softmax_kappa = 1, Stereo_a = 0, Sobolev_vk2 = c(0, 0, 1),
                      ncps = 0, verbose = TRUE, N = 320, x_tail = NULL, ...) {

  weights_dfs <- weights_dfs_Sobolev(p = p, K_max = K_max, thre = thre,
                                     type = type, Rothman_t = Rothman_t,
                                     Pycke_q = Pycke_q, Riesz_s = Riesz_s,
                                     Poisson_rho = Poisson_rho,
                                     Softmax_kappa = Softmax_kappa,
                                     Stereo_a = Stereo_a,
                                     Sobolev_vk2 = Sobolev_vk2,
                                     verbose = verbose, Gauss = TRUE, N = N,
                                     x_tail = x_tail)
  q_wschisq(u = u, weights = weights_dfs$weights, dfs = weights_dfs$dfs,
            ncps = ncps, method = method, ...)

}

#' @rdname Sobolev
#' @export
null_var_Sobolev <- function(n, p, type, lambda_grid, K_max = 1e3,
                             verbose = TRUE) {

  # Check existence of Stereo depending on dimension
  if (type == "Stereo") {

    if (p == 2) {

      stop(paste0("'Stereo' statistic is not defined when p = 2."))

    } else if (p == 3) {

      stop(paste0("Variance of \'Stereo\' statistic under uniformity ",
                  "(H_0) is not finite when p = 3. If using it as input of ",
                  "unif_test_cv(), set \"rep(1, length(lambda_grid))\" as ",
                  "null_variance in order not to account it in the ",
                  "q-score function."))

    }
  }

  b0 <- numeric(length(lambda_grid))
  b0_sq <- numeric(length(lambda_grid))
  for (i in seq_along(lambda_grid)) {

    lambda <- lambda_grid[i]

    # Gegenbauer coefficients
    bk <- vk2_to_bk(weights_dfs_Sobolev(p = p, K_max = K_max, K_start = 0,
                                        type = type, Poisson_rho = lambda,
                                        Softmax_kappa = lambda,
                                        Stereo_a = lambda,
                                        verbose = verbose)$weights,
                    p = p, K_start = 0)

    # Normalizing constants (required for both cases)
    c_kp <- Gegen_coefs(k = seq_along(bk) - 1, p = p, only_const = TRUE)

    # Kernel bias coefficient
    b0[i] <- bk[1]

    # Bias of squared kernel coefficient
    b0_sq[i] <- Gegen_norm(coefs = bk, k = seq(0, K_max, 1), p = p,
                           c_kp = c_kp)^2 * (rotasym::w_p(p - 1)
                                             / rotasym::w_p(p))

  }

  return(2 * (n - 1) / n * (b0_sq - b0^2))

}


#' @title Finite Sobolev statistics for testing (hyper)spherical uniformity
#'
#' @description Computes the finite Sobolev statistic \deqn{
#' S_{n, p}(\{b_{k, p}\}_{k=1}^K) = \sum_{i, j = 1}^n
#' \sum_{k = 1}^K b_{k, p}C_k^(p / 2 - 1)(\cos^{-1}({\bf X}_i'{\bf X}_j)),}
#' for a sequence  \eqn{\{b_{k, p}\}_{k = 1}^K} of non-negative weights. For
#' \eqn{p = 2}, the Gegenbauer polynomials are replaced by Chebyshev ones.
#' @inheritParams sph_stat
#' @inheritParams cir_stat
#' @param vk2 weights for the finite Sobolev test. A non-negative vector or
#' matrix. Defaults to \code{c(0, 0, 1)}.
#' @return A matrix of size \code{c(M, ncol(vk2))} containing the statistics for
#' each of the \code{M} samples.
#' @export
sph_stat_Sobolev <- function(X, Psi_in_X = FALSE, p = 0, vk2 = c(0, 0, 1)) {

  # Compute Psi matrix with angles between pairs?
  if (Psi_in_X) {

    n <- n_from_dist_vector(nrow(X))
    if (p == 0) {

      stop("p >= 2 must be specified if Psi_in_X = TRUE.")

    }
    X <- X[, , 1, drop = FALSE]
    dim(X) <- dim(X)[1:2]
    M <- ncol(X)

  } else {

    n <- nrow(X)
    p <- ncol(X)
    M <- dim(X)[3]
    X <- Psi_mat(data = X)

  }

  # Statistics for k with vk2 != 0
  vk2 <- rbind(vk2)
  nonzero_vk2 <- which(apply(vk2 != 0, 2, any))
  Tnk <- matrix(0, nrow = M, ncol = ncol(vk2))
  for (j in seq_len(M)) {

      Tnk[j, nonzero_vk2] <- rowSums(Gegen_polyn(theta = X[, j],
                                                 k = nonzero_vk2, p = p))

  }

  # Get Gegenbauer coefficients
  bk <- vk2_to_bk(vk2 = vk2, p = p)

  # Construct statistic, a matrix of size c(M, length(bk))
  Tn <- (2 / n) * (Tnk %*% t(bk))

  # Add diagonal bias
  bias <- drop(bk[, nonzero_vk2] %*% Gegen_polyn(theta = 0, k = nonzero_vk2,
                                                 p = p))
  Tn <- t(t(Tn) + bias)
  return(unname(Tn))

}


#' @rdname sph_stat_Sobolev
#' @export
cir_stat_Sobolev <- function(Theta, Psi_in_Theta = FALSE, vk2 = c(0, 0, 1)) {

  if (Psi_in_Theta) {

    if (length(dim(Theta)) < 3) {

      dim(Theta) <- c(dim(Theta), 1)

    }
    return(sph_stat_Sobolev(X = Theta, Psi_in_X = TRUE, p = 2, vk2 = vk2))

  } else {

    return(sph_stat_Sobolev(X = Theta_to_X(Theta), Psi_in_X = FALSE, vk2 = vk2))

  }

}


#' @title Transformation between different coefficients in Sobolev statistics
#'
#' @description Given a Sobolev statistic
#' \deqn{S_{n, p} = \sum_{i, j = 1}^n \psi(\cos^{-1}({\bf X}_i'{\bf X}_j)),}{
#' S_{n, p} = \sum_{i, j = 1}^n \psi(\cos^{-1}(X_i'X_j)),}
#' for a sample \eqn{{\bf X}_1, \ldots, {\bf X}_n \in S^{p - 1} := \{{\bf x}
#' \in R^p : ||{\bf x}|| = 1\}}{X_1, \ldots, X_n \in S^{p - 1} :=
#' \{x \in R^p : ||x|| = 1\}}, \eqn{p\ge 2}, three important sequences
#' are related to \eqn{S_{n, p}}.
#' \itemize{
#' \item \link[=Gegen_coefs]{Gegenbauer coefficients} \eqn{\{b_{k, p}\}} of
#' \eqn{\psi_p} (see, e.g., the \link[=Pn]{projected-ecdf statistics}), given
#' by
#' \deqn{b_{k, p} := \frac{1}{c_{k, p}}\int_0^\pi \psi_p(\theta)
#' C_k^{p / 2 - 1}(\cos\theta)\,\mathrm{d}\theta.}{
#' b_{k, p} := \frac{1}{c_{k, p}} \int_0^\pi \psi_p(\theta)
#' C_k^(p / 2 - 1)(\cos\theta) d\theta.}
#' \item Weights \eqn{\{v_{k, p}^2\}} of the
#' \link[=Sobolev]{asymptotic distribution} of the Sobolev statistic,
#' \eqn{\sum_{k = 1}^\infty v_k^2 \chi^2_{d_{p, k}}}, given by
#' \deqn{v_{k, p}^2 = \left(1 + \frac{2k}{p - 2}\right)^{-1} b_{k, p},
#' \quad p \ge 3.}{v_{k, p}^2 = (1 + 2k / (p - 2))^{-1} b_{k, p}, p \ge 3.}
#' \item Gegenbauer coefficients \eqn{\{u_{k, p}\}} of the
#' \link[=locdev]{local projected alternative} associated to \eqn{S_{n, p}},
#' given by
#' \deqn{u_{k, p} = \left(1 + \frac{2k}{p - 2}\right) v_{k, p},
#' \quad p \ge 3.}{u_{k, p} = (1 + 2k / (p - 2)) b_{k, p}, p \ge 3.}
#' }
#' For \eqn{p = 2}, the factor \eqn{(1 + 2k / (p - 2))} is replaced by \eqn{2}.
#'
#' @param bk coefficients \eqn{b_{k, p}} associated to the indexes
#' \code{1:length(bk)} (if \code{K_start = 1}) or \code{0:(length(bk)-1)} (if
#' \code{K_start = 0}), a vector.
#' @param vk2 \bold{squared} coefficients \eqn{v_{k, p}^2} associated to the
#' indexes  \code{1:length(vk2)} (if \code{K_start = 1}) or
#' \code{0:(length(vk2)-1)} (if \code{K_start = 0}), a vector.
#' @param uk coefficients \eqn{u_{k, p}} associated to the indexes
#' \code{1:length(uk)} (if \code{K_start = 1}) or \code{0:(length(uk)-1)}
#' (if \code{K_start = 0}), a vector.
#' @inheritParams r_unif_sph
#' @param signs signs of the coefficients \eqn{u_{k, p}}, a vector of the
#' same size as \code{vk2} or \code{bk}, or a scalar. Defaults to \code{1}.
#' @param log do operations in log scale (log-in, log-out)? Defaults to
#' \code{FALSE}.
#' @param K_start either \code{0} or \code{1}, giving the order \eqn{k} of the
#'  first coefficient in \eqn{b_{k,p}}, \eqn{v^2_{k,p}}, and \eqn{u_{k,p}}.
#' Defaults to \code{1}.
#' @return The corresponding vectors of coefficients \code{vk2}, \code{bk}, or
#' \code{uk}, depending on the call.
#' @details
#' See more details in Prentice (1978) and García-Portugués et al. (2023). The
#' adequate signs of \code{uk} for the \code{"PRt"} \link[=Pn]{Rothman test}
#' can be retrieved with \code{\link{akx}} and \code{sqr = TRUE}, see the
#' examples.
#' @references
#' García-Portugués, E., Navarro-Esteban, P., Cuesta-Albertos, J. A. (2023)
#' On a projection-based class of uniformity tests on the hypersphere.
#' \emph{Bernoulli}, 29(1):181--204. \doi{10.3150/21-BEJ1454}.
#'
#' Prentice, M. J. (1978). On invariant tests of uniformity for directions and
#' orientations. \emph{The Annals of Statistics}, 6(1):169--176.
#' \doi{10.1214/aos/1176344075}
#' @examples
#' # bk, vk2, and uk for the PCvM test in p = 3
#' (bk <- Gegen_coefs_Pn(k = 1:5, type = "PCvM", p = 3))
#' (vk2 <- bk_to_vk2(bk = bk, p = 3))
#' (uk <- bk_to_uk(bk = bk, p = 3))
#'
#' # vk2 is the same as
#' weights_dfs_Sobolev(K_max = 10, thre = 0, p = 3, type = "PCvM")$weights
#'
#' # bk and uk for the Rothman test in p = 3, with adequate signs
#' t <- 1 / 3
#' (bk <- Gegen_coefs_Pn(k = 1:5, type = "PRt", p = 3, Rothman_t = t))
#' (ak <- akx(x = drop(q_proj_unif(t, p = 3)), p = 3, k = 1:5, sqr = TRUE))
#' (uk <- bk_to_uk(bk = bk, p = 3, signs = ak))
#' @name Sobolev_coefs


#' @rdname Sobolev_coefs
#' @export
bk_to_vk2 <- function(bk, p, log = FALSE, K_start = 1) {

  # Check dimension
  p <- as.integer(p)
  stopifnot(p >= 2)

  # Check K_start
  stopifnot(K_start %in% c(0, 1))

  # Compute k
  if (is.matrix(bk)) {

    if (K_start == 1) {

      k <- matrix(seq_len(ncol(bk)), nrow = nrow(bk), ncol = ncol(bk),
                byrow = TRUE)

    } else {

      k <- matrix(seq_len(ncol(bk)) - 1, nrow = nrow(bk), ncol = ncol(bk),
                  byrow = TRUE)

    }

  } else {

    if (K_start == 1) {

      k <- seq_along(bk)

    } else {

      k <- seq_along(bk) - 1

    }

  }

  # Add factor
  if (log) {

    if (p == 2) {

      return(bk - log(2))

    } else {

      return(bk - log(1 + 2 * k / (p - 2)))

    }

  } else {

    if (p == 2) {

      return(bk / 2)

    } else {

      return(bk / (1 + 2 * k / (p - 2)))

    }

  }

}


#' @rdname Sobolev_coefs
#' @export
bk_to_uk <- function(bk, p, signs = 1, K_start = 1) {

  # Check dimension
  p <- as.integer(p)
  stopifnot(p >= 2)

  # Check signs
  stopifnot(length(signs) %in% c(1, length(bk)))

  # Check K_start
  stopifnot(K_start %in% c(0, 1))

  # Compute k
  if (is.matrix(bk)) {

    if (K_start == 1) {

      k <- matrix(seq_len(ncol(bk)), nrow = nrow(bk), ncol = ncol(bk),
                  byrow = TRUE)

    } else {

      k <- matrix(seq_len(ncol(bk)) - 1, nrow = nrow(bk), ncol = ncol(bk),
                  byrow = TRUE)

    }

  } else {

    if (K_start == 1) {

      k <- seq_along(bk)

    } else {

      k <- seq_along(bk) - 1

    }

  }

  # Add factor
  if (p == 2) {

    return(sign(signs) * sqrt(2 * bk))

  } else {

    return(sign(signs) * sqrt((1 + 2 * k / (p - 2)) * bk))

  }

}


#' @rdname Sobolev_coefs
#' @export
vk2_to_bk <- function(vk2, p, log = FALSE, K_start = 1) {

  # Check dimension
  p <- as.integer(p)
  stopifnot(p >= 2)

  # Check K_start
  stopifnot(K_start %in% c(0, 1))

  # Compute k
  if (is.matrix(vk2)) {

    if (K_start == 1) {

      k <- matrix(seq_len(ncol(vk2)), nrow = nrow(vk2), ncol = ncol(vk2),
                  byrow = TRUE)

    } else {

      k <- matrix(seq_len(ncol(vk2)) - 1, nrow = nrow(vk2), ncol = ncol(vk2),
                  byrow = TRUE)

    }

  } else {

    if (K_start == 1) {

      k <- seq_along(vk2)

    } else {

      k <- seq_along(vk2) - 1

    }

  }

  # Add factor
  if (log) {

    if (p == 2) {

      return(vk2 + log(2))

    } else {

      return(vk2 + log(1 + 2 * k / (p - 2)))

    }

  } else {

    if (p == 2) {

      return(2 * vk2)

    } else {

      return((1 + 2 * k / (p - 2)) * vk2)

    }

  }

}


#' @rdname Sobolev_coefs
#' @export
vk2_to_uk <- function(vk2, p, signs = 1, K_start = 1) {

  # Check dimension
  p <- as.integer(p)
  stopifnot(p >= 2)

  # Check signs
  stopifnot(length(signs) %in% c(1, length(vk2)))

  # Check K_start
  stopifnot(K_start %in% c(0, 1))

  # Compute k
  if (is.matrix(vk2)) {

    if (K_start == 1) {

      k <- matrix(seq_len(ncol(vk2)), nrow = nrow(vk2), ncol = ncol(vk2),
                  byrow = TRUE)

    } else {

      k <- matrix(seq_len(ncol(vk2)) - 1, nrow = nrow(vk2), ncol = ncol(vk2),
                  byrow = TRUE)

    }

  } else {

    if (K_start == 1) {

      k <- seq_along(vk2)

    } else {

      k <- seq_along(vk2) - 1

    }

  }

  # Add factor
  if (p == 2) {

    return(2 * sign(signs) * sqrt(vk2))

  } else {

    return((1 + 2 * k / (p - 2)) * sign(signs) * sqrt(vk2))

  }

}


#' @rdname Sobolev_coefs
#' @export
uk_to_vk2 <- function(uk, p, K_start = 1) {

  # Check dimension
  p <- as.integer(p)
  stopifnot(p >= 2)

  # Check K_start
  stopifnot(K_start %in% c(0, 1))

  # Compute k
  if (is.matrix(uk)) {

    if (K_start == 1) {

      k <- matrix(seq_len(ncol(uk)), nrow = nrow(uk), ncol = ncol(uk),
                  byrow = TRUE)

    } else {

      k <- matrix(seq_len(ncol(uk)) - 1, nrow = nrow(uk), ncol = ncol(uk),
                  byrow = TRUE)

    }

  } else {

    if (K_start == 1) {

      k <- seq_along(uk)

    } else {

      k <- seq_along(uk) - 1

    }

  }

  # Add factor
  if (p == 2) {

    return((uk / 2)^2)

  } else {

    return((uk / (1 + 2 * k / (p - 2)))^2)

  }

}


#' @rdname Sobolev_coefs
#' @export
uk_to_bk <- function(uk, p, K_start = 1) {

  # Check dimension
  p <- as.integer(p)
  stopifnot(p >= 2)

  # Check K_start
  stopifnot(K_start %in% c(0, 1))

  # Compute k
  if (is.matrix(uk)) {

    if (K_start == 1) {

      k <- matrix(seq_len(ncol(uk)), nrow = nrow(uk), ncol = ncol(uk),
                  byrow = TRUE)

    } else {

      k <- matrix(seq_len(ncol(uk)) - 1, nrow = nrow(uk), ncol = ncol(uk),
                  byrow = TRUE)

    }

  } else {

    if (K_start == 1) {

      k <- seq_along(uk)

    } else {

      k <- seq_along(uk) - 1

    }

  }

  # Add factor
  if (p == 2) {

    return(uk^2 / 2)

  } else {

    return(uk^2 / (1 + 2 * k / (p - 2)))

  }

}
