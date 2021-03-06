

#' @title Asymptotic distributions of Sobolev statistics of spherical uniformity
#'
#' @description Approximated density, distribution, and quantile functions for
#' the asymptotic null distributions of Sobolev statistics of uniformity
#' on \eqn{S^{p-1}:=\{{\bf x}\in R^p:||{\bf x}||=1\}}{S^{p-1}:=
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
#' @param K_max integer giving the truncation of the series that compute the
#' asymptotic p-value of a Sobolev test. Defaults to \code{1e3}.
#' @param thre error threshold for the tail probability given by the
#' the first terms of the truncated series of a Sobolev test. Defaults to
#' \code{1e-3}.
#' @param type Sobolev statistic. For \eqn{p = 2}, either \code{"Rothman"},
#' \code{"Pycke_q"}, or \code{"Hermans_Rasson"}. For \eqn{p \ge 2},
#' \code{"Ajne"}, \code{"Gine_Gn"}, \code{"Gine_Fn"}, \code{"Bakshaev"},
#' \code{"PCvM"}, \code{"PAD"} or \code{"PRt"}.
#' @param log compute the logarithm of \eqn{d_{p,k}}? Defaults to
#' \code{FALSE}.
#' @param verbose output information about the truncation? Defaults to
#' \code{TRUE}.
#' @inheritParams unif_stat
#' @inheritParams wschisq
#' @inheritParams Gegenbauer
#' @inheritParams Pn
#' @param force_positive set negative
#' @param ... further parameters passed to \code{*_\link{wschisq}}.
#' @return
#' \itemize{
#'   \item \code{d_p_k}: a vector of length \code{length(k)} with the
#'   evaluation of \eqn{d_{p,k}}.
#'   \item \code{weights_dfs_Sobolev}: a list with entries \code{weights} and
#'   \code{dfs}, automatically truncated according to \code{K_max} and
#'   \code{thre} (see details).
#'   \item \code{d_Sobolev}: density function evaluated at \code{x}, a vector.
#'   \item \code{p_Sobolev}: distribution function evaluated at \code{x},
#'   a vector.
#'   \item \code{q_Sobolev}: quantile function evaluated at \code{u}, a vector.
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
#' curve(p_Sobolev(x = x, p = 2, type = "Rothman"),
#'       n = 2e2, ylab = "Distribution", main = "Rothman")
#' curve(p_Sobolev(x = x, p = 2, type = "Pycke_q"), to = 10,
#'       n = 2e2, ylab = "Distribution", main = "Pycke_q")
#' curve(p_Sobolev(x = x, p = 2, type = "Hermans_Rasson"),
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
#' }
#' @name Sobolev


#' @rdname Sobolev
#' @export
d_p_k <- function(p, k, log = FALSE) {

  # log(nu_{p, k})
  p <- p - 2
  log_dfs <- lchoose(n = p + k - 1, k = p) + log(2 + p / k)

  # nu_{p, k}
  if (!log) {

    log_dfs <- exp(log_dfs)

  }
  return(log_dfs)

}


#' @rdname Sobolev
#' @export
weights_dfs_Sobolev <- function(p, K_max = 1e3, thre = 1e-3, type,
                                Rothman_t = 1 / 3, Pycke_q = 0.5, Riesz_s = 1,
                                log = FALSE, verbose = TRUE, Gauss = TRUE,
                                N = 320, tol = 1e-6, force_positive = TRUE,
                                x_tail = NULL) {

  # alpha
  alpha <- 0.5 * p - 1

  # Sobolev weights and dfs
  if (p == 2 & type %in% c("Rothman", "Hermans_Rasson", "Pycke_q")) {

    if (type == "Rothman") {

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

      # Sequence of indexes
      k <- 1:K_max

      # log(v_{2 * k - 1}^2)
      log_v2km12 <- log(2 / pi) - log((2 * k - 1)^2)

      # log(v_{2 * k}^2)
      log_v2k2 <- log(5.79 / pi) - log((2 * k)^2 - 1)

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

      # log(v_k^2)
      log_vk2 <- switch((p > 2) + 1, log(8) - log(pi) - log(4 * k^2 - 1),
                        (p - 3) * log(2) + log(2 * k + p - 2) +
                          lgamma(k - 1 / 2) + lgamma(p / 2 - 1) +
                          lgamma(p / 2) - log(pi) - lgamma(k + p - 1 / 2))

      # Divide by 2 if p = 2 and (1 + k / alpha) if p > 2
      if (p == 2) {

        log_vk2 <- log_vk2 - log(2)

      } else {

        log_vk2 <- log_vk2 - log(1 + k / alpha)

      }

      # log(d_{p, k})
      log_dk <- d_p_k(p = p, k = k, log = TRUE)

      # Log weights and dfs
      log_weights <- log_vk2
      log_dfs <- log_dk

    } else if (type == "Riesz") {

      # Sequence of indexes
      k <- 1:K_max

      # tau_{s, k, p}
      if (p == 2) {

        tau <- 2^(Riesz_s + 1) / (1 + (k == 0))

      } else {

        tau <- 2^(p + Riesz_s - 3) * (p + 2 * k - 2) * gamma((p - 2) / 2)

      }

      # log(v_k^2) = log(b_{s, k, p})
      log_vk2 <- log(tau / sqrt(pi)) + lgamma((p - 1 + Riesz_s) / 2) +
        lgamma(-Riesz_s / 2 + k) - lgamma(p - 1 + Riesz_s / 2 + k) -
        lgamma(-Riesz_s / 2)

      # Divide by 2 if p = 2 and (1 + k / alpha) if p > 2
      if (p == 2) {

        log_vk2 <- log_vk2 - log(2)

      } else {

        log_vk2 <- log_vk2 - log(1 + k / alpha)

      }

      # log(d_{p, k})
      log_dk <- d_p_k(p = p, k = k, log = TRUE)

      # Log weights and dfs
      log_weights <- log_vk2
      log_dfs <- log_dk

    } else if (type == "PCvM") {

      # Sequence of indexes
      k <- 1:K_max

      # log(v_k^2)
      log_vk2 <- log(Gegen_coefs_Pn(k = k, p = p, type = "PCvM",
                                    Gauss = Gauss, N = N, tol = tol))

      # Divide by 2 if p = 2 and (1 + k / alpha) if p > 2
      if (p == 2) {

        log_vk2 <- log_vk2 - log(2)

      } else {

        log_vk2 <- log_vk2 - log(1 + k / alpha)

      }

      # log(d_{p, k})
      log_dk <- d_p_k(p = p, k = k, log = TRUE)

      # Log weights and dfs
      log_weights <- log_vk2
      log_dfs <- log_dk

    } else if (type == "PAD") {

      # Sequence of indexes
      k <- 1:K_max

      # log(v_k^2)
      log_vk2 <- log(Gegen_coefs_Pn(k = k, p = p, type = "PAD", 
                                    Gauss = Gauss, N = N, tol = tol))

      # Divide by 2 if p = 2 and (1 + k / alpha) if p > 2
      if (p == 2) {

        log_vk2 <- log_vk2 - log(2)

      } else {

        log_vk2 <- log_vk2 - log(1 + k / alpha)

      }

      # log(d_{p, k})
      log_dk <- d_p_k(p = p, k = k, log = TRUE)

      # Log weights and dfs
      log_weights <- log_vk2
      log_dfs <- log_dk

    } else if (type == "PRt") {

      # Sequence of indexes
      k <- 1:K_max

      # log(v_k^2)
      log_vk2 <- log(Gegen_coefs_Pn(k = k, p = p, type = "PRt",
                                    Rothman_t = Rothman_t, Gauss = Gauss,
                                    N = N, tol = tol))

      # Divide by 2 if p = 2 and (1 + k / alpha) if p > 2
      if (p == 2) {

        log_vk2 <- log_vk2 - log(2)

      } else {

        log_vk2 <- log_vk2 - log(1 + k / alpha)

      }

      # log(d_{p, k})
      log_dk <- d_p_k(p = p, k = k, log = TRUE)

      # Log weights and dfs
      log_weights <- log_vk2
      log_dfs <- log_dk

    } else {

      stop("Incompatible choice of p and type")

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
            " last weight = ", sprintf("%.3e", exp(log_weights[cutoff])), ")")

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
                      Riesz_s = 1, verbose = TRUE, N = 320, x_tail = NULL,
                      ...) {

  weights_dfs <- weights_dfs_Sobolev(p = p, K_max = K_max, thre = thre,
                                     type = type, Rothman_t = Rothman_t,
                                     Pycke_q = Pycke_q, Riesz_s = Riesz_s, 
                                     verbose = verbose, Gauss = TRUE, N = N,
                                     x_tail = x_tail)
  d_wschisq(x = x, weights = weights_dfs$weights, dfs = weights_dfs$dfs,
            method = method, ...)

}


#' @rdname Sobolev
#' @export
p_Sobolev <- function(x, p, type, method = c("I", "SW", "HBE", "MC")[1],
                      K_max = 1e3, thre = 1e-3, Rothman_t = 1 / 3,
                      Pycke_q = 0.5, Riesz_s = 1, verbose = TRUE,
                      N = 320, x_tail = NULL, ...) {

  weights_dfs <- weights_dfs_Sobolev(p = p, K_max = K_max, thre = thre,
                                     type = type, Rothman_t = Rothman_t,
                                     Pycke_q = Pycke_q, Riesz_s = Riesz_s,
                                     verbose = verbose, Gauss = TRUE,
                                     N = N, x_tail = x_tail)
  p_wschisq(x = x, weights = weights_dfs$weights, dfs = weights_dfs$dfs,
            method = method, ...)

}


#' @rdname Sobolev
#' @export
q_Sobolev <- function(u, p, type, method = c("I", "SW", "HBE", "MC")[1],
                      K_max = 1e3, thre = 1e-3, Rothman_t = 1 / 3,
                      Pycke_q = 0.5, Riesz_s = 1, verbose = TRUE, N = 320,
                      x_tail = NULL, ...) {

  weights_dfs <- weights_dfs_Sobolev(p = p, K_max = K_max, thre = thre,
                                     type = type, Rothman_t = Rothman_t,
                                     Pycke_q = Pycke_q, Riesz_s = Riesz_s,
                                     verbose = verbose, Gauss = TRUE, N = N,
                                     x_tail = x_tail)
  q_wschisq(u = u, weights = weights_dfs$weights, dfs = weights_dfs$dfs,
            method = method, ...)

}

