

#' @title Utilities for projected-ecdf statistics of spherical uniformity
#'
#' @description Computation of the kernels
#' \deqn{\psi_p^W(\theta) := \int_{-1}^1 A_x(\theta)\,\mathrm{d}W(F_p(x)),}{
#' \psi_p^W(\theta) := \int_{-1}^1 A_x(\theta) dW(F_p(x)),
#' }
#' where \eqn{A_x(\theta)} is the proportion of area surface of
#' \eqn{S^{p - 1}} covered by the
#' \link[=A_theta_x]{intersection of two hyperspherical caps} with common solid
#' angle \eqn{\pi - \cos^{-1}(x)} and centers separated by
#' an angle \eqn{\theta \in [0, \pi]}, \eqn{F_p} is the distribution function
#' of the \link[=p_proj_unif]{projected spherical uniform distribution},
#' and \eqn{W} is a measure on \eqn{[0, 1]}.
#'
#' Also, computation of the \link[=Gegen_coefs]{Gegenbauer coefficients} of
#' \eqn{\psi_p^W}:
#' \deqn{b_{k, p}^W := \frac{1}{c_{k, p}}\int_0^\pi \psi_p^W(\theta)
#' C_k^{p / 2 - 1}(\cos\theta)\,\mathrm{d}\theta.}{
#' b_{k, p}^W := \frac{1}{c_{k, p}} \int_0^\pi \psi_p^W(\theta)
#' C_k^(p / 2 - 1)(\cos\theta) d\theta.}
#' These coefficients can also be computed via
#' \deqn{b_{k, p}^W = \int_{-1}^1 a_{k, p}^x\,\mathrm{d}W(F_p(x))}{
#' b_{k, p}^W = \int_{-1}^1 a_{k, p}^x dW(F_p(x))}
#' for a certain function \eqn{x\rightarrow a_{k, p}^x}. They serve to define
#' \link[=locdev]{projected alternatives to uniformity}.
#' @inheritParams Gegenbauer
#' @param q integer giving the dimension of the sphere \eqn{S^q}.
#' @param type type of projected-ecdf test statistic. Must be either
#' \code{"PCvM"} (Cramér--von Mises), \code{"PAD"} (Anderson--Darling), or
#' \code{"PRt"} (Rothman).
#' @inheritParams unif_stat
#' @param tilde include the constant and bias term? Defaults to \code{FALSE}.
#' @param psi_Gauss use a \link[=Gauss_Legen_nodes]{Gauss--Legendre quadrature}
#' rule with \code{psi_N} nodes in the computation of the kernel function?
#' Defaults to \code{TRUE}.
#' @param psi_N number of points used in the Gauss--Legendre quadrature for
#' computing the kernel function. Defaults to \code{320}.
#' @param verbose flag to print informative messages. Defaults to \code{FALSE}.
#' @param x evaluation points for \eqn{a_{k, p}^x}, a vector with values in
#' \eqn{[-1, 1]}.
#' @param sqr return the \emph{signed} square root of \eqn{a_{k, p}^x}?
#' Defaults to \code{FALSE}.
#' @param K number of equispaced points on \eqn{[-1, 1]} used for evaluating
#' \eqn{f} and then interpolating. Defaults to \code{1e3}.
#' @inheritParams locdev
#' @return
#' \itemize{
#'   \item \code{psi_Pn}: a vector of length \code{length(theta)} with the
#'   evaluation of \eqn{\psi}.
#'   \item \code{Gegen_coefs_Pn}: a vector of size \code{length(k)} containing
#'   the coefficients \eqn{b_{k, p}^W}.
#'   \item \code{akx}: a matrix of size \code{c(length(x), length(k))}
#'   containing the coefficients \eqn{a_{k, p}^x}.
#'   \item \code{f_locdev_Pn}: the projected alternative \eqn{f} as a function
#'   ready to be evaluated.
#' }
#' @author Eduardo García-Portugués and Paula Navarro-Esteban.
#' @details
#' The evaluation of \eqn{\psi_p^W} and \eqn{b_{k, p}^W} depends on the type of
#' projected-ecdf statistic:
#' \itemize{
#' \item PCvM: closed-form expressions for \eqn{\psi_p^W} and \eqn{b_{k, p}^W}
#' with \eqn{p = 2, 3, 4}, numerical integration required for \eqn{p \ge 5}.
#' \item PAD: closed-form expressions for \eqn{\psi_2^W} and \eqn{b_{k, 3}^W},
#' numerical integration required for \eqn{\psi_p^W} with \eqn{p \ge 3} and
#' \eqn{b_{k, p}^W} with \eqn{p = 2} and \eqn{p \ge 4}.
#' \item PRt: closed-form expressions for \eqn{\psi_p^W} and \eqn{b_{k, p}^W}
#' for any \eqn{p \ge 2}.
#' }
#' See García-Portugués et al. (2020) for more details.
#' @references
#' García-Portugués, E., Navarro-Esteban, P., Cuesta-Albertos, J. A. (2020)
#' On a projection-based class of uniformity tests on the hypersphere.
#' \emph{arXiv:2008.09897}. \url{https://arxiv.org/abs/2008.09897}
#' @examples
#' # Kernels in the projected-ecdf test statistics
#' k <- 0:10
#' coefs <- list()
#' (coefs$PCvM <- t(sapply(2:5, function(p)
#'   Gegen_coefs_Pn(k = k, p = p, type = "PCvM"))))
#' (coefs$PAD <- t(sapply(2:5, function(p)
#'   Gegen_coefs_Pn(k = k, p = p, type = "PAD"))))
#' (coefs$PRt <- t(sapply(2:5, function(p)
#'   Gegen_coefs_Pn(k = k, p = p, type = "PRt"))))
#'
#' # Gegenbauer expansion
#' th <- seq(0, pi, length.out = 501)[-501]
#' old_par <- par(mfrow = c(3, 4))
#' for (type in c("PCvM", "PAD", "PRt")) {
#'
#'   for (p in 2:5) {
#'
#'     plot(th, psi_Pn(theta = th, q = p - 1, type = type), type = "l",
#'          main = paste0(type, ", p = ", p), xlab = expression(theta),
#'          ylab = expression(psi(theta)), axes = FALSE, ylim = c(-1.5, 1))
#'     axis(1, at = c(0, pi / 4, pi / 2, 3 * pi / 4, pi),
#'          labels = expression(0, pi / 4, pi / 2, 3 * pi / 4, pi))
#'     axis(2); box()
#'     lines(th, Gegen_series(theta = th, coefs = coefs[[type]][p - 1, ],
#'                            k = k, p = p), col = 2)
#'
#'   }
#'
#' }
#' par(old_par)
#'
#' # Analytical coefficients vs. numerical integration
#' test_coef <- function(type, p, k = 0:20) {
#'
#'   plot(k, log1p(abs(Gegen_coefs_Pn(k = k, p = p, type = type))),
#'        ylab = "Coefficients", main = paste0(type, ", p = ", p))
#'   points(k, log1p(abs(Gegen_coefs(k = k, p = p, psi = psi_Pn, type = type,
#'                                   q = p - 1))), col = 2)
#'   legend("topright", legend = c("log(1 + Gegen_coefs_Pn))",
#'                                 "log(1 + Gegen_coefs(psi_Pn))"),
#'          lwd = 2, col = 1:2)
#'
#' }
#'
#' # PCvM statistic
#' old_par <- par(mfrow = c(2, 2))
#' for (p in 2:5) test_coef(type = "PCvM", p = p)
#' par(old_par)
#'
#' # PAD statistic
#' old_par <- par(mfrow = c(2, 2))
#' for (p in 2:5) test_coef(type = "PAD", p = p)
#' par(old_par)
#'
#' # PRt statistic
#' old_par <- par(mfrow = c(2, 2))
#' for (p in 2:5) test_coef(type = "PRt", p = p)
#' par(old_par)
#'
#' # akx
#' akx(x = seq(-1, 1, l = 5), k = 1:4, p = 2)
#' akx(x = 0, k = 1:4, p = 3)
#'
#' # PRt alternative to uniformity
#' z <- seq(-1, 1, l = 1e3)
#' p <- c(2:5, 10, 15, 17)
#' col <- viridisLite::viridis(length(p))
#' plot(z, f_locdev_Pn(p = p[1], type = "PRt")(z), type = "s",
#'      col = col[1], ylim = c(0, 0.6), ylab = expression(f[Rt](z)))
#' for (k in 2:length(p)) {
#'   lines(z, f_locdev_Pn(p = p[k], type = "PRt")(z), type = "s", col = col[k])
#' }
#' legend("topleft", legend = paste("p =", p), col = col, lwd = 2)
#' @name Pn


#' @rdname Pn
#' @export
psi_Pn <- function(theta, q, type, Rothman_t = 1 / 3, tilde = FALSE,
                   psi_Gauss = TRUE, psi_N = 320, tol = 1e-6) {

  # Check dimension
  stopifnot(q >= 1)

  # Check theta
  stopifnot(all(0 <= theta & theta <= pi))

  # Type of statistic
  if (type == "PCvM") {

    # Compute psi_q
    ps <-
      switch(as.character(q),
             "1" = { # psi_1

               1 / 2 + theta * (theta - 2 * pi)  / (4 * pi^2)

             },
             "2" = { # psi_2

               1 / 2 - sin(theta / 2) / 4

             },
             "3" = { # psi_3

               psi_Pn(theta = theta, q = 1, type = "PCvM") +
                 ((pi - theta) * tan(theta / 2) - 2 * sin(theta / 2)^2) /
                 (4 * pi^2)

             },
             { # psi_q, q >= 4

               # Explicit part
               psi_q <- -3/4 + theta / (2 * pi) +
                 2 * drop(p_proj_unif(x = cos(theta / 2), p = q + 1))^2

               # Integral
               psi_q <- psi_q -
                 4 * sapply(theta, function(th) {

                   # Gauss--Legendre quadrature? Otherwise, call integrate
                   if (psi_Gauss) {

                     # Nodes and weights
                     t_k <- Gauss_Legen_nodes(a = 0, b = cos(th / 2), N = psi_N)
                     w_k <- Gauss_Legen_weights(a = 0, b = cos(th / 2),
                                                N = psi_N)

                     # Integrand evaluation
                     f_k <- exp(p_proj_unif(x = t_k, p = q + 1, log = TRUE) +
                                  p_proj_unif(x = tan(th / 2) * t_k /
                                                sqrt(1 - t_k^2),
                                              p = q, log = TRUE) +
                                  d_proj_unif(x = t_k, p = q + 1, log = TRUE))

                     # Remove non-finite data
                     f_k[!is.finite(f_k)] <- NA

                     # Integration
                     sum(w_k * f_k, na.rm = TRUE)

                   } else {

                     # Integrand
                     f <- function(t) {

                       f_k <- exp(p_proj_unif(x = t, p = q + 1, log = TRUE) +
                                    p_proj_unif(x = tan(th / 2) * t /
                                                  sqrt(1 - t^2),
                                                p = q, log = TRUE) +
                                    d_proj_unif(x = t, p = q + 1, log = TRUE))
                       f_k[!is.finite(f_k)] <- NA
                       return(f_k)

                     }

                     # Integrate
                     integrate(f = f, lower = 0, upper = cos(th / 2),
                               subdivisions = 1e3, rel.tol = tol, abs.tol = tol,
                               stop.on.error = FALSE)$value

                   }

                 })

             }
      )

    # Add constant?
    if (tilde) {

      ps <- ps - 1 / 3

    }

  } else if (type == "PAD") {

    # Compute psi_q
    ps <-
      switch(as.character(q),
             "1" = { # psi_1

               -2 * log(2 * pi) + (theta * log(theta) + (2 * pi - theta) *
                                     log(2 * pi - theta)) / pi

             },
             "2" = { # psi_2

               -log(4) + 2 / pi * sapply(theta, function(th) {

                 # Gauss--Legendre quadrature? Otherwise, call integrate
                 if (psi_Gauss) {

                   # Nodes and weights
                   t_k <- Gauss_Legen_nodes(a = 0, b = cos(th / 2), N = psi_N)
                   w_k <- Gauss_Legen_weights(a = 0, b = cos(th / 2), N = psi_N)

                   # Integrand evaluation
                   f_k <- log((1 + t_k) / (1 - t_k)) *
                     acos(tan(th / 2) * t_k / sqrt(1 - t_k^2))

                   # Remove non-finite data
                   f_k[!is.finite(f_k)] <- NA

                   # Integration
                   sum(w_k * f_k, na.rm = TRUE)

                 } else {

                   # Integrand
                   f <- function(t) {

                     f_k <- log((1 + t) / (1 - t)) *
                       acos(tan(th / 2) * t / sqrt(1 - t^2))
                     f_k[!is.finite(f_k)] <- NA
                     return(f_k)

                   }

                   # Integrate
                   integrate(f = f, lower = 0, upper = cos(th / 2),
                             subdivisions = 1e3, rel.tol = tol, abs.tol = tol,
                             stop.on.error = FALSE)$value

                 }

               })

             },
             "3" = { # psi_3

               s_theta <- theta - sin(theta)
               -2 * log(2 * pi) + (s_theta * log(s_theta) + (2 * pi - s_theta) *
                                     log(2 * pi - s_theta)) / pi -
                 4 / pi * tan(theta / 2) * sapply(theta, function(th) {

                 # Gauss--Legendre quadrature? Otherwise, call integrate
                 if (psi_Gauss) {

                   # Nodes and weights
                   t_k <- Gauss_Legen_nodes(a = 0, b = cos(th / 2), N = psi_N)
                   w_k <- Gauss_Legen_weights(a = 0, b = cos(th / 2), N = psi_N)

                   # Integrand evaluation
                   f_k <- t_k *
                     log(pi / (acos(t_k) - t_k * sqrt(1 - t_k^2)) - 1)

                   # Remove non-finite data
                   f_k[!is.finite(f_k)] <- NA

                   # Integration
                   sum(w_k * f_k, na.rm = TRUE)

                 } else {

                   # Integrand
                   f <- function(t) {

                     f_k <- t * log(pi / (acos(t) - t * sqrt(1 - t^2)) - 1)
                     f_k[!is.finite(f_k)] <- NA
                     return(f_k)

                   }

                   # Integrate
                   integrate(f = f, lower = 0, upper = cos(th / 2),
                             subdivisions = 1e3, rel.tol = tol, abs.tol = tol,
                             stop.on.error = FALSE)$value

                 }

               })

             },
             { # psi_q, q >= 4

               -log(4) + 4 * sapply(theta, function(th) {

                 # Gauss--Legendre quadrature? Otherwise, call integrate
                 if (psi_Gauss) {

                   # Nodes and weights
                   t_k <- Gauss_Legen_nodes(a = 0, b = cos(th / 2), N = psi_N)
                   w_k <- Gauss_Legen_weights(a = 0, b = cos(th / 2), N = psi_N)

                   # Integrand evaluation
                   log_F_k <- p_proj_unif(x = t_k, p = q + 1, log = TRUE) -
                     p_proj_unif(x = -t_k, p = q + 1, log = TRUE)
                   f_k <- log_F_k *
                     exp(p_proj_unif(x = -tan(th / 2) * t_k / sqrt(1 - t_k^2),
                                     p = q, log = TRUE) +
                           d_proj_unif(x = t_k, p = q + 1, log = TRUE))

                   # Remove non-finite data
                   f_k[!is.finite(f_k)] <- NA

                   # Integration
                   sum(w_k * f_k, na.rm = TRUE)

                 } else {

                   # Integrand
                   f <- function(t) {

                     log_F_k <- p_proj_unif(x = t, p = q + 1, log = TRUE) -
                       p_proj_unif(x = -t, p = q + 1, log = TRUE)
                     f_k <- log_F_k *
                       exp(p_proj_unif(x = -tan(th / 2) * t / sqrt(1 - t^2),
                                       p = q, log = TRUE) +
                             d_proj_unif(x = t, p = q + 1, log = TRUE))
                     f_k[!is.finite(f_k)] <- NA
                     return(f_k)

                   }

                   # Integrate
                   integrate(f = f, lower = 0, upper = cos(th / 2),
                             subdivisions = 1e3, rel.tol = tol, abs.tol = tol,
                             stop.on.error = FALSE)$value

                 }

               })

             }
      )

    # Set NaNs to 0, as they correspond to theta = 0
    ps[!is.finite(ps)] <- 0

    # Add constant?
    if (tilde) {

      ps <- ps + 1

    }

  } else if (type == "PRt") {

    # Separation angle
    t_m <- min(Rothman_t, 1 - Rothman_t)
    theta_t_m <- 2 * acos(drop(q_proj_unif(u = 1 - t_m, p = q + 1)))

    # Two cases
    psi_q <- numeric(length(theta))
    ind <- (theta >= theta_t_m)

    # Explicit part: theta >= theta_t_m
    psi_q[ind] <- 1 / 2 - t_m

    # Compute psi_q
    ps <-
      switch(as.character(q),
             "1" = { # psi_1

               -pmin(theta / (2 * pi) - t_m * (1 - t_m), t_m^2) +
                 1 / 2 - t_m * (1 - t_m)

             },
             "2" = { # psi_2

               # Longer part
               psi_q[!ind] <- -t_m + 1 / 2 - (1 - 2 * t_m) / pi *
                 acos((1 / 2 - t_m) / sqrt(t_m * (1 - t_m)) *
                        tan(theta[!ind] / 2)) +
                 atan(sqrt(cos(theta[!ind] / 2)^2 - (1 - 2 * t_m)^2) /
                        sin(theta[!ind] / 2)) / pi
               psi_q

             },
             "3" = { # psi_3

               # Longer part
               psi_q[!ind] <- 1/2 + t_m - (theta[!ind] + theta_t_m) / (2 * pi) +
                 (sin(theta_t_m) / 2 + tan(theta[!ind] / 2) *
                    cos(theta_t_m / 2)^2) / pi
               psi_q

             },
             { # psi_q, q >= 4

               # Integral part
               psi_q[!ind] <-
                 t_m - theta[!ind] / (2 * pi) +
                 2 * sapply(theta[!ind], function(th) {

                   # Gauss--Legendre quadrature? Otherwise, call integrate
                   if (psi_Gauss) {

                     # Nodes and weights
                     u_k <- Gauss_Legen_nodes(a = 0, b = cos(theta_t_m / 2),
                                              N = psi_N)
                     w_k <- Gauss_Legen_weights(a = 0, b = cos(theta_t_m / 2),
                                                N = psi_N)

                     # Integrand
                     f_k <- p_proj_unif(x = (u_k / sqrt(1 - u_k^2)) *
                                          tan(th / 2), p = q) *
                       d_proj_unif(x = u_k, p = q + 1)

                     # Remove non-finite data
                     f_k[!is.finite(f_k)] <- NA

                     # Integration
                     sum(w_k * f_k, na.rm = TRUE)

                   } else {

                     # Integrand
                     f <- function(u) {

                       f_k <- p_proj_unif(x = (u / sqrt(1 - u^2)) *
                                            tan(th / 2), p = q) *
                         d_proj_unif(x = u, p = q + 1)
                       f_k[!is.finite(f_k)] <- NA
                       return(f_k)

                     }

                     # Integrate
                     integrate(f = f, lower = 0, upper = cos(theta_t_m / 2),
                               subdivisions = 1e3, rel.tol = tol, abs.tol = tol,
                               stop.on.error = FALSE)$value

                   }

                 })
               psi_q

             }
      )

    # Add constant?
    if (tilde) {

      ps <- ps + t_m * (1 - t_m) - 1 / 2

    }

  } else {

    stop("type must be either \"PCvM\", \"PAD\", or \"PRt\".")

  }

  return(ps)

}


#' @rdname Pn
#' @export
Gegen_coefs_Pn <- function(k, p, type, Rothman_t = 1 / 3, Gauss = TRUE,
                           N = 320, tol = 1e-6, verbose = FALSE) {

  # Check dimension
  stopifnot(p >= 2)

  # Check orders
  stopifnot(all(k >= 0))

  # Coefficients
  coefs <- numeric(length(k))

  # Positive indexes
  ind1 <- k > 0
  k1 <- k[ind1]
  dok1 <- any(ind1)

  # Type of statistic
  if (type == "PCvM") {

    # k = 0
    coefs[!ind1] <- 1 / 3

    # k > 0
    if (dok1) {

      if (p == 2) {

        coefs[ind1] <- 1 / (k1 * pi)^2

      } else if (p == 3) {

        coefs[ind1] <- 1 / (8 * (k1 + 3/2) * (k1 - 1/2))

      } else if (p == 4) {

        coefs[ind1] <- (k1 == 1) * (35 / (72 * pi^2)) +
          (k1 > 1) / (2 * pi^2) * exp(log(3 * k1^2 + 6 * k1 + 4) - 2 * log(k1) -
                                        log(k1 + 1) - 2 * log(k1 + 2))

      } else {

        # Gauss--Legendre quadrature?
        if (Gauss) {

          # Nodes and weights
          w_k <- drop(Gauss_Legen_weights(a = -1, b = 1, N = N))
          x_k <- drop(Gauss_Legen_nodes(a = -1, b = 1, N = N))

          # k > 0
          coefs[ind1] <- colSums(w_k * akx(x = x_k, p = p, k = k1, sqr = FALSE) *
                                   drop(d_proj_unif(x = x_k, p = p)),
                                 na.rm = TRUE)

        } else {

          # k > 0
          coefs[ind1] <- sapply(k1, function(l) {
            f <- function(x) drop(akx(x = x, p = p, k = l, sqr = FALSE) *
                                    d_proj_unif(x = x, p = p))
            integrate(f = f, lower = -1, upper = 1,
                      subdivisions = 1e3, abs.tol = tol, rel.tol = tol,
                      stop.on.error = FALSE)$value
          })

        }
        if (verbose) {

          message("Gegenbauer coefficients of the ", type,
                  " statistic with p = ", p, " unknown in analytical form, ",
                  "using numerical integration of akx.")

        }

      }

    }

  } else if (type == "PAD") {

    # k = 0
    coefs[!ind1] <- -1

    # k > 0
    if (dok1) {

      if (p == 2) {

        # Gauss--Legendre quadrature?
        if (Gauss) {

          # Nodes and weights
          th_k <- drop(Gauss_Legen_nodes(a = 0, b = pi, N = N))
          w_k <- drop(Gauss_Legen_weights(a = 0, b = pi, N = N))

          # Integrand evaluation
          f_k <- (1 - cos(2 * th_k %o% k1)) / ((pi - th_k) * th_k)

          # k > 0
          coefs[ind1] <- colSums(w_k * f_k) / (pi * k1^2)

        } else {

          # k > 0
          coefs[ind1] <- sapply(k1, function(m) {
            f <- function(th) (1 - cos(2 * th * m)) / ((pi - th) * th)
            integrate(f = f, lower = 0 + tol, upper = pi - tol,
                      subdivisions = 1e3, abs.tol = tol, rel.tol = tol,
                      stop.on.error = FALSE)$value
          }) / (pi * k1^2)

        }

      } else if (p == 3) {

        # k > 0
        coefs[ind1] <- 1 / ((k1 + 1) * k1)

      } else {

        # Gauss--Legendre quadrature?
        if (Gauss) {

          # Nodes and weights
          stopifnot(N > 5)
          w_k <- drop(Gauss_Legen_weights(a = -1, b = 1, N = N))
          x_k <- drop(Gauss_Legen_nodes(a = -1, b = 1, N = N))

          # Weight evaluation
          dWF_k <- drop(p_proj_unif(x = x_k, p = p, log = TRUE))
          dWF_k <- rep(dWF_k[1:(N / 2)] + dWF_k[(N / 2 + 1):N], 2)
          dWF_k <- exp(drop(d_proj_unif(x = x_k, p = p, log = TRUE)) - dWF_k)
          dWF_k[!is.finite(dWF_k)] <- NA

          # k > 0
          coefs[ind1] <-
            colSums(w_k * akx(x = x_k, p = p, k = k1, sqr = FALSE) * dWF_k,
                    na.rm = TRUE)

        } else {

          # k > 0
          coefs[ind1] <- sapply(k1, function(l) {
            f <- function(x) {
              dWF_k <- p_proj_unif(x = x, p = p, log = TRUE) +
                p_proj_unif(x = -x, p = p, log = TRUE)
              dWF_k <- drop(exp(d_proj_unif(x = x, p = p, log = TRUE) - dWF_k))
              dWF_k[!is.finite(dWF_k)] <- NA
              drop(akx(x = x, p = p, k = l, sqr = FALSE) * dWF_k)
            }
            integrate(f = f, lower = -1, upper = 1, subdivisions = 1e3,
                      abs.tol = tol, rel.tol = tol, stop.on.error = FALSE)$value
          })

        }
        if (verbose) {

          message("Gegenbauer coefficients of the ", type,
                  " statistic with p = ", p, " unknown in analytical form, ",
                  "using numerical integration of akx.")

        }

      }

    }

  } else if (type == "PRt") {

    # k = 0
    t_m <- min(Rothman_t, 1 - Rothman_t)
    coefs[!ind1] <- 0.5 - t_m * (1 - t_m)

    # k > 0
    if (dok1) {

      # Call akx
      coefs[ind1] <- ifelse(p == 2, 2, 1) *
        drop(akx(x = drop(q_proj_unif(u = t_m, p = p)), p = p, k = k1,
                 sqr = FALSE))

    }

  } else {

    stop("type must be either \"PCvM\", \"PAD\", or \"PRt\".")

  }

  return(coefs)

}


#' @rdname Pn
#' @export
akx <- function(x, p, k, sqr = FALSE) {

  # Check dimension
  stopifnot(p >= 2)

  # Check orders
  stopifnot(all(k >= 0))

  # Coeffiecients
  a <- matrix(nrow = length(x), ncol = length(k))

  # Positive orders
  ind1 <- k > 0
  k1 <- k[ind1]
  dok1 <- any(ind1)

  # k = 0
  a[, !ind1] <- drop(p_proj_unif(x = x, p = p))

  # k > 0
  if (dok1) {

    if (p == 2) {

      a[, ind1] <- t(sin(k1 %o% acos(x)) / (k1 * pi))

    } else {

      q <- p - 1
      a[, ind1] <- t(Gegen_polyn(theta = acos(x), k = k1 - 1, p = p + 2))
      a[, ind1] <- a[, ind1] *
        t(exp(0.5 * log(1 + 2 * k1 / (q - 1)) +
                ((q - 1) * log(2) + 2 * lgamma((q + 1) / 2) +
                   lgamma(k1) - log(pi) - lgamma(k1 + q)) +
                rep(q / 2, length(k1)) %o% log(1 - x^2)))

    }

  }

  # Return sqrt(a) or a
  if (!sqr) {

    a <- a^2

  }
  return(a)

}


#' @rdname Pn
#' @export
f_locdev_Pn <- function(p, type, K = 1e3, N = 320, K_max = 1e4, thre = 1e-3,
                        Rothman_t = 1 / 3, verbose = FALSE) {

  # Exact f for PRt family
  if (type %in% c("PRt", "Rothman", "Ajne")) {

    t <- ifelse(type == "Ajne", 0.5, Rothman_t)
    x_t <- drop(q_proj_unif(u = t, p = p))
    f <- function(x) ((x >= x_t) + t) / rotasym::w_p(p = p)

  # Exact f for PCvM and p = 2
  } else if (type == "PCvM" & p == 2) {

    f <- function(x) (1 - log(2 * (1 - x)) * sqrt(2) / (2 * pi)) / (2 * pi)

  # Series expansion otherwise
  } else {

    uk <- cutoff_locdev(p = p, K_max = K_max, type = type, thre = thre,
                        Rothman_t = Rothman_t, verbose = verbose,
                        Gauss = TRUE, N = N)
    z <- seq(-1, 1, l = 1e4)
    f <- splinefun(x = z, y = pmax(f_locdev(z = z, p = p, uk = uk), 0))

  }
  return(f)

}

