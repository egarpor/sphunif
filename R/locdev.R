

#' @title Local projected alternatives to uniformity
#'
#' @description Density and random generation for local projected alternatives
#' to uniformity with densities
#' \deqn{f_{\kappa, \boldsymbol{\mu}}({\bf x}): =
#' \frac{1 - \kappa}{\omega_p} + \kappa f({\bf x}'\boldsymbol{\mu})}{
#' f_{\kappa, \mu}(x) = (1 - \kappa) / \omega_p + \kappa f(x'\mu)}
#' where
#' \deqn{f(z) = \frac{1}{\omega_p}\left\{1 + \sum_{k = 1}^\infty u_{k, p}
#' C_k^{p / 2 - 1}(z)\right\}}{f(x) = (1 / \omega_p)
#' \{1 + \sum_{k = 1}^\infty u_{k, p} C_k^(p / 2 - 1)(z)\}}
#' is the \emph{angular function} controlling the local alternative in a
#' \link[=Gegenbauer]{Gegenbauer series}, \eqn{0\le \kappa \le 1},
#' \eqn{\boldsymbol{\mu}}{\mu} is a direction on \eqn{S^{p - 1}}, and
#' \eqn{\omega_p} is the surface area of \eqn{S^{p - 1}}. The sequence
#' \eqn{\{u_{k, p}\}} is typically such that
#' \eqn{u_{k, p} = \left(1 + \frac{2k}{p - 2}\right) b_{k, p}}{
#' u_{k, p} = (1 + 2k / (p - 2)) b_{k, p}} for the Gegenbauer coefficients
#' \eqn{\{b_{k, p}\}} of the kernel function of a Sobolev statistic (see the
#' \link[=Sobolev_coefs]{transformation} between the coefficients \eqn{u_{k, p}}
#' and \eqn{b_{k, p}}).
#'
#' Also, automatic truncation of the series \eqn{\sum_{k = 1}^\infty u_{k, p}
#' C_k^{p / 2 - 1}(z)}{\sum_{k = 1}^\infty u_{k, p} C_k^(p / 2 - 1)(z)}
#' according to the proportion of \link[=Gegenbauer]{"Gegenbauer norm"}
#' explained.
#'
#' @param z projected evaluation points for \eqn{f}, a vector with entries on
#' \eqn{[-1, 1]}.
#' @inheritParams Sobolev_coefs
#' @inheritParams rotasym::d_tang_norm
#' @param mu a unit norm vector of size \code{p} giving the axis of rotational
#' symmetry.
#' @param f angular function defined on \eqn{[-1, 1]}. Must be vectorized.
#' @param kappa the strength of the local alternative, between \code{0}
#' and \code{1}.
#' @inheritParams r_unif
#' @param F_inv quantile function associated to \eqn{f}. Computed by
#' \code{\link{F_inv_from_f}} if \code{NULL} (default).
#' @inheritParams Gegenbauer
#' @param ... further parameters passed to \code{\link{F_inv_from_f}}.
#' @param K_max integer giving the truncation of the series. Defaults to
#' \code{1e4}.
#' @param thre proportion of norm \emph{not} explained by the first terms of the
#' truncated series. Defaults to \code{1e-3}.
#' @inheritParams Sobolev
#' @param verbose output information about the truncation (\code{TRUE} or
#' \code{1}) and a diagnostic plot (\code{2})? Defaults to \code{FALSE}.
#' @return
#' \itemize{
#'   \item \code{f_locdev}: angular function evaluated at \code{x}, a vector.
#'   \item \code{con_f}: normalizing constant \eqn{c_f} of \eqn{f}, a scalar.
#'   \item \code{d_locdev}: density function evaluated at \code{x}, a vector.
#'   \item \code{r_locdev}: a matrix of size \code{c(n, p)} containing a random
#'   sample from the density \eqn{f_{\kappa, \boldsymbol{\mu}}}{
#'   f_{\kappa, \mu}}.
#'   \item \code{cutoff_locdev}: vector of coefficients \eqn{\{u_{k, p}\}}
#'   automatically truncated according to \code{K_max} and \code{thre}
#'   (see details).
#' }
#' @details
#' See the definitions of local alternatives in Prentice (1978) and in
#' García-Portugués et al. (2023).
#'
#' The truncation of \eqn{\sum_{k = 1}^\infty u_{k, p} C_k^{p / 2 - 1}(z)}{
#' \sum_{k = 1}^\infty u_{k, p} C_k^(p / 2 - 1)(z)} is done to the first
#' \code{K_max} terms and then up to the index such that the first terms
#' leave unexplained the proportion \code{thre} of the norm of the whole series.
#' Setting \code{thre = 0} truncates to \code{K_max} terms exactly. If the
#' series only contains odd or even non-zero terms, then only \code{K_max / 2}
#' addends are \emph{effectively} taken into account in the first truncation.
#' @references
#' García-Portugués, E., Navarro-Esteban, P., Cuesta-Albertos, J. A. (2023)
#' On a projection-based class of uniformity tests on the hypersphere.
#' \emph{Bernoulli}, 29(1):181--204. \doi{10.3150/21-BEJ1454}.
#'
#' Prentice, M. J. (1978). On invariant tests of uniformity for directions and
#' orientations. \emph{The Annals of Statistics}, 6(1):169--176.
#' \doi{10.1214/aos/1176344075}
#' @examples
#' ## Local alternatives diagnostics
#'
#' loc_alt_diagnostic  <- function(p, type, thre = 1e-3, K_max = 1e3) {
#'
#'   # Coefficients of the alternative
#'   uk <- cutoff_locdev(K_max = K_max, p = p, type = type, thre = thre,
#'                       N = 640)
#'
#'   old_par <- par(mfrow = c(2, 2))
#'
#'   # Construction of f
#'   z <- seq(-1, 1, l = 1e3)
#'   f <- function(z) f_locdev(z = z, p = p, uk = uk)
#'   plot(z, f(z), type = "l", xlab = expression(z), ylab = expression(f(z)),
#'        main = paste0("Local alternative f, ", type, ", p = ", p), log = "y")
#'
#'   # Projected density on [-1, 1]
#'   f_proj <- function(z) rotasym::w_p(p = p - 1) * f(z) *
#'     (1 - z^2)^((p - 3) / 2)
#'   plot(z, f_proj(z), type = "l", xlab = expression(z),
#'        ylab = expression(omega[p - 1] * f(z) * (1 - z^2)^{(p - 3) / 2}),
#'        main = paste0("Projected density, ", type, ", p = ", p), log = "y",
#'        sub = paste("Integral:", round(con_f(f = f, p = p), 4)))
#'
#'   # Quantile function for projected density
#'   mu <- c(rep(0, p - 1), 1)
#'   F_inv <- F_inv_from_f(f = f, p = p, K = 5e2)
#'   plot(F_inv, xlab = expression(x), ylab = expression(F^{-1}*(x)),
#'        main = paste0("Quantile function, ", type, ", p = ", p))
#'
#'   # Sample from the alternative and plot the projected sample
#'   n <- 5e4
#'   samp <- r_locdev(n = n, mu = mu, f = f, kappa = 1, F_inv = F_inv)
#'   plot(z, f_proj(z), col = 2, type = "l",
#'        main = paste0("Simulated projected data, ", type, ", p = ", p),
#'        ylim = c(0, 1.75))
#'   hist(samp %*% mu, freq = FALSE, breaks = seq(-1, 1, l = 50), add = TRUE)
#'
#'   par(old_par)
#'
#' }
#' \donttest{
#' ## Local alternatives for the PCvM test
#'
#' loc_alt_diagnostic(p = 2, type = "PCvM")
#' loc_alt_diagnostic(p = 3, type = "PCvM")
#' loc_alt_diagnostic(p = 4, type = "PCvM")
#' loc_alt_diagnostic(p = 5, type = "PCvM")
#' loc_alt_diagnostic(p = 11, type = "PCvM")
#'
#' ## Local alternatives for the PAD test
#'
#' loc_alt_diagnostic(p = 2, type = "PAD")
#' loc_alt_diagnostic(p = 3, type = "PAD")
#' loc_alt_diagnostic(p = 4, type = "PAD")
#' loc_alt_diagnostic(p = 5, type = "PAD")
#' loc_alt_diagnostic(p = 11, type = "PAD")
#'
#' ## Local alternatives for the PRt test
#'
#' loc_alt_diagnostic(p = 2, type = "PRt")
#' loc_alt_diagnostic(p = 3, type = "PRt")
#' loc_alt_diagnostic(p = 4, type = "PRt")
#' loc_alt_diagnostic(p = 5, type = "PRt")
#' loc_alt_diagnostic(p = 11, type = "PRt")
#' }
#' @name locdev


#' @rdname locdev
#' @export
f_locdev <- function(z, p, uk) {

  # Check dimension
  stopifnot(p >= 2)

  # Unnormalized local alternative
  f <- 1 + Gegen_series(theta = acos(z), coefs = uk, p = p, k = seq_along(uk))

  # Normalize by \omega_p such that
  # \omega_{p - 1} * f(z) * (1 - z^2)^((p - 3) / 2) integrates one in [-1, 1]
  f <- f / rotasym::w_p(p = p)
  return(f)

}


#' @rdname locdev
#' @export
con_f <- function(f, p, N = 320) {

  # Gauss--Legendre nodes and weights
  th_k <- drop(Gauss_Legen_nodes(a = 0, b = pi, N = N))
  w_k <- drop(Gauss_Legen_weights(a = 0, b = pi, N = N))

  # Integral of w_{p - 1} * f(z) * (1 - z^2)^((p - 3) / 2) in [-1, 1],
  # using z = cos(theta) for theta in [0, pi]
  int <- rotasym::w_p(p = p - 1) *
    sum(w_k * f(cos(th_k)) * sin(th_k)^(p - 2), na.rm = TRUE)
  return(1 / int)

}


#' @rdname locdev
#' @export
d_locdev <- function(x, mu, f, kappa) {

  # Check dimension
  if (is.null(dim(x))) {

    x <- rbind(x)

  }
  stopifnot(ncol(x) == length(mu))

  # Check kappa
  stopifnot(0 <= kappa & kappa <= 1)

  # Alternative density
  if (kappa > 0) {

    f1 <- rotasym::d_tang_norm(x = x, theta = mu,
                               d_U = rotasym::d_unif_sphere,
                               g_scaled = function(z, log = TRUE) log(f(z)))

  } else {

    f1 <- 0

  }

  # Uniform density
  f0 <- rotasym::d_unif_sphere(x = x)

  # Merge densities
  return((1 - kappa) * f0 + kappa * f1)

}


#' @rdname locdev
#' @export
r_locdev <- function(n, mu, f, kappa, F_inv = NULL, ...) {

  # Dimension
  p <- length(mu)

  # Check kappa
  stopifnot(0 <= kappa & kappa <= 1)
  if (kappa == 0) {

    return(r_unif_sph(n = n, p = p, M = 1)[, , 1])

  }

  # Compute the inverse of the distribution function F?
  if (is.null(F_inv)) {

    F_inv <- F_inv_from_f(f = f, p = p, ...)

  }

  # Sample object
  samp <- matrix(0, nrow = n, ncol = p)
  ind_1 <- runif(n = n) <= kappa
  n_1 <- sum(ind_1)

  # Sample under the alternative
  if (n_1 > 0) {

    r_V <- function(n) F_inv(runif(n = n))
    r_U <- function(n) r_unif_sph(n = n, p = p - 1, M = 1)[, , 1]
    samp[ind_1, ] <- rotasym::r_tang_norm(n = n_1, theta = mu,
                                          r_V = r_V, r_U = r_U)

  }

  # Sample under the null
  if (n_1 < n) {

    samp[!ind_1, ] <- r_unif_sph(n = n - n_1, p = p, M = 1)[, , 1]

  }

  # Sample
  return(samp)

}


#' @rdname locdev
#' @export
cutoff_locdev <- function(p, K_max = 1e4, thre = 1e-3, type, Rothman_t = 1 / 3,
                          Pycke_q = 0.5, verbose = FALSE, Gauss = TRUE, N = 320,
                          tol = 1e-6) {

  # vk2
  vk2 <- weights_dfs_Sobolev(p, K_max = K_max, thre = 0, type = type,
                             Rothman_t = Rothman_t, Pycke_q = Pycke_q,
                             log = FALSE, verbose = FALSE, Gauss = Gauss,
                             N = N, tol = tol)$weights
  K_max_new <- length(vk2)

  # Signs
  if (type %in% c("PRt", "Rothman", "Ajne")) {

    x_t <- drop(q_proj_unif(u = ifelse(type == "Ajne", 0.5, Rothman_t), p = p))
    signs <- akx(x = x_t, p = p, k = seq_len(K_max_new), sqr = TRUE)

  } else {

    if (verbose > 1) {

      message("Signs unknown for the ", type,
              " statistic, using positive signs experimentally.")

    }
    signs <- 1

  }

  # uk
  uk <- vk2_to_uk(vk2 = vk2, p = p, signs = signs)

  # Cutoff based on the explained squared norm
  cum_norm <- Gegen_norm(coefs = uk, k = seq_len(K_max_new), p = p,
                         cumulative = TRUE)^2
  cum_norm <- cum_norm / cum_norm[K_max_new]
  cutoff <- which(cum_norm >= 1 - thre)[1]

  # Truncate displaying optional information
  uk_cutoff <- uk[1:cutoff]
  if (verbose) {

    message("Series truncated from ", K_max_new, " to ", cutoff,
            " terms (", 100 * (1 - thre),
            "% of cumulated norm; last coefficient = ",
            sprintf("%.3e", uk_cutoff[cutoff]), ").")

    # Diagnostic plots
    if (verbose > 1) {

      old_par <- par(mfrow = c(1, 2), mar = c(5, 5.5, 4, 2) + 0.1)

      # Cumulated norm
      plot(seq_len(K_max), 100 * c(cum_norm, rep(1, K_max - K_max_new)),
           xlab = "k", ylab = "Percentage of cumulated squared norm",
           type = "s", log = "x")
      segments(x0 = cutoff, y0 = par()$usr[3],
               x1 = cutoff, y1 = 100 * cum_norm[cutoff], col = 3)
      segments(x0 = 1, y0 = 100 * (1 - thre),
               x1 = cutoff, y1 = 100 * (1 - thre), col = 2)
      abline(v = K_max_new, col = "gray", lty = 2)

      # Function and truncation
      z <- seq(-1, 1, l = 1e3)[-c(1, 1e3)]
      th <- acos(z)
      G1 <- Gegen_series(theta = th, coefs = c(1, uk),
                         k = c(0, seq_along(uk)), p = p)
      G2 <- Gegen_series(theta = th, coefs = c(1, uk_cutoff),
                         k = c(0, seq_along(uk_cutoff)), p = p)
      e <- expression(f(z) == 1 + sum(u[k] * C[k]^(p / 2 - 1) * (z), k == 1, K))
      plot(z, G1, ylim = c(1e-3, max(c(G1, G2))), xlab = expression(z),
           ylab = e, type = "l", log = "y")
      lines(z, G2, col = 2)
      legend("top", legend = paste("K =", c(K_max_new, cutoff)),
             col = 1:2, lwd = 2)
      par(old_par)

    }

  }
  return(uk_cutoff)

}


#' @title Distribution and quantile functions from angular function
#'
#' @description Numerical computation of the distribution function \eqn{F} and
#' the quantile function \eqn{F^{-1}} for an \link[=locdev]{angular function}
#' \eqn{f} in a \link[rotasym:d_tang_norm]{tangent-normal decomposition}.
#' \eqn{F^{-1}(x)} results from the inversion of
#' \deqn{F(x) = \int_{-1}^x \omega_{p - 1}c_f f(z) (1 - z^2)^{(p - 3) / 2}
#' \,\mathrm{d}z}{F(x) = \int_{-1}^x \omega_{p - 1}c_f f(z)
#' (1 - z^2)^{(p - 3) / 2} dz}
#' for \eqn{x\in [-1, 1]}, where \eqn{c_f} is a normalizing constant and
#' \eqn{\omega_{p - 1}} is the surface area of \eqn{S^{p - 2}}.
#'
#' @inheritParams locdev
#' @inheritParams r_unif
#' @param Gauss use a \link[=Gauss_Legen_nodes]{Gauss--Legendre quadrature}
#' rule to integrate \eqn{f} with \code{N} nodes? Otherwise, rely on
#' \code{\link{integrate}} Defaults to \code{TRUE}.
#' @param N number of points used in the Gauss--Legendre quadrature. Defaults
#' to \code{320}.
#' @param K number of equispaced points on \eqn{[-1, 1]} used for evaluating
#' \eqn{F^{-1}} and then interpolating. Defaults to \code{1e3}.
#' @param tol tolerance passed to \code{\link{uniroot}} for the inversion of
#' \eqn{F}. Also, passed to \code{\link{integrate}}'s \code{rel.tol} and
#' \code{abs.tol} if \code{Gauss = FALSE}. Defaults to \code{1e-6}.
#' @param ... further parameters passed to \code{f}.
#' @details
#' The normalizing constant \eqn{c_f} is such that \eqn{F(1) = 1}. It does not
#' need to be part of \code{f} as it is computed internally.
#'
#' Interpolation is performed by a monotone cubic spline. \code{Gauss = TRUE}
#' yields more accurate results, at expenses of a heavier computation.
#'
#' If \code{f} yields negative values, these are silently truncated to zero.
#' @return A \code{\link{splinefun}} object ready to evaluate \eqn{F} or
#' \eqn{F^{-1}}, as specified.
#' @examples
#' f <- function(x) rep(1, length(x))
#' plot(F_from_f(f = f, p = 4, Gauss = TRUE), ylab = "F(x)", xlim = c(-1, 1))
#' plot(F_from_f(f = f, p = 4, Gauss = FALSE), col = 2, add = TRUE,
#'      xlim = c(-1, 1))
#' curve(p_proj_unif(x = x, p = 4), col = 3, add = TRUE, n = 300)
#' plot(F_inv_from_f(f = f, p = 4, Gauss = TRUE), ylab = "F^{-1}(x)")
#' plot(F_inv_from_f(f = f, p = 4, Gauss = FALSE), col = 2, add = TRUE)
#' curve(q_proj_unif(u = x, p = 4), col = 3, add = TRUE, n = 300)
#' @name F_from_f


#' @rdname F_from_f
#' @export
F_from_f <- function(f, p, Gauss = TRUE, N = 320, K = 1e3, tol = 1e-6, ...) {

  # Integration grid
  z <- seq(-1, 1, length.out = K)

  # Integrated \omega_{p - 1} * f(z) * (1 - z^2)^{(p - 3) / 2}
  if (Gauss) {

    # Using Gauss--Legendre quadrature but ensuring monotonicity?
    F_grid <- rotasym::w_p(p = p - 1) * sapply(z, function(t) {
      z_k <- drop(Gauss_Legen_nodes(a = -1, b = t, N = N))
      w_k <- drop(Gauss_Legen_weights(a = -1, b = t, N = N))
      sum(w_k * exp(log(pmax(f(z_k, ...), 0)) + ((p - 3) / 2) * log1p(-z_k^2)),
          na.rm = TRUE)
    })

    # Normalize f (the normalizing constant may not be included in f)
    c_f <- F_grid[length(F_grid)]
    F_grid <- F_grid / c_f

  } else {

    # Using integrate
    g <- function(t) exp(log(pmax(f(t, ...), 0)) + ((p - 3) / 2) * log1p(-t^2))
    F_grid <- sapply(z[-1], function(u) rotasym::w_p(p = p - 1) *
                       integrate(f = g, lower = -1, upper = u,
                                 subdivisions = 1e3, rel.tol = tol,
                                 abs.tol = tol, stop.on.error = FALSE)$value)

    # Normalize f (the normalizing constant may not be included in f)
    c_f <- F_grid[length(F_grid)]
    F_grid <- F_grid / c_f

    # Add 0
    F_grid <- c(0, F_grid)

  }

  # Use method = "hyman" for monotone interpolations if possible
  if (anyNA(F_grid)) stop("Numerical error (NAs) in F_grid")
  F_appf <- switch(is.unsorted(F_grid) + 1,
                   splinefun(x = z, y = F_grid, method = "hyman"),
                   approxfun(x = z, y = F_grid, method = "linear", rule = 2))
  return(F_appf)

}


#' @rdname F_from_f
#' @export
F_inv_from_f <- function(f, p, Gauss = TRUE, N = 320, K = 1e3, tol = 1e-6,
                         ...) {

  # Approximate F
  F_appf <- F_from_f(f = f, p = p, Gauss = Gauss, N = N, K = K, tol = tol, ...)

  # Inversion of F
  u <- seq(0, 1, length.out = K)
  F_inv_grid <- sapply(u[-c(1, K)], function(v) {
    uniroot(f = function(x) F_appf(x) - v, lower = -1, upper = 1,
            tol = tol)$root
  })
  F_inv_grid <- c(-1, F_inv_grid, 1)

  # Use method = "hyman" for monotone interpolations if possible
  stopifnot(!anyNA(F_inv_grid))
  F_inv <- switch(is.unsorted(F_inv_grid) + 1,
                  splinefun(x = u, y = F_inv_grid, method = "hyman"),
                  approxfun(x = u, y = F_inv_grid, method = "linear",
                            rule = 2))
  return(F_inv)

}
