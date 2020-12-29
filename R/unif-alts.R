

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
#'   \item \code{cte_f}: normalizing constant \eqn{c_f} of \eqn{f}, a scalar.
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
#' García-Portugués et al. (2020).
#'
#' The truncation of \eqn{\sum_{k = 1}^\infty u_{k, p} C_k^{p / 2 - 1}(z)}{
#' \sum_{k = 1}^\infty u_{k, p} C_k^(p / 2 - 1)(z)} is done to the first
#' \code{K_max} terms and then up to the index such that the first terms
#' leave unexplained the proportion \code{thre} of the norm of the whole series.
#' Setting \code{thre = 0} truncates to \code{K_max} terms exactly. If the
#' series only contains odd or even non-zero terms, then only \code{K_max / 2}
#' addends are \emph{effectively} taken into account in the first truncation.
#' @references
#' García-Portugués, E., Navarro-Esteban, P., Cuesta-Albertos, J. A. (2020)
#' On a projection-based class of uniformity tests on the hypersphere.
#' \emph{arXiv:2008.09897}. \url{https://arxiv.org/abs/2008.09897}
#'
#' Prentice, M. J. (1978). On invariant tests of uniformity for directions and
#' orientations. \emph{The Annals of Statistics}, 6(1):169--176.
#' \url{https://doi.org/10.1214/aos/1176344075}
#' @examples
#' ## Local alternatives diagnostics
#'
#' loc_alt_diagnostic  <- function(p, type, thre = 1e-3, K_max = 1e3) {
#'
#'   # Coefficients of the alternative
#'   uk <- cutoff_locdev(K_max = K_max, p = p, type = type, thre = thre,
#'                       N = 640)
#'
#'   par_old <- par(mfrow = c(2, 2))
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
#'        sub = paste("Integral:", round(cte_f(f = f, p = p), 4)))
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
#'   par(par_old)
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
cte_f <- function(f, p, N = 320) {

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

      message("signs unknown for the ", type,
              " statistic, using positive signs experimentally.")

    }
    signs <- 1

  }

  # uk
  uk <- vk2_to_uk(vk2 = vk2, p = p, signs = signs)

  # Cutoff based on the explained squared norm
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
            sprintf("%.3e", uk_cutoff[cutoff]), ")")

    # Diagnostic plots
    if (verbose > 1) {

      par_old <- par(mfrow = c(1, 2), mar = c(5, 5.5, 4, 2) + 0.1)

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
      plot(z, G1, ylim = c(1e-3, max(c(G1, G2))), xlab = expression(z),
           ylab = expression(f(z) == 1 + sum(u[k] * C[k]^{(p / 2 - 1)} * (z),
                                             k == 1, K)), type = "l", log = "y")
      lines(z, G2, col = 2)
      legend("top", legend = paste("K =", c(K_max_new, cutoff)),
             col = 1:2, lwd = 2)
      par(par_old)

    }

  }
  return(uk_cutoff)

}


#' @title Quantile function from angular function
#'
#' @description Numerical computation of the quantile function \eqn{F^{-1}}
#' for an \link[=locdev]{angular function} \eqn{f} in a
#' \link[=tang-norm-decomp]{tangent-normal decomposition}. \eqn{F^{-1}(x)}
#' results from the inversion of
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
#' @return A \code{\link{splinefun}} object ready to evaluate \eqn{F^{-1}}.
#' @examples
#' f <- function(x) rep(1, length(x))
#' plot(F_inv_from_f(f = f, p = 4, Gauss = TRUE), ylab = "F^{-1}(x)")
#' plot(F_inv_from_f(f = f, p = 4, Gauss = FALSE), col = 2, add = TRUE)
#' curve(q_proj_unif(u = x, p = 4), col = 3, add = TRUE, n = 300)
#' @export
F_inv_from_f <- function(f, p, Gauss = TRUE, N = 320, K = 1e3, tol = 1e-6,
                         ...) {

  # Integration grid
  z <- seq(-1, 1, length.out = K)

  # Integrated \omega_{p - 1} * f(z) * (1 - z^2)^{(p - 3) / 2}
  if (Gauss) {

    # Using Gauss--Legendre quadrature but ensuring monotnonicity?
    y <- rotasym::w_p(p = p - 1) * sapply(z, function(t) {
      z_k <- drop(Gauss_Legen_nodes(a = -1, b = t, N = N))
      w_k <- drop(Gauss_Legen_weights(a = -1, b = t, N = N))
      sum(w_k * pmax(f(z_k, ...), 0) * (1 - z_k^2)^((p - 3) / 2), na.rm = TRUE)
    })

    # Normalize f (the normalizing constant may not be included in f)
    c_f <- y[length(y)]
    y <- y / c_f

  } else {

    # Using integrate
    g <- function(t) pmax(f(t, ...), 0) * (1 - t^2)^((p - 3) / 2)
    y <- sapply(z[-1], function(u) rotasym::w_p(p = p - 1) *
                  integrate(f = g, lower = -1, upper = u, subdivisions = 1e3,
                            rel.tol = tol, abs.tol = tol,
                            stop.on.error = FALSE)$value)

    # Normalize f (the normalizing constant may not be included in f)
    c_f <- y[length(y)]
    y <- y / c_f

    # Add 0
    y <- c(0, y)

  }

  # Use method = "hyman" for monotone interpolations if possible
  F_appf <- switch(is.unsorted(y) + 1,
                   splinefun(x = z, y = y, method = "hyman"),
                   approxfun(x = z, y = y, method = "linear", rule = 2))

  # Inversion of F
  u <- seq(0, 1, length.out = K)
  F_inv_grid <- sapply(u[-c(1, K)], function(v) {
    uniroot(f = function(x) F_appf(x) - v, lower = -1, upper = 1,
            tol = tol)$root
  })
  F_inv_grid <- c(-1, F_inv_grid, 1)

  # Use method = "hyman" for monotone interpolations if possible
  F_inv <- switch(is.unsorted(F_inv_grid) + 1,
                  splinefun(x = u, y = F_inv_grid, method = "hyman"),
                  approxfun(x = u, y = F_inv_grid, method = "linear",
                            rule = 2))
  return(F_inv)

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
#' \code{1:length(bk)}, a vector.
#' @param vk2 \bold{squared} coefficients \eqn{v_{k, p}^2} associated to the
#' indexes \code{1:length(vk2)}, a vector.
#' @param uk coefficients \eqn{u_{k, p}} associated to the indexes
#' \code{1:length(uk)}, a vector.
#' @inheritParams r_unif_sph
#' @param signs signs of the coefficients \eqn{u_{k, p}}, a vector of the
#' same size as \code{vk2} or \code{bk}, or a scalar. Defaults to \code{1}.
#' @return
#' The corresponding vectors of coefficients \code{vk2}, \code{bk} or \code{uk},
#' depending on the call.
#' @details
#' See more details in Prentice (1978) and García-Portugués et al. (2020). The
#' adequate signs of \code{uk} for the \code{"PRt"} \link[=Pn]{Rothman test}
#' can be retrieved with \code{\link{akx}} and \code{sqr = TRUE}, see the
#' examples.
#' @references
#' García-Portugués, E., Navarro-Esteban, P., Cuesta-Albertos, J. A. (2020)
#' On a projection-based class of uniformity tests on the hypersphere.
#' \emph{arXiv:2008.09897}. \url{https://arxiv.org/abs/2008.09897}
#'
#' Prentice, M. J. (1978). On invariant tests of uniformity for directions and
#' orientations. \emph{The Annals of Statistics}, 6(1):169--176.
#' \url{https://doi.org/10.1214/aos/1176344075}
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
bk_to_vk2 <- function(bk, p) {

  # Check dimension
  p <- as.integer(p)
  stopifnot(p >= 2)

  # Add factor
  if (p == 2) {

    return(bk / 2)

  } else {

    return(bk / (1 + 2 * seq_along(bk) / (p - 2)))

  }

}


#' @rdname Sobolev_coefs
#' @export
bk_to_uk <- function(bk, p, signs = 1) {

  # Check dimension
  p <- as.integer(p)
  stopifnot(p >= 2)

  # Check signs
  stopifnot(length(signs) %in% c(1, length(bk)))

  # Add factor
  if (p == 2) {

    return(sign(signs) * sqrt(2 * bk))

  } else {

    return(sign(signs) * sqrt((1 + 2 * seq_along(bk) / (p - 2)) * bk))

  }

}


#' @rdname Sobolev_coefs
#' @export
vk2_to_bk <- function(vk2, p) {

  # Check dimension
  p <- as.integer(p)
  stopifnot(p >= 2)

  # Add factor
  if (p == 2) {

    return(2 * vk2)

  } else {

    return((1 + 2 * seq_along(vk2) / (p - 2)) * vk2)

  }

}


#' @rdname Sobolev_coefs
#' @export
vk2_to_uk <- function(vk2, p, signs = 1) {

  # Check dimension
  p <- as.integer(p)
  stopifnot(p >= 2)

  # Check signs
  stopifnot(length(signs) %in% c(1, length(vk2)))

  # Add factor
  if (p == 2) {

    return(2 * sign(signs) * sqrt(vk2))

  } else {

    return((1 + 2 * seq_along(vk2) / (p - 2)) * sign(signs) * sqrt(vk2))

  }

}


#' @rdname Sobolev_coefs
#' @export
uk_to_vk2 <- function(uk, p) {

  # Check dimension
  p <- as.integer(p)
  stopifnot(p >= 2)

  # Add factor
  if (p == 2) {

    return((uk / 2)^2)

  } else {

    return((uk / (1 + 2 * seq_along(uk) / (p - 2)))^2)

  }

}


#' @rdname Sobolev_coefs
#' @export
uk_to_bk <- function(uk, p) {

  # Check dimension
  p <- as.integer(p)
  stopifnot(p >= 2)

  # Add factor
  if (p == 2) {

    return(uk^2 / 2)

  } else {

    return(uk^2 / (1 + 2 * seq_along(uk) / (p - 2)))

  }

}


#' @title Sample non-uniformly distributed spherical data
#'
#' @description Simple simulation of prespecified non-uniform spherical
#' distributions: von Mises--Fisher (vMF), Mixture of vMF (MvMF),
#' Angular Central Gaussian (ACG), Small Circle (SC) or Watson (W).
#'
#' @inheritParams r_unif
#' @param scenario simulation scenario, must be \code{"vMF"}, \code{"MvMF"},
#' \code{"ACG"}, \code{"SC"} or \code{"W"}. See details below.
#' @param kappa non-negative parameter measuring the strength of the deviation
#' with respect to uniformity.
#' @param nu projection along \eqn{{\bf e}_p}{e_p} controlling the modal
#' strip of the small circle distribution. Must be in (-1, 1). Defaults to
#' \code{0.5}.
#' @param F_inv quantile function returned by \code{\link{F_inv_from_f}}. Used
#' for \code{"SC"} and \code{"W"}. Computed by internally if \code{NULL}
#' (default).
#' @inheritParams F_inv_from_f
#' @details
#' The parameter \code{kappa} is used as \eqn{\kappa} in the following
#' distributions:
#' \itemize{
#'   \item \code{"vMF"}: von Mises--Fisher distribution with concentration
#'   \eqn{\kappa} and directional mean \eqn{{\bf e}_p = (0, 0, \ldots, 1)}{
#'   e_p = (0, 0, \ldots, 1)}.
#'   \item \code{"MvMF"}: equally-weighted mixture of \eqn{p} von Mises--Fisher
#'   distributions with common concentration \eqn{\kappa} and directional means
#'   \eqn{{\bf e}_1, \ldots, {\bf e}_p}{e_1, \ldots, e_p}.
#'   \item \code{"ACG"}: Angular Central Gaussian distribution with diagonal
#'   shape matrix with diagonal given by
#'   \deqn{(1, \ldots, 1, 1 + \kappa) / (p + \kappa).}
#'   \item \code{"SC"}: Small Circle distribution with axis mean
#'   \eqn{{\bf e}_p = (0, 0, \ldots, 1)}{e_p = (0, 0, \ldots, 1)} and
#'   concentration \eqn{\kappa} about the projection along the mean, \eqn{\nu}.
#'   \item \code{"W"}: Watson distribution with axis mean
#'   \eqn{{\bf e}_p = (0, 0, \ldots, 1)}{e_p = (0, 0, \ldots, 1)} and
#'   concentration \eqn{\kappa}. The Watson distribution is a particular case
#'   of the Bingham distribution.
#' }
#' @details
#' Much faster sampling for \code{"SC"} and \code{"W"} is achieved providing
#' \code{F_inv}, see examples.
#' @examples
#' ## Simulation with p = 2
#'
#' p <- 2
#' n <- 200
#' kappa <- 20
#' nu <- 0.5
#' F_inv_SC_2 <- F_inv_from_f(f = function(z) exp(-kappa * (z - nu)^2), p = 2)
#' F_inv_W_2 <- F_inv_from_f(f = function(z) exp(kappa * z^2), p = 2)
#' x1 <- r_alt(n = n, p = p, scenario = "vMF", kappa = kappa)[, , 1]
#' x2 <- r_alt(n = n, p = p, scenario = "MvMF", kappa = kappa)[, , 1]
#' x3 <- r_alt(n = n, p = p, scenario = "ACG", kappa = kappa)[, , 1]
#' x4 <- r_alt(n = n, p = p, scenario = "SC", F_inv = F_inv_SC_2)[, , 1]
#' x5 <- r_alt(n = n, p = p, scenario = "W", F_inv = F_inv_W_2)[, , 1]
#' r <- runif(n, 0.95, 1.05) # Radius perturbation to improve visualization
#' plot(r * x1, pch = 16, xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1), col = 1)
#' points(r * x2, pch = 16, col = 2)
#' points(r * x3, pch = 16, col = 3)
#' points(r * x4, pch = 16, col = 4)
#' points(r * x5, pch = 16, col = 5)
#'
#' ## Simulation with p = 3
#'
#' n <- 200
#' p <- 3
#' kappa <- 20
#' nu <- 0.5
#' F_inv_SC_3 <- F_inv_from_f(f = function(z) exp(-kappa * (z - nu)^2), p = 3)
#' F_inv_W_3 <- F_inv_from_f(f = function(z) exp(kappa * z^2), p = 3)
#' x1 <- r_alt(n = n, p = p, scenario = "vMF", kappa = kappa)[, , 1]
#' x2 <- r_alt(n = n, p = p, scenario = "MvMF", kappa = kappa)[, , 1]
#' x3 <- r_alt(n = n, p = p, scenario = "ACG", kappa = kappa)[, , 1]
#' x4 <- r_alt(n = n, p = p, scenario = "SC", F_inv = F_inv_SC_3)[, , 1]
#' x5 <- r_alt(n = n, p = p, scenario = "W", F_inv = F_inv_W_3)[, , 1]
#' if (requireNamespace("rgl"))
#'   rgl::plot3d(x1, size = 5, xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1),
#'               zlim = c(-1.1, 1.1), col = 1)
#'   rgl::points3d(x2, size = 5, col = 2)
#'   rgl::points3d(x3, size = 5, col = 3)
#'   rgl::points3d(x4, size = 5, col = 4)
#'   rgl::points3d(x5, size = 5, col = 5)
#' }
#' @export
r_alt <- function(n, p, M = 1, scenario = "vMF", kappa = 1, nu = 0.5,
                  F_inv = NULL, K = 1e3) {

  # Common mean (North pole)
  mu <- c(rep(0, p - 1), 1)

  # Choose scenario
  if (scenario == "vMF") {

    long_samp <- rotasym::r_vMF(n = n * M, mu = mu, kappa = kappa)

  } else if (scenario == "MvMF") {

    # Mixture components
    j <- sample(x = 1:p, size = n * M, replace = TRUE)
    nM_j <- table(j)
    mu_j <- diag(1, nrow = p, ncol = p)

    # Sample components
    long_samp <- matrix(nrow = n * M, ncol = p)
    for (k in 1:p) {

      long_samp[j == k, ] <- rotasym::r_vMF(n = nM_j[k], mu = mu_j[k, ],
                                            kappa = kappa)

    }

    # Add plus and minus means
    long_samp <- sample(x = c(-1, 1), size = n * M, replace = TRUE) * long_samp

    # Shuffle data
    long_samp <- long_samp[sample(x = n * M), ]

  } else if (scenario == "ACG") {

    Lambda <- diag(c(rep(1 / (p + kappa), p - 1),
                     (1 + kappa) / (p + kappa)), nrow = p, ncol = p)
    long_samp <- rotasym::r_ACG(n = n * M, Lambda = Lambda)

  } else if (scenario == "SC") {

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

  } else if (scenario == "W") {

    # Compute the inverse of the distribution function F?
    if (is.null(F_inv)) {

      f <- function(z) exp(kappa * z^2)
      F_inv <- F_inv_from_f(f = f, p = p, K = K)

    }

    # Sample the small circle distribution
    r_U <- function(n) r_unif_sph(n = n, p = p - 1, M = 1)[, , 1]
    r_V <- function(n) F_inv(runif(n = n))
    long_samp <- rotasym::r_tang_norm(n = n * M, theta = mu,
                                      r_U = r_U, r_V = r_V)

  } else {

    stop(paste("Wrong scenario; must be \"vMF\", \"MvMF\", \"Bing\"",
               "\"ACG\", \"SC\" or \"W\"."))

  }

  # As an array
  samp <- array(dim = c(n, p, M))
  for (j in 1:M) {

    samp[, , j] <- long_samp[(1 + (j - 1) * n):(j * n), , drop = FALSE]

  }
  return(samp)

}
