

#' @title Gegenbauer polynomials and coefficients
#'
#' @description The \href{https://dlmf.nist.gov/18.3}{Gegenbauer polynomials}
#' \eqn{\{C_k^{(\lambda)}(x)\}_{k = 0}^\infty}{
#' {C_k^(\lambda)(x)}_{k = 0}^\infty}
#' form a family of orthogonal polynomials on the interval \eqn{[-1, 1]}
#' with respect to the weight function \eqn{(1 - x^2)^{\lambda - 1/2}},
#' for \eqn{\lambda > -1/2}, \eqn{\lambda \neq 0}. They usually appear
#' when dealing with functions defined on
#' \eqn{S^{p-1} := \{{\bf x} \in R^p : ||{\bf x}|| = 1\}}{
#' S^{p-1} := \{x \in R^p : ||x|| = 1\}} with index \eqn{\lambda = p / 2 - 1}.
#'
#' The Gegenbauer polynomials are somehow simpler to evaluate for
#' \eqn{x = \cos(\theta)}, with \eqn{\theta \in [0, \pi]}. This simplifies
#' also the connection with the Chebyshev polynomials
#' \eqn{\{T_k(x)\}_{k = 0}^\infty}{{T_k(x)}_{k = 0}^\infty}, which admit
#' the \href{https://dlmf.nist.gov/18.5.E1}{explicit expression}
#' \eqn{T_k(\cos(\theta)) = \cos(k\theta)}. The Chebyshev polynomials
#' appear as the limit of the Gegenbauer polynomials
#' (divided by \eqn{\lambda}) when \eqn{\lambda} goes to \eqn{0}, so they
#' can be regarded as the extension by continuity of
#' \eqn{\{C_k^{(p/2 - 1)}(x)\}_{k = 0}^\infty}{
#' {C_k^(p/2 - 1)(x)}_{k = 0}^\infty} to the case \eqn{p = 2}.
#'
#' For a \href{https://dlmf.nist.gov/18.18.i}{reasonably smooth} function
#' \eqn{\psi} defined on \eqn{[0, \pi]},
#' \deqn{\psi(\theta) = \sum_{k = 0}^\infty b_{k, p}
#' C_k^{(p/2 - 1)}(\cos(\theta)),}{\psi(\theta) = \sum_{k = 0}^\infty b_{k, p}
#' C_k^(p/2 - 1)(\cos(\theta)),}
#' provided that the coefficients
#' \deqn{b_{k, p} := \frac{1}{c_{k, p}} \int_0^\pi \psi(\theta)
#' C_k^{(p/2 - 1)}(\cos(\theta)) (\sin(\theta))^{p - 2}\,\mathrm{d}\theta}{
#' b_{k, p} := \frac{1}{c_{k, p}} \int_0^\pi \psi(\theta)
#' C_k^(p/2 - 1)(\cos(\theta)) (\sin(\theta))^{p - 2} d\theta}
#' are finite, where the normalizing constants are
#' \deqn{c_{k, p} := \int_0^\pi (C_k^{(p/2 - 1)}(\cos(\theta)))^2
#' (\sin(\theta))^{p - 2} \,\mathrm{d}\theta.}{
#' c_{k, p} := \int_0^\pi (C_k^(p/2 - 1)(\cos(\theta)))^2
#' (\sin(\theta))^{p - 2} d\theta.}
#' The (squared) "Gegenbauer norm" of \eqn{\psi} is
#' \deqn{\|\psi\|_{G, p}^2 := \int_0^\pi \psi(\theta)^2
#' C_k^{(p/2 - 1)}(\cos(\theta)) (\sin(\theta))^{p - 2}\,\mathrm{d}\theta.}{
#' ||\psi||_{G, p}^2 := \int_0^\pi \psi(\theta)^2
#' C_k^(p/2 - 1)(\cos(\theta)) (\sin(\theta))^{p - 2} d\theta.}
#'
#' The previous expansion can be generalized for a 2-dimensional function
#' \eqn{\psi} defined on \eqn{[0, \pi] \times [0, \pi]}:
#' \deqn{\psi(\theta_1, \theta_2) = \sum_{k = 0}^\infty \sum_{m = 0}^\infty
#' b_{k, m, p} C_k^{(p/2 - 1)}(\cos(\theta_1))
#' C_k^{(p/2 - 1)}(\cos(\theta_2)),}{
#' \psi(\theta_1, \theta_2) = \sum_{k = 0}^\infty \sum_{m = 0}^\infty
#' b_{k, m, p} C_k^(p/2 - 1)(\cos(\theta_1)) C_k^(p/2 - 1)(\cos(\theta_2))}
#' with coefficients
#' \deqn{b_{k, m, p} := \frac{1}{c_{k, p} c_{m, p}} \int_0^\pi\int_0^\pi
#' \psi(\theta_1, \theta_2) C_k^{(p/2 - 1)}(\cos(\theta_1))
#' C_k^{(p/2 - 1)}(\cos(\theta_2)) (\sin(\theta_1))^{p - 2}
#' (\sin(\theta_2))^{p - 2}\,\mathrm{d}\theta_1\,\mathrm{d}\theta_2.}{
#' b_{k, m, p} := \frac{1}{c_{k, p} c_{m, p}} \int_0^\pi \int_0^\pi
#' \psi(\theta_1, \theta_2) C_k^(p/2 - 1)(\cos(\theta_1))
#' C_k^(p/2 - 1)(\cos(\theta_2)) (\sin(\theta_1))^{p - 2}
#' (\sin(\theta_2))^{p - 2} d\theta_1 d\theta_2.}
#' The (squared) "Gegenbauer norm" of \eqn{\psi} is
#' \deqn{\|\psi\|_{G, p}^2 := \int_0^\pi\int_0^\pi \psi(\theta_1, \theta_2)^2
#' C_k^{(p/2 - 1)}(\cos(\theta_1)) C_k^{(p/2 - 1)}(\cos(\theta_2))
#' (\sin(\theta_1))^{p - 2} (\sin(\theta_2))^{p - 2}
#' \,\mathrm{d}\theta_1\,\mathrm{d}\theta_2.}{
#' ||\psi||_{G, p}^2 := \int_0^\pi\int_0^\pi \psi(\theta_1, \theta_2)^2
#' C_k^(p/2 - 1)(\cos(\theta_1)) C_k^(p/2 - 1)(\cos(\theta_2))
#' (\sin(\theta_1))^{p - 2} (\sin(\theta_2))^{p - 2} d\theta_1 d\theta_2.}
#'
#' @param theta,theta_1,theta_2 vectors with values in \eqn{[0, \pi]}.
#' @param k,m vectors with the orders of the Gegenbauer polynomials. Must
#' be integers larger or equal than \code{0}.
#' @inheritParams r_unif
#' @param psi function defined in \eqn{[0, \pi]} and whose Gegenbauer
#' coefficients are to be computed. Must be vectorized. For
#' \code{Gegen_coefs_2d}, it must return a matrix of size
#' \code{c(length(theta_1), length(theta_2))}.
#' @param coefs for \code{Gegen_series} and \code{Gegen_norm}, a vector of
#' coefficients \eqn{b_{k, p}} with length \code{length(k)}. For
#' \code{Gegen_series_2d} and \code{Gegen_norm_2d}, a matrix of coefficients
#' \eqn{b_{k, m, p}} with size \code{c(length(k), length(m))}. The
#' order of the coefficients is given by \code{k} and \code{m}.
#' @param Gauss use a Gauss--Legendre quadrature rule of \code{N} nodes
#' in the computation of the Gegenbauer coefficients? Otherwise, call
#' \code{\link{integrate}}. Defaults to \code{TRUE}.
#' @param N number of points used in the \link[=Gauss_Legen_nodes]{
#' Gauss--Legendre quadrature} for computing the Gegenbauer coefficients.
#' Defaults to \code{320}.
#' @param normalize consider normalized coefficients (divided by
#' \eqn{c_{k, p}})? Defaults to \code{TRUE}.
#' @param only_const return only the normalizing constants \eqn{c_{k, p}}?
#' Defaults to \code{FALSE}.
#' @param tol tolerance passed to \code{\link{integrate}}'s \code{rel.tol} and
#' \code{abs.tol} if \code{Gauss = FALSE}. Defaults to \code{1e-6}.
#' @param ... further arguments to be passed to \code{psi}.
#' @param cumulative return the cumulative norm for increasing truncation of
#' the series? Defaults to \code{FALSE}.
#' @return
#' \itemize{
#'   \item \code{Gegen_polyn}: a matrix of size
#'   \code{c(length(theta), length(k))} containing the evaluation of the
#'   \code{length(k)} Gegenbauer polynomials at \code{theta}.
#'   \item \code{Gegen_coefs}: a vector of size \code{length(k)} containing
#'   the coefficients \eqn{b_{k, p}}.
#'   \item \code{Gegen_series}: the evaluation of the truncated series
#'   expansion, a vector of size \code{length(theta)}.
#'   \item \code{Gegen_norm}: the Gegenbauer norm of the truncated series,
#'   a scalar if \code{cumulative = FALSE}, otherwise a vector of size
#'   \code{length(k)}.
#'   \item \code{Gegen_polyn_2d}: a 4-dimensional array of size
#'   \code{c(length(theta_1), length(theta_2), length(k), length(m))}
#'   containing the evaluation of the \code{length(k) * length(m)}
#'   2-dimensional Gegenbauer polynomials at the bivariate grid
#'   spanned by \code{theta_1} and \code{theta_2}.
#'   \item \code{Gegen_coefs_2d}: a matrix of size
#'   \code{c(length(k), length(m))} containing the coefficients
#'   \eqn{b_{k, m, p}}.
#'   \item \code{Gegen_series_2d}: the evaluation of the truncated series
#'   expansion, a matrix of size \code{c(length(theta_1), length(theta_2))}.
#'   \item \code{Gegen_norm_2d}: the 2-dimensional Gegenbauer norm of the
#'   truncated series, a scalar.
#' }
#' @details
#' The \code{Gegen_polyn} function is a wrapper to the functions
#' \link[gsl:Gegenbauer]{gegenpoly_n} and \link[gsl:Gegenbauer]{gegenpoly_array}
#' in the \link[gsl]{gsl-package}, which they interface the functions
#' defined in the header file \code{gsl_sf_gegenbauer.h} (documented
#' \href{https://www.gnu.org/software/gsl/doc/html/specfunc.html#gegenbauer-functions}{
#' here}) of the \href{https://www.gnu.org/software/gsl/}{
#' GNU Scientific Library}.
#'
#' Note that the function \code{Gegen_polyn} computes the regular
#' \emph{unnormalized} Gegenbauer polynomials.
#'
#' For the case \eqn{p = 2}, the Chebyshev polynomials are considered.
#' @references
#' Galassi, M., Davies, J., Theiler, J., Gough, B., Jungman, G., Alken, P.,
#' Booth, M., and Rossi, F. (2009) \emph{GNU Scientific Library Reference
#' Manual}. Network Theory Ltd. \url{http://www.gnu.org/software/gsl/}
#'
#' \emph{NIST Digital Library of Mathematical Functions}. Release
#' 1.0.20 of 2018-09-15. F. W. J. Olver, A. B. Olde Daalhuis, D. W. Lozier,
#' B. I. Schneider, R. F. Boisvert, C. W. Clark, B. R. Miller,
#' and B. V. Saunders, eds. \url{https://dlmf.nist.gov/}
#' @examples
#' ## Representation of Gegenbauer polynomials (Chebyshev polynomials for p = 2)
#'
#' th <- seq(0, pi, l = 500)
#' k <- 0:3
#' old_par <- par(mfrow = c(2, 2))
#' for (p in 2:5) {
#'   matplot(th, t(Gegen_polyn(theta = th, k = k, p = p)), lty = 1,
#'           type = "l", main = substitute(p == d, list(d = p)),
#'           axes = FALSE, xlab = expression(theta), ylab = "")
#'   axis(1, at = c(0, pi / 4, pi / 2, 3 * pi / 4, pi),
#'        labels = expression(0, pi / 4, pi / 2, 3 * pi / 4, pi))
#'   axis(2); box()
#'   mtext(text = expression({C[k]^{p/2 - 1}}(cos(theta))), side = 2,
#'         line = 2, cex = 0.75)
#'   legend("bottomleft", legend = paste("k =", k), lwd = 2, col = seq_along(k))
#' }
#' par(old_par)
#'
#' ## Coefficients and series in p = 2
#'
#' # Function in [0, pi] to be projected in Chebyshev polynomials
#' psi <- function(th) -sin(th / 2)
#'
#' # Coefficients
#' p <- 2
#' k <- 0:4
#' (coefs <- Gegen_coefs(k = k, p = p, psi = psi))
#'
#' # Series
#' plot(th, psi(th), type = "l", axes = FALSE, xlab = expression(theta),
#'       ylab = "", ylim = c(-1.25, 0))
#' axis(1, at = c(0, pi / 4, pi / 2, 3 * pi / 4, pi),
#'      labels = expression(0, pi / 4, pi / 2, 3 * pi / 4, pi))
#' axis(2); box()
#' col <- viridisLite::viridis(length(coefs))
#' for (i in seq_along(coefs)) {
#'   lines(th, Gegen_series(theta = th, coefs = coefs[1:(i + 1)], k = 0:i,
#'                          p = p), col = col[i])
#' }
#' lines(th, psi(th), lwd = 2)
#'
#' ## Coefficients and series in p = 3
#'
#' # Function in [0, pi] to be projected in Gegenbauer polynomials
#' psi <- function(th) tan(th / 3)
#'
#' # Coefficients
#' p <- 3
#' k <- 0:10
#' (coefs <- Gegen_coefs(k = k, p = p, psi = psi))
#'
#' # Series
#' plot(th, psi(th), type = "l", axes = FALSE, xlab = expression(theta),
#'       ylab = "", ylim = c(0, 2))
#' axis(1, at = c(0, pi / 4, pi / 2, 3 * pi / 4, pi),
#'      labels = expression(0, pi / 4, pi / 2, 3 * pi / 4, pi))
#' axis(2); box()
#' col <- viridisLite::viridis(length(coefs))
#' for (i in seq_along(coefs)) {
#'   lines(th, Gegen_series(theta = th, coefs = coefs[1:(i + 1)], k = 0:i,
#'                          p = p), col = col[i])
#' }
#' lines(th, psi(th), lwd = 2)
#'
#' ## Surface representation
#'
#' # Surface in [0, pi]^2 to be projected in Gegenbauer polynomials
#' p <- 3
#' psi <- function(th_1, th_2) A_theta_x(theta = th_1, x = cos(th_2),
#'                                       p = p, as_matrix = TRUE)
#'
#' # Coefficients
#' k <- 0:20
#' m <- 0:10
#' coefs <- Gegen_coefs_2d(k = k, m = m, p = p, psi = psi)
#'
#' # Series
#' th <- seq(0, pi, l = 100)
#' col <- viridisLite::viridis(20)
#' old_par <- par(mfrow = c(2, 2))
#' image(th, th, A_theta_x(theta = th, x = cos(th), p = p), axes = FALSE,
#'       col = col, zlim = c(0, 1), xlab = expression(theta[1]),
#'       ylab = expression(theta[2]), main = "Original")
#' axis(1, at = c(0, pi / 4, pi / 2, 3 * pi / 4, pi),
#'      labels = expression(0, pi / 4, pi / 2, 3 * pi / 4, pi))
#' axis(2, at = c(0, pi / 4, pi / 2, 3 * pi / 4, pi),
#'      labels = expression(0, pi / 4, pi / 2, 3 * pi / 4, pi))
#' box()
#' for(K in c(5, 10, 20)) {
#'   A <- Gegen_series_2d(theta_1 = th, theta_2 = th,
#'                        coefs = coefs[1:(K + 1), ], k = 0:K, m = m, p = p)
#'   image(th, th, A, axes = FALSE, col = col, zlim = c(0, 1),
#'         xlab = expression(theta[1]), ylab = expression(theta[2]),
#'         main = paste(K, "x", m[length(m)], "coefficients"))
#'   axis(1, at = c(0, pi / 4, pi / 2, 3 * pi / 4, pi),
#'        labels = expression(0, pi / 4, pi / 2, 3 * pi / 4, pi))
#'   axis(2, at = c(0, pi / 4, pi / 2, 3 * pi / 4, pi),
#'        labels = expression(0, pi / 4, pi / 2, 3 * pi / 4, pi))
#'   box()
#' }
#' par(old_par)
#' @name Gegenbauer


#' @rdname Gegenbauer
#' @export
Gegen_polyn <- function(theta, k, p) {

  # Check orders
  stopifnot(all(k >= 0))

  # Check dimension
  stopifnot(p >= 2)

  # Check that theta is in range
  stopifnot(all(0 <= theta & theta <= pi))

  # Handle the special case p = 2 (Chebyshev series) separately
  p <- as.integer(p)
  if (p == 2) {

    return(cos(k %o% theta))

  } else {

    # Vectorization on k or not
    if (length(k) == 1) {

      return(rbind(gsl::gegenpoly_n(n = k, lambda = 0.5 * p - 1,
                                    x = cos(theta))))

    } else {

      return(gsl::gegenpoly_array(nmax = max(k), lambda = 0.5 * p - 1,
                                  x = cos(theta))[k + 1, , drop = FALSE])

    }

  }

}


#' @rdname Gegenbauer
#' @export
Gegen_coefs <- function(k, p, psi, Gauss = TRUE, N = 320, normalize = TRUE,
                        only_const = FALSE, tol = 1e-6, ...) {

  # Check orders
  stopifnot(all(k >= 0))

  # Check dimension
  stopifnot(p >= 2)

  # Log-normalizing constants for each projection
  if (normalize) {

    if (p == 2) {

      # Beware if p = 2: the normalizing constant for the zeroth-order Chebyshev
      # term (constant term) is pi / 2, whereas for the rest is pi
      # Check table https://dlmf.nist.gov/18.3 in NIST
      log_const <- log(ifelse(k > 0, pi / 2, pi))

    } else {

      # Normalizing constants worked out from table https://dlmf.nist.gov/18.3
      # in NIST (with lambda = p / 2 - 1)
      log_const <- (4 - p) * log(2) + log(pi) + lgamma(p + k - 2) -
        (log(p + 2 * (k - 1)) + lfactorial(k) + 2 * lgamma(p / 2 - 1))

    }

  } else {

    log_const <- 0

  }

  # Return only the normalizing constants?
  if (only_const) {

    return(exp(log_const))

  }

  # Compute projections
  if (Gauss) {

    # Gauss--Legendre nodes and weigths
    th_k <- drop(Gauss_Legen_nodes(a = 0, b = pi, N = N))
    w_k <- drop(Gauss_Legen_weights(a = 0, b = pi, N = N))

    # Integral approximation, using implicit column recycling for
    # multiplication by vectors
    projs <- colSums((w_k * psi(th_k, ...) * sin(th_k)^(p - 2)) *
                       t(Gegen_polyn(theta = th_k, k = k, p = p)),
                     na.rm = TRUE)

  } else {

    # Use intergrate()
    projs <- sapply(k, function(k_i) {
      f <- function(th) psi(th, ...) * sin(th)^(p - 2) *
        drop(Gegen_polyn(theta = th, k = k_i, p = p))
      integrate(f = f, lower = 0, upper = pi, subdivisions = 1e3,
                rel.tol = tol, abs.tol = tol,
                stop.on.error = FALSE)$value
    })

  }

  # Normalization of the projections
  return(sign(projs) * exp(log(abs(projs)) - log_const))

}


#' @rdname Gegenbauer
#' @export
Gegen_series <- function(theta, coefs, k, p, normalize = TRUE) {

  # Check orders
  stopifnot(all(k >= 0))

  # Check dimension
  stopifnot(p >= 2)

  # Check leength of coefficients
  stopifnot(length(coefs) == length(k))

  # Normalize unnormalized coefficients?
  if (!normalize) {

    coefs <- coefs / Gegen_coefs(k = k, p = p, only_const = TRUE)

  }

  # Series
  return(drop(coefs %*% Gegen_polyn(theta = theta, k = k, p = p)))

}


#' @rdname Gegenbauer
#' @export
Gegen_norm <- function(coefs, k, p, normalize = TRUE, cumulative = FALSE) {

  # Normalizing constants (required for both cases)
  c_kp <- Gegen_coefs(k = k, p = p, only_const = TRUE)

  # Norm
  op <- ifelse(cumulative, cumsum, sum)
  if (normalize) {

    return(sqrt(op(coefs^2 * c_kp)))

  } else {

    return(sqrt(op(coefs^2 / c_kp)))

  }

}


#' @rdname Gegenbauer
#' @export
Gegen_polyn_2d <- function(theta_1, theta_2, k, m, p) {

  # Check orders
  stopifnot(all(k >= 0))
  stopifnot(all(m >= 0))

  # Check dimension
  stopifnot(p >= 2)

  # Check that theta_1 and theta_2 are in range
  stopifnot(all(0 <= theta_1 & theta_1 <= pi))
  stopifnot(all(0 <= theta_2 & theta_2 <= pi))

  # Compute a 4-dimensional array of size length(theta_1) x length(theta_2) x
  # length(k) x length(m) with the k x m surfaces of the 2-dimensional
  # Gegenbauer polynomials
  polyn_2d <- aperm(a = Gegen_polyn(theta = theta_1, k = k, p = p) %o%
                      Gegen_polyn(theta = theta_2, k = m, p = p),
                    perm = c(2, 4, 1, 3))
  return(polyn_2d)

}


#' @rdname Gegenbauer
#' @export
Gegen_coefs_2d <- function(k, m, p, psi, Gauss = TRUE, N = 320,
                           normalize = TRUE, only_const = FALSE, tol = 1e-6,
                           ...) {

  # Check orders
  stopifnot(all(k >= 0))
  stopifnot(all(m >= 0))

  # Check dimension
  stopifnot(p >= 2)

  # Log-normalizing constants for each projection
  if (normalize) {

    if (p == 2) {

      # Beware if p = 2: the normalizing constant for the zeroth-order Chebyshev
      # term (constant term) is pi / 2, whereas for the rest is pi
      # Check table https://dlmf.nist.gov/18.3 in NIST
      log_const_k <- log(ifelse(k > 0, pi / 2, pi))
      log_const_m <- log(ifelse(m > 0, pi / 2, pi))

    } else {

      # Normalizing constants worked out from table https://dlmf.nist.gov/18.3
      # in NIST (with lambda = p / 2 - 1)
      log_const_k <- (4 - p) * log(2) + log(pi) + lgamma(p + k - 2) -
        (log(p + 2 * (k - 1)) + lfactorial(k) + 2 * lgamma(p / 2 - 1))
      log_const_m <- (4 - p) * log(2) + log(pi) + lgamma(p + m - 2) -
        (log(p + 2 * (m - 1)) + lfactorial(m) + 2 * lgamma(p / 2 - 1))

    }

    # Products of univariate normalizing constants
    log_const <- outer(log_const_k, log_const_m, `+`)

  } else {

    log_const <- matrix(0, nrow = length(k), ncol = length(m))

  }

  # Return only the normalizing constants?
  if (only_const) {

    return(exp(log_const))

  }

  # Compute projections
  if (Gauss) {

    # Gauss--Legendre nodes and weigths
    th_k <- drop(Gauss_Legen_nodes(a = 0, b = pi, N = N))
    w_k <- drop(Gauss_Legen_weights(a = 0, b = pi, N = N))

    # Compute polynomials
    polyn_2d <- Gegen_polyn_2d(theta_1 = th_k, theta_2 = th_k, k = k, m = m,
                               p = p)

    # Jacobian
    sin_q <- tcrossprod(sin(th_k)^(p - 2))

    # Expand weights
    w_k <- tcrossprod(w_k)

    # psi evaluated on the bivariate grid
    psi_q <- matrix(psi(th_k, th_k, ...), nrow = N, ncol = N)

    # Integral approximation
    projs <- apply(polyn_2d, c(3, 4),
                   function(C_km) sum(w_k * psi_q * C_km * sin_q, na.rm = TRUE))

  } else {

    # Use intergrate()
    projs <- sapply(m, function(m_i) {
      sapply(k, function(k_i) {

        # Integrand on th_1
        f_1 <- function(th_1) {
          sapply(th_1, function(th) {

            # Integrand on th_2
            f_2 <- function(th_2) {
              drop(psi(th, th_2, ...)) * sin(th)^(p - 2) * sin(th_2)^(p - 2) *
                drop(Gegen_polyn_2d(theta_1 = th, theta_2 = th_2, k = k_i,
                                    m = m_i, p = p))
            }

            # Integral on th_2
            integrate(f = f_2, lower = 0, upper = pi, subdivisions = 1e3,
                      rel.tol = tol, abs.tol = tol, stop.on.error = FALSE)$value

          })
        }

        # Integral on th_1
        integrate(f = f_1, lower = 0, upper = pi, subdivisions = 1e3,
                  rel.tol = tol, abs.tol = tol, stop.on.error = FALSE)$value
      })
    })

  }

  # Normalization of the projections
  return(sign(projs) * exp(log(abs(projs)) - log_const))

}


#' @rdname Gegenbauer
#' @export
Gegen_series_2d <- function(theta_1, theta_2, coefs, k, m, p,
                            normalize = TRUE) {

  # Check orders
  stopifnot(all(k >= 0))
  stopifnot(all(m >= 0))

  # Check dimension
  stopifnot(p >= 2)

  # Check size of coefficients
  stopifnot(nrow(coefs) == length(k) & ncol(coefs) == length(m))

  # Normalize unnormalized coefficients?
  if (!normalize) {

    coefs <- coefs / Gegen_coefs_2d(k = k, m = m, p = p, only_const = TRUE)

  }

  # Compute polynomials
  polyn_2d <- Gegen_polyn_2d(theta_1 = theta_1, theta_2 = theta_2,
                             k = k, m = m, p = p)

  # Add them into the series
  return(apply(apply(polyn_2d, c(1, 2),
                     function(C_km) c(coefs * C_km)), 2:3, sum))

}


#' @rdname Gegenbauer
#' @export
Gegen_norm_2d <- function(coefs, k, m, p, normalize = TRUE) {

  # Normalizing constants (required for both cases)
  c_kmp <- Gegen_coefs_2d(k = k, m = m, p = p, only_const = TRUE)

  # Norm
  if (normalize) {

    return(sqrt(sum(coefs^2 * c_kmp)))

  } else {

    return(sqrt(sum(coefs^2 / c_kmp)))

  }

}
