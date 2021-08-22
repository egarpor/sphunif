

#' @title Conversion between angular and Cartesian coordinates of the
#' (hyper)sphere
#'
#' @description Transforms the angles \eqn{(\theta_1,\ldots,\theta_{p-1})} in
#' \eqn{[0,\pi)^{p-2}\times[-\pi,\pi)} into the Cartesian coordinates
#' \deqn{(\cos(x_1),\sin(x_1)\cos(x_2),\sin(x_1)\sin(x_2)\cos(x_3),\ldots,
#' \sin(x_1)\cdots\sin(x_{p-2})\cos(x_{p-1}),
#' \sin(x_1)\cdots\sin(x_{p-2})\sin(x_{p-1}))}
#' of \eqn{S^{p-1}}, and viceversa.
#'
#' @param theta matrix of size \code{c(n, p - 1)} with the angles.
#' @param x matrix of size \code{c(n, p)} with the Cartesian coordinates.
#' Assumed to be of unit norm by rows.
#' @return For \code{angles_to_sphere}, the matrix \code{x}. For
#' \code{sphere_to_angles}, the matrix \code{theta}.
#' @examples
#' # Check changes of coordinates
#' sphere_to_angles(angles_to_sphere(c(pi / 2, 0, pi)))
#' sphere_to_angles(angles_to_sphere(rbind(c(pi / 2, 0, pi), c(pi, pi / 2, 0))))
#' angles_to_sphere(sphere_to_angles(c(0, sqrt(0.5), sqrt(0.1), sqrt(0.4))))
#' angles_to_sphere(sphere_to_angles(
#'   rbind(c(0, sqrt(0.5), sqrt(0.1), sqrt(0.4)),
#'         c(0, sqrt(0.5), sqrt(0.5), 0),
#'         c(0, 1, 0, 0),
#'         c(0, 0, 0, -1),
#'         c(0, 0, 1, 0))))
#'
#' # Circle
#' sphere_to_angles(angles_to_sphere(0))
#' sphere_to_angles(angles_to_sphere(cbind(0:3)))
#' angles_to_sphere(cbind(sphere_to_angles(rbind(c(0, 1), c(1, 0)))))
#' angles_to_sphere(cbind(sphere_to_angles(rbind(c(0, 1)))))
#' @export
angles_to_sphere <- function(theta) {

  # Convert to matrix
  if (!is.matrix(theta)) {

    theta <- rbind(theta)

  }

  # Apply transformation
  q <- ncol(theta)
  x <- matrix(0, nrow = nrow(theta), ncol = q + 1)
  cos_theta <- cos(theta)
  sin_theta <- sin(theta)
  if (q == 1) {

    x <- cbind(cos_theta, sin_theta)

  } else {

    prod_sin_theta <- t(apply(sin_theta, 1, cumprod))
    x <- cbind(cos_theta[, 1, drop = FALSE],
               prod_sin_theta[, -q, drop = FALSE] * 
                 cos_theta[, -1, drop = FALSE],
               prod_sin_theta[, q, drop = FALSE])

  }
  return(unname(x))

}


#' @rdname angles_to_sphere
#' @export
sphere_to_angles <- function(x) {

  # Convert to matrix
  if (!is.matrix(x)) {

    x <- rbind(x)

  }

  # Check unit norm
  x <- rotasym::check_unit_norm(x)

  # Apply transformation
  p <- ncol(x)
  i_norm <- t(apply(x^2, 1, function(x) rev(sqrt(cumsum(rev(x))))))[, -p]
  theta <- x[, -p, drop = FALSE] / i_norm
  theta[is.nan(theta)] <- 1 # Avoid NaNs
  theta <- acos(theta)
  xp <- (x[, p] < 0)
  theta[, p - 1] <- (2 * pi) * xp + (1 - 2 * xp) * theta[, p - 1]
  return(unname(theta))

}


#' @title (Hyper)spherical harmonics
#'
#' @description Computation of a certain explicit representation of
#' (hyper)spherical harmonics on
#' \eqn{S^{p-1}:=\{{\bf x}\in R^p:||{\bf x}||=1\}}{
#' S^{p-1}:=\{x\in R^p:||x||=1\}}, \eqn{p\ge 2}. Details are available in
#' García-Portugués et al. (2021).
#'
#' @param x locations in \eqn{S^{p-1}} to evaluate \eqn{g_{i,k}}. Either a
#' matrix of size \code{c(nx, p)} or a vector of length \code{p}. Normalized
#' internally if required (with a \code{warning} message).
#' @param i,k alternative indexing to refer to the \code{i}-th (hyper)spherical
#' harmonic of order \code{k}. \code{i} is a positive integer smaller than
#' \code{\link[=Sobolev]{d_p_k}} and \code{k} is a non-negative integer.
#'
#' @param m (hyper)spherical harmonic index, as used in Proposition 2.1. The
#' index is computed internally from \code{i} and \code{k}. Defaults to
#' \code{NULL}.
#' @param show_m flag to print \code{m} if computed internally when
#' \code{m = NULL}.
#' @return A vector of size \code{nrow(x)}.
#' @details
#' The implementation uses Proposition 2.1 in García-Portugués et al. (2021),
#' which adapts Theorem 1.5.1 in Dai and Xu (2013) with the correction of
#' typos in the normalizing constant \eqn{h_\alpha} and in the definition of
#' the function \eqn{g_\alpha} of the latter theorem.
#' @references 
#' Dai, F. and Xu, Y. (2013). \emph{Approximation Theory and Harmonic Analysis
#' on Spheres and Balls}. Springer, New York. \doi{10.1007/978-1-4614-6660-4}
#' 
#' García-Portugués, E., Paindaveine, D., and Verdebout, T. (2021). On the
#' power of Sobolev tests for isotropy under local rotationally symmetric
#' alternatives. \emph{arXiv:2108.XXXXX}. \url{https://arxiv.org/abs/2108.XXXXX}
#' @examples
#' n <- 5e3
#' old_par <- par(mfrow = c(2, 3))
#' k <- 2
#' for (i in 1:d_p_k(p = 3, k = k)) {
#'   X <- r_unif_sph(n = n, p = 3, M = 1)[, , 1]
#'   col <- rainbow(n)[rank(g_i_k(x = X, k = k, i = i, show_m = TRUE))]
#'   scatterplot3d::scatterplot3d(X[, 1], X[, 2], X[, 3], color = col,
#'                                axis = FALSE)
#' }
#' for (k in 0:5) {
#'   X <- r_unif_sph(n = n, p = 3, M = 1)[, , 1]
#'   col <- rainbow(n)[rank(g_i_k(x = X, k = k, i = 1, show_m = TRUE))]
#'   scatterplot3d::scatterplot3d(X[, 1], X[, 2], X[, 3], color = col,
#'                                axis = FALSE)
#' }
#' par(old_par)
#' @export
#' @name harmonics


#' @rdname harmonics
g_i_k <- function(x, i = 1, k = 1, m = NULL, show_m = FALSE) {

  # Set x as a matrix
  if (!is.matrix(x)) {

    x <- rbind(x)

  }

  # Check unit norm
  x <- rotasym::check_unit_norm(x)

  # Get dimension
  p <- ncol(x)

  # Obtain m if not supplied
  if (is.null(m)) {

    # Initial check
    dpk <- d_p_k(p = p, k = k)
    stopifnot(1 <= i && i <= dpk + 0.1)

    # Possible m's for k
    ms <- as.matrix(do.call(expand.grid, c(rep(list(0:k), p - 1), list(0:1))))
    ms <- ms[rowSums(ms) == k, , drop = FALSE]
    stopifnot(abs(dpk - nrow(ms)) < 0.1)

    # Selected m
    m <- ms[i, ]
    if (show_m) message(paste0("m = (", paste0(m, collapse = ", "), ")"))

  } else {

    # Retrieve k for p = 2
    k <- sum(m)

    # Initial checks
    stopifnot(p == length(m))
    stopifnot(all(m >= 0) && m[p] %in% c(0, 1))

  }

  # Obtain angles
  theta <- sphere_to_angles(x = x)[, (p - 1):1, drop = FALSE]

  # Treat differently the p = 2 and p >= 3 cases
  if (p == 2) {

    if (m[2] == 0) {

      g <- cos(k * theta[, 1])

    } else if (m[2] == 1) {

      g <- sin(k * theta[, 1])

    }
    g <- sqrt(2) * g

  } else if (p >= 3) {

    # zeta_m
    if (m[p] == 0) {

      zeta_m <- cos(m[p - 1] * theta[, 1])

    } else if (m[p] == 1) {

      zeta_m <- sin((m[p - 1] + 1) * theta[, 1])

    }

    # b_m
    b_m <- ifelse(m[p - 1] + m[p] > 0, 2, 1)

    # Special notation
    abs_j <- function(m, j) sum(m[j:length(m)])
    lpoch <- function(a, n) lgamma(a + n) - lgamma(a)
    lamb_j <- function(m, j, p) abs_j(m = m, j = j + 1) + (p - j - 1) / 2

    # B_m
    B_m <- 0
    for (j in 1:(p - 2)) {

      la_j <- lamb_j(m = m, j = j, p = p)
      m_j1 <- abs_j(m = m, j = j + 1)
      B_m <- B_m + lfactorial(m[j]) + lpoch(a = (p - j + 1) / 2, n = m_j1) +
        log(m[j] + la_j) - lpoch(a = 2 * la_j, n = m[j]) -
        lpoch(a = (p - j) / 2, n = m_j1) - log(la_j)

    }
    B_m <- exp(log(b_m) + B_m)

    # g
    g <- 1
    for (j in 1:(p - 2)) {

      g <- g * sin(theta[, p - j])^abs_j(m = m, j = j + 1) *
        gsl::gegenpoly_n(n = m[j], lambda = lamb_j(m = m, j = j, p = p),
                         x = cos(theta[, p - j]))

    }
    g <- g * sqrt(B_m) * zeta_m

  }

  return(unname(g))

}
