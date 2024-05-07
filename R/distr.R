

#' @rdname proj_unif
#' @export
r_proj_unif <- function(n, p) {

  # Check p
  if (p <= 1) {

    stop("p must be >= 2.")

  }

  # The projection is the square of a B(1 / 2, (p - 1)/2)
  return(sqrt(rbeta(n = n, shape1 = 0.5, shape2 = (p - 1) / 2)) *
           sample(c(-1, 1), size = n, replace = TRUE))

}


#' @title Uniform cap distribution
#'
#' @description Density and simulation of a uniform distribution within a
#' spherical cap of radius \eqn{0\leq r\leq \pi} around a direction
#' \eqn{\boldsymbol{\mu}}{\mu} on \eqn{S^{p-1}:=\{\mathbf{x}\in
#' R^p:||\mathbf{x}||=1\}}{S^{p-1}:=\{x\in R^p:||x||=1\}}, \eqn{p > 1}.
#' The density at \eqn{\mathbf{x} \in S^{p-1}}{x \in S^{p-1}} is given by
#' \deqn{c_{p,r}1_{[1-\tilde{r}, 1]}(\mathbf{x}' \boldsymbol{\mu})
#' \quad\mathrm{with}\quad c_{p,r}:=\omega_{p}\left[1-F_{p}(1-\tilde{r})
#' \right]}{c_{p,r}1_{[1-\tilde{r}, 1]}(x' \mu) with \quad c_{p,r}:=\omega_{p}
#' \left[1-F_{p}(1-\tilde{r})\right]} where \eqn{\boldsymbol{\mu}\in S^{p-1}}{
#' \mu\inS^{p-1}} is the directional mean, \eqn{\tilde{r}=\cos(r)} is the
#' projected radius of the spherical cap about \eqn{\boldsymbol{\mu}}{\mu},
#' \eqn{\omega_{p}} is the surface area of \eqn{S^{p-1}}{S^{p-1}}, and
#' \eqn{F_p(x)} is the cdf of a \eqn{B(\frac{1}{2},\frac{p-1}{2})}
#' beta distribution.
#'
#' The angular function of the uniform cap distribution is
#' \eqn{g(t):=1_{[1-\tilde{r}, 1]}(t)}{g(t):=I_{[1-\tilde{r}, 1]}(t)}. The
#' associated projected density is \eqn{\tilde{g}(v):=\omega_{p-1}c_{p,r}
#' (1-v^2)^{\frac{p-3}{2}}1_{[1-\tilde{r}, 1]}(v)}{\tilde{g}(v):=\omega_{p-1}
#' c_{p,r}(1-v^2)^{\frac{p-3}{2}}I_{[1-\tilde{r}, 1]}(v)}.
#'
#' @inheritParams rotasym::d_tang_norm
#' @inheritParams r_unif
#' @inheritParams proj_unif
#' @param x locations to evaluate the density or distribution. For
#' \code{d_unif_cap}, positions on \eqn{S^{p - 1}} given either as a matrix of
#' size \code{c(nx, p)} or a vector of length \code{p}. Normalized internally if
#' required (with a \code{warning} message). For \code{d_unif_proj_cap} and
#' \code{p_unif_proj_cap}, a vector with values in \eqn{[-1, 1]}.
#' @param mu the directional mean \eqn{\boldsymbol{\mu}}{\mu} of the
#' distribution. A unit-norm vector of length \code{p}.
#' @param r angle of the spherical cap around the antipodal observation in
#' radians. A nonnegative scalar in \eqn{[0, \pi]}. Defaults to \eqn{\pi / 10}.
#' @param scaled whether to scale the angular function by the normalizing
#' constant. Defaults to \code{TRUE}.
#' @return
#' Depending on the function:
#' \itemize{
#' \item \code{d_unif_cap}: a vector of length \code{nx} or \code{1} with the
#' evaluated density at \code{x}.
#' \item \code{r_unif_cap}: a matrix of size \code{c(n, p)} with the random
#' sample.
#' \item \code{c_unif_cap}: the normalizing constant.
#' \item \code{q_proj_unif_cap}: a vector of length \code{x} with the evaluated
#' distribution function at \code{x}.
#' \item \code{q_proj_unif_cap}: a vector of length \code{u} with the evaluated
#' quantile function at \code{u}.
#' \item \code{d_proj_unif_cap}: a vector of size \code{length(t)} with the
#' evaluated angular function.
#' \item \code{r_proj_unif_cap}: a vector of length \code{n} containing simulated
#' values from the cosines density associated to the angular function.
#'
#' }
#' @author Alberto Fernández-de-Marcos and Eduardo García-Portugués.
#' @examples
#' # Simulation and density evaluation for p = 2
#' mu <- c(0, 1)
#' r <- pi / 5
#' n <- 1e2
#' x <- r_unif_cap(n = n, mu = mu, r = r)
#' col <- viridisLite::viridis(n)
#' r_noise <- runif(n, 0.95, 1.05) # Radius perturbation to improve visualization
#' plot(r_noise * x, pch = 16,
#'      col = col[rank(d_unif_cap(x = x, mu = mu, r = r))],
#'      xlim = c(-1,1), ylim = c(-1,1))
#'
#' # Simulation and density evaluation for p = 3
#' mu <- c(0, 0, 1)
#' r <- pi / 5
#' x <- r_unif_cap(n = n, mu = mu, r = r)
#' scatterplot3d::scatterplot3d(x, size = 5, xlim = c(-1, 1), ylim = c(-1, 1),
#'                              zlim = c(-1, 1), color =
#'                               col[rank(d_unif_cap(x = x, mu = mu, r = r))])
#'
#' # Simulated data from the cosines density
#' n <- 1e4
#' p <- 3
#' r <- pi/3
#' hist(r_proj_unif_cap(n = n, p = p, r = r), breaks = seq(cos(r), 1, l = 10),
#'      probability = TRUE, main = "Simulated data from d_proj_unif_cap",
#'      xlab = "t", xlim = c(-1, 1))
#' t <- seq(-1, 1, by = 0.01)
#' lines(t, d_proj_unif_cap(x = t, p = p, r = r), col = "red")
#' @name unif_cap


#' @rdname unif_cap
#' @export
d_unif_cap <- function(x, mu, r = pi / 10) {

  if (is.null(dim(x))) {

    x <- rbind(x)

  }
  x <- rotasym::check_unit_norm(x = x, warnings = TRUE)
  mu <- rotasym::check_unit_norm(x = mu, warnings = TRUE)
  p <- ncol(x)
  if (p != length(mu)) {

    stop("x and mu do not have the same dimension.")

  }
  c_norm <- c_unif_cap(p, r)
  r <- 1 - cos(r)
  return(c_norm * (x %*% mu >= 1 - r))

}


#' @rdname unif_cap
#' @export
c_unif_cap <- function(p, r) {

  # Check r
  if (r < 0) {

    stop("r must be non-negative.")

  }
  if (r > pi) {

    stop("r must be lower than pi.")

  }
  r <- 1 - cos(r)
  return(1 / (rotasym::w_p(p = p - 1) *
    drop(1 - p_proj_unif(x = 1 - r, p = p))))

}


#' @rdname unif_cap
#' @export
r_unif_cap <- function(n, mu, r) {

  mu <- rotasym::check_unit_norm(x = mu, warnings = TRUE)
  p <- length(mu)
  if (r == 0) {

    return(mu)

  } else {

    r_V <- function(n) r_proj_unif_cap(n = n, p = p, r = r)
    r_U <- function(n) rotasym::r_unif_sphere(n = n, p = p - 1)
    samp <- rotasym::r_tang_norm(n = n, theta = mu, r_V = r_V, r_U = r_U)

  }
  return(samp)

}


#' @rdname unif_cap
#' @export
p_proj_unif_cap <- function(x, p, r) {

  if (r < 0) {

    stop("r must be non-negative.")

  }
  if (r > pi) {

    stop("r must be lower than pi.")

  }
  r <- 1 - cos(r)
  p1r <- drop(p_proj_unif(x = 1 - r, p = p))
  return(drop((p_proj_unif(x = x, p = p) - p1r) / (1 - p1r)))

}


#' @rdname unif_cap
#' @export
q_proj_unif_cap <- function(u, p, r) {

  if (r < 0) {

    stop("r must be non-negative.")

  }
  if (r > pi) {

    stop("r must be lower than pi.")

  }
  r <- 1 - cos(r)
  p1r <- drop(p_proj_unif(x = 1 - r, p = p))
  u <- (1 - p1r) * u + p1r
  u <- pmax(pmin(u, 1), 0)
  return(drop(q_proj_unif(u = u, p = p)))

}


#' @rdname unif_cap
#' @export
d_proj_unif_cap <- function(x, p, r, scaled = TRUE) {

  if (r < 0) {

    stop("r must be non-negative.")

  }
  if (r > pi) {

    stop("r must be lower than pi.")

  }
  # return(ifelse(scaled, c_unif_cap(p = p, r = r), 1) *
  #          ((t >= (1 - (1 - cos(r)))) * (t <= 1)))
  return(ifelse(scaled, c_unif_cap(p = p, r = r), 1) *
           ((x >= (1 + cos(r))) * (x <= 1)))

}


#' @rdname unif_cap
#' @export
r_proj_unif_cap <- function(n, p, r) {

  U <- runif(n = n)
  return(q_proj_unif_cap(U, p = p, r = r))

}

