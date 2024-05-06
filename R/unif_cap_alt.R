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
#' @param x locations in \eqn{S^{p-1}} to evaluate the density.
#' Either a matrix of size \code{c(nx, p)} or a vector of length \code{p}.
#' Normalized internally if required (with a warning message).
#' @param n sample size, a positive integer.
#' @param p dimension of the ambient space \eqn{R^p} that contains
#' \eqn{S^{p-1}}. A positive integer.
#' @param u vector of probabilities.
#' @param mu the directional mean \eqn{\boldsymbol{\mu}}{\mu} of the
#' distribution. A unit-norm vector of length \code{p}.
#' @param r radius of the spherical cap in radians.
#' A nonnegative scalar in \eqn{[0, \pi]}.
#' @param t a vector with values in \eqn{[-1, 1]}.
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
#' \item \code{q_unif_cap}: a vector of length \code{u} with the evaluated
#' quantile function at \code{u}.
#' \item \code{g_unif_cap}: a vector of size \code{length(t)} with the
#' evaluated angular function.
#' \item \code{r_g_unif_cap}: a vector of length \code{n} containing simulated
#' values from the cosines density associated to the angular function.
#'
#' }
#' @author Alberto Fernández-de-Marcos and Eduardo García-Portugués.
#' @examples
#' # Simulation and density evaluation for p = 2
#' mu <- c(0, 1)
#' r <- pi/5
#' n <- 1e2
#' x <- r_unif_cap(n = n, mu = mu, r = r)
#' col <- viridisLite::viridis(n)
#' r_noise <- runif(n, 0.95, 1.05)# Radius perturbation to improve visualization
#' plot(r_noise * x, pch = 16,
#'      col = col[rank(d_unif_cap(x = x, mu = mu, r = r))],
#'      xlim = c(-1,1), ylim = c(-1,1))
#'
#' # Simulation and density evaluation for p = 3
#' mu <- c(0, 0, 1)
#' r <- pi/5
#' x <- r_unif_cap(n = n, mu = mu, r = r)
#' if (requireNamespace("rgl")) {
#'   rgl::plot3d(x, col = col[rank(d_unif_cap(x = x, mu = mu, r = r))],
#'               size = 5, xlim = c(-1,1), ylim = c(-1,1), zlim = c(-1,1))
#' }
#'
#' # Cosines density
#' g_tilde <- function(t, p, r) {
#' r <- 1 - cos(r)
#' q <- p - 1
#' return((1-t^2)^(q/2 - 1)/drop((beta(1/2, q/2) * (
#'   1 - p_proj_unif(x = 1 - r, p = p)))) * ((t >= 1 - r)*(t <= 1)))
#' }
#'
#' # Simulated data from the cosines density
#' n <- 1e4
#' p <- 3
#' r <- pi/3
#' hist(r_g_unif_cap(n = n, p = p, r = r), breaks = seq(cos(r), 1, l = 10),
#'      probability = TRUE, main = "Simulated data from g_unif_cap", xlab = "t",
#'      xlim = c(-1,1))
#' t <- seq(-1, 1, by = 0.01)
#' lines(t, g_tilde(t = t, p = p, r = r), col = "red")
#'
#' # Cosine density as a function of the dimension
#' M <- 100
#' col <- viridisLite::viridis(M)
#' plot(t, g_tilde(t = t, p = 2, r = r), col = col[2], type = "l",
#'      ylab = "Density")
#' for (p in 3:M) {
#'   lines(t, g_tilde(t = t, p = p, r = r), col = col[p])
#' }

#' @name unif_cap_alt

#' @rdname unif_cap_alt
#' @export
d_unif_cap <- function(x, mu, r) {
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
  return(c_norm * (x %*% mu >= 1 - r) * (x %*% mu <= 1))
}

#' @rdname unif_cap_alt
#' @export
c_unif_cap <- function(p, r) {
  if (r < 0) {
    stop("radius must be non-negative.")
  }
  if (r > pi) {
    stop("radius must be lower than pi.")
  }
  r <- 1 - cos(r)
  1 / (rotasym::w_p(p = p - 1) * drop((1 - p_proj_unif(x = 1 - r, p = p))))
}

#' @rdname unif_cap_alt
#' @export
r_unif_cap <- function(n, mu, r) {

  mu <- rotasym::check_unit_norm(x = mu, warnings = TRUE)
  p <- length(mu)
  if (r == 0) {
    return(mu)
  } else {
    r_V <- function(n) r_g_unif_cap(n = n, p = p, r = r)
    r_U <- function(n) rotasym::r_unif_sphere(n = n, p = p - 1)
    samp <- rotasym::r_tang_norm(n = n, theta = mu, r_V = r_V,
                                 r_U = r_U)
  }
  return(samp)
}

#' @rdname unif_cap_alt
#' @export
q_unif_cap <- function(u, p, r) {
  if (r < 0) {
    stop("radius must be non-negative.")
  }
  if (r > pi) {
    stop("radius must be lower than pi.")
  }
  r <- 1 - cos(r)
  drop(q_proj_unif(u = drop(1 - p_proj_unif(x = 1 - r,
                                            p = p)
                            ) * u + drop(p_proj_unif(x = 1 - r, p = p)), p = p))
}

#' @rdname unif_cap_alt
#' @export
g_unif_cap <- function(t, p, r, scaled = TRUE) {
  if (r < 0) {
    stop("radius must be non-negative.")
  }
  if (r > pi) {
    stop("radius must be lower than pi.")
  }
  return(ifelse(scaled, c_unif_cap(p = p, r = r),
                1) * ((t >= 1 - (1 - cos(r))) * (t <= 1)))
}

#' @rdname unif_cap_alt
#' @export
r_g_unif_cap <- function(n, p, r) {
  U <- runif(n = n)
  return(q_unif_cap(U, p = p, r = r))
}
