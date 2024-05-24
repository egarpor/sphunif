

#' @title Statistics for testing (hyper)spherical uniformity
#'
#' @description Low-level implementation of several statistics for assessing
#' uniformity on the (hyper)sphere
#' \eqn{S^{p-1}:=\{{\bf x}\in R^p:||{\bf x}||=1\}}{S^{p-1}:=
#' \{x\in R^p:||x||=1\}}, \eqn{p\ge 2}.
#'
#' @param X an \bold{array} of size \code{c(n, p, M)} containing the Cartesian
#' coordinates of \code{M} samples of size \code{n} of directions on
#' \eqn{S^{p-1}}. Must not contain \code{NA}'s.
#' @param Psi_in_X does \code{X} contain the shortest angles matrix
#' \eqn{\boldsymbol\Psi}{\Psi} that is obtained with \code{\link{Psi_mat}(X)}?
#' If \code{FALSE} (default), \eqn{\boldsymbol\Psi}{\Psi} is computed
#' internally.
#' @inheritParams unif_stat_distr
#' @inheritParams unif_stat
#' @inheritParams sph_stat_distr
#' @inheritParams cir_stat
#' @param dirs a matrix of size \code{c(n_proj, p)} containing \code{n_proj}
#' random directions (in Cartesian coordinates) on \eqn{S^{p-1}} to perform
#' the CCF09 test.
#' @param N number of points used in the
#' \link[=Gauss_Legen_nodes]{Gauss-Legendre quadrature}. Defaults to
#' \code{160}.
#' @param L number of discretization points to interpolate angular functions
#' that require evaluating an integral. Defaults to \code{1e3}.
#' @return A matrix of size \code{c(M, 1)} containing the statistics for each
#' of the \code{M} samples.
#' @section Warning:
#' Be careful on avoiding the next bad usages of the functions, which will
#' produce spurious results:
#' \itemize{
#'   \item The directions in \code{X} do \emph{not} have unit norm.
#'   \item \code{X} does \emph{not} contain \code{Psi_mat(X)} when
#'   \code{X_in_Theta = TRUE}.
#'   \item The parameter \code{p} does \emph{not} match with the dimension of
#'   \eqn{R^p}.
#'   \item \emph{Not} passing the scalar products matrix to \code{sph_stat_CJ12}
#'   when \code{Psi_in_X = TRUE}.
#'   \item The directions in \code{dirs} do \emph{not} have unit norm.
#' }
#' @details
#' Detailed descriptions and references of the statistics are available
#' in García-Portugués and Verdebout (2018).
#'
#' The Pycke and CJ12 statistics employ the \emph{scalar products} matrix,
#' rather than the shortest angles matrix, when \code{Psi_in_X = TRUE}. This
#' matrix is obtained by setting \code{scalar_prod = TRUE} in
#' \code{\link{Psi_mat}}.
#' @references
#' García-Portugués, E. and Verdebout, T. (2018) An overview of uniformity
#' tests on the hypersphere. \emph{arXiv:1804.00286}.
#' \doi{10.48550/arXiv.1804.00286}.
#' @examples
#' ## Sample uniform spherical data
#'
#' M <- 2
#' n <- 100
#' p <- 3
#' set.seed(123456789)
#' X <- r_unif_sph(n = n, p = p, M = M)
#'
#' ## Sobolev tests
#'
#' # Rayleigh
#' sph_stat_Rayleigh(X)
#'
#' # Bingham
#' sph_stat_Bingham(X)
#'
#' # Ajne
#' Psi <- Psi_mat(X)
#' dim(Psi) <- c(dim(Psi), 1)
#' sph_stat_Ajne(X)
#' sph_stat_Ajne(Psi, Psi_in_X = TRUE)
#'
#' # Gine Gn
#' sph_stat_Gine_Gn(X)
#' sph_stat_Gine_Gn(Psi, Psi_in_X = TRUE, p = p)
#'
#' # Gine Fn
#' sph_stat_Gine_Fn(X)
#' sph_stat_Gine_Fn(Psi, Psi_in_X = TRUE, p = p)
#'
#' # Pycke
#' sph_stat_Pycke(X)
#' sph_stat_Pycke(Psi, Psi_in_X = TRUE, p = p)
#'
#' # Bakshaev
#' sph_stat_Bakshaev(X)
#' sph_stat_Bakshaev(Psi, Psi_in_X = TRUE, p = p)
#'
#' # Riesz
#' sph_stat_Riesz(X, s = 1)
#' sph_stat_Riesz(Psi, Psi_in_X = TRUE, p = p, s = 1)
#'
#' # Projected Cramér-von Mises
#' sph_stat_PCvM(X)
#' sph_stat_PCvM(Psi, Psi_in_X = TRUE, p = p)
#'
#' # Projected Rothman
#' sph_stat_PRt(X)
#' sph_stat_PRt(Psi, Psi_in_X = TRUE, p = p)
#'
#' # Projected Anderson-Darling
#' sph_stat_PAD(X)
#' sph_stat_PAD(Psi, Psi_in_X = TRUE, p = p)
#'
#' ## Other tests
#'
#' # CCF09
#' dirs <- r_unif_sph(n = 3, p = p, M = 1)[, , 1]
#' sph_stat_CCF09(X, dirs = dirs)
#'
#' ## High-dimensional tests
#'
#' # Rayleigh HD-Standardized
#' sph_stat_Rayleigh_HD(X)
#'
#' # CJ12
#' sph_stat_CJ12(X, regime = 1)
#' sph_stat_CJ12(Psi, regime = 1, Psi_in_X = TRUE, p = p)
#' sph_stat_CJ12(X, regime = 2)
#' sph_stat_CJ12(Psi, regime = 2, Psi_in_X = TRUE, p = p)
#' sph_stat_CJ12(X, regime = 3)
#' sph_stat_CJ12(Psi, regime = 3, Psi_in_X = TRUE, p = p)
#' @name sph_stat
NULL
