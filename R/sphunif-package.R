

#' @title \code{sphunif}: Uniformity Tests on the Circle, Sphere, and
#' Hypersphere
#'
#' @description Implementation of uniformity tests on the circle and
#' (hyper)sphere. The main function of the package is \code{\link{unif_test}},
#' which conveniently collects more than 35 tests for assessing uniformity on
#' \eqn{S^{p-1}=\{{\bf x}\in R^p:||{\bf x}||=1\}}{
#' S^{p-1}=\{x\in R^p:||x||=1\}}, \eqn{p\ge 2}. The test statistics are
#' implemented in the \code{\link{unif_stat}} function, which allows computing
#' several statistics for different samples within a single call, thus
#' facilitating Monte Carlo experiments. Furthermore, the
#' \code{\link{unif_stat_MC}} function allows parallelizing them in
#' a simple way. The asymptotic null distributions of the statistics are
#' available through the function \code{\link{unif_stat_distr}}. The core of
#' \code{\link{sphunif-package}} is coded in C++ by relying on the
#' \code{\link[Rcpp]{Rcpp-package}}. The package also provides several
#' novel datasets and gives the replicability for the data applications/
#' simulations in García-Portugués et al. (2021)
#' <doi:10.1007/978-3-030-69944-4_12>, García-Portugués et al. (2023)
#' <doi:10.3150/21-BEJ1454>, García-Portugués et al. (2024)
#' <arXiv:2108.09874v2>, and Fernández-de-Marcos and García-Portugués (2024)
#' <arXiv:2405.13531>.
#'
#' @author Eduardo García-Portugués and Thomas Verdebout.
#' @references
#' Fernández-de-Marcos, A. and García-Portugués, E. (2024) A stereographic test
#' of spherical uniformity. \emph{arXiv:2405.13531}.
#' \url{https://arxiv.org/abs/2405.13531}.
#'
#' García-Portugués, E. and Verdebout, T. (2018) An overview of uniformity
#' tests on the hypersphere. \emph{arXiv:1804.00286}.
#' \url{https://arxiv.org/abs/1804.00286}.
#'
#' García-Portugués, E., Navarro-Esteban, P., Cuesta-Albertos, J. A. (2023)
#' On a projection-based class of uniformity tests on the hypersphere.
#' \emph{Bernoulli}, 29(1):181--204. \doi{10.3150/21-BEJ1454}.
#'
#' García-Portugués, E., Navarro-Esteban, P., and Cuesta-Albertos, J. A. (2021).
#' A Cramér–von Mises test of uniformity on the hypersphere. In Balzano, S.,
#' Porzio, G. C., Salvatore, R., Vistocco, D., and Vichi, M. (Eds.), \emph{
#' Statistical Learning and Modeling in Data Analysis}, Studies in
#' Classification, Data Analysis and Knowledge Organization, pp. 107--116.
#' Springer, Cham. \doi{10.1007/978-3-030-69944-4_12}.
#'
#' García-Portugués, E., Paindaveine, D., and Verdebout, T. (2024). On a class
#' of Sobolev tests for symmetry of directions, their detection thresholds, and
#' asymptotic powers. \emph{arXiv:2108.09874v2}.
#' \url{https://arxiv.org/abs/2108.09874v2}
#' @docType package
#' @name sphunif-package
#' @import Rcpp
#' @importFrom graphics abline legend lines par segments
#' @importFrom stats approxfun dchisq dgamma integrate na.omit nlm pchisq pgamma
#' qchisq qgamma quantile rbeta rchisq runif splinefun uniroot
#' @importFrom utils capture.output
#' @useDynLib sphunif
#' @aliases sphunif sphunif-package
NULL
