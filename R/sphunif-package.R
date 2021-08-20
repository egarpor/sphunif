

#' @title \code{sphunif} -- Uniformity Tests on the Circle, Sphere, and
#' Hypersphere
#'
#' @description Implementation of uniformity tests on the circle and
#' (hyper)sphere. The main function of the package is \code{\link{unif_test}},
#' which conveniently collects more than 30 tests for assessing uniformity on
#' \eqn{S^{p-1}=\{{\bf x}\in R^p:||{\bf x}||=1\}}{
#' S^{p-1}=\{x\in R^p:||x||=1\}}, \eqn{p\ge 2}. The test statistics are
#' implemented in the \code{\link{unif_stat}} function, which allows to compute
#' several statistics to different samples within a single call, facilitating
#' thus Monte Carlo experiments. Furthermore, the
#' \code{\link{unif_stat_MC}} function allows to parallelize them in
#' a simple way. The asymptotic null distributions of the statistics are
#' available through the function \code{\link{unif_stat_distr}}. The core of
#' the \code{sphunif} is coded in C++ by relying on the
#' \code{\link[Rcpp]{Rcpp-package}}. The package allows the replication of
#' the data application in García-Portugués, Navarro-Esteban and
#' Cuesta-Albertos (2020) <arXiv:2008.09897>.
#'
#' @author Eduardo García-Portugués and Thomas Verdebout.
#' @references
#' García-Portugués, E. and Verdebout, T. (2018) An overview of uniformity
#' tests on the hypersphere. \emph{arXiv:1804.00286}.
#' \url{https://arxiv.org/abs/1804.00286}.
#'
#' García-Portugués, E., Navarro-Esteban, P., Cuesta-Albertos, J. A. (2020)
#' On a projection-based class of uniformity tests on the hypersphere.
#' \emph{arXiv:2008.09897}. \url{https://arxiv.org/abs/2008.09897}
#'
#' García-Portugués, E., Paindaveine, D., and Verdebout, T. (2021). On the
#' power of Sobolev tests for isotropy under local rotationally symmetric
#' alternatives. \emph{Submitted}
#' @docType package
#' @name sphunif-package
#' @import graphics stats foreach parallel doSNOW Rcpp
#' @useDynLib sphunif
#' @aliases sphunif sphunif-package
NULL
