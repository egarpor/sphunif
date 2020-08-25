

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
#' \code{\link[Rcpp]{Rcpp-package}}.
#'
#' @author Eduardo García-Portugués and Thomas Verdebout.
#' @references
#' García-Portugués, E. and Verdebout, T. (2018) An overview of uniformity
#' tests on the hypersphere. \emph{arXiv:1804.00286}.
#' \url{https://arxiv.org/abs/1804.00286}.
#' @docType package
#' @name sphunif-package
#' @import graphics stats foreach parallel doSNOW Rcpp
#' @useDynLib sphunif
#' @aliases sphunif sphunif-package
NULL
