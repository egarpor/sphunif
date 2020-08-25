

#' @title Statistics for testing circular uniformity
#'
#' @description Low-level implementation of several statistics for assessing
#' circular uniformity on \eqn{[0, 2\pi)} or, equivalently,
#' \eqn{S^1:=\{{\bf x}\in R^2:||{\bf x}||=1\}}{S^1:=\{x\in R^2:||x||=1\}}.
#'
#' @param Theta a \bold{matrix} of size \code{c(n, M)} with \code{M} samples
#' of size \code{n} of circular data on \eqn{[0, 2\pi)}. Must not contain
#' \code{NA}'s.
#' @param gaps_in_Theta does \code{Theta} contain the matrix of
#' \emph{circular gaps} that is obtained with \cr\code{\link{cir_gaps}(Theta)}?
#' If \code{FALSE} (default), the circular gaps are computed internally.
#' @inheritParams cir_gaps
#' @param Psi_in_Theta does \code{Theta} contain the shortest angles matrix
#' \eqn{\boldsymbol\Psi}{\Psi} that is obtained with
#' \cr\code{\link{Psi_mat}(array(Theta, dim = c(n, 1, M)))}? If \code{FALSE}
#' (default), \eqn{\boldsymbol\Psi}{\Psi} is computed internally.
#' @param KS compute the Kolmogorov-Smirnov statistic (which is
#' \emph{not} invariant under origin shifts) instead of the Kuiper statistic?
#' Defaults to \code{FALSE}.
#' @param CvM compute the Cramér-von Mises statistic (which is \emph{not}
#' invariant under origin shifts) instead of the Watson statistic? Defaults to
#' \code{FALSE}.
#' @param m integer \eqn{m} for the \eqn{m}-modal Rayleigh test. Defaults to
#' \code{m = 1} (the standard Rayleigh test).
#' @param max_gap compute the maximum gap for the range statistic? If
#' \code{TRUE} (default), rejection happens for \emph{large} values of
#' the statistic, which is consistent with the rest of tests. Otherwise,
#' the minimum gap is computed and rejection happens for \emph{low} values.
#' @param asymp_std normalize the Hodges-Ajne statistic in terms of its
#' asymptotic distribution? Defaults to \code{FALSE}.
#' @param use_Cressie compute the Hodges-Ajne statistic as a particular case
#' of the Cressie statistic? Defaults to \code{TRUE} as it is more efficient.
#' If \code{FALSE}, the geometric construction in Ajne (1968) is employed.
#' @param minus compute the invariant \eqn{D_n^-} instead of \eqn{D_n^+}?
#' Defaults to \code{FALSE}.
#' @param a \eqn{a_n = a / n} parameter used in the length of the arcs
#' of the coverage-based tests. Must be positive. Defaults to \code{2 * pi}.
#' @param t \eqn{t} parameter for the Rothman and Cressie tests, a real in
#' \eqn{(0, 1)}. Defaults to \code{1 / 3}.
#' @param q \eqn{q} parameter for the Pycke "\eqn{q}-test", a real in
#' \eqn{(0, 1)}. Defaults to \code{1 / 2}.
#' @param abs_val return the absolute value of the Darling's log gaps
#' statistic? If \code{TRUE} (default), rejection happens for \emph{large}
#' values of the statistic, which is consistent with the rest of tests.
#' Otherwise, the signed statistic is computed and rejection happens for
#' large \emph{absolute} values.
#' @param minus_val return the negative value of the (standardized) number of
#' uncovered spacings? If \code{TRUE} (default), rejection happens for
#' \emph{large} values of the statistic, which is consistent with the rest of
#' tests. Otherwise, rejection happens for \emph{low} values.
#' @param rand_dirs a matrix of size \code{c(n_proj, 2)} containing
#' \code{n_proj} random directions (in Cartesian coordinates) on \eqn{S^1} to
#' perform the Cuesta-Albertos test.
#' @param K_Cuesta_Albertos integer giving the truncation of the series present
#' in the asymptotic distribution of the Kolmogorov-Smirnov statistic. Defaults
#' to \code{25}.
#' @param original return the Cuesta-Albertos statistic as originally defined?
#' If \code{FALSE} (default), a faster and equivalent statistic is computed,
#' and rejection happens for \emph{large} values of the statistic, which is
#' consistent with the rest of tests. Otherwise, rejection happens for
#' \emph{low} values.
#' @param Stephens compute Stephens (1970) modification so that the null
#' distribution of the is less dependent on the sample size? The modification
#' does not alter the test decision.
#' @return A matrix of size \code{c(M, 1)} containing the statistics for each
#' of the \code{M} samples.
#' @section Warning:
#' Be careful on avoiding the next bad usages of the functions, which will
#' produce spurious results:
#' \itemize{
#'   \item The entries of \code{Theta} are \emph{not} in \eqn{[0, 2\pi)}.
#'   \item \code{Theta} does \emph{not} contain the circular gaps when
#'   \code{gaps_in_Theta = TRUE}.
#'   \item \code{Theta} is \emph{not} sorted increasingly when
#'   \code{data_sorted = TRUE}.
#'   \item \code{Theta} does \emph{not} contain
#'   \code{\link{Psi_mat}(array(Theta, dim = c(n, 1, M)))} when
#'   \cr\code{Psi_in_Theta = TRUE}.
#'   \item The directions in \code{rand_dirs} do \emph{not} have unit norm.
#' }
#' @details
#' Detailed descriptions and references of the statistics are available
#' in García-Portugués and Verdebout (2018).
#'
#' The statistics \code{cir_stat_PCvM} and \code{cir_stat_PRt} are provided
#' for the sake of completion, but they equal the more efficiently-implemented
#' statistics \code{2 * cir_stat_Watson} and \code{cir_stat_Rothman},
#' respectively.
#' @references
#' García-Portugués, E. and Verdebout, T. (2018) An overview of uniformity
#' tests on the hypersphere. \emph{arXiv:1804.00286}.
#' \url{https://arxiv.org/abs/1804.00286}.
#' @examples
#' ## Sample uniform circular data
#'
#' M <- 2
#' n <- 100
#' set.seed(987202226)
#' Theta <- r_unif_cir(n = n, M = M)
#'
#' ## Tests based on the empirical cumulative distribution function
#'
#' # Kuiper
#' cir_stat_Kuiper(Theta)
#' cir_stat_Kuiper(Theta, Stephens = TRUE)
#'
#' # Watson
#' cir_stat_Watson(Theta)
#' cir_stat_Watson(Theta, Stephens = TRUE)
#'
#' # Watson (1976)
#' cir_stat_Watson_1976(Theta)
#'
#' ## Partition-based tests
#'
#' # Ajne
#' Theta_array <- Theta
#' dim(Theta_array) <- c(nrow(Theta), 1, ncol(Theta))
#' Psi <- Psi_mat(Theta_array)
#' cir_stat_Ajne(Theta)
#' cir_stat_Ajne(Psi, Psi_in_Theta = TRUE)
#'
#' # Rothman
#' cir_stat_Rothman(Theta, t = 0.5)
#' cir_stat_Rothman(Theta)
#' cir_stat_Rothman(Psi, Psi_in_Theta = TRUE)
#'
#' # Hodges-Ajne
#' cir_stat_Hodges_Ajne(Theta)
#' cir_stat_Hodges_Ajne(Theta, use_Cressie = FALSE)
#'
#' # Cressie
#' cir_stat_Cressie(Theta, t = 0.5)
#' cir_stat_Cressie(Theta)
#'
#' # Feltz-Goldin
#' cir_stat_Feltz_Goldin(Theta)
#'
#' ## Spacings-based tests
#'
#' # Range
#' cir_stat_Range(Theta)
#'
#' # Rao
#' cir_stat_Rao(Theta)
#'
#' # Greenwood
#' cir_stat_Greenwood(Theta)
#'
#' # Log gaps
#' cir_stat_Log_gaps(Theta)
#'
#' # Vacancy
#' cir_stat_Vacancy(Theta)
#'
#' # Maximum uncovered spacing
#' cir_stat_Max_uncover(Theta)
#'
#' # Number of uncovered spacings
#' cir_stat_Num_uncover(Theta)
#'
#' # Gini mean difference
#' cir_stat_Gini(Theta)
#'
#' # Gini mean squared difference
#' cir_stat_Gini_squared(Theta)
#'
#' ## Sobolev tests
#'
#' # Rayleigh
#' cir_stat_Rayleigh(Theta)
#' cir_stat_Rayleigh(Theta, m = 2)
#'
#' # Bingham
#' cir_stat_Bingham(Theta)
#'
#' # Hermans-Rasson
#' cir_stat_Hermans_Rasson(Theta)
#' cir_stat_Hermans_Rasson(Psi, Psi_in_Theta = TRUE)
#'
#' # Gine Fn
#' cir_stat_Gine_Fn(Theta)
#' cir_stat_Gine_Fn(Psi, Psi_in_Theta = TRUE)
#'
#' # Gine Gn
#' cir_stat_Gine_Gn(Theta)
#' cir_stat_Gine_Gn(Psi, Psi_in_Theta = TRUE)
#'
#' # Pycke
#' cir_stat_Pycke(Theta)
#' cir_stat_Pycke(Psi, Psi_in_Theta = TRUE)
#'
#' # Pycke q
#' cir_stat_Pycke_q(Theta)
#' cir_stat_Pycke_q(Psi, Psi_in_Theta = TRUE)
#'
#' # Bakshaev
#' cir_stat_Bakshaev(Theta)
#' cir_stat_Bakshaev(Psi, Psi_in_Theta = TRUE)
#'
#' # Projected Cramér-von Mises
#' cir_stat_PCvM(Theta)
#' cir_stat_PCvM(Psi, Psi_in_Theta = TRUE)
#'
#' # Projected Rothman
#' cir_stat_PRt(Theta, t = 0.5)
#' cir_stat_PRt(Theta)
#' cir_stat_PRt(Psi, Psi_in_Theta = TRUE)
#'
#' # Projected Anderson-Darling
#' cir_stat_PAD(Theta)
#' cir_stat_PAD(Psi, Psi_in_Theta = TRUE)
#'
#' ## Other tests
#'
#' # Cuesta-Albertos
#' rand_dirs <- r_unif_sph(n = 3, p = 2, M = 1)[, , 1]
#' cir_stat_Cuesta_Albertos(Theta, rand_dirs = rand_dirs)
#'
#' ## Connection of Kuiper and Watson statistics with KS and CvM, respectively
#'
#' # Rotate sample for KS and CvM
#' alpha <- seq(0, 2 * pi, l = 1e4)
#' KS_alpha <- sapply(alpha, function(a) {
#'   cir_stat_Kuiper((Theta[, 2, drop = FALSE] + a) %% (2 * pi), KS = TRUE)
#' })
#' CvM_alpha <- sapply(alpha, function(a) {
#'   cir_stat_Watson((Theta[, 2, drop = FALSE] + a) %% (2 * pi), CvM = TRUE)
#' })
#'
#' # Kuiper is the maximum rotated KS
#' plot(alpha, KS_alpha, type = "l")
#' abline(h = cir_stat_Kuiper(Theta[, 2, drop = FALSE]), col = 2)
#' points(alpha[which.max(KS_alpha)], max(KS_alpha), col = 2, pch = 16)
#'
#' # Watson is the minimum rotated CvM
#' plot(alpha, CvM_alpha, type = "l")
#' abline(h = cir_stat_Watson(Theta[, 2, drop = FALSE]), col = 2)
#' points(alpha[which.min(CvM_alpha)], min(CvM_alpha), col = 2, pch = 16)
#' @name cir_stat
NULL
