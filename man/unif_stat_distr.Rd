% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unif_stat_distr.R
\name{unif_stat_distr}
\alias{unif_stat_distr}
\title{Null distributions for circular and (hyper)spherical uniformity
statistics}
\usage{
unif_stat_distr(x, type, p, n, approx = "asymp", M = 10000,
  stats_MC = NULL, K_max = 10000, method = "I", Stephens = FALSE,
  CCF09_dirs = NULL, CJ12_beta = 0, CJ12_reg = 3, cov_a = 2 * pi,
  Cressie_t = 1/3, K_Ajne = 500, K_CCF09 = 25, K_Kuiper = 25,
  K_Watson = 25, K_Watson_1976 = 5, Poisson_rho = 0.5, Pycke_q = 0.5,
  Rayleigh_m = 1, Riesz_s = 1, Rothman_t = 1/3, Sobolev_vk2 = c(0, 0,
  1), Softmax_kappa = 1, Stein_K = 10, Stein_cf = FALSE, Stereo_a = 0,
  ...)
}
\arguments{
\item{x}{evaluation points for the null distribution(s). Either a vector of
size \code{nx}, if the evaluation points are common for the tests in
\code{type}, or a matrix of size \code{c(nx, length(type))} with columns
containing the evaluation points for each test. Must not contain \code{NA}'s.}

\item{type}{type of test to be applied. A character vector containing any of
the following types of tests, depending on the dimension \eqn{p}:
\itemize{
  \item Circular data: any of the names available at object
  \code{\link{avail_cir_tests}}.
  \item (Hyper)spherical data: any of the names available at object
  \code{\link{avail_sph_tests}}.
}
If \code{type = "all"} (default), then \code{type} is set as
\code{avail_cir_tests} or \code{avail_sph_tests}, depending on the value of
\eqn{p}.}

\item{p}{integer giving the dimension of the ambient space \eqn{R^p} that
contains \eqn{S^{p-1}}.}

\item{n}{sample size employed for computing the statistic.}

\item{approx}{type of approximation to the null distribution, either
\code{"asymp"} (default) for employing the asymptotic null distribution, if
available, or \code{"MC"}, for employing the Monte Carlo approximation of
the exact null distribution.}

\item{M}{number of Monte Carlo replications for approximating the null
distribution when \code{approx = "MC"}. Also, number of Monte Carlo samples
for approximating the asymptotic distributions based on weighted sums of chi
squared random variables. Defaults to \code{1e4}.}

\item{stats_MC}{a data frame of size \code{c(M, length(type))}, with column
names containing the character vector \code{type}, that results from
extracting \code{$stats_MC} from a call to \code{\link{unif_stat_MC}}. If
provided, the computation of Monte Carlo statistics when \code{approx = "MC"}
is skipped. \code{stats_MC} is checked internally to see if it is sorted.
Internally computed if \code{NULL} (default).}

\item{K_max}{integer giving the truncation of the series that compute the
asymptotic p-value of a Sobolev test. Defaults to \code{1e4}.}

\item{method}{method for approximating the density, distribution, or
quantile function of the weighted sum of chi squared random variables. Must
be \code{"I"} (Imhof), \code{"SW"} (Satterthwaite--Welch), \code{"HBE"}
(Hall--Buckley--Eagleson), or \code{"MC"} (Monte Carlo; only for distribution
or quantile functions). Defaults to \code{"I"}.}

\item{Stephens}{compute Stephens (1970) modification so that the null
distribution of the is less dependent on the sample size? The modification
does not alter the test decision.}

\item{CCF09_dirs}{a matrix of size \code{c(n_proj, p)} containing
\code{n_proj} random directions (in Cartesian coordinates) on \eqn{S^{p-1}}
to perform the CCF09 test. If \code{NULL} (default), a sample of size
\code{n_proj = 50} directions is computed internally.}

\item{CJ12_beta}{\eqn{\beta} parameter in the exponential regime of CJ12
test, a positive real.}

\item{CJ12_reg}{type of asymptotic regime for CJ12 test, either \code{1}
(sub-exponential regime), \code{2} (exponential), or \code{3}
(super-exponential; default).}

\item{cov_a}{\eqn{a_n = a / n} parameter used in the length of the arcs
of the coverage-based tests. Must be positive. Defaults to \code{2 * pi}.}

\item{Cressie_t}{\eqn{t} parameter for the Cressie test, a real in
\eqn{(0, 1)}. Defaults to \code{1 / 3}.}

\item{K_CCF09}{integer giving the truncation of the series present in the
asymptotic distribution of the Kolmogorov-Smirnov statistic. Defaults to
\code{25}.}

\item{K_Kuiper, K_Watson, K_Watson_1976, K_Ajne}{integer giving the truncation
of the series present in the null asymptotic distributions. For the
Kolmogorov-Smirnov-related series defaults to \code{25}.}

\item{Poisson_rho}{\eqn{\rho} parameter for the Poisson test, a real in
\eqn{[0, 1)}. Defaults to \code{0.5}.}

\item{Pycke_q}{\eqn{q} parameter for the Pycke "\eqn{q}-test", a real in
\eqn{(0, 1)}. Defaults to \code{1 / 2}.}

\item{Rayleigh_m}{integer \eqn{m} for the \eqn{m}-modal Rayleigh test.
Defaults to \code{m = 1} (the standard Rayleigh test).}

\item{Riesz_s}{\eqn{s} parameter for the \eqn{s}-Riesz test, a real in
\eqn{(0, 2)}. Defaults to \code{1}.}

\item{Rothman_t}{\eqn{t} parameter for the Rothman test, a real in
\eqn{(0, 1)}. Defaults to \code{1 / 3}.}

\item{Sobolev_vk2}{weights for the finite Sobolev test. A non-negative
vector or matrix. Defaults to \code{c(0, 0, 1)}.}

\item{Softmax_kappa}{\eqn{\kappa} parameter for the Softmax test, a
non-negative real. Defaults to \code{1}.}

\item{Stein_K}{truncation \eqn{K} parameter for the Stein test, a positive
integer. Defaults to \code{10}.}

\item{Stein_cf}{logical indicating whether to use the characteristic
function in the Stein test. Defaults to \code{FALSE} (moment generating
function).}

\item{Stereo_a}{\eqn{a} parameter for the Stereo test, a real in
\eqn{[-1, 1]}. Defaults to \code{0}.}

\item{...}{if \code{approx = "MC"}, optional performance parameters to be
passed to \cr\code{\link{unif_stat_MC}}: \code{chunks}, \code{cores},
and \code{seed}.}
}
\value{
A data frame of size \code{c(nx, length(type))}, with column names
given by \code{type}, that contains the values of the null distributions of
the statistics evaluated at \code{x}.
}
\description{
Approximate computation of the null distributions of several
statistics for assessing uniformity on the (hyper)sphere
\eqn{S^{p-1}:=\{{\bf x}\in R^p:||{\bf x}||=1\}}{S^{p-1}:=
\{x\in R^p:||x||=1\}}, \eqn{p\ge 2}. The approximation is done either by
means of the asymptotic distribution or by Monte Carlo.
}
\details{
When \code{approx = "asymp"}, statistics that do not have an implemented or
known asymptotic are omitted, and a warning is generated.

For Sobolev tests, \code{K_max = 1e4} produces probabilities uniformly
accurate with three digits for the \code{"PCvM"}, \code{"PAD"}, and
\code{"PRt"} tests, for dimensions \eqn{p \le 11}. With \code{K_max = 5e4},
these probabilities are uniformly accurate in the fourth digit. With
\code{K_max = 1e3}, only two-digit uniform accuracy is obtained. Uniform
accuracy deteriorates when \eqn{p} increases, e.g., a digit accuracy is lost
when \eqn{p = 51}.

Descriptions and references on most of the asymptotic distributions
are available in García-Portugués and Verdebout (2018).
}
\examples{
## Asymptotic distribution

# Circular statistics
x <- seq(0, 1, l = 5)
unif_stat_distr(x = x, type = "Kuiper", p = 2, n = 10)
unif_stat_distr(x = x, type = c("Ajne", "Kuiper"), p = 2, n = 10)
unif_stat_distr(x = x, type = c("Ajne", "Kuiper"), p = 2, n = 10, K_Ajne = 5)
\donttest{
# All circular statistics
unif_stat_distr(x = x, type = avail_cir_tests, p = 2, n = 10, K_max = 1e3)
}
# Spherical statistics
unif_stat_distr(x = cbind(x, x + 1), type = c("Rayleigh", "Bingham"),
                p = 3, n = 10)
unif_stat_distr(x = cbind(x, x + 1), type = c("Rayleigh", "Bingham"),
                p = 3, n = 10, M = 100)
\donttest{
# All spherical statistics
unif_stat_distr(x = x, type = avail_sph_tests, p = 3, n = 10, K_max = 1e3)

## Monte Carlo distribution

# Circular statistics
x <- seq(0, 5, l = 10)
unif_stat_distr(x = x, type = avail_cir_tests, p = 2, n = 10, approx = "MC")
unif_stat_distr(x = x, type = "Kuiper", p = 2, n = 10, approx = "MC")
unif_stat_distr(x = x, type = c("Ajne", "Kuiper"), p = 2, n = 10,
                approx = "MC")

# Spherical statistics
unif_stat_distr(x = x, type = avail_sph_tests, p = 3, n = 10,
                approx = "MC")
unif_stat_distr(x = cbind(x, x + 1), type = c("Rayleigh", "Bingham"),
                p = 3, n = 10, approx = "MC")
unif_stat_distr(x = cbind(x, x + 1), type = c("Rayleigh", "Bingham"),
                p = 3, n = 10, approx = "MC")

## Specific arguments

# Rothman
unif_stat_distr(x = x, type = "Rothman", p = 2, n = 10, Rothman_t = 0.5,
                approx = "MC")

# CCF09
dirs <- r_unif_sph(n = 5, p = 3, M = 1)[, , 1]
x <- seq(0, 1, l = 10)
unif_stat_distr(x = x, type = "CCF09", p = 3, n = 10, approx = "MC",
                CCF09_dirs = dirs)
unif_stat_distr(x = x, type = "CCF09", p = 3, n = 10, approx = "MC")

# CJ12
unif_stat_distr(x = x, type = "CJ12", p = 3, n = 100, CJ12_reg = 3)
unif_stat_distr(x = x, type = "CJ12", p = 3, n = 100, CJ12_reg = 2,
               CJ12_beta = 0.01)
unif_stat_distr(x = x, type = "CJ12", p = 3, n = 100, CJ12_reg = 1)

## Sobolev

x <- seq(0, 1, l = 10)
vk2 <- diag(1, nrow = 3)
unif_stat_distr(x = x, type = "Sobolev", approx = "asymp", p = 3, n = 100,
                Sobolev_vk2 = vk2)
sapply(1:3, function(i)
  unif_stat_distr(x = x, type = "Sobolev", approx = "asymp", p = 3, n = 100,
                  Sobolev_vk2 = vk2[i, ])$Sobolev)
sapply(1:3, function(i)
  unif_stat_distr(x = x, type = "Sobolev", approx = "MC", p = 3, n = 100,
                  Sobolev_vk2 = vk2[i, ], M = 1e3)$Sobolev)
unif_stat_distr(x = x, type = "Sobolev", approx = "MC", p = 3, n = 100,
                Sobolev_vk2 = vk2, M = 1e3)
}
}
\references{
García-Portugués, E. and Verdebout, T. (2018) An overview of uniformity
tests on the hypersphere. \emph{arXiv:1804.00286}.
\doi{10.48550/arXiv.1804.00286}.
}
