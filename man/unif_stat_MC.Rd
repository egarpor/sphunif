% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unif_stat_MC.R
\name{unif_stat_MC}
\alias{unif_stat_MC}
\title{Monte Carlo simulation of circular and (hyper)spherical uniformity
statistics}
\usage{
unif_stat_MC(n, type = "all", p, M = 10000, r_H1 = NULL,
  crit_val = NULL, alpha = c(0.1, 0.05, 0.01), return_stats = TRUE,
  stats_sorted = FALSE, chunks = ceiling((n * M)/1e+05), cores = 1,
  seeds = NULL, CCF09_dirs = NULL, CJ12_reg = 3, cov_a = 2 * pi,
  Cressie_t = 1/3, K_CCF09 = 25, Poisson_rho = 0.5, Pycke_q = 0.5,
  Rayleigh_m = 1, Riesz_s = 1, Rothman_t = 1/3, Sobolev_vk2 = c(0, 0,
  1), Softmax_kappa = 1, Stein_K = 10, Stein_cf = FALSE, Stereo_a = 0,
  ...)
}
\arguments{
\item{n}{sample size.}

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

\item{M}{number of Monte Carlo replications. Defaults to \code{1e4}.}

\item{r_H1}{if provided, the computation of empirical powers is
carried out for the alternative hypothesis sampled with \code{r_H1}.
This must be a function with the same arguments and value as
\code{\link{r_unif_sph}} (see examples). Defaults to \code{NULL}, indicating
that the critical values are estimated from samples of \code{r_unif_sph}.}

\item{crit_val}{if provided, must be the critical values as returned by
\code{$stats_MC} in a call to \code{unif_stat_MC}. They are used for
computing the empirical powers of the tests present in \code{type}.
Defaults to \code{NULL}, which means that no power computation is done.}

\item{alpha}{vector with significance levels. Defaults to
\code{c(0.10, 0.05, 0.01)}.}

\item{return_stats}{return the Monte Carlo statistics? If only the critical
values or powers are desired, \code{FALSE} saves memory in the returned
object. Defaults to \code{TRUE}.}

\item{stats_sorted}{sort the returned Monte Carlo statistics? If
\code{TRUE}, this is useful for evaluating faster the empirical cumulative
distribution function when approximating the distribution in
\code{\link{unif_stat_distr}}. Defaults to \code{FALSE}.}

\item{chunks}{number of chunks to split the \code{M} Monte Carlo
replications. Useful for parallelizing the simulation study in \code{chunks}
tasks containing \code{ceiling(M / chunks)} replications. Useful also for
avoiding memory bottlenecks when \code{M} is large. Defaults to
\cr\code{ceiling((n * M) / 1e5)}.}

\item{cores}{number of cores to perform the simulation. Defaults to \code{1}.}

\item{seeds}{if provided, a vector of size \code{chunks} for fixing the
seeds on each of the simulation chunks (useful for reproducing parallel
simulations). Specifically, for \code{k in 1:chunks}, seeds are
set as \code{set.seed(seeds[k], kind = "Mersenne-Twister")} in each chunk.
Defaults to \code{NULL} (no seed setting is done).}

\item{CCF09_dirs}{a matrix of size \code{c(n_proj, p)} containing
\code{n_proj} random directions (in Cartesian coordinates) on \eqn{S^{p-1}}
to perform the CCF09 test. If \code{NULL} (default), a sample of size
\code{n_proj = 50} directions is computed internally.}

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

\item{...}{optional arguments to be passed to the \code{r_H1} sampler or to
\code{\link[foreach]{foreach}} (for example, \code{.export} to export global
variables or other functions to the \code{foreach} environment).}
}
\value{
A list with the following entries:
\itemize{
  \item \code{crit_val_MC}: a data frame of size
  \code{c(length(alpha), length(type))}, with column names given by
  \code{type} and rows corresponding to the significance levels \code{alpha},
  that contains the estimated critical values of the tests.
  \item \code{power_MC}: a data frame of size
  \code{c(nrow(crit_val), length(type))}, with column names given by
  \code{type} and rows corresponding to the significance levels of
  \code{crit_val}, that contains the empirical powers of the tests. \code{NA}
  if \code{crit_val = NULL}.
  \item \code{stats_MC}: a data frame of size \code{c(M, length(type))}, with
  column names given by \code{type}, that contains the Monte Carlo
  statistics.
}
}
\description{
Utility for performing Monte Carlo simulation of several
statistics for assessing uniformity on the (hyper)sphere
\eqn{S^{p-1}:=\{{\bf x}\in R^p:||{\bf x}||=1\}}{
S^{p-1}:=\{x\in R^p:||x||=1\}}, \eqn{p\ge 2}.

\code{unif_stat_MC} provides a convenient wrapper for parallel
evaluation of \code{unif_stat}, the estimation of critical values under the
null distribution, and the computation of empirical powers under the
alternative.
}
\details{
It is possible to have a progress bar if \code{unif_stat_MC} is wrapped with
\code{\link[progressr:with_progress]{progressr::with_progress}} or if
\code{progressr::handlers(global = TRUE)} is invoked (once) by the user.
See the examples below. The progress bar is updated with the number of
finished chunks.

All the tests reject for large values of the test statistic
(\code{max_gap = TRUE} is assumed for the Range test), so the critical
values for the significance levels \code{alpha} correspond to the
\code{alpha}-upper quantiles of the null distribution of the test statistic.

The Monte Carlo simulation for the CCF09 test is made conditionally
on the choice of \code{CCF09_dirs}. That is, all the Monte Carlo statistics
share the same random directions.

Except for \code{CCF09_dirs}, \code{K_CCF09}, and \code{CJ12_reg}, all the
test-specific parameters are vectorized.
}
\examples{
## Critical values

# Single statistic, specific alpha
cir <- unif_stat_MC(n = 10, M = 1e2, type = "Ajne", p = 2, alpha = 0.15)
summary(cir$stats_MC)
cir$crit_val_MC

# All circular statistics
cir <- unif_stat_MC(n = 10, M = 1e2, p = 2)
head(cir$stats_MC)
cir$crit_val_MC

# All spherical statistics
sph <- unif_stat_MC(n = 10, M = 1e2, p = 3)
head(sph$stats_MC)
sph$crit_val_MC

## Using a progress bar

# Define a progress bar
require(progress)
require(progressr)
handlers(handler_progress(
  format = paste("(:spin) [:bar] :percent Iter: :current/:total Rate:",
                 ":tick_rate iter/sec ETA: :eta Elapsed: :elapsedfull"),
  clear = FALSE))

# Call unif_stat_MC() within with_progress()
with_progress(unif_stat_MC(n = 10, M = 1e2, p = 3, chunks = 10))

# With several cores
with_progress(unif_stat_MC(n = 10, M = 1e2, p = 3, chunks = 10, cores = 2))

# Instead of using with_progress() each time, it is more practical to run
# handlers(global = TRUE)
# once to activate progress bars in your R session

## Power computation

# Single statistic
cir_pow <- unif_stat_MC(n = 10, M = 1e2, type = "Ajne", p = 2,
                        crit_val = cir$crit_val_MC)
cir_pow$crit_val_MC
cir_pow$power_MC

# All circular statistics
cir_pow <- unif_stat_MC(n = 10, M = 1e2, p = 2, crit_val = cir$crit_val_MC)
cir_pow$crit_val_MC
cir_pow$power_MC

# All spherical statistics
sph_pow <- unif_stat_MC(n = 10, M = 1e2, p = 3, crit_val = sph$crit_val_MC)
sph_pow$crit_val_MC
sph_pow$power_MC
\donttest{
## Custom r_H1

# Circular
r_H1 <- function(n, p, M, l = 0.05) {

  stopifnot(p == 2)
  Theta_to_X(matrix(runif(n * M, 0, (2 - l) * pi), n, M))

}
dirs <- r_unif_sph(n = 5, p = 2, M = 1)[, , 1]
cir <- unif_stat_MC(n = 50, M = 1e2, p = 2, CCF09_dirs = dirs)
cir_pow <- unif_stat_MC(n = 50, M = 1e2, p = 2, r_H1 = r_H1, l = 0.10,
                        crit_val = cir$crit_val_MC, CCF09_dirs = dirs)
cir_pow$crit_val_MC
cir_pow$power_MC

# Spherical
r_H1 <- function(n, p, M, l = 0.5) {

  samp <- array(dim = c(n, p, M))
  for (j in 1:M) {

    samp[, , j] <- mvtnorm::rmvnorm(n = n, mean = c(l, rep(0, p - 1)),
                                    sigma = diag(rep(1, p)))
    samp[, , j] <- samp[, , j] / sqrt(rowSums(samp[, , j]^2))

  }
  return(samp)

}
dirs <- r_unif_sph(n = 5, p = 3, M = 1)[, , 1]
sph <- unif_stat_MC(n = 50, M = 1e2, p = 3, CCF09_dirs = dirs)
sph_pow <- unif_stat_MC(n = 50, M = 1e2, p = 3, r_H1 = r_H1, l = 0.5,
                       crit_val = sph$crit_val_MC, CCF09_dirs = dirs)
sph_pow$power_MC

## Pre-built r_H1

# Circular
dirs <- r_unif_sph(n = 5, p = 2, M = 1)[, , 1]
cir_pow <- unif_stat_MC(n = 50, M = 1e2, p = 2, r_H1 = r_alt, alt = "vMF",
                        kappa = 1, crit_val = cir$crit_val_MC,
                        CCF09_dirs = dirs)
cir_pow$power_MC

# Spherical
dirs <- r_unif_sph(n = 5, p = 3, M = 1)[, , 1]
sph_pow <- unif_stat_MC(n = 50, M = 1e2, p = 3, r_H1 = r_alt, alt = "vMF",
                        kappa = 1, crit_val = sph$crit_val_MC,
                        CCF09_dirs = dirs)
sph_pow$power_MC
}
}
