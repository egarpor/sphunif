
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::plugins("cpp11")]]
#include <RcppArmadillo.h>
#include <Rmath.h>
#include <R.h>
using namespace Rcpp;

// Declaration for functions
arma::vec ecdf_bin(arma::vec data, arma::vec sorted_x, bool data_sorted,
                   bool efic, bool divide_n);
arma::vec beta_inc(arma::vec x, double a, double b, bool lower_tail,
                   bool log);
arma::vec beta_inc_inv(arma::vec u, double a, double b, bool lower_tail,
                       bool log);

// Constants
const double log_two = std::log(2.0);
const double two_M_PI = 2.0 * M_PI;
const double inv_M_PI = 1.0 / M_PI;


//' @title Projection of the spherical uniform distribution
//'
//' @description Density, distribution, and quantile functions of the
//' projection of the spherical uniform random variable on an arbitrary
//' direction, that is, the random variable
//' \eqn{\boldsymbol{\gamma}'{\bf X}}{\gamma'X}, where \eqn{{\bf X}}{X}
//' is uniformly distributed on the (hyper)sphere
//' \eqn{S^{p-1}:=\{{\bf x}\in R^p:||{\bf x}||=1\}}{S^{p-1}:=
//' \{x\in R^p:||x||=1\}}, \eqn{p\ge 2}, and
//' \eqn{\boldsymbol{\gamma}\in S^{p-1}}{\gamma\in S^{p-1}} is an
//' \emph{arbitrary} projection direction. Note that the distribution is
//' invariant to the choice of \eqn{\boldsymbol{\gamma}}{\gamma}. Also,
//' efficient simulation of \eqn{\boldsymbol{\gamma}'{\bf X}}{\gamma'X}.
//'
//' @inheritParams cir_stat_distr
//' @inheritParams r_unif
//' @param u vector of probabilities.
//' @param log compute the logarithm of the density or distribution?
//' @return A matrix of size \code{c(nx, 1)} with the evaluation of the
//' density, distribution, or quantile function at \code{x} or \code{u}.
//' For \code{r_proj_unif}, a random vector of size \code{n}.
//' @author Eduardo García-Portugués and Paula Navarro-Esteban.
//' @examples
//' # Density function
//' curve(d_proj_unif(x, p = 2), from = -2, to = 2, n = 2e2, ylim = c(0, 2))
//' curve(d_proj_unif(x, p = 3), n = 2e2, col = 2, add = TRUE)
//' curve(d_proj_unif(x, p = 4), n = 2e2, col = 3, add = TRUE)
//' curve(d_proj_unif(x, p = 5), n = 2e2, col = 4, add = TRUE)
//' curve(d_proj_unif(x, p = 6), n = 2e2, col = 5, add = TRUE)
//'
//' # Distribution function
//' curve(p_proj_unif(x, p = 2), from = -2, to = 2, n = 2e2, ylim = c(0, 1))
//' curve(p_proj_unif(x, p = 3), n = 2e2, col = 2, add = TRUE)
//' curve(p_proj_unif(x, p = 4), n = 2e2, col = 3, add = TRUE)
//' curve(p_proj_unif(x, p = 5), n = 2e2, col = 4, add = TRUE)
//' curve(p_proj_unif(x, p = 6), n = 2e2, col = 5, add = TRUE)
//'
//' # Quantile function
//' curve(q_proj_unif(u = x, p = 2), from = 0, to = 1, n = 2e2, ylim = c(-1, 1))
//' curve(q_proj_unif(u = x, p = 3), n = 2e2, col = 2, add = TRUE)
//' curve(q_proj_unif(u = x, p = 4), n = 2e2, col = 3, add = TRUE)
//' curve(q_proj_unif(u = x, p = 5), n = 2e2, col = 4, add = TRUE)
//' curve(q_proj_unif(u = x, p = 6), n = 2e2, col = 5, add = TRUE)
//'
//' # Sampling
//' hist(r_proj_unif(n = 1e4, p = 4), freq = FALSE, breaks = 50)
//' curve(d_proj_unif(x, p = 4), n = 2e2, col = 3, add = TRUE)
//' @name proj_unif


//' @rdname proj_unif
//' @export
// [[Rcpp::export]]
arma::vec d_proj_unif(arma::vec x, arma::uword p, bool log = false) {

  // Check p
  if (p <= 1) {

    stop("p must be >= 2.");

  }

  // Set to log(0) outside (-1, 1)
  arma::vec f0 = arma::vec(x.n_elem).fill(-arma::datum::inf);

  // f0 inside (-1, 1)
  arma::uvec ind = arma::find(x > -1 && x < 1);
  f0.elem(ind) = (0.5 * p - 1.5) * arma::log1p(-arma::square(x.elem(ind))) -
    R::lbeta(0.5, 0.5 * (p - 1));

  // Return log-density or density
  if (!log) {

    f0 = arma::exp(f0);

  }
  return f0;

}


//' @rdname proj_unif
//' @export
// [[Rcpp::export]]
arma::vec p_proj_unif(arma::vec x, arma::uword p, bool log = false) {

  // Check p
  if (p <= 1) {

    stop("p must be >= 2.");

  }

  // Set to log(0) for (-\infty, -1]
  arma::vec F0 = arma::vec(x.n_elem).fill(-arma::datum::inf);

  // Special cases
  if (p == 2 || p == 3 || p == 4 || p == 5) {

    // F0 inside (-1, 1)
    arma::uvec ind = arma::find(x > -1 && x < 1);

    // Cases
    if (p == 2) {

      F0.elem(ind) = arma::log1p(-arma::acos(x.elem(ind)) * inv_M_PI);

    } else if (p == 3) {

      F0.elem(ind) = arma::log1p(x.elem(ind)) - log_two;

    } else if (p == 4) {

      F0.elem(ind) = arma::log1p((x.elem(ind) %
        arma::sqrt(1 - arma::square(x.elem(ind))) -
        arma::acos(x.elem(ind))) * inv_M_PI);

    } else if (p == 5) {

      F0.elem(ind) = arma::log1p(x.elem(ind)) +
        arma::log1p(0.5 * x.elem(ind) % (1 - x.elem(ind))) - log_two;

    }

    // Set to log(1) for [1, +\infty)
    ind = arma::find(x >= 1);
    F0.elem(ind).zeros();

  // General case
  } else {

    F0 = -log_two + arma::log1p(arma::sign(x) %
      beta_inc(arma::square(x), 0.5, 0.5 * (p - 1), true, false));

  }

  // Return log-probability or probability
  if (!log) {

    F0 = arma::exp(F0);

  }
  return F0;

}


//' @rdname proj_unif
//' @export
// [[Rcpp::export]]
arma::vec q_proj_unif(arma::vec u, arma::uword p) {

  // Check p
  if (p <= 1) {

    stop("p must be >= 2.");

  }

  // Special cases
  if (p == 2 || p == 3) {

    // Set to NaN the quantiles of unproper probabilities
    arma::vec x = arma::vec(u.n_elem).fill(arma::datum::nan);

    // u inside (-1, 1)
    arma::uvec ind = arma::find(u >= 0 && u <= 1);

    // Cases
    if (p == 2) {

      x.elem(ind) = arma::cos(M_PI * (1 - u.elem(ind)));

    } else if (p == 3) {

      x.elem(ind) = 2 * u.elem(ind) - 1;

    }

    // Return quantile
    return x;

  // General case
  } else {

    u = 2 * u - 1;
    return arma::sign(u) %
      arma::sqrt(beta_inc_inv(arma::abs(u), 0.5, 0.5 * (p - 1), true, false));

  }

}


//' @title Sample uniformly distributed circular and spherical data
//'
//' @description Simulation of the uniform distribution on \eqn{[0, 2\pi)} and
//' \eqn{S^{p-1}:=\{{\bf x}\in R^p:||{\bf x}||=1\}}{
//' S^{p-1}:=\{x\in R^p:||x||=1\}}, \eqn{p\ge 2}.
//'
//' @param n sample size.
//' @param M number of samples of size \code{n}. Defaults to \code{1}.
//' @param p integer giving the dimension of the ambient space \eqn{R^p} that
//' contains \eqn{S^{p-1}}.
//' @param sorted return each circular sample sorted? Defaults to \code{FALSE}.
//' @return
//' \itemize{
//'   \item \code{r_unif_cir}: a \bold{matrix} of size \code{c(n, M)} with
//'   \code{M} random samples of size \code{n} of uniformly-generated circular
//'   data on \eqn{[0, 2\pi)}.
//'   \item \code{r_unif_sph}: an \bold{array} of size \code{c(n, p, M)} with
//'   \code{M} random samples of size \code{n} of uniformly-generated
//'   directions on \eqn{S^{p-1}}.
//' }
//' @examples
//' # A sample on [0, 2*pi)
//' n <- 5
//' r_unif_cir(n = n)
//'
//' # A sample on S^1
//' p <- 2
//' samp <- r_unif_sph(n = n, p = p)
//' samp
//' rowSums(samp^2)
//'
//' # A sample on S^2
//' p <- 3
//' samp <- r_unif_sph(n = n, p = p)
//' samp
//' rowSums(samp^2)
//' @name r_unif


//' @rdname r_unif
//' @export
// [[Rcpp::export]]
arma::mat r_unif_cir(arma::uword n, arma::uword M = 1, bool sorted = false) {

  arma::mat sample = arma::randu(n, M) * two_M_PI;
  if (sorted) {

    sample = arma::sort(sample);

  }
  return sample;

}


//' @rdname r_unif
//' @export
// [[Rcpp::export]]
arma::cube r_unif_sph(arma::uword n, arma::uword p, arma::uword M = 1) {

  arma::cube sample = arma::randn(n, p, M);
  for (arma::uword k = 0; k < M; k++) {

    sample.slice(k) = arma::normalise(sample.slice(k), 2, 1);

  }
  return sample;

}


//' @title Utilities for weighted sums of non-central chi squared random
//' variables
//'
//' @description Simulation from a weighted sum of non-central chi squared
//' random variables and Monte Carlo approximation of its distribution function.
//'
//' @inheritParams r_unif
//' @inheritParams wschisq
//' @param ncps non-negative non-centrality parameters. A vector with the same
//' length as \code{weights}.
//' @param M number of Monte Carlo samples for approximating the distribution.
//' Defaults to \code{1e4}.
//' @param sample if \code{use_sample = TRUE}, the Monte Carlo sample to
//' approximate the distribution. If not, it is computed internally. Defaults
//' to \code{1e4}.
//' @param use_sample use the already computed \code{sample}? If \code{FALSE}
//' (default), \code{sample} is computed internally.
//' @return
//' \itemize{
//'   \item \code{r_wschisq_Cpp}: a matrix of size \code{c(n, 1)}
//'   containing a random sample.
//'   \item \code{p_wschisq_MC}: a matrix of size \code{c(nx, 1)}
//'   with the evaluation of the distribution function at \code{x}.
//' }
//' @examples
//' x <- seq(0, 50, l = 1e3)
//' weights <- c(2, 1, 0.5)
//' dfs <- c(3, 6, 12)
//' ncps <- c(0, 0, 1)
//' samp <- sphunif:::r_wschisq_Cpp(n = 5e2, weights = weights, dfs = dfs,
//'                                 ncps = ncps)
//' plot(ecdf(samp), main = "")
//' lines(x, sphunif:::p_wschisq_MC(x, weights = weights, dfs = dfs,
//'                                 ncps = ncps),
//'       type = "s", col = 2)
//' @name wschisq_utils


//' @rdname wschisq_utils
//' @keywords internal
// [[Rcpp::export]]
arma::vec r_wschisq_Cpp(arma::uword n, arma::vec weights, arma::vec dfs,
                        arma::vec ncps) {

  // Sample and number of elements in the sum
  arma::vec samp = arma::zeros(n);
  arma::uword p = weights.n_elem;

  // Loop on the sum
  for (arma::uword k = 0; k < p; k++) {

    // Loop on sample elements
    arma::vec samp_k = arma::zeros(n);
    for (arma::uword i = 0; i < n; i++) {

      samp_k(i) = R::rnchisq(dfs(k), ncps(k));

    }

    // Weighted sum
    samp += weights(k) * samp_k;

  }

  return(samp);

}


//' @rdname wschisq_utils
//' @keywords internal
// [[Rcpp::export]]
arma::vec p_wschisq_MC(arma::vec x, arma::vec weights = 0, arma::vec dfs = 0,
                       arma::vec ncps = 0, arma::uword M = 1e4,
                       arma::vec sample = 0, bool use_sample = false,
                       bool x_sorted = false) {

  // Sample
  if (!use_sample) {

    sample.set_size(M);
    sample = r_wschisq_Cpp(M, weights, dfs, ncps);

  }

  // Estimated cdf
  if (x_sorted || (x.n_elem == 1)) {

    // Call ecdf_bin directly
    x = ecdf_bin(sample, x, false, true, true);

  } else {

    // Sort x for calling ecdf_bin and unsort afterwards
    arma::uvec ind = arma::sort_index(x);
    x = ecdf_bin(sample, x(ind), false, true, true);
    x = x(arma::sort_index(ind));

  }

  return x;

}


//' @title Density and distribution of a chi squared
//'
//' @description Computation of the density and distribution functions of a chi
//' squared.
//'
//' @inheritParams cir_stat_distr
//' @param df degrees of freedom.
//' @param ncp non-centrality parameter.
//' @return A matrix of size \code{c(nx, 1)} with the evaluation of the
//' density or distribution function at \code{x}.
//' @examples
//' curve(sphunif:::d_chisq(x, df = 2), from = 0, to = 5, n = 200)
//' curve(sphunif:::p_chisq(x, df = 2), add = TRUE, col = 2)
//' @name chisq


//' @rdname chisq
//' @keywords internal
// [[Rcpp::export]]
arma::vec d_chisq(arma::vec x, arma::uword df, arma::uword ncp = 0) {

  // Vectorize by a lambda function -- requires C++ 11
  x.transform([df, ncp](double y) {
    return R::dnchisq(y, df, ncp, false);
  });
  return x;

}


//' @rdname chisq
//' @keywords internal
// [[Rcpp::export]]
arma::vec p_chisq(arma::vec x, arma::uword df, arma::uword ncp = 0) {

  // Vectorize by a lambda function -- requires C++ 11
  x.transform([df, ncp](double y) {
    return R::pnchisq(y, df, ncp, true, false);
  });
  return x;

}

