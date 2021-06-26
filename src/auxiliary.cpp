
// // [[Rcpp::depends(BH)]]
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::plugins("cpp11")]]
#include <RcppArmadillo.h>
#include <Rmath.h>
#include <R.h>
// #include <boost/math/special_functions/hypergeometric_pFq.hpp>
// #include <boost/math/special_functions/expm1.hpp>
// #include <boost/math/special_functions/gamma.hpp>
// #include <boost/multiprecision/mpfr.hpp>
// #include <boost/chrono/include.hpp>
// #include <boost/chrono/chrono.hpp>
// #include <boost/chrono/chrono_io.hpp>
// #include <boost/chrono/duration.hpp>
// #include <boost/chrono/time_point.hpp>
// #include <boost/chrono/system_clocks.hpp>
// typedef boost::multiprecision::mpfr_float mp_type;
using namespace Rcpp;

// Constants
const double inv_two_M_PI = 0.5 / M_PI;
const double two_M_PI = 2.0 * M_PI;

//' @title Transforming between polar and Cartesian coordinates
//'
//' @description Transformation between a matrix \code{Theta} containing
//' \code{M} circular samples of size \code{n} on \eqn{[0, 2\pi)} and an array
//' \code{X} containing the associated Cartesian coordinates on
//' \eqn{S^1:=\{{\bf x}\in R^2:||{\bf x}||=1\}}{S^1:=\{x\in R^2:||x||=1\}}.
//'
//' @inheritParams cir_stat
//' @param X an \bold{array} of size \code{c(n, 2, M)} containing the Cartesian
//' coordinates of \code{M} samples of size \code{n} of directions on
//' \eqn{S^{1}}. Must not contain \code{NA}'s.
//' @return
//' \itemize{
//'   \item \code{Theta_to_X}: the corresponding \code{X}.
//'   \item \code{X_to_Theta}: the corresponding \code{Theta}.
//' }
//' @examples
//' # Sample
//' Theta <- r_unif_cir(n = 10, M = 2)
//' X <- r_unif_sph(n = 10, p = 2, M = 2)
//'
//' # Check equality
//' sum(abs(X - Theta_to_X(X_to_Theta(X))))
//' sum(abs(Theta - X_to_Theta(Theta_to_X(Theta))))
//' @name cir_coord_conv


//' @rdname cir_coord_conv
//' @export
// [[Rcpp::export]]
arma::cube Theta_to_X(arma::mat Theta) {

  // Create X
  arma::cube X(Theta.n_rows, 2, Theta.n_cols);

  // Fill X by columns
  X(arma::span::all, arma::span(0), arma::span::all) = arma::cos(Theta);
  X(arma::span::all, arma::span(1), arma::span::all) = arma::sin(Theta);

  return X;

}


//' @rdname cir_coord_conv
//' @export
// [[Rcpp::export]]
arma::mat X_to_Theta(arma::cube X) {

  // Check dimension
  if (X.n_cols != 2) {

    stop("The number of columns in X must be 2");

  }

  // Create Theta
  arma::mat Theta = arma::atan2(
    X(arma::span::all, arma::span(1), arma::span::all),
    X(arma::span::all, arma::span(0), arma::span::all)
    );

  // Convert to [0, 2 * pi)
  Theta -= two_M_PI * arma::floor(Theta * inv_two_M_PI);
  return Theta;

}


//' @title Circular gaps
//'
//' @description Computation of the circular gaps of an angular sample
//' \eqn{\Theta_1,\ldots,\Theta_n} on \eqn{[0, 2\pi)}, defined as
//' \deqn{\Theta_{(2)} - \Theta_{(1)},\ldots,\Theta_{(n)} - \Theta_{(n - 1)},
//' 2\pi - \Theta_{(n)} - \Theta_{(1)},}
//' where
//' \deqn{0 \le \Theta_{(1)} \le \Theta_{(2)} \le \ldots \le
//' \Theta_{(n)} \le 2\pi.}
//'
//' @inheritParams cir_stat
//' @param sorted are the columns of \code{Theta} sorted increasingly? If
//' \code{TRUE}, performance is improved. If \code{FALSE} (default), each
//' column of \code{Theta} is sorted internally.
//' @return A matrix of size \code{c(n, M)} containing the \code{n} circular
//' gaps for each of the \code{M} circular samples.
//' @section Warning:
//' Be careful on avoiding the next bad usages of \code{cir_gaps}, which will
//' produce spurious results:
//' \itemize{
//'   \item The entries of \code{Theta} are \emph{not} in \eqn{[0, 2\pi)}.
//'   \item \code{Theta} is \emph{not} sorted increasingly when
//'   \code{data_sorted = TRUE}.
//' }
//' @examples
//' Theta <- cbind(c(pi, 0, 3 * pi / 2), c(0, 3 * pi / 2, pi), c(5, 3, 1))
//' cir_gaps(Theta)
//' @export
// [[Rcpp::export]]
arma::mat cir_gaps(arma::mat Theta, bool sorted = false) {

  // Sample size
  arma::uword n = Theta.n_rows;

  // Gaps
  arma::mat Di = arma::zeros(n, Theta.n_cols);

  // Sort data on each column and scale
  if (!sorted) {

    Theta = arma::sort(Theta);

  }

  // Gaps for each column
  Di.head_rows(n - 1) = arma::diff(Theta);
  Di.row(n - 1) = two_M_PI - (Theta.row(n - 1) - Theta.row(0));
  return Di;

}


//' @title Efficient evaluation of the empirical cumulative distribution
//' function
//'
//' @description Evaluates the empirical cumulative distribution function
//' (ecdf) of a sample \code{data} at the evaluation points \code{sorted_x}.
//' This is done through binary search.
//'
//' @param data a vector or column matrix containing the sample.
//' @param sorted_x a vector or column matrix with the evaluation points
//' \bold{sorted increasingly}.
//' @param data_sorted is \code{data} is already sorted increasingly?
//' This avoids sorting the data internally.
//' @param efic use the more efficient version of the ecdf evaluation? Set to
//' \code{FALSE} only for debugging purposes.
//' @param divide_n if \code{FALSE}, returns the absolute frequencies instead
//' of the relative frequencies. Defaults to \code{TRUE}.
//' @return The ecdf evaluated at \code{sorted_x}.
//' @author Original code from Douglas Bates'
//' \url{https://github.com/dmbates/ecdfExample}. Minor adaptations by Eduardo
//' García-Portugués.
//' @section Warning:
//' Be careful on avoiding the next bad usages of the function, which will
//' produce spurious results:
//' \itemize{
//'   \item \code{sorted_x} is not sorted increasingly.
//'   \item \code{data} is not sorted increasingly when
//'   \code{data_sorted = TRUE}-
//' }
//' @examples
//' set.seed(1234567)
//' samp <- rnorm(200)
//' x <- seq(-1, 1, l = 10)
//' ecdf(samp)(x)
//' sphunif:::ecdf_bin(samp, x, efic = FALSE)
//' sphunif:::ecdf_bin(samp, x, efic = TRUE)
//' @keywords internal
// [[Rcpp::export]]
arma::vec ecdf_bin(arma::vec data, arma::vec sorted_x, bool data_sorted = false,
                   bool efic = true, bool divide_n = true) {

  // Sizes
  arma::uword nx = sorted_x.n_elem;
  arma::uword n = data.n_elem;

  // Sort the sample
  if (!data_sorted) {

    data = arma::sort(data);

  }

  // Do binary search
  arma::vec res = arma::zeros(nx);
  if (efic) {

    // Code borrowed from Douglas Bates'
    // https://github.com/dmbates/ecdfExample#using-a-sorted-sample-in-c
    for (arma::uword i = 0, j = 0; i < nx; ++i) {

      while (j < n && data(j) <= sorted_x(i)) ++j;
      res(i) = j; // j is the 1-based index of the lower bound

    }

  } else {

    // Code borrowed from Douglas Bates'
    // https://github.com/dmbates/ecdfExample#using-stdlower_bound-in-c
    for (arma::uword i = 0; i < nx; ++i) {

      res(i) = std::lower_bound(data.begin(), data.end(), sorted_x(i)) -
        data.begin();

    }

  }

  // Divide by n the indexes
  if (divide_n) {

    res /= n;

  }
  return res;

}


//' @title The incomplete beta function and its inverse
//'
//' @description Computes the incomplete beta function
//' \deqn{I_x(a,b):=\int_0^x u^{a-1}(1-u)^{b-1}\,d\mathrm{u},\quad a,b>0}{
//' I_x(a,b):=\int_0^x u^{a-1}(1-u)^{b-1}du, a,b>0}
//' and its inverse function.
//'
//' @inheritParams cir_stat_distr
//' @param u a vector of probabilities of size \code{nu} or a matrix of size
//' \code{c(nu, 1)}.
//' @param a,b scalars giving the parameters of the beta function.
//' @param lower_tail accumulate the probability from the lower tail? If
//' \code{FALSE}, the probability is accumulated from the \emph{upper} tail.
//' Defaults to \code{FALSE}.
//' @param log use log-scale? If \code{TRUE}, returns the logarithm of the
//' incomplete beta function and uses log-scale for \code{u} in \code{beta_inc}.
//' Defaults to \code{FALSE}.
//' @return
//' \itemize{
//'   \item \code{beta_inc}: a matrix of size \code{c(nx, 1)} with the
//'   evaluation of the incomplete beta function at \code{x}.
//'   \item \code{beta_inc_inv}: a matrix of size \code{c(nu, 1)} with the
//'   evaluation of the inverse incomplete beta function at \code{u}.
//' }
//' @details
//' The functions are mere wrappers to R's internal \code{pbeta} and
//' \code{qbeta} functions.
//' @examples
//' # Comparison with R
//' old_par <- par(mfrow = c(1, 2))
//' x <- seq(-1, 2, l = 1e3)
//' plot(x, sphunif:::beta_inc(x, 0.75, 2), type = "l")
//' lines(x, pbeta(x, 0.75, 2), col = 2)
//' u <- seq(0, 1, l = 1e3)
//' plot(u, sphunif:::beta_inc_inv(u, 0.75, 2), type = "l")
//' lines(u, qbeta(u, 0.75, 2), col = 2)
//' par(old_par)
//' @keywords internal
// [[Rcpp::export]]
arma::vec beta_inc(arma::vec x, double a, double b, bool lower_tail = true,
                   bool log = false) {

  // Vectorize by a lambda function -- requires C++ 11
  x.transform([a, b, lower_tail, log](double y) {
    return R::pbeta(y, a, b, lower_tail, log);
  });
  return x;

}


//' @rdname beta_inc
// [[Rcpp::export]]
arma::vec beta_inc_inv(arma::vec u, double a, double b, bool lower_tail = true,
                       bool log = false) {

  // Vectorize by a lambda function -- requires C++ 11
  u.transform([a, b, lower_tail, log](double y) {
    return R::qbeta(y, a, b, lower_tail, log);
  });
  return u;

}


// //' @title The generalized hypergeometric function
// //'
// //' @description Computes the generalized hypergeometric function
// //' \deqn{{}_pF_q(\mathbf{a}; \mathbf{b}; 1)}{pFq(a; b; 1)}.
// //'
// //' @param a,b vectors of parameters for \eqn{{}_pF_q(\mathbf{a}; \mathbf{b}; 1)
// //' }{pFq(a; b; 1)}.
// //' @return TODO
// //' @details
// //' The function is a mere wrapper to \href{boost}{https://www.boost.org/}'s
// //' \code{hypergeometric_pFq} function (documented
// //' \href{https://www.boost.org/doc/libs/1_72_0/libs/math/doc/html/math_toolkit/hypergeometric/hypergeometric_pfq.html}{
// //' here}).
// //' @references
// //' Maddock, J., Bristow, P., Holin, H., and Zhang, X. \emph{Math/Special
// //' Functions} in Boost 1.72.0 Library Documentation.
// //' \url{https://www.boost.org/doc/libs/1_72_0/libs/math/doc/html/special.html}
// //' @examples
// //' # 2F1(c(1, 2), c(4)) = gamma(4) * gamma(1) / (gamma(3) * gamma(2)) = 3
// //' pFq(c(1, 2), c(4))
// //'
// //' # 4F3(c(1, 2, 5, 5), c(4, 5, 5); 1) = 3
// //' hypergeo::genhypergeo(U = c(1, 2, 5, 5), L = c(4, 5, 5), z = 0.99)
// //' pFq(c(1, 2, 5, 5), c(4, 5, 5))
// //' @keywords internal
// // [[Rcpp::export]]
// double pFq(double a, double b, double c, double d,
//            double e, double f, double g) {
//
//   mp_type result = boost::math::hypergeometric_pFq_precision(
//     {mp_type(a), mp_type(b), mp_type(c), mp_type(d)},
//     {mp_type(e), mp_type(f), mp_type(g)}, mp_type(1.0), 10, 1.0);
//   double d_result = static_cast<double>(result);
//   return d_result;
//
// }


//' @title Low-level utilities for \pkg{sphunif}
//'
//' @description Internal and undocumented low-level utilities for
//' \pkg{sphunif}.
//'
//' @param n_dist a positive integer \eqn{(n - 1) * n / 2} for which \eqn{n}
//' is to be recovered.
//' @param t a vector to evaluate \eqn{t / \sqrt{1 - t^2}}.
//' @examples
//' n <- 10
//' sphunif:::n_from_dist_vector((n - 1) * n / 2)
//' sphunif:::n_from_dist_vector((n - 1) * n / 2 + 1) # Invalid input
//' @name utils


//' @rdname utils
//' @keywords internal
// [[Rcpp::export]]
arma::uword n_from_dist_vector(arma::uword n_dist) {

  return(0.5 * (std::sqrt(8.0 * n_dist + 1) + 1));

}


//' @rdname utils
//' @keywords internal
// [[Rcpp::export]]
arma::vec t_inv_sqrt_one(arma::vec t) {

  return t / arma::sqrt(1 - arma::square(t));

}

