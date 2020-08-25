
// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <Rmath.h>
#include <R.h>
using namespace Rcpp;

// Declaration for functions
arma::vec d_chisq(arma::vec, arma::uword df);
arma::vec p_chisq(arma::vec, arma::uword df);

// Constants
const double eight_PI = 8.0 * PI;
const double two_PI = 2.0 * PI;
const double inv_sqrt_eight_PI = 1.0 / std::sqrt(eight_PI);
const double inv_sqrt_two_PI = 1.0 / std::sqrt(two_PI);

//' @rdname sph_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec p_sph_stat_Rayleigh(arma::vec x, arma::uword p) {

  return p_chisq(x, p);

}


//' @rdname sph_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec d_sph_stat_Rayleigh(arma::vec x, arma::uword p) {

  return d_chisq(x, p);

}


//' @rdname sph_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec p_sph_stat_Bingham(arma::vec x, arma::uword p) {

  return p_chisq(x, 0.5 * (p - 1) * (p + 2));

}


//' @rdname sph_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec d_sph_stat_Bingham(arma::vec x, arma::uword p) {

  return d_chisq(x, 0.5 * (p - 1) * (p + 2));

}


//' @rdname sph_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec p_sph_stat_Rayleigh_HD(arma::vec x, arma::uword p) {

  return arma::normcdf(x, 0, 1);

}


//' @rdname sph_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec d_sph_stat_Rayleigh_HD(arma::vec x, arma::uword p) {

  return arma::normpdf(x, 0, 1);

}


//' @rdname sph_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec p_sph_stat_Cai(arma::vec x, arma::uword regime = 1, double beta = 0) {

  // Avoid degeneracy
  if (beta == 0 && regime == 2) {

    regime = 1;

  }

  // Super-exponential regime
  double fact = 1;
  if (regime == 3) {

    fact *= inv_sqrt_two_PI;

  // Exponential regime
  } else if (regime == 2) {

    fact = std::sqrt(beta / (two_PI * (1 - std::exp(-4 * beta))));

  // Sub-exponential regime
  } else {

    fact *= inv_sqrt_eight_PI;

  }

  // Distribution
  return 1 - arma::exp(-fact * arma::exp(0.5 * (x + 8 * beta)));

}


//' @rdname sph_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec d_sph_stat_Cai(arma::vec x, arma::uword regime = 3, double beta = 0) {

  // Avoid degeneracy
  if (beta == 0 && regime == 2) {

    regime = 1;

  }

  // Super-exponential regime
  double fact = 1;
  if (regime == 3) {

    fact *= inv_sqrt_two_PI;

  // Exponential regime
  } else if (regime == 2) {

    fact = std::sqrt(beta / (2 * PI * (1 - std::exp(-4 * beta))));

  // Sub-exponential regime
  } else {

    fact *= inv_sqrt_eight_PI;

  }

  // Density
  x = 0.5 * (x + 8 * beta);
  return 0.5 * fact * arma::exp(-fact * arma::exp(x) + x);

}

