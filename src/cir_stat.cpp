
// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// Declaration for functions
arma::mat cir_gaps(arma::mat Theta, bool sorted);
arma::cube Theta_to_X(arma::mat Theta);
arma::mat Psi_mat(arma::cube data, arma::uvec ind_tri, bool use_ind_tri,
                  bool scalar_prod, bool angles_diff);
arma::uvec upper_tri_ind(arma::uword n);
arma::uword n_from_dist_vector(arma::uword n_dist);
arma::vec ecdf_bin(arma::vec data, arma::vec sorted_x, bool data_sorted,
                   bool efic, bool divide_n);
arma::vec cir_stat_An_Psi(arma::mat Psi, arma::uword n);
arma::vec cir_stat_Rothman_Psi(arma::mat Psi, double t_m2, double t_min2,
                               arma::uword n);
arma::vec cir_stat_Cressie(arma::mat Theta, double t, bool sorted);
arma::vec sph_stat_Bingham(arma::cube X);
arma::vec cir_stat_Hermans_Rasson_Psi(arma::mat Psi, arma::uword n);
arma::vec sph_stat_Gine_Gn(arma::cube X, bool Psi_in_X, arma::uword p);
arma::vec sph_stat_Gine_Fn(arma::cube X, bool Psi_in_X, arma::uword p);
arma::vec cir_stat_Pycke_Psi(arma::mat Psi, arma::uword n);
arma::vec cir_stat_Pycke_q_Psi(arma::mat Psi, arma::uword n, double q);
arma::vec sph_stat_Riesz(arma::cube X, bool Psi_in_X, arma::uword p, double s);
arma::vec sph_stat_PCvM(arma::cube X, bool Psi_in_X, arma::uword p,
                       arma::uword N, arma::uword L);
arma::vec sph_stat_PRt(arma::cube X, double t,  bool Psi_in_X, arma::uword p,
                       arma::uword N, arma::uword L);
arma::vec sph_stat_PAD(arma::cube X, bool Psi_in_X, arma::uword p,
                       arma::uword N, arma::uword L);
arma::vec sph_stat_Cuesta_Albertos(arma::cube X, arma::mat rand_dirs,
                                   arma::uword K_Cuesta_Albertos,
                                   bool original = false);

// Constants
const double pi = PI;
const double inv_PI = 1.0 / pi; // To force including pi, needed for defaults
const double inv_two_PI = 0.5 / PI;
const double sq_PI = PI * PI;
const double inv_two_PI_sq = 0.25 / sq_PI;
const double two_sq_PI = 2 *sq_PI;
const double two_PI = 2.0 * PI;
const double half_PI = 0.5 * PI;
const double two_inv_PI = 2.0 / PI;
const double beta2 = (sq_PI / 36) / (0.5 - 4 / sq_PI);
const double const_Hn = half_PI + beta2 * two_inv_PI;


/*
 * Tests based on the empirical cumulative distribution function
 */


//' @rdname cir_stat
//' @export
// [[Rcpp::export]]
arma::vec cir_stat_Kuiper(arma::mat Theta, bool sorted = false,
                          bool KS = false, bool Stephens = false) {

  // Sample size
  arma::uword n = Theta.n_rows;

  // Sort data on each column
  if (!sorted) {

    Theta = arma::sort(Theta);

  }

  // Scale
  Theta *= inv_two_PI;

  // Substract i / n to columns, Theta becomes Dn^+ when taking the -min on
  // each of the columns
  arma::vec i = arma::linspace(1.0 / n, 1.0, n);
  Theta.each_col() -= i;

  // Statistic for each column
  arma::vec Vn = arma::zeros(Theta.n_cols);
  if (KS) {

    // K-S statistic max(Dn^+, Dn^-) for each column
    Vn = std::sqrt(n) * arma::max(arma::join_vert(-arma::min(Theta),
                             arma::max(Theta) + 1.0 / n)).t();

    // Add Stephens (1970) modification?
    if (Stephens) {

      Vn *= (1 + 0.12 / std::sqrt(n) + 0.21 / n);

    }

  } else {

    // Kuiper statistic Dn^- + Dn^+ for each column
    Vn = std::sqrt(n) * (arma::max(Theta) + 1.0 / n - arma::min(Theta)).t();

    // Add Stephens (1970) modification?
    if (Stephens) {

      Vn *= (1 + 0.155 / std::sqrt(n) + 0.24 / n);

    }

  }
  return Vn;

}


//' @rdname cir_stat
//' @export
// [[Rcpp::export]]
arma::vec cir_stat_Watson(arma::mat Theta, bool sorted = false,
                          bool CvM = false, bool Stephens = false) {

  // Sample size
  arma::uword n = Theta.n_rows;

  // Sort data on each column
  if (!sorted) {

    Theta = arma::sort(Theta);

  }

  // Scale
  Theta *= inv_two_PI;

  // Mean
  arma::rowvec Theta_bar = arma::mean(Theta);

  // Substract (i - 0.5) / n to columns
  arma::vec i = arma::linspace(0.5 / n, 1.0 - 0.5 / n, n);
  Theta.each_col() -= i;

  // Substract mean to rows - the only difference between the Watson and CvM
  // statistics
  if (!CvM) {

    Theta.each_row() -= Theta_bar - 0.5;

  }

  // Sum of squares plus bias
  arma::vec Un2 = arma::sum(arma::square(Theta), 0).t() + 1.0 / (12.0 * n);

  // Add Stephens (1970) modification?
  if (Stephens) {

    if (CvM) {

      Un2 += -0.4 / n + 0.6 / (n * n);
      Un2 *= 1 + 1.0 / n;

    } else {

      Un2 += -0.1 / n + 0.1 / (n * n);
      Un2 *= 1 + 0.8 / n;

    }

  }

  // Statistic for each column
  return Un2;

}


//' @rdname cir_stat
//' @export
// [[Rcpp::export]]
arma::vec cir_stat_Watson_1976(arma::mat Theta, bool sorted = false,
                               bool minus = false) {

  // Sample size
  arma::uword n = Theta.n_rows;

  // Sort data on each column
  if (!sorted) {

    Theta = arma::sort(Theta);

  }

  // Scale
  Theta *= inv_two_PI;

  // Mean
  arma::rowvec Theta_bar = arma::mean(Theta);

  // Substract i / n to columns
  arma::vec i = arma::linspace(1.0 / n, 1.0, n);
  Theta.each_col() -= i;

  // Statistic for each column: invariant Dn^+ or Dn^-?
  arma::rowvec Mn = arma::zeros(1, Theta.n_cols);
  if (minus) {

    // Dn^- minus mean
    Mn = arma::max(Theta) + 1.0 / n;
    Mn -= Theta_bar - 0.5;

  } else {

    // Dn^+ plus mean
    Mn = -arma::min(Theta);
    Mn += Theta_bar - 0.5;

  }
  return std::sqrt(n) * Mn.t();

}


/*
 * Spacings-based tests
 */


//' @rdname cir_stat
//' @export
// [[Rcpp::export]]
arma::vec cir_stat_Range(arma::mat Theta, bool sorted = false,
                         bool gaps_in_Theta = false, bool max_gap = true) {

  // Are the gaps stored in Theta?
  if (!gaps_in_Theta) {

    Theta = cir_gaps(Theta, sorted);

  }

  // Statistic for each column
  arma::vec Tn = arma::max(Theta, 0).t();
  if (!max_gap) {

    Tn = two_PI - Tn;

  }
  return Tn;

}


//' @rdname cir_stat
//' @export
// [[Rcpp::export]]
arma::vec cir_stat_Rao(arma::mat Theta, bool sorted = false,
                       bool gaps_in_Theta = false) {

  // Sample size
  arma::uword n = Theta.n_rows;

  // Are the gaps stored in Theta?
  if (!gaps_in_Theta) {

    Theta = cir_gaps(Theta, sorted);

  }

  // Statistic for each column
  return std::sqrt(n) * (0.5 * arma::sum(arma::abs(Theta - two_PI / n), 0).t() -
                   two_PI / arma::datum::e);

}


//' @rdname cir_stat
//' @export
// [[Rcpp::export]]
arma::vec cir_stat_Greenwood(arma::mat Theta, bool sorted = false,
                             bool gaps_in_Theta = false) {

  // Sample size
  arma::uword n = Theta.n_rows;

  // Are the gaps stored in Theta?
  if (!gaps_in_Theta) {

    Theta = cir_gaps(Theta, sorted);

  }

  // Statistic for each column
  return std::sqrt(n) *
    (n * sum(arma::square(Theta), 0).t() * inv_two_PI_sq - 2.0);

}


//' @rdname cir_stat
//' @export
// [[Rcpp::export]]
arma::vec cir_stat_Log_gaps(arma::mat Theta, bool sorted = false,
                            bool gaps_in_Theta = false, bool abs_val = true) {

  // Sample size
  arma::uword n = Theta.n_rows;

  // Are the gaps stored in Theta?
  if (!gaps_in_Theta) {

    Theta = cir_gaps(Theta, sorted);

  }

  // Statistic for each column (minus for changing the rejection to
  // large values instead of small values)
  arma::vec In = std::sqrt(n) * (std::log(two_PI / n) -
    arma::mean(arma::log(Theta), 0).t() - arma::datum::euler);

  // Return the absolute value? Useful because rejection happens for large
  // absolute values of In
  if (abs_val) {

    In = arma::abs(In);

  }

  return In;

}


//' @rdname cir_stat
//' @export
// [[Rcpp::export]]
arma::vec cir_stat_Vacancy(arma::mat Theta, double a = 2 * pi,
                           bool sorted = false, bool gaps_in_Theta = false) {

  // Sample size
  arma::uword n = Theta.n_rows;

  // Are the gaps stored in Theta?
  if (!gaps_in_Theta) {

    Theta = cir_gaps(Theta, sorted);

  }

  // Statistic for each column
  Theta = arma::clamp(Theta - a / n, 0, arma::datum::inf);
  arma::vec Yn = arma::sum(Theta, 0).t()* inv_two_PI;

  // Set a as in a unit-length circle
  a *= inv_two_PI;
  double exp_a = std::exp(-a);

  // Standardize
  return std::sqrt(n) * (Yn - exp_a) /
    std::sqrt(2 * exp_a * (1 - exp_a * (1 + a * (1 + 0.5 * a))));

}


//' @rdname cir_stat
//' @export
// [[Rcpp::export]]
arma::vec cir_stat_Max_uncover(arma::mat Theta, double a = 2 * pi,
                               bool sorted = false,
                               bool gaps_in_Theta = false) {

  // Sample size
  arma::uword n = Theta.n_rows;

  // Are the gaps stored in Theta?
  if (!gaps_in_Theta) {

    Theta = cir_gaps(Theta, sorted);

  }

  // Statistic for each column
  arma::vec Xn = arma::max(Theta, 0).t();
  Xn = arma::clamp(Xn - a / n, 0, arma::datum::inf) * inv_two_PI;

  // Set a as in a unit-length circle
  a *= inv_two_PI;

  // Standardize
  return n * Xn - std::log(n) + a;

}


//' @rdname cir_stat
//' @export
// [[Rcpp::export]]
arma::vec cir_stat_Num_uncover(arma::mat Theta, double a = 2 * pi,
                               bool sorted = false, bool gaps_in_Theta = false,
                               bool minus_val = true) {

  // Sample size
  arma::uword n = Theta.n_rows;

  // Are the gaps stored in Theta?
  if (!gaps_in_Theta) {

    Theta = cir_gaps(Theta, sorted);

  }

  // Statistic for each column
  arma::vec Ln =
    arma::sum(arma::conv_to<arma::mat>::from(Theta > (a / n)), 0).t();

  // Set a as in a unit-length circle
  a *= inv_two_PI;
  double exp_a = std::exp(-a);

  // Standardize
  double sign = minus_val ? -1.0 : 1.0;
  return sign * (-Ln + n * exp_a) /
    std::sqrt(n * exp_a * (1 - exp_a * (1 + a * a)));

}


//' @rdname cir_stat
//' @export
// [[Rcpp::export]]
arma::vec cir_stat_Gini(arma::mat Theta, bool sorted = false,
                        bool gaps_in_Theta = false) {

  // Sample size
  arma::uword n = Theta.n_rows;

  // Are the gaps stored in Theta?
  if (!gaps_in_Theta) {

    Theta = cir_gaps(Theta, sorted);

  }

  // Theta as cube for Psi_mat
  arma::cube Theta_cube(n, 1, Theta.n_cols);
  Theta_cube.col(0) = Theta;

  // Statistic for each column
  arma::uvec ind_tri = upper_tri_ind(n);
  arma::vec Fn = arma::sum(arma::abs(Psi_mat(Theta_cube, ind_tri,
                                             false, false, true)), 0).t();

  // Standardize
  return std::sqrt(n) * (Fn / (PI * (n - 1.0)) - 1.0);

}


//' @rdname cir_stat
//' @export
// [[Rcpp::export]]
arma::vec cir_stat_Gini_squared(arma::mat Theta, bool sorted = false,
                                bool gaps_in_Theta = false) {

  // Sample size
  arma::uword n = Theta.n_rows;

  // Are the gaps stored in Theta?
  if (!gaps_in_Theta) {

    Theta = cir_gaps(Theta, sorted);

  }

  // Theta as cube for Psi_mat
  arma::cube Theta_cube(n, 1, Theta.n_cols);
  Theta_cube.col(0) = Theta;

  // Statistic for each column
  arma::uvec ind_tri = upper_tri_ind(n);
  arma::vec Qn = arma::sum(arma::square(Psi_mat(Theta_cube, ind_tri,
                                                false, false, true)), 0).t();

  // Standardize
  return std::sqrt(n) * (n / ((n - 1.0) * two_sq_PI) * Qn - 2.0);

}


/*
 * Partition-based tests
 */


//' @rdname cir_stat
//' @export
// [[Rcpp::export]]
arma::vec cir_stat_Ajne(arma::mat Theta, bool Psi_in_Theta = false) {

  // Sample size
  arma::uword n = Psi_in_Theta ?
  n_from_dist_vector(Theta.n_rows) : Theta.n_rows;

  // Number of samples
  arma::uword M = Theta.n_cols;

  // Compute statistic using precomputed Psi matrix?
  if (Psi_in_Theta) {

    // Compute statistic
    return cir_stat_An_Psi(Theta, n);

  } else {

    // Statistic for each column
    arma::vec An = arma::zeros(M);
    arma::uvec ind_tri = upper_tri_ind(n);
    arma::cube Theta_cube(n, 1, 1);
    for (arma::uword k = 0; k < M; k++) {

      // Compute Psi matrix
      Theta_cube.slice(0) = Theta.col(k);
      arma::mat Psi = Psi_mat(Theta_cube, ind_tri, true, false, false);

      // Compute statistic
      An(k) = arma::as_scalar(cir_stat_An_Psi(Psi, n));

    }

    return An;

  }

}


//' @keywords internal
// [[Rcpp::export]]
arma::vec cir_stat_An_Psi(arma::mat Psi, arma::uword n) {

  // Statistic
  arma::vec An = arma::sum(Psi, 0).t();

  // Factors statistic
  An *= -inv_PI / n;
  An += 0.25 * n;
  return An;

}


//' @rdname cir_stat
//' @export
// [[Rcpp::export]]
arma::vec cir_stat_Rothman(arma::mat Theta, double t = 1.0 / 3.0,
                           bool Psi_in_Theta = false) {

  // Sample size
  arma::uword n = Psi_in_Theta ?
  n_from_dist_vector(Theta.n_rows) : Theta.n_rows;

  // Number of samples
  arma::uword M = Theta.n_cols;

  // t's
  t = std::min(t, 1.0 - t);
  double t_m2 = t * (1.0 - t);
  double t_min2 = t * t;

  // Compute statistic using precomputed Psi matrix?
  if (Psi_in_Theta) {

    // Compute statistic
    return cir_stat_Rothman_Psi(Theta, t_m2, t_min2, n);

  } else {

    // Statistic for each column
    arma::vec Ant = arma::zeros(M);
    arma::uvec ind_tri = upper_tri_ind(n);
    arma::cube Theta_cube(n, 1, 1);
    for (arma::uword k = 0; k < M; k++) {

      // Compute Psi matrix
      Theta_cube.slice(0) = Theta.col(k);
      arma::mat Psi = Psi_mat(Theta_cube, ind_tri, true, false, false);

      // Compute statistic
      Ant(k) = arma::as_scalar(cir_stat_Rothman_Psi(Psi, t_m2, t_min2, n));

    }

    return Ant;

  }

}


//' @keywords internal
// [[Rcpp::export]]
arma::vec cir_stat_Rothman_Psi(arma::mat Psi, double t_m2, double t_min2,
                               arma::uword n) {

  // Statistic
  arma::vec Ant = -arma::sum(arma::clamp(Psi * inv_two_PI - t_m2,
                                         -arma::datum::inf, t_min2), 0).t();

  // Factors statistic
  Ant *= 2.0 / n;
  Ant += t_m2;
  return Ant;

}


//' @rdname cir_stat
//' @export
// [[Rcpp::export]]
arma::vec cir_stat_Hodges_Ajne(arma::mat Theta, bool asymp_std = false,
                               bool sorted = false, bool use_Cressie = true) {

  // Sample size
  arma::uword n = Theta.n_rows;

  // Number of samples
  arma::uword M = Theta.n_cols;

  // Sort data on each column
  if (!sorted) {

    Theta = arma::sort(Theta);

  }

  // Computation of Hn
  arma::vec Hn = arma::zeros(M);
  if (use_Cressie) {

    // Compute the statistic as a particular case of the statistic in
    // "On some properties of the scan statistic on the circle and the
    // line" (Cressie, 1977)
    Hn = cir_stat_Cressie(Theta, 0.5, true);

  } else {

    // Compute sup_\alpha N(\alpha) as in "A simple test for uniformity of a
    // circular distribution" (Ajne, 1968)
    arma::vec counts = arma::ones(2 * n);
    counts.tail(n).fill(-1);
    for (arma::uword k = 0; k < M; k++) {

      // Theta + pi
      arma::vec Theta_shift = Theta.col(k) + PI;
      Theta_shift -= two_PI * arma::floor(Theta_shift * inv_two_PI);

      // Sort augmented sample (only the first n elements are necessary)
      arma::uvec sort_ind = arma::sort_index(
        arma::join_horiz(Theta.col(k), Theta_shift)
      );
      arma::vec S = arma::cumsum(counts.elem(sort_ind.head(n)));

      // Compute K_i, i = 1, ..., n
      arma::vec K = arma::zeros(n);
      K(0) = 0.5 * (n - S(n - 1));
      K.tail(n - 1) = S.head(n - 1) + K(0);

      // sup_\alpha N(\alpha) = max(K_1, \ldots, K_{2n})
      Hn(k) = arma::max(arma::max(K, n - K));

    }

  }

  // Standardize in terms of the asymptotic distribution?
  if (asymp_std) {

    Hn = (2 * Hn - n) / std::sqrt(n);

  }
  return Hn;

}


//' @rdname cir_stat
//' @export
// [[Rcpp::export]]
arma::vec cir_stat_Cressie(arma::mat Theta, double t = 1.0 / 3.0,
                           bool sorted = false) {

  // Sample size
  arma::uword n = Theta.n_rows;

  // Sort data on each column
  if (!sorted) {

    Theta = arma::sort(Theta);

  }

  // Number of samples
  arma::uword M = Theta.n_cols;

  // Evaluation of Fn(Theta_{(i)}^-)
  arma::vec i = arma::regspace(0, 1, n - 1.0);

  // For computing Nnt we use that N((Theta_{(i)} + \pi * t) mod 2 * \pi, t) =
  // #{Theta_j \in arc(Theta_{(i)}, Theta_{(i)} + 2\pi t]} and then compute
  // the counts using ecdf_bin() for a better efficiency. The number of points
  // in the arc arc(Theta_{(i)}, Theta_{(i)} + 2\pi t] is
  // F_i - (i - 1) + n * wind_i
  // where
  // F_i = n * Fn((Theta_{(i)} + 2 * \pi * t) mod 2 * \pi)
  // wind_i = floor(((Theta_{(i)} + 2 * \pi * t) mod 2 * \pi) / (2 * \pi))
  // This computation actually is faster for t =  0.5 than the original in
  // cir_stat_Hodges_Ajne that was based on "A simple test for uniformity of a
  // circular distribution" (Ajne, 1968)
  arma::mat Nnt = arma::vec(M);
  double two_PIt = two_PI * (t - arma::datum::eps);
  for (arma::uword k = 0; k < M; k++) {

    // Winding number
    arma::vec wind = arma::floor((Theta.col(k) + two_PIt) * inv_two_PI);

    // Wrapping
    arma::vec wrap = (Theta.col(k) + two_PIt) - two_PI * wind;

    // Sort the evaluation points for ecdf_bin, keeping track of the sorting
    // that was made for later sorting -(i - 1) + n * wind_i
    arma::uvec sort_wrap_ind = arma::sort_index(wrap);
    wrap = wrap(sort_wrap_ind);

    // Call ecdf_bin
    arma::vec F = ecdf_bin(Theta.col(k), wrap, true, true, false);

    // Build - (i - 1) + n * wind_i
    wind *= n;
    wind -= i;

    // Take the max
    Nnt(k) = arma::max(F + wind(sort_wrap_ind));

  }

  return Nnt;

}


//' @rdname cir_stat
//' @export
// [[Rcpp::export]]
arma::vec cir_stat_Feltz_Goldin(arma::mat Theta, bool sorted = false) {

  // Sample size
  arma::uword n = Theta.n_rows;

  // Number of samples
  arma::uword M = Theta.n_cols;

  // Sort data on each column
  if (!sorted) {

    Theta = arma::sort(Theta);

  }

  // Scale
  Theta *= inv_two_PI;

  // Mean
  arma::rowvec Theta_bar = arma::mean(Theta);

  // Substract (i - 0.5) / n to columns
  arma::mat i = arma::repmat(arma::linspace(0.5 / n, 1.0 - 0.5 / n, n), 1, M);

  // Cn1 / 6
  arma::rowvec Cn4 = (arma::sum(arma::square(Theta - i)) +
    1.0 / (12.0 * n)) / 6.0;

  // sigma(j - i)
  arma::mat sigma = arma::ones(n, n);
  sigma.elem(arma::find(arma::trimatl(sigma, -1))).fill(-1);

  // sigma_sum = \sum_i v_i \sum_j v_j^2 sigma(j - i)
  arma::rowvec sigma_sum = arma::sum((sigma * arma::square(Theta)) % Theta);

  // xi
  i = arma::repmat(arma::regspace(2, 1, n + 1), 1, M);
  arma::rowvec xi = n / 15.0 - n * Theta_bar / 3.0 + sigma_sum / n +
    0.5 * arma::sum(arma::square(Theta)) -
    2 * arma::sum(i % arma::pow(Theta, 3)) / (3.0 * n)  +
    arma::sum(arma::pow(Theta, 4)) / 6.0;

  // Rest of Cn4
  Cn4 += -0.5 * n * arma::square(0.5 - Theta_bar) +
    2 * xi - 4 * (0.5 - Theta_bar) %
    (n / 12.0 - arma::sum(arma::pow(Theta, 3), 0) / 3.0) +
    n * arma::square(1.0 / 3.0 - arma::sum(arma::square(Theta)) / n);

  // Statistic for each column
  return Cn4.t();

}


/*
 * Sobolev tests
 */


//' @rdname cir_stat
//' @export
// [[Rcpp::export]]
arma::vec cir_stat_Rayleigh(arma::mat Theta, arma::uword m = 1) {

  // Sample size
  arma::uword n = Theta.n_rows;

  // Statistic for each column
  arma::rowvec Rn = 2 * (arma::square(arma::sum(arma::cos(m * Theta))) +
    arma::square(arma::sum(arma::sin(m * Theta)))) / n;
  return Rn.t();

}


//' @rdname cir_stat
//' @export
// [[Rcpp::export]]
arma::vec cir_stat_Bingham(arma::mat Theta) {

  return sph_stat_Bingham(Theta_to_X(Theta));

}


//' @rdname cir_stat
//' @export
// [[Rcpp::export]]
arma::vec cir_stat_Hermans_Rasson(arma::mat Theta, bool Psi_in_Theta = false) {

  // Sample size
  arma::uword n = Psi_in_Theta ?
    n_from_dist_vector(Theta.n_rows) : Theta.n_rows;

  // Number of samples
  arma::uword M = Theta.n_cols;

  // Compute statistic using precomputed Psi matrix?
  if (Psi_in_Theta) {

    // Compute statistic
    return cir_stat_Hermans_Rasson_Psi(Theta, n);

  } else {

    // Statistic for each column
    arma::vec Hn = arma::zeros(M);
    arma::uvec ind_tri = upper_tri_ind(n);
    arma::cube Theta_cube(n, 1, 1);
    for (arma::uword k = 0; k < M; k++) {

      // Compute Psi matrix
      Theta_cube.slice(0) = Theta.col(k);
      arma::mat Psi = Psi_mat(Theta_cube, ind_tri, true, false, false);

      // Compute statistic
      Hn(k) = arma::as_scalar(cir_stat_Hermans_Rasson_Psi(Psi, n));

    }

    return Hn;

  }

}


//' @keywords internal
// [[Rcpp::export]]
arma::vec cir_stat_Hermans_Rasson_Psi(arma::mat Psi, arma::uword n) {

  // Statistic as h_3 in Pycke (2010)
  arma::vec Hn = -arma::sum(Psi + beta2 * arma::sin(Psi), 0).t();

  // Factors statistic
  Hn *= 2.0 / n;
  Hn += n * const_Hn;
  return Hn;

}


//' @rdname cir_stat
//' @export
// [[Rcpp::export]]
arma::vec cir_stat_Gine_Gn(arma::mat Theta, bool Psi_in_Theta = false) {

  if (Psi_in_Theta) {

    arma::cube Theta_cube(Theta.n_rows, Theta.n_cols, 1);
    Theta_cube.slice(0) = Theta;
    return sph_stat_Gine_Gn(Theta_cube, true, 2);

  } else {

    return sph_stat_Gine_Gn(Theta_to_X(Theta), false, 2);

  }

}


//' @rdname cir_stat
//' @export
// [[Rcpp::export]]
arma::vec cir_stat_Gine_Fn(arma::mat Theta, bool Psi_in_Theta = false) {

  if (Psi_in_Theta) {

    arma::cube Theta_cube(Theta.n_rows, Theta.n_cols, 1);
    Theta_cube.slice(0) = Theta;
    return sph_stat_Gine_Fn(Theta_cube, true, 2);

  } else {

    return sph_stat_Gine_Fn(Theta_to_X(Theta), false, 2);

  }

}


//' @rdname cir_stat
//' @export
// [[Rcpp::export]]
arma::vec cir_stat_Pycke(arma::mat Theta, bool Psi_in_Theta = false) {

  // Sample size
  arma::uword n = Psi_in_Theta ?
  n_from_dist_vector(Theta.n_rows) : Theta.n_rows;

  // Number of samples
  arma::uword M = Theta.n_cols;

  // Compute statistic using precomputed Psi matrix?
  if (Psi_in_Theta) {

    // Compute statistic
    return cir_stat_Pycke_Psi(Theta, n);

  } else {

    // Statistic for each column
    arma::vec Pn = arma::zeros(M);
    arma::uvec ind_tri = upper_tri_ind(n);
    arma::cube Theta_cube(n, 1, 1);
    for (arma::uword k = 0; k < M; k++) {

      // Compute Psi matrix
      Theta_cube.slice(0) = Theta.col(k);
      arma::mat Psi = Psi_mat(Theta_cube, ind_tri, true, false, false);

      // Compute statistic
      Pn(k) = arma::as_scalar(cir_stat_Pycke_Psi(Psi, n));

    }

    return Pn;

  }

}


//' @keywords internal
// [[Rcpp::export]]
arma::vec cir_stat_Pycke_Psi(arma::mat Psi, arma::uword n) {

  // Statistic
  arma::vec Pn = arma::sum(-2 * arma::log(2 * (1 - arma::cos(Psi))), 0).t();
  return Pn / (n - 1.0);

}


//' @rdname cir_stat
//' @export
// [[Rcpp::export]]
arma::vec cir_stat_Pycke_q(arma::mat Theta, bool Psi_in_Theta = false,
                           double q = 0.5) {

  // Sample size
  arma::uword n = Psi_in_Theta ?
  n_from_dist_vector(Theta.n_rows) : Theta.n_rows;

  // Number of samples
  arma::uword M = Theta.n_cols;

  // Compute statistic using precomputed Psi matrix?
  if (Psi_in_Theta) {

    // Compute statistic
    return cir_stat_Pycke_q_Psi(Theta, n, q);

  } else {

    // Statistic for each column
    arma::vec Qn = arma::zeros(M);
    arma::uvec ind_tri = upper_tri_ind(n);
    arma::cube Theta_cube(n, 1, 1);
    for (arma::uword k = 0; k < M; k++) {

      // Compute Psi matrix
      Theta_cube.slice(0) = Theta.col(k);
      arma::mat Psi = Psi_mat(Theta_cube, ind_tri, true, false, false);

      // Compute statistic
      Qn(k) = arma::as_scalar(cir_stat_Pycke_q_Psi(Psi, n, q));

    }

    return Qn;

  }

}


//' @keywords internal
// [[Rcpp::export]]
arma::vec cir_stat_Pycke_q_Psi(arma::mat Psi, arma::uword n, double q = 0.5) {

  // Statistic
  Psi = arma::cos(Psi);
  arma::vec Qn = arma::sum((Psi - q) / (1 - 2 * q * Psi + q * q), 0).t();

  // Factors statistic
  Qn *= 4.0 / n;
  Qn += 2.0 / (1.0 - q);
  return Qn;

}


//' @rdname cir_stat
//' @export
// [[Rcpp::export]]
arma::vec cir_stat_Bakshaev(arma::mat Theta, bool Psi_in_Theta = false) {

  if (Psi_in_Theta) {

    arma::cube Theta_cube(Theta.n_rows, Theta.n_cols, 1);
    Theta_cube.slice(0) = Theta;
    return sph_stat_Riesz(Theta_cube, true, 2, 1.0);

  } else {

    return sph_stat_Riesz(Theta_to_X(Theta), false, 2, 1.0);

  }

}


//' @rdname cir_stat
//' @export
// [[Rcpp::export]]
arma::vec cir_stat_Riesz(arma::mat Theta, bool Psi_in_Theta = false, 
                         double s = 1.0) {

  if (Psi_in_Theta) {

    arma::cube Theta_cube(Theta.n_rows, Theta.n_cols, 1);
    Theta_cube.slice(0) = Theta;
    return sph_stat_Riesz(Theta_cube, true, 2, s);

  } else {

    return sph_stat_Riesz(Theta_to_X(Theta), false, 2, s);

  }

}


//' @rdname cir_stat
//' @export
// [[Rcpp::export]]
arma::vec cir_stat_PCvM(arma::mat Theta, bool Psi_in_Theta = false) {

  if (Psi_in_Theta) {

    arma::cube Theta_cube(Theta.n_rows, Theta.n_cols, 1);
    Theta_cube.slice(0) = Theta;
    return sph_stat_PCvM(Theta_cube, true, 2, 0, 0);

  } else {

    return sph_stat_PCvM(Theta_to_X(Theta), false, 2, 0, 0);

  }

}


//' @rdname cir_stat
//' @export
// [[Rcpp::export]]
arma::vec cir_stat_PRt(arma::mat Theta, double t = 1.0 / 3.0,
                       bool Psi_in_Theta = false) {

  if (Psi_in_Theta) {

    arma::cube Theta_cube(Theta.n_rows, Theta.n_cols, 1);
    Theta_cube.slice(0) = Theta;
    return sph_stat_PRt(Theta_cube, t, true, 2, 0, 0);

  } else {

    return sph_stat_PRt(Theta_to_X(Theta), t, Psi_in_Theta, 2, 0, 0);

  }

}


//' @rdname cir_stat
//' @export
// [[Rcpp::export]]
arma::vec cir_stat_PAD(arma::mat Theta, bool Psi_in_Theta = false) {

  if (Psi_in_Theta) {

    arma::cube Theta_cube(Theta.n_rows, Theta.n_cols, 1);
    Theta_cube.slice(0) = Theta;
    return sph_stat_PAD(Theta_cube, true, 2, 0, 0);

  } else {

    return sph_stat_PAD(Theta_to_X(Theta), false, 2, 0, 0);

  }

}


/*
 * Other tests
 */


//' @rdname cir_stat
//' @export
// [[Rcpp::export]]
arma::vec cir_stat_Cuesta_Albertos(arma::mat Theta, arma::mat rand_dirs,
                                   arma::uword K_Cuesta_Albertos = 25,
                                   bool original = false) {

  return sph_stat_Cuesta_Albertos(Theta_to_X(Theta), rand_dirs,
                                  K_Cuesta_Albertos, original);

}

