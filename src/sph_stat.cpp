
// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <Rmath.h>
#include <R.h>
using namespace Rcpp;

// Declaration for functions
arma::vec sph_stat_Rayleigh(arma::cube X);
arma::mat Psi_mat(arma::cube data, arma::uvec ind_tri, bool use_ind_tri,
                  bool scalar_prod, bool angles_diff);
arma::uvec upper_tri_ind(arma::uword n);
arma::uword n_from_dist_vector(arma::uword n_dist);
arma::vec t_inv_sqrt_one(arma::vec t);
arma::mat sort_each_col(arma::mat A);
arma::vec d_proj_unif(arma::vec x, arma::uword p, bool log);
arma::vec p_proj_unif(arma::vec x, arma::uword p, bool log);
arma::vec q_proj_unif(arma::vec u, arma::uword p);
arma::vec p_Kolmogorov(arma::vec, arma::uword K_Kolmogorov, bool alternating);
arma::vec Gauss_Legen_nodes(double a, double b, arma::uword N);
arma::vec Gauss_Legen_weights(double a, double b, arma::uword N);
arma::vec cir_stat_An_Psi(arma::mat Psi, arma::uword n);
arma::vec sph_stat_Gine_Gn_Psi(arma::mat Psi, arma::uword n, arma::uword p);
arma::vec sph_stat_Gine_Fn_Psi(arma::mat Psi, arma::uword n, arma::uword p);
arma::vec sph_stat_Pycke_Psi(arma::mat Psi, arma::uword n, arma::uword p);
arma::vec sph_stat_Riesz(arma::cube X, bool Psi_in_X, arma::uword p, double s);
arma::vec sph_stat_Riesz_Psi(arma::mat Psi, arma::uword n, double s);
arma::vec sph_stat_PCvM_Psi(arma::mat Psi, arma::uword n, arma::uword p,
                            arma::vec th_grid, arma::vec int_grid);
arma::vec sph_stat_PRt_Psi(arma::mat Psi, double t_m, double theta_t_m,
                           arma::uword n, arma::uword p, arma::vec th_grid,
                           arma::vec int_grid);
arma::vec sph_stat_PAD_Psi(arma::mat Psi, arma::uword n, arma::uword p,
                           arma::vec th_grid, arma::vec int_grid);
arma::vec sph_stat_Cai_Psi(arma::mat Psi, arma::uword n, arma::uword p);

// Constants
const double inv_M_PI = 1.0 / M_PI;
const double inv_two_M_PI = 0.5 / M_PI;
const double inv_two_M_PI_sq = 0.25 / (M_PI * M_PI);
const double two_M_PI = 2.0 * M_PI;
const double sqrt_M_PI = std::sqrt(M_PI);
const double log_two_M_PI = std::log(2 * M_PI);
const double log_two = std::log(2.0);
const double inv_four_M_PI = 0.25 / M_PI;
const double const_Pycke = (1.0 - std::log(4.0)) * inv_four_M_PI;

/*
 * Sobolev tests
 */


//' @rdname sph_stat
//' @export
// [[Rcpp::export]]
arma::vec sph_stat_Rayleigh(arma::cube X) {

  // Sample size
  arma::uword n = X.n_rows;

  // Dimension
  arma::uword p = X.n_cols;

  // Statistic applied to each slice
  arma::vec Rn = n * p * sum(arma::square(arma::mean(X)), 1);
  return Rn;

}


//' @rdname sph_stat
//' @export
// [[Rcpp::export]]
arma::vec sph_stat_Bingham(arma::cube X) {

  // Sample size
  arma::uword n = X.n_rows;

  // Dimension
  arma::uword p = X.n_cols;

  // Number of samples
  arma::uword M = X.n_slices;

  // Statistic applied to each slice
  arma::vec Bn = arma::zeros(M);
  for (arma::uword k = 0; k < M; k++) {

    // Second moment without 1 / n
    arma::mat S = X.slice(k).t() * X.slice(k);

    // Trace of the squared second moment
    Bn(k) = arma::trace(S * S) / (n * n);

  }

  // Factors Bn
  Bn -= 1.0 / p;
  Bn *= 0.5 * n * p * (p + 2);
  return Bn;

}


//' @rdname sph_stat
//' @export
// [[Rcpp::export]]
arma::vec sph_stat_Ajne(arma::cube X, bool Psi_in_X = false) {

  // Sample size
  arma::uword n = Psi_in_X ? n_from_dist_vector(X.n_rows) : X.n_rows;

  // Number of samples
  arma::uword M = Psi_in_X ? X.n_cols : X.n_slices;

  // Compute statistic using precomputed Psi matrix?
  if (Psi_in_X) {

    // Compute statistic
    return cir_stat_An_Psi(X.slice(0), n);

  } else {

    // Statistic for each slice
    arma::vec An = arma::zeros(M);
    arma::uvec ind_tri = upper_tri_ind(n);
    for (arma::uword k = 0; k < M; k++) {

      // Compute Psi matrix
      arma::mat Psi = Psi_mat(X(arma::span::all, arma::span::all,
                                arma::span(k)), ind_tri, true, false, false);

      // Compute statistic
      An(k) = arma::as_scalar(cir_stat_An_Psi(Psi, n));

    }

    return An;

  }

}


//' @rdname sph_stat
//' @export
// [[Rcpp::export]]
arma::vec sph_stat_Gine_Gn(arma::cube X, bool Psi_in_X = false,
                           arma::uword p = 0) {

  // Sample size
  arma::uword n = Psi_in_X ? n_from_dist_vector(X.n_rows) : X.n_rows;

  // Dimension
  p = Psi_in_X ? p : X.n_cols;
  if (Psi_in_X && (p == 0)) {

    stop("p >= 2 must be specified if Psi_in_X = TRUE.");

  }

  // Number of samples
  arma::uword M = Psi_in_X ? X.n_cols : X.n_slices;

  // Compute statistic using precomputed Psi matrix?
  if (Psi_in_X) {

    // Compute statistic
    return sph_stat_Gine_Gn_Psi(X.slice(0), n, p);

  } else {

    // Statistic for each slice
    arma::vec Gn = arma::zeros(M);
    arma::uvec ind_tri = upper_tri_ind(n);
    for (arma::uword k = 0; k < M; k++) {

      // Compute Psi matrix
      arma::mat Psi = Psi_mat(X(arma::span::all, arma::span::all,
                                arma::span(k)), ind_tri, true, false, false);

      // Compute statistic
      Gn(k) = arma::as_scalar(sph_stat_Gine_Gn_Psi(Psi, n, p));

    }

    return Gn;

  }

}


//' @keywords internal
// [[Rcpp::export]]
arma::vec sph_stat_Gine_Gn_Psi(arma::mat Psi, arma::uword n, arma::uword p) {

  // Statistic
  arma::vec Gn = arma::sum(arma::sin(Psi), 0).t();

  // Factors statistic
  Gn *= -(p - 1.0) / (2.0 * n) *
    std::pow(R::gammafn(0.5 * (p - 1)) / R::gammafn(0.5 * p), 2);
  Gn += 0.5 * n;
  return Gn;

}


//' @rdname sph_stat
//' @export
// [[Rcpp::export]]
arma::vec sph_stat_Gine_Fn(arma::cube X, bool Psi_in_X = false,
                           arma::uword p = 0) {

  // Sample size
  arma::uword n = Psi_in_X ? n_from_dist_vector(X.n_rows) : X.n_rows;

  // Dimension
  p = Psi_in_X ? p : X.n_cols;
  if (Psi_in_X && (p == 0)) {

    stop("p >= 2 must be specified if Psi_in_X = TRUE.");

  }

  // Number of samples
  arma::uword M = Psi_in_X ? X.n_cols : X.n_slices;

  // Compute statistic using precomputed Psi matrix?
  if (Psi_in_X) {

    // Compute statistic
    return sph_stat_Gine_Fn_Psi(X.slice(0), n, p);

  } else {

    // Statistic for each slice
    arma::vec Fn = arma::zeros(M);
    arma::uvec ind_tri = upper_tri_ind(n);
    for (arma::uword k = 0; k < M; k++) {

      // Compute Psi matrix
      arma::mat Psi = Psi_mat(X(arma::span::all, arma::span::all,
                                arma::span(k)), ind_tri, true, false, false);

      // Compute statistic
      Fn(k) = arma::as_scalar(sph_stat_Gine_Fn_Psi(Psi, n, p));

    }

    return Fn;

  }

}


//' @keywords internal
// [[Rcpp::export]]
arma::vec sph_stat_Gine_Fn_Psi(arma::mat Psi, arma::uword n, arma::uword p) {

  // Statistics An and Gn
  arma::vec An = arma::sum(Psi, 0).t();
  arma::vec Gn = arma::sum(arma::sin(Psi), 0).t();

  // Factors An
  An *= -inv_M_PI / n;
  An += 0.25 * n;

  // Factors Gn
  Gn *= -(p - 1.0) * pow(R::gammafn(0.5 * (p - 1)), 2) /
    (2 * n * std::pow(R::gammafn(0.5 * p), 2));
  Gn += 0.5 * n;

  // Fn
  return 4.0 * An + Gn;

}


//' @rdname sph_stat
//' @export
// [[Rcpp::export]]
arma::vec sph_stat_Pycke(arma::cube X, bool Psi_in_X = false,
                         arma::uword p = 0) {

  // Sample size
  arma::uword n = Psi_in_X ? n_from_dist_vector(X.n_rows) : X.n_rows;

  // Dimension
  p = Psi_in_X ? p : X.n_cols;
  if (Psi_in_X && (p == 0)) {

    stop("p >= 2 must be specified if Psi_in_X = TRUE.");

  }

  // Compute Riesz statistic?
  if (p > 3) {

    // Undo the cosine of shortest angles if Psi is given
    if (Psi_in_X) {

      X = arma::acos(X);

    }

    // Compute Riesz statistic
    Rcpp::warning("Pycke statistic is only defined for p = 2,3. Using Riesz statistic with s = 0 instead, which behaves consistently across dimensions.");
    arma::vec Gamman = sph_stat_Riesz(X, Psi_in_X, p, 0);
    return Gamman;

  }

  // Number of samples
  arma::uword M = Psi_in_X ? X.n_cols : X.n_slices;

  // Compute statistic using precomputed Psi matrix?
  arma::vec Gamman = arma::zeros(M);
  if (Psi_in_X) {

    // Compute statistic
    Gamman = sph_stat_Pycke_Psi(X.slice(0), n, p);

  } else {

    // Statistic for each slice
    arma::uvec ind_tri = upper_tri_ind(n);
    for (arma::uword k = 0; k < M; k++) {

      // Compute Psi matrix
      arma::mat Psi = Psi_mat(X(arma::span::all, arma::span::all,
                                arma::span(k)), ind_tri, true, true, false);

      // Compute statistic
      Gamman(k) = arma::as_scalar(sph_stat_Pycke_Psi(Psi, n, p));

    }

  }

  return Gamman;

}


//' @keywords internal
// [[Rcpp::export]]
arma::vec sph_stat_Pycke_Psi(arma::mat Psi, arma::uword n, arma::uword p) {

  // log(Gamman) without n-constants: 2 \sum_{i < j} \log(sqrt(1 - \Psi_{ij}))
  // Psi contains scalar products!
  arma::vec Gamman = arma::sum(arma::log1p(-Psi), 0).t();

  // Factors statistic
  if (p == 2) {

    Gamman *= -2.0 / (n - 1.0); // Factor 2 of log(Gamman^2)
    Gamman += -log_two * n; // -n * log(2) / 2 * 2 where the
    // last 2 is because of factor 2 of log(Gamman^2)

  } else if (p == 3) {

    Gamman *= -inv_two_M_PI / (n - 1.0); // Included factor 2 of log(Gamman^2)
    // and 1 / (4 * pi)
    Gamman += -(log_two * inv_four_M_PI + const_Pycke) * n;

  } else {

    stop("Pycke statistic is only defined for p = 2,3.");

  }
  return Gamman;

}


//' @rdname sph_stat
//' @export
// [[Rcpp::export]]
arma::vec sph_stat_Bakshaev(arma::cube X, bool Psi_in_X = false,
                            arma::uword p = 0) {

  return sph_stat_Riesz(X, Psi_in_X, p, 1.0);

}


//' @rdname sph_stat
//' @export
// [[Rcpp::export]]
arma::vec sph_stat_Riesz(arma::cube X, bool Psi_in_X = false,
                         arma::uword p = 0, double s = 1.0) {

  // Sample size
  arma::uword n = Psi_in_X ? n_from_dist_vector(X.n_rows) : X.n_rows;

  // Dimension
  p = Psi_in_X ? p : X.n_cols;
  if (Psi_in_X && (p == 0)) {

    stop("p >= 2 must be specified if Psi_in_X = TRUE.");

  }

  // Number of samples
  arma::uword M = Psi_in_X ? X.n_cols : X.n_slices;

  // Compute statistic using precomputed Psi matrix?
  arma::vec Rn = arma::zeros(M);
  if (Psi_in_X) {

    // Compute statistic
    Rn = sph_stat_Riesz_Psi(X.slice(0), n, s);

  } else {

    // Statistic for each slice
    arma::uvec ind_tri = upper_tri_ind(n);
    for (arma::uword k = 0; k < M; k++) {

      // Compute Psi matrix
      arma::mat Psi = Psi_mat(X(arma::span::all, arma::span::all,
                                arma::span(k)), ind_tri, true, false, false);

      // Compute statistic
      Rn(k) = arma::as_scalar(sph_stat_Riesz_Psi(Psi, n, s));

    }

  }

  // Compute bias integral
  double tau = 0;
  if (s == 0) {

    if (p == 2) {

      tau = 0;

    } else if (p == 3) {

      tau = log_two - 0.5;

    } else if (p == 4) {

      tau = 0.25;

    } else {

      arma::vec m = arma::regspace(1, 1, p - 2);
      arma::vec signs = 2 * m - 4 * arma::ceil(0.5 * m) + 1;
      tau = log_two +
        std::pow(-1, p - 1) * (log_two + arma::accu(signs / m));
      tau = 0.5 * tau;

    }
    Rn += n * tau;

  } else {

    if (p == 2) {

      tau = std::pow(2, s);

    } else {

      tau = std::pow(2, p + s - 3) * (p - 2) * R::gammafn(0.5 * p - 1);

    }
    Rn += n * tau * R::gammafn(0.5 * (p - 1 + s)) /
      (sqrt_M_PI * R::gammafn(p - 1 + 0.5 * s));

  }
  return Rn;

}


//' @keywords internal
// [[Rcpp::export]]
arma::vec sph_stat_Riesz_Psi(arma::mat Psi, arma::uword n, double s) {

  // Statistic
  arma::vec Rn = arma::zeros(Psi.n_cols);
  if (s == 0) {

    // log(Gamman) without n-constants: 2 \sum_{i < j} \log(sqrt(1 - \Psi_{ij}))
    Rn = arma::sum(arma::log1p(-arma::cos(Psi)), 0).t();

    // Divide by n and add log(sqrt(2))
    Rn *= -1.0 / n;
    Rn += -0.5 * log_two * (n - 1.0); // -n * log(2) / 2

  } else {

    Rn = -std::pow(2, s + 1) *
      arma::sum(arma::pow(arma::sin(0.5 * Psi), s), 0).t() / n;

  }
  return Rn;

}


//' @rdname sph_stat
//' @export
// [[Rcpp::export]]
arma::vec sph_stat_PCvM(arma::cube X, bool Psi_in_X = false, arma::uword p = 0,
                        arma::uword N = 160, arma::uword L = 1e3) {

  // Sample size
  arma::uword n = Psi_in_X ? n_from_dist_vector(X.n_rows) : X.n_rows;

  // Dimension
  p = Psi_in_X ? p : X.n_cols;
  if (Psi_in_X && (p == 0)) {

    stop("p >= 2 must be specified if Psi_in_X = TRUE.");

  }

  // Number of samples
  arma::uword M = Psi_in_X ? X.n_cols : X.n_slices;

  // Compute integral on a grid if p > 4
  arma::vec th_grid = arma::linspace(0, M_PI, L);
  arma::vec int_grid = arma::zeros(L);
  if (p > 4) {

    for (arma::uword k = 0; k < L; k++) {

      // Integration nodes and weights
      double theta = th_grid(k);
      double cos_theta = std::cos(0.5 * theta);
      arma::vec t_k = Gauss_Legen_nodes(0, cos_theta, N);
      arma::vec w_k = Gauss_Legen_weights(0, cos_theta, N);

      // Integral
      int_grid(k) = arma::accu(w_k % d_proj_unif(t_k, p, false) %
                               p_proj_unif(t_k, p, false) %
                               p_proj_unif(std::tan(0.5 * theta) *
                                 t_inv_sqrt_one(t_k), p - 1, false));

    }

    // Debug int_grid -- see check-NaNs-PAD.R
    // int_grid.print();

  }

  // Compute statistic using precomputed Psi matrix?
  if (Psi_in_X) {

    // Compute statistic
    return sph_stat_PCvM_Psi(X.slice(0), n, p, th_grid, int_grid);

  } else {

    // Statistic for each slice
    arma::vec PCvMn = arma::zeros(M);
    arma::uvec ind_tri = upper_tri_ind(n);
    for (arma::uword k = 0; k < M; k++) {

      // Compute Psi matrix
      arma::mat Psi = Psi_mat(X(arma::span::all, arma::span::all,
                                arma::span(k)), ind_tri, true, false, false);

      // Compute statistic
      PCvMn(k) = arma::as_scalar(sph_stat_PCvM_Psi(Psi, n, p, th_grid,
                                                   int_grid));

    }

    return PCvMn;

  }

}


//' @keywords internal
// [[Rcpp::export]]
arma::vec sph_stat_PCvM_Psi(arma::mat Psi, arma::uword n, arma::uword p,
                            arma::vec th_grid, arma::vec int_grid) {

  // Create returned statistic
  arma::vec PCvMn = arma::zeros(Psi.n_cols);

  // Compute statistic from Psi matrix. Observe that a constant addend c inside
  // the sum becomes n * (n - 1) / 2 * c
  if (p == 2) {

    // Addends
    Psi %= Psi - two_M_PI;

    // Sum
    PCvMn = 0.25 * n * (n - 1) + arma::sum(Psi, 0).t() * inv_two_M_PI_sq;

  } else if (p == 3) {

    // Addends
    Psi = arma::sin(0.5 * Psi);

    // Sum
    PCvMn = 0.25 * n * (n - 1) - 0.25 * arma::sum(Psi, 0).t();

  } else if (p == 4) {

    // Addends
    Psi = Psi % (Psi - two_M_PI) + (M_PI - Psi) % arma::tan(0.5 * Psi) -
      2 * arma::square(arma::sin(0.5 * Psi));

    // Sum
    PCvMn = 0.25 * n * (n - 1) + arma::sum(Psi, 0).t() * inv_two_M_PI_sq;

  } else if (p > 4) {

    // Flatten Psi
    arma::uword n_rows_Psi = Psi.n_rows;
    arma::uword n_cols_Psi = Psi.n_cols;
    Psi.reshape(n_rows_Psi * n_cols_Psi, 1);

    // Perform interpolation for integral
    arma::vec int_interp = arma::zeros(Psi.n_elem);
    arma::interp1(th_grid, int_grid, Psi, int_interp);

    // Addends
    Psi = Psi * inv_two_M_PI +
      2 * arma::square(p_proj_unif(arma::cos(0.5 * Psi), p, false)) -
      4 * int_interp;
    Psi.reshape(n_rows_Psi, n_cols_Psi);

    // Sum
    PCvMn = -0.375 * n * (n - 1) + arma::sum(Psi, 0).t();

  } else {

    stop("p must be >= 2.");

  }

  // Factors statistic
  PCvMn *= 2.0 / n;
  PCvMn += 0.5 - n / 3.0;
  return PCvMn;

}


//' @rdname sph_stat
//' @export
// [[Rcpp::export]]
arma::vec sph_stat_PRt(arma::cube X, double t = 1.0 / 3.0,
                       bool Psi_in_X = false, arma::uword p = 0,
                       arma::uword N = 160, arma::uword L = 1e3) {

  // Sample size
  arma::uword n = Psi_in_X ? n_from_dist_vector(X.n_rows) : X.n_rows;

  // Dimension
  p = Psi_in_X ? p : X.n_cols;
  if (Psi_in_X && (p == 0)) {

    stop("p >= 2 must be specified if Psi_in_X = TRUE.");

  }

  // Number of samples
  arma::uword M = Psi_in_X ? X.n_cols : X.n_slices;

  // t's quantities
  double t_m = std::min(t, 1.0 - t);
  double theta_t_m = 2.0 * std::acos(-arma::as_scalar(q_proj_unif({t_m}, p)));

  // Compute integral on a grid if p > 4
  arma::vec th_grid = arma::linspace(0, M_PI, L);
  arma::vec int_grid = arma::zeros(L);
  if (p > 4) {

    // Integration nodes and weights
    double cos_theta_t_m = std::cos(0.5 * theta_t_m);
    arma::vec t_k = Gauss_Legen_nodes(0, cos_theta_t_m, N);
    arma::vec w_k = Gauss_Legen_weights(0, cos_theta_t_m, N);

    // Integral
    for (arma::uword k = 0; k < L; k++) {

      int_grid(k) = arma::accu(w_k % d_proj_unif(t_k, p, false) %
                               p_proj_unif(std::tan(0.5 * th_grid(k)) *
                                 t_inv_sqrt_one(t_k), p - 1, false));

    }

    // Debug int_grid -- see check-NaNs-PAD.R
    // int_grid.print();

  }

  // Compute statistic using precomputed Psi matrix?
  if (Psi_in_X) {

    // Compute statistic
    return sph_stat_PRt_Psi(X.slice(0), t_m, theta_t_m, n, p, th_grid,
                            int_grid);

  } else {

    // Statistic for each slice
    arma::vec PRtn = arma::zeros(M);
    arma::uvec ind_tri = upper_tri_ind(n);
    for (arma::uword k = 0; k < M; k++) {

      // Compute Psi matrix
      arma::mat Psi = Psi_mat(X(arma::span::all, arma::span::all,
                                arma::span(k)), ind_tri, true, false, false);

      // Compute statistic
      PRtn(k) = arma::as_scalar(sph_stat_PRt_Psi(Psi, t_m, theta_t_m, n, p,
                                                 th_grid, int_grid));

    }

    return PRtn;

  }

}


//' @keywords internal
// [[Rcpp::export]]
arma::vec sph_stat_PRt_Psi(arma::mat Psi, double t_m, double theta_t_m,
                           arma::uword n, arma::uword p, arma::vec th_grid,
                           arma::vec int_grid) {

  // Create returned statistic
  arma::vec PRtn = arma::zeros(Psi.n_cols);

  // Compute statistic from Psi vector. Observe that a constant addend c
  // inside the sum becomes 0.5 * n * (n - 1) * c
  double t_m2 = t_m * (1 - t_m);
  if (p == 2) {

    Psi = Psi * inv_two_M_PI - t_m2;
    PRtn = 0.5 * (0.5 - t_m2) * n * (n - 1) -
      arma::sum(arma::clamp(Psi, -arma::datum::inf, t_m * t_m), 0).t();

  } else {

    // Flatten Psi
    arma::uword n_rows_Psi = Psi.n_rows;
    arma::uword n_cols_Psi = Psi.n_cols;
    Psi.reshape(n_rows_Psi * n_cols_Psi, 1);

    // Index for the variable part
    arma::uvec ind_var = arma::find(Psi < theta_t_m);

    // Constant part
    double half_minus_t_m = 0.5 - t_m;
    Psi.transform([theta_t_m, half_minus_t_m](double theta) {
      return (theta >= theta_t_m) ? half_minus_t_m : theta;
    });

    // Variable part
    if (p == 3) {

      Psi.elem(ind_var) *= 0.5;
      Psi.elem(ind_var) = -t_m + 0.5 +
        inv_M_PI * ((2 * t_m - 1) *
          arma::acos(((0.5 - t_m) / std::sqrt(t_m * (1 - t_m))) *
            arma::tan(Psi.elem(ind_var))) +
          arma::atan(arma::sqrt(arma::square(arma::cos(Psi.elem(ind_var))) -
            std::pow(1 - 2 * t_m, 2)) / arma::sin(Psi.elem(ind_var))));

    } else if (p == 4) {

      Psi.elem(ind_var) = 0.5 + t_m +
        inv_M_PI * (-0.5 * (Psi.elem(ind_var) + theta_t_m) +
          0.5 * std::sin(theta_t_m) + arma::tan(0.5 * Psi.elem(ind_var)) *
          std::pow(std::cos(0.5 * theta_t_m), 2));

    } else if (p > 4) {

      // Perform interpolation for integral for theta < theta_t_m
      arma::vec int_interp = arma::zeros(ind_var.n_elem);
      arma::interp1(th_grid, int_grid, Psi.elem(ind_var), int_interp);

      // Integral
      Psi.elem(ind_var) = t_m - Psi.elem(ind_var) * inv_two_M_PI +
        2 * int_interp;

    } else {

      stop("p must be >= 2.");

    }

    // Sum
    Psi.reshape(n_rows_Psi, n_cols_Psi);
    PRtn = arma::sum(Psi, 0).t();

  }

  // Factors statistic
  PRtn *= 2.0 / n;
  PRtn += 0.5 * (1.0 - n) + n * t_m2;
  return PRtn;

}


//' @rdname sph_stat
//' @export
// [[Rcpp::export]]
arma::vec sph_stat_PAD(arma::cube X, bool Psi_in_X = false, arma::uword p = 0,
                       arma::uword N = 160, arma::uword L = 1e3) {

  // Sample size
  arma::uword n = Psi_in_X ? n_from_dist_vector(X.n_rows) : X.n_rows;

  // Dimension
  p = Psi_in_X ? p : X.n_cols;
  if (Psi_in_X && (p == 0)) {

    stop("p >= 2 must be specified if Psi_in_X = TRUE.");

  }

  // Number of samples
  arma::uword M = Psi_in_X ? X.n_cols : X.n_slices;

  // Compute integral on a grid if p > 2
  arma::vec th_grid = arma::linspace(0, M_PI, L);
  arma::vec int_grid = arma::zeros(L);
  if (p > 2) {

    for (arma::uword k = 0; k < L; k++) {

      // Integration nodes and weights
      double theta = th_grid(k);
      double cos_theta = std::cos(0.5 * theta);
      arma::vec t_k = Gauss_Legen_nodes(0, cos_theta, N);
      arma::vec w_k = Gauss_Legen_weights(0, cos_theta, N);

      // Integral
      if (p == 3) {

        int_grid(k) = arma::accu(w_k % arma::log((1 + t_k) / (1 - t_k)) %
                                 arma::acos(std::tan(0.5 * theta) *
                                   t_inv_sqrt_one(t_k)));

      } else if (p == 4) {

        int_grid(k) = arma::accu(w_k % t_k % arma::log(M_PI / (arma::acos(t_k) -
                                 t_k % arma::sqrt(1 - arma::square(t_k))) - 1));

      } else if (p > 4){

        // Problematic difference of logs with overflow
        arma::vec logs = -p_proj_unif(-t_k, p, TRUE);
        logs.elem(arma::find_nonfinite(logs)).zeros();
        logs += p_proj_unif(t_k, p, TRUE);

        // Integral
        int_grid(k) = arma::accu(w_k % d_proj_unif(t_k, p, false) % logs %
          p_proj_unif(-std::tan(0.5 * theta) * t_inv_sqrt_one(t_k),
                      p - 1, false));

      }

    }

    // Debug int_grid -- see check-NaNs-PAD.R
    // int_grid.print();

  }

  // Compute statistic using precomputed Psi matrix?
  if (Psi_in_X) {

    // Compute statistic
    return sph_stat_PAD_Psi(X.slice(0), n, p, th_grid, int_grid);

  } else {

    // Statistic for each slice
    arma::vec PADn = arma::zeros(M);
    arma::uvec ind_tri = upper_tri_ind(n);
    for (arma::uword k = 0; k < M; k++) {

      // Compute Psi matrix
      arma::mat Psi = Psi_mat(X(arma::span::all, arma::span::all,
                                arma::span(k)), ind_tri, true, false, false);

      // Compute statistic
      PADn(k) = arma::as_scalar(sph_stat_PAD_Psi(Psi, n, p, th_grid, int_grid));

    }

    return PADn;

  }

}


//' @keywords internal
// [[Rcpp::export]]
arma::vec sph_stat_PAD_Psi(arma::mat Psi, arma::uword n, arma::uword p,
                           arma::vec th_grid, arma::vec int_grid) {

  // Create returned statistic
  arma::vec PADn = arma::zeros(Psi.n_cols);

  // Compute statistic from Psi vector. Observe that a constant addend c
  // inside the sum becomes 0.5 * n * (n - 1) * c
  if (p == 2) {

    // Addends
    Psi = Psi % arma::log(Psi) + (two_M_PI - Psi) % arma::log(two_M_PI - Psi);

    // Replace NaNs (created by theta = 0) with 0
    Psi.replace(arma::datum::nan, 0);

    // Sum
    PADn = -log_two_M_PI * n * (n - 1) + arma::sum(Psi, 0).t() * inv_M_PI;

  } else {

    // Flatten Psi
    arma::uword n_rows_Psi = Psi.n_rows;
    arma::uword n_cols_Psi = Psi.n_cols;
    Psi.reshape(n_rows_Psi * n_cols_Psi, 1);

    // Perform interpolation for integral
    arma::mat int_interp = arma::zeros(Psi.n_elem, 1);
    arma::interp1(th_grid, int_grid, Psi, int_interp);

    if (p == 3) {

      // Sum
      int_interp.reshape(n_rows_Psi, n_cols_Psi);
      PADn = -log_two * n * (n - 1) + 2 * inv_M_PI * arma::sum(int_interp, 0).t();

    } else if (p == 4) {

      // Addends
      arma::vec s_theta = Psi - arma::sin(Psi);
      s_theta = s_theta % arma::log(s_theta) +
        (two_M_PI - s_theta) % arma::log(two_M_PI - s_theta);
      Psi = s_theta - 4 * arma::tan(0.5 * Psi) % int_interp;

      // Replace NaNs (created by theta = 0) with 0
      Psi.replace(arma::datum::nan, 0);

      // Sum
      Psi.reshape(n_rows_Psi, n_cols_Psi);
      PADn = -log_two_M_PI * n * (n - 1) + inv_M_PI * arma::sum(Psi, 0).t();

    } else if (p > 4){

      // Sum
      int_interp.reshape(n_rows_Psi, n_cols_Psi);
      PADn = -log_two * n * (n - 1) + 4 * arma::sum(int_interp, 0).t();

    } else {

      stop("p must be >= 2.");

    }

  }

  // Factors statistic
  PADn *= 2.0 / n;
  PADn += n;
  return PADn;

}


/*
 * Other tests
 */


//' @rdname sph_stat
//' @export
// [[Rcpp::export]]
arma::vec sph_stat_Cuesta_Albertos(arma::cube X, arma::mat rand_dirs,
                                   arma::uword K_Cuesta_Albertos = 25,
                                   bool original = false) {

  // Sample size
  arma::uword n = X.n_rows;

  // Number of samples
  arma::uword M = X.n_slices;

  // Check dimension
  arma::uword p = X.n_cols;
  if (p != rand_dirs.n_cols) {

    stop("The dimension of the directions of X and rand_dirs are incompatible.");

  }

  // Number of projections
  arma::uword n_proj = rand_dirs.n_rows;

  // Project data
  arma::cube X_proj = arma::zeros(n, n_proj, M);
  arma::mat rand_dirs_t = rand_dirs.t();
  for (arma::uword k = 0; k < M; k++) {

    X_proj.slice(k) = sort_each_col(X.slice(k) * rand_dirs_t);

  }

  // Apply F0
  for (arma::uword k = 0; k < M; k++) {

    X_proj.slice(k) = arma::reshape(
      p_proj_unif(arma::vectorise(X_proj.slice(k)), p, false), n, n_proj);

  }

  // Subtract i / n to slices and columns, X becomes Dn^+ when taking the
  // -min on each of the slices and columns
  arma::vec i = arma::linspace(1.0 / n, 1.0, n);
  X_proj.each_slice() -= arma::repmat(i, 1, n_proj);

  // K-S statistics
  arma::mat KSn = arma::zeros(n_proj, M);
  for (arma::uword k = 0; k < M; k++) {

    // max(Dn^+, Dn^-) for each column
    KSn.col(k) = arma::max(arma::join_vert(-arma::min(X_proj.slice(k)),
                           arma::max(X_proj.slice(k)) + 1.0 / n)).t();

  }
  KSn *= std::sqrt(n);

  // The statistic is min(pval_1, ..., pval_{n_proj}), but we compute
  // 1 - K(max(KS_1, ..., KS_{n_proj})) which is equivalent and faster

  // The maximum of the statistics
  arma::vec CAn = arma::max(KSn, 0).t();

  // Original definition of the statistic or just the maximum of statistics?
  if (original) {

    CAn = 1.0 - p_Kolmogorov(CAn, K_Cuesta_Albertos, true);

  }
  return CAn;

}


/*
 * High-dimensional tests
 */


//' @rdname sph_stat
//' @export
// [[Rcpp::export]]
arma::vec sph_stat_Rayleigh_HD(arma::cube X) {

  // Unstandardized statistic
  arma::vec Rn = sph_stat_Rayleigh(X);

  // Dimension
  arma::uword p = X.n_cols;

  // HD standardization
  Rn -= p;
  Rn /= std::sqrt(2.0 * p);
  return Rn;

}


//' @rdname sph_stat
//' @export
// [[Rcpp::export]]
arma::vec sph_stat_Cai(arma::cube X, arma::uword regime = 3,
                       bool Psi_in_X = false, arma::uword p = 0) {

  // Sample size
  arma::uword n = Psi_in_X ? n_from_dist_vector(X.n_rows) : X.n_rows;

  // Dimension
  p = Psi_in_X ? p : X.n_cols;
  if (Psi_in_X && (p == 0)) {

    stop("p >= 2 must be specified if Psi_in_X = TRUE.");

  }

  // Number of samples
  arma::uword M = Psi_in_X ? X.n_cols : X.n_slices;

  // Compute statistic using precomputed Psi matrix?
  arma::vec Cn = arma::zeros(M);
  if (Psi_in_X) {

    // Compute statistic
    Cn = sph_stat_Cai_Psi(X.slice(0), n, p);

  } else {

    // Statistic for each slice
    arma::uvec ind_tri = upper_tri_ind(n);
    for (arma::uword k = 0; k < M; k++) {

      // Compute Psi matrix
      arma::mat Psi = Psi_mat(X(arma::span::all, arma::span::all,
                                arma::span(k)), ind_tri, true, true, false);

      // Compute statistic
      Cn(k) = arma::as_scalar(sph_stat_Cai_Psi(Psi, n, p));

    }

  }

  // Sub-exponential regime
  if (regime == 1 | regime == 2) {

    Cn = p * Cn + 4 * std::log(n) - std::log(std::log(n));

  // Super-exponential regime
  } else if (regime == 3) {

    Cn = p * Cn + 4.0 * p / (p - 2.0) * std::log(n) - std::log(p);

  }

  return Cn;

}


//' @keywords internal
// [[Rcpp::export]]
arma::vec sph_stat_Cai_Psi(arma::mat Psi, arma::uword n, arma::uword p) {

  // Statistic
  arma::vec Cn = arma::max(arma::abs(Psi), 0).t();
  Cn = arma::log1p(-arma::square(Cn));

  return Cn;

}

