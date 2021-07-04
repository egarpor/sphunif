
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::plugins("cpp11")]]
#include <RcppArmadillo.h>
#include <Rmath.h>
#include <R.h>
using namespace Rcpp;

// Declaration for functions
arma::vec p_Kolmogorov(arma::vec, arma::uword K_Kolmogorov, bool alternating);
arma::vec d_Kolmogorov(arma::vec, arma::uword K_Kolmogorov, bool alternating);
arma::vec Gauss_Legen_nodes(double a, double b, arma::uword N);
arma::vec Gauss_Legen_weights(double a, double b, arma::uword N);
arma::vec d_chisq(arma::vec, arma::uword df);
arma::vec p_chisq(arma::vec, arma::uword df);

// Constants
const double half_M_PI = 0.5 * M_PI;
const double two_M_PI = 2.0 * M_PI;
const double inv_two_M_PI = 0.5 / M_PI;
const double sqrt_two_M_PI = std::sqrt(2 * M_PI);
const double sqrt_M_PI = std::sqrt(M_PI);
const double sqrt_eight = std::sqrt(8.0);
const double log_two = std::log(2.0);
const double const_1976_1 = 8.0 / (3 * sqrt_M_PI);
const double const_1976_2 = 8 * sqrt_eight / (9 * sqrt_M_PI);
const double two_thirds = 2.0 / 3.0;
const double four_thirds = 4.0 / 3.0;
const double three_inv_sqrt_eight = 3.0 / sqrt_eight;
const double sqrt_M_PI_sixth_one = std::sqrt(M_PI * M_PI / 6.0 - 1);
const double sqrt_inv_three = 1 / std::sqrt(3.0);
const double sd_Rao = two_M_PI *
  std::sqrt((2 - 5.0 / arma::datum::e) / arma::datum::e);


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec p_Kolmogorov(arma::vec x, arma::uword K_Kolmogorov = 25,
                       bool alternating = true) {
  
  // For 0 <= x <= 0.16, it happens that p_Kolmogorov(x) < 1e-19
  // The alternating series is adequate for evaluating upper tail
  // probabilities, but for 0.16 < x < 0.5 it may yield artificial increments
  // of the probability (if the truncation K_Kolmogorov is small and even,
  // e.g., K_Kolmogorov = 12) or negative values (K_Kolmogorov small and odd,
  // e.g, K_Kolmogorov = 11). The non-alternating version behaves better in
  // the lower tail for small K (e.g., K_Kolmogorov = 5), but converges slower
  // for the upper tail than the alternating series. With the truncation for
  // 0 <= x <= 0.16 and K_Kolmogorov = 25, both series yield an accuracy of
  // 3.6e-15 uniformly over [0, 10] discretized with step 1e-4. The
  // alternating implementation is slightly faster and precise, so it is the
  // default
  
  // Zero for the unfeasible values
  arma::vec F = arma::zeros(x.n_elem);
  arma::uvec ind = arma::find(x > 0.16);
  
  // Feasible values?
  if (ind.n_elem > 0) {
    
    if (alternating) {
      
      // Signs
      arma::rowvec m = arma::regspace(1, 1, K_Kolmogorov).t();
      arma::rowvec signs = -2 * m + 4 * arma::ceil(0.5 * m) - 1;
      
      // Series
      arma::mat S = arma::exp(-2 * arma::square(x.elem(ind) * m));
      S.each_row() %= signs;
      
      // Cdf F(x)
      F.elem(ind) = 1 - 2 * arma::sum(S, 1);
      
    } else {
      
      // (2 * m - 1)^2 * pi^2 / 8
      arma::rowvec m_odd_pi_sq_8 =
        arma::linspace(M_PI, (2 * K_Kolmogorov - 1) * M_PI, K_Kolmogorov).t();
      m_odd_pi_sq_8 = arma::square(m_odd_pi_sq_8) * 0.125;
      
      // Series
      arma::vec x1 = 1 / x.elem(ind);
      arma::mat S = arma::square(x1) * m_odd_pi_sq_8;
      
      // Cdf F(x)
      F.elem(ind) = sqrt_two_M_PI * arma::sum(exp(-S), 1) % x1;
      
    }
    
  }
  
  // F(x)
  return F;
  
}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec d_Kolmogorov(arma::vec x, arma::uword K_Kolmogorov = 25,
                       bool alternating = true) {
  
  // Same comment as for p_Kolmogorov holds
  
  // Zero for the unfeasible values
  arma::vec f = arma::zeros(x.n_elem);
  arma::uvec ind = arma::find(x > 0.16);
  
  // Feasible values?
  if (ind.n_elem > 0) {
    
    if (alternating) {
      
      // Signs
      arma::rowvec m = arma::regspace(1, 1, K_Kolmogorov).t();
      arma::rowvec signs = -2 * m + 4 * arma::ceil(0.5 * m) - 1;
      
      // Series
      m = arma::square(m);
      arma::mat S = arma::exp(-2 * arma::square(x.elem(ind)) * m) %
        (x.elem(ind) * m);
      S.each_row() %= signs;
      
      // Pdf f(x)
      f.elem(ind) = 8 * arma::sum(S, 1);
      
    } else {
      
      // (2 * m - 1)^2 * pi^2 / 8
      arma::rowvec m_odd_pi_sq_8 =
        arma::linspace(M_PI, (2 * K_Kolmogorov - 1) * M_PI, K_Kolmogorov).t();
      m_odd_pi_sq_8 = arma::square(m_odd_pi_sq_8) * 0.125;
      
      // Series
      arma::vec x2 = arma::square(1.0 / x.elem(ind));
      arma::mat S = x2 * m_odd_pi_sq_8;
      S = arma::exp(-S) % (2 * x2 * m_odd_pi_sq_8 - 1);
      S.each_col() %= x2;
      
      // Pdf f(x)
      f.elem(ind) = sqrt_two_M_PI * arma::sum(S, 1);
      
    }
    
  }
  
  // f(x)
  return f;
  
}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec p_cir_stat_Ajne(arma::vec x, arma::uword K_Ajne = 15) {
  
  // Zero for the unfeasible values
  arma::vec F = arma::zeros(x.n_elem);
  arma::uvec ind = arma::find(x > 0);
  
  // Feasible values?
  if (ind.n_elem > 0) {
    
    // Signs
    arma::rowvec m = arma::regspace(1, 1, K_Ajne).t();
    arma::rowvec signs = -2 * m + 4 * arma::ceil(0.5 * m) - 1;
    
    // Series
    m = M_PI * (2 * m - 1);
    arma::mat S = arma::exp(-0.5 * x.elem(ind) * arma::square(m));
    S.each_row() %= signs / m;
    
    // Cdf F(x)
    F.elem(ind) = 1 - 4 * arma::sum(S, 1);
    
    // Avoid potential problematic cases
    F = arma::clamp(F, 0, 1);
    
  }
  
  // F(x)
  return F;
  
}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec d_cir_stat_Ajne(arma::vec x, arma::uword K_Ajne = 15) {
  
  // Zero for the unfeasible values
  arma::vec f = arma::zeros(x.n_elem);
  arma::uvec ind = arma::find(x > 0);
  
  // Feasible values?
  if (ind.n_elem > 0) {
    
    // Signs
    arma::rowvec m = arma::regspace(1, 1, K_Ajne).t();
    arma::rowvec signs = -2 * m + 4 * arma::ceil(0.5 * m) - 1;
    
    // Series
    m = M_PI * (2 * m - 1);
    arma::mat S = arma::exp(-0.5 * x.elem(ind) * arma::square(m));
    S.each_row() %= signs % m;
    
    // Pdf f(x)
    f.elem(ind) = 2 * arma::sum(S, 1);
    
    // Avoid potential problematic cases
    f.elem(arma::find(f < 0)).zeros();
    
  }
  
  // f(x)
  return f;
  
}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec p_cir_stat_Bingham(arma::vec x) {
  
  return p_chisq(x, 2);
  
}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec d_cir_stat_Bingham(arma::vec x) {
  
  return d_chisq(x, 2);
  
}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec p_cir_stat_Greenwood(arma::vec x) {
  
  return arma::normcdf(x, 0, 2);
  
}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec d_cir_stat_Greenwood(arma::vec x) {
  
  return arma::normpdf(x, 0, 2);
  
}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec p_cir_stat_Gini(arma::vec x) {
  
  return arma::normcdf(x, 0, sqrt_inv_three);
  
}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec d_cir_stat_Gini(arma::vec x) {
  
  return arma::normpdf(x, 0, sqrt_inv_three);
  
}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec p_cir_stat_Gini_squared(arma::vec x) {
  
  return arma::normcdf(x, 0, 4);
  
}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec d_cir_stat_Gini_squared(arma::vec x) {
  
  return arma::normpdf(x, 0, 4);
  
}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec p_cir_stat_Hodges_Ajne2(arma::vec x, arma::uword n,
                                  bool asymp_std = false) {
  
  if (asymp_std) {
    
    return 1 - p_Kolmogorov(half_M_PI / x);
    
  } else {
    
    // j's
    x += 1;
    arma::uword M = arma::max(arma::floor((n - x) / (2.0 * x - n)));
    arma::rowvec j = arma::linspace(0, M, M + 1).t();
    
    // Binomial coefficients -- requires C++ 11
    arma::mat log_coef = (2.0 * x - n) * j;
    log_coef.each_col() += x;
    log_coef.transform([n](int m) {
      return R::lchoose(n, m);
    });
    
    // Add factor
    log_coef.each_col() += arma::log(2.0 * x - n) - log_two * (n - 1);
    return 1 - arma::sum(arma::exp(log_coef), 1);
    
  }
  
}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec p_cir_stat_Hodges_Ajne(arma::vec x, arma::uword n, bool exact = true,
                                 bool asymp_std = false) {
  
  // TODO: way slower than p_cir_stat_Hodges_Ajne2
  
  // Unstandardize evaluation points
  if (asymp_std) {
    
    x = 0.5 * (std::sqrt(n) * x + n);
    
  }
  
  // Zero and one for the unfeasible values
  arma::vec F = arma::zeros(x.n_elem);
  arma::uvec ind = arma::find(x >= floor(0.5 * n) && x < n);
  F.elem(arma::find(x >= n)).ones();
  
  // Feasible values?
  if (ind.n_elem > 0) {
    
    // Exact (discrete) or asymptotic distribution?
    if (exact) {
      
      // The + 1 because Ajne (1986) gives P[sup_\alpha N(\alpha) >= k], but
      // P[sup_\alpha N(\alpha) <= k] = 1 - P[sup_\alpha N(\alpha) >= k + 1]
      x.elem(ind) += 1;
      
      // Determine the minimum atom necessary to compute (from k_min to n)
      arma::uvec x_floor = arma::conv_to<arma::uvec>::from(
        arma::floor(x.elem(ind)));
      arma::uword k_min = arma::min(x_floor);
      arma::uword k_max = arma::max(x_floor);
      arma::vec k = arma::regspace(k_min, 1, k_max);
      
      // Largest upper limit in the sum (obtained for k_min)
      arma::uword j_max = std::floor((n - k_min) / (2.0 * k_min - n));
      
      // j's in the sum
      arma::rowvec j = arma::regspace(0, 1, j_max).t();
      
      // Log-binomial coefficients, a matrix of size k.n_elem x (j_max + 1)
      arma::mat log_coef = (2.0 * k - n) * j;
      log_coef.each_col() += k;
      log_coef.transform([n](int m) {
        return R::lchoose(n, m);
      });
      
      // Add factor
      log_coef.each_col() += arma::log(2.0 * k - n) - log_two * (n - 1);
      
      // Accumulated probability for k
      arma::vec F_k = 1 - arma::sum(arma::exp(log_coef), 1);
      
      // Obtain the corresponding values of F_k for x
      x_floor -= k_min;
      F.elem(ind) = F_k.elem(x_floor);
      
    } else {
      
      F.elem(ind) = 1 -
        p_Kolmogorov(half_M_PI * std::sqrt(n) / (2.0 * (x.elem(ind) + 1) - n));
      
    }
    
    // Avoid potential problematic cases
    F = arma::clamp(F, 0, 1);
    
  }
  
  // F(x)
  return F;
  
}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec d_cir_stat_Hodges_Ajne(arma::vec x, arma::uword n, bool exact = true,
                                 bool asymp_std = false) {
  
  // Zero for the unfeasible values
  arma::vec f = arma::zeros(x.n_elem);
  arma::uvec ind = arma::find(x >= floor(0.5 * n) && x <= n);
  f.elem(arma::find(x > n)).zeros();
  
  // Feasible values?
  if (ind.n_elem > 0) {
    
    if (exact) {
      
      // f(x) = F(x) - F(x^-)
      arma::vec F_dif = p_cir_stat_Hodges_Ajne(
        arma::join_vert(x.elem(ind), x.elem(ind) - 1.0), n, true, asymp_std);
      f.elem(ind) = F_dif.head(ind.n_elem) - F_dif.tail(ind.n_elem);
      
    } else {
      
      // Standardize evaluation points
      double c = half_M_PI;
      if (!asymp_std) {
        
        double inv_sqrt_n = 1.0 / std::sqrt(n);
        x.elem(ind) = (2.0 * x.elem(ind) - n) * inv_sqrt_n;
        c = M_PI * inv_sqrt_n;
        
      }
      
      // f(x)
      return c * d_Kolmogorov(half_M_PI / x.elem(ind)) /
        arma::square(x.elem(ind));
      
    }
    
  }
  
  // f(x)
  return f;
  
}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec p_cir_stat_Kuiper(arma::vec x, arma::uword n,
                            arma::uword K_Kuiper = 12, bool second_term = true,
                            bool Stephens = false) {

  // For 0 <= x <= 0.32 and 2 <= n <= 1e6, it happens that p_Kuiper(x) < 1e-15

  // Zero for the unfeasible values
  arma::vec F = arma::zeros(x.n_elem);
  arma::uvec ind = arma::find(x > 0.32);

  // Feasible values?
  if (ind.n_elem > 0) {

    // Add Stephens (1970) modification?
    if (Stephens) {

      x /= (1 + 0.155 / std::sqrt(n) + 0.24 / n);

    }

    // Series
    arma::rowvec m = arma::regspace(1, 1, K_Kuiper).t();
    arma::mat S = arma::square(x.elem(ind) * m);
    if (second_term) {

      S = (4 * S - 1 - four_thirds / std::sqrt(n) *
        (x.elem(ind) * arma::square(m)) % (4 * S - 3)) % arma::exp(-2 * S);

    } else {

      S = (4 * S - 1) % arma::exp(-2 * S);

    }

    // Cdf F(x)
    F.elem(ind) = 1 - 2 * arma::sum(S, 1);

    // Avoid potential problematic cases
    F = arma::clamp(F, 0, 1);

  }

  // F(x)
  return F;

}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec d_cir_stat_Kuiper(arma::vec x, arma::uword n,
                            arma::uword K_Kuiper = 12, bool second_term = true,
                            bool Stephens = false) {

  // Same comment as for p_cir_stat_Kuiper holds

  // Zero for the unfeasible values
  arma::vec f = arma::zeros(x.n_elem);
  arma::uvec ind = arma::find(x > 0.32);

  // Feasible values?
  if (ind.n_elem > 0) {

    // Add Stephens (1970) modification?
    double c = 1.0;
    if (Stephens) {

      c = 1.0 / (1 + 0.155 / std::sqrt(n) + 0.24 / n);
      x *= c;

    }

    // Series
    arma::rowvec m = arma::regspace(1, 1, K_Kuiper).t();
    arma::mat S = arma::square(x.elem(ind) * m);
    m = arma::square(m);
    if (second_term) {

      S = (-16 * arma::pow(x.elem(ind), 3) * arma::square(m) +
        12 * x.elem(ind) * m + four_thirds / std::sqrt(n) *
        arma::repmat(m, ind.n_elem, 1) %
        (16 * arma::square(S) - 24 * S + 3)) % arma::exp(-2 * S);

    } else {

      S = (-16 * arma::pow(x.elem(ind), 3) * arma::square(m) +
        12 * x.elem(ind) * m) % arma::exp(-2 * S);

    }

    // Pdf f(x)
    f.elem(ind) = -2 * arma::sum(S, 1) * c;

    // Avoid potential problematic cases
    f.elem(arma::find(f < 0)).zeros();

  }

  // f(x)
  return f;

}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec p_cir_stat_Log_gaps(arma::vec x, bool abs_val = true) {
  
  if (abs_val) {
    
    arma::vec F = 2 * arma::normcdf(x, 0, sqrt_M_PI_sixth_one) - 1;
    F.elem(arma::find(x < 0)).zeros();
    return F;
    
  } else {
    
    return arma::normcdf(x, 0, sqrt_M_PI_sixth_one);
    
  }
  
}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec d_cir_stat_Log_gaps(arma::vec x, bool abs_val = true) {
  
  if (abs_val) {
    
    arma::vec f = 2 * arma::normpdf(x, 0, sqrt_M_PI_sixth_one);
    f.elem(arma::find(x < 0)).zeros();
    return f;
    
  } else {
    
    return arma::normpdf(x, 0, sqrt_M_PI_sixth_one);
    
  }
  
}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec p_cir_stat_Max_uncover(arma::vec x) {
  
  return arma::exp(-arma::exp(-x));
  
}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec d_cir_stat_Max_uncover(arma::vec x) {
  
  return arma::exp(-(x + arma::exp(-x)));
  
}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec p_cir_stat_Num_uncover(arma::vec x) {
  
  return arma::normcdf(x, 0, 1);
  
}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec d_cir_stat_Num_uncover(arma::vec x) {
  
  return arma::normpdf(x, 0, 1);
  
}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec p_cir_stat_Pycke(arma::vec x) {
  
  return arma::exp(-arma::exp(-0.5 * (x + 2 * arma::datum::euler)));
  
}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec d_cir_stat_Pycke(arma::vec x) {
  
  x = 0.5 * (x + 2 * arma::datum::euler);
  return 0.5 * arma::exp(-(x + arma::exp(-x)));
  
}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec p_cir_stat_Vacancy(arma::vec x) {
  
  return arma::normcdf(x, 0, 1);
  
}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec d_cir_stat_Vacancy(arma::vec x) {
  
  return arma::normpdf(x, 0, 1);
  
}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec p_cir_stat_Watson(arma::vec x, arma::uword n = 0,
                            arma::uword K_Watson = 25, bool Stephens = false) {

  // Zero for the unfeasible values
  arma::vec F = arma::zeros(x.n_elem);
  arma::uvec ind = arma::find(x > 0);

  // Feasible values?
  if (ind.n_elem > 0) {

    // Add Stephens (1970) modification?
    if (Stephens) {

      x /= (1 + 0.155 / std::sqrt(n) + 0.24 / n);

    }

    // Cdf F(x)
    F.elem(ind) = p_Kolmogorov(arma::sqrt(x.elem(ind)) * M_PI, K_Watson, true);

  }

  // F(x)
  return F;

}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec d_cir_stat_Watson(arma::vec x, arma::uword n = 0,
                            arma::uword K_Watson = 25, bool Stephens = false) {

  // Zero for the unfeasible values
  arma::vec f = arma::zeros(x.n_elem);
  arma::uvec ind = arma::find(x > 0);

  // Feasible values?
  if (ind.n_elem > 0) {

    // Add Stephens (1970) modification?
    double c = 1.0;
    if (Stephens) {

      c /= 1 + 0.8 / n;
      x *= c;
      x -= -0.1 / n + 0.1 / (n * n);

    }

    // Pdf f(x)
    x = arma::sqrt(x);
    f.elem(ind) = d_Kolmogorov(x.elem(ind) * M_PI, K_Watson, true) /
      x.elem(ind) * (half_M_PI * c);

  }

  // f(x)
  return f;

}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec p_cir_stat_Watson_1976(arma::vec x, arma::uword K_Watson_1976 = 8,
                                 arma::uword N = 40) {

  // Zero and one for the unfeasible values
  arma::vec F = arma::zeros(x.n_elem);
  arma::uvec ind = arma::find(x > 0 && x < 2.4);
  F.elem(arma::find(x >= 2.4)).ones();

  // Feasible values?
  if (ind.n_elem > 0) {

    // Zeros of Ai(-x) = 0
    arma::vec alpha = {2.338107410459767038489, 4.087949444130970616637,
                       5.520559828095551059128, 6.786708090071758998780,
                       7.944133587120853123138, 9.022650853340980380158,
                       10.04017434155808593059, 11.00852430373326289324,
                       11.93601556323626251701, 12.82877675286575720041,
                       13.69148903521071792830, 14.52782995177533498207,
                       15.34075513597799685715, 16.13268515694577143935,
                       16.90563399742994262704, 17.66130010569705750925,
                       18.40113259920711541586, 19.12638047424695214412,
                       19.83812989172149970095, 20.53733290767756635998,
                       21.22482994364209695520, 21.90136759558513070741,
                       22.56761291749650283146, 23.22416500112168106132,
                       23.87156445553591856712};
    if (K_Watson_1976 > 25) {

      stop("Only K_Watson_1976 smaller than 25 are implemented");

    } else {

      alpha = alpha.head(K_Watson_1976);

    }

    // Zeros of J(x, 1/3) + J(x, -1/3) = Ai(-(3/2 * x)^(2/3)) = 0
    alpha = two_thirds * arma::pow(alpha, 1.5);

    // Gauss--Legendre numerical integration
    arma::vec x_k = Gauss_Legen_nodes(0, M_PI, N);
    arma::vec w_k = Gauss_Legen_weights(0, M_PI, N);

    // g(theta) for the integrand of v
    arma::rowvec g_theta = (arma::square(arma::sin(two_thirds * x_k)) %
      arma::sin(x_k / 3.0) / arma::pow(arma::sin(x_k), 3)).t();
    arma::mat log_g_theta = arma::repmat(arma::log(g_theta), ind.n_elem, 1);

    // Loop on the sum, vectorize on the evaluation of the integrand
    arma::vec x_m_inv = arma::zeros(ind.n_elem);
    arma::vec x_inv = three_inv_sqrt_eight / x.elem(ind);
    arma::vec v = arma::zeros(ind.n_elem);
    for (arma::uword m = 0; m < K_Watson_1976; m++) {

      // Evaluate integrand
      x_m_inv = alpha(m) * x_inv;
      arma::mat integrand = arma::exp(-arma::square(x_m_inv) * g_theta +
        log_g_theta);

      // Update sum by integral
      v += (integrand * w_k) % arma::pow(x_m_inv, 3) / alpha(m);

    }

    // Cdf F(x)
    F.elem(ind) = const_1976_1 * v;

    // Avoid potential problematic cases
    F = arma::clamp(F, 0, 1);

  }

  // F(x)
  return F;

}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec d_cir_stat_Watson_1976(arma::vec x, arma::uword K_Watson_1976 = 8) {

  // Zero for the unfeasible values
  arma::vec f = arma::zeros(x.n_elem);
  arma::uvec ind = arma::find(x > 0 && x < 2.4);

  // Feasible values?
  if (ind.n_elem > 0) {

    // Zeros of Ai(-x) = 0
    arma::vec alpha = {2.338107410459767038489, 4.087949444130970616637,
                       5.520559828095551059128, 6.786708090071758998780,
                       7.944133587120853123138, 9.022650853340980380158,
                       10.04017434155808593059, 11.00852430373326289324,
                       11.93601556323626251701, 12.82877675286575720041,
                       13.69148903521071792830, 14.52782995177533498207,
                       15.34075513597799685715, 16.13268515694577143935,
                       16.90563399742994262704, 17.66130010569705750925,
                       18.40113259920711541586, 19.12638047424695214412,
                       19.83812989172149970095, 20.53733290767756635998,
                       21.22482994364209695520, 21.90136759558513070741,
                       22.56761291749650283146, 23.22416500112168106132,
                       23.87156445553591856712};
    if (K_Watson_1976 > 25) {

      stop("Only K_Watson_1976 smaller than 25 are implemented");

    } else {

      alpha = alpha.head(K_Watson_1976);

    }

    // Zeros of J(x, 1/3) + J(x, -1/3) = Ai(-(3/2 * x)^(2/3)) = 0
    alpha = two_thirds * arma::pow(alpha, 1.5);

    // Gauss--Legendre numerical integration
    arma::vec x_k = Gauss_Legen_nodes(0, M_PI, 40);
    arma::vec w_k = Gauss_Legen_weights(0, M_PI, 40);

    // g(theta) for the integrand of v
    arma::rowvec g_theta = (arma::square(arma::sin(two_thirds * x_k)) %
      arma::sin(x_k / 3.0) / arma::pow(arma::sin(x_k), 3)).t();
    arma::mat log_g_theta = arma::repmat(arma::log(g_theta), ind.n_elem, 1);

    // Loop on the sum, vectorize on the evaluation of the integrand
    arma::vec x_m_inv_2 = arma::zeros(ind.n_elem);
    arma::vec v = arma::zeros(ind.n_elem);
    for (arma::uword m = 0; m < K_Watson_1976; m++) {

      // Evaluate integrand
      x_m_inv_2 = arma::square(3 * alpha(m) / (sqrt_eight * x.elem(ind)));
      arma::mat integrand = arma::exp(-x_m_inv_2 * g_theta + log_g_theta) %
        (-3 + x_m_inv_2 * g_theta * 2);

      // Update sum by integral
      v += (integrand * w_k) % arma::square(x_m_inv_2) / (alpha(m) * alpha(m));

    }

    // Pdf f(x)
    f.elem(ind) = const_1976_2 * v;

    // Avoid potential problematic cases
    f.elem(arma::find(f < 0)).zeros();

  }

  // f(x)
  return f;

}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec p_cir_stat_Range(arma::vec x, arma::uword n, bool max_gap = true) {

  // Standardize x
  x *= inv_two_M_PI;

  // Transform to minimum gap
  if (max_gap) {

    x = 1 - x;

  }

  // The maximum gap can not be smaller than 2 * pi / n (best separation).
  // Then, the minimum gap can not be larger than 2 * pi * (n - 1) / n

  // Zero and one for the unfeasible values
  arma::vec F = arma::zeros(x.n_elem);
  arma::uvec ind = arma::find(x > 0 && x < (1 - 1.0 / n));
  F.elem(arma::find(x >= (1 - 1.0 / n))).ones();

  // Feasible values?
  if (ind.n_elem > 0) {

    // m's
    arma::uword M = arma::max(arma::floor(1.0 / (1 - x.elem(ind))));
    M = std::min(M, n);
    arma::rowvec m = arma::regspace(1, 1, M).t();

    // Signs
    arma::rowvec signs = -2 * m + 4 * arma::ceil(0.5 * m) - 1;

    // Log-binomial coefficients -- requires C++ 11
    arma::rowvec log_coef = m;
    log_coef.transform([n](int m) {
      return R::lchoose(n, m);
    });

    // We compute the log-terms of the sum with killing of negative factors
    // (1 - m * (1 - t / (2 * pi))) that appear for m larger than
    // min(floor(2 * pi / (2 * pi - t)), n) (because M is not dependent on x).
    // The killing is done by excluding the non-finite terms, as those are
    // the logarithms of non-positive values

    // Log-terms of the sum
    arma::mat log_S = (n - 1) * arma::log(1 - (1 - x.elem(ind)) * m);
    log_S.each_row() += log_coef;

    // Kill terms later with exp(log_S) = 0
    log_S.elem(arma::find_nonfinite(log_S)).fill(-arma::datum::inf);

    // Improve the numerical stability of the sum by shifting the
    // exponents with the maximum log-term
    arma::vec log_S_max = arma::max(log_S, 1);
    log_S.each_col() -= log_S_max;

    // Exponentiate the log-terms and add signs
    log_S = arma::exp(log_S);
    log_S.each_row() %= signs;

    // Sum and multiply reshift exponentials
    F.elem(ind) = arma::sum(log_S, 1) % arma::exp(log_S_max);

    // Transform for maximum gap
    if (max_gap) {

      F = 1 - F;

    }

    // Avoid potential problematic cases
    F = arma::clamp(F, 0, 1);

  }

  // F(x)
  return F;

}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec d_cir_stat_Range(arma::vec x, arma::uword n, bool max_gap = true) {

  // Standardize x
  x *= inv_two_M_PI;

  // Transform to minimum gap
  if (max_gap) {

    x = 1 - x;

  }

  // Same comment as for p_cir_stat_Range holds

  // Zero for the unfeasible values
  arma::vec f = arma::zeros(x.n_elem);
  arma::uvec ind = arma::find(x > 0 && x < (1 - 1.0 / n));

  // Feasible values?
  if (ind.n_elem > 0) {

    // m's
    arma::uword M = arma::max(arma::floor(1 / (1 - x.elem(ind))));
    M = std::min(M, n);
    arma::rowvec m = arma::regspace(1, 1, M).t();

    // Signs
    arma::rowvec signs = -2 * m + 4 * arma::ceil(0.5 * m) - 1;

    // Log-binomial coefficients -- requires C++ 11
    arma::rowvec log_coef = m;
    log_coef.transform([n](int m) {
      return R::lchoose(n, m);
    });

    // Log-terms of the sum
    arma::mat log_S = (n - 2) * arma::log(1 - (1 - x.elem(ind)) * m);
    log_S.each_row() += log_coef + arma::log((n - 1) * m * inv_two_M_PI);

    // Kill terms later with exp(log_S) = 0
    log_S.elem(arma::find_nonfinite(log_S)).fill(-arma::datum::inf);

    // Improve the numerical stability of the sum by shifting the
    // exponents with the maximum log-term
    arma::vec log_S_max = arma::max(log_S, 1);
    log_S.each_col() -= log_S_max;

    // Exponentiate the log-terms and add signs
    log_S = arma::exp(log_S);
    log_S.each_row() %= signs;

    // Sum and multiply reshift exponentials
    f.elem(ind) = arma::sum(log_S, 1) % arma::exp(log_S_max);

    // Avoid potential problematic cases
    f.elem(arma::find(f < 0)).zeros();

  }

  // f(x)
  return f;

}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec p_cir_stat_Rao(arma::vec x) {

  return arma::normcdf(x, 0, sd_Rao);

}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec d_cir_stat_Rao(arma::vec x) {

  return arma::normpdf(x, 0, sd_Rao);

}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec p_cir_stat_Rayleigh(arma::vec x) {

  return p_chisq(x, 2);

}


//' @rdname cir_stat_distr
//' @export
// [[Rcpp::export]]
arma::vec d_cir_stat_Rayleigh(arma::vec x) {

  return d_chisq(x, 2);

}

