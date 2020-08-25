
// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// Declaration for functions
arma::vec Gauss_Legen_nodes(double a, double b, arma::uword N);
arma::vec Gauss_Legen_weights(double a, double b, arma::uword N);
arma::vec d_proj_unif(arma::vec x, arma::uword p, bool log);
arma::vec p_proj_unif(arma::vec x, arma::uword p, bool log);
arma::vec t_inv_sqrt_one(arma::vec t);

// Constants
const double inv_PI = 1.0 / PI;
const double inv_two_PI = 0.5 / PI;

//' @title Surface area of the intersection of two hyperspherical caps
//'
//' @description Computation of
//' \deqn{A_x(\theta_{ij}) := \frac{1}{\omega_p}
//' \int_{S^{p - 1}} 1_{\{{\bf X}_i'\boldsymbol\gamma \le x,
//' {\bf X}_j'\boldsymbol\gamma \le x\}}\,\mathrm{d}\boldsymbol\gamma,}{
//' A_x(\theta_{ij}) := \frac{1}{\omega_p} \int_{S^{p - 1}}
//' 1_{X_i'\gamma \le x, X_j'\gamma \le x} d\gamma,}
//' where \eqn{\theta_{ij} := \cos^{-1}({\bf X}_i'{\bf X}_j)
//' \in [0, \pi]}{\theta_{ij} := \cos^{-1}(X_i'X_j) \in [0, \pi]},
//' \eqn{x \in [-1, 1]}, and \eqn{\omega_{p}} is the surface area of
//' \eqn{S^{p - 1}}. \eqn{A_x(\theta_{ij})} is the proportion of surface area
//' of \eqn{S^{p - 1}} covered by the intersection of two hyperspherical caps
//' centered at \eqn{{\bf X}_i}{X_i} and \eqn{{\bf X}_j}{X_j} and with
//' common solid angle \eqn{\pi - \cos^{-1}(x)}.
//'
//' @inheritParams Gegenbauer
//' @param x vector with values in \eqn{[-1, 1]}.
//' @inheritParams r_unif
//' @param N number of points used in the
//' \link[=Gauss_Legen_nodes]{Gauss-Legendre quadrature}. Defaults to
//' \code{160}.
//' @param as_matrix return a matrix with the values of \eqn{A_x(\theta)} on
//' the grid formed by \code{theta} and \code{x}? If \code{FALSE},
//' \eqn{A_x(\theta)} is evaluated on \code{theta} and \code{x} if they equal
//' in size. Defaults to \code{TRUE}.
//' @return A matrix of size \code{c(length(theta), length(x))} containing the
//' evaluation of \eqn{A_x(\theta)} if \code{as_matrix = TRUE}. Otherwise,
//' a vector of size \code{c(length(theta)} if \code{theta} and \code{x} equal
//' in size.
//' @details
//' See García-Portugués et al. (2020) for more details about the
//' \eqn{A_x(\theta)} function.
//' @references
//' García-Portugués, E., Navarro-Esteban, P., Cuesta-Albertos, J. A. (2020)
//' On a projection-based class of uniformity tests on the hypersphere.
//' \emph{arXiv:2008.09897}. \url{https://arxiv.org/abs/2008.09897}
//' @examples
//' # Plot A_x(theta) for several dimensions and x's
//' A_lines <- function(x, th = seq(0, pi, l = 200)) {
//'
//'   plot(th, A_theta_x(theta = th, x = x, p = 2), type = "l",
//'        col = 1, ylim = c(0, 1.25), main = paste("x =", x),
//'        ylab = expression(A[x](theta)),
//'        xlab = expression(theta), axes = FALSE)
//'   axis(1, at = c(0, pi / 4, pi / 2, 3 * pi / 4, pi),
//'        labels = expression(0, pi / 4, pi / 2, 3 * pi / 4, pi))
//'   axis(2); box()
//'   abline(h = c(0, 1), lty = 2)
//'   lines(th, A_theta_x(theta = th, x = x, p = 3), col = 2)
//'   lines(th, A_theta_x(theta = th, x = x, p = 4), col = 3)
//'   lines(th, A_theta_x(theta = th, x = x, p = 5), col = 4)
//'   legend("top", lwd = 2, legend = paste("p =", 2:5),
//'          col = 1:4, cex = 0.75, horiz = TRUE)
//'
//' }
//' old_par <- par(mfrow = c(2, 3))
//' A_lines(x = -0.75)
//' A_lines(x = -0.25)
//' A_lines(x = 0)
//' A_lines(x = 0.25)
//' A_lines(x = 0.5)
//' A_lines(x = 0.75)
//' par(old_par)
//'
//' # As surface of (theta, x) for several dimensions
//' A_surf <- function(p, x = seq(-1, 1, l = 201), th = seq(0, pi, l = 201)) {
//'
//'   col <- c("white", viridisLite::viridis(20))
//'   breaks <- c(-1, seq(1e-15, 1, l = 21))
//'   A <- A_theta_x(theta = th, x = x, p = p)
//'   image(th, x, A, main = paste("p =", p), col = col, breaks = breaks,
//'         xlab = expression(theta), axes = FALSE)
//'   axis(1, at = c(0, pi / 4, pi / 2, 3 * pi / 4, pi),
//'        labels = expression(0, pi / 4, pi / 2, 3 * pi / 4, pi))
//'   axis(2); box()
//'   contour(th, x, A, levels = breaks, add = TRUE)
//'
//' }
//' old_par <- par(mfrow = c(2, 2))
//' A_surf(p = 2)
//' A_surf(p = 3)
//' A_surf(p = 4)
//' A_surf(p = 5)
//' par(old_par)
//'
//' # No matrix return
//' th <- seq(0, pi, l = 5)
//' x <- seq(-1, 1, l = 5)
//' diag(A_theta_x(theta = th, x = x, p = 2))
//' A_theta_x(theta = th, x = x, p = 2, as_matrix = FALSE)
//' @export
// [[Rcpp::export]]
arma::mat A_theta_x(arma::vec theta, arma::vec x, arma::uword p,
                    arma::uword N = 160, bool as_matrix = true) {

  // Check if compatible lengths
  if (!as_matrix && (theta.n_elem != x.n_elem)) {

    stop("Incompatible lenghts of theta and x, consider as_matrix = TRUE");

  }

  // Identify negative x's
  arma::uvec x_neg = find(x < 0);

  // Overwrite x by its absolute value since signs are no longer needed
  x = abs(x);

  // Compute the cheap explicit part
  arma::vec exp_part = 2 * p_proj_unif(x, p, false) - 1;

  // Result matrix
  arma::mat A(theta.n_elem, as_matrix ? x.n_elem : 1);

  // Evaluate A_x(theta) on the grid of theta and x?
  if (!as_matrix) {

    // Initialize A with the explicit part
    A = exp_part;

    // Special case for dimension p = 2
    if (p == 2) {

      // Explicit form
      arma::vec acos = arma::acos(x);
      A = 1 - (acos + arma::min(0.5 * theta, acos)) * inv_PI;

    } else {

      // Integral part
      arma::uvec ind_int = arma::find((0 <= theta) %
                                      (theta < 2 * arma::acos(x)));
      if (ind_int.n_elem > 0) {

        // Loop on data, the kind of computation does not allow
        // vectorized benefits
        for (arma::uword j = 0; j < ind_int.n_elem; j++) {

          // Data index
          arma::uword i = ind_int(j);

          // i-th datums
          double x_i = x(i);
          double th_i = theta(i);

          // Nodes and weights for Gauss--Legendre quadrature
          arma::vec t_k = Gauss_Legen_nodes(0, x_i, N);
          arma::vec w_k = Gauss_Legen_weights(0, x_i, N);

          // Integrand multiplied by weights
          arma::vec integrand = p_proj_unif(t_inv_sqrt_one(t_k) *
                                            std::tan(0.5 * th_i), p - 1, false);
          integrand %= w_k % d_proj_unif(t_k, p, false);

          // Approximate integral
          A(i) = 0.5 - th_i * inv_two_PI;
          A(i) += 2 * arma::accu(integrand);

        }

      }

    }

    // Negative x? If so, apply symmetry
    if (x_neg.n_elem > 0) {

      A(x_neg) -= exp_part(x_neg);

    }

  } else {

    // Initialize A with the explicit part
    A = arma::repmat(exp_part.t(), theta.n_elem, 1);

    // Special case for dimension p = 2
    if (p == 2) {

      // Explicit form
      arma::mat acos = arma::repmat(arma::acos(x.as_row()), theta.n_elem, 1);
      A = 1 - (acos +
        arma::min(arma::repmat(0.5 * theta, 1, x.n_elem), acos)) * inv_PI;

    } else {

      // Loop on x
      for (arma::uword i = 0; i < x.n_elem; i++) {

        // i-th datum in x
        double x_i = x(i);

        // i-th column of A
        arma::vec A_i = A.col(i);

        // Integral part
        arma::uvec ind_int = arma::find((0 <= theta) %
                                        (theta < 2 * std::acos(x_i)));
        if (ind_int.n_elem > 0) {

          // Nodes and weights for Gauss--Legendre quadrature
          arma::vec t_k = Gauss_Legen_nodes(0, x_i, N);
          arma::vec w_k = Gauss_Legen_weights(0, x_i, N);

          // Integrand multiplied by weights
          arma::mat integrand = arma::reshape(
            p_proj_unif(arma::vectorise(t_inv_sqrt_one(t_k) *
              arma::tan(0.5 * theta(ind_int).t())), p - 1, false),
              N, ind_int.n_elem);
          integrand.each_col() %= w_k % d_proj_unif(t_k, p, false);

          // Approximate integral
          A_i(ind_int) = 0.5 - theta(ind_int) * inv_two_PI;
          A_i(ind_int) += 2 * arma::sum(integrand, 0);

        }

        // Save the i-th column
        A.col(i) = A_i;

      }

    }

    // Negative x? If so, apply symmetry
    if (x_neg.n_elem > 0) {

      A.cols(x_neg) -= arma::repmat(exp_part(x_neg).t(), theta.n_elem, 1);

    }

  }

  // Value A_x(theta)
  return A;

}
