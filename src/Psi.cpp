
// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// Declaration for functions
arma::uvec upper_tri_ind(arma::uword n);


//' @title Shortest angles matrix
//'
//' @description Efficient computation of the shortest angles matrix
//' \eqn{\boldsymbol\Psi}{\Psi}, defined as
//' \deqn{\Psi_{ij}:=\cos^{-1}({\bf X}_i'{\bf X}_j),\quad
//' i,j=1,\ldots,n,}{\Psi_{ij} = \cos^{-1}(X_i'X_j), i, j = 1, \ldots, n,}
//' for a sample \eqn{{\bf X}_1,\ldots,{\bf X}_n\in S^{p-1}:=\{{\bf x}\in
//' R^p:||{\bf x}||=1\}}{X_1, \ldots, X_n \in
//' S^{p - 1} := \{x \in R^p : ||x|| = 1\}}, \eqn{p\ge 2}.
//'
//' For a circular sample \eqn{\Theta_1, \ldots, \Theta_n \in [0, 2\pi)},
//' \eqn{\boldsymbol\Psi}{\Psi} can be expressed as
//' \deqn{\Psi_{ij}=\pi-|\pi-|\Theta_i-\Theta_j||,\quad
//' i,j=1,\ldots,n.}{\Psi_{ij}=\pi-|\pi-|\Theta_i-\Theta_j||, i,j=1,\ldots,n.}
//'
//' @param data an array of size \code{c(n, p, M)} containing the Cartesian
//' coordinates of \code{M} samples of size \code{n} of directions on
//' \eqn{S^{p-1}}. Alternatively if \code{p = 2}, an array of size
//' \code{c(n, 1, M)} containing the angles on \eqn{[0, 2\pi)} of the \code{M}
//' circular samples of size \code{n} on \eqn{S^{1}}. Must not contain
//' \code{NA}'s.
//' @param ind_tri if \code{use_ind_tri = TRUE}, the vector of 0-based indexes
//' provided by \code{upper_tri_ind(n)}, which allows to extract the upper
//' triangular part of the matrix \eqn{\boldsymbol\Psi}{\Psi}. See the examples.
//' @param use_ind_tri use the already computed vector index \code{ind_tri}? If
//' \code{FALSE} (default), \code{ind_tri} is computed internally.
//' @param scalar_prod return the scalar products
//' \eqn{{\bf X}_i'{\bf X}}{X_i'X_j} instead of the shortest angles? Only taken
//' into account for data in \emph{Cartesian} form. Defaults to
//' \code{FALSE}.
//' @param angles_diff return the (unwrapped) angles difference
//' \eqn{\Theta_i-\Theta_j} instead of the shortest angles? Only taken into
//' account for data in \emph{angular} form. Defaults to \code{FALSE}.
//' @param n sample size, used to determine the index vector that gives the
//' upper triangular part of \eqn{\boldsymbol\Psi}{\Psi}.
//' @return
//'
//' \itemize{
//'   \item \code{Psi_mat}: a matrix of size
//'   \code{c(n * (n - 1) / 2, M)} containing, for each column, the vector
//'   half of \eqn{\boldsymbol\Psi}{\Psi} for each of the \code{M} samples.
//'   \item \code{upper_tri_ind}: a matrix of size \code{n * (n - 1) / 2}
//'   containing the 0-based linear indexes for extracting the upper triangular
//'   matrix of a matrix of size \code{c(n, n)}, diagonal excluded, assuming
//'   column-major order.
//' }
//' @section Warning:
//' Be careful on avoiding the next bad usages of \code{Psi_mat}, which will
//' produce spurious results:
//' \itemize{
//'   \item The directions in \code{data} do \emph{not} have unit norm when
//'   Cartesian coordinates are employed.
//'   \item The entries of \code{data} are \emph{not} in \eqn{[0, 2\pi)} when
//'   polar coordinates are employed.
//'   \item \code{ind_tri} is a vector of size \code{n * (n - 1) / 2} that
//'   does \emph{not} contain the indexes produced by \code{upper_tri_ind(n)}.
//' }
//' @examples
//' # Shortest angles
//' n <- 5
//' X <- r_unif_sph(n = n, p = 2, M = 2)
//' Theta <- X_to_Theta(X)
//' dim(Theta) <- c(n, 1, 2)
//' Psi_mat(X)
//' Psi_mat(Theta)
//'
//' # Precompute ind_tri
//' ind_tri <- upper_tri_ind(n)
//' Psi_mat(X, ind_tri = ind_tri, use_ind_tri = TRUE)
//'
//' # Compare with R
//' A <- acos(tcrossprod(X[, , 1]))
//' ind <- upper.tri(A)
//' A[ind]
//'
//' # Reconstruct matrix
//' Psi_vec <- Psi_mat(Theta[, , 1, drop = FALSE])
//' Psi <- matrix(0, nrow = n, ncol = n)
//' Psi[upper.tri(Psi)] <- Psi_vec
//' Psi <- Psi + t(Psi)
//' @name Psi


//' @rdname Psi
//' @export
// [[Rcpp::export]]
arma::mat Psi_mat(arma::cube data, arma::uvec ind_tri = 0,
                  bool use_ind_tri = false, bool scalar_prod = false,
                  bool angles_diff = false) {

  // Create vec_Psi matrix
  arma::uword n = data.n_rows;
  arma::uword p = data.n_cols;
  arma::uword M = data.n_slices;
  arma::mat vec_Psi = arma::zeros(0.5 * n * (n - 1), M);

  // Upper triangular matrix index
  if (!use_ind_tri) {

    ind_tri.set_size(0.5 * n * (n - 1));
    ind_tri = upper_tri_ind(n);

  }

  // Angles
  if (p == 1) {

    // Process first slice
    arma::mat Psi = arma::repmat(data.slice(0), 1, n);
    Psi -= Psi.t();
    vec_Psi.col(0) = Psi.elem(ind_tri);

    // Loop on the remaining slices reusing ind
    for (arma::uword j = 1; j < M; j++) {

      Psi = arma::repmat(data.slice(j), 1, n);
      Psi -= Psi.t();
      vec_Psi.col(j) = Psi.elem(ind_tri);

    }

    // Angles difference or shortest angle?
    if (angles_diff) {

      return vec_Psi;

    } else {

      return M_PI - arma::abs(M_PI - arma::abs(vec_Psi));

    }

  // Cartesian coordinates
  } else {

    // Process first slice
    arma::mat Psi = arma::trimatu(data.slice(0) * data.slice(0).t(), 1);
    vec_Psi.col(0) = Psi.elem(ind_tri);

    // Loop on the remaining slices reusing ind
    for (arma::uword j = 1; j < M; j++) {

      Psi = arma::trimatu(data.slice(j) * data.slice(j).t(), 1);
      vec_Psi.col(j) = Psi.elem(ind_tri);

    }

    // Return scalar product
    if (scalar_prod) {

      return vec_Psi;

    // Return shortest angles
    } else {

      return arma::acos(arma::clamp(vec_Psi, -1.0, 1.0));

    }

  }

}


//' @rdname Psi
//' @export
// [[Rcpp::export]]
arma::uvec upper_tri_ind(arma::uword n) {

  return arma::find(arma::trimatu(arma::ones(n, n), 1));

}


//' @title Sort the columns of a matrix
//'
//' @description Convenience functions to sort the columns of a matrix in an
//' increasing way.
//'
//' @param A a matrix of arbitrary dimensions.
//' @return A matrix with the same dimensions as \code{A} such that each
//' column is sorted increasingly (for \code{sort_each_col}) or contains the
//' sorting indexes (for \code{sort_index_each_col}).
//' \code{sort_index_each_col}.
//' @keywords internal
// [[Rcpp::export]]
arma::mat sort_each_col(arma::mat A) {

  return arma::sort(A);

}


//' @rdname sort_each_col
//' @keywords internal
// [[Rcpp::export]]
arma::umat sort_index_each_col(arma::mat A) {

  // Loop on columns and extract sort_index
  arma::umat I(arma::size(A));
  for (arma::uword j = 0; j < A.n_cols; j ++) {

    I.col(j) = arma::sort_index(A.col(j));

  }
  return(I + 1);

}
