

#' @title Type-7 quantiles of a vector
#'
#' @description Computation of the type-7 quantiles of a vector, avoiding
#' sorting if observations are already sorted.
#'
#' @param x_sorted a vector \bold{sorted increasingly} and \bold{free of}
#' \code{NA}s.
#' @param x a vector, possibly unsorted or containing \code{NA}s.
#' @param probs a vector of probabilities in \eqn{[0, 1]}.
#' @return A vector of size \code{length(probs)} with the type-7 quantiles.
#' @section Warning:
#' \code{quantile_sorted} does \bold{not} check its assumptions on
#' \code{x_sorted}, and it returns spurious results if \code{x_sorted} is
#' unsorted or contains \code{NA}s. Use \code{quantile_col} if these
#' assumptions cannot be guaranteed.
#' @details
#' \code{quantile_sorted} reproduces \code{quantile(x_sorted, probs,
#' type = 7, names = FALSE)} bit-for-bit, as it mirrors the interpolation in
#' \code{stats:::quantile.default} (including its tie guard), but without
#' re-sorting the vector.
#'
#' \code{quantile_col} takes the \code{quantile_sorted} shortcut only when
#' \code{x} is sorted increasingly and free of \code{NA}s, and delegates to
#' \code{\link[stats]{quantile}} otherwise.
#' @keywords internal
quantile_sorted <- function(x_sorted, probs) {

  n <- length(x_sorted)
  index <- 1 + max(n - 1, 0) * probs
  lo <- floor(index)
  hi <- ceiling(index)
  qs <- x_sorted[lo]
  i <- which(index > lo & x_sorted[hi] != qs)
  h <- (index - lo)[i]
  qs[i] <- (1 - h) * qs[i] + h * x_sorted[hi[i]]
  qs

}


#' @rdname quantile_sorted
quantile_col <- function(x, probs) {

  # anyNA() must be checked first, as is.unsorted() returns NA if there are NAs
  if (!anyNA(x) && !is.unsorted(x)) {

    quantile_sorted(x_sorted = x, probs = probs)

  } else {

    quantile(x, probs = probs, na.rm = TRUE, names = FALSE)

  }

}
