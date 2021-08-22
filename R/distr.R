

#' @rdname proj_unif
#' @export
r_proj_unif <- function(n, p) {

  # Check p
  if (p <= 1) {

    stop("p must be >= 2.")

  }

  # The projection is the square of a B(1 / 2, (p - 1)/2)
  return(sqrt(rbeta(n = n, shape1 = 0.5, shape2 = (p - 1) / 2)) *
           sample(c(-1, 1), size = n, replace = TRUE))

}
