
#' @title Evaluation of \eqn{\log(e^{-x} \mathcal{I}_{\nu}(x))} with asymptotics
#'
#' @description Computes the logarithm of the scaled modified Bessel function
#' of the first kind with automatic handling of asymptotic cases.
#'
#' @param nu a scalar with the order of the Bessel function.
#' @param x vector with evaluation points for the Bessel function.
#' @details For \code{x} larger than \code{1e4}, an asymptotic approximation
#' with \code{\link[Bessel]{besselIasym}} is used. For \code{nu} larger than
#' \code{100}, an asymptotic approximation with
#' \code{\link[Bessel]{besselI.nuAsym}} is used.
#' @return A vector of size \code{length(x)} with the evaluated function.
#' @examples
#' curve(log_besselI_scaled_asymp(nu = 0.5, x = x), from = 0.01, to = 10)
#' @export
log_besselI_scaled_asymp <- function(nu, x) {

  # Check length of nu and see if vectorized
  if (!(length(x) == 1 || length(nu) == 1 || length(nu) == length(x))) {

    stop("nu and must have equal lengths or one being a single number.")

  } else {

    single_nu <- length(nu) == 1

  }

  # Remove NAs in x and propagate to nu
  x_na <- is.na(x)
  x <- x[!x_na]
  if (anyNA(nu)) {

    stop("nu must not contain NAs.")

  }
  if (all(x_na)) {

    return(as.numeric(rep(NA, max(length(x_na), max(length(nu))))))

  }
  if (!single_nu) {

    nu <- nu[!x_na]

  }

  # Desired result
  res <- rep(NA, max(length(x), length(nu)))

  # Needed a x-asymptotic approximation?
  ind_x_asymp <- x >= 10000
  all_x_asymp <- all(ind_x_asymp)
  any_x_asymp <- all_x_asymp || any(ind_x_asymp)

  # Needed a nu-asymptotic approximation?
  ind_nu_asymp <- nu > 100
  all_nu_asymp <- all(ind_nu_asymp)
  any_nu_asymp <- all_nu_asymp || any(ind_nu_asymp)

  # Asymptotic approximation
  if (any_nu_asymp) {

    # Single nu?
    if (single_nu) {

      # https://dlmf.nist.gov/10.41#E3
      res <- Bessel::besselI.nuAsym(x = x, nu = nu, k.max = 5,
                                    expon.scaled = TRUE, log = TRUE)

    # Vectorized nu
    } else {

      if (length(x) == 1) {

        x <- rep(x, length.out = length(nu))

      }
      if (all_nu_asymp) {

        # https://dlmf.nist.gov/10.41#E3
        res <- Bessel::besselI.nuAsym(x = x, nu = nu, k.max = 5,
                                      expon.scaled = TRUE, log = TRUE)

      } else {

        # Regular Bessel + https://dlmf.nist.gov/10.41#E3
        res[!ind_nu_asymp] <- log(besselI(x = x[!ind_nu_asymp],
                                          nu = nu[!ind_nu_asymp],
                                          expon.scaled = TRUE))
        res[ind_nu_asymp] <- Bessel::besselI.nuAsym(x = x[ind_nu_asymp],
                                                    nu = nu[ind_nu_asymp],
                                                    k.max = 5,
                                                    expon.scaled = TRUE,
                                                    log = TRUE)

      }

    }

  # Standard algorithm
  } else {

    res <- log(besselI(x = x, nu = nu, expon.scaled = TRUE))

  }

  # Add NAs to x
  if (single_nu) {

    res_na <- rep(NA, length(x_na))
    res_na[!x_na] <- res

  } else {

    res_na <- rep(NA, length(nu))
    res_na[!x_na] <- res

  }
  return(res_na)

}