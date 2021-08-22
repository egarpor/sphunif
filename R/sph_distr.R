

#' @title Asymptotic distributions for spherical uniformity statistics
#'
#' @description Computation of the asymptotic null distributions of
#' spherical uniformity statistics.
#'
#' @inheritParams cir_stat_distr
#' @inheritParams Sobolev
#' @inheritParams r_unif
#' @inheritParams wschisq
#' @inheritParams unif_stat_distr
#' @param regime type of asymptotic regime for the CJ12 test, either \code{1}
#' (sub-exponential regime), \code{2} (exponential), or \code{3}
#' (super-exponential; default).
#' @param beta \eqn{\beta} parameter in the exponential regime of the CJ12
#' test, a nonnegative real. Defaults to \code{0}.
#' @inheritParams cir_stat
#' @param ... further parameters passed to \code{\link{p_Sobolev}} or
#' \code{\link{d_Sobolev}} (such as \code{x_tail}).
#' @return
#' \itemize{
#'   \item \code{r_sph_stat_*}: a matrix of size \code{c(n, 1)} containing
#'   the sample.
#'   \item \code{p_sph_stat_*}, \code{d_sph_stat_*}: a matrix of size
#'   \code{c(nx, 1)} with the evaluation of the distribution or density
#'   functions at \code{x}.
#' }
#' @details
#' Descriptions and references on most of the asymptotic distributions
#' are available in García-Portugués and Verdebout (2018).
#' @examples
#' # Ajne
#' curve(d_sph_stat_Ajne(x, p = 3, method = "HBE"), n = 2e2, ylim = c(0, 4))
#' curve(p_sph_stat_Ajne(x, p = 3, method = "HBE"), n = 2e2, col = 2,
#'       add = TRUE)
#'
#' # Bakshaev
#' curve(d_sph_stat_Bakshaev(x, p = 3, method = "HBE"), to = 5, n = 2e2,
#'       ylim = c(0, 2))
#' curve(p_sph_stat_Bakshaev(x, p = 3, method = "HBE"), n = 2e2, col = 2,
#'       add = TRUE)
#'
#' # Bingham
#' curve(d_sph_stat_Bingham(x, p = 3), to = 20, n = 2e2, ylim = c(0, 1))
#' curve(p_sph_stat_Bingham(x, p = 3), n = 2e2, col = 2, add = TRUE)
#'
#' # CJ12
#' curve(d_sph_stat_CJ12(x, regime = 1), from = -10, to = 10, n = 2e2,
#'       ylim = c(0, 1))
#' curve(d_sph_stat_CJ12(x, regime = 2, beta = 0.1), n = 2e2, col = 2,
#'       add = TRUE)
#' curve(d_sph_stat_CJ12(x, regime = 3), n = 2e2, col = 3, add = TRUE)
#' curve(p_sph_stat_CJ12(x, regime = 1), n = 2e2, col = 1, add = TRUE)
#' curve(p_sph_stat_CJ12(x, regime = 2, beta = 0.1), n = 2e2, col = 2,
#'       add = TRUE)
#' curve(p_sph_stat_CJ12(x, regime = 3), col = 3, add = TRUE)
#'
#' # Gine Fn
#' curve(d_sph_stat_Gine_Fn(x, p = 3, method = "HBE"), to = 2, n = 2e2,
#'       ylim = c(0, 2))
#' curve(p_sph_stat_Gine_Fn(x, p = 3, method = "HBE"), n = 2e2, col = 2,
#'       add = TRUE)
#'
#' # Gine Gn
#' curve(d_sph_stat_Gine_Gn(x, p = 3, method = "HBE"), to = 1.5, n = 2e2,
#'       ylim = c(0, 2.5))
#' curve(p_sph_stat_Gine_Gn(x, p = 3, method = "HBE"), n = 2e2, col = 2,
#'       add = TRUE)
#'
#' # PAD
#' curve(d_sph_stat_PAD(x, p = 3, method = "HBE"), to = 3, n = 2e2,
#'       ylim = c(0, 1.5))
#' curve(p_sph_stat_PAD(x, p = 3, method = "HBE"), n = 2e2, col = 2,
#'       add = TRUE)
#'
#' # PCvM
#' curve(d_sph_stat_PCvM(x, p = 3, method = "HBE"), to = 0.6, n = 2e2,
#'       ylim = c(0, 7))
#' curve(p_sph_stat_PCvM(x, p = 3, method = "HBE"), n = 2e2, col = 2,
#'       add = TRUE)
#'
#' # PRt
#' curve(d_sph_stat_PRt(x, p = 3, method = "HBE"), n = 2e2, ylim = c(0, 5))
#' curve(p_sph_stat_PRt(x, p = 3, method = "HBE"), n = 2e2, col = 2, add = TRUE)
#'
#' # Rayleigh
#' curve(d_sph_stat_Rayleigh(x, p = 3), to = 15, n = 2e2, ylim = c(0, 1))
#' curve(p_sph_stat_Rayleigh(x, p = 3), n = 2e2, col = 2, add = TRUE)
#'
#' # HD-standardized Rayleigh
#' curve(d_sph_stat_Rayleigh_HD(x, p = 3), from = -4, to = 4, n = 2e2,
#'       ylim = c(0, 1))
#' curve(p_sph_stat_Rayleigh_HD(x, p = 3), n = 2e2, col = 2, add = TRUE)
#'
#' # Riesz
#' curve(d_sph_stat_Riesz(x, p = 3, method = "HBE"), n = 2e2, from = 0, to = 5,
#'       ylim = c(0, 2))
#' curve(p_sph_stat_Riesz(x, p = 3, method = "HBE"), n = 2e2, col = 2,
#'       add = TRUE)
#' @name sph_stat_distr
NULL


#' @rdname sph_stat_distr
#' @export
p_sph_stat_Ajne <- function(x, p, K_max = 1e3, thre = 0, ...) {

  cbind(p_Sobolev(x = x, p = p, type = "Ajne", K_max = K_max, thre = thre, ...))

}


#' @rdname sph_stat_distr
#' @export
d_sph_stat_Ajne <- function(x, p, K_max = 1e3, thre = 0, ...) {

  cbind(d_Sobolev(x = x, p = p, type = "Ajne", K_max = K_max, thre = thre, ...))

}


#' @rdname sph_stat_distr
#' @export
p_sph_stat_Bakshaev <- function(x, p, K_max = 1e3, thre = 0, ...) {

  cbind(p_Sobolev(x = x, p = p, type = "Bakshaev", K_max = K_max, thre = thre,
                  ...))

}


#' @rdname sph_stat_distr
#' @export
d_sph_stat_Bakshaev <- function(x, p, K_max = 1e3, thre = 0, ...) {

  cbind(d_Sobolev(x = x, p = p, type = "Bakshaev", K_max = K_max, thre = thre,
                  ...))

}


#' @rdname sph_stat_distr
#' @export
p_sph_stat_Gine_Fn <- function(x, p, K_max = 1e3, thre = 0, ...) {

  cbind(p_Sobolev(x = x, p = p, type = "Gine_Fn", K_max = K_max, thre = thre,
                  ...))

}


#' @rdname sph_stat_distr
#' @export
d_sph_stat_Gine_Fn <- function(x, p, K_max = 1e3, thre = 0, ...) {

  cbind(d_Sobolev(x = x, p = p, type = "Gine_Fn", K_max = K_max, thre = thre,
                  ...))

}


#' @rdname sph_stat_distr
#' @export
p_sph_stat_Gine_Gn <- function(x, p, K_max = 1e3, thre = 0, ...) {

  cbind(p_Sobolev(x = x, p = p, type = "Gine_Gn", K_max = K_max, thre = thre,
                  ...))

}


#' @rdname sph_stat_distr
#' @export
d_sph_stat_Gine_Gn <- function(x, p, K_max = 1e3, thre = 0, ...) {

  cbind(d_Sobolev(x = x, p = p, type = "Gine_Gn", K_max = K_max, thre = thre,
                  ...))

}


#' @rdname sph_stat_distr
#' @export
p_sph_stat_PAD <- function(x, p, K_max = 1e3, thre = 0, ...) {

  cbind(p_Sobolev(x = x, p = p, type = "PAD", K_max = K_max, thre = thre,
                  ...))

}


#' @rdname sph_stat_distr
#' @export
d_sph_stat_PAD <- function(x, p, K_max = 1e3, thre = 0, ...) {

  cbind(d_Sobolev(x = x, p = p, type = "PAD", K_max = K_max, thre = thre,
                  ...))

}


#' @rdname sph_stat_distr
#' @export
p_sph_stat_PCvM <- function(x, p, K_max = 1e3, thre = 0, ...) {

  cbind(p_Sobolev(x = x, p = p, type = "PCvM", K_max = K_max, thre = thre,
                  ...))

}


#' @rdname sph_stat_distr
#' @export
d_sph_stat_PCvM <- function(x, p, K_max = 1e3, thre = 0, ...) {

  cbind(d_Sobolev(x = x, p = p, type = "PCvM", K_max = K_max, thre = thre,
                  ...))

}


#' @rdname sph_stat_distr
#' @export
p_sph_stat_PRt <- function(x, p, t = 1 / 3, K_max = 1e3, thre = 0, ...) {

  cbind(p_Sobolev(x = x, p = p, type = "PRt", Rothman_t = t,
                  K_max = K_max, thre = thre, ...))

}


#' @rdname sph_stat_distr
#' @export
d_sph_stat_PRt <- function(x, p, t = 1 / 3, K_max = 1e3, thre = 0, ...) {

  cbind(d_Sobolev(x = x, p = p, type = "PRt", Rothman_t = t,
                  K_max = K_max, thre = thre, ...))

}


#' @rdname sph_stat_distr
#' @export
p_sph_stat_Riesz <- function(x, p, s = 1, K_max = 1e3, thre = 0, ...) {
  
  cbind(p_Sobolev(x = x, p = p, type = "Riesz", Riesz_s = s, K_max = K_max,
                  thre = thre, ...))
  
}


#' @rdname sph_stat_distr
#' @export
d_sph_stat_Riesz <- function(x, p, s = 1, K_max = 1e3, thre = 0, ...) {

  cbind(d_Sobolev(x = x, p = p, type = "Riesz", Riesz_s = s, K_max = K_max,
                  thre = thre, ...))

}

