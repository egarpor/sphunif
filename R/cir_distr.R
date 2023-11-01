

#' @title Asymptotic distributions for circular uniformity statistics
#'
#' @description Computation of the asymptotic null distributions of circular
#' uniformity statistics.
#'
#' @param x a vector of size \code{nx} or a matrix of size \code{c(nx, 1)}.
#' @inheritParams Sobolev
#' @param thre error threshold for the tail probability given by the
#' the first terms of the truncated series of a Sobolev test. Defaults to
#' \code{0} (no further truncation).
#' @param K_Kolmogorov,K_Kuiper,K_Watson,K_Watson_1976,K_Ajne integer giving
#' the truncation of the series present in the null asymptotic distributions.
#' For the Kolmogorov-Smirnov-related series defaults to \code{25}; for the
#' others series defaults to a smaller number.
#' @inheritParams unif_stat_distr
#' @inheritParams cir_stat
#' @param alternating use the alternating series expansion for the distribution
#' of the Kolmogorov-Smirnov statistic? Defaults to \code{TRUE}.
#' @param second_term use the second-order series expansion for the
#' distribution of the Kuiper statistic? Defaults to \code{TRUE}.
#' @param N number of points used in the
#' \link[=Gauss_Legen_nodes]{Gauss-Legendre quadrature}. Defaults to \code{40}.
#' @param exact use the exact distribution for the Hodges-Ajne statistic?
#' Defaults to \code{TRUE}.
#' @param asymp_std compute the distribution associated to the normalized
#' Hodges-Ajne statistic? Defaults to \code{FALSE}.
#' @param max_gap compute the distribution associated to the maximum gap for
#' the range statistic? Defaults to \code{TRUE}.
#' @param abs_val compute the distribution associated to the absolute value of
#' the Darling's log gaps statistic? Defaults to \code{TRUE}.
#' @inheritParams sph_stat_distr
#' @param ... further parameters passed to \code{\link{p_Sobolev}} or
#' \code{\link{d_Sobolev}} (such as \code{x_tail}).
#' @return A matrix of size \code{c(nx, 1)} with the evaluation of the
#' distribution or density function at \code{x}.
#' @references García-Portugués, E. and Verdebout, T. (2018) An overview of
#' uniformity tests on the hypersphere. \emph{arXiv:1804.00286}.
#' \url{https://arxiv.org/abs/1804.00286}.
#' @details
#' Descriptions and references for most of the tests are available
#' in García-Portugués and Verdebout (2018).
#' @examples
#' # Ajne
#' curve(d_cir_stat_Ajne(x), to = 1.5, n = 2e2, ylim = c(0, 4))
#' curve(p_cir_stat_Ajne(x), n = 2e2, col = 2, add = TRUE)
#'
#' # Bakshaev
#' curve(d_cir_stat_Bakshaev(x, method = "HBE"), to = 6, n = 2e2,
#'       ylim = c(0, 1))
#' curve(p_cir_stat_Bakshaev(x, method = "HBE"), n = 2e2, add = TRUE, col = 2)
#'
#' # Bingham
#' curve(d_cir_stat_Bingham(x), to = 12, n = 2e2, ylim = c(0, 1))
#' curve(p_cir_stat_Bingham(x), n = 2e2, col = 2, add = TRUE)
#
#' # Greenwood
#' curve(d_cir_stat_Greenwood(x), from = -6, to = 6, n = 2e2, ylim = c(0, 1))
#' curve(p_cir_stat_Greenwood(x), n = 2e2, col = 2, add = TRUE)
#'
#' # Hermans-Rasson
#' curve(p_cir_stat_Hermans_Rasson(x, method = "HBE"), to = 10, n = 2e2,
#'       ylim = c(0, 1))
#' curve(d_cir_stat_Hermans_Rasson(x, method = "HBE"), n = 2e2, add = TRUE,
#'       col = 2)
#'
#' # Hodges-Ajne
#' plot(25:45, d_cir_stat_Hodges_Ajne(cbind(25:45), n = 50), type = "h",
#'      lwd = 2, ylim = c(0, 1))
#' lines(25:45, p_cir_stat_Hodges_Ajne(cbind(25:45), n = 50), type = "s",
#'       col = 2)
#'
#' # Kolmogorov-Smirnov
#' curve(d_Kolmogorov(x), to = 3, n = 2e2, ylim = c(0, 2))
#' curve(p_Kolmogorov(x), n = 2e2, col = 2, add = TRUE)
#'
#' # Kuiper
#' curve(d_cir_stat_Kuiper(x, n = 50), to = 3, n = 2e2, ylim = c(0, 2))
#' curve(p_cir_stat_Kuiper(x, n = 50), n = 2e2, col = 2, add = TRUE)
#'
#' # Kuiper and Watson with Stephens modification
#' curve(d_cir_stat_Kuiper(x, n = 8, Stephens = TRUE), to = 2.5, n = 2e2,
#'       ylim = c(0, 10))
#' curve(d_cir_stat_Watson(x, n = 8, Stephens = TRUE), n = 2e2, lty = 2,
#'       add = TRUE)
#' n <- c(10, 20, 30, 40, 50, 100, 500)
#' col <- rainbow(length(n))
#' for (i in seq_along(n)) {
#'   curve(d_cir_stat_Kuiper(x, n = n[i], Stephens = TRUE), n = 2e2,
#'         col = col[i], add = TRUE)
#'   curve(d_cir_stat_Watson(x, n = n[i], Stephens = TRUE), n = 2e2,
#'         col = col[i], lty = 2, add = TRUE)
#' }
#'
#' # Maximum uncovered spacing
#' curve(d_cir_stat_Max_uncover(x), from = -3, to = 6, n = 2e2, ylim = c(0, 1))
#' curve(p_cir_stat_Max_uncover(x), n = 2e2, col = 2, add = TRUE)
#'
#' # Number of uncovered spacing
#' curve(d_cir_stat_Num_uncover(x), from = -4, to = 4, n = 2e2, ylim = c(0, 1))
#' curve(p_cir_stat_Num_uncover(x), n = 2e2, col = 2, add = TRUE)
#'
#' # Log gaps
#' curve(d_cir_stat_Log_gaps(x), from = -1, to = 4, n = 2e2, ylim = c(0, 1))
#' curve(p_cir_stat_Log_gaps(x), n = 2e2, col = 2, add = TRUE)
#'
#' # Gine Fn
#' curve(d_cir_stat_Gine_Fn(x, method = "HBE"), to = 2.5, n = 2e2,
#'       ylim = c(0, 2))
#' curve(p_cir_stat_Gine_Fn(x, method = "HBE"), n = 2e2, add = TRUE, col = 2)
#'
#' # Gine Gn
#' curve(d_cir_stat_Gine_Gn(x, method = "HBE"), to = 2.5, n = 2e2,
#'       ylim = c(0, 2))
#' curve(p_cir_stat_Gine_Gn(x, method = "HBE"), n = 2e2, add = TRUE, col = 2)
#'
#' # Gini mean difference
#' curve(d_cir_stat_Gini(x), from = -4, to = 4, n = 2e2, ylim = c(0, 1))
#' curve(p_cir_stat_Gini(x), n = 2e2, col = 2, add = TRUE)
#'
#' # Gini mean squared difference
#' curve(d_cir_stat_Gini_squared(x), from = -10, to = 10, n = 2e2,
#'       ylim = c(0, 1))
#' curve(p_cir_stat_Gini_squared(x), n = 2e2, col = 2, add = TRUE)
#'
#' # PAD
#' curve(d_cir_stat_PAD(x, method = "HBE"), to = 3, n = 2e2, ylim = c(0, 1.5))
#' curve(p_cir_stat_PAD(x, method = "HBE"), n = 2e2, add = TRUE, col = 2)
#'
#' # PCvM
#' curve(d_cir_stat_PCvM(x, method = "HBE"), to = 4, n = 2e2, ylim = c(0, 2))
#' curve(p_cir_stat_PCvM(x, method = "HBE"), n = 2e2, add = TRUE, col = 2)
#'
#' # PRt
#' curve(d_cir_stat_PRt(x, method = "HBE"), n = 2e2, ylim = c(0, 5))
#' curve(p_cir_stat_PRt(x, method = "HBE"), n = 2e2, add = TRUE, col = 2)
#'
#' # Pycke
#' curve(d_cir_stat_Pycke(x), from = -5, to = 10, n = 2e2, ylim = c(0, 1))
#' curve(p_cir_stat_Pycke(x), n = 2e2, col = 2, add = TRUE)
#'
#' # Pycke q
#' curve(d_cir_stat_Pycke_q(x, method = "HBE"), to = 15, n = 2e2,
#'       ylim = c(0, 1))
#' curve(p_cir_stat_Pycke_q(x, method = "HBE"), n = 2e2, add = TRUE, col = 2)
#'
#' # Range
#' curve(d_cir_stat_Range(x, n = 50), to = 2, n = 2e2, ylim = c(0, 4))
#' curve(p_cir_stat_Range(x, n = 50), n = 2e2, col = 2, add = TRUE)
#'
#' # Rao
#' curve(d_cir_stat_Rao(x), from = -6, to = 6, n = 2e2, ylim = c(0, 1))
#' curve(p_cir_stat_Rao(x), n = 2e2, col = 2, add = TRUE)
#'
#' # Rayleigh
#' curve(d_cir_stat_Rayleigh(x), to = 12, n = 2e2, ylim = c(0, 1))
#' curve(p_cir_stat_Rayleigh(x), n = 2e2, col = 2, add = TRUE)
#'
#' # Riesz
#' curve(d_cir_stat_Riesz(x, method = "HBE"), to = 6, n = 2e2,
#'       ylim = c(0, 1))
#' curve(p_cir_stat_Riesz(x, method = "HBE"), n = 2e2, add = TRUE, col = 2)
#'
#' # Rothman
#' curve(d_cir_stat_Rothman(x, method = "HBE"), n = 2e2, ylim = c(0, 5))
#' curve(p_cir_stat_Rothman(x, method = "HBE"), n = 2e2, add = TRUE, col = 2)
#'
#' # Vacancy
#' curve(d_cir_stat_Vacancy(x), from = -4, to = 4, n = 2e2, ylim = c(0, 1))
#' curve(p_cir_stat_Vacancy(x), n = 2e2, col = 2, add = TRUE)
#'
#' # Watson
#' curve(d_cir_stat_Watson(x), to = 0.5, n = 2e2, ylim = c(0, 15))
#' curve(p_cir_stat_Watson(x), n = 2e2, col = 2, add = TRUE)
#'
#' # Watson (1976)
#' curve(d_cir_stat_Watson_1976(x), to = 1.5, n = 2e2, ylim = c(0, 3))
#' curve(p_cir_stat_Watson_1976(x), n = 2e2, col = 2, add = TRUE)
#'
#' # Sobolev
#' vk2 <- c(0.5, 0)
#' curve(d_cir_stat_Sobolev(x = x, vk2 = vk2), to = 3, n = 2e2, ylim = c(0, 2))
#' curve(p_cir_stat_Sobolev(x = x, vk2 = vk2), n = 2e2, col = 2, add = TRUE)
#' @name cir_stat_distr
NULL


#' @rdname cir_stat_distr
#' @export
p_cir_stat_Bakshaev <- function(x, K_max = 1e3, thre = 0, ...) {

  p_sph_stat_Bakshaev(x = x, p = 2, K_max = K_max, thre = thre, ...)

}


#' @rdname cir_stat_distr
#' @export
d_cir_stat_Bakshaev <- function(x, K_max = 1e3, thre = 0, ...) {

  d_sph_stat_Bakshaev(x = x, p = 2, K_max = K_max, thre = thre, ...)

}


#' @rdname cir_stat_distr
#' @export
p_cir_stat_Gine_Fn <- function(x, K_max = 1e3, thre = 0, ...) {

  p_sph_stat_Gine_Fn(x = x, p = 2, K_max = K_max, thre = thre, ...)

}


#' @rdname cir_stat_distr
#' @export
d_cir_stat_Gine_Fn <- function(x, K_max = 1e3, thre = 0, ...) {

  d_sph_stat_Gine_Fn(x = x, p = 2, K_max = K_max, thre = thre, ...)

}


#' @rdname cir_stat_distr
#' @export
p_cir_stat_Gine_Gn <- function(x, K_max = 1e3, thre = 0, ...) {

  p_sph_stat_Gine_Gn(x = x, p = 2, K_max = K_max, thre = thre, ...)

}


#' @rdname cir_stat_distr
#' @export
d_cir_stat_Gine_Gn <- function(x, K_max = 1e3, thre = 0, ...) {

  d_sph_stat_Gine_Gn(x = x, p = 2, K_max = K_max, thre = thre, ...)

}



#' @rdname cir_stat_distr
#' @export
p_cir_stat_Hermans_Rasson <- function(x, K_max = 1e3, thre = 0, ...) {

  cbind(p_Sobolev(x = x, p = 2, type = "Hermans_Rasson", K_max = K_max,
                  thre = thre, ...))

}


#' @rdname cir_stat_distr
#' @export
d_cir_stat_Hermans_Rasson <- function(x, K_max = 1e3, thre = 0, ...) {

  cbind(d_Sobolev(x = x, p = 2, type = "Hermans_Rasson", K_max = K_max,
                  thre = thre, ...))

}


#' @rdname cir_stat_distr
#' @export
p_cir_stat_PAD <- function(x, K_max = 1e3, thre = 0, ...) {

  p_sph_stat_PAD(x = x, p = 2, K_max = K_max, thre = thre, ...)

}


#' @rdname cir_stat_distr
#' @export
d_cir_stat_PAD <- function(x, K_max = 1e3, thre = 0, ...) {

  d_sph_stat_PAD(x = x, p = 2, K_max = K_max, thre = thre, ...)

}


#' @rdname cir_stat_distr
#' @export
p_cir_stat_PCvM <- function(x, K_max = 1e3, thre = 0, ...) {

  p_sph_stat_PCvM(x = x, p = 2, K_max = K_max, thre = thre, ...)

}


#' @rdname cir_stat_distr
#' @export
d_cir_stat_PCvM <- function(x, K_max = 1e3, thre = 0, ...) {

  d_sph_stat_PCvM(x = x, p = 2, K_max = K_max, thre = thre, ...)

}


#' @rdname cir_stat_distr
#' @export
p_cir_stat_PRt <- function(x, t = 1 / 3, K_max = 1e3, thre = 0, ...) {

  p_sph_stat_PRt(x = x, p = 2, t = t, K_max = K_max, thre = thre, ...)

}


#' @rdname cir_stat_distr
#' @export
d_cir_stat_PRt <- function(x, t = 1 / 3, K_max = 1e3, thre = 0, ...) {

  d_sph_stat_PRt(x = x, p = 2, t = t, K_max = K_max, thre = thre, ...)

}


#' @rdname cir_stat_distr
#' @export
p_cir_stat_Pycke_q <- function(x, q = 0.5, K_max = 1e3, thre = 0, ...) {

  cbind(p_Sobolev(x = x, p = 2, type = "Pycke_q", Pycke_q = q,
                  K_max = K_max, thre = thre, ...))

}


#' @rdname cir_stat_distr
#' @export
d_cir_stat_Pycke_q <- function(x, q = 0.5, K_max = 1e3, thre = 0, ...) {

  cbind(d_Sobolev(x = x, p = 2, type = "Pycke_q", Pycke_q = q,
                  K_max = K_max, thre = thre, ...))

}


#' @rdname cir_stat_distr
#' @export
p_cir_stat_Rothman <- function(x, t = 1 / 3, K_max = 1e3, thre = 0, ...) {

  cbind(p_Sobolev(x = x, p = 2, type = "Rothman", Rothman_t = t,
                  K_max = K_max, thre = thre, ...))

}


#' @rdname cir_stat_distr
#' @export
d_cir_stat_Rothman <- function(x, t = 1 / 3, K_max = 1e3, thre = 0, ...) {

  cbind(d_Sobolev(x = x, p = 2, type = "Rothman", Rothman_t = t,
                  K_max = K_max, thre = thre, ...))

}


#' @rdname cir_stat_distr
#' @export
p_cir_stat_Riesz <- function(x, s = 1, K_max = 1e3, thre = 0, ...) {

  p_sph_stat_Riesz(x = x, p = 2, s = s, K_max = K_max, thre = thre, ...)

}


#' @rdname cir_stat_distr
#' @export
d_cir_stat_Riesz <- function(x, s = 1, K_max = 1e3, thre = 0, ...) {

  d_sph_stat_Riesz(x = x, p = 2, s = s, K_max = K_max, thre = thre, ...)

}


#' @rdname cir_stat_distr
#' @export
p_cir_stat_Sobolev <- function(x, vk2 = c(0, 0, 1), ...) {

  p_sph_stat_Sobolev(x = x, p = 2, vk2 = vk2, ...)

}


#' @rdname cir_stat_distr
#' @export
d_cir_stat_Sobolev <- function(x, vk2 = c(0, 0, 1), ...) {

  d_sph_stat_Sobolev(x = x, p = 2, vk2 = vk2, ...)

}