

#' @title Null distributions for circular and (hyper)spherical uniformity
#' statistics
#'
#' @description Approximate computation of the null distributions of several
#' statistics for assessing uniformity on the (hyper)sphere
#' \eqn{S^{p-1}:=\{{\bf x}\in R^p:||{\bf x}||=1\}}{S^{p-1}:=
#' \{x\in R^p:||x||=1\}}, \eqn{p\ge 2}. The approximation is done either by
#' means of the asymptotic distribution or by Monte Carlo.
#'
#' @param x evaluation points for the null distribution(s). Either a vector of
#' size \code{nx}, if the evaluation points are common for the tests in
#' \code{type}, or a matrix of size \code{c(nx, length(type))} with columns
#' containing the evaluation points for each test. Must not contain \code{NA}'s.
#' @inheritParams unif_test
#' @inheritParams r_unif
#' @param n sample size employed for computing the statistic.
#' @param approx type of approximation to the null distribution, either
#' \code{"asymp"} (default) for employing the asymptotic null distribution, if
#' available, or \code{"MC"}, for employing the Monte Carlo approximation of
#' the exact null distribution.
#' @param M number of Monte Carlo replications for approximating the null
#' distribution when \code{approx = "MC"}. Also, number of Monte Carlo samples
#' for approximating the asymptotic distributions based on weighted sums of chi
#' squared random variables. Defaults to \code{1e4}.
#' @param stats_MC a data frame of size \code{c(M, length(type))}, with column
#' names containing the character vector \code{type}, that results from
#' extracting \code{$stats_MC} from a call to \code{\link{unif_stat_MC}}. If
#' provided, the computation of Monte Carlo statistics when \code{approx = "MC"}
#' is skipped. \code{stats_MC} is checked internally to see if it is sorted.
#' Internally computed if \code{NULL} (default).
#' @inheritParams unif_stat
#' @inheritParams cir_stat
#' @inheritParams cir_stat_distr
#' @inheritParams sph_stat_distr
#' @param CJ12_beta \eqn{\beta} parameter in the exponential regime of CJ12
#' test, a positive real.
#' @param K_Kuiper,K_Watson,K_Watson_1976,K_Ajne integer giving the truncation
#' of the series present in the null asymptotic distributions. For the
#' Kolmogorov-Smirnov-related series defaults to \code{25}.
#' @param K_max integer giving the truncation of the series that compute the
#' asymptotic p-value of a Sobolev test. Defaults to \code{1e4}.
#' @param ... if \code{approx = "MC"}, optional performance parameters to be
#' passed to \cr\code{\link{unif_stat_MC}}: \code{chunks}, \code{cores},
#' and \code{seed}.
#' @return A data frame of size \code{c(nx, length(type))}, with column names
#' given by \code{type}, that contains the values of the null distributions of
#' the statistics evaluated at \code{x}.
#' @details
#' When \code{approx = "asymp"}, statistics that do not have an implemented or
#' known asymptotic are omitted, and a warning is generated.
#'
#' For Sobolev tests, \code{K_max = 1e4} produces probabilities uniformly
#' accurate with three digits for the \code{"PCvM"}, \code{"PAD"}, and
#' \code{"PRt"} tests, for dimensions \eqn{p \le 11}. With \code{K_max = 5e4},
#' these probabilities are uniformly accurate in the fourth digit. With
#' \code{K_max = 1e3}, only two-digit uniform accuracy is obtained. Uniform
#' accuracy deteriorates when \eqn{p} increases, e.g., a digit accuracy is lost
#' when \eqn{p = 51}.
#'
#' Descriptions and references on most of the asymptotic distributions
#' are available in García-Portugués and Verdebout (2018).
#' @references
#' García-Portugués, E. and Verdebout, T. (2018) An overview of uniformity
#' tests on the hypersphere. \emph{arXiv:1804.00286}.
#' \url{https://arxiv.org/abs/1804.00286}.
#' @examples
#' ## Asymptotic distribution
#'
#' # Circular statistics
#' x <- seq(0, 1, l = 5)
#' unif_stat_distr(x = x, type = "Kuiper", p = 2, n = 10)
#' unif_stat_distr(x = x, type = c("Ajne", "Kuiper"), p = 2, n = 10)
#' unif_stat_distr(x = x, type = c("Ajne", "Kuiper"), p = 2, n = 10, K_Ajne = 5)
#'\donttest{
#' # All circular statistics
#' unif_stat_distr(x = x, type = avail_cir_tests, p = 2, n = 10, K_max = 1e3)
#' }
#' # Spherical statistics
#' unif_stat_distr(x = cbind(x, x + 1), type = c("Rayleigh", "Bingham"),
#'                 p = 3, n = 10)
#' unif_stat_distr(x = cbind(x, x + 1), type = c("Rayleigh", "Bingham"),
#'                 p = 3, n = 10, M = 100)
#'\donttest{
#' # All spherical statistics
#' unif_stat_distr(x = x, type = avail_sph_tests, p = 3, n = 10, K_max = 1e3)
#'
#' ## Monte Carlo distribution
#'
#' # Circular statistics
#' x <- seq(0, 5, l = 10)
#' unif_stat_distr(x = x, type = avail_cir_tests, p = 2, n = 10, approx = "MC")
#' unif_stat_distr(x = x, type = "Kuiper", p = 2, n = 10, approx = "MC")
#' unif_stat_distr(x = x, type = c("Ajne", "Kuiper"), p = 2, n = 10,
#'                 approx = "MC")
#'
#' # Spherical statistics
#' unif_stat_distr(x = x, type = avail_sph_tests, p = 3, n = 10,
#'                 approx = "MC")
#' unif_stat_distr(x = cbind(x, x + 1), type = c("Rayleigh", "Bingham"),
#'                 p = 3, n = 10, approx = "MC")
#' unif_stat_distr(x = cbind(x, x + 1), type = c("Rayleigh", "Bingham"),
#'                 p = 3, n = 10, approx = "MC")
#'
#' ## Specific arguments
#'
#' # Rothman
#' unif_stat_distr(x = x, type = "Rothman", p = 2, n = 10, Rothman_t = 0.5,
#'                 approx = "MC")
#'
#' # CCF09
#' dirs <- r_unif_sph(n = 5, p = 3, M = 1)[, , 1]
#' x <- seq(0, 1, l = 10)
#' unif_stat_distr(x = x, type = "CCF09", p = 3, n = 10, approx = "MC",
#'                 CCF09_dirs = dirs)
#' unif_stat_distr(x = x, type = "CCF09", p = 3, n = 10, approx = "MC")
#'
#' # CJ12
#' unif_stat_distr(x = x, type = "CJ12", p = 3, n = 100, CJ12_reg = 3)
#' unif_stat_distr(x = x, type = "CJ12", p = 3, n = 100, CJ12_reg = 2,
#'                CJ12_beta = 0.01)
#' unif_stat_distr(x = x, type = "CJ12", p = 3, n = 100, CJ12_reg = 1)
#' }
#' @export
unif_stat_distr <- function(x, type, p, n, approx = "asymp", M = 1e4,
                            stats_MC = NULL, Rothman_t = 1 / 3, Pycke_q = 0.5,
                            Riesz_s = 1, CCF09_dirs = NULL, CJ12_reg = 3,
                            CJ12_beta = 0, Stephens = FALSE, K_Kuiper = 25,
                            K_Watson = 25, K_Watson_1976 = 5, K_Ajne = 5e2,
                            K_CCF09 = 25, K_max = 1e4, ...) {

  # Stop if NA's
  if (anyNA(x)) {

    stop("NAs present in x, please remove them.")

  }

  # As a matrix
  if (is.vector(x)) {

    x <- matrix(x, ncol = 1)

  }

  # Dimension
  if (p == 2) {

    avail_stats <- avail_cir_tests
    prefix_distr <- "p_cir_stat_"

  } else {

    avail_stats <- avail_sph_tests
    prefix_distr <- "p_sph_stat_"

  }

  # Get the type of statistics
  if (is.character(type)) {

    type <- gsub("-", "_", type)
    type <- unique(type)
    if ("all" %in% type) {

      stats_type <- avail_stats

    } else {

      stats_type <- match.arg(arg = type, choices = avail_stats,
                              several.ok = TRUE)

    }

  } else if (is.numeric(type)) {

    type <- unique(type)
    stats_type <- avail_stats[type]

  } else {

    stop("type must be a character or a numeric vector.")

  }

  # Omit statistics that do not have asymptotic distribution
  if (approx == "asymp") {

    ind_asymp <- sapply(paste0(prefix_distr, stats_type),
                        exists, where = asNamespace("sphunif"))
    if (!all(ind_asymp)) {

      warning(paste0(paste("Omitting the following statistics with not",
                           "implemented or known asymptotic distributions: "),
                     paste(stats_type[!ind_asymp], collapse = ", "), "."))
      stats_type <- stats_type[ind_asymp]

    }

  }

  # Number of statistics
  n_stats <- length(stats_type)

  # Check dimension
  n_col_x <- ncol(x)
  if (n_col_x != n_stats) {

    if (n_col_x == 1) {

      x <- matrix(x, nrow = nrow(x), ncol = n_stats, byrow = FALSE)

    } else {

      stop(paste("Check x and type, ncol(x) must equal the number of valid",
                 "type values."))
    }

  }

  # Put x as a data frame
  x <- as.data.frame(x)
  names(x) <- stats_type

  # Asymptotic or Monte Carlo distribution?
  if (approx == "asymp") {

    # Optional arguments
    args <- list("t" = Rothman_t, "q" = Pycke_q, "s" = Riesz_s,
                 "dirs" = CCF09_dirs, "regime" = CJ12_reg, "beta" = CJ12_beta,
                 "Stephens" = Stephens, "K_Kuiper" = K_Kuiper,
                 "K_Watson" = K_Watson, "K_Watson_1976" = K_Watson_1976,
                 "K_Ajne" = K_Ajne, "K_CCF09" = K_CCF09, "K_max" = K_max,
                 "thre" = 0, "n" = n, "p" = p)
    names_args <- names(args)

    # Evaluate distributions
    distrs <- sapply(stats_type, function(distr) {

      # Additional arguments
      name_distr <- paste0(prefix_distr, distr)
      distr_args <- args[names_args %in% names(formals(name_distr))]

      # Call distribution
      do.call(what = name_distr, args = c(list(x = x[[distr]]), distr_args))

    }, simplify = FALSE)

  } else if (approx == "MC") {

    # Sample statistics under the null
    if (is.null(stats_MC)) {

      stats_MC <- unif_stat_MC(n = n, type = stats_type, p = p, M = M,
                               r_H1 = NULL, crit_val = NULL,
                               alpha = c(0.10, 0.05, 0.01),
                               return_stats = TRUE, stats_sorted = TRUE,
                               Rothman_t = Rothman_t, Pycke_q = Pycke_q,
                               CCF09_dirs = CCF09_dirs, CJ12_reg = CJ12_reg,
                               K_CCF09 = K_CCF09, ...)$stats_MC

    } else {

      # Check if there is any missing statistic in stats_MC
      names_stats_MC <- names(stats_MC)
      missing_stats_MC <- !(stats_type %in% names_stats_MC)
      if (any(missing_stats_MC)) {

        stop(paste0("Missing statistics in stats_MC: \"",
                    paste(stats_type[missing_stats_MC],
                          collapse = "\", \""), "\"."))

      }

      # Sort for faster comparisons
      if (any(apply(stats_MC, 2, is.unsorted))) {

        stats_MC <- as.data.frame(sort_each_col(as.matrix(stats_MC)))
        names(stats_MC) <- names_stats_MC

      }

      # Match columns of stats_MC with stats_type
      stats_MC <- stats_MC[, pmatch(x = stats_type, table = names_stats_MC),
                           drop = FALSE]

    }

    # Is it required to sort x?
    if (any(apply(x, 2, is.unsorted))) {

      # Indexes for sorting and then unsorting
      ind_1 <- sort_index_each_col(x)
      ind_2 <- sort_index_each_col(ind_1)

      # Approximate distributions
      distrs <- sapply(stats_type, function(distr) {
        ecdf_bin(data = stats_MC[[distr]], sorted_x = x[[distr]][ind_1],
                 data_sorted = TRUE, efic = TRUE, divide_n = TRUE)
      }, simplify = FALSE)
      distrs <- distrs[ind_2, ]

    } else {

      # Approximate distributions
      distrs <- sapply(stats_type, function(distr) {
        ecdf_bin(data = stats_MC[[distr]], sorted_x = x[[distr]],
                 data_sorted = TRUE, efic = TRUE, divide_n = TRUE)
      }, simplify = FALSE)

    }

  } else {

    stop("Wrong choice for approx.")

  }

  # Return data frame
  distrs <- as.data.frame(lapply(lapply(distrs, pmin, 1), pmax, 0))
  names(distrs) <- stats_type
  return(distrs)

}
