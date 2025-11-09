

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
#' @inheritParams Sobolev method
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
#' \doi{10.48550/arXiv.1804.00286}.
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
#'
#' ## Sobolev
#'
#' x <- seq(0, 1, l = 10)
#' vk2 <- diag(1, nrow = 3)
#' unif_stat_distr(x = x, type = "Sobolev", approx = "asymp", p = 3, n = 100,
#'                 Sobolev_vk2 = vk2)
#' sapply(1:3, function(i)
#'   unif_stat_distr(x = x, type = "Sobolev", approx = "asymp", p = 3, n = 100,
#'                   Sobolev_vk2 = vk2[i, ])$Sobolev)
#' sapply(1:3, function(i)
#'   unif_stat_distr(x = x, type = "Sobolev", approx = "MC", p = 3, n = 100,
#'                   Sobolev_vk2 = vk2[i, ], M = 1e3)$Sobolev)
#' unif_stat_distr(x = x, type = "Sobolev", approx = "MC", p = 3, n = 100,
#'                 Sobolev_vk2 = vk2, M = 1e3)
#' }
#' @export
unif_stat_distr <- function(x, type, p, n, approx = "asymp", M = 1e4,
                            stats_MC = NULL, K_max = 1e4, method = "I",
                            Stephens = FALSE, CCF09_dirs = NULL, CJ12_beta = 0,
                            CJ12_reg = 3, cov_a = 2 * pi, Cressie_t = 1 / 3,
                            K_Ajne = 5e2, K_CCF09 = 25, K_Kuiper = 25,
                            K_Watson = 25, K_Watson_1976 = 5, Poisson_rho = 0.5,
                            Pycke_q = 0.5, Rayleigh_m = 1, Riesz_s = 1,
                            Rothman_t = 1 / 3, Sobolev_vk2 = c(0, 0, 1),
                            Softmax_kappa = 1, Stein_K = 10, Stein_cf = FALSE,
                            Stereo_a = 0, ...) {


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

      stats_type <- try(match.arg(arg = type, choices = avail_stats,
                                  several.ok = TRUE), silent = TRUE)
      if (inherits(stats_type, "try-error")) {

        stop(paste(strwrap(paste0(
          "Statistic with type = \"", type, "\" is unsupported. ",
          "Must be one of the following tests: \"",
          paste(avail_stats, collapse = "\", \""), "\"."),
          width = 80, indent = 0, exdent = 2), collapse = "\n"))

      }

    }

  } else if (is.numeric(type)) {

    type <- unique(type)
    if (type > length(avail_stats)) {

      stop("type must be a numeric vector with values between 1 and ",
           length(avail_stats), ".")

    } else {

        stats_type <- avail_stats[type]

    }

  } else {

    stop("type must be a character or a numeric vector.")

  }

  # Check if n is missing
  if (missing(n)) {

    if (approx == "MC") {

      stop("n must be specified for approx = \"MC\".")

    } else if (approx == "asymp") {

      if ("Watson" %in% type) {

        warning(paste("n is not specified; set to n = 0. The asymptotic",
                      "distribution of the Watson statistic accepts n to",
                      "improve its accuracy with Stephens = TRUE."))

      }
      if ("Kuiper" %in% type) {

        stop(paste("n is not specified. The asymptotic distribution of the",
                   "Kuiper statistic requires n to improve its accuracy",
                   "(and Stephens = TRUE can also be set)."))

      }
      if (any(c("Hodges_Ajne", "Range") %in% type)) {

        stop(paste("n is not specified. The distributions of the Hodges-Ajne",
                   "and Range statistics require to specify n."))

      }
      n <- 0

    }

  }

  # Omit statistics that do not have asymptotic distribution
  if (approx == "asymp") {

    # Exclude vectorizations. TODO: Allow for vectorizations
    param_vectorized <- c("cov_a", "Cressie_t", "Poisson_rho", "Pycke_q",
                          "Rayleigh_m", "Riesz_s", "Rothman_t",
                          "Softmax_kappa", "Stein_K", "Stein_cf",
                          "Stereo_a") # "Sobolev_vk2"
    n_param_vectorized <- sapply(param_vectorized, function(par) {
      obj <- get(x = par)
      ifelse(is.null(dim(obj)), length(obj), nrow(obj))
    })
    if (any(n_param_vectorized > 1)) {

      stop(paste("No vectorization is implemented for",
                 paste(param_vectorized[n_param_vectorized > 1],
                       collapse = ", "),
                 "with approx = \"asymp\". Either use approx =",
                 "\"MC\" or unvectorize"))

    }

    # Search for p_*_stat_*
    ind_asymp <- sapply(paste0(prefix_distr, stats_type),
                        exists, where = asNamespace("sphunif"))

    # Exclude Stereo if p <= 3
    if (("Stereo" %in% stats_type) && (p <= 3)) {

      ind_asymp[which(stats_type == "Stereo")] <- FALSE

    }

    # Exclude statistics and throw warning
    if (!all(ind_asymp)) {

      warning(paste0(paste0("Omitting the following statistics with not ",
                            "implemented or known asymptotic distributions ",
                            "for p = ", p, ": "),
                     paste(stats_type[!ind_asymp], collapse = ", "), "."))
      stats_type <- stats_type[ind_asymp]

    }

  }

  # Number of statistics, counting statistics with vectorized parameters
  Sobolev_vk2 <- rbind(Sobolev_vk2)
  n_Sobolev_vk2 <- nrow(Sobolev_vk2)
  n_stats <- length(stats_type) +
    ifelse("Sobolev" %in% stats_type, n_Sobolev_vk2 - 1, 0)

  # Expand the names of statistics including the ones with vectorized parameters
  if ("Sobolev" %in% stats_type && n_Sobolev_vk2 > 1) {

    stats_type_vec <-
      strsplit(gsub(x = paste(stats_type, collapse = " "), pattern = "Sobolev",
                    replacement = paste0("Sobolev.", seq_len(n_Sobolev_vk2),
                                         collapse = " ")), split = " ")[[1]]

  } else {

    stats_type_vec <- stats_type

  }

  # Check dimension
  n_col_x <- ncol(x)
  if (n_col_x != n_stats) {

    if (n_col_x == 1) {

      x <- matrix(x, nrow = nrow(x), ncol = n_stats, byrow = FALSE)

    } else {

      stop(paste(strwrap(paste0(
        "Incompatible number of columns in x (", ncol(x), ") and tests in type",
        " (", n_stats, "; after expanding tests with vectorized arguments). ",
        "The tests being considered are: \"",
        paste(stats_type_vec, collapse = "\", \""), "\"."),
        width = 80, indent = 0, exdent = 2), collapse = "\n"))

    }

  }

  # Put x as a data frame
  x <- as.data.frame(x)
  names(x) <- stats_type_vec

  # Asymptotic or Monte Carlo distribution?
  if (approx == "asymp") {

    # Optional arguments
    args <- list("n" = n, "p" = p, "K_max" = K_max, "method" = method,
                 "thre" = 0, "Stephens" = Stephens, "dirs" = CCF09_dirs,
                 "regime" = CJ12_reg, "beta" = CJ12_beta, "K_Ajne" = K_Ajne,
                 "K_CCF09" = K_CCF09, "K_Kuiper" = K_Kuiper,
                 "K_Watson" = K_Watson, "K_Watson_1976" = K_Watson_1976,
                 "t" = Rothman_t, "rho" = Poisson_rho, "q" = Pycke_q,
                 "s" = Riesz_s, "vk2" = Sobolev_vk2, "kappa" = Softmax_kappa,
                 "Stein_K" = Stein_K, "Stein_cf" = Stein_cf, "a" = Stereo_a)
    names_args <- names(args)

    # Evaluate distributions
    distrs <- sapply(stats_type, function(distr) {

      # Additional arguments for the distribution
      name_distr <- paste0(prefix_distr, distr)
      distr_args <- args[names_args %in% names(formals(name_distr))]

      # Where to evaluate the distribution, depending on whether the
      # test is vectorized
      if (distr %in% c("Sobolev")) {

        x_distr <- as.matrix(x[, grep(pattern = distr, x = stats_type_vec,
                                      value = TRUE)])

      } else {

        x_distr <- x[[distr]]

      }

      # Call distribution
      do.call(what = name_distr, args = c(list(x = x_distr), distr_args))

    }, simplify = FALSE)

  } else if (approx == "MC") {

    # Sample statistics under the null
    if (is.null(stats_MC)) {

      stats_MC <- unif_stat_MC(n = n, type = stats_type, p = p, M = M,
                               r_H1 = NULL, crit_val = NULL,
                               return_stats = TRUE, stats_sorted = TRUE,
                               CCF09_dirs = CCF09_dirs, CJ12_reg = CJ12_reg,
                               cov_a = cov_a, Cressie_t = Cressie_t,
                               K_CCF09 = K_CCF09, Poisson_rho = Poisson_rho,
                               Pycke_q = Pycke_q, Rayleigh_m = Rayleigh_m,
                               Riesz_s = Riesz_s, Rothman_t = Rothman_t,
                               Sobolev_vk2 = Sobolev_vk2,
                               Softmax_kappa = Softmax_kappa,
                               Stein_K = Stein_K, Stein_cf = Stein_cf,
                               Stereo_a = Stereo_a, ...)$stats_MC

    } else {

      # Check if there is any missing statistic in stats_MC
      names_stats_MC <- names(stats_MC)
      missing_stats_MC <- !(stats_type_vec %in% names_stats_MC)
      if (any(missing_stats_MC)) {

        stop(paste0("Missing statistics in stats_MC: \"",
                    paste(stats_type_vec[missing_stats_MC],
                          collapse = "\", \""), "\"."))

      }

      # Sort for faster comparisons
      if (any(apply(stats_MC, 2, is.unsorted))) {

        stats_MC <- as.data.frame(sort_each_col(as.matrix(stats_MC)))
        names(stats_MC) <- names_stats_MC

      }

      # Match columns of stats_MC with stats_type
      stats_MC <- stats_MC[, pmatch(x = stats_type_vec, table = names_stats_MC),
                           drop = FALSE]

    }

    # Is it required to sort x?
    if (any(apply(x, 2, is.unsorted))) {

      # Indexes for sorting and then unsorting
      ind_1 <- sort_index_each_col(x)
      ind_2 <- sort_index_each_col(ind_1)

      # Approximate distributions
      distrs <- sapply(stats_type_vec, function(distr) {
        ecdf_bin(data = stats_MC[[distr]], sorted_x = x[[distr]][ind_1],
                 data_sorted = TRUE, efic = TRUE, divide_n = TRUE)
      }, simplify = FALSE)
      distrs <- distrs[ind_2, ]

    } else {

      # Approximate distributions
      distrs <- sapply(stats_type_vec, function(distr) {
        ecdf_bin(data = stats_MC[[distr]], sorted_x = x[[distr]],
                 data_sorted = TRUE, efic = TRUE, divide_n = TRUE)
      }, simplify = FALSE)

    }

  } else {

    stop("Wrong choice for approx.")

  }

  # Return data frame
  distrs <- as.data.frame(lapply(lapply(distrs, pmin, 1), pmax, 0))
  names(distrs) <- stats_type_vec
  return(distrs)

}
