

#' @title Circular and (hyper)spherical uniformity tests
#'
#' @description Implementation of several uniformity tests on the (hyper)sphere
#' \eqn{S^{p-1}:=\{{\bf x}\in R^p:||{\bf x}||=1\}}{
#' S^{p-1}:=\{x\in R^p:||x||=1\}}, \eqn{p\ge 2}, with calibration either in
#' terms of their asymptotic/exact distributions, if available, or Monte Carlo.
#'
#' \code{unif_test} receives a sample of directions
#' \eqn{{\bf X}_1,\ldots,{\bf X}_n\in S^{p-1}}{X_1,\ldots,X_n\in S^{p-1}} in
#' \emph{Cartesian coordinates}, except for the circular case (\eqn{p=2}) in
#' which the sample can be represented in terms of \emph{angles}
#' \eqn{\Theta_1,\ldots,\Theta_n\in [0, 2\pi)}.
#'
#' \code{unif_test} allows to perform several tests within a single call,
#' facilitating thus the exploration of a dataset by applying several tests.
#'
#' @param data sample to perform the test. A matrix of size \code{c(n, p)}
#' containing a sample of size \code{n} of directions (in Cartesian
#' coordinates) on \eqn{S^{p-1}}. Alternatively if \code{p = 2}, a matrix of
#' size \code{c(n, 1)} containing the \code{n} angles on \eqn{[0, 2\pi)} of the
#' circular sample on \eqn{S^{1}}. Other objects accepted are an array of size
#' \code{c(n, p, 1)} with directions (in Cartesian coordinates), or a vector of
#' size \code{n} or an array of size \code{c(n, 1, 1)} with angular data.
#' Must not contain \code{NA}'s.
#' @param type type of test to be applied. A character vector containing any of
#' the following types of tests, depending on the dimension \eqn{p}:
#' \itemize{
#'   \item Circular data: any of the names available at object
#'   \code{\link{avail_cir_tests}}.
#'   \item (Hyper)spherical data: any of the names available at object
#'   \code{\link{avail_sph_tests}}.
#' }
#' If \code{type = "all"} (default), then \code{type} is set as
#' \code{avail_cir_tests} or \code{avail_sph_tests}, depending on the value of
#' \eqn{p}.
#' @param p_value type of \eqn{p}-value computation. Either \code{"MC"} for
#' employing the approximation by Monte Carlo of the exact null distribution,
#' \code{"asymp"} (default) for the use of the asymptotic/exact null
#' distribution (if available), or \code{"crit_val"} for approximation by means
#' of the table of critical values \code{crit_val}.
#' @param alpha vector with significance levels. Defaults to
#' \code{c(0.10, 0.05, 0.01)}.
#' @inheritParams unif_stat_distr
#' @inheritParams unif_stat_MC
#' @param crit_val table with critical values for the tests, to be used if
#' \code{p_value = "crit_val"}. A data frame, with column names containing the
#' character vector \code{type} and rows corresponding to the significance
#' levels \code{alpha}, that results from extracting \code{$crit_val_MC} from
#' a call to \code{\link{unif_stat_MC}}. Internally computed if
#' \code{NULL} (default).
#' @inheritParams unif_stat
#' @param ... If \code{p_value = "MC"} or \code{p_value = "crit_val"}, optional
#' performance parameters to be passed to \code{\link{unif_stat_MC}}:
#' \code{chunks}, \code{cores}, and \code{seed}. If \code{p_value = "MC"},
#' additional parameters to \code{\link{unif_stat_distr}}.
#' @return If only a \bold{single test} is performed, a list with class
#' \code{htest} containing the following components:
#' \itemize{
#'   \item \code{statistic}: the value of the test statistic.
#'   \item \code{p.value}: the p-value of the test. If
#'   \code{p_value = "crit_val"}, an \code{NA}.
#'   \item \code{alternative}: a character string describing the alternative
#'   hypothesis.
#'   \item \code{method}: a character string indicating what type of test was
#'   performed.
#'   \item \code{data.name}: a character string giving the name of the data.
#'   \item \code{reject}: the rejection decision for the levels of significance
#'   \code{alpha}.
#'   \item \code{crit_val}: a vector with the critical values for the
#'   significance levels \code{alpha} used with \code{p_value = "MC"} or
#'   \code{p_value = "asymp"}.
#'   \item \code{param}: parameter(s) used in the test (if any).
#' }
#' If \bold{several tests} are performed, a \code{type}-named list with
#' entries for each test given by the above list.
#' @details
#' All the tests reject for large values of the test statistic, so the critical
#' values for the significance levels \code{alpha} correspond to the
#' \code{alpha}-upper quantiles of the null distribution of the test statistic.
#'
#' When \code{p_value = "asymp"}, tests that do not have an implemented or
#' known asymptotic are omitted, and a warning is generated.
#'
#' When \code{p_value = "MC"}, it is possible to have a progress bar indicating
#' the Monte Carlo simulation progress if \code{unif_test} is wrapped with
#' \code{\link[progressr:with_progress]{progressr::with_progress}} or if
#' \code{progressr::handlers(global = TRUE)} is invoked (once) by the user.
#' See the examples below. The progress bar is updated with the number of
#' finished chunks.
#'
#' All the statistics are continuous random variables except the
#' Hodges--Ajne statistic (\code{"Hodges_Ajne"}), the Cressie statistic
#' (\code{"Cressie"}), and the number of (different) uncovered spacings
#' (\code{"Num_uncover"}). These three statistics are discrete random variables.
#'
#' The Monte Carlo calibration for the CCF09 test is made conditionally
#' on the choice of \cr\code{CCF09_dirs}. That is, all the Monte
#' Carlo statistics share the same random directions.
#'
#' Except for \code{CCF09_dirs}, \code{K_CCF09}, and \code{CJ12_reg}, all the
#' test-specific parameters are vectorized.
#'
#' Descriptions and references for most of the tests are available
#' in García-Portugués and Verdebout (2018).
#' @references
#' García-Portugués, E. and Verdebout, T. (2018) An overview of uniformity
#' tests on the hypersphere. \emph{arXiv:1804.00286}.
#' \url{https://arxiv.org/abs/1804.00286}.
#' @examples
#' ## Asymptotic distribution
#'
#' # Circular data
#' n <- 10
#' samp_cir <- r_unif_cir(n = n)
#'
#' # Matrix
#' unif_test(data = samp_cir, type = "Ajne", p_value = "asymp")
#'
#' # Vector
#' unif_test(data = samp_cir[, 1], type = "Ajne", p_value = "asymp")
#'
#' # Array
#' unif_test(data = array(samp_cir, dim = c(n, 1, 1)), type = "Ajne",
#'           p_value = "asymp")
#' \donttest{
#' # Several tests
#' unif_test(data = samp_cir, type = avail_cir_tests, p_value = "asymp")
#' }
#' # Spherical data
#' n <- 10
#' samp_sph <- r_unif_sph(n = n, p = 3)
#'
#' # Array
#' unif_test(data = samp_sph, type = "Bingham", p_value = "asymp")
#'
#' # Matrix
#' unif_test(data = samp_sph[, , 1], type = "Bingham", p_value = "asymp")
#' \donttest{
#' # Several tests
#' unif_test(data = samp_sph, type = avail_sph_tests, p_value = "asymp")
#'
#' ## Monte Carlo
#'
#' # Circular data
#' unif_test(data = samp_cir, type = "Ajne", p_value = "MC")
#' unif_test(data = samp_cir, type = avail_cir_tests, p_value = "MC")
#'
#' # Spherical data
#' unif_test(data = samp_sph, type = "Bingham", p_value = "MC")
#' unif_test(data = samp_sph, type = avail_sph_tests, p_value = "MC")
#'
#' # Caching stats_MC
#' stats_MC_cir <- unif_stat_MC(n = nrow(samp_cir), p = 2)$stats_MC
#' stats_MC_sph <- unif_stat_MC(n = nrow(samp_sph), p = 3)$stats_MC
#' unif_test(data = samp_cir, type = avail_cir_tests,
#'           p_value = "MC", stats_MC = stats_MC_cir)
#' unif_test(data = samp_sph, type = avail_sph_tests, p_value = "MC",
#'           stats_MC = stats_MC_sph)
#'
#' ## Critical values
#'
#' # Circular data
#' unif_test(data = samp_cir, type = avail_cir_tests, p_value = "crit_val")
#'
#' # Spherical data
#' unif_test(data = samp_sph, type = avail_sph_tests, p_value = "crit_val")
#'
#' # Caching crit_val
#' crit_val_cir <- unif_stat_MC(n = n, p = 2)$crit_val_MC
#' crit_val_sph <- unif_stat_MC(n = n, p = 3)$crit_val_MC
#' unif_test(data = samp_cir, type = avail_cir_tests,
#'           p_value = "crit_val", crit_val = crit_val_cir)
#' unif_test(data = samp_sph, type = avail_sph_tests, p_value = "crit_val",
#'           crit_val = crit_val_sph)
#'
#' ## Specific arguments
#'
#' # Rothman
#' unif_test(data = samp_cir, type = "Rothman", Rothman_t = 0.5)
#'
#' # CCF09
#' unif_test(data = samp_sph, type = "CCF09", p_value = "MC",
#'           CCF09_dirs = samp_sph[1:2, , 1])
#' unif_test(data = samp_sph, type = "CCF09", p_value = "MC",
#'           CCF09_dirs = samp_sph[3:4, , 1])
#'
#' ## Using a progress bar when p_value = "MC"
#'
#' # Define a progress bar
#' require(progress)
#' require(progressr)
#' handlers(handler_progress(
#'   format = paste("(:spin) [:bar] :percent Iter: :current/:total Rate:",
#'                  ":tick_rate iter/sec ETA: :eta Elapsed: :elapsedfull"),
#'   clear = FALSE))
#'
#' # Call unif_test() within with_progress()
#' with_progress(
#'   unif_test(data = samp_sph, type = avail_sph_tests, p_value = "MC",
#'             chunks = 10, M = 1e3)
#' )
#'
#' # With several cores
#' with_progress(
#'   unif_test(data = samp_sph, type = avail_sph_tests, p_value = "MC",
#'             cores = 2, chunks = 10, M = 1e3)
#' )
#'
#' # Instead of using with_progress() each time, it is more practical to run
#' # handlers(global = TRUE)
#' # once to activate progress bars in your R session
#' }
#' @export
unif_test <- function(data, type = "all", p_value = "asymp",
                      alpha = c(0.10, 0.05, 0.01), M = 1e4, stats_MC = NULL,
                      crit_val = NULL, data_sorted = FALSE, K_max = 1e4,
                      method = "I", CCF09_dirs = NULL, CJ12_beta = 0,
                      CJ12_reg = 3, cov_a = 2 * pi, Cressie_t = 1 / 3,
                      K_CCF09 = 25, Poisson_rho = 0.5, Pycke_q = 0.5,
                      Rayleigh_m = 1, Riesz_s = 1, Rothman_t = 1 / 3,
                      Sobolev_vk2 = c(0, 0, 1), Softmax_kappa = 1,
                      Stereo_a = 0, ...) {

  # Read data's name
  data_name <- deparse(substitute(data))

  # Stop if NA's
  if (anyNA(data)) {

    stop("NAs present in data, please remove them.")

  }

  # If data is a vector, transform it to matrix
  if (is.vector(data)) {

    data <- matrix(data, ncol = 1)

  }

  # If data is an array, transform it to matrix
  d <- dim(data)
  l <- length(d)
  if (l == 3) {

    # As matrix
    data <- matrix(data[, , 1], nrow = d[1], ncol = d[2])

    # First slice only
    if (d[3] != 1) {

      message(paste("data is an array with more than one slice,",
                    "only the first one is employed."))

    }

  } else if (l > 3) {

    stop("data must be a vector, matrix, or a 3-dimensional array.")

  }

  # Sample size and dimension
  n <- nrow(data)
  d <- ncol(data)

  # Circular or spherical data?
  if (d == 1 || d == 2) {

    avail_stats <- avail_cir_tests
    prefix_stat <- "cir_stat_"
    p <- 2

    # As polar coordinates
    if (d == 2) {

      dim(data) <- c(n, d, 1)
      data <- X_to_Theta(X = data)

    }

  } else {

    avail_stats <- avail_sph_tests
    prefix_stat <- "sph_stat_"
    p <- d
    dim(data) <- c(n, d, 1) # As an array

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
          "Test with type = \"", type, "\" is unsupported. ",
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

  # Omit statistics that do not have asymptotic distribution
  if (p_value == "asymp") {

    # Search for p_*_stat_*
    ind_asymp <- sapply(paste0("p_", prefix_stat, stats_type),
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
      if (length(stats_type) == 0) {

        stop("No remaining statistics to use.")

      }

    }

  }

  # Number of statistics
  n_stats <- length(stats_type)

  ## Statistics

  # Sample random directions for the CCF09 test
  if ("CCF09" %in% stats_type && is.null(CCF09_dirs)) {

    CCF09_dirs <- r_unif_sph(n = 50, p = p, M = 1)[, , 1]

  }

  # Compute statistics
  stat <- unif_stat(data = data, type = stats_type, data_sorted = data_sorted,
                    CCF09_dirs = CCF09_dirs, CJ12_reg = CJ12_reg,
                    cov_a = cov_a, Cressie_t = Cressie_t, K_CCF09 = K_CCF09,
                    Poisson_rho = Poisson_rho, Pycke_q = Pycke_q,
                    Rayleigh_m = Rayleigh_m, Riesz_s = Riesz_s,
                    Rothman_t = Rothman_t, Sobolev_vk2 = Sobolev_vk2,
                    Softmax_kappa = Softmax_kappa, Stereo_a = Stereo_a)
  stats_type_vec <- names(stat) # We can have Sobolev.1, Sobolev.2, etc.

  # Update the number of statistics (to count those Sobolev.1, Sobolev.2, etc.)
  n_stats <- length(stats_type_vec)

  ## Calibration

  if (p_value == "crit_val") {

    # Get the critical values
    if (is.null(crit_val)) {

      crit_val <- unif_stat_MC(n = n, type = stats_type, p = p, M = M,
                               r_H1 = NULL, crit_val = NULL, alpha = alpha,
                               return_stats = FALSE, stats_sorted = FALSE,
                               CCF09_dirs = CCF09_dirs, CJ12_reg = CJ12_reg,
                               cov_a = cov_a, Cressie_t = Cressie_t,
                               K_CCF09 = K_CCF09, Poisson_rho = Poisson_rho,
                               Pycke_q = Pycke_q, Rayleigh_m = Rayleigh_m,
                               Riesz_s = Riesz_s, Rothman_t = Rothman_t,
                               Sobolev_vk2 = Sobolev_vk2, Softmax_kappa =
                                 Softmax_kappa, Stereo_a = Stereo_a,
                               ...)$crit_val_MC

    } else {

      # Check if there is any missing statistic in crit_val
      names_crit_val <- names(crit_val)
      missing_crit_val <- !(stats_type_vec %in% names_crit_val)
      if (any(missing_crit_val)) {

        stop(paste0("Missing statistics in crit_val: \"",
                    paste(stats_type_vec[missing_crit_val],
                          collapse = "\", \""), "\"."))

      }

      # Match columns of crit_val with stats_type_vec
      crit_val <- crit_val[, pmatch(x = stats_type_vec, table = names_crit_val),
                           drop = FALSE]

    }

    # Rejection?
    reject <- rbind(apply(crit_val, 1, function(x) x > stat))

    # p-values
    p_val <- matrix(NA, nrow = 1, ncol = n_stats)

  } else if (p_value == "MC") {

    # Get the stats_MC
    if (is.null(stats_MC)) {

      stats_MC <- unif_stat_MC(n = n, type = stats_type, p = p, M = M,
                               r_H1 = NULL, crit_val = NULL, alpha = alpha,
                               return_stats = TRUE, stats_sorted = TRUE,
                               CCF09_dirs = CCF09_dirs, CJ12_reg = CJ12_reg,
                               cov_a = cov_a, Cressie_t = Cressie_t,
                               K_CCF09 = K_CCF09, Poisson_rho = Poisson_rho,
                               Pycke_q = Pycke_q, Rayleigh_m = Rayleigh_m,
                               Riesz_s = Riesz_s, Rothman_t = Rothman_t,
                               Sobolev_vk2 = Sobolev_vk2, Softmax_kappa =
                                 Softmax_kappa, Stereo_a = Stereo_a,
                               ...)$stats_MC

    }

    # p-values
    p_val <- 1 - as.data.frame(sapply(stats_type_vec, function(distr) {
      ecdf_bin(data = stats_MC[[distr]], sorted_x = stat[[distr]],
               data_sorted = TRUE, efic = TRUE, divide_n = TRUE)
    }, simplify = FALSE))

    # Critical values
    crit_val <- apply(stats_MC, 2, quantile, probs = 1 - alpha, na.rm = TRUE)
    rownames(crit_val) <- alpha

    # Rejection?
    reject <- rbind(apply(crit_val, 1, function(x) stat >= x))

  } else if (p_value == "asymp") {

    # p-values
    p_val <- 1 - unif_stat_distr(x = stat, type = stats_type, p = p, n = n,
                                 approx = "asymp", stats_MC = NULL, M = M,
                                 K_max = K_max, method = method, thre = 0,
                                 CCF09_dirs = CCF09_dirs, CJ12_reg = CJ12_reg,
                                 CJ12_beta = CJ12_beta, cov_a = cov_a,
                                 Cressie_t = Cressie_t, K_CCF09 = K_CCF09,
                                 Poisson_rho = Poisson_rho, Pycke_q = Pycke_q,
                                 Rayleigh_m = Rayleigh_m, Riesz_s = Riesz_s,
                                 Rothman_t = Rothman_t, Sobolev_vk2 =
                                   Sobolev_vk2, Softmax_kappa = Softmax_kappa,
                                 Stereo_a = Stereo_a, ...)

    # Critical values
    crit_val <- as.data.frame(matrix(NA, nrow = length(alpha), ncol = n_stats))
    rownames(crit_val) <- alpha

    # Rejection?
    reject <- rbind(sapply(alpha, function(a) p_val < a))
    colnames(reject) <- alpha

  } else {

    stop(paste("Wrong choice for calibration, must be \"MC\", \"asymp\",",
               "or \"crit_val\"."))

  }

  ## Prepare return object

  # Create list of htest objects
  test <- vector(mode = "list", length = n_stats)
  names(test) <- stats_type_vec
  stats_type_rep <- gsub(x = stats_type_vec, pattern = "[.][0-9]+",
                         replacement = "")
  Sobolev_vk2 <- rbind(Sobolev_vk2)
  for (i in seq_along(stats_type_vec)) {

    # Parameter for statistics with vectorized arguments
    par_i <- as.numeric(strsplit(stats_type_vec[i], split = ".",
                                 fixed = TRUE)[[1]][2])

    # Type of test
    if (p == 2) {

      param <- switch(stats_type_rep[i],
         "Ajne" = NA,
         "Bakshaev" = NA,
         "Bingham" = NA,
         "CCF09" = CCF09_dirs,
         "Cressie" = c("Cressie_t" = ifelse(length(Cressie_t) == 1, Cressie_t,
                                            Cressie_t[par_i])),
         "FG01" = NA,
         "Gine_Fn" = NA,
         "Gine_Gn" = NA,
         "Gini" = NA,
         "Gini_squared" = NA,
         "Greenwood" = NA,
         "Hermans_Rasson" = NA,
         "Hodges_Ajne" = NA,
         "Kuiper" = NA,
         "Log_gaps" = NA,
         "Max_uncover" = c("cov_a" = ifelse(length(cov_a) == 1, cov_a,
                                            cov_a[par_i])),
         "Num_uncover" = c("cov_a" = ifelse(length(cov_a) == 1, cov_a,
                                            cov_a[par_i])),
         "PAD" = NA,
         "PCvM" = NA,
         "Poisson" = c("Poisson_rho" = ifelse(length(Poisson_rho) == 1,
                                              Poisson_rho, Poisson_rho[par_i])),
         "PRt" = c("Rothman_t" = ifelse(length(Rothman_t) == 1, Rothman_t,
                                        Rothman_t[par_i])),
         "Pycke" = NA,
         "Pycke_q" = c("Pycke_q" = ifelse(length(Pycke_q) == 1, Pycke_q,
                                          Pycke_q[par_i])),
         "Range" = NA,
         "Rao" = NA,
         "Rayleigh" = c("Rayleigh_m" = ifelse(length(Rayleigh_m) == 1,
                                              Rayleigh_m, Rayleigh_m[par_i])),
         "Riesz" = c("Riesz_s" = ifelse(length(Riesz_s) == 1, Riesz_s,
                                        Riesz_s[par_i])),
         "Rothman" = c("Rothman_t" = ifelse(length(Rothman_t) == 1, Rothman_t,
                                            Rothman_t[par_i])),
         "Sobolev" = c("Sobolev_vk2" = Sobolev_vk2[
           ifelse(nrow(Sobolev_vk2) == 1, 1, par_i), ]),
         "Softmax" = c("Softmax_kappa" = ifelse(length(Softmax_kappa) == 1,
                                                Softmax_kappa,
                                                Softmax_kappa[par_i])),
         "Vacancy" = c("cov_a" = ifelse(length(cov_a) == 1, cov_a,
                                        cov_a[par_i])),
         "Watson" = NA,
         "Watson_1976" = NA
      )

      method <- switch(stats_type_rep[i],
         "Ajne" = "Ajne test of circular uniformity",
         "Bakshaev" = "Bakshaev (2010) test of circular uniformity",
         "Bingham" = "Bingham test of circular uniformity",
         "CCF09" = paste("Cuesta-Albertos et al. (2009) test of circular",
                         "uniformity with k =", nrow(param)),
         "Cressie" = paste("Cressie test of circular uniformity with t =",
                           round(param, 3)),
         "FG01" = paste("Cramer-von Mises 4-point test of Feltz and",
                        "Goldin (2001)"),
         "Gine_Fn" = "Gine's Fn test of circular uniformity",
         "Gine_Gn" = "Gine's Gn test of circular uniformity",
         "Gini" = "Gini mean difference test of circular uniformity",
         "Gini_squared" = paste("Gini squared mean difference test of",
                                "circular uniformity"),
         "Greenwood" = "Greenwood spacings test of circular uniformity",
         "Hermans_Rasson" = "Hermans-Rasson test of circular uniformity",
         "Hodges_Ajne" = "Hodges-Ajne test of circular uniformity",
         "Kuiper" = "Kuiper test of circular uniformity",
         "Log_gaps" = paste("Darling\'s log-gaps spacings test of",
                            "circular uniformity"),
         "Max_uncover" = paste("Maximum uncovered spacing test of",
                               "circular uniformity with a =", round(param, 3)),
         "Num_uncover" = paste("Number of uncovered spacings test of",
                               "circular uniformity with a =", round(param, 3)),
         "PAD" = paste("Projected Anderson-Darling test of",
                       "circular uniformity"),
         "PCvM" = paste("Projected Cramer-von Mises test of",
                        "circular uniformity"),
         "Poisson" = paste("Poisson-kernel test of circular uniformity with",
                           "rho =", round(param, 3)),
         "PRt" = paste("Projected Rothman test of circular uniformity",
                       "with t =", round(param, 3)),
         "Pycke" = "Pycke test of circular uniformity",
         "Pycke_q" = paste("Pycke \"q-test\" of spherical uniformity with q =",
                           round(param, 3)),
         "Range" = "Range test of circular uniformity",
         "Rao" = "Rao\'s spacings test of circular uniformity",
         "Rayleigh" = paste0("Rayleigh test of circular uniformity with m =",
                             round(param, 3)),
         "Riesz" = paste("Warning! This is an experimental test not meant to",
                         "be used with s =", round(param, 3)),
         "Rothman" = paste("Rothman test of circular uniformity with t =",
                           round(param, 3)),
         "Sobolev" = paste("Finite Sobolev test of circular uniformity with",
                           "vk2 =", capture.output(dput(unname(param)))),
         "Softmax" = paste("Softmax test of circular uniformit with kappa =",
                           round(param, 3)),
         "Vacancy" = paste("Vacancy test of circular uniformity with a =",
                           round(param, 3)),
         "Watson" = "Watson test of circular uniformity",
         "Watson_1976" = "Watson (1976) test of circular uniformity"
      )

      alternative <- switch(stats_type_rep[i],
         "Ajne" = "any non-axial alternative to circular uniformity",
         "Bakshaev" = "any alternative to circular uniformity",
         "Bingham" = "scatter matrix different from constant",
         "Cressie" = paste("any alternative to circular uniformity",
                           "if t is irrational (conjectured)"),
         "CCF09" = "any alternative to circular uniformity",
         "FG01" = "any alternative to circular uniformity",
         "Gine_Fn" = "any alternative to circular uniformity",
         "Gine_Gn" = "any axial alternative to circular uniformity",
         "Gini" = "any alternative to circular uniformity",
         "Gini_squared" = "any alternative to circular uniformity",
         "Greenwood" = "any alternative to circular uniformity",
         "Hermans_Rasson" = "any alternative to circular uniformity",
         "Hodges_Ajne" = "any non-axial alternative to circular uniformity",
         "Kuiper" = "any alternative to circular uniformity",
         "Log_gaps" = "any alternative to circular uniformity",
         "Max_uncover" = "any alternative to circular uniformity",
         "Num_uncover" = paste("any alternative to circular uniformity",
                               "if a \u2264 2\u03c0"),
         "PAD" = "any alternative to circular uniformity",
         "PCvM" = "any alternative to circular uniformity",
         "Poisson" = "any alternative to circular uniformity for rho > 0",
         "PRt" = paste("any alternative to circular uniformity",
                       "if t is irrational"),
         "Pycke" = "any alternative to circular uniformity",
         "Pycke_q" = "any alternative to circular uniformity",
         "Range" = "any alternative to circular uniformity",
         "Rao" = "any alternative to circular uniformity",
         "Rayleigh" = "mean direction different from zero",
         "Rothman" = paste("any alternative to circular uniformity",
                           "if t is irrational"),
         "Riesz" = "unclear, experimental test",
         "Sobolev" = paste("alternatives in the Fourier subspace",
                           "with vk2 \u2260 0"),
         "Softmax" = "any alternative to circular uniformity for kappa > 0",
         "Vacancy" = "any alternative to circular uniformity",
         "Watson" = "any alternative to circular uniformity",
         "Watson_1976" = "unclear consistency"
      )

    } else {

      param <- switch(stats_type_rep[i],
         "Ajne" = NA,
         "Bakshaev" = NA,
         "Bingham" = NA,
         "CJ12" = NA,
         "CCF09" = CCF09_dirs,
         "Gine_Fn" = NA,
         "Gine_Gn" = NA,
         "PAD" = NA,
         "PCvM" = NA,
         "Poisson" = c("Poisson_rho" = ifelse(length(Poisson_rho) == 1,
                                              Poisson_rho, Poisson_rho[par_i])),
         "PRt" = c("Rothman_t" = ifelse(length(Rothman_t) == 1, Rothman_t,
                                        Rothman_t[par_i])),
         "Pycke" = NA,
         "Rayleigh" = NA,
         "Rayleigh_HD" = NA,
         "Riesz" = c("Riesz_s" = ifelse(length(Riesz_s) == 1, Riesz_s,
                                        Riesz_s[par_i])),
         "Sobolev" = c("Sobolev_vk2" = Sobolev_vk2[
           ifelse(nrow(Sobolev_vk2) == 1, 1, par_i), ]),
         "Softmax" = c("Softmax_kappa" = ifelse(length(Softmax_kappa) == 1,
                                                Softmax_kappa,
                                                Softmax_kappa[par_i])),
         "Stereo" = c("Stereo_a" = ifelse(length(Stereo_a) == 1,
                                           Stereo_a, Stereo_a[par_i]))
      )

      method <- switch(stats_type_rep[i],
         "Ajne" = "Ajne test of spherical uniformity",
         "Bakshaev" = "Bakshaev (2010) test of spherical uniformity",
         "Bingham" = "Bingham test of spherical uniformity",
         "CJ12" = "Cai and Jiang (2012) test of spherical uniformity",
         "CCF09" = paste("Cuesta-Albertos et al. (2009) test of spherical",
                         "uniformity with k =", nrow(param)),
         "Gine_Fn" = "Gine's Fn test of spherical uniformity",
         "Gine_Gn" = "Gine's Gn test of spherical uniformity",
         "PAD" = paste("Projected Anderson-Darling test of",
                        "spherical uniformity"),
         "PCvM" = paste("Projected Cramer-von Mises test of",
                        "spherical uniformity"),
         "Poisson" = paste("Poisson-kernel test of spherical uniformity with",
                           "rho =", round(param, 3)),
         "PRt" = paste("Projected Rothman test of spherical uniformity",
                       "with t =", round(param, 3)),
         "Pycke" = "Pycke test of spherical uniformity",
         "Rayleigh" = "Rayleigh test of spherical uniformity",
         "Rayleigh_HD" = paste("HD-standardized Rayleigh test of",
                               "spherical uniformity"),
         "Riesz" = paste("Warning! This is an experimental test not meant to",
                         "be used with s =", round(param, 3)),
         "Sobolev" = paste("Finite Sobolev test of spherical uniformity with",
                           "vk2 =", capture.output(dput(unname(param)))),
         "Softmax" = paste("Softmax test of spherical uniformit with kappa =",
                           round(param, 3)),
         "Stereo" = paste("Stereographic projection test of spherical",
                           "uniformity with a =", round(param, 3))
      )

      alternative <- switch(stats_type_rep[i],
         "Ajne" = "any non-axial alternative to spherical uniformity",
         "Bakshaev" = "any alternative to spherical uniformity",
         "Bingham" = "scatter matrix different from constant",
         "CJ12" = "unclear consistency",
         "CCF09" = "any alternative to spherical uniformity",
         "Gine_Fn" = "any alternative to spherical uniformity",
         "Gine_Gn" = "any axial alternative to spherical uniformity",
         "PAD" = "any alternative to spherical uniformity",
         "PCvM" = "any alternative to spherical uniformity",
         "Poisson" = "any alternative to spherical uniformity for rho > 0",
         "PRt" = paste("any alternative to spherical uniformity",
                       "if t is irrational"),
         "Pycke" = "any alternative to spherical uniformity",
         "Rayleigh" = "mean direction different from zero",
         "Rayleigh_HD" = "mean direction different from zero",
         "Riesz" = "unclear, experimental test",
         "Sobolev" = paste("alternatives in the spherical harmonics subspace",
                           "with vk2 \u2260 0"),
         "Softmax" = "any alternative to spherical uniformity for kappa > 0",
         "Stereo" = "any alternative to spherical uniformity for |a| < 1"
      )

    }

    # htest object
    test[[i]] <- list(statistic = c("statistic" = stat[1, i]),
                      p.value = p_val[1, i], alternative = alternative,
                      method = method, data.name = data_name,
                      reject = reject[i, ], crit_val = crit_val[, i],
                      param = param)
    class(test[[i]]) <- "htest"

  }

  # If there is only one test, return htest directly
  if (n_stats == 1) {

    test <- test[[1]]

  }
  return(test)

}


#' @title Available circular and (hyper)spherical uniformity tests
#'
#' @description Listing of the tests implemented in the \code{\link{sphunif}}
#' package.
#'
#' @return A character vector whose elements are valid inputs for the
#' \code{type} argument in \code{\link{unif_test}}, \code{\link{unif_stat}},
#' \code{\link{unif_stat_distr}}, and \code{\link{unif_stat_MC}}.
#' \code{avail_cir_tests} provides the available circular tests and
#' \code{avail_sph_tests} the (hyper)spherical tests.
#' @examples
#' # Circular tests
#' avail_cir_tests
#'
#' # Spherical tests
#' avail_sph_tests
#' @name avail_tests


#' @rdname avail_tests
#' @export
avail_cir_tests <- c("Ajne",
                     "Bakshaev",
                     "Bingham",
                     "Cressie",
                     "CCF09",
                     # "CJ12",
                     "FG01",
                     "Gine_Fn",
                     "Gine_Gn",
                     "Gini",
                     "Gini_squared",
                     "Greenwood",
                     "Hermans_Rasson",
                     "Hodges_Ajne",
                     "Kuiper",
                     "Log_gaps",
                     "Max_uncover",
                     "Num_uncover",
                     "PAD",
                     "PCvM",
                     "Poisson",
                     "PRt",
                     "Pycke",
                     "Pycke_q",
                     "Range",
                     "Rao",
                     "Rayleigh",
                     "Riesz",
                     "Rothman",
                     "Sobolev",
                     "Softmax",
                     "Vacancy",
                     "Watson",
                     "Watson_1976")


#' @rdname avail_tests
#' @export
avail_sph_tests <- c("Ajne",
                     "Bakshaev",
                     "Bingham",
                     "CCF09",
                     "CJ12",
                     "Gine_Fn",
                     "Gine_Gn",
                     "PAD",
                     "PCvM",
                     "Poisson",
                     "PRt",
                     "Pycke",
                     "Sobolev",
                     "Softmax",
                     "Stereo",
                     "Rayleigh",
                     "Rayleigh_HD",
                     "Riesz")
