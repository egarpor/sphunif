

#' @title Circular and (hyper)spherical uniformity statistics
#'
#' @description Implementation of several statistics for assessing uniformity
#' on the (hyper)sphere
#' \eqn{S^{p-1} := \{{\bf x} \in R^p : ||{\bf x}|| = 1\}}{
#' S^{p-1} := \{x \in R^p : ||x|| = 1\}}, \eqn{p\ge 2}, for a sample
#' \eqn{{\bf X}_1,\ldots,{\bf X}_n\in S^{p-1}}{X_1,\ldots,X_n\in S^{p-1}}.
#'
#' \code{unif_stat} receives a (several) sample(s) of directions in
#' \emph{Cartesian coordinates}, except for the circular case (\eqn{p=2}) in
#' which the sample(s) can be \emph{angles}
#' \eqn{\Theta_1,\ldots,\Theta_n\in [0, 2\pi)}.
#'
#' \code{unif_stat} allows to compute several statistics to several samples
#' within a single call, facilitating thus Monte Carlo experiments.
#'
#' @param data sample to compute the test statistic. An \bold{array} of size
#' \code{c(n, p, M)} containing \code{M} samples of size \code{n} of directions
#' (in Cartesian coordinates) on \eqn{S^{p-1}}. Alternatively, a
#' \bold{matrix} of size \code{c(n, M)} with the angles on \eqn{[0, 2\pi)} of
#' the \code{M} circular samples of size \code{n} on \eqn{S^{1}}. Other objects
#' accepted are an array of size \code{c(n, 1, M)} or a vector of size
#' \code{n} with angular data. Must not contain \code{NA}'s.
#' @inheritParams unif_test
#' @param data_sorted is the circular data sorted? If \code{TRUE}, certain
#' statistics are faster to compute. Defaults to \code{FALSE}.
#' @param Rayleigh_m integer \eqn{m} for the \eqn{m}-modal Rayleigh test.
#' Defaults to \code{m = 1} (the standard Rayleigh test).
#' @param cov_a \eqn{a_n = a / n} parameter used in the length of the arcs
#' of the coverage-based tests. Must be positive. Defaults to \code{2 * pi}.
#' @param Rothman_t \eqn{t} parameter for the Rothman test, a real in
#' \eqn{(0, 1)}. Defaults to \code{1 / 3}.
#' @param Cressie_t \eqn{t} parameter for the Cressie test, a real in
#' \eqn{(0, 1)}. Defaults to \code{1 / 3}.
#' @param Pycke_q \eqn{q} parameter for the Pycke "\eqn{q}-test", a real in
#' \eqn{(0, 1)}. Defaults to \code{1 / 2}.
#' @param Riesz_s \eqn{s} parameter for the \eqn{s}-Riesz test, a real in
#' \eqn{(0, 2)}. Defaults to \code{1}.
#' @param CCF09_dirs a matrix of size \code{c(n_proj, p)} containing
#' \code{n_proj} random directions (in Cartesian coordinates) on \eqn{S^{p-1}}
#' to perform the CCF09 test. If \code{NULL} (default), a sample of size
#' \code{n_proj = 50} directions is computed internally.
#' @param K_CCF09 integer giving the truncation of the series present in the
#' asymptotic distribution of the Kolmogorov-Smirnov statistic. Defaults to
#' \code{25}.
#' @param CJ12_reg type of asymptotic regime for CJ12 test, either \code{1}
#' (sub-exponential regime), \code{2} (exponential), or \code{3}
#' (super-exponential; default).
#' @param Poisson_rho \eqn{\rho} parameter for the Poisson test, a real in
#' \eqn{[0, 1)}. Defaults to \code{0.5}.
#' @param Softmax_kappa \eqn{\kappa} parameter for the Softmax test, a
#' non-negative real. Defaults to \code{1}.
#' @param Stein_K truncation \eqn{K} parameter for the Stein test, a positive
#' integer. Defaults to \code{10}.
#' @param Stein_cf logical indicating whether to use the characteristic
#' function in the Stein test. Defaults to \code{FALSE} (moment generating
#' function).
#' @param Stereo_a \eqn{a} parameter for the Stereo test, a real in
#' \eqn{[-1, 1]}. Defaults to \code{0}.
#' @param Sobolev_vk2 weights for the finite Sobolev test. A non-negative
#' vector or matrix. Defaults to \code{c(0, 0, 1)}.
#' @return A data frame of size \code{c(M, length(type))}, with column names
#' given by \code{type}, that contains the values of the test statistics.
#' @details
#' Except for \code{CCF09_dirs}, \code{K_CCF09}, and \code{CJ12_reg}, all the
#' test-specific parameters are vectorized.
#'
#' Descriptions and references for most of the statistics are available
#' in García-Portugués and Verdebout (2018).
#' @references
#' García-Portugués, E. and Verdebout, T. (2018) An overview of uniformity
#' tests on the hypersphere. \emph{arXiv:1804.00286}.
#' \doi{10.48550/arXiv.1804.00286}.
#' @examples
#' ## Circular data
#'
#' # Sample
#' n <- 10
#' M <- 2
#' Theta <- r_unif_cir(n = n, M = M)
#'
#' # Matrix
#' unif_stat(data = Theta, type = "all")
#'
#' # Array
#' unif_stat(data = array(Theta, dim = c(n, 1, M)), type = "all")
#'
#' # Vector
#' unif_stat(data = Theta[, 1], type = "all")
#'
#' ## Spherical data
#'
#' # Circular sample in Cartesian coordinates
#' n <- 10
#' M <- 2
#' X <- array(dim = c(n, 2, M))
#' for (i in 1:M) X[, , i] <- cbind(cos(Theta[, i]), sin(Theta[, i]))
#'
#' # Array
#' unif_stat(data = X, type = "all")
#'
#' # High-dimensional data
#' X <- r_unif_sph(n = n, p = 3, M = M)
#' unif_stat(data = X, type = "all")
#'
#' ## Specific arguments
#'
#' # Rothman
#' unif_stat(data = Theta, type = "Rothman", Rothman_t = 0.5)
#'
#' # CCF09
#' unif_stat(data = X, type = "CCF09", CCF09_dirs = X[, , 1])
#' unif_stat(data = X, type = "CCF09", CCF09_dirs = X[, , 1], K_CCF09 = 1)
#'
#' # CJ12
#' unif_stat(data = X, type = "CJ12", CJ12_reg = 3)
#' unif_stat(data = X, type = "CJ12", CJ12_reg = 1)
#' @export
unif_stat <- function(data, type = "all", data_sorted = FALSE,
                      CCF09_dirs = NULL, CJ12_reg = 3, cov_a = 2 * pi,
                      Cressie_t = 1 / 3, K_CCF09 = 25, Poisson_rho = 0.5,
                      Pycke_q = 0.5, Rayleigh_m = 1, Riesz_s = 1,
                      Rothman_t = 1 / 3, Sobolev_vk2 = c(0, 0, 1),
                      Softmax_kappa = 1, Stein_K = 10, Stein_cf = FALSE,
                      Stereo_a = 0) {

  # Stop if NA's
  if (anyNA(data)) {

    stop("NAs present in data, please remove them.")

  }

  # If data is a vector, transform it to matrix
  if (is.vector(data)) {

    data <- matrix(data, ncol = 1)

  }

  # Decode the type of data: vector, matrix, or array
  d <- dim(data)

  # Check sample size
  n <- d[1]
  if (n < 2) {

    stop("Sample size is one!")

  }

  # Array or matrix?
  if (length(d) == 3) {

    # Dimension and number of samples
    p <- d[2]
    M <- d[3]

    # Convert to polar coordinates?
    if (p == 2) {

      data <- X_to_Theta(X = data)

    } else if (p == 1) {

      data <- data[, 1, ]
      p <- 2

    }

  } else {

    # Dimension and number of samples
    p <- 2
    M <- d[2]

  }

  # Available statistics
  if (p == 2) {

    avail_stats <- avail_cir_tests

  } else {

    avail_stats <- avail_sph_tests

  }

  # Get the type of statistics
  if (is.character(type)) {

    type <- gsub("-", "_", type)
    type <- unique(type)
    if ("all" %in% type) {

      stats_type <- avail_stats

    } else {

      # Allow KS, CvM, and AD in unif_stat() (not in unif_test())
      if (p == 2) {

        avail_stats <- c(avail_stats, "KS", "CvM", "AD")

      }
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

  # Create a data frame for the statistics
  n_stats <- length(stats_type)
  stats <- matrix(NA, nrow = M, ncol = n_stats)
  colnames(stats) <- stats_type
  stats <- as.data.frame(stats)

  # Compute efficiently several statistics
  if (p == 2) { # p == 2

    # Statistics using gaps and sorted data
    stats_using_gaps <- c("Gini", "Gini_squared", "Greenwood", "Log_gaps",
                          "Max_uncover", "Num_uncover", "Range", "Rao",
                          "Vacancy")
    stats_using_sorted_data <- c("Cressie", "FG01", "Hodges_Ajne", "Kuiper",
                                 "Watson", "Watson_1976", "KS", "CvM", "AD",
                                 stats_using_gaps)

    # Statistics using the shortest angles matrix Psi
    stats_using_Psi <- c("Ajne", "Bakshaev", "Gine_Fn", "Gine_Gn",
                         "Hermans_Rasson", "PAD", "PCvM", "Poisson", "PRt",
                         "Pycke", "Pycke_q", "Rothman", "Riesz", "Sobolev",
                         "Softmax", "Stereo")

    # Statistics with vectorized parameters
    stats_vectorized <- c("Cressie", c("Max_uncover", "Num_uncover", "Vacancy"),
                          "Rayleigh", "Riesz", c("Rothman", "PRt"), "Poisson",
                          "Pycke_q", "Sobolev", "Softmax")
    param_vectorized <- c("Cressie_t", rep("cov_a", 3), "Rayleigh_m", "Riesz_s",
                          rep("Rothman_t", 2), "Poisson_rho", "Pycke_q",
                          "Sobolev_vk2", "Softmax_kappa")

    # Evaluate which statistics to apply
    run_test <- as.list(c(avail_cir_tests, "KS", "CvM", "AD") %in% stats_type)
    names(run_test) <- c(avail_cir_tests, "KS", "CvM", "AD")

    ## Sorted-data statistics

    # Is it worth to sort Theta? Only if there is more than one statistic
    # using sorted Theta and if it has not already been done
    # CAUTION: replacement of data with sorted_data!
    n_stats_type_sorted <- sum(stats_type %in% stats_using_sorted_data)
    if (!data_sorted && (n_stats_type_sorted > 1)) {

      data <- sort_each_col(data)
      data_sorted <- TRUE

    }

    # Statistics
    if (run_test$Cressie) {

      Cressie <- matrix(NA, nrow = M, ncol = length(Cressie_t))
      for (i in seq_along(Cressie_t)) {

        Cressie[, i] <- cir_stat_Cressie(Theta = data, t = Cressie_t[i],
                                         sorted = data_sorted)

      }
      stats$Cressie <- Cressie

    }
    if (run_test$FG01) {

      stats$FG01 <- cir_stat_FG01(Theta = data, sorted = data_sorted)

    }
    if (run_test$Hodges_Ajne) {

      stats$Hodges_Ajne <- cir_stat_Hodges_Ajne(Theta = data, asymp_std = FALSE,
                                                sorted = data_sorted,
                                                use_Cressie = TRUE)

    }
    if (run_test$Kuiper) {

      stats$Kuiper <- cir_stat_Kuiper(Theta = data, sorted = data_sorted,
                                      KS = FALSE, Stephens = FALSE)

    }
    if (run_test$Watson) {

      stats$Watson <- cir_stat_Watson(Theta = data, sorted = data_sorted,
                                      CvM = FALSE, Stephens = FALSE)

    }
    if (run_test$KS) {

      stats$KS <- cir_stat_Kuiper(Theta = data, sorted = data_sorted,
                                  KS = TRUE, Stephens = FALSE)

    }
    if (run_test$CvM) {

      stats$CvM <- cir_stat_Watson(Theta = data, sorted = data_sorted,
                                   CvM = TRUE, Stephens = FALSE)

    }
    if (run_test$AD) {

      stats$AD <- cir_stat_PAD(Theta = data, sorted = data_sorted, AD = TRUE)

    }
    if (run_test$Watson_1976) {

      stats$Watson_1976 <- cir_stat_Watson_1976(Theta = data,
                                                sorted = data_sorted,
                                                minus = FALSE)

    }

    ## Gaps-based statistics

    # Is it worth to compute the circular gaps? Only if there is more than
    # one statistic using them
    n_stats_type_gaps <- sum(stats_type %in% stats_using_gaps)
    if (n_stats_type_gaps > 1) {

      gaps <- cir_gaps(Theta = data, sorted = data_sorted)
      gaps_in_Theta <- TRUE

    } else {

      gaps <- data
      gaps_in_Theta <- FALSE

    }
    if (run_test$Range) {

      stats$Range <- cir_stat_Range(Theta = gaps, gaps_in_Theta = gaps_in_Theta,
                                    max_gap = TRUE, sorted = data_sorted)

    }
    if (run_test$Rao) {

      stats$Rao <- cir_stat_Rao(Theta = gaps, gaps_in_Theta = gaps_in_Theta,
                                sorted = data_sorted)

    }
    if (run_test$Gini) {

      stats$Gini <- cir_stat_Gini(Theta = gaps, gaps_in_Theta = gaps_in_Theta,
                                  sorted = data_sorted)

    }
    if (run_test$Gini_squared) {

      stats$Gini_squared <- cir_stat_Gini_squared(Theta = gaps,
                                                  gaps_in_Theta = gaps_in_Theta,
                                                  sorted = data_sorted)

    }
    if (run_test$Greenwood) {

      stats$Greenwood <- cir_stat_Greenwood(Theta = gaps,
                                            gaps_in_Theta = gaps_in_Theta,
                                            sorted = data_sorted)

    }
    if (run_test$Log_gaps) {

      stats$Log_gaps <- cir_stat_Log_gaps(Theta = gaps,
                                          gaps_in_Theta = gaps_in_Theta,
                                          sorted = data_sorted)

    }
    if (run_test$Max_uncover) {

      Max_uncover <- matrix(NA, nrow = M, ncol = length(cov_a))
      for (i in seq_along(cov_a)) {

        Max_uncover[, i] <- cir_stat_Max_uncover(Theta = gaps, a = cov_a[i],
                                                 gaps_in_Theta = gaps_in_Theta,
                                                 sorted = data_sorted)

      }
      stats$Max_uncover <- Max_uncover

    }
    if (run_test$Num_uncover) {

      Num_uncover <- matrix(NA, nrow = M, ncol = length(cov_a))
      for (i in seq_along(cov_a)) {

        Num_uncover[, i] <- cir_stat_Num_uncover(Theta = gaps, a = cov_a[i],
                                                 gaps_in_Theta = gaps_in_Theta,
                                                 sorted = data_sorted)

      }
      stats$Num_uncover <- Num_uncover

    }
    if (run_test$Vacancy) {

      Vacancy <- matrix(NA, nrow = M, ncol = length(cov_a))
      for (i in seq_along(cov_a)) {

        Vacancy[, i] <- cir_stat_Vacancy(Theta = gaps, a = cov_a[i],
                                         gaps_in_Theta = gaps_in_Theta,
                                         sorted = data_sorted)

      }
      stats$Vacancy <- Vacancy

    }

    ## Other statistics

    if (run_test$Rayleigh) {

      Rayleigh <- matrix(NA, nrow = M, ncol = length(Rayleigh_m))
      for (i in seq_along(Rayleigh_m)) {

        Rayleigh[, i] <- cir_stat_Rayleigh(Theta = data, m = Rayleigh_m[i])

      }
      stats$Rayleigh <- Rayleigh

    }
    if (run_test$Bingham) {

      stats$Bingham <- cir_stat_Bingham(Theta = data)

    }
    if (run_test$CCF09) {

      # Sample random directions
      if (is.null(CCF09_dirs)) {

        CCF09_dirs <- r_unif_sph(n = 50, p = 2, M = 1)[, , 1]

      }
      stats$CCF09 <- cir_stat_CCF09(Theta = data, dirs = CCF09_dirs,
                                    K_CCF09 = K_CCF09, original = FALSE)

    }

    ## Psi-based statistics

    # Is it worth to compute Psi? Only if there is more than one statistic
    # using it (after discounting equivalent tests if they are simultaneously
    # present) AND if 0.5 * n * (n - 1) * M is not too large.
    # CAUTION: replacement of data with Psi!

    # Number of statistics
    n_stats_type_Psi <- sum(stats_type %in% stats_using_Psi)

    # Add statistics with vectorized parameters
    n_stats_type_Psi <- n_stats_type_Psi + (
      run_test$Riesz * (length(Riesz_s) - 1) +
      (run_test$PRt || run_test$Rothman) * (length(Rothman_t) - 1) +
      run_test$Poisson * (length(Poisson_rho) - 1) +
      run_test$Softmax * (length(Softmax_kappa) - 1) +
      run_test$Sobolev * (nrow(rbind(Sobolev_vk2)) - 1))

    # Remove equivalent tests
    n_stats_type_Psi <- n_stats_type_Psi - (
      (run_test$Watson && run_test$PCvM) +
      (run_test$Rothman && run_test$PRt) +
      (any(Riesz_s == 2) && run_test$Riesz) +
      (any(Riesz_s == 0) && run_test$Pycke && run_test$Riesz) +
      (any(Riesz_s == 1) && run_test$Bakshaev && run_test$Riesz) +
      (any(Poisson_rho == 0) && run_test$Poisson) +
      (any(Softmax_kappa == 0) && run_test$Softmax) +
      (run_test$Gine_Fn && run_test$Gine_Gn && run_test$Ajne))

    # Compute Psi
    if (n_stats_type_Psi > 1 && 0.5 * n * (n - 1) * M <= 1e8) {

      dim(data) <- c(n, 1, M)
      data <- Psi_mat(data = data)
      Psi_in_Theta <- TRUE

    } else {

      Psi_in_Theta <- FALSE

    }

    # Statistics
    if (run_test$Ajne) {

      stats$Ajne <- cir_stat_Ajne(Theta = data, Psi_in_Theta = Psi_in_Theta)

    }
    if (run_test$Bakshaev) {

      stats$Bakshaev <- cir_stat_Riesz(Theta = data,
                                       Psi_in_Theta = Psi_in_Theta, s = 1)

    }
    if (run_test$Riesz) {

      Riesz <- matrix(NA, nrow = M, ncol = length(Riesz_s))
      for (i in seq_along(Riesz_s)) {

        if (Riesz_s[i] == 1 && run_test$Bakshaev) {

          Riesz[, i] <- stats$Bakshaev

        } else if (Riesz_s[i] == 2 &&
                   ((run_test$Rayleigh && any(Rayleigh_m == 1)) ||
                    !Psi_in_Theta)) {

          if (run_test$Rayleigh) {

            Riesz[, i] <- stats$Rayleigh[, which(Rayleigh_m == 1)]

          } else {

            Riesz[, i] <- cir_stat_Rayleigh(Theta = data, m = 1)

          }

        } else {

          Riesz[, i] <- (1 - 2 * (Riesz_s[i] < 0)) *
            cir_stat_Riesz(Theta = data, Psi_in_Theta = Psi_in_Theta,
                           s = Riesz_s[i])

        }

      }
      stats$Riesz <- Riesz

    }
    if (run_test$Gine_Gn) {

      stats$Gine_Gn <- cir_stat_Gine_Gn(Theta = data,
                                        Psi_in_Theta = Psi_in_Theta)

    }
    if (run_test$Gine_Fn) {

      if (run_test$Ajne && run_test$Gine_Gn) {

        stats$Gine_Fn <- 4 * stats$Ajne + stats$Gine_Gn

      } else {

        stats$Gine_Fn <- cir_stat_Gine_Fn(Theta = data,
                                          Psi_in_Theta = Psi_in_Theta)

      }

    }
    if (run_test$Hermans_Rasson) {

      stats$Hermans_Rasson <-
        cir_stat_Hermans_Rasson(Theta = data, Psi_in_Theta = Psi_in_Theta)

    }
    if (run_test$PAD) {

      stats$PAD <- cir_stat_PAD(Theta = data, Psi_in_Theta = Psi_in_Theta,
                                AD = FALSE, sorted = FALSE)

    }
    if (run_test$PCvM) {

      if (run_test$Watson) {

        stats$PCvM <- 2 * stats$Watson

      } else {

        stats$PCvM <- cir_stat_PCvM(Theta = data, Psi_in_Theta = Psi_in_Theta)

      }

    }
    if (run_test$Rothman) {

      Rothman <- matrix(NA, nrow = M, ncol = length(Rothman_t))
      for (i in seq_along(Rothman_t)) {

        Rothman[, i] <- cir_stat_Rothman(Theta = data,
                                         Psi_in_Theta = Psi_in_Theta,
                                         t = Rothman_t[i])

      }
      stats$Rothman <- Rothman

    }
    if (run_test$PRt) {

      if (run_test$Rothman) {

        stats$PRt <- stats$Rothman

      } else {

        PRt <- matrix(NA, nrow = M, ncol = length(Rothman_t))
        for (i in seq_along(Rothman_t)) {

          PRt[, i] <- cir_stat_PRt(Theta = data, Psi_in_Theta = Psi_in_Theta,
                                   t = Rothman_t[i])

        }
        stats$PRt <- PRt

      }

    }
    if (run_test$Pycke) {

      if (run_test$Riesz && any(Riesz_s == 0)) {

        stats$Pycke <- (2 * n) / (n - 1) * stats$Riesz[, which(Riesz_s == 0)]

      } else {

        stats$Pycke <- cir_stat_Pycke(Theta = data, Psi_in_Theta = Psi_in_Theta)

      }

    }
    if (run_test$Pycke_q) {

      Pycke_q_stat <- matrix(NA, nrow = M, ncol = length(Pycke_q))
      for (i in seq_along(Pycke_q)) {

        Pycke_q_stat[, i] <- cir_stat_Pycke_q(Theta = data,
                                              Psi_in_Theta = Psi_in_Theta,
                                              q = Pycke_q[i])

      }
      stats$Pycke_q <- Pycke_q_stat

    }
    if (run_test$Poisson) {

      Poisson <- matrix(NA, nrow = M, ncol = length(Poisson_rho))
      for (i in seq_along(Poisson_rho)) {

        if (Poisson_rho[i] == 0) {

          if (run_test$Rayleigh && any(Rayleigh_m == 1)) {

            Poisson[, i] <- stats$Rayleigh[, which(Rayleigh_m == 1)]

          } else {

            Poisson[, i] <- cir_stat_Rayleigh(Theta = data, m = 1)

          }

        } else {

          Poisson[, i] <- cir_stat_Poisson(Theta = data,
                                           Psi_in_Theta = Psi_in_Theta,
                                           rho = Poisson_rho[i])

        }

      }
      stats$Poisson <- Poisson

    }
    if (run_test$Softmax) {

      Softmax <- matrix(NA, nrow = M, ncol = length(Softmax_kappa))
      for (i in seq_along(Softmax_kappa)) {

        if (Softmax_kappa[i] == 0) {

          if (run_test$Rayleigh && any(Rayleigh_m == 1)) {

            Softmax[, i] <- stats$Rayleigh[, which(Rayleigh_m == 1)]

          } else {

            Softmax[, i] <- cir_stat_Rayleigh(Theta = data, m = 1)

          }

        } else {

            Softmax[, i] <- cir_stat_Softmax(Theta = data,
                                             Psi_in_Theta = Psi_in_Theta,
                                             kappa = Softmax_kappa[i])

        }

      }
      stats$Softmax <- Softmax

    }
    if (run_test$Sobolev) {

      stats$Sobolev <- cir_stat_Sobolev(Theta = data,
                                        Psi_in_Theta = Psi_in_Theta,
                                        vk2 = Sobolev_vk2)

    }
    if (run_test$Stein) {

      Stein_vk2 <- weights_dfs_Sobolev(p = 2, K_max = Stein_K, thre = 0,
                                       type = "Stein", verbose = FALSE,
                                       Stein_cf = Stein_cf)$weights
      stats$Stein <- cir_stat_Sobolev(Theta = data, Psi_in_Theta = Psi_in_Theta,
                                      vk2 = Stein_vk2)

    }

  } else { # p >= 3

    # Statistics using the shortest angles matrix Psi
    stats_using_Psi <- c("Ajne", "Bakshaev", "CJ12", "Gine_Fn", "Gine_Gn",
                         "PAD", "PCvM", "PRt", "Poisson", "Pycke", "Riesz",
                         "Sobolev", "Softmax", "Stereo")

    # Statistics with vectorized parameters
    stats_vectorized <- c("Poisson", "PRt", "Riesz", "Sobolev", "Softmax",
                          "Stereo")
    param_vectorized <- c("Poisson_rho", "Rothman_t", "Riesz_s",
                          "Sobolev_vk2", "Softmax_kappa", "Stereo_a")

    # Evaluate which statistics to apply
    run_test <- as.list(avail_sph_tests %in% stats_type)
    names(run_test) <- avail_sph_tests

    ## Other statistics

    if (run_test$Bingham) {

      stats$Bingham <- sph_stat_Bingham(X = data)

    }
    if (run_test$CCF09) {

      # Sample random directions
      if (is.null(CCF09_dirs)) {

        CCF09_dirs <- r_unif_sph(n = 50, p = p, M = 1)[, , 1]

      }
      stats$CCF09 <- sph_stat_CCF09(X = data, dirs = CCF09_dirs,
                                    K_CCF09 = K_CCF09, original = FALSE)

    }
    if (run_test$Rayleigh) {

      stats$Rayleigh <- sph_stat_Rayleigh(X = data)

    }
    if (run_test$Rayleigh_HD) {

      if (run_test$Rayleigh) {

        stats$Rayleigh_HD <- (stats$Rayleigh - p) / sqrt(2 * p)

      } else {

        stats$Rayleigh_HD <- sph_stat_Rayleigh_HD(X = data)

      }

    }

    ## Psi-based statistics

    # Is it worth to compute Psi? Only if there is more than one statistic
    # using it (after discounting equivalent tests if they are simultaneously
    # present) AND if 0.5 * n * (n - 1) * M is not too large.
    # CAUTION: replacement of data with Psi!

    # Number of statistics
    n_stats_type_Psi <- sum(stats_type %in% stats_using_Psi)

    # Add statistics with vectorized parameters
    n_stats_type_Psi <- n_stats_type_Psi + (
      run_test$Riesz * (length(Riesz_s) - 1) +
      run_test$PRt * (length(Rothman_t) - 1) +
      run_test$Poisson * (length(Poisson_rho) - 1) +
      run_test$Softmax * (length(Softmax_kappa) - 1) +
      run_test$Stereo * (length(Stereo_a) - 1) +
      run_test$Sobolev * (nrow(rbind(Sobolev_vk2)) - 1))

    # Remove equivalent tests
    n_stats_type_Psi <- n_stats_type_Psi - (
      (any(Riesz_s == 2) && run_test$Riesz) +
      (any(Riesz_s == 0) && run_test$Pycke && run_test$Riesz) +
      (any(Riesz_s == 1) && run_test$Bakshaev && run_test$Riesz) +
      (any(Poisson_rho == 0) && run_test$Poisson) +
      (any(Softmax_kappa == 0) && run_test$Softmax) +
      (p == 3 && run_test$PCvM && run_test$Bakshaev) +
      (run_test$Gine_Fn && run_test$Gine_Gn && run_test$Ajne))

    # Compute Psi
    if (n_stats_type_Psi > 1 && 0.5 * n * (n - 1) * M <= 1e8) {

      data <- Psi_mat(data = data)
      dim(data) <- c(dim(data), 1)
      Psi_in_X <- TRUE

    } else {

      Psi_in_X <- FALSE

    }

    # Statistics
    if (run_test$Ajne) {

      stats$Ajne <- sph_stat_Ajne(X = data, Psi_in_X = Psi_in_X)

    }
    if (run_test$Bakshaev) {

      stats$Bakshaev <- sph_stat_Riesz(X = data, Psi_in_X = Psi_in_X, p = p,
                                       s = 1)

    }
    if (run_test$Riesz) {

      Riesz <- matrix(NA, nrow = M, ncol = length(Riesz_s))
      for (i in seq_along(Riesz_s)) {

        if (Riesz_s[i] == 1 && run_test$Bakshaev) {

          Riesz[, i] <- stats$Bakshaev

        } else if (Riesz_s[i] == 2 && (run_test$Rayleigh || !Psi_in_X)) {

          if (run_test$Rayleigh) {

            Riesz[, i] <- (2 / p) * stats$Rayleigh

          } else {

            Riesz[, i] <- (2 / p) * sph_stat_Rayleigh(X = data)

          }

        } else {

          Riesz[, i] <- (1 - 2 * (Riesz_s[i] < 0)) *
            sph_stat_Riesz(X = data, Psi_in_X = Psi_in_X, p = p, s = Riesz_s[i])

        }

      }
      stats$Riesz <- Riesz

    }
    if (run_test$CJ12) {

      if (Psi_in_X) {

        stats$CJ12 <- sph_stat_CJ12(X = data, Psi_in_X = TRUE, p = p,
                                    regime = CJ12_reg)

      } else {

        stats$CJ12 <- sph_stat_CJ12(X = data, Psi_in_X = FALSE, p = p,
                                    regime = CJ12_reg)

      }

    }
    if (run_test$Gine_Gn) {

      stats$Gine_Gn <- sph_stat_Gine_Gn(X = data, Psi_in_X = Psi_in_X, p = p)

    }
    if (run_test$Gine_Fn) {

      if (run_test$Ajne && run_test$Gine_Gn) {

        stats$Gine_Fn <- 4 * stats$Ajne + stats$Gine_Gn

      } else {

        stats$Gine_Fn <- sph_stat_Gine_Fn(X = data, Psi_in_X = Psi_in_X, p = p)

      }

    }
    if (run_test$PAD) {

      stats$PAD <- sph_stat_PAD(X = data, Psi_in_X = Psi_in_X, p = p)

    }
    if (run_test$PCvM) {

      if (run_test$Bakshaev && p == 3) {

        stats$PCvM <- 0.125 * stats$Bakshaev

      } else {

        stats$PCvM <- sph_stat_PCvM(X = data, Psi_in_X = Psi_in_X, p = p)

      }

    }
    if (run_test$PRt) {

      PRt <- matrix(NA, nrow = M, ncol = length(Rothman_t))
      for (i in seq_along(Rothman_t)) {

        PRt[, i] <- sph_stat_PRt(X = data, Psi_in_X = Psi_in_X, p = p,
                                 t = Rothman_t[i])

      }
      stats$PRt <- PRt

    }
    if (run_test$Pycke) {

      if (run_test$Riesz && any(Riesz_s == 0)) {

        if (p == 3) {

          stats$Pycke <- n / (2 * pi * (n - 1)) *
            (stats$Riesz[, which(Riesz_s == 0)] - (log(4) - 1) / 2)

        } else {

          warning(paste("Pycke statistic is only defined for p = 2,3.",
                        "Using Riesz statistic with s = 0 instead,",
                        "which behaves consistently across dimensions."))
          stats$Pycke <- stats$Riesz

        }

      } else {

        if (Psi_in_X) {

          stats$Pycke <- sph_stat_Pycke(X = data, Psi_in_X = TRUE, p = p)

        } else {

          stats$Pycke <- sph_stat_Pycke(X = data, Psi_in_X = FALSE, p = p)

        }

      }

    }
    if (run_test$Stereo) {

      Stereo <- matrix(NA, nrow = M, ncol = length(Stereo_a))
      for (i in seq_along(Stereo_a)) {

        Stereo[, i] <- sph_stat_Stereo(X = data, Psi_in_X = Psi_in_X, p = p,
                                       a = Stereo_a[i])

      }
      stats$Stereo <- Stereo

    }
    if (run_test$Poisson) {

      Poisson <- matrix(NA, nrow = M, ncol = length(Poisson_rho))
      for (i in seq_along(Poisson_rho)) {

        if (Poisson_rho[i] == 0) {

          if (run_test$Rayleigh) {

            Poisson[, i] <- stats$Rayleigh

          } else {

            Poisson[, i] <- sph_stat_Rayleigh(X = data)

          }

        } else {

          Poisson[, i] <- sph_stat_Poisson(X = data, Psi_in_X = Psi_in_X, p = p,
                                           rho = Poisson_rho[i])

        }

      }
      stats$Poisson <- Poisson

    }
    if (run_test$Softmax) {

      Softmax <- matrix(NA, nrow = M, ncol = length(Softmax_kappa))
      for (i in seq_along(Softmax_kappa)) {

        if (Softmax_kappa[i] == 0) {

          if (run_test$Rayleigh) {

            Softmax[, i] <- stats$Rayleigh

          } else {

            Softmax[, i] <- sph_stat_Rayleigh(X = data)

          }

        } else {

          Softmax[, i] <- sph_stat_Softmax(X = data, Psi_in_X = Psi_in_X,
                                           p = p, kappa = Softmax_kappa[i])

        }

      }
      stats$Softmax <- Softmax

    }
    if (run_test$Sobolev) {

      stats$Sobolev <- sph_stat_Sobolev(X = data, Psi_in_X = Psi_in_X, p = p,
                                        vk2 = Sobolev_vk2)

    }
    if (run_test$Stein) {

      Stein_vk2 <- weights_dfs_Sobolev(p = p, K_max = Stein_K, thre = 0,
                                       type = "Stein", Stein_cf = Stein_cf,
                                       verbose = FALSE,)$weights
      stats$Stein <- sph_stat_Sobolev(X = data, Psi_in_X = Psi_in_X, p = p,
                                      vk2 = Stein_vk2)

    }

  }

  # Avoid returning matrices in variables if there are vectorized tests.
  # Instead, return a data frame with Sobolev.1, Sobolev.2, etc. variables
  n_param_vectorized <- sapply(param_vectorized, function(par) {

    obj <- get(x = par)
    return(ifelse(is.null(dim(obj)), length(obj), nrow(obj)))

  })
  if (any(stats_vectorized %in% type) && any(n_param_vectorized > 1)) {

    stats <- do.call(data.frame, stats)

  }

  return(stats)

}
