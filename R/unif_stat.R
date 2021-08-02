

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
#' accepted are an array of size \code{c(n, 1, M)} or a vector of length
#' \code{n} with angular data. Must not contain \code{NA}'s.
#' @inheritParams unif_test
#' @param data_sorted is the circular data sorted? If \code{TRUE}, certain
#' statistics are faster to compute. Defaults to \code{FALSE}.
#' @param Rayleigh_m integer \eqn{m} for the \eqn{m}-modal Rayleigh test.
#' Defaults to \code{m = 1} (the standard Rayleigh test).
#' @param coverage_a \eqn{a_n = a / n} parameter used in the length of the arcs
#' of the coverage-based tests. Must be positive. Defaults to \code{2 * pi}.
#' @param Rothman_t \eqn{t} parameter for the Rothman test, a real in
#' \eqn{(0, 1)}. Defaults to \code{1 / 3}.
#' @param Cressie_t \eqn{t} parameter for the Cressie test, a real in
#' \eqn{(0, 1)}. Defaults to \code{1 / 3}.
#' @param Pycke_q \eqn{q} parameter for the Pycke "\eqn{q}-test", a real in
#' \eqn{(0, 1)}. Defaults to \code{1 / 2}.
#' @param Riesz_s \eqn{s} parameter for the \eqn{s}-Riesz test, a real in
#' \eqn{(0, 2)}. Defaults to \code{1}.
#' @param Cuesta_Albertos_rand_dirs a matrix of size \code{c(n_proj, p)}
#' containing \code{n_proj} random directions (in Cartesian coordinates) on
#' \eqn{S^{p-1}} to perform the Cuesta-Albertos test. If \code{NULL} (default),
#' a sample of size \code{n_proj = 50} directions is computed internally.
#' @param K_Cuesta_Albertos integer giving the truncation of the series
#' present in the asymptotic distribution of the Kolmogorov-Smirnov statistic.
#' Defaults to \code{5e2}.
#' @param Cai_regime type of asymptotic regime for Cai test, either \code{1}
#' (sub-exponential regime), \code{2} (exponential), or \code{3}
#' (super-exponential; default).
#' @return A data frame of size \code{c(M, length(type))}, with column names
#' given by \code{type}, that contains the values of the test statistics.
#' @details
#' Detailed descriptions and references of the statistics are available
#' in García-Portugués and Verdebout (2018).
#' @references
#' García-Portugués, E. and Verdebout, T. (2018) An overview of uniformity
#' tests on the hypersphere. \emph{arXiv:1804.00286}.
#' \url{https://arxiv.org/abs/1804.00286}.
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
#' # Cuesta-Albertos
#' unif_stat(data = X, type = "Cuesta_Albertos",
#'           Cuesta_Albertos_rand_dirs = X[, , 1])
#' unif_stat(data = X, type = "Cuesta_Albertos",
#'           Cuesta_Albertos_rand_dirs = X[, , 1], K_Cuesta_Albertos = 1)
#'
#' # Cai
#' unif_stat(data = X, type = "Cai", Cai_regime = 3)
#' unif_stat(data = X, type = "Cai", Cai_regime = 1)
#' @export
unif_stat <- function(data, type = "all", data_sorted = FALSE,
                      Rayleigh_m = 1, coverage_a = 2 * pi, Rothman_t = 1 / 3,
                      Cressie_t = 1 / 3, Pycke_q = 0.5, Riesz_s = 1, 
                      Cuesta_Albertos_rand_dirs = NULL, K_Cuesta_Albertos = 25,
                      Cai_regime = 3) {

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
      stats_type <- match.arg(arg = type, choices = avail_stats,
                              several.ok = TRUE)

    }

  } else if (is.numeric(type)) {

    type <- unique(type)
    stats_type <- avail_stats[type]

  } else {

    stop("type must be a character or a numeric vector.")

  }

  # Create a data frame for the statistics
  n_stats <- length(stats_type)
  stats <- matrix(0, nrow = M, ncol = n_stats)
  colnames(stats) <- stats_type
  stats <- as.data.frame(stats)

  # Compute efficiently several statistics
  if (p == 2) {

    # Statistics using gaps and sorted data
    stats_using_gaps <- c("Gini", "Gini_squared", "Greenwood", "Log_gaps",
                          "Max_uncover", "Num_uncover", "Range", "Rao",
                          "Vacancy")
    stats_using_sorted_data <- c("Cressie", "Feltz_Goldin", "Hodges_Ajne",
                                 "Kuiper", "Watson", "Watson_1976",
                                 "KS", "CvM", "AD", stats_using_gaps)

    # Statistics using the shortest angles matrix Psi
    stats_using_Psi <- c("Ajne", "Bakshaev", "Gine_Fn", "Gine_Gn",
                         "Hermans_Rasson", "PAD", "PCvM", "PRt", "Pycke",
                         "Pycke_q", "Rothman", "Riesz")

    # Evaluate which statistics to apply
    run_test <- as.list(c(avail_cir_tests, "KS", "CvM", "AD") %in% stats_type)
    names(run_test) <- c(avail_cir_tests, "KS", "CvM", "AD")

    ## Sorted-data statistics

    # Is it worth to sort Theta? Only if there is more than one statistic
    # using sorted Theta and if it has not already been done
    # CAUTION: replacement of data with sorted_data!
    n_stats_type_sorted <- sum(stats_type %in% stats_using_sorted_data)
    if (!data_sorted & (n_stats_type_sorted > 1)) {

      data <- sort_each_col(data)
      data_sorted <- TRUE

    }
    if (run_test$Cressie) {

      stats$Cressie <- cir_stat_Cressie(Theta = data, t = Cressie_t,
                                        sorted = data_sorted)

    }
    if (run_test$Feltz_Goldin) {

      stats$Feltz_Goldin <- cir_stat_Feltz_Goldin(Theta = data,
                                                  sorted = data_sorted)

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

      stats$AD <- cir_stat_PAD(Theta = data, sorted = data_sorted,
                               AD = TRUE)

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

      stats$Max_uncover <- cir_stat_Max_uncover(Theta = gaps, a = coverage_a,
                                                gaps_in_Theta = gaps_in_Theta,
                                                sorted = data_sorted)

    }
    if (run_test$Num_uncover) {

      stats$Num_uncover <- cir_stat_Num_uncover(Theta = gaps, a = coverage_a,
                                                gaps_in_Theta = gaps_in_Theta,
                                                sorted = data_sorted)

    }
    if (run_test$Vacancy) {

      stats$Vacancy <- cir_stat_Vacancy(Theta = gaps, a = coverage_a,
                                        gaps_in_Theta = gaps_in_Theta,
                                        sorted = data_sorted)

    }

    ## Other statistics

    if (run_test$Rayleigh) {

      stats$Rayleigh <- cir_stat_Rayleigh(Theta = data, m = Rayleigh_m)

    }
    if (run_test$Bingham) {

      stats$Bingham <- cir_stat_Bingham(Theta = data)

    }
    if (run_test$Cuesta_Albertos) {

      # Sample random directions
      if (is.null(Cuesta_Albertos_rand_dirs)) {

        Cuesta_Albertos_rand_dirs <- r_unif_sph(n = 50, p = 2, M = 1)[, , 1]

      }
      stats$Cuesta_Albertos <-
        cir_stat_Cuesta_Albertos(Theta = data,
                                 rand_dirs = Cuesta_Albertos_rand_dirs,
                                 K_Cuesta_Albertos = K_Cuesta_Albertos,
                                 original = FALSE)

    }

    ## Psi-based statistics

    # Is it worth to compute Psi? Only if there is more than one statistic
    # using it (after discounting equivalent tests if they are simultaneously
    # present) AND if 0.5 * n * (n - 1) * M is not too large.
    # CAUTION: replacement of data with Psi!
    n_stats_type_Psi <- sum(stats_type %in% stats_using_Psi)
    n_stats_type_Psi <- n_stats_type_Psi - ((run_test$Watson && run_test$PCvM) +
         (run_test$Rothman && run_test$PRt) + (Riesz_s == 2) +
         (Riesz_s == 0 && run_test$Pycke && run_test$Riesz) + 
         (Riesz_s == 1 && run_test$Bakshaev && run_test$Riesz) +
         (run_test$Gine_Fn && run_test$Gine_Gn && run_test$Ajne))
    if (n_stats_type_Psi > 1 & 0.5 * n * (n - 1) * M <= 1e8) {

      dim(data) <- c(n, 1, M)
      data <- Psi_mat(data = data)
      Psi_in_Theta <- TRUE

    } else {

      Psi_in_Theta <- FALSE

    }
    if (run_test$Ajne) {

      stats$Ajne <- cir_stat_Ajne(Theta = data, Psi_in_Theta = Psi_in_Theta)

    }
    if (run_test$Bakshaev) {

      stats$Bakshaev <- cir_stat_Riesz(Theta = data,
                                       Psi_in_Theta = Psi_in_Theta, s = 1)

    }
    if (run_test$Riesz) {

      if (Riesz_s == 1 && run_test$Bakshaev) {

        stats$Riesz <- stats$Bakshaev

      } else if (Riesz_s == 2 && (run_test$Rayleigh | !Psi_in_Theta)) {

        if (run_test$Rayleigh) {

          stats$Riesz <- stats$Rayleigh

        } else {

          stats$Riesz <- cir_stat_Rayleigh(Theta = data, m = 1)

        }

      } else {

        stats$Riesz <- cir_stat_Riesz(Theta = data, Psi_in_Theta = Psi_in_Theta,
                                      s = Riesz_s)

      }

    }
    if (run_test$Gine_Gn) {

      stats$Gine_Gn <- cir_stat_Gine_Gn(Theta = data,
                                        Psi_in_Theta = Psi_in_Theta)

    }
    if (run_test$Gine_Fn) {

      if (run_test$Ajne & run_test$Gine_Gn) {

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

      stats$Rothman <- cir_stat_Rothman(Theta = data, t = Rothman_t,
                                        Psi_in_Theta = Psi_in_Theta)

    }
    if (run_test$PRt) {

      if (run_test$Rothman) {

        stats$PRt <- stats$Rothman

      } else {

        stats$PRt <- cir_stat_PRt(Theta = data, Psi_in_Theta = Psi_in_Theta,
                                  t = Rothman_t)

      }

    }
    if (run_test$Pycke) {

      if (Riesz_s == 0 && run_test$Riesz) {

        stats$Pycke <- (2 * n) / (n - 1) * stats$Riesz

      } else {

        stats$Pycke <- cir_stat_Pycke(Theta = data, Psi_in_Theta = Psi_in_Theta)

      }

    }
    if (run_test$Pycke_q) {

      stats$Pycke_q <- cir_stat_Pycke_q(Theta = data,
                                        Psi_in_Theta = Psi_in_Theta,
                                        q = Pycke_q)

    }

  } else {

    # Statistics using the shortest angles matrix Psi
    stats_using_Psi <- c("Ajne", "Bakshaev", "Cai", "Gine_Fn", "Gine_Gn",
                         "PAD", "PCvM", "PRt", "Pycke", "Riesz")

    # Evaluate which statistics to apply
    run_test <- as.list(avail_sph_tests %in% stats_type)
    names(run_test) <- avail_sph_tests

    ## Other statistics

    if (run_test$Bingham) {

      stats$Bingham <- sph_stat_Bingham(X = data)

    }
    if (run_test$Cuesta_Albertos) {

      # Sample random directions
      if (is.null(Cuesta_Albertos_rand_dirs)) {

        Cuesta_Albertos_rand_dirs <- r_unif_sph(n = 50, p = p, M = 1)[, , 1]

      }
      stats$Cuesta_Albertos <-
        sph_stat_Cuesta_Albertos(X = data,
                                 rand_dirs = Cuesta_Albertos_rand_dirs,
                                 K_Cuesta_Albertos = K_Cuesta_Albertos,
                                 original = FALSE)

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
    n_stats_type_Psi <- sum(stats_type %in% stats_using_Psi)
    n_stats_type_Psi <- n_stats_type_Psi - ((Riesz_s == 2) +
         (Riesz_s == 0 && run_test$Pycke && run_test$Riesz) +
         (Riesz_s == 1 && run_test$Bakshaev && run_test$Riesz) +
         (run_test$PCvM && run_test$Bakshaev && p == 3) +
         (run_test$Gine_Fn && run_test$Gine_Gn && run_test$Ajne))
    if (n_stats_type_Psi > 1 & 0.5 * n * (n - 1) * M <= 1e8) {

      data <- Psi_mat(data = data)
      dim(data) <- c(dim(data), 1)
      Psi_in_X <- TRUE

    } else {

      Psi_in_X <- FALSE

    }
    if (run_test$Ajne) {

      stats$Ajne <- sph_stat_Ajne(X = data, Psi_in_X = Psi_in_X)

    }
    if (run_test$Bakshaev) {

      stats$Bakshaev <- sph_stat_Riesz(X = data, Psi_in_X = Psi_in_X, p = p,
                                       s = 1)

    }
    if (run_test$Riesz) {

      if (Riesz_s == 1 && run_test$Bakshaev) {

        stats$Riesz <- stats$Bakshaev

      } else if (Riesz_s == 2 && (run_test$Rayleigh | !Psi_in_X)) {

        if (run_test$Rayleigh) {

          stats$Riesz <- (2 / p) * stats$Rayleigh

        } else {

          stats$Riesz <- (2 / p) * sph_stat_Rayleigh(X = data) 

        }

      } else {

        stats$Riesz <- sph_stat_Riesz(X = data, Psi_in_X = Psi_in_X, p = p,
                                      s = Riesz_s)

      }

    }
    if (run_test$Cai) {

      if (Psi_in_X) {

        stats$Cai <- sph_stat_Cai(X = cos(data), Psi_in_X = TRUE, p = p,
                                  regime = Cai_regime)

      } else {

        stats$Cai <- sph_stat_Cai(X = data, Psi_in_X = FALSE, p = p,
                                  regime = Cai_regime)

      }

    }
    if (run_test$Gine_Gn) {

      stats$Gine_Gn <- sph_stat_Gine_Gn(X = data, Psi_in_X = Psi_in_X, p = p)

    }
    if (run_test$Gine_Fn) {

      if (run_test$Ajne & run_test$Gine_Gn) {

        stats$Gine_Fn <- 4 * stats$Ajne + stats$Gine_Gn

      } else {

        stats$Gine_Fn <- sph_stat_Gine_Fn(X = data, Psi_in_X = Psi_in_X, p = p)

      }

    }
    if (run_test$PAD) {

      stats$PAD <- sph_stat_PAD(X = data, Psi_in_X = Psi_in_X, p = p)

    }
    if (run_test$PCvM) {

      if (run_test$Bakshaev & p == 3) {

        stats$PCvM <- 0.125 * stats$Bakshaev

      } else {

        stats$PCvM <- sph_stat_PCvM(X = data, Psi_in_X = Psi_in_X, p = p)

      }

    }
    if (run_test$PRt) {

      stats$PRt <- sph_stat_PRt(X = data, t = Rothman_t, Psi_in_X = Psi_in_X,
                                p = p)

    }
    if (run_test$Pycke) {

      if (Riesz_s == 0 && run_test$Riesz) {

        if (p == 3) {

          stats$Pycke <- n / (2 * pi * (n - 1)) * 
            (stats$Riesz - (log(4) - 1) / 2)

        } else {

          warning("Pycke statistic is only defined for p = 2,3. Using Riesz statistic with s = 0 instead, which behaves consistently across dimensions.")
          stats$Pycke <- stats$Riesz

        }

      } else {

        if (Psi_in_X) {

          stats$Pycke <- sph_stat_Pycke(X = cos(data), Psi_in_X = TRUE, p = p)

        } else {

          stats$Pycke <- sph_stat_Pycke(X = data, Psi_in_X = FALSE, p = p)

        }

      }

    }

  }

  return(stats)

}
