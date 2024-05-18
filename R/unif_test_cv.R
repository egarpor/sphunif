#' @title Circular and (hyper)spherical uniformity K-fold tests
#'
#' @description Implementation of parameter-dependent uniformity tests
#' on the (hyper)sphere \eqn{S^{p-1}:=\{{\bf x}\in R^p:||{\bf x}||=1\}}{
#' S^{p-1}:=\{x\in R^p:||x||=1\}}, \eqn{p\ge 2} following the \eqn{K}-fold
#' cross-validation procedure, with asymptotically-exact Harmonic Mean P-value
#' calibration either in terms of their asymptotic distributions, if available,
#' or Monte Carlo.
#'
#' \code{unif_test_cv} receives a sample of directions
#' \eqn{{\bf X}_1,\ldots,{\bf X}_n\in S^{p-1}}{X_1,\ldots,X_n\in S^{p-1}} in
#' \emph{Cartesian coordinates}, except for the circular case (\eqn{p=2}) in
#' which the sample can be represented in terms of \emph{angles}
#' \eqn{\Theta_1,\ldots,\Theta_n\in [0, 2\pi)}.
#'
#' \code{unif_test_cv} allows to perform several tests within a single call,
#' facilitating thus the exploration of a dataset by applying several tests.
#'
#' \code{unif_test_cv} needs a grid of parameters to find the one that maximizes
#' the power proxy. The grids are specified for each statistic parameter.
#'
#' \code{null_var} computes the exact-\code{n} variance of the statistic for
#' a set of parameters \code{lambda_grid}.
#'
#' @inheritParams unif_test
#' @inheritParams unif_stat_distr
#' @param type type of test to be applied. A character vector containing any of
#' the following types of tests, depending on the dimension \eqn{p}:
#' \itemize{
#'   \item Circular data: any of the names available at object
#'   \code{\link{avail_cir_cv_tests}}.
#'   \item (Hyper)spherical data: any of the names available at object
#'   \code{\link{avail_sph_cv_tests}}.
#' }
#' If \code{type = "all"} (default), then \code{type} is set as
#' \code{avail_cir_cv_tests} or \code{avail_sph_cv_tests}, depending on the
#' value of \eqn{p}.
#' @param K Number of folds of (roughly) equal sizes to split \code{data}.
#' @param p_value type of \eqn{p}-value computation. Either \code{"MC"} for
#' employing the approximation by Monte Carlo of the exact null distribution or
#' \code{"asymp"} (default) for the use of the asymptotic null distribution
#'  (if available).
#' @param null_variance TODO
#' @param rel.tol TODO
#' @param seed_fold an integer that fixes the seed for splitting data into
#' \code{K} folds.
#' @param lambda_grid vector with parameters to compute null variance of the
#' statistic
#' @param ... If \code{p_value = "MC"}, optional performance parameters to
#' be passed to \code{\link{unif_stat_MC}}: \code{chunks},
#' \code{cores}, and \code{seed}.
#'
#' @return If only a \bold{single test} is performed, a list with class
#' \code{htest} containing the following components:
#' \itemize{
#'   \item \code{fold_statistics}: the value of the test statistic for
#'   each fold.
#'   \item \code{fold_params}: the value of the optimal parameter for
#'   each fold.
#'   \item \code{fold_p.values}: the p-values of the test for each fold.
#'   \item \code{p.value}: the HMP-aggregated p-value of the test.
#'   \item \code{alternative}: a character string describing the alternative
#'   hypothesis.
#'   \item \code{method}: a character string indicating what type of test was
#'   performed.
#'   \item \code{data.name}: a character string giving the name of the data.
#'   \item \code{reject}: the rejection decision for the levels of significance
#'   \code{alpha}.
#   \item \code{crit_val}: a vector with the critical values for the
#   significance levels \code{alpha} used with \code{p_value = "MC"} or
#   \code{p_value = "asymp"}.
#' }
#' If \bold{several tests} are performed, a \code{type}-named list with
#' entries for each test given by the above list.
#'
#' \code{null_var} returns a vector of length \code{length(lambda_grid)}
#' that contains the values of the exact-n variance of the statistic
#' \code{type} under the null hypothesis.
#' @details
#' All the tests reject for large values of the test statistic, so the critical
#' values for the significance levels \code{alpha} correspond to the
#' \code{alpha}-upper quantiles of the null distribution of the test statistic.
#'
#  TODO:
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
#' Description of the \eqn{K}-fold cross-validation procedure is available
#' in Fernández-de-Marcos and García-Portugués (2023).
#' Descriptions and references for most of the tests are available
#' in García-Portugués and Verdebout (2018).
#' @references
#' Fernández-de-Marcos, A. and García-Portugués, E. (2023) On new omnibus tests
#' of uniformity on the hypersphere. \emph{Test}, 32(4):1508-–1529.
#' \doi{10.1007/s11749-023-00882-x}
#'
#' García-Portugués, E. and Verdebout, T. (2018) An overview of uniformity
#' tests on the hypersphere. \emph{arXiv:1804.00286}.
#' \url{https://arxiv.org/abs/1804.00286}.
#' @examples
#' ## Asymptotic distribution
#'
#' seed <- 12345
#'
#' # Circular data
#' n <- 50
#' samp_cir <- r_unif_cir(n = n)
#'
#' # Matrix
#' unif_test_cv(data = samp_cir, type = "all", K = 3, p_value = "asymp",
#'              seed_fold = seed)
#'
#' # Vector
#' unif_test_cv(data = samp_cir[, 1], type = "all", K = 3, p_value = "asymp",
#'              seed_fold = seed)
#'
#' # Array
#' unif_test_cv(data = array(samp_cir, dim = c(n, 1, 1)), type = "all", K = 3,
#'              p_value = "asymp", seed_fold = seed)
#'
#' # Spherical data
#' n <- 50
#' samp_sph <- r_unif_sph(n = n, p = 3)
#'
#' # Array
#' unif_test_cv(data = samp_sph, type = c("Poisson", "Softmax"), K = 3,
#'              p_value = "asymp", seed_fold = seed)
#'
#' # Matrix
#' unif_test_cv(data = samp_sph[, , 1], type = c("Poisson", "Softmax"), K = 3,
#'              p_value = "asymp", seed_fold = seed)
#'
#' ## Monte Carlo
#'
#' # Circular data
#' unif_test_cv(data = samp_cir, type = "all", K = 3, p_value = "MC", M = 1e3,
#'              seed_fold = seed)
#'
#' # Spherical data
#' unif_test_cv(data = samp_sph, type = c("Poisson", "Softmax"), K = 3, M = 1e3,
#'              p_value = "MC", seed_fold = seed)
#'
#' # Caching stats_MC
#' stats_MC_cir <- unif_stat_MC(n = nrow(samp_cir), type = avail_cir_cv_tests,
#'                              p = 2, M = 1e3, r_H1 = NULL, crit_val = NULL,
#'                              return_stats = TRUE, stats_sorted = TRUE,
#'                              Poisson_rho = seq(0.1, 0.9, 0.1),
#'                              Softmax_kappa = seq(0.1, 20, 1),
#'                              Stereo_a = seq(-1, 1, 0.25))$stats_MC
#' stats_MC_sph <- unif_stat_MC(n = nrow(samp_sph), type = avail_sph_cv_tests,
#'                              p = 3, M = 1e3, r_H1 = NULL, crit_val = NULL,
#'                              return_stats = TRUE, stats_sorted = TRUE,
#'                              Poisson_rho = seq(0.1, 0.9, 0.1),
#'                              Softmax_kappa = seq(0.1, 20, 1),
#'                              Stereo_a = seq(-1, 1, 0.25))$stats_MC
#' unif_test_cv(data = samp_cir, type = avail_cir_tests, K = 3, p_value = "MC",
#'              stats_MC = stats_MC_cir, seed_fold = seed)
#' unif_test_cv(data = samp_sph, type = c("Poisson", "Softmax"), K = 3,
#'              p_value = "MC", stats_MC = stats_MC_sph, seed_fold = seed)
#' @name unif_test_cv

#' @rdname unif_test_cv
#' @export
unif_test_cv <- function(data, type = "all", K = 10, p_value = "asymp",
                         alpha = c(0.10, 0.05, 0.01), M = 1e4, stats_MC = NULL,
                         K_max = 1e4, method = "I", null_variance = NULL,
                         rel.tol = 1e-10, Poisson_rho = seq(0.1, 0.9, 0.1),
                         Softmax_kappa = seq(0.1, 20, 1),
                         Stereo_a = seq(-1, 1, 0.25),
                         seed_fold = NULL, ...) {

  # Check K > 1
  if (K < 2) {

    stop(paste0("The number of folds (K = ", K, ") must be at least 2."))

  }

  # Read data's name
  data_name <- deparse(substitute(data))

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

    avail_stats <- avail_cir_cv_tests
    p <- 2

    # As polar coordinates
    if (d == 2) {

      dim(data) <- c(n, d, 1)
      data <- X_to_Theta(X = data)

    }

  } else {

    avail_stats <- avail_sph_cv_tests
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

      # TODO: Improvement. Warn the user that some statistics are being omitted.
      stats_type <- try(match.arg(arg = type, choices = avail_stats,
                                  several.ok = TRUE), silent = TRUE)
      if (inherits(stats_type, "try-error")) {

        stop(
          paste(
            strwrap(
              paste0("CV test with type = \"", type, "\" is unsupported. ",
                     "Must be one of the following tests: \"",
                     paste(avail_stats, collapse = "\", \""), "\"."),
              width = 80, indent = 0, exdent = 2),
            collapse = "\n")
          )

      }

    }

  } else {

    stop("type must be a character vector")

  }

  # Check if there is any missing statistic in stats_MC
  if (!is.null(stats_MC)) {

    # Dummy stats
    check_stat <- unif_stat(data = r_unif_sph(n = 2, p = p, M = 1),
                            type = stats_type, Poisson_rho = Poisson_rho,
                            Softmax_kappa = Softmax_kappa, Stereo_a = Stereo_a)

    # Names check
    checks <- names(check_stat) %in% colnames(stats_MC)
    if (any(!checks)) {

      stop(paste("stats_MC must be a data.frame with colnames containing",
                 "the tests names returned by unif_stat(...). stats_MC misses",
                 paste(paste0("\"", colnames(check_stat)[!checks], "\""),
                       collapse = ", "), ". Check grids of parameters."))

    }

  }

  # Number of statistics
  n_stats <- length(stats_type)

  # Get null variance of statistic. If null_variance is given, check it is
  # the same size of lambda_grid. Otherwise, compute asymptotic null variance.
  if (is.null(null_variance)) {

    null_variance <- sapply(stats_type, function(stat_type) {

      lambda_grid <- switch(stat_type,
                            "Poisson" = Poisson_rho,
                            "Softmax" = Softmax_kappa,
                            "Stereo" = Stereo_a)

      return(null_var(n = round(n / K), p = p, type = stat_type,
                      lambda_grid = lambda_grid,
                      rel.tol = rel.tol))

    })


  } else {

    # TODO: Include null variance as an input, whether it has been done
    # asymptotically or exactly.

  }

  # Define parameter names for unif_stat_ calls
  param_name <- list("Poisson" = "Poisson_rho",
                     "Softmax" = "Softmax_kappa",
                     "Stereo" = "Stereo_a")
  param_args_name <- unname(sapply(stats_type, function(stat_type) {

    param_name[[stat_type]]

  }))

  # Split data into K disjoint subsamples of (roughly) same sizes
  folds <- k_fold_split(n = n, K = K, seed = seed_fold)

  # Create summary data.frames for statistic, p.values and optimal parameters.
  p_val <- vector("list", length = n_stats)
  names(p_val) <- stats_type
  p_val <- as.data.frame(p_val)

  lambda_hat <- vector("list", length = n_stats)
  names(lambda_hat) <- stats_type
  lambda_hat <- as.data.frame(lambda_hat)

  stat <- vector("list", length = n_stats)
  names(stat) <- stats_type
  stat <- as.data.frame(stat)

  # K-fold testing
  for (k in 1:K){

    # Compute estimator of approximate oracle parameter \hat{lambda}(S_k)
    if (p == 2) {
      Sk <- data[folds[[k]], ]
    } else {
      Sk <- data[folds[[k]], , ]
      dim(Sk) <- c(dim(Sk), 1)
    }

    # Compute power-approximate score in grid
    stat_k <- unif_stat(Sk, type = stats_type,
                        Poisson_rho = Poisson_rho,
                        Softmax_kappa = Softmax_kappa,
                        Stereo_a = Stereo_a)

    lambda_hat_k <- sapply(stats_type, function(stat_type) {

      lambda_grid <- switch(stat_type,
                            "Poisson" = Poisson_rho,
                            "Softmax" = Softmax_kappa,
                            "Stereo" = Stereo_a)

      stat_k_cols <- names(stat_k)
      idx_stat_type <- sapply(strsplit(stat_k_cols, "\\."),
                              function(x) x[1] == stat_type)
      stat_type_k <- stat_k_cols[idx_stat_type]
      # TODO: Generalize to other statistics: Score in case of V-statistic must
      # be E_H1 - H_H0
      q <- stat_k[stat_type_k] / sqrt(null_variance[[stat_type]])

      return(lambda_grid[which.max(q)])

    })

    lambda_hat <- rbind(lambda_hat, lambda_hat_k)

    # Perform test based on T(\hat{\lambda}; S\S_k) on the remaining subsamples
    if (p == 2) {
      S_notk <- data[!((1:n) %in% folds[[k]]), ]
    } else {
      S_notk <- data[!((1:n) %in% folds[[k]]), , ]
      dim(S_notk) <- c(dim(S_notk), 1)
    }

    unif_stat_args <- list(data = S_notk, type = stats_type)
    specific_args <- sapply(stats_type, function(stat_type) {

      return(lambda_hat_k[[stat_type]])

    })

    # Statistic on K - 1 folds
    names(specific_args) <- param_args_name
    stat_k <- do.call(what = unif_stat, args = c(unif_stat_args, specific_args))

    stat <- rbind(stat, stat_k)

    # Calibration
    if (p_value == "MC") {

      # Get the stats_MC
      if (is.null(stats_MC)) {

        unif_stat_MC_args <- list(n = n, type = stats_type, p = p, M = M,
                                  r_H1 = NULL, crit_val = NULL, alpha = alpha,
                                  return_stats = TRUE, stats_sorted = TRUE)
        stats_MC_k <- do.call(what = unif_stat_MC,
                              args = c(unif_stat_MC_args,
                                       specific_args))$stats_MC

        # p-values
        p_val_k <- 1 - as.data.frame(sapply(stats_type, function(distr) {
          ecdf_bin(data = stats_MC_k[[distr]], sorted_x = stat_k[[distr]],
                   data_sorted = TRUE,
                   efic = TRUE, divide_n = TRUE)
        }, simplify = FALSE))

      } else {

        # Get the index of distr.idx for lambda_hat_k
        idx_lambda_hat <- sapply(stats_type, function(stat_type) {

          lambda_grid <- switch(stat_type,
                                "Poisson" = Poisson_rho,
                                "Softmax" = Softmax_kappa,
                                "Stereo" = Stereo_a)

          return(which(lambda_grid == lambda_hat_k[[stat_type]]))

        })

        # p-values
        p_val_k <- 1 - as.data.frame(sapply(stats_type, function(distr) {
          ecdf_bin(data = stats_MC[[paste(distr, idx_lambda_hat[[distr]],
                                          sep = ".")]],
                   sorted_x = stat_k[[distr]], data_sorted = TRUE,
                   efic = TRUE, divide_n = TRUE)
        }, simplify = FALSE))

      }

      p_val <- rbind(p_val, p_val_k)

    } else if (p_value == "asymp") {

      # TODO: Include verbose
      unif_stat_distr_args <- list(x = stat_k, type = stats_type, p = p,
                                   n = n, approx = "asymp", stats_MC = NULL,
                                   M = M, K_max = K_max, method = method)
      p_val_k <- 1 - do.call(what = unif_stat_distr,
                             args = c(unif_stat_distr_args, specific_args))

      p_val <- rbind(p_val, p_val_k)

    } else {

      stop("Wrong choice for calibration, must be \"MC\" or \"asymp\".")

    }

  }

  # Set a minimum value for p_values when they are 0 to deal with HMP.
  # TODO: When p_value = "asymp", seems a little odd.
  p_val <- apply(p_val, 2, function(p) pmax(p, 1 / M))

  # Compute asymptotically-exact HMP
  p_val_hmp <- apply(p_val, 2, function(p) {

    harmonicmeanp::p.hmp(p, L = length(p))

  })

  # Rejection?
  reject <- rbind(sapply(alpha, function(a) p_val_hmp < a))
  colnames(reject) <- alpha

  ## Prepare return object

  # Create list of htest objects
  test <- vector(mode = "list", length = n_stats)
  names(test) <- stats_type

  for (i in seq_along(stats_type)) {

    # Type of test
    if (p == 2) {

      method <- switch(stats_type[i],
                       "Poisson" = "Poisson-kernel test of circular uniformity (CV)",
                       "Softmax" = "Softmax test of circular uniformity (CV)"
      )

      alternative <- switch(stats_type[i],
                            "Poisson" = "any alternative to circular uniformity for rho > 0",
                            "Softmax" = "any alternative to circular uniformity for kappa > 0"
      )

    } else {

      method <- switch(stats_type[i],
                       "Poisson" = "Poisson-kernel test of spherical uniformity (CV)",
                       "Softmax" = "Softmax test of spherical uniformit (CV)",
                       "Stereo" = "Stereographic projection test of spherical uniformity (CV)",
      )

      alternative <- switch(stats_type[i],
                            "Poisson" = "any alternative to spherical uniformity for rho > 0",
                            "Softmax" = "any alternative to spherical uniformity for kappa > 0",
                            "Stereo" = "any alternative to spherical uniformity for |a| < 1"
      )

    }

    # htest object
    test[[i]] <- list(fold_statistics = stat[, i],
                      fold_params = lambda_hat[, i],
                      fold_p.values = p_val[, i], p.value = p_val_hmp[i],
                      alternative = alternative,
                      method = method, data.name = data_name,
                      reject = reject[i, ])
    class(test[[i]]) <- "htest"

  }

  # If there is only one test, return htest directly
  if (n_stats == 1) {

    test <- test[[1]]

  }

  return(test)

}

k_fold_split <- function(n, K, seed = NULL) {

  # Divide into (roughly) same sizes K folds
  if (K > (n / 2)) {
    stop(paste0("Number of folds, K = ", K, ", must be lower or equal ",
                "to half of the sample size (n/2) = ", floor(n / 2), "."))
  }
  quo <- as.integer(n / K)
  rem <- n / K - quo

  # Random ordering setting seed
  if (!is.null(seed)) {

    set.seed(seed)

  }

  rand_order <- sample.int(n)

  # Get exact number of elements for each partition from random ordered index
  idx <- vector("list", length = K)
  for (k in 1:K){

    idx[[k]] <- rand_order[(((k - 1) * quo + as.integer(rem * k)) :
                              (k * quo + as.integer(rem * (k + 1)) - 1)) + 1]

  }

  return(idx)

}

#' @rdname unif_test_cv
#' @export
null_var <- function(n, p, type, lambda_grid, rel.tol = 1e-10) {
  # TODO: Expand the number of statistics available

  alpha <- 0.5 * p - 1

  if (type == "Poisson") {

    psi <- function(th, lambda) {

      exp(log(1 - lambda^2) - (p / 2) * log(1 - 2 * lambda * cos(th)
                                            + lambda^2)
          + (p - 1) * log(1 - lambda) - log(1 + lambda))

    }
    # TODO: Check if rel.tol affects anything else.
    # SOLVED: Throws an error when p = 4, because of the finiteness of
    # the integral.


    if (p == 2) {

      log_b0 <- log(1 - lambda_grid) - log(1 + lambda_grid)

    } else {

      log_b0 <- (p - 1) * log(1 - lambda_grid) - log(1 + lambda_grid)

    }

  } else if (type == "Softmax") {

    psi <- function(th, lambda) exp(lambda * (cos(th) - 1))

    if (p == 2) {

      log_b0 <- log(besselI(x = lambda_grid, nu = 0, expon.scaled = TRUE))

    } else {

      log_b0 <- alpha * log(2 / lambda_grid) + lgamma(alpha) +
        log(alpha) + log(besselI(x = lambda_grid,
                                 nu = alpha, expon.scaled = TRUE))

    }

  } else if (type == "Stereo") {

    if (p == 2) {

      stop(paste0("\'Stereo\' statistic is not defined when p = 2."))

    } else if (p == 3) {

      stop(paste0("Variance of \'Stereo\' statistic under uniformity ",
                  "(H_0) is not finite when p = 3."))

    }

    psi <- function(th, lambda) 1 / tan(th / 2) + lambda * tan(th / 2)

    log_b0 <- log(1 + lambda_grid) + log(alpha) + 2 * (
      lgamma(alpha) - lgamma((p - 1) / 2))

  } else {

    stop("Incompatible choice of p and type.")

  }

  b0 <- exp(log_b0)
  b0_sq <- numeric(length(lambda_grid))
  for (i in seq_along(lambda_grid)) {

    lambda <- lambda_grid[i]
    b0_sq[i] <- rotasym::w_p(p - 1) / rotasym::w_p(p) * integrate(
      function(x) psi(acos(x), lambda)^2 * (1 - x^2)^((p - 3) / 2),
      lower = -1, upper = 1, rel.tol = rel.tol)$value

  }

  return(2 * (n - 1) / n * (b0_sq - b0^2))

}

#' @rdname unif_test_cv
#' @export
avail_cir_cv_tests <- c("Poisson", "Softmax")

#' @rdname unif_test_cv
#' @export
avail_sph_cv_tests <- c("Poisson", "Softmax", "Stereo")
