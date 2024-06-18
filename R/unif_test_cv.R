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
#' the power proxy. The grids are specified for each test statistic parameter.
#'
#' @param type type of test to be applied. A character vector or a single string
#'  containing any of the following types of tests, depending on the dimension
#'  \eqn{p}:
#' \itemize{
#'   \item Circular data: any of the names available at object
#'   \code{\link{avail_cir_cv_tests}}.
#'   \item (Hyper)spherical data: any of the names available at object
#'   \code{\link{avail_sph_cv_tests}}.
#' }
#' In \code{unif_test_cv}, if \code{type = "all"} (default), then \code{type}
#' is set as \code{avail_cir_cv_tests} or \code{avail_sph_cv_tests}, depending
#' on the value of \eqn{p}.
#' @param K Number of folds of (roughly) equal sizes to split \code{data}.
#' @param folds list of length \code{K}, each element containing a vector with
#' the indexes of observations included in each fold, as returned by a call
#' to \code{k_fold_split}. Optional in \code{k_fold_stat}, defaults to
#' \code{NULL}. If \code{NULL}, folds are computed internally calling
#' \code{k_fold_split}, but \code{K} must be specified.
#' @param lambda_grid list of size \code{length(type)} and names equal to
#' \code{type} with a vector of parameters to compute the optimal
#' \code{lambda_hat} (see \code{Poisson_rho}, \code{Softmax_kappa},
#' or \code{Stereo_a} for additional requirements).
#' @param p_vals matrix of size \code{K x length(type)} containing p_values
#' computed for each \code{K} fold and test in \code{type}.
#' @param p_value type of \eqn{p}-value computation. Either \code{"MC"} for
#' employing the approximation by Monte Carlo of the exact null distribution or
#' \code{"asymp"} (default) for the use of the asymptotic null distribution
#'  (if available).
#' @param null_variance list of length \code{length(type)} and names equal
#' to \code{type} that contains the null variance values returned by
#' \code{null_var_Sobolev()} for the required grid of parameters. If \code{NULL}
#' (default), it is computed internally.
#' @param seed_fold an integer that fixes the seed for splitting data into
#' \code{K} folds. It also applies to \code{seed} in \code{k_fold_split()} and
#' \code{k_fold_stat()}.
#' @param ... If \code{p_value = "MC"}, optional performance parameters to
#' be passed to \code{\link{unif_stat_MC}}: \code{chunks},
#' \code{cores}, and \code{seed}.
#' @inheritParams Sobolev
#' @inheritParams unif_stat_distr
#' @inheritParams unif_test
#'
#' @return \code{unif_test_cv}: If only a \bold{single test} is performed, a
#' list with class \code{htest} containing the following components:
#' \itemize{
#'   \item \code{fold_statistics}: the value of the test statistic for
#'   each fold, a vector of length \code{K}.
#'   \item \code{fold_params}: the value of the optimal parameter for
#'   each fold, a vector of length \code{K}.
#'   \item \code{fold_p.values}: the p-values of the test for each fold, a
#'   vector of length \code{K}.
#'   \item \code{p.value}: the HMP-aggregated p-value of the test.
#'   \item \code{alternative}: a character string describing the alternative
#'   hypothesis.
#'   \item \code{method}: a character string indicating what type of test was
#'   performed.
#'   \item \code{data.name}: a character string giving the name of the data.
#'   \item \code{reject}: the rejection decision for the levels of significance
#'   \code{alpha}.
#' }
#' If \bold{several tests} are performed, a \code{type}-named list with
#' entries for each test given by the above list.
#'
#' Intermediate functions return:
#' \itemize{
#'   \item \code{lambda_hat}: a matrix of size \code{K x length(type)}
#'   containing the optimal parameters for each of the \code{K} folds.
#'   \item \code{k_fold_stat}: a named list with \code{statistic} and
#'   \code{lambda_hat} with a data.frame and a matrix, respectively,
#'   of size \code{K x length(type)} containing the statistic and the optimal
#'   parameters for each of the \code{K} folds.
#'   \item \code{p_val_hmp}: a matrix of size \code{K x length(type)}
#'   containing the aggregated Harmonic Mean P-value for each test in \code{type}.
#'   \item \code{k_fold_split}: a list of length \code{K}, each element
#'   containing a vector with the indexes of observations included in each fold.
#'   \item \code{prepare_test_data}: a named list containing the processed data
#'   for testing uniformity in \code{data}, the sample size in \code{n}, the
#'   ambient dimension in \code{p}, and a list of the available CV tests
#'   depending on \code{p} in \code{avail_stats}.
#' }
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
#' in Fernández-de-Marcos and García-Portugués (2023). Descriptions and
#' references for most of the tests are available in
#' García-Portugués andVerdebout (2018).
#' @references
#' Fernández-de-Marcos, A. and García-Portugués, E. (2023) On new omnibus tests
#' of uniformity on the hypersphere. \emph{Test}, 32(4):1508-–1529.
#' \doi{10.1007/s11749-023-00882-x}
#'
#' García-Portugués, E. and Verdebout, T. (2018) An overview of uniformity
#' tests on the hypersphere. \emph{arXiv:1804.00286}.
#' \doi{10.48550/arXiv.1804.00286}.
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
#'              seed_fold = seed, K_max = 1e3, verbose = FALSE)
#'
#' # Vector
#' unif_test_cv(data = samp_cir[, 1], type = "all", K = 3, p_value = "asymp",
#'              seed_fold = seed, K_max = 1e3, verbose = FALSE)
#'
#' # Array
#' unif_test_cv(data = array(samp_cir, dim = c(n, 1, 1)), type = "all", K = 3,
#'              p_value = "asymp", seed_fold = seed, K_max = 1e3,
#'              verbose = FALSE)
#'
#' # Spherical data
#' n <- 50
#' samp_sph <- r_unif_sph(n = n, p = 3)
#'
#' # Array
#' unif_test_cv(data = samp_sph, type = c("Poisson", "Softmax"), K = 3,
#'              p_value = "asymp", seed_fold = seed, K_max = 1e3,
#'              verbose = FALSE)
#'
#' # Matrix
#' unif_test_cv(data = samp_sph[, , 1], type = c("Poisson", "Softmax"), K = 3,
#'              p_value = "asymp", seed_fold = seed, K_max = 1e3,
#'              verbose = FALSE)
#'
#' ## Monte Carlo
#'
#' # Circular data
#' unif_test_cv(data = samp_cir, type = "all", K = 3, p_value = "MC", M = 1e3,
#'              seed_fold = seed, K_max = 1e3, verbose = FALSE)
#'
#' # Spherical data
#' unif_test_cv(data = samp_sph, type = c("Poisson", "Softmax"), K = 3, M = 1e3,
#'              p_value = "MC", seed_fold = seed, K_max = 1e3, verbose = FALSE)
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
#'              stats_MC = stats_MC_cir, seed_fold = seed, K_max = 1e3,
#'              verbose = FALSE)
#' unif_test_cv(data = samp_sph, type = c("Poisson", "Softmax"), K = 3,
#'              p_value = "MC", stats_MC = stats_MC_sph, seed_fold = seed,
#'              K_max = 1e3, verbose = FALSE)
#'
#' ## Pre-specifying null_variance
#' Poisson_grid <- c(0.1, 0.5, 0.7)
#' Softmax_grid <- c(0.1, 0.5, 1, 5, 10)
#' null_variance <- sapply(c("Poisson", "Softmax"), function(stat_type) {
#'
#'   lambda_grid <- switch(stat_type,
#'                         "Poisson" = Poisson_grid,
#'                         "Softmax" = Softmax_grid)
#'
#'   return(null_var_Sobolev(n = round(n / 3), p = 3, type = stat_type,
#'                   lambda_grid = lambda_grid, verbose = FALSE))
#'
#' })
#' unif_test_cv(data = samp_sph, type = c("Poisson", "Softmax"), K = 3,
#'              p_value = "MC", M = 1e3, null_variance = null_variance,
#'              seed_fold = seed, Poisson_rho = Poisson_grid,
#'              Softmax_kappa = Softmax_grid, K_max = 1e3, verbose = FALSE)
#'
#' ## Using a progress bar when p_value = "MC"
#'
#' # Define a progress bar
#' require(progress)
#' require(progressr)
#' handlers(handler_progress(
#'   format = ":spin [:bar] :percent Total: :elapsedfull End \u2248 :eta",
#'   clear = FALSE))
#'
#' # Call unif_test() within with_progress()
#' with_progress(
#'   unif_test_cv(data = samp_sph, type = c("Poisson", "Softmax"), K = 3,
#'                p_value = "MC", , M = 1e3, seed_fold = seed,
#'                K_max = 1e3, verbose = FALSE, chunks = 10)
#' )
#' @name unif_test_cv
#'
# TODO: Update simulations-stereo when function is uploaded as part of the
#       package
# TODO: Remove unif_test_cv from simulations

#' @rdname unif_test_cv
#' @export
lambda_hat <- function(data, type, lambda_grid, folds, K_max = 1e3,
                       verbose = TRUE) {

  # Get data parameters
  d <- prepare_test_data(data)
  data <- d[["data"]]
  n <- d[["n"]]
  p <- d[["p"]]
  K <- length(folds)

  # Compute estimator of E_{H_1}[T_n(S_k)]
  stat <- data.frame(t(sapply(seq_along(folds), function(k) {

    # Subsample S_k
    if (p == 2) {
      Sk <- data[folds[[k]], ]
    } else {
      Sk <- data[folds[[k]], , ]
      dim(Sk) <- c(dim(Sk), 1)
    }

    # Compute statistics on S_k with lambda_grid parameters
    unif_stat_args <- list(data = Sk, type = type)
    stat_args <- lapply(type, function(t) lambda_grid[[t]])
    names(stat_args) <- lapply(type, function(t) stat_args_name[[t]])
    stat_k <- do.call(what = unif_stat, args = c(unif_stat_args, stat_args))
    names_stat <- names(stat_k)

    stat_k <- as.double(stat_k)
    names(stat_k) <- names_stat

    return(stat_k)

  })))

  # Compute null variance
  null_variance <- lapply(type, function(t) {

    return(null_var_Sobolev(n = round(n / K), p = p, type = t,
                    lambda_grid = lambda_grid[[t]], K_max = K_max,
                    verbose = verbose))

  })
  names(null_variance) <- type

  stat_cols <- names(stat)
  lambda_optimal <- sapply(type, function(t) {

    idx_type <- sapply(strsplit(stat_cols, "\\."),
                            function(x) x[1] == t)
    type_k <- stat_cols[idx_type]

    # TODO: Generalize to other statistics: Score in case of V-statistic must
    # be E_H1 - H_H0
    q <- stat[type_k] / sqrt(null_variance[[t]])

    return(lambda_grid[[t]][apply(q, 1, which.max)])

  })
  lambda_optimal <- rbind(lambda_optimal)

  return(lambda_optimal)

}

#' @rdname unif_test_cv
#' @export
k_fold_stat <- function(data, type, lambda_grid,
                        folds = NULL, K = NULL, seed = NULL, K_max = 1e3,
                        verbose = TRUE) {

  # Get data parameters
  d <- prepare_test_data(data)
  data <- d[["data"]]
  n <- d[["n"]]
  p <- d[["p"]]

  # If folds are not given, compute K folds internally
  if (is.null(folds)) {

    if (is.null(K)) {

      stop("Number of folds, K, must be given in order to split the data.")

    } else {

      folds <- k_fold_split(n = n, K = K, seed = seed)

    }

  }

  # Compute power-approximate score in grid
  lambda_opt <- lambda_hat(data = data, type = type,
                           lambda_grid = lambda_grid, folds = folds,
                           K_max = K_max, verbose = verbose)

  # Perform test based on T(\hat{\lambda}; S\S_k) on the remaining subsamples
  stat <- sapply(seq_along(folds), function(k) {

    if (p == 2) {
      S_notk <- data[!((1:n) %in% folds[[k]]), ]
    } else {
      S_notk <- data[!((1:n) %in% folds[[k]]), , ]
      dim(S_notk) <- c(dim(S_notk), 1)
    }

    unif_stat_args <- list(data = S_notk, type = type)
    stat_args <- lapply(type, function(t) lambda_opt[k, t])
    names(stat_args) <- lapply(type, function(t) stat_args_name[[t]])

    # Statistic on K - 1 folds
    stat_k <- as.matrix(do.call(what = unif_stat,
                                args = c(unif_stat_args, stat_args)))

    return(stat_k)

  })

  stat <- t(rbind(stat))
  colnames(stat) <- type
  stat <- as.data.frame(stat)

  return(list("statistic" = stat, "lambda_hat" = lambda_opt))

}

#' @rdname unif_test_cv
#' @export
p_val_hmp <- function(p_vals, M) {

  # Set a minimum value for p_values when they are 0 to deal with HMP.
  p_val <- apply(p_vals, 2, function(p) pmax(p, 1 / M))

  # Compute asymptotically-exact HMP
  p_val <- apply(p_val, 2, function(p) {

    harmonicmeanp::p.hmp(p, L = length(p))

  })

  p_val <- rbind(p_val)

  return(p_val)

}

#' @rdname unif_test_cv
#' @export
unif_test_cv <- function(data, type, K, p_value = "asymp",
                         alpha = c(0.10, 0.05, 0.01), M = 1e4, stats_MC = NULL,
                         K_max = 1e4, method = "I", null_variance = NULL,
                         Poisson_rho = seq(0.1, 0.9, 0.1),
                         Softmax_kappa = c(0.1, seq(1, 5, 1), seq(10, 30, 5)),
                         Stereo_a = seq(-1, 1, 0.25),
                         seed_fold = NULL, verbose = TRUE, ...) {

  # TODO: Include null variance as parameter to avoid calculation.

  # Read data's name
  data_name <- deparse(substitute(data))

  # Prepare data
  d <- prepare_test_data(data)
  data <- d[["data"]]
  n <- d[["n"]]
  p <- d[["p"]]
  avail_stats <- d[["avail_stats"]]

  # Check K > 1
  if (K < 2) {

    stop(paste0("The number of folds (K = ", K, ") must be at least 2."))

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

  # If null_variance is given, check it is the same size of lambda_grid.
  if (!is.null(null_variance)) {

    # Check that all needed statistics are given
    checks <- stats_type %in% names(null_variance)
    if (any(!checks)) {

      stop(paste("null_variance must be a list with names containing",
                 "the tests names required in type with the values from ",
                 "null_var(...). null_variance misses",
                 paste(paste0("\"", stats_type[!checks], "\""),
                       collapse = ", "), "."))

    } else {

      # Check all parameters in grids for each specific stat_type are given
      check_grid_size <- sapply(stats_type, function(t) {

        length_grid <- switch(t,
                              "Poisson" = length(Poisson_rho),
                              "Softmax" = length(Softmax_kappa),
                              "Stereo" = length(Stereo_a))

        return(length(null_variance[[t]]) == length_grid)

      })

      if (any(!check_grid_size)) {

        stop(paste("null_variance contains the statistics",
                   paste(paste0("\"", stats_type[!check_grid_size], "\""),
                         collapse = ", "), "that misses some of the parameters",
                   "required for the grid. Check grids of parameters."))

      }

    }

  }

  # Number of statistics
  n_stats <- length(stats_type)

  # Build lambda_grid
  lambda_grid <- list("Poisson" = Poisson_rho,
                      "Softmax" = Softmax_kappa,
                      "Stereo" = Stereo_a)

  # Split data in K folds
  folds <- k_fold_split(n = n, K = K, seed = seed_fold)

  # Compute statistic
  stat <- k_fold_stat(data = data, type = stats_type, lambda_grid = lambda_grid,
                      folds = folds, K_max = K_max, verbose = verbose)

  # Calibration
  if (p_value == "MC") {

    # Get the stats_MC
    if (is.null(stats_MC)) {

      p_values <- sapply(seq_along(folds), function (k) {

        # Sample size for K - 1 folds
        # Subsample S_k
        n_not_k <- if (p == 2) {
          length(data[!((1:n) %in% folds[[k]]), ])
        } else {
          nrow(data[!((1:n) %in% folds[[k]]), , ])
        }
        # Monte Carlo exact-n null distribution
        unif_stat_MC_args <- list(n = n_not_k, type = stats_type, p = p,
                                  M = M, r_H1 = NULL, crit_val = NULL,
                                  alpha = alpha, return_stats = TRUE,
                                  stats_sorted = TRUE)
        stat_args <- lapply(stats_type, function(t) stat$lambda_hat[k, t])
        names(stat_args) <- lapply(stats_type, function(t) stat_args_name[[t]])
        stats_MC_k <- do.call(what = unif_stat_MC,
                              args = c(unif_stat_MC_args,
                                       stat_args))$stats_MC

        # p-values
        p_val <- as.double(1 - as.data.frame(sapply(stats_type, function(t) {
          sphunif:::ecdf_bin(data = stats_MC_k[[t]],
                             sorted_x = stat$statistic[k, t],
                             data_sorted = TRUE,
                             efic = TRUE, divide_n = TRUE)
        }, simplify = FALSE)))

        return(p_val)

      })


    } else {

      # Get the index of distr.idx for lambda_hat_k
      idx_lambda_hat <- sapply(stats_type, function(t) {

        return(lapply(stat$lambda_hat[, t],
                      function (lambda) which(lambda_grid[[t]] == lambda)))

      })

      # p-values
      p_values <- sapply(seq_along(folds), function(k) {

        p_val <- as.double(1 - as.data.frame(sapply(stats_type, function(t) {
          sphunif:::ecdf_bin(data = stats_MC[[paste(t,
                                                    idx_lambda_hat[k, t],
                                                    sep = ".")]],
                             sorted_x = stat$statistic[k, t],
                             data_sorted = TRUE, efic = TRUE, divide_n = TRUE)
        }, simplify = FALSE)))

        return(p_val)

      })

    }

  } else if (p_value == "asymp") {

    # TODO: Include verbose
    p_values <- sapply(seq_along(folds), function(k) {

      unif_stat_distr_args <- list(x = stat$statistic[k,], type = stats_type,
                                   p = p, n = n, approx = "asymp",
                                   stats_MC = NULL, M = M, K_max = K_max,
                                   method = method)
      stat_args <- lapply(stats_type, function(t) stat$lambda_hat[k, t])
      names(stat_args) <- lapply(stats_type, function(t) stat_args_name[[t]])
      p_val <- as.matrix(1 - do.call(what = unif_stat_distr,
                             args = c(unif_stat_distr_args, stat_args)))
      return(p_val)

    })

  } else {

    stop("Wrong choice for calibration, must be \"MC\" or \"asymp\".")

  }

  # Aggregate p_values
  p_values <- t(rbind(p_values))
  colnames(p_values) <- stats_type
  p_val <- p_val_hmp(p_vals = p_values, M = M)

  # Rejection?
  reject <- rbind(sapply(alpha, function(a) p_val < a))
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
    test[[i]] <- list(fold_statistics = stat$statistic[, i],
                      fold_params = stat$lambda_hat[, i],
                      fold_p.values = p_values[, i], p.value = p_val[i],
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

#' @rdname unif_test_cv
#' @export
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
#' @keywords internal
prepare_test_data <- function(data) {

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

  return(list("data" = data, "n" = n, "p" = p, "avail_stats" = avail_stats))
}

stat_args_name <- list("Poisson" = "Poisson_rho",
                       "Softmax" = "Softmax_kappa",
                       "Stereo" = "Stereo_a")

avail_cir_cv_tests <- c("Poisson", "Softmax")

avail_sph_cv_tests <- c("Poisson", "Softmax", "Stereo")
