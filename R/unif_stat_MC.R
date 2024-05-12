

#' @title Monte Carlo simulation of circular and (hyper)spherical uniformity
#' statistics
#'
#' @description Utility for performing Monte Carlo simulation of several
#' statistics for assessing uniformity on the (hyper)sphere
#' \eqn{S^{p-1}:=\{{\bf x}\in R^p:||{\bf x}||=1\}}{
#' S^{p-1}:=\{x\in R^p:||x||=1\}}, \eqn{p\ge 2}.
#'
#' \code{unif_stat_MC} provides a convenient wrapper for parallel
#' evaluation of \code{unif_stat}, the estimation of critical values under the
#' null distribution, and the computation of empirical powers under the
#' alternative.
#'
#' @inheritParams unif_test
#' @inheritParams r_unif
#' @param M number of Monte Carlo replications. Defaults to \code{1e4}.
#' @param r_H1 if provided, the computation of empirical powers is
#' carried out for the alternative hypothesis sampled with \code{r_H1}.
#' This must be a function with the same arguments and value as
#' \code{\link{r_unif_sph}} (see examples). Defaults to \code{NULL}, indicating
#' that the critical values are estimated from samples of \code{r_unif_sph}.
#' @param crit_val if provided, must be the critical values as returned by
#' \code{$stats_MC} in a call to \code{unif_stat_MC}. They are used for
#' computing the empirical powers of the tests present in \code{type}.
#' Defaults to \code{NULL}, which means that no power computation is done.
#' @param return_stats return the Monte Carlo statistics? If only the critical
#' values or powers are desired, \code{FALSE} saves memory in the returned
#' object. Defaults to \code{TRUE}.
#' @param stats_sorted sort the returned Monte Carlo statistics? If
#' \code{TRUE}, this is useful for evaluating faster the empirical cumulative
#' distribution function when approximating the distribution in
#' \code{\link{unif_stat_distr}}. Defaults to \code{FALSE}.
#' @param chunks number of chunks to split the \code{M} Monte Carlo
#' replications. Useful for parallelizing the simulation study in \code{chunks}
#' tasks containing \code{ceiling(M / chunks)} replications. Useful also for
#' avoiding memory bottlenecks when \code{M} is large. Defaults to
#' \cr\code{ceiling((n * M) / 1e5)}.
#' @param cores number of cores to perform the simulation. Defaults to \code{1}.
#' @param seeds if provided, a vector of size \code{chunks} for fixing the
#' seeds on each of the simulation chunks (useful for reproducing parallel
#' simulations). Specifically, for \code{k in 1:chunks}, seeds are
#' set as \code{set.seed(seeds[k], kind = "Mersenne-Twister")} in each chunk.
#' Defaults to \code{NULL} (no seed setting is done).
#' @inheritParams unif_stat
#' @param ... optional arguments to be passed to the \code{r_H1} sampler or to
#' \code{\link[foreach]{foreach}} (for example, \code{.export} to export global
#' variables or other functions to the \code{foreach} environment).
#' @return A list with the following entries:
#' \itemize{
#'   \item \code{crit_val_MC}: a data frame of size
#'   \code{c(length(alpha), length(type))}, with column names given by
#'   \code{type} and rows corresponding to the significance levels \code{alpha},
#'   that contains the estimated critical values of the tests.
#'   \item \code{power_MC}: a data frame of size
#'   \code{c(nrow(crit_val), length(type))}, with column names given by
#'   \code{type} and rows corresponding to the significance levels of
#'   \code{crit_val}, that contains the empirical powers of the tests. \code{NA}
#'   if \code{crit_val = NULL}.
#'   \item \code{stats_MC}: a data frame of size \code{c(M, length(type))}, with
#'   column names given by \code{type}, that contains the Monte Carlo
#'   statistics.
#' }
#' @details
#' It is possible to have a progress bar if \code{unif_stat_MC} is wrapped with
#' \code{\link[progressr:with_progress]{progressr::with_progress}} or if
#' \code{progressr::handlers(global = TRUE)} is invoked (once) by the user.
#' See the examples below. The progress bar is updated with the number of
#' finished chunks.
#'
#' All the tests reject for large values of the test statistic
#' (\code{max_gap = TRUE} is assumed for the Range test), so the critical
#' values for the significance levels \code{alpha} correspond to the
#' \code{alpha}-upper quantiles of the null distribution of the test statistic.
#'
#' The Monte Carlo simulation for the CCF09 test is made conditionally
#' on the choice of \code{CCF09_dirs}. That is, all the Monte Carlo statistics
#' share the same random directions.
#'
#' Except for \code{CCF09_dirs}, \code{K_CCF09}, and \code{CJ12_reg}, all the
#' test-specific parameters are vectorized.
#' @examples
#' ## Critical values
#'
#' # Single statistic, specific alpha
#' cir <- unif_stat_MC(n = 10, M = 1e2, type = "Ajne", p = 2, alpha = 0.15)
#' summary(cir$stats_MC)
#' cir$crit_val_MC
#'
#' # All circular statistics
#' cir <- unif_stat_MC(n = 10, M = 1e2, p = 2)
#' head(cir$stats_MC)
#' cir$crit_val_MC
#'
#' # All spherical statistics
#' sph <- unif_stat_MC(n = 10, M = 1e2, p = 3)
#' head(sph$stats_MC)
#' sph$crit_val_MC
#'
#' ## Using a progress bar
#'
#' # Define a progress bar
#' require(progress)
#' require(progressr)
#' handlers(handler_progress(
#'   format = ":spin [:bar] :percent Total: :elapsedfull End \u2248 :eta",
#'   clear = FALSE))
#'
#' # Call unif_stat_MC() within with_progress()
#' with_progress(unif_stat_MC(n = 10, M = 1e2, p = 3, chunks = 10))
#'
#' # With several cores
#' with_progress(unif_stat_MC(n = 10, M = 1e2, p = 3, chunks = 10, cores = 2))
#'
#' # Instead of using with_progress() each time, it is more practical to run
#' # handlers(global = TRUE)
#' # once to activate progress bars in your R session
#'
#' ## Power computation
#'
#' # Single statistic
#' cir_pow <- unif_stat_MC(n = 10, M = 1e2, type = "Ajne", p = 2,
#'                         crit_val = cir$crit_val_MC)
#' cir_pow$crit_val_MC
#' cir_pow$power_MC
#'
#' # All circular statistics
#' cir_pow <- unif_stat_MC(n = 10, M = 1e2, p = 2, crit_val = cir$crit_val_MC)
#' cir_pow$crit_val_MC
#' cir_pow$power_MC
#'
#' # All spherical statistics
#' sph_pow <- unif_stat_MC(n = 10, M = 1e2, p = 3, crit_val = sph$crit_val_MC)
#' sph_pow$crit_val_MC
#' sph_pow$power_MC
#' \donttest{
#' ## Custom r_H1
#'
#' # Circular
#' r_H1 <- function(n, p, M, l = 0.05) {
#'
#'   stopifnot(p == 2)
#'   Theta_to_X(matrix(runif(n * M, 0, (2 - l) * pi), n, M))
#'
#' }
#' dirs <- r_unif_sph(n = 5, p = 2, M = 1)[, , 1]
#' cir <- unif_stat_MC(n = 50, M = 1e2, p = 2, CCF09_dirs = dirs)
#' cir_pow <- unif_stat_MC(n = 50, M = 1e2, p = 2, r_H1 = r_H1, l = 0.10,
#'                         crit_val = cir$crit_val_MC, CCF09_dirs = dirs)
#' cir_pow$crit_val_MC
#' cir_pow$power_MC
#'
#' # Spherical
#' r_H1 <- function(n, p, M, l = 0.5) {
#'
#'   samp <- array(dim = c(n, p, M))
#'   for (j in 1:M) {
#'
#'     samp[, , j] <- mvtnorm::rmvnorm(n = n, mean = c(l, rep(0, p - 1)),
#'                                     sigma = diag(rep(1, p)))
#'     samp[, , j] <- samp[, , j] / sqrt(rowSums(samp[, , j]^2))
#'
#'   }
#'   return(samp)
#'
#' }
#' dirs <- r_unif_sph(n = 5, p = 3, M = 1)[, , 1]
#' sph <- unif_stat_MC(n = 50, M = 1e2, p = 3, CCF09_dirs = dirs)
#' sph_pow <- unif_stat_MC(n = 50, M = 1e2, p = 3, r_H1 = r_H1, l = 0.5,
#'                        crit_val = sph$crit_val_MC, CCF09_dirs = dirs)
#' sph_pow$power_MC
#'
#' ## Pre-built r_H1
#'
#' # Circular
#' dirs <- r_unif_sph(n = 5, p = 2, M = 1)[, , 1]
#' cir_pow <- unif_stat_MC(n = 50, M = 1e2, p = 2, r_H1 = r_alt, alt = "vMF",
#'                         kappa = 1, crit_val = cir$crit_val_MC,
#'                         CCF09_dirs = dirs)
#' cir_pow$power_MC
#'
#' # Spherical
#' dirs <- r_unif_sph(n = 5, p = 3, M = 1)[, , 1]
#' sph_pow <- unif_stat_MC(n = 50, M = 1e2, p = 3, r_H1 = r_alt, alt = "vMF",
#'                         kappa = 1, crit_val = sph$crit_val_MC,
#'                         CCF09_dirs = dirs)
#' sph_pow$power_MC
#' }
#' @export
unif_stat_MC <- function(n, type = "all", p, M = 1e4, r_H1 = NULL,
                         crit_val = NULL, alpha = c(0.10, 0.05, 0.01),
                         return_stats = TRUE, stats_sorted = FALSE,
                         chunks = ceiling((n * M) / 1e5), cores = 1,
                         seeds = NULL, CCF09_dirs = NULL, CJ12_reg = 3,
                         cov_a = 2 * pi, Cressie_t = 1 / 3, K_CCF09 = 25,
                         Poisson_rho = 0.5, Pycke_q = 0.5, Rayleigh_m = 1,
                         Riesz_s = 1, Rothman_t = 1 / 3,
                         Sobolev_vk2 = c(0, 0, 1), Softmax_kappa = 1,
                         Stereo_a = 0, ...) {

  # Check dimension
  if (p < 2) {

    stop("Dimension p must be p >= 2.")

  }

  # Chunk large n * M to avoid memory issues
  small_M <- M %/% chunks

  # Fix projections for the CCF09 test
  if ("CCF09" %in% type && is.null(CCF09_dirs)) {

    CCF09_dirs <- r_unif_sph(n = 50, p = p, M = 1)[, , 1]

  }

  # If the uniformity is simulated
  if (is.null(r_H1)) {

    r_H1 <- r_unif_sph

  }

  # Check if crit_val is a compatible data.frame with the output by unif_stat()
  if (!is.null(crit_val)) {

    # Dummy stats
    check_stat <- unif_stat(data = r_unif_sph(n = 2, p = p, M = 1), type = type,
                            CCF09_dirs = CCF09_dirs, CJ12_reg = CJ12_reg,
                            cov_a = cov_a, Cressie_t = Cressie_t,
                            K_CCF09 = K_CCF09, Poisson_rho = Poisson_rho,
                            Pycke_q = Pycke_q, Rayleigh_m = Rayleigh_m,
                            Riesz_s = Riesz_s, Rothman_t = Rothman_t,
                            Sobolev_vk2 = Sobolev_vk2, Softmax_kappa =
                              Softmax_kappa, Stereo_a = Stereo_a)


    # Names check
    checks <- names(check_stat) %in% colnames(crit_val)
    if (any(!checks)) {

      stop(paste("crit_val must be a data.frame with colnames containing",
                 "the tests names returned by unif_stat(...). crit_val misses",
                 paste(paste0("\"", colnames(check_stat)[!checks], "\""),
                       collapse = ", "), "."))

    }

  }

  # Extra arguments for foreach::foreach and r_H1
  dots <- list(...)
  foreach_args <- dots[names(dots) %in% names(formals(foreach::foreach))]
  r_H1_args <- dots[names(dots) %in% names(formals(r_H1))]

  # Check seeds
  if (!is.null(seeds) && length(seeds) != chunks) {

    warning(paste("seeds and chunks do not have the same length:",
                  "seeds are ignored."))
    seeds <- NULL

  }

  # Parallel backend
  old_dopar <- doFuture::registerDoFuture()
  old_plan <- future::plan(future::multisession(), workers = cores)
  on.exit({

    with(old_dopar, foreach::setDoPar(fun = fun, data = data, info = info))
    future::plan(old_plan)

  })
  `%op%` <- doRNG::`%dorng%`

  # Measure progress?
  if (requireNamespace("progressr", quietly = TRUE)) {

    prog <- progressr::progressor(along = 1:chunks)

  }

  # Monte Carlo
  k <- 0
  stats <- do.call(what = foreach::foreach,
                   args = c(list(k = 1:chunks, .combine = rbind,
                                 .inorder = TRUE, .multicombine = TRUE,
                                 .maxcombine = 100, .packages = "sphunif"),
                            foreach_args)) %op% {

    # Samples
    if (!is.null(seeds)) {

      set.seed(seeds[k], kind = "Mersenne-Twister")

    }
    X <- do.call(what = r_H1, args = c(list(n = n, p = p, M = small_M),
                                            r_H1_args))

    # Statistics
    stats <- unif_stat(data = X, type = type, CCF09_dirs = CCF09_dirs,
                       CJ12_reg = CJ12_reg, cov_a = cov_a,
                       Cressie_t = Cressie_t, K_CCF09 = K_CCF09,
                       Poisson_rho = Poisson_rho, Pycke_q = Pycke_q,
                       Rayleigh_m = Rayleigh_m, Riesz_s = Riesz_s,
                       Rothman_t = Rothman_t, Sobolev_vk2 = Sobolev_vk2,
                       Softmax_kappa =
                         Softmax_kappa, Stereo_a = Stereo_a)

    # Remove X
    rm(X)

    # Signal progress
    if (requireNamespace("progressr", quietly = TRUE)) {

      prog()

    }

    # Return stats
    stats

  }

  # Sort statistics
  if (stats_sorted && return_stats) {

    names <- names(stats)
    stats <- as.data.frame(sort_each_col(as.matrix(stats)))
    names(stats) <- names

  }

  # Build tables
  if (is.null(crit_val)) {

    # Critical values
    crit_val <- rbind(apply(stats, 2, quantile, probs = 1 - alpha,
                            na.rm = TRUE))
    crit_val <- as.data.frame(crit_val)
    rownames(crit_val) <- alpha

    # Power
    power <- NA

  } else {

    # In case crit_val contains more tests than the ones considered in type
    # and in a potentially different order, match them accordingly with stats
    names_stats <- colnames(stats)
    ind_tests <- match(x = names_stats, table = colnames(crit_val))
    if (any(is.na(ind_tests))) {

      stop(paste0(
        "No critical values in \"crit_val\" associated with the test(s) ",
        paste(paste0("\"", names_stats[is.na(ind_tests)], "\""),
              collapse = ", "), "."))

    }
    crit_val <- crit_val[, ind_tests, drop = FALSE]

    # Empirical rejection level
    power <- t(apply(crit_val, 1,
                     function(c) rowMeans(t(stats) > c, na.rm = TRUE)))
    if (ncol(crit_val) == 1) {

      power <- t(power)
      colnames(power) <- names_stats

    }
    power <- as.data.frame(power)

  }

  # Empty stats if not desired
  if (!return_stats) {

    stats <- NA

  }

  return(list("crit_val_MC" = crit_val, "power_MC" = power, "stats_MC" = stats))

}
