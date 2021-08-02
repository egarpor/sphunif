

#' @title Monte Carlo integration of functions on the (hyper)sphere
#'
#' @description Monte Carlo approximation of the integral
#' \deqn{\int_{S^{p-1}}f(x)\,\mathrm{d}x}{\int_{S^{p-1}} f(x) dx}
#' of a function \eqn{f:S^{p-1} \rightarrow R} defined on the (hyper)sphere
#' \eqn{S^{p-1}:=\{{\bf x}\in R^p:||{\bf x}||=1\}}{
#' S^{p-1}:=\{x\in R^p:||x||=1\}}, \eqn{p\ge 2}.
#'
#' @param f function to be integrated. Its first argument must be the
#' (hyper)sphere position. Must be vectorized and return a vector of size
#' \code{nrow(x)} for a matrix input \code{x}. See examples.
#' @inheritParams r_unif
#' @param M number of Monte Carlo samples. Defaults to \code{1e4}.
#' @param chunks number of chunks to split the \code{M} Monte Carlo
#' samples. Useful for parallelizing the integration in \code{chunks}
#' tasks containing \code{ceiling(M / chunks)} replications. Useful also for
#' avoiding memory bottlenecks when \code{M} is large. Defaults to
#' \cr\code{ceiling(M / 1e3)}.
#' @param cores number of cores to perform the integration. Defaults to
#' \code{1}.
#' @inheritParams unif_stat_MC
#' @param ... optional arguments to be passed to \code{f} or to
#' \code{\link[foreach]{foreach}} (for example, \code{.export} to export global
#' variables or other functions to the \code{foreach} environment).
#' @return A scalar with the approximate integral.
#' @examples
#' # Vectorized functions to be integrated
#' x1 <- function(x) x[, 1]
#' quad <- function(x, a = 0) a + rowSums(x^4)
#'
#' # Approximate \int_{S^{p-1}} x_1 dx = 0
#' int_sph_MC(f = x1, p = 3, cores = 1, M = 1e4, chunks = 2)
#' int_sph_MC(f = x1, p = 3, cores = 2, M = 1e5, chunks = 2)
#'
#' # Approximate \int_{S^{p-1}} (a + \sum_i x_i^4) dx
#' int_sph_MC(f = quad, p = 2, cores = 1, M = 1e4, a = 0, chunks = 2)
#' int_sph_MC(f = quad, p = 2, cores = 1, M = 1e4, a = 1, chunks = 2)
#'
#' # Compare with Gauss--Legendre integration on S^2
#' th_k <- Gauss_Legen_nodes(a = 0, b = 2 * pi, N = 40)
#' w_k <- Gauss_Legen_weights(a = 0, b = 2 * pi, N = 40)
#' sum(w_k * quad(cbind(cos(th_k), sin(th_k)), a = 1))
#' @export
int_sph_MC <- function(f, p, M = 1e4, cores = 1, chunks = ceiling(M / 1e3),
                       verbose = TRUE, seeds = NULL, ...) {

  # Parallel backend
  cl <- parallel::makeCluster(cores)
  doSNOW::registerDoSNOW(cl = cl)
  `%op%` <- foreach::`%dopar%`

  # Show progress?
  if (verbose) {

    pb <- utils::txtProgressBar(max = chunks, style = 3)
    progress <- function(n) utils::setTxtProgressBar(pb = pb, value = n)
    opts <- list(progress = progress)

  } else {

    opts <- list()

  }

  # Extra arguments for foreach::foreach and f
  dots <- list(...)
  foreach_args <- dots[names(dots) %in% names(formals(foreach::foreach))]
  f_args <- dots[names(dots) %in% names(formals(f))]
  name_arg1_f <- names(formals(f))[1]
  names_args_f <- names(f_args)

  # Check seeds
  if (!is.null(seeds) & length(seeds) != chunks) {

    warning(paste("seeds and chunks do not have the same length,",
                  "seeds are ignored."))
    seeds <- NULL

  }

  # Monte Carlo integration
  k <- 0
  small_M <- M / chunks
  int <- do.call(what = foreach::foreach,
                 args = c(list(k = 1:chunks, .combine = "+",
                               .options.snow = opts, .inorder	= FALSE,
                               .packages = "sphunif"),
                          foreach_args)) %op% {

    # Sample uniform data
    if (!is.null(seeds)) {

      set.seed(seeds[k], kind = "Mersenne-Twister")

    }
    X <- r_unif_sph(n = small_M, p = p, M = 1)[, , 1]

    # Evaluate f
    args <- c(list(X), f_args)
    names(args) <- c(name_arg1_f, names_args_f)
    sum(do.call(what = f, args = args))

  }

  # Close loop
  parallel::stopCluster(cl)
  if (verbose) {

    close(pb)

  }

  # Return integral approximation
  int <- rotasym::w_p(p = p) * int / M
  return(int)

}

