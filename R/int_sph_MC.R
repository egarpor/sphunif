

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
#' @details
#' It is possible to have a progress bar if \code{int_sph_MC} is wrapped with
#' \code{\link[progressr:with_progress]{progressr::with_progress}} or if
#' \code{progressr::handlers(global = TRUE)} is invoked (once) by the user.
#' See the examples below. The progress bar is updated with the number of
#' finished chunks.
#' @examples
#' ## Sequential simulation
#'
#' # Vectorized functions to be integrated
#' x1 <- function(x) x[, 1]
#' quad <- function(x, a = 0) a + rowSums(x^4)
#'
#' # Approximate \int_{S^{p-1}} x_1 dx = 0
#' int_sph_MC(f = x1, p = 3, M = 1e4, chunks = 2)
#'
#' # Approximate \int_{S^{p-1}} (a + \sum_i x_i^4) dx
#' int_sph_MC(f = quad, p = 2, M = 1e4, a = 0, chunks = 2)
#'
#' # Compare with Gauss--Legendre integration on S^2
#' th_k <- Gauss_Legen_nodes(a = 0, b = 2 * pi, N = 40)
#' w_k <- Gauss_Legen_weights(a = 0, b = 2 * pi, N = 40)
#' sum(w_k * quad(cbind(cos(th_k), sin(th_k)), a = 1))
#'
#' ## Parallel simulation with a progress bar
#' \donttest{
#' # Define a progress bar
#' require(progress)
#' require(progressr)
#' handlers(handler_progress(
#'   format = ":spin [:bar] :percent Total: :elapsedfull End \u2248 :eta",
#'   clear = FALSE))
#'
#' # Call int_sph_MC() within with_progress()
#' with_progress(int_sph_MC(f = x1, p = 3, cores = 2, M = 1e5, chunks = 100))
#'
#' # Instead of using with_progress() each time, it is more practical to run
#' # handlers(global = TRUE)
#' # once to activate progress bars in your R session
#' }
#' @export
int_sph_MC <- function(f, p, M = 1e4, cores = 1, chunks = ceiling(M / 1e3),
                       seeds = NULL, ...) {

  # Chunk large n * M to avoid memory issues
  small_M <- M %/% chunks

  # Extra arguments for foreach::foreach and f
  dots <- list(...)
  foreach_args <- dots[names(dots) %in% names(formals(foreach::foreach))]
  f_args <- dots[names(dots) %in% names(formals(f))]
  name_arg1_f <- names(formals(f))[1]
  names_args_f <- names(f_args)

  # Check seeds
  if (!is.null(seeds) & length(seeds) != chunks) {

    warning(paste("seeds and chunks do not have the same length:",
                  "seeds are ignored."))
    seeds <- NULL

  }

  # Parallel backend
  old_dopar <- doFuture::registerDoFuture()
  old_plan <- future::plan(future::multisession(), workers = cores)
  options(future.rng.onMisuse = "ignore")
  on.exit({

    with(old_dopar, foreach::setDoPar(fun = fun, data = data, info = info))
    future::plan(old_plan)
    options(future.rng.onMisuse = NULL)

  })
  `%op%` <- foreach::`%dopar%`

  # Measure progress?
  if (requireNamespace("progressr", quietly = TRUE)) {

    prog <- progressr::progressor(along = 1:chunks)

  }

  # Monte Carlo
  k <- 0
  int <- do.call(what = foreach::foreach,
                 args = c(list(k = 1:chunks, .combine = "sum",
                               .inorder = FALSE, .multicombine = TRUE,
                               .packages = "sphunif"),
                          foreach_args)) %op% {

    # Sample uniform data
    if (!is.null(seeds)) {

      set.seed(seeds[k], kind = "Mersenne-Twister")

    }
    X <- rbind(r_unif_sph(n = small_M, p = p, M = 1)[, , 1])

    # Evaluate f
    args <- c(list(X), f_args)
    names(args) <- c(name_arg1_f, names_args_f)
    sf <- sum(do.call(what = f, args = args))

    # Remove X and clean memory
    rm(X)
    gc()

    # Signal progress
    if (requireNamespace("progressr", quietly = TRUE)) {

      prog()

    }

    # Return sum
    sf

  }

  # Return integral approximation
  int <- rotasym::w_p(p = p) * int / M
  return(int)

}
