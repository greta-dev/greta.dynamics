#' @name iterate_dynamic_matrix
#'
#' @title iterate dynamic transition matrices
#'
#' @description Calculate the stable population size for a stage-structured
#'   dynamical system, encoded by a transition matrix, the value of which
#'   changes at each iteration, given by function of the previous state:
#'   \code{state[t] = f(state[t-1]) \%*\% state[t-1]}.
#'
#' @details Because \code{iterate_matrix} iterates with a static transition
#'   matrix, it converges to a stable \emph{growth rate} and \emph{relative}
#'   population sizes for a dynamical system. \code{iterate_dynamic_matrix}
#'   instead uses a matrix which changes at each iteration, and can be dependent
#'   on the population sizes after the previous iteration, and the iteration
#'   number. Because this can encode density-dependence, the dynamics can
#'   converge to \emph{absolute} population sizes. The convergence criterion is
#'   therefore based on growth rates conveerging on 0.
#'
#'   As in \code{iterate_matrix}, the greta array returned by
#'   \code{matrix_function} can either be a square matrix, or a 3D array
#'   representing (on the first dimension) \emph{n} different matrices.
#'   \code{initial_state} should be shaped accordingly, as detailed in
#'   \code{iterate_matrix}.
#'
#'   To ensure the matrix is iterated for a specific number of iterations, you
#'   can set that number as \code{niter}, and set \code{tol} to 0 or a negative
#'   number to ensure that the iterations are not stopped early.
#'
#' @param matrix_function a function taking in the previous population state and
#'   the current iteration (and possibly other greta arrays) and returning a
#'   transition matrix to use for this iteration. The first two arguments must
#'   be named 'state' and 'iter', the state vector and scalar iteration number
#'   respectively. The remaining parameters must be named arguments representing
#'   (temporally static) model parameters. Variables and distributions cannot be
#'   defined inside the function.
#' @param initial_state either a column vector (with m elements) or a 3D array
#'   (with dimensions n x m x 1) giving one or more initial states from which to
#'   iterate the matrix
#' @param niter a positive integer giving the maximum number of times to iterate
#'   the matrix
#' @param tol a scalar giving a numerical tolerance, below which the algorithm
#'   is determined to have converged to a stable population size in all stages
#' @param ... optional named arguments to \code{matrix_function}, giving greta
#'   arrays for additional parameters
#' @param state_limits a numeric vector of length 2 giving minimum and maximum
#'   values at which to clamp the values of state after each iteration to
#'   prevent numerical under/overflow; i.e. elements with values below the
#'   minimum (maximum) will be set to the minimum (maximum).
#'
#' @return a named list with four greta arrays:
#' \itemize{
#'   \item{\code{stable_state}} {a vector or matrix (with the same dimensions as
#'   \code{initial_state}) giving the state after the final iteration.}
#'   \item{\code{all_states}} {an n x m x niter matrix of the state values at
#'   each iteration. This will be 0 for all entries after \code{iterations}.}
#'   \item{\code{converged}} {an integer scalar indicating whether \emph{all}
#'   the matrix iterations converged to a tolerance less than \code{tol} (1 if
#'   so, 0 if not) before the algorithm finished.}
#'   \item{\code{iterations}} {a scalar of the maximum number of iterations
#'   completed before the algorithm terminated. This should match \code{niter}
#'   if \code{converged} is \code{FALSE}.}
#' }
#'
#' @note because greta vectorises across both MCMC chains and the calculation of
#'   greta array values, the algorithm is run until all chains (or posterior
#'   samples), sites and stages have converged to stable growth. So a single
#'   value of both \code{converged} and \code{iterations} is returned, and the
#'   value of this will always have the same value in an `mcmc.list` object. So
#'   inspecting the MCMC trace of these parameters will only tell you whether
#'   the iteration converged in \emph{all} posterior samples, and the maximum
#'   number of iterations required to do so across all these samples
#'
#' @export
#'
iterate_dynamic_matrix <- function(
  matrix_function,
  initial_state,
  niter,
  tol,
  ...,
  state_limits = c(-Inf, Inf)
) {

  # generalise checking of inputs from iterate_matrix into functions
  niter <- as.integer(niter)
  state <- as.greta_array(initial_state)

  # check input dimensions
  state_dim <- dim(state)
  state_n_dim <- length(state_dim)

  check_state_is_2d_or_3d(state_n_dim)

  check_initial_state_col_vec_or_3d_dim_1(state)

  # if this is multisite
  state_multisite <- state_n_dim == 3

  # create a tensorflow function from the matrix function
  dots <- list(...)
  dots <- lapply(dots, as.greta_array)
  tf_matrix_function <- as_tf_matrix_function(
    matrix_function,
    state,
    iter = as_data(1),
    dots
    )

  # op returning a fake greta array which is actually a list containing both
  # values and states
  results <- op("iterate_dynamic_matrix",
                state,
                ...,
                operation_args = list(
                  tf_matrix_function = tf_matrix_function,
                  niter = niter,
                  tol = tol,
                  state_limits = state_limits
                ),
                tf_operation = "tf_iterate_dynamic_matrix",
                dim = c(1, 1))

  # ops to extract the components
  stable_population <- op("stable_population",
                          results,
                          tf_operation = "tf_extract_stable_population",
                          dim = state_dim)

  all_states_dim <- state_dim
  if (state_multisite) {
    all_states_dim[3] <- niter
  } else {
    all_states_dim[2] <- niter
  }

  all_states <- op("all_states",
                   results,
                   tf_operation = "tf_extract_states",
                   dim = all_states_dim)

  converged <- op("converged",
                  results,
                  tf_operation = "tf_extract_converged",
                  dim = c(1, 1))

  iterations <- op("iterations",
                   results,
                   tf_operation = "tf_extract_iterations",
                   dim = c(1, 1))

  list(stable_population = stable_population,
       all_states = all_states,
       converged = converged,
       iterations = iterations)
}

# given a greta/R function derivative function, and greta arrays for the inputs,
# return a tensorflow function taking tensors for y and t and returning a tensor
# for dydt
as_tf_matrix_function <- function (matrix_function, state, iter, dots) {

  # create a function acting on the full set of inputs, as tensors
  args <- list(r_fun = matrix_function, state = state, iter = iter)
  tf_fun <- do.call(as_tf_function, c(args, dots))

  # for CRAN's benefit
  tf_dots <- NULL

  # return a function acting only on tensors y and t, to feed to the ode solver
  function (state, iter) {

    # t will be dimensionless when used in the ode solver, need to expand out t
    # to have same dim as a scalar constant so that it can be used in the same
    # way as the greta array in the R function
    iter <- tf$reshape(iter, shape = shape(1, 1, 1))

    # tf_dots will have been added to this environment by
    # tf_iterate_dynamic_matrix
    args <- list(state = state, iter = iter)
    do.call(tf_fun, c(args, tf_dots))

  }

}

# tensorflow code
# iterate matrix tensor `matrix` `niter` times, each time using and updating vector
# tensor `state`, and return lambda for the final iteration
tf_iterate_dynamic_matrix <- function (state, ..., tf_matrix_function, niter, tol,
                                       state_limits) {

  # assign the dots (as tensors) to the matrix function's environment
  assign("tf_dots", list(...),
         environment(tf_matrix_function))

  # use a tensorflow while loop to do the recursion:
  body <- function(old_state, t_all_states, growth_rates, converged, iter, maxiter) {

    # look up the batch size from old_state and put it in the greta stash, so
    # greta can use it to appropriately create tensors for constants defined in
    # the user-provided function. Note we need to do this inside body (not
    # before), because a new control flow graph is created by tf_while_loop and
    # it otherwise becomes an unknown and causes shape variance
    batch_size <- tf$shape(old_state)[[0]]

    # note we need to access the greta stash directly here, rather than including
    # it in internals.R, because otherwise it makes a copy of the environment
    # instead and the contents can't be accessed by greta:::as_tf_function()
    assign("batch_size",
           batch_size,
           envir = greta::.internals$greta_stash)

    # create matrix (dots have been inserted into its environment, since TF
    # while loops are treacherous things)
    matrix <- tf_matrix_function(old_state, iter)

    # do matrix multiplication
    new_state <- tf$matmul(matrix, old_state, transpose_a = FALSE)

    # clamp to max and min
    new_state <- tf$clip_by_value(new_state,
                                  clip_value_min = tf_min,
                                  clip_value_max = tf_max)

    # store new state object
    t_all_states <- tf$tensor_scatter_nd_update(
      tensor = t_all_states,
      indices = iter,
      updates = tf$transpose(new_state)
    )

    # get the growth rate, and whether it has converged
    growth_rates <- new_state / old_state
    converged <- state_converged(growth_rates, tf_tol)

    list(
      new_state,
      t_all_states,
      growth_rates,
      converged,
      iter + 1L,
      maxiter
    )

  }

  # create a matrix of zeros to store all the states, but use *the transpose* so
  # it's easier to update

  # iter needs to have rank 2 for slice updating; make niter the same shape
  iter <- tf$constant(0L, shape = shape(1, 1))
  tf_niter <- tf$constant(as.integer(niter), shape = shape(1, 1))

  # add convergence tolerance and indicator
  tf_tol <- tf$constant(tol, dtype = tf_float())
  converged <- tf$constant(FALSE, dtype = tf$bool)

  # coerce limits to max and min
  tf_min <- tf$constant(state_limits[1],
                        dtype = tf_float())
  tf_max <- tf$constant(state_limits[2],
                        dtype = tf_float())

  # make a single slice (w.r.t. batch dimension) and tile along batch dimension
  state_dim <- dim(state)[-1]
  n_dim <- length(state_dim)

  shp <- to_shape(c(niter, rev(state_dim[-n_dim]), 1))
  t_all_states_slice <- tf$zeros(shp, dtype = tf_float())
  batch_size <- tf$shape(state)[[0]]
  ndim <- length(dim(t_all_states_slice))
  t_all_states <- tf$tile(t_all_states_slice, c(rep(1L, ndim - 1), batch_size))

  # create an initial growth rate, and expand its batches
  shp <- to_shape(c(1, state_dim))
  growth_rates_slice <- tf$zeros(shp, dtype = tf_float())
  growth_rates <- expand_to_batch(growth_rates_slice, state)

  # add tolerance next
  values <- list(
    state,
    t_all_states,
    growth_rates,
    converged,
    iter,
    tf_niter
  )

  # add tolerance next
  cond <- function(
    new_state,
    t_all_states,
    growth_rates,
    converged,
    iter,
    maxiter
  ) {
    tf$squeeze(tf$less(iter, maxiter)) & tf$logical_not(converged)
  }

  # iterate
  out <- tf$while_loop(cond, body, values)

  # return some elements: the transposed tensor of all the states
  list(
    state = out[[1]],
    t_all_states = out[[2]],
    growth_rates = out[[3]],
    converged = out[[4]],
    iterations = out[[5]]
  )

}

# assess convergence of the iterations to a stable growth rate. If it has
# converged, the rate of change will be the same for all stages.
state_converged <- function (growth_rates, tol) {

  # calculate the largest deviation across all stages, sites, and the batch
  # dimension and determiine whether it's acceptatble
  error <- tf$reduce_max(tf$abs(growth_rates - fl(1)))
  error < tol

}

# pull out the final population sizes
tf_extract_stable_population <- function (results) {

  state <- results$state

  # reshape if needed
  ndim <- length(dim(state))
  if (ndim > 3) {
    state <- tf$squeeze(state, ndim - 1L)
  }

  state

}

