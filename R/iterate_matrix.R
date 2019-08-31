#' @name iterate_matrix
#'
#' @title iterate transition matrices
#'
#' @description calculate the intrinsic growth rate(s) and stable stage
#'   distribution(s) for a stage-structured dynamical system, encoded as a
#'   transition matrix
#'
#' @details \code{iterate_matrix} can either act on single transition matrix and
#'   initial state (if \code{matrix} is 2D and \code{initial_state} is a column
#'   vector), or it can simultaneously act on \emph{n} different matrices and/or
#'   \emph{n} different initial states (if \code{matrix} and
#'   \code{initial_state} are 3D arrays). In the latter case, the first
#'   dimension of both objects should be the batch dimension \emph{n}.
#'
#' @param matrix either a square 2D transition matrix (with dimensions m x m),
#'   or a 3D array (with dimensions n x m x m), giving one or more transition
#'   matrices to iterate
#' @param initial_state either a column vector (with m elements) or a 3D array
#'   (with dimensions n x m x 1) giving one or more initial states from which to
#'   iterate the matrix
#' @param niter a positive integer giving the number of times to iterate the
#'   matrix
#'
#' @return a named list with three greta arrays: \code{lambda} a scalar or
#'   vector giving the ratio of the first stage values between the final two
#'   iterations, \code{stable_state} a vector or matrix (with the same
#'   dimensions as \code{initial_state}) giving the state after the final
#'   iteration, normalised so that the values for all stages sum to one, and
#'   \code{all_states} an n x m x niter matrix of the state values at each
#'   iteration. If the system has converged in \code{niter} iterations,
#'   \code{lambda} and \code{stable_state} correspond to the asymptotic growth
#'   rate and stable stage distribution respectively.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # simulate from a probabilistic 4-stage transition matrix model
#' k <- 4
#'
#' # component variables
#' # survival probability for all stages
#' survival <- uniform(0, 1, dim = k)
#' # conditional (on survival) probability of staying in a stage
#' stasis <- c(uniform(0, 1, dim = k - 1), 1)
#' # marginal probability of staying/progressing
#' stay <- survival * stasis
#' progress <- (survival * (1 - stay))[1:(k - 1)]
#' # recruitment rate for the largest two stages
#' recruit <- exponential(c(3, 5))
#'
#' # combine into a matrix:
#' tmat <- zeros(k, k)
#' diag(tmat) <- stay
#' progress_idx <- row(tmat) - col(tmat) == 1
#' tmat[progress_idx] <- progress
#' tmat[1, k - (1:0)] <- recruit
#'
#' # analyse this to get the intrinsic growth rate and stable state
#' iterations <- iterate_matrix(tmat)
#' iterations$lambda
#' iterations$stable_distribution
#' iterations$all_states
#'
#' # Can also do this simultaneously for a collection of transition matrices
#' k <- 2
#' n <- 10
#' survival <- uniform(0, 1, dim = c(n, k))
#' stasis <- cbind(uniform(0, 1, dim = n), rep(1, n))
#' stay <- survival * stasis
#' progress <- (survival * (1 - stasis))[, 1]
#' recruit_rate <- 1 / seq(0.1, 5, length.out = n)
#' recruit <- exponential(recruit_rate, dim = n)
#' tmats <- zeros(10, 2, 2)
#' tmats[, 1, 1] <- stasis[, 1]
#' tmats[, 2, 2] <- stasis[, 2]
#' tmats[, 2, 1] <- progress
#' tmats[, 1, 2] <- recruit
#'
#' iterations <- iterate_matrix(tmats)
#' iterations$lambda
#' iterations$stable_distribution
#' iterations$all_states
#'
#' }
iterate_matrix <- function(matrix,
                           initial_state = rep(1, ncol(matrix)),
                           niter = 100) {

  niter <- as.integer(niter)
  matrix <- as.greta_array(matrix)
  state <- as.greta_array(initial_state)

  # check input dimensions
  matrix_dim <- dim(matrix)
  state_dim <- dim(state)
  matrix_n_dim <- length(matrix_dim)
  state_n_dim <- length(state_dim)

  if (!matrix_n_dim %in% 2:3 | !state_n_dim %in% 2:3) {
    stop ("matrix and state must be either two- or three-dimensional",
          call. = FALSE)
  }

  # ensure the last dim of state is 1
  if (state_dim[state_n_dim] != 1) {
    stop ("initial_state must be either a column vector, ",
          "or a 3D array with final dimension 1",
          call. = FALSE)
  }

  # if this is multisite
  matrix_multisite <- matrix_n_dim == 3
  state_multisite <- state_n_dim == 3
  multisite <- matrix_multisite | state_multisite

  # ensure the site dimension matches
  if (multisite) {

    if (!state_multisite) {

      n <- matrix_dim[1]

      # expand state
      state_list <- replicate(n, state, simplify = FALSE)
      state <- t(do.call(cbind, state_list))

      # add trailing dimension back on
      dim(state) <- c(dim(state), 1)

      state_dim <- dim(state)
      state_n_dim <- 3

    }

    if (!matrix_multisite) {

      n <- state_dim[1]

      # expand matrix
      dim(matrix) <- c(1, dim(matrix))
      matrix_list <- replicate(n, matrix, simplify = FALSE)
      matrix <- do.call(abind, c(matrix_list, list(along = 1)))

      matrix_dim <- dim(matrix)
      matrix_n_dim <- 3

    }

    if (matrix_multisite & state_multisite) {
      n <- matrix_dim[1]
      if (state_dim[1] != n) {
        stop ("if matrix is 3D and initial_state is a matrix",
              "the batch dimension (n) must match",
              call. = FALSE)
      }
    }

  } else {
    n <- 1
  }

  # check the number of stages matches
  m <- matrix_dim[2]
  if (multisite) {
    matrix_raw_dim <- matrix_dim[2:3]
    state_raw_dim <- state_dim[2]
  } else {
    matrix_raw_dim <- matrix_dim
    state_raw_dim <- state_dim[1]
  }

  if (!all(matrix_raw_dim == m)) {
    stop ("each matrix must be a two-dimensional square greta array ",
          call. = FALSE)
  }

  if (state_raw_dim != m) {
    stop ("length of each initial_state must match the dimension of matrix",
          call. = FALSE)
  }

  # op returning a fake greta array which is actually a list containing both
  # values and states
  states <- op('iterate_matrix',
               matrix,
               state,
               operation_args = list(niter = niter),
               tf_operation = "tf_iterate_matrix",
               dim = c(1, 1))

  # ops to extract the components
  lambda <- op('lambda',
               states,
               tf_operation = "tf_extract_lambda",
               dim = c(n, 1))

  stable_distribution <- op('stable_distribution',
                            states,
                            tf_operation = "tf_extract_stable_distribution",
                            dim = state_dim)

  all_states_dim <- state_dim
  if (multisite) {
    all_states_dim[3] <- niter
  } else {
    all_states_dim[2] <- niter
  }

  all_states <- op('all_states',
               states,
               tf_operation = "tf_extract_states",
               dim = all_states_dim)

  list(lambda = lambda,
       stable_distribution = stable_distribution,
       all_states = all_states)

}

op <- greta::.internals$nodes$constructors$op
as.greta_array <- greta::.internals$greta_arrays$as.greta_array
tf_sweep <- greta::.internals$tensors$tf_sweep
tf_colsums <- greta::.internals$tensors$tf_colsums

# tensorflow code
# iterate matrix tensor `mat` `niter` times, each time using and updating vector
# tensor `state`, and return lambda for the final iteration
tf_iterate_matrix_old <- function (mat, state, niter) {

  states <- list(state)

  # iterate the matrix, storing states in a list
  for (i in seq_len(niter))
    states[[i + 1]] <-  tf$matmul(mat, states[[i]], transpose_a = TRUE)

  # return all the states
  states

}

# tensorflow code
# iterate matrix tensor `matrix` `niter` times, each time using and updating vector
# tensor `state`, and return lambda for the final iteration
tf_iterate_matrix <- function (matrix, state, niter) {

  # handle variable input dimensions?

  # use a tensorflow while loop to do the recursion:
  body <- function(matrix, old_state, t_all_states, iter, maxiter) {

    # do matrix multiplication
    new_state <- tf$matmul(matrix, old_state, transpose_a = TRUE)

    # store new state object
    t_all_states <- tf$tensor_scatter_nd_update(
      tensor = t_all_states,
      indices = iter,
      updates = tf$transpose(new_state)
    )

    list(matrix, new_state, t_all_states, iter + 1L, maxiter)
  }

  # create a matrix of zeros to store all the states, but use *the transpose* so
  # it's easier to update

  # iter needs to have rank 2 for slice updating; make niter the same shape
  iter <- tf$constant(1L, shape = shape(1, 1))
  tf_niter <- tf$constant(as.integer(niter), shape = shape(1, 1))

  # make a single slice (w.r.t. batch dimension) and tile along batch dimension
  state_dim <- dim(state)[-1]
  n_dim <- length(state_dim)

  shp <- to_shape(c(niter, rev(state_dim[-n_dim]), 1))
  t_all_states_slice <- tf$zeros(shp, dtype = tf_float())
  batch_size <- tf$shape(matrix)[[0]]
  ndim <- length(dim(t_all_states_slice))
  t_all_states <- tf$tile(t_all_states_slice, c(rep(1L, ndim - 1), batch_size))

  # add tolerance next
  values <- list(matrix,
                 state,
                 t_all_states,
                 iter,
                 tf_niter)

  # add tolerance next
  cond <- function(matrix, new_state, t_all_states, iter, maxiter) {
    tf$squeeze(tf$less(iter, maxiter))
  }

  # iterate
  out <- tf$while_loop(cond,
                       body,
                       values)

  # return the tensor of all the states
  all_states <- tf$transpose(out[[3]])
  all_states

}



# return the ratio of the first stage values for the last two states, which
# should be the intrinsic growth rate if the iteration has converged
tf_extract_lambda_old <- function (states) {
  niter <- length(states)

  # handle possible site dimension, and squeeze additional dimension from this
  if (length(dim(states[[1]])) == 4) {
    before <- states[[niter - 1]][, , 0, 0, drop = FALSE]
    after <- states[[niter]][, , 0, 0, drop = FALSE]
    lambda <- after / before
    lambda <- tf$squeeze(lambda, axis = 3L)
  } else {
    before <-  states[[niter - 1]][, 0, 0, drop = FALSE]
    after <- states[[niter]][, 0, 0, drop = FALSE]
    lambda <- after / before
  }

  lambda

}

tf_extract_lambda <- function (states) {

  niter <- dim(states)[[4]]

  # handle possible site dimension, and squeeze additional dimension from this
  if (length(dim(states[[1]])) == 4) {
    before <- states[[niter - 1]][, , 0, 0, drop = FALSE]
    after <- states[[niter]][, , 0, 0, drop = FALSE]
    lambda <- after / before
    lambda <- tf$squeeze(lambda, axis = 3L)
  } else {
    before <-  states[[niter - 1]][, 0, 0, drop = FALSE]
    after <- states[[niter]][, 0, 0, drop = FALSE]
    lambda <- after / before
  }

  lambda

}

# return the final state from matrix iteration (should have stabilised)
tf_extract_stable_distribution <- function (states) {
  niter <- length(states)
  state <- states[[niter]]
  axis <- length(dim(state)) - 2L
  sums <- tf$reduce_sum(state, axis = axis, keepdims = TRUE)
  tf$truediv(state, sums)
}

# return the final state from matrix iteration (should have stabilised)
tf_extract_states <- function (states) {
  axis <- length(dim(states[[1]])) - 1L
  tf$concat(states[-1], axis = axis)
}
