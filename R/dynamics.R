# methods for using stage-structured (matrix) population models using greta
# dynamics module

op <- greta::.internals$nodes$constructors$op

# tensorflow code
# iterate matrix tensor `mat` `niter` times, each time using and updating vector
# tensor `state`, and return lambda for the final iteration
tf_iterate_lambda <- function (mat, state, niter) {

  # store states (can't overwrite since we need to maintain the chain of nodes)
  states <- list(state)

  # iterate the matrix
  for (i in seq_len(niter))
    states[[i + 1]] <-  tf$matmul(mat, states[[i]], transpose_a = TRUE)

  # return the final growth rate (should be same for all states at convergence)
  lambda <- states[[niter + 1]][1] / states[[niter]][1]
  tf$reshape(lambda, shape = c(1L, 1L))

}

# iterate matrix tensor `mat` `max(niter)` times, each time using and updating vector
# tensor `state`, and return states corresponding to niter
tf_iterate_state <- function (mat, state,
                              dens_param,
                              niter,
                              dens_form) {
  
  # store states (can't overwrite since we need to maintain the chain of nodes)
  states <- list(state)
  
  # iterate the matrix
  for (i in seq_len(max(niter))) {
    
    # include density dependence
    nkm1 <- tf$reduce_sum(states[[i]])
    scale_factor <- switch(dens_form,
                           "bh" = tf$divide(nkm1, tf$add(tf$constant(1, dtype = tf$float32),
                                                         tf$multiply(dens_param, nkm1))),
                           "ricker" = tf$multiply(nkm1,
                                                  tf$exp(tf$multiply(tf$multiply(tf$constant(-1, dtype = tf$float32),
                                                                                 dens_param),
                                                                     nkm1))),
                           tf$constant(1, dtype = tf$float32))
    
    # update state
    mat_tmp <- tf$multiply(scale_factor, mat)
    states[[i + 1]] <-  tf$matmul(mat_tmp, states[[i]], transpose_a = TRUE)
  }
  
  # return the final state
  do.call(greta::.internals$tensors$tf_cbind, states[niter + 1])
  
}

# apply iterate_lambda to a series of n matrices of dimension m, stored as an n
# x m^2 matrix, each row being unpacked *rowwise*
tf_iterate_lambda_vectorised <- function (mat, state, n, m, niter) {

  # create indices for a block-diagonal sparse matrix
  # column wise offset to move the matrices along one horizontally
  offset <- rep(seq_len(n) - 1, each = m ^ 2) * m
  idx_cols <- rep(seq_len(m), m * n) + offset
  idx_rows <- rep(rep(seq_len(m), each = m), n) + offset

  # set up sparse tensor
  indices <- tf$constant(cbind(idx_rows, idx_cols) - 1, dtype = tf$int64)
  values <- tf$reshape(mat, shape = shape(n * m ^ 2))
  shape <- tf$constant(as.integer(rep(n * m, 2)), dtype = tf$int64)
  full_mat <- tf$SparseTensor(indices = indices,
                              values = values,
                              dense_shape = shape)

  # replicate state for all patches
  state <- tf$tile(state,
                   tf$constant(c(n, 1L),
                               dtype = tf$int32,
                               shape = shape(2)))

  # store states (can't overwrite since we need to maintain the chain of nodes)
  states <- list(state)

  # iterate the matrix
  for (i in seq_len(niter))
    states[[i + 1]] <-  tf$sparse_tensor_dense_matmul(full_mat, states[[i]])

  # indices to subset the first state for each patch
  idx_begin <- tf$constant(c(0L, 0L), dtype = tf$int32)
  idx_end <- tf$constant(c(m * n, 1L), dtype = tf$int32)
  idx_strides <- tf$constant(c(m, 1L), dtype = tf$int32)

  # subset these iterated states to get the first state per sub-matrix
  before <- tf$strided_slice(states[[niter]],
                             begin = idx_begin,
                             end = idx_end,
                             strides = idx_strides)

  after <- tf$strided_slice(states[[niter + 1]],
                            begin = idx_begin,
                            end = idx_end,
                            strides = idx_strides)

  # get growth rate
  lambdas <- after / before
  lambdas

}

#' @name gretaDynamics
#' @title iterate transition matrices
#'
#' @description greta functions to calculate the intrinsic growth rate or stable
#'   stage distribution for stage-structured dynamical systems, encoded as greta
#'   arrays representing transition matrices
#'
#' @details \code{iterate_lambda} iterates a matrix a certain number of
#'   times and returns, as a scalar greta array, the terminal growth rate for
#'   the first element of the state vector. \code{iterate_state} carries out the
#'   same procedure, but returns the final state vector.
#'   \code{iterate_lambda_vectorised} is a vectorised version of
#'   \code{iterate_lambda} for iterating over multiple matrices, returning a
#'   vector of growth rates.
#'
#' @import greta
#' @importFrom tensorflow tf shape
NULL

#' @name iterate_state
#' @rdname gretaDynamics
#'
#' @param matrix a square, two-dimensional (i.e. matrix-like) greta array
#'   representing transition probabilities between states
#' @param state a column vector greta array representing the initial state from
#'   which to iterate the matrix
#' @param dens_param details
#' @param niter a positive integer giving the number of times to iterate the
#'   matrix
#' @param dens_form details
#'
#' @export
iterate_state <- function(matrix, state,
                          dens_param,
                          niter,
                          dens_form) {
  
  niter <- as.integer(niter)
  
  dimfun <- function(elem_list) {
    
    # input dimensions
    matrix_dim <- dim(elem_list[[1]])
    state_dim <- dim(elem_list[[2]])
    
    if (length(state_dim) != 2 | state_dim[2] != 1)
      stop ('state must be a column vector greta array',
            call. = FALSE)
    
    if (length(matrix_dim) != 2 | matrix_dim[1] != matrix_dim[2])
      stop ('matrix must be a two-dimensional square greta array',
            call. = FALSE)
    
    if (matrix_dim[2] != state_dim[1])
      stop ('number of elements in state must match the dimension of matrix',
            call. = FALSE)
    
    # output dimensions
    c(state_dim[1], length(niter))
  }
  
  op('iterate_state',
     matrix,
     state,
     dens_param,
     operation_args = list(niter = niter,
                           dens_form = dens_form),
     tf_operation = tf_iterate_state,
     dimfun = dimfun)
  
}

#' @name iterate_lambda
#' @rdname gretaDynamics
#' @export
iterate_lambda <- function(matrix, state, niter) {

  niter <- as.integer(niter)

  dimfun <- function(elem_list) {

    # input dimensions
    matrix_dim <- dim(elem_list[[1]])
    state_dim <- dim(elem_list[[2]])

    if (length(state_dim) != 2 | state_dim[2] != 1)
      stop ('state must be a column vector greta array',
            call. = FALSE)

    if (length(matrix_dim) != 2 | matrix_dim[1] != matrix_dim[2])
      stop ('matrix must be a two-dimensional square greta array',
            call. = FALSE)

    if (matrix_dim[2] != state_dim[1])
      stop ('number of elements in state must match the dimension of matrix',
            call. = FALSE)

    # output dimensions
    c(1, 1)
  }

  op('iterate_lambda',
     matrix,
     state,
     operation_args = list(niter = niter),
     tf_operation = tf_iterate_lambda,
     dimfun = dimfun)

}

#' @name iterate_lambda_vectorised
#' @rdname gretaDynamics
#'
#' @param matrices a rectangular two-dimensional greta array of dimension n x
#'   m^2, each row of which gives the rowwise elements of a different m x m
#'   matrix to iterate
#' @param n the number of m x m matrices to be iterated (first dimensions of
#'   \code{matrices})
#' @param m the dimension of each matrix to be iterated
#'
#' @export
iterate_lambda_vectorised <- function(matrices, state, n, m, niter) {

  n <- as.integer(n)
  m <- as.integer(m)
  niter <- as.integer(niter)

  dimfun <- function(elem_list) {

    # input dimensions
    matrices_dim <- dim(elem_list[[1]])
    state_dim <- dim(elem_list[[2]])

    if (length(state_dim) != 2 | state_dim[2] != 1)
      stop ('state must be a column vector greta array',
            call. = FALSE)

    if (m != state_dim[1])
      stop ('number of elements in state must match the dimension of matrix',
            call. = FALSE)

    if (length(matrices_dim) != 2 | matrices_dim[2] != (m ^ 2) | matrices_dim[1] != n)
      stop ('matrix must be a rectangular greta array with dimensions n x m^2',
            call. = FALSE)

    # output dimensions
    c(n, 1)

  }

  op('iterate_lambda_vectorised',
     matrices,
     state,
     operation_args = list(n = n,
                           m = m,
                           niter = niter),
     tf_operation = tf_iterate_lambda_vectorised,
     dimfun = dimfun)

}
