# test functions

library(tensorflow)

# set the seed before running tests
set.seed(2017-05-01)

dag_class <- greta::.internals$inference$dag_class

randu <- function (...) {
  dim <- c(...)
  array(runif(prod(dim)), dim = dim)
}

# a baseline transition matrix
base_matrix <- function(n) {
  mat <- matrix(0, n, n)
  diag(mat) <- 0.5
  sub_diag <- (row(mat) - col(mat)) == 1
  mat[sub_diag] <- 0.4
  mat[1, n] <- 2
  mat
}

# R versions of dynamics module methods
r_iterate_matrix <- function (matrix, initial_state, niter = 100, tol = 1e-6) {

  states <- list(initial_state)

  i <- 0L
  diff <- Inf

  while(i < niter & diff > tol) {
    i <- i + 1
    states[[i + 1]] <- matrix %*% states[[i]]
    growth <- states[[i + 1]] / states[[i]]
    diffs <- growth - mean(growth)
    diff <- max(abs(diffs))
  }

  lambda <- states[[i]][1] / states[[i - 1]][1]
  stable_distribution <- states[[i]]
  stable_distribution <- stable_distribution / sum(stable_distribution)
  all_states <- matrix(0, ncol(matrix), niter)
  states_keep <- states[-1]
  all_states[, seq_along(states_keep)] <- t(do.call(rbind, states_keep))

  list(lambda = lambda,
       stable_distribution = stable_distribution,
       all_states = all_states,
       converged = as.integer(diff < tol),
       max_iter = i)
}

# R versions of dynamics module methods
r_iterate_dynamic_matrix <- function (matrix_function, initial_state, niter = 100, tol = 1e-6, ...) {

  states <- list(initial_state)

  i <- 0L
  diff <- Inf

  while(i < niter & diff > tol) {
    i <- i + 1L
    matrix <- matrix_function(states[[i]], i, ...)
    states[[i + 1]] <- matrix %*% states[[i]]
    growth <- states[[i + 1]] / states[[i]]
    diffs <- growth - 1
    diff <- max(abs(diffs))
  }

  all_states <- matrix(0, ncol(matrix), niter)
  states_keep <- states[-1]
  all_states[, seq_along(states_keep)] <- t(do.call(rbind, states_keep))

  list(stable_state = states[[i]],
       all_states = all_states,
       converged = as.integer(diff < tol),
       max_iter = i)
}

# a midpoint solver for use in deSolve, from the vignette p8
rk_midpoint <- deSolve::rkMethod(ID = "midpoint",
                                 varstep = FALSE,
                                 A = c(0, 1/2),
                                 b1 = c(0, 1),
                                 c = c(0, 1/2),
                                 stage = 2,
                                 Qerr = 1)
