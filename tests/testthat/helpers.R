# test functions

library(tensorflow)

# set the seed before running tests
set.seed(2017-05-01)

dag_class <- greta::.internals$inference$dag_class

randu <- function (...) {
  dim <- c(...)
  array(runif(prod(dim)), dim = dim)
}

# R versions of dynamics module methods
r_iterate_matrix <- function (matrix, state, niter) {
  states <- list(state)
  for (i in seq_len(niter))
    states[[i + 1]] <- states[[i]] %*% matrix

  lambda <- states[[niter + 1]][1] / states[[niter]][1]
  stable_distribution <- t(states[[niter + 1]])
  stable_distribution <- stable_distribution / sum(stable_distribution)
  all_states <- t(do.call(rbind, states[-1]))

  list(lambda = lambda,
       stable_distribution = stable_distribution,
       all_states = all_states)
}

