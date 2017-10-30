# test functions

library(tensorflow)

# set the seed before running tests
set.seed(2017-05-01)

dag_class <- greta::.internals$inference$dag_class

# evaluate a greta_array, node, or tensor
grab <- function (x) {

  if (inherits(x, "node"))
    x <- as.greta_array(x)

  if (inherits(x, "greta_array")) {
    dag <- dag_class$new(list(x))
    x$node$define_tf(dag)
    x <- get(dag$tf_name(x$node),
             envir = dag$tf_environment)
  }

  tf$Session()$run(x)

}

randu <- function (...) {
  dim <- c(...)
  array(runif(prod(dim)), dim = dim)
}

# R versions of dynamics module methods
it_lambda <- function (matrix, state, niter) {
  states <- list(state)
  for (i in seq_len(niter))
    states[[i + 1]] <- states[[i]] %*% matrix
  states[[niter + 1]][1] / states[[niter]][1]
}

it_state <- function (matrix, state,
                      dens_param,
                      niter,
                      dens_form) {
  for (i in seq_len(niter)) {
    # include density dependence
    nkm1 <- sum(state)
    scale_factor <- switch(dens_form,
                           "bh" = (nkm1 / (1 + dens_param * nkm1)),
                           "ricker" = (nkm1 * exp(-dens_param * nkm1)),
                           1)

    # update state
    mat_tmp <- scale_factor * matrix
    state <- state %*% mat_tmp
  }
  state[1, ]
}
