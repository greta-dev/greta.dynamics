context("dynamic function iteration")

test_that("single iteration works", {

  skip_if_not(greta:::check_tf_version())
  source ("helpers.R")

  n <- 4
  init <- rep(1, n)
  niter <- 100
  tol <- 1e-8
  test_tol <- tol * 100

  fun <- function(state, iter) {

    # make fecundity a Ricker-like function of the total population, by
    # pro-rating down the fecundity
    Nt <- sum(state)
    K <- 100
    ratio <- exp(1 - Nt / K)
    multiplier <- 1 + (ratio - 1)
    state * multiplier

  }

  # r version
  r_iterates <- r_iterate_dynamic_function(
    transition_function = fun,
    initial_state = init,
    niter = niter,
    tol = tol
  )

  target_stable <- r_iterates$stable_state
  target_states <- r_iterates$all_states

  # greta version
  iterates <- iterate_dynamic_function(
    transition_function = fun,
    initial_state = init,
    niter = niter,
    tol = tol
  )

  stable <- iterates$stable_population
  states <- iterates$all_states
  converged <- iterates$converged
  iterations <- iterates$iterations

  greta_stable <- calculate(stable)[[1]]
  difference <- abs(greta_stable - target_stable)
  expect_true(all(difference < test_tol))

  greta_states <- calculate(states)[[1]]
  difference <- abs(greta_states - target_states)
  expect_true(all(difference < test_tol))

  greta_converged <- calculate(converged)[[1]]
  expect_true(greta_converged == 1)

  greta_iterations <- calculate(iterations)[[1]]
  expect_lt(greta_iterations, niter)

})



test_that("iteration works with time-varying parameters", {

  skip_if_not(greta:::check_tf_version())
  source ("helpers.R")

  n <- 4
  init <- runif(n)
  niter <- 100
  tol <- 0
  test_tol <- 1e-06

  # time-varying covariate
  x <- matrix(rnorm(niter * n), niter, n)

  fun <- function(state, iter, x) {

    # make fecundity a Ricker-like function of the total population, with random
    # fluctuations on each state
    Nt <- sum(state)
    K <- 100
    ratio <- exp(1 - Nt / K)
    multiplier <- 1 + (ratio - 1)
    state * multiplier * exp(x)

  }

  # r version
  r_iterates <- r_iterate_dynamic_function(
    transition_function = fun,
    initial_state = init,
    niter = niter,
    tol = tol,
    x = x,
    parameter_is_time_varying = "x"
  )

  target_states <- r_iterates$all_states

  # greta version
  iterates <- iterate_dynamic_function(
    transition_function = fun,
    initial_state = init,
    niter = niter,
    tol = tol,
    x = x,
    parameter_is_time_varying = "x"
  )

  states <- iterates$all_states

  greta_states <- calculate(states)[[1]]
  difference <- abs(greta_states - target_states)
  expect_true(all(difference < test_tol))

})
