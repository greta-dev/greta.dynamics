context('iteration functions')

test_that('single iteration works', {

  skip_if_not(greta:::check_tf_version())
  source ('helpers.R')

  n <- 10
  mat <- randu(n, n)
  init <- rep(1, n)
  niter <- 50
  tol <- 1e-6
  test_tol <- tol * 10

  # r version
  r_iterations <- r_iterate_matrix(matrix = mat,
                             state = init,
                             niter = niter,
                             tol = tol)

  target_lambda <- r_iterations$lambda
  target_stable <- r_iterations$stable_distribution
  target_states <- r_iterations$all_states

  # greta version
  iterations <- iterate_matrix(matrix = mat,
                               initial_state = init,
                               niter = niter,
                               tol = tol)

  lambda <- iterations$lambda
  stable <- iterations$stable_distribution
  states <- iterations$all_states

  greta_lambda <- calculate(lambda)
  difference <- abs(greta_lambda - target_lambda)
  expect_true(difference < test_tol)

  greta_stable <- calculate(stable)
  difference <- abs(greta_stable - target_stable)
  expect_true(all(difference < test_tol))

  greta_states <- calculate(states)
  difference <- abs(greta_states - target_states)
  expect_true(all(difference < test_tol))

})

test_that('vectorised matrix iteration works', {

  skip_if_not(greta:::check_tf_version())
  source('helpers.R')

  n <- 10
  n_mat <- 20
  mat <- randu(n_mat, n, n)
  init <- rep(1, n)
  niter <- 50
  tol <- 1e-6
  test_tol <- tol * 10
  mat_list <- lapply(seq_len(n_mat), function (i) mat[i, , ])

  # r version
  target_iterations <- lapply(mat_list,
                              r_iterate_matrix,
                              state = init,
                              niter = niter,
                              tol = tol)

  target_lambdas <- sapply(target_iterations, function(x) x$lambda)
  target_stable <- t(sapply(target_iterations, function(x) x$stable_distribution))

  iterations <- iterate_matrix(mat,
                               initial_state = init,
                               niter = niter,
                               tol = tol)
  greta_lambdas <- calculate(iterations$lambda)
  greta_stable <- calculate(iterations$stable_distribution)
  dim(greta_stable) <- dim(greta_stable)[1:2]

  difference <- abs(greta_lambdas - target_lambdas)
  expect_true(all(difference < test_tol))

  difference <- abs(greta_stable - target_stable)
  expect_true(all(difference < test_tol))

})

test_that('vectorised initial_state iteration works', {

  skip_if_not(greta:::check_tf_version())
  source('helpers.R')

  n <- 10
  n_mat <- 20
  mat <- randu(n, n)
  init <- randu(n_mat, n, 1)
  niter <- 50
  tol <- 1e-6
  test_tol <- tol * 10
  init_list <- lapply(seq_len(n_mat), function (i) init[i, , ])

  # r version
  target_iterations <- lapply(init_list,
                              function (init) {
                                r_iterate_matrix(
                                  mat,
                                  init,
                                  niter = niter,
                                  tol = tol
                                )
                              })

  target_lambdas <- sapply(target_iterations, function(x) x$lambda)
  target_stable <- t(sapply(target_iterations, function(x) x$stable_distribution))

  iterations <- iterate_matrix(mat,
                               initial_state = init,
                               niter = niter,
                               tol = tol)
  greta_lambdas <- calculate(iterations$lambda)
  greta_stable <- calculate(iterations$stable_distribution)
  dim(greta_stable) <- dim(greta_stable)[1:2]

  difference <- abs(greta_lambdas - target_lambdas)
  expect_true(all(difference < test_tol))

  difference <- abs(greta_stable - target_stable)
  expect_true(all(difference < test_tol))

})

test_that('dynamics module errors informatively', {

  source ('helpers.R')

  n <- 10
  m <- 3

  bad_mat <- randu(m, m + 1)
  bad_state <- randu(m, 2)
  bad_matrices1 <- randu(n, m, m, 1)
  bad_matrices2 <- randu(n, m, m - 1)

  good_mat <- randu(m, m)
  good_state <- randu(m, 1)
  good_states <- randu(n, m, 1)
  good_matrices <- randu(n, m, m)

  mismatched_state <- randu(m + 1, 1)

  # wrongly shaped matrix
  expect_error(iterate_matrix(matrix = bad_mat,
                              initial_state = good_state),
               'matrix must be a two-dimensional square greta array')

  expect_error(iterate_matrix(matrix = bad_matrices1,
                              initial_state = good_state),
               '^matrix and state must be either two- or three-dimensional')

  expect_error(iterate_matrix(matrix = bad_matrices1,
                              initial_state = good_states),
               '^matrix and state must be either two- or three-dimensional')

  expect_error(iterate_matrix(matrix = bad_matrices2,
                              initial_state = good_state),
               '^each matrix must be a two-dimensional square greta array')

  expect_error(iterate_matrix(matrix = bad_matrices2,
                              initial_state = good_states),
               '^each matrix must be a two-dimensional square greta array')

  # wrongly shaped state
  expect_error(iterate_matrix(matrix = good_mat,
                              initial_state = bad_state),
               'initial_state must be either a column vector, or a 3D array')

  expect_error(iterate_matrix(matrix = good_matrices,
                              initial_state = bad_state),
               'initial_state must be either a column vector, or a 3D array')

  # mismatched matrix and state
  expect_error(iterate_matrix(matrix = good_mat,
                              initial_state = mismatched_state),
               'length of each initial_state must match the dimension')

  expect_error(iterate_matrix(matrix = good_matrices,
                              initial_state = mismatched_state),
               'length of each initial_state must match the dimension')

})

test_that('convergence tolerance works', {

  skip_if_not(greta:::check_tf_version())
  source ('helpers.R')

  n <- 10
  niter <- 100

  # with an identity matrix it should converge instantly
  iterations <- iterate_matrix(matrix = diag(n))
  converged <- calculate(iterations$converged)
  expect_true(converged, 0L)

  # with tolerance of 0, it should time out
  iterations <- iterate_matrix(matrix = randu(n, n), tol = 0)
  converged <- calculate(iterations$converged)
  expect_false(converged, 1L)

})

test_that("iteration works in mcmc", {

  skip_if_not(greta:::check_tf_version())
  source ('helpers.R')

  n <- 10
  n_site <- 30

  # non-batched case
  mat <- uniform(0, 1, dim = c(n, n))
  iterations <- iterate_matrix(matrix = mat)
  lambda <- iterations$lambda
  m <- model(lambda)
  draws <- mcmc(m, warmup = 100, n_samples = 100, verbose = FALSE)
  expect_s3_class(draws, "mcmc.list")

  # batched case
  mat <- uniform(0, 1, dim = c(n_site, n, n))
  iterations <- iterate_matrix(matrix = mat)
  lambda <- iterations$lambda
  m <- model(lambda)
  draws <- mcmc(m, warmup = 100, n_samples = 100, verbose = FALSE)
  expect_s3_class(draws, "mcmc.list")

})
