context('iteratation functions')

test_that('single iteration works', {

  skip_if_not(greta:::check_tf_version())
  source ('helpers.R')

  n <- 10
  mat <- randu(n, n)
  init <- rep(1, n)
  niter <- 50

  # r version
  r_iterations <- r_iterate_matrix(matrix = mat,
                             state = init,
                             niter = niter)

  target_lambda <- r_iterations$lambda
  target_state <- r_iterations$final_state
  target_state <- target_state / sum(target_state)

  # greta version
  iterations <- iterate_matrix(matrix = mat,
                               initial_state = init,
                               niter = niter)

  lambda <- iterations$lambda
  final_state <- iterations$final_state
  final_state <- final_state / sum(final_state)

  greta_lambda <- calculate(lambda)
  difference <- abs(greta_lambda - target_lambda)
  expect_true(difference < 1e-12)

  greta_state <- calculate(final_state)
  difference <- abs(greta_state - target_state)
  expect_true(all(difference < 1e-12))

})

test_that('vectorised matrix iteration works', {

  skip_if_not(greta:::check_tf_version())
  source('helpers.R')

  n <- 10
  n_mat <- 20
  mat <- randu(n_mat, n, n)
  init <- rep(1, n)
  niter <- 50
  mat_list <- lapply(seq_len(n_mat), function (i) mat[i, , ])

  # r version
  target_iterations <- lapply(mat_list,
                              r_iterate_matrix,
                              state = init,
                              niter = niter)

  target_lambdas <- sapply(target_iterations, function(x) x$lambda)
  target_states <- t(sapply(target_iterations, function(x) x$final_state))

  iterations <- iterate_matrix(mat,
                               initial_state = init,
                               niter = niter)
  greta_lambdas <- calculate(iterations$lambda)
  greta_states <- calculate(iterations$final_state)
  dim(greta_states) <- dim(greta_states)[1:2]

  difference <- abs(greta_lambdas - target_lambdas)
  expect_true(all(difference < 1e-12))

  difference <- abs(greta_states - target_states)
  expect_true(all(difference < 1e-12))

})

test_that('vectorised initial_state iteration works', {

  skip_if_not(greta:::check_tf_version())
  source('helpers.R')

  n <- 10
  n_mat <- 20
  mat <- randu(n, n)
  init <- randu(n_mat, n, 1)
  niter <- 50
  init_list <- lapply(seq_len(n_mat), function (i) init[i, , ])

  # r version
  target_iterations <- lapply(init_list,
                              function (init) {
                                r_iterate_matrix(mat, init, niter = niter)
                              })

  target_lambdas <- sapply(target_iterations, function(x) x$lambda)
  target_states <- t(sapply(target_iterations, function(x) x$final_state))

  iterations <- iterate_matrix(mat,
                               initial_state = init,
                               niter = niter)
  greta_lambdas <- calculate(iterations$lambda)
  greta_states <- calculate(iterations$final_state)
  dim(greta_states) <- dim(greta_states)[1:2]

  difference <- abs(greta_lambdas - target_lambdas)
  expect_true(all(difference < 1e-12))

  difference <- abs(greta_states - target_states)
  expect_true(all(difference < 1e-12))

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
