# iterate dynamic transition matrices

Calculate the stable population size for a stage-structured dynamical
system, encoded by a transition matrix, the value of which changes at
each iteration, given by function of the previous state:
`state[t] = f(state[t-1]) %*% state[t-1]`.

## Usage

``` r
iterate_dynamic_matrix(
  matrix_function,
  initial_state,
  niter,
  tol,
  ...,
  state_limits = c(-Inf, Inf)
)
```

## Arguments

- matrix_function:

  a function taking in the previous population state and the current
  iteration (and possibly other greta arrays) and returning a transition
  matrix to use for this iteration. The first two arguments must be
  named 'state' and 'iter', the state vector and scalar iteration number
  respectively. The remaining parameters must be named arguments
  representing (temporally static) model parameters. Variables and
  distributions cannot be defined inside the function.

- initial_state:

  either a column vector (with m elements) or a 3D array (with
  dimensions n x m x 1) giving one or more initial states from which to
  iterate the matrix

- niter:

  a positive integer giving the maximum number of times to iterate the
  matrix

- tol:

  a scalar giving a numerical tolerance, below which the algorithm is
  determined to have converged to a stable population size in all stages

- ...:

  optional named arguments to `matrix_function`, giving greta arrays for
  additional parameters

- state_limits:

  a numeric vector of length 2 giving minimum and maximum values at
  which to clamp the values of state after each iteration to prevent
  numerical under/overflow; i.e. elements with values below the minimum
  (maximum) will be set to the minimum (maximum).

## Value

a named list with four greta arrays:

- `stable_population` a vector or matrix (with the same dimensions as
  `initial_state`) giving the state after the final iteration.

- `all_states` an n x m x niter matrix of the state values at each
  iteration. This will be 0 for all entries after `iterations`.

- `converged` an integer scalar indicating whether *all* the matrix
  iterations converged to a tolerance less than `tol` (1 if so, 0 if
  not) before the algorithm finished.

- `iterations` a scalar of the maximum number of iterations completed
  before the algorithm terminated. This should match `niter` if
  `converged` is `FALSE`

## Details

Because `iterate_matrix` iterates with a static transition matrix, it
converges to a stable *growth rate* and *relative* population sizes for
a dynamical system. `iterate_dynamic_matrix` instead uses a matrix which
changes at each iteration, and can be dependent on the population sizes
after the previous iteration, and the iteration number. Because this can
encode density-dependence, the dynamics can converge to *absolute*
population sizes. The convergence criterion is therefore based on growth
rates converging on 0.

As in `iterate_matrix`, the greta array returned by `matrix_function`
can either be a square matrix, or a 3D array representing (on the first
dimension) *n* different matrices. `initial_state` should be shaped
accordingly, as detailed in `iterate_matrix`.

To ensure the matrix is iterated for a specific number of iterations,
you can set that number as `niter`, and set `tol` to 0 or a negative
number to ensure that the iterations are not stopped early.

## Note

because greta vectorises across both MCMC chains and the calculation of
greta array values, the algorithm is run until all chains (or posterior
samples), sites and stages have converged to stable growth. So a single
value of both `converged` and `iterations` is returned, and the value of
this will always have the same value in an `mcmc.list` object. So
inspecting the MCMC trace of these parameters will only tell you whether
the iteration converged in *all* posterior samples, and the maximum
number of iterations required to do so across all these samples

## Examples

``` r
if (FALSE) { # \dontrun{
# a transition matrix that is constant across iterations. Its columns sum to
# one, so the total population size is conserved and the two stages converge
# on a stable ratio of 2:1
matrix_function <- function(state, iter) {
  mat <- zeros(2, 2)
  mat[1, 1] <- 0.9
  mat[1, 2] <- 0.2
  mat[2, 1] <- 0.1
  mat[2, 2] <- 0.8
  mat
}

# start away from the stable state, so the iteration has work to do
initial_state <- as_data(matrix(c(2, 20), nrow = 2, ncol = 1))

results <- iterate_dynamic_matrix(
  matrix_function = matrix_function,
  initial_state = initial_state,
  niter = 100,
  tol = 1e-6
)

# the population size after the final iteration, and whether the iteration
# converged before hitting `niter`
calculate(results$stable_population, results$converged, results$iterations)
} # }
```
