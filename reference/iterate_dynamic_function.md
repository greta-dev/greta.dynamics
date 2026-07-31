# iterate dynamic transition functions

Calculate the stable population size for a stage-structured dynamical
system, encoded by a transition function, the value of which changes at
each iteration, given by function of the previous state:
`state[t] = f(state[t-1])`.

## Usage

``` r
iterate_dynamic_function(
  transition_function,
  initial_state,
  niter,
  tol,
  ...,
  parameter_is_time_varying = c(),
  state_limits = c(-Inf, Inf)
)
```

## Arguments

- transition_function:

  a function taking in the previous population state and the current
  iteration (and possibly other greta arrays) and returning the
  population state at the next iteration. The first two arguments must
  be named 'state' and 'iter', the state vector and scalar iteration
  number respectively. The remaining parameters must be named arguments
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

- parameter_is_time_varying:

  a character vector naming the parameters (ie. the named arguments of
  the function that are passed via `...`) that should be considered to
  be time-varying. That is, at each iteration only the corresponding
  slice from the first dimension of the object passed in should be used
  at that iteration.

- state_limits:

  a numeric vector of length 2 giving minimum and maximum values at
  which to clamp the values of state after each iteration to prevent
  numerical under/overflow; i.e. elements with values below the minimum
  (maximum) will be set to the minimum (maximum).

## Value

a named list with four greta arrays:

- `stable_state` a vector or matrix (with the same dimensions as
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

Like `iterate_dynamic_matrix` this converges to *absolute* population
sizes. The convergence criterion is therefore based on growth rates
converging on 0.

The greta array returned by `transition_function` must have the same
dimension as the state input and `initial_state` should be shaped
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
