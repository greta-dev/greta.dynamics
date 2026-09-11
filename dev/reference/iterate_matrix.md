# iterate transition matrices

Calculate the intrinsic growth rate(s) and stable stage distribution(s)
for a stage-structured dynamical system, encoded as
`state_t = matrix \%*\% state_tm1`.

## Usage

``` r
iterate_matrix(
  matrix,
  initial_state = rep(1, ncol(matrix)),
  niter = 100,
  tol = 1e-06
)
```

## Arguments

- matrix:

  either a square 2D transition matrix (with dimensions m x m), or a 3D
  array (with dimensions n x m x m), giving one or more transition
  matrices to iterate

- initial_state:

  either a column vector (with m elements) or a 3D array (with
  dimensions n x m x 1) giving one or more initial states from which to
  iterate the matrix

- niter:

  a positive integer giving the maximum number of times to iterate the
  matrix

- tol:

  a scalar giving a numerical tolerance, below which the algorithm is
  determined to have converged to the same growth rate in all stages

## Value

a named list with five greta arrays:

- `lambda` a scalar or vector giving the ratio of the first stage values
  between the final two iterations.

- `stable_distribution` a vector or matrix (with the same dimensions as
  `initial_state`) giving the state after the final iteration,
  normalised so that the values for all stages sum to one.

- `all_states` an n x m x niter matrix of the state values at each
  iteration. This will be 0 for all entries after `iterations`.

- `converged` an integer scalar or vector indicating whether the
  iterations for each matrix have converged to a tolerance less than
  `tol` (1 if so, 0 if not) before the algorithm finished.

- `iterations` a scalar of the maximum number of iterations completed
  before the algorithm terminated. This should match `niter` if
  `converged` is `FALSE`.

## Details

`iterate_matrix` can either act on a single transition matrix and
initial state (if `matrix` is 2D and `initial_state` is a column
vector), or it can simultaneously act on *n* different matrices and/or
*n* different initial states (if `matrix` and `initial_state` are 3D
arrays). In the latter case, the first dimension of both objects should
be the batch dimension *n*.

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
# simulate from a probabilistic 4-stage transition matrix model
k <- 4

# component variables
# survival probability for all stages
survival <- uniform(0, 1, dim = k)
# conditional (on survival) probability of staying in a stage
stasis <- c(uniform(0, 1, dim = k - 1), 1)
# marginal probability of staying/progressing
stay <- survival * stasis
progress <- (survival * (1 - stay))[1:(k - 1)]
# recruitment rate for the largest two stages
recruit <- exponential(c(3, 5))

# combine into a matrix:
tmat <- zeros(k, k)
diag(tmat) <- stay
progress_idx <- row(tmat) - col(tmat) == 1
tmat[progress_idx] <- progress
tmat[1, k - (1:0)] <- recruit

# analyse this to get the intrinsic growth rate and stable state
iterations <- iterate_matrix(tmat)
iterations$lambda
iterations$stable_distribution
iterations$all_states

# Can also do this simultaneously for a collection of transition matrices
k <- 2
n <- 10
survival <- uniform(0, 1, dim = c(n, k))
stasis <- cbind(uniform(0, 1, dim = n), rep(1, n))
stay <- survival * stasis
progress <- (survival * (1 - stasis))[, 1]
recruit_rate <- 1 / seq(0.1, 5, length.out = n)
recruit <- exponential(recruit_rate, dim = n)
tmats <- zeros(10, 2, 2)
tmats[, 1, 1] <- stasis[, 1]
tmats[, 2, 2] <- stasis[, 2]
tmats[, 2, 1] <- progress
tmats[, 1, 2] <- recruit

iterations <- iterate_matrix(tmats)
iterations$lambda
iterations$stable_distribution
iterations$all_states
} # }
```
