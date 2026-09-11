# greta.dynamics (development version)

## Bug fixes

- Resolved an issue where `iterate_dynamic_function()` and 
  `iterate_dynamic_matrix()` could not be used in a model sampled with `mcmc()`. 
  This was a result of greta using `.batch_size` over `batch_size` - see https://github.com/greta-dev/greta/issues/634. We now use `.batch_size`. Tests
  Added to catch this in the future.

# greta.dynamics 0.2.3

This release restores `greta.dynamics` to CRAN. It was archived in September
2025 as a consequence of its dependency, `greta`, being archived. `greta` 0.6.0
returned to CRAN in July 2026, which unblocks this resubmission.

## Documentation fixes

- The `@return` documentation for `iterate_matrix()`, `iterate_dynamic_matrix()`
  and `iterate_dynamic_function()` named a list element `stable_state`, which
  none of these functions actually return. The documentation now uses the real
  names: `stable_distribution` for `iterate_matrix()`, and `stable_population`
  for `iterate_dynamic_matrix()` and `iterate_dynamic_function()`.
- Added worked examples to `iterate_dynamic_matrix()` and
  `iterate_dynamic_function()`, which previously had none.

## Internal changes

- Now requires `greta` >= 0.6.0.
- Removed `LazyData` from `DESCRIPTION`, as the package ships no data.

# greta.dynamics 0.2.2

## New features

- Added `iterate_dynamic_function()` and `iterate_dynamic_matrix()`

## Breaking Changes

- Updated internal ODE interface to match new tensorflow probability API. This 
  involves now only supporting two solvers, "dp", and "bdf". The Default is 
  "dp", which is similar to deSolve's "ode45".  The "dp" solver is 
  Dormand-Prince explicit solver for non-stiff ODEs. The "bdf" solver is
  Backward Differentiation Formula (BDF) solver for stiff ODEs. Currently no
  arguments for "bdf" or "dp" are able to be specified.
  
## Internal changes

- import rlang, cli, use latest version of greta, 0.5.0
- use internal checking functions
- use snapshot testing for checking error messages

# greta.dynamics 0.2.1

- Use sentinel "_PACKAGE"

# greta.dynamics 0.2.0

- Added a `NEWS.md` file to track changes to the package.
