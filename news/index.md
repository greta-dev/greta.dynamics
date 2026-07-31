# Changelog

## greta.dynamics 0.2.3

CRAN release: 2026-07-22

This release restores `greta.dynamics` to CRAN. It was archived in
September 2025 as a consequence of its dependency, `greta`, being
archived. `greta` 0.6.0 returned to CRAN in July 2026, which unblocks
this resubmission.

### Documentation fixes

- The `@return` documentation for
  [`iterate_matrix()`](https://greta-dev.github.io/greta.dynamics/reference/iterate_matrix.md),
  [`iterate_dynamic_matrix()`](https://greta-dev.github.io/greta.dynamics/reference/iterate_dynamic_matrix.md)
  and
  [`iterate_dynamic_function()`](https://greta-dev.github.io/greta.dynamics/reference/iterate_dynamic_function.md)
  named a list element `stable_state`, which none of these functions
  actually return. The documentation now uses the real names:
  `stable_distribution` for
  [`iterate_matrix()`](https://greta-dev.github.io/greta.dynamics/reference/iterate_matrix.md),
  and `stable_population` for
  [`iterate_dynamic_matrix()`](https://greta-dev.github.io/greta.dynamics/reference/iterate_dynamic_matrix.md)
  and
  [`iterate_dynamic_function()`](https://greta-dev.github.io/greta.dynamics/reference/iterate_dynamic_function.md).
- Added worked examples to
  [`iterate_dynamic_matrix()`](https://greta-dev.github.io/greta.dynamics/reference/iterate_dynamic_matrix.md)
  and
  [`iterate_dynamic_function()`](https://greta-dev.github.io/greta.dynamics/reference/iterate_dynamic_function.md),
  which previously had none.

### Internal changes

- Now requires `greta` \>= 0.6.0.
- Removed `LazyData` from `DESCRIPTION`, as the package ships no data.

## greta.dynamics 0.2.2

CRAN release: 2024-11-14

### New features

- Added
  [`iterate_dynamic_function()`](https://greta-dev.github.io/greta.dynamics/reference/iterate_dynamic_function.md)
  and
  [`iterate_dynamic_matrix()`](https://greta-dev.github.io/greta.dynamics/reference/iterate_dynamic_matrix.md)

### Breaking Changes

- Updated internal ODE interface to match new tensorflow probability
  API. This involves now only supporting two solvers, “dp”, and “bdf”.
  The Default is “dp”, which is similar to deSolve’s “ode45”. The “dp”
  solver is Dormand-Prince explicit solver for non-stiff ODEs. The “bdf”
  solver is Backward Differentiation Formula (BDF) solver for stiff
  ODEs. Currently no arguments for “bdf” or “dp” are able to be
  specified.

### Internal changes

- import rlang, cli, use latest version of greta, 0.5.0
- use internal checking functions
- use snapshot testing for checking error messages

## greta.dynamics 0.2.1

- Use sentinel “\_PACKAGE”

## greta.dynamics 0.2.0

CRAN release: 2022-09-05

- Added a `NEWS.md` file to track changes to the package.
