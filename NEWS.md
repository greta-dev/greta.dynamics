# greta.dynamics (development version)

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
