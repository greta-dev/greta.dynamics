#' greta.dynamics: a greta extension for modelling dynamical systems
#' @name greta.dynamics
#'
#' @description an extension to [greta](https://greta-stats.org/)
#'   with functions for simulating dynamical systems, defined by of ordinary
#'   differential equations (see [ode_solve()]) or transition matrices
#'   ([iterate_matrix()]).
#'
#' @importFrom greta .internals abind
#' @importFrom tensorflow tf shape
#'
"_PACKAGE"


# crate the node list object whenever the package is loaded
.onLoad <- function(libname, pkgname) {

  # add a numerical message for the ODE solver to the stash
  greta_stash <- greta::.internals$greta_stash
  greta_stash$numerical_messages <- c(
    greta_stash$numerical_messages,
    "underflow in dt",
    "non-finite values in state",
    "max_num_steps exceeded"
  )
}


## usethis namespace: start
## usethis namespace: end
NULL

globalVariables(
  c(
    "as_data"
  )
)
