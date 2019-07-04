# package file

#' greta.dynamics: a greta extension for modelling dynamical systems
#' @name greta.dynamics
#'
#' @description an extension to \href{https://greta-dev.github.io/greta}{greta}
#'   with functions for simulating dynamical systems, defined by of ordinary
#'   differential equations (see \code{\link{ode_solve}}) or transition matrices
#'   (\code{\link{iterate_matrix}}).
#'
#' @importFrom greta .internals
#' @importFrom tensorflow tf shape
#'
#' @docType package
NULL

# crate the node list object whenever the package is loaded
.onLoad <- function (libname, pkgname) {

  # add a numerical message for the ODE solver to the stash
  greta_stash <- greta::.internals$greta_stash
  greta_stash$numerical_messages <- c(greta_stash$numerical_messages,
                                      "underflow in dt",
                                      "non-finite values in state",
                                      "max_num_steps exceeded")

}

