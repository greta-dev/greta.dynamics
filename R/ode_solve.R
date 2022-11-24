#' @name ode_solve
#'
#' @title solve ODEs
#'
#' @description Solve a system of ordinary differential equations.
#'
#' @param derivative a derivative function. The first two arguments must be 'y'
#'   and 't', the state parameter and scalar timestep respectively. The
#'   remaining parameters must be named arguments representing (temporally
#'   static) model parameters. Variables and distributions cannot be defined in
#'   the function.
#' @param y0 a greta array for the value of the state parameter y at time 0
#' @param times a column vector of times at which to evaluate y
#' @param ... named arguments giving greta arrays for the additional (fixed)
#'   parameters
#' @param method which solver to use. `"ode45"` uses adaptive step
#'   sizes, whilst `"rk4"` and `"midpoint"` use the fixed grid defined
#'   by `times`; they may be faster but less accurate than `"ode45"`.
#'
#' @return greta array
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # replicate the Lotka-Volterra example from deSolve
#' library(deSolve)
#' LVmod <- function(Time, State, Pars) {
#'   with(as.list(c(State, Pars)), {
#'     Ingestion <- rIng * Prey * Predator
#'     GrowthPrey <- rGrow * Prey * (1 - Prey / K)
#'     MortPredator <- rMort * Predator
#'
#'     dPrey <- GrowthPrey - Ingestion
#'     dPredator <- Ingestion * assEff - MortPredator
#'
#'     return(list(c(dPrey, dPredator)))
#'   })
#' }
#'
#' pars <- c(
#'   rIng = 0.2, # /day, rate of ingestion
#'   rGrow = 1.0, # /day, growth rate of prey
#'   rMort = 0.2, # /day, mortality rate of predator
#'   assEff = 0.5, # -, assimilation efficiency
#'   K = 10
#' ) # mmol/m3, carrying capacity
#'
#' yini <- c(Prey = 1, Predator = 2)
#' times <- seq(0, 30, by = 1)
#' out <- ode(yini, times, LVmod, pars)
#'
#' # simulate observations
#' jitter <- rnorm(2 * length(times), 0, 0.1)
#' y_obs <- out[, -1] + matrix(jitter, ncol = 2)
#'
#' # ~~~~~~~~~
#' # fit a greta model to infer the parameters from this simulated data
#'
#' # greta version of the function
#' lotka_volterra <- function(y, t, rIng, rGrow, rMort, assEff, K) {
#'   Prey <- y[1, 1]
#'   Predator <- y[1, 2]
#'
#'   Ingestion <- rIng * Prey * Predator
#'   GrowthPrey <- rGrow * Prey * (1 - Prey / K)
#'   MortPredator <- rMort * Predator
#'
#'   dPrey <- GrowthPrey - Ingestion
#'   dPredator <- Ingestion * assEff - MortPredator
#'
#'   cbind(dPrey, dPredator)
#' }
#'
#' # priors for the parameters
#' rIng <- uniform(0, 2) # /day, rate of ingestion
#' rGrow <- uniform(0, 3) # /day, growth rate of prey
#' rMort <- uniform(0, 1) # /day, mortality rate of predator
#' assEff <- uniform(0, 1) # -, assimilation efficiency
#' K <- uniform(0, 30) # mmol/m3, carrying capacity
#'
#' # initial values and observation error
#' y0 <- uniform(0, 5, dim = c(1, 2))
#' obs_sd <- uniform(0, 1)
#'
#' # solution to the ODE
#' y <- ode_solve(lotka_volterra, y0, times, rIng, rGrow, rMort, assEff, K)
#'
#' # sampling statement/observation model
#' distribution(y_obs) <- normal(y, obs_sd)
#'
#' # we can use greta to solve directly, for a fixed set of parameters (the true
#' # ones in this case)
#' values <- c(
#'   list(y0 = t(1:2)),
#'   as.list(pars)
#' )
#' vals <- calculate(y, values = values)[[1]]
#' plot(vals[, 1] ~ times, type = "l", ylim = range(vals))
#' lines(vals[, 2] ~ times, lty = 2)
#' points(y_obs[, 1] ~ times)
#' points(y_obs[, 2] ~ times, pch = 2)
#'
#' # or we can do inference on the parameters:
#'
#' # build the model (takes a few seconds to define the tensorflow graph)
#' m <- model(rIng, rGrow, rMort, assEff, K, obs_sd)
#'
#' # compute MAP estimate
#' o <- opt(m)
#' o
#' }
ode_solve <- function(derivative, y0, times, ...,
                      method = c("ode45", "rk4", "midpoint")) {
  y0 <- as.greta_array(y0)
  times <- as.greta_array(times)
  method <- match.arg(method)

  # check times is a column vector
  t_dim <- dim(times)
  if (length(t_dim != 2) && t_dim[2] != 1) {
    stop("",
      call. = FALSE
    )
  }

  dots <- list(...)
  dots <- lapply(dots, as.greta_array)
  # check all additional parameters are named and there are no extras

  # check derivative is a function

  # check the arguments of derivative are valid and match dots

  # create a tensorflow version of the function
  tf_derivative <- as_tf_derivative(derivative, y0, times[1], dots)

  # the dimensions should be the dimensions of the y0, duplicated along times
  n_time <- dim(times)[1]
  y0_dims <- dim(y0)

  # drop the first element if it's a one
  if (y0_dims[1] == 1) {
    y0_dims <- y0_dims[-1]
  }

  dims <- c(n_time, y0_dims)

  op("ode", y0, times, ...,
    dim = dims,
    tf_operation = "tf_ode_solve",
    operation_args = list(
      tf_derivative = tf_derivative,
      method = method
    )
  )
}

# internal tf function wrapping the core TF method
# return a tensor for the integral of derivative evaluated at times, given
# starting state y0 and other parameters dots
tf_ode_solve <- function(y0, times, ..., tf_derivative, method) {

  # drop the columns and batch dimension in times
  times <- tf_flatten(times)
  times <- tf$slice(times, c(0L, 0L), c(1L, -1L))
  times <- tf$squeeze(times, 0L)

  # assign the dots (as tensors) to the function environment
  assign(
    "tf_dots", list(...),
    environment(tf_derivative)
  )

  # integrate - need to run this with the dag's TF graph as default, so that
  # tf_derivative creates tensors correctly
  dag <- parent.frame()$dag

  tf_int <- tf$contrib$integrate

  integrator <- switch(method,
    ode45 = tf_int$odeint,
    rk4 = function(...) {
      tf_int$odeint_fixed(..., method = "rk4")
    },
    midpoint = function(...) {
      tf_int$odeint_fixed(..., method = "midpoint")
    }
  )

  dag$on_graph(integral <- integrator(tf_derivative, y0, times))

  # reshape to put batch dimension first
  permutation <- seq_along(dim(integral)) - 1L
  permutation[1:2] <- permutation[2:1]
  integral <- tf$transpose(integral, perm = permutation)

  # if the first (non-batch) dimension of y0 was 1, drop it in the results
  if (dim(y0)[[2]] == 1) {
    integral <- tf$squeeze(integral, 2L)
  }

  integral
}

# given a greta/R function derivative function, and greta arrays for the inputs,
# return a tensorflow function taking tensors for y and t and returning a tensor
# for dydt
as_tf_derivative <- function(derivative, y, t, dots) {

  # create a function acting on the full set of inputs, as tensors
  args <- list(r_fun = derivative, y = y, t = t)
  tf_fun <- do.call(as_tf_function, c(args, dots))

  # for CRAN's benefit
  tf_dots <- NULL

  # return a function acting only on tensors y and t, to feed to the ode solver
  function(y, t) {

    # t will be dimensionless when used in the ode solver, need to expand out t
    # to have same dim as a scalar constant so that it can be used in the same
    # way as the greta array in the R function
    t <- tf$reshape(t, shape = shape(1, 1, 1))

    # tf_dots will have been added to this environment by tf_ode_solve
    args <- list(y = y, t = t)
    do.call(tf_fun, c(args, tf_dots))
  }
}
