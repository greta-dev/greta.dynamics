test_that("ode_solve works like deSolve::ode", {
  set.seed(2017 - 05 - 01)
  skip_if_not(check_tf_version())

  # deSolve version of the Lotka Volterra model
  LVmod <- function(Time, State, Pars) {
    with(as.list(c(State, Pars)), {
      Ingestion <- rIng * Prey * Predator
      GrowthPrey <- rGrow * Prey * (1 - Prey / K)
      MortPredator <- rMort * Predator

      dPrey <- GrowthPrey - Ingestion
      dPredator <- Ingestion * assEff - MortPredator

      return(list(c(dPrey, dPredator)))
    })
  }

  # greta version of the model
  lotka_volterra <- function(y, t, rIng, rGrow, rMort, assEff, K) {
    Prey <- y[1, 1]
    Predator <- y[1, 2]

    Ingestion <- rIng * Prey * Predator
    GrowthPrey <- rGrow * Prey * (1 - Prey / K)
    MortPredator <- rMort * Predator

    dPrey <- GrowthPrey - Ingestion
    dPredator <- Ingestion * assEff - MortPredator

    cbind(dPrey, dPredator)
  }

  pars <- c(
    rIng = 0.2, # /day, rate of ingestion
    rGrow = 1.0, # /day, growth rate of prey
    rMort = 0.2, # /day, mortality rate of predator
    assEff = 0.5, # -, assimilation efficiency
    K = 10
  ) # mmol/m3, carrying capacity

  yini <- c(Prey = 1, Predator = 2)
  times <- seq(0, 50, by = 1)

  # loop through the solvers (ode45 should be similar to the dopri5 method in TF)
  methods <- c("bdf", "dp")

  desolve_ode_fun <- function(method){
    out <- deSolve::ode(yini, times, LVmod, pars, method = method)
    out
  }

  greta_ode_fun <- function(method){
    y <- ode_solve(lotka_volterra,
                   y0 = t(yini),
                   times,
                   rIng = pars["rIng"],
                   rGrow = pars["rGrow"],
                   rMort = pars["rMort"],
                   assEff = pars["assEff"],
                   K = pars["K"],
                   method = method
    )
    g_out <- cbind(times, y)

    greta_out <- calculate(g_out)[[1]]

    greta_out
  }

  desolve_bdf <- desolve_ode_fun("bdf")
  greta_bdf <- greta_ode_fun("bdf")
  # desolve equivalent to dp is ode45
  desolve_dp <- desolve_ode_fun("ode45")
  greta_dp <- greta_ode_fun("dp")

  difference_bdf <- abs(greta_bdf - desolve_bdf)
  difference_dp <- abs(greta_dp - desolve_dp)

  # these aren't a great match (during regions of rapid change), apparently due
  # to hard-coded differences in implementation between deSolve and TFP
  expect_true(all(difference_bdf < 1e-2))
  expect_true(all(difference_dp < 1e-2))

})

test_that("inference works with ode_solve", {
  skip_if_not(check_tf_version())
  set.seed(2017 - 05 - 01)

  lotka_volterra <- function(y, t, rIng, rGrow, rMort, assEff, K) {
    Prey <- y[1, 1]
    Predator <- y[1, 2]

    Ingestion <- rIng * Prey * Predator
    GrowthPrey <- rGrow * Prey * (1 - Prey / K)
    MortPredator <- rMort * Predator

    dPrey <- GrowthPrey - Ingestion
    dPredator <- Ingestion * assEff - MortPredator

    cbind(dPrey, dPredator)
  }

  rIng <- uniform(0, 2) # /day, rate of ingestion
  rGrow <- uniform(0, 3) # /day, growth rate of prey
  rMort <- uniform(0, 1) # /day, mortality rate of predator
  assEff <- uniform(0, 1) # -, assimilation efficiency
  K <- uniform(0, 30) # mmol/m3, carrying capacity

  yini <- c(Prey = 1, Predator = 2)
  times <- seq(0, 50, by = 1)

  y <- ode_solve(lotka_volterra,
    y0 = t(yini),
    times,
    rIng = rIng,
    rGrow = rGrow,
    rMort = rMort,
    assEff = assEff,
    K = K,
    method = "dp"
  )

  # simulate some data and fit to it
  y_true <- calculate(y, nsim = 1)[[1]][1, , ]
  y_obs <- y_true + rnorm(prod(dim(y_true)), 0, 0.1)
  distribution(y_obs) <- normal(y_true, 0.1)

  m <- model(rIng)

  # should be fine in opt() and mcmc()
  o <- opt(m)
  expect_true(is.list(o))
  expect_identical(names(o), c("par", "value", "iterations", "convergence"))

  draws <- mcmc(m, chains = 2, warmup = 100, n_samples = 100, verbose = FALSE)
  expect_s3_class(draws, "greta_mcmc_list")
})
