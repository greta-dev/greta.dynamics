# ODE solve example

``` r

library(greta.dynamics)
#> Loading required package: greta
#> 
#> Attaching package: 'greta'
#> The following objects are masked from 'package:stats':
#> 
#>     binomial, cov2cor, poisson
#> The following objects are masked from 'package:base':
#> 
#>     %*%, %o%, apply, backsolve, beta, chol2inv, colMeans, colSums,
#>     diag, eigen, forwardsolve, gamma, identity, outer, rowMeans,
#>     rowSums, sweep, tapply
greta_sitrep()
#> ℹ checking if python available
#> ✔ python (v3.11) available
#> 
#> ℹ checking if TensorFlow available
#> ✔ TensorFlow (v2.15.1) available
#> 
#> ℹ checking if TensorFlow Probability available
#> ✔ TensorFlow Probability (v0.23.0) available
#> 
#> ℹ greta conda environment: not used (managed (uv) environment active)
#> • backend: "managed (uv) environment"
#> • selected via: RETICULATE_PYTHON environment variable
#> ✖ `UV_OFFLINE`=1 is set but the uv cache is not yet populated, so the next
#> start may fail to resolve dependencies
#> ℹ See the installation vignette: `vignette(greta::installation)`.
#> ℹ All dependencies available; greta is ready to use!
```

This example is taken from the `ode_solve.R` helpfile, and run here to
demonstrate the outputs from this in documentation.

``` r

# replicate the Lotka-Volterra example from deSolve
library(deSolve)
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

pars <- c(
  rIng = 0.2, # /day, rate of ingestion
  rGrow = 1.0, # /day, growth rate of prey
  rMort = 0.2, # /day, mortality rate of predator
  assEff = 0.5, # -, assimilation efficiency
  K = 10
) # mmol/m3, carrying capacity

yini <- c(Prey = 1, Predator = 2)
times <- seq(0, 30, by = 1)
out <- ode(y = yini, 
           times = times, 
           func = LVmod, 
           parms = pars)

# simulate observations
jitter <- rnorm(2 * length(times), 0, 0.1)
y_obs <- out[, -1] + matrix(jitter, ncol = 2)
```

``` r

# fit a greta model to infer the parameters from this simulated data

# greta version of the function
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

# priors for the parameters
rIng <- uniform(0, 2) # /day, rate of ingestion
rGrow <- uniform(0, 3) # /day, growth rate of prey
rMort <- uniform(0, 1) # /day, mortality rate of predator
assEff <- uniform(0, 1) # -, assimilation efficiency
K <- uniform(0, 30) # mmol/m3, carrying capacity

# initial values and observation error
y0 <- uniform(0, 5, dim = c(1, 2))
obs_sd <- uniform(0, 1)

# solution to the ODE
y <- ode_solve(
  derivative = lotka_volterra, 
  y0 = y0, 
  times = times, 
  rIng, rGrow, rMort, assEff, K,
  method = "dp"
  )

# sampling statement/observation model
distribution(y_obs) <- normal(y, obs_sd)
```

``` r

# we can use greta to solve directly, for a fixed set of parameters (the true
# ones in this case)
values <- c(
  list(y0 = t(1:2)),
  as.list(pars)
)
vals <- calculate(y, values = values)[[1]]
plot(vals[, 1] ~ times, type = "l", ylim = range(vals, na.rm = TRUE))
lines(vals[, 2] ~ times, lty = 2)
points(y_obs[, 1] ~ times)
points(y_obs[, 2] ~ times, pch = 2)
```

![](ode-solve-example_files/figure-html/greta-solve-directly-1.png)

``` r

# We can also do inference on the parameters
# priors for the parameters
rIng <- uniform(0, 2) # /day, rate of ingestion
rGrow <- uniform(0, 3) # /day, growth rate of prey
rMort <- uniform(0, 1) # /day, mortality rate of predator
assEff <- uniform(0, 1) # -, assimilation efficiency
K <- uniform(0, 30) # mmol/m3, carrying capacity
obs_sd <- uniform(0, 1)

# build the model
m <- model(rIng, rGrow, rMort, assEff, K, obs_sd)

# compute MAP estimate
o <- opt(
  model = m,
  initial_values = initials(
    rIng = 0.2,
    rGrow = 1.0,
    rMort = 0.2,
    assEff = 0.5,
    K = 10.0
  )
  )
o
#> $par
#> $par$rIng
#> [1] 1
#> 
#> $par$rGrow
#> [1] 1.5
#> 
#> $par$rMort
#> [1] 0.5
#> 
#> $par$assEff
#> [1] 0.5
#> 
#> $par$K
#> [1] 15
#> 
#> $par$obs_sd
#> [1] 0.5
#> 
#> 
#> $value
#> [1] 8.317766
#> 
#> $iterations
#> [1] 5
#> 
#> $convergence
#> [1] 0
```
