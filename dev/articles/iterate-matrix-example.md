# iterate-matrix-example

``` r

library(greta.dynamics)
```

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

``` r

greta_sitrep()
```

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

    #> ✖ will need internet on next start: uv cache not yet populated

    #> ℹ See the installation vignette: `vignette(greta::installation)`.
    #> ℹ All dependencies available; greta is ready to use!

``` r

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
```

``` r

# combine into a matrix:
tmat <- zeros(k, k)
diag(tmat) <- stay
progress_idx <- row(tmat) - col(tmat) == 1
tmat[progress_idx] <- progress
tmat[1, k - (1:0)] <- recruit
```

``` r

# analyse this to get the intrinsic growth rate and stable state
iterations <- iterate_matrix(tmat)
iterations$lambda
```

    #> greta array <operation>

    #> 

    #>      [,1]
    #> [1,]  ?

    #> 

``` r

iterations$stable_distribution
```

    #> greta array <operation>

    #>      [,1]
    #> [1,]  ?  
    #> [2,]  ?  
    #> [3,]  ?  
    #> [4,]  ?

    #> 

``` r

iterations$all_states
```

    #> greta array <operation>

    #>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
    #> [1,]  ?    ?    ?    ?    ?    ?    ?    ?    ?    ?     ?     ?     ?     ?   
    #> [2,]  ?    ?    ?    ?    ?    ?    ?    ?    ?    ?     ?     ?     ?     ?   
    #> [3,]  ?    ?    ?    ?    ?    ?    ?    ?    ?    ?     ?     ?     ?     ?   
    #> [4,]  ?    ?    ?    ?    ?    ?    ?    ?    ?    ?     ?     ?     ?     ?   
    #>      [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24] [,25] [,26]
    #> [1,]  ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?   
    #> [2,]  ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?   
    #> [3,]  ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?   
    #> [4,]  ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?   
    #>      [,27] [,28] [,29] [,30] [,31] [,32] [,33] [,34] [,35] [,36] [,37] [,38]
    #> [1,]  ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?   
    #> [2,]  ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?   
    #> [3,]  ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?   
    #> [4,]  ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?   
    #>      [,39] [,40] [,41] [,42] [,43] [,44] [,45] [,46] [,47] [,48] [,49] [,50]
    #> [1,]  ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?   
    #> [2,]  ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?   
    #> [3,]  ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?   
    #> [4,]  ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?   
    #>      [,51] [,52] [,53] [,54] [,55] [,56] [,57] [,58] [,59] [,60] [,61] [,62]
    #> [1,]  ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?   
    #> [2,]  ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?   
    #> [3,]  ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?   
    #> [4,]  ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?   
    #>      [,63] [,64] [,65] [,66] [,67] [,68] [,69] [,70] [,71] [,72] [,73] [,74]
    #> [1,]  ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?   
    #> [2,]  ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?   
    #> [3,]  ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?   
    #> [4,]  ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?   
    #>      [,75] [,76] [,77] [,78] [,79] [,80] [,81] [,82] [,83] [,84] [,85] [,86]
    #> [1,]  ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?   
    #> [2,]  ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?   
    #> [3,]  ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?   
    #> [4,]  ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?   
    #>      [,87] [,88] [,89] [,90] [,91] [,92] [,93] [,94] [,95] [,96] [,97] [,98]
    #> [1,]  ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?   
    #> [2,]  ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?   
    #> [3,]  ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?   
    #> [4,]  ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?     ?   
    #>      [,99] [,100]
    #> [1,]  ?     ?    
    #> [2,]  ?     ?    
    #> [3,]  ?     ?    
    #> [4,]  ?     ?

    #> 

    #> ℹ 390 more values
    #> Use `print(n = ...)` to see more values

``` r

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
```

``` r

iterations <- iterate_matrix(tmats)
iterations$lambda
```

    #> greta array <operation>

    #> 

    #>       [,1]
    #>  [1,]  ?  
    #>  [2,]  ?  
    #>  [3,]  ?  
    #>  [4,]  ?  
    #>  [5,]  ?  
    #>  [6,]  ?  
    #>  [7,]  ?  
    #>  [8,]  ?  
    #>  [9,]  ?  
    #> [10,]  ?

    #> 

``` r

iterations$stable_distribution
```

    #> greta array <operation>

    #> , , 1
    #> 
    #>       [,1] [,2]
    #>  [1,]  ?    ?  
    #>  [2,]  ?    ?  
    #>  [3,]  ?    ?  
    #>  [4,]  ?    ?  
    #>  [5,]  ?    ?  
    #>  [6,]  ?    ?  
    #>  [7,]  ?    ?  
    #>  [8,]  ?    ?  
    #>  [9,]  ?    ?  
    #> [10,]  ?    ?

    #> 

    #> ℹ 10 more values
    #> Use `print(n = ...)` to see more values

``` r

dim(iterations$all_states)
```

    #> [1]  10   2 100

``` r

iterations$all_states[, , 1]
```

    #> greta array <operation>

    #> 

    #> , , 1
    #> 
    #>       [,1] [,2]
    #>  [1,]  ?    ?  
    #>  [2,]  ?    ?  
    #>  [3,]  ?    ?  
    #>  [4,]  ?    ?  
    #>  [5,]  ?    ?  
    #>  [6,]  ?    ?  
    #>  [7,]  ?    ?  
    #>  [8,]  ?    ?  
    #>  [9,]  ?    ?  
    #> [10,]  ?    ?

    #> 

    #> ℹ 10 more values
    #> Use `print(n = ...)` to see more values

``` r

iterations$all_states[, , 100]
```

    #> greta array <operation>

    #> 

    #> , , 1
    #> 
    #>       [,1] [,2]
    #>  [1,]  ?    ?  
    #>  [2,]  ?    ?  
    #>  [3,]  ?    ?  
    #>  [4,]  ?    ?  
    #>  [5,]  ?    ?  
    #>  [6,]  ?    ?  
    #>  [7,]  ?    ?  
    #>  [8,]  ?    ?  
    #>  [9,]  ?    ?  
    #> [10,]  ?    ?

    #> 

    #> ℹ 10 more values
    #> Use `print(n = ...)` to see more values
