globpso
=======

Features
--------
This package collects the powerful variants of particle swarm optimization (PSO) and differential evolution (DE) algorithms for R users. 

Currently it supports five types of PSO:
- Particle Swarm Optimization (PSO, Eberhart and Kennedy, 1995)
- Quantum-behaved particle Swarm Optimization (QPSO, Sun et al., 2004a,b)
- Locally convergent rotationally invariant particle
swarm optimization (LcRiPSO, Bonyadi and Michalewicz, 2014)
- Competitive Swarm Optimizer (CSO, Cheng and Jin, 2015)
- Double exponential particle swarm optimization (DExPSO, Stehlik et al., 2024+)

and six types of DE (Storn, R. and Price, K., 1997):
- DE/rand/1: Mutation operation on the current position with one random direction.
- DE/rand/2: Mutation operation on the current position with two random directions.
- DE/best/1: Mutation operation on the best position with one random direction.
- DE/best/2: Mutation operation on the best position with two random directions.
- DE/rand_to-best/1: Mutation operation on the current position with direction to the best and one random direction.
- DE/rand-to-best/2: Mutation operation on the current position with direction to the best and two random directions.


Installation
------------
Install the stable version of `globpso` on CRAN via

``` r
install.packages("globpso")
```

(Not recommended) Or, install the latest under-development version from github with
``` r
install.packages("devtools")
devtools::install_github("PingYangChen/globpso")
```

Examples of PSO
---------------
For setting the PSO type, please find the instruction through the command.
``` r
?getPSOInfo
```

``` r
library(globpso)

# Optimize the 3-dimensional quadratic objective function with location shift
objf <- function(x, loc) {
 val <- 0
 for (i in 1:length(x)) val <- val + (x[i] - loc)^2
 return(val)
}

# The search domain is [-5, 5]^3
upp_bound <- rep(5, 3)
low_bound <- rep(-5, 3)

# Define the location shift to be 1
loc_shift <- 1

# Run PSO for this optimization problem
# Also input the enviorment variable, the location shift 'loc_shift'
res <- globpso(objFunc = objf, lower = low_bound, upper = upp_bound, 
               loc = loc_shift)
res$par
res$val

# Set initial values to run PSO
initSwarm <- c(2, 2, 2)
# initSwarm <- rbind(c(2, 2, 2), c(1, 0.1, 1)) # for more than 2 initial points
res_i <- globpso(objFunc = objf, lower = low_bound, upper = upp_bound, 
                 init = initSwarm, loc = loc_shift)

# One can also write C++ objective function to further accelerate the computation
library(Rcpp)
library(RcppArmadillo)
objf_c <- cppFunction('double objf_c(SEXP x, SEXP loc) {
   double val = 0;
   double loc_c = (double)Rcpp::as<double>(loc);
   arma::rowvec x_c = (arma::rowvec)Rcpp::as<arma::rowvec>(x);
   for (arma::uword i = 0; i < x_c.n_elem; i++) {
     val += (x_c(i) - loc_c)*(x_c(i) - loc_c);
   }
   return val;
 }', depends = "RcppArmadillo")
res_c <- globpso(objFunc = objf_c, lower = low_bound, upper = upp_bound, 
                 loc = loc_shift)
res_c$par
res_c$val

# Use getPSOInfo() to change the PSO options
alg_setting <- getPSOInfo(nSwarm = 64, maxIter = 200, psoType = "quantum")
res_c_large <- globpso(objFunc = objf_c, lower = low_bound, upper = upp_bound, PSO_INFO = alg_setting, loc = loc_shift)
res_c_large$history
```

Examples of DE
--------------
For setting the DE type, please find the instruction through the command.
``` r
?getDEInfo
```

``` r
library(globpso)

# Optimize the 3-dimensional quadratic objective function with location shift
objf <- function(x, loc) {
  val <- 0
  for (i in 1:length(x)) val <- val + (x[i] - loc)^2
  return(val)
}

# The search domain is [-5, 5]^3
upp_bound <- rep(5, 3)
low_bound <- rep(-5, 3)

# Define the location shift to be 1
loc_shift <- 1

# Use getDEInfo() to setup the DE options
de_setting <- getDEInfo(nPop = 64, maxIter = 200, deType = "best-2", sf = 0.5, cr = 0.1)

# Run DE for this optimization problem
# Also input the enviorment variable, the location shift 'loc_shift'
res <- diffevo(objFunc = objf, lower = low_bound, upper = upp_bound, 
               DE_INFO = de_setting, verbose = FALSE, loc = loc_shift)
res$par
res$val

# Simialr to PSO, 
# One can also write C++ objective function to further accelerate the computation
```

Reference
---------
- Bonyadi, M. R. and Michalewicz, Z. (2014). A locally convergent rotationally invariant particle swarm optimization algorithm. Swarm Intelligence, 8(3):159-198.
- Cheng, R. and Jin, Y. (2015). A competitive swarm optimizer for large scale optimization. IEEE transactions on cybernetics, 45(2):191-204.
- Eberhart, R. and Kennedy, J. (1995). A new optimizer using particle swarm theory. In The 6th International Symposium on Micro Machine and Human Science, pages 39-43. IEEE.
- Stehlík, M., Chen, P. Y., Wong, W. K., and Kiseľák, J. (2024). A double exponential particle swarm optimization with non-uniform variates as stochastic tuning and guaranteed convergence to a global optimum with sample applications to finding optimal exact designs in biostatistics. Applied Soft Computing, 163, 111913.
- Storn, R., & Price, K. (1997). Differential evolution-a simple and efficient heuristic for global optimization over continuous spaces. Journal of global optimization, 11, 341-359.
- Sun, J., Feng, B., and Xu, W. (2004a). Particle swarm optimization with particles having quantum behavior. In Evolutionary Computation, 2004. CEC2004. Congress on, volume 1, pages 325-331. IEEE.
- Sun, J., Xu, W., and Feng, B. (2004b). A global search strategy of quantum-behaved particle swarm optimization. In Cybernetics and Intelligent Systems, 2004 IEEE Conference on, volume 1, pages 111-116. IEEE.


If you encounter a bug, please file a reproducible example on [github](https://github.com/PingYangChen/globpso/issues).
