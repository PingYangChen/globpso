globpso
=======
This package collects the powerful variants of particle swarm optimization algorithm.  Currently it supports the basic PSO algorithm and the Quantum-behaved PSO algorithm.  

Installation
------------
Please install the latest development version from github with

``` r
install.packages("devtools")
devtools::install_github("PingYangChen/globpso")
```

If you encounter a bug, please file a reproducible example on [github](https://github.com/PingYangChen/globpso/issues).

Examples
--------
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