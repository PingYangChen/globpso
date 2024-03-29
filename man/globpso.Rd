% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/globpso.R
\name{globpso}
\alias{globpso}
\title{Particle Swarm Optimization Algorithms for Minimization Problems}
\usage{
globpso(
  objFunc,
  lower,
  upper,
  init = NULL,
  fixed = NULL,
  PSO_INFO = NULL,
  seed = NULL,
  verbose = TRUE,
  environment,
  ...
)
}
\arguments{
\item{objFunc}{The R or Rcpp compiled objective function. See the example.}

\item{lower}{The vector of finite lower bounds of the search domain.
Its length should be equal to the dimension of the domain space.}

\item{upper}{The vector of finite upper bounds of the search domain.
Its length should be equal to the dimension of the domain space.}

\item{init}{The vector of initial swarm.
Its length should be equal to the dimension of the domain space.
When there are more than one initial vectors, specify \code{init} as
a matrix. Each row vector represents one initial point.
The default for \code{init} is \code{NULL}.}

\item{fixed}{The vector of real values and NA values that controls PSO to search only 
for the NA-valued components.}

\item{PSO_INFO}{The list of PSO parameters generated by \code{getPSOInfo()}.}

\item{seed}{The random seed that controls initial swarm of PSO. The default is \code{NULL}.}

\item{verbose}{The logical value controls if PSO would reports the updating progress. The default is \code{TRUE}.}
}
\value{
An List.
\itemize{
\item{par}{ the global best particle.}
\item{val}{ the objective function value of the global best particle.}
\item{history}{ a vector of objective function values of the global best particle in PSO search history.}
\item{cputime}{ the computational time in seconds.}
}
}
\description{
Particle Swarm Optimization Algorithms for Minimization Problems
}
\details{
TBD
}
\examples{
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
res <- globpso(objFunc = objf, lower = low_bound, upper = upp_bound, loc = loc_shift)
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
res_c <- globpso(objFunc = objf_c, lower = low_bound, upper = upp_bound, loc = loc_shift)
res_c$par
res_c$val
# Use getPSOInfo() to change the PSO options
alg_setting <- getPSOInfo(nSwarm = 64, maxIter = 200, psoType = "quantum")
res_c_large <- globpso(objFunc = objf_c, lower = low_bound, upper = upp_bound, PSO_INFO = alg_setting, loc = loc_shift)
res_c_large$history
}
\references{
Shi, Y., & Eberhart, R. (1998, May). A modified particle swarm optimizer. In Evolutionary Computation Proceedings, 1998. IEEE World Congress on Computational Intelligence., The 1998 IEEE International Conference on (pp. 69-73). IEEE.

Sun, J., Wu, X., Palade, V., Fang, W., Lai, C.-H., and Xu, W. (2012). Convergence analysis and improvements of quantum-behaved particle swarm optimization. Information Sciences, 193:81-103.
}
