#' The globpso package: Particle Swarm Optimization Algorithms and Differential Evolution for Minimization Problems
#'
#' Particle Swarm Optimization Algorithms and Differential Evolution for Minimization Problems
#'
#' @references 
#' \enumerate{
#'   \item Bonyadi, M. R., & Michalewicz, Z. (2014). A locally convergent rotationally invariant particle swarm optimization algorithm. Swarm intelligence, 8(3), 159-198.
#'   \item Cheng, R., & Jin, Y. (2014). A competitive swarm optimizer for large scale optimization. IEEE transactions on cybernetics, 45(2), 191-204.
#'   \item Shi, Y., & Eberhart, R. (1998, May). A modified particle swarm optimizer. In Evolutionary Computation Proceedings, 1998. IEEE World Congress on Computational Intelligence., The 1998 IEEE International Conference on (pp. 69-73). IEEE.
#'   \item Sun, J., Wu, X., Palade, V., Fang, W., Lai, C.-H., and Xu, W. (2012). Convergence analysis and improvements of quantum-behaved particle swarm optimization. Information Sciences, 193:81-103.
#'   \item Storn, R., & Price, K. (1997). Differential evolution-a simple and efficient heuristic for global optimization over continuous spaces. Journal of global optimization, 11, 341-359.
#' }
#'
#' @docType package
#' @author Ping-Yang Chen <pychen.ping@gmail.com>
#' @name globpso-package
#' @useDynLib globpso
#' @importFrom Rcpp evalCpp RcppArmadillo cppFunction sourceCpp
#' @importFrom utils globalVariables
