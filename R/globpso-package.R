#' **globpso**: Particle Swarm Optimization Algorithms and Differential Evolution for Minimization Problems
#'
#' A fast and flexible general-purpose implementation of Particle Swarm Optimization (PSO) and Differential Evolution (DE) for solving global minimization problems is provided. It is designed to handle complex optimization tasks with nonlinear, non-differentiable, and multi-modal objective functions defined by users. There are five types of PSO variants: Particle Swarm Optimization (PSO, Eberhart & Kennedy, 1995) <doi:10.1109/MHS.1995.494215>, Quantum-behaved particle Swarm Optimization (QPSO, Sun et al., 2004) <doi:10.1109/CEC.2004.1330875>, Locally convergent rotationally invariant particle swarm optimization (LcRiPSO, Bonyadi & Michalewicz, 2014) <doi:10.1007/s11721-014-0095-1>, Competitive Swarm Optimizer (CSO, Cheng & Jin, 2015) <doi:10.1109/TCYB.2014.2322602> and Double exponential particle swarm optimization (DExPSO, Stehlik et al., 2024) <doi:10.1016/j.asoc.2024.111913>. For the DE algorithm, six types in Storn, R. & Price, K. (1997) <doi:10.1023/A:1008202821328> are included: DE/rand/1, DE/rand/2, DE/best/1, DE/best/2, DE/rand_to-best/1 and DE/rand_to-best/2.
#'
#' @references 
#' \enumerate{
#'   \item Bonyadi, M. R., & Michalewicz, Z. (2014). A locally convergent rotationally invariant particle swarm optimization algorithm. Swarm intelligence, 8(3), 159-198.
#'   \item Cheng, R., & Jin, Y. (2014). A competitive swarm optimizer for large scale optimization. IEEE transactions on cybernetics, 45(2), 191-204.
#   \item Eberhart, R. & Kennedy, J. (1995). A new optimizer using particle swarm theory. In The 6th International Symposium on Micro Machine and Human Science, pages 39-43. IEEE.
#'   \item Shi, Y., & Eberhart, R. (1998, May). A modified particle swarm optimizer. In Evolutionary Computation Proceedings, 1998. IEEE World Congress on Computational Intelligence., The 1998 IEEE International Conference on (pp. 69-73). IEEE.
#'   \item Stehlík, M., Chen, P. Y., Wong, W. K., and Kiseľák, J. (2024). A double exponential particle swarm optimization with non-uniform variates as stochastic tuning and guaranteed convergence to a global optimum with sample applications to finding optimal exact designs in biostatistics. Applied Soft Computing, 163, 111913.
#'   \item Storn, R., & Price, K. (1997). Differential evolution-a simple and efficient heuristic for global optimization over continuous spaces. Journal of global optimization, 11, 341-359.
#'   \item Sun, J., Feng, B., and Xu, W. (2004a). Particle swarm optimization with particles having quantum behavior. In Evolutionary Computation, 2004. CEC2004. Congress on, volume 1, pages 325-331. IEEE.
#' }
#'
#' @docType package
#' @author Ping-Yang Chen <pychen.ping@gmail.com>
#' @name globpso-package
#' @useDynLib globpso, .registration=TRUE
#' @import Rcpp 
#' @import RcppArmadillo
#' @importFrom Rcpp evalCpp cppFunction sourceCpp
#' @importFrom utils globalVariables
"_PACKAGE"