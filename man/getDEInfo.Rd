% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/de.R
\name{getDEInfo}
\alias{getDEInfo}
\title{Generation function of DE parameter settings}
\usage{
getDEInfo(
  nPop = 32,
  maxIter = 100,
  deType = "rand-1",
  freeRun = 1,
  tol = 1e-06,
  sf = 0.5,
  cr = 0.1
)
}
\arguments{
\item{nPop}{A integer number of population size in DE algorithm.}

\item{maxIter}{A integer number of maximal DE iterations.}

\item{deType}{string. The type of DE. This package current supports the following types:
\itemize{
\item{"rand-1"}{ Mutation operation on the current position with one random direction}
\item{"rand-2"}{ Mutation operation on the current position with two random directions}
\item{"best-1"}{ Mutation operation on the best position with one random direction}
\item{"best-2"}{ Mutation operation on the best position with two random directions}
\item{"rand_to-best-1"}{ Mutation operation on the current position with direction to the best and one random direction}
\item{"rand-to-best-2"}{ Mutation operation on the current position with direction to the best and two random directions}
}
The default type is `rand-1`.}

\item{freeRun}{A number between \eqn{[0,1]} that controls the percentage of DE iterations which are free from examining the
stopping criterion, \eqn{|f'-f|<\varepsilon}
where \eqn{f'} and \eqn{f} are the objective function values in the previous and current iterations, respectively.
The default is 1.0 implying the DE will completely ignore the stopping criterion.
Otherwise, the DE checks the stopping criterion after free iterations.}

\item{tol}{A small value for the tolerance, \eqn{\varepsilon}, in the stopping criterion.
For \code{freeRun} smaller than 1.0, the default is \code{1e-6}. Otherwise, this value would not affect the algorithm.}

\item{sf}{The value of scaling factor in DE updating procedure. The default is 0.5.}

\item{cr}{The value of crossover rate in DE updating procedure. The default is 0.1.}

\item{checkConv}{A logical value which controls whether DE checks the stopping criterion during updating procedure.
Specify \code{TRUE} for DE to compute the stopping criterion \eqn{|f'-f|<\varepsilon}
where \eqn{f'} and \eqn{f} are the objective function values in the previous and current iterations, respectively.
The default is \code{FALSE}.}
}
\value{
A list of DE parameter settings.
}
\description{
Create a list with DE parameters for Minimization.
}
\examples{
# Get default settings with specified population size and maximal number of iterations.
DE_INFO <- getDEInfo(nPop = 32, maxIter = 100)
}
