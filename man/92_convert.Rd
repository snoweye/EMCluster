\name{Conversion}
\alias{LTSigma2variance}
\alias{variance2LTSigma}
\alias{LTsigma2var}
\alias{var2LTsigma}
\alias{class2Gamma}
\alias{Gamma2class}
\title{Convert Matrices in Different Format}
\description{
  These utility functions are to convert matrices in different formats.
}
\usage{
LTSigma2variance(x)
variance2LTSigma(x)
LTsigma2var(x1, p = NULL)
var2LTsigma(x1)
class2Gamma(class)
Gamma2class(Gamma)
}
\arguments{
  \item{x}{a matrix/array to be converted, the dimension could be
           \eqn{K\times p(p+1)/2}{K * p(p+1)/2} or
           \eqn{p\times p\times K}{p * p * K}.}
  \item{x1}{a vector/matrix to be converted, the length and dimension could be
            \eqn{p(p+1)/2}{p(p+1)/2} and \eqn{p\times p}{p * p}.}
  \item{p}{dimension of matrix.}
  \item{class}{id of clusters for each observation, length \eqn{n}{n}.}
  \item{Gamma}{containing posterior probabilities if normalized,
               otherwise containing component densities weighted by
               mixing proportion, dimension \eqn{n\times K}{n * K}.}
}
\details{
  \code{LTSigma2variance} converts \code{LTSigma} format to 3D array, and
  \code{variance2LTSigma} is the inversion function. 

  \code{LTsigma2var} converts \code{LTsigma} format to a matrix, and
  \code{var2LTsigma} is the inversion function.
  Note that \code{LTsigma} is one component of \code{LTSigma}.

  \code{class2Gamma} converts id to a Gamma matrix where with probability 1
  for the cluster where the observation belongs to, and \code{Gamma2class}
  converts posterior to cluster id where largest posterior is picked for
  each observation.
}
\value{
  A vector/matrix/array is returned.
}
\references{
  \url{https://www.stat.iastate.edu/people/ranjan-maitra}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com} and Ranjan Maitra
}
\seealso{
  \code{\link{init.EM}}, \code{\link{emcluster}}.
}
\examples{
\dontrun{
library(EMCluster, quietly = TRUE)
x <- da2$LTSigma
class <- da2$class

y <- LTSigma2variance(x)
str(y)
y <- variance2LTSigma(y)
str(y)
sum(x != y)

Gamma <- class2Gamma(class)
class.new <- Gamma2class(Gamma)
sum(class != class.new)
}
}
\keyword{utility}
