\name{Assign Class}
\alias{assign.class}
\title{Assign Class Id}
\description{
  This function assigns cluster id to each observation in \code{x}
  according to the desired model \code{emobj} or specified
  parameters \code{pi}, \code{Mu}, and \code{LTSigma}. 
}
\usage{
  assign.class(x, emobj = NULL, pi = NULL, Mu = NULL, LTSigma = NULL,
               lab = NULL, return.all = TRUE)
}
\arguments{
  \item{x}{the data matrix, dimension \eqn{n\times p}{n * p}.}
  \item{emobj}{the desired model which is a list mainly contains \code{pi},
               \code{Mu}, and \code{LTSigma}, usually a returned object from
               \code{init.EM}.}
  \item{pi}{the mixing proportion, length \eqn{K}{K}.}
  \item{Mu}{the centers of clusters, dimension \eqn{K\times p}{K * p}.}
  \item{LTSigma}{the lower triangular matrices of dispersion, dimension
                 \eqn{K\times p(p+1)/2}{K * p(p+1)/2}.}
  \item{lab}{labeled data for semi-supervised clustering,
             length \eqn{n}{n}.}
  \item{return.all}{if returning with a whole \code{emobj} object.}
}
\details{
  This function are based either an input \code{emobj} or inputs \code{pi},
  \code{Mu}, and \code{LTSigma} to assign class id to each observation of
  \code{x}.

  If \code{lab} is submitted, then the observation with label id greater 0
  will not be assigned new class.
}
\value{
  This function returns a list containing mainly two new variables:
  \code{nc} (length \eqn{K}{K} numbers of observations in each class) and
  \code{class} (length \eqn{n}{n} class id).
}
\references{
  \url{https://www.stat.iastate.edu/people/ranjan-maitra/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com} and Ranjan Maitra.
}
\seealso{
  \code{\link{init.EM}}, \code{\link{emcluster}}.
}
\examples{
%\dontrun{
library(EMCluster, quietly = TRUE)
set.seed(1234)
x2 <- da2$da

ret <- init.EM(x2, nclass = 2)
ret.new <- assign.class(x2, ret, return.all = FALSE)
str(ret.new)
%}
}
\keyword{tool}
