\name{MVN}
\alias{dmvn}
\alias{dlmvn}
\alias{dmixmvn}
\alias{logL}
\title{Density of (Mixture) Multivariate Normal Distribution}
\description{
  These functions are tools for compute density of (mixture) multivariate
  Gaussian distribution with unstructured dispersion.
}
\usage{
dmvn(x, mu, LTsigma, log = FALSE)
dlmvn(x, mu, LTsigma, log = TRUE)
dmixmvn(x, emobj = NULL, pi = NULL, Mu = NULL, LTSigma = NULL, log = FALSE)
logL(x, emobj = NULL, pi = NULL, Mu = NULL, LTSigma = NULL)
}
\arguments{
  \item{x}{the data matrix, dimension \eqn{n\times p}{n * p}.}
  \item{mu}{the centers of clusters, length \eqn{p}{p}.}
  \item{LTsigma}{the lower triangular matrices of dispersion, length
                 \eqn{p(p+1)/2}{p(p+1)/2}.}
  \item{log}{if logarithm returned.}

  \item{emobj}{the desired model which is a list mainly contains \code{pi},
               \code{Mu}, and \code{LTSigma}, usually a returned object from
               \code{init.EM}.}
  \item{pi}{the mixing proportion, length \eqn{K}{K}.}
  \item{Mu}{the centers of clusters, dimension \eqn{K\times p}{K * p}.}
  \item{LTSigma}{the lower triangular matrices of dispersion,
                 \eqn{K\times p(p+1)/2}{K * p(p+1)/2}.}
}
\details{
  The \code{dmvn} and \code{dlmvn} compute density and log density of
  multivariate distribution.

  The \code{dmixmvn} computes density of mixture multivariate distribution
  and is based either an input \code{emobj} or inputs \code{pi},
  \code{Mu}, and \code{LTSigma} to assign class id to each observation of
  \code{x}.

  The \code{logL} returns the value of the observed log likelihood function
  of the parameters at the current values of the parameters \code{pi},
  \code{Mu}, and \code{LTSigma}, with the suplied data matrix \code{x}.
}
\value{
  A density value is returned.
}
\references{
  \url{https://maitra.public.iastate.edu/}
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
x2 <- da2$da
x3 <- da3$da

emobj2 <- list(pi = da2$pi, Mu = da2$Mu, LTSigma = da2$LTSigma)
emobj3 <- list(pi = da3$pi, Mu = da3$Mu, LTSigma = da3$LTSigma)

logL(x2, emobj = emobj2)
logL(x3, emobj = emobj3)

dmixmvn2 <- dmixmvn(x2, emobj2)
dmixmvn3 <- dmixmvn(x3, emobj3)

dlmvn(da2$da[1,], da2$Mu[1,], da2$LTSigma[1,])
log(dmvn(da2$da[1,], da2$Mu[1,], da2$LTSigma[1,]))
%}
}
\keyword{utility}
