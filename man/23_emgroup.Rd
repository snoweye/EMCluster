\name{Other Initializations}
\alias{starts.via.svd}
\alias{emgroup}
\title{Other Initializations}
\description{
  Two more functions with different initialization method.
}
\usage{
starts.via.svd(x, nclass = 1, method = c("em", "kmeans"),
               EMC = .EMC)
emgroup(x, nclass = 1, EMC = .EMC)
}
\arguments{
  \item{x}{the data matrix, dimension \eqn{n\times p}{n * p}.}
  \item{nclass}{the desired number of clusters, \eqn{K}{K}.}
  \item{method}{method with the svd initializations.}
  \item{EMC}{the control for the EM iterations.}
}
\details{
  The \code{starts.via.svd} utilizes SVD to initial parameters,
  and the \code{emgroup} runs the EM algorithm starting from the
  initial.
}
\value{
  The \code{starts.via.svd} returns an object with class \code{svd},
  and the \code{emgroup} returns and object \code{emobj} with class
  \code{emret}.
}
\references{
  \url{http://maitra.public.iastate.edu/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com} and Ranjan Maitra.
}
\seealso{
  \code{\link{init.EM}}, \code{\link{.EMControl}}.
}
\examples{
%\dontrun{
library(EMCluster, quietly = TRUE)
set.seed(1234)
x1 <- da1$da

emobj <- emgroup(x1, nclass = 10)
summary(emobj)

ret.0 <- starts.via.svd(x1, nclass = 10, method = "kmeans")
summary(ret.0)
%}
}
\keyword{programming}
