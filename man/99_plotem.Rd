\name{Plot EM Results}
\alias{plotem}
\alias{plot2d}
\title{Plot Two Dimensional Data with clusters}
\description{
  The functions plot two dimensional data for clusters.
}
\usage{
plotem(emobj, x, main = NULL, xlab = NULL, ylab = NULL,
       ...)
plot2d(x, emobj = NULL, k = NULL, color.pch = 1,
       append.BN = TRUE, ...)
}
\arguments{
  \item{emobj}{the desired model which is a list mainly contains \code{pi},
               \code{Mu}, and \code{LTSigma}, usually a returned object from
               \code{init.EM}.}
  \item{x}{the data matrix, dimension \eqn{n\times p}{n * p}.}
  \item{main}{title of plot.}
  \item{xlab}{label of x-axis.}
  \item{ylab}{label of y-axis.}
  \item{\dots}{other parameters to the plot.}

  \item{k}{index for symbols.}
  \item{color.pch}{color and style for symbols.}
  \item{append.BN}{if appending bivariate normal ellipsoid.}
}
\details{
  This a simple x-y lot.
}
\value{
  A plot is returned.
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
\dontrun{
library(EMCluster, quietly = TRUE)
x1 <- da1$da

ret.1 <- starts.via.svd(x1, nclass = 10, method = "em")
summary(ret.1)

plotem(ret.1, x1)
}
}
\keyword{utility}
