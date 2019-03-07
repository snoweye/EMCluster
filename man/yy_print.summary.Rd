\name{Print and Summary}
\alias{print.emret}
\alias{summary.emret}
\alias{summary.svd}
\title{Functions for Printing or Summarizing Objects According to Classes}
\description{
  Several classes are declared in \pkg{EMCluster}, and these are functions
  to print and summary objects.
}
\usage{
\method{print}{emret}(x, digits = max(4, getOption("digits") - 3), ...)
\method{summary}{emret}(object, ...)
\method{summary}{svd}(object, ...)
}
\arguments{
  \item{x}{an object with the class attributes.}
  \item{digits}{for printing out numbers.}
  \item{object}{an object with the class attributes.}
  \item{\dots}{other possible options.}
}
\details{
  These are useful functions for summarizing and debugging.
}
\value{
  The results will cat or print on the STDOUT by default.
}
\references{
  \url{https://www.stat.iastate.edu/people/ranjan-maitra/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com} and Ranjan Maitra.
}
\seealso{
  \code{\link{init.EM}}, \code{\link{emcluster}}, \code{\link{starts.via.svd}}.
}
\examples{
\dontrun{
library(EMCluster, quietly = TRUE)
x2 <- da2$da

emobj <- list(pi = da2$pi, Mu = da2$Mu, LTSigma = da2$LTSigma)
eobj <- e.step(x2, emobj = emobj)
emobj <- m.step(x2, emobj = eobj)
summary(emobj)

ret <- starts.via.svd(x2, nclass = 10, method = "kmeans")
summary(ret)
}
}
\keyword{summary}
