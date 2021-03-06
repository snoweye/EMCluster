\name{Dataset}
\alias{da1}
\alias{da2}
\alias{da3}
\docType{data}
\title{Dataset for demonstrations}
\description{
  There are four small datasets to test and demonstrate \pkg{EMCluster}.
}
\usage{
da1
da2
da3
}
\format{
  \code{da1}, \code{da2}, \code{da3} are in \code{list}.
}
\details{
  \code{da1} has 500 observations in two dimensions
  \code{da1$da$x} and \code{da1$da$y},
  and they are in 10 clusters given in \code{da1$class}.

  \code{da2} has 2,500 observations in two dimensions, too. The true
  parameters are given in \code{da1$pi}, \code{da1$Mu}, and \code{da1$LTSigma}.
  There are 40 clusters given in \code{da1$class} for this dataset.

  \code{da3} is similar to \code{da2}, but with lower overlaps between
  clusters.
}
\references{
  \url{https://www.stat.iastate.edu/people/ranjan-maitra}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com} and Ranjan Maitra.
}
\keyword{datasets}
