\name{Rand Index}
\alias{RRand}
\title{Rand Index and Adjusted Rand Index}
\description{
  This function returns the Rand index and the adjusted Rand index for
  given true class ids and predicted class ids.
}
\usage{
RRand(trcl, prcl, lab = NULL)
}
\arguments{
  \item{trcl}{true class ids.}
  \item{prcl}{predicted class ids.}
  \item{lab}{known ids for semi-supervised clustering.}
}
\details{
  All ids, \code{trcl} and \code{prcl}, should be positive integers and
  started from 1 to K, and the maximums are allowed to be different.

  \code{lab} used in semi-supervised clustering contains the labels which
  are known before clustering. It should be positive integer and
  started from 1 for labeled data and 0 for unlabeled data.
}
\value{
  Return a \code{Class RRand} contains Rand index and adjusted
  Rand index.
}
\references{
  \url{http://maitra.public.iastate.edu/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com} and Ranjan Maitra.
}
\examples{
%\dontrun{
library(EMCluster, quietly = TRUE)

true.id <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
pred.id <- c(2, 1, 2, 1, 1, 1, 2, 1, 1)
label   <- c(0, 0, 0, 0, 1, 0, 2, 0, 0)

RRand(true.id, pred.id)
RRand(true.id, pred.id, lab = label)
%}
}
\keyword{summary}
