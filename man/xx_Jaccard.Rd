\name{Jaccard Index}
\alias{Jaccard.Index}
\title{Jaccard Index}
\description{
  This function returns the Jaccard index for binary ids.
}
\usage{
Jaccard.Index(x, y)
}
\arguments{
  \item{x}{true binary ids, 0 or 1.}
  \item{y}{predicted binary ids, 0 or 1.}
}
\details{
  All ids, \code{x} and \code{y}, should be either 0 (not active) or 1 (active).
  Any value other than 1 will be converted to 0.
}
\value{
  Return the value of Jaccard index.
}
\references{
  \url{https://www.stat.iastate.edu/people/ranjan-maitra}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com} and Ranjan Maitra.
}
\examples{
%\dontrun{
library(EMCluster, quietly = TRUE)

x.id <- c(1, 1, 1, 0, 0, 0, 3, 3, 3)
y.id <- c(0, 1, 0, 1, 1, 1, 0, 1, 1)

Jaccard.Index(x.id, y.id)
%}
}
\keyword{summary}
