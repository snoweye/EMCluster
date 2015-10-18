\name{Recolor Classification IDs}
\alias{recolor}
\title{Recolor Classification IDs}
\description{
  This function returns new classification IDs.
}
\usage{
recolor(id.class, id.target)
}
\arguments{
  \item{id.class}{original class ids.}
  \item{id.target}{target class ids.}
}
\details{
  This function colors id.target in accordance with the most likely
  candidate in id.class.

  All ids should be positive integers and
  started from 1 to K, and the maximums are allowed to be different.
}
\value{
  Return a list contains two factors, id.trcl and id.prcl, which are
  factor versions of id.class and id.target.
}
\references{
  \url{http://maitra.public.iastate.edu/}
}
\author{
  Ranjan Maitra.
}
\examples{
\dontrun{
library(EMCluster, quietly = TRUE)

true.id <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
pred.id <- c(2, 1, 2, 1, 1, 1, 2, 1, 1)

recolor(pred.id, true.id)
}
}
\keyword{summary}
