\name{Recolor Classification IDs}
\alias{rematch}
\alias{recode}
\alias{recolor}
\title{Recolor Classification IDs}
\description{
  These functions return new classification IDs.
}
\usage{
recolor(id.target, id.class, scatter.class = NULL, scatter.target = NULL)
rematch(tg.id, cl.id)
recode(id)
}
\arguments{
  \item{id.target}{target class ids.}
  \item{id.class}{original class ids.}
  \item{scatter.class}{scatter class ids.}
  \item{scatter.target}{scatter target class ids.}
  \item{id}{class ids.}
  \item{tg.id}{target class ids.}
  \item{cl.id}{class ids.}
}
\details{
  The function \code{recolor} colors \code{id.target} in accordance with the
  most likely candidate in \code{id.class}.
  Note that if scatter is present, then the class given by 0 is represented
  as scatter and it is assumed to be the same for both classifications.

  The function \code{rematch} returns a list as \code{id.trcl} and
  \code{id.prcl}. It is the heart of the recolor function
  and is usually called from recolor.

  The function \code{recode} reoders classes to eliminate group ids without
  any members. It is assumed that the group ids are integers.
}
\value{
  See Details.
}
\references{
  \url{https://www.stat.iastate.edu/people/ranjan-maitra}
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
