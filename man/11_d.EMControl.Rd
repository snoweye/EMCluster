\name{EM Control}
\alias{.EMControl}
\alias{.EMC}
\alias{.EMC.Rnd}
\alias{.EMC.Rndp}
\title{EM Control Generator and Controller}
\description{
  The \code{.EMControl} generates an EM control (\code{.EMC})
  controlling the options and conditions of EM algorithms,
  i.e. this function generate a default template.
  One can either modify \code{.EMC} or employ this function to
  control EM algorithms.
  By default, \code{.EMC}, \code{.EMC.Rnd}, and \code{.EC.Rndp} are
  three native controllers as the \pkg{EMCluster} is loaded.
}
\usage{
.EMControl(alpha = 0.99, short.iter = 200, short.eps = 1e-2,
    fixed.iter = 1, n.candidate = 3,
    em.iter = 1000, em.eps = 1e-6, exhaust.iter = 5)
.EMC
.EMC.Rnd
.EMC.Rndp
}
\arguments{
  \item{alpha}{only used in \code{emgroup} for "SVD" initialization.}
  \item{short.iter}{number of short-EM steps, default = 200.}
  \item{short.eps}{tolerance of short-EM steps, default = 1e-2.}
  \item{fixed.iter}{fixed iterations of EM for "RndEM" initialization,
                    default = 1.}
  \item{n.candidate}{reserved for other initialization methods (unimplemented).}
  \item{em.iter}{maximum number of long-EM steps, default = 1000.}
  \item{em.eps}{tolerance of long-EM steps, default = 1e-6.}
  \item{exhaust.iter}{number of iterations for "exhaustEM" initialization,
                      default = 5.}
}
\details{
  \code{exhaust.iter} and \code{fixed.iter} are used to control the
  iterations of initialization procedures.

  \code{short.iter} and \code{short.eps} are used to control the
  short-EM iterations.

  \code{em.iter} and \code{em.eps} are used to control the long-EM iterations.

  Moeover, \code{short.eps} and \code{em.eps} are for checking convergence of
  the iterations.
}
\value{
  This function returns a list as \code{.EMC} by default.

  The \code{.EMC.Rnd} is equal to \code{.EMControl(short.eps = Inf)} and
  usually used by the \code{rand.EM} method.

  The \code{.EMC.Rndp} is equal to \code{.EMControl(fixed.iter = 5)} where
  each random initials run 5 EM iterations in the \code{rand.EM} method.
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

.EMC <- .EMControl()
.EMC.Rnd <- .EMControl(short.eps = Inf)
.EMC.Rndp <- .EMControl(fixed.iter = 5)
%}
}
\keyword{programming}
