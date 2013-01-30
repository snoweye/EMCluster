\name{Information Criteria}
\alias{em.ic}
\alias{em.aic}
\alias{em.bic}
\alias{em.clc}
\alias{em.icl}
\alias{em.icl.bic}
\title{Information Criteria for Model-Based Clustering}
\description{
  These functions are tools for compute information criteria
  for the fitted models.
}
\usage{
em.ic(x, emobj = NULL, pi = NULL, Mu = NULL, LTSigma = NULL,
      llhdval = NULL)
em.aic(x, emobj = NULL, pi = NULL, Mu = NULL, LTSigma = NULL)
em.bic(x, emobj = NULL, pi = NULL, Mu = NULL, LTSigma = NULL)
em.clc(x, emobj = NULL, pi = NULL, Mu = NULL, LTSigma = NULL)
em.icl(x, emobj = NULL, pi = NULL, Mu = NULL, LTSigma = NULL)
em.icl.bic(x, emobj = NULL, pi = NULL, Mu = NULL, LTSigma = NULL)
}
\arguments{
  \item{x}{the data matrix, dimension \eqn{n\times p}{n * p}.}
  \item{emobj}{the desired model which is a list mainly contains \code{pi},
               \code{Mu}, and \code{LTSigma}, usually a returned object from
               \code{init.EM}.}
  \item{pi}{the mixing proportion, length \eqn{K}{K}.}
  \item{Mu}{the centers of clusters, dimension \eqn{K\times p}{K * p}.}
  \item{LTSigma}{the lower triangular matrices of dispersion,
                 \eqn{K\times p(p+1)/2}{K * p(p+1)/2}.}
  \item{llhdval}{the total log likelihood value of \code{x} given \code{emobj}.}
}
\details{
  The \code{em.ic} calls all other functions to compute AIC (\code{em.aic}),
  BIC (\code{em.bic}), CLC (\code{em.clc}), ICL (\code{em.icl}), and
  ICL.BIC (\code{em.icl.bic}). All are useful information criteria for
  model selections, mainly choosing number of cluster.
}
\value{
  \code{em.ic} returns a list containing all other information criteria
  for given the data \code{x} and the desired model \code{emobj}.
}
\references{
  \url{http://maitra.public.iastate.edu/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}, Ranjan Maitra, and
  Volodymyr Melnykov.
}
\seealso{
  \code{\link{init.EM}}.
}
\examples{
\dontrun{
library(EMCluster, quiet = TRUE)
x2 <- da2$da

emobj <- list(pi = da2$pi, Mu = da2$Mu, LTSigma = da2$LTSigma)
em.ic(x2, emobj = emobj)
}
}
\keyword{tool}
