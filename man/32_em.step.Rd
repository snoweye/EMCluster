\name{Single Step}
\alias{e.step}
\alias{m.step}
\title{Single E- and M-step}
\description{
  These functions are single E- and M-step of EM algorithm for
  model-based clustering of finite mixture multivariate Gaussian
  distribution with unstructured dispersion.
}
\usage{
e.step(x, emobj = NULL, pi = NULL, Mu = NULL, LTSigma = NULL,
          norm = TRUE)
m.step(x, emobj = NULL, Gamma = NULL, assign.class = FALSE)
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
  \item{norm}{if returning normalized \code{Gamma}.}
  \item{Gamma}{containing posterior probabilities if normalized,
               otherwise containing component densities weighted by
               mixing proportion, dimension \eqn{n\times K}{n * K}.}
  \item{assign.class}{if assigning class id.}
}
\details{
  These two functions are mainly used in debugging for development and
  post process after model fitting.
}
\value{
  The \code{e.step} returns a list contains \code{Gamma}, the posterior
  probabilities if \code{norm=TRUE}, otherwise it contains component densities.
  This is one E-step and \code{Gamma} is used to update \code{emobj} in
  the M-step next.

  The \code{m.step} returns a new \code{emobj} according to the \code{Gamma}
  from the E-step above.
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
eobj <- e.step(x2, emobj = emobj)
emobj <- m.step(x2, emobj = eobj)
emobj
}
}
\keyword{tool}
