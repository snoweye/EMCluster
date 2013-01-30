\name{EM Algorithm}
\alias{emcluster}
\alias{shortemcluster}
\alias{simple.init}
\alias{emret class}
\alias{emobj object}
\title{EM Algorithm for model-based clustering}
\description{
  These are core functions of \pkg{EMCluster} performing EM algorithm
  for model-based clustering of finite mixture multivariate Gaussian
  distribution with unstructured dispersion.
}
\usage{
emcluster(x, emobj = NULL, pi = NULL, Mu = NULL, LTSigma = NULL,
          lab = NULL, EMC = .EMC, assign.class = FALSE)
shortemcluster(x, emobj = NULL, pi = NULL, Mu = NULL,
               LTSigma = NULL, maxiter = 100, eps = 1e-2)
simple.init(x, nclass = 1)
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
  \item{lab}{labeled data for semi-supervised clustering,
             length \eqn{n}{n}.}
  \item{EMC}{the control for the EM iterations.}
  \item{assign.class}{if assigning class id.}

  \item{maxiter}{maximum number of iterations.}
  \item{eps}{convergent tolerance.}

  \item{nclass}{the desired number of clusters, \eqn{K}{K}.}
}
\details{
  The \code{emcluster} mainly performs EM iterations starting from the given
  parameters \code{emobj} without other initializations.

  The \code{shortemcluster} performs short-EM iterations as described in
  \code{init.EM}.
}
\value{
  The \code{emcluster} returns an object \code{emobj} with class \code{emret}
  which can be used in post-process or other functions such as
  \code{e.step}, \code{m.step}, \code{assign.class}, \code{em.ic},
  and \code{dmixmvn}.

  The \code{shortemcluster} also returns an object \code{emobj} with class
  \code{emret} which is the best of several random initializations.

  The \code{simple.init} utilizes \code{rand.EM} to obtain a simple initial.
}
\references{
  \url{http://maitra.public.iastate.edu/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}, Ranjan Maitra, and
  Volodymyr Melnykov.
}
\seealso{
  \code{\link{init.EM}}, \code{\link{e.step}}, \code{\link{m.step}},
  \code{\link{.EMControl}}.
}
\examples{
\dontrun{
library(EMCluster, quiet = TRUE)
set.seed(1234)
x1 <- da1$da

emobj <- simple.init(x1, nclass = 10)
emobj <- shortemcluster(x1, emobj)
summary(emobj)

ret <- emcluster(x1, emobj, assign.class = TRUE)
summary(ret)
}
}
\keyword{programming}
