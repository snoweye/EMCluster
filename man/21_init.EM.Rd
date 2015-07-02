\name{Initialization and EM}
\alias{init.EM}
\alias{em.EM}
\alias{rand.EM}
\alias{exhaust.EM}
\title{Initialization and EM Algorithm}
\description{
  These functions perform initializations (including em.EM and RndEM)
  followed by the EM iterations for model-based clustering of finite mixture
  multivariate Gaussian distribution with unstructured dispersion
  in both of unsupervised and semi-supervised clusterings.
}
\usage{
init.EM(x, nclass = 1, lab = NULL, EMC = .EMC,
        stable.solution = TRUE, min.n = NULL, min.n.iter = 10,
        method = c("em.EM", "Rnd.EM"))
em.EM(x, nclass = 1, lab = NULL, EMC = .EMC,
      stable.solution = TRUE, min.n = NULL, min.n.iter = 10)
rand.EM(x, nclass = 1, lab = NULL, EMC = .EMC.Rnd,
        stable.solution = TRUE, min.n = NULL, min.n.iter = 10)
exhaust.EM(x, nclass = 1, lab = NULL,
           EMC = .EMControl(short.iter = 1, short.eps = Inf),
           method = c("em.EM", "Rnd.EM"),
           stable.solution = TRUE, min.n = NULL, min.n.iter = 10);
}
\arguments{
  \item{x}{the data matrix, dimension \eqn{n\times p}{n * p}.}
  \item{nclass}{the desired number of clusters, \eqn{K}{K}.}
  \item{lab}{labeled data for semi-supervised clustering,
             length \eqn{n}{n}.}
  \item{EMC}{the control for the EM iterations.}
  \item{stable.solution}{if returning a stable solution.}
  \item{min.n}{restriction for a stable solution, the minimum number
               of observations for every final clusters.}
  \item{min.n.iter}{restriction for a stable solution, the minimum number
                    of iterations for trying a stable solution.}
  \item{method}{an initialization method.}
}
\details{
  The \code{init.EM} calls either \code{em.EM} if \code{method="em.EM"} or
  \code{rand.EM} if \code{method="Rnd.EM"}.

  The \code{em.EM} has two steps: short-EM has loose convergent
  tolerance controlled by \code{.EMC$short.eps} and try several random
  initializations controlled by \code{.EMC$short.iter}, while long-EM
  starts from the best short-EM result (in terms of log likelihood) and
  run to convergence with a tight tolerance controlled by \code{.EMC$EM.eps}.

  The \code{rand.EM} also has two steps: first randomly pick several
  random initializations controlled by \code{.EMC$short.iter}, and
  second starts from the best of the random result
  (in terms of log likelihood) and run to convergence.

  The \code{lab} is only for the semi-supervised clustering, and it contains
  pre-labeled indices between 1 and \eqn{K}{K} for labeled observations.
  Observations with index 0 is non-labeled and has to be clustered by
  the EM algorithm. Indices will be assigned by the results of the EM
  algorithm. See \code{demo(allinit_ss,'EMCluster')} for details.

  The \code{exhaust.EM} also calls the \code{init.EM} with different
  \code{EMC} and perform \code{exhaust.iter} times of EM algorithm
  with different initials. The best result is returned.
}
\value{
  These functions return an object \code{emobj} with class \code{emret}
  which can be used in post-process or other functions such as
  \code{e.step}, \code{m.step}, \code{assign.class}, \code{em.ic},
  and \code{dmixmvn}.
}
\references{
  \url{http://maitra.public.iastate.edu/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com} and Ranjan Maitra.
}
\seealso{
  \code{\link{emcluster}}, \code{\link{.EMControl}}.
}
\examples{
\dontrun{
library(EMCluster, quietly = TRUE)
set.seed(1234)
x <- da1$da

ret.em <- init.EM(x, nclass = 10, method = "em.EM")
ret.Rnd <- init.EM(x, nclass = 10, method = "Rnd.EM", EMC = .EMC.Rnd)

emobj <- simple.init(x, nclass = 10)
ret.init <- emcluster(x, emobj, assign.class = TRUE)

par(mfrow = c(2, 2))
plotem(ret.em, x)
plotem(ret.Rnd, x)
plotem(ret.init, x)
}
}
\keyword{programming}
