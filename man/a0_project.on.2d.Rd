\name{Projection On 2D}
\alias{project.on.2d}
\title{Produce Projection on 2D}
\description{
  The function projects multivariate data on 2D plane which can be displayed
  by \code{\link{plotppcontour}()} later.
}
\usage{
project.on.2d(x, emobj = NULL, pi = NULL, Mu = NULL,
    LTSigma = NULL, class = NULL, method = c("PP", "SVD"))
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
  \item{class}{id of classifications, length \eqn{n}{n}.}
  \item{method}{either projection pursuit or singular value decomposition.}
}
\details{
  This function produces projection outputs of \code{x} and \code{emobj}.
}
\value{
  A projection is returned which is a list contains
  \itemize{
    \item \code{da} is a \eqn{n \times 2}{n * 2} projected matrix of \code{x}.
    \item \code{Pi} is the original proportion \code{emobj$pi} of
                    length \eqn{K}{K}.
    \item \code{Mu} is a \eqn{K \times 2}{K * 2} projected matrix of
                    \code{emboj$Mu}.
    \item \code{S} is a \eqn{2 \times 2 \times K}{2 * 2 * K} projected array of
                    \code{emboj$LTSigma}.
    \item \code{class} is the original class id \code{emobj$class}.
    \item \code{proj.mat} is the projection matrix of dimension
                          \eqn{p \time 2}{p * 2}.
  }
}
\references{
  \url{https://www.stat.iastate.edu/people/ranjan-maitra/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com} and Ranjan Maitra.
}
\seealso{
  \code{\link{project.on.2d}()}.
}
\examples{
\dontrun{
library(EMCluster, quietly = TRUE)
set.seed(1234)

### Iris.
x <- as.matrix(iris[, 1:4])
ret <- init.EM(x, nclass = 3, min.n = 30)
ret.proj <- project.on.2d(x, ret)

### Plot.
pdf("iris_ppcontour.pdf", height = 5, width = 5)
plotppcontour(ret.proj$da, ret.proj$Pi, ret.proj$Mu, ret.proj$S,
              ret.proj$class, main = "Iris K = 3")
dev.off()
}
}
\keyword{utility}
