\name{Plot Projection and Contour}
\alias{plotppcontour}
\title{Plot Contour}
\description{
  The function plots multivariate data on 2D plane with contour.
  Typically, the contour is built via projection pursuit or SVD
  algorithms, such as \code{\link{project.on.2d}()}.
}
\usage{
plotppcontour(da, Pi, Mu, S, class, class.true = NULL, n.grid = 128,
    angle = 0, xlab = "", ylab = "", main = "")
}
\arguments{
  \item{da}{a projected data matrix, dimension \eqn{n\times 2}{n * 2}.}
  \item{Pi}{proportion, length \eqn{K}{K}.}
  \item{Mu}{the projected centers of cluster, dimension \eqn{K\times 2}{K * 2}.}
  \item{S}{projected matrices of dispersion, dimension
           \eqn{p\times p \times K}{2 * 2 * K}.}
  \item{class}{id of classifications, length \eqn{n}{n}.}
  \item{class.true}{ture id of classifications if available, length \eqn{n}{n}.}
  \item{n.grid}{number of grid points.}
  \item{angle}{a rotation angle (\eqn{0}{0} to \eqn{2\pi}{2 pi}).}
  \item{xlab}{an option for \code{plot()} function.}
  \item{ylab}{an option for \code{plot()} function.}
  \item{main}{an option for \code{plot()} function.}
}
\details{
  This function plots projection output of \code{project.on.2d()}.

  \code{da}, \code{Mu}, and \code{S} are projected by some projection matrices
  obtained via SVD or projection pursuit algorithms. The projection is made
  on a 2D plane in the direction in which clusters of data \code{x}
  are most distinguishable to visualize.
}
\value{
  A 2D projection plot is returned.
}
\note{
  Only distinguishable for up to 7 clusters due to the limited color schemes.
}
\references{
  \url{https://www.stat.iastate.edu/people/ranjan-maitra}
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
library(MASS, quietly = TRUE)
set.seed(1234)

### Crabs.
x <- as.matrix(crabs[, 4:8])
ret <- init.EM(x, nclass = 4, min.n = 20)
ret.proj <- project.on.2d(x, ret)

### Plot.
pdf("crabs_ppcontour.pdf", height = 5, width = 5)
plotppcontour(ret.proj$da, ret.proj$Pi, ret.proj$Mu, ret.proj$S,
              ret.proj$class, angle = pi/6, main = "Crabs K = 4")
dev.off()
}
}
\keyword{utility}
