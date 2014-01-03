\name{Likelihood Mixture Tests}
\alias{lmt}
\title{Likelihood Mixture Tests}
\description{
  This function test two mixture Gaussian models with unstructured
  covariance matrix and different numbers of clusters.
}
\usage{
lmt(emobj.0, emobj.a, x, tau = 0.5, n.mc = 1000, verbose = FALSE)
}
\arguments{
  \item{emobj.0}{a \code{emret} object for the null hypothesis.}
  \item{emobj.a}{a \code{emret} object for the alternative hypothesis.}
  \item{x}{the data matrix, dimension \eqn{n\times p}{n * p}.}
  \item{tau}{proportion of null and alternative hypotheses.}
  \item{n.mc}{number of Monte Carlo simulations.}
  \item{verbose}{if verbose.}
}
\details{
  This function calls several subroutines to compute information,
  likelihood ratio statistics, degrees of freedom, non-centrality
  of chi-squared distributions \dots etc. Based on Monte Carlo methods
  to estimate parameters of likelihood mixture tests, this function
  return a p-value for testing H0: \code{emobj.0} v.s. Ha: \code{emboj.a}.
}
\value{
  A list of class \code{lmt} are returned.
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
library(EMCluster, quietly = TRUE)
set.seed(1234)

x <- as.matrix(iris[, 1:4])
p <- ncol(x)
min.n <- p * (p + 1) / 2
.EMC$short.iter <- 200

ret.2 <- init.EM(x, nclass = 2, min.n = min.n, method = "Rnd.EM")
ret.3 <- init.EM(x, nclass = 3, min.n = min.n, method = "Rnd.EM")
ret.4 <- init.EM(x, nclass = 4, min.n = min.n, method = "Rnd.EM")

(lmt.23 <- lmt(ret.2, ret.3, x))
(lmt.34 <- lmt(ret.3, ret.4, x))
(lmt.24 <- lmt(ret.2, ret.4, x))
}
}
\keyword{programming}
