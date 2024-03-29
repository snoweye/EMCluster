\name{EMCluster-package}
\alias{EMCluster-package}
\docType{package}
\title{EM Algorithm for Model-Based Clustering of Finite Mixture Gaussian
       Distribution}
\description{
  EMCluster provides EM algorithms and several efficient
  initialization methods for model-based clustering of finite
  mixture Gaussian distribution with unstructured dispersion
  in both of unsupervised and semi-supervised clustering.
}
\details{
  The install command is simply as \cr
  \cr
  \code{> R CMD INSTALL EMCluster_0.2-0.tar.gz } \cr
  \cr
  from a command mode or \cr
  \cr
  \code{R> install.packages("EMCluster") } \cr
  \cr
  inside an R session.
}
\references{
  \url{https://www.stat.iastate.edu/people/ranjan-maitra}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com} and Ranjan Maitra.
}
\seealso{
  \code{\link{init.EM}}, \code{\link{emcluster}}.
}
\examples{
\dontrun{
demo(allinit, 'EMCluster', ask = F, echo = F)
demo(allinit_ss, 'EMCluster', ask = F, echo = F)
}
}
\keyword{package}

