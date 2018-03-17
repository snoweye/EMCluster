### This function utilize init.EM to obtain an initial.

simple.init <- function(x, nclass = 1){
  .emc <- .EMControl(short.iter = 1, em.iter = 0)
  tmp <- rand.EM(x, nclass, EMC = .emc, stable.solution = FALSE,
                 min.n.iter = 1)
  ret <- list(pi = tmp$pi, Mu = tmp$Mu, LTSigma = tmp$LTSigma)
  class(ret) <- "emret"
  ret
}
