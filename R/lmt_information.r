
# This function creates G matrix needed for taking
# derivatives of symmetric matrices
# Namely, vec(X) = G %*% vech(X)
Gmat <- function(p){
  M <- p * (p + 1) / 2
  G <- matrix(NA, ncol = M, nrow = p * p)
  n <- 1

  Z <- rep(0, M)
  for(a in 1:p){
    for(b in 1:p){
      Zn <- Z
      i1 <- max(a, b)      
      i2 <- min(a, b)
      ind <- M - (p - i2 + 1) * (p - i2 + 2) / 2 + i1 - i2 + 1
      Zn[ind] <- 1
      G[n,] <- Zn
      n <- n + 1
    }
  }
  return(G)
} # End of Gmat().


### partial logL. Return a matrix with dimension M * N.
partial.logL <- function(x, PI, MU, S, t){
  K <- length(PI)
  N <- dim(x)[1]
  p <- dim(x)[2]

  mp <- p * (p + 1) / 2
  M <- K - 1 + K * p + K * p * (p + 1) / 2
  G <- Gmat(p)
  Id <- diag(p)

  invS <- apply(S, 3, solve)
  dim(invS) <- c(p, p, K)

  SS <- lapply(1:N, function(i){
    Sbuf.1 <- NULL
    if(K != 1){
      tmp <- t[i, ] / PI
      Sbuf.1 <- tmp[1:(K-1)] - tmp[K]
    }

    Sbuf.2 <- list()
    Sbuf.3 <- list()
    for(j in 1:K){
      tmp <- x[i,] - MU[j,]
      tmp.1 <- t[i,j] * invS[,,j]
      Sbuf.2[[j]] <- tmp.1 %*% tmp
      Smat <-  tmp.1 %*% (tmp %*% t(tmp) %*% invS[,,j] - Id) / 2
      Sbuf.3[[j]] <- as.vector(Smat) %*% G
    }

    c(Sbuf.1, do.call("c", Sbuf.2), do.call("c", Sbuf.3))
  })

  do.call("cbind", SS)
} # End of partial.logL().


### Back compartiable with Volodymyr's code.
Iy <- function(x, PI, MU, S, t){
  SS <- partial.logL(x, PI, MU, S, t)
  SS %*% t(SS)
} # End of Iy().

Iy2 <- function(x, PI0, MU0, S0, t0, PIa, MUa, Sa, ta){
  SS <- rbind(partial.logL(x, PI0, MU0, S0, t0),
              partial.logL(x, PIa, MUa, Sa, ta))
  SS %*% t(SS)
} # End of Iy2().


### Compute posterior.
postPI <- function(x, emobj){
  e.step(x, emobj)$Gamma
} # End of postPI().


### Generate dataset.
GenDataSet <- function(N, PI, MU, S){
  K <- dim(MU)[1]
  p <- dim(MU)[2]
  Nk <- drop(rmultinom(1, N, PI))

  id <- rep(1:K, Nk)
  Sample <- lapply(1:K, function(i){
    ### To avoid degeneration.
    if(Nk[i] > 0){
      ret <- mvrnorm(n = Nk[i], mu = MU[i,], Sigma = S[,,i])
    } else{
      ret <- NULL
    }
    ret
  })
  x <- do.call("rbind", Sample)

  list(x = x, id = id)
} # End of GenDataSet().

GenMixDataSet <- function(N, PI0, MU0, S0, PIa, MUa, Sa, tau = 0.5){
  N0 <- rbinom(1, N, c(tau, 1-tau))
  Na <- N - N0

  ret0 <- GenDataSet(N0, PI0, MU0, S0)
  reta <- GenDataSet(Na, PIa, MUa, Sa)

  list(x = rbind(ret0$x, reta$x), id = c(ret0$id, reta$id),
       hid = c(rep(0, N0), rep(1, Na)))
} # End of GenMixDataSet().


### Observed Information for DataSet.
w.2 <- function(x, emobj0, emobja, tau = 0.5){
  f0 <- sum(log(dmixmvn(x, emobj0))) + log(tau)
  fa <- sum(log(dmixmvn(x, emobja))) + log(1 - tau)
  ### Numerical unstable.
  # g <- exp(f0) + exp(fa)
  # c(exp(f0) / g, exp(fa) / g)
  pi0 <- 1 / (exp(fa - f0) + 1)
  c(pi0, 1 - pi0)
} # End of w.2().


### Obtain parameters.
get.E.chi2 <- function(x, emobj0, emobja, given = c("0", "a"), tau = 0.5,
    n.mc = 1000, verbose = TRUE){
  N <- nrow(x)
  p <- ncol(x)

  K0 <- emobj0$nclass
  PI0 <- emobj0$pi
  MU0 <- emobj0$Mu
  S0 <- LTSigma2variance(emobj0$LTSigma) 

  Ka <- emobja$nclass
  PIa <- emobja$pi
  MUa <- emobja$Mu
  Sa <- LTSigma2variance(emobja$LTSigma) 

  if(given[1] == "0"){
    PI <- PI0
    MU <- MU0
    S <- S0
  } else if(given[1] == "a"){
    PI <- PIa
    MU <- MUa
    S <- Sa
  } else{
    stop("given should be '0' or 'a'.")
  }

  ### Obtain nabla logL via Monte Carlo.
  x.new <- GenDataSet(n.mc, PI, MU, S)$x
  t0 <- postPI(x.new, emobj0)
  ta <- postPI(x.new, emobja)
  pl0 <- partial.logL(x.new, PI0, MU0, S0, t0)
  pla <- partial.logL(x.new, PIa, MUa, Sa, ta)

  ### Expected degrees of freedom.
  nl <- rbind(pl0, pla)
  mu <- rowMeans(nl)
  nl <- nl - mu
  J <- nl %*% t(nl) / n.mc
  par.df <- eigen(J, TRUE, only.values = TRUE)$values
  par.df <- par.df[par.df > 0]
  par.df <- sum((par.df / par.df[1] > 1e-10) |
                (cumsum(par.df) / sum(par.df) < 0.90))

  ### Expected noncenteriality
  M0 <- nrow(pl0)
  Ma <- nrow(pla)
  if(given == "0"){
    id <- M0 + (1:Ma)
  } else{
    id <- 1:M0
  }
  J.nc <- matrix(J[id, id], nrow = length(id))
  nu <- matrix(mu[id], nrow = length(id))
  ### Numerical unstable.
  # par.nc <- t(nu) %*% solve(J.nc) %*% nu
  tmp <- eigen(J.nc, TRUE)
  tmp.nc <- tmp$values[tmp$values > 0]
  tmp.nc <- sum((tmp.nc / tmp.nc[1] > 1e-8) |
                (cumsum(tmp.nc) / sum(tmp.nc) < 0.90))
  nu <- t(nu) %*% tmp$vectors[, 1:tmp.nc]
  par.nc <- nu %*% diag(1 / tmp$values[1:tmp.nc], tmp.nc, tmp.nc) %*% t(nu)

  ### For returns.
  ret <- c(par.df, par.nc)
  if(verbose){
    cat("K0=", K0, ", M0=", M0, " v.s. Ka=", Ka, ", Ma=", Ma,
        " | given=", given, " : df=", ret[1], ", nc=", ret[2],
        ".\n", sep = "")
  }
  ret
} # End of get.E.chi2().

get.E.delta <- function(x, emobj0, emobja, tau = 0.5, n.mc = 1000){
  N <- nrow(x)

  PI0 <- emobj0$pi
  MU0 <- emobj0$Mu
  S0 <- LTSigma2variance(emobj0$LTSigma) 

  PIa <- emobja$pi
  MUa <- emobja$Mu
  Sa <- LTSigma2variance(emobja$LTSigma) 

  ### This is inaccurate and incorrect!!!
  # E.delta <- lapply(1:n.mc, function(i){
  #   x.new <- GenMixDataSet(N, PI0, MU0, S0, PIa, MUa, Sa, tau = tau)$x
  #   logL(x.new, emobja) - logL(x.new, emobj0)
  # })

  n.mc.0 <- rbinom(1, n.mc, c(tau, 1-tau))
  n.mc.a <- n.mc - n.mc.0

  E.delta <- NULL
  if(n.mc.0 > 0){
    tmp <- lapply(1:n.mc.0, function(i){
      x.new <- GenDataSet(N, PI0, MU0, S0)$x
      logL(x.new, emobja) - logL(x.new, emobj0)
    })
    E.delta <- c(E.delta, tmp)
  }
  if(n.mc.a > 0){
    tmp <- lapply(1:n.mc.a, function(i){
      x.new <- GenDataSet(N, PIa, MUa, Sa)$x
      logL(x.new, emobja) - logL(x.new, emobj0)
    })
    E.delta <- c(E.delta, tmp)
  }

  do.call("sum", E.delta) / n.mc
} # End of get.E.delta().


pchisq.my <- function(q, df, ncp = 0, lower.tail = TRUE){
  if(ncp == Inf){
    if(lower.tail){
      ret <- 0
    } else{
      ret <- 1
    }
  } else{
    ret <- pchisq(q, df, ncp, lower.tail)
  }
  ret
} # End of pchisq.my().
