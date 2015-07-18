Jaccard.Index <- function(x, y){
  xx <- as.vector(x)
  yy <- as.vector(y)

  if(length(xx) != length(yy)){
    return("Error: x and y have different numbers of elements")
  }

  id <- (!is.na(xx)) & (!is.na(yy))
  xx <- xx[id]
  yy <- yy[id]
  tmpx <- xx
  tmpy <- yy

  xx[tmpx == 1] <- 0
  xx[tmpx != 1] <- 1

  yy[tmpy == 1] <- 0
  yy[tmpy != 1] <- 1

  sum(xx == yy) / (length(xx) - sum((1 - xx) * (1 - yy)))
} # End of Jaccard.Index().
