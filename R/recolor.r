rematch <- function(tg.id, cl.id) {
  ##
  ## This function remaps the colors in the classification vector cl.id to the
  ## target class vector given by tg.id. It is the heart of the recolor function
  ## and is usually called from recolor. However, it can easily be used as a
  ## standalone function.
  ##
  ## written Ranjan Maitra, Ames, IA 50011-1210, 2015/10/26

  which.max.matrix <- function(mat) (which(x = mat == max(mat), arr.ind=T))
  
  cl.id.tmp <- recode(cl.id) - min(cl.id) + 1
  tg.id.tmp <- recode(tg.id) - min(tg.id) + 1
  orig.cl.id <- sort(unique(cl.id))
  orig.tg.id <- sort(unique(tg.id))
  xtabs <- table(tg.id.tmp, cl.id.tmp)
  ncl <- max(cl.id.tmp)
  ntg <- max(tg.id.tmp)
  id.tg <- NULL
  id.cl <- NULL
  tg.idx <- 1:ntg 
  cl.idx <- 1:ncl 
  for (i in 1:min(c(ncl,ntg))) {
    xtab <- xtabs[tg.idx, cl.idx]
    if (!is.null(dim(xtab))) {
      ind <- which.max.matrix(xtab)[1,]
      id.tg <- c(id.tg, tg.idx[ind[1]])
      id.cl <- c(id.cl, cl.idx[ind[2]])
      tg.idx <- setdiff(tg.idx, tg.idx[ind[1]])
      cl.idx <- setdiff(cl.idx, cl.idx[ind[2]])
    }
    else {
      if (ncl == ntg) {
        id.cl <- c(id.cl, cl.idx)
        tg.cl <- c(tg.cl, tg.idx)
      }
      else {
        ind <- which.max(xtab)
        if (ncl > ntg) {
          id.cl <- c(id.cl,cl.idx[ind])
          id.tg <- c(id.tg,tg.idx)
          cl.idx <- setdiff(cl.idx, cl.idx[ind])
          id.cl <- c(id.cl, cl.idx)
        }
        else {
          id.cl <- c(id.cl,cl.idx)
          id.tg <- c(id.tg,tg.idx[ind])
          tg.idx <- setdiff(tg.idx, tg.idx[ind])
          id.tg <- c(id.tg, tg.idx)
        }
      }
    }
  }
  return(list(id.class = orig.cl.id[id.cl], id.target = orig.tg.id[id.tg]))
}

recode <- function(id) {
  ##
  ## This function reoders classes to eliminate group ids without any members.
  ## It is assumed that the group ids are integers
  ## Written Ranjan Maitra, Ames, IA 50011-1210, 2015/10/26
  ##

  cl.sort <- sort(unique(id))
  if (length(cl.sort) < 1 + diff(range(cl.sort))) {
    j <- min(id)
    for (i in 1:length(cl.sort)) {
      id[id==cl.sort[i]] <- j
      j <- j + 1
    }
  }
  return(id)
}

recolor <- function(id.target, id.class, scatter.class = NULL, scatter.target = NULL) {
  ##
  ## This function colors id.target in accordance with the most likely candidate
  ## in id.class. It returns a list as id.trcl (which is a factor version of 
  ## id.class) and id.prcl (which is a factored version of the colored id.target)
  ## Note that if scatter is present, then the class given by 0 is represented
  ## as scatter and it is assumed to be the same for both classifications.
  ##
  ## written Ranjan Maitra, Ames, IA 50011-1210, 2015/10/19
  ##

  ## first erase missing classes
  id.cl <- recode(id.class)
  id.tg <- recode(id.target)
  tg.id <- id.tg
  cl.id <- id.cl
  if (!is.null(scatter.target) | !is.null(scatter.class)) {
    tg.id <- id.tg[(id.tg != scatter.target) & (id.cl != scatter.class)]
    cl.id <- cl.id[(id.tg != scatter.target) & (id.cl != scatter.class)]
  }
  cls <- rematch(tg.id, cl.id)
  tg.ids <- cls$id.target
  cl.ids <- cls$id.class
  for (i in 1:min(c(length(tg.ids),length(cl.ids))))  {
    id.cl[id.cl == tg.ids[i]] <- -1
    id.cl[id.cl == cl.ids[i]] <- tg.ids[i]
    id.cl[id.cl == -1] <- cl.ids[i]
    j <- i+1
    while (j <= min(c(length(tg.ids),length(cl.ids)))) {
      if (cl.ids[j] == tg.ids[i]) {
        cl.ids[j] <- cl.ids[i]
        j <-  min(c(length(tg.ids),length(cl.ids)))
      }
      else
        j <- j + 1
    }
  }
  return(id.cl)
}
