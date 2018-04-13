seqCluster2 <- function(x=NULL, diss=NULL, k0,
                     subsample=TRUE, beta, top.can = 5, remain.n = 30, k.min = 3,
                     k.max=k0+10, verbose=TRUE, subsampleArgs=NULL,
                     mainClusterArgs=NULL, checkDiss=TRUE) {


  ####Checks
  ####
  checkOut <- .checkSubsampleClusterDArgs(x=x, diss=diss,
                                          subsample=subsample, sequential=TRUE,
                                          mainClusterArgs=mainClusterArgs,
                                          subsampleArgs=subsampleArgs,
                                          checkDiss=checkDiss)

  if(is.character(checkOut)) {
    stop(checkOut)
  } else {
    mainClusterArgs <- checkOut$mainClusterArgs
    subsampleArgs <- checkOut$subsampleArgs
    input <- checkOut$inputClusterD
  }

  ###The following is legacy of tight.clust. They originally had programmed ability to look across more than 2 at each step to determing the stability of a cluster. This was not what they described in paper, and function is hard-coded at 2, but I have left code here in case we ever wanted to reconsider this issue.
  seq.num <- 2
  kReturn <- "last"
  ## when look at stability, return stable as first or last? For seq.num=2, not really matter, take last like paper
  kReturn <- match.arg(kReturn, c("last","first"))
  betaNum <- "all"
  betaNum <- match.arg(betaNum,c("all","last","first"))
  #This makes all combinations of 1:top.can, seq.num times (could be simplified if seq.num=2):
  #a ncombinations x seq.num matrix -- each row gives a combination of clusters to compare stability
  index.m <- as.matrix(expand.grid(lapply(seq(0, seq.num - 1), function(x) seq(0, top.can - 1))))
  whReturn <- switch(kReturn, "last"=seq.num-1, "first"=0) #way to index which one gets returned.

  if(input %in% c("X")) {
    N <- dim(x)[2]
  }

  if(input=="diss") {
    N <- dim(diss)[2]
  }

  if(verbose){
    if(input %in% c("X")) cat(paste("Number of points:", N, "\tDimension:", dim(x)[1], "\n"))
    else cat(paste("Number of points:", N,"\n"))
  }

  if(input %in% c("X")) {
    colnames(x) <- as.character(1:N)
  }
  if(input %in% c("diss")) {
    colnames(diss) <- rownames(diss) <- as.character(1:N)
  }

  updateClustering<-function(newk, x, diss){

    if(verbose) cat(paste("k =", newk,"\n"))
    if(subsample){
      tempArgs<-subsampleArgs
      tempArgs[["clusterArgs"]]<-c(list(k=newk), subsampleArgs[["clusterArgs"]]) #set k
      #also set the k for the mainClustering to be the same as in subsampling.
      tempClusterDArgs<-mainClusterArgs
      tempClusterDArgs[["clusterArgs"]] <- c(list(k=newk), mainClusterArgs[["clusterArgs"]])

      res <- .clusterWrapper(x=x, diss=diss,subsample=subsample,  subsampleArgs=tempArgs, mainClusterArgs=tempClusterDArgs)$results
    }
    else{
      tempArgs<-mainClusterArgs
      tempArgs[["clusterArgs"]]<-c(list(k=newk), mainClusterArgs[["clusterArgs"]]) #set k
      res <- .clusterWrapper(x=x, diss=diss, subsample=subsample,  subsampleArgs=subsampleArgs, mainClusterArgs=tempArgs)$results

    }
    return(res)
  }

  # if(is.null(diss)) {
  #   diss <- matrix(0)
  # }
  #
  # if(is.null(x)) {
  #   x <- matrix(0)
  # }

dd <- as.matrix(dist(t(x)))
  retval <- do_seq_cluster(x, dd, N, k0, remain.n, seq.num, k.min, k.max,
                           top.can, beta, index.m, betaNum, whReturn,
                           input, updateClustering)

  if(retval$whyStop == 1) {
    whyStop <- "Stopped in midst of searching for cluster because no clusters meeting criteria found in last iteration and previous clusters not similar enough."
  } else {
    if(retval$remain< retval$remain.n) whyStop<-"Ran out of samples"
    if(!retval$found & retval$k>retval$k.max) whyStop<-paste("Went past k.max=",retval$k.max,"in looking for cluster with similarity to previous.")

  }

  clusterVector<-.convertClusterListToVector(retval$tclust,N)
  if(all(clusterVector==-1) & length(retval$tclust)>0) stop("coding error")
  if(retval$nfound>0){
    size <- sapply(retval$tclust, length)
    sizeMat<-cbind(size=size,
                   kStart=retval$kstart,
                   kEnd=retval$kend,
                   nIter=retval$kend-retval$kstart)
    res <- list(clustering = clusterVector, clusterInfo = sizeMat, whyStop=whyStop)
    if(verbose) cat(paste("Stopped because:", whyStop),"")
    return(res)
  } else{
    if(verbose) cat("No tight clusters could be found with given parameters")
    return(list(clustering = clusterVector, whyStop=whyStop))
  }

}
