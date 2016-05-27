#' Program for sequentially clustering, removing cluster, and starting again.
#'
#' Given a data matrix, this function will call clustering
#' routines, and sequentially remove best clusters, and iterate to find
#' clusters.
#'
#' @param x \code{p x n} data matrix on which to run the clustering (samples in
#'   columns).
#' @param diss \code{n x n} data matrix of dissimilarities between the samples
#'   on which to run the clustering
#' @param k0 the value of K at the first iteration of sequential algorithm, see
#'   details below or vignette.
#' @param clusterFunction passed to clusterDMat option 'clusterFunction' to
#'   indicate method of clustering, see \code{\link{clusterD}}.
#' @param subsample logical as to whether to subsample via 
#'   \code{\link{subsampleClustering}} to get the distance matrix at each 
#'   iteration; otherwise the distance matrix is set by arguments to
#'   \code{\link{clusterD}}.
#' @param beta value between 0 and 1 to decide how stable clustership membership
#'   has to be before 'finding' and removing the cluster.
#' @param top.can only the top.can clusters from \code{\link{clusterD}} (ranked
#'   by 'orderBy' argument given to \code{\link{clusterD}}) will be compared
#'   pairwise for stability. Making this very big will effectively remove this
#'   parameter and all pairwise comparisons of all clusters found will be
#'   considered. This might result in smaller clusters being found. Current
#'   default is fairly large, so probably will have little effect.
#' @param remain.n when only this number of samples are left (i.e. not yet
#'   clustered) then algorithm will stop.
#' @param k.min each iteration of sequential detection of clustering will
#'   decrease the beginning K of subsampling, but not lower than k.min.
#' @param k.max algorithm will stop if K in iteration is increased beyond this
#'   point.
#' @param verbose whether the algorithm should print out information as to its
#'   progress.
#' @param subsampleArgs list of arguments to be passed to
#'   \code{\link{subsampleClustering}}.
#' @param clusterDArgs list of arguments to be passed to
#'   \code{\link{clusterD}}(which can include arguments to be passed to
#'   \code{\link{cluster01}} or \code{\link{clusterK}}).
#'
#' @details Each iteration of the algorithm will cluster the current set of
#'   samples. Depending on the method, the number of clusters resulting from
#'   \code{\link{clusterD}} may not be equal to the K used in the clustering of
#'   the (subsampled) data. The resulting clusters will then be compared to
#'   clusters found in the previous iteration that set the subsampling
#'   clustering to K-1. For computational (and other?) convenience, only the
#'   first top.can clusters of each iteration will be compared to the first
#'   top.can clusters of previous iteration for similarity (where top.can
#'   currently refers to ordering by size, so first top.can largest clusters).
#'
#' @details If there is a cluster in the current iteration that has overlap
#'   similarity > beta to a cluster in the previous iteration, then the cluster
#'   with the largest such similarity will be identified as a 'final' cluster
#'   and the samples in it will be removed for future iterations. The algorithm
#'   will then continue to the next iteration, but without these samples.
#'   Furthermore, in this case K for the next iteration will NOT be set to K+1,
#'   but will be reset to kinit-1, where kinit was the first K used after the
#'   previous 'final' cluster was removed. If kinit-1<k.min, then K will be set
#'   to k.min.
#'
#' @details If there is no cluster of the first top.can in the current iteration
#'   that has overlap similarity > beta to any in the previous iteration, then
#'   the algorithm will move to the next iteration (i.e. redo after increasing K
#'   to K+1).
#'
#' @details If there are less than remain.n samples left after finding a cluster
#'   and removing its samples, the algorithm will stop, as subsampling is deamed
#'   to no longer be appropriate. If the K has to be increased to beyond k.max
#'   without finding any pair of clusters with overlap > beta, then the
#'   algorithm will stop. Any samples not found as part of a 'final' cluster
#'   after the algorithm stops, will be classified as unclustered (given a value
#'   of -1)
#'
#' @details 'subsample' controls what is the D (distance) matrix used for
#'   clustering at each iteration. If subsample=TRUE, D is given via
#'   \code{\link{subsampleClustering}} function with k=K (with additional
#'   arguments passed via subsampleArgs). If subsample=FALSE, D is dist(x), for
#'   the samples currently considered in the iteration and clusterFunction must
#'   be of the 'K' type (e.g. "pam", see \code{\link{clusterD}}) or an error
#'   will be produced. The nsample x nsample matrix D is then clustered via
#'   \code{\link{clusterD}} to find clusters. The option 'clusterFunction' is
#'   passed to the argument 'clusterFunction' of \code{\link{clusterD}} to
#'   control what method is used to cluster D.
#'
#' @details If clusterFunction is of type 'K' (e.g. "pam", see
#'   \code{\link{clusterD}}) the 'k' argument of \code{\link{clusterK}} called
#'   by \code{\link{clusterD}} is set to the current iteration of K by the
#'   sequential iteration, so setting 'k=' in the list given to clusterDArgs
#'   will not do anything and will produce a warning to that effect.
#'
#' @details Similarly, the current K of the iteration also determines the 'k'
#'   argument passed to \code{\link{subsampleClustering}}  so setting 'k=' in
#'   the list given to the subsampleArgs will not do anything and will produce a
#'   warning to that effect.
#'
#' @details If subsample=FALSE and 'findBestK=FALSE' is passed to clusterDArgs,
#'   then each iteration will run the clustering given by clusterFunction on
#'   dist(x) iterating over k. However, if subsample=FALSE, you should not set
#'   'findBestK=TRUE' (otherwise clustering dist(x) will be essentially the same
#'   for iterating over different k and there is no method implemented to change
#'   the choice of how to remove a cluster other than similarity as you change
#'   k); an error message will be given if this combination of options are set.
#'
#' @details However, if clusterFunction="pam" (or is of type 'K') and
#'   subsample=TRUE passing either 'findBestK=TRUE' or 'findBestK=FALSE' will
#'   function as expected. In particular, the iteration over K will set the
#'   number of clusters for clustering of each subsample. If findBestK=FALSE,
#'   that same K will be used for clustering of DMat. If findBestK=TRUE, then
#'   \code{\link{clusterD}} will search for best k; note that the default
#'   'kRange' over which \code{\link{clusterD}} searches when findBestK=TRUE
#'   depends on the input value of 'k' (you can change this to a fixed set of
#'   values by setting 'kRange' explicitly in the clusterDArgs list).
#'
#' @return A list with values
#' \itemize{
#'
#' \item{\code{clustering}}{ a vector of length equal to nrows(x) giving the
#' integer-valued cluster ids for each sample. The integer values are assigned
#' in the order that the clusters were found. "-1" indicates the sample was not
#' clustered.}
#'
#' \item{\code{clusterInfo}}{ if clusters were successfully found, a matrix of
#' information regarding the algorithm behavior for each cluster (the starting
#' and stopping K for each cluster, and the number of iterations for each
#' cluster).}
#'
#' \item{\code{whyStop}}{ a character string explaining what triggered the
#' algorithm to stop.}
#' }
#' @references Tseng and Wong (2005), "Tight Clustering: A Resampling-Based
#'   Approach for Identifying Stable and Tight Patterns in Data", Biometrics,
#'   61:10-16.
#' 
#' @examples
#' \dontrun{
#' data(simData)
#'
#' set.seed(12908)
#'
#' clustSeqHier <- seqCluster(t(simData), k0=5, subsample=TRUE,
#' clusterFunction="hierarchical01", beta=0.8, subsampleArgs=list(resamp.n=100,
#' samp.p=0.7, clusterFunction="kmeans", clusterArgs=list(nstart=10)),
#' clusterDArgs=list(minSize=5))
#' }
#' @export
seqCluster<-function (x=NULL, diss=NULL, k0, clusterFunction=c("tight","hierarchical01","pam","hierarchicalK"), subsample=TRUE,beta = 0.7, top.can = 15, remain.n = 30, k.min = 2, k.max=k0+10,verbose=TRUE, subsampleArgs=NULL,clusterDArgs=NULL)
{
  input<-.checkXDissInput(x,diss)
    #for now, if use pam for subsampleClusterMethod, just use given k.
  if(!is.function(clusterFunction)){
    clusterFunction<-match.arg(clusterFunction)
    if(!is.function(clusterFunction)) typeAlg<-.checkAlgType(clusterFunction)
  }
  else{
    if(! "typeAlg" %in% clusterDArgs) stop("if you provide your own clustering algorithm to be passed to clusterD, then you must specify 'typeAlg' in clusterDArgs")
    else typeAlg<-clusterDArgs[["typeAlg"]]
  }
  if(typeAlg == "K"){
    if("findBestK" %in% names(clusterDArgs) & !subsample){
      if(clusterDArgs[["findBestK"]]) stop("Cannot do sequential clustering where subsample=FALSE and 'findBestK=TRUE' is passed via clusterDArgs. See help documentation.")
    }
    
  }
  ################
  ################
  ###The following is legacy of tight.clust. They originally had programmed ability to look across more than 2 at each step to determing the stability of a cluster. This was not what they described in paper, and function is hard-coded at 2, but I have left code here in case we ever wanted to reconsider this issue.
  seq.num<-2
  kReturn<-"last" # when look at stability, return stable as first or last? For seq.num=2, not really matter, take last like paper
  kReturn<-match.arg(kReturn,c("last","first"))
  betaNum<-"all"
  betaNum<-match.arg(betaNum,c("all","last","first"))
  #This makes all combinations of 1:top.can, seq.num times (could be simplified if seq.num=2):
  #a ncombinations x seq.num matrix -- each row gives a combination of clusters to compare stability
  index.m <- as.matrix(expand.grid(lapply(1:seq.num, function(x) 1:top.can)))
  whReturn<-switch(kReturn,"last"=seq.num,"first"=1) #way to index which one gets returned.
  ################
  ################
  if(input %in% c("X","both")) N <- dim(x)[2]
  if(input=="diss") N<-dim(diss)[2]
  if(verbose){
    if(input %in% c("X","both")) cat(paste("Number of points:", N, "\tDimension:", dim(x)[1], "\n"))
    else cat(paste("Number of points:", N,"\n"))
  }
#   if(input %in% c("X","both")){
#     original.data <- x
#     colnames(x) <- as.character(1:N)
#     id <- colnames(x)
#   }
#   else{
#     original.data <- diss
#     id<-colnames(diss)
#   }
  if(input %in% c("X","both"))  colnames(x) <- as.character(1:N)
  if(input %in% c("diss","both")) colnames(diss)<-rownames(diss)<-as.character(1:N)
  
  #iterative setup
  remain <- N #keep track of how many samples not yet clustered (stop when less than remain.n)
  nfound <- 0 #keep track of how many clusters found/removed so far
  found <- TRUE #has a cluster been found/removed in last iteration
  k.start <- k0 #the starting k for the next cluster
  k <- k0
  
  candidates <- list() #list of length seq.num of possible clusters found for each k to be compared
  tclust <- list() #list of final cluster identifications (indices of rows of x)
  kstart<-c() #the starting k for the cluster
  kend<-c() #the ending k for the cluster
  whyStop<-NULL
  if("k" %in% names(subsampleArgs)){
    #remove predefined versions of k from both.
    whK<-which(names(subsampleArgs)=="k")
    warning("Setting 'k' in subsampleArgs when the seqCluster is called will have no effect.")
    subsampleArgs<-subsampleArgs[-whK]
  }
  if("k" %in% names(clusterDArgs)){
    whK<-which(names(clusterDArgs)=="k")
    warning("Setting 'k' in clusterDArgs when the seqCluster is called will have no effect.")
    clusterDArgs<-clusterDArgs[-whK]
  }
  while (remain >= remain.n && (found || k <= k.max)) {
    if (found) { #i.e. start finding new cluster
      if(verbose) cat(paste("Looking for cluster", nfound + 1, "...\n"))
      k <- k.start
      currentStart<-k.start #will add this to kstart if successful in finding cluster
      
      #find clusters for K,K+1
      for (i in 1:seq.num) {
        if(verbose) cat(paste("k =", k + i - 1,"\n"))
        if(subsample){
          tempArgs<-c(list(k=k + i - 1),subsampleArgs) #set k
          res <- .clusterWrapper(x=x, subsample=subsample, clusterFunction=clusterFunction, subsampleArgs=tempArgs, clusterDArgs=clusterDArgs,typeAlg=typeAlg)$results
        }
        else{
          tempArgs<-c(list(k=k + i - 1),clusterDArgs) #set k
          res <- .clusterWrapper(x=x, diss=diss, subsample=subsample, clusterFunction=clusterFunction, subsampleArgs=subsampleArgs, clusterDArgs=tempArgs,typeAlg=typeAlg)$results
          
        }
        # if(length(res)==0) {
        # 					cat(paste("Found",paste(nClusterPerK,collapse=","),"clusters for k=",paste(k+1:seq.num-1,collapse=","),". Stopping because zero-length cluster.\n"))
        # 								whyStop<-paste("Stopped in midst of searching for cluster",nfound+1," because no clusters meeting criteria found for iteration k=",k+i-1,"and previous clusters not similar enough.")
        # 				}
        if(length(res)>0) res <- res[1:min(top.can,length(res))]
        candidates[[i]]<-res
      }
    }
    else { #need to go increase to K+2,K+3, etc.
      candidates <- candidates[-1] #remove old k
      if(verbose) cat(paste("k =", k + seq.num - 1, "\n"))
      #add new k (because always list o)
      if(subsample){
        tempArgs<-c(list(k=k + seq.num - 1),subsampleArgs)  #set k
        res <- .clusterWrapper(x=x, diss=diss, subsample=subsample, clusterFunction=clusterFunction, subsampleArgs=tempArgs, clusterDArgs=clusterDArgs,typeAlg=typeAlg)$results
      }
      else{
        tempArgs<-c(list(k=k + seq.num - 1),clusterDArgs) #set k
        res <- .clusterWrapper(x=x, diss=diss, subsample=subsample, clusterFunction=clusterFunction, subsampleArgs=subsampleArgs, clusterDArgs=tempArgs,typeAlg=typeAlg)$results
        
      }
      if(length(res)>0) res <- res[1:min(top.can,length(res))]
      candidates[[seq.num]] <- res
    }
    ##################
    #check whether all got top.can values for each -- could be less.
    #find which rows of index.m define cluster combinations that don't exist
    ##################
    nClusterPerK<-sapply(candidates,length) #number of clusters found per k sequence
    whInvalid<-unique(unlist(lapply(1:ncol(index.m),function(i){which(index.m[,i] > nClusterPerK[i])})))
    if(length(whInvalid)==nrow(index.m)){
      #all invalid -- probably means that for some k there were no candidates found. So should stop.
      if(verbose) cat(paste("Found ",paste(nClusterPerK,collapse=","),"clusters for k=",paste(k+1:seq.num-1,collapse=","),", respectively. Stopping iterating because zero-length cluster.\n"))
      whyStop<-paste("Stopped in midst of searching for cluster",nfound+1," because no clusters meeting criteria found for iteration k=",k+i-1,"and previous clusters not similar enough.")
      #browser()
      break
    }
    if(length(whInvalid)>0){
      if(verbose) cat("Did not find", top.can,"clusters: ")
      if(verbose) cat(paste("found",paste(nClusterPerK,collapse=","),"clusters for k=",paste(k+1:seq.num-1,collapse=","),", respectively\n"))
      
      tempIndex<-index.m[-whInvalid,,drop=FALSE]
    }
    else tempIndex<-index.m
    
    ##################
    #Calculate the stability pairwise between all of cluster combinations
    ##################
    #function to calculate the stability across the given combination
    #	y is a combination (row of index.m) giving clusters to compare stability from k, k+1
    calc.beta <- function(y) {
      #written generally enough to deal with seq.num>2; could be a lot simpler with seq.num=2.
      temp <- lapply(1:seq.num, function(z) candidates[[z]][[y[z]]]) #each
      i.temp <- temp[[1]]
      if(betaNum %in% c("all","first")) u.temp <- temp[[1]] ###eap: changed here
      if(betaNum == "last") u.temp<-temp[[seq.num]]
      for (j in 2:seq.num) {
        i.temp <- intersect(i.temp, temp[[j]])
        if(betaNum=="all") u.temp <- union(u.temp, temp[[j]])
      }
      out<-length(i.temp)/length(u.temp)
      #if(is.na(out)) stop("coding error: invalid similarity calculation")
      #else
      return(out)
      
    }
    beta.temp <- apply(tempIndex, 1, calc.beta) #original code had unlist. I removed it...might cause problems, but if so, should figure them out!
    if (any(beta.temp >= beta)){
      found <- TRUE
      nfound <- nfound + 1
      if(verbose) cat(paste("Cluster",nfound,"found."), "")
      if (k.start > k.min) k.start <- k.start - 1 #decrease
      found.temp <- candidates[[whReturn]][[tempIndex[which.max(beta.temp)[1], whReturn]]]
      kend[[nfound]]<-k+seq.num-1 #just assuming returning last here!
      kstart[[nfound]]<-currentStart
      if(input %in% c("X","both")) tclust[[nfound]] <- colnames(x)[found.temp] #need to do rownames, because remove rows from x
      else tclust[[nfound]] <- colnames(diss)[found.temp] #need to do rownames, because remove rows from x
      mode(tclust[[nfound]]) <- "numeric"
      if(input %in% c("X","both")) x <- x[-found.temp, ] 
      if(input %in% c("diss","both")) diss<-diss[-found.temp,-found.temp]
      remain <- remain - length(tclust[[nfound]])
      if(verbose) cat(paste("Cluster size:", length(tclust[[nfound]]),
                            "\tRemaining number of points:", remain, "\n"),
                      "")
    }
    else {
      found = FALSE
      k = k + 1
    }
  }
  if(is.null(whyStop)){
    if(remain< remain.n) whyStop<-"Ran out of samples"
    if(!found & k>k.max) whyStop<-paste("Went past k.max=",k.max,"in looking for cluster with similarity to previous.")
  }
  #browser()
  clusterVector<-.convertClusterListToVector(tclust,N)
  if(all(clusterVector==-1) & length(tclust)>0) stop("coding error")
  if(nfound>0){
    size <- sapply(tclust, length)
    sizeMat<-cbind(size=size,kStart=kstart,kEnd=kend,nIter=kend-kstart)
    res <- list(clustering = clusterVector, clusterInfo = sizeMat, whyStop=whyStop)
    if(verbose) cat(paste("Stopped because:", whyStop),"")
    return(res)
  }
  else{
    if(verbose) cat("No tight clusters could be found with given parameters")
    return(list(clustering = clusterVector, whyStop=whyStop))
  }
  
}
