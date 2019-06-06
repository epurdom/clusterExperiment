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
#' @param subsample logical as to whether to subsample via 
#'   \code{\link{subsampleClustering}} to get the distance matrix at each 
#'   iteration; otherwise the distance matrix is set by arguments to
#'   \code{\link{mainClustering}}.
#' @param beta value between 0 and 1 to decide how stable clustership membership
#'   has to be before 'finding' and removing the cluster.
#' @param top.can only the top.can clusters from \code{\link{mainClustering}} (ranked
#'   by 'orderBy' argument given to \code{\link{mainClustering}}) will be compared
#'   pairwise for stability. Making this very big will effectively remove this
#'   parameter and all pairwise comparisons of all clusters found will be
#'   considered. This might result in smaller clusters being found. The current
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
#' @param mainClusterArgs list of arguments to be passed to
#'   \code{\link{mainClustering}}).
#' @inheritParams clusterSingle
#' @details \code{seqCluster} is not meant to be called by the user. It is only
#'   an exported function so as to be able to clearly document the arguments for
#'   \code{seqCluster} which can be passed via the argument \code{seqArgs} in
#'   functions like \code{\link{clusterSingle}} and \code{\link{clusterMany}}.
#' @details This code is adapted from the sequential protion of the code of the
#'   tightClust package of Tseng and Wong. At each iteration of the algorithm it
#'   finds a set of samples that constitute a homogeneous cluster and remove
#'   them, and iterate again to find the next set of samples that form a
#'   cluster.
#' @details In each iteration, to determine the next set of homogeneous set of
#'   samples, the algorithm will iteratively cluster the current set of samples
#'   for a series of increasing values of the parameter $K$, starting at a value
#'   \code{kinit} and increasing by 1 at each iteration, until a sufficiently
#'   homogeneous set of clusters is found. For the first set of homogeneous
#'   samples, \code{kinit} is set to the argument $k0$, and for iteration,
#'   \code{kinit} is increased internally.
#' @details Depending on the value of \code{subsample} how the value of $K$ is
#'   used differs. If \code{subsample=TRUE}, $K$ is the \code{k} sent to the
#'   cluster function \code{clusterFunction} sent to 
#'   \code{\link{subsampleClustering}} via \code{subsampleArgs}; then
#'   \code{\link{mainClustering}} is run on the result of the co-occurance matrix from
#'   \code{\link{subsampleClustering}} with the \code{ClusterFunction} object
#'   defined in the argument \code{clusterFunction} set via \code{mainClusterArgs}.
#'   The number of clusters actually resulting from this run of
#'   \code{\link{mainClustering}} may not be equal to the $K$ sent to  the clustering
#'   done in \code{\link{subsampleClustering}}. If \code{subsample=FALSE},
#'   \code{\link{mainClustering}} is called directly on the data to determine the
#'   clusters and $K$ set by \code{seqCluster} for this iteration determines the
#'   parameter of the clustering done by \code{\link{mainClustering}}. Specifically,
#'   the argument \code{clusterFunction} defines the clustering of the
#'   \code{\link{mainClustering}} step and \code{k} is sent to that
#'   \code{ClusterFunction} object. This means that if \code{subsample=FALSE},
#'   the \code{clusterFunction} must be of \code{algorithmType} "K".
#' @details In either setting of \code{subsample}, the resulting clusters from
#'   \code{\link{mainClustering}} for a particular $K$ will be compared to clusters
#'   found in the previous iteration of $K-1$. For computational (and other?)
#'   convenience, only the first \code{top.can} clusters of each iteration will
#'   be compared to the first \code{top.can} clusters of previous iteration for
#'   similarity (where \code{top.can} currently refers to ordering by size, so
#'   first \code{top.can} largest clusters.
#' @details If there is no cluster of the first \code{top.can} in the current
#'   iteration $K$ that has overlap similarity > \code{beta} to any in the
#'   previous iteration, then the algorithm will move to the next iteration,
#'   increasing to $K+1$.
#'   
#' @details If, however, of these clusters there is a cluster in the current
#'   iteration $K$ that has overlap similarity > beta to a cluster in the
#'   previous iteration $K-1$, then the cluster with the largest such similarity
#'   will be identified as a homogenous set of samples and the samples in it
#'   will be removed and designated as such. The algorithm will then start again
#'   to determine the next set of homogenous samples, but without these samples.
#'   Furthermore, in this case (i.e. a cluster was found and removed), the value
#'   of \code{kinit} will be be reset to \code{kinit-1}; i.e. the range of
#'   increasing $K$ that will be iterated over to find a set of homogenous
#'   samples will start off one value less than was the case for the previous
#'   set of homogeneous samples. If \code{kinit-1}<\code{k.min}, then
#'   \code{kinit} will be set to \code{k.min}.
#'   
#'   
#' @details If there are less than \code{remain.n} samples left after finding a
#'   cluster and removing its samples, the algorithm will stop, as subsampling
#'   is deamed to no longer be appropriate. If the K has to be increased to
#'   beyond \code{k.max} without finding any pair of clusters with overlap >
#'   beta, then the algorithm will stop. Any samples not found as part of a
#'   homogenous set of clusters at that point will be classified as unclustered
#'   (given a value of -1)
#' @details Certain combinations of inputs to \code{mainClusterArgs} and
#'   \code{subsampleArgs} are not allowed. See \code{\link{clusterSingle}} for
#'   these explanations.
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
#' @seealso tight.clust,
#'   \code{\link{clusterSingle}},\code{\link{mainClustering}},\code{\link{subsampleClustering}}
#'   
#' @examples
#' \dontrun{
#' data(simData)
#'
#' set.seed(12908)
#' clustSeqHier <- seqCluster(simData, k0=5, subsample=TRUE,
#' beta=0.8, subsampleArgs=list(resamp.n=100,
#' samp.p=0.7, clusterFunction="kmeans", clusterArgs=list(nstart=10)),
#' mainClusterArgs=list(minSize=5,clusterFunction="hierarchical01",clusterArgs=list(alpha=0.1)))
#' }
#' @export
#' @rdname seqCluster
#' @export
seqCluster<-function(x=NULL, diss=NULL, cat=NULL, k0,  
                     subsample=TRUE, beta, top.can = 5, remain.n = 30, k.min = 3, 
                     k.max=k0+10,verbose=TRUE, subsampleArgs=NULL,mainClusterArgs=NULL,checkDiss=FALSE)
{
  ########
  ####Checks
  ########
  
  checkOut<-.checkArgs(x=x,diss=diss,cat=cat, subsample=subsample,sequential=TRUE,mainClusterArgs=mainClusterArgs,subsampleArgs=subsampleArgs,checkDiss=checkDiss)
  if(is.character(checkOut)) stop(checkOut)
  else {
    mainClusterArgs<-checkOut$mainClusterArgs
    subsampleArgs<-checkOut$subsampleArgs
    input<-checkOut$inputClusterD
  }		
  #Note, for the purpose of this
  N <- switch(input, "X"=dim(x)[2], "diss"=dim(diss)[2],"cat"=dim(cat)[2])
  if(input %in% c("X"))  colnames(x) <- as.character(seq_len(N))
  if(input %in% c("cat"))  colnames(cat) <- as.character(seq_len(N))
  if(input %in% c("diss")) colnames(diss)<-rownames(diss)<-as.character(seq_len(N))
  
  #########
  #Function that does the actual clustering steps
  # (Called at each iteration)
  #########
  updateClustering<-function(newk){
    if(verbose) cat(paste("k =", newk,"\n"))
    if(subsample){
      tempArgs<-subsampleArgs
      tempArgs[["clusterArgs"]]<-c(list(k=newk), subsampleArgs[["clusterArgs"]]) #set k  
      #also set the k for the mainClustering to be the same as in subsampling.
      tempClusterDArgs<-mainClusterArgs
      tempClusterDArgs[["clusterArgs"]] <- c(list(k=newk), mainClusterArgs[["clusterArgs"]])
      
      res <- .clusterWrapper(x=x, diss=diss, cat=cat, subsample=subsample,  subsampleArgs=tempArgs, mainClusterArgs=tempClusterDArgs)$results
    }
    else{
      tempArgs<-mainClusterArgs
	  #set k:
	  tempArgs[["clusterArgs"]]<-c(list(k=newk), mainClusterArgs[["clusterArgs"]])       res <- .clusterWrapper(x=x, diss=diss, cat=cat, subsample=subsample,  subsampleArgs=subsampleArgs, mainClusterArgs=tempArgs)$results
      
    }
    return(res)
  }
  
  
  #########
  ## Iteration code
  #########
  
  ###---------
  ###---------
  ###The following is legacy of tight.clust. They originally had programmed ability to look across more than 2 at each step to determing the stability of a cluster. This was not what they described in paper, and function is hard-coded at 2, but I have left code here in case we ever wanted to reconsider this issue.
  ###---------
  seq.num<-2
  kReturn<-"last" # when look at stability, return stable as first or last? For seq.num=2, not really matter, take last like paper
  kReturn<-match.arg(kReturn,c("last","first"))
  betaNum<-"all"
  betaNum<-match.arg(betaNum,c("all","last","first"))
  #This makes all combinations of 1:top.can, seq.num times (could be simplified if seq.num=2):
  #a ncombinations x seq.num matrix -- each row gives a combination of clusters to compare stability
  index.m <- as.matrix(expand.grid(lapply(seq_len(seq.num), function(x) seq_len(top.can))))
  whReturn<-switch(kReturn,"last"=seq.num,"first"=1) #way to index which one gets returned.
  ###---------
  ###---------
  
  #---------
  #iterative setup
  #-------
  remain <- N #keep track of how many samples not yet clustered (stop when less than remain.n)
  nfound <- 0 #keep track of how many clusters found/removed so far
  found <- TRUE #has a cluster been found/removed in last iteration
  k.start <- k0 #the starting k for the next cluster
  k <- k0
  ### Blank values that will be filled in.
  candidates <- list() #list of length seq.num of possible clusters found for each k to be compared
  tclust <- list() #list of final cluster identifications (indices of rows of x)
  kstart<-c() #the starting k for the cluster
  kend<-c() #the ending k for the cluster
  whyStop<-NULL
  #---------
  #iterative loop
  #-------
  while (remain >= remain.n && (found || k <= k.max)) {
    if (found) { 
	  #i.e. previous iteration found cluster; 
	  #need to start finding new cluster
      if(verbose) cat(paste("Looking for cluster", nfound + 1, "...\n"))
      k <- k.start
      currentStart<-k.start #will add this to kstart if successful in finding cluster
      #find clusters for K,K+1
      for (i in seq_len(seq.num)) {
        newk<-k + i - 1
        res<-updateClustering(newk)
        if(length(res)>0) res <- res[seq_len(min(top.can,length(res)))]
        candidates[[i]]<-res
      }
    }
    else { 
	  #previous iteration clustering wasn't good enough
	  #need to go increase to K+2,K+3, etc.
      candidates <- candidates[-1] #remove old k
      newk<-k + seq.num - 1
      if(verbose) cat(paste("k =", newk, "\n"))
      #add new k (because always list o)
      res<-updateClustering(newk)
      if(length(res)>0) res <- res[seq_len(min(top.can,length(res)))]
      candidates[[seq.num]] <- res
    }
    ##################
    #check whether all got top.can values for each -- could be less.
    #find which rows of index.m define cluster combinations that don't exist
    ##################
    nClusterPerK<-sapply(candidates,length) #number of clusters found per k sequence
    whInvalid<-unique(unlist(lapply(seq_len(ncol(index.m)), 
		function(i){ which(index.m[,i] > nClusterPerK[i])}
	)))
    if(length(whInvalid)==nrow(index.m)){
      #all invalid -- probably means that for some k there were no candidates found. So should stop.
      if(verbose) cat(paste("Found ",paste(nClusterPerK,collapse=","),"clusters for k=",paste(k+seq_len(seq.num)-1,collapse=","),", respectively. Stopping iterating because zero-length cluster.\n"))
      whyStop<-paste("Stopped in midst of searching for cluster",nfound+1," because no clusters meeting criteria found for iteration k=",k+i-1,"and previous clusters not similar enough.")
      break
    }
    if(length(whInvalid)>0){
      if(verbose) cat("Did not find", top.can,"clusters: ")
      if(verbose) cat(paste("found",paste(nClusterPerK,collapse=","),"clusters for k=",paste(k+seq_len(seq.num)-1,collapse=","),", respectively\n"))
      
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
      temp <- lapply(seq_len(seq.num), function(z) candidates[[z]][[y[z]]]) #each
      i.temp <- temp[[1]]
      if(betaNum %in% c("all","first")) u.temp <- temp[[1]] ###eap: changed here
      if(betaNum == "last") u.temp<-temp[[seq.num]]
      for (j in 2:seq.num) {
        i.temp <- intersect(i.temp, temp[[j]])
        if(betaNum=="all") u.temp <- union(u.temp, temp[[j]])
      }
      out<-length(i.temp)/length(u.temp)
      return(out)
      
    }
    beta.temp <- apply(tempIndex, 1, calc.beta) 
    if (any(beta.temp >= beta)){
      found <- TRUE
      nfound <- nfound + 1
      if(verbose) cat(paste("Cluster",nfound,"found."), "")
      if (k.start > k.min) k.start <- k.start - 1 #decrease
      found.temp <- candidates[[whReturn]][[tempIndex[which.max(beta.temp)[1], whReturn]]]
      kend[[nfound]]<-k+seq.num-1 #just assuming returning last here!
      kstart[[nfound]]<-currentStart
	  
	  tclust[[nfound]]<-switch(input,
		  "X"=colnames(x),
		  "diss"=colnames(diss),
		  "cat"=colnames(cat))[found.temp]
      tclust[[nfound]]<-as.numeric(tclust[[nfound]])
      if(input %in% c("X")) x <- x[,-found.temp] 
      if(input %in% c("diss")) diss<-diss[-found.temp,-found.temp]
	  if(input %in% c("cat")) cat<-cat[,-found.temp]
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
  ###############
  # Clean up and return results
  ###############
  if(is.null(whyStop)){
    if(remain< remain.n) whyStop<-"Ran out of samples"
    if(!found & k>k.max) whyStop<-paste("Went past k.max=",k.max,"in looking for cluster with similarity to previous.")
  }
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
