#' @title Cluster distance matrix from subsampling
#'
#' @description Given a \code{n x n} matrix of distances, these functions will
#'   try to find the clusters based on the given clustering function. cluster01
#'   and clusterK are internal functions and clusterD is a wrapper around these
#'   two functions for easier user interface. cluster01 and clusterK are not
#'   expected to be called directly by the user, except for ease in debugging
#'   user-defined clustering functions.
#'
#'
#' @param x \code{p x n} data matrix on which to run the clustering (samples in
#'   columns).
#' @param diss \code{n x n} data matrix of dissimilarities between the samples
#'   on which to run the clustering
##' @param clusterFunction clusterFunction a function that clusters a nxn matrix
#'   of dissimilarities/distances. Can also be given character values to
#'   indicate use of internal wrapper functions for default methods. See Details
#'   for the format of what the function must take as arguments and what format
#'   the function must return.
#' @param typeAlg character value of either '01' or 'K' determining whether the
#'   function given in clusterFunction should be called by clusterK or
#'   cluster01. Only used if clusterFunction is a user-defined function.
#'   Otherwise, for methods provided by the package (i.e. by user setting
#'   clusterFunction to a character value) clusterD will determine the
#'   appropriate input for 'typeAlg' and will ignore user input.
#' @param distFunction a distance function to be applied to \code{D}. Only relevant if
#'   input \code{D} is a matrix of data, rather than a distance. See details.
#' @param minSize the minimum number of samples in a cluster. Clusters found
#'   below this size will be discarded and samples in the cluster will be given
#'   a cluster assignment of "-1" to indicate that they were not clustered.
#' @param orderBy how to order the cluster (either by size or by maximum alpha
#'   value).
#' @param format whether to return a list of indices in a cluster or a vector of
#'   clustering assignments. List is mainly for compatibility with sequential
#'   part.
#' @param clusterArgs arguments to be passed directly to the clusterFunction,
#'   beyond the required input.
#' @param alpha a cutoff value of how much similarity needed for drawing blocks
#'   (lower values more strict).
#' @param findBestK logical, whether should find best K based on average
#'   silhouette width (only used if clusterFunction of type "K").
#' @param k single value to be used to determine how many clusters to find, if
#'   findBestK=FALSE (only used if clusterFunction of type "K").
#' @param kRange vector of integers. If findBestK=TRUE, this gives the range of
#'   k's to look over. Default is k-2 to k+20, subject to those values being
#'   greater than 2. Note that default values depend on the input k, so running
#'   for different choices of k and findBestK=TRUE can give different answers
#'   unless kRange is set to be the same.
#' @param silCutoff Requirement on minimum silhouette width to be included in
#'   cluster (only if removeSil=TRUE).
#' @param removeSil logical as to whether remove when silhouette < silCutoff
#'   (only used if clusterFunction of type "K")
#' @param checkArgs logical as to whether should give warning if arguments given
#'   that don't match clustering choices given. Otherwise, inapplicable
#'   arguments will be ignored without warning.
#' @param returnData logical as to whether to return the D matrix in output.
#' @param ... arguments given to clusterD to be passed to cluster01 or clusterK
#'   (depending on the value of typeAlg). Examples include 'k' for clusterK or
#'   'alpha' for cluster01. These should not be the arguments needed by
#'   clusterFunction (which should be passed via the argument 'clusterArgs') but
#'   the actual arguments of cluster01 or clusterK.
#' @details To provide a distance matrix via the argument \code{distFunction},
#'     the function must be defined to take the distance of the rows of a matrix
#'     (internally, the function will call \code{distFunction(t(x))}. This is to
#'     be compatible with the input for the \code{dist} function.
#'     \code{as.matrix} will be performed on the output of \code{distFunction},
#'     so if the object returned has a \code{as.matrix} method that will convert
#'     the output into a symmetric matrix of distances, this is fine (for
#'     example the class \code{dist} for objects returned by \code{dist} have
#'     such a method). If \code{distFunction=NA}, then a default distance will 
#'     be calculated based on the type of clustering algorithm of 
#'     \code{clusterFunction}. For type "K" the default is to take \code{dist}
#'     as the distance function. For type "01", the default is to take the
#'     (1-cor(x))/2.
#'
#' @details cluster01 and
#'   clusterK are given as separate functions in order to allow the user to
#'   provide different clustering functions that expect different types of input
#'   and for us to provide different shared processing of the results that is
#'   different for these different types of clustering methods (for example,
#'   removing low silhouette values is appropriate for clusterK clustering
#'   functions rather than cluster01 functions). 
#' @return clusterD returns a vector of cluster assignments (if format="vector")
#'   or a list of indices for each cluster (if format="list"). Clusters less
#'   than minSize are removed. If orderBy="size" the clusters are reordered by
#'   the size of the cluster, instead of by the internal ordering of the
#'   clusterFunction.
#'
#' @return cluster01 and clusterK return a list of indices of the clusters found,
#'   which each element of the list corresponding to a cluster and the elements
#'   of that list a vector of indices giving the indices of the samples assigned
#'   to that cluster. Indices not included in any list are assumed to have not
#'   been clustered. The list is assumed to be ordered in terms of the `best'
#'   cluster (as defined by the clusterFunction for cluster01 or by average
#'   silhoute for clusterK), for example in terms of most internal similarity of
#'   the elements, or average silhouette width.
#'
#' @examples
#' data(simData)
#' cl1<-clusterD(simData,clusterFunction="pam",k=3)
#' cl2<-clusterD(simData,clusterFunction="hierarchical01")
#' cl3<-clusterD(simData,clusterFunction="tight")
#' #change distance to manhattan distance
#' cl4<-clusterD(simData,clusterFunction="pam",k=3,
#'      distFunction=function(x){dist(x,method="manhattan")})
#' 
#' #run hierarchical method for finding blocks, with method of evaluating
#' #coherence of block set to evalClusterMethod="average", and the hierarchical
#' #clustering using single linkage:
#' clustSubHier <- clusterD(simData, clusterFunction="hierarchical01", alpha=0.1,
#' minSize=5, clusterArgs=list(evalClusterMethod="average", method="single"))
#'
#' #do tight
#' clustSubTight <- clusterD(simData, clusterFunction="tight", alpha=0.1,
#' minSize=5)
#'
#' #two twists to pam
#' clustSubPamK <- clusterD(simData, clusterFunction="pam", silCutoff=0, minSize=5,
#' removeSil=TRUE, k=3)
#' clustSubPamBestK <- clusterD(simData, clusterFunction="pam", silCutoff=0,
#' minSize=5, removeSil=TRUE, findBestK=TRUE, kRange=2:10)
#'
#' # note that passing the wrong arguments for an algorithm results in warnings
#' # (which can be turned off with checkArgs=FALSE)
#' clustSubTight_test <- clusterD(simData, clusterFunction="tight", alpha=0.1,
#' minSize=5, removeSil=TRUE)
#' clustSubTight_test2 <- clusterD(simData, clusterFunction="tight", alpha=0.1,
#' clusterArgs=list(evalClusterMethod="average"))
#' @importFrom cluster daisy silhouette pam
#' @export
setMethod(
  f = "clusterD",
  signature = signature(clusterFunction = "character"),
  definition = function(clusterFunction,...){
  	clusterD(getBuiltInClusterFunction(clusterFunction),...)
	  
  }
 )
#' @rdname clusterD
#' @export
setMethod(
   f = "clusterD",
   signature = signature(clusterFunction = "ClusterFunction"),
definition=function(clusterFunction,x=NULL, diss=NULL,
                   distFunction=NA,clusterArgs=NULL,minSize=1, orderBy=c("size","best"),
                   format=c("vector","list"),checkArgs=TRUE,checkDiss=TRUE,returnData=FALSE,...){
	orderBy<-match.arg(orderBy)
	format<-match.arg(format)
	postProcessArgs<-list(...)
	if(length(postProcessArgs)>0){
	#get rid of wrong args passed because of user confusion between the two
		whRightArgs<-which(names(postProcessArgs) %in% getPostProcessingArgs(clusterFunction))
		if(length(whRightArgs)!=length(postProcessArgs) & checkArgs) warning("Some arguments passed via '...' in clusterD do not match the algorithmType of the given ClusterFunction object")
		if(length(whRightArgs)>0) postProcessArgs<-postProcessArgs[whRightArgs]
		else postProcessArgs<-NULL
	}
	#######################
	### Check input and Create distance if needed, and check it.
	#######################
	input<-.checkXDissInput(x,diss,inputType=clusterFunction@inputType,algType=clusterFunction@algorithmType,checkDiss=checkDiss)
	if(input=="X" & clusterFunction@inputType=="diss"){
		diss<-.makeDiss(x,distFunction=distFunction,checkDiss=checkDiss,algType=clusterFunction@algorithmType)
		input<-"diss"
	}
	#-----
	# Other Checks
	#-----
	 reqArgs<-requiredArgs(clusterFunction)
	 #remove required args not needed if certain postProcessArgs are given:
	 if(algorithmType(clusterFunction)=="K" & "findBestK" %in% names(postProcessArgs)){
		 if(postProcessArgs[["findBestK"]]) reqArgs<-reqArgs[-which(reqArgs=="k")]
	 }
	 if(length(reqArgs)>0 & !all(reqArgs %in% names(clusterArgs))) stop(paste("For this clusterFunction algorithm type ('",algorithmType(clusterFunction),"') must supply arguments",reqArgs,"as elements of the list of 'clusterArgs'"))
	 if(input %in% c("X","both")) N <- dim(x)[2] else N<-dim(diss)[2]
	
	#######################
	####Run clustering:
	#######################
	
	if(algorithmType(clusterFunction)=="01") {
	 	argsClusterList<-switch(input,"X"=list(x=x), "diss"=list(diss=diss))
	 	argsClusterList<-c(argsClusterList,list("checkArgs"=checkArgs,"cluster.only"=TRUE))
	    result<-do.call(clusterFunction@clusterFUN,c(argsClusterList,clusterArgs))
	}
	if(algorithmType(clusterFunction)=="K") {
	 	argsClusterList<-switch(input,"X"=list(x=x), "diss"=list(diss=diss))
	 	argsClusterList<-c(argsClusterList,list("checkArgs"=checkArgs,"cluster.only"=TRUE))
		res<-do.call(".postProcessClusterK",c(list(clusterFunction=clusterFunction,clusterArgs=argsClusterList,N=N,orderBy=orderBy),passedArgs))
		###Note to self: .postProcessClusterK returns clusters in list form.
	}

	#######################
	#Now format into desired output, order
	#######################
	#this is perhaps not efficient. For now will do this, then consider going back and only converting when, where needed.
	if(clusterFunction@outputType=="vector" & algorithmType(clusterFunction)!="K"){
		res<-.clusterVectorToList(res)
	}
	clusterSize<-sapply(res, length)
    if(length(res)>0) res <- res[clusterSize>=minSize]
	if(length(res)!=0 & orderBy=="size"){ #i.e. there exist clusters found that passed minSize
		  clusterSize<-sapply(res, length) #redo because dropped small clusters earlier
		  res <- res[order(clusterSize,decreasing=TRUE)]
	}
	if(format=="vector"){
			res<-.clusterListToVector(res,N)
			names(res)<-if(input=="X") colnames(X) else rownames(diss)
	}
	if(!returnData) return(res)
	else return(list(result=res,diss=diss,x=x))
}



#' @rdname clusterD
#' @aliases getPostProcessingArgs
#' @export
setMethod(
  f = "getPostProcessingArgs",
  signature = c("character"),
  definition = function(clusterFunction) {
  	switch(algorithmType(clusterFunction),"01"=.argsPostCluster01,"K"=.argsPostClusterK)
)

.argsPostCluster01<-c("")
.argsPostClusterK<-c("findBestK","kRange","removeSil","silCutoff")

.postProcessClusterK<-function(clusterFunction,findBestK=FALSE,  kRange,removeSil=FALSE,silCutoff=0,clusterArgs,N,orderBy)
{
  k<-clusterArgs[["k"]]
  D<-diss
  if(!findBestK && is.null(k)) stop("If findBestK=FALSE, must provide k")
  if(findBestK){
    if(missing(kRange)){
      if(!is.null(k)) kRange<-(k-2):(k+20)
      else kRange<-2:20
    }
    if(any(kRange<2)){
      kRange<-kRange[kRange>=2]
      if(length(kRange)==0) stop("Undefined values for kRange; must be greater than or equal to 2")
    }
  }
  if(findBestK) ks<-kRange else ks<-k
  if(any(ks>= nrow(D))) ks<-ks[ks<nrow(D)]
  	clusters<-lapply(ks,FUN=function(currk){
		cl<-do.call(clusterFunction@clusterFUN,c(list(k=currk),clusterArgs))
		if(clusterFunction@outputType=="list") cl<-.clusterListToVector(cl,N=N)
		return(cl)
	})
  silClusters<-lapply(clusters,function(cl){
    silhouette(cl,dmatrix=D)
  })
  if(length(ks)>1){
    whichBest<-which.max(sapply(silClusters, mean))
    finalCluster<-clusters[[whichBest]]
    sil<-silClusters[[whichBest]][,"sil_width"]
  }
  else{
    finalCluster<-clusters[[1]]
    sil<-silClusters[[1]][,"sil_width"]
  }
  if(removeSil){
    cl<-as.numeric(sil>silCutoff)
    cl[cl==0]<- -1
    cl[cl>0]<-finalCluster[cl>0]
    sil[cl == -1] <- -Inf #make the -1 cluster the last one in order
  }
  else{
    cl<-finalCluster
  }
  
  #make list of indices and put in order of silhouette width (of positive)
  clList<-tapply(1:length(cl),cl,function(x){x},simplify=FALSE)
  if(orderBy=="best"){
	  clAveWidth<-tapply(sil,cl,mean,na.rm=TRUE)
	  clList[order(clAveWidth,decreasing=TRUE)]
  }
  #remove -1 group
  if(removeSil){
    whNotAssign<-which(sapply(clList,function(x){all(cl[x]== -1)}))
    if(length(whNotAssign)>1) stop("Coding error in removing unclustered samples")
    if(length(whNotAssign)>0) clList<-clList[-whNotAssign]
  }
  return(clList)
  
}


