#' @title Cluster distance matrix from subsampling
#'
#' @description Given input data, this function will
#'   try to find the clusters based on the given ClusterFunction object. 

#' @name clusterD
#' @aliases clusterD-character-method
#'
#' @param x \code{p x n} data matrix on which to run the clustering (samples in
#'   columns).
#' @param diss \code{n x n} data matrix of dissimilarities between the samples
#'   on which to run the clustering
#' @param distFunction a distance function to be applied to \code{D}. Only relevant if
#'   input \code{D} is a matrix of data, rather than a distance. See details.
#' @param minSize the minimum number of samples in a cluster. Clusters found
#'   below this size will be discarded and samples in the cluster will be given
#'   a cluster assignment of "-1" to indicate that they were not clustered.
#' @param orderBy how to order the cluster (either by size or by maximum alpha
#'   value). If orderBy="size" the numbering of the clusters are reordered by
#'   the size of the cluster, instead of by the internal ordering of the
#'   \code{clusterFUN} defined in the \code{ClusterFunction} object (an internal ordering is only possible if slot \code{outputType} of the \code{ClusterFunction} is \code{"list"}).
#' @param format whether to return a list of indices in a cluster or a vector of
#'   clustering assignments. List is mainly for compatibility with sequential
#'   part.
#' @param clusterArgs arguments to be passed directly to the \code{clusterFUN} slot of the \code{ClusterFunction} object
#' @param checkArgs logical as to whether should give warning if arguments given
#'   that don't match clustering choices given. Otherwise, inapplicable
#'   arguments will be ignored without warning.
#' @param returnData logical as to whether to return the \code{diss} or \code{x} matrix in the output. If \code{FALSE} only the clustering vector is returned.
#' @param ... arguments passed to the post-processing steps of the clustering. The available post-processing arguments for a \code{ClusterFunction} object depend on it's algorithm type and can be found by calling \code{getPostProcessingArgs}. See details below for documentation.
#' @inherits subsampleClustering
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
#' @return clusterD returns a vector of cluster assignments (if format="vector")
#'   or a list of indices for each cluster (if format="list"). Clusters less
#'   than minSize are removed. 
#'
#' @examples
#' data(simData)
#' cl1<-clusterD(x=simData,clusterFunction="pam",clusterArgs=list(k=3))
#' cl2<-clusterD(simData,clusterFunction="hierarchical01",clusterArgs=list(alpha=.1))
#' cl3<-clusterD(simData,clusterFunction="tight",clusterArgs=list(alpha=.1))
#' #change distance to manhattan distance
#' cl4<-clusterD(simData,clusterFunction="pam",clusterArgs=list(k=3),
#'      distFunction=function(x){dist(x,method="manhattan")})
#' 
#' #run hierarchical method for finding blocks, with method of evaluating
#' #coherence of block set to evalClusterMethod="average", and the hierarchical
#' #clustering using single linkage:
#' clustSubHier <- clusterD(simData, clusterFunction="hierarchical01",
#' minSize=5, clusterArgs=list(alpha=0.1,evalClusterMethod="average", method="single"))
#'
#' #do tight
#' clustSubTight <- clusterD(simData, clusterFunction="tight", clusterArgs=list(alpha=0.1),
#' minSize=5)
#'
#' #two twists to pam
#' clustSubPamK <- clusterD(simData, clusterFunction="pam", silCutoff=0, minSize=5,
#' removeSil=TRUE, clusterArgs=list(k=3))
#' clustSubPamBestK <- clusterD(simData, clusterFunction="pam", silCutoff=0,
#' minSize=5, removeSil=TRUE, findBestK=TRUE, kRange=2:10)
#'
#' # note that passing the wrong arguments for an algorithm results in warnings
#' # (which can be turned off with checkArgs=FALSE)
#' clustSubTight_test <- clusterD(simData, clusterFunction="tight",
#' clusterArgs=list(alpha=0.1), minSize=5, removeSil=TRUE)
#' clustSubTight_test2 <- clusterD(simData, clusterFunction="tight",
#' clusterArgs=list(alpha=0.1,evalClusterMethod="average"))
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
	 # don't need define 'k' if choose 'findBestK=TRUE'
	 if(algorithmType(clusterFunction)=="K" & "findBestK" %in% names(postProcessArgs)){
		 if(postProcessArgs[["findBestK"]]) reqArgs<-reqArgs[-which(reqArgs=="k")]
	 }
	 if(length(reqArgs)>0 & !all(reqArgs %in% names(clusterArgs))) stop(paste("For this clusterFunction algorithm type ('",algorithmType(clusterFunction),"') must supply arguments",reqArgs,"as elements of the list of 'clusterArgs'"))
	 if(input %in% c("X","both")) N <- dim(x)[2] else N<-dim(diss)[2]
	
	#######################
	####Run clustering:
	#######################
	
	argsClusterList<-.makeDataArgs(dataInput=input,funInput=clusterFunction@inputType, xData=x, dissData=diss)
	argsClusterList<-c(argsClusterList, clusterArgs, list("checkArgs"=checkArgs, "cluster.only"=TRUE))
	if(algorithmType(clusterFunction)=="01") {
	    res<-do.call(clusterFunction@clusterFUN,argsClusterList)
	}
	if(algorithmType(clusterFunction)=="K") {
		res<-do.call(".postProcessClusterK",c(list(clusterFunction=clusterFunction,clusterArgs=argsClusterList,N=N,orderBy=orderBy,diss=diss),postProcessArgs))
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
			names(res)<-if(input=="X") colnames(x) else rownames(diss)
	}
	if(!returnData) return(res)
	else return(list(result=res,diss=diss,x=x))
}
)



#' @rdname clusterD
#' @aliases getPostProcessingArgs
#' @details Post-processing Arguments: For post-processing the clustering, currently only type 'K' algorithms have a defined post-processing. Specifically
#' \itemize{
#'  \item{"findBestK"}{logical, whether should find best K based on average
#'   silhouette width (only used if clusterFunction of type "K").}
#'  \item{"kRange"}{vector of integers to try for k values if findBestK=TRUE. If \code{k} is given in \code{clusterArgs}, then default is k-2 to k+20, subject to those values being
#'   greater than 2; if not the default is \code{2:20}. Note that default values depend on the input k, so running
#'   for different choices of k and findBestK=TRUE can give different answers
#'   unless kRange is set to be the same.}
#'  \item{"removeSil"}{logical as to whether remove the assignment of a sample to a cluster when the sample's silhouette value is less than \code{silCutoff}}
#'  \item{"silCutoff"}{Cutoff on the minimum silhouette width to be included in
#'   cluster (only used if removeSil=TRUE).}
#' }
#' @export
setMethod(
  f = "getPostProcessingArgs",
  signature = c("ClusterFunction"),
  definition = function(clusterFunction) {
  	switch(algorithmType(clusterFunction),"01"=.argsPostCluster01,"K"=.argsPostClusterK)
}
)

.argsPostCluster01<-c("")
.argsPostClusterK<-c("findBestK","kRange","removeSil","silCutoff")

#' @importFrom cluster silhouette
.postProcessClusterK<-function(clusterFunction,findBestK=FALSE,  kRange,removeSil=FALSE,silCutoff=0,clusterArgs,N,orderBy,diss=NULL)
{
  doPostProcess<-findBestK | removeSil #whether will calculate silhouette or not; if not, speeds up the function... 
  k<-clusterArgs[["k"]]
  if(!findBestK && is.null(k)) stop("If findBestK=FALSE, must provide k")
  if(!is.null(k)) clusterArgs<-clusterArgs[-which(names(clusterArgs)=="k")]
  if(findBestK){
    if(missing(kRange)){
      if(!is.null(k)) kRange<-(k-2):(k+20)
      else kRange<-2:20
    }
    if(any(kRange<2)){
      kRange<-kRange[kRange>=2]
      if(length(kRange)==0) stop("Undefined values for kRange; must be greater than or equal to 2")
    }
	ks<-kRange 
  }
  else ks<-k
  if(any(ks>= N)) ks<-ks[ks<N]
  clusters<-lapply(ks,FUN=function(currk){
		cl<-do.call(clusterFunction@clusterFUN,c(list(k=currk),clusterArgs))
		if(clusterFunction@outputType=="list") cl<-.clusterListToVector(cl,N=N)
		return(cl)
	})
	if(doPostProcess){
	    silClusters<-lapply(clusters,function(cl){
	      cluster::silhouette(cl,dmatrix=diss)
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
	}
	else{
		cl<-clusters[[1]]
	}
  
  
  #make list of indices and put in order of silhouette width (of positive)
  clList<-tapply(1:length(cl),cl,function(x){x},simplify=FALSE)
  if(doPostProcess){
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
  }

  return(clList)
  
}


