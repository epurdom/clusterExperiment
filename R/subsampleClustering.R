#' Cluster subsamples of the data
#'
#' Given a data matrix, this function will subsample the rows
#' (samples), cluster the subsamples, and return a \code{n x n} matrix with the
#' probability of co-occurance.
#'
#' @param x the data on which to run the clustering (samples in columns).
#' @param k number of clusters to find for each clustering of a subsample
#'   (passed to clusterFunction).
#' @param clusterFunction a \code{\link{ClusterFunction}} object that clusters defines clustering routine. Can also be given character values to indicate use
#'   of internal wrapper functions. Must accept arguments 'x' and 'k' (whether
#'   uses them or not). See Details for format of what must return.
#' @param clusterArgs a list of parameter arguments to be passed to
#'   clusterFunction.
#' @param resamp.num the number of subsamples to draw.
#' @param samp.p the proportion of samples to sample for each subsample.
#' @param classifyMethod method for determining which samples should be used in
#'   the co-occurance matrix. "All"= all samples, "OutOfSample"= those not
#'   subsampled, and "InSample"=those in the subsample.  "All" and "OutOfSample"
#'   require that you provide classifyFunction to define how to classify those
#'   samples not in the subsample into a cluster. If "All" is chosen, all
#'   samples will be classified into clusters via the classifyFunctions, not
#'   just those that are out-of-sample. Note if not choose 'All' possible to get
#'   NAs in resulting D matrix (particularly if not enough subsamples taken).
#' @param classifyFunction a function which, given the output of clusterFunction
#'   and new data points, will classify the new data points into a cluster.
#' @param ncores integer giving the number of cores. If ncores>1, mclapply will
#'   be called.
#' @param ... arguments passed to mclapply (if ncores>1).
#'
#' @details The \code{clusterFunction} must be a function that takes as an
#'   argument 'x' which is a \code{p x n} matrix  of data and integer 'k'. It
#'   minimally must return a list with element named 'clustering' giving the
#'   vector of cluster ids. To be incorporated with the larger hierarchy, it
#'   should be list with elements of a partition object, just as is returned by
#'   \code{\link[cluster]{pam}}. Generally, the user will need to write a
#'   wrapper function to do this. In the case of pam or kmeans, the user can
#'   identify clusterFunction as "pam" or "kmeans", and the package functions
#'   will use internally written wrappers for the clusterFunction and
#'   classifyFunction arguments. Additional arguments should be supplied via
#'   clusterArgs.
#'
#' @details The classifyFunction should take as an object a data matrix 'x' with
#'   samples on the columns, and the output of the clusterFunction. Note that the
#'   function should assume that the input 'x' is not the same samples that were
#'   input to the clusterFunction (but can assume that it is the same number of
#'   features/columns).
#'
#' @return A \code{n x n} matrix of co-occurances.
#'
#' @examples
#' data(simData)
#'
#' subD <- subsampleClustering(t(simData), clusterArgs=list(k=3), clusterFunction="kmeans",
#' clusterArgs=list(nstart=10), resamp.n=100, samp.p=0.7)
#'
#' heatmap(subD)
#' @export
setMethod(
  f = "subsampleClustering",
  signature = signature(clusterFunction = "character"),
  definition = function(clusterFunction,...){
  	subsampleClustering(getBuiltInClusterFunction(clusterFunction),...)
	  
  }
 )
 
#' @rdname subsampleClustering
#' @export
setMethod(
   f = "subsampleClustering",
   signature = signature(clusterFunction = "ClusterFunction"),
definition=function(clusterFunction, x,diss,clusterArgs=NULL, 
                              classifyMethod=c("InSample","OutOfSample","All"),
                              resamp.num = 100, samp.p = 0.7,ncores=1,checkDiss=TRUE,... )
{
	#-----
	# Checks
	#-----
  classifyMethod<-match.arg(classifyMethod)
  input<-.checkXDissInput(x,diss,inputType=clusterFunction@inputType,checkDiss=checkDiss)
  if(classifyMethod %in% c("All","OutOfSample") && is.null(clusterFunction@classifyFUN)){
    classifyMethod<-"InSample" #silently change it...
  }
  else{
	  inputClassify<-.checkXDissInput(x, diss, inputType=clusterFunction@inputClassifyType, checkDiss=checkDiss)  	
  }
  reqArgs<-requiredArgs(clusterFunction)
  if(!all(reqArgs %in% clusterArgs)) stop("For this clusterFunction algorithm type (",algorithmType(clusterFunction),") must supply arguments",reqArgs,"as list in 'clusterArgs'")
  
#-----
# Basic parameters, subsamples
#-----
  if(input %in% c("X","both")) N <- dim(x)[2] else N<-dim(diss)[2]
  subSize <- round(samp.p * N)
  idx<-replicate(resamp.num,sample(1:N,size=subSize)) #each column a set of indices for the subsample.
  #-----
  # Function that calls the clustering for each subsample
  #-----
  perSample<-function(ids){
	  ##----
	  ##Cluster part of subsample
	  ##----
	 argsClusterList<-switch(input,"X"=list(x=x[,ids,drop=FALSE]), "diss"=list(diss=diss[ids,ids,drop=FALSE]))
	 argsClusterList<-c(argsClusterList,list("checkArgs"=TRUE,"cluster.only"=FALSE))
    result<-do.call(clusterFunction@clusterFUN,c(argsClusterList,clusterArgs))

	  ##----
	  ##Classify part of subsample
	  ##----
    if(classifyMethod=="All"){
	    argsClassifyList<-switch(inputClassify,"X"=list(x=x), "diss"=list(diss=diss))
		classX<-do.call(clusterFunction@classifyFUN,c(argsClassifyList,list(result=result)))
	}
    if(classifyMethod=="OutOfSample"){
	    argsClassifyList<-switch(inputClassify,"X"=list(x=x[,-ids,drop=FALSE]), "diss"=list(diss=diss[-ids,-ids,drop=FALSE]))
		classElse<-do.call(clusterFunction@classifyFUN,c(argsClassifyList, list(result=result)))
      classX<-rep(NA,N)
      classX[-ids]<-classElse
    }
    if(classifyMethod=="InSample"){
      classX<-rep(NA,N)
      classX[ids]<-result$clustering
    }
    D <- outer(classX, classX, function(a, b) a == b)
    Dinclude<-matrix(1,N,N)
    whNA<-which(is.na(classX))
    if(length(whNA)>0){
      Dinclude[whNA,]<-0 #don't add them to the denominator either
      Dinclude[,whNA]<-0
      D[whNA,]<-0 #don't add to sum
      D[,whNA]<-0
    }
    return(list(D=D,Dinclude=Dinclude))
  }
  if(ncores==1){
    DList<-apply(idx,2,perSample)
  }
  else{
    DList<-parallel::mclapply(1:ncol(idx), function(nc){ perSample(idx[,nc]) }, mc.cores=ncores,...)
  }
  DDenom<-Reduce("+",lapply(DList,function(y){y$Dinclude}))
  DNum<-Reduce("+",lapply(DList,function(y){y$D}))
  Dbar = DNum/DDenom
  if(input %in% c("X")) rownames(Dbar)<-colnames(Dbar)<-colnames(x)
  else rownames(Dbar)<-colnames(Dbar)<-colnames(diss)
  rownames(Dbar)<-colnames(Dbar)<-colnames(x)
  return(Dbar)
}
