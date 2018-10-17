#' @name getClusterIndex
#' @title getClusterIndex
#' @description Finds index of clustering in clusterMatrix slot of object based on descriptions of clusters.
#' @param object a ClusterExperiment object
#' @param silentlyRemove logical as to whether to silently remove mismatches. Otherwise values that do not match are given NA values, unless all values are NA in which an error is returned.
#' @param whichClusters argument that can be either numeric or character vector
#'   indicating the clusterings to be used. See details of \code{\link{getClusterIndex}}.
#' @param whichCluster argument that can be a single numeric or character value
#'   indicating the \emph{single} clustering to be used. Giving values that result in more than one clustering will result in an error. See details of \code{\link{getClusterIndex}}.
#' @param passedArgs other arguments passed to the function (only used internally)
#' @details The function \code{getClusterIndex} is largely used internally to parse the argument \code{whichClusters} which is used as an argument extensively across functions in this package. Note that some functions require the match return a single clustering, in which case those functions use the function \code{getSingleClusterIndex} with the singular argument \code{whichCluster} and returns an error if it indicates more than one clustering. Furthermore \code{getSingleClusterIndex} does not allow for any mismatches (\code{noMatch="throwError"}. Otherwise the parsing of the two arguments \code{whichClusters} and \code{whichCluster} is the same, and is described in what follows.
#' @details If \code{whichClusters} is numeric, then the function just returns the 
#'  numeric values of \code{whichClusters}, after checking that they are valid. If any are invalid, they are silently removed if \code{silentlyRemove=TRUE}. The values will be returned  \emph{in the order given}, so this argument can also be used to defined by functions to give an
#'   ordering for the clusterings (as relevant).
#' @details If \code{whichClusters} is a character value, then it the function attempts to use the character value to identify the clustering. The value of the argument is first matched against a set of "special" values: "workflow","all","none","primaryCluster","dendro" using the argument \code{\link{match.arg}}, which does partial matching. If whichClusters is a vector of values, only the first value of the vector is matched against these values and if it matches, the remaining values are ignored.  If it matches one of these values, then the cluster indices are given as follows.
##' \itemize{
##'  \item{"workflow"}{all clusterings in the current workflow (see \code{\link{workflowClusters}})}
##'  \item{"all"}{all clusterings, with the primary clustering put first.}
##'  \item{"none"}{no clusterings}
##'  \item{"primaryCluster"}{the primary clustering index as given by \code{\link{primaryClusterIndex}}}
##'  \item{"dendro"}{the index of the clustering given to create the cluster dendrogram, if it exists}
##' }
#'  @details If \code{whichClusters} is a character value, but its first element does not match these predesignated values, then all the values of \code{whichClusters} are attempted to be matched to the \code{\link{clusterTypes}} of the object. Note that there may be more than one clustering that matches a given type. For any entries that do not match a value in  \code{clusterTypes(object)} are then matched based on the value of \code{\link{clusterLabels}} of the object. 
#' @export
#' @return \code{getClusterIndex} returns a vector of all numeric indices that are indicated by the requested \code{whichClusters}. Note that there is not a one-to-one match between input values and returned values since there may be more than one value for a given value of \code{whichClusters} or no value at all.  
setMethod(
	f="getClusterIndex",
	signature="ClusterExperiment",
	definition=function(object,whichClusters,noMatch=c("silentlyRemove","throwError")){
	noMatch<-match.arg(noMatch)
	.createMismatch<-function(errorMessage){
	  if(noMatch=="silentlyRemove")
		  return(vector("integer",length=0))
	  else if(noMatch=="throwError")
		  stop(errorMessage)
	  else if(noMatch=="NA")
		  return(NA)  
	}
	if(is.numeric(whichClusters)) wh<-whichClusters
	else{
	  test<-try( match.arg(whichClusters[1], c("workflow","all","none","primaryCluster","dendro")), silent=TRUE)
	  if(!inherits(test,"try-error")){
	    if(test=="workflow"){
	      ppIndex<-workflowClusterDetails(object)
	      if(!is.null(ppIndex) && sum(ppIndex[,"iteration"]==0)>0){
	        wh<-unlist(lapply(.workflowValues,function(tt){
	          ppIndex[ppIndex[,"iteration"]==0 & ppIndex[,"type"]==tt,"index"]
	        }))
	      }
	      else{
			  wh<-.createMismatch("There are no workflow clusterings")
		  }
	    }
	    if(test=="all"){
	      #put primary cluster first
	      ppcl<-primaryClusterIndex(object)
	      wh<-c(ppcl,c(seq_len(nClusterings(object)))[-ppcl])
	    }
	    if(test=="none") wh<-.createMismatch("No clusterings requested")
	    if(test=="primaryCluster") wh<-primaryClusterIndex(object)
	    if(test=="dendro"){
	      wh<-dendroClusterIndex(object)
	      if(is.na(wh)) wh<-.createMismatch("There is no dendrogram, cannot return clustering that created one")
	    }
	  }
	  else{
	    #first match to clusterTypes  
	    mClType<-match(whichClusters,clusterTypes(object))  
	    mClLabel<-match(whichClusters,clusterLabels(object))  
	    totalMatch<-mapply(whichClusters,mClType,mClLabel,FUN=function(cl,type,lab){
	      if(is.na(type) & !is.na(lab)) return(lab)
	      if(is.na(type) & is.na(lab)) return(NA)
	      if(!is.na(type)){
	        return(which(clusterTypes(object) %in% cl)) #prioritize clusterType and get ALL of them, not just first match
	      }
	    },SIMPLIFY=FALSE)
	    wh<-unlist(totalMatch,use.names=FALSE)
	  }
	} 
	if(any(!wh %in% seq_len(nClusterings(object)))) {
		if(noMatch=="throwError") 
			stop("Invalid value for 'whichCluster'. Must be integer between 1 and ", nClusterings(object))
		if(noMatch=="silentlyRemove") wh[wh<=nClusterings(object) & wh>0]<-NA
	}
	if(all(is.na(wh))){
		if(noMatch=="throwError") 
			stop("whichCluster(s) did not match ANY clustering in the object")
		if(noMatch=="silentlyRemove") wh<-vector("integer",length=0)	
	} 
	else{
	  if(any(is.na(wh))){
	  	  if(noMatch=="throwError") 
	  		stop("Not all values in whichCluster(s) matched a clustering in the object")
		  if(noMatch=="silentlyRemove") wh<-na.omit(wh) #silently ignore things that don't match.
	  }
	}
	return(wh)



})
	
#' @rdname getClusterIndex
#' @export
setMethod(
	f="getSingleClusterIndex",
	signature="ClusterExperiment",
	definition=function(object,whichCluster,passedArgs=NULL){
	if(!is.null(passedArgs) && any(c("whichClusters") %in% names(passedArgs))){
		stop("The argument of this function is 'whichCluster' (singular) not 'whichClusters' indicating only a single clustering can be used for this cluster")
	}
  whCl<-getClusterIndex(object,whichClusters=whichCluster,throwError=TRUE)
  if(length(whCl)!=1) stop("Invalid value for 'whichCluster'. Current value of the argument identifies ",length(whCl)," clusterings, but this function requires that it must identify only a single clustering (a singular 'whichCluster' argument).")

  return(whCl)
})

