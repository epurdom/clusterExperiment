#' @name addClusterings
#' @title Add clusterings to RSECClass object
#' @description Function for adding new clusterings in form of vector (single
#'   cluster) or matrix (multiple clusterings) to an existing RSECClass
#'   object
#' @param x a RSECClass object
#' @param y additional clusters to add to x. Can be a ClusterExperiment object,
#'   RSECClass object, or a matrix/vector of clusters.
#' @param clusterLabels label(s) for the clusters being added. If \code{y} a
#'   matrix, the column names of that matrix will be used by default, if
#'   \code{clusterLabels} is not given.
#' @param clusterLegend a list giving the cluster legend for the clusters added.
#' @aliases addClusterings removeClusterings
#'   addClusterings,ClusterExperiment,matrix-method
#' @inheritParams ClusterExperiment-class
#' @return A \code{ClusterExperiment} object.
#' @details addClusterings adds y to x, and is thus not symmetric in the two
#'   arguments. In particular, the \code{primaryCluster}, all of the dendrogram
#'   information, the merge information, \code{coClustering}, and 
#'   \code{orderSamples} are all kept from
#'   the x object, even if y is a ClusterExperiment.
#' @export
#' @examples
#' data(simData)
#'
#' cl1 <- clusterSingle(simData, subsample=FALSE,
#' sequential=FALSE, mainClusterArgs=list(clusterArgs=list(k=3), 
#' clusterFunction="pam"))

#' cl2 <- clusterSingle(simData, subsample=FALSE,
#' sequential=FALSE, mainClusterArgs=list(clusterArgs=list(k=3), 
#' clusterFunction="pam"))
#'
#' addClusterings(cl1, cl2)
setMethod(
  f = "addClusterings",
  signature = signature("RSECClass", "ANY"),
  definition = function(x, y, transferFrom=c("x","y"),mergeObjects,...) {
    retval <- callNextMethod()
    if(transferFrom=="y"){
        retval@merge_index<-y@merge_index+nClusterings(x) #update index to where merge from
        retval@merge_dendrocluster_index <- 
            y@merge_dendrocluster_index + nClusterings(x)
        if(.typeOfCoClustering(y)=="indices")
            retval@coClustering<-y@coClustering+nClusterings(x)
    }
		if(mergeObjects){
        ### If missing in transfer object, pull from the other one
        if(transferFrom=="x"){
            addedObj<-y
            defaultObj<-x
        }
        else{
            addedObj<-x
            defaultObj<-y
        }
        if(is.na(retval@merge_index) & !is.na(addedObj@merge_index)){
            if(transferFrom=="x") 
                retval@merge_index<-addedObj@merge_index+nClusterings(defaultObj) 
            retval@merge_nodeMerge<-addedObj@merge_nodeMerge
            retval@merge_cutoff<-addedObj@merge_cutoff
            retval@merge_method<-addedObj@merge_method
            retval@merge_demethod<-addedObj@merge_demethod
        }
        if(is.null(retval@merge_nodeProp) & !is.null(addedObj@merge_nodeProp)){
            retval@merge_nodeProp<-addedObj@merge_nodeProp
            retval@merge_dendrocluster_index<-
                addedObj@merge_dendrocluster_index+nClusterings(defaultObj)                }
        if(is.null(retval@coClustering)){
            if(transferFrom=="x" & .typeOfCoClustering(addedObj)=="indices")    
                retval@coClustering<-addedObj@coClustering + nClusterings(defaultObj)
            else retval@coClustering<-addedObj@coClustering
        }
    }

  }
)

