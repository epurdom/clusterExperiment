#' @name addClusterings
#' @title Add clusterings to ClusterExperiment object
#' @description Function for adding new clusterings in form of vector (single
#'   cluster) or matrix (multiple clusterings) to an existing ClusterExperiment
#'   object
#' @param x a ClusterExperiment object
#' @param y additional clusters to add to x. Can be a ClusterExperiment object
#'   or a matrix/vector of clusters.
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
  signature = signature("ClusterExperiment", "matrix"),
  definition = function(x, y, clusterTypes="User",clusterLabels=NULL,clusterLegend=NULL) {
    if(!is.null(clusterLabels)){
      if(length(clusterLabels)!=ncol(y)) stop("clusterLabels must vector of length equal to the number of clusterings (columns of y)")
      colnames(y)<-clusterLabels
    }
    ccObj<-ClusterExperiment(x,
                             clusters=y,
                             transformation=transformation(x),
                             clusterTypes=clusterTypes,
                             checkTransformAndAssay=FALSE,
                             clusterLegend=clusterLegend)
    addClusterings(x,ccObj)
  }
)

#' @rdname addClusterings
#' @param transferFrom If x and y are both \code{ClusterExperiment} objects
#'   indicates from which object the clustering info should be taken (regarding merging, dendrogram, etc). Does not affect the order of the clusterings, which will always be the clusterings of x, followed by those of y (along with slots `clusterType`, `clusterInfo`, `clusterLegend`)
#' @param mergeCEObjects logical If x and y are both \code{ClusterExperiment} objects indicates as to whether should try to grab in the information missing from x from y (or vice versa if transferFrom=y).

#' @export
setMethod(
    f = "addClusterings",
    signature = signature("ClusterExperiment", "ClusterExperiment"),
    definition = function(x, y, transferFrom=c("x","y"),mergeCEObjects=FALSE) 
{
    transferFrom<-match.arg(transferFrom)
    if(!all(dim(assay(y)) == dim(assay(x))) || !all(assay(y) == assay(x))) {
      stop("Cannot merge clusters from different data.")
    }
    ## FIXME: does make copy -- could be much more memory than before. 
    if(transferFrom=="x") retval<-x
    else retval<-y
    retval@clusterMatrix <- cbind(x@clusterMatrix, y@clusterMatrix)
    retval@clusterTypes <- c(x@clusterTypes, y@clusterTypes)
    retval@clusterInfo<-c(x@clusterInfo,y@clusterInfo)
    retval@clusterLegend<-c(x@clusterLegend,y@clusterLegend)
    if(any(duplicated(colnames(retval@clusterMatrix)))){
      colnames(retval@clusterMatrix)<-
          make.names(colnames(retval@clusterMatrix),unique=TRUE)
    }
    if(transferFrom=="y"){
        retval@dendro_index<-y@dendro_index+nClusterings(x)
        retval@merge_index<-y@merge_index+nClusterings(x) #update index to where merge from
        retval@merge_dendrocluster_index <- 
            y@merge_dendrocluster_index + nClusterings(x)
        if(.typeOfCoClustering(y)=="indices")
            retval@coClustering<-y@coClustering+nClusterings(x)
    }

    if(mergeCEObjects){
        ### If missing in transfer object, pull from the other one
        if(transferFrom=="x"){
            addedObj<-y
            defaultObj<-x
        }
        else{
            addedObj<-x
            defaultObj<-y
        }
        if(is.na(retval@dendro_index) & !is.na(addedObj@dendro_index)){
            retval@dendro_samples<-addedObj@dendro_samples
            retval@dendro_clusters<-addedObj@dendro_clusters
            if(transferFrom=="x") 
                retval@dendro_index <- 
                    addedObj@dendro_index+ nClusterings(defaultObj) 
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
        #put back orderSamples
        if(all(retval@orderSamples==seq_len(nSamples(retval))) & 
            !all(addedObj@orderSamples==seq_len(nSamples(retval)))) 
                retval@orderSamples<-addedObj@orderSamples
        if(is.null(retval@coClustering)){
            if(transferFrom=="x" & .typeOfCoClustering(addedObj)=="indices")    
                retval@coClustering<-addedObj@coClustering + nClusterings(defaultObj)
            else retval@coClustering<-addedObj@coClustering
        }
    }
    retval<-.unnameClusterSlots(retval) #just gets rid of the names of objects that shouldn't have them

    ch<-.checkClusterMatrix(retval)
    if(!is.logical(ch)) stop(ch)
    ch<-.checkClusterTypes(retval)
    if(!is.logical(ch)) stop(ch)
    ch<-.checkClusterLegend(retval)
    if(!is.logical(ch)) stop(ch)
    #would it be less memory to do a call to "new"? What is difference versus having to check dendrogram, coClustering, etc if they exists? Should do checks on large data.
    return(retval)
  }
)

#' @rdname addClusterings
#' @export
#' @param ... For \code{addClusterings}, passed to signature
#'   \code{ClusterExperiment,matrix}. For \code{[} (subsetting), passed to
#'   \code{SingleCellExperiment} subsetting function.
#' @param makePrimary whether to make the added cluster the primary cluster
#'   (only relevant if \code{y} is a vector)
setMethod(
  f = "addClusterings",
  signature = signature("ClusterExperiment", "vector"),
  definition = function(x, y, makePrimary=FALSE,...) {
    mat<-matrix(y,ncol=1)
    x<-addClusterings(x,mat,...)
    if(makePrimary){
      x@primaryIndex<-ncol(clusterMatrix(x))
    }
    return(x)
  }
)
