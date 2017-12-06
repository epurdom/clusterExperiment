#' Functions to add/remove clusters to ClusterExperiment
#'
#' These functions are used to add or remove clusters to a
#' \code{\link{ClusterExperiment}} object.
#'
#' @param x a ClusterExperiment object.
#' @param y additional clusters to add to x. Can be a ClusterExperiment object
#'   or a matrix/vector of clusters.
#' @param clusterLabel label(s) for the clusters being added.
#' @inheritParams ClusterExperiment
#' @details addClusterings adds y to x, and is thus not symmetric in the two 
#'   arguments. In particular, the \code{primaryCluster}, all of the dendrogram
#'   information, \code{coClustering}, and \code{orderSamples} are all kept from
#'   the x object, even if y is a ClusterExperiment.
#'
#' @return A \code{\link{ClusterExperiment}} object with the added clusters.
#'
#' @rdname addClusterings
#' @aliases addClusterings removeClusterings
#' @export
#' @examples
#' data(simData)
#'
#' cl1 <- clusterSingle(simData, subsample=FALSE,
#' sequential=FALSE, mainClusterArgs=list(clusterArgs=list(k=3), clusterFunction="pam"))

#' cl2 <- clusterSingle(simData, subsample=FALSE,
#' sequential=FALSE, mainClusterArgs=list(clusterArgs=list(k=3), clusterFunction="pam"))
#'
#' addClusterings(cl1, cl2)
setMethod(
  f = "addClusterings",
  signature = signature("ClusterExperiment", "matrix"),
  definition = function(x, y, clusterTypes="User") {
    ccObj<-ClusterExperiment(assay(x),y,transformation=transformation(x),clusterTypes=clusterTypes,checkTransformAndAssay=FALSE)
    addClusterings(x,ccObj)
  }
)

#' @rdname addClusterings
#' @export
setMethod(
  f = "addClusterings",
  signature = signature("ClusterExperiment", "ClusterExperiment"),
  definition = function(x, y) {
    if(!all(dim(assay(y)) == dim(assay(x))) || !all(assay(y) == assay(x))) {
      stop("Cannot merge clusters from different data.")
    }
    x@clusterMatrix <- cbind(x@clusterMatrix, y@clusterMatrix)
    x@clusterTypes <- c(x@clusterTypes, y@clusterTypes)
    x@clusterInfo<-c(x@clusterInfo,y@clusterInfo)
    x@clusterLegend<-c(x@clusterLegend,y@clusterLegend)
    if(any(duplicated(colnames(x@clusterMatrix)))){
      colnames(x@clusterMatrix)<-make.names(colnames(x@clusterMatrix),unique=TRUE)
    }
    x<-.unnameClusterSlots(x) #just gets rid of the names of objects that shouldn't have them
	ch<-.checkClusterMatrix(x)
	if(!is.logical(ch)) stop(ch)
	ch<-.checkClusterTypes(x)
	if(!is.logical(ch)) stop(ch)
	ch<-.checkClusterLegend(x)
	if(!is.logical(ch)) stop(ch)
		#would it be less memory to do a call to "new"? What is difference versus having to check dendrogram, coClustering, etc if they exists? Should do checks on large data.
    return(x)
  }
)

#' @rdname addClusterings
#' @export
#' @param ... Passed to signature \code{ClusterExperiment,matrix}.
setMethod(
  f = "addClusterings",
  signature = signature("ClusterExperiment", "numeric"),
  definition = function(x, y, clusterLabel=NULL,...) {
    mat<-matrix(y,ncol=1)
    if(!is.null(clusterLabel)) colnames(mat)<-clusterLabel
    addClusterings(x,mat,...)
  }
)

#' @rdname addClusterings
#' @export
setMethod(
  f = "removeClusterings",
  signature = signature("ClusterExperiment","character"),
  definition = function(x, whichClusters,...) {
	  whichClusters<-.TypeIntoIndices(x,whichClusters)
	  removeClusterings(x,whichClusters,...)
  }
)

#' @inheritParams ClusterExperiment-methods
#'
#'@details \code{removeClusterings} removes the clusters given by
#'  \code{whichClusters}. If the
#'  \code{primaryCluster} is one of the clusters removed, the
#'  \code{primaryClusterIndex} is set to 1 and the dendrogram and coclustering
#'  matrix are discarded and orderSamples is set to \code{1:NCOL(x)}.
#' @return \code{removeClusterings} returns a \code{ClusterExperiment} object, 
#'  unless all clusters are removed, in which case it returns a
#'  \code{\link{SingleCellExperiment}} object.
#' @rdname addClusterings
#' @export
setMethod(
  f = "removeClusterings",
  signature = signature("ClusterExperiment","numeric"),
  definition = function(x, whichClusters) {
    if(any(whichClusters>NCOL(clusterMatrix(x)))) stop("invalid indices -- must be between 1 and",NCOL(clusterMatrix(x)))
    if(length(whichClusters)==NCOL(clusterMatrix(x))){
      warning("All clusters have been removed. Will return just a Summarized Experiment Object")
      #make it Summarized Experiment
      return(as(x,"SingleCellExperiment"))
    }
    
    newClLabels<-clusterMatrix(x)[,-whichClusters,drop=FALSE]
    newClusterInfo<-clusteringInfo(x)[-whichClusters]
    newClusterType<-clusterTypes(x)[-whichClusters]
    newClusterColors<-clusterLegend(x)[-whichClusters]
    dend_samples <- x@dendro_samples
    dend_cl <- x@dendro_clusters
    dend_ind<-dendroClusterIndex(x)
    dend_out<-x@dendro_outbranch
    coMat<-x@coClustering
    orderSamples<-orderSamples(x)
    if(primaryClusterIndex(x) %in% whichClusters) pIndex<-1
    else pIndex<-match(primaryClusterIndex(x),(1:NCOL(clusterMatrix(x)))[-whichClusters])
    if(dendroClusterIndex(x) %in% whichClusters){
        dend_cl<-NULL
        dend_samples<-NULL
        dend_ind<-NA_real_
        dend_out<-NA
    }
    else{
      dend_ind<-match(dend_ind,(1:NCOL(clusterMatrix(x)))[-whichClusters])
    }
    
    retval<-ClusterExperiment(as(x,"SingleCellExperiment"),newClLabels,transformation(x),
                              clusterTypes=newClusterType,
                              clusterInfo<-newClusterInfo,
                              primaryIndex=pIndex,
                              dendro_samples=dend_samples,
                              dendro_clusters=dend_cl,
                            dendro_index=dend_ind,
                            dendro_outbranch=dend_out,
                            coClustering=coMat,
                            orderSamples=orderSamples,
							checkTransformAndAssay=FALSE
                              )
    clusterLegend(retval)<-newClusterColors
    return(retval)
  }
)

#' @details \code{removeUnclustered} removes all samples that are unclustered
#'   (i.e. -1 or -2 assignment) in the \code{primaryCluster} of x (so they may
#'   be unclustered in other clusters found in \code{clusterMatrix(x)}).
#' @rdname addClusterings
#' @aliases removeUnclustered
#' @export
setMethod(
  f = "removeUnclustered",
  signature = "ClusterExperiment",
  definition = function(x) {
    return(x[,primaryCluster(x) >= 0])
  }
)


#' @details \code{unassignSamples} unassigns samples in \code{clustersToRemove} and assigns them to -1 (unassigned) 
#' @rdname addClusterings
#' @aliases unassignSamples
#' @export
setMethod(
  f = "unassignSamples",
  signature = c("ClusterExperiment","numeric"),
  definition = function(x,whichClusters,clustersToRemove) {
	  if(length(whichClusters)!=1) stop("whichClusters should identify a single clustering.")
	  cl<-clusterMatrix(x)[,whichClusters]
	  if(is.numeric(clustersToRemove)){
		  cl[cl %in% clustersToRemove]<- -1
	  }
	  else if(is.character(clustersToRemove)){
	  	
	  }
	return(x[,primaryCluster(x) >= 0])
	
  }
)
#' @rdname addClusterings
#' @export
setMethod(
  f = "unassignSamples",
  signature = signature("ClusterExperiment","character"),
  definition = function(x, whichClusters,...) {
	  whichClusters<-.TypeIntoIndices(x,whichClusters)
	  removeClusters(x,whichClusters,...)
  }
)