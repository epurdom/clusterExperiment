#' Functions to add/remove clusters to ClusterExperiment
#'
#' These functions are used to add or remove clusters to a
#' \code{\link{ClusterExperiment}} object.
#'
#' @param x a ClusterExperiment object.
#' @param y additional clusters to add to x. Can be a ClusterExperiment object
#'   or a matrix/vector of clusters.
#' @inheritParams clusterExperiment
#' @details addClusters adds y to x, and is thus not symmetric in the two
#'   arguments. In particular, the \code{primaryCluster} and all of its
#'   supporting information (dendrogram, coClustering, and orderSamples) are all
#'   kept from the x object, even if y is a ClusterExperiment.
#' @rdname addClusters
#' @aliases addClusters removeClusters
#' @export
#' @examples
#' data(simData)
#'
#' cl1 <- clusterSingle(simData, clusterFunction="pam", subsample=FALSE,
#' sequential=FALSE, clusterDArgs=list(k=3))
#'
#' cl2 <- clusterSingle(simData, clusterFunction="pam", subsample=FALSE,
#' sequential=FALSE, clusterDArgs=list(k=5))
#'
#' addClusters(cl1, cl2)
setMethod(
  f = "addClusters",
  signature = signature("ClusterExperiment", "matrix"),
  definition = function(x, y, clusterType="User") {
    ccObj<-clusterExperiment(assay(x),y,transformation=transformation(x),clusterType=clusterType)
    addClusters(x,ccObj)
  }
)

#' @rdname addClusters
#' @export
setMethod(
  f = "addClusters",
  signature = signature("ClusterExperiment", "ClusterExperiment"),
  definition = function(x, y) {
    if(!all(assay(y) == assay(x))) {
      stop("Cannot merge clusters from different data.")
    }
    x@clusterMatrix <- cbind(x@clusterMatrix, y@clusterMatrix)
    #browser()
    x@clusterType <- c(x@clusterType, y@clusterType)
    x@clusterInfo<-c(x@clusterInfo,y@clusterInfo)
    x@clusterLegend<-c(x@clusterLegend,y@clusterLegend)
    if(any(duplicated(colnames(x@clusterMatrix)))){
      colnames(x@clusterMatrix)<-make.names(colnames(x@clusterMatrix),unique=TRUE)
    }
    x<-.unnameClusterSlots(x)
    validObject(x)
    return(x)
  }
)

#' @rdname addClusters
#' @export
#' @param ... Passed to signature \code{ClusterExperiment,matrix}.
setMethod(
  f = "addClusters",
  signature = signature("ClusterExperiment", "numeric"),
  definition = function(x, y, ...) {
    addClusters(x,matrix(y,ncol=1),...)
  }
)

#' @rdname addClusters
#' @export
setMethod(
  f = "removeClusters",
  signature = signature("ClusterExperiment","character"),
  definition = function(x, whichRemove,exactMatch=TRUE) {
    if(exactMatch) wh<-which(clusterType(x) %in% whichRemove)
    else{
      sapply(whichRemove,grep, clusterType(x))
    }
    removeClusters(x,wh)
  }
)

#' @param exactMatch logical. Whether \code{whichRemove} must exactly match a
#'   value of \code{clusterType(x)}. Only relevant if whichRemove is character.
#' @param whichRemove which clusters to remove. Can be numeric or character. If
#'   numeric, must give indices of \code{clusterMatrix(x)} to remove. If
#'   character, should match a \code{clusterType} of x.
#'
#'@details \code{removeClusters} removes the clusters given by
#'  \code{whichRemove}. If all clusters are implied, then returns a
#'  \code{\link[SummarizedExperiment]{SummarizedExperiment}} object. If the
#'  \code{primaryCluster} is one of the clusters removed, the
#'  \code{primaryClusterIndex} is set to 1 and the dendrogram and cooccurance
#'  matrix are discarded and orderSamples is set to \code{1:NCOL(x)}.
#' @rdname addClusters
#' @export
setMethod(
  f = "removeClusters",
  signature = signature("ClusterExperiment","numeric"),
  definition = function(x, whichRemove) {
    #browser()
    if(any(whichRemove>NCOL(clusterMatrix(x)))) stop("invalid indices -- must be between 1 and",NCOL(clusterMatrix(x)))
    if(length(whichRemove)==NCOL(clusterMatrix(x))){
      warning("All clusters have been removed. Will return just a Summarized Experiment Object")
      #make it Summarized Experiment
    }
    newClLabels<-clusterMatrix(x)[,-whichRemove,drop=FALSE]
    newClusterInfo<-clusterInfo(x)[-whichRemove]
    newClusterType<-clusterType(x)[-whichRemove]
    newClusterColors<-clusterLegend(x)[-whichRemove]
    if(primaryClusterIndex(x) %in% whichRemove){
      pIndex<-1
      dend_samples<-NULL
      dend_cl <- NULL
      coMat<-new("matrix")
      orderSamples<-1:NCOL(x)
    }
    else{
      pIndex<-match(primaryClusterIndex(x),1:NCOL(clusterMatrix(x))[-whichRemove])
      dend_samples <- x@dendro_samples
      dend_cl <- x@dendro_clusters
      coMat<-x@coClustering
      orderSamples<-orderSamples(x)
    }
    retval<-clusterExperiment(as(x,"SummarizedExperiment"),newClLabels,transformation(x),clusterType=newClusterType,clusterInfo<-newClusterInfo)
    retval@coClustering<-coMat
    validObject(retval)
    clusterLegend(retval)<-newClusterColors
    primaryClusterIndex(retval)<-pIndex #Note can only set it on valid object so put it here...
    retval@dendro_samples <- dend_samples
    retval@dendro_clusters <- dend_cl
    orderSamples(retval)<-orderSamples
    return(retval)
  }
)

#' @details \code{removeUnclustered} removes all samples that are unclustered
#'   (i.e. -1 or -2 assignment) in the \code{primaryCluster} of x (so they may
#'   be unclustered in other clusters found in \code{clusterMatrix(x)}).
#' @rdname addClusters
#' @export
setMethod(
  f = "removeUnclustered",
  signature = "ClusterExperiment",
  definition = function(x) {
    return(x[,primaryCluster(x) >= 0])
  }
)
