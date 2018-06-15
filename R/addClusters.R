#' Functions to add/remove clusters to ClusterExperiment
#'
#' These functions are used to add or remove clusters to a
#' \code{\link{ClusterExperiment}} object.
#' @name addClusterings
#' @param x a ClusterExperiment object.
#' @param y additional clusters to add to x. Can be a ClusterExperiment object
#'   or a matrix/vector of clusters.
#' @param clusterLabels label(s) for the clusters being added. If \code{y} a
#'   matrix, the column names of that matrix will be used by default, if
#'   \code{clusterLabels} is not given.
#' @param clusterLegend a list giving the cluster legend for the clusters added.
#' @inheritParams ClusterExperiment
#' @details addClusterings adds y to x, and is thus not symmetric in the two
#'   arguments. In particular, the \code{primaryCluster}, all of the dendrogram
#'   information, \code{coClustering}, and \code{orderSamples} are all kept from
#'   the x object, even if y is a ClusterExperiment.
#'
#' @return A \code{\link{ClusterExperiment}} object with the added clusters.
#'
#' @rdname addClusterings
#' @aliases addClusterings removeClusterings addClusterings,ClusterExperiment,matrix-method
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
#' @param makePrimary whether to make the added cluster the primary cluster (only relevant if \code{y} is a vector)
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
    coMat<-x@coClustering
    orderSamples<-orderSamples(x)
    if(primaryClusterIndex(x) %in% whichClusters) pIndex<-1
    else 
			pIndex<-match(primaryClusterIndex(x),seq_len(NCOL(clusterMatrix(x)))[-whichClusters])
		
		#fix dendro info
	  dend_samples <- x@dendro_samples
    dend_cl <- x@dendro_clusters
    dend_ind<-dendroClusterIndex(x)
    dend_out<-x@dendro_outbranch
    if(dendroClusterIndex(x) %in% whichClusters){
      dend_cl<-NULL
      dend_samples<-NULL
      dend_ind<-NA_real_
      dend_out<-NA
    }
    else{
      dend_ind<-match(dend_ind,seq_len(NCOL(clusterMatrix(x)))[-whichClusters])
    }
		#fix merge info:
		#erase merge info if either dendro or merge index deleted.
	  if(mergeClusterIndex(x) %in% whichClusters | x@merge_dendrocluster_index %in% whichClusters){
      merge_index=NA_real_
      merge_cutoff=NA_real_
      merge_dendrocluster_index=NA_real_
      merge_nodeProp=NULL
      merge_nodeMerge=NULL
      merge_method=NA_character_
	  }
		else{
      merge_index<-match(x@merge_index,seq_len(NCOL(clusterMatrix(x)))[-whichClusters])
      merge_dendrocluster_index<-match(x@merge_dendrocluster_index, seq_len(NCOL(clusterMatrix(x)))[-whichClusters])
      merge_cutoff=x@merge_cutoff
      merge_nodeProp=x@merge_nodeProp
      merge_nodeMerge=x@merge_nodeMerge
      merge_method=x@merge_method
		}

    retval<-ClusterExperiment(as(x,"SingleCellExperiment"),
                              clusters=newClLabels,
                              transformation=transformation(x),
                              clusterTypes=newClusterType,
                              clusterInfo<-newClusterInfo,
                              primaryIndex=pIndex,
                              dendro_samples=dend_samples,
                              dendro_clusters=dend_cl,
                              dendro_index=dend_ind,
                              dendro_outbranch=dend_out,
															merge_index=merge_index,
															merge_dendrocluster_index=merge_dendrocluster_index,
												      merge_cutoff=merge_cutoff,
												      merge_nodeProp=merge_nodeProp,
												      merge_nodeMerge=merge_nodeMerge,
												      merge_method=merge_method,
                              coClustering=coMat,
                              orderSamples=orderSamples,
                              clusterLegend=newClusterColors,
                              checkTransformAndAssay=FALSE
    )
    #    clusterLegend(retval)<-newClusterColors
    return(retval)
  }
)



#' @details \code{removeClusters} creates a new cluster that unassigns samples in cluster \code{clustersToRemove} (in the clustering defined by \code{whichClusters}) and assigns them to -1 (unassigned)
#' @param clustersToRemove numeric vector identifying the clusters to remove (whose samples will be reassigned to -1 value).
#' @param whichCluster Clustering from which to remove clusters for
#'  \code{removeCluster}. Note that it is a singular cluster.
#' @rdname addClusterings
#' @aliases removeClusters
#' @export
setMethod(
  f = "removeClusters",
  signature = c("ClusterExperiment","numeric"),
  definition = function(x,whichCluster,clustersToRemove,makePrimary=FALSE,clusterLabels=NULL) {
    whCl<-.convertSingleWhichCluster(x,whichCluster)
    cl<-clusterMatrix(x)[,whCl]
    leg<-clusterLegend(x)[[whCl]]
    if(is.character(clustersToRemove)){
      m<- match(clustersToRemove,leg[,"name"] )
      if(any(is.na(m)))
        stop("invalid names of clusters in 'clustersToRemove'")
      clustersToRemove<-as.numeric(leg[m,"clusterIds"])
    }
    if(is.numeric(clustersToRemove)){
      if(any(!clustersToRemove %in% cl)) stop("invalid clusterIds in 'clustersToRemove'")
      if(any(clustersToRemove== -1)) stop("cannot remove -1 clusters using this function. See 'assignUnassigned' to assign unassigned samples.")
      cl[cl %in% clustersToRemove]<- -1
    }
    else stop("clustersToRemove must be either character or numeric")
    if(is.null(clusterLabels)){
      currlabel<-clusterLabels(x)[whCl]
      clusterLabels<-paste0(currlabel,"_unassignClusters")
    }
    if(clusterLabels %in% clusterLabels(x))
      stop("must give a 'clusterLabels' value that is not already assigned to a clustering")
    newleg<-leg
    if(!"-1" %in% leg[,"clusterIds"] & any(cl== -1)){
      newleg<-rbind(newleg,c("-1","white","-1"))
    }
    whRm<-which(as.numeric(newleg[,"clusterIds"]) %in% clustersToRemove )
    if(length(whRm)>0){
      newleg<-newleg[-whRm,,drop=FALSE]
    }
		newCl<-list(newleg)
		#names(newCl)<-clusterLabels
    return(addClusterings(x, cl,  clusterLabels = clusterLabels,clusterLegend=newCl,makePrimary=makePrimary))


  }
)
#' @rdname addClusterings
#' @export
setMethod(
  f = "removeClusters",
  signature = signature("ClusterExperiment","character"),
  definition = function(x, whichCluster,...) {
    whichCluster<-.TypeIntoIndices(x,whichCluster)
    removeClusters(x,whichCluster,...)
  }
)
