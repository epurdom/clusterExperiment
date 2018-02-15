#' Functions to add/remove clusters to ClusterExperiment
#'
#' These functions are used to add or remove clusters to a
#' \code{\link{ClusterExperiment}} object.
#'
#' @param x a ClusterExperiment object.
#' @param y additional clusters to add to x. Can be a ClusterExperiment object
#'   or a matrix/vector of clusters.
#' @param clusterLabels label(s) for the clusters being added. If \code{y} a matrix, the column names of that matrix will be used by default, if \code{clusterLabels} is not given. 
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
  definition = function(x, y, clusterTypes="User",clusterLabels=NULL,clusterLegend=NULL) {
   if(!is.null(clusterLabels)){
	   if(length(clusterLabels)!=ncol(y)) stop("clusterLabels must vector of length equal to the number of clusterings (columns of y)")
	   colnames(y)<-clusterLabels
   }
    ccObj<-ClusterExperiment(assay(x),
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
		coClustering=coMat,
		orderSamples=orderSamples,
		clusterLegend=newClusterColors,
		checkTransformAndAssay=FALSE
     )
#    clusterLegend(retval)<-newClusterColors
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


#' @details \code{removeClusters} creates a new cluster that unassigns samples in cluster \code{clustersToRemove} (in the clustering defined by \code{whichClusters}) and assigns them to -1 (unassigned) 
#' @param clustersToRemove numeric vector identifying the clusters to remove (whose samples will be reassigned to -1 value). 
#' @rdname addClusterings
#' @aliases removeClusters
#' @export
setMethod(
  f = "removeClusters",
  signature = c("ClusterExperiment","numeric"),
  definition = function(x,whichClusters,clustersToRemove,clusterLabels) {
	  if(length(whichClusters)!=1) stop("whichClusters should identify a single clustering.")
		 makePrimary<-whichClusters==x@primaryIndex
	  cl<-clusterMatrix(x)[,whichClusters]
	  leg<-clusterLegend(x)[[whichClusters]]
	  if(is.character(clustersToRemove)){
	  		 m<- match(clustersToRemove,leg[,"name"] )
	  		 if(any(is.na(m))) 
	  			 stop("invalid names of clusters in 'clustersToRemove'")
	  		 clustersToRemove<-as.numeric(leg[m,"clusterIds"])
	  	  }
	  if(is.numeric(clustersToRemove)){
		  if(any(!clustersToRemove %in% cl)) stop("invalid clusterIds in 'clustersToRemove'")
		  if(any(clustersToRemove== -1)) stop("cannot remove -1 clusters using this function")
		  cl[cl %in% clustersToRemove]<- -1
	  }
	  else stop("clustersToRemove must be either character or numeric")
	  if(missing(clusterLabels)){
		  currlabel<-clusterLabels(x)[whichClusters]
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
	  return(addClusterings(x, cl,  clusterLabels = clusterLabels,clusterLegend=list(newleg),makePrimary=makePrimary))
	  
	 
  }
)
#' @rdname addClusterings
#' @export
setMethod(
  f = "removeClusters",
  signature = signature("ClusterExperiment","character"),
  definition = function(x, whichClusters,...) {
	  whichClusters<-.TypeIntoIndices(x,whichClusters)
	  removeClusters(x,whichClusters,...)
  }
)