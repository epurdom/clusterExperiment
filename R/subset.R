#' @title Functions to subset ClusterExperiment Objects
#' @description These functions are used to subset ClusterExperiment objects, either by removing samples, genes, or clusterings
#' @name subset
#' @param x a ClusterExperiment object.
#' @inheritParams ClusterExperiment-class 
#' @inheritParams ClusterExperiment-methods
#' @return A \code{\link{ClusterExperiment}} object.
#' @aliases removeClusterings,ClusterExperiment,character-method
#' @details \code{removeClusterings} removes the clusters given by
#'  \code{whichClusters}. If the
#'  \code{primaryCluster} is one of the clusters removed, the
#'  \code{primaryClusterIndex} is set to 1 and the dendrogram and coclustering
#'  matrix are discarded and orderSamples is set to \code{1:NCOL(x)}.
#' @return \code{removeClusterings} returns a \code{ClusterExperiment} object,
#'  unless all clusters are removed, in which case it returns a
#'  \code{\link{SingleCellExperiment}} object.
#' @export
setMethod(
  f = "removeClusterings",
  signature = signature("ClusterExperiment","character"),
  definition = function(x, whichClusters,...) {
    whichClusters<-.TypeIntoIndices(x,whichClusters)
    removeClusterings(x,whichClusters,...)
  }
)

#' @rdname subset
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
    if(primaryClusterIndex(x) %in% whichClusters) 
		pIndex<-1
    else 
		pIndex<-match(primaryClusterIndex(x),seq_len(NCOL(clusterMatrix(x)))[-whichClusters])
		
	
	#fix merge info:
	#erase merge info if either dendro or merge index deleted.
	if(mergeClusterIndex(x) %in% whichClusters | x@merge_dendrocluster_index %in% whichClusters){
      x<-.eraseMerge(x)
	  merge_index<-x@merge_index
	  merge_dendrocluster_index<-x@merge_dendrocluster_index
	}
	else{
      merge_index<-match(x@merge_index,seq_len(NCOL(clusterMatrix(x)))[-whichClusters])
      merge_dendrocluster_index<-match(x@merge_dendrocluster_index, seq_len(NCOL(clusterMatrix(x)))[-whichClusters])
      
	}
	#fix dendro info
	dend_samples <- x@dendro_samples
    dend_cl <- x@dendro_clusters
    dend_ind<-dendroClusterIndex(x)
    if(dendroClusterIndex(x) %in% whichClusters){
      dend_cl<-NULL
      dend_samples<-NULL
      dend_ind<-NA_real_
    }
    else{
      dend_ind<-match(dend_ind,seq_len(NCOL(clusterMatrix(x)))[-whichClusters])
    }
    retval<-ClusterExperiment(
		as(x,"SingleCellExperiment"),
		clusters=newClLabels,
		transformation=transformation(x),
		clusterTypes=newClusterType,
		clusterInfo<-newClusterInfo,
		primaryIndex=pIndex,
		dendro_samples=dend_samples,
		dendro_clusters=dend_cl,
		dendro_index=dend_ind,
		merge_index=merge_index,
		merge_dendrocluster_index=merge_dendrocluster_index,
		merge_cutoff=x@merge_cutoff,
		merge_nodeProp=x@merge_nodeProp,
		merge_nodeMerge=x@merge_nodeMerge,
		merge_method=x@merge_method,
		merge_demethod=x@merge_demethod,                              
		coClustering=coMat,
		orderSamples=orderSamples,
		clusterLegend=newClusterColors,
		checkTransformAndAssay=FALSE
    )
    return(retval)
  }
)



#' @details \code{removeClusters} creates a new cluster that unassigns samples in cluster \code{clustersToRemove} (in the clustering defined by \code{whichClusters}) and assigns them to -1 (unassigned)
#' @param clustersToRemove numeric vector identifying the clusters to remove (whose samples will be reassigned to -1 value).
#' @param whichCluster Clustering from which to remove clusters for
#'  \code{removeCluster}. Note that it is a singular cluster.
#' @rdname subset
#' @inheritParams addClusterings
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
#' @rdname subset
#' @export
setMethod(
  f = "removeClusters",
  signature = signature("ClusterExperiment","character"),
  definition = function(x, whichCluster,...) {
    whichCluster<-.TypeIntoIndices(x,whichCluster)
    removeClusters(x,whichCluster,...)
  }
)

#' @details Note that when subsetting the data, the dendrogram information and
#' the co-clustering matrix are lost.
#' @aliases [,ClusterExperiment,ANY,ANY,ANY-method [,ClusterExperiment,ANY,character,ANY-method
#' @param i,j A vector of logical or integer subscripts, indicating the rows and columns to be subsetted for \code{i} and \code{j}, respectively.
#' @param drop A logical scalar that is ignored.
#' @rdname subset
#' @export
setMethod(
  f = "[",
  signature = c("ClusterExperiment", "ANY", "character"),
  definition = function(x, i, j, ..., drop=TRUE) {
    j<-match(j, colnames(x))
    callGeneric()
    
  }
)
#' @rdname subset
#' @export
setMethod(
  f = "[",
  signature = c("ClusterExperiment", "ANY", "logical"),
  definition = function(x, i, j, ..., drop=TRUE) {
    j<-which(j)
    callGeneric()
  }
)
#' @rdname subset
#' @export
setMethod(
  f = "[",
  signature = c("ClusterExperiment", "ANY", "numeric"),
  definition = function(x, i, j, ..., drop=TRUE) {
    # #out <- callNextMethod() #doesn't work once I added the logical and character choices.
    # out<-selectMethod("[",c("SingleCellExperiment","ANY","numeric"))(x,i,j) #have to explicitly give the inherintence... not great.
    ###Note: Could fix subsetting, so that if subset on genes, but same set of samples, doesn't do any of this...
    #Following Martin Morgan advice, do "new" rather than @<- to create changed object
    #need to subset cluster matrix and convert to consecutive integer valued clusters:
		
		#pull names out so can match it to the clusterLegend. 
		subMat<-clusterMatrixNamed(x)[j, ,drop=FALSE]
		
		#danger if not unique names
		whNotUniqueNames<-vapply(clusterLegend(x),FUN=function(mat){length(unique(mat[,"name"]))!=nrow(mat)},FUN.VALUE=TRUE)
		if(any(whNotUniqueNames)){
			warning("Some clusterings do not have unique names; information in clusterLegend will not be transferred to subset.")
			subMatInt<-x@clusterMatrix[j, whNotUniqueNames,drop=FALSE]
			subMat[,whNotUniqueNames]<-subMatInt
		}
    nms<-colnames(subMat)
    ##Fix clusterLegend slot, in case now lost a level and to match new integer values
		#shouldn't need give colors, but function needs argument
    if(nrow(subMat)>0){
			out<-.makeColors(clMat=subMat, distinctColors=FALSE,colors=massivePalette,                           matchClusterLegend=clusterLegend(x),matchTo="name") 
			newMat<-out$numClusters
	    colnames(newMat)<-nms
	    newClLegend<-out$colorList
	    #fix order of samples so same
	    newOrder<-rank(x@orderSamples[j])
	    #
    	
    }
		else{
			newClLegend<-list()
			newOrder<-NA_real_
			newMat<-subMat
		}
    out<- ClusterExperiment(
      object=as(selectMethod("[",c("SingleCellExperiment","ANY","numeric"))(x,i,j),"SingleCellExperiment"),#have to explicitly give the inherintence... not great.
      clusters = newMat,
      transformation=x@transformation,
      primaryIndex = x@primaryIndex,
      clusterTypes = x@clusterTypes,
      clusterInfo=x@clusterInfo,
      orderSamples=newOrder,
      clusterLegend=newClLegend,
      checkTransformAndAssay=FALSE
    )
    #	clusterLegend(out)<-newClLegend
    return(out)
  }
)
