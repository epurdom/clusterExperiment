#' Helper methods for the ClusterExperiment class
#'
#' This is a collection of helper methods for the ClusterExperiment class.
#' @name ClusterExperiment-methods
#' @aliases ClusterExperiment-methods [,ClusterExperiment,ANY,ANY,ANY-method [,ClusterExperiment,ANY,character,ANY-method
#' @details Note that when subsetting the data, the dendrogram information and
#' the co-clustering matrix are lost.
#' @export
#' @param ...,i,j,drop Forwarded to the
#'   \code{\link{SingleCellExperiment}} method.
#' @param value The value to be substituted in the corresponding slot. See the
#'   slot descriptions in \code{\link{ClusterExperiment}} for details on what
#'   objects may be passed to these functions.
setMethod(
  f = "[",
  signature = c("ClusterExperiment", "ANY", "character"),
  definition = function(x, i, j, ..., drop=TRUE) {
    j<-match(j, colnames(x))
    callGeneric()
    
  }
)
#' @rdname ClusterExperiment-methods
#' @export
setMethod(
  f = "[",
  signature = c("ClusterExperiment", "ANY", "logical"),
  definition = function(x, i, j, ..., drop=TRUE) {
    j<-which(j)
    callGeneric()
  }
)
#' @rdname ClusterExperiment-methods
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
    out<-.makeColors(clMat=subMat, distinctColors=FALSE,colors=massivePalette,                           matchClusterLegend=clusterLegend(x),matchTo="name") 
		newMat<-out$numClusters
    colnames(newMat)<-nms
    newClLegend<-out$colorList
    #fix order of samples so same
    newOrder<-rank(x@orderSamples[j])
    #
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

## show
#' @rdname ClusterExperiment-methods
#' @export
setMethod(
  f = "show",
  signature = "ClusterExperiment",
  definition = function(object) {
    cat("class:", class(object), "\n")
    cat("dim:", dim(object), "\n")
    cat("reducedDimNames:",if(anyValidReducedDims(object)) reducedDimNames(object) else "no reduced dims stored","\n")
    cat("filterStats:",if(anyValidFilterStats(object)) filterNames(object) else "no valid filtering stats stored","\n")
    cat("-----------\n")
    cat("Primary cluster type:", clusterTypes(object)[primaryClusterIndex(object)],"\n")
    cat("Primary cluster label:", clusterLabels(object)[primaryClusterIndex(object)],"\n")
    cat("Table of clusters (of primary clustering):")
    print(table(primaryClusterNamed(object)))
    cat("Total number of clusterings:", NCOL(clusterMatrix(object)),"\n")
    if(!is.na(dendroClusterIndex(object)) ) cat("Dendrogram run on '",clusterLabels(object)[dendroClusterIndex(object)],"' (cluster index: ", dendroClusterIndex(object),")\n",sep="") else cat("No dendrogram present\n")
    cat("-----------\n")
    cat("Workflow progress:\n")
    typeTab<-names(table(clusterTypes(object)))
    cat("clusterMany run?",if("clusterMany" %in% typeTab) "Yes" else "No","\n")
    cat("makeConsensus run?",if("makeConsensus" %in% typeTab) "Yes" else "No","\n")
    cat("makeDendrogram run?",if(!is.null(object@dendro_samples) & !is.null(object@dendro_clusters) ) "Yes" else "No","\n")
    cat("mergeClusters run?",if("mergeClusters" %in% typeTab) "Yes" else "No","\n")
  }
)




#' @rdname ClusterExperiment-methods
#' @return \code{transformation} prints the function used to transform the data
#' prior to clustering.
#' @export
#' @aliases transformation
setMethod(
  f = "transformation",
  signature = "ClusterExperiment",
  definition = function(x) {
    return(x@transformation)
  }
)

#' @rdname ClusterExperiment-methods
#' @export
#' @details Note that redefining the transformation function via
#'   \code{transformation(x)<-} will check the validity of the transformation on
#'   the data assay. If the assay is large, this may be time consuming. Consider
#'   using a call to ClusterExperiment, which has the option as to whether to
#'   check the validity of the transformation.
#' @aliases transformation<-
setReplaceMethod(
  f = "transformation",
  signature = signature("ClusterExperiment", "function"),
  definition = function(object, value) {
	checkValidity=TRUE
    object@transformation <- value
    if(checkValidity){
		ch<-.checkTransform(object)
    	if(ch) return(object) else stop(ch)
	}
	else return(object)
  }
)

#' @rdname ClusterExperiment-methods
#' @return \code{nClusterings} returns the number of clusterings (i.e., ncol of
#' clusterMatrix).
#' @export
#' @aliases nClusterings
setMethod(
  f = "nClusterings",
  signature = "ClusterExperiment",
  definition = function(x){
    return(NCOL(clusterMatrix(x)))
  }
)


#' @rdname ClusterExperiment-methods
#' @return \code{nClusters} returns the number of clusters per clustering
#' @param ignoreUnassigned logical. If true, ignore the clusters with -1 or -2 assignments in calculating the number of clusters per clustering. 
#' @export
#' @aliases nClusters
setMethod(
  f = "nClusters",
  signature = "ClusterExperiment",
  definition = function(x,ignoreUnassigned=TRUE){
	  if(ignoreUnassigned){
		  return(apply(clusterMatrix(x),2,function(x){length(unique(x[x>0]))}))
	  }
	  else return(apply(clusterMatrix(x),2,function(x){length(unique(x))}))
  }
)

#' @rdname ClusterExperiment-methods
#' @return \code{nFeatures} returns the number of features (same as `nrow`).
#' @aliases nFeatures
#' @export
setMethod(
  f = "nFeatures",
  signature =  "ClusterExperiment",
  definition = function(x){
    return(NROW(assay(x)))
  }
)

#' @rdname ClusterExperiment-methods
#' @return \code{nSamples} returns the number of samples (same as `ncol`).
#' @aliases nSamples
#' @export
setMethod(
  f = "nSamples",
  signature = "ClusterExperiment",
  definition = function(x){
    return(NCOL(assay(x)))
  }
)

#' @rdname ClusterExperiment-methods
#' @return \code{clusterMatrixNamed} returns a matrix with cluster labels.
#' @export
#' @aliases clusterMatrixNamed
#' @param x,object a ClusterExperiment object.
setMethod(
  f = "clusterMatrixNamed",
  signature = "ClusterExperiment",
  definition = function(x, whichClusters="all") {
    convertClusterLegend(x,output="matrixNames",whichClusters=whichClusters)
  }
)

#' @rdname ClusterExperiment-methods
#' @return \code{clusterMatrixColors} returns the matrix with all the clusterings, using the internally stored colors for each cluster
#' @export
#' @aliases clusterMatrixColors
setMethod(
  f = "clusterMatrixColors",
  signature = c("ClusterExperiment"),
  definition = function(x,whichClusters) {
    convertClusterLegend(x,output="matrixColors",whichClusters=whichClusters)
  }
)

#' @rdname ClusterExperiment-methods
#' @param whichClusters argument that can be either numeric or
#'   character value. If numeric, gives the indices of the \code{clusterMatrix}
#'   to return; this can also be used to defined an ordering for the
#'   clusterings. \code{whichClusters} can be a character value identifying the 
#'   \code{clusterTypes} to be used, or if not matching \code{clusterTypes} then
#'   \code{clusterLabels}; alternatively \code{whichClusters} can be either 
#'   'all' or 'workflow' to indicate choosing all clusters or choosing all 
#'   \code{\link{workflowClusters}}. If missing, the entire matrix of all
#'   clusterings is returned.
#' @return \code{clusterMatrix} returns the matrix with all the clusterings.
#' @export
#' @aliases clusterMatrix
setMethod(
  f = "clusterMatrix",
  signature = c("ClusterExperiment","missing"),
  definition = function(x,whichClusters) {
    wh<-seq_len(ncol(x@clusterMatrix))
    return(clusterMatrix(x,whichClusters=wh))
  }
)
#' @rdname ClusterExperiment-methods
#' @return \code{clusterMatrix} returns the matrix with all the clusterings.
#' @export
#' @aliases clusterMatrix
setMethod(
  f = "clusterMatrix",
  signature = c("ClusterExperiment","numeric"),
  definition = function(x,whichClusters) {
	  mat<-x@clusterMatrix[,whichClusters,drop=FALSE]
	  rownames(mat)<-colnames(x)
    return(mat)
  }
)
#' @rdname ClusterExperiment-methods
#' @return \code{clusterMatrix} returns the matrix with all the clusterings.
#' @export
#' @aliases clusterMatrix
setMethod(
  f = "clusterMatrix",
  signature = c("ClusterExperiment","character"),
  definition = function(x,whichClusters) {
	  wh<-.TypeIntoIndices(x,whClusters=whichClusters)
	  return(clusterMatrix(x,whichClusters=wh))
  }
)


#' @rdname ClusterExperiment-methods
#' @return \code{primaryCluster} returns the primary clustering (as numeric).
#' @export
#' @aliases primaryCluster
setMethod(
  f = "primaryCluster",
  signature = "ClusterExperiment",
  definition = function(x) {
    return(x@clusterMatrix[,primaryClusterIndex(x)])
  }
)

#' @rdname ClusterExperiment-methods
#' @return \code{primaryClusterIndex} returns/sets the primary clustering index
#' (i.e., which column of clusterMatrix corresponds to the primary clustering).
#' @export
#' @aliases primaryClusterIndex
setMethod(
  f = "primaryClusterIndex",
  signature = "ClusterExperiment",
  definition = function(x) {
    return(x@primaryIndex)
  }
)

#' @rdname ClusterExperiment-methods
#' @return \code{primaryClusterIndex} returns/sets the primary clustering index
#' (i.e., which column of clusterMatrix corresponds to the primary clustering).
#' @export
#' @aliases primaryClusterLabel
setMethod(
  f = "primaryClusterLabel",
  signature = "ClusterExperiment",
  definition = function(x) {
    return(clusterLabels(x)[primaryClusterIndex(x)])
  }
)
#' @rdname ClusterExperiment-methods
#' @return \code{primaryClusterNamed} returns the primary cluster (using cluster
#' labels).
#' @export
#' @aliases primaryClusterNamed
setMethod(
  f = "primaryClusterNamed",
  signature = "ClusterExperiment",
  definition = function(x) {
    as.vector(clusterMatrixNamed(x,whichCluster="primary"))
  })

#' @rdname ClusterExperiment-methods
#' @return \code{primaryClusterIndex} returns/sets the primary clustering index
#' (i.e., which column of clusterMatrix corresponds to the primary clustering).
#' @export
#' @aliases primaryClusterType
setMethod(
  f = "primaryClusterType",
  signature = "ClusterExperiment",
  definition = function(x) {
    return(clusterTypes(x)[primaryClusterIndex(x)])
  }
)


#' @rdname ClusterExperiment-methods
#' @return \code{subsetByCluster} subsets the object by clusters in a clustering
#' and returns a ClusterExperiment object with only those samples
#' @param clusterValue values of the cluster to match to for subsetting
#' @param matchTo for subsetting, whether to match to the cluster name
#'   (\code{"name"}) or internal cluster id (\code{"clusterIds"})
#' @export
#' @aliases subsetByCluster
setMethod(
  f = "subsetByCluster",
  signature = "ClusterExperiment",
  definition = function(x,clusterValue,whichCluster="primary",matchTo=c("name","clusterIds")) {
    
		whCl<-.convertSingleWhichCluster(x,whichCluster)
		matchTo<-match.arg(matchTo)
		if(matchTo=="name"){
			cl<-clusterMatrixNamed(x)[,whCl]
		}
		else cl<-clusterMatrix(x)[,whCl]
		return(x[,which(cl %in% clusterValue)])
  }
)


#' @rdname ClusterExperiment-methods
#' @export
#' @aliases primaryClusterIndex<-
setReplaceMethod(
  f = "primaryClusterIndex",
  signature = signature("ClusterExperiment", "numeric"),
  definition = function(object, value) {
    object@primaryIndex <- value
    ch<-.checkPrimaryIndex(object)
    if(is.logical(ch) && ch) return(object) else stop(ch)
  }
)
#' @rdname ClusterExperiment-methods
#' @return \code{dendroClusterIndex} returns/sets the clustering index 
#' of the clusters used to create dendrogram
#' (i.e., which column of clusterMatrix corresponds to the clustering).
#' @export
#' @aliases dendroClusterIndex
setMethod(
  f = "dendroClusterIndex",
  signature = "ClusterExperiment",
  definition = function(x) {
    return(x@dendro_index)
  }
)



#' @rdname ClusterExperiment-methods
#' @return \code{coClustering} returns/sets the co-clustering matrix.
#' @export
#' @aliases coClustering
setMethod(
  f = "coClustering",
  signature = "ClusterExperiment",
  definition = function(x) {
    return(x@coClustering)
  }
)

#' @rdname ClusterExperiment-methods
#' @export
#' @aliases coClustering<-
setReplaceMethod(
  f = "coClustering",
  signature = signature(object="ClusterExperiment", value="matrix"),
  definition = function(object, value) {
    object@coClustering <- value
    ch<-.checkCoClustering(object)
    if(is.logical(ch) && ch) return(object) else stop(ch)
  }
)

#' @rdname ClusterExperiment-methods
#' @return \code{clusterTypes} returns/sets the clusterTypes slot.
#' @export
#' @aliases clusterTypes
setMethod(
  f = "clusterTypes",
  signature = "ClusterExperiment",
  definition = function(x) {
    out<-x@clusterTypes
    #names(out)<-clusterLabels(x)
    return(out)
  }
)

#' @rdname ClusterExperiment-methods
#' @return \code{clusteringInfo} returns the clusterInfo slot.
#' @aliases clusteringInfo
#' @export
setMethod(
  f = "clusteringInfo",
  signature = "ClusterExperiment",
  definition = function(x) {
    out<-x@clusterInfo
    names(out)<-clusterLabels(x)
    return(out)
  }
)


#' @rdname ClusterExperiment-methods
#' @return \code{clusterLabels} returns/sets the column names of the clusterMatrix slot.
#' @export
#' @aliases clusterLabels
setMethod(
  f = "clusterLabels",
  signature = signature(x = "ClusterExperiment"),
  definition = function(x){
    labels<-colnames(clusterMatrix(x))
    if(is.null(labels)) cat("No labels found for clusterings\n")
    return(labels)

  }
)
#' @export
#' @rdname ClusterExperiment-methods
#' @aliases clusterLabels<-
setReplaceMethod( 
  f = "clusterLabels",
  signature = signature(object="ClusterExperiment", value="character"),
  definition = function(object, value) {
    if(length(value)!=NCOL(clusterMatrix(object))) stop("value must be a vector of length equal to NCOL(clusterMatrix(object)):",NCOL(clusterMatrix(object)))
    colnames(object@clusterMatrix) <- value
    ch<-.checkClusterLabels(object)
    if(is.logical(ch) && ch) return(object) else stop(ch)
  }
)
#' @rdname ClusterExperiment-methods
#' @return \code{clusterLegend} returns/sets the clusterLegend slot.
#' @export
#' @aliases clusterLegend
setMethod(
    f = "clusterLegend",
    signature = "ClusterExperiment",
    definition = function(x) {
      out<-x@clusterLegend
      names(out)<-clusterLabels(x)
      return(out)
    }
)

#' @rdname ClusterExperiment-methods
#' @export
#' @aliases clusterLegend<-
setReplaceMethod( 
  f = "clusterLegend",
  signature = signature(object="ClusterExperiment", value="list"),
  definition = function(object, value) {
    object@clusterLegend<-unname(value)
    ch<-.checkClusterLegend(object)
    if(is.logical(ch) && ch) return(object) else stop(ch)
  }
)

#' @rdname ClusterExperiment-methods
#' @return \code{renameClusters} changes the names assigned to clusters within a clustering
#' @param whichCluster argument to identify cluster, taking input like 
#' \code{whichClusters}, only requires that only a single clustering can be identified. 
#' @export
#' @aliases renameClusters
setMethod( 
  f = "renameClusters",
  signature = signature(object="ClusterExperiment", value="character"),
  definition = function(object, value,whichCluster="primary") {
		whCl<-.convertSingleWhichCluster(object,whichCluster)
		mat<-clusterLegend(object)[[whCl]]
		clVals<-as.numeric(mat[,"clusterIds"])
		if(is.null(names(value))){
			
			if(length(value)== nrow(mat)) names(value)<-mat[,"clusterIds"]
			else if(length(value)==length(clVals[clVals>0])) names(value)<-mat[clVals>0,"clusterIds"]
			else stop("length of argument 'value' not equal to number of clusters, nor does it have names to identify it to 'clusterIds' of this clustering.")
		} 
		if(!all(names(value) %in% mat[,"clusterIds"])) stop("'value' must be vector with names that match the 'clusterIds' column of the requested clusterLegend")
			
			m<-match(names(value),mat[,"clusterIds"])
		mat[m,"name"]<-value
		clusterLegend(object)[[whCl]]<-mat
		
    ch<-.checkClusterLegend(object)
    if(is.logical(ch) && ch) return(object) else stop(ch)
  }
)

#' @rdname ClusterExperiment-methods
#' @return \code{recolorClusters} changes the colors assigned to clusters within a clustering
#' @export
#' @aliases recolorClusters
setMethod( 
  f = "recolorClusters",
  signature = signature(object="ClusterExperiment", value="character"),
  definition = function(object, value,whichCluster="primary") {
		whCl<-.convertSingleWhichCluster(object,whichCluster)
		mat<-clusterLegend(object)[[whCl]]
		
		if(is.null(names(value)) || !all(names(value) %in% mat[,"clusterIds"])) stop("'value' must be vector with names matching the 'clusterIds' column of the requested clusterLegend")
			
			m<-match(names(value),mat[,"clusterIds"])
		mat[m,"color"]<-value
		clusterLegend(object)[[whCl]]<-mat
		
    ch<-.checkClusterLegend(object)
    if(is.logical(ch) && ch) return(object) else stop(ch)
  }
)

#' @rdname ClusterExperiment-methods
#' @return \code{orderSamples} returns/sets the orderSamples slot.
#' @export
#' @aliases orderSamples
setMethod(
  f = "orderSamples",
  signature = "ClusterExperiment",
  definition = function(x) {
    return(x@orderSamples)
  }
)
#' @rdname ClusterExperiment-methods
#' @export
#' @aliases orderSamples<-
setReplaceMethod( 
  f = "orderSamples",
  signature = signature(object="ClusterExperiment", value="numeric"),
  definition = function(object, value) {
    object@orderSamples<-value
    ch<-.checkOrderSamples(object) 
    if(is.logical(ch) && ch) return(object) else stop(ch)
    
  }
)

#' @rdname ClusterExperiment-methods
#' @export
#' @aliases clusterTypes<-
setReplaceMethod( 
  f = "clusterTypes",
  signature = signature(object="ClusterExperiment", value="character"),
  definition = function(object,value) {
    object@clusterTypes<-value
    object<-.unnameClusterSlots(object)
    ch<-.checkClusterTypes(object)
    if(is.logical(ch) && ch) return(object) else stop(ch)
    
  }
)

