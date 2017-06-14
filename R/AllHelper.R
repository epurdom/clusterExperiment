#' Helper methods for the ClusterExperiment class
#'
#' This is a collection of helper methods for the ClusterExperiment class.
#' @name ClusterExperiment-methods
#' @aliases ClusterExperiment-methods [,ClusterExperiment,ANY,ANY,ANY-method [,ClusterExperiment,ANY,character,ANY-method
#' @details Note that when subsetting the data, the dendrogram information and
#' the co-clustering matrix are lost.
#' @export
#' @param ...,i,j,drop Forwarded to the
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}} method.
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
    origN<-NCOL(x)
    #out <- callNextMethod() #doesn't work once I added the logical and character choices.
    out<-selectMethod("[",c("SummarizedExperiment","ANY","numeric"))(x,i,j) #have to explicitly give the inherintence... not great.
    #browser()
    out@clusterMatrix <- as.matrix(x@clusterMatrix[j, ,drop=FALSE])
    out@coClustering <- NULL
    out@dendro_samples <- NULL
    out@dendro_clusters <- NULL
    out@dendro_index <- NA_real_
   # browser()
    #out@orderSamples<-match(out@orderSamples[j],c(1:origN)[j])
	out@orderSamples <- rank(x@orderSamples[j])
	
    #need to convert to consecutive integer valued clusters:
    newMat<-.makeIntegerClusters(out@clusterMatrix)
    colnames(newMat)<-colnames(out@clusterMatrix)
    ##Fix clusterLegend slot, in case now lost a level and to match new integer values
    newClLegend<-lapply(1:NCOL(out@clusterMatrix),function(ii){
        colMat<-out@clusterLegend[[ii]]
        newCl<-newMat[,ii]
        cl<-out@clusterMatrix[,ii]
        #remove (possible) levels lost
        whRm<-which(!colMat[,"clusterIds"] %in% as.character(cl))
        if(length(whRm)>0){
            colMat<-colMat[-whRm,,drop=FALSE]
        }
        #convert
        oldNew<-unique(cbind(old=cl,new=newCl))
        if(nrow(oldNew)!=nrow(colMat)) stop("error in converting colorLegend")
        m<-match(colMat[,"clusterIds"],oldNew[,"old"])
        colMat[,"clusterIds"]<-oldNew[m,"new"]
        return(colMat)
    })
    out@clusterMatrix<-newMat
    out@clusterLegend<-newClLegend
    #browser()
    validObject(out)
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
    cat("combineMany run?",if("combineMany" %in% typeTab) "Yes" else "No","\n")
    cat("makeDendrogram run?",if(!is.null(object@dendro_samples) & !is.null(object@dendro_clusters) ) "Yes" else "No","\n")
    cat("mergeClusters run?",if("mergeClusters" %in% typeTab) "Yes" else "No","\n")
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
  definition = function(x) {
    clMat<-clusterMatrix(x)
    out<-do.call("cbind",lapply(1:ncol(clMat),function(ii){
      cl<-clMat[,ii]
      leg<-clusterLegend(x)[[ii]]
      leg[,"name"][match(cl,leg[,"clusterIds"])]
    }))
    colnames(out)<-colnames(clMat)
    rownames(out)<-NULL
    return(out)
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
    clusterMatrixNamed(x)[,primaryClusterIndex(x)]
  })

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
#' @aliases transformation<-
setReplaceMethod(
  f = "transformation",
  signature = signature("ClusterExperiment", "function"),
  definition = function(object, value) {
    object@transformation <- value
    validObject(object)
    return(object)
  }
)

#' @rdname ClusterExperiment-methods
#' @return \code{nClusters} returns the number of clusterings (i.e., ncol of
#' clusterMatrix).
#' @export
#' @aliases nClusters
setMethod(
  f = "nClusters",
  signature = "ClusterExperiment",
  definition = function(x){
    return(NCOL(clusterMatrix(x)))
  }
)
#' @rdname ClusterExperiment-methods
#' @return \code{nFeatures} returns the number of features (same as `nrow`).
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
#' @export
setMethod(
  f = "nSamples",
  signature = "ClusterExperiment",
  definition = function(x){
    return(NCOL(assay(x)))
  }
)

#' @rdname ClusterExperiment-methods
#' @param whichClusters optional argument that can be either numeric or
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
    return(x@clusterMatrix)
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
    return(x@clusterMatrix[,whichClusters,drop=FALSE])
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
#' @export
#' @aliases primaryClusterIndex<-
setReplaceMethod(
  f = "primaryClusterIndex",
  signature = signature("ClusterExperiment", "numeric"),
  definition = function(object, value) {
    object@primaryIndex <- value
    validObject(object)
    return(object)
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
    validObject(object)
    return(object)
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
#' @return \code{clusterInfo} returns the clusterInfo slot.
#' @aliases clusterInfo
#' @export
setMethod(
  f = "clusterInfo",
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
    if(any(duplicated(value))) stop("cannot have duplicated clusterLabels")
    colnames(object@clusterMatrix) <- value
    validObject(object)
    return(object)
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
        validObject(object)
        return(object)
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
        validObject(object)
        return(object)
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
        validObject(object)
        return(object)
    }
)

# # Need to implement: wrapper to get a nice summary of the parameters choosen, similar to that of paramMatrix of clusterMany (and compatible with it)
# #' @rdname ClusterExperiment-class
# setMethod(
#   f= "paramValues",
#   signature = "ClusterExperiment",
#   definition=function(x,type){
#     whCC<-which(clusterTypes(x)==type)
#     if(length(wwCC)==0) stop("No clusterings of type equal to ",type,"are found")
#     if(type=="clusterMany"){
#       #recreate the paramMatrix return value
#       paramMatrix<-do.call("rbind",lapply(wwCC,function(ii){
#         data.frame(index=ii,clusterInfo(x)[[ii]]["choicesParam"])
#       }))
#
#     }
#     else if(type=="clusterSingle"){
#
#     }
#     else{
#       return(clusterInfo(x)[whCC])
#     }
#   }
# )
