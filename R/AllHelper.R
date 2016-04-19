#' Helper methods for the ClusterExperiment class
#'
#' This is a collection of helper methods for the ClusterExperiment class.
#' @name ClusterExperiment-methods
#' @aliases ClusterExperiment-methods
#' @details Note that when subsetting the data, the dendrogram information and
#' the co-clustering matrix are lost.
#' @export
setMethod(
  f = "[",
  signature = c("ClusterExperiment", "ANY", "ANY"),
  definition = function(x, i, j, ..., drop=TRUE) {
    origN<-NCOL(x)
    out <- callNextMethod()
    out@clusterMatrix <- as.matrix(x@clusterMatrix[j, ,drop=FALSE])
    out@coClustering <- new("matrix") ###Need to think about this
    out@dendro_samples <- NULL
    out@dendro_clusters <- NULL

   # browser()
    out@orderSamples<-match(out@orderSamples[j],c(1:origN)[j])

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
    validObject(out)
    return(out)
  }
)

## show
#' @rdname ClusterExperiment-methods
#' @export
#' @param object a ClusterExperiment object.
setMethod(
  f = "show",
  signature = "ClusterExperiment",
  definition = function(object) {
    cat("class:", class(object), "\n")
    cat("dim:", dim(object), "\n")
     cat("Primary cluster type:", clusterType(object)[primaryClusterIndex(object)],"\n")
    cat("Primary cluster label:", clusterLabels(object)[primaryClusterIndex(object)],"\n")
    cat("Table of clusters (of primary clustering):")
    print(table(primaryClusterNamed(object)))
    cat("Total number of clusterings:", NCOL(clusterMatrix(object)),"\n")
    typeTab<-names(table(clusterType(object)))
    cat("clusterMany run?",if("clusterMany" %in% typeTab) "Yes" else "No","\n")
    cat("combineMany run?",if("combineMany" %in% typeTab) "Yes" else "No","\n")
    cat("mergeClusters run?",if("mergeClusters" %in% typeTab) "Yes" else "No","\n")
  }
)

#' @rdname ClusterExperiment-methods
#' @details \code{clusterMatrixNamed} returns a matrix with cluster labels.
#' @export
#' @param x a ClusterExperiment object.
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
#' @details \code{primaryClusterNamed} returns the primary cluster (using cluster
#' labels).
#' @export
setMethod(
  f = "primaryClusterNamed",
  signature = "ClusterExperiment",
  definition = function(x) {
    clusterMatrixNamed(x)[,primaryClusterIndex(x)]
  })

#' @rdname ClusterExperiment-methods
#' @details \code{transformation} prints the function used to transform the data
#' prior to clustering.
#' @export
setMethod(
  f = "transformation",
  signature = "ClusterExperiment",
  definition = function(x) {
    return(x@transformation)
  }
)

#' @export
#' @rdname clusterLabels
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

#' Methods to set and return the clustering labels
#'
#' These methods will return or set the clustering labels (i.e., the
#' colnames of the clusterMatrix slot).
#'
#' @param `whichClusters` controls which labels to be returned. It is either
#' numeric, in which case gives the indices of the clusters, or character, in
#' which case it matches to clusterType(x) to find the indices of the clusters.
#' @export
#' @rdname clusterLabels
setMethod(
  f = "clusterLabels",
  signature = signature(x = "ClusterExperiment",whichClusters="numeric"),
  definition = function(x, whichClusters){
    if(!all(whichClusters %in% 1:NCOL(clusterMatrix(x)))) stop("Invalid indices for clusterLabels")
    labels<-colnames(clusterMatrix(x))[whichClusters]
    if(is.null(labels)) cat("No labels found for clusterings\n")
    return(labels)
  }
)

#' @rdname clusterLabels
#' @export
setMethod(
  f = "clusterLabels",
  signature = signature(x = "ClusterExperiment", whichClusters ="character"),
  definition = function(x, whichClusters="all"){
    wh<-.TypeIntoIndices(x,whClusters=whichClusters)
    return(clusterLabels(x,wh))
  }
)
#' @rdname clusterLabels
#' @export
setMethod(
  f = "clusterLabels",
  signature = signature(x = "ClusterExperiment",whichClusters="missing"),
  definition = function(x, whichClusters){
    return(clusterLabels(x,whichClusters="all"))
  }
)
#' @rdname ClusterExperiment-methods
#' @details \code{nClusters} returns the number of clusterings (i.e., ncol of
#' clusterMatrix).
#' @export
setMethod(
  f = "nClusters",
  signature = "ClusterExperiment",
  definition = function(x){
    return(NCOL(clusterMatrix(x)))
  }
)
#' @rdname ClusterExperiment-methods
#' @details \code{nFeatures} returns the number of features (same as `nrow`).
#' @export
setMethod(
  f = "nFeatures",
  signature =  "ClusterExperiment",
  definition = function(x){
    return(NROW(assay(x)))
  }
)

#' @rdname ClusterExperiment-methods
#' @details \code{nSamples} returns the number of samples (same as `ncol`).
#' @export
setMethod(
  f = "nSamples",
  signature = "ClusterExperiment",
  definition = function(x){
    return(NCOL(assay(x)))
  }
)

#' @rdname ClusterExperiment-methods
#' @details \code{clusterMatrix} returns the matrix with all the clusterings.
#' @export
setMethod(
  f = "clusterMatrix",
  signature = "ClusterExperiment",
  definition = function(x) {
    return(x@clusterMatrix)
  }
)

#' @rdname ClusterExperiment-methods
#' @details \code{primaryCluster} returns the primary clustering (as numeric).
#' @export
setMethod(
  f = "primaryCluster",
  signature = "ClusterExperiment",
  definition = function(x) {
    return(x@clusterMatrix[,primaryClusterIndex(x)])
  }
)

#' @rdname ClusterExperiment-methods
#' @details \code{primaryClusterIndex} returns/sets the primary clustering index
#' (i.e., which column of clusterMatrix corresponds to the primary clustering).
#' @export
setMethod(
  f = "primaryClusterIndex",
  signature = "ClusterExperiment",
  definition = function(x) {
    return(x@primaryIndex)
  }
)

#' @rdname ClusterExperiment-methods
#' @export
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
#' @details \code{coClustering} returns/sets the co-clustering matrix.
#' @export
setMethod(
  f = "coClustering",
  signature = "ClusterExperiment",
  definition = function(x) {
    return(x@coClustering)
  }
)

#' @rdname ClusterExperiment-methods
#' @export
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
#' @details \code{clusterType} returns/sets the clusterType slot.
#' @export
setMethod(
  f = "clusterType",
  signature = "ClusterExperiment",
  definition = function(x) {
    out<-x@clusterType
    names(out)<-clusterLabels(x)
    return(out)
  }
)


#Update here if change pipeline values. Also defines the order of them.
.pipelineValues<-c("final","mergeClusters","combineMany","clusterMany")
#' @rdname pipelineClusters
#' @export
setMethod(
  f = "pipelineClusterDetails",
  signature = signature("ClusterExperiment"),
  definition = function(x) {

    if(length(clusterType(x))!=NCOL(clusterMatrix(x))) stop("Invalid ClusterExperiment object")
    #check if old iterations already exist; note assumes won't have previous iteration unless have current one.
    existingOld<-lapply(.pipelineValues,function(ch){
      regex<-paste(ch,"_",sep="")
      grep(regex,clusterType(x))

    })
    st<-strsplit(clusterType(x)[unlist(existingOld)],"_")
    oldValues<-data.frame(index=unlist(existingOld),type=sapply(st,.subset2,1),iteration=as.numeric(sapply(st,.subset2,2)),stringsAsFactors=FALSE)

    wh<-which(clusterType(x) %in% .pipelineValues) #current iteration
    if(length(wh)>0){
      existingValues<-data.frame(index=wh,type=clusterType(x)[wh], iteration=0,stringsAsFactors=FALSE) #0 indicates current iteration
      if(nrow(oldValues)>0) existingValues<-rbind(oldValues,existingValues)
    }
    else{
      if(nrow(oldValues)>0) existingValues<-oldValues
      else   return(NULL)
    }

    return(existingValues)

  }
)
#' @rdname pipelineClusters
#' @export
setMethod(
  f = "pipelineClusterTable",
  signature = signature("ClusterExperiment"),
  definition = function(x){
    ppIndex<-pipelineClusterDetails(x)
    table(Type=factor(ppIndex[,"type"],levels=.pipelineValues),Iteration=factor(ppIndex[,"iteration"]))
  }
)
#' @rdname pipelineClusters
#' @title Methods for pipeline clusters
#' @name pipelineClusters
#' @aliases pipelineClusters pipelineClusterTable pipelineClusterDetails
#' @description  The main pipeline of the package is made of \code{\link{clusterMany}},
#' \code{\link{combineMany}}, and \code{\link{mergeClusters}}.
#' The clusterings from these functions (and not those obtained in a different
#' way) can be obtained with the functions documented here.
#' @export
setMethod(
  f = "pipelineClusters",
  signature = signature("ClusterExperiment"),
  definition = function(x,iteration=0) {
    ppIndex<-pipelineClusterDetails(x)
    if(is.na(iteration)) iteration<-unique(ppIndex[,"iteration"])
    if(!is.null(ppIndex)){
      whIteration<-which(ppIndex[,"iteration"]%in%iteration)
      if(length(whIteration)>0){
        index<-ppIndex[whIteration,"index"]
        return(clusterMatrix(x)[,index,drop=FALSE])
      }
      else return(NULL)
    }
    else return(NULL)
}
)


#' @rdname ClusterExperiment-methods
#' @details \code{clusterInfo} returns the clusterInfo slot.
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
#' @details \code{clusterLegend} returns/sets the clusterLegend slot.
#' @export
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
#' @details \code{orderSamples} returns/sets the orderSamples slot.
#' @export
setMethod(
    f = "orderSamples",
    signature = "ClusterExperiment",
    definition = function(x) {
        return(x@orderSamples)
    }
)
#' @rdname ClusterExperiment-methods
#' @export
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
setReplaceMethod(
    f = "clusterType",
    signature = signature(object="ClusterExperiment", value="character"),
    definition = function(object,value) {
        object@clusterType<-value
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
#     whCC<-which(clusterType(x)==type)
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
