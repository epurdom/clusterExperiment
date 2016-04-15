## subsetting
#' @rdname ClusterExperiment-class
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
#' @details removeUnclustered removes all samples that are unclustered
#' (i.e. -1 or -2 assignment) in the primaryCluster of x (so they may be
#' unclustered in other clusters found in clusterMatrix(x))
#' @rdname addClusters
setMethod(
    f = "removeUnclustered",
    signature = "ClusterExperiment",
    definition = function(x) {
        return(x[,primaryCluster(x) >= 0])
    }
)
## show
#' @rdname ClusterExperiment-class
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

#' @rdname ClusterExperiment-class
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

#' @rdname ClusterExperiment-class
setMethod(
  f = "primaryClusterNamed",
  signature = "ClusterExperiment",
  definition = function(x) {
    clusterMatrixNamed(x)[,primaryClusterIndex(x)]
  })

#' @rdname ClusterExperiment-class
setMethod(
  f = "transformation",
  signature = "ClusterExperiment",
  definition = function(x) {
    return(x@transformation)
  }
)

#' @rdname ClusterExperiment-class
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

#' @rdname ClusterExperiment-class
#' @param whichClusters either numeric, in which case gives the indices of the clusters, or character, in which case it matches to clusterType(x) to find the indices of the clusters
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
#' @rdname ClusterExperiment-class
setMethod(
  f = "clusterLabels",
  signature = signature(x = "ClusterExperiment", whichClusters ="character"),
  definition = function(x, whichClusters="all"){
    wh<-.TypeIntoIndices(x,whClusters=whichClusters)
    return(clusterLabels(x,wh))
  }
)
#' @rdname ClusterExperiment-class
setMethod(
  f = "clusterLabels",
  signature = signature(x = "ClusterExperiment",whichClusters="missing"),
  definition = function(x, whichClusters){
    return(clusterLabels(x,whichClusters="all"))
  }
)
#' @rdname ClusterExperiment-class
setMethod(
  f = "nClusters",
  signature = "ClusterExperiment",
  definition = function(x){
    return(NCOL(clusterMatrix(x)))
  }
)
#' @rdname ClusterExperiment-class
setMethod(
  f = "nFeatures",
  signature =  "ClusterExperiment",
  definition = function(x){
    return(NROW(assay(x)))
  }
)
#' @rdname ClusterExperiment-class
setMethod(
  f = "nSamples",
  signature = "ClusterExperiment",
  definition = function(x){
    return(NCOL(assay(x)))
  }
)
#' @rdname ClusterExperiment-class
setMethod(
  f = "clusterMatrix",
  signature = "ClusterExperiment",
  definition = function(x) {
    return(x@clusterMatrix)
  }
)

#' @rdname ClusterExperiment-class
setMethod(
  f = "primaryCluster",
  signature = "ClusterExperiment",
  definition = function(x) {
    return(x@clusterMatrix[,primaryClusterIndex(x)])
  }
)

#' @rdname ClusterExperiment-class
setMethod(
  f = "primaryClusterIndex",
  signature = "ClusterExperiment",
  definition = function(x) {
    return(x@primaryIndex)
  }
)

#' @rdname ClusterExperiment-class
setReplaceMethod(
  f = "primaryClusterIndex",
  signature = signature("ClusterExperiment", "numeric"),
  definition = function(object, value) {
    object@primaryIndex <- value
    validObject(object)
    return(object)
  }
)

#' @rdname ClusterExperiment-class
setMethod(
  f = "coClustering",
  signature = "ClusterExperiment",
  definition = function(x) {
    return(x@coClustering)
  }
)

#' @rdname ClusterExperiment-class
setReplaceMethod(
  f = "coClustering",
  signature = signature(object="ClusterExperiment", value="matrix"),
  definition = function(object, value) {
    object@coClustering <- value
    validObject(object)
    return(object)
  }
)

#' @rdname ClusterExperiment-class
setMethod(
  f = "clusterType",
  signature = "ClusterExperiment",
  definition = function(x) {
    out<-x@clusterType
    names(out)<-clusterLabels(x)
    return(out)
  }
)

#' @rdname addClusters
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
#' @param exactMatch logical. Whether the whichRemove must exactly match a value
#' of clusterType(x). Only relevant if whichRemove is character.
#' @param whichRemove which clusters to remove. Can be numeric or character. If
#' numeric, must give indices of clusterMatrix(x) to remove. If character,
#' should match a clusterType of x
#' @details removeClusters removes the clusters given by whichRemove. If all
#' clusters are implied, then returns a SummarizedExperiment Object. If the
#' primaryCluster is one of the clusters removed, the primaryClusterIndex is set
#' to 1 and the dendrogram and cooccurance matrix are discarded and orderSamples
#'is set to 1:NCOL(x).
#' @rdname addClusters
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
    retval<-clusterExperiment(assay(x),newClLabels,transformation(x),clusterType=newClusterType,clusterInfo<-newClusterInfo)
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



#Update here if change pipeline values. Also defines the order of them.
.pipelineValues<-c("final","mergeClusters","combineMany","clusterMany")
#' @rdname pipelineClusters
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
#'
#' The main pipeline of the package is made of \code{\link{clusterMany}},
#' \code{\link{combineMany}}, and \code{\link{mergeClusters}}.
#' The clusterings from these functions (and not those obtained in a different
#' way) can be obtained with the functions documented here.
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

#' @param x a ClusterExperiment Object
#' @param y additional clusters to add to x. Can be ClusterExperiment Object or
#' a matrix/vector of clusters
#' @details addClusters adds y to x, and is thus not symmetric in the two
#' arguments. In particular, the primaryCluster and all of its supporting
#' information (dendrogram, coClustering, and orderSamples) are all kept from
#' the x object, even if y is a ClusterExperiment.
#' @rdname addClusters
setMethod(
  f = "addClusters",
  signature = signature("ClusterExperiment", "ClusterExperiment"),
  definition = function(x, y) {
    if(!all(assay(y) == assay(x))) {
      stop("Cannot merge clusters from different data.")
    }
    x@clusterMatrix <- cbind(x@clusterMatrix, y@clusterMatrix)
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
#' @title Function to add clusters to an existing ClusterExperiment object
#' @name addClusters
#' @aliases addClusters removeClusters
setMethod(
    f = "addClusters",
    signature = signature("ClusterExperiment", "matrix"),
    definition = function(x, y, type="User") {
        ccObj<-clusterExperiment(assay(x),y,transformation=transformation(x),clusterType=type)
        addClusters(x,ccObj)
    }
)
#' @rdname addClusters
setMethod(
  f = "addClusters",
  signature = signature("ClusterExperiment", "numeric"),
  definition = function(x, y, ...) {
    addClusters(x,matrix(y,ncol=1),...)
  }
)





#' @rdname ClusterExperiment-class
setMethod(
  f = "clusterInfo",
  signature = "ClusterExperiment",
  definition = function(x) {
    out<-x@clusterInfo
    names(out)<-clusterLabels(x)
    return(out)
  }
)

#' @rdname ClusterExperiment-class
setMethod(
    f = "clusterLegend",
    signature = "ClusterExperiment",
    definition = function(x) {
      out<-x@clusterLegend
      names(out)<-clusterLabels(x)
      return(out)
    }
)
#' @rdname ClusterExperiment-class
setReplaceMethod(
    f = "clusterLegend",
    signature = signature(object="ClusterExperiment", value="list"),
    definition = function(object, value) {
        object@clusterLegend<-unname(value)
        validObject(object)
        return(object)
    }
)

#' @rdname ClusterExperiment-class
setMethod(
    f = "orderSamples",
    signature = "ClusterExperiment",
    definition = function(x) {
        return(x@orderSamples)
    }
)
#' @rdname ClusterExperiment-class
setReplaceMethod(
    f = "orderSamples",
    signature = signature(object="ClusterExperiment", value="numeric"),
    definition = function(object, value) {
        object@orderSamples<-value
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
