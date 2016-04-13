## subsetting
#' @rdname ClusterExperiment-class
setMethod(
  f = "[",
  signature = c("ClusterExperiment", "ANY", "ANY"),
  definition = function(x, i, j, ..., drop=TRUE) {
    out <- callNextMethod()
    out@clusterMatrix <- as.matrix(x@clusterMatrix[j,])
    out@coClustering <- new("matrix")
    out@dendrogram <- list()
    return(out)
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
    cat("Table of clusters (of primary clustering):")
    print(table(allClusters(object)[,primaryClusterIndex(object)]))
    cat("Primary cluster type:", clusterType(object)[primaryClusterIndex(object)],"\n")
    cat("Total number of clusterings:", NCOL(allClusters(object)),"\n")
    typeTab<-names(table(clusterType(object)))
        cat("clusterMany run?",if("clusterMany" %in% typeTab) "Yes" else "No","\n")
        cat("findSharedClusters run?",if("findSharedClusters" %in% typeTab) "Yes" else "No","\n")
        cat("mergeClusters run?",if("mergeClusters" %in% typeTab) "Yes" else "No","\n")
  }
)

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
    if(length(value)!=NCOL(allClusters(object))) stop("value must be a vector of length equal to NCOL(allClusters(object)):",NCOL(allClusters(object)))
#note, don't currently require unique labels. Probably best, since mainly used for plotting
    colnames(object@clusterMatrix) <- value
    validObject(object)
    return(object)
  }
)

#' @rdname ClusterExperiment-class
setMethod(
  f = "clusterLabels",
  signature = signature(x = "ClusterExperiment",whichClusters="numeric"),
  definition = function(x, whichClusters){
    if(!all(whichClusters %in% 1:NCOL(allClusters(x)))) stop("Invalid indices for clusterLabels")
    labels<-colnames(allClusters(x))[whichClusters]
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
    return(NCOL(allClusters(x)))
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
  f = "allClusters",
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
  f = "dendrogram",
  signature = "ClusterExperiment",
  definition = function(x) {
    return(x@dendrogram)
  }
)

#' @rdname ClusterExperiment-class
setMethod(
  f = "clusterType",
  signature = "ClusterExperiment",
  definition = function(x) {
    return(x@clusterType)
  }
)

#' @rdname ClusterExperiment-class
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
#' @rdname ClusterExperiment-class
setMethod(
  f = "removeClusters",
  signature = signature("ClusterExperiment","numeric"),
  definition = function(x, whichRemove) {
   #browser()
    if(any(whichRemove>NCOL(allClusters(x)))) stop("invalid indices -- must be between 1 and",NCOL(allClusters(x)))
    if(length(whichRemove)==NCOL(allClusters(x))){
      warning("All clusters have been removed. Will return just a Summarized Experiment Object")
      #make it Summarized Experiment
    }
    newClLabels<-allClusters(x)[,-whichRemove,drop=FALSE]
    newClusterInfo<-clusterInfo(x)[-whichRemove]
    newClusterType<-clusterType(x)[-whichRemove]
    newClusterColors<-clusterColors(x)[-whichRemove]
    if(primaryClusterIndex(x) %in% whichRemove) pIndex<-1
    else pIndex<-match(primaryClusterIndex(x),1:NCOL(allClusters(x))[-whichRemove])
    retval<-clusterExperiment(assay(x),newClLabels,transformation(x),clusterType=newClusterType,clusterInfo<-newClusterInfo)
    validObject(retval)
    clusterColors(retval)<-newClusterColors
    primaryClusterIndex(retval)<-pIndex #Note can only set it on valid object so put it here...
    return(retval)
  }
)



#Update here if change pipeline values. Also defines the order of them.
.pipelineValues<-c("final","mergeClusters","findSharedClusters","clusterMany")
#' @rdname pipelineClusters
setMethod(
  f = "pipelineClusterDetails",
  signature = signature("ClusterExperiment"),
  definition = function(x) {

    if(length(clusterType(x))!=NCOL(allClusters(x))) stop("Invalid ClusterExperiment object")
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
        return(allClusters(x)[,index,drop=FALSE])
      }
      else return(NULL)
    }
    else return(NULL)
}
)

#' @rdname ClusterExperiment-class
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
        x@clusterColors<-c(x@clusterColors,y@clusterColors)
        validObject(x)
        return(x)
    }
)
#' @rdname ClusterExperiment-class
setMethod(
    f = "addClusters",
    signature = signature("ClusterExperiment", "matrix"),
    definition = function(x, y, type="User") {
        ccObj<-clusterExperiment(assay(x),y,transformation=transformation(x),clusterType=type)
        addClusters(x,ccObj)
    }
)
#' @rdname ClusterExperiment-class
setMethod(
  f = "addClusters",
  signature = signature("ClusterExperiment", "numeric"),
  definition = function(x, y, ...) {
    addClusters(x,matrix(y,ncol=1),...)
  }
)



#' @rdname ClusterExperiment-class
setMethod(
  f = "removeUnclustered",
  signature = "ClusterExperiment",
  definition = function(x) {
    return(x[,primaryCluster(x) > 0])
  }
)

#' @rdname ClusterExperiment-class
setMethod(
  f = "clusterInfo",
  signature = "ClusterExperiment",
  definition = function(x) {
    return(x@clusterInfo)
  }
)
#need to implement replacement functions for these clusterInfo, clusterType, clusterColors

#' @rdname ClusterExperiment-class
setMethod(
    f = "clusterColors",
    signature = "ClusterExperiment",
    definition = function(x) {
        return(x@clusterColors)
    }
)
#' @rdname ClusterExperiment-class
setReplaceMethod(
    f = "clusterColors",
    signature = signature(object="ClusterExperiment", value="list"),
    definition = function(object, value) {
        object@clusterColors<-value
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
#     else if(type=="clusterAll"){
#
#     }
#     else{
#       return(clusterInfo(x)[whCC])
#     }
#   }
# )
