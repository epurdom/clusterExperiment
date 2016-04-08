## subsetting
setMethod(
  f = "[",
  signature = c("ClusterCells", "ANY", "ANY"),
  definition = function(x, i, j, ..., drop=TRUE) {
    out <- callNextMethod()
    out@clusterLabels <- as.matrix(x@clusterLabels[j,])
    out@coClustering <- new("matrix")
    out@dendrogram <- list()
    return(out)
  }
)

## show
#' @rdname clusterCells
setMethod(
  f = "show",
  signature = "ClusterCells",
  definition = function(object) {
    cat("class:", class(object), "\n")
    cat("dim:", dim(object), "\n")
    cat("Table of clusters (of primary clustering):")
    print(table(object@clusterLabels[,object@primaryIndex]))
    cat("Primary cluster type:", object@clusterType[object@primaryIndex],"\n")
    cat("Total number of clusterings:", NCOL(object@clusterLabels),"\n")
    typeTab<-names(table(clusterType(object)))
        cat("compareChoices run?",if("compareChoices" %in% typeTab) "Yes" else "No","\n")
        cat("findSharedClusters run?",if("findSharedClusters" %in% typeTab) "Yes" else "No","\n")
        cat("mergeClusters run?",if("mergeClusters" %in% typeTab) "Yes" else "No","\n")
  }
)

# #' @rdname clusterCells
# setMethod(
#   f = "isLog",
#   signature = "ClusterCells",
#   definition = function(x) {
#     return(x@isLog)
#   }
# )


#' @rdname clusterCells
setMethod(
  f = "transformation",
  signature = "ClusterCells",
  definition = function(x) {
    return(x@transformation)
  }
)



#' @rdname clusterCells
setMethod(
  f = "allClusters",
  signature = "ClusterCells",
  definition = function(x) {
    return(x@clusterLabels)
  }
)

#' @rdname clusterCells
setMethod(
  f = "primaryCluster",
  signature = "ClusterCells",
  definition = function(x) {
    return(x@clusterLabels[,x@primaryIndex])
  }
)

#' @rdname clusterCells
setMethod(
  f = "primaryClusterIndex",
  signature = "ClusterCells",
  definition = function(x) {
    return(x@primaryIndex)
  }
)

#' @rdname clusterCells
setReplaceMethod(
  f = "primaryClusterIndex",
  signature = signature("ClusterCells", "numeric"),
  definition = function(object, value) {
    object@primaryIndex <- value
    validObject(object)
    return(object)
  }
)

#' @rdname clusterCells
setMethod(
  f = "coClustering",
  signature = "ClusterCells",
  definition = function(x) {
    return(x@coClustering)
  }
)

#' @rdname clusterCells
setMethod(
  f = "dendrogram",
  signature = "ClusterCells",
  definition = function(x) {
    return(x@dendrogram)
  }
)

#' @rdname clusterCells
setMethod(
  f = "clusterType",
  signature = "ClusterCells",
  definition = function(x) {
    return(x@clusterType)
  }
)

#' @rdname clusterCells
setMethod(
  f = "removeClusters",
  signature = signature("ClusterCells","character"),
  definition = function(x, whichRemove,exactMatch=TRUE) {
    if(exactMatch) wh<-which(clusterType(x) %in% whichRemove)
    else{
      sapply(whichRemove,grep, clusterType(x))
    }
    removeClusters(x,wh)
  }
)
#' @rdname clusterCells
setMethod(
  f = "removeClusters",
  signature = signature("ClusterCells","numeric"),
  definition = function(x, whichRemove) {
   
    if(any(whichRemove>NCOL(allClusters(x)))) stop("invalid indices -- must be between 1 and",NCOL(allClusters(x)))
    if(length(whichRemove)==NCOL(allClusters(x))){ 
      warning("All clusters have been removed. Will return just a Summarized Experiment Object")
      #make it Summarized Experiment
    }
    newClLabels<-allClusters(x)[,-whichRemove,drop=FALSE]
    newClusterInfo<-clusterInfo(x)[-whichRemove,drop=FALSE]
    newClusterType<-clusterType(x)[-whichRemove,drop=FALSE]
    if(primaryClusterIndex(x) %in% whichRemove) primaryIndex<-1
    else primaryIndex<-match(primaryClusterIndex(x),1:NCOL(allClusters(x))[-whichRemove])
    retval<-clusterCells(assay(x),newClLabels[,1],transformation(x))
    retval@clusterLabels<-newClLabels
    retval@primaryIndex<-primaryIndex
    retval@clusterInfo<-newClusterInfo
    retval@clusterType<-newClusterType
    return(retval)
  }
)

#' @rdname clusterCells
setMethod(
  f = "addClusters",
  signature = signature("ClusterCells", "matrix"),
  definition = function(x, y, type="User") {
    if(!(NROW(y) == NCOL(x))) {
      stop("Incompatible dimensions.")
    }
    x@clusterLabels <- cbind(x@clusterLabels, y)
    if(length(type)==1) type<-rep(type, NCOL(y))
    x@clusterType <- c(x@clusterType, type)
    yClusterInfo<-rep(list(NULL),NCOL(y))
    x@clusterInfo<-c(x@clusterInfo,yClusterInfo)
    return(x)
  }
)

#' @rdname clusterCells
setMethod(
  f = "pipelineClusterIndex",
  signature = signature("ClusterCells"),
  definition = function(x, printTable=TRUE) {
    pipelineValues<-c("compareChoices","findSharedClusters","mergeClusters")
    if(length(clusterType(x))!=NCOL(allClusters(x))) stop("Invalid ClusterCells object")
    #check if old iterations already exist; note assumes won't have previous iteration unless have current one. 
    existingOld<-lapply(pipelineValues,function(ch){
      regex<-paste(ch,"_",sep="")
      grep(regex,clusterType(x))
      
    })
    st<-strsplit(clusterType(x)[unlist(existingOld)],"_")
    oldValues<-data.frame(index=unlist(existingOld),type=sapply(st,.subset2,1),iteration=as.numeric(sapply(st,.subset2,2)),stringsAsFactors=FALSE)
    
    wh<-which(clusterType(x) %in% pipelineValues) #current iteration
    if(length(wh)>0){
      existingValues<-data.frame(index=wh,type=clusterType(x)[wh], iteration=0,stringsAsFactors=FALSE) #0 indicates current iteration
      if(nrow(oldValues)>0) existingValues<-rbind(oldValues,existingValues)
    }
    else{
      if(nrow(oldValues)>0) existingValues<-oldValues
      else   return(NULL)
    }
    
    if(printTable) print(table(Type=factor(existingValues[,"type"],levels=pipelineValues),Iteration=factor(existingValues[,"iteration"])))
    invisible(existingValues)
    
  }
)
#' @rdname clusterCells
setMethod(
  f = "pipelineClusters",
  signature = signature("ClusterCells"),
  definition = function(x,iteration=0) {
    ppIndex<-pipelineClusterIndex(x,print=FALSE)
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


#' @rdname clusterCells
setMethod(
  f = "addClusters",
  signature = signature("ClusterCells", "numeric"),
  definition = function(x, y, type="User") {
    if(!(length(y) == NCOL(x))) {
      stop("Incompatible dimensions.")
    }
    x@clusterLabels <- cbind(x@clusterLabels, y)
    if(length(type)==1) type<-rep(type, 1)
    yClusterInfo<-rep(list(NULL),1)
    x@clusterInfo<-c(x@clusterInfo,yClusterInfo)
    x@clusterType <- c(x@clusterType, type)
    return(x)
  }
)

#' @rdname clusterCells
setMethod(
  f = "addClusters",
  signature = signature("ClusterCells", "ClusterCells"),
  definition = function(x, y) {
    if(!all(assay(y) == assay(x))) {
      stop("Cannot merge clusters from different data.")
    }
    x@clusterLabels <- cbind(x@clusterLabels, y@clusterLabels)
    x@clusterType <- c(x@clusterType, y@clusterType)
    x@clusterInfo<-c(x@clusterInfo,y@clusterInfo)
    return(x)
  }
)

#' @rdname clusterCells
setMethod(
  f = "removeUnclustered",
  signature = "ClusterCells",
  definition = function(x) {
    return(x[,primaryCluster(x) > 0])
  }
)

#' @rdname clusterCells
setMethod(
  f = "clusterInfo",
  signature = "ClusterCells",
  definition = function(x) {
    return(x@clusterInfo)
  }
)
