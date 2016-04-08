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
#' @rdname ClusterCells-class
setMethod(
  f = "show",
  signature = "ClusterCells",
  definition = function(object) {
    cat("class:", class(object), "\n")
    cat("dim:", dim(object), "\n")
    cat("Table of clusters (of primary clustering):")
    print(table(allClusters(object)[,primaryClusterIndex(object)]))
    cat("Primary cluster type:", clusterType(object)[primaryClusterIndex(object)],"\n")
    cat("Total number of clusterings:", NCOL(allClusters(object)),"\n")
    typeTab<-names(table(clusterType(object)))
        cat("compareChoices run?",if("compareChoices" %in% typeTab) "Yes" else "No","\n")
        cat("findSharedClusters run?",if("findSharedClusters" %in% typeTab) "Yes" else "No","\n")
        cat("mergeClusters run?",if("mergeClusters" %in% typeTab) "Yes" else "No","\n")
  }
)

# #' @rdname ClusterCells-class
# setMethod(
#   f = "isLog",
#   signature = "ClusterCells",
#   definition = function(x) {
#     return(x@isLog)
#   }
# )


#' @rdname ClusterCells-class
setMethod(
  f = "transformation",
  signature = "ClusterCells",
  definition = function(x) {
    return(x@transformation)
  }
)



#' @rdname ClusterCells-class
setMethod(
  f = "allClusters",
  signature = "ClusterCells",
  definition = function(x) {
    return(x@clusterLabels)
  }
)

#' @rdname ClusterCells-class
setMethod(
  f = "primaryCluster",
  signature = "ClusterCells",
  definition = function(x) {
    return(x@clusterLabels[,primaryClusterIndex(x)])
  }
)

#' @rdname ClusterCells-class
setMethod(
  f = "primaryClusterIndex",
  signature = "ClusterCells",
  definition = function(x) {
    return(x@primaryIndex)
  }
)

#' @rdname ClusterCells-class
setReplaceMethod(
  f = "primaryClusterIndex",
  signature = signature("ClusterCells", "numeric"),
  definition = function(object, value) {
    primaryClusterIndex(object) <- value
    validObject(object)
    return(object)
  }
)

#' @rdname ClusterCells-class
setMethod(
  f = "coClustering",
  signature = "ClusterCells",
  definition = function(x) {
    return(x@coClustering)
  }
)

#' @rdname ClusterCells-class
setMethod(
  f = "dendrogram",
  signature = "ClusterCells",
  definition = function(x) {
    return(x@dendrogram)
  }
)

#' @rdname ClusterCells-class
setMethod(
  f = "clusterType",
  signature = "ClusterCells",
  definition = function(x) {
    return(x@clusterType)
  }
)

#' @rdname ClusterCells-class
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
#' @rdname ClusterCells-class
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
    primaryClusterIndex(retval)<-primaryIndex
    retval@clusterInfo<-newClusterInfo
    retval@clusterType<-newClusterType
    return(retval)
  }
)

#' @rdname ClusterCells-class
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
    validObject(x)
    return(x)
  }
)

#' @rdname ClusterCells-class
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
#' @rdname ClusterCells-class
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


#' @rdname ClusterCells-class
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
    validObject(x)
    return(x)
  }
)

#' @rdname ClusterCells-class
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
    validObject(x)
    return(x)
  }
)

#' @rdname ClusterCells-class
setMethod(
  f = "removeUnclustered",
  signature = "ClusterCells",
  definition = function(x) {
    return(x[,primaryCluster(x) > 0])
  }
)

#' @rdname ClusterCells-class
setMethod(
  f = "clusterInfo",
  signature = "ClusterCells",
  definition = function(x) {
    return(x@clusterInfo)
  }
)

# # Need to implement: wrapper to get a nice summary of the parameters choosen, similar to that of paramMatrix of compareChoices (and compatible with it)
# #' @rdname ClusterCells-class
# setMethod(
#   f= "paramValues",
#   signature = "ClusterCells",
#   definition=function(x,type){
#     whCC<-which(clusterType(x)==type)
#     if(length(wwCC)==0) stop("No clusterings of type equal to ",type,"are found")
#     if(type=="compareChoices"){
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