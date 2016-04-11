## subsetting
setMethod(
  f = "[",
  signature = c("ClusterCells", "ANY", "ANY"),
  definition = function(x, i, j, ..., drop=TRUE) {
    out <- callNextMethod()
    out@clusterMatrix <- as.matrix(x@clusterMatrix[j,])
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
        cat("clusterMany run?",if("clusterMany" %in% typeTab) "Yes" else "No","\n")
        cat("findSharedClusters run?",if("findSharedClusters" %in% typeTab) "Yes" else "No","\n")
        cat("mergeClusters run?",if("mergeClusters" %in% typeTab) "Yes" else "No","\n")
  }
)

#' @rdname ClusterCells-class
setMethod(
  f = "transformation",
  signature = "ClusterCells",
  definition = function(x) {
    return(x@transformation)
  }
)

#' @rdname ClusterCells-class
setReplaceMethod(
  f = "clusterLabels",
  signature = signature(object="ClusterCells", value="character"),
  definition = function(object, value) {
    if(length(value)!=NCOL(allClusters(object))) stop("value must be a vector of length equal to NCOL(allClusters(object)):",NCOL(allClusters(object)))
#note, don't currently require unique labels. Probably best, since mainly used for plotting
    colnames(object@clusterMatrix) <- value
    validObject(object)
    return(object)
  }
)

#' @rdname ClusterCells-class
setMethod(
  f = "clusterLabels",
  signature = signature(x = "ClusterCells",whichClusters="numeric"),
  definition = function(x, whichClusters){
    if(!all(whichClusters %in% 1:NCOL(allClusters(x)))) stop("Invalid indices for clusterLabels")
    labels<-colnames(allClusters(x))[whichClusters]
    if(is.null(labels)) cat("No labels found for clusterings\n")
    return(labels)
  }
)
#' @rdname ClusterCells-class
setMethod(
  f = "clusterLabels",
  signature = signature(x = "ClusterCells", whichClusters ="character"),
  definition = function(x, whichClusters="all"){
    wh<-.TypeIntoIndices(x,whClusters=whichClusters)
    return(clusterLabels(x,wh))
  }
)
#' @rdname ClusterCells-class
setMethod(
  f = "clusterLabels",
  signature = signature(x = "ClusterCells",whichClusters="missing"),
  definition = function(x, whichClusters){
    return(clusterLabels(x,whichClusters="all"))
  }
)
#' @rdname ClusterCells-class
setMethod(
  f = "nClusters",
  signature = "ClusterCells",
  definition = function(x){
    return(NCOL(allClusters(x)))
  }
)
#' @rdname ClusterCells-class
setMethod(
  f = "nFeatures",
  signature =  "ClusterCells",
  definition = function(x){
    return(NROW(assay(x)))
  }
)
#' @rdname ClusterCells-class
setMethod(
  f = "nSamples",
  signature = "ClusterCells",
  definition = function(x){
    return(NCOL(assay(x)))
  }
)
#' @rdname ClusterCells-class
setMethod(
  f = "allClusters",
  signature = "ClusterCells",
  definition = function(x) {
    return(x@clusterMatrix)
  }
)

#' @rdname ClusterCells-class
setMethod(
  f = "primaryCluster",
  signature = "ClusterCells",
  definition = function(x) {
    return(x@clusterMatrix[,primaryClusterIndex(x)])
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
    object@primaryIndex <- value
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
   #browser()
    if(any(whichRemove>NCOL(allClusters(x)))) stop("invalid indices -- must be between 1 and",NCOL(allClusters(x)))
    if(length(whichRemove)==NCOL(allClusters(x))){ 
      warning("All clusters have been removed. Will return just a Summarized Experiment Object")
      #make it Summarized Experiment
    }
    newClLabels<-allClusters(x)[,-whichRemove,drop=FALSE]
    newClusterInfo<-clusterInfo(x)[-whichRemove,drop=FALSE]
    newClusterType<-clusterType(x)[-whichRemove,drop=FALSE]
    if(primaryClusterIndex(x) %in% whichRemove) pIndex<-1
    else pIndex<-match(primaryClusterIndex(x),1:NCOL(allClusters(x))[-whichRemove])
    retval<-clusterCells(assay(x),newClLabels[,1],transformation(x))
    retval@clusterMatrix<-newClLabels
    retval@clusterInfo<-newClusterInfo
    retval@clusterType<-newClusterType
    validObject(retval)
    primaryClusterIndex(retval)<-pIndex #Note can only set it on valid object so put it here...
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
    x@clusterMatrix <- cbind(x@clusterMatrix, y)
    if(length(type)==1) type<-rep(type, NCOL(y))
    x@clusterType <- c(x@clusterType, type)
    yClusterInfo<-rep(list(NULL),NCOL(y))
    x@clusterInfo<-c(x@clusterInfo,yClusterInfo)
    validObject(x)
    return(x)
  }
)

#Update here if change pipeline values. Also defines the order of them.
.pipelineValues<-c("final","mergeClusters","findSharedClusters","clusterMany")
#' @rdname pipelineClusters
setMethod(
  f = "pipelineClusterDetails",
  signature = signature("ClusterCells"),
  definition = function(x) {
    
    if(length(clusterType(x))!=NCOL(allClusters(x))) stop("Invalid ClusterCells object")
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
  signature = signature("ClusterCells"),
  definition = function(x){
    ppIndex<-pipelineClusterDetails(x)
    table(Type=factor(ppIndex[,"type"],levels=.pipelineValues),Iteration=factor(ppIndex[,"iteration"]))
  }
)
#' @rdname pipelineClusters
setMethod(
  f = "pipelineClusters",
  signature = signature("ClusterCells"),
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


#' @rdname ClusterCells-class
setMethod(
  f = "addClusters",
  signature = signature("ClusterCells", "numeric"),
  definition = function(x, y, type="User") {
    if(!(length(y) == NCOL(x))) {
      stop("Incompatible dimensions.")
    }
    x@clusterMatrix <- cbind(x@clusterMatrix, y)
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
    x@clusterMatrix <- cbind(x@clusterMatrix, y@clusterMatrix)
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

# # Need to implement: wrapper to get a nice summary of the parameters choosen, similar to that of paramMatrix of clusterMany (and compatible with it)
# #' @rdname ClusterCells-class
# setMethod(
#   f= "paramValues",
#   signature = "ClusterCells",
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