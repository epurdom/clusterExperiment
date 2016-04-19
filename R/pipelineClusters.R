#Update here if change pipeline values. Also defines the order of them.
.pipelineValues<-c("final","mergeClusters","combineMany","clusterMany")

#' @title Methods for pipeline clusters
#'
#' @description The main pipeline of the package is made of
#'   \code{\link{clusterMany}}, \code{\link{combineMany}}, and
#'   \code{\link{mergeClusters}}. The clusterings from these functions (and not
#'   those obtained in a different way) can be obtained with the functions
#'   documented here.
#' @param x a \code{\link{ClusterExperiment}} object.
#'
#' @return \code{pipelineClusters} returns a matrix consisting of the
#'   appropriate columns of the \code{clusterMatrix} slot.
#' @name pipelineClusters
#' @aliases pipelineClusters pipelineClusterTable pipelineClusterDetails
#' @export
#' @examples
#' data(simData)
#'
#' cl <- clusterMany(simData,nPCADims=c(5,10,50),  dimReduce="PCA",
#' clusterFunction="pam", ks=2:4, findBestK=c(FALSE), removeSil=TRUE,
#' subsample=FALSE)
#'
#' clCommon <- combineMany(cl, whichClusters="pipeline", proportion=0.7,
#' minSize=10)
#'
#' clCommon <- makeDendrogram(clCommon)
#'
#' clMerged <- mergeClusters(clCommon)
#'
#' head(pipelineClusters(clMerged))
#' pipelineClusterDetails(clMerged)
#' pipelineClusterTable(clMerged)
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


#' @rdname pipelineClusters
#' @return \code{pipelineClusterDetails} returns a \code{data.frame} with some
#'   details on the clusterings, such as the type (e.g., `clusterMany`,
#'   `combineMany`) and iteration.
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
#' @return \code{pipelineClusterTable} returns a table of how many of the
#'   clusterings belong to each of the following possible values: `final`,
#'   `mergeClusters`, `combineMany` and `clusterMany`.
#' @export
setMethod(
  f = "pipelineClusterTable",
  signature = signature("ClusterExperiment"),
  definition = function(x){
    ppIndex<-pipelineClusterDetails(x)
    table(Type=factor(ppIndex[,"type"],levels=.pipelineValues),Iteration=factor(ppIndex[,"iteration"]))
  }
)
