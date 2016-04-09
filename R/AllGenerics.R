setGeneric(
  name = "clusterAll",
  def = function(x,  subsample=TRUE, sequential=FALSE,
                 clusterFunction=c("tight", "hierarchical", "pam", "kmeans"),
                 clusterDArgs=NULL, subsampleArgs=NULL, seqArgs=NULL,
                 isCount=FALSE,transFun=NULL,dimReduce=c("none","PCA","mostVar"),ndims=NA) {
    standardGeneric("clusterAll")
  }
)
#' @rdname ClusterCells-class
setGeneric(
  name = "pipelineClusters",
  def = function(x,iteration=0) {
    standardGeneric("pipelineClusters")
  }
)
#' @rdname ClusterCells-class
setGeneric(
  name = "pipelineClusterIndex",
  def = function(x, printTable=TRUE) {
    standardGeneric("pipelineClusterIndex")
  }
)
#' @rdname ClusterCells-class
setGeneric(
  name = "paramValues",
  def = function(x,type) {
    standardGeneric("paramValues")
  }
)

#' @rdname compareChoices
setGeneric(
  name = "compareChoices",
  def = compareChoices <- function(x, ks=3:5, clusterMethod="pam", alphas=0.1, findBestK=FALSE,sequential=FALSE,
                                   removeSil=FALSE, subsample=FALSE,silCutoff=0,
                                   dimReduce="none",nVarDims=NA,nPCADims=NA,
                                   clusterDArgs=list(minSize=5),
                                   subsampleArgs=list(resamp.num=50),
                                   seqArgs=list(beta=0.9,k.min=3, verbose=FALSE),
                                   transFun=NULL,isCount=FALSE,
                                   ncores=1,random.seed=NULL,run=TRUE,paramMatrix=NULL,eraseOld=FALSE,...
                                   ) {
    standardGeneric("compareChoices")
  }
)

#' @rdname plotTracking
setGeneric(
  name="plotTracking",
  def=function(clusters, orderClusters=NULL,
               orderSamples=NULL,index=NULL,reuseColors=FALSE,matchToTop=FALSE,
               plot=TRUE,unassignedColor="white",missingColor="grey",
               minRequireColor=0.3,startNewColors=FALSE,
               colPalette=bigPalette,input=c("clusters","colors"),
               clNames=colnames(clusters),add=FALSE,xCoord=NULL,
               ylim=NULL,tick=FALSE,ylab="",xlab="",axisLine=0,box=FALSE,...)
    {
    standardGeneric("plotTracking")
  }
)
#' @rdname ClusterCells-class
setGeneric(
  name = "transformation",
  def = function(x) {
    standardGeneric("transformation")
  }
)
#' @rdname ClusterCells-class
  setGeneric(
  name = "transform",
  def = function(x,nPCADims=NA,nVarDims=NA,dimReduce="none") {
    standardGeneric("transform")
  }
)

#' @rdname ClusterCells-class
setGeneric(
  name = "allClusters",
  def = function(x) {
    standardGeneric("allClusters")
  }
)

setGeneric(
  name = "primaryCluster",
  def = function(x) {
    standardGeneric("primaryCluster")
  }
)

setGeneric(
  name = "primaryClusterIndex",
  def = function(x) {
    standardGeneric("primaryClusterIndex")
  }
)

setGeneric(
  name = "primaryClusterIndex<-",
  def = function(object, value) {
    standardGeneric("primaryClusterIndex<-")
  }
)

setGeneric(
  name = "coClustering",
  def = function(x) {
    standardGeneric("coClustering")
  }
)

setGeneric(
  name = "dendrogram",
  def = function(x) {
    standardGeneric("dendrogram")
  }
)

setGeneric(
  name = "clusterType",
  def = function(x) {
    standardGeneric("clusterType")
  }
)

setGeneric(
  name = "addClusters",
  def = function(x, y,...) {
    standardGeneric("addClusters")
  }
)
#' @rdname ClusterCells-class
setGeneric(
  name = "removeClusters",
  def = function(x, whichRemove,exactMatch=TRUE) {
    standardGeneric("removeClusters")
  }
)
setGeneric(
  name = "removeUnclustered",
  def = function(x) {
    standardGeneric("removeUnclustered")
  }
)

setGeneric(
  name = "clusterInfo",
  def = function(x) {
    standardGeneric("clusterInfo")
  }
)
