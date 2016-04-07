setGeneric(
  name = "clusterAll",
  def = function(x,  subsample=TRUE, sequential=FALSE,
                 clusterFunction=c("tight", "hierarchical", "pam", "kmeans"),
                 clusterDArgs=NULL, subsampleArgs=NULL, seqArgs=NULL,
                 isCount=FALSE,transFun=NULL,dimReduce=c("none","PCA","mostVar"),ndims=NA) {
    standardGeneric("clusterAll")
  }
)

setGeneric(
  name = "compareChoices",
  def = compareChoices <- function(x, ks=3:5, clusterMethod="pam", alphas=0.1, findBestK=FALSE,sequential=FALSE,
                                   removeSil=FALSE, subsample=FALSE,silCutoff=0,
                                   dimReduce="none",nVarDims=NA,nPCADims=NA,
                                   clusterDArgs=list(minSize=5),
                                   subsampleArgs=list(resamp.num=50),
                                   seqArgs=list(beta=0.9,k.min=3, verbose=FALSE),
                                   transFun=NULL,isCount=FALSE,
                                   ncores=1,random.seed=NULL,run=TRUE,paramMatrix=NULL,...
                                   ) {
    standardGeneric("compareChoices")
  }
)

setGeneric(
  name = "transformation",
  def = function(x) {
    standardGeneric("transformation")
  }
)
setGeneric(
  name = "transform",
  def = function(x,nPCADims=NA,nVarDims=NA,dimReduce="none") {
    standardGeneric("transform")
  }
)

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
