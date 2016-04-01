setGeneric(
  name = "clusterAll",
  def = function(x,  subsample=TRUE, sequential=FALSE,
                 clusterFunction=c("tight", "hierarchical", "pam", "kmeans"),
                 clusterDArgs=NULL, subsampleArgs=NULL, seqArgs=NULL,
                 isLog=TRUE) {
    standardGeneric("clusterAll")
  }
)

setGeneric(
  name = "isLog",
  def = function(x) {
    standardGeneric("isLog")
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
  def = function(x, y) {
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
