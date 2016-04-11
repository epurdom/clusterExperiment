#' @rdname clusterAll
setGeneric(
  name = "clusterAll",
  def = function(x,  ...) {
    standardGeneric("clusterAll")
  }
)


#' @rdname ClusterExperiments-class
setGeneric(
    name = "clusterColors",
    def = function(x) {
        standardGeneric("clusterColors")
    }
)
#' @rdname ClusterExperiments-class
setGeneric(
    name = "clusterColors<-",
    def = function(object, value) {
        standardGeneric("clusterColors<-")
    }
)
#' @rdname ClusterExperiments-class
setGeneric(
    name = "orderSamples",
    def = function(x) {
        standardGeneric("orderSamples")
    }
)
#' @rdname ClusterExperiments-class
setGeneric(
    name = "orderSamples<-",
    def = function(object, value) {
        standardGeneric("orderSamples<-")
    }
)
#' @rdname ClusterExperiments-class
setGeneric(
  name = "clusterLabels",
  def = function(x, whichClusters) {
    standardGeneric("clusterLabels")
  }
)
#' @rdname ClusterExperiments-class
setGeneric(
    name = "clusterLabels<-",
    def = function(object, value) {
        standardGeneric("clusterLabels<-")
    }
)
#' @rdname ClusterExperiments-class
setGeneric(
  name = "nClusters",
  def = function(x) {
    standardGeneric("nClusters")
  }
)
#' @rdname ClusterExperiments-class
setGeneric(
  name = "nFeatures",
  def = function(x) {
    standardGeneric("nFeatures")
  }
)
#' @rdname ClusterExperiments-class
setGeneric(
  name = "nSamples",
  def = function(x) {
    standardGeneric("nSamples")
  }
)


#' @rdname ClusterExperiments-class
setGeneric(
  name = "pipelineClusters",
  def = function(x,iteration=0) {
    standardGeneric("pipelineClusters")
  }
)
#' @rdname ClusterExperiments-class
setGeneric(
  name = "pipelineClusterDetails",
  def = function(x) {
    standardGeneric("pipelineClusterDetails")
  }
)
#' @rdname ClusterExperiments-class
setGeneric(
  name = "pipelineClusterTable",
  def = function(x) {
    standardGeneric("pipelineClusterTable")
  }
)
#' @rdname ClusterExperiments-class
setGeneric(
  name = "paramValues",
  def = function(x,type) {
    standardGeneric("paramValues")
  }
)

#' @rdname clusterMany
setGeneric(
  name = "clusterMany",
  def = clusterMany <- function(x, ... ) {
    standardGeneric("clusterMany")
  }
)

#' @rdname plotClusters
setGeneric(
  name="plotClusters",
  def=function(clusters, orderClusters,...)
    {
    standardGeneric("plotClusters")
  }
)
#' @rdname plotHeatmap
setGeneric(
  name="plotHeatmap",
  def=function(x,...)
  {
    standardGeneric("plotHeatmap")
  }
)
#' @rdname ClusterExperiments-class
setGeneric(
  name = "transformation",
  def = function(x) {
    standardGeneric("transformation")
  }
)
#' @rdname ClusterExperiments-class
  setGeneric(
  name = "transform",
  def = function(x,nPCADims=NA,nVarDims=NA,dimReduce="none") {
    standardGeneric("transform")
  }
)

#' @rdname ClusterExperiments-class
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
#' @rdname ClusterExperiments-class
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
