#' @rdname clusterSingle
setGeneric(
  name = "clusterSingle",
  def = function(x,  ...) {
    standardGeneric("clusterSingle")
  }
)

#' @rdname plotHeatmap
setGeneric(
  name = "plotCoClustering",
  def = function(data,  ...) {
    standardGeneric("plotCoClustering")
  }
)

#' @rdname makeDendrogram
setGeneric(
  name = "makeDendrogram",
  def = function(x,  ...) {
    standardGeneric("makeDendrogram")
  }
)

#' @rdname ClusterExperiment-class
setGeneric(
    name = "clusterLegend",
    def = function(x) {
        standardGeneric("clusterLegend")
    }
)
#' @rdname ClusterExperiment-class
setGeneric(
    name = "clusterLegend<-",
    def = function(object, value) {
        standardGeneric("clusterLegend<-")
    }
)
#' @rdname ClusterExperiment-class
setGeneric(
    name = "orderSamples",
    def = function(x) {
        standardGeneric("orderSamples")
    }
)
#' @rdname ClusterExperiment-class
setGeneric(
    name = "orderSamples<-",
    def = function(object, value) {
        standardGeneric("orderSamples<-")
    }
)
#' @rdname ClusterExperiment-class
setGeneric(
  name = "clusterLabels",
  def = function(x, whichClusters) {
    standardGeneric("clusterLabels")
  }
)
#' @rdname ClusterExperiment-class
setGeneric(
    name = "clusterLabels<-",
    def = function(object, value) {
        standardGeneric("clusterLabels<-")
    }
)
#' @rdname ClusterExperiment-class
setGeneric(
  name = "nClusters",
  def = function(x) {
    standardGeneric("nClusters")
  }
)
#' @rdname ClusterExperiment-class
setGeneric(
  name = "nFeatures",
  def = function(x) {
    standardGeneric("nFeatures")
  }
)
#' @rdname ClusterExperiment-class
setGeneric(
  name = "nSamples",
  def = function(x) {
    standardGeneric("nSamples")
  }
)


#' @rdname ClusterExperiment-class
setGeneric(
  name = "pipelineClusters",
  def = function(x,iteration=0) {
    standardGeneric("pipelineClusters")
  }
)
#' @rdname ClusterExperiment-class
setGeneric(
  name = "pipelineClusterDetails",
  def = function(x) {
    standardGeneric("pipelineClusterDetails")
  }
)
#' @rdname ClusterExperiment-class
setGeneric(
  name = "pipelineClusterTable",
  def = function(x) {
    standardGeneric("pipelineClusterTable")
  }
)
#' @rdname ClusterExperiment-class
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
#' @rdname plotDendrogram
setGeneric(
  name="plotDendrogram",
  def=function(x,...)
    {
    standardGeneric("plotDendrogram")
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
  def=function(data,...)
  {
    standardGeneric("plotHeatmap")
  }
)
#' @rdname ClusterExperiment-class
setGeneric(
  name = "transformation",
  def = function(x) {
    standardGeneric("transformation")
  }
)
#' @rdname ClusterExperiment-class
  setGeneric(
  name = "transform",
  def = function(x,...) {
    standardGeneric("transform")
  }
)

#' @rdname ClusterExperiment-class
setGeneric(
  name = "clusterMatrix",
  def = function(x) {
    standardGeneric("clusterMatrix")
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
  name = "coClustering<-",
  def = function(object, value) {
    standardGeneric("coClustering<-")
  }
)

setGeneric(
  name = "clusterType",
  def = function(x) {
    standardGeneric("clusterType")
  }
)
setGeneric(
    name = "clusterType<-",
    def = function(object,value) {
        standardGeneric("clusterType<-")
    }
)
setGeneric(
  name = "addClusters",
  def = function(x, y,...) {
    standardGeneric("addClusters")
  }
)
#' @rdname ClusterExperiment-class
setGeneric(
  name = "removeClusters",
  def = function(x, whichRemove,exactMatch=TRUE) {
    standardGeneric("removeClusters")
  }
)
#' @rdname ClusterExperiment-class
setGeneric(
  name = "removeUnclustered",
  def = function(x) {
    standardGeneric("removeUnclustered")
  }
)

#' @rdname ClusterExperiment-class
setGeneric(
  name = "clusterInfo",
  def = function(x) {
    standardGeneric("clusterInfo")
  }
)

#' @rdname ClusterExperiment-class
setGeneric(
    name = "convertClusterColors",
    def = function(x) {
        standardGeneric("convertClusterColors")
    }
)

#' @rdname combineMany
setGeneric(
  name = "combineMany",
  def = function(x, whichClusters, ...) {
    standardGeneric("combineMany")
  }
)

#' @rdname clusterMatrixNamed
setGeneric(
  name = "clusterMatrixNamed",
  def = function(x,  ...) {
    standardGeneric("clusterMatrixNamed")
  }
)
#' @rdname primaryClusterNamed
setGeneric(
  name = "primaryClusterNamed",
  def = function(x,  ...) {
    standardGeneric("primaryClusterNamed")
  }
)
#' @rdname getBestFeatures
setGeneric(
  name = "getBestFeatures",
  def = function(x, ...) {
    standardGeneric("getBestFeatures")
  }
)
#' @rdname mergeClusters
setGeneric(
  name = "mergeClusters",
  def = function(x, ...) {
    standardGeneric("mergeClusters")
  }
)
