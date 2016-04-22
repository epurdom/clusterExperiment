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

#' @rdname convertClusterLegend
setGeneric(
  name = "convertClusterLegend",
  def = function(object,  ...) {
    standardGeneric("convertClusterLegend")
  }
)

#' @rdname makeDendrogram
setGeneric(
  name = "makeDendrogram",
  def = function(x,  ...) {
    standardGeneric("makeDendrogram")
  }
)

#' @rdname ClusterExperiment-methods
#' @name ClusterExperiment-methods
#' @aliases clusterLegend
setGeneric(
    name = "clusterLegend",
    def = function(x) {
        standardGeneric("clusterLegend")
    }
)

#' @rdname ClusterExperiment-methods
setGeneric(
    name = "clusterLegend<-",
    def = function(object, value) {
        standardGeneric("clusterLegend<-")
    }
)

#' @rdname ClusterExperiment-methods
setGeneric(
    name = "orderSamples",
    def = function(x) {
        standardGeneric("orderSamples")
    }
)

#' @rdname ClusterExperiment-methods
setGeneric(
    name = "orderSamples<-",
    def = function(object, value) {
        standardGeneric("orderSamples<-")
    }
)

#' @rdname clusterLabels
setGeneric(
  name = "clusterLabels",
  def = function(x, whichClusters) {
    standardGeneric("clusterLabels")
  }
)

#' @rdname clusterLabels
setGeneric(
    name = "clusterLabels<-",
    def = function(object, value) {
        standardGeneric("clusterLabels<-")
    }
)

#' @rdname ClusterExperiment-methods
setGeneric(
  name = "nClusters",
  def = function(x) {
    standardGeneric("nClusters")
  }
)

#' @rdname ClusterExperiment-methods
setGeneric(
  name = "nFeatures",
  def = function(x) {
    standardGeneric("nFeatures")
  }
)

#' @rdname ClusterExperiment-methods
setGeneric(
  name = "nSamples",
  def = function(x) {
    standardGeneric("nSamples")
  }
)

#' @rdname pipelineClusters
setGeneric(
  name = "pipelineClusters",
  def = function(x,iteration=0) {
    standardGeneric("pipelineClusters")
  }
)

#' @rdname pipelineClusters
setGeneric(
  name = "pipelineClusterDetails",
  def = function(x) {
    standardGeneric("pipelineClusterDetails")
  }
)

#' @rdname pipelineClusters
setGeneric(
  name = "pipelineClusterTable",
  def = function(x) {
    standardGeneric("pipelineClusterTable")
  }
)

# setGeneric(
#   name = "paramValues",
#   def = function(x,type) {
#     standardGeneric("paramValues")
#   }
# )

#' @rdname clusterMany
setGeneric(
  name = "clusterMany",
  def = clusterMany <- function(x, ... ) {
    standardGeneric("clusterMany")
  }
)

#' @rdname makeDendrogram
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
    def=function(clusters, whichClusters,...)
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

#' @rdname ClusterExperiment-methods
setGeneric(
  name = "transformation",
  def = function(x) {
    standardGeneric("transformation")
  }
)

#' @rdname transform
setGeneric(
  name = "transform",
  def = function(x,...) {
    standardGeneric("transform")
  }
)

#' @rdname ClusterExperiment-methods
setGeneric(
  name = "clusterMatrix",
  def = function(x) {
    standardGeneric("clusterMatrix")
  }
)

#' @rdname ClusterExperiment-methods
setGeneric(
  name = "primaryCluster",
  def = function(x) {
    standardGeneric("primaryCluster")
  }
)

#' @rdname ClusterExperiment-methods
setGeneric(
  name = "primaryClusterIndex",
  def = function(x) {
    standardGeneric("primaryClusterIndex")
  }
)

#' @rdname ClusterExperiment-methods
setGeneric(
  name = "primaryClusterIndex<-",
  def = function(object, value) {
    standardGeneric("primaryClusterIndex<-")
  }
)

#' @rdname ClusterExperiment-methods
setGeneric(
  name = "coClustering",
  def = function(x) {
    standardGeneric("coClustering")
  }
)

#' @rdname ClusterExperiment-methods
setGeneric(
  name = "coClustering<-",
  def = function(object, value) {
    standardGeneric("coClustering<-")
  }
)

#' @rdname ClusterExperiment-methods
setGeneric(
  name = "clusterType",
  def = function(x) {
    standardGeneric("clusterType")
  }
)

#' @rdname ClusterExperiment-methods
setGeneric(
    name = "clusterType<-",
    def = function(object,value) {
        standardGeneric("clusterType<-")
    }
)

#' @rdname addClusters
setGeneric(
  name = "addClusters",
  def = function(x, y,...) {
    standardGeneric("addClusters")
  }
)

#' @rdname addClusters
setGeneric(
  name = "removeClusters",
  def = function(x, whichRemove,exactMatch=TRUE) {
    standardGeneric("removeClusters")
  }
)

#' @rdname addClusters
setGeneric(
  name = "removeUnclustered",
  def = function(x) {
    standardGeneric("removeUnclustered")
  }
)

#' @rdname ClusterExperiment-methods
setGeneric(
  name = "clusterInfo",
  def = function(x) {
    standardGeneric("clusterInfo")
  }
)

# setGeneric(
#     name = "convertClusterColors",
#     def = function(x) {
#         standardGeneric("convertClusterColors")
#     }
# )

#' @rdname combineMany
setGeneric(
  name = "combineMany",
  def = function(x, whichClusters, ...) {
    standardGeneric("combineMany")
  }
)

#' @rdname ClusterExperiment-methods
setGeneric(
  name = "clusterMatrixNamed",
  def = function(x) {
    standardGeneric("clusterMatrixNamed")
  }
)

#' @rdname ClusterExperiment-methods
setGeneric(
  name = "primaryClusterNamed",
  def = function(x) {
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
