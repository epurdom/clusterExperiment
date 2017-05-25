setGeneric(
    name = "RSEC",
    def = function(x, ...) {
        standardGeneric("RSEC")
    }
)

setGeneric(
  name = "clusterSingle",
  def = function(x, diss,  ...) {
    standardGeneric("clusterSingle")
  }
)
setGeneric(
  name = "setToCurrent",
  def = function(x,  ...) {
    standardGeneric("setToCurrent")
  }
)
setGeneric(
  name = "setToFinal",
  def = function(x,  ...) {
    standardGeneric("setToFinal")
  }
)
setGeneric(
  name = "plotCoClustering",
  def = function(data,  ...) {
    standardGeneric("plotCoClustering")
  }
)
setGeneric(
    name="clusterContrasts",
    def=function(cluster,...){
        standardGeneric("clusterContrasts")
    }
)
setGeneric(
  name = "convertClusterLegend",
  def = function(object,  ...) {
    standardGeneric("convertClusterLegend")
  }
)

setGeneric(
  name = "makeDendrogram",
  def = function(x,  ...) {
    standardGeneric("makeDendrogram")
  }
)

setGeneric(
    name = "clusterLegend",
    def = function(x) {
        standardGeneric("clusterLegend")
    }
)

setGeneric(
    name = "clusterLegend<-",
    def = function(object, value) {
        standardGeneric("clusterLegend<-")
    }
)

setGeneric(
    name = "orderSamples",
    def = function(x) {
        standardGeneric("orderSamples")
    }
)

setGeneric(
    name = "orderSamples<-",
    def = function(object, value) {
        standardGeneric("orderSamples<-")
    }
)

setGeneric(
  name = "clusterLabels",
  def = function(x) {
    standardGeneric("clusterLabels")
  }
)

setGeneric(
    name = "clusterLabels<-",
    def = function(object, value) {
        standardGeneric("clusterLabels<-")
    }
)

setGeneric(
  name = "nClusters",
  def = function(x) {
    standardGeneric("nClusters")
  }
)

#' Generic function that returns the number of features
#'
#' Given an object that describes a dataset or a model, it returns the number of
#' features.
#' @param x an object that describes a dataset or a model.
#' @return the number of features.
#' @export
setGeneric(
  name = "nFeatures",
  def = function(x) {
    standardGeneric("nFeatures")
  }
)

#' Generic function that returns the number of samples
#'
#' Given an object that describes a model or a dataset, it returns the number of
#' samples.
#' @param x an object that describes a dataset or a model.
#' @return the number of samples.
#' @export
setGeneric(
  name = "nSamples",
  def = function(x) {
    standardGeneric("nSamples")
  }
)

setGeneric(
  name = "workflowClusters",
  def = function(x,iteration=0) {
    standardGeneric("workflowClusters")
  }
)

setGeneric(
  name = "workflowClusterDetails",
  def = function(x) {
    standardGeneric("workflowClusterDetails")
  }
)

setGeneric(
  name = "workflowClusterTable",
  def = function(x) {
    standardGeneric("workflowClusterTable")
  }
)

setGeneric(
  name = "clusterMany",
  def = clusterMany <- function(x, ... ) {
    standardGeneric("clusterMany")
  }
)

setGeneric(
  name="plotDendrogram",
  def=function(x,...)
    {
    standardGeneric("plotDendrogram")
  }
)

setGeneric(
    name="plotClusters",
    def=function(clusters, whichClusters,...)
    {
        standardGeneric("plotClusters")
    }
)

setGeneric(
    name="plotBarplot",
    def=function(clusters, whichClusters,...)
    {
        standardGeneric("plotBarplot")
    }
)

setGeneric(
  name="plotHeatmap",
  def=function(data,...)
  {
    standardGeneric("plotHeatmap")
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
  def = function(x,...) {
    standardGeneric("transform")
  }
)

setGeneric(
  name = "clusterMatrix",
  def = function(x,whichClusters) {
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
  name = "clusterTypes",
  def = function(x) {
    standardGeneric("clusterTypes")
  }
)

setGeneric(
    name = "clusterTypes<-",
    def = function(object,value) {
        standardGeneric("clusterTypes<-")
    }
)

setGeneric(
  name = "addClusters",
  def = function(x, y,...) {
    standardGeneric("addClusters")
  }
)

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

setGeneric(
  name = "combineMany",
  def = function(x, whichClusters, ...) {
    standardGeneric("combineMany")
  }
)

setGeneric(
  name = "clusterMatrixNamed",
  def = function(x) {
    standardGeneric("clusterMatrixNamed")
  }
)

setGeneric(
  name = "primaryClusterNamed",
  def = function(x) {
    standardGeneric("primaryClusterNamed")
  }
)

setGeneric(
  name = "getBestFeatures",
  def = function(x, ...) {
    standardGeneric("getBestFeatures")
  }
)

setGeneric(
  name = "mergeClusters",
  def = function(x, ...) {
    standardGeneric("mergeClusters")
  }
)
