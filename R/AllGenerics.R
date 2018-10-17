setGeneric(".matchToStats", function(object,...) { standardGeneric(".matchToStats")})

setGeneric("filterNames", function(object,...) { standardGeneric("filterNames")})
setGeneric("filterData", function(object,...) { standardGeneric("filterData")})
setGeneric("defaultNDims",function(object,...){standardGeneric("defaultNDims")})
setGeneric(name = "makeFilterStats", function(object,...){ standardGeneric("makeFilterStats")})
setGeneric(name="makeReducedDims", function(object,...){ standardGeneric("makeReducedDims")})

setGeneric(name="plotClusterLegend", function(object,...){standardGeneric("plotClusterLegend")})
setGeneric("plotFeatureBoxplot", function(object,feature,...) { standardGeneric("plotFeatureBoxplot")})
setGeneric("plotFeatureScatter", function(object,features,...) { standardGeneric("plotFeatureScatter")})

setGeneric( "tableClusters", function(object,...){ standardGeneric("tableClusters") })
setGeneric( "plotClustersTable", function(object,...){ standardGeneric("plotClustersTable") })
setGeneric("bubblePlot", function(propTable,sizeTable,...) { standardGeneric("bubblePlot")})

setGeneric(name = "assignUnassigned", def=function(object,...){ standardGeneric("assignUnassigned") })
setGeneric( name = "removeUnassigned", def = function(object,...) {
	standardGeneric("removeUnassigned")})
setGeneric(name="subsetByCluster", def = function(x,...) {
	standardGeneric("subsetByCluster")})
setGeneric(name = "getReducedData", def=function(object,...){ standardGeneric("getReducedData") })

###Merge Info
setGeneric("getMergeCorrespond", def=function(x,...){ standardGeneric("getMergeCorrespond")})
setGeneric(name = "nodeMergeInfo", def=function(x,...){ standardGeneric("nodeMergeInfo") })
setGeneric(name = "mergeClusterIndex",def=function(x,...){ standardGeneric("mergeClusterIndex")})
setGeneric(name = "mergeMethod",def=function(x,...){standardGeneric("mergeMethod")})
setGeneric(name = "mergeCutoff", def=function(x,...){ standardGeneric("mergeCutoff")})
setGeneric(name="eraseMergeInfo",def=function(x,...){ standardGeneric("eraseMergeInfo")})
setGeneric(name = "nClusters", function(x,...){ standardGeneric("nClusters")})

setGeneric("plotClustersWorkflow", function(object,...){ standardGeneric("plotClustersWorkflow")}
)
setGeneric(
	name="plotContrastHeatmap",
	def=function(object,...){
		standardGeneric("plotContrastHeatmap")
	}
)
setGeneric(
    name = "RSEC",
    def = function(x, ...) {
        standardGeneric("RSEC")
    }
)
setGeneric(
  name = "subsampleClustering",
  def = function(clusterFunction, ...) {
    standardGeneric("subsampleClustering")
  }
)
setGeneric(
  name = "mainClustering",
  def = function(clusterFunction, ...) {
    standardGeneric("mainClustering")
  }
)
setGeneric(
  name = "seqCluster",
  def = function(clusterFunction, ...) {
    standardGeneric("seqCluster")
  }
)
setGeneric(
  name = "clusterSingle",
  def = function(x, diss, ...) {
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

setGeneric("clusterLegend", function(x) { standardGeneric("clusterLegend")})
setGeneric("clusterLegend<-", function(object, value){standardGeneric("clusterLegend<-") })
setGeneric("renameClusters", function(object,value,...) { standardGeneric("renameClusters")})
setGeneric("recolorClusters", function(object,value,...) { standardGeneric("recolorClusters")})
setGeneric("sampleDendrogram", function(x) { standardGeneric("sampleDendrogram")})
setGeneric("clusterDendrogram", function(x) { standardGeneric("clusterDendrogram")})
setGeneric("convertToDendrogram",function(x){standardGeneric("convertToDendrogram")})
setGeneric("nodeIds",function(x,type){standardGeneric("nodeIds")})
setGeneric("checkDendrogram",function(x,dendroCluster,dendroSample,...){standardGeneric("checkDendrogram")})
setGeneric("nInternalNodes",function(x){standardGeneric("nInternalNodes")})
# #Have to recreate these generics, because phylobase's are nonstandardGenericFunctions because used braces.
# setGeneric("nodeLabels",function(x){standardGeneric("nodeLabels")})
# setGeneric("nodeLabels<-",function(x,...,value){standardGeneric("nodeLabels<-")})
# setGeneric("tipLabels",function(x){standardGeneric("tipLabels")})
# setGeneric("nNodes",function(x){standardGeneric("nNodes")})
# setGeneric("nTips",function(x){standardGeneric("nTips")})

setGeneric("getClusterIndex",function(object,...){standardGeneric("getClusterIndex")})
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
	name="getClusterManyParams",
	def=function(x,...){
		standardGeneric("getClusterManyParams")
	})
setGeneric(
  name = "nClusterings",
  def = function(x) {
    standardGeneric("nClusterings")
  }
)

setGeneric(
  name = "nFeatures",
  def = function(x) {
    standardGeneric("nFeatures")
  }
)

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

setGeneric( "plotClusters", function(object, ...) { standardGeneric("plotClusters")})
setGeneric( name="plotBarplot", def=function(object, ...) { standardGeneric("plotBarplot") })
setGeneric("plotReducedDims",function(object, whichCluster,...){ standardGeneric("plotReducedDims")})

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
  name = "transformation<-",
  def = function(object, value) {
    standardGeneric("transformation<-")
  }
)
setGeneric(
  name = "transformData",
  def = function(object,...) {
    standardGeneric("transformData")
  }
)
setGeneric("colDataClusters",function(object,...){standardGeneric("colDataClusters")})
setGeneric("addToColData",function(object,...){standardGeneric("addToColData")})
setGeneric("clusterMatrix", function(x,...) { standardGeneric("clusterMatrix")})
setGeneric("clusterMatrixNamed", function(x,...) { standardGeneric("clusterMatrixNamed")})
setGeneric("clusterMatrixColors", function(x,...) { standardGeneric("clusterMatrixColors")})
setGeneric("primaryCluster", function(x) { standardGeneric("primaryCluster")})
setGeneric("primaryClusterIndex", function(x) {standardGeneric("primaryClusterIndex")})
setGeneric("primaryClusterLabel", function(x) {standardGeneric("primaryClusterLabel")})
setGeneric("primaryClusterType", function(x) {standardGeneric("primaryClusterType")})
setGeneric("primaryClusterNamed",function(x) { standardGeneric("primaryClusterNamed")})
setGeneric(name = "primaryClusterIndex<-", def = function(object, value) {standardGeneric("primaryClusterIndex<-")})
setGeneric(name = "dendroClusterIndex",def = function(x) {standardGeneric("dendroClusterIndex")})
setGeneric( name = "coClustering", def = function(x) { standardGeneric("coClustering")})
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


setGeneric( name = "removeClusters", def = function(x, whichCluster,...) {
    standardGeneric("removeClusters")})
setGeneric( name = "addClusterings", def = function(x, y,...) {
    standardGeneric("addClusterings")})
setGeneric(name = "removeClusterings", def = function(x, whichClusters,...) {
    standardGeneric("removeClusterings")})
setGeneric( name = "clusteringInfo", def = function(x) { standardGeneric("clusteringInfo")})
setGeneric( name = "makeConsensus", def = function(x, ...) {
    standardGeneric("makeConsensus")})




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
setGeneric(
  name = "getBuiltInFunction",
  def = function(object, ...) {
    standardGeneric("getBuiltInFunction")
  }
)

 
setGeneric(
  name = "requiredArgs",
  def = function(object, ...) {
    standardGeneric("requiredArgs")
  }
)
setGeneric(
  name = "algorithmType",
  def = function(object, ...) {
    standardGeneric("algorithmType")
  }
)
setGeneric(
  name = "inputType",
  def = function(object, ...) {
    standardGeneric("inputType")
  }
)
setGeneric(
  name = "getPostProcessingArgs",
  def = function(clusterFunction, ...) {
    standardGeneric("getPostProcessingArgs")
  }
)
