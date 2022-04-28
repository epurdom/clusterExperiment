## Additional Object
setGeneric( name = "coClustering", def = function(x) { standardGeneric("coClustering")})
setGeneric(name = "coClustering<-",def = function(object, value) {standardGeneric("coClustering<-")})



###Merge Info
setGeneric("getMergeCorrespond", def=function(x,...){ standardGeneric("getMergeCorrespond")})
setGeneric(name = "nodeMergeInfo", def=function(x,...){ standardGeneric("nodeMergeInfo") })
setGeneric(name = "mergeClusterIndex",def=function(x,...){ standardGeneric("mergeClusterIndex")})
setGeneric(name = "mergeMethod",def=function(x,...){standardGeneric("mergeMethod")})
setGeneric(name = "mergeCutoff", def=function(x,...){ standardGeneric("mergeCutoff")})
setGeneric(name="eraseMergeInfo",def=function(x,...){ standardGeneric("eraseMergeInfo")})
setGeneric(name = "nClusters", function(x,...){ standardGeneric("nClusters")})

# Clustering/RSEC
setGeneric(name = "RSEC",def = function(x, ...) {standardGeneric("RSEC")})
setGeneric(name = "subsampleClustering",def = function(clusterFunction, ...) {standardGeneric("subsampleClustering")})
setGeneric(name = "mainClustering",def = function(clusterFunction, ...) {standardGeneric("mainClustering")})
setGeneric(name = "seqCluster",def = function(clusterFunction, ...) {standardGeneric("seqCluster")})
setGeneric( name = "clusterSingle",def = function(inputMatrix, ...) { standardGeneric("clusterSingle") })
setGeneric(name="getClusterManyParams",def=function(x,...){standardGeneric("getClusterManyParams")})
setGeneric(name = "clusterMany",def = clusterMany <- function(x, ... ) {standardGeneric("clusterMany")})

# Merge/Coclustering
setGeneric(name = "plotCoClustering",def = function(data,  ...) {standardGeneric("plotCoClustering")})
setGeneric( name = "makeConsensus", def = function(x, ...) {
    standardGeneric("makeConsensus")})
setGeneric(name = "mergeClusters",def = function(x, ...) {standardGeneric("mergeClusters")})


# ClusterFunction
setGeneric(name = "getBuiltInFunction",def = function(object, ...) {standardGeneric("getBuiltInFunction")})
setGeneric(name = "requiredArgs",def = function(object, ...) {standardGeneric("requiredArgs")})
setGeneric(name = "algorithmType",def = function(object, ...) {standardGeneric("algorithmType")})
setGeneric(name = "inputType",def = function(object, ...) {standardGeneric("inputType")})
setGeneric(name = "getPostProcessingArgs",def = function(clusterFunction, ...) {standardGeneric("getPostProcessingArgs")})
