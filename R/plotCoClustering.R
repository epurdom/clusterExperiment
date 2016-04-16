#' @rdname plotHeatmap
setMethod(
  f = "plotCoClustering",
  signature = "ClusterExperiment",
  definition = function(data, ...){
    if(nrow(data@coClustering)==0) stop("coClustering slot is empty")
    fakeCE<-clusterExperiment(data@coClustering,
                              clusterMatrix(data),
                              transformation=transformation(data),
                              clusterInfo=clusterInfo(data),
                              clusterType=clusterType(data)
                              )
    clusterLegend(fakeCE)<-clusterLegend(data)
    orderSamples(fakeCE)<-orderSamples(data)
    colData(fakeCE)<-colData(data)
    fakeCE@dendro_samples<-data@dendro_samples
    fakeCE@primaryIndex<-data@primaryIndex
    validObject(fakeCE) #just in case screwed something up
    plotHeatmap(fakeCE,...)
  })