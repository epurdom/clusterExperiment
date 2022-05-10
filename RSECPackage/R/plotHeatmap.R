#' @name plotHeatmap
#' @title Extend plotHeatmap to RSECClass
#
#' @rdname plotHeatmap
#' @aliases plotHeatmap,RSECClass-method
setMethod(
  f = "plotHeatmap",
  signature = signature(data = "RSECClass"),
  definition = function(data,
                        whichClusters= c("primary","workflow","all","none"),
                        ...
  ){
    whCl<-getClusterIndex(data,whichClusters=whichClusters,noMatch="silentlyRemove")
		callNextMethod(data,whichClusters=whCl)
		
})