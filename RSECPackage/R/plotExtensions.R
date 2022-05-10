#' @details If \code{whichClusters = "workflow"}, then the workflow clusterings
#'   will be plotted in the following order: final, mergeClusters, makeConsensus,
#'   clusterMany.
#' @export
#' @rdname plotClusters
#' @aliases plotClusters,RSECClass-method

setMethod(
  f = "plotClusters",
  signature = signature(object = "RSECClass"),
  definition = function(object, whichClusters="workflow", 
  	existingColors=c("ignore","all","firstOnly"),
    resetNames=FALSE, resetColors=FALSE, resetOrderSamples=FALSE, colData=NULL, clusterLabels=NULL,...)
  {
    whichClusters<-getClusterIndex(object,
      whichClusters=whichClusters,noMatch="throwError")
    callNextMethod(object,whichClusters=whichClusters)
	}
)

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