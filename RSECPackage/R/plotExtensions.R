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