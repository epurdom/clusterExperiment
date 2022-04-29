#' @name ClusterExperiment-methods
#' @title Helper methods for the ClusterExperiment class
#'
#' @description This is a collection of helper methods for the ClusterExperiment
#'   class.
#' @param ... For \code{addToColData}, arguments passed to
#'   \code{colDataClusters}.
#' @param value The value to be substituted in the corresponding slot. See the
#'   slot descriptions in \code{\link{ClusterExperiment}} for details on what
#'   objects may be passed to these functions.
#' @rdname ClusterExperiment-methods
#' @aliases show show,ClusterExperiment-method
#' @examples
#' # load data:
#' data(rsecFluidigm)
#' show(rsecFluidigm)
#' #Number of clusterings
#' nClusterings(rsecFluidigm)
#' # Number of clusters per clustering
#' nClusters(rsecFluidigm)
#' # Number of features/samples
#' nSamples(rsecFluidigm)
#' nFeatures(rsecFluidigm)
#' # retrieve all clustering assignments
#' # (either as cluster ids, cluster names or cluster colors)
#' head(clusterMatrix(rsecFluidigm)[,1:5])
#' head(clusterMatrixNamed(rsecFluidigm)[,1:5])
#' head(clusterMatrixColors(rsecFluidigm)[,1:5])
#' # clustering Types/Labels
#' clusterTypes(rsecFluidigm)
#' clusterLabels(rsecFluidigm)
#' # Add a clustering assignment to the colData of the object
#' # (useful if working with function that relies on colData)
#' colData(rsecFluidigm)
#' test<-addToColData(rsecFluidigm,whichCluster="primary")
#' colData(test)
#' @export
setMethod(
  f = "show",
  signature = "ClusterExperiment",
  definition = function(object) {
		callNextMethod(object)
    cat("reducedDimNames:",if(anyValidReducedDims(object)) reducedDimNames(object) else "no reduced dims stored","\n")
    cat("filterStats:",if(anyValidFilterStats(object)) filterNames(object) else "no valid filtering stats stored","\n")
    cat("Workflow progress:\n")
    typeTab<-names(table(clusterTypes(object)))
    cat("clusterMany run?",if("clusterMany" %in% typeTab) "Yes" else "No","\n")
    cat("makeConsensus run?",if("makeConsensus" %in% typeTab) "Yes" else "No","\n")
    cat("makeDendrogram run?",if(!is.null(object@dendro_samples) & !is.null(object@dendro_clusters) ) "Yes" else "No","\n")
    cat("mergeClusters run?",if("mergeClusters" %in% typeTab) "Yes" else "No","\n")
  }
)


#' @rdname ClusterExperiment-methods
#' @return \code{coClustering} returns/sets the co-clustering matrix.
#' @export
#' @aliases coClustering
setMethod(
  f = "coClustering",
  signature = "ClusterExperiment",
  definition = function(x) {
    return(x@coClustering)
  }
)
#' @rdname ClusterExperiment-methods
#' @export
#' @aliases coClustering<-
setReplaceMethod(
  f = "coClustering",
  signature = signature(object="ClusterExperiment", value="matrix"),
  definition = function(object, value) {
    #This appears to detect if matrix is symmetric without any further effort
    object@coClustering <- Matrix::Matrix(value,sparse=TRUE)        
    ch<-.checkCoClustering(object)
    if(is.logical(ch) && ch) return(object) else stop(ch)
  }
)
#' @rdname ClusterExperiment-methods
#' @export
#' @aliases coClustering<-
setReplaceMethod(
  f = "coClustering",
  signature = signature(object="ClusterExperiment", value="dsCMatrix"),
  definition = function(object, value) {
    object@coClustering <- value
    ch<-.checkCoClustering(object)
    if(is.logical(ch) && ch) return(object) else stop(ch)
  }
)

#' @rdname ClusterExperiment-methods
#' @export
#' @aliases coClustering<-
setReplaceMethod(
  f = "coClustering",
  signature = signature(object="ClusterExperiment", value="numeric"),
  definition = function(object, value) {
    object@coClustering <- value
    ch<-.checkCoClustering(object)
    if(is.logical(ch) && ch) return(object) else stop(ch)
  }
)

