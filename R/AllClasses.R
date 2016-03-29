#' @rdname clusterCells
#' @export
setClass(
  Class = "ClusterCells",
  contains = "SummarizedExperiment",
  slots = list(isLog = "logical", labels = "matrix", dendrogram = "list",
                        coClustering = "matrix")
)

## this class extends SummarizedExperiment (I'm on an old version of R now, should be changed from SummarizedExperiment0)
## the idea is to store the data into the assay slot, having a primary cluster label as a column of colData
## in principle, we could just use the additional slots as element of the metadata list, but we loose the ability to check for the right class.
## One question is how to extend the "[" method, i.e., how do we subset the co-occurance matrix and the dendrogram?
## For now, if subsetting, these are lost, but perhaps we can do something smarter?
## Q: do we want to enforce the cluster labels to be numeric? Or factor?

setValidity("ClusterCells", function(object) {
  if(length(assays(object)) != 1) {
    return("There must be one assay slot.")
  }
  if(!is.numeric(assay(object))) {
    return("The data must be numeric.")
  }
  if(any(is.na(assay(object)))) {
    return("NA values are not allowed.")
  }
  if(!("clusterLabels" %in% names(colData(object)))) {
    return("colData must contain a column named `clusterLabels`.")
  }
  if(NROW(object@labels) > 0 & !(NROW(object@labels) == NCOL(object))) {
    return("If present, `labels` must have as many row as cells.")
  }
  if(!is.numeric(object@labels)) {
    return("`labels` must be a numeric matrix.")
  }
  if(length(object@dendrogram) > 0) {
    if(class(object@dendrogram) != "dendrogram") {
      return("`dendrogram` must be of class dendrogram.")
    }
  }
  if(NROW(object@coClustering)>0 &
       (NROW(object@coClustering) != NCOL(object@coClustering)
       | NCOL(object@coClustering) != NCOL(object))) {
    return("`coClustering` must be a sample by sample matrix.")
  }

  return(TRUE)
})

#' ClusterCells object and constructor
#'
#' \code{ClusterCells} is a class that extends \code{SummarizedExperiment},
#' used to store the single-cell RNA-seq data and clustering information.
#' In addition to the slots of the \code{SummarizedExperiment} class, the
#' \code{ClusterCells} object has the following additional slots:
#' \itemize{
#' \item isLog: logical. Whether the data are in the linear or log scale.
#' \item labels: matrix. A matrix of cluster labels, useful for consensus
#' clustering.
#' \item dendrogram: dendrogram. A dendrogram containing the cluster
#' relationship.
#' \item coClustering: matrix. A matrix with the cluster co-occurrence
#' information; this can either be based on subsampling or on co-clustering
#' across parameter sets (see \code{compareChoices}).
#' }
#'
#' The constructor \code{clusterCells} creates an object of the class
#' \code{ClusterCells}. However, the typical way of creating these objects is
#' the result of a call to \code{compareChoices} or \code{clusterAll}.
#'
#' Note that when subsetting the data, the co-clustering and dendrogram
#' information are lost.
#'
#'@param se a matrix or \code{SummarizedExperiment} containing the clustered
#'data.
#'@param labels a vector with cluster labels.
#'@param isLog logical. Whether the data are in log (TRUE) or linear (FALSE)
#'scale
#'
#'@return A \code{ClusterCells} object.
#'
#'@aliases ClusterCells ClusterCells-class
#'
#'@docType class
#'
#'@examples
#'
#'se <- matrix(data=rnorm(200), ncol=10)
#'labels <- gl(5, 2)
#'
#'cc <- clusterCells(se, as.numeric(labels), isLog = TRUE)
#'
clusterCells <- function(se, labels, isLog) {
  if(NCOL(se) != length(labels)) {
    stop("`labels` must be a vector of length equal to the number of samples.")
  }
  if(is(se, "SummarizedExperiment")) {
    cd <- colData(se)
    cd$clusterLabels <- labels
    out <- new("ClusterCells",
               assays = Assays(assays(se)),
               elementMetadata = mcols(se),
               colData = cd,
               isLog = isLog,
               labels = matrix(data=labels, ncol=1))
  } else if(is.matrix(se)) {
    out <- new("ClusterCells",
               assays = Assays(se),
               elementMetadata = new("DataFrame", nrows=nrow(se)),
               isLog = isLog,
               labels = matrix(data=labels, ncol=1),
               colData = DataFrame(clusterLabels=labels)
               )
  } else {
    stop("`se` must be a matrix or SummarizedExperiment object.")
  }

  return(out)
}

## subsetting
setMethod("[", c("ClusterCells", "ANY", "ANY"),
          function(x, i, j, ..., drop=TRUE) {
            out <- callNextMethod()
            out@labels <- as.matrix(out@labels[j,])
            out@coClustering <- new("matrix")
            out@dendrogram <- list()
            return(out)
          }
)

