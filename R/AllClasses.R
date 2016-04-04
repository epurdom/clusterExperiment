#' @rdname clusterCells
#' @export
setClass(
  Class = "ClusterCells",
  contains = "SummarizedExperiment",
  slots = list(
    transformation="function",
              clusterLabels = "matrix",
               primaryIndex = "numeric",
               clusterInfo = "list",
               clusterType = "character",
               dendrogram = "list",
               coClustering = "matrix")
)

## One question is how to extend the "[" method, i.e., how do we subset the co-occurance matrix and the dendrogram?
## For now, if subsetting, these are lost, but perhaps we can do something smarter?

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
  tX<-try(transform(object))
  if(inherits(tX, "try-error")){
    stop(paste("User-supplied `transformation` produces error on the input data matrix:\n",x))
  }
  if(any(is.na(tX))) {
    return("NA values after transforming data matrix are not allowed.")
  }
  
  if(!all(is.na((object@clusterLabels))) &
     !(NROW(object@clusterLabels) == NCOL(object))) {
    return("If present, `clusterLabels` must have as many row as cells.")
  }
  if(!is.numeric(object@clusterLabels)) {
    return("`clusterLabels` must be a numeric matrix.")
  }
  
  if(NCOL(object@clusterLabels)!= length(object@clusterType)) return("length of clusterType must be same as NCOL of the clusterLabels")
  
  if(NCOL(object@clusterLabels)!= length(object@clusterInfo)) return("length of clusterInfo must be same as NCOL of the clusterLabels")

  if(length(object@dendrogram) > 0) {
    if(class(object@dendrogram) != "dendrogram") {
      return("`dendrogram` must be of class dendrogram.")
    }
  }
  if(!all(is.na(object@coClustering)>0) &
       (NROW(object@coClustering) != NCOL(object@coClustering)
       | NCOL(object@coClustering) != NCOL(object))) {
    return("`coClustering` must be a sample by sample matrix.")
  }
  if(!all(is.na(object@clusterLabels)) & length(object@primaryIndex) != 1) {
    return("If more than one set of cluster labels, a primary cluster must
           be specified.")
  }
  if(!all(is.na(object@clusterLabels)) &
     (object@primaryIndex > NCOL(object@clusterLabels) |
      object@primaryIndex < 1)) {
    return("`primaryIndex` out of bounds.")
  }
  if(!all(is.na(object@clusterLabels)) &
     NCOL(object@clusterLabels) != length(object@clusterType)) {
    return("`clusterType` must be the same length of `clusterLabels`.")
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
#' \item transformation: function. Function to transform the data by when methods that assume normal-like data (e.g. log)
#' \item clusterLabels: matrix. A matrix giving the
#' integer-valued cluster ids for each sample. The rows of the matrix correspond to clusterings and columns to samples. 
#' The integer values are assigned
#' in the order that the clusters were found, if found by setting sequential=TRUE in clusterAll. "-1" indicates
#' the sample was not clustered.
#' \item primaryIndex: numeric. An index that specifies the primary set of
#' labels.
#' \item clusterInfo: list. A list with info about the clustering.
#' \itam clusterType: character vector with the origin of each column of
#' clusterLabels.
#' \item dendrogram: dendrogram. A dendrogram containing the cluster
#' relationship.
#' \item coClustering: matrix. A matrix with the cluster co-occurrence
#' information; this can either be based on subsampling or on co-clustering
#' across parameter sets (see \code{compareChoices}). The matrix is a square matrix with number of rows/columns equal to the number of samples.
#' }
#'
#' If created from \code{\link{clusterAll}}, clusterInfo will include the
#' parameter used for the call, and the call itself. If \code{sequential = TRUE}
#' it will also include the following components.
#'
#' \itemize{
#' \item{\code{clusterInfo}}{if sequential=TRUE and clusters were successfully
#' found, a matrix of information regarding the algorithm behavior for each
#' cluster (the starting and stopping K for each cluster, and the number of
#' iterations for each cluster).}
#' \item{\code{whyStop}}{if sequential=TRUE and clusters were successfully
#' found, a character string explaining what triggered the algorithm to stop.}
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
#'@param transformation function. A function to transform the data before performing steps that assume normal-like (i.e. constant variance), such as the log
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
#'cc <- clusterCells(se, as.numeric(labels), transformation = function(x){x})
#'
clusterCells <- function(se, labels, transformation,clusterType="User",clusterInfo=list(NULL)) {
  if(NCOL(se) != length(labels)) {
    stop("`labels` must be a vector of length equal to the number of samples.")
  }
  if(!is(se, "SummarizedExperiment") & !is.matrix(se)) {
    stop("`se` must be a matrix or SummarizedExperiment object.")
  }

  if(is.factor(labels)) {
    labels <- as.numeric(labels)
    warning("The factor `labels` was coerced to numeric.")
  }

  if(is.matrix(se)) {
    se <- SummarizedExperiment(se)
  }

  out <- new("ClusterCells",
             assays = Assays(assays(se)),
             elementMetadata = mcols(se),
             colData = colData(se),
             transformation=transformation,
             clusterLabels = matrix(data=labels, ncol=1),
             primaryIndex = 1,
             clusterType = clusterType,
             clusterInfo=clusterInfo
             )

  return(out)
}
