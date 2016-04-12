
#' @title Class ClusterExperiment
#'
#' @description \code{ClusterExperiment} is a class that extends
#' \code{SummarizedExperiment} and is used to store the single-cell RNA-seq data
#' and clustering information.
#'
#' @docType class
#' @aliases ClusterExperiment ClusterExperiment-class clusterExperiment
#'
#' @description In addition to the slots of the \code{SummarizedExperiment}
#' class, the \code{ClusterExperiment} object has the following additional
#' slots:
#' @slot transformation function. Function to transform the data by when methods
#' that assume normal-like data (e.g. log)
#' @slot clusterMatrix matrix. A matrix giving the integer-valued cluster ids
#' for each sample. The rows of the matrix correspond to clusterings and columns
#' to samples. The integer values are assigned in the order that the clusters
#' were found, if found by setting sequential=TRUE in clusterAll. "-1" indicates
#' the sample was not clustered.
#' @slot primaryIndex: numeric. An index that specifies the primary set of
#' labels.
#' @slot clusterInfo: list. A list with info about the clustering.
#' If created from \code{\link{clusterAll}}, clusterInfo will include the
#' parameter used for the call, and the call itself. If \code{sequential = TRUE}
#' it will also include the following components.
#' \itemize{
#' \item{\code{clusterInfo}}{if sequential=TRUE and clusters were successfully
#' found, a matrix of information regarding the algorithm behavior for each
#' cluster (the starting and stopping K for each cluster, and the number of
#' iterations for each cluster).}
#' \item{\code{whyStop}}{if sequential=TRUE and clusters were successfully
#' found, a character string explaining what triggered the algorithm to stop.}
#' }
#' @slot clusterType character vector with the origin of each column of
#' clusterMatrix.
#' @slot dendrogram dendrogram. A dendrogram containing the cluster
#' relationship.
#' @slot coClustering matrix. A matrix with the cluster co-occurrence
#' information; this can either be based on subsampling or on co-clustering
#' across parameter sets (see \code{clusterMany}). The matrix is a square matrix
#' with number of rows/columns equal to the number of samples.
#' @slot clusterColors a list, one per cluster in \code{clusterMatrix}. Each
#' element of the list is a matrix with nrows equal to the number of different
#' clusters in the clustering, and consisting of at least two columns with the
#' following column names: "clusterId" and "color".
#' @slot orderSamples a numeric vector (of integers) defining the order of
#' samples to be used for plotting of samples. Usually set internally by other
#' functions.
#' @name ClusterExperiment-class
#' @rdname ClusterExperiment-class
#' @exportClass
#'
setClass(
  Class = "ClusterExperiment",
  contains = "SummarizedExperiment",
  slots = list(
    transformation="function",
              clusterMatrix = "matrix",
               primaryIndex = "numeric",
               clusterInfo = "list",
               clusterType = "character",
               dendrogram = "list",
               coClustering = "matrix",
              clusterColors="list",
              orderSamples="numeric"
    )
)

## One question is how to extend the "[" method, i.e., how do we subset the co-occurance matrix and the dendrogram?
## For now, if subsetting, these are lost, but perhaps we can do something smarter?

setValidity("ClusterExperiment", function(object) {
  if(length(assays(object)) != 1) {
    return("There must be one assay slot.")
  }
  if(!is.numeric(assay(object))) {
    return("The data must be numeric.")
  }
  if(any(is.na(assay(object)))) {
    return("NA values are not allowed.")
  }
  tX <- try(transform(object),silent=TRUE)
  if(inherits(tX, "try-error")){
    stop(paste("User-supplied `transformation` produces error on the input data
               matrix:\n",x))
  }
  if(any(is.na(tX))) {
    return("NA values after transforming data matrix are not allowed.")
  }

  if(!all(is.na((object@clusterMatrix))) &
     !(NROW(object@clusterMatrix) == NCOL(object))) {
    return("If present, `clusterMatrix` must have as many row as cells.")
  }
  if(!is.numeric(object@clusterMatrix)) {
    return("`clusterMatrix` must be a numeric matrix.")
  }

  if(NCOL(object@clusterMatrix)!= length(object@clusterType)) {
    return("length of clusterType must be same as NCOL of the clusterMatrix")
  }

  if(NCOL(object@clusterMatrix)!= length(object@clusterInfo)) {
    return("length of clusterInfo must be same as NCOL of the clusterMatrix")
  }

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
  if(!all(is.na(object@clusterMatrix))){ #what does this mean, how can they be all NA?
      if(length(object@primaryIndex) != 1) {
        return("If more than one set of cluster labels, a primary cluster must
               be specified.")
      }
      if(object@primaryIndex > NCOL(object@clusterMatrix) |
         object@primaryIndex < 1) {
        return("`primaryIndex` out of bounds.")
      }
      if(NCOL(object@clusterMatrix) != length(object@clusterType)) {
        return("`clusterType` must be the same length as NCOL of
               `clusterMatrix`.")
      }
      #test that @clusterColors is proper form
      if(length(object@clusterColors) != NCOL(object@clusterMatrix)) {
        return("`clusterColors` must be list of same length as NCOL of
               `clusterMatrix`")
      }
      testIsMatrix <- sapply(object@clusterColors,
                             function(x) {!is.null(dim(x))})
      if(!all(testIsMatrix)) {
        return("Each element of `clusterColors` list must be a matrix")
      }
      testColorRows <- sapply(object@clusterColors, function(x){nrow(x)})
      testClusterMat <- apply(object@clusterMatrix, 2, function(x) {
        length(unique(x))})
      if(!all(testColorRows == testClusterMat)) {
        return("each element of `clusterColors` must be matrix with number of
               rows equal to the number of clusters (including -1 or -2 values)
               in `clusterMatrix`")
      }
      testColorCols1 <- sapply(object@clusterColors, function(x) {
        "color" %in% colnames(x)})
      testColorCols2<-sapply(object@clusterColors,function(x){"clusterIds" %in% colnames(x)})
    if(!all(testColorCols1) || !all(testColorCols2)) return("each element of 'clusterColors' must be matrix with at least 2 columns, and at least 2 columns have names 'clusterIds' and 'color'")
    testColorCols1<-sapply(object@clusterColors,function(x){is.character(x)})
    if(!all(testColorCols1)) return("each element of 'clusterColors' must be matrix of character values")
    testColorCols1<-sapply(1:length(object@clusterColors),function(ii){
        col<-object@clusterColors[[ii]]
        x<-object@clusterMatrix[,ii]
        y<-as.numeric(col[,"clusterIds"])
        all(y %in% x)
    })
    if(!all(testColorCols1)) return("each element of 'clusterColors' must be matrix with column 'clusterIds' matching the corresponding integer valued clusterMatrix values")

  }
  if(length(object@orderSamples)!=NCOL(assay(object))) return("'orderSamples' must be of same length as number of samples (NCOL(assay(object)))")
  if(any(!object@orderSamples %in% 1:NCOL(assay(object)))) return("'orderSamples' must be values between 1 and the number of samples.")
  return(TRUE)
})

#' @description The constructor \code{clusterExperiment} creates an object of
#' the class \code{ClusterExperiment}. However, the typical way of creating
#' these objects is the result of a call to \code{clusterMany} or
#' \code{clusterAll}.
#'
#' @description Note that when subsetting the data, the co-clustering and
#' dendrogram information are lost.
#'
#'@param se a matrix or \code{SummarizedExperiment} containing the clustered
#'data.
#'@param labels a vector with cluster labels.
#'@param transformation function. A function to transform the data before
#'performing steps that assume normal-like (i.e. constant variance), such as
#'the log
#'
#'@return A \code{ClusterExperiment} object.
#'
#'@examples
#'
#'se <- matrix(data=rnorm(200), ncol=10)
#'labels <- gl(5, 2)
#'
#'cc <- clusterExperiment(se, as.numeric(labels),
#'                         transformation = function(x){x})
#'
#' @rdname ClusterExperiment-class
setGeneric(
  name = "clusterExperiment",
  def = function(se,  clusters,...) {
    standardGeneric("clusterExperiment")
  }
)
#' @rdname ClusterExperiment-class
setMethod(
  f = "clusterExperiment",
  signature = signature("matrix","ANY"),
  definition = function(se, clusters, ...){
    clusterExperiment(SummarizedExperiment(se), clusters, ...)
  })
#' @rdname ClusterExperiment-class
setMethod(
  f = "clusterExperiment",
  signature = signature("SummarizedExperiment", "numeric"),
  definition = function(se, clusters, ...){
    if(NCOL(se) != length(clusters)) {
      stop("`clusters` must be a vector of length equal to the number of
           samples.")
    }
  clusterExperiment(se,matrix(data=clusters, ncol=1),...)
})
#' @rdname ClusterExperiment-class
setMethod(
  f = "clusterExperiment",
  signature = signature("SummarizedExperiment","character"),
  definition = function(se, clusters,...){
    clusters <- as.numeric(factor(clusters))
    warning("The character vector `clusters` was coerced to integer values (one
            per cluster)")
    clusterExperiment(se,clusters,...)
    })
#' @rdname ClusterExperiment-class
setMethod(
  f = "clusterExperiment",
  signature = signature("SummarizedExperiment","factor"),
  definition = function(se, clusters,...){
    clusters <- as.numeric(clusters)
    warning("The factor `clusters` was coerced to numeric.")
    clusterExperiment(se,clusters,...)
  })
#' @rdname ClusterExperiment-class
setMethod(
  f = "clusterExperiment",
  signature = signature("SummarizedExperiment","matrix"),
  definition = function(se, clusters, transformation,clusterType="User",
                        clusterInfo=NULL){
    if(NCOL(se) != nrow(clusters)) {
      stop("`clusters` must be a matrix of rows equal to the number of
           samples.")
    }
    if(length(clusterType)==1) {
      clusterType <- rep(clusterType, length=NCOL(clusters))
    }
    if(is.null(clusterInfo)) {
      clusterInfo<-rep(list(NULL),length=NCOL(clusters))
    }
    if(length(clusterType)!=NCOL(clusters)) {
      stop("clusterType must be of length equal to number of clusters in
           `clusters`")
    }
    if(length(clusterType) == 1) {
        clusterType <- rep(clusterType, length=NCOL(clusters))
    }
    if(is.null(clusterInfo)) {
        clusterInfo <- rep(list(NULL), length=NCOL(clusters))
    }
    out <- new("ClusterExperiment",
               assays = Assays(assays(se)),
               elementMetadata = mcols(se),
               colData = colData(se),
               transformation=transformation,
               clusterMatrix = clusters,
               primaryIndex = 1,
               clusterType = clusterType,
               clusterInfo=clusterInfo,
               clusterColors=.makeColors(clusters, colors=bigPalette)$colorList,
               orderSamples=1:ncol(se)
    )
    validObject(out)
    return(out)
  })
