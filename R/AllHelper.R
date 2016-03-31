## subsetting
#' @rdname clusterCells
setMethod(
  f = "[",
  signature = c("ClusterCells", "ANY", "ANY"),
  definition = function(x, i, j, ..., drop=TRUE) {
    out <- callNextMethod()
    out@clusterLabels <- as.matrix(x@clusterLabels[j,])
    out@coClustering <- new("matrix")
    out@dendrogram <- list()
    return(out)
  }
)

## show
#' @rdname clusterCells
setMethod(
  f = "show",
  signature = "ClusterCells",
  definition = function(object) {
    cat("class:", class(object), "\n")
    cat("dim:", dim(object), "\n")
    cat("primary clustering:")
    print(table(object@clusterLabels[,object@primaryIndex]))
    cat("cluster type:", object@clusterType[object@primaryIndex])
  }
)

#' @rdname clusterCells
setMethod(
  f = "isLog",
  signature = "ClusterCells",
  definition = function(x) {
    return(x@isLog)
  }
)

#' @rdname clusterCells
setMethod(
  f = "allClusters",
  signature = "ClusterCells",
  definition = function(x) {
    return(x@clusterLabels)
  }
)

#' @rdname clusterCells
setMethod(
  f = "primaryCluster",
  signature = "ClusterCells",
  definition = function(x) {
    return(x@clusterLabels[,x@primaryIndex])
  }
)

#' @rdname clusterCells
setMethod(
  f = "primaryClusterIndex",
  signature = "ClusterCells",
  definition = function(x) {
    return(x@primaryIndex)
  }
)

#' @rdname clusterCells
setReplaceMethod(
  f = "primaryClusterIndex",
  signature = signature("ClusterCells", "numeric"),
  definition = function(object, value) {
    object@primaryIndex <- value
    validObject(object)
    return(object)
  }
)

#' @rdname clusterCells
setMethod(
  f = "coClustering",
  signature = "ClusterCells",
  definition = function(x) {
    return(x@coClustering)
  }
)

#' @rdname clusterCells
setMethod(
  f = "dendrogram",
  signature = "ClusterCells",
  definition = function(x) {
    return(x@dendrogram)
  }
)

#' @rdname clusterCells
setMethod(
  f = "clusterType",
  signature = "ClusterCells",
  definition = function(x) {
    return(x@clusterType)
  }
)

#' @rdname clusterCells
setMethod(
  f = "addClusters",
  signature = signature("ClusterCells", "matrix"),
  definition = function(x, y) {
    if(!(NROW(y) == NCOL(x))) {
      stop("Incompatible dimensions.")
    }
    x@clusterLabels <- cbind(x@clusterLabels, y)
    x@clusterType <- c(x@clusterType, rep("User", NCOL(y)))
    return(x)
  }
)

#' @rdname clusterCells
setMethod(
  f = "addClusters",
  signature = signature("ClusterCells", "numeric"),
  definition = function(x, y) {
    if(!(length(y) == NCOL(x))) {
      stop("Incompatible dimensions.")
    }
    x@clusterLabels <- cbind(x@clusterLabels, y)
    x@clusterType <- c(x@clusterType, "User")
    return(x)
  }
)

#' @rdname clusterCells
setMethod(
  f = "addClusters",
  signature = signature("ClusterCells", "ClusterCells"),
  definition = function(x, y) {
    if(!all(assay(y) == assay(x))) {
      stop("Cannot merge clusters from different data.")
    }
    x@clusterLabels <- cbind(x@clusterLabels, y@clusterLabels)
    x@clusterType <- c(x@clusterType, y@clusterType)
    return(x)
  }
)

#' @rdname clusterCells
setMethod(
  f = "removeUnclustered",
  signature = "ClusterCells",
  definition = function(x) {
    return(x[,primaryCluster(x) > 0])
  }
)
