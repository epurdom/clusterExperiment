##########
###This code is now defunct. Only saving it in case want to resurrect it
##########

# #' Methods to set and return the clustering labels
# #'
# #' These methods will return or set the clustering labels (i.e., the
# #' colnames of the clusterMatrix slot).
# #'
# #' @aliases clusterLabels
# #'
# #' @param x,object a \code{\link{ClusterExperiment}} object.
# #' @param whichClusters controls which labels to be returned. It is either
# #'   numeric, in which case gives the indices of the clusters, or character, in
# #'   which case it matches to \code{clusterTypes(x)} to find the indices of the
# #'   clusters.
# #'
# #' @return A vector with the clustering labels.
# #' @export
# #' @rdname clusterLabels
# #' @examples
# #' data(simData)
# #'
# #' cl <- clusterMany(simData,nReduceDims=c(5,10,50),  reduceMethod="PCA",
# #' clusterFunction="pam", ks=2:4, findBestK=c(FALSE), removeSil=TRUE,
# #' subsample=FALSE)
# #'
# #' clusterLabels(cl)
# #' clusterLabels(cl, whichClusters=1:3)
# #'
# #' clusterTypes(cl)
# #' clusterLabels(cl, whichClusters="clusterMany")
# #' @rdname clusterLabels
# #' @export
# setMethod(
#   f = "clusterLabels",
#   signature = signature(x = "ClusterExperiment"),
#   definition = function(x){
#     labels<-colnames(clusterMatrix(x))
#     if(is.null(labels)) cat("No labels found for clusterings\n")
#     return(labels)
#     
#   }
# )
# 
# setMethod(
#   f = "clusterLabels",
#   signature = signature(x = "ClusterExperiment",whichClusters="numeric"),
#   definition = function(x, whichClusters){
#     if(!all(whichClusters %in% 1:NCOL(clusterMatrix(x)))) stop("Invalid indices for clusterLabels")
#     labels<-colnames(clusterMatrix(x))[whichClusters]
#     if(is.null(labels)) cat("No labels found for clusterings\n")
#     return(labels)
#   }
# )
# setMethod(
#   f = "clusterLabels",
#   signature = signature(x = "ClusterExperiment",whichClusters="missing"),
#   definition = function(x, whichClusters){
#     clusterLabels(x,whichClusters=1:nClusterings(x))
#   }
# )
# #' @export
# #' @rdname clusterLabels
# #' @aliases clusterLabels<-
# #' @param value character. A vector of length equal to
# #'   \code{NCOL(clusterMatrix(object))}.
# setReplaceMethod(
#   f = "clusterLabels",
#   signature = signature(object="ClusterExperiment", value="character"),
#   definition = function(object, value) {
#     if(length(value)!=NCOL(clusterMatrix(object))) stop("value must be a vector of length equal to NCOL(clusterMatrix(object)):",NCOL(clusterMatrix(object)))
#     if(any(duplicated(value))) stop("cannot have duplicated clusterLabels")
#     colnames(object@clusterMatrix) <- value
#     validObject(object)
#     return(object)
#   }
# )
# 
# #' @rdname clusterLabels
# #' @export
# setMethod(
#   f = "clusterLabels",
#   signature = signature(x = "ClusterExperiment", whichClusters ="character"),
#   definition = function(x, whichClusters="all"){
#     wh<-.TypeIntoIndices(x,whClusters=whichClusters)
#     return(clusterLabels(x,wh))
#   }
# )
# 
# #' @rdname clusterLabels
# #' @export
# setMethod(
#   f = "clusterLabels",
#   signature = signature(x = "ClusterExperiment",whichClusters="missing"),
#   definition = function(x, whichClusters){
#     return(clusterLabels(x,whichClusters="all"))
#   }
# )
