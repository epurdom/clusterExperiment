#' @name addClusterings
#' @title Add clusterings to ClusterExperiment object
#' @description Function for adding new clusterings in form of vector (single cluster) or matrix (multiple clusterings) to an existing ClusterExperiment object
#' @aliases addClusterings removeClusterings addClusterings,ClusterExperiment,matrix-method
#' @return A \code{ClusterExperiment} object.
#' @export
#' @examples
#' data(simData)
#'
#' cl1 <- clusterSingle(simData, subsample=FALSE,
#' sequential=FALSE, mainClusterArgs=list(clusterArgs=list(k=3), clusterFunction="pam"))

#' cl2 <- clusterSingle(simData, subsample=FALSE,
#' sequential=FALSE, mainClusterArgs=list(clusterArgs=list(k=3), clusterFunction="pam"))
#'
#' addClusterings(cl1, cl2)
setMethod(
  f = "addClusterings",
  signature = signature("ClusterExperiment", "matrix"),
  definition = function(x, y, clusterTypes="User",clusterLabels=NULL,clusterLegend=NULL) {
    if(!is.null(clusterLabels)){
      if(length(clusterLabels)!=ncol(y)) stop("clusterLabels must vector of length equal to the number of clusterings (columns of y)")
      colnames(y)<-clusterLabels
    }
    ccObj<-ClusterExperiment(x,
                             clusters=y,
                             transformation=transformation(x),
                             clusterTypes=clusterTypes,
                             checkTransformAndAssay=FALSE,
                             clusterLegend=clusterLegend)
    addClusterings(x,ccObj)
  }
)

#' @rdname addClusterings
#' @export
setMethod(
  f = "addClusterings",
  signature = signature("ClusterExperiment", "ClusterExperiment"),
  definition = function(x, y) {
    if(!all(dim(assay(y)) == dim(assay(x))) || !all(assay(y) == assay(x))) {
      stop("Cannot merge clusters from different data.")
    }
    x@clusterMatrix <- cbind(x@clusterMatrix, y@clusterMatrix)
    x@clusterTypes <- c(x@clusterTypes, y@clusterTypes)
    x@clusterInfo<-c(x@clusterInfo,y@clusterInfo)
    x@clusterLegend<-c(x@clusterLegend,y@clusterLegend)
    if(any(duplicated(colnames(x@clusterMatrix)))){
      colnames(x@clusterMatrix)<-make.names(colnames(x@clusterMatrix),unique=TRUE)
    }
    x<-.unnameClusterSlots(x) #just gets rid of the names of objects that shouldn't have them
    ch<-.checkClusterMatrix(x)
    if(!is.logical(ch)) stop(ch)
    ch<-.checkClusterTypes(x)
    if(!is.logical(ch)) stop(ch)
    ch<-.checkClusterLegend(x)
    if(!is.logical(ch)) stop(ch)
    #would it be less memory to do a call to "new"? What is difference versus having to check dendrogram, coClustering, etc if they exists? Should do checks on large data.
    return(x)
  }
)

#' @rdname addClusterings
#' @export
#' @param ... For \code{addClusterings}, passed to signature 
#' \code{ClusterExperiment,matrix}. For \code{[} (subsetting), passed to 
#' \code{SingleCellExperiment} subsetting function.
#' @param makePrimary whether to make the added cluster the primary cluster (only relevant if \code{y} is a vector)
setMethod(
  f = "addClusterings",
  signature = signature("ClusterExperiment", "vector"),
  definition = function(x, y, makePrimary=FALSE,...) {
    mat<-matrix(y,ncol=1)
    x<-addClusterings(x,mat,...)
    if(makePrimary){
      x@primaryIndex<-ncol(clusterMatrix(x))
    }
    return(x)
  }
)
