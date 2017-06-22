setOldClass("dendrogram")
setClassUnion("dendrogramOrNULL",members=c("dendrogram", "NULL"))
setClassUnion("matrixOrNULL",members=c("matrix", "NULL"))
setClassUnion("matrixOrMissing",members=c("matrix", "missing"))
setClassUnion("functionOrNULL",members=c("function", "NULL"))
#' @title Class ClusterExperiment
#'
#' @description \code{ClusterExperiment} is a class that extends
#' \code{SummarizedExperiment} and is used to store the data
#' and clustering information.
#'
#' @docType class
#' @aliases ClusterExperiment ClusterExperiment-class clusterExperiment
#'
#' @description In addition to the slots of the \code{SummarizedExperiment}
#' class, the \code{ClusterExperiment} object has the additional slots described
#' in the Slots section.
#'
#' @description There are several methods implemented for this class. The most
#' important methods (e.g., \code{\link{clusterMany}}, \code{\link{combineMany}},
#' ...) have their own help page. Simple helper methods are described in the
#' Methods section. For a comprehensive list of methods specific to this class
#' see the Reference Manual.
#'
#' @slot transformation function. Function to transform the data by when methods
#' that assume normal-like data (e.g. log)
#' @slot clusterMatrix matrix. A matrix giving the integer-valued cluster ids
#' for each sample. The rows of the matrix correspond to clusterings and columns
#' to samples. The integer values are assigned in the order that the clusters
#' were found, if found by setting sequential=TRUE in clusterSingle. "-1" indicates
#' the sample was not clustered.
#' @slot primaryIndex numeric. An index that specifies the primary set of
#' labels.
#' @slot clusterInfo list. A list with info about the clustering.
#' If created from \code{\link{clusterSingle}}, clusterInfo will include the
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
#' @slot clusterTypes character vector with the origin of each column of
#' clusterMatrix.
#' @slot dendro_samples dendrogram. A dendrogram containing the cluster
#' relationship (leaves are samples; see \code{\link{makeDendrogram}} for
#' details).
#' @slot dendro_clusters dendrogram. A dendrogram containing the cluster
#' relationship (leaves are clusters; see \code{\link{makeDendrogram}} for
#' details).
#' @slot dendro_index numeric. An integer giving the cluster that was used to
#'   make the dendrograms. NA_real_ value if no dendrograms are saved.
#' @slot dendro_outbranch logical. Whether the dendro_samples dendrogram put 
#' missing/non-clustered samples in an outbranch, or intermixed in the dendrogram.
#' @slot coClustering matrix. A matrix with the cluster co-occurrence
#' information; this can either be based on subsampling or on co-clustering
#' across parameter sets (see \code{clusterMany}). The matrix is a square matrix
#' with number of rows/columns equal to the number of samples.
#' @slot clusterLegend a list, one per cluster in \code{clusterMatrix}. Each
#' element of the list is a matrix with nrows equal to the number of different
#' clusters in the clustering, and consisting of at least two columns with the
#' following column names: "clusterId" and "color".
#' @slot orderSamples a numeric vector (of integers) defining the order of
#' samples to be used for plotting of samples. Usually set internally by other
#' functions.
#'
#' @name ClusterExperiment-class
#' @aliases ClusterExperiment
#' @rdname ClusterExperiment-class
#' @import SummarizedExperiment
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom dendextend as.phylo.dendrogram
#' @export
#'
setClass(
  Class = "ClusterExperiment",
  contains = "SummarizedExperiment",
  slots = list(
    transformation="function",
    clusterMatrix = "matrix",
    primaryIndex = "numeric",
    clusterInfo = "list",
    clusterTypes = "character",
    dendro_samples = "dendrogramOrNULL",
    dendro_clusters = "dendrogramOrNULL",
    dendro_index = "numeric",
	dendro_outbranch = "logical",
    coClustering = "matrixOrNULL",
    clusterLegend="list",
    orderSamples="numeric"
    )
)

## One question is how to extend the "[" method, i.e., how do we subset the co-occurance matrix and the dendrogram?
## For now, if subsetting, these are lost, but perhaps we can do something smarter?

setValidity("ClusterExperiment", function(object) {
  #browser()
  if(length(assays(object)) < 1) {
    return("There must be at least one assay slot.")
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

  if(NCOL(object@clusterMatrix)!= length(object@clusterTypes)) {
    return("length of clusterTypes must be same as NCOL of the clusterMatrix")
  }

  if(NCOL(object@clusterMatrix)!= length(object@clusterInfo)) {
    return("length of clusterInfo must be same as NCOL of the clusterMatrix")
  }
  ############
  ##Check dendrogram slotNames
  ############
  #browser()
  if(!is.null(object@dendro_samples)){
    if(nobs(object@dendro_samples) != NCOL(object)) {
      return("dendro_samples must have the same number of leaves as the number of samples")
    }
	if(is.na(object@dendro_outbranch)) return("if dendro_samples is defined, must also define dendro_outbranch")
  }
  else{
    if(!is.null(object@dendro_clusters)) return("dendro_samples should not be null if dendro_clusters is non-null")
	if(!is.na(object@dendro_outbranch)) return("dendro_samples should not be null if dendro_outbranch is not NA")
  }
  if(!is.null(object@dendro_clusters)){
    if(is.na(dendroClusterIndex(object))) return("if dendrogram slots are filled, must have corresponding dendro_index defined.")
    dcluster<-clusterMatrix(object)[,dendroClusterIndex(object)]
    if(nobs(object@dendro_clusters) != max(dcluster)) {
      return("dendro_clusters must have the same number of leaves as the number of (non-negative) clusters")
    }
  }
  else{
    if(!is.null(object@dendro_samples)) return("dendro_clusters should not be null if dendro_samples is non-null")
  }
  ## Check co-clustering
  if(!is.null(object@coClustering) &&
     (NROW(object@coClustering) != NCOL(object@coClustering)
      | NCOL(object@coClustering) != NCOL(object))) {
    return("`coClustering` must be a sample by sample matrix.")
  }
  ## If have a cluster matrix
  if(!all(is.na(object@clusterMatrix))){ #what does this mean, how can they be all NA?
    #check primary index
    if(length(object@primaryIndex) != 1) {
        if(length(object@primaryIndex) == 0) return("If more than one set of clusterings, a primary cluster must
               be specified.")
        if(length(object@primaryIndex) > 0) return("Only a single primary index may be specified")
    }
    if(object@primaryIndex > NCOL(object@clusterMatrix) |
       object@primaryIndex < 1) {
      return("`primaryIndex` out of bounds.")
    }
    #check clusterTypes
    if(NCOL(object@clusterMatrix) != length(object@clusterTypes)) {
      return("`clusterTypes` must be the same length as NCOL of
               `clusterMatrix`.")
    }
    #check internally stored as integers
    testConsecIntegers<-apply(object@clusterMatrix,2,function(x){
      whCl<-which(!x %in% c(-1,-2))
      uniqVals<-unique(x[whCl])
      return(all(sort(uniqVals)==1:length(uniqVals)))
    })
    #browser()
    if(!all(testConsecIntegers)) return("the cluster ids in clusterMatrix must be stored internally as consecutive integer values")

    ####
    #test that colnames of clusterMatrix appropriately aligns with everything else
    ####
    if(is.null(colnames(object@clusterMatrix))) return("clusterMatrix must have column names")
    if(any(duplicated(colnames(object@clusterMatrix)))) return("clusterMatrix must have unique column names")
    if(!is.null(names(object@clusterTypes))) return("clusterTypes should not have names")
    if(!is.null(names(object@clusterInfo))) return("clusterInfo should not have names")
    if(!is.null(names(object@clusterLegend))) return("clusterLegend should not have names")
    ####
    #test that @clusterLegend is proper form
    ####
    if(length(object@clusterLegend) != NCOL(object@clusterMatrix)) {
      return("`clusterLegend` must be list of same length as NCOL of
               `clusterMatrix`")
    }
    testIsMatrix <- sapply(object@clusterLegend,
                           function(x) {!is.null(dim(x))})
    if(!all(testIsMatrix)) {
      return("Each element of `clusterLegend` list must be a matrix")
    }
    testColorRows <- sapply(object@clusterLegend, function(x){nrow(x)})
    testClusterMat <- apply(object@clusterMatrix, 2, function(x) {
      length(unique(x))})
    if(!all(testColorRows == testClusterMat)) {
      return("each element of `clusterLegend` must be matrix with number of
               rows equal to the number of clusters (including -1 or -2 values)
               in `clusterMatrix`")
    }
    testColorCols1 <- sapply(object@clusterLegend, function(x) {
      "color" %in% colnames(x)})
    testColorCols2 <- sapply(object@clusterLegend, function(x) {
      "clusterIds" %in% colnames(x)})
    testColorCols3 <- sapply(object@clusterLegend, function(x) {
      "name" %in% colnames(x)})
    if(!all(testColorCols1) || !all(testColorCols2) || !all(testColorCols3)) {
      return("each element of `clusterLegend` must be matrix with at least 3
             columns, and at least 3 columns have names `clusterIds`,
             `color` and `name`")
    }
#     testUniqueName <- sapply(object@clusterLegend, function(x) {
#       any(duplicated(x[,"name"]))})
#     if(any(testUniqueName)) return("the column")
    testColorCols1 <- sapply(object@clusterLegend, function(x){is.character(x)})
    if(!all(testColorCols1)) {
      return("each element of `clusterLegend` must be matrix of character
             values")
    }
    testColorCols1 <- sapply(1:length(object@clusterLegend), function(ii){
      col<-object@clusterLegend[[ii]]
      x<-object@clusterMatrix[,ii]
      y<-as.numeric(col[,"clusterIds"])
      all(y %in% x)
    })
    if(!all(testColorCols1)) {
      return("each element of `clusterLegend` must be matrix with column
             `clusterIds` matching the corresponding integer valued
             clusterMatrix values")
    }
  }
  if(length(object@orderSamples)!=NCOL(assay(object))) {
    return("`orderSamples` must be of same length as number of samples
         (NCOL(assay(object)))")
  }
  if(any(!object@orderSamples %in% 1:NCOL(assay(object)))) {
    return("`orderSamples` must be values between 1 and the number of samples.")
  }
  return(TRUE)
})

#' @description The constructor \code{clusterExperiment} creates an object of
#' the class \code{ClusterExperiment}. However, the typical way of creating
#' these objects is the result of a call to \code{\link{clusterMany}} or
#' \code{\link{clusterSingle}}.
#'
#' @description Note that when subsetting the data, the co-clustering and
#' dendrogram information are lost.
#'
#'@param se a matrix or \code{SummarizedExperiment} containing the data to be
#'clustered.
#'@param clusters can be either a numeric or character vector, a factor, or a
#'numeric matrix, containing the cluster labels.
#'@param transformation function. A function to transform the data before
#'performing steps that assume normal-like data (i.e. constant variance), such
#'as the log.
#'@param ... The arguments \code{transformation}, \code{clusterTypes} and
#'  \code{clusterInfo} to be passed to the constructor for signature
#'  \code{SummarizedExperiment,matrix}.
#'
#'@return A \code{ClusterExperiment} object.
#'
#'@examples
#'
#'se <- matrix(data=rnorm(200), ncol=10)
#'labels <- gl(5, 2)
#'
#'cc <- clusterExperiment(se, as.numeric(labels), transformation =
#'function(x){x})
#'
#' @rdname ClusterExperiment-class
#' @export
setGeneric(
  name = "clusterExperiment",
  def = function(se,  clusters,...) {
    standardGeneric("clusterExperiment")
  }
)
#' @rdname ClusterExperiment-class
#' @export
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
      stop("`clusters` must be a vector of length equal to the number of samples.")
    }
  clusterExperiment(se,matrix(clusters, ncol=1),...)
})
#' @rdname ClusterExperiment-class
setMethod(
  f = "clusterExperiment",
  signature = signature("SummarizedExperiment","character"),
  definition = function(se, clusters,...){
    clusterExperiment(se,matrix(clusters,ncol=1),...)
    })
#' @rdname ClusterExperiment-class
setMethod(
  f = "clusterExperiment",
  signature = signature("SummarizedExperiment","factor"),
  definition = function(se, clusters,...){
    clusters <- as.character(clusters)
    clusterExperiment(se,clusters,...)
  })
#'@rdname ClusterExperiment-class
#'@param clusterTypes a string describing the nature of the clustering. The
#'  values `clusterSingle`, `clusterMany`, `mergeClusters`, `combineMany` are
#'  reserved for the clustering coming from the package workflow and should not
#'  be used when creating a new object with the constructor.
#'@param clusterInfo a list with information on the clustering (see Slots).
#'@param primaryIndex integer. Sets the `primaryIndex` slot (see Slots).
#'@param orderSamples a vector of integers. Sets the `orderSamples` slot (see
#'  Slots).
#'@param dendro_samples dendrogram. Sets the `dendro_samples` slot (see Slots).
#'@param dendro_clusters dendrogram. Sets the `dendro_clusters` slot (see
#'  Slots).
#'@param dendro_index numeric. Sets the dendro_index slot (see Slots).
#' @param dendro_outbranch logical. Sets the dendro_outbranch slot (see Slots).
#'@param coClustering matrix. Sets the `coClustering` slot (see Slots).
#'@details The \code{clusterExperiment} constructor function gives clusterLabels
#'  based on the column names of the input matrix/SummarizedExperiment. If
#'  missing, will assign labels "cluster1","cluster2", etc.
setMethod(
  f = "clusterExperiment",
  signature = signature("SummarizedExperiment","matrix"),
  definition = function(se, clusters,
            transformation,
            primaryIndex=1,
            clusterTypes="User",
            clusterInfo=NULL,
            orderSamples=1:ncol(se),
            dendro_samples=NULL,
            dendro_index=NA_real_,
            dendro_clusters=NULL,
			dendro_outbranch=NA,
            coClustering=NULL
            ){
    if(NCOL(se) != nrow(clusters)) {
      stop("`clusters` must be a matrix of rows equal to the number of
           samples.")
    }
    if(length(clusterTypes)==1) {
      clusterTypes <- rep(clusterTypes, length=NCOL(clusters))
    }
    if(is.null(clusterInfo)) {
      clusterInfo<-rep(list(NULL),length=NCOL(clusters))
    }
    if(length(clusterTypes)!=NCOL(clusters)) {
      stop("clusterTypes must be of length equal to number of clusters in
           `clusters`")
    }
    #fix up names of clusters and match
    #browser()
    if(is.null(colnames(clusters))){
      colnames(clusters)<-paste("cluster",1:NCOL(clusters),sep="")
    }
    if(any(duplicated(colnames(clusters)))){#probably not possible
      colnames(clusters)<-make.names(colnames(clusters),unique=TRUE)
    }
    if(length(clusterTypes) == 1) {
        clusterTypes <- rep(clusterTypes, length=NCOL(clusters))
    }
    if(is.null(clusterInfo)) {
        clusterInfo <- rep(list(NULL), length=NCOL(clusters))
    }
    #make clusters consecutive integer valued:
    tmp<-.makeColors(clusters, colors=bigPalette)
    clusterLegend<-tmp$colorList
    clustersNum<-tmp$numClusters
    colnames(clustersNum)<-colnames(clusters)
    #can just give se in constructor, and then don't loose any information!
    out <- new("ClusterExperiment",
               se,
               transformation=transformation,
               clusterMatrix = clustersNum,
               primaryIndex = primaryIndex,
               clusterTypes = unname(clusterTypes),
               clusterInfo=unname(clusterInfo),
               clusterLegend=unname(clusterLegend),
               orderSamples=1:ncol(se),
               dendro_samples=dendro_samples,
               dendro_clusters=dendro_clusters,
               dendro_index=dendro_index,
               dendro_outbranch=dendro_outbranch,
               coClustering=coClustering
    )
    validObject(out)
    return(out)
  })


################ clusterFunction class

#' @title Class ClusterFunction
#'   
#' @description \code{ClusterFunction} is a class for holding functions that can
#'   be used for clustering in the clustering algorithms in this package.
#'   
#' @docType class
#' @aliases ClusterFunction ClusterFunction-class clusterFunction
#' @slot clusterFUN a function defining the clustering function. See details for
#'   required arguments.
#' @slot inputType a character defining what type of input \code{clusterFUN}
#'   takes. Must be one of either "diss","X", or "either"
#' @slot algorithmType a character defining what type of clustering algorithm
#'   \code{clusterFUN} is. Must be one of either "01" or "K". \code{clusterFUN}
#'   must take the corresponding required arguments (see details below).
#' @slot classifyFUN a function that takes as input new data and the output of
#'   \code{clusterFUN} (when \code{cluster.only=FALSE} and results in cluster
#'   assignments of the new data.  Note that the function should assume that the
#'   input 'x' is not the same samples that were input to the clusterFunction
#'   (but can assume that it is the same number of features/columns). Used in
#'   subsampling clustering. If given value \code{NULL} then subsampling can
#'   only be \code{"InSample"}, see \code{\link{subsampleClustering}}.
#' @slot inputClassifyType the input type for the classification function (if
#'   not NULL); like \code{inputType}, must be one of "diss","X", or "either"
#' @slot outputType the type of output given by \code{clusterFUN}. Must either
#'   be "vector" or "list". If "vector" then the output should be a vector of
#'   length equal to the number of observations   with integer-valued elements
#'   identifying them to different clusters; the vector assignments should be in
#'   the same order as the original input of the data. Samples that are not
#'   assigned to any cluster should be given a '-1' value.  If "list", then it
#'   must be a list equal to the length of the number of clusters, and the
#'   elements of the list contain the indices of the samples in that cluster.
#'   Any indices not in any of the list elements are assumed to be -1. The main
#'   advantage of "list" is that it can preserve the order of the clusters if
#'   the \code{clusterFUN} desires to do so. In which case the \code{orderBy}
#'   argument of \code{\link{mainClustering}} can preserve this ordering (default is
#'   to order by size).
#' @slot requiredArgs Any additional required arguments for \code{clusterFUN}
#'   (beyond those required of all \code{clusterFUN}, described in details).
#' @slot checkFunctions logical. If TRUE, the validity check of the
#'   \code{ClusterFunction} object will check the \code{clusterFUN} with simple
#'   toy data using the function \code{internalFunctionCheck}.
#' @details Required arguments for \code{clusterFUN}: \itemize{ \item{"x or
#'   diss"}{either \code{x} and/or \code{diss} depending on \code{inputType}. If
#'   \code{x}, then \code{x} is assumed to be nfeatures x nsamples (like
#'   assay(CEObj) would give)} \item{"checkArgs"}{logical argument. If
#'   \code{checkArgs=TRUE}, the \code{clusterFUN} should check if the arguments
#'   passed in \code{...} are valid and return an error if not; otherwise, no
#'   error will be given, but the check should be done and only valid arguments
#'   in \code{...} passed along. This is necessary for the function to work with
#'   \code{clusterMany} which passes all arguments to all functions without
#'   checking. } \item{"cluster.only"}{logical argument. If
#'   \code{cluster.only=TRUE}, then \code{clusterFUN} should return only the
#'   vector of cluster assignments (or list if \code{outputType="list"}). If
#'   \code{cluster.only=FALSE} then the \code{clusterFUN} should return a named
#'   list where one of the elements entitled \code{clustering} contains the
#'   vector described above (no list!); anything else needed by the
#'   \code{classifyFUN} to classify new data should be contained in the output
#'   list as well. \code{cluster.only} is set internally depending on whether
#'   \code{classifyFUN} will be used by subsampling or only for clustering the
#'   final product.} \item{"..."}{Any additional arguments specific to the
#'   algorithm used by \code{clusterFUN} should be passed via \code{...} and NOT
#'   passed via arguments to \code{clusterFUN}} \item{"Other required
#'   arguments"}{\code{clusterFUN} must also accept arguments required for its
#'   \code{algorithmType} (see Details below).} }
#'   
#'   
#' @details \code{algorithmType}: Type "01" is for clustering functions that 
#'   expect as an input a dissimilarity matrix that takes on 0-1 values (e.g.
#'   from subclustering) with 1 indicating more dissimilarity between samples.
#'   "01" algorithm types must also have \code{inputType} equal to
#'   \code{"diss"}. It is also generally expected that "01" algorithms use the
#'   0-1 nature of the input to set criteria as to where to find clusters. "01"
#'   functions must take as an argument \code{alpha} between 0 and 1 to
#'   determine the clusters, where larger values of \code{alpha} require less
#'   similarity between samples in the same cluster. "K" is for clustering
#'   functions that require an argument \code{k} (the number of clusters), but
#'   arbitrary \code{inputType}.  On the other hand, "K" algorithms are assumed
#'   to need a predetermined 'k' and are also assumed to cluster all samples to
#'   a cluster. If not, the post-processing steps in \code{\link{mainClustering}} such
#'   as \code{findBestK} and \code{removeSil} may not operate correctly since
#'   they rely on silhouette distances.
#' @name ClusterFunction-class
#' @aliases ClusterFunction
#' @rdname ClusterFunction-class
#' @export
#'
setClass(
	Class = "ClusterFunction",
	slots = list(
  	  	clusterFUN="function",
  		inputType = "character",
  		algorithmType = "character",
  		inputClassifyType = "character",
		classifyFUN="functionOrNULL",
		outputType = "character",
		requiredArgs= "character",
		checkFunctions="logical"
  	)
)
.inputTypes<-c("X","diss","either")
.algTypes<-c("01","K")
.required01Args<-c("alpha")
.requiredKArgs<-c("k")
.outputTypes<-c("vector","list")

#' @rdname ClusterFunction-class
#' @export
#' @aliases internalFunctionCheck
#' @examples
#' #Use internalFunctionCheck to check possible function
#' goodFUN<-function(x,diss,k,checkArgs,cluster.only,...){cluster::pam(x=t(x),k=k,cluster.only=cluster.only)}
#' #passes internal check
#' internalFunctionCheck(goodFUN,inputType="X",algorithmType="K",outputType="vector")
#' #Note it doesn't pass if inputType="either" because no catches for x=NULL
#' internalFunctionCheck(goodFUN, inputType="either",algorithmType="K",outputType="vector")
#' myCF<-clusterFunction(clusterFUN=goodFUN, inputType="X",algorithmType="K", outputType="vector")
#' badFUN<-function(x,diss,k,checkArgs,cluster.only,...){cluster::pam(x=x,k=k)}
#' internalFunctionCheck(badFUN,inputType="X",algorithmType="K",outputType="vector")
#' @details \code{internalFunctionCheck} is the function that is called by the 
#'   validity check of the \code{clusterFunction} constructor (if 
#'   \code{checkFunctions=TRUE}). It is available as an S3 function for the user
#'   to be able to test their functions and debug them, which is difficult to do
#'   with a S4 validity function.
internalFunctionCheck<-function(clusterFUN,inputType,algorithmType,outputType){
	#--- Make small data
	N<-20
	set.seed(2851)
	x<-matrix(rnorm(N*3),ncol=N,nrow=3)
	set.seed(2851)
	diss<-matrix(runif(N^2,min=0,max=0.5),ncol=N,nrow=N)
	diss<-diss + t(diss)
	diag(diss)<-0
	#--- Set parameters
	if(algorithmType=="01") argList<-list(alpha=.5)	
	if(algorithmType=="K") argList<-list(k=2)	
	argList<-c(argList,list(cluster.only=TRUE,checkArgs=FALSE))
	#--- Run function on small data
	if(inputType %in% c("X")){
		test<-try(do.call(clusterFUN,c(list(x=x),argList)),silent=TRUE)
		if(inherits(test,"try-error")) return(paste("function test fails with input X",test[1]))
	}
	if(inputType %in% c("diss")){
		test<-try(do.call(clusterFUN,c(list(diss=diss),argList)),silent=TRUE)
		if(inherits(test,"try-error")) return(paste("function test fails with input diss",test[1]))
	}
	if(inputType %in% c("either")){
		test1<-try(do.call(clusterFUN,c(list(x=x,diss=NULL),argList)),silent=TRUE)
		if(inherits(test1,"try-error")) return(paste("function test fails with input x and NULL diss",test1[1]))
		test2<-try(do.call(clusterFUN,c(list(x=NULL,diss=diss),argList)),silent=TRUE)
		if(inherits(test2,"try-error")){
			return(paste("function test fails with input diss and NULL x",test2[1]))
		}
		test3<-try(do.call(clusterFUN,c(list(x=x,diss=diss),argList)),silent=TRUE)
		if(inherits(test3,"try-error")) return(paste("function test fails both diss and x input",test3[1]))
		if(outputType=="vector" & length(test1)!=N || length(test2)!=N || length(test3)!=N) return("clusterFUN does not return a vector equal to the number of observations")
	}
	else{
		if(outputType=="vector"){
			if(length(test)!=N) return("clusterFUN does not return a vector equal to the number of observations")
		}
	}
	return(TRUE)
}


.checkHasArgs<-function(FUN,requiredArgs){
    funArgs<-names(as.list(args(FUN)))
	all(requiredArgs %in% funArgs)
}
setValidity("ClusterFunction", function(object) {
    if(is.na(object@outputType)) {
      return("Must define outputType.")
    }
	if(!object@outputType%in%.outputTypes) return(paste("outputType must be one of",paste(.outputTypes,collapse=",")))
    #----
	# inputType
	#----
    if(is.na(object@inputType)) {
      return("Must define inputType.")
    }
	if(!object@inputType%in%.inputTypes) return(paste("inputType must be one of",paste(.inputTypes,collapse=",")))
	if(is.null(object@classifyFUN)& !is.na(object@inputClassifyType)) return("should not define inputClassifyType if classifyFUN is not defined")
    if(!is.null(object@classifyFUN) & is.na(object@inputClassifyType)) {
      return("Must define inputClassifyType if define classifyFUN.")
    }
	if(!is.null(object@classifyFUN) & !object@inputClassifyType%in%.inputTypes) return(paste("inputClassifyType must be one of",paste(.inputTypes,collapse=",")))
    #----
	# algorithmType
	#----
	if(is.na(object@algorithmType)) return("Must define algorithmType")
	if(!object@algorithmType%in%.algTypes) return(paste("algorithmType must be one of",paste(.algTypes,collapse=",")))
	### Add additional checks that 'k' and '01' work as expected... in particular that take particular arguments, etc. that 'k' and '01' are expected to take. 
		
		
	#----
	# function arguments are as needed
	#----
	if(object@inputType%in%c("X","either") & !.checkHasArgs(FUN=object@clusterFUN,requiredArgs="x")) return("inputType is either 'X' or 'either' but arguments to clusterFunction doesn't contain 'x'")
		if(object@inputType%in%c("diss","either") & !.checkHasArgs(FUN=object@clusterFUN,requiredArgs="diss")) return("inputType is either 'diss' or 'either' but arguments to clusterFunction doesn't contain 'diss'")	
	if(object@algorithmType=="K" & !.checkHasArgs(FUN=object@clusterFUN,requiredArgs=.requiredKArgs)) return("algorithmType is 'K' but arguments to clusterFunction doesn't contain",paste(.requiredKArgs,collapse=","))
	if(object@algorithmType=="01" & !.checkHasArgs(FUN=object@clusterFUN, requiredArgs=.required01Args)) return("algorithmType is '01' but arguments to clusterFunction doesn't contain", paste(.required01Args,collapse=","))
	
	
	if(object@checkFunctions){ #user can skip the check.
		out<-internalFunctionCheck(object@clusterFUN,object@inputType,object@algorithmType,object@outputType)
		if(!is.logical(out) || !out) return(out)
		
	}
	return(TRUE)
  })

#'@description The constructor \code{clusterFunction} creates an object of the
#'  class \code{ClusterFunction}.
#'  
#'@param clusterFUN function bassed to slot \code{clusterFUN}.
#'@param inputType character for slot \code{inputType}
#'@param algorithmType character for slot \code{inputType}
#'@param classifyFUN function for slot \code{classifyFUN}
#'@param outputType character for slot \code{outputType}
#'@param inputClassifyType character for slot \code{inputClassifyType}
#'@param requiredArgs character for slot \code{requiredArgs}
#'@param checkFunctions logical for whether to check the input functions with
#'  \code{internalFunctionsCheck}
#'@param ... arguments passed to different methods of \code{clusterFunction}  
#'@return A \code{ClusterFunction} object.
#'
#' @aliases clusterFunction
#' @rdname ClusterFunction-class
#' @export
setGeneric(
	name = "clusterFunction",
	def = function(clusterFUN,...) {
	  standardGeneric("clusterFunction")
	}
)
#' @rdname ClusterFunction-class
#' @export
setMethod(
	f = "clusterFunction",
	signature = signature("function"),
	definition = function(clusterFUN, inputType,outputType,algorithmType,inputClassifyType=NA_character_,requiredArgs=NA_character_,classifyFUN=NULL,checkFunctions=TRUE){
		out <- new("ClusterFunction",
	         clusterFUN=clusterFUN,
	         inputType=inputType,
	         algorithmType = algorithmType,
			 inputClassifyType=inputClassifyType,
			 classifyFUN=classifyFUN,
			 outputType=outputType,
			 requiredArgs=requiredArgs,
			 checkFunctions=checkFunctions
			 )
		validObject(out)
		return(out)
	}
)
