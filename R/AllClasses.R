#' @include AllChecks.R
setOldClass("dendrogram")
setClassUnion("matrixOrMissing",members=c("matrix", "missing"))
setClassUnion("dendrogramOrNULL",members=c("dendrogram", "NULL"))
setClassUnion("matrixOrNULL",members=c("matrix", "NULL"))
setClassUnion("listOrNULL",members=c("list", "NULL"))
setClassUnion("functionOrNULL",members=c("function", "NULL"))
setClassUnion("data.frameOrNULL",members=c("data.frame", "NULL"))
setClassUnion("matrixOrHDF5",members=c("matrix", "HDF5Matrix"))
setClassUnion("matrixOrHDF5OrNULL",members=c("matrix", "HDF5Matrix","NULL"))

#############################################################
#############################################################
############### ClusterExperiment Class #####################
#############################################################
#############################################################
#' @title Class ClusterExperiment
#'
#' @description \code{ClusterExperiment} is a class that extends
#' \code{SingleCellExperiment} and is used to store the data
#' and clustering information.
#'
#' @docType class
#' @aliases ClusterExperiment ClusterExperiment-class ClusterExperiment
#'
#' @description In addition to the slots of the \code{SingleCellExperiment}
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
#' @slot merge_index index of the current merged cluster
#' @slot merge_cutoff value for the cutoff used to determine whether to merge
#'   clusters
#' @slot merge_dendrocluster_index index of the cluster merged with the current
#'   merge
#' @slot merge_nodeMerge data.frame of information about nodes merged in the
#'   current merge
#' @slot merge_nodeProp data.frame of information of proportion estimated
#'   non-null at each node of dendrogram
#' @slot merge_method character indicating method used for merging
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
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import methods
#' @importFrom dendextend as.phylo.dendrogram
#' @export
#'
setClass(
  Class = "ClusterExperiment",
  contains = "SingleCellExperiment",
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
    orderSamples="numeric",
	merge_index="numeric",
	merge_dendrocluster_index="numeric",
	merge_method="character",
	merge_cutoff="numeric",
	merge_nodeProp="data.frameOrNULL",
	merge_nodeMerge="data.frameOrNULL"

    )
)

setValidity("ClusterExperiment", function(object) {
    ####
    #test that clusterInfo not have names
    ####
    if(!is.null(names(object@clusterInfo))) return("clusterInfo should not have names")

	ch<-.checkClusterMatrix(object)
	if(!is.logical(ch))  return(ch)
	ch<-.checkPrimaryIndex(object)
	if(!is.logical(ch))  return(ch)
	ch<-.checkClusterTypes(object)
	if(!is.logical(ch))  return(ch)
	ch<-.checkClusterLegend(object)
	if(!is.logical(ch))  return(ch)
	ch<-.checkOrderSamples(object)
	if(!is.logical(ch))  return(ch)
	ch<-.checkClusterLabels(object)
	if(!is.logical(ch))  return(ch)
	ch<-.checkMerge(object)
	if(!is.logical(ch))  return(ch)
	ch<-.checkDendrogram(object)
	if(!is.logical(ch))  return(ch)
	ch<-.checkCoClustering(object)
	if(!is.logical(ch))  return(ch)
	return(TRUE)
})

#'@description The constructor \code{ClusterExperiment} creates an object of the
#'  class \code{ClusterExperiment}. However, the typical way of creating these
#'  objects is the result of a call to \code{\link{clusterMany}} or
#'  \code{\link{clusterSingle}}.
#'
#'@description Note that when subsetting the data, the co-clustering and
#'  dendrogram information are lost.
#'
#'@param object a matrix or \code{SummarizedExperiment} or
#'  \code{SingleCellExperiment} containing the data that was clustered.
#'@param clusters can be either a numeric or character vector, a factor, or a
#'  numeric matrix, containing the cluster labels.
#'@param transformation function. A function to transform the data before
#'  performing steps that assume normal-like data (i.e. constant variance), such
#'  as the log.
#'@param ... The arguments \code{transformation}, \code{clusterTypes} and
#'  \code{clusterInfo} to be passed to the constructor for signature
#'  \code{SingleCellExperiment,matrix}.
#'
#'@return A \code{ClusterExperiment} object.
#'
#'@examples
#'
#'sce <- matrix(data=rnorm(200), ncol=10)
#'labels <- gl(5, 2)
#'
#'cc <- ClusterExperiment(sce, as.numeric(labels), transformation =
#'function(x){x})
#'
#' @rdname ClusterExperiment-class
#' @export
setGeneric(
  name = "ClusterExperiment",
  def = function(object,  clusters,...) {
    standardGeneric("ClusterExperiment")
  }
)
#' @rdname ClusterExperiment-class
#' @export
setMethod(
  f = "ClusterExperiment",
  signature = signature("matrixOrHDF5","ANY"),
  definition = function(object, clusters, ...){
    ClusterExperiment(SummarizedExperiment(object), clusters, ...)
  })
#' @rdname ClusterExperiment-class
setMethod(
  f = "ClusterExperiment",
  signature = signature("SummarizedExperiment", "ANY"),
  definition = function(object, clusters, ...){
  ClusterExperiment(as(object, "SingleCellExperiment"),clusters,...)
})
#' @rdname ClusterExperiment-class
setMethod(
  f = "ClusterExperiment",
  signature = signature("SingleCellExperiment", "numeric"),
  definition = function(object, clusters, ...){
  ClusterExperiment(object,matrix(clusters, ncol=1),...)
})
#' @rdname ClusterExperiment-class
setMethod(
  f = "ClusterExperiment",
  signature = signature("SingleCellExperiment","character"),
  definition = function(object, clusters,...){
    ClusterExperiment(object,matrix(clusters,ncol=1),...)
    })
#' @rdname ClusterExperiment-class
setMethod(
  f = "ClusterExperiment",
  signature = signature("SingleCellExperiment","factor"),
  definition = function(object, clusters,...){
    clusters <- as.character(clusters)
    ClusterExperiment(object,clusters,...)
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
#'@param dendro_index numeric. Sets the \code{dendro_index} slot (see Slots).
#'@param dendro_outbranch logical. Sets the \code{dendro_outbranch} slot (see
#'  Slots).
#'@param coClustering matrix. Sets the \code{coClustering} slot (see Slots).
#'@param checkTransformAndAssay logical. Whether to check the content of the
#'  assay and given transformation function for whether they are valid.
#'@param merge_index integer. Sets the \code{merge_index} slot (see Slots)
#'@param merge_cutoff numeric. Sets the \code{merge_cutoff} slot (see Slots)
#'@param merge_dendrocluster_index integer. Sets the
#'  \code{merge_dendrocluster_index} slot (see Slots)
#'@param merge_nodeMerge data.frame. Sets the \code{merge_nodeMerge} slot (see
#'  Slots)
#'@param merge_nodeProp data.frame. Sets the \code{merge_nodeProp} slot (see
#'  Slots)
#'@param merge_method character, Sets the \code{merge_method} slot (see Slots)
#'@param clusterLegend list, Sets the \code{clusterLegend} slot (see Slots)
#' @details The \code{ClusterExperiment} constructor function gives
#'   clusterLabels based on the column names of the input
#'   matrix/SingleCellExperiment. If missing, will assign labels
#'   "cluster1","cluster2", etc.
#' @details Note that the validity check when creating a new
#'   \code{ClusterExperiment} object with \code{new} is less extensive than when
#'   using \code{ClusterExperiment} function with
#'   \code{checkTransformAndAssay=TRUE} (the default). Users are advised to use
#'   \code{ClusterExperiment} to create new \code{ClusterExperiment} objects.
setMethod(
  f = "ClusterExperiment",
  signature = signature("SingleCellExperiment","matrix"),
  definition = function(object, clusters,
                        transformation=function(x){x},
                        primaryIndex=1,
                        clusterTypes="User",
                        clusterInfo=NULL,
                        orderSamples=1:ncol(object),
                        dendro_samples=NULL,
                        dendro_index=NA_real_,
                        dendro_clusters=NULL,
                        dendro_outbranch=NA,
                        coClustering=NULL,
                        merge_index=NA_real_,
						merge_cutoff=NA_real_,
                        merge_dendrocluster_index=NA_real_,
                        merge_nodeProp=NULL,
                        merge_nodeMerge=NULL,
                        merge_method=NA_character_,
						clusterLegend=NULL,
                        checkTransformAndAssay=TRUE
  ){
    if(NCOL(object) != nrow(clusters)) {
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
    tmp<-.makeColors(clusters, colors=massivePalette)
    if(is.null(clusterLegend)) clusterLegend<-tmp$colorList
	else{
		clusterLegend<-unname(clusterLegend)
		ch<-.checkIndClusterLegend(clusters,clusterLegend)
		if(!is.logical(ch)) stop(ch)
		#need to grab colors/names in given clusterLegend
		autoLegend<-tmp$colorList
		clusterLegend<-mapply(clusterLegend,autoLegend,FUN=function(orig,auto){
			m<-match(orig[,"clusterIds"],auto[,"name"])
			if(any(is.na(m))) stop("coding error -- do not have all of original clusters in new clusterLegend") #shouldn't happen!
			orig[,"clusterIds"]<-auto[m,"clusterIds"]
			return(orig)

		},SIMPLIFY=FALSE)
	}
    clustersNum<-tmp$numClusters
    colnames(clustersNum)<-colnames(clusters)
    #can just give object in constructor, and then don't loose any information!
    out <- new("ClusterExperiment",
               object,
               transformation=transformation,
               clusterMatrix = clustersNum,
               primaryIndex = primaryIndex,
               clusterTypes = unname(clusterTypes),
               clusterInfo=unname(clusterInfo),
               clusterLegend=unname(clusterLegend),
               orderSamples=1:ncol(object),
               dendro_samples=dendro_samples,
               dendro_clusters=dendro_clusters,
               dendro_index=dendro_index,
               dendro_outbranch=dendro_outbranch,
               merge_index=merge_index,
			   merge_cutoff=merge_cutoff,
               merge_dendrocluster_index=merge_dendrocluster_index,
               merge_nodeProp=merge_nodeProp,
               merge_nodeMerge=merge_nodeMerge,
               merge_method=merge_method,
               coClustering=coClustering

    )
    if(checkTransformAndAssay){
      chass<-.checkAssays(out)
      if(is.logical(chass) && !chass) stop(chass)
      chtrans<-.checkTransform(out)
      if(is.logical(chtrans) && !chtrans) stop(chtrans)
    }
    return(out)
  })


#############################################################
#############################################################
############### ClusterFunction Class #####################
#############################################################
#############################################################

#' @title Class ClusterFunction
#'
#' @description \code{ClusterFunction} is a class for holding functions that can
#'   be used for clustering in the clustering algorithms in this package.
#'
#' @docType class
#' @aliases ClusterFunction ClusterFunction-class ClusterFunction
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
#'   input 'x' is not the same samples that were input to the ClusterFunction
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
#'   argument of \code{\link{mainClustering}} can preserve this ordering
#'   (default is to order by size).
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
#'   a cluster. If not, the post-processing steps in
#'   \code{\link{mainClustering}} such as \code{findBestK} and \code{removeSil}
#'   may not operate correctly since they rely on silhouette distances.
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
	if(object@inputType%in%c("X","either") & !.checkHasArgs(FUN=object@clusterFUN,requiredArgs="x")) return("inputType is either 'X' or 'either' but arguments to ClusterFunction doesn't contain 'x'")
		if(object@inputType%in%c("diss","either") & !.checkHasArgs(FUN=object@clusterFUN,requiredArgs="diss")) return("inputType is either 'diss' or 'either' but arguments to ClusterFunction doesn't contain 'diss'")
	if(object@algorithmType=="K" & !.checkHasArgs(FUN=object@clusterFUN,requiredArgs=.requiredKArgs)) return("algorithmType is 'K' but arguments to ClusterFunction doesn't contain",paste(.requiredKArgs,collapse=","))
	if(object@algorithmType=="01" & !.checkHasArgs(FUN=object@clusterFUN, requiredArgs=.required01Args)) return("algorithmType is '01' but arguments to ClusterFunction doesn't contain", paste(.required01Args,collapse=","))


	if(object@checkFunctions){ #user can skip the check.
		out<-internalFunctionCheck(object@clusterFUN,object@inputType,object@algorithmType,object@outputType)
		if(!is.logical(out) || !out) return(out)

	}
	return(TRUE)
  })

#'@description The constructor \code{ClusterFunction} creates an object of the
#'  class \code{ClusterFunction}.
#'
#'@param clusterFUN function passed to slot \code{clusterFUN}.
#'@param inputType character for slot \code{inputType}
#'@param algorithmType character for slot \code{inputType}
#'@param classifyFUN function for slot \code{classifyFUN}
#'@param outputType character for slot \code{outputType}
#'@param inputClassifyType character for slot \code{inputClassifyType}
#'@param requiredArgs character for slot \code{requiredArgs}
#'@param checkFunctions logical for whether to check the input functions with
#'  \code{internalFunctionsCheck}
#'@param ... arguments passed to different methods of \code{ClusterFunction}
#'@return A \code{ClusterFunction} object.
#'
#' @aliases ClusterFunction
#' @rdname ClusterFunction-class
#' @export
setGeneric(
	name = "ClusterFunction",
	def = function(clusterFUN,...) {
	  standardGeneric("ClusterFunction")
	}
)
#' @rdname ClusterFunction-class
#' @export
setMethod(
	f = "ClusterFunction",
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
		return(out)
	}
)




