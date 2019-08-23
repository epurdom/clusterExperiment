#' @include AllChecks.R
#' @importClassesFrom HDF5Array HDF5Matrix
#' @importClassesFrom DelayedArray DelayedArray DelayedArray1 DelayedMatrix RleArray RleMatrix
#' @importClassesFrom phylobase phylo4 phylo4d
#' @rawNamespace import(phylobase, except = plot)
#' @import Matrix
setClassUnion("matrixOrMissing",members=c("matrix", "missing"))
setClassUnion("phylo4OrNULL",members=c("phylo4d", "NULL"))
setClassUnion("matrixOrNULL",members=c("matrix", "NULL"))
setClassUnion("listOrNULL",members=c("list", "NULL"))
setClassUnion("functionOrNULL",members=c("function", "NULL"))
setClassUnion("data.frameOrNULL",members=c("data.frame", "NULL"))
setClassUnion("matrixOrHDF5",members=c("matrix", "DelayedArray","HDF5Matrix")) #Sometimes it appears necessary to have HDF5Matrix listed separately, not sure why, but otherwise not finding it. 
setClassUnion("matrixOrHDF5OrNULL",members=c("matrix","DelayedArray","HDF5Matrix","NULL"))
setClassUnion("sparseOrHDF5OrNULL",members=c("sparseMatrix","DelayedArray","HDF5Matrix","NULL"))

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
#' important methods (e.g., \code{\link{clusterMany}}, \code{\link{makeConsensus}},
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
#' @slot merge_index index of the current merged cluster
#' @slot merge_cutoff value for the cutoff used to determine whether to merge
#'   clusters
#' @slot merge_dendrocluster_index index of the cluster merged with the current
#'   merge
#' @slot merge_nodeMerge data.frame of information about nodes merged in the
#'   current merge. See \code{\link{mergeClusters}}
#' @slot merge_nodeProp data.frame of information of proportion estimated
#'   non-null at each node of dendrogram. See \code{\link{mergeClusters}}
#' @slot merge_method character indicating method used for merging. See 
#'   \code{\link{mergeClusters}}
#' @slot merge_demethod character indicating the DE method used for merging. See 
#'   \code{\link{mergeClusters}}
#' @slot clusterTypes character vector with the origin of each column of
#' clusterMatrix.
#' @slot dendro_samples \code{\link[phylobase]{phylo4d}} object. A dendrogram
#'   containing the cluster relationship (leaves are samples; see
#'   \code{\link{clusterDendrogram}} for details).
#' @slot dendro_clusters \code{\link[phylobase]{phylo4d}} object. A dendrogram
#'   containing the cluster relationship (leaves are clusters; see see
#'   \code{\link{sampleDendrogram}} for details).
#' @slot dendro_index numeric. An integer giving the cluster that was used to
#'   make the dendrograms. NA_real_ value if no dendrograms are saved.
#' @slot coClustering \code{\link[Matrix]{sparseMatrix}} object. A sparse 
#' representation of the matrix with the cluster co-occurrence
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
#' @seealso \code{\link[Matrix]{sparseMatrix}} \code{\link[phylobase]{phylo4d}} 
#' @name ClusterExperiment-class
#' @aliases ClusterExperiment
#' @rdname ClusterExperiment-class
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import methods
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
    dendro_samples = "phylo4OrNULL",
    dendro_clusters = "phylo4OrNULL",
    dendro_index = "numeric",
    coClustering = "sparseOrHDF5OrNULL",
    clusterLegend="list",
    orderSamples="numeric",
	merge_index="numeric",
	merge_dendrocluster_index="numeric",
	merge_method="character",
	merge_demethod="character",
	merge_cutoff="numeric",
	merge_nodeProp="data.frameOrNULL",
	merge_nodeMerge="data.frameOrNULL"

    ),
    prototype=(    
        dendro_samples = NULL,
        dendro_clusters = NULL,
        dendro_index=NA_real_, 
        merge_index=NA_real_,
    	merge_dendrocluster_index=NA_real,
    	merge_method=NA_character,
    	merge_demethod=NA_character,
    	merge_cutoff=NA_real,
        merge_nodeProp=NULL,
    	merge_nodeMerge=NULL
        
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
#'  values `clusterSingle`, `clusterMany`, `mergeClusters`, `makeConsensus` are
#'  reserved for the clustering coming from the package workflow and should not
#'  be used when creating a new object with the constructor.
#'@param clusterInfo a list with information on the clustering (see Slots).
#'@param primaryIndex integer. Sets the `primaryIndex` slot (see Slots).
#'@param orderSamples a vector of integers. Sets the `orderSamples` slot (see
#'  Slots).
#'@param dendro_samples phylo4 object. Sets the `dendro_samples` slot (see Slots).
#'@param dendro_clusters phylo4 object. Sets the `dendro_clusters` slot (see
#'  Slots).
#'@param dendro_index numeric. Sets the \code{dendro_index} slot (see Slots).
#'@param coClustering matrix. Sets the \code{coClustering} slot (see Slots).
#'@param checkTransformAndAssay logical. Whether to check the content of the
#'  assay and given transformation function for whether they are valid.
#'@param merge_index integer. Sets the \code{merge_index} slot (see Slots)
#'@param merge_cutoff numeric. Sets the \code{merge_cutoff} slot (see Slots)
#'@param merge_dendrocluster_index integer. Sets the
#'  \code{merge_dendrocluster_index} slot (see Slots)
#'@param merge_demethod character, Sets the
#'  \code{merge_demethod} slot (see Slots)
#'@param merge_nodeMerge data.frame. Sets the \code{merge_nodeMerge} slot (see
#'  Slots)
#'@param merge_nodeProp data.frame. Sets the \code{merge_nodeProp} slot (see
#'  Slots)
#'@param merge_method character, Sets the \code{merge_method} slot (see Slots)
#'@param clusterLegend list, Sets the \code{clusterLegend} slot (see details).

#' @details The \code{clusterLegend} argument to \code{ClusterExperiment} 
#'  must be a valid clusterLegend format and match the values in \code{clusters}, 
#' in that the "clusterIds" column must matches the value in the clustering matrix
#'  \code{clusters}. If \code{names(clusterLegend)==NULL}, it is assumed that the 
#'  entries of \code{clusterLegend} are in the same order as the columns of 
#'  \code{clusters}. Generally, this is not a good way for users to set the
#' clusterLegend slot.
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
                        orderSamples=seq_len(ncol(object)),
                        dendro_samples=NULL,
                        dendro_index=NA_real_,
                        dendro_clusters=NULL,
                        coClustering=NULL,
                        merge_index=NA_real_,
                        merge_cutoff=NA_real_,
                        merge_dendrocluster_index=NA_real_,
                        merge_nodeProp=NULL,
                        merge_nodeMerge=NULL,
                        merge_method=NA_character_,
												merge_demethod=NA_character_,
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
      colnames(clusters)<-paste("cluster",seq_len(NCOL(clusters)),sep="")
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
    if(nrow(clusters)>0){
	    tmp<-.makeColors(clusters, colors=massivePalette,matchClusterLegend=clusterLegend,matchTo="givenIds")
	    if(is.null(clusterLegend)){
	    	clusterLegend<-tmp$colorList
	    } 
	    else{
	      #need to check matches the clusters, which .makeColors doesn't do.
	      clusterLegend<-unname(clusterLegend)
	      ch<-.checkClustersWithClusterLegend(clusters,clusterLegend)
	      if(!is.logical(ch)) stop(ch)
	 			clusterLegend<-tmp$colorList      
	    }
	    clustersNum<-tmp$numClusters
	    colnames(clustersNum)<-colnames(clusters)
		
    }
	else{
		clustersNum<-clusters
		clusterLegend<-lapply(seq_len(ncol(clusters)),function(ii){
			out<-matrix(nrow=0,ncol=3)
			colnames(out)<-c("clusterIds","color","name")
			return(out)
		})
	}
    #can just give object in constructor, and then don't loose any information!
    out <- new("ClusterExperiment",
               object,
               transformation=transformation,
               clusterMatrix = clustersNum,
               primaryIndex = primaryIndex,
               clusterTypes = unname(clusterTypes),
               clusterInfo=unname(clusterInfo),
               clusterLegend=unname(clusterLegend),
               orderSamples=seq_len(ncol(object)),
               dendro_samples=dendro_samples,
               dendro_clusters=dendro_clusters,
               dendro_index=dendro_index,
               merge_index=merge_index,
               merge_cutoff=merge_cutoff,
               merge_dendrocluster_index=merge_dendrocluster_index,
               merge_nodeProp=merge_nodeProp,
               merge_nodeMerge=merge_nodeMerge,
               merge_method=merge_method,
			   merge_demethod=merge_demethod,
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







