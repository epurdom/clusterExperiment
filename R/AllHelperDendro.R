#' @name clusterDendrogram
#' @title Accessing and manipulating the dendrograms 
#' @rdname clusterDendrogram
#' @description These functions are for accessing and manipulating the dendrograms 
#' stored in a \code{ClusterExperiment} object. We also document the required format of 
#' these dendrograms here.
#' @param x a ClusterExperiment object
#' @section The Stored Dendrograms
#' @details Two dendrograms are stored in a \code{ClusterExperiment} object. One
#'   is a dendrogram that describes the hierarchy between the clusters
#'   (\code{@dendro_clusters}), and the other is a dendrogram that extends that
#'   hierarchy to include the clusters (\code{@dendro_samples}). The clustering
#'   that is used to make these hierarchies is saved in as well
#'   (\code{@dendro_index})
#' @details The dendrograms stored in a \code{ClusterExperiment} object are
#'   required to be a \code{\link[phylobase]{phylo4d-class}} from the package
#'   \code{phylobase} (which uses the basic format of the S3 class
#'   \code{\link[ape]{phylo}} in the \code{ape} package to store the edges;
#'   \code{phylobase} makes it a S4 class with some useful helpers). This class
#'   allows storage of a data.frame of information corresponding to information
#'   on each node (see \code{\link[phylobase]{tdata}}).
#' @details Additional requirements are made of these dendrograms to be a valid
#'   for the slots of the \code{ClusterExperiment} class, described below,
#'   regarding the data that must be stored with it and the labels which can be
#'   assigned. Possible dendrograms can be checked for validity with the
#'   function \code{checkDendrogram}. The reason for the restrictions on the
#'   labels is so as to not duplicate storage of the names, see below
#'   descriptions for where to save user-defined names.
#' @section Cluster Hierarchy
#' @details \itemize{
#' \item{Labels}{The cluster dendrogram can only have labels on the
#' \emph{internal} nodes. Labels on the internal nodes of the cluster dendrogram
#' can be set by the user (the function \code{nodeLabels<-} is defined to work
#' on a \code{ClusterExperiment} object to make this easy). The tips of the
#' cluster dendrogram, corresponding to the clusters, cannot have labels; users
#' can set the labels (e.g. for plotting, etc) in the
#' \code{\link{clusterLegend}} slot using the function
#' \code{\link{renameClusters}}.}
#' \item{Data}{ 
#' The cluster hierarchy must have data stored with it that has the following
#' columns (additional ones are allowed):
#'  \itemize{
#'   \item{NodeId}{The permanent node id for the node. Must be of the format "NodeIdX" where "X" is a integer.}
#'   \item{Position}{The type of node, in terms of its position. The internal nodes should have the values "cluster hierarchy node" while the tips should have "cluster hierarchy tip".}
#'   \item{ClusterIdDendro}{Only for tips of dendrogram, should have the id that corresponds to its cluster in the clustering of the @dendro_index. Of the form "ClusterIdX", where "X" is the internal cluster id (see \code{\link{clusterLegend}}). Internal nodes should have NA values.}
#'   \item{ClusterIdMerge}{The id that corresponds to the cluster in the clustering of the @merge_index, if it exists. Of the form "ClusterIdX", where "X" is the internal cluster id (see \code{\link{clusterLegend}}}
#'  }
#' }
#' } 
#' @section Sample Hierarchy
#' @details \itemize{
#' \item{Labels}{The sample dendrogram is not allowed to have ANY labels. The
#' names for those nodes that correspond to the cluster hierarchy will be pulled
#' from the names in the cluster hierarchy for plotting, etc. and should be set
#' there (see above). Sample names for the tips of the tree will be pulled from
#' \code{colnames} of the object and should be set there. }
#' \item{Data}{ 
#' The cluster hierarchy must have data stored with it that has the following
#' columns (additional ones are allowed):
#'  \itemize{
#'   \item{NodeId}{For those nodes that correspond to a node in the cluster
#'   hierarchy, should have its permanent node id in this column. Other nodes
#'   should be NA.}
#'   \item{Position}{The type of node, in terms of its position. The internal
#'   nodes should have the values "cluster hierarchy node" while the tips should
#'   have "cluster hierarchy tip".}
#'   \item{SampleIndex}{Only for tips of dendrogram, the index of the sample at
#'   that tip to the samples in the object.}
#'  }
#' } 
#' }
#' @section Helper Functions
#' @seealso \code{\link{makeDendrogram}}, \code{\link[phylobase]{phylo4d-class}}, \code{\link[ape]{phylo}}
#' @return \code{clusterDendrogram} returns the dendrogram describing 
#' the clustering hierarchy.
#' @export
#' @aliases clusterDendrogram,ClusterExperiment-method
setMethod(
  f = "clusterDendrogram",
  signature = "ClusterExperiment",
  definition = function(x) {
    return(x@dendro_clusters)
  }
)

#' @rdname clusterDendrogram
#' @return \code{sampleDendrogram} returns the dendrogram that expands the
#'   cluster hierarchy to the samples.
#' @export
#' @aliases sampleDendrogram
setMethod(
  f = "sampleDendrogram",
  signature = "ClusterExperiment",
  definition = function(x) {
    return(x@dendro_samples)
  }
)


#' @importMethodsFrom phylobase nNodes nTips
#' @rdname clusterDendrogram
#' @return \code{nInternalNodes} returns the number of \emph{internal} nodes of
#'   the cluster hierarchy.
#' @export
#' @aliases nInternalNodes
setMethod(
  f = "nInternalNodes",
  signature = "ClusterExperiment",
  definition = function(x) {
    return(length(phylobase::nodeLabels(x@dendro_clusters)))
  }
)

#' @rdname clusterDendrogram
#' @return \code{nTips} returns the number of tips of the cluster hierarchy
#'   (same as number of non-negative clusters in the dendrogram clustering)
#' @export
#' @aliases nTips
setMethod(
  f = "nTips",
  signature = "ClusterExperiment",
  definition = function(x) {
    return(phylobase::nTips(x@dendro_clusters))
  }
)

#' @rdname clusterDendrogram
#' @return \code{nNodes} returns the number of total nodes of the cluster hierarchy
#' @export
#' @aliases nNodes
setMethod(
  f = "nNodes",
  signature = "ClusterExperiment",
  definition = function(x) {
    return(phylobase::nNodes(x@dendro_clusters))
  }
)

#' @rdname clusterDendrogram
#' @return \code{nodeLabels<-} sets the node labels of the \emph{cluster}
#' dendrogram
#' @export
#' @details For setting the node labels of the cluster dendrogram via 
#' \code{nodeLabels<-}, the replacement value has to have names that match the 
#' internal ids of the cluster dendrogram (the \code{NodeId} column). 
#' @aliases nodeLabels<-
#' @param value replacement value for \code{nodeLabels}. See details.
#' @param ... additional options passed to \code{nodeLabels<-} (ignored)
#' @importMethodsFrom phylobase nodeLabels<- nodeLabels
setMethod(
  f = "nodeLabels<-",
  signature = signature(x="ClusterExperiment", value="ANY"),
  definition = function(x, ..., value) {
	if(!is.null(names(value))){
		#match them
		m<-try(.matchToDendroData(inputValue=names(value),dendro=x@dendro_clusters,matchColumn="NodeId",returnColumn="NodeIndex"),silent=TRUE)
		if(inherits(m,"try-error")) stop("names of replace vector do not match in internal node id:",m)
		phylobase::labels(x@dendro_clusters)[m]<-value
		ch<-.checkDendrogram(x)
	    if(!is.logical(ch)) stop(ch)
		return(x)
	}
	else{stop("the replacement value needs to have names that match the internal node ids")}
  }
)

#' @rdname clusterDendrogram
#' @return \code{checkClusterDendrogram} checks if a \code{phylo4d} objects are
#'   valid for the cluster and sample dendrogram slots of the given
#'   \code{ClusterExperiment} object. Returns \code{TRUE} if there are no
#'   problems. Otherwise creates error.
#' @param dendroCluster a \code{phylo4d} to be check as for being cluster hierarchy
#' @param dendroSample a \code{phylo4d} to be check as for being cluster hierarchy
#' @param whichCluster which cluster are the dendrograms clustering.
#' @export
#' @aliases checkDendrogram
setMethod(
  f = "checkDendrogram",
  signature = signature(x="ClusterExperiment",dendroCluster="phylo4d",dendroSample="phylo4d"),
  definition = function(x,dendroCluster,dendroSample,whichCluster="dendro") {
	  whCl<-.convertSingleWhichCluster(x,whichCluster)
	  x@dendro_clusters<-dendroCluster
	  x@dendro_samples<-dendroSample
	  x@dendro_index<-whCl
	  ch<-.checkDendrogram(x)
      if(!is.logical(ch)) stop(ch)
	  else return(TRUE)
  }
)

#' @rdname clusterDendrogram
#' @return \code{nodeLabels} returns the node labels of the \emph{cluster}
#' dendrogram
#' @export
#' @aliases nodeLabels
setMethod(
  f = "nodeLabels",
  signature = signature(x="ClusterExperiment"),
  definition = function(x) {
	phylobase::nodeLabels(x@dendro_clusters)
  }
)

#' @rdname clusterDendrogram
#' @return \code{nodeIds} returns the internal (permanent) node ids of the
#'   \emph{cluster} dendrogram
#' @export
#' @aliases nodeIds
setMethod(
  f = "nodeIds",
  signature = signature(x="ClusterExperiment"),
  definition = function(x) {
	  phylobase::tdata(x@dendro_clusters,type="internal")[,"NodeId"]
  }
)

#' @rdname clusterDendrogram
#' @return \code{convertToDendrogram} returns the sample dendrogram converted to
#'   a \code{\link[stats]{dendrogram}} class.
#' @seealso \code{\link[stats]{dendrogram}}
#' @export
#' @aliases convertToDendrogram
setMethod(
	f="convertToDendrogram",
	signature = "ClusterExperiment",
	definition=function(x){
		dendroId <- .setNodeLabels(x,labelType="id")$dendro_samples
    return(.convertToDendrogram(dendroId,tipNames=colnames(x)))
		
	})
#' @importFrom ape as.hclust.phylo
#' @importFrom stats as.dendrogram
.convertToDendrogram<-function(x,tipNames=NULL){
	if(inherits(x,"dendrogram") )return(x)
	if(inherits(x,"phylo4")){
		x<-.convertToPhyClasses(x,"phylo")
	} 
	if(inherits(x,"phylo")){
		if(any(is.na(as.numeric(x$tip.label)))) stop("can only convert to dendrogram if sample indices are the tip.labels")
		x<-try(ape::as.hclust.phylo(x),FALSE)
		if(inherits(x, "try-error")) stop("coding error -- could not convert to hclust object. Reported error:",x)
		#need to convert integers in $order and (negative) entries in $merge into the values of $labels
		newOrder<-as.integer(as.numeric(x$labels))
		if(!all(sort(newOrder)==1:length(newOrder))) stop("could not convert hclust object because dendrogram does not have all consecutive sample indices as tip labels.")
		m<-match(abs(x$merge[x$merge<0]),x$order)
		newMerge<-x$merge
		newMerge[x$merge<0]<-newOrder[m]
		newMerge[x$merge<0]<-newMerge[x$merge<0]*sign(x$merge[x$merge<0])
		newMerge<-matrix(as.integer(newMerge),ncol=2)
		x$order<-newOrder
		x$merge<-newMerge
		if(is.null(tipNames)) x$labels<-NULL
		else x$labels<-tipNames[x$order]
	}
	if(inherits(x,"hclust")){
		x<-try(stats::as.dendrogram(x),FALSE)
		if(inherits(x, "try-error")) stop("could not convert from hclust to dendrogram object. Reported error:",x)
		if(!is.integer(unlist(x))) stop("coding error -- did not create integer valued dendrogram")
		return(x)
	}
	else{ stop("input x is not of hclust, phylo4 or phylo class")}
}