#' @name clusterDendrogram
#' @title Accessing and manipulating the dendrograms 
#' @rdname clusterDendrogram
#' @description These functions are for accessing and manipulating the dendrograms 
#' stored in a \code{ClusterExperiment} object. We also document the required format of 
#' these dendrograms here.
#' @details Two dendrograms are stored in a \code{ClusterExperiment} object. One is a dendrogram that describes the hierarchy between the clusters (\code{@dendro_clusters}), and the other is a dendrogram that extends that hierarchy to include the clusters (\code{@dendro_samples}). The clustering that is used to make these hierarchies is saved in as well (\code{@dendro_index})
#' @details The dendrograms stored in a \code{ClusterExperiment} object
#' are required to be a \code{\link[package=phylobase]{phylo4d-class}} from the package \code{phylobase} (which 
#' uses the basic format of the S3 class \code{\link[package=ape]{phylo}} in the \code{ape} package to store 
#' the edges, \code{phylobase} makes it a S4 class with some useful helpers)
#' @details Additional requirements are made of the dendrograms to be a valid \code{ClusterExperiment} class. The cluster hierarchy must have data stored with it that has column names \code{NodeId},\code{Position},\code{ClusterIdDendro},\code{ClusterIdMerge}. The sample hierarchy must have data stored with it with column names \code{SampleIndex},\code{NodeId}, \code{Position}. These columns allow for the linking of the two dendrograms, as well as linking the dendrograms to the clustering and must not be changed by the user.
#' Furthermore there are a restriction as to setting the node labels of the dendrograms. The sample dendrogram is not allowed to have any labels. The cluster dendrogram can only have labels on the \emph{internal} nodes. Labels on the internal nodes of the cluster dendrogram can be set by the user (the function \code{nodeLabels<-} is defined to work on a \code{ClusterExperiment} object to make this easy); but the tips of the cluster dendrogram cannot have labels. The reason for these restrictions is so as to not duplicate storage of the names: tips of the cluster dendrogram are clusters by definition (so their names are stored in the \code{clusterLegend} slot), the internal nodes of the sample dendrogram are the same as the cluster dendrogram (so can be pulled from there), and the tips of the sample dendrogram are the individual samples (so the names are pulled from colnames of the object).
#' @seealso \code{\link{makeDendrogram}}, \code{\link[package=phylobase]{phylo4d-class}}, \code{\link[package=ape]{phylo}}
#' @return \code{clusterDendrogram} returns the dendrogram describing 
#' the clustering hierarchy.
#' @export
#' @aliases clusterDendrogram
setMethod(
  f = "clusterDendrogram",
  signature = "ClusterExperiment",
  definition = function(x) {
    return(x@dendro_clusters)
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
		.checkDendrogram(x)
	    return(x)
	}
	else{stop("the replacement value needs to have names that match the internal node ids")}
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
#' @export
#' @aliases sampleDendrogram
setMethod(
  f = "sampleDendrogram",
  signature = "ClusterExperiment",
  definition = function(x) {
    return(x@dendro_samples)
  }
)
