#' @title Make hierarchy of set of clusters
#' @description Extends the makeDendrogram function from clusterExperiment package to the output of RSEC.
#' @name makeDendrogram
#' @rdname makeDendrogram
#' @param reduceMethod character A character identifying what type of
#'   dimensionality reduction to perform before clustering. Can be either a
#'   value stored in either of reducedDims or filterStats slot or a built-in
#'   diminsionality reduction/filtering. The option "coCluster" will use the
#'   co-Clustering matrix stored in the 'coClustering' slot of the
#'   \code{RSECClass} object
#' @details The makeDendrogram function for the \code{RSECClass} differs from that of \code{Cluster} Experiment only in that it allows for \code{makeDendrogram} to be performed on the output of the coClustering matrix from \code{\link{makeConsensus}}
#' @seealso \code{\link{makeConsensus}}
#' @export
#' @aliases makeDendrogram,RSECClass-method
#' @name makeDendrogram
#' @rdname makeDendrogram
setMethod(
    f = "makeDendrogram",
    signature = "RSECClass",
    definition = function(x, whichCluster="primaryCluster",reduceMethod="mad",...)
    {
        passedArgs<-list(...)
        
        unassignedSamples<-match.arg(unassignedSamples)
        whCl<-getSingleClusterIndex(x,whichCluster,passedArgs)
        cl<-clusterMatrix(x)[,whCl]
        if(!is.na(mergeClusterIndex(x)) || !is.na(x@merge_dendrocluster_index)) x<-.eraseMerge(x)

	    if(length(reduceMethod)>1) stop('makeDendrogram only takes one choice of "reduceMethod" as argument')
       
		if(reduceMethod=="coCluster"){
          if(is.null(x@coClustering)) stop("Cannot choose 'coCluster' if 'coClustering' slot is empty. Run makeConsensus before running 'makeDendrogram' or choose another option for 'reduceMethod'")
          if(is.null(dimnames(x@coClustering))) stop("This ClusterExperiment object was made with an old version of clusterExperiment and did not give dimnames to the coClustering slot.")
          outlist<-do.call("makeDendrogram",c(list(
              x=as.dist(1-x@coClustering),
              cluster=cl,calculateSample=TRUE),
              passedArgs)) 
	        #Add clusterNames as Ids to cluster and sample dendrogram. Code replicated from ClusterExperiment. 
	 				x@dendro_clusters <- outlist$clusters
					phylobase::tipLabels(x@dendro_clusters)<-NA #erase any labels of the tips, internal nodes already have the defaults.
					x@dendro_samples <- outlist$samples #labels should have been erased already
	        x@dendro_index<-whCl
	        ch<-.checkDendrogram(x)
	        if(!is.logical(ch)) stop(ch)
					return(x)
		}
		else{
			#call ClusterExperiment version
			return(callNextMethod())
		}

        
	}
)
