#' @title Make hierarchy of set of clusters
#'
#' @aliases makeDendrogram,RSECClass-method
#' @name makeDendrogram
#' @rdname makeDendrogram
setMethod(
    f = "makeDendrogram",
    signature = "RSECClass",
    definition = function(x, whichCluster="primaryCluster",reduceMethod="mad",...)
    {
        passedArgs<-list(...)
        
        checkIgnore<-.depricateArgument(passedArgs=passedArgs,"filterIgnoresUnassigned","ignoreUnassignedVar") #06/2018 added in BioC 3.8
        if(!is.null(checkIgnore)){
            passedArgs<-checkIgnore$passedArgs
            filterIgnoresUnassigned<-checkIgnore$val
        }
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
              cluster=cl,calculateSample=TRUE,
              unassignedSamples=unassignedSamples),
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