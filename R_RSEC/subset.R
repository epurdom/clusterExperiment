#' @title Functions to subset ClusterExperiment Objects
#' @description These functions are used to subset ClusterExperiment objects,
#'   either by removing samples, genes, or clusterings
#' @name subset
#' @param x a ClusterExperiment object.
#' @inheritParams ClusterExperiment-class
#' @inheritParams getClusterIndex
#' @return A \code{\link{ClusterExperiment}} object.
#' @details \code{removeClusterings} removes the clusters given by
#'   \code{whichClusters}. If the \code{primaryCluster} is one of the clusters
#'   removed, the \code{primaryClusterIndex} is set to 1 and the dendrogram and
#'   coclustering matrix are discarded and orderSamples is set to
#'   \code{1:NCOL(x)}.
#' @return \code{removeClusterings} returns a \code{ClusterExperiment} object,
#'   unless all clusters are removed, in which case it returns a
#'   \code{\link{SingleCellExperiment}} object.
#' @examples
#' #load CE object
#' data(rsecFluidigm)
#' # remove the mergeClusters step from the object
#' clusterLabels(rsecFluidigm)
#' test<-removeClusterings(rsecFluidigm,whichClusters="mergeClusters")
#' clusterLabels(test)
#' tableClusters(rsecFluidigm)
#' test<-removeClusters(rsecFluidigm,whichCluster="mergeClusters",clustersToRemove=c("m01","m04"))
#' tableClusters(test,whichCluster="mergeClusters")
#' @export
#' @aliases removeClusterings,RSECClass-method
setMethod(
    f = "removeClusterings",
    signature = signature("ClusterExperiment"),
    definition = function(x, whichClusters) {
        whichClusters<-getClusterIndex(object=x,whichClusters=whichClusters,noMatch="throwError")
				retval<-callNextMethod(x=x, whichClusters=whichClusters)
        
        ## Fix CoClustering information
        ## Erase if any are part of clusters to remove
        coMat<-coClustering(x)
        typeCoCl<-.typeOfCoClustering(x)
        if(typeCoCl=="indices"){
            if(any(coMat %in% whichClusters)){
                warning("removing clusterings that were used in makeConsensus (i.e. stored in CoClustering slot). Will delete the coClustering slot")
                coMat<-NULL
            }
            else{
                #Fix so indexes the right clustering...
                coMat<-match(coMat,
                     seq_len(NCOL(clusterMatrix(x)))[-whichClusters])
            }
            
        }
        
        #fix merge info:
        #erase merge info if either dendro or merge index deleted.
        if(mergeClusterIndex(x) %in% whichClusters | x@merge_dendrocluster_index %in% whichClusters){
            x<-.eraseMerge(x)
            merge_index<-x@merge_index
            merge_dendrocluster_index<-x@merge_dendrocluster_index
        }
        else{
            merge_index<-match(x@merge_index, seq_len(NCOL(clusterMatrix(x)))[-whichClusters])
            merge_dendrocluster_index<-match(x@merge_dendrocluster_index,  seq_len(NCOL(clusterMatrix(x)))[-whichClusters])
            
        }
				
        retval<-RSECClass(
            retval,
            merge_index=merge_index,
            merge_dendrocluster_index=merge_dendrocluster_index,
            merge_cutoff=x@merge_cutoff,
            merge_nodeProp=x@merge_nodeProp,
            merge_nodeMerge=x@merge_nodeMerge,
            merge_method=x@merge_method,
            merge_demethod=x@merge_demethod,                              
            coClustering=coMat
        )
        return(retval)
    }
)