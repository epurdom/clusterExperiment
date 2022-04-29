#' @rdname plotCoClustering
#' @aliases plotCoClustering
#' @param invert logical determining whether the coClustering matrix should be
#'   inverted to be 1-coClustering for plotting. By default, if the diagonal
#'   elements are all zero, invert=TRUE, and otherwise invert=FALSE. If
#'   coClustering matrix is not a 0-1 matrix (e.g. if equal to a distance matrix
#'   output from \code{\link{clusterSingle}}, then the user should manually set
#'   this parameter to FALSE.)
#' @param saveDistance logical. When the \code{coClustering} slot contains
#'   indices of the clusterings or a NxB set of clusterings, the hamming
#'   distance will be calculated before running the plot. This argument
#'   determines whether the \code{ClusterExperiment} object with that distance
#'   in \code{coClustering} slot should be returned (so as to avoid
#'   re-calculating it in the future) or not.
#' @details \code{plotCoClustering} is a convenience function to plot the
#'   heatmap of the co-clustering distance matrix from the \code{coClustering}
#'   slot of a \code{ClusterExperiment} object (either by calculating the
#'   hamming distance of the clusterings stored in the \code{coClustering} slot,
#'   or the distance stored in the \code{coClustering} slot if it has already
#'   been calculated.
#' @export
setMethod(
  f = "plotCoClustering",
  signature = "ClusterExperiment",
  definition = function(data, invert, saveDistance=FALSE,...){
    typeCoCl<-.typeOfCoClustering(data)
    if(typeCoCl=="null") stop("coClustering slot is empty")
    if(typeCoCl=="indices"){
        #calculate the distance
        coClustering(data)<-.clustersHammingDistance(
            t(clusterMatrix(data,whichClusters=data@coClustering)))        
    }
    else{
        #Need to catch if not a symmetric matrix
        #(case of subsampling matrix)
        if(typeCoCl=="clusterings"){
            #calculate the distance
            coClustering(data)<-.clustersHammingDistance(
                t(coClustering(data)))        
            
        }
    }
    if(missing(invert)) 
        invert<-ifelse(all(diag(data@coClustering)==0),TRUE,FALSE)
    if(invert){
        #Make it a similarity matrix (better anyway for the sparse representation)
        coClustering(data) <- 
            as(as(1-data@coClustering,"sparseMatrix"),"symmetricMatrix")
    }
    	
    # Do all this so don't have to erase merge info from data (so can return calculated distance to user)
    fakeCE<-ClusterExperiment(as(data@coClustering,"matrix"),
                              clusterMatrix(data),
                              transformation=function(x){x},
                              checkTransformAndAssay=FALSE


    )
    for(sName in c('clusterMatrix', 'primaryIndex', 'clusterInfo', 
        'clusterTypes', 'dendro_samples', 'dendro_clusters', 'dendro_index', 
        'clusterLegend', 'orderSamples', 'merge_index', 
        'merge_dendrocluster_index', 'merge_method', 'merge_demethod', 
        'merge_cutoff', 'merge_nodeProp', 'merge_nodeMerge','colData')){
        slot(fakeCE, sName)<-slot(data,sName)
    }   
    plotHeatmap(fakeCE,isSymmetric=TRUE,clusterFeaturesData="all",...)
    if(saveDistance) return(data)
    else invisible()
  })
