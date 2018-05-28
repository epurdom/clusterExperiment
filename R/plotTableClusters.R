#' @name plotTableClusters
#' @title Plot heatmap of cross-tabs of 2 clusterings
#' @description Plot heatmap of cross-tabs of two clusterings
#' @aliases plotTableClusters
#' @aliases plotTableClusters,ClusterExperiment-method
#' @param object ClusterExperiment object (or matrix with table result)
#' @param ignoreUnassigned logical as to whether to ignore unassigned clusters in the plotting.
#' @param margin if NA, the counts from \code{tableClusters} will be plotted. Otherwise, \code{\link[package=base]{prop.table}} will be called and the argument \code{\margin} will be passed to \code{prop.table} to determine how proportions should be calculated.
#' @rdname plotTableClusters
#' @seealso \code{\link[package=base]{prop.table}}
#' @export
setMethod( 
  f = "plotTableClusters",
  signature = signature(object = "ClusterExperiment"),
  definition = function(object, whichClusters,ignoreUnassigned=FALSE,margin=NA,...){
    whCl<-.TypeIntoIndices(object,whClusters=whichClusters,cluster=FALSE)
    if(length(whCl)!=2) stop("invalid choice of 'whichClusters' -- must be exactly 2 clusterings chosen.")
		tableAll<-tableClusters(object,whichClusters=whCl,useNames=TRUE)
		if(!is.na(margin)){
			tableAll<-prop.table(tableAll,margin=margin)
		}
		plotTableClusters(tableAll,...)
}
)
		
		
#' @rdname plotTableClusters
#' @export
setMethod( 
  f = "plotTableClusters",
  signature = signature(object = "matrix"),
  definition = function(object, ...){
		tableAll<-object
		rNms<-rownames(tableAll)
		cNms<-colnames(tableAll)

		if(ignoreUnassigned){
			if("-1" %in% rNms || "-2" %in% rNms){
				wh<-which(! rNms %in% c("-1","-2"))
				if(length(wh)>0) tableAll<-tableAll[wh, ,drop=FALSE]
					else stop("All of the first clustering are unassigned, cannot use ignoreUnassigned=TRUE")
			}
			if("-1" %in% cNms || "-2" %in% cNms){
				wh<-which(! cNms %in% c("-1","-2"))
				if(length(wh)>0) tableAll<-tableAll[, wh,drop=FALSE]
					else stop("All of the second clustering are unassigned, cannot use ignoreUnassigned=TRUE")
			}
			
		}
		rankValues<-rank(sapply(1:ncol(tableAll),FUN=function(ii){
			whMax<-which.max(tableAll[,ii])
	
		}),ties.method="first")
		order2<-order(rankValues)
 	 	rData<-data.frame(colnames(tableAll))
		names(rData)<-clusterLabels(ceObj)[whCl[1]]

 	 	cData<-data.frame(colnames(tableAll)[order2])
		names(cData)<-clusterLabels(ceObj)[whCl[2]]
		
		#shouldn't need unclass if table method works!
		plotHeatmap(unclass(tableAll[,order2]),
			colData=cData, annRow=rData,
			clusterLegend=clusterLegend(ceObj)[whCl],
			clusterFeatures=cluster,clusterSamples=cluster
			)
		invisible(tableAll[,order2])
  	
  }
	)

#' @aliases tableClusters
#' @rdname plotTableClusters
#' @export
setMethod( 
  f = "tableClusters",
  signature = signature(object= "ClusterExperiment",whichClusters="character"),
  definition = function(object, whichClusters,...)
  {
    wh<-.TypeIntoIndices(object,whClusters=whichClusters)
    if(length(wh)==0) stop("invalid choice of 'whichClusters'")
    return(tableClusters(object,whichClusters=wh,...))
    
  })

#' @rdname plotTableClusters
#' @export
setMethod( 
  f = "tableClusters",
  signature = signature(object = "ClusterExperiment",whichClusters="missing"),
  definition = function(object, whichClusters,...)
  {
    tableClusters(object,whichClusters="primaryCluster")
    
  })

#' @rdname plotTableClusters
#' @param useNames for \code{tableClusters}, whether the output should be tabled
#'   with names (\code{useNames=TRUE}) or ids (\code{useNames=FALSE})
#' @export
setMethod( 
  f = "tableClusters",
  signature = signature(object = "ClusterExperiment",whichClusters="numeric"),
  definition = function(object, whichClusters, useNames=TRUE,...)
  { 
    if(useNames) numCluster<-clusterMatrixNames(object,whichClusters=whichClusters)
    else numCluster<-clusterMatrix(object)[,whichClusters]
    table(data.frame(numCluster))
  })
