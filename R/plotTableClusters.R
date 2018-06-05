#' @name plotTableClusters
#' @title Plot heatmap of cross-tabs of 2 clusterings
#' @description Plot heatmap of cross-tabulations of two clusterings
#' @aliases plotTableClusters,ClusterExperiment-method
#' @param object ClusterExperiment object (or matrix with table result)
#' @param ignoreUnassigned logical as to whether to ignore unassigned clusters in the 
#' plotting. This means they will also be ignored in the calculations of the proportions 
#' (if \code{margin} not NA).
#' @param margin if NA, the counts from \code{tableClusters} will be plotted. Otherwise, 
#' \code{\link[package=base]{prop.table}} will be called and the argument \code{margin} #' will be passed to \code{prop.table} to determine how proportions should be calculated.
#' @param whichClusters which clusters to tabulate. For \code{plotTableClusters} should
#'  be 2 clusters, for \code{tableClusters} can indicate arbitrary number.
#' @rdname plotTableClusters
#' @seealso \code{\link[package=base]{prop.table}}
#' @examples
#' #clustering using pam: try using different dimensions of pca and different k
#' data(simData)
#'
#' cl <- clusterMany(simData, nReducedDims=c(5, 10, 50), reducedDim="PCA",
#' clusterFunction="pam", ks=2:4, findBestK=c(TRUE,FALSE),
#' removeSil=c(TRUE,FALSE))
#' #give arbitrary names to clusters for demonstration
#' cl<-renameClusters(cl,value=letters[1:nClusters(cl)[1]],whichCluster=1)
#' tableClusters(cl,whichClusters=1:2)
#' #heatmap of the counts in each entry of table:
#' plotTableClusters(cl,whichClusters=1:2, ignoreUnassigned=TRUE)
#' @export
setMethod( 
  f = "plotTableClusters",
  signature = signature(object = "ClusterExperiment"),
  definition = function(object, whichClusters,ignoreUnassigned=FALSE,margin=NA,...){
    whCl<-.TypeIntoIndices(object,whClusters=whichClusters)
    if(length(whCl)!=2) stop("invalid choice of 'whichClusters' -- must be exactly 2 clusterings chosen.")
		tableAll<-tableClusters(object,whichClusters=whCl,useNames=TRUE)
		cL<-clusterLegend(object)[whCl]
		if(ignoreUnassigned){
			rNms<-rownames(tableAll)
			cNms<-colnames(tableAll)
			mat1<-cL[[whCl[1]]]
			mat2<-cL[[whCl[2]]]
			#convert names to clusterIds to check for missing
			rNms<-mat1[match(rNms,mat1[,"name"]),"clusterIds"]
			cNms<-mat2[match(cNms,mat2[,"name"]),"clusterIds"]
			
			if("-1" %in% rNms || "-2" %in% rNms){
				wh<-which(! rNms %in% c("-1","-2"))
				if(length(wh)>0){
					tableAll<-tableAll[wh, ,drop=FALSE]
					#deal with fact that plotHeatmap doesn't fix problem with rowData!
					mat1<-mat1[mat1[,"clusterIds"] %in% rNms[wh], ,drop=FALSE]
					cL[[1]]<-mat1
				}
				else stop("All of the first clustering are unassigned, cannot use ignoreUnassigned=TRUE")
				
			}
			if("-1" %in% cNms || "-2" %in% cNms){
				wh<-which(! cNms %in% c("-1","-2"))
				if(length(wh)>0) tableAll<-tableAll[, wh,drop=FALSE]
				else stop("All of the second clustering are unassigned, cannot use ignoreUnassigned=TRUE")
			}
			
		}
		
		if(!is.na(margin)){
			tableAll<-prop.table(tableAll,margin=margin)
		}
		plotTableClusters(tableAll,clusterLegend=cL,...)
}
)

		
#' @rdname plotTableClusters
#' @param cluster logical, whether to cluster the rows and columns of the table. Passed
#'  to arguments \code{clusterFeatures} AND \code{clusterSamples} of \code{plotHeatmap}.
#' @param clusterLegend list in \code{clusterLegend} format that gives colors for the
#'  clusters tabulated.
#' @param ... arguments passed on to \code{plotHeatmap}
#' @seealso \code{\link{plotHeatmap}}
#' @details Note that the cluster labels in \code{plotTableClusters} and 
#' \code{tableClusters} are converted to "proper" R names via \code{make.names}. This is 
#' because \code{tableClusters} calls the R function \code{table}, which makes this 
#' conversion
#' @inheritParams plotHeatmap
#' conversion.
#' @seealso \code{\link[base]{table}}
#' @export
setMethod( 
  f = "plotTableClusters",
  signature = signature(object = "table"),
  definition = function(object,clusterLegend=NULL,cluster=FALSE, ...){
		tableAll<-object
		varNames<-make.names(names(dimnames(tableAll)))
 	 	
		if(!cluster){
			rankValues<-rank(sapply(1:ncol(tableAll),FUN=function(ii){
				whMax<-which.max(tableAll[,ii])
	
			}),ties.method="first")
			order2<-order(rankValues)
			
		}
		else order2<-1:ncol(tableAll)
			
		rData<-data.frame(rownames(tableAll))
		names(rData)<-varNames[1]

 	 	cData<-data.frame(colnames(tableAll)[order2])
		names(cData)<-varNames[2]
		if(!is.null(clusterLegend)) names(clusterLegend)<-make.names(names(clusterLegend))
		plotHeatmap(tableAll[,order2],
			colData=cData, rowData=rData,
			clusterLegend=clusterLegend,
			clusterFeatures=cluster,clusterSamples=cluster,...
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
    if(useNames) numCluster<-clusterMatrixNamed(object,whichClusters=whichClusters)
    else numCluster<-clusterMatrix(object)[,whichClusters]
    table(data.frame(numCluster))
  })
