#' @title Plot 2-dimensionsal representation with clusters
#' @description Plot a 2-dimensional representation of the data, color-code by a clustering.
#' @rdname plot2D
#' @export
setMethod(
  f = "plot2D",
  signature = signature(object = "ClusterExperiment",whichClusters="character"),
  definition = function(object, whichClusters,...)
  {
	wh<-.TypeIntoIndices(object,whClusters=whichClusters)
	if(length(wh)==0) stop("invalid choice of 'whichClusters'")
	wh<-head(wh,2) #limit it to 2
    return(plot2D(object,whichClusters=wh,...))

  })
  
#' @rdname plot2D
#' @export
setMethod(
  f = "plot2D",
  signature = signature(object = "ClusterExperiment",whichClusters="missing"),
  definition = function(object, whichClusters,...)
  {
    plot2D(object,whichClusters="primaryCluster")

  })
  
#' @param object a ClusterExperiment object
#' @param whichClusters which clusters to show on the plot
#' @param dimReduce What dimensionality reduction method to use (currently only PCA). See \code{\link{transform}} for how this is implemented.
#' @param whichDims vector of length 2 giving the indices of which dimensions to show. The first value goes on the x-axis and the second on the y-axis. 
#' @param clusterLegend matrix with three columns and colnames 'clusterIds','name', and 'color' that give the color and name of the clusters in whichClusters. If NULL, pulls the information from \code{object}.
#' @param pch the point type, passed to \code{plot.default}
#' @param ... arguments passed to \code{\link{plot.default}}
#' @seealso \code{\link{plot.default}}, \code{\link{transform}}
#' @rdname plot2D
#' @examples
#' #clustering using pam: try using different dimensions of pca and different k
#' data(simData)
#'
#' cl <- clusterMany(simData, nPCADims=c(5, 10, 50), dimReduce="PCA",
#' clusterFunction="pam", ks=2:4, findBestK=c(TRUE,FALSE),
#' removeSil=c(TRUE,FALSE))
#'
#' plot2D(cl)
#' plot2D(cl,whichClusters=1:2)

#' @export
setMethod(
f = "plot2D",
signature = signature(object = "ClusterExperiment",whichClusters="numeric"),
definition = function(object, whichClusters,dimReduce=c("PCA"),whichDims,clusterLegend=NULL,pch=19,...)
{
	dimReduce<-match.arg(dimReduce)
	if(length(whichDims)!=2) stop("whichDims must be a vector of length 2 giving the which dimensions of the dimensionality reduction to plot")
	if(length(whichClusters)!=1) stop("whichClusters must identify a single clustering.")
    
	cluster<-clusterMatrix(object)[,whichClusters]
	if("col" %in% names(list(...))) stop("plotting parameter 'col' may not be passed to plot.default. Colors must be set via 'clusterLegend' argument.")
	if(is.null(clusterLegend)){
		clusterLegend<-clusterLegend(object)[[whichClusters]]
	}
	else{
		if(is.null(dim(clusterLegend)) || ncol(clusterLegend)!=3 || !all(c("clusterIds","name","color") %in% colnames(clusterLegend))) stop("clusterLegend must be a matrix with three columns and names 'clusterIds','name', and 'color'")
		if(!all(as.character(unique(cluster)) %in% clusterLegend[,"clusterIds"])) stop("'clusterIds' in 'clusterLegend' do not match clustering values")
	}
	clFactor<-factor(as.character(cluster),levels=clusterLegend[,"clusterIds"], labels=clusterLegend[,"name"])
	clColor<-clusterLegend[,"color"]
	names(clColor)<-clusterLegend[,"name"]

	transObj <- .transData(assay(object), nPCADims=max(whichDims), nVarDims=NA,
                           dimReduce=dimReduce, transFun=transformation(object),clustering=dimReduceCl)
	dat <- transObj$x[,whichDims]
	plot(dat,col=clColor[clFactor],pch=pch,...)
	
	
	
})
