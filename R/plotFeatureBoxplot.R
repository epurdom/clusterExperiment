#' @name plotFeatureBoxplot
#' @title Plot boxplot of feature values by cluster
#' @description Plot a boxplot of the (transformed) values for a particular gene, separated by cluster
#' @aliases plotFeatureBoxplot plotFeatureBoxplot,ClusterExperiment,character,ANY-method
#' @rdname plotFeatureBoxplot
#' @export
setMethod(
  f = "plotFeatureBoxplot",
  signature = signature(object = "ClusterExperiment",whichCluster="character",feature="ANY"),
  definition = function(object, whichCluster,feature,...)
  {
  	wh<-.TypeIntoIndices(object,whClusters=whichCluster)
  	if(length(wh)==0) stop("invalid choice of 'whichCluster'")
  	if(length(wh)>1) stop("only a single clustering can be shown'whichCluster'")
     invisible(plotFeatureBoxplot(object,whichCluster=wh,feature=feature,...))

  })
  
#' @rdname plotFeatureBoxplot
#' @export
setMethod(
  f = "plotFeatureBoxplot",
  signature = signature(object = "ClusterExperiment",whichCluster="missing",feature="ANY"),
  definition = function(object, whichCluster,feature,...)
  {
    invisible(plotFeatureBoxplot(object,whichCluster="primaryCluster",feature=feature,...))

  })


#' @rdname plotFeatureBoxplot
#' @export
setMethod(
f = "plotFeatureBoxplot",
signature = signature(object = "ClusterExperiment",whichCluster="numeric",feature="character"),
definition = function(object, whichCluster,feature,...)
{
	m<-match(feature,rownames(object))
	if(is.na(m)) stop("feature does not match one of the rownames of this object")
	else invisible(plotFeatureBoxplot(object,whichCluster=whichCluster,feature=m,...))
})


  
#' @param feature identification of feature to plot, either row name or index
#' @param unassignedColor If not NULL, should be character value giving the
#'   color for unassigned (-2) samples (overrides \code{clusterLegend}) default.
#' @param missingColor If not NULL, should be character value giving the color
#'   for missing (-2) samples (overrides \code{clusterLegend}) default.
#' @param main title of plot. If NULL, given default title.
#' @param plotUnassigned whether to plot the unassigned samples as a cluster (either -1 or -2)
#' @param ... arguments passed to \code{\link{boxplot}}
#' @inheritParams plotReducedDims
#' @seealso \code{\link{boxplot}}
#' @rdname plotFeatureBoxplot
#' @return A plot is created. The output of boxplot is returned 
#' @examples
#' #clustering using pam: try using different dimensions of pca and different k
#' data(simData)
#'
#' cl <- clusterMany(simData, nReducedDims=c(5, 10, 50), reducedDim="PCA",
#' clusterFunction="pam", ks=2:4, findBestK=c(TRUE,FALSE),
#' removeSil=c(TRUE,FALSE))
#' clusterLegend(cl)[[1]][,"name"]<-letters[1:nClusters(cl,ignoreUnassigned = FALSE)[1]]
#' plotFeatureBoxplot(cl,feature=1)
#' @export
setMethod(
f = "plotFeatureBoxplot",
signature = signature(object = "ClusterExperiment",whichCluster="numeric",feature="numeric"),
definition = function(object, whichCluster, feature,plotUnassigned=FALSE,unassignedColor=NULL,missingColor=NULL,main=NULL,assay=NULL,...)
{
	#get data:
	if(is.null(assay)) dat<-transformData(object)[feature,]
	else{
		if(!is.matrix(assay)) assay<-assays(object)[[assay]]
		if(!all(dim(assay)==dim(object))) stop("input assay does not match the dimensions of object")
		if(is.null(assay)) stop("assay value does not match an assay in object")
		else dat<-assay[feature,]
	}
	
	clLegend<-clusterLegend(object)[[whichCluster]]
	uniqueNames<-length(unique(clLegend[,"name"]))==nrow(clLegend)
	#put in alpha order by cluster name / id:
	whNotMissing<-which(as.numeric(clLegend[,"clusterIds"])>0)
	if(length(whNotMissing)>0){
		orderedLegend<-clLegend[whNotMissing,]
		if(uniqueNames) orderedLegend<-orderedLegend[order(orderedLegend[,"name"]),]
		else orderedLegend<-orderedLegend[order(orderedLegend[,"name"],orderedLegend[,"clusterIds"]),]
	}
	#add missing if exist to end.
	whMissing<-which(as.numeric(clLegend[,"clusterIds"])<0)
	if(length(whMissing)>0 & !plotUnassigned){
		if(length(whNotMissing)>0) orderedLegend<-rbind(orderedLegend,clLegend[whMissing,])
		else orderedLegend<-clLegend[whMissing,]
		if(!is.null(unassignedColor) & any(orderedLegend[,"clusterIds"]=="-1")) 
			orderedLegend[orderedLegend[,"clusterIds"]=="-1","color"]<-unassignedColor
		if(!is.null(missingColor) & any(orderedLegend[,"clusterIds"]=="-2")) 
			orderedLegend[orderedLegend[,"clusterIds"]=="-2","color"]<-missingColor
	}
	clLegend<-orderedLegend
	if(uniqueNames){
		cl<-convertClusterLegend(object,output="matrixNames",whichClusters=whichCluster)
		cl<-factor(cl,levels=orderedLegend[,"name"])
	}
	else{
		warning("Non-unique names for the clusters. Will order them by internal cluster ids")
		cl<-clusterMatrix(object)[,whichCluster,drop=FALSE]
	}
	if(!is.null(dim(col))){
		if(ncol(cl)>1) stop("only a single cluster may be used in whichCluster")
		else cl<-cl[,1]
	}
	if(is.null(main)){
		if(!is.null(rownames(object))) main<-paste("Gene expression of",rownames(object)[feature])
		else paste("Gene expression of feature number",feature)
	}
	cl<-factor(cl,levels=if(uniqueNames) clLegend[,"name"] else clLegend[,"clusterIds"])
	invisible(boxplot(dat ~ cl, names=clLegend[,"name"],main=main,col=clLegend[,"color"],...))
	})