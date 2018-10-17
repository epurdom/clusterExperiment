#' @name plotFeatureBoxplot
#' @title Plot boxplot of feature values by cluster
#' @description Plot a boxplot of the (transformed) values for a particular
#'   gene, separated by cluster
#' @aliases plotFeatureBoxplot
#'   plotFeatureBoxplot,ClusterExperiment,character-method
#' @rdname plotFeatureBoxplot
#' @export
setMethod(
  f = "plotFeatureBoxplot",
  signature = signature(object = "ClusterExperiment",feature="character"),
  definition = function(object, feature,whichCluster="primary",...)
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
#' @param plotUnassigned whether to plot the unassigned samples as a cluster
#'   (either -1 or -2)
#' @param whichAssay numeric or character specifying which assay to use. See
#'   \code{\link[SummarizedExperiment]{assay}} for details.
#' @param ... arguments passed to \code{\link{boxplot}}
#' @inheritParams plotReducedDims
#' @seealso \code{\link{boxplot}}
#' @rdname plotFeatureBoxplot
#' @return A plot is created. The output of boxplot is returned (see
#' \code{\link{boxplot}}), with additional elements \code{colors} and 
#'  \code{clusterIds} that gives the colors and internal ids that match 
#'  each boxplot (pulled from \code{clusterLegend} but in the order of plot)
#' @examples
#' data(simData)
#' #Create a ClusterExperiment object
#' cl <- clusterMany(simData, nReducedDims=c(5, 10, 50), reducedDim="PCA",
#'    clusterFunction="pam", ks=2:4, findBestK=c(TRUE,FALSE),
#'    removeSil=c(TRUE,FALSE))
#' #give names to the clusters
#' cl<-renameClusters(cl, whichCluster=1, 
#'    value=letters[1:nClusters(cl)[1]])
#' plotFeatureBoxplot(cl,feature=1)
#' @export
setMethod(
  f = "plotFeatureBoxplot",
  signature = signature(object = "ClusterExperiment",feature="numeric"),
  definition = function(object,  feature,whichCluster="primary",plotUnassigned=FALSE,unassignedColor=NULL,missingColor=NULL,main=NULL,whichAssay=1,...)
  {
    whCl<-getSingleClusterIndex(object,whichCluster,list(...))
    #get data:
    dat<-transformData(object, whichAssay=whichAssay)[feature,]
    clLegend<-clusterLegend(object)[[whCl]]
    uniqueNames<-length(unique(clLegend[,"name"]))==nrow(clLegend)
    #-----
		#put in alpha order by cluster name / id (except for the missing ones)
    #-----
    whNotMissing<-which(as.numeric(clLegend[,"clusterIds"])>0)
    if(length(whNotMissing)>0){
      orderedLegend<-clLegend[whNotMissing,,drop=FALSE]
      if(uniqueNames) orderedLegend<-orderedLegend[order(orderedLegend[,"name"]), ,drop=FALSE]
      else orderedLegend<-orderedLegend[order(orderedLegend[,"name"],orderedLegend[,"clusterIds"]), ,drop=FALSE]
    }
		else orderedLegend<-clLegend
    #-----
		#if plot unassigned, then add unassigned back to the end of group
    #-----
    whMissing<-which(as.numeric(clLegend[,"clusterIds"])<0)
    if(length(whMissing)>0 & plotUnassigned){
      if(length(whNotMissing)>0) orderedLegend<-rbind(orderedLegend,clLegend[whMissing, ,drop=FALSE])
      else orderedLegend<-clLegend[whMissing, ,drop=FALSE]
      if(!is.null(unassignedColor) & any(orderedLegend[,"clusterIds"]=="-1"))
        orderedLegend[orderedLegend[,"clusterIds"]=="-1","color"]<-unassignedColor
      if(!is.null(missingColor) & any(orderedLegend[,"clusterIds"]=="-2"))
        orderedLegend[orderedLegend[,"clusterIds"]=="-2","color"]<-missingColor
    }
    clLegend<-orderedLegend
    if(uniqueNames){
      cl<-clusterMatrixNamed(object)[,whCl]
      cl<-factor(cl,levels=orderedLegend[,"name"])
    }
    else{
      warning("Non-unique names for the clusters. Will order them by internal cluster ids")
      cl<-clusterMatrix(object)[,whCl,drop=FALSE]
    }
    if(!is.null(dim(cl))){
      if(ncol(cl)>1) stop("only a single cluster may be used in whichCluster")
      else cl<-cl[,1]
    }
    if(is.null(main)){
      if(!is.null(rownames(object))) main<-paste("Gene expression of",rownames(object)[feature])
      else paste("Gene expression of feature, index number",feature)
    }
    cl<-factor(cl,levels=if(uniqueNames) clLegend[,"name"] else clLegend[,"clusterIds"])
    bxpOut<-boxplot(as.vector(dat) ~ cl, names=clLegend[,"name"],main=main,col=clLegend[,"color"],...)
		bxpOut<-c(bxpOut,list(colors=clLegend[,"color"],clusterIds=clLegend[,"clusterIds"]))
		invisible(bxpOut)
  })
