#' @name plotReducedDims
#' @title Plot 2-dimensionsal representation with clusters
#' @description Plot a 2-dimensional representation of the data, color-code by a clustering.
#' @aliases plotReducedDims plotReducedDims,ClusterExperiment,character-method
#' @export
setMethod(
  f = "plotReducedDims",
  signature = signature(object = "ClusterExperiment",whichCluster="character"),
  definition = function(object, whichCluster,...)
  {
	wh<-.TypeIntoIndices(object,whClusters=whichCluster)
	if(length(wh)==0) stop("invalid choice of 'whichCluster'")
	if(length(wh)>1) stop("only a single clustering can be shown'whichCluster'")
    return(plotReducedDims(object,whichCluster=wh,...))

  })
  
#' @rdname plotReducedDims
#' @export
setMethod(
  f = "plotReducedDims",
  signature = signature(object = "ClusterExperiment",whichCluster="missing"),
  definition = function(object, whichCluster,...)
  {
    plotReducedDims(object,whichCluster="primaryCluster",...)

  })
  
#' @param object a ClusterExperiment object
#' @param whichCluster which clusters to show on the plot
#' @param reducedDim What dimensionality reduction method to use. Should match
#'   either a value in \code{reducedDimNames(object)} or one of the built-in 
#'   functions of \code{\link{listBuiltInReducedDims}()}
#' @param whichDims vector of length 2 giving the indices of which dimensions to
#'   show. The first value goes on the x-axis and the second on the y-axis.
#' @param clusterLegend matrix with three columns and colnames
#'   'clusterIds','name', and 'color' that give the color and name of the
#'   clusters in whichCluster. If NULL, pulls the information from
#'   \code{object}.
#' @param legend either logical, indicating whether to plot legend, or character
#'   giving the location of the legend (passed to \code{\link{legend}})
#' @param xlab Label for x axis
#' @param ylab Label for y axis
#' @param unassignedColor If not NULL, should be character value giving the
#'   color for unassigned (-2) samples (overrides \code{clusterLegend}) default.
#' @param missingColor If not NULL, should be character value giving the color
#'   for missing (-2) samples (overrides \code{clusterLegend}) default.
#' @param plotUnassigned logical as to whether unassigned (either -1 or -2
#'   cluster values) should be plotted. 
#' @param pch the point type, passed to \code{plot.default}
#' @param legendTitle character value giving title for the legend. If NULL, uses
#'   the clusterLabels value for clustering.
#' @param ... arguments passed to \code{\link{plot.default}}
#' @details If \code{plotUnassigned=TRUE}, and the color for -1 or -2 is set to
#'   "white", will be coerced to "lightgrey" regardless of user input to 
#'   \code{missingColor} and \code{unassignedColor}. If \code{plotUnassigned=FALSE}, 
#'   the samples with -1/-2 will not be plotted, nor will the category show up in the 
#'   legend.
#' @seealso \code{\link{plot.default}}, \code{\link{makeReducedDims}}, \code{\link{listBuiltInReducedDims}()}
#' @rdname plotReducedDims
#' @return A plot is created. Nothing is returned. 
#' @examples
#' #clustering using pam: try using different dimensions of pca and different k
#' data(simData)
#'
#' cl <- clusterMany(simData, nReducedDims=c(5, 10, 50), reducedDim="PCA",
#' clusterFunction="pam", ks=2:4, findBestK=c(TRUE,FALSE),
#' removeSil=c(TRUE,FALSE))
#'
#' plotReducedDims(cl,legend="bottomright")
#' @export
setMethod(
f = "plotReducedDims",
signature = signature(object = "ClusterExperiment",whichCluster="numeric"),
definition = function(object, whichCluster,
	reducedDim="PCA",whichDims=c(1:2),plotUnassigned=TRUE,legend=TRUE,legendTitle="",
	clusterLegend=NULL,unassignedColor=NULL,missingColor=NULL,pch=19,xlab=NULL,ylab=NULL,...)
{
	if(length(whichCluster)!=1) stop("whichCluster must identify a single clustering.")
	cluster<-clusterMatrix(object)[,whichCluster]
	if("col" %in% names(list(...))) stop("plotting parameter 'col' may not be passed to plot.default. Colors must be set via 'clusterLegend' argument.")
	if(is.null(clusterLegend)){
		clusterLegend<-clusterLegend(object)[[whichCluster]]
	}
	else{
		if(is.null(dim(clusterLegend)) || ncol(clusterLegend)!=3 || !all(c("clusterIds","name","color") %in% colnames(clusterLegend))) stop("clusterLegend must be a matrix with three columns and names 'clusterIds','name', and 'color'")
		if(!all(as.character(unique(cluster)) %in% clusterLegend[,"clusterIds"])) stop("'clusterIds' in 'clusterLegend' do not match clustering values")
	}
	clColor<-clusterLegend[,"color"]
	names(clColor)<-clusterLegend[,"name"]
	wh1<-which(clusterLegend[,"clusterIds"]=="-1")
	wh2<-which(clusterLegend[,"clusterIds"]=="-2")
	if(!is.null(unassignedColor) && is.character(unassignedColor)){
		if(length(wh1)>0) clColor[wh1]<-unassignedColor
	}
	if(!is.null(missingColor) && is.character(missingColor)){
		if(length(wh2)>0) clColor[wh2]<-missingColor
	}
	if(length(wh1)>0){
		if(plotUnassigned && clColor[wh1]=="white") clColor[wh1]<-"lightgrey"
		else if(!plotUnassigned) clColor[wh1]<-NA
		names(clColor)[wh1]<-"Unassigned"
		clusterLegend[wh1,"name"]<-"Unassigned"
		clColor<-c(clColor[-wh1],clColor[wh1])
	} 
	if(length(wh2)>0){
		if(plotUnassigned && clColor[wh2]=="white") clColor[wh2]<-"lightgrey"
		else if(!plotUnassigned) clColor[wh2]<-NA
		names(clColor)[wh2]<-"Missing"
		clusterLegend[wh2,"name"]<-"Missing"
		clColor<-c(clColor[-wh2],clColor[wh2])
	}
	
	clFactor<-factor(as.character(cluster),levels=clusterLegend[,"clusterIds"], labels=clusterLegend[,"name"])
	
	#################
	####Dim reduction stuff:
	#################
	if(length(whichDims)<2) 
		stop("whichDims must be a vector of length at least 2 giving the which dimensions of the dimensionality reduction to plot")
	redoDim<-FALSE
	if(!reducedDim %in% reducedDimNames(object) & reducedDim %in% listBuiltInReducedDims()) redoDim<-TRUE
	if(reducedDim %in% reducedDimNames(object)){ 
		#check if ask for higher dim than available
		if(max(whichDims)>NROW(object) || max(whichDims)>NCOL(object)) stop("Invalid value for whichDims: larger than row or column")
		if(max(whichDims)>ncol(reducedDim(object,type=reducedDim))) redoDim<-TRUE
	}	
		
	if(redoDim) object<-makeReducedDims(object,reducedDims=reducedDim,maxDims=max(whichDims))
	if(reducedDim %in% reducedDimNames(object)){
		
		dat<-reducedDim(object,type=reducedDim)[,whichDims]
	}
	else stop("'reducedDim' does not match saved dimensionality reduction nor built in methods.")

	if(length(whichDims)==2){
		if(is.null(ylab)) ylab<-paste("Dimension",whichDims[2])
		if(is.null(xlab)) xlab<-paste("Dimension",whichDims[1])
		plot(dat,col=clColor[as.character(clFactor)],pch=pch,xlab=xlab,ylab=ylab,...)
		doLegend<-FALSE
		if(is.logical(legend) && legend){
		  doLegend<-TRUE
		  legend<-"topright"
		}
		else{
		  if(is.character(legend)) doLegend<-TRUE
		  else stop("legend must either be logical or character value.")
		}
		if(doLegend){
			if(!plotUnassigned){
				wh1<-which(names(clColor)=="Unassigned")
				wh2<-which(names(clColor)=="Missing")
				if(length(wh1)>0) clColor<-clColor[-wh1]
				if(length(wh2)>0) clColor<-clColor[-wh2]
			}
		 	if(is.null(legendTitle)) legendTitle<-clusterLabels(object)[whichCluster]
			legend(x=legend,legend=names(clColor),fill=clColor,title=legendTitle)
		}
		
	}
	else if(length(whichDims)>2){
		colnames(dat)<-paste("Dim.",whichDims,sep="")
		pairs(data.frame(dat),col=clColor[as.character(clFactor)],pch=pch,...)
	}

	invisible(object)
	
})
