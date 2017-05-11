#' Barplot of 1 or 2 clusterings
#'
#' Make a barplot of sample's assignments to clusters for single clustering, or
#' cross comparison for two clusterings.
#'
#' @aliases plotBarplot
#' @docType methods
#' @param clusters A matrix of with each column corresponding to a clustering
#'   and each row a sample or a \code{\link{ClusterExperiment}} object. If a
#'   matrix, the function will plot the clusterings in order of this matrix, and
#'   their order influences the plot greatly.
#' @param whichClusters If numeric, a predefined order for the clusterings in
#'   the plot. If x is a \code{\link{ClusterExperiment}} object,
#'   \code{whichClusters} can be a character value identifying the
#'   \code{clusterTypess} to be used; alternatively \code{whichClusters}
#'   can be either 'all' or 'workflow' to indicate choosing all clusters or
#'   choosing all \code{\link{workflowClusters}}.
#' @param unassignedColor If ``-1'' in \code{clusters}, will be given this color
#'   (meant for samples not assigned to cluster).
#' @param missingColor If ``-2'' in clusters, will be given this color (meant
#'   for samples that were missing from the clustering, mainly when comparing
#'   clusterings run on different sets of samples)
#' @param colPalette a vector of colors used for the different clusters. Must be
#'   as long as the maximum number of clusters found in any single
#'   clustering/column given in \code{clusters} or will otherwise return an
#'   error.
#' @param xNames names for the first clusters (on x-axis). By default uses
#'   values in 1st cluster of clusters matrix
#' @param legNames names for the first clusters (in legend). By default uses
#'   values in 2nd cluster of clusters matrix
#' @param legend whether to plot the legend
#' @param xlab label for x-axis. By default or if equal NULL the column name of
#'   the 1st cluster of clusters matrix
#' @param legend.title label for legend. By default or if equal NULL the column
#'   name of the 2st cluster of clusters matrix
#' @param labels if clusters is a clusterExperiment object, then labels defines
#'   whether the clusters will be identified by their names values in
#'   clusterLegend (labels="names", the default) or by their clusterIds value in
#'   clusterLegend (labels="ids").
#' @param ... for \code{plotBarplot} arguments passed either to the method
#'   of \code{plotBarplot} for matrices or ultimately to \code{\link{barplot}}.
#' @details The first column of the cluster matrix will be on the x-axis and the
#'   second column will separate the groups of the first column.
#' @details All arguments of the matrix version can be passed to the
#'   \code{ClusterExperiment} version. As noted above, however, some arguments
#'   have different interpretations.
#' @details If \code{whichClusters = "workflow"}, then the most recent two
#'   clusters of the workflow will be chosen where recent is based on the
#'   following order (most recent first): final, mergeClusters, combineMany,
#'   clusterMany.
#'
#' @author Elizabeth Purdom

#' @export
#'
#' @examples
#' #clustering using pam: try using different dimensions of pca and different k
#' data(simData)
#'
#' cl <- clusterMany(simData, nPCADims=c(5, 10, 50), dimReduce="PCA",
#' clusterFunction="pam", ks=2:4, findBestK=c(TRUE,FALSE),
#' removeSil=c(TRUE,FALSE))
#'
#' plotBarplot(cl)
#' plotBarplot(cl,whichClusters=1:2)
#'
#' @rdname plotBarplot
setMethod(
  f = "plotBarplot",
  signature = signature(clusters = "ClusterExperiment",whichClusters="character"),
  definition = function(clusters, whichClusters,...)
  {
    wh<-head(.TypeIntoIndices(clusters,whClusters=whichClusters),2)
    return(plotBarplot(clusters,whichClusters=wh,...))

  })

#' @rdname plotBarplot
#' @export
setMethod(
    f = "plotBarplot",
    signature = signature(clusters = "ClusterExperiment",whichClusters="missing"),
    definition = function(clusters, whichClusters,...)
    {
      plotBarplot(clusters,whichClusters="primaryCluster")

    })

#' @rdname plotBarplot
#' @export
setMethod(
  f = "plotBarplot",
  signature = signature(clusters = "ClusterExperiment",whichClusters="numeric"),
  definition = function(clusters, whichClusters,labels=c("names","ids"),...)
  { 
  	labels<-match.arg(labels)
	legend<-clusterLegend(clusters)[[tail(whichClusters,1)]]
	colPalette<-legend[,"color"]
	numClusterMat<-clusterMatrix(clusters,whichClusters=whichClusters)
	if(labels=="names"){
		clusterMat<-convertClusterLegend(clusters,output="matrixNames")[,whichClusters]
		names(colPalette)<-legend[,"name"]
		#make sure "-1" stays "-1" 
		clusterMat[numClusterMat== -1]<- "-1"
		clusterMat[numClusterMat== -2]<- "-2"
		if(any(legend[,"clusterIds"]== "-1")){
			names(colPalette)[which(legend[,"clusterIds"]== "-1")]<-"-1"
		}
		if(any(legend[,"clusterIds"]== "-2")){
			names(colPalette)[which(legend[,"clusterIds"]== "-2")]<-"-2"
		}
	}
	else{
		clusterMat<-numClusterMat
		names(colPalette)<-legend[,"clusterIds"]
	}
	args<-list(...)
	if(!"unassignedColor" %in% names(args) & any(legend[,"clusterIds"]== "-1")){
		args$unassignedColor<-legend[legend[,"clusterIds"]== "-1","color"]
	}
	if(!"missingColor" %in% names(args) & any(legend[,"clusterIds"]== "-2")){
		args$missingColor<-legend[legend[,"clusterIds"]== "-2","color"]
	}
	#browser()
	do.call("plotBarplot",c(list(clusters=clusterMat,colPalette=colPalette),args))

  })

#' @rdname plotBarplot
setMethod(
  f = "plotBarplot",
  signature = signature(clusters = "ClusterExperiment",whichClusters="missing"),
  definition = function(clusters, whichClusters,...)
  {
    plotBarplot(clusters,whichClusters="primaryCluster",...)
  })



#' @rdname plotBarplot
setMethod(
  f = "plotBarplot",
  signature = signature(clusters = "vector",whichClusters="missing"),
  definition = function(clusters, whichClusters, ...){
	  plotBarplot(matrix(clusters,ncol=1),...)
  })

#' @rdname plotBarplot
setMethod(
  f = "plotBarplot",
  signature = signature(clusters = "matrix",whichClusters="missing"),
  definition = function(clusters, whichClusters, xNames=NULL, legNames=NULL, legend=TRUE, xlab=NULL, legend.title=NULL, unassignedColor="white", missingColor="grey", colPalette=bigPalette,...){
	if(ncol(clusters)>2) stop("clusters must at most 2 clusters (i.e. 2 columns)")
	clLeg<-clusters[,1]
	if(is.null(xlab)) xlab<-colnames(clusters)[1]
    if(ncol(clusters)==2){
		pair<-TRUE
		clX<-clusters[,2]
	    x<-t(table(clLeg,clX)) #references is on the columns, alt on rows
		if(is.null(legend.title)) legend.title<-colnames(clusters)[2]	
		#browser()
	   	
	    if(is.null(names(colPalette))) colPalette<-rep(colPalette,length=nrow(x))   
		else colPalette<-colPalette[rownames(x)]
		#change name and color of missing/unassigned
	    whAltNotAssigned<-which(row.names(x)=="-1")
	    whAltMissing<-which(row.names(x)=="-2")
	    whRefNotAssigned<-which(colnames(x)=="-1")
	    whRefMissing<-which(colnames(x)=="-2")
		if(length(whAltNotAssigned)>0){
			row.names(x)[whAltNotAssigned]<-"Not Assigned"
			colPalette[whRefNotAssigned]<-unassignedColor
		}
		if(length(whAltMissing)>0){
			row.names(x)[whAltMissing]<-"Not Included in Clustering"
			colPalette[whRefMissing]<-missingColor
		}
		if(length(whRefNotAssigned)>0){
			colnames(x)[whRefNotAssigned]<-"Not Assigned"
		}
		if(length(whRefMissing)>0){
			colnames(x)[whRefMissing]<-"Not Included in Clustering"
		}
		#change order so those are last
		if(any(length(whAltNotAssigned)>0 | length(whAltMissing)>0)){
			nm<-row.names(x)
			wh<-c(whAltNotAssigned,whAltMissing)
			x<-rbind(x[-wh,,drop=FALSE],x[wh,,drop=FALSE])
			rownames(x)<-c(nm[-wh],nm[wh]) #annoying, but otherwise still loose the names
		}
		if(any(length(whRefNotAssigned)>0 | length(whRefMissing)>0)){
			nm<-colnames(x)
			wh<-c(whRefNotAssigned,whRefMissing)
			x<-cbind(x[,-wh,drop=FALSE],x[,wh,drop=FALSE])
			colPalette<-c(colPalette[-wh],colPalette[wh])
			colnames(x)<-c(nm[-wh],nm[wh]) #annoying, but otherwise still loose the names
		}
		 if(is.null(legNames)){
	            legNames<-colnames(x)
	            names(legNames)<-colnames(x)
	            labs<-legNames
	    }
	    else{
	            if(is.null(names(legNames))) stop("must give names to legNames that match values of reference cluster")
	            if(length(legNames)!=ncol(x)) stop("Invalid reference cluster names -- not same length as number of reference clusters")
	            if(!all(sort(names(legNames))==sort(colnames(x)))) stop("Invalid names for reference cluster names -- not match names of reference clusters")
	            #put in same order
	            legNames<-legNames[colnames(x)]
	            labs<-paste(legNames," (",colnames(x),")",sep="")
	    }
	}
	else{
		x<-table(clLeg)
	    if(is.null(names(colPalette))) colPalette<-rep(colPalette,length=length(x))   
		else colPalette<-colPalette[names(x)]
	    if(is.null(legNames)){
	            legNames<-names(x)
	            names(legNames)<-names(x)
	            labs<-legNames
	    }
	    else{
	            if(is.null(names(legNames))) stop("must give names to legNames that match values of reference cluster")
	            if(length(legNames)!=ncol(x)) stop("Invalid reference cluster names -- not same length as number of reference clusters")
	            if(!all(sort(names(legNames))==sort(names(x)))) stop("Invalid names for reference cluster names -- not match names of reference clusters")
	            #put in same order
	            legNames<-legNames[names(x)]
	            labs<-paste(legNames," (",names(x),")",sep="")
	    }
	
	}
    par(mar=c(9.1,4.1,4.1,1.1),las=2)
    bp<-barplot(x,col=colPalette,legend=legend,args.legend=list(title=legend.title), names.arg=rep("",length(labs)),xlab="",...)
    xsize<-diff(par("usr")[3:4])
    text(bp, par("usr")[3]+.0*xsize, labels=labs, srt=45, adj=c(1,2), xpd=TRUE)
    title(xlab=xlab,line=7)

})


