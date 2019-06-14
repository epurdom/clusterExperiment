#' Barplot of 1 or 2 clusterings
#' 
#' Make a barplot of sample's assignments to clusters for single clustering, or 
#' cross comparison for two clusterings.
#' 
#' @aliases plotBarplot
#' @docType methods
#' @param object A matrix of with each column corresponding to a clustering and
#'   each row a sample or a \code{\link{ClusterExperiment}} object.
#' @param colPalette a vector of colors used for the different clusters. See
#'   details.
#' @param xNames names for the clusters on x-axis (i.e. clustering given 1st). 
#'   By default uses names of the 1st column of clusters matrix. See details.
#' @param legNames names for the clusters dividing up the 1st clusters (will 
#'   appear in legend). By default uses names of the 2nd cluster of clusters
#'   matrix. If only one clustering, \code{xNames} and \code{legNames} refer to
#'   the same clustering. See details.
#' @param legend whether to plot the legend
#' @param xlab label for x-axis. By default or if equal NULL the column name of 
#'   the 1st cluster of clusters matrix
#' @param legend.title label for legend. By default or if equal NULL the column 
#'   name of the 2st cluster of clusters matrix
#' @param labels if object is a ClusterExperiment object, then labels defines 
#'   whether the clusters will be identified by their names values in 
#'   clusterLegend (labels="names", the default) or by their clusterIds value in
#'   clusterLegend (labels="ids").
#' @param ... for \code{plotBarplot} arguments passed either to the method of
#'   \code{plotBarplot} for matrices or ultimately to \code{\link{barplot}}.
#' @details The first column of the cluster matrix will be on the x-axis and the
#'   second column (if present) will separate the groups of the first column.
#' @details All arguments of the matrix version can be passed to the 
#'   \code{ClusterExperiment} version. As noted above, however, some arguments 
#'   have different interpretations.
#' @details If \code{whichClusters = "workflow"}, then the most recent two 
#'   clusters of the workflow will be chosen where recent is based on the 
#'   following order (most recent first): final, mergeClusters, makeConsensus, 
#'   clusterMany.
#' @details \code{xNames}, \code{legNames} and \code{colPalette} should all be
#'   named vectors, with the names referring to the clusters they should match
#'   to (for \code{ClusterExperiment} objects, it is determined by the argument
#'   \code{labels} as to whether the names should match the cluster names or the
#'   clusterIds). \code{colPalette} and \code{legNames} must be same length of
#'   the number of clusters found in the second clustering, or if  only a single
#'   clustering, the same length as the number of clusters in that clustering.
#' @return A plot is produced, nothing is returned
#' @author Elizabeth Purdom
#' @inheritParams plotClusters
#' @inheritParams getClusterIndex
#' @examples
#' #clustering using pam: try using different dimensions of pca and different k
#' data(simData)
#'
#' cl <- clusterMany(simData, nReducedDims=c(5, 10, 50), reduceMethod="PCA",
#' clusterFunction="pam", ks=2:4, findBestK=c(TRUE,FALSE),
#' removeSil=c(TRUE,FALSE), makeMissingDiss=TRUE)
#'
#' plotBarplot(cl)
#' plotBarplot(cl,whichClusters=1:2)
#'
#' @rdname plotBarplot
#' @export
setMethod(
  f = "plotBarplot",
  signature = signature(object = "ClusterExperiment"),
  definition = function(object, whichClusters="primary",labels=c("names","ids"),...)
  { 
	whichClusters<-getClusterIndex(object,whichClusters=whichClusters,noMatch="throwError")
	whichClusters<-head(whichClusters,2) #limit it to 2
	labels<-match.arg(labels)
    legend<-clusterLegend(object)[[tail(whichClusters,1)]]
    args<-list(...)
    if(!"colPalette" %in% names(args)){
      useBuiltInColors<-TRUE
      colPalette<-legend[,"color"]
    }
    else{
      colPalette<-args$colPalette
      useBuiltInColors<-FALSE
      args<-args[-grep("colPalette",names(args))]
    }
    numClusterMat<-clusterMatrix(object,whichClusters=whichClusters)
    if(labels=="names"){
	  clusterMat<-convertClusterLegend(object,output="matrixNames")[,whichClusters]
      if(useBuiltInColors) names(colPalette)<-legend[,"name"]
      #make sure "-1" stays "-1" 
      clusterMat[numClusterMat== -1]<- "-1"
      clusterMat[numClusterMat== -2]<- "-2"
      if(any(legend[,"clusterIds"]== "-1") & useBuiltInColors){
        names(colPalette)[which(legend[,"clusterIds"]== "-1")]<-"-1"
      }
      if(any(legend[,"clusterIds"]== "-2") & useBuiltInColors){
        names(colPalette)[which(legend[,"clusterIds"]== "-2")]<-"-2"
      }
    }
    else{
      clusterMat<-numClusterMat
      if(useBuiltInColors) names(colPalette)<-legend[,"clusterIds"]
    }
    if(!"unassignedColor" %in% names(args) & any(legend[,"clusterIds"]== "-1")){
      args$unassignedColor<-legend[legend[,"clusterIds"]== "-1","color"]
    }
    if(!"missingColor" %in% names(args) & any(legend[,"clusterIds"]== "-2")){
      args$missingColor<-legend[legend[,"clusterIds"]== "-2","color"]
    }
    
    do.call("plotBarplot",c(list(object=clusterMat,colPalette=colPalette),args))
    
  })


#' @rdname plotBarplot
#' @export
setMethod(
  f = "plotBarplot",
  signature = signature(object = "vector"),
  definition = function(object, ...){
    plotBarplot(matrix(object,ncol=1),...)
  })

#' @rdname plotBarplot
#' @importFrom graphics barplot
#' @export
setMethod(
  f = "plotBarplot",
  signature = signature(object = "matrix"),
  definition = function(object,  xNames=NULL, legNames=NULL, legend=ifelse(ncol(object)==2,TRUE,FALSE), xlab=NULL, legend.title=NULL, unassignedColor="white", missingColor="grey", colPalette=NULL,...){
    if(ncol(object)>2) stop("if 'object' a matrix, must contain at most 2 clusters (i.e. 2 columns)")
    clLeg<-object[,1]
    if(is.null(xlab)) xlab<-colnames(object)[1]
    if(ncol(object)==2){
      clX<-object[,2]
      x<-t(table(clLeg,clX))
      nX<-nrow(x)
      xnames<-rownames(x)
    }
    else{
      x<-table(clLeg)
      nX<-length(x)
      xnames<-names(x)
    }
    #-----------
    ###Check colors given:
    #-----------
    if(is.null(colPalette)){
      colPalette<-bigPalette[seq_len(nX)]
      names(colPalette)<-xnames
    }
    if(is.null(names(colPalette)) & length(colPalette)>1) stop("must give names to colPalette") 
    if(length(colPalette)==1){
      if(ncol(object)==1){
        colPalette<-rep(colPalette,length=nX)
        names(colPalette)<-xnames
      }
      else stop("cannot give a single color to colPalette if comparing 2 clusters.")
    } 
    if(!all(xnames %in% names(colPalette))) stop("invalid names for colPalette")
    colPalette<-colPalette[rownames(x)]
    
    if(ncol(object)==2){
      pair<-TRUE
      
      #references is on the columns, alt on rows
      if(is.null(legend.title)) legend.title<-colnames(object)[2]	
      ###cols of x go along the x-axis ("Ref")
      ###rows of x divide up the cols of x ("Alt")
      ###Note: colPalette is of length of rows of x ("Alt")
      #change name and color of missing/unassigned
      whAltNotAssigned<-which(row.names(x)=="-1")
      whAltMissing<-which(row.names(x)=="-2")
      whRefNotAssigned<-which(colnames(x)=="-1")
      whRefMissing<-which(colnames(x)=="-2")
      if(length(whAltNotAssigned)>0){
        row.names(x)[whAltNotAssigned]<-"Not Assigned"
        colPalette[whAltNotAssigned]<-unassignedColor
      }
      if(length(whAltMissing)>0){
        row.names(x)[whAltMissing]<-"Not Included in Clustering"
        colPalette[whAltMissing]<-missingColor
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
        if(!identical(sort(names(legNames)),sort(colnames(x)))) stop("Invalid names for reference cluster names -- not match names of reference clusters")
        #put in same order
        legNames<-legNames[colnames(x)]
        labs<-paste(legNames," (",colnames(x),")",sep="")
      }
    }
    else{
      if(is.null(legNames)){
        if("-1" %in% names(x)){
          names(x)[names(x)=="-1"]<-"Not Assigned"
        }
        if("-2" %in% names(x)){
          names(x)[names(x)=="-2"]<-"Not Included in Clustering"
        }
        legNames<-names(x)
        names(legNames)<-names(x)
        labs<-legNames
      }
      else{
        if(is.null(names(legNames))) stop("must give names to legNames that match values of reference cluster")
        if(length(legNames)!=ncol(x)) stop("Invalid reference cluster names -- not same length as number of reference clusters")
        if(!identical(sort(names(legNames)),sort(names(x)))) stop("Invalid names for reference cluster names -- not match names of reference clusters")
        #put in same order
        legNames<-legNames[names(x)]
        labs<-paste(legNames," (",names(x),")",sep="")
      }
      
    }
    par(mar=c(9.1,4.1,4.1,1.1),las=2)
    bp<-graphics::barplot(x,col=colPalette,legend=legend,args.legend=list(title=legend.title), names.arg=rep("",length(labs)),xlab="",...)
    xsize<-diff(par("usr")[3:4])
    text(bp, par("usr")[3]+.0*xsize, labels=labs, srt=45, adj=c(1,2), xpd=TRUE)
    title(xlab=xlab,line=7)
    
  })


