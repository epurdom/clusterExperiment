#' @rdname plotClustersWorkflow
#' @aliases plotClustersWorkflow
#' @title A plot of clusterings specific for clusterMany and workflow
#'   visualization
#' @description A realization of \code{\link{plotClusters}} call specific to
#'   separating out the results of \code{clusterMany} and other clustering
#'   results.
#' @param object A \code{ClusterExperiment} object on which
#'   \code{\link{clusterMany}} has been run
#' @param whichClusterMany indicate which clusterings to plot in the bulk of the
#'   plot, see argument \code{whichClusters} of \code{\link{getClusterIndex}}
#'   for description of format allowed.
#' @param whichClusters which clusterings to "highlight", i.e draw separately
#'   from the bulk of the plot, see argument \code{whichClusters} of
#'   \code{\link{getClusterIndex}} for description of format allowed.
#' @param nBlankLines the number of blank (i.e. white) rows to add between the
#'   clusterMany clusterings and the highlighted clusterings.
#' @param nSizeResult the number of rows each highlighted clustering should take up.
#'   Increasing the number increases the thickness of the rectangles
#'   representing the highlighted clusterings.
#' @param clusterManyLabels either logical, indicating whether to plot the
#'   labels for the clusterings from clusterMany identified in the
#'   \code{whichClusterMany}, or a character vector of labels to use.
#' @param clusterLabels either logical, indicating whether to plot the labels
#'   for the clusterings identified to be highlighted in the
#'   \code{whichClusters} argument, or a character vector of labels to use.
#' @param sortBy how to align the clusters. If "highlighted" then the
#'   highlighted clusters indicated in the argument \code{whichClusters} are
#'   first in the alignment done by \code{plotClusters}. If "clusterMany", then
#'   the clusterMany results are first in the alignment. (Note this does not
#'   determine where they will be plotted, but how they are ordered in the
#'   aligning step done by \code{plotClusters})
#' @param highlightOnTop logical. Whether the highlighted clusters should be
#'   plotted on the top of clusterMany results or underneath.
#' @param existingColors one of "ignore","all","highlightOnly". Whether the plot
#'   should use the stored colors in the \code{ClusterExperiment} object given.
#'   "highlightOnly" means only the highlighted clusters will use the stored
#'   colors, not the clusterMany clusterings.
#' @param ... arguments passed to the matrix version of
#'     \code{\link{plotClusters}}
#' @details This plot is solely intended to make it easier to use the 
#'     \code{\link{plotClusters}} visualization when there are a large number of
#'     clusterings from a call to \code{\link{clusterMany}}. This plot separates
#'     out the \code{\link{clusterMany}} results from a designated clustering of
#'     interest, as indicated by the \code{whichClusters} argument
#'     (by default clusterings from a call to \code{\link{makeConsensus}} or 
#'     \code{\link{mergeClusters}}). In addition the highlighted clusters are
#'     made bigger so that they can be easily seen.
#' @seealso \code{\link{plotClusters}}, \code{\link{clusterMany}}
#' @return A plot is produced, nothing is returned.
#' @examples
#' #clustering using pam: try using different dimensions of pca and different k
#' data(simData)
#'
#' cl <- clusterMany(simData, nReducedDims=c(5, 10, 50), reduceMethod="PCA",
#' clusterFunction="pam", ks=2:4, findBestK=c(TRUE,FALSE),
#' removeSil=c(TRUE,FALSE))
#' cl <- makeConsensus(cl, proportion=0.7)
#' plotClustersWorkflow(cl)
#' @export
setMethod(
  f = "plotClustersWorkflow",
  signature = signature(object = "ClusterExperiment"),
  definition = function(object, whichClusters=c("mergeClusters","makeConsensus"), whichClusterMany=NULL, nBlankLines=ceiling(nClusterings(object)*.05), existingColors=c("ignore","all","highlightOnly"),
                        nSizeResult=ceiling(nClusterings(object)*.02), clusterLabels=TRUE, clusterManyLabels=TRUE, sortBy=c("highlighted","clusterMany"), highlightOnTop=TRUE,...)
  {
    sortBy<-match.arg(sortBy)
    existingColors<-match.arg(existingColors)
    allClusterMany<-which(clusterTypes(object)=="clusterMany")
    if("colData" %in% names(list(...))) stop("this function does not (yet) allow the designation of 'colData' argument. You must use plotClusters for this option.")
    if(is.null(whichClusterMany)){
      whichClusterMany<-allClusterMany
    }
    #convert to indices
	whichClusters<-getClusterIndex(object,whichClusters=whichClusters,noMatch="removeSilently")
    if(length(whichClusters)==0) stop("invalid identification of clusters for whichClusters argument")
	whichClusters<-getClusterIndex(object,whichClusters=whichClusterMany,noMatch="removeSilently")
	if(length(whichClusterMany)==0) stop("invalid identification of clusters for whichClusters argument")
    
    
    #result labels (yaxis):
    if(is.logical(clusterLabels)){
      if(clusterLabels) clusterLabels<-clusterLabels(object)[whichClusters]
      else clusterLabels<-rep("",length(whichClusters))  
    }
    else{
      if(length(clusterLabels)!=length(whichClusters) & !is.null(clusterLabels)){
        stop("number of cluster labels given in clusterLabels must be equal to the number of clusterings in 'whichClusters'")
        
      }
    }
    #clusterMany labels (yaxis):
    if(is.logical(clusterManyLabels)){
      if(clusterManyLabels) clusterManyLabels<-clusterLabels(object)[whichClusterMany]
      else clusterManyLabels<-rep("",length(whichClusterMany))
    }
    else{
      if(length(clusterManyLabels)!=length(whichClusterMany) & !is.null(clusterMany)){
        stop("number of cluster labels given in clusterManyLabels must be equal to the number of clusterings in 'whichClusterMany'")
        
      }
    }
    ###Get the sorted index using the matrix version of plotClusters
    ### out is the result of plotClusters
    if(sortBy=="highlighted"){
      orderOfClusters<-c(whichClusters,whichClusterMany)
      highIndex<-c(seq_along(whichClusters))
      cmIndex<-seq_along(orderOfClusters)[-highIndex]
    }
    else{
      orderOfClusters<-c(whichClusterMany,whichClusters)
      cmIndex<-c(seq_along(whichClusterMany))
      highIndex<-seq_along(orderOfClusters)[-cmIndex]
      
    }
    tempClusters<-clusterMatrix(object)[,orderOfClusters,drop=FALSE]
    out<-plotClusters(tempClusters,plot=FALSE)	 	
      
    ### Create color matrix
    ### resM is the highlighted clusters (columns the clusters)
    ### cmM is the clusterMany clusters (columns the clusters)
    
    if(existingColors!="ignore") 
      existingColorMat<-convertClusterLegend(object, whichClusters=orderOfClusters, output="matrixColors")
    
    if(existingColors %in% c("all","highlightOnly")){
      resM<-existingColorMat[,highIndex,drop=FALSE]
    }
    else{
      resM<-out$colors[,highIndex,drop=FALSE]
      
    }
    if(existingColors =="all"){
      cmM<-existingColorMat[,cmIndex,drop=FALSE]
    }
    else{
      cmM<-out$colors[,cmIndex,drop=FALSE] 
    }
    
    # make replication of results
    repResults<-lapply(seq_len(ncol(resM)),function(ii){
      x<-resM[,ii]
      mat<-matrix(x,nrow=length(x),ncol=nSizeResult,byrow=FALSE)
      colnames(mat)<-rep("",ncol(mat))
      colnames(mat)[ceiling(ncol(mat)/2)]<-clusterLabels[ii]
      return(mat)
    })
    repResults<-do.call("cbind",repResults)
    ##Add blanks
    if(highlightOnTop){
      bd<-makeBlankData(data=t(cbind(resM,cmM)),
						groupsOfFeatures=list("Results"=seq_along(whichClusters),"ClusterMany"=(length(whichClusters)+1):(length(whichClusters)+length(whichClusterMany))),
						nBlankFeatures=nBlankLines
						)
      whNotRes<-(length(whichClusters)+1):nrow(bd$dataWBlanks) #includes blanks
      whCM<-whNotRes[-c(seq_len(nBlankLines))] #no blanks
    } 	
    else{
      bd<-makeBlankData(data=t(cbind(cmM,resM)), groupsOfFeatures= list("ClusterMany"=seq_along(whichClusterMany), "Results"=(length(whichClusterMany)+1):(length(whichClusterMany)+length(whichClusters))),nBlankFeatures=nBlankLines)
      whNotRes<-  seq_len(length(whichClusterMany)+nBlankLines) #includes blanks
      whCM<-  seq_along(whichClusterMany) #no blanks
    }  
    test<-t(bd$dataWBlanks)
    test[is.na(test)]<-"white"
    
    ###ClusterLabels
    
    colnames(test)<-rep("",ncol(test))
    colnames(test)[whCM]<-clusterManyLabels
    
    if(highlightOnTop) test<-cbind(repResults,test[,whNotRes])
    else test<-cbind(test[,whNotRes],repResults)
    
    plotClusters(test[out$orderSamples,],input="colors",...)
    
  })