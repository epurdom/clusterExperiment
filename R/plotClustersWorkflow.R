#' @rdname plotClustersWorkflow
#' @aliases plotClustersWorkflow
#' @title A plot of clusterings specific for clusterMany and workflow
#'   visualization
#' @description A realization of \code{\link{plotClusters}} call specific to
#'   separating out the results of \code{clusterMany} and other clustering
#'   results.
#' @param object A \code{ClusterExperiment} object on which
#'   \code{\link{clusterMany}} has been run
#' @param whichClusterMany numeric indices of which of the clusterMany 
#'   clusterings to plot (if NULL, defaults to all). Unlike
#'   \code{whichClusters}, these must be numeric indices. They must also refer
#'   to clusterings of clusterType \code{clusterMany}.
#' @param whichClusters which clusterings to "highlight", i.e draw separately 
#'   from the \code{clusterMany} results. Can be numeric or character vector, 
#'   indicating the indices or clusterLabels/clusterTypes of the clusterings of 
#'   interest, respectively.
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
#' @param existingColors logical. If logical, whether the highlighted clusters
#'   should use colors matched to the clusterMany results, or should the stored
#'   colors in \code{ClusterExperiment} object be used. This argument has no
#'   effect on the colors of the clusterMany results, whose colors will be
#'   chosen based on the alignment of plotClusters.
#'   @param ... arguments passed to the matrix version of
#'     \code{\link{plotClusters}}
#'   @details This plot is solely intended to make it easier to use the 
#'     \code{\link{plotClusters}} visualization when there are a large number of
#'     clusterings from a call to \code{\link{clusterMany}}. This plot separates
#'     out the \code{\link{clusterMany}} results from a designated clustering of
#'     interest, as indicated by the \code{whichClusters} argument
#'     (by default clusterings from a call to \code{\link{combineMany}} or 
#'     \code{\link{mergeClusters}}). In addition the highlighted clusters are
#'     made bigger so that they can be easily seen.
#' @seealso \code{\link{plotClusters}}, \code{\link{clusterMany}}
#' @export
setMethod(
  f = "plotClustersWorkflow",
  signature = signature(object = "ClusterExperiment"),
  definition = function(object, whichClusters=c("mergeClusters","combineMany"), whichClusterMany=NULL, nBlankLines=ceiling(nClusters(object)*.05), existingColors=FALSE,
  nSizeResult=ceiling(nClusters(object)*.02), clusterLabels=TRUE, clusterManyLabels=TRUE, sortBy=c("highlighted","clusterMany"), highlightOnTop=TRUE,...)
  {
	  sortBy<-match.arg(sortBy)
	  allClusterMany<-which(clusterTypes(object)=="clusterMany")
	  if("sampleData" %in% names(list(...))) stop("this function does not (yet) allow the designation of 'sampleData' argument. You must use plotClusters for this option.")
	 if(is.null(whichClusterMany)){
		 whichClusterMany<-allClusterMany
	 }
	 if(!is.numeric(whichClusterMany)) stop("'whichClusterMany' must give numeric indices of clusters of the ClusterExperiment object")
	 if(any(!whichClusterMany %in% allClusterMany)) stop("input to `whichClusterMany` must be indices to clusters of type 'clusterMany' ")
		 #convert to indices
 	if(is.character(whichClusters)){
 		whichClusters<- .TypeIntoIndices(object,whClusters=whichClusters)
 		if(length(whichClusters)==0) stop("invalid identification of clusters for whichClusters argument")
 	}
	 #result labels:
       if(is.logical(clusterLabels)){
 		  if(clusterLabels) clusterLabels<-clusterLabels(object)[whichClusters]
		  else clusterLabels<-rep("",length(whichClusters))  
 	  }
 	  else{
 	       if(length(clusterLabels)!=length(whichClusters) & !is.null(clusterLabels)){
 	   			stop("number of cluster labels given in clusterLabels must be equal to the number of clusterings in 'whichClusters'")

 	   		}
 	  }
 	 #clusterMany labels:
        if(is.logical(clusterManyLabels)){
  		  if(clusterManyLabels) clusterManyLabels<-clusterLabels(object)[whichClusterMany]
			  else clusterManyLabels<-rep("",length(whichClusterMany))
  	  }
  	  else{
  	       if(length(clusterManyLabels)!=length(whichClusterMany) & !is.null(clusterMany)){
  	   			stop("number of cluster labels given in clusterManyLabels must be equal to the number of clusterings in 'whichClusterMany'")

  	   		}
  	  }


	 if(sortBy=="highlighted"){
		 tempClusters<-clusterMatrix(object)[,c(whichClusters,whichClusterMany),drop=FALSE]
		 out<-plotClusters(tempClusters,plot=FALSE)	 	
		 cmM<-out$colors[,-c(1:length(whichClusters)),drop=FALSE]
		 resM<-out$colors[,c(1:length(whichClusters)),drop=FALSE]
	 }
	 else{
		 tempClusters<-clusterMatrix(object)[,c(whichClusterMany,whichClusters),drop=FALSE]
		 out<-plotClusters(tempClusters,plot=FALSE)
		 cmM<-out$colors[,c(1:length(whichClusterMany)),drop=FALSE]
		 resM<-out$colors[,-c(1:length(whichClusterMany)),drop=FALSE]
	 }
	 
	 if(existingColors){
		 #don't use resM for the colors, but make one
		 resMatch<-lapply(1:length(whichClusters),function(ii){
			 x<-clusterMatrix(object)[,ii]
			 col<-clusterLegend(object)[[ii]]
			 m<-match(as.character(x),col[,"clusterIds"])
			 return(col[m,"color"])
		 })
		 resM<-do.call("cbind",resMatch)
	 }
	 
 	# make replication of results
 	 repResults<-lapply(1:ncol(resM),function(ii){
 		 x<-resM[,ii]
 		 mat<-matrix(x,nrow=length(x),ncol=nSizeResult,byrow=FALSE)
 	 	 colnames(mat)<-rep("",ncol(mat))
 		 colnames(mat)[ceiling(ncol(mat)/2)]<-clusterLabels[ii]
 		 return(mat)
 	 })
 	 repResults<-do.call("cbind",repResults)
	 ##Add blanks
	 if(highlightOnTop){
    		bd<-makeBlankData(t(cbind(resM,cmM)), list("Results"=1:length(whichClusters),"ClusterMany"=(length(whichClusters)+1):(length(whichClusters)+length(whichClusterMany))),nBlankLines=nBlankLines)
			whNotRes<-(length(whichClusters)+1):nrow(bd$dataWBlanks) #includes blanks
			whCM<-whNotRes[-c(1:nBlankLines)] #no blanks
	 } 	
	 else{
   		bd<-makeBlankData(t(cbind(cmM,resM)), list("ClusterMany"=1:length(whichClusterMany), "Results"=(length(whichClusterMany)+1):(length(whichClusterMany)+length(whichClusters))),nBlank=nBlankLines)
   	  	whNotRes<-  1:(length(whichClusterMany)+nBlankLines) #includes blanks
   	  	whCM<-  1:(length(whichClusterMany)) #no blanks
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