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
#'   clusterings to plot (if NULL, defaults to all)
#' @param whichResults which clusterings to use as the results.
#' @param nBlankLines the number of blank (i.e. white) rows to add between the
#'   clusterMany clusterings and the results
#' @param nSizeResult the number of rows each result clustering should take up
#' @param clusterManyLabels either logical, whether to plot the labels for the
#'   clusterings from clusterMany, or a character vector of labels to use
#' @param resultLabels either logical, whether to plot the labels for the
#'   clusterings identified in the results , or a character vector of labels to
#'   use.
#' @param sortBy how to align the clusters. If "results" then the results are in
#'   the top of the alignment done by plotClusters. If "clusterMany", then the
#'   clusterMany results are in the top. (Note this does not determine where
#'   they will be plotted, but how they are ordered in the aligning step done by
#'   \code{plotClusters})
#' @param resultsOnTop logical. Whether the results should be plotted on the top
#'   of clusterMany results or the bottom.
#' @seealso \code{\link{plotClusters}}, \code{\link{clusterMany}}
#' @export
setMethod(
  f = "plotClustersWorkflow",
  signature = signature(object = "ClusterExperiment"),
  definition = function(object, whichClusterMany=NULL, whichResults=c("mergeClusters","combineMany"),nBlankLines=ceiling(nClusters(object)*.05), 
  nSizeResult=ceiling(nClusters(object)*.02), clusterManyLabels=TRUE, resultLabels=TRUE, sortBy=c("results","clusterMany"), resultsOnTop=TRUE,...)
  {
	  sortBy<-match.arg(sortBy)
	  allClusterMany<-which(clusterTypes(object)=="clusterMany")
	 if(is.null(whichClusterMany)){
		 whichClusterMany<-allClusterMany
	 }
	 if(!is.numeric(whichClusterMany)) stop("'whichClusterMany' must give numeric indices of clusters of the ClusterExperiment object")
	 if(any(!whichClusterMany %in% allClusterMany)) stop("input to `whichClusterMany` must be indices to clusters of type 'clusterMany' ")
		 #convert to indices
 	if(is.character(whichResults)){
 		whichResults<- .TypeIntoIndices(object,whClusters=whichResults)
 		if(length(whichResults)==0) stop("invalid identification of clusters for whichResults argument")
 	}
	 #result labels:
       if(is.logical(resultLabels)){
 		  if(resultLabels) resultLabels<-clusterLabels(object)[whichResults]
		  else resultLabels<-rep("",length(whichResults))  
 	  }
 	  else{
 	       if(length(resultLabels)!=length(whichResults) & !is.null(resultLabels)){
 	   			stop("number of cluster labels given in resultLabels must be equal to the number of clusterings in 'whichResults'")

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


	 if(sortBy=="results"){
		 tempClusters<-clusterMatrix(object)[,c(whichResults,whichClusterMany),drop=FALSE]
		 out<-plotClusters(tempClusters,plot=FALSE)	 	
		 cmM<-out$colors[,-c(1:length(whichResults)),drop=FALSE]
		 resM<-out$colors[,c(1:length(whichResults)),drop=FALSE]
	 }
	 else{
		 tempClusters<-clusterMatrix(object)[,c(whichClusterMany,whichResults),drop=FALSE]
		 out<-plotClusters(tempClusters,plot=FALSE)
		 cmM<-out$colors[,c(1:length(whichClusterMany)),drop=FALSE]
		 resM<-out$colors[,-c(1:length(whichClusterMany)),drop=FALSE]
	 }
	 

	 
 	# make replication of results
 	 repResults<-lapply(1:ncol(resM),function(ii){
 		 x<-resM[,ii]
 		 mat<-matrix(x,nrow=length(x),ncol=nSizeResult,byrow=FALSE)
 	 	 colnames(mat)<-rep("",ncol(mat))
 		 colnames(mat)[ceiling(ncol(mat)/2)]<-resultLabels[ii]
 		 return(mat)
 	 })
 	 repResults<-do.call("cbind",repResults)
	 ##Add blanks
	 if(resultsOnTop){
    		bd<-makeBlankData(t(cbind(resM,cmM)), list("Results"=1:length(whichResults),"ClusterMany"=(length(whichResults)+1):(length(whichResults)+length(whichClusterMany))),nBlankLines=nBlankLines)
			whNotRes<-(length(whichResults)+1):nrow(bd$dataWBlanks) #includes blanks
			whCM<-whNotRes[-c(1:nBlankLines)] #no blanks
	 } 	
	 else{
   		bd<-makeBlankData(t(cbind(cmM,resM)), list("ClusterMany"=1:length(whichClusterMany), "Results"=(length(whichClusterMany)+1):(length(whichClusterMany)+length(whichResults))),nBlank=nBlankLines)
   	  	whNotRes<-  1:(length(whichClusterMany)+nBlankLines) #includes blanks
   	  	whCM<-  1:(length(whichClusterMany)) #no blanks
	 }  
	 test<-t(bd$dataWBlanks)
  	 test[is.na(test)]<-"white"
	 
	 ###ClusterLabels

	colnames(test)<-rep("",ncol(test))
	colnames(test)[whCM]<-clusterManyLabels
 
	if(resultsOnTop) test<-cbind(repResults,test[,whNotRes])
	else test<-cbind(test[,whNotRes],repResults)
	
	plotClusters(test[out$orderSamples,],input="colors",...)

})