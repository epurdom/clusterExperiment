#' @rdname plotContrastHeatmap
#' @title Plot heatmaps showing significant genes per contrast
#' @description Plots a heatmap of the data, with the genes grouped based on the contrast for which they were significant. 
#' @param object ClusterExperiment object on which biomarkers were found
#' @param signifTable A \code{data.frame} in format of the result of \code{\link{getBestFeatures}}. It must minimally contain columns 'Contrast' and 'IndexInOriginal' giving the grouping and original index of the features in the \code{assay(object)}
#' @param nBlankLines integer:  many blank lines to add between contrast groups.
#' @param ... Arguments passed to \code{\link{plotHeatmap}}
#' @seealso \code{\link{plotHeatmap}}, \code{\link{makeBlankData}}, \code{\link{getBestFeatures}}
#' @export
#'
#' @examples
#' data(simData)
#'
#' cl <- clusterSingle(simData, subsample=FALSE,
#' sequential=FALSE, mainClusterArgs=list(clusterFunction="pam", clusterArgs=list(k=8)))
#'
#' #Do all pairwise, only return significant, try different adjustments:
#' pairsPerC <- getBestFeatures(cl, contrastType="Pairs", number=5,
#' p.value=0.05, isCount=FALSE)
#' plotContrastHeatmap(cl,pairsPerC)
setMethod(
  f = "plotContrastHeatmap",
  signature = "ClusterExperiment",
  definition = function(object,signifTable,...) {
	  if(!all(c("IndexInOriginal","Contrast") %in% colnames(signifTable ))) stop("signifTable must have columns 'IndexInOriginal' and 'Contrast'")
	if(!is.numeric(signifTable$IndexInOriginal)) stop("Column 'IndexInOriginal' Must consist of numeric values")
	 if(!all(signifTable$IndexInOriginal %in% 1:nrow(object))) stop("Column 'IndexInOriginal' must consist of indices that match the row indices of 'object'")
	 	  #divide by contrast
	  geneByContrast<-by(signifTable,signifTable$Contrast,function(x){x$IndexInOriginal})
	  if("ContrastName" %in% colnames(signifTable)){
		  gpnames<-unique(signifTable[,c("Contrast","ContrastName")])
		  if(nrow(gpnames)==length(geneByContrast)){
			  m<-match(names(geneByContrast),gpnames[,"Contrast"])
			  names(geneByContrast)<-gpnames[m,"ContrastName"]
		  }
	  }
	 plotHeatmap(object,clusterFeaturesData=geneByContrast,...)

}
)