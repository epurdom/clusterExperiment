#' @name plotFeatureScatter
#' @title Plot scatter plot of feature values colored by cluster
#' @description Plot a scatter plot of the (transformed) values for a set of
#'   gene expression values, colored by cluster
#' @aliases plotFeatureScatter plotFeatureScatter,ClusterExperiment,character-method
#' @inheritParams plotFeatureBoxplot
#' @rdname plotFeatureScatter
#' @examples
#' data(simData)
#' #Create a ClusterExperiment object
#' cl <- clusterMany(simData, nReducedDims=c(5, 10, 50), reducedDim="PCA",
#'    clusterFunction="pam", ks=2:4, findBestK=c(TRUE,FALSE),
#'    removeSil=c(TRUE,FALSE))
#' #give names to the clusters
#' cl<-renameClusters(cl, whichCluster=1, 
#'    value=letters[1:nClusters(cl)[1]])
#' plotFeatureScatter(cl,feature=1:2,pch=19,legendLocation="topright")
#' @export
setMethod(
  f = "plotFeatureScatter",
  signature = signature(object = "ClusterExperiment",features="character"),
  definition = function(object, features,...)
  {
    m<-match(features,rownames(object))
    if(any(is.na(m))) stop("not all of features match one of the rownames of this object")
    else invisible(plotFeatureScatter(object,features=m,...))
  })

#' @param legendLocation character value passed to \code{location} argument of 
#'  \code{plotClusterLegend} indicating where to put the legend. If NA, legend 
#'  will not be plotted.
#' @param features the indices of the features (either numeric or character 
#'  matching rownames of object) to be plotted.
#' @param jitterFactor numeric. If NA, no jittering is done. Otherwise, passed to 
#' \code{factor} of function \code{\link[base]{jitter}}
#'  (useful for low counts)
#' @rdname plotFeatureScatter
#' @export
setMethod(
  f = "plotFeatureScatter",
  signature = signature(object = "ClusterExperiment",features="numeric"),
  definition = function(object, features, whichCluster="primary", plotUnassigned=TRUE,unassignedColor="grey",missingColor="white",whichAssay=1,legendLocation=NA,jitterFactor=NA,...)
  {
		if(length(features)<2) stop("plotFeatureScatter requires at least 2 features")
    whCl<-getSingleClusterIndex(object,whichCluster,list(...))
 
    #get data:
    dat<-transformData(object, whichAssay=whichAssay)[features,]
		if(!is.na(jitterFactor)) dat<-jitter(dat,factor=jitterFactor)
		if(is.null(rownames(dat))){
			rownames(dat)<-paste("Feature",features,sep=" ")
		}
    clLegend<-clusterLegend(object)[[whCl]]
		cl<-clusterMatrix(object)[,whCl]
		whLegUn<-which(as.numeric(clLegend[,"clusterIds"])<0)
		if(!plotUnassigned){
			whUnassign<-which(cl<0)
			if(length(whUnassign)>0){
				cl<-cl[-whUnassign]
				dat<-dat[,-whUnassign]
				clLegend<-clLegend[-whLegUn,,drop=FALSE]

			}	
		}
		else{
			if(length(whLegUn)>0){
				wh1<-which(as.numeric(clLegend[,"clusterIds"])== -1)
				if(length(wh1)>0 & !is.null(unassignedColor))
					clLegend[wh1,"color"]<-unassignedColor
				wh2<-which(as.numeric(clLegend[,"clusterIds"])== -2)
				if(length(wh2)>0 & !is.null(missingColor))
					clLegend[wh2,"color"]<-missingColor
				
			}
		}
    uniqueNames<-length(unique(clLegend[,"name"]))==nrow(clLegend)


		m<-match(as.character(cl),clLegend[,"clusterIds"])
		cols<-clLegend[m,"color"]
		if(length(features)==2){
			invisible(plot(t(dat),col=cols,bg=cols,...))
			if(!is.na(legendLocation)) 
				legend(x=legendLocation,legend=clLegend[,"name"],fill=clLegend[,"color"])
		}
		else{
			invisible(pairs(t(data.matrix(dat)),col=cols,bg=cols,...))
		}
  })