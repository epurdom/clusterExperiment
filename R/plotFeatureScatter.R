#' @name plotFeatureScatter
#' @title Plot scatter plot of feature values colored by cluster
#' @description Plot a scatter plot of the (transformed) values for a set of
#'   gene expression values, colored by cluster
#' @aliases plotFeatureScatter plotFeatureScatter,ClusterExperiment-method
#' @inheritParams plotFeatureBoxplot
#' @rdname plotFeatureScatter
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
setMethod(
  f = "plotFeatureScatter",
  signature = signature(object = "ClusterExperiment",features="numeric"),
  definition = function(object, features, whichCluster, plotUnassigned=TRUE,unassignedColor="grey",missingColor="white",whichAssay=1,legendLocation=NA,...)
  {
		if(length(features)<2) stop("plotFeatureScatter requires at least 2 features")
    whCl<-.convertSingleWhichCluster(object,whichCluster,list(...))
 
    #get data:
    dat<-transformData(object, whichAssay=whichAssay)[features,]
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