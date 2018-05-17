#' @title Assign unassigned samples to nearest cluster
#' @description Assigns the unassigned samples in a cluster to the nearest cluster based on distance to the medians of the clusters.
#' @name assignUnassigned
#' @rdname assignUnassigned
#' @aliases assignUnassigned,ClusterExperiment-method
#' @param whichCluster in which clustering to assign the unassigned samples
#' @param whichAssay which assay to use to calculate the median per cluster
#' @param reduceMethod which reduction method to apply to the data before taking the median per cluster
#' @param clusterLabel if missing, the current cluster label of the cluster will be appended with the string "_AllAssigned".
#' @param nDims the number of dimensions to keep if using a reduceMethod
#' @details The function calculates the median values of each variable for each cluster, and then calculates the euclidean distance of each unassigned sample to the median of each cluster. Each unassigned sample is assigned to the cluster for which it closest to the median.
#' @details All unassigned samples in the cluster are given a clustering, regardless of whether they are classified as -1 or -2.
#' @inheritParams addClusterings
#' @export
setMethod(
  f = "assignUnassigned",
  signature = signature("ClusterExperiment"),
  definition = function(x,whichCluster="primary",clusterLabel,makePrimary=TRUE,whichAssay=1,reduceMethod="none",nDims=defaultNDims(x,reduceMethod)){
		
		if(is.character(whichCluster)) whCl<-.TypeIntoIndices(x,whClusters=whichCluster) else whCl<-whichCluster
    if(length(whCl)!=1) stop("Invalid value for 'whichCluster'. Current value identifies ",length(whCl)," clusterings, but 'whichCluster' must identify only a single clustering.")
    if(!whCl %in% seq_len(nClusterings(x))) stop("Invalid value for 'whichCluster'. Must be integer between 1 and ", nClusterings(x))
    cl<-clusterMatrix(x)[,whCl]
		if(missing(clusterLabel)) clusterLabel<-paste0(clusterLabels(x)[whCl],"_AllAssigned")
		whichUnassigned<-which(cl<0)
		if(length(whichUnassigned)>0){
				if(length(whichUnassigned)< length(cl)){
			    ########
			    ##Transform the data
			    ########
			    if(length(reduceMethod)>1) stop('assignUnassigned only takes one choice of "reduceMethod" as argument')
			    if(length(nDims) > 1) {
			      stop("makeDendrogram only handles one choice of dimensions.")
			    }
			    if(!is.na(nDims) & reduceMethod=="none") {
			      warning("specifying nDims has no effect if reduceMethod==`none`")
			    }
			    ###Calculate filters/reduceMethod if needed...
			    if(!isReducedDims(x,reduceMethod) & isBuiltInReducedDims(reduceMethod)){
			      x<-makeReducedDims(x,reducedDims=reduceMethod,maxDims=nDims,whichAssay=whichAssay)
			    }
			    else if(!isFilterStats(x,reduceMethod) & isBuiltInFilterStats(reduceMethod)){
			      x<-makeFilterStats(x,filterStat=reduceMethod, whichAssay=whichAssay,
			                         whichClusterIgnoreUnassigned=if(ignoreUnassignedVar) whCl else NULL)
      
			    }
					if(reduceMethod=="none")
			      dat<-transformData(x, whichAssay=whichAssay)
			    else if(isReducedDims(x,reduceMethod))
			      dat<-t(reducedDim(x,type=reduceMethod)[,seq_len(nDims)])
			    else if(isFilterStats(x,reduceMethod))
			      dat<-transformData(filterData(x,filterStats=reduceMethod,percentile=nDims), whichAssay=whichAssay)
			    else
			      stop("'x' does not contain the given 'reduceMethod' value nor does 'reduceMethod' value match any built-in filters or dimensionality reduction options.")

					###############
					#find centers of clusters based on assigned samples:
					###############
			    clFactor <- factor(cl[-whichUnassigned])
					medoids <- do.call("rbind", by(t(dat[,-whichUnassigned]), clFactor, function(z){apply(z, 2, median)}))
			    rownames(medoids) <- levels(clFactor)
			
					classif<-.genericClassify(x=dat[,whichUnassigned],centers=medoids)
			
					###############
					#check reasonable results:
					###############
					if(any(!unique(classif)%in%unique(cl[cl>0]))) stop("programming error in classifying -- clusters not seen before given as assignment")
					if(	length(which(cl<0))!=length(classif)) stop("programming error in classifying -- not right length")
					###############
					# Assign cluster to object
					###############
					cl[cl<0]<-classif
					x<-addClusterings(x,cl, clusterLabel= clusterLabel,makePrimary=makePrimary)
					return(x)
			}
			else{
				stop("All cells are unassigned, cannot assign them")
			}
		}
		else stop("No cells are unassigned in the designated cluster")
		
		


  }
)


