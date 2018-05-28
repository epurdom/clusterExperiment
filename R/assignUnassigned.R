#' @title Assign unassigned samples to nearest cluster
#' @description Assigns the unassigned samples in a cluster to the nearest
#'   cluster based on distance to the medians of the clusters.
#' @name assignUnassigned
#' @rdname assignUnassigned
#' @aliases assignUnassigned,ClusterExperiment-method
#' @param object A Cluster Experiment object
#' @param whichCluster in which clustering to get cluster assignments and assign
#'   the unassigned samples
#' @param whichAssay which assay to use to calculate the median per cluster and
#'   take dimensionality reduction (if requested)
#' @param clusterLabel if missing, the current cluster label of the cluster will
#'   be appended with the string "_AllAssigned".
#' @param ... arguments passed to \code{\link{getReducedData}} specifying the
#'   dimensionality reduction (if any) to be taken of the data for calculating
#'   the medians of the clusters
#' @details The function \code{assignUnassigned} calculates the median values of each variable for each
#'   cluster, and then calculates the euclidean distance of each unassigned
#'   sample to the median of each cluster. Each unassigned sample is assigned to
#'   the cluster for which it closest to the median.
#' @details All unassigned samples in the cluster are given a clustering,
#'   regardless of whether they are classified as -1 or -2.
#' @return The function \code{assignUnassigned} returns a \code{ClusterExperiment}
#' object with the unassigned samples assigned to one of the existing clusters. 
#' @seealso \code{\link{getReducedData}}
#' @examples
#' #load CE object
#' data(rsecFluidigm)
#' smallCE<-rsecFluidigm[,1:50]
#' #assign the unassigned samples
#' assignUnassigned(smallCE, makePrimary=TRUE)
#' 
#' #note how samples are REMOVED:
#' removeUnassigned(smallCE)
#' @inheritParams addClusterings 
#' @inheritParams reduceFunctions
#' @export
setMethod(
  f = "assignUnassigned",
  signature = signature("ClusterExperiment"),
  definition = function(object,whichCluster="primary",clusterLabel,
                        makePrimary=TRUE,whichAssay=1,reduceMethod="none",...){
		
    whCl<-.convertSingleWhichCluster(object,whichCluster,list(...))
    cl<-clusterMatrix(object)[,whCl]
		if(missing(clusterLabel)) clusterLabel<-paste0(clusterLabels(object)[whCl],"_AllAssigned")
		whichUnassigned<-which(cl<0)
		if(length(whichUnassigned)>0){
				if(length(whichUnassigned)< length(cl)){
			    ########
			    ##Transform the data
			    ########
			    datList<-getReducedData(object,reduceMethod=reduceMethod,returnValue="list",...)
			    object<-datList$objectUpdate
			    dat<-datList$dat
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
					object<-addClusterings(object,cl, clusterLabel= clusterLabel,makePrimary=makePrimary)
					return(object)
			}
			else{
				stop("All cells are unassigned, cannot assign them")
			}
		}
		else stop("No cells are unassigned in the designated cluster")

  }
)


#' @rdname assignUnassigned
#' @aliases removeUnassigned
#' @details \code{removeUnclustered} removes all samples that are unclustered
#'   (i.e. -1 or -2 assignment) in the designated cluster of \code{object} (so they may
#'   be unclustered in other clusters found in \code{clusterMatrix(object)}).
#' @return The function \code{removeUnassigned} returns a \code{ClusterExperiment}
#' object with the unassigned samples removed. 
#' @export
setMethod(
  f = "removeUnassigned",
  signature = "ClusterExperiment",
  definition = function(object,whichCluster="primary") {
    whCl<-.convertSingleWhichCluster(object,whichCluster)
		cl<-clusterMatrix(object)[,whCl]
		return(object[,cl>= 0])
  }
)
