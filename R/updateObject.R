#' @title Update old ClusterExperiment object to current class definition
#'
#' @description This function updates ClusterExperiment objects from previous versions of package into the current definition
#' @param object a \code{ClusterExperiment} (or \code{clusterExperiment} from older 
#' versions). Must have at a minimum a slot \code{clusterMatrix}.
#' @inheritParams ClusterExperiment-class
#' @details The function creates a valid \code{ClusterExperiment} object by adding the default values of missing slots. It does so by calling the \code{\link{ClusterExperiment}} function, which imputs default (empty) values for missing slots.
#' @details The object is required to have minimal components to be updated. Specifically, it must have all the required elements of a Summarized Experiment as well as the basic slots of a ClusterExperiment object which have not changed over time. These are:   \code{clusterMatrix},\code{primaryIndex},\code{clusterInfo},\code{transformation}, \code{clusterTypes},\code{clusterLegend},\code{orderSamples}.
#' @return A valid \code{ClusterExperiment} object based on the current definition of 
#' ClusterExperiment.
#' @seealso \code{\link{ClusterExperiment}}
#' @export
updateObject<-function(object,checkTransformAndAssay=FALSE) {
	snames<-slotNames(object)
		#list names of all current required slots
		ceSlots<-getSlots("ClusterExperiment")
		#id the ones I create
		myslots<- c("transformation",
    	"primaryIndex", "clusterInfo",
    "clusterTypes", "dendro_samples", "dendro_clusters",
    "dendro_index", "dendro_outbranch", "coClustering",
    "clusterLegend", "orderSamples", "merge_index",
	"merge_dendrocluster_index",
	"merge_method", "merge_demethod", "merge_cutoff",
	"merge_nodeProp", "merge_nodeMerge")
	
	#--------
	#check has at least the required slots of SummarizedExperiment class
	#--------
	if(!all(names(getSlots("SummarizedExperiment")) %in% snames)){
		missSE<-which(!names(getSlots("SummarizedExperiment")) %in% snames)
		stop("given object does not have the basic slots of SummarizedExperiment, cannot be updated (missing:",paste(names(getSlots("SummarizedExperiment"))[missSE],collapse=","),"). To construct a ClusterExperiment object from its original components, use the function 'ClusterExperiment'")
	}
	
	#--------
	#check has minimal required slots of a clusterExperiment object of any version
	#--------
  requiredSlots<-c("clusterMatrix","primaryIndex", "clusterInfo","transformation","clusterTypes","clusterLegend","orderSamples")
	if(!all(requiredSlots %in% snames)){
		missCE<-which(!requiredSlots %in% snames)
		stop("given object does not have the basic slots of ClusterExperiment, cannot be updated (missing:",paste(requiredSlots[missCE],collapse=","),"). To construct a ClusterExperiment object from its original components, use the function 'ClusterExperiment'")
	}

	#--------
	#extract either SE or SCE object
	#--------
	if(!all(names(getSlots("SingleCellExperiment")) %in% snames)){
		#if object was from before SCE requirement (2.0.0)
		se<-as(object,"SummarizedExperiment")
	}
	else{
		se<-as(object,"SingleCellExperiment")
	}
	
	#--------
	# Get included slots into list
	#--------
	objectList<-lapply(myslots,function(slotName){
		if(slotName %in% snames) return(object@slotName)
			else return(NULL)
	})
	return(do.call("ClusterExperiment",c(list(object=se,clusters=object@clusterMatrix,checkTransformAndAssay=checkTransformAndAssay),objectList)))
}


