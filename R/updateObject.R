#' @title Update old ClusterExperiment object to current class definition
#'
#' @description This function updates ClusterExperiment objects from previous versions of package into the current definition
#' @param object a \code{ClusterExperiment} (or \code{clusterExperiment} from older 
#' versions). Must have at a minimum a slot \code{clusterMatrix}.
#' @inheritParams BiocGenerics::updateObject
#' @details The function creates a valid \code{ClusterExperiment} object by adding the default values of missing slots. It does so by calling the \code{\link{ClusterExperiment}} function, which imputs default (empty) values for missing slots.
#' @details The object is required to have minimal components to be updated. Specifically, it must have all the required elements of a Summarized Experiment as well as the basic slots of a ClusterExperiment object which have not changed over time. These are:   \code{clusterMatrix},\code{primaryIndex},\code{clusterInfo},\code{transformation}, \code{clusterTypes},\code{clusterLegend},\code{orderSamples}.
#' @return A valid \code{ClusterExperiment} object based on the current definition of 
#' ClusterExperiment.
#' @seealso \code{\link{ClusterExperiment}}
#' @export
#' @importFrom BiocGenerics updateObject
setMethod(
  f = "updateObject",
  signature = signature(object = "ClusterExperiment"),
  definition = function(object#, ..., verbose=FALSE
  ){
		#--------
		#check has at least the required slots of SummarizedExperiment class
		#--------
		if(!all(slotNames("SummarizedExperiment") %in% snames)){
			missSE<-which(!slotNames("SummarizedExperiment") %in% snames)
			stop("given object does not have the basic slots of SummarizedExperiment, cannot be updated (missing:",paste(slotNames("SummarizedExperiment")[missSE],collapse=","),"). To construct a ClusterExperiment object from its original components, use the function 'ClusterExperiment'")
		}
			

		#list names of all current required slots
		ceSlots<-slotNames(object)
		testSnames<-sapply(ceSlots,.hasSlot,object=object)
		snames<-ceSlots[testSnames]



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
	if(canCoerce(object,"SummarizedExperiment")) se<-as(object,"SummarizedExperiment")
	# if(canCoerce(object,"SingleCellExperiment")){
	# 	#if object was from before SCE requirement (2.0.0)
	# 	se<-as(object,"SingleCellExperiment")
	# }
	# else{
	# 	if(canCoerce(object,"SummarizedExperiment")) se<-as(object,"SummarizedExperiment")
	# 	else stop("cannot coerce object to SummarizedExperiment")
	# }

	#--------
	# Ignore slots that have to be in groups, with warnings
	#--------

	dendroSlots<-c("dendro_samples", "dendro_clusters",
  "dendro_index", "dendro_outbranch")

	mergeSlots<-c("merge_index",
"merge_dendrocluster_index",
"merge_method", "merge_demethod", "merge_cutoff",
"merge_nodeProp", "merge_nodeMerge")

	if(any(!dendroSlots %in% snames)& any(dendroSlots %in% snames)){
		warning("'object' does not contain ALL required slots saving the dendro-related information. Will remove all dendro AND merge related slots")
		snames<-snames[-which(snames %in% dendroSlots)]
		snames<-snames[-which(snames %in% mergeSlots)]
	}

	if(any(!mergeSlots %in% snames) & any(mergeSlots %in% snames)){
		warning("'object' does not contain ALL required slots saving the merge-related information.  Will remove all merge related slots")
		snames<-snames[-which(snames %in% mergeSlots)]
	}

	#--------
	# Get included slots of ones I create 
	#--------
	myslots<- c("transformation",
  	"primaryIndex", "clusterInfo",
  "clusterTypes", "dendro_samples", "dendro_clusters",
  "dendro_index", "dendro_outbranch", "coClustering",
  "clusterLegend", "orderSamples", "merge_index",
"merge_dendrocluster_index",
"merge_method", "merge_demethod", "merge_cutoff",
"merge_nodeProp", "merge_nodeMerge")
	snames<-snames[snames %in% myslots]


	object<-try(do.call("ClusterExperiment",c(list(object=se,clusters=object@clusterMatrix,checkTransformAndAssay=checkTransformAndAssay),attributes(object)[snames])),silent=TRUE)
	if(!inherits(object,"try-error")){
		object<-callNextMethod()
		return(object)
	} 
		else stop("Attempt to convert did not result in valid object. Here is the error from 'ClusterExperiment':\n",out)
}
)

