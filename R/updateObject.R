#' @title Update old ClusterExperiment object to current class definition
#' @name updateObject
#' @description This function updates ClusterExperiment objects from previous
#'   versions of package into the current definition
#' @param object a \code{ClusterExperiment} (or \code{clusterExperiment} from 
#'   older versions). Must have at a minimum a slot \code{clusterMatrix}.
#' @inheritParams BiocGenerics::updateObject 
#' @inheritParams ClusterExperiment-class
#' @details The function creates a valid \code{ClusterExperiment} object by 
#'   adding the default values of missing slots. It does so by calling the 
#'   \code{\link{ClusterExperiment}} function, which imputs default (empty) 
#'   values for missing slots.
#' @details The object is required to have minimal components to be updated. 
#'   Specifically, it must have all the required elements of a Summarized 
#'   Experiment as well as the basic slots of a ClusterExperiment object which 
#'   have not changed over time. These are: \code{clusterMatrix},
#'   \code{primaryIndex}, \code{clusterInfo}, \code{transformation}, 
#'   \code{clusterTypes}, \code{clusterLegend}, \code{orderSamples}.
#' @details If \emph{any} of the dendrogram-related slots are missing, ALL of 
#'   the dendrogram \emph{and} merge related slots will be cleared to default 
#'   values. Similarly, if \emph{any} of the merge-related slots are missing, 
#'   ALL of the merge-related slots will be cleared to the default values.
#' @details The function currently only works for object of 
#'   \code{ClusterExperiment}, not the older name \code{clusterExperiment}.
#' @return A valid \code{ClusterExperiment} object based on the current 
#'   definition of ClusterExperiment.
#' @seealso \code{\link{ClusterExperiment}}
#' @aliases updateObject,ClusterExperiment-method
#' @rdname updateObject
#' @export
#' @importFrom BiocGenerics updateObject
setMethod(
  f = "updateObject",
  signature = signature(object = "ClusterExperiment"),
  definition = function(object, checkTransformAndAssay=FALSE,...,verbose=FALSE){
		#create snames, which is the slots the object actually has
		#and will eventually be narrowed down to only those slots will pass to `ClusterExperiment`
		#list names of all current required slots
		ceSlots<-slotNames(object)
		testSnames<-sapply(ceSlots,.hasSlot,object=object)
		snames<-ceSlots[testSnames]

		#--------
		#check has at least the required slots of SummarizedExperiment class
		#--------
		if(!all(slotNames("SummarizedExperiment") %in% snames)){
			missSE<-which(!slotNames("SummarizedExperiment") %in% snames)
			stop("given object does not have the basic slots of SummarizedExperiment, cannot be updated (missing:",paste(slotNames("SummarizedExperiment")[missSE],collapse=","),"). To construct a ClusterExperiment object from its original components, use the function 'ClusterExperiment'")
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
	# if(canCoerce(object,"SummarizedExperiment")) se<-updateObject(as(object,"SummarizedExperiment"))
	if(canCoerce(object,"SingleCellExperiment")){
		#if object was from before SCE requirement (2.0.0)
		se<-updateObject(as(object,"SingleCellExperiment"))
	}
	else{
		if(canCoerce(object,"SummarizedExperiment")) se<-updateObject(as(object,"SummarizedExperiment"))
		else stop("cannot coerce object to SummarizedExperiment")
	}

	#--------
	# Ignore slots that have to come together, with warnings
	#--------

	dendroSlots<-c("dendro_samples", "dendro_clusters",
  "dendro_index", "dendro_outbranch")

	mergeSlots<-c("merge_index",
"merge_dendrocluster_index",
"merge_method", "merge_demethod", "merge_cutoff",
"merge_nodeProp", "merge_nodeMerge")

	if(any(!dendroSlots %in% snames)& any(dendroSlots %in% snames)){
		warning("'object' does not contain ALL required slots saving the dendro-related information. Updated object will remove all dendro AND merge related slots")
		snames<-snames[-which(snames %in% dendroSlots)]
		snames<-snames[-which(snames %in% mergeSlots)]
	}

	if(any(!mergeSlots %in% snames) & any(mergeSlots %in% snames)){
		warning("'object' does not contain ALL required slots saving the merge-related information.  Updated object will remove all merge related slots")
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

	## Fix class of the dendrograms:
	if("dendro_samples" %in% snames){
		if(class(object@dendro_clusters)=="dendrogram"){
			newPhylo<-.makePhylobaseTree(object@dendro_clusters,isSamples=FALSE)
			#Add necessary information: Position, NodeId, ClusterIdDendro, ClusterIdMerge
			
			#make tip labels match cluster names
		}
		if(class(object@dendro_samples)=="dendrogram"){
			newPhylo<-.makePhylobaseTree(object@dendro_samples,isSamples=TRUE,outbranch=object@dendro_outbranch)
			#Add necessary information: Position, NodeId, SampleIndex
			
			
			#make node labels match cluster names and sample names
		}
	}

	object<-try(do.call("ClusterExperiment",c(list(object=se,clusters=object@clusterMatrix,checkTransformAndAssay=checkTransformAndAssay),attributes(object)[snames])),silent=TRUE)
	if(!inherits(object,"try-error")){
		return(object)
	} 
		else stop("Attempt to convert did not result in valid object. Here is the error from 'ClusterExperiment':\n",object)
}
)




####
#Convert dendrogram class slots to class used by phylobase (phylo4) so can navigate easily. Does so by first converting to class of ape (phylo)
#' @importFrom phylobase edgeLength rootNode descendants nodeLabels
#' @importFrom ape as.phylo
#' @importFrom stats as.hclust
#' @importClassesFrom phylobase phylo4 
.makePhylobaseTree<-function(x,isSamples=FALSE,outbranch=FALSE, returnOnlyPhylo=FALSE){
  type<-class(x)
  if(type=="dendrogram"){
	  x<-try(stats::as.hclust(x))
	  if(inherits(x, "try-error")) stop("the dendrogram object cannot be converted to a hclust object class with the methods of the 'stats' package. Reported error from stats package:",x)
  }
  type<-class(x)
  if(type=="hclust"){
    #first into phylo from ape package
    tempPhylo<-try(ape::as.phylo(x),FALSE)
    if(inherits(tempPhylo, "try-error")) stop("the hclust object cannot be converted to a phylo class with the methods of the 'ape' package. Reported error from ape package:",tempPhylo)
  }
  if(type=="phylo") tempPhylo<-x

	  #put this in because some zero edges are becoming very small negative values...
  if(any(tempPhylo$edge.length<0)){
	  if(all(tempPhylo$edge.length> -1e-5)) tempPhylo$edge.length[tempPhylo$edge.length<0]<-0
	  else stop("coding error -- given object results in negative edge lengths")
  }
  if(returnOnlyPhylo ) return(tempPhylo) #don't do the rest of fixing up...

  phylo4Obj<-try(as(tempPhylo,"phylo4"),FALSE) 
  if(inherits(phylo4Obj, "try-error")) stop(paste("the internally created phylo object cannot be converted to a phylo4 class. Check that you gave simple hierarchy of clusters, and not one with fake data per sample. Reported error from phylobase package:",tempPhylo))
  
  if(isSamples){
	  ###Adds (in uniform way) the node names:
	  ### "Root"
	  ### "NodeX" where X is a number *to those NOT in outgroup*
	  ### "MissingNodeX" where X is a number *to those IN outgroup*
    #NOTE: clusterNodes are determined by those with descendants with zero edge lengths, but non-zero edge-length between them and their decendents (i.e. fake hierarchy of the cluster must have edge length zero to root)
	#Moreover, they are distinguished from the -1/-2 outgroup by the fact that the tips of the -1/-2 DO HAVE non-zero edges to them. 
	#WHAT IF HAVE CLUSTER WITH ONE SAMPLE???? (just doesn't happen in practice, but often in our tests...)

	whNonZeroEdges<-which(sapply(phylobase::edgeLength(phylo4Obj),function(x){!isTRUE(all.equal(x,0)) & !is.na(x)}))
    nonZeroEdges<-phylobase::edgeLength(phylo4Obj)[ whNonZeroEdges ] #doesn't include root
    
	
	#all nodes where edge going *into* node is >0 -- excludes the root
	#this also picks up the outbranch between -1,-2 and all the internal nodes/tips there
	trueInternal<-sort(unique(as.numeric(sapply(strsplit(names(nonZeroEdges),"-"),.subset2,2)))) 
	#add root to trueInternal
    rootNode<-phylobase::rootNode(phylo4Obj)
    trueInternal<-c(rootNode,trueInternal)
	
	#make sure no tips included (single sample clusters)
	trueInternal<-trueInternal[!trueInternal%in%phylobase::nodeId(phylo4Obj,"tip")]
	
    if(outbranch){ 
      #######
      #remove root from labeling schema if there exists -1 outbranch
      #######
      trueInternal<-trueInternal[!trueInternal%in%rootNode]
      
      #######
      #find the -1/-2 internal node (if it exists)
      #determine it as the one without 0-length tip edges.
      #assumes all tips in the non-outbranch have 0-length (so max value is zero)
      #######
      rootChild<-phylobase::descendants(phylo4Obj,node=rootNode,type="children")
      #find tip descendants of each of these:
      rootChildDesc<-lapply(rootChild,phylobase::descendants,phy=phylo4Obj,type="tip")
      rootChildLeng<-lapply(rootChildDesc,phylobase::edgeLength,x=phylo4Obj)
      
      #Problem here!!! if there is single sample in a cluster, then could be a tip with length not equal to zero. Need to ignore these....how? If take the min, then a single zero length in outbranch will result in both having zeros...
      #maybe should change function so have to provide a single name of a sample that is in outbranch so as to identify it that way. 
      #for now, lets hope that never happens! i.e. that BOTH a single sample in a cluster and that outbranch has a zero length
      rootChildNum<-sapply(rootChildLeng,min) #minimum length 
      
      #indicator of which child node is the 
      whPos<-sapply(rootChildNum,function(x){isTRUE(all.equal(x,0))}) #just incase not *exactly* 0
      if(sum(whPos)!=1){
        #if both sides have a zero, then use max instead. 
        rootChildNum<-sapply(rootChildLeng,max) #maximum length 
        whPos<-sapply(rootChildNum,function(x){isTRUE(all.equal(x,0))}) 
      }
      if(sum(whPos)!=1) stop("Internal coding error in finding which is the outbranch in the dendro_samples slot. Please report to git repository!")
      outbranchNode<-rootChild[!whPos]
      
      if(outbranchNode %in% trueInternal){
        outbranchIsInternal<-TRUE
        outbranchNodeDesc<-phylobase::descendants(phylo4Obj,node=outbranchNode,type="ALL") #includes itself
        trueInternal<-trueInternal[!trueInternal%in%outbranchNodeDesc]
        outbranchNodeDesc<-outbranchNodeDesc[outbranchNodeDesc %in% phylobase::getNode(phylo4Obj,type="internal")]
      }
      else outbranchIsInternal<-FALSE
      
    }
    
    phylobase::nodeLabels(phylo4Obj)[as.character(trueInternal)] <- paste("Node", seq_along(trueInternal), sep="")
    #add new label for root 
    if(outbranch){
      phylobase::nodeLabels(phylo4Obj)[as.character(rootNode)] <- "Root"
      if(outbranchIsInternal) phylobase::nodeLabels(phylo4Obj)[ as.character(outbranchNodeDesc) ] <- paste("MissingNode", seq_along(outbranchNodeDesc),sep="")
    }
  }
  else phylobase::nodeLabels(phylo4Obj)<-paste("Node",seq_len(phylobase::nNodes(phylo4Obj)),sep="")
  
  return(phylo4Obj)
}