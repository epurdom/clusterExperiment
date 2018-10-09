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
			newPhyloCluster<-.makePhylobaseTree(object@dendro_clusters,isSamples=FALSE)
			#-----
			#Add necessary information: Position, NodeId, ClusterIdDendro, ClusterIdMerge
			#-----
			nTotal<-length(phylobase::getNode(newPhyloCluster,type="all"))
			clIds<-phylobase::tipLabels(newPhyloCluster)
			clIdsDendro<-rep(NA,length= nTotal)
			clIdsDendro[phylobase::getNode(newPhyloCluster,type="tip")]<-paste("ClusterIdDendro",clIds,sep="")
			clIdsMerge<-rep(NA,length= nTotal)
			mergeCorr<-getMergeCorrespond(object,by="original")
			clIdsMerge[phylobase::getNode(newPhyloCluster,type="tip")]<-paste("ClusterIdMerge",mergeCorr[clIds],sep="")
			maxInternal<-max(as.numeric(gsub("Node","",phylobase::nodeLabels(newPhyloCluster))))
			tipIds<-paste("NodeId",seq(from=maxInternal+1,length=nTips(newPhyloCluster)),sep="")
			phylobase::tipLabels(newPhyloCluster)<-tipIds
			phylobase::nodeLabels(newPhyloCluster)<-gsub("Node","NodeId",phylobase::nodeLabels(newPhyloCluster))
			labs<-phylobase::labels(newPhyloCluster,type="all")
			pos<-factor(rep(NA,length= nTotal),levels=.positionLevels)
			pos[phylobase::getNode(newPhyloCluster,type="tip")]<-"cluster hierarchy tip"
			pos[phylobase::getNode(newPhyloCluster,type="internal")]<-"cluster hierarchy node"
			data.cl<-data.frame(NodeId=labs, ClusterIdDendro=clIdsDendro, ClusterIdMerge=clIdsMerge,              Position=pos,stringsAsFactors=FALSE)
			#get rid of node labels of tips
			phylobase::tipLabels(newPhyloCluster)<-NA
			newPhyloCluster<-phylobase::phylo4d(x=newPhyloCluster, all.data = data.cl)
			ch<-.checkDendroClusterFormat(newPhyloCluster,checkLabels=TRUE)
			if(!is.logical(ch)) stop("Error in converting to new format for dendro_clusters slot:",ch)
			#make tip labels match cluster names
		}
		if(class(object@dendro_samples)=="dendrogram"){
			dataNode<-phylobase::tdata(newPhyloCluster,type="internal")
			data.cl<-phylobase::tdata(newPhyloCluster,type="all")
			newPhyloSample<-.makePhylobaseTree(object@dendro_samples,isSamples=TRUE,outbranch=object@dendro_outbranch)
			
			#-------
			#Add necessary information: Position, NodeId, SampleIndex
			#-------
			origLabs<-phylobase::labels(newPhyloSample,type="all")
			tipNodes<-phylobase::getNode(newPhyloSample,type="tip")
			nTotal<-length(origLabs)
			## SampleIndex
			#Note: Previous naming convention of samples was to give names if not null, and if null, give values "1","2", etc. corresponding to index.
			sampIndex<-rep(NA,length=nTotal)
			if(!is.null(colnames(object))){
				mSample<-match(origLabs,colnames(object))
				if(any(is.na(mSample[tipNodes]))) stop("error in converting tip labels in the dendro_samples slot: the tip labels do not all match colnames of object")
				sampIndex<-mSample
			}
			else{
				tipLabs<-as.numeric(phylobase::tipLabels(newPhyloSample))
				if(any(is.na(tipLabs))) stop("error in converting tip labels in the dendro_samples slot: the object has no colnames and tip labels are not numeric index")
				sampIndex[tipNodes]<-tipLabs
			}
			
			## NodeId
			mNode<-match(origLabs,gsub("Id","",dataNode$NodeId))
			newNodeLabs<-rep(NA,length=length(mNode))
			newNodeLabs[!is.na(mNode)]<-as.character(dataNode$NodeId)[mNode[!is.na(mNode)]]
			#Now need to match nodes that correspond to cluster id... they don't have any reasonable naming scheme...
			whNodes<-phylobase::getNode(newPhyloSample,type="internal")
			origNodes<-phylobase::nodeLabels(newPhyloSample)
			clusterNodes<-origNodes[!is.na(origNodes) & is.na(mNode[whNodes]) ]
			clusterNodes<-clusterNodes[-grep("Missing",clusterNodes)]
			clusterNodes<-clusterNodes[-grep("Root",clusterNodes)]
			#find descendants of each node
			clMat<-clusterLegend(object)[[object@dendro_index]]
			desc<-sapply(clusterNodes,function(x){
				dd<-phylobase::descendants(phy=newPhyloSample,node=x,type="tips")
				
				if(any(is.na(sampIndex[dd]))) stop("coding error -- NA in samples")
				clusterVal<-unique(clusterMatrix(object)[sampIndex[dd],object@dendro_index])
				#if(length(clusterVal)!=1) stop("invalid dendro_samples slot -- multiple clusters in descendant of one cluster node")
				return(paste("ClusterId",clusterVal,sep=""))
			})
			
			
			
			## Postion
			# .positionLevels<-c("cluster hierarchy node","cluster hierarchy tip","tip hierarchy","assigned tip","outbranch hierarchy node","unassigned tip","outbranch root")
			#
			pos<-factor(rep(NA,length= nTotal),levels=.positionLevels)
			pos[grep("Missing",origLabs)]<-"outbranch hierarchy node"
			pos[grep("Root",origLabs)]<-"outbranch root"
			pos[!is.na(match(newNodeLabs,data.cl$NodeId))]<-"cluster hierarchy tip"
			pos[!is.na(mNode)]<-"cluster hierarchy node"
			
			
			
		}
	}

	object<-try(do.call("ClusterExperiment",c(list(object=se,clusters=object@clusterMatrix,checkTransformAndAssay=checkTransformAndAssay),attributes(object)[snames])),silent=TRUE)
	if(!inherits(object,"try-error")){
		return(object)
	} 
		else stop("Attempt to convert did not result in valid object. Here is the error from 'ClusterExperiment':\n",object)
}
)


#---
#remove the outbranch from the dendrogram and from cl
#(note this is using phylo4 obj)
#---
.removeOutbranch<-function(phylo4Obj){
  rootNode<-phylobase::rootNode(phylo4Obj)
  rootChild<-phylobase::descendants(phylo4Obj,node=rootNode,type="children")
  tips<-phylobase::getNode(phylo4Obj,type="tip")
  whMissingNode<-grep("MissingNode",names(rootChild))
  if(length(whMissingNode)==0){
    #check not a single -1 sample from root:
    if(any(rootChild %in% tips)){
      #which ever rootChild is in tips must be single missing sample 
      #because can't make dendrogram with only 1 cluster so couldn't run plot or mergeClusters. 
	#Note only true because outbranch=TRUE
      clusterNode<-rootChild[!rootChild %in% tips]
      #stop("Internal coding error: need to fix .plotDendro to deal with when single missing sample")
    }
    else stop("Internal coding error: no outbranch nodes")	
  } 
  else clusterNode<-rootChild[-whMissingNode]
  if(length(clusterNode)!=1) stop("Internal coding error: removing missing node does not leave exactly 1 descendent of root")
  clusterTips<-phylobase::descendants(phylo4Obj,node=clusterNode,type="tip")
  if(length(clusterTips)==0) stop("Internal coding error: no none missing samples in tree")
  namesClusterTips<-names(clusterTips)
  #
  
  if(is.matrix(cl)) cl<-cl[namesClusterTips,] else cl<-cl[namesClusterTips]
  phylo4Obj<-phylobase::subset(phylo4Obj, node.subtree=clusterNode)
  return(phylo4Obj)
}



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