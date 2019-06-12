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
#' @details Cluster and Sample dendrograms of the class \code{dendrogram} will
#'   be updated to the \code{\link[phylobase]{phylo4d}} class now used in
#'   \code{ClusterExperiment} objects; the merge information on these nodes will
#'   be updated to have the correct format (i.e. match to the internal node id
#'   names in the new dendrogram). The previous identification of nodes that was
#'   previously created internally by plotDendrogram and the merging (labels in
#'   the form of 'Node1','Node2'), will be kept as
#'   \code{\link[phylobase]{nodeLabels}} in the new dendrogram class.
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
  "dendro_index")

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
  "dendro_index", "coClustering",
  "clusterLegend", "orderSamples", "merge_index",
"merge_dendrocluster_index",
"merge_method", "merge_demethod", "merge_cutoff",
"merge_nodeProp", "merge_nodeMerge")
	snames<-snames[snames %in% myslots]

	## Fix class of the dendrograms:
	success<-TRUE
	successMerge<-TRUE
	if("dendro_samples" %in% snames){
		if(is(object@dendro_clusters,"dendrogram")){
			newPhyloCluster<-.makePhylobaseTree(object@dendro_clusters,type="dendro",isSamples=FALSE)
			#-----
			#Add necessary information: Position, NodeId, ClusterIdDendro, ClusterIdMerge
			#-----
			nTotal<-length(phylobase::getNode(newPhyloCluster,type="all"))
			clIds<-phylobase::tipLabels(newPhyloCluster)
			clIdsDendro<-rep(NA,length= nTotal)
			clIdsDendro[phylobase::getNode(newPhyloCluster,type="tip")]<-paste("ClusterId",clIds,sep="")
			clIdsMerge<-rep(NA,length= nTotal)
			mergeCorr<-getMergeCorrespond(object,by="original")
			clIdsMerge[phylobase::getNode(newPhyloCluster,type="tip")]<-paste("ClusterId",mergeCorr[clIds],sep="")
			maxInternal<-max(as.numeric(gsub("Node","",phylobase::nodeLabels(newPhyloCluster))))
			tipIds<-paste("NodeId",seq(from=maxInternal+1,length=nTips(newPhyloCluster)),sep="")
			phylobase::tipLabels(newPhyloCluster)<-tipIds
			labs<-phylobase::labels(newPhyloCluster,type="all")
			labs<-gsub("Node","NodeId",labs)
			pos<-factor(rep(NA,length= nTotal),levels=.positionLevels)
			pos[phylobase::getNode(newPhyloCluster,type="tip")]<-"cluster hierarchy tip"
			pos[phylobase::getNode(newPhyloCluster,type="internal")]<-"cluster hierarchy node"
			data.cl<-data.frame(NodeId=labs, ClusterIdDendro=clIdsDendro, ClusterIdMerge=clIdsMerge,              Position=pos,stringsAsFactors=FALSE)
			#get rid of node labels of tips (keep internal node labels)
			phylobase::tipLabels(newPhyloCluster)<-NA
			newPhyloCluster<-phylobase::phylo4d(x=newPhyloCluster, all.data = data.cl)
			ch<-.checkDendroClusterFormat(newPhyloCluster,checkLabels=TRUE)
			if(!is.logical(ch)) success<-success<-paste("Error in converting to new format for dendro_clusters slot:",ch)
			#make tip labels match cluster names
			object@dendro_clusters<-newPhyloCluster
			
			###Need to match up the merge information with NodeIds:
			if(!is.na(object@merge_nodeMerge) && !"NodeId" %in% colnames(object@merge_nodeMerge) && "Node" %in% colnames(object@merge_nodeMerge)){
				nodes<-object@merge_nodeMerge[,"Node"]
				m<-match(nodes,phylobase::nodeLabels(object@dendro_clusters))
				if(any(is.na(m))) successMerge<-"Error in converting merge_nodeMerge information -- doesn't match node labels of dendro_clusters"
				nodeId<-phylobase::tdata(object@dendro_clusters,type="internal")[m,"NodeId"]
				object@merge_nodeMerge[,"Node"]<-nodeId
				colnames(object@merge_nodeMerge)[grep("Node",object@merge_nodeMerge)]<-"NodeId"
			}
			if(!is.na(object@merge_nodeProp)&& !"NodeId" %in% colnames(object@merge_nodeProp) && "Node" %in% colnames(object@merge_nodeProp)){
				nodes<-object@merge_nodeProp[,"Node"]
				m<-match(nodes,phylobase::nodeLabels(object@dendro_clusters))
				if(any(is.na(m))) successMerge<-"Error in converting merge_nodeProp information -- doesn't match node labels of dendro_clusters"
				nodeId<-phylobase::tdata(object@dendro_clusters,type="internal")[m,"NodeId"]
				object@merge_nodeProp[,"Node"]<-nodeId
				colnames(object@merge_nodeProp)[grep("Node",object@merge_nodeProp)]<-"NodeId"
				
			}
		}
		if(is(object@dendro_samples,"dendrogram")){
			#requires old slot name @dendro_outbranch too.
			dataNode<-phylobase::tdata(newPhyloCluster,type="internal")
			data.cl<-phylobase::tdata(newPhyloCluster,type="all")
			newPhyloSample<-.makePhylobaseTree(object@dendro_samples,type="dendro",isSamples=TRUE,outbranch=object@dendro_outbranch)
			
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
				if(any(is.na(mSample[tipNodes]))) success<-"error in converting tip labels in the dendro_samples slot: the tip labels do not all match colnames of object"
				sampIndex<-mSample
			}
			else{
				tipLabs<-as.numeric(phylobase::tipLabels(newPhyloSample))
				if(any(is.na(tipLabs))) success<-"error in converting tip labels in the dendro_samples slot: the object has no colnames and tip labels are not numeric index"
				sampIndex[tipNodes]<-tipLabs
			}
			
			## NodeId
			mNode<-match(origLabs,gsub("Id","",dataNode$NodeId))
			newNodeLabs<-rep(NA,length=length(mNode))
			newNodeLabs[!is.na(mNode)]<-as.character(dataNode$NodeId)[mNode[!is.na(mNode)]]
			#Now need to match nodes that correspond to cluster id... they don't have any reasonable naming scheme...
			clusterNodes<-phylobase::getNode(newPhyloSample,type="all")
			clusterNodes<-clusterNodes[which(!is.na(newNodeLabs))]
			getClusterDescendants<-function(node){
				dd<-phylobase::descendants(phy=newPhyloSample,node=node,type="tip")
				if(any(is.na(sampIndex[dd]))) success<<-"Not all descendants have sample ids; unable to identify clusters to cluster ids in sample dendrogram"
				else{
					cl<-clusterMatrix(object)[sampIndex[dd],object@dendro_index]
					tab<-table(clusterMatrix(object)[,object@dendro_index])
					tabNode<-table(cl)
					tab<-tab[names(tabNode)]
					clusterVal<-as.numeric(names(tabNode))
					allCluster<-tabNode==tab	
					return(list(clusters=clusterVal,completeCluster=allCluster))
				}
			}
			#of cluster hierarchy nodes, figure out which of them have child that == cluster
			#returns a character matrix where column "Node" is numberic node index and "Cluster" is numeric
			clusterLink<-do.call("rbind",lapply(clusterNodes,function(x){
				cc<-phylobase::descendants(phy=newPhyloSample,node=x,type="children")
				if(!all(cc %in% clusterNodes)){
					cc<-cc[!cc%in%clusterNodes]
					out<-sapply(cc,function(node){
						ret<-NA
						clVal<-getClusterDescendants(node)
						if(length(clVal$clusters)==1){
							if(clVal$completeCluster) ret<-clVal$clusters
						}
						return(ret)
					})
					names(out)<-cc
					out<-out[!is.na(out)]
					return(cbind("Node"=as.numeric(names(out)),"Cluster"=as.numeric( out)))
				}
				else return(NULL)
				})
			)
			#check got all of them
			if(any(is.na(match(na.omit(data.cl$ClusterIdDendro),paste("ClusterId",clusterLink[,"Cluster"],sep=""))))) success<-"error in updating sample dendrogram -- not all clusters found in cluster dendro find match to node in samples dendro"
			nodeIdMatch<-match(paste("ClusterId",clusterLink[,"Cluster"],sep=""),data.cl$ClusterIdDendro)
			newNodeLabs[clusterLink[,"Node"]]<-data.cl$NodeId[nodeIdMatch]
			## Postion
			# .positionLevels<-c("cluster hierarchy node","cluster hierarchy tip","tip hierarchy","assigned tip","outbranch hierarchy node","unassigned tip","outbranch root")
			#
			pos<-factor(rep(NA,length= nTotal),levels=.positionLevels)
			pos[grep("Missing",origLabs)]<-"outbranch hierarchy node"
			pos[grep("Root",origLabs)]<-"outbranch root"
			tipUnassigned<-unlist(phylobase::descendants(phy=newPhyloSample,node=grep("Missing",origLabs),type="tip"))
			pos[tipUnassigned]<-"unassigned tip"
			pos[!is.na(match(newNodeLabs,data.cl$NodeId))]<-"cluster hierarchy tip"
			pos[!is.na(mNode)]<-"cluster hierarchy node"
			tipAssigned<-unlist(phylobase::descendants(phy=newPhyloSample,node=which(pos=="cluster hierarchy tip"),type="tip"))
			pos[tipAssigned]<-"assigned tip"
			pos[is.na(pos)]<-"tip hierarchy"
			sample.data<-data.frame(SampleIndex=sampIndex,NodeId=newNodeLabs,Position=pos,stringsAsFactors=FALSE)
			#get rid of labels
			phylobase::labels(newPhyloSample)<-NA
			newPhyloSample<-phylobase::phylo4d(x=newPhyloSample, all.data = sample.data)
			ch<-.checkDendroSamplesFormat(newPhyloSample,checkLabels=TRUE)
			if(!is.logical(ch)) success<-paste("Error in converting to new format for dendro_samples slot:",ch)
			object@dendro_samples<-newPhyloSample
			
		}
		if(!is.logical(success)){
			warning("Could not successfully convert dendrogram to class 'phylo4d'. Updated object will remove all dendro AND merge related slots. Error:",success)
			snames<-snames[-which(snames %in% dendroSlots)]
			snames<-snames[-which(snames %in% mergeSlots)]
		}
		else if(!is.logical(successMerge)){
			warning("Could not successfully convert merge information to match the nodes in updated dendrogram object. Updated object will remove ALL merge related slots. Error:",successMerge)
			snames<-snames[-which(snames %in% mergeSlots)]
		}
	}

	##Fix class of coClustering -- will save it as sparse
	if(!is.null(object@coClustering) &&
        inherits(object@coClustering,"matrix")){
		coClustering(object)<-object@coClustering
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
#' @importFrom ape as.phylo
#' @importFrom stats as.hclust
#' @importClassesFrom phylobase phylo4 
####
#Convert to object used by phylobase so can navigate easily 
.makePhylobaseTree<-function(x,type,isSamples=FALSE,outbranch=FALSE){
  type<-match.arg(type,c("hclust","dendro"))
  if(type=="dendro"){
    tempPhylo<-try(stats::as.hclust(x),FALSE)
    if(inherits(tempPhylo, "try-error")) stop("the dendrogram object cannot be converted to a hclust class with 'as.hclust.dendrogram'. Check that you gave simple hierarchy of clusters, and not one with fake data per sample")
  }
  if(type=="hclust"){
    #first into phylo from ape package
    tempPhylo<-try(ape::as.phylo(x),FALSE)
    if(inherits(tempPhylo, "try-error")) stop("the hclust object cannot be converted to a phylo class with the methods of the 'ape' package.")
  }

  phylo4Obj<-try(as(tempPhylo,"phylo4"),FALSE) 
  if(inherits(phylo4Obj, "try-error")) stop("the internally created phylo object cannot be converted to a phylo4 class. Check that you gave simple hierarchy of clusters, and not one with fake data per sample")
  
  if(isSamples){
    #NOTE: clusterNodes are found by those with non-zero edge-length between them and their decendents
    nonZeroEdges<-phylobase::edgeLength(phylo4Obj)[which(phylobase::edgeLength(phylo4Obj)>0)] #doesn't include root
    trueInternal<-sort(unique(as.numeric(sapply(strsplit(names(nonZeroEdges),"-"),.subset2,1)))) #this also picks up the outbranch between -1,-2
    #old way of doing it:
    #clusterNodes<-sort(unique(unlist(phylobase::ancestors(phylo4Obj,node=phylobase::getNode(phylo4Obj,type="tip"),type="parent"),recursive=FALSE,use.names=FALSE)))
    if(outbranch){#remove root from labeling if -1 outbranch
      #######
      #remove root
      #######
      rootNode<-phylobase::rootNode(phylo4Obj)
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
      whKeep<-sapply(rootChildNum,function(x){isTRUE(all.equal(x,0))}) #just incase not *exactly* 0
      if(sum(whKeep)!=1){
        #if both sides have a zero, then use max instead. 
        rootChildNum<-sapply(rootChildLeng,max) #maximum length 
        whKeep<-sapply(rootChildNum,function(x){isTRUE(all.equal(x,0))}) 
      }
      if(sum(whKeep)!=1) stop("Internal coding error in finding which is the outbranch in the dendro_samples slot. Please report to git repository!")
      outbranchNode<-rootChild[!whKeep]
      
      if(outbranchNode %in% trueInternal){
        outbranchIsInternal<-TRUE
        outbranchNodeDesc<-phylobase::descendants(phylo4Obj,node=outbranchNode,type="ALL") #includes itself
        trueInternal<-trueInternal[!trueInternal%in%outbranchNodeDesc]
        outbranchNodeDesc<-outbranchNodeDesc[outbranchNodeDesc %in% phylobase::getNode(phylo4Obj,type="internal")]
      }
      else outbranchIsInternal<-FALSE
      
    }
    #trueInternal<-allInternal[!allInternal%in%clusterNodes]
    
    phylobase::nodeLabels(phylo4Obj)[as.character(trueInternal)]<-paste("Node",seq_along(trueInternal),sep="")
    #add new label for root 
    if(outbranch){
      phylobase::nodeLabels(phylo4Obj)[as.character(rootNode)]<-"Root"
      if(outbranchIsInternal) phylobase::nodeLabels(phylo4Obj)[as.character(outbranchNodeDesc)]<-paste("MissingNode",seq_along(outbranchNodeDesc),sep="")
    }
  }
  else phylobase::nodeLabels(phylo4Obj)<-paste("Node",seq_len(phylobase::nNodes(phylo4Obj)),sep="")
  
  return(phylo4Obj)
}
