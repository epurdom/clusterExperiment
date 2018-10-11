.makeSampleNames<-function(x){paste("Sample",as.character(x))}

.clusterDendroColumns<-c("Position","ClusterIdDendro","ClusterIdMerge","NodeId")
.clusterSampleColumns<-c("Position","NodeId","SampleIndex")

.positionLevels<-c("cluster hierarchy node","cluster hierarchy tip","tip hierarchy","assigned tip","outbranch hierarchy node","unassigned tip","outbranch root")

.hasOutBranch<-function(object){
	"outbranch root" %in% phylobase::tdata(object@dendro_samples)$Position
}

#' @param inputValue a vector of values to match to dendro data
#' @param dendro a phylo4d object
#' @param matchColumn a character value giving the name of a column in the \code{tdata} of the dendro to match to. If "NodeIndex", means that will match to the row number (i.e. inputValue is a integer giving the row index); in this case it is the same as doing tdata(dendro)[inputValue,returnColumn]
#' @param returnColumn a character value giving the name of a column to return -- i.e. after matching to the matchColumn, return the corresponding values in the returnColumn. If returnColumn="NodeIndex", then just the index of the match is returned.
#' @param returnColumn
#' @noRd
.matchToDendroData<-function(inputValue,dendro,matchColumn="NodeId",returnColumn){
	#Note that matchColumn="NodeIndex" means inputs give rows of the df. returnColumn="NodeIndex" means to return the index that matches (i.e. nodeIndex too)
	df<-phylobase::tdata(dendro,type="all")
	if(!is.character(returnColumn) || !is.character(matchColumn)) stop("coding error -- returnColumn and matchColumn must be character valued")
	if(!returnColumn=="NodeIndex" && !returnColumn %in% names(df)) stop("coding error -- invalid value of returnColumn ")
	if(!matchColumn=="NodeIndex" && !matchColumn %in% names(df)) stop("coding error -- invalid value of matchColumn")
	if(!matchColumn=="NodeIndex") m<-match(inputValue, df[,matchColumn]) else m<-inputValue
	if(any(is.na(m))) stop("coding error -- invalid value in inputValue (not in dendro)")
	if(returnColumn=="NodeIndex") return(m)
	else return(df[m,returnColumn])
}

.getClusterIds<-function(tipIndex,clusterDendro,returnValue=c("ClusterIdDendro","ClusterIdMerge")){
	returnValue<-match.arg(returnValue)
  tipNames <- .matchToDendroData(tipIndex,clusterDendro,matchColumn="NodeIndex",returnColumn=returnColumn)
  tipNames<-gsub("ClusterId","",as.character(tipNames))
  return(tipNames)
}


#' @description phylo4d objects for cluster and sample dendrograms with pretty (and matched!) node and tip labels.
#' @param object a ClusterExperiment object
#' @param labelType. If 'id' returns the NodeId value for all the internal nodes that match the cluster dendrogram (cluster and sample dendrograms) and the cluster id for the tips (cluster dendrogram) or sample index for the tips (sample dendrogram). If "name", then the cluster ids and sample indices are converted into the cluster names (in the clusterMat) and the sample names (in colnames)
#' @param useMergeClusters if TRUE and there is an active merge, will remove the dendrogram cluster ids, and instead use merge cluster ids (which means will nave no label for dendrogram cluster merged)
#' @return Returns list of the two dendrograms with nodes that have been updated. Note they do not match requirement of the clusterExperiment object because have labels for nodes they "shouldn't"
#' @details  Different from convertToPhyClasses, which is trying to get the needed info into the phylo or phylo4 class that doesn't have tdata. Calls that function internally
#' @noRd
.setNodeLabels<-function(object,labelType=c("name","id"),useMergeClusters=FALSE,overrideExistingNode=FALSE,singletonCluster=c("sample","cluster"),...){
	labelType<-match.arg(labelType)
	singletonCluster<-match.arg(singletonCluster)
	 if(!inherits(object,"ClusterExperiment")) stop("coding error -- function should be used for ClusterExperiment objects")
	 if(is.na(object@merge_index)|| all(is.na(phylobase::tdata(object@dendro_clusters)$ClusterIdMerge))) useMergeClusters<-FALSE
	 if(is.null(object@dendro_clusters)) return(NULL)
	 #internal option returns same as the .convertPhyClasses with both =TRUE
	 internalCluster<-.convertToPhyClasses(object@dendro_clusters, returnClass=c("phylo4d"), convertNodes=TRUE, convertTips=TRUE)
	 internalSamples<-.convertToPhyClasses(object@dendro_samples, returnClass=c("phylo4d"), convertNodes=TRUE, convertTips=TRUE)
	 
	 
	 if(useMergeClusters){
		 #erase all the dendro cluster Ids, which are only on the tips
		 whConvert<-which(!is.na(phylobase::tdata(internalCluster,type="tip")$ClusterIdDendro))
		 phylobase::tipLabels(internalCluster)[whConvert]<-NA 
		 #now add merge cluster ids.
		 whConvert<-which(!is.na(phylobase::tdata(internalCluster,type="all")$ClusterIdMerge))
		 phylobase::labels(internalCluster)[whConvert] <- phylobase::tdata(internalCluster,type="all")$ClusterIdMerge[whConvert]
		 
	 }
	 if(!overrideExistingNode){
		 #located here so merge labels won't be displayed if user changed that node name manually.
		 existingNodeLabels<-phylobase::nodeLabels(object@dendro_clusters)
		 whKeep<-which(!is.na(existingNodeLabels) & existingNodeLabels!=phylobase::tdata(internalCluster,type="internal")$NodeId)
		 phylobase::nodeLabels(internalCluster)[whKeep]<-existingNodeLabels[whKeep]
	 }
	 	 
	 
	 if(labelType=="name"){
		 ###Convert to the name versions
		 if(useMergeClusters){
			 legMat<-clusterLegend(object)[[object@merge_index]]
			 labs<-phylobase::labels(internalCluster)[whConvert]
			 labs<-as.character(as.numeric(gsub("ClusterId","",labs)))
			 m<-match(labs,legMat[,"clusterIds"])	
			 phylobase::labels(internalCluster)[whConvert]<-legMat[m,"name"]
		 }
		 else{
			 legMat<-clusterLegend(object)[[object@dendro_index]]
		 	 labs<-phylobase::tipLabels(internalCluster)
			 labs<-as.character(as.numeric(gsub("ClusterId","",labs)))
			 m<-match(labs,legMat[,"clusterIds"])	
			 phylobase::tipLabels(internalCluster)<-legMat[m,"name"]
		 }
		 
		 nms<-if(!is.null(colnames(object))) colnames(object) else .makeSampleNames(seq_len(ncol(object)))
		 phylobase::tipLabels(internalSamples) <- nms[as.numeric(phylobase::tipLabels(internalSamples))]

			 
	 }
	 #Make sample node names match
	 type<-switch(singletonCluster,"sample"="internal","cluster"="all")
	 vals<-phylobase::tdata(internalSamples,type=type)$NodeId
	 wh<-which(!is.na(vals))
	 mToSample <- .matchToDendroData(inputValue=vals[wh], dendro=internalCluster, matchColumn="NodeId", returnColumn="NodeIndex")
	 phylobase::labels(internalSamples,type=type)[wh]<- phylobase::labels(internalCluster)[mToSample]
	 
	 
	 return(list(dendro_clusters=internalCluster,dendro_samples=internalSamples))
}


#' @importFrom ape as.hclust.phylo
#' @importFrom stats as.dendrogram
.convertToDendrogram<-function(x){
	if(inherits(x,"dendrogram") )return(x)
	if(inherits(x,"phylo4")){
		x<-.convertToPhyClasses(x,"phylo")
	} 
	if(inherits(x,"phylo")){
		x<-try(ape::as.hclust.phylo(x),FALSE)
		if(inherits(x, "try-error")) stop("coding error -- could not convert to hclust object. Reported error:",x)
	}
	if(inherits(x,"hclust")){
		x<-try(stats::as.dendrogram(x),FALSE)
		if(inherits(x, "try-error")) stop("coding error -- could not convert from hclust to dendrogram object. Reported error:",x)
		return(x)
	}
	else{ stop("input x is not of hclust, phylo4 or phylo class")}
}



#' @param convertNodes logical. If true, the returned dendrogram will have the node labels change to be NodeId
#' @param convertTips logical. If true, the returned dendrogram will have the tip labels changed. If 'ClusterIdDendro' is in tdata (i.e. its a cluster dendrogram) then they are converted to Dendro cluster id (if no NAs) or if there are NAs, then merge cluster id (if dendros are NA in tips). If 'SampleIndex' is column name in tdata (i.e. is a samples dendrogram), then returns tip labels that are the SampleIndex value
#' @importFrom stats as.hclust
#' @importFrom ape as.phylo.hclust
#' @importClassesFrom phylobase phylo4 
.convertToPhyClasses<-function(x,returnClass=c("phylo4","phylo","phylo4d"),convertNodes=FALSE,convertTips=FALSE){
	returnClass<-match.arg(returnClass)
	if(inherits(x,"phylo4d") ){
		if(returnClass %in% c("phylo","phylo4","phylo4d")){
			#------
			#Before convert,
			#make internal node and cluster ids 
			#the node and tip labels (i.e. erase existing)
			#------
			if(convertNodes) phylobase::nodeLabels(x)<-as.character(phylobase::tdata(x,type="internal")$NodeId)
	  
			if(convertTips){
				if("ClusterIdDendro" %in% names(phylobase::tdata(x,type="all"))){
					phylobase::tipLabels(x)<-as.character(phylobase::tdata(x,type="tip")$ClusterIdDendro)
					if(any(is.na(phylobase::tipLabels(x)))){
						#this should mean its the cluster dendrogram limited to the merge ids.
						phylobase::tipLabels(x)<-as.character(phylobase::tdata(x,type="tip")$ClusterIdMerge)
					}	
					if(any(is.na(phylobase::tipLabels(x)))) stop("coding error -- even after accounting for merge id, still have NA values")
				}
				else if("SampleIndex" %in% names(phylobase::tdata(x,type="all")))
					phylobase::tipLabels(x)<-as.character(phylobase::tdata(x,type="tip")$SampleIndex)
				else stop("coding error -- tree should have either 'ClusterIdDendro' or 'SampleIndex' as names in tdata(x)")
			} 			
		}
		if(returnClass %in% c("phylo4","phylo4d")) return(x)
	}
	if(!inherits(x,"phylo4d") & returnClass=="phylo4d") stop("coding error -- can't convert other classes to phylo4d at this time and still retain all correct information")
	if(inherits(x,"dendrogram")){
		x<-try(stats::as.hclust(x),FALSE)
		if(inherits(x, "try-error")) stop("coding error -- could not convert from dendrogram to hclust object. Reported error:",x)
	}
	if(inherits(x,"hclust")){
		x<-try(ape::as.phylo.hclust(x),FALSE)
		if(inherits(x, "try-error")) stop("coding error -- could not convert from hclust to phylo object. Reported error:",x)
	}

	if(inherits(x,"phylo")){
		if(returnClass=="phylo") return(x)
		else{
			x<-try(as(x,"phylo4"),FALSE)
			if(inherits(x, "try-error")) stop("coding error -- could not convert from phylo to phylo4 object. Reported error:",x)
			return(x)
		}
	}
	if(inherits(x,"phylo4")){
		if(returnClass=="phylo4") return(x)
		else{
			#phylobase warnings that trees with unknown edge order may be unsafe. This is probably because of this problem they ran into: http://lists.r-forge.r-project.org/pipermail/phylobase-devl/2009-January/000353.html
			#but the problem has probably been fixed by now in ape!
			x<-try(suppressWarnings(as(x,"phylo")),FALSE)
			if(inherits(x, "try-error")) stop("coding error -- could not convert from phylo4 to phylo object. Reported error:",x)
			return(x)
		}
	}
	else{ stop("input x is not of hclust, dendrogram, phylo or phylo4 class")}
	
}

.pruneToNodes<-function(phylo4,nodesPruned){
	##nodesPruned should be *the node numerical index* of those nodes that need "pruning" -- i.e. drop their children; other nodes kept the same.
	tipsPruned<-unlist(phylobase::descendants(phylo4,node=nodesPruned,type="tip"),use.names=FALSE)
	tipsKept<-phylobase::getNode(phylo4,type="tip")
	tipsKept<-tipsKept[!tipsKept %in% tipsPruned]
	nodesKeep<-c(nodesPruned,tipsKept)
	
	allNodes<-phylobase::getNode(phylo4,type="all")
	if(!all(nodesKeep %in% allNodes)) 
		stop("nodes specified are not valid")
	nodesKeep<-unique(unlist(phylobase::ancestors(phylo4,nodesKeep,"ALL"),use.names=FALSE)) #ALL means include self
	names(nodesKeep)<-names(allNodes)[nodesKeep]
	whEdgesKeep<-which(phylobase::edges(phylo4)[,2] %in% nodesKeep) 
	
	newEdge<-phylobase::edges(phylo4)[whEdgesKeep,]
	
	#need to get rid of singleton nodes -- i.e. nodes with only 1 descendant
	#this won't happen I think in my merge tree example, so make it error...
	singletonNode<-if(length(which(table(newEdge[-which(newEdge[,1]==0),1])==1))>0) stop("coding error -- removing these nodes creates internal nodes that are singletons. Should be only choice of nodes that 'prunes' back the tree, not actually subsetting.")
	
	#----------
	#have to renumber everything so consecutive... technically documentation doesn't say that, but actually do (error otherwise)
	#----------
	#get the current numbers of tips, internal nodes, and root
	#this is so know how to match the names, etc. back to them
	newTips<-newEdge[which(!newEdge[,2] %in% newEdge[,1]),2]
	newNodes<-nodesKeep
	root<-newEdge[newEdge[,1] ==0,2]
	whRoot<-which(newEdge[,1] ==0)

	#tips: make them 1:(#tips)
	#tips should be in edge matrix once:
	whEdgeTips<-match(newTips,newEdge) #find edges that contain tips
	consecutiveTips<-match(newTips,sort(newTips)) #give them consecutive numbers in same order as previous
	newEdge[whEdgeTips]<-consecutiveTips
	whNodeTips<-match(newTips,newNodes) #this is for reordering
	newNodes[whNodeTips]<-consecutiveTips
	
	#internal nodes: make them consecutive but in same order as previously, then add root number + nTips
	consecutiveNodes<-match(newNodes[-whNodeTips],sort(newNodes[-whNodeTips])) #doesn't include the 0 edge to root, so root is #1
	mNodesToEdge<-match(newEdge[-whEdgeTips],newNodes[-whNodeTips]) #get one NA from 0
	if(sum(is.na(mNodesToEdge))!=1) stop("coding error -- should have one na because of root")
	newEdge[-whEdgeTips]<-consecutiveNodes[mNodesToEdge]+length(newTips)
	newNodes[-whNodeTips]<-consecutiveNodes+length(newTips)	
	newEdge[whRoot,1]<-0
	if(any(sort(newNodes)!=1:length(newNodes))) stop("coding error -- did not result in consecutive node numbers")
		
	#now need to make sure everything in right order with new edges:
	newdata<-phylo4@data[nodesKeep,]
	row.names(newdata)<-as.character(newNodes)
	newdata<-newdata[order(newNodes),]
	
	tipLabels<-names(newNodes)[whNodeTips]
	tipLabels<-tipLabels[order(newNodes[whNodeTips])]
	nodeLabels<-names(newNodes)[-whNodeTips]
	nodeLabels<-nodeLabels[order(newNodes[-whNodeTips])]

	##edge length slots: edge, edge.length, edge.label
	##node length slots: data (rows), label
	##Don't change: order, metadata, annote
	return(phylobase::phylo4d(x=newEdge, 
		edge.length = phylo4@edge.length[whEdgesKeep], 
		tip.label = tipLabels,
		node.label = nodeLabels, 
		edge.label = phylo4@edge.label[whEdgesKeep], 
		order = phylo4@order,
		annote = phylo4@annote,
		all.data=newdata,
		metadata=phylo4@metadata)
	)	
	
}

.safePhyloSubset<-function(phylo4,tipsRemove,nodeName){
  if(length(phylobase::tipLabels(phylo4))-length(tipsRemove)<2){
    ###Check that would have >1 tips left after remove (otherwise gives an error, not sure why with trim.internal=FALSE; should report it)
    ###Remove all but 1 tip seems to work -- collapse down desptie trim.internal=FALSE. Very weird.
    keptTip<-TRUE
    tipKeep<-names(tipsRemove)[1] #label of the tip removed (tipsRemove has internal names as value)
    tipsRemove<-tipsRemove[-1] 
  }
  else keptTip<-FALSE
  phylo4<-phylobase::subset(phylo4,tips.exclude=tipsRemove,trim.internal =FALSE)
  #have to give that 
  if(keptTip){
    labs<-phylobase::tipLabels(phylo4)
    wh<-which(labs==tipKeep)
    labs[wh]<-nodeName
    phylobase::tipLabels(phylo4)<-labs
  }
  return(phylo4)
}


###Note, cluster nodes DEFINED as those whose descendants are of length>0. So edge from root of these fake binary trees needs to be > 0, and rest =0 
### But also use this to make fake binary when n<5 samples and need it to stay ultrametric. Here, the edgeLength>0 will make it not ultrametric! Sigh. I fix that issue after returned (i.e. not here)
.makeFakeBinary<-function(tipNames,rootEdgeLength=0,edgeLength=0){
	newPhylo<-list()
	n<-length(tipNames)
	if(n>1){
		if(n>2){
			newPhylo$edge<-cbind(seq(from=n+1,to=n+(n-1)-1,by=1),seq(from=n+2,to=n+(n-1),by=1))
			newPhylo$edge<-rbind(newPhylo$edge,cbind(n+1:(n-1),1:(n-1)))
			newPhylo$edge<-rbind(newPhylo$edge,c(n+n-1,n))
		}
		else if(n==2){
			newPhylo$edge<-rbind(c(3,1),c(3,2))		
		}
		newPhylo$tip.label<-tipNames
		newPhylo$Nnode<-(n-1)
		newPhylo$edge.length<-rep(edgeLength,length=nrow(newPhylo$edge))
		whRoot<-which(newPhylo$edge==n+1,arr.ind = TRUE)
		if(nrow(whRoot)>2) stop("coding error -- found more than two descendants of root")
		if(any(whRoot[,"col"]==2)) stop("coding error -- found root as a descendant")
		newPhylo$edge.length[whRoot[,"row"]]<-rootEdgeLength
	}
	else{
		if(n==1){
			newPhylo<-list()
			newPhylo$edge<-matrix(NA,nrow=0,ncol=2)
			newPhylo$tip.label<-tipNames
			newPhylo$edge.length<-rootEdgeLength
			newPhylo$Nnode<-0
		}
		
		else stop("coding error -- zero or less length tipNames")
	}
	newPhylo$edge<-.makeIntegerMatrix(newPhylo$edge)
	
	class(newPhylo)<-"phylo"
	return(newPhylo)
}


.makeIntegerMatrix<-function(mat){
	#make the entries integer valued.
	matrix(as.integer(mat),ncol=ncol(mat),byrow=FALSE)
}


#' @importFrom ape node.depth.edgelength
.mergePhylo<-function(tree1,tree2,mergeEdgeLength,balanceLength=TRUE){
	#if balanceLength==TRUE, want to rebalance so that the two trees have same node height -- without adding to 0-length edges!
	n1<-length(tree1$tip.label)
	m1<-tree1$Nnode
	isTree<-is.list(tree2) & class(tree2)=="phylo" #whether tree2 merging is singleton or actual tree (different protocal from .addTipsToTrees...)
	mergeEdgeLength1<-mergeEdgeLength2<-mergeEdgeLength
	depth1<-max(ape::node.depth.edgelength(tree1))	
	if(isTree){
		if(is.null(tree1$edge.length) || is.null(tree1$edge.length) ) stop("must have edge length on both trees")
		n2<-length(tree2$tip.label)
		m2<-tree2$Nnode
		if(balanceLength==TRUE){
			depth2<-max(ape::node.depth.edgelength(tree2))	
			finalDepth<-max(depth1,depth2)
			
			if(depth1==0){
				depth1<-1
				mergeEdgeLength1<-finalDepth+mergeEdgeLength
			}
			if(depth2==0){
				depth2<-1
				mergeEdgeLength2<-finalDepth+mergeEdgeLength	
			}
			tree1$edge.length<-tree1$edge.length/depth1*finalDepth
			tree2$edge.length<-tree2$edge.length/depth2*finalDepth
		}
	}else{
		if(is.null(tree1$edge.length)) stop("must have edge length to merge to tree")
		if(!is.vector(tree2) || length(tree2)!=1) stop("coding error -- trying to merge add more than 1 tip to tree")
		n2<-1
		m2<-0
		mergeEdgeLength2<-depth1+mergeEdgeLength
	}
	newPhylo<-list()
	
	#1) Tree 1 inner nodes: add n2+1 (because now will have new root)
	whinner1<-which(tree1$edge[,2] > n1) #including root of tree 1
	if(length(whinner1)>0) tree1$edge[whinner1,2]<-tree1$edge[whinner1,2]+n2+1
	tree1$edge[,1]<-tree1$edge[,1]+n2+1
	#2) Number of internal nodes
	newPhylo$Nnode<-m1+m2+1
	#3) Rescue node labels in tree1:
	#were 1:tree1$Nnode
	#in edge were +n1
	#moved to +n2+1 in new edge matrix
	#in node index, subtract n1+n2 -> +1
	newPhylo$node.label<-rep(NA,times=newPhylo$Nnode)
	if(!is.null(tree1$node.label)) newPhylo$node.label[1:m1+1]<-tree1$node.label

	if(isTree){
		#4) Tree 2 tips: add n1
		whtip2<-which(tree2$edge[,2]<=n2)
		tree2$edge[whtip2,2]<-tree2$edge[whtip2,2]+n1
	
		#5) Tree 2 inner nodes: -n2 to get to 1:... then add n2+n1+m1+1 -> +n1+m1+1
		if(nrow(tree2$edge)!=length(whtip2)) tree2$edge[-whtip2,2]<-tree2$edge[-whtip2,2]+n1+m1+1
		tree2$edge[,1]<-tree2$edge[,1]+n1+m1+1
	
		#6) Add root edge (n2+n1+1) to (n1+n2+2) and (n2+m1+n2+2) to bottom of edge matrix:
		newPhylo$edge<-rbind(tree1$edge,tree2$edge)
		newPhylo$edge<-rbind(newPhylo$edge, rbind(c(n2+n1+1,n1+n2+2),c(n2+n1+1,n1+m1+n2+2)))
		#7) tip.label
		newPhylo$tip.label<-c(tree1$tip.label,tree2$tip.label)
		#8) 		#should there be check that there is edge length?
		newPhylo$edge.length<-c(tree1$edge.length,tree2$edge.length,c(mergeEdgeLength1,mergeEdgeLength2))
		#9) Rescue node labels in tree 2:
		#were 1:tree2$Nnode
		#in edge were +n2
		#moved to +n1+m1+1 in new edge matrix
		#in node index, subtract n1+n2 -> +m1+1
			if(!is.null(tree2$node.label)) newPhylo$node.label[1:m2+m1+1]<-tree2$node.label

	}
	else{
		#6) Add root edge (n2+n1+1) to tree 1 (n1+n2+2) and root to tip to bottom of edge matrix:
		newPhylo$edge<-tree1$edge
		newPhylo$edge<-rbind(newPhylo$edge, rbind(c(n2+n1+1,n1+n2+2),c(n2+n1+1,n1+n2)))
		#7) tip.label: c(tree1$tip.label,tree2$tip.label)
		newPhylo$tip.label<-c(tree1$tip.label,tree2)		
		#8) 
		newPhylo$edge.length<-c(tree1$edge.length,c(mergeEdgeLength1,mergeEdgeLength2))
	}
	newPhylo$edge<-.makeIntegerMatrix(newPhylo$edge)
	class(newPhylo)<-"phylo"
	return(newPhylo)
	
}


#Input:
#tipTrees should be list of trees in same order as tips in mainTree
#assumes rooted tree
#single values to be added to tree should have $edge=0-row matrix, $Nnode=0, $tip.label name of cluster, $edge.length=value to be added to length of cluster node
#Output:

.addTreesToTips<-function(mainTree,tipTrees){
	N<-length(mainTree$tip.label)
	M<-mainTree$Nnode
	
	if(length(tipTrees)!=N) stop("coding error -- length of list of trees must match number of tips in main Tree")
	
	#check valid phylo objects
	isPhylo<-sapply(tipTrees,class)=="phylo"
	if(!all(isPhylo)) stop("coding error -- each of the new trees must be of class phylo")
	# check<-sapply(tipTrees,.testPhyloObject,verbose=FALSE)
	# check<-.testPhyloObject(mainTree)
	
	
	nVec<-sapply(tipTrees,function(x){length(x$tip.label)})
	mVec<-sapply(tipTrees,function(x){x$Nnode})
	isSingle<-sapply(tipTrees,function(x){nrow(x$edge)==0})
	cumN<-c(0,cumsum(nVec))
	cumM<-c(0,cumsum(pmax(mVec-1,0))) #cumulative number of NON-ROOT internal
	cumSingle<-c(0,cumsum(as.numeric(isSingle))) #counts how many clusters/tips on main tree were singleton before i
	#don't need the last entry.
	cumN<-head(cumN,-1) 
	cumM<-head(cumM,-1) 
	cumSingle<-head(cumSingle,-1) 
	
	##############
	#change mainTree edge matrix
	##############
	whRoot<-which(mainTree$edge==N+1)
	whInternal<-which(mainTree$edge>N)
	#1) inner add sum(nVec)-N (including root)
	#note root becomes sum(nVec)+1, and last internal edge should be sum(nVec)+M
	mainTree$edge[whInternal]<-mainTree$edge[whInternal]-N+sum(nVec)
	#2) tips add sum(nVec)+M so now internal edges  -- except those that will be continue being tip because singleton added!
	whTips<-which(mainTree$edge<=N & !mainTree$edge %in% which(isSingle))
	# main tree tips that are not single get consecutive numbers
	mainTree$edge[whTips]<-mainTree$edge[whTips]+sum(nVec)+M-cumSingle[!isSingle]
	if(sum(isSingle)>0){
		whSingle<-match(which(isSingle),mainTree$edge[,2])
		mainTree$edge[whSingle,2]<-cumN[isSingle]+1
	
	#Note that singletons, unlike others, short value of rootEdgeLength to their edge length, comparatively, when made with fakeBinaryTrees function. (and rootEdgeLength was added to deal with identification of cluster problem)
	#Need to add to those
		addValues<-sapply(tipTrees[isSingle],.subset2,"edge.length")
		mainTree$edge.length[whSingle]<-mainTree$edge.length[whSingle]+addValues
	}

	newPhylo<-list()
	newPhylo$Nnode<-mainTree$Nnode+sum(sapply(tipTrees,.subset2,"Nnode"))
	
	#------
	#Add back original node labels of the tree:
	#------
	newPhylo$node.label<-rep(NA,length=newPhylo$Nnode)
	#previously internal in main were 1:mainTree$Nnode+N in edge branch; 
	#Moved to -N +sum(nVec)
	#Then to move to nodes index by subtracting off length(newPhylo$tip.label)=sum(nVec)
	#so same index...
	newPhylo$node.label[1:mainTree$Nnode]<-mainTree$node.label

	# main tree tips that are not single become nodes -- same formula as before
	# were c(1:N)[!isSingle] in edge matrix (and tip.label)
	# Then moved by +sum(nVec)+M-cumSingle[!isSingle] in edge matrix
	# Then subtract of sum(nVec) to get to node index
	notSingleTips<-c(1:N)[!isSingle]
	newPhylo$node.label[notSingleTips+M-cumSingle[!isSingle]]<-mainTree$tip.label[notSingleTips]
	##############
	#change treeList edgeMatrix
	##############
	tipTrees<-lapply(1:length(tipTrees),function(ii){
		n<-nVec[[ii]]
		x<-tipTrees[[ii]]
		if(nrow(x$edge)>0){
			whInternal<-which(x$edge>n+1)
			whRoot<-which(x$edge==n+1)
			whTips<-which(x$edge<=n)
			x$edge[whRoot]<-ii+sum(nVec)+M-cumSingle[ii]
			# give new numbers to internal:
			# 1) subtract n[i]-1 to non-root internal nodes go from 1:(m[i]-1)
			# 2) + no. new tips (sum(nVec))
			# 3) + no. internal of main tree (M) 
			# 4) + total no. tips of main that become internal (i.e. not singletons) (sum(as.numeric(!isSingle)))
			# 5) + no. non-root internal of the ii-1 trees already added (cumM[ii])
			x$edge[whInternal]<-x$edge[whInternal]-(n+1)+ sum(nVec) + M  + sum(as.numeric(!isSingle)) + cumM[ii]
			x$edge[whTips]<-x$edge[whTips]+ cumN[ii]
		}
		return(x)
	})
	
	newPhylo$edge<-do.call("rbind",c(lapply(tipTrees[!isSingle],.subset2,"edge"),list(mainTree$edge)))
	newPhylo$tip.label<-do.call("c",lapply(tipTrees,.subset2,"tip.label"))
	newPhylo$edge.length<-do.call("c",c(lapply(tipTrees[!isSingle],.subset2,"edge.length"),list(mainTree$edge.length)))
	newPhylo$edge<-.makeIntegerMatrix(newPhylo$edge)

	#singleton clusters stay tips... make sure they get the right label (namely cluster, not sample -- sort of arbitrary choice)
	# were tip.label/edge values: which(isSingle)
	# move to = cumN[isSingle]+1 (still tip values don't subtract off
	newPhylo$tip.label[cumN[isSingle]+1]<-mainTree$tip.label[which(isSingle)]
	
	class(newPhylo)<-"phylo"
	return(newPhylo)	
}

.testPhyloObject<-function(phyloObj,verbose=FALSE){
	if(verbose) cat("Tests (should all be yes):\n")
	.testPrint("Required elements with right name:",all(c("tip.label","edge","Nnode") %in% names(phyloObj)),verbose=verbose)
	ntip<-length(phyloObj$tip.label)
	ninterior<-phyloObj$Nnode
	allNodes<-unique(as.numeric(phyloObj$edge))

	.testPrint("Edge matrix has two columns",dim(phyloObj$edge)[2]==2,verbose=verbose)

	#no gaps in series
	.testPrint("No gaps", all(sort(allNodes)==1:length(allNodes) ),verbose=verbose)

	#matches expected range
	.testPrint("Range of series matches expected", length(allNodes)== ntip+ninterior ,verbose=verbose)
	
	#All elements, except the root
	#n+1, appear once in the second column.
	rootVal<-ntip+1
	tab<-table(phyloObj$edge[,2])
	.testPrint("All elements except root once in second column",all(tab[names(tab)!= as.character(rootVal)]==1),verbose=verbose)
	
	if("edge.length" %in% names(phyloObj)){
		.testPrint("edge length matches number edges",nrow(phyloObj$edge)==length(phyloObj$edge.length),verbose=verbose)
	}
	if(verbose){
		cat("-----------\n")
		cat("Number of tips:",ntip,"\n")
		cat("Number of interior:",ninterior,"\n")
		cat("Tips labels:",paste(phyloObj$tip.label,collapse=","),"\n")
		cat("Range of edge matrix values:",paste(range(allNodes),collapse=","),"\n")		
	}

}
.testPrint<-function(name,logic,verbose=FALSE){
	if(!is.logical(logic)) stop("logic input is not logical")
	if(!verbose){
		if(!logic)stop(paste("coding error -- the following check for phylo class did not pass",name))
	}
	else{
		cat(paste0(name,":"))
		if(logic) cat("Yes\n") else cat("No\n")		
	}
	
}

.force.ultrametric<-function(tree){
	##From http://blog.phytools.org/2017/03/forceultrametric-method-for-ultrametric.html
	#calculates, per tip, the amount missing and adds it to the tips.
	if(inherits(tree,"phylo")){
		depth<-ape::node.depth.edgelength(tree) #in order of all nodes
		ntips<-length(tree$tip.label)
		maxD<-max(depth[1:ntips])
		addValue<-maxD-depth[1:ntips]
		whTips<-match(1:ntips,tree$edge[,2])
		tree$edge.length[whTips]<-tree$edge.length[whTips]+addValue
		return(tree)
	}
	else if(!inherits(tree,"phylo4")){
		stop("tree must be of class phylo4")
	} else{
		allTips<-phylobase::tipLabels(tree)
		depthToTips<-phylobase::nodeHeight(tree,allTips,from="root")
		maxD<-max(depthToTips)
		addValue<-maxD-depthToTips
		allLen<-phylobase::edgeLength(tree)
	  edgeMat<-phylobase::edges(tree)
  
	  #add 'addValue' to the tips so ultrametric
	  tipIds<-as.numeric(names(allTips))
	  m<-match(tipIds,edgeMat[,2])
	  edgeIds<-paste(edgeMat[m,1],edgeMat[m,2],sep="-")

	  #check didn't do something stupid:
	  checkTipEdges<-phylobase::edgeId(tree,type="tip")
	  if(!identical(sort(unname(checkTipEdges)),sort(unname(edgeIds)))) stop("coding error -- didn't correctly get edge ids for tips")

	  #replace with new edges:
		allLen[edgeIds]<-allLen[edgeIds]+addValue
		phylobase::edgeLength(tree)<-allLen
		return(tree) #returns phylo4 tree
	}

}



# #---------
# 		  #old code that made non-binary tree by hand.
# #---------
# 		  tab<-table(cl[whPos])
# 		  newPhylo<-list()
#
# 		  if(any(tab==1)){
# 	cluster1<-names(tab)[tab==1]
# 	nCluster1<-length(cluster1)
#   mClToTip1<-match(cluster1,phyloObj$tip.label) #find their tip number in matrix
# 	tipEdges1<-phyloObj$edge[phyloObj$edge[,2]==mClToTip1, ,drop=FALSE] #get the edge matrix for them
# 	mToPositiveCl<-match(tipEdges1[,2],cl[whPos])
# 	tipEdges1[,2]<-mToPositiveCl #change their tip number to their index in positive cluster vector
#
# }
# else{
# 	nCluster1<-0
# 	tipEdges1<-matrix(0,nrow=0,ncol=2)
# }
# if(any(tab>1)){
# 	clusterGr1<-names(tab)[tab>1]
# 	mGr1<-which(as.character(cl[whPos]) %in% clusterGr1)
#   newPhylo$edge<-cbind(cl[whPos][mGr1],mGr1)
#
#   mClToTip<-match(as.character(newPhylo$edge[,1]),phyloObj$tip.label) #gives for each cluster, tip number in old cluster tree
# 	if(any(is.na(mClToTip))) stop("coding error -- some of tip labels in phylo object do not match cluster names")
#   newPhylo$edge[,1]<-mClToTip
# 	newPhylo$edge<-rbind(newPhylo$edge,tipEdges1)
# }
# #---------
# ##Some tests -- might drop them for speed later
# #---------
# if(any(newPhylo$edge[,2]>length(cl[whPos]))) stop("coding error -- gave tip values greater than number of assigned samples")
# if(nrow(newPhylo$edge)!=length(cl[whPos])) stop("coding error -- did not assign all assigned samples to tip of tree")
# if(length(newPhylo$edge[,2])!=length(unique(newPhylo$edge[,2]))) stop("coding error -- did not give unique ids to tips of tree")
#
# #-------
# #make tips have values larger than root
# #-------
# internalEdges<-phyloObj$edge
# nOrigTips<-length(phyloObj$tip.label)
# maxOrig<-nOrigTips+phyloObj$Nnode
# whTipInternal<-internalEdges[,2]<=nOrigTips
# internalEdges[whTipInternal,2]<-internalEdges[whTipInternal,2]+maxOrig
# #fix up the newPhylo ones to have same number (again, for speed if fix to begin with?)
# whTipNew<-newPhylo$edge[,1]<=nOrigTips
# newPhylo$edge[whTipNew,1]<-newPhylo$edge[whTipNew,1]+maxOrig
#
#
# if(any(tab==1)){
# 	#have to remove those edges that went to clusters of size 1, because won't be internal edges.
# 	internalEdges<-internalEdges[ phyloObj$edge[,2]!=mClToTip1, ,drop=FALSE]
# }
# #also need to renumber everything, because the tip lost will leave gap... Need to keep the same order.
# oldValues<-sort( unique( as.numeric( internalEdges )))
# newValues<-1:length(oldValues)
# m1<-match(internalEdges[,1],oldValues)
# m2<-match(internalEdges[,2],oldValues)
# internalEdges<-cbind(newValues[m1],newValues[m2])
#
# #now fix up the edge matrix too; speed wise, is it faster to do it before? More complicated to code it before...
# newPhylo$edge[,1]<-newValues[match(newPhylo$edge[,1],oldValues)]
# internalEdges<-internalEdges+length(cl[whPos])
# newPhylo$edge[,1]<-newPhylo$edge[,1]+length(cl[whPos])
# newPhylo$edge<-rbind(internalEdges,newPhylo$edge)
#
# 		  if(!is.null(sampleNames)) newPhylo$tip.label<-sampleNames[whPos]
# 		  else{
# 		      newPhylo$tip.label<-paste("Sample",whPos)
# 		    }
# 		  newPhylo$Nnode<-phyloObj$Nnode+length(phyloObj$tip.label)-nCluster1
# if("edge.length" %in% names(phyloObj)){
# 	oldEdgeLength<-phyloObj$edge.length
# 	if( any(tab==1) ) oldEdgeLength<-oldEdgeLength[phyloObj$edge[,2]!=mClToTip1]
# 	newPhylo$edge.length<-rep(sampleEdgeLength,length=length(newPhylo$tip.label))
# 	if(any(tab==1)){ #need get back length to single-sample cluster
# 		oldEdgeLength1<-phyloObj$edge.length[phyloObj$edge[,2]==mClToTip1]
# 		newPhylo$edge.length[(length(newPhylo$edge.length)-nCluster1+1):length(newPhylo$edge.length)]<-oldEdgeLength1
# 	}
# 	newPhylo$edge.length<-c(oldEdgeLength,newPhylo$edge.length)
# }
# 		  class(newPhylo)<-"phylo"
