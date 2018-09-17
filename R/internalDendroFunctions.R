####
#Convert dendrogram class slots to class used by phylobase (phylo4) so can navigate easily. Does so by first converting to class of ape (phylo)

#' @importFrom phylobase edgeLength rootNode descendants nodeLabels
#' @importFrom dendextend as.phylo.dendrogram
#' @importFrom ape as.phylo
.makePhylobaseTree<-function(x,isSamples=FALSE,outbranch=FALSE, returnOnlyPhylo=FALSE){
  type<-class(x)#match.arg(type,c("hclust","dendro"))
  if(type=="hclust"){
    #first into phylo from ape package
    tempPhylo<-try(ape::as.phylo(x),FALSE)
    if(inherits(tempPhylo, "try-error")) stop("the hclust object cannot be converted to a phylo class with the methods of the 'ape' package. Reported error from ape package:",tempPhylo)
  }
  if(type=="dendrogram"){
    tempPhylo<-try(dendextend::as.phylo.dendrogram(x),FALSE)
    if(inherits(tempPhylo, "try-error")) stop(paste("the dendrogram object cannot be converted to a phylo class with the methods of 'dendextend' package. Check that you gave simple hierarchy of clusters, and not one with fake data per sample. Reported error from dendextend package:",tempPhylo))
  }
  if(type=="phylo") tempPhylo<-x

	  #put this in because some zero edges are becoming very small negative values...
  if(any(tempPhylo$edge.length<0)){
	  if(all(tempPhylo$edge.length> -1e-5)) tempPhylo$edge.length[tempPhylo$edge.length<0]<-0
	  else stop("coding error -- given object results in negative edge lengths")
  }
  if(returnOnlyPhylo ) return(tempPhylo) #don't do the rest of fixing up...

  phylo4Obj<-try(as(tempPhylo,"phylo4"),FALSE) 
  if(inherits(phylo4Obj, "try-error")) stop(paste("the internally created phylo object cannot be converted to a phylo4 class. Check that you gave simple hierarchy of clusters, and not one with fake data per sample. Reported error from dendextend package:",tempPhylo))
  
  if(isSamples){
	  ###Adds (in uniform way) the node names:
	  ### "Root"
	  ### "NodeX" where X is a number *to those NOT in outgroup*
	  ### "MissingNodeX" where X is a number *to those IN outgroup*
#	  browser()
    #NOTE: clusterNodes are found by those with non-zero edge-length between them and their decendents (i.e. fake hierarchy of the cluster should have edge length zero)
	#Moreover, they are distinguished from the -1/-2 outgroup by the fact that the tips of the -1/-2 DO HAVE non-zero edges to them. 
	#WHAT IF HAVE CLUSTER WITH ONE SAMPLE???? (just doesn't happen in practice, but often in our tests...)
    nonZeroEdges<-phylobase::edgeLength(phylo4Obj)[ which(phylobase::edgeLength(phylo4Obj)>0) ] #doesn't include root
    trueInternal<-sort(unique(as.numeric(sapply(strsplit(names(nonZeroEdges),"-"),.subset2,1)))) #this also picks up the outbranch between -1,-2
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

# clTree<-.makePhylobaseTree(clustWithDendro@dendro_clusters,"dendro")
# sampTree<-.makePhylobaseTree(clustWithDendro@dendro_samples,"dendro",isSamples=TRUE,outbranch=FALSE)

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




###Note, cluster nodes DEFINED as those whose descendants are of length>0. So root of these fake binary trees needs to be > 0
.makeFakeBinary<-function(tipNames,rootEdgeLength=1,edgeLength=0){
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
		mergeEdgeLength2<-depth1+mergeEdgeLength
	}
	newPhylo<-list()
	#3) Tree 1 inner nodes: add n2+1
	whinner1<-which(tree1$edge[,2] > n1)
	if(length(whinner1)>0) tree1$edge[whinner1,2]<-tree1$edge[whinner1,2]+n2+1
	tree1$edge[,1]<-tree1$edge[,1]+n2+1

	if(isTree){
		#1) Tree 2 tips: add n1
		whtip2<-which(tree2$edge[,2]<=n2)
		tree2$edge[whtip2,2]<-tree2$edge[whtip2,2]+n1
	
		#2) Tree 2 inner nodes: add n1+m1+1
		if(nrow(tree2$edge)!=length(whtip2)) tree2$edge[-whtip2,2]<-tree2$edge[-whtip2,2]+n1+m1+1
		tree2$edge[,1]<-tree2$edge[,1]+n1+m1+1
	
		#4) Add root edge (n2+n1+1) to (n1+n2+2) and (n2+m1+n2+2) to bottom of edge matrix:
		newPhylo$edge<-rbind(tree1$edge,tree2$edge)
		newPhylo$edge<-rbind(newPhylo$edge, rbind(c(n2+n1+1,n1+n2+2),c(n2+n1+1,n1+m1+n2+2)))
		#5) tip.label: c(tree1$tip.label,tree2$tip.label)
		newPhylo$tip.label<-c(tree1$tip.label,tree2$tip.label)
		#6) edge.length<-c(tree1$edge.length,tree2$edge.length,rep(mergeEdgeLength,2))
		#should there be check that there is edge length?
		
		newPhylo$edge.length<-c(tree1$edge.length,tree2$edge.length,c(mergeEdgeLength1,mergeEdgeLength2))
		newPhylo$Nnode<-m1+m2+1
	}
	else{
		#4) Add root edge (n2+n1+1) to tree 1 (n1+n2+2) and root to tip to bottom of edge matrix:
		newPhylo$edge<-tree1$edge
		newPhylo$edge<-rbind(newPhylo$edge, rbind(c(n2+n1+1,n1+n2+2),c(n2+n1+1,n1+n2)))
		
		#5) tip.label: c(tree1$tip.label,tree2$tip.label)
		newPhylo$tip.label<-c(tree1$tip.label,tree2)
		
		#6) 
		newPhylo$edge.length<-c(tree1$edge.length,c(mergeEdgeLength1,mergeEdgeLength2))
		newPhylo$Nnode<-m1+1
	}
	newPhylo$edge<-.makeIntegerMatrix(newPhylo$edge)
	class(newPhylo)<-"phylo"
	return(newPhylo)
	
}



#tipTrees should be list of trees in same order as tips in mainTree
#assumes rooted tree
#single values to be added to tree should have $edge=0-row matrix, $Nnode=0, $tip.label name of cluster, $edge.length=value to be added to length of cluster node
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
	whSingle<-match(which(isSingle),mainTree$edge[,2])
	if(sum(isSingle)>0){
		mainTree$edge[whSingle,2]<-cumN[isSingle]+1
	
	#Note that singletons, unlike others, short value of rootEdgeLength to their edge length, comparatively, when made with fakeBinaryTrees function. (and rootEdgeLength was added to deal with identification of cluster problem)
	#Need to add to those
		addValues<-sapply(tipTrees[isSingle],.subset2,"edge.length")
		mainTree$edge.length[whSingle]<-mainTree$edge.length[whSingle]+addValues
	}
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
	
	newPhylo<-list()
	newPhylo$edge<-do.call("rbind",c(lapply(tipTrees[!isSingle],.subset2,"edge"),list(mainTree$edge)))
	newPhylo$tip.label<-do.call("c",lapply(tipTrees,.subset2,"tip.label"))
	newPhylo$edge.length<-do.call("c",c(lapply(tipTrees[!isSingle],.subset2,"edge.length"),list(mainTree$edge.length)))
	newPhylo$Nnode<-mainTree$Nnode+sum(sapply(tipTrees,.subset2,"Nnode"))
	newPhylo$edge<-.makeIntegerMatrix(newPhylo$edge)
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

#' @importFrom phylobase nodeHeight tipLabels edgeLength edges edgeId
.force.ultrametric<-function(tree){
	##From http://blog.phytools.org/2017/03/forceultrametric-method-for-ultrametric.html
	if(!inherits(tree,"phylo4")) stop("tree must be of class phylo4")

	allTips<-phylobase::tipLabels(tree)
	depthToTips<-phylobase::nodeHeight(tree,allTips,from="root")
	maxD<-max(depthToTips)
	addValue<-maxD-depthToTips
	allLen<-phylobase::edgeLength(tree)
  edgeMat<-phylobase::edges(tree)
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