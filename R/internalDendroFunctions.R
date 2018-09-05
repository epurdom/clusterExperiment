####
#Convert dendrogram class slots to class used by phylobase (phylo4) so can navigate easily. Does so by first converting to class of ape (phylo)

#' @importFrom phylobase edgeLength rootNode descendants nodeLabels
#' @importFrom dendextend as.phylo.dendrogram
#' @importFrom ape as.phylo
.makePhylobaseTree<-function(x,type,isSamples=FALSE,outbranch=FALSE, returnOnlyPhylo=FALSE){
  type<-match.arg(type,c("hclust","dendro"))
  if(type=="hclust"){
    #first into phylo from ape package
    tempPhylo<-try(ape::as.phylo(x),FALSE)
    if(inherits(tempPhylo, "try-error")) stop("the hclust object cannot be converted to a phylo class with the methods of the 'ape' package.")
  }
  if(type=="dendro"){
    tempPhylo<-try(dendextend::as.phylo.dendrogram(x),FALSE)
    if(inherits(tempPhylo, "try-error")) stop(paste("the dendrogram object cannot be converted to a phylo class with the methods of 'dendextend' package. Check that you gave simple hierarchy of clusters, and not one with fake data per sample. Reported error from dendextend package:",tempPhylo))
  }
  if(returnOnlyPhylo) return(tempPhylo) #don't do the rest of fixing up...
  phylo4Obj<-try(as(tempPhylo,"phylo4"),FALSE) 
  if(inherits(phylo4Obj, "try-error")) stop(paste("the internally created phylo object cannot be converted to a phylo4 class. Check that you gave simple hierarchy of clusters, and not one with fake data per sample. Reported error from dendextend package:",tempPhylo))
  
  if(isSamples){
    #NOTE: clusterNodes are found by those with non-zero edge-length between them and their decendents
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


#Note, in long term, want to be able to subset the CE object without tossing the dendrogram! How does this work if you completely lose a cluster in a subset?

## Large things it makes:
# 2 vectors roughly size of n
# phylo tree: n(n-1)/2 x 2 matrix (edge matrix), n length vector of characters (tip labels), n(n-1)/2 length vector (edge lengths)
.makeSampleDendro<-function(clusterDendro, cl, unassignedSamples=c("remove","outgroup"),sampleEdgeLength=0, sampleNames=NULL, outbranchLength,outlierHclust=NULL){
	unassignedSamples<-match.arg(unassignedSamples)
	#right now outlierHclust is expected to be a hclust object.
	#need check that outlierHclust is right number of observations!
	#if no outbranchHist given, assume only a few observations (e.g. 1!) and do arbitrary order from single node, like a cluster.
	
	
  whPos<-which(cl>0) #this is copy close to length of n
	if(!is.null(sampleNames) && length(sampleNames)!=length(cl)) stop("sampleNames must be same length as cluster vector")
  #loses internal node names. Don't think that matters.
  phyloObj <- .makePhylobaseTree(x=clusterDendro, "dendro",isSamples=FALSE,returnOnlyPhylo = TRUE)

	#########
  #make new edge matrix
	#########
  tab<-table(cl[whPos])
  newPhylo<-list()
	
  if(any(tab==1)){
		cluster1<-names(tab)[tab==1]
		nCluster1<-length(cluster1)
	  mClToTip1<-match(cluster1,phyloObj$tip.label) #find their tip number in matrix
		tipEdges1<-phyloObj$edge[phyloObj$edge[,2]==mClToTip1, ,drop=FALSE] #get the edge matrix for them
		mToPositiveCl<-match(tipEdges1[,2],cl[whPos])
		tipEdges1[,2]<-mToPositiveCl #change their tip number to their index in positive cluster vector
		
	}
	else{
		nCluster1<-0
		tipEdges1<-matrix(0,nrow=0,ncol=2)
	}
	if(any(tab>1)){
		clusterGr1<-names(tab)[tab>1]
		mGr1<-which(as.character(cl[whPos]) %in% clusterGr1)
	  newPhylo$edge<-cbind(cl[whPos][mGr1],mGr1)
		
	  mClToTip<-match(as.character(newPhylo$edge[,1]),phyloObj$tip.label) #gives for each cluster, tip number in old cluster tree
		if(any(is.na(mClToTip))) stop("coding error -- some of tip labels in phylo object do not match cluster names")
	  newPhylo$edge[,1]<-mClToTip
		newPhylo$edge<-rbind(newPhylo$edge,tipEdges1)		
	}
	
	#########
	##Some tests -- might drop them for speed later
	#########
	if(any(newPhylo$edge[,2]>length(cl[whPos]))) stop("coding error -- gave tip values greater than number of assigned samples")
	if(nrow(newPhylo$edge)!=length(cl[whPos])) stop("coding error -- did not assign all assigned samples to tip of tree")
	if(length(newPhylo$edge[,2])!=length(unique(newPhylo$edge[,2]))) stop("coding error -- did not give unique ids to tips of tree")	
	
	#-------
	#make tips have values larger than root
	#-------
	internalEdges<-phyloObj$edge
	nOrigTips<-length(phyloObj$tip.label)
	maxOrig<-nOrigTips+phyloObj$Nnode
	whTipInternal<-internalEdges[,2]<=nOrigTips
	internalEdges[whTipInternal,2]<-internalEdges[whTipInternal,2]+maxOrig
	#fix up the newPhylo ones to have same number (again, for speed if fix to begin with?)
	whTipNew<-newPhylo$edge[,1]<=nOrigTips
	newPhylo$edge[whTipNew,1]<-newPhylo$edge[whTipNew,1]+maxOrig


	if(any(tab==1)){
		#have to remove those edges that went to clusters of size 1, because won't be internal edges.
		internalEdges<-internalEdges[ phyloObj$edge[,2]!=mClToTip1, ,drop=FALSE]		
	}
	#also need to renumber everything, because the tip lost will leave gap... Need to keep the same order. 
	oldValues<-sort( unique( as.numeric( internalEdges )))
	newValues<-1:length(oldValues)
	m1<-match(internalEdges[,1],oldValues)
	m2<-match(internalEdges[,2],oldValues)
	internalEdges<-cbind(newValues[m1],newValues[m2])
	#now fix up the edge matrix too; speed wise, is it faster to do it before? More complicated to code...
	newPhylo$edge[,1]<-newValues[match(newPhylo$edge[,1],oldValues)]

	internalEdges<-internalEdges+length(cl[whPos])
	newPhylo$edge[,1]<-newPhylo$edge[,1]+length(cl[whPos])
	
	newPhylo$edge<-rbind(internalEdges,newPhylo$edge)
	
  if(!is.null(sampleNames)) newPhylo$tip.label<-sampleNames[whPos]
  else{
      newPhylo$tip.label<-paste("Sample",1:length(cl))[whPos]
    }
  newPhylo$Nnode<-phyloObj$Nnode+length(phyloObj$tip.label)-nCluster1

	if("edge.length" %in% names(phyloObj)){
		oldEdgeLength<-phyloObj$edge.length
		if( any(tab==1) ) oldEdgeLength<-oldEdgeLength[phyloObj$edge[,2]!=mClToTip1]
		newPhylo$edge.length<-rep(sampleEdgeLength,length=length(newPhylo$tip.label))
		if(any(tab==1)){ #need get back length to single-sample cluster
			oldEdgeLength1<-phyloObj$edge.length[phyloObj$edge[,2]==mClToTip1]
			newPhylo$edge.length[(length(newPhylo$edge.length)-nCluster1+1):length(newPhylo$edge.length)]<-oldEdgeLength1
		}
		newPhylo$edge.length<-c(oldEdgeLength,newPhylo$edge.length)
	}
	
	
	if(unassignedSamples=="outgroup"){
		#make a tree out of the unassigned samples with hclust:
		
		
		
	}
		
		
	
	
  class(newPhylo)<-"phylo"
  return(newPhylo)

}

#only for merging actual treeds
.mergePhylo<-function(tree1,tree2,mergeEdgeLength){
	n1<-length(tree1$tip.label)
	n2<-length(tree2$tip.label)
	m1<-tree1$Nnode
	m2<-tree2$Nnode
	
	newPhylo<-list()
	#1) Tree 2 tips: add n1
	whtip2<-tree2$edge[,2]<=n2
	tree2$edge[whtip2,2]<-tree2$edge[whtip2,2]+n1
	#2) Tree 2 inner nodes: add n1+m1+1
	tree2$edge[-whtip2,2]<-tree2$edge[-whtip2,2]+n1+m1+1
	tree2$edge[,1]<-tree2$edge[,1]+n1+m1+1
	
	#3) Tree 1 inner nodes: add n1
	whinner1<-tree1$edge[,2]>=n1
	tree1$edge[whinner1,2]<-tree1$edge[whinner1,2]+n1
	tree1$edge[,1]<-tree1$edge[,1]+n1
	
	#4) Add root edge (n2+n1+1) to (n1+n2+2) and (n2+m1+n2+2) to bottom of edge matrix:
	newPhylo$edge<-rbind(tree1$edge,tree2$edge)
	newPhylo$edge<-rbind(newPhylo$edge, rbind(c(n2+n1+1,n1+n2+2),c(n2+n1+1,n2+m1+n2+2)))

	#5) tip.label: c(tree1$tip.label,tree2$tip.label)
	newPhylo$tip.label<-c(tree1$tip.label,tree2$tip.label)
	
	#6) edge.length<-c(tree1$edge.length,tree2$edge.length,rep(mergeEdgeLength,2))
	#should there be check that there is edge length?
	if(is.null(tree1$edge.length) || is.null(tree1$edge.length) ) stop("must have edge length")
	newPhylo$edge.length<-c(tree1$edge.length,tree2$edge.length,rep(mergeEdgeLength,2))
	return(newPhylo)
	
}

.mergeTip2Phylo<-function(tree1, tipName, mergeEdgeLength){

}