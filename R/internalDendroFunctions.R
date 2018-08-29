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

###Potential problems:
## What if only 1 sample in a cluster? Need to make it take the place of the cluster node...
.makeSampleDendro<-function(ceObj,unassignedSamples=c("ignore","outgroup", "cluster")){
	unassignedSamples<-match.arg(unassignedSamples)
  cl<-ceObj@clusterMatrix[,ceObj@dendro_index]
  whPos<-which(cl>0)
  clPos<-cl[whPos]
    
  #Note, perhaps should make a "phylo" only option in .makePhylobase so don't do this twice.
  #Would loose internal node names. Need check if that matters.
  phyloObj <- .makePhylobaseTree(x=ceObj@dendro_clusters, "dendro",isSamples=FALSE,returnOnlyPhylo = TRUE)

  #make edge matrix per cluster
	nClusterGr1<-
  tipEdges<-cbind(clPos,1:length(clPos))
  internalEdges<- phyloObj$edge+length(clPos)
  mClToTip<-match(as.character(clPos),phyloObj$tip.label) #gives for each cluster, tip number in cluster tree
	if(any(is.na(mClToTip))) stop("coding error -- some of tip labels in phylo object do not match cluster names")
  tipEdges[,1]<-mClToTip+length(clPos)
  newPhylo<-list()
  newPhylo$edge<-rbind(internalEdges,tipEdges)
  if(!is.null(colnames(ceObj))) newPhylo$tip.label<-colnames(ceObj)[whPos]
  else{
      newPhylo$tip.label<-paste("Sample",1:NCOL(ceObj))[whPos]
    }
  newPhylo$Nnode<-phyloObj$Nnode+length(phyloObj$tip.label)
  class(newPhylo)<-"phylo"
  return(newPhylo)

}

