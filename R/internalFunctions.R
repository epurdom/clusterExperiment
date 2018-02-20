.mynote<-function(x){
	message(paste("Note:",x))
}

.eraseMerge<-function(x){
  x@merge_index<-NA_real_
  x@merge_dendrocluster_index<-NA_real_
  x@merge_method<-NA_character_
  x@merge_cutoff<-NA_real_
  x@merge_nodeProp<-NULL
  x@merge_nodeMerge<-NULL
  ch<-.checkMerge(x)      
  if(!is.logical(ch)) stop(ch)
  else return(x)
}
.addPrefixToClusterNames<-function(ceObj,prefix,whCluster){
    ceLegend<-clusterLegend(ceObj)[[whCluster]]
    whPos<-which(as.numeric(ceLegend[,"clusterIds"]) >0)
    if(length(whPos)>0) ceLegend[whPos,"name"]<-paste(prefix,ceLegend[whPos,"clusterIds"],sep="")
    clusterLegend(ceObj)[[whCluster]]<-ceLegend
    return(ceObj)
}

.addNewResult<-function(newObj,oldObj){
    retval<-addClusterings(newObj,oldObj) #want most recent addition on top of clusterMatrix
    #erases dendrogram from oldObj -- only keeps newObj -- so need to put it back if wasn't already there
    if(is.na(retval@dendro_index) & !is.na(newObj@dendro_index)) stop("Coding error -- addClusterings lost dendro_index")
    if(is.na(retval@merge_index) & !is.na(newObj@merge_index)) stop("Coding error -- addClusterings lost merge_index")
    if(is.na(retval@dendro_index) & !is.na(oldObj@dendro_index)){
      retval@dendro_samples<-oldObj@dendro_samples
      retval@dendro_clusters<-oldObj@dendro_clusters
      retval@dendro_outbranch<-oldObj@dendro_outbranch
      retval@dendro_index<-oldObj@dendro_index+nClusterings(newObj) #update index to where dendrogram from
    }
	if(is.na(retval@merge_index) & !is.na(oldObj@merge_index)){
      retval@merge_index<-oldObj@merge_index+nClusterings(newObj) #update index to where merge from
      retval@merge_nodeMerge<-oldObj@merge_nodeMerge
      retval@merge_cutoff<-oldObj@merge_cutoff
      retval@merge_method<-oldObj@merge_method
    }
    if(is.null(retval@merge_nodeProp) & !is.null(oldObj@merge_nodeProp)){
      retval@merge_nodeProp<-oldObj@merge_nodeProp
      retval@merge_dendrocluster_index<-oldObj@merge_dendrocluster_index+nClusterings(newObj) #update index to where merge from
    }
    #put back orderSamples, coClustering
    if(all(retval@orderSamples==1:nSamples(retval)) & !all(oldObj@orderSamples==1:nSamples(retval))) retval@orderSamples<-oldObj@orderSamples
    if(is.null(retval@coClustering)) retval@coClustering<-oldObj@coClustering
    retval<-.addBackSEInfo(newObj=retval,oldObj=oldObj) #make sure keeps SE info
#   Note: .addBackSEInfo calls ClusterExperiment (i.e. validates)
    return(retval)
}

.addBackSEInfo<-function(newObj,oldObj){
  retval<-ClusterExperiment(as(oldObj,"SingleCellExperiment"),
                            clusters=clusterMatrix(newObj),
                            transformation=transformation(newObj),
                            clusterTypes=clusterTypes(newObj),
                            clusterInfo=clusteringInfo(newObj),
                            orderSamples=orderSamples(newObj),
                            coClustering=coClustering(newObj),
                            dendro_samples=newObj@dendro_samples,
                            dendro_outbranch=newObj@dendro_outbranch,
                            dendro_clusters=newObj@dendro_clusters,
                            dendro_index=newObj@dendro_index,
                            merge_index=newObj@merge_index,
                            merge_cutoff=newObj@merge_cutoff,
                            merge_dendrocluster_index=newObj@merge_dendrocluster_index,
                            merge_nodeProp=newObj@merge_nodeProp,
                            merge_nodeMerge=newObj@merge_nodeMerge,
                            merge_method=newObj@merge_method,
                            primaryIndex=primaryClusterIndex(newObj),
							clusterLegend=clusterLegend(newObj),
                            checkTransformAndAssay=FALSE
  )
#  clusterLegend(retval)<-clusterLegend(newObj)
  return(retval)
}
#Returns NULL if no sample data
.pullSampleData<-function(ce,wh,fixNA=c("keepNA","unassigned","missing")){
	fixNA<-match.arg(fixNA)
  if(!is.null(wh)){
    sData<-colData(ce)
	if(!is.logical(wh)){

	    if(NCOL(sData)==0) stop("no colData for object data, so cannot pull sampleData")
	    if(is.character(wh)){
	      if(all(wh=="all")) wh<-1:NCOL(sData)
	      else{
	        if(!all(wh %in% colnames(sData))) stop("Invalid names for pulling sampleData (some do not match names of colData)")
	        else wh<-match(wh,colnames(sData))
	      }
	    }
	    else if(is.numeric(wh)){
	      if(!all(wh %in% 1:NCOL(sData))) stop("Invalid indices for for pulling sampleData (some indices are not in 1:NCOL(colData)")
	    }
	    else stop("invalid values for pulling sampleData from colData of object")
		sData<-as.data.frame(sData[,wh,drop=FALSE])
	}
	else{ #if 
		if(wh) sData<- colData(ce)
		else sData<-NULL
	}
}
  else sData<-NULL
  if(!is.null(sData) && fixNA!="keepNA"){
  		newValue<-switch(fixNA,"unassigned"=-1,"missing"=-2)
  		fixNAFunction<-function(x,newValue){
  			if(is.factor(x)){ #change to character
  				waslevels<-levels(x)
  				wasFactor<-TRUE
  				x<-as.character(x)
  			}
  			else wasFactor<-FALSE
  			if(is.character(x)){
  				x[which(is.na(x))]<-as.character(newValue)
  				if(wasFactor) x<-factor(x,levels=c(waslevels,as.character(newValue))) #keeps order of previous factors
  			}
  			else x[which(is.na(x))]<- newValue #assume numeric if not character/factor
  			return(x)
  		}
  		#have to do this; otherwise makes them all characters if use apply...
  		cnames<-colnames(sData)
  		sData<-do.call("data.frame",lapply(1:ncol(sData),function(ii){fixNAFunction(sData[,ii],newValue=newValue)}))
  		colnames(sData)<-cnames
    	
      }
  return(sData)
}

.unnameClusterSlots<-function(ce){
    names(ce@clusterLegend)<-names(ce@clusterInfo)<-names(ce@clusterTypes)<-NULL
    return(ce)
}

.convertToNum<-function(x){
	nms<-names(x)
    if(is.factor(x)){
        x<-as.character(x)
    }
    if(is.character(x)){
        op <- options(warn=2)
        #if character values convert to numeric, will use that, otherwise convert to factor first
        test <- try(as.numeric(x) ,silent=TRUE)
        if(inherits(test,"try-error")) x<-as.numeric(factor(x))
        else x<-test
        options(op)
        
    }
	names(x)<-nms
	return(x)
}
##Universal way to convert matrix of clusters (of any value) into integers, preserving -1, -2 values
.makeIntegerClusters<-function(clMat){
    if(!is.matrix(clMat)) stop("must give matrix input")
    fun<-function(x){ #make it numbers from 1:length(x), except for -1,-2
      #id special values of -1,-2
        isChar<-as.character(is.character(x))
        wh1<-switch(isChar,"TRUE"=x =="-1","FALSE"=x ==-1)
        wh2<-switch(isChar,"TRUE"=x =="-2","FALSE"=x ==-2)
      wh<-wh1 | wh2
      vals<-unique(x[!wh])
      y<-match(x,vals)
      y[wh1]<- -1
      y[wh2]<- -2
      return(y)
    }
    
    if(!is.null(dim(clMat)) && ncol(clMat)>1){
        x<-apply(clMat,2,fun)
        if(is.null(dim(x))) x<-matrix(x,nrow=1) #in case clMat was matrix with 1 row
        return(x  )
    } 
    else{
        if(is.matrix(clMat)) clMat<-clMat[,1]
        return(matrix(fun(clMat),ncol=1))
    }
}
##Universal way to convert matrix of clusters into colors
.makeColors<-function(clMat,colors,unassignedColor="white",missingColor="grey",makeIntegers=TRUE){
    if(any(apply(clMat,2,function(x){length(unique(x))})>length(colors))) warning("too many clusters to have unique color assignments")
    if(any(apply(clMat,2,function(x){any(is.na(x))}))) stop("clusters should not have 'NA' values; non-clustered samples should get a '-1' or '-2' value depending on why they are not clustered.")
    cNames<-colnames(clMat)
    origClMat<-clMat
    if(makeIntegers) clMat<-.makeIntegerClusters(clMat) #don't use when call from some plots where very carefully already chosen
    
    if(ncol(clMat)>1){
        colorMat<-apply(clMat,2,function(x){
            y<-vector("character",length(x))
            currcolors<-rep(colors,length=length(unique(x[x>0]))) #just duplicate colors if more than in existing palate
            y[x>0]<-currcolors[x[x>0]]
            return(y)
        })
		if(nrow(clMat)==1) colorMat<-matrix(colorMat,nrow=1) #in case only 1 sample!
    }
    else{
        if(is.matrix(clMat)) x<-clMat[,1] else x<-clMat
        y<-vector("character",length(x))
        currcolors<-rep(colors,length=length(unique(x[x>0]))) #just duplicate colors if more than in existing palate
        y[x>=0]<-currcolors[x[x>=0]]
        colorMat<-matrix(y,ncol=1)
    }
    colorMat[clMat== -1]<-unassignedColor
    colorMat[clMat== -2]<-missingColor

    #convert ids into list of matrices:
    colorList<-lapply(1:ncol(clMat),function(ii){
mat<-unique(cbind("clusterIds"=clMat[,ii],"color"=colorMat[,ii],"name"=origClMat[,ii]))
        rownames(mat)<-NULL
        return(mat)
    })
    names(colorList)<-cNames
    colnames(colorMat)<-cNames
    return(list(colorList=colorList,convertedToColor=colorMat,numClusters=clMat))
}

##Universal way to change character indication of clusterTypes into integer indices.
##If no match, returns vector length 0
.TypeIntoIndices<-function(x,whClusters){
 if(is.numeric(whClusters)) wh<-whClusters
	 else{
		 test<-try(match.arg(whClusters[1],c("workflow","all","none","primaryCluster","dendro")),silent=TRUE)
		   if(!inherits(test,"try-error")){
		     if(test=="workflow"){
		       ppIndex<-workflowClusterDetails(x)
		       if(!is.null(ppIndex) && sum(ppIndex[,"iteration"]==0)>0){
		         wh<-unlist(lapply(.workflowValues,function(tt){
		           ppIndex[ppIndex[,"iteration"]==0 & ppIndex[,"type"]==tt,"index"]
		         }))
		       }
		       else wh<-vector("integer",length=0)
		     }
		     if(test=="all"){
		       #put primary cluster first
		       ppcl<-primaryClusterIndex(x)
		       wh<-c(ppcl,c(1:nClusterings(x))[-ppcl])
		     }
		     if(test=="none") wh<-vector("integer",length=0)
		     if(test=="primaryCluster") wh<-primaryClusterIndex(x)
		 	if(test=="dendro"){
		 		wh<-dendroClusterIndex(x)
		 		if(is.na(wh)) wh<-vector("integer",length=0)
		 	}
		   }
		   else{
		     #first match to clusterTypes  
		     mClType<-match(whClusters,clusterTypes(x))  
		     mClLabel<-match(whClusters,clusterLabels(x))  
		     totalMatch<-mapply(whClusters,mClType,mClLabel,FUN=function(cl,type,lab){
		         if(is.na(type) & !is.na(lab)) return(lab)
		         if(is.na(type) & is.na(lab)) return(NA)
		         if(!is.na(type)){
		             return(which(clusterTypes(x) %in% cl)) #prioritize clusterType and get ALL of them, not just first match
		         }
		     },SIMPLIFY=FALSE)
		     totalMatch<-unlist(totalMatch,use.names=FALSE)
    
		     if(all(is.na(totalMatch))) wh<-vector("integer",length=0)
		     else wh<-na.omit(totalMatch) #silently ignore things that don't match.
		   }
	 } 
  	 if(any(wh>nClusterings(x) | wh<1)){
		 wh<-wh[wh<=nClusterings(x) & wh>0]
		 
  	 }
#	 if(length(wh)>0) wh<-wh[is.integer(wh)]
  return(wh)
}



####
#Convert to object used by phylobase so can navigate easily 
.makePhylobaseTree<-function(x,type,isSamples=FALSE,outbranch=FALSE){
    type<-match.arg(type,c("hclust","dendro"))
    if(type=="hclust"){
        #first into phylo from ape package
        tempPhylo<-try(ape::as.phylo(x),FALSE)
        if(inherits(tempPhylo, "try-error")) stop("the hclust object cannot be converted to a phylo class with the methods of the 'ape' package.")
    }
    if(type=="dendro"){
        tempPhylo<-try(dendextend::as.phylo.dendrogram(x),FALSE)
        if(inherits(tempPhylo, "try-error")) stop("the dendrogram object cannot be converted to a phylo class with the methods of 'dendextend' package. Check that you gave simple hierarchy of clusters, and not one with fake data per sample")
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
		
		phylobase::nodeLabels(phylo4Obj)[as.character(trueInternal)]<-paste("Node",1:length(trueInternal),sep="")
		#add new label for root 
		if(outbranch){
			phylobase::nodeLabels(phylo4Obj)[as.character(rootNode)]<-"Root"
			if(outbranchIsInternal) phylobase::nodeLabels(phylo4Obj)[as.character(outbranchNodeDesc)]<-paste("MissingNode",1:length(outbranchNodeDesc),sep="")
		}
	}
	else phylobase::nodeLabels(phylo4Obj)<-paste("Node",1:phylobase::nNodes(phylo4Obj),sep="")
    
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