.checkXDissInput<-function(x,diss){
  if(is.null(x) & is.null(diss)) stop("must give either x or diss argument")
  #  if(!is.null(x) & !is.null(diss)) stop("cannot give both x and diss argument")
  if(!is.null(x) & is.null(diss)) input<-"X"
  if(!is.null(x) & !is.null(diss)) input<-"both"
  if(is.null(x) & !is.null(diss)) input<-"diss"
  if(input %in% c("diss","both")) .checkDistFunction(diss)
  if(input == "both" && ncol(x)!=ncol(diss)) stop("ncol(x)!=ncol(diss): if both x and diss then must have compatible dimensions.") 
  return(input)
}
.checkDistFunction<-function(D){
  if(any(is.na(as.vector(D)))) stop("NA values found in D (could be from too small of subsampling if classifyMethod!='All', see documentation of subsampleClustering)")
  if(any(is.na(D) | is.nan(D) | is.infinite(D))) stop("D matrix contains either NAs, NANs or Infinite values.")
  if(any(D<0)) stop("distance function must give strictly positive values")
  if(any(diag(D)!=0)) stop("distance function must have zero values on the diagonal of the distance matrix")
}

.addPrefixToClusterNames<-function(ceObj,prefix,whCluster){
    ceLegend<-clusterLegend(ceObj)[[whCluster]]
    whPos<-which(ceLegend[,"clusterIds"] >0)
    if(length(whPos)>0) ceLegend[whPos,"name"]<-paste(prefix,ceLegend[whPos,"clusterIds"],sep="")
    clusterLegend(ceObj)[[whCluster]]<-ceLegend
    return(ceObj)
}

.addNewResult<-function(newObj,oldObj){
    retval<-addClusters(newObj,oldObj) #want most recent addition on top of clusterMatrix
    #erases dendrogram so need to put it back if wasn't already there
    if(is.na(retval@dendro_index) & !is.na(oldObj@dendro_index)){
        retval@dendro_samples<-oldObj@dendro_samples
        retval@dendro_clusters<-oldObj@dendro_clusters
        retval@dendro_index<-oldObj@dendro_index+nClusters(newObj) #update index to where dendrogram from
    }
    #put back orderSamples, coClustering
    if(all(retval@orderSamples==1:nSamples(retval)) & !all(oldObj@orderSamples==1:nSamples(retval))) retval@orderSamples<-oldObj@orderSamples
    if(is.null(retval@coClustering)) retval@coClustering<-oldObj@coClustering
    retval<-.addBackSEInfo(newObj=retval,oldObj=oldObj) #make sure keeps SE info
    validObject(retval)
    return(retval)
}

.addBackSEInfo<-function(newObj,oldObj){
  retval<-clusterExperiment(oldObj,
                            clusters=clusterMatrix(newObj),
                            transformation=transformation(newObj),
                            clusterTypes=clusterTypes(newObj),
                            clusterInfo=clusterInfo(newObj),
                            orderSamples=orderSamples(newObj),
                            coClustering=coClustering(newObj),
                            dendro_samples=newObj@dendro_samples,
                            dendro_clusters=newObj@dendro_clusters,
                            dendro_index=newObj@dendro_index)
  clusterLegend(retval)<-clusterLegend(newObj)
  return(retval)
}
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
    #browser()
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
#         col<-colorMat[,ii]
#         ids<-clMat[,ii]
#         origids<-origClMat[,ii]
#         uniqueIds<-unique(ids)
#         mIds<-match(uniqueIds,ids)
#         uniqueCols<-col[mIds]
#         mat<-cbind("clusterIds"=uniqueIds,"color"=uniqueCols)
        mat<-unique(cbind("clusterIds"=clMat[,ii],"color"=colorMat[,ii],"name"=origClMat[,ii]))
        rownames(mat)<-mat[,"clusterIds"]
        return(mat)
    })
    names(colorList)<-cNames
    colnames(colorMat)<-cNames
    return(list(colorList=colorList,convertedToColor=colorMat,numClusters=clMat))
}

##Universal way to change character indication of clusterTypes into indices.
.TypeIntoIndices<-function(x,whClusters){
  test<-try(match.arg(whClusters[1],c("workflow","all","none","primaryCluster")),silent=TRUE)
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
      wh<-c(ppcl,c(1:nClusters(x))[-ppcl])
    }
    if(test=="none") wh<-vector("integer",length=0)
    if(test=="primaryCluster") wh<-primaryClusterIndex(x)
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
    #browser()
    if(all(is.na(totalMatch))) wh<-vector("integer",length=0)
    else wh<-na.omit(totalMatch) #silently ignore things that don't match.
#     
#         if(!any(whClusters %in% clusterTypes(x))){
#         if(!any(whClusters %in% clusterLabels(x))) wh<-vector("integer",length=0)
#         else{
#             wh<-which(clusterLabels(x) %in% whClusters)
#         }
#     }
#     else{
#       #if(!all(whClusters %in% clusterTypes(x))) warning("not all indicated clusters match a clusterTypes")
#       wh<-which(clusterTypes(x) %in% whClusters)
#     }
  }
  return(wh)
}

#######
#Internal algorithms for clustering
#######
#check what type
.checkAlgType<-function(clusterFunction){
	##These return lists of indices of clusters satisifying alpha criteria
	if(clusterFunction=="tight") type<-"01"
	if(clusterFunction=="hierarchical01") type<-"01"
	if(clusterFunction=="hierarchicalK") type<-"K"
	if(clusterFunction=="pam") type<-"K"
	return(type)
}



#convert list output into cluster vector.
.convertClusterListToVector<-function(clusterList,N)
{
    clust.id <- rep(-1, N)
    nfound<-length(clusterList)
    if(nfound>0){
        #make cluster ids in order of when found
        for (i in 1:length(clusterList)) clust.id[clusterList[[i]]] <- i 
    }
    return(clust.id)
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
	#browser()
	if(isSamples){
		#NOTE: clusterNodes are found by those with non-zero edge-length between them and their decendents
		nonZeroEdges<-edgeLength(phylo4Obj)[which(edgeLength(phylo4Obj)>0)] #doesn't include root
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
			#######
			rootChild<-phylobase::descendants(phylo4Obj,node=rootNode,type="children")
			#find node descendants of these:
			rootChildDesc<-lapply(rootChild,phylobase::descendants,phy=phylo4Obj,type="all")
			rootChildNum<-sapply(rootChildDesc,function(x){length(x[x%in%trueInternal])})
			outbranchNode<-rootChild[rootChildNum<=1]
			if(outbranchNode %in% trueInternal){
				outbranchIsInternal<-TRUE
				trueInternal<-trueInternal[!trueInternal%in%outbranchNode]
			}
			
		}
		#trueInternal<-allInternal[!allInternal%in%clusterNodes]
		
		phylobase::nodeLabels(phylo4Obj)[as.character(trueInternal)]<-paste("Node",1:length(trueInternal),sep="")
		#add new label for root 
		if(outbranch){
			phylobase::nodeLabels(phylo4Obj)[as.character(rootNode)]<-"Root"
			if(outbranchIsInternal) phylobase::nodeLabels(phylo4Obj)[as.character(outbranchNode)]<-"MissingSamples"
		}
	}
	else phylobase::nodeLabels(phylo4Obj)<-paste("Node",1:phylobase::nNodes(phylo4Obj),sep="")
    
	return(phylo4Obj)
}

# clTree<-.makePhylobaseTree(clustWithDendro@dendro_clusters,"dendro")
# sampTree<-.makePhylobaseTree(clustWithDendro@dendro_samples,"dendro",isSamples=TRUE,outbranch=FALSE)

