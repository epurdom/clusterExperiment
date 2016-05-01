.addBackSEInfo<-function(newObj,oldObj){
  retval<-clusterExperiment(oldObj,
                            clusters=clusterMatrix(newObj),
                            transformation=transformation(newObj),
                            clusterType=clusterType(newObj),
                            clusterInfo=clusterInfo(newObj))
  clusterLegend(retval)<-clusterLegend(newObj)
  orderSamples(retval)<-orderSamples(newObj)
  coClustering(retval)<-coClustering(newObj)
  return(retval)
}
.pullSampleData<-function(ce,wh){
  if(!is.null(wh)){
    sData<-colData(ce)
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
  else sData<-NULL
  return(sData)
}

.unnameClusterSlots<-function(ce){
    names(ce@clusterLegend)<-names(ce@clusterInfo)<-names(ce@clusterType)<-NULL
    return(ce)
}

.convertToNum<-function(x){
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
        return(x)
    }
    else return(x)
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

##Universal way to change character indication of clusterType into indices.
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
    if(!any(whClusters %in% clusterType(x))){
      #warning("none of indicated clusters match a clusterType")
      wh<-vector("integer",length=0)
    }
    else{
      #if(!all(whClusters %in% clusterType(x))) warning("not all indicated clusters match a clusterType")
      wh<-which(clusterType(x) %in% whClusters)
    }
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
	if(clusterFunction=="pam") type<-"K"
	return(type)
}

#wrapper that calls the clusterSampling and clusterD routines in reasonable order.
.clusterWrapper <- function(x, subsample, clusterFunction,clusterDArgs=NULL,
    subsampleArgs=NULL,typeAlg) 
{
	if(subsample){
		if(is.null(subsampleArgs) || !"k" %in% names(subsampleArgs)) stop("must provide k in 'subsampleArgs' (or if sequential should have been set by sequential strategy)")
		Dbar<-do.call("subsampleClustering",c(list(x=x),subsampleArgs))
		if(typeAlg=="K"){
			if(is.null(clusterDArgs)) clusterDArgs<-list(k=subsampleArgs[["k"]])
			else if(!"k" %in% names(clusterDArgs)) clusterDArgs[["k"]]<-subsampleArgs[["k"]] #either sequential sets this value, or get error in subsampleClustering, so always defined.
		}
    subDbar<-Dbar
	}
	else{ #create distance matrix if not subsampling
	  subDbar<-NULL
	}
  if(typeAlg=="K"){
    findBestK<-FALSE	
    if(!is.null(clusterDArgs) && "findBestK" %in% names(clusterDArgs)){
      findBestK<-clusterDArgs[["findBestK"]]
    }
    if(is.null(clusterDArgs) || (!"k" %in% names(clusterDArgs) && !findBestK)) stop("if not type 'K' algorithm, must give k in 'clusterDArgs' (or if sequential should have been set by sequential strategy)")
  }
	
	res<-do.call("clusterD",c(list(D=Dbar,format="list", clusterFunction=clusterFunction),clusterDArgs)) 
	return(list(results=res,subsampleCocluster=subDbar)) 
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
.makePhylobaseTree<-function(x,type){
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
    phylobase::nodeLabels(phylo4Obj)<-paste("Node",1:phylobase::nNodes(phylo4Obj),sep="")
    return(phylo4Obj)
}


