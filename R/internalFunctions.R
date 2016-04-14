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
    if(!is.null(dim(clMat)) && ncol(clMat)>1){
        return(apply(clMat,2,function(x){ #make it numbers from 1:length(x), except for -1,-2
            vals<-unique(x[!x<0])
            y<-match(x,vals)
            y[x%in%c(-1,-2)]<-x[x<0]
            return(y)
        })  )
    } 
    else{
        if(is.matrix(clMat)) clMat<-clMat[,1]
        vals<-unique(clMat[clMat>=0])
        y<-match(clMat,vals)
        y[clMat%in%c(-1,-2)]<-clMat[clMat<0]
        return(matrix(y,ncol=1))
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
        mat<-unique(cbind("clusterIds"=origClMat[,ii],"color"=colorMat[,ii],"name"=origClMat[,ii]))
        rownames(mat)<-mat[,"clusterIds"]
        return(mat)
    })
    names(colorList)<-cNames
    colnames(colorMat)<-cNames
    return(list(colorList=colorList,convertedToColor=colorMat))
}

##Universal way to change character indication of clusterType into indices.
.TypeIntoIndices<-function(x,whClusters){
  test<-try(match.arg(whClusters[1],c("pipeline","all","none","primary")),silent=TRUE)
  if(!inherits(test,"try-error")){
    if(test=="pipeline"){
      ppIndex<-pipelineClusterDetails(x)
      if(!is.null(ppIndex) && sum(ppIndex[,"iteration"]==0)>0){
        wh<-unlist(lapply(.pipelineValues,function(tt){
          ppIndex[ppIndex[,"iteration"]==0 & ppIndex[,"type"]==tt,"index"]
        }))
      }
      else wh<-vector("integer",length=0)
    }
    if(test=="all") wh<-1:ncol(clusterMatrix(x))
    if(test=="none") wh<-vector("integer",length=0)
    if(test=="primary") wh<-primaryClusterIndex(x)
  }
  else{
    if(!any(whClusters %in% clusterType(x))) stop("none of orderClusters match a clusterType of x")
    if(!all(whClusters %in% clusterType(x))) warning("not all of orderClusters match a clusterType of x")
    wh<-which(clusterType(x) %in% whClusters)
  }
  return(wh)
}

#change current pipeline to old iteration 
# add number to it if eraseOld=FALSE
# delete ALL pipeline if eraseOld=TRUE (not just the current iteration)
.updateCurrentPipeline<-function(x,eraseOld){
  ppIndex<-pipelineClusterDetails(x)
  if(!is.null(ppIndex)){ #need to change the clusterType values (or erase them) before get new ones
    if(eraseOld){ #removes all of them, not just current
      newX<-removeClusters(x,ppIndex[,"index"]) 
    }
    else{
      if(0 %in% ppIndex[,"iteration"]){
        newIteration<-max(ppIndex[,"iteration"])+1
        whCurrent<-ppIndex[ppIndex[,"iteration"]==0,"index"]
        updateCluster<-clusterType(x)
        updateCluster[whCurrent]<-paste(updateCluster[whCurrent],newIteration,sep="_")
        newX<-x
        newX@clusterType<-updateCluster          
      }
    }
    
  }
  else newX<-x
  validObject(newX)
  return(newX)
}

#check what type
.checkAlgType<-function(clusterFunction){
	##These return lists of indices of clusters satisifying alpha criteria
	if(clusterFunction=="tight") type<-"01"
	if(clusterFunction=="hierarchical") type<-"01"
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
	else{
		if(typeAlg!="K") stop("currently, if not subsampling, must use 'pam' or a clusterFunction defined as typeAlg='K' as clusterMethod")
		Dbar<-as.matrix(dist(x)	)	
		findBestK<-FALSE	
		if(!is.null(clusterDArgs) && "findBestK" %in% names(clusterDArgs)){
				findBestK<-clusterDArgs[["findBestK"]]
			}
		if(is.null(clusterDArgs) || (!"k" %in% names(clusterDArgs) && !findBestK)) stop("if not subsampling, must give k in 'clusterDArgs' (or if sequential should have been set by sequential strategy)")
	  subDbar<-NULL
	}
	if(any(is.na(as.vector(Dbar)))) stop("NA values found in Dbar (could be from too small of subsampling if classifyMethod!='All', see documentation of subsampleClustering)")
	
	res<-do.call("clusterD",c(list(D=Dbar,format="list", clusterFunction=clusterFunction),clusterDArgs)) 
	return(list(results=res,subsampleCocluster=subDbar)) #nothing found
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




