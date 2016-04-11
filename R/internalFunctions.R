##Universal way to change character indication of clusterType into indices.
.TypeIntoIndices<-function(x,whClusters=c("pipeline","all")){
  test<-try(match.arg(whClusters),silent=TRUE)
  if(!inherits(test,"try-error")){
    if(test=="pipeline"){
      ppIndex<-pipelineClusterDetails(x)
      if(!is.null(ppIndex) && sum(ppIndex[,"iteration"]==0)>0){
        wh<-unlist(lapply(.pipelineValues,function(tt){
          ppIndex[ppIndex[,"iteration"]==0 & ppIndex[,"type"]==tt,"index"]
        }))
      }
      else stop("There are no (current) pipeline clusters in the ClusterExperiments object")
    }
    if(test=="all") wh<-1:ncol(allClusters(x))
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
      newX<-removeClusters(x,ppIndex[,"index"]) ###Getting error: Error: evaluation nested too deeply: infinite recursion / options(expressions=)?
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




