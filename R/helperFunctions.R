
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




