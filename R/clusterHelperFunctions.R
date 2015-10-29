#check what type
.checkAlgType<-function(clusterFunction){
	##These return lists of indices of clusters satisifying alpha criteria
	if(clusterFunction=="tight") type<-"01"
	if(clusterFunction=="hierarchical") type<-"01"
	if(clusterFunction=="pam") type<-"K"
	return(type)
}

#wrapper that calls the clusterSampling and clusterD routines in reasonable order.
.clusterWrapper <- function(x, subsample, clusterFunction,clusterArgs=NULL,
    subsampleArgs=NULL,typeAlg) 
{
	if(subsample){
		if(is.null(subsampleArgs) || !"k" %in% names(subsampleArgs)) stop("must provide k in 'subsampleArgs' (or if sequential should have been set by sequential strategy)")
		Dbar<-do.call("subsampleClustering",c(list(x=x),subsampleArgs))
		if(typeAlg=="K"){
			if(is.null(DclusterArgs)) DclusterArgs<-list(k=subsampleArgs[["k"]])
			else if(!"k" %in% names(DclusterArgs)) DclusterArgs[["k"]]<-subsampleArgs[["k"]] #either sequential sets this value, or get error in subsampleClustering, so always defined.
		}

	}
	else{
		if(typeAlg!="K") stop("currently, if not subsampling, must use 'pam' or a clusterFunction defined as typeAlg='K' as clusterMethod")
		Dbar<-as.matrix(dist(x)	)		
		if(is.null(DclusterArgs) || !"k" %in% names(DclusterArgs)) stop("if not subsampling, must give k in 'DclusterArgs' (or if sequential should have been set by sequential strategy)")
	}
	if(any(is.na(as.vector(Dbar)))) stop("NA values found in Dbar (could be from too small of subsampling if classifyMethod!='All', see documentation of subsampleClustering)")
	
	res<-do.call("clusterD",c(list(D=Dbar,format="list", clusterFunction=clusterFunction),DclusterArgs)) 
	return(res) #nothing found
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



#make partition object same form as pam output
kmeansPartitionObject<-function(x,kmeansObj){
	dissE<-(cluster:::daisy(x))^2
	silObj<-cluster:::silhouette(kmeansObj$cl,dissE^2)
	silinfo<-list(widths=silObj, clus.avg.widths=summary(silObj)$clus.avg.widths, ave.width=summary(silObj)$avg.width)
	return(list(mediods=kmeansObj$centers,clustering=kmeansObj$cluster,call=NA,silinfo=silinfo,objective=NA,diss=dissE,data=x))
}

##Internal wrapper functions for kmeans and pam
.kmeansClassify <- function(x, clusterResult) { 
       centers <- clusterResult$mediods
	   suppressWarnings(stats:::kmeans(x, centers, iter.max = 1, algorithm = "Lloyd")$cluster) #probably uses this so always classifies points to centers
    } 
.kmeansCluster <- function(x,k, ...) { 
	out<-stats:::kmeans(x,centers=k,...)
	out<-kmeansPartitionObject(x,out) #make it a partition object like pam.
	#out$clustering<-out$cluster #stupid difference in naming...
	return(out)
 } 
.pamClassify <- function(x, clusterResult) {
	center<-clusterResult$medoids
	innerProd<-tcrossprod(x,center) #a n x k matrix of inner-products between them
	distMat<-as.matrix(dist(rbind(x,center)))
	distMat<-distMat[1:nrow(x),(nrow(x)+1):ncol(distMat)]
	apply(distMat,1,which.min)
} 
.pamCluster <- function(x,k, ...) { cluster:::pam(x=x,k=k,...) }
	
	
.hierCluster<-function(x,k,...){
	argList<-list(...)
	hout<-do.call("hclust",c(list(dist(x)),argList))
	stats:::cutree(tree, k = k)
	
}
.hierClassify<-function(x,clusterResult){
	
}