#' Create a matrix of clustering across values of k
#'
#' Given a range of k's, this funciton will return a matrix with the clustering of the samples
#' across the range, which can be passed to \code{plotTracking} for visualization.
#'
#' @param data the data on which to run the clustering (samples in rows)
#' @param ks the range of k values (see details for meaning for different choices)
#' @param alphas values of alpha to be tried. Only used for subsampleClusterMethod either 'tight' or 'hierarchical'; otherwise (i.e. pam) alpha=0 for removal of negative silhouette values [can't adjust this at this time.]
#' @param findBestK values of findBestK to be tried (logical).
#' @param sequential values of sequential to be tried (logical)
#' @param removeSil values of removeSil to be tried (logical)
#' @param clusterMethod if subsampleClusterMethod!="none", method used in clustering of subsampled data; if ="none", method for clustering dist(x) [currently must be "pam]
#' @param subsampleClusterMethod if "none", no subsampling is done, and clustering on dist(x) with clusterMethod technique. Otherwise, passed to clusterDMat option 'method' to indicate method of clustering of the co-occurrence matrix from subsampleClustering.
#' @param clusterArgs list of arguments to be passed to subsampleClustering's argument 'clusterArgs' that control clustering of subsampled data
#' @param subsampleArgs list of arguments to be passed to subsamplingClustering
#' @param subsampleClusterArgs list of arguments to be passed to clusterDMat
#' @param ncores the number of threads
#' @param random.seed a value to set seed before each run of clusterAll (so that all of the runs are run on the same subsample of the data)
#' @param ... to be passed on to mclapply (if ncores>1)
#'
#' @details Currently there is no functionality to run over the different subsampleClusterMethods, in the sense that there is no way to reuse the same subsampling matrix and try different subsampleClusterMethods on it. This is because if sequential=TRUE, different subsampleClusterMethods will create different sets of data to subsample (and if sequential=FALSE, we have not implemented functionality for this reuse). Therefore, only one value of subsampleClusterMethod is allowed
#'
#' @return A list with the following objects:
#' \itemize{
#' \item{\code{clMat}}{a matrix of with each row corresponding to a clustering and each column a sample.}
#' \item{\code{clusterInfo}}{a list with information regarding clustering result for those clusterings with sequential=TRUE}
#' }
#' @examples
#' data(simData)
#' system.time(clusterTrack<-compareChoices(simData, ks=2:15, alphas=c(0.1,0.2,0.3), #'	findBestK=c(TRUE,FALSE),sequential=c(FALSE),subsample=c(FALSE),removeSil=c(TRUE),
#'	clusterMethod="pam",
#'	clusterArgs = list(min.size = 5,kRange=2:15),
#'	ncores=1,random.seed=48120))

compareChoices <- function(data, ks=2:15, alphas=c(0.1,0.2,0.3), findBestK=c(TRUE,FALSE),sequential=c(TRUE,FALSE),
removeSil=c(TRUE,FALSE), subsample=c(TRUE,FALSE),
	clusterMethod=c("tight","hierarchical","pam"),
	subsampleClusterArgs=list(nstart=10),
    clusterArgs=list(min.size=5),
	subsampleArgs=list(resamp.num=50),
	seqArgs=list(beta=0.9,k.min=3),
	ncores=1,random.seed=NULL,...
	)
{
	if(clusterMethod %in% c("pam")){ #alpha values not used -- assume cutoff silhouette at negative.
		alphas<-0
	}
	param<-expand.grid(k=ks,alpha=alphas,findBestK=findBestK,sequential=sequential,removeSil=removeSil,subsample=subsample,clusterMethod=clusterMethod)
	param<-unique(param) #just in case the user gave duplicated values of something by mistake. 
   	whInvalid<-which(!param[,"subsample"] & param[,"sequential"] & param[,"findBestK"])
	if(length(whInvalid)>0) param<-param[-whInvalid,]
	whVary<-which(apply(param,2,function(x){length(unique(x))>1}))
	if(length(whVary)>0) cnames<-apply(param[,whVary,drop=FALSE],1,function(x){
		paste(colnames(param)[whVary],x,sep="=",collapse=",")})
	else stop("set of parameters imply only 1 combination")
	cat(nrow(param),"parameter combinations,",sum(param[,"sequential"]),"use sequential method.\n")
	paramFun<-function(i){
		par<-param[i,]
		#make them logical values... otherwise adds a space before the TRUE and doesn't recognize.
		#well, sometimes. Maybe don't need this?
		removeSil<-as.logical(gsub(" ","",par["removeSil"]))
		sequential<-as.logical(gsub(" ","",par["sequential"]))
		subsample<-as.logical(gsub(" ","",par["subsample"]))
		findBestK<-as.logical(gsub(" ","",par["findBestK"]))
		clusterMethod<-as.character(par[["clusterMethod"]])
		if(sequential) seqArgs[["k0"]]<-par[["k"]] 
		else{
			#to be safe, set both in case user set one. 
			subsampleClusterArgs[["k"]]<-par[["k"]]
			clusterArgs[["k"]]<-par[["k"]]
		}
		clusterArgs[["alpha"]]<-par[["alpha"]]
		clusterArgs[["findBestK"]]<-findBestK
		clusterArgs[["removeSil"]]<-removeSil
		if(!is.null(random.seed)) set.seed(random.seed)
		clusterAll(x=data,  subsample=subsample,clusterMethod=clusterMethod,  clusterArgs=clusterArgs,subsampleArgs=subsampleArgs,subsampleClusterArgs=subsampleClusterArgs,
			seqArgs=seqArgs, sequential=sequential) 
	}
	if(ncores>1) out<-mclapply(1:nrow(param),FUN=paramFun,mc.cores=ncores,...)
	else out<-lapply(1:nrow(param),FUN=paramFun)

	clMat<-sapply(out,function(x){x$clustering})
	colnames(clMat)<-cnames
	
	#   #just return matrices of clusters
	#   out<-lapply(allTracking,function(trackTightAlpha){
	#     sapply(trackTightAlpha,function(x){x$cluster})
	#   })
	clusterInfo<-lapply(out,function(x){
		# if(all(c("clusterInfo","whyStop") %in% names(x))) return(x[c("clusterInfo","whyStop")])
# 		else return(NULL)
		return(x[c("clusterInfo","whyStop")])
	})
	 names(clusterInfo)<-cnames
	
	return(list(clMat=clMat,clusterInfo=clusterInfo))
}

# .trackPam<-function(data, ks){
#   d <- dist(data)
#   pamRes <- lapply(ks, function(z) pam(d,k=z))
#   clMat<-sapply(pamRes,function(x){x$clustering})
#   colnames(clMat)<-paste("K=",ks,sep="")
#   return(list(clMat=clMat))
# }
#
# .trackTight<-function(data, ks=2:15, pickLargestBeta=TRUE, nstart=10, resamp.num=50,
#                       k.stop=3, alphas=c(0.1,0.2,0.3), target=100, beta=0.9, standardize.gene=FALSE,
#                       random.seed=NULL, ncores=1, verbose=FALSE, ...){
#   tightFun<-function(k,a){
#     if(verbose) cat("--------- k0=",k,"----------")
#     tight.clust.fixed(data,target=target,k.min=k,alpha=a,beta=beta,standardize.gene=standardize.gene,random.seed =random.seed,pickLargestBeta=pickLargestBeta,nstart=nstart,resamp.num=resamp.num,k.stop=k.stop,verbose=verbose)
#   }
#   alphaFun<-function(a){
#     if(verbose) cat("######## alpha=",a,"########")
#     trackTightAlpha<-lapply(ks,tightFun,a=a)
#     names(trackTightAlpha)<-paste("K=",ks,sep="")
#     return(trackTightAlpha)
#   }
#   if(ncores>1){
#     allTracking<-mclapply(alphas,FUN=alphaFun,mc.cores=ncores,...)
#   }
#   else{
#     allTracking<-lapply(alphas,alphaFun)
#   }
#   names(allTracking)<-paste("alpha=",alphas,sep="")
#
#   #just return matrices of clusters
#   out<-lapply(allTracking,function(trackTightAlpha){
#     sapply(trackTightAlpha,function(x){x$cluster})
#   })
#   clusterInfo<-lapply(allTracking,function(trackTightAlpha){
#     lapply(trackTightAlpha,function(x){x$clusterInfo})
#   })
#   names(out)<-names(allTracking)
#   return(list(clMat=out,clusterInfo=clusterInfo))
# }