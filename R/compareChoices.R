#' Create a matrix of clustering across values of k
#'
#' Given a range of k's, this funciton will return a matrix with the clustering of the samples
#' across the range, which can be passed to \code{plotTracking} for visualization.
#'
#' @param data the data on which to run the clustering (samples in rows)
#' @param ks the range of k values (see details for meaning for different choices)
#' @param alphas values of alpha to be tried. Only used for subsampleClusterMethod either 'tight' or 'hierarchical'.
#' @param findBestK values of findBestK to be tried (logical) (only for 'pam').
#' @param sequential values of sequential to be tried (logical) (only for 'pam')
#' @param removeSil values of removeSil to be tried (logical) (only for 'pam')
#' @param silCutoff values of silCutoff to be tried (only for 'pam')
#' @param clusterMethod method used in clustering of subsampled data passed to argument 'cluserFunction' of \code{\link{clusterD}}. Note that unlike other functions of this package, this must be a character vector of pre-defined clustering techniques provided by the package, and can not be user-defined.
#' @param DclusterArgs list of arguments to be passed to \code{\link{clusterD}}
#' @param subsampleArgs list of arguments to be passed to \code{\link{subsamplingClustering}}
#' @param seqArgs list of arguments to be passed to \code{\link{seqCluster}}
#' @param ncores the number of threads
#' @param random.seed a value to set seed before each run of clusterAll (so that all of the runs are run on the same subsample of the data)
#' @param ... to be passed on to mclapply (if ncores>1)
#'
#' @details While the function allows for different clusterMethod, they do not reuse the same subsampling matrix and try different clusterMethods on it. If sequential=TRUE, different subsampleClusterMethods will create different sets of data to subsample so it is not possible; if sequential=FALSE, we have not implemented functionality for this reuse. Setting the random.seed value, however, should mean that the subsampled matrix is the same, but there is no gain in computational complexity (i.e. each subsampled co-occurence matrix is recalculated for each set of parameters). 
#'
#' @return A list with the following objects:
#' \itemize{
#' \item{\code{clMat}}{a matrix of with each row corresponding to a clustering and each column a sample.}
#' \item{\code{clusterInfo}}{a list with information regarding clustering result (for those clusterings with sequential=TRUE)}
#' }
#' @examples
#' data(simData)
#' \dontrun{
#'	#following code takes around 1+ minutes to run:
#'	system.time(clusterTrack<-compareChoices(simData, ks=2:15, alphas=c(0.1,0.2,0.3), findBestK=c(TRUE,FALSE),sequential=c(FALSE),subsample=c(FALSE),removeSil=c(TRUE), clusterMethod="pam", clusterArgs = list(min.size = 5,kRange=2:15),ncores=1,random.seed=48120))
#' }


compareChoices <- function(data, ks=2:15, alphas=c(0.1,0.2,0.3), findBestK=c(TRUE,FALSE),sequential=c(TRUE,FALSE),
removeSil=c(TRUE,FALSE), subsample=c(TRUE,FALSE),silCutoff=0,
	clusterMethod=c("tight","hierarchical","pam"),
    DclusterArgs=list(min.size=5),
	subsampleArgs=list(resamp.num=50),
	seqArgs=list(beta=0.9,k.min=3),
	ncores=1,random.seed=NULL,...
	)
{
	# if(clusterMethod %in% c("pam")){ #alpha values not used -- assume cutoff silhouette at negative.
	# 	alphas<-0
	# }
	#browser()
	param<-expand.grid(k=ks,alpha=alphas,findBestK=findBestK,sequential=sequential,removeSil=removeSil,subsample=subsample,clusterMethod=clusterMethod,silCutoff=silCutoff)
	#don't vary them across ones that don't matter (i.e. 0-1 versus K); 
	#code sets to single value and then will do unique
	#also deals with just in case the user gave duplicated values of something by mistake.
	typeK<-which(param[,"clusterMethod"] %in% c("pam"))
	if(length(typeK)>0){
		param[typeK,"alpha"]<-0.01
	}
	type01<-which(param[,"clusterMethod"] %in% c("hierarchical","tight"))
	if(length(type01)>0){
		param[type01,"findBestK"]<-FALSE
		param[type01,"removeSil"]<-FALSE
		param[type01,"silCutoff"]<-0
	}
	param<-unique(param)  
	
	#deal with those that are invalid combinations:
	#Error in clusterAll(x = data, subsample = subsample, clusterFunction = clusterMethod,  : 
#  Cannot do sequential clustering where subsample=FALSE and 'findBestK=TRUE' is passed via DclusterArgs. See help documentation.

   	whInvalid<-which(!param[,"subsample"] & param[,"sequential"] & param[,"findBestK"])
	if(length(whInvalid)>0) param<-param[-whInvalid,]

	if(nrow(param)<=1) stop("set of parameters imply only 1 combination")

	#find names of the parameter combinations.
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
			subsampleArgs[["k"]]<-par[["k"]]
			DclusterArgs[["k"]]<-par[["k"]]
		}
		DclusterArgs[["alpha"]]<-par[["alpha"]]
		DclusterArgs[["findBestK"]]<-findBestK
		DclusterArgs[["removeSil"]]<-removeSil
		DclusterArgs[["silCutoff"]]<-par[["silCutoff"]]
		if(!is.null(random.seed)) set.seed(random.seed)
		clusterAll(x=data,  subsample=subsample,clusterFunction=clusterMethod,  DclusterArgs=DclusterArgs,subsampleArgs=subsampleArgs,
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