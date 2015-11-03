

#' @title Cluster distance matrix from subsampling
#' 
#' @description Given a nxn matrix of distances, these functions will try to find the clusters based on the given clustering function. cluster01 and clusterK are internal functions and clusterD is a wrapper around these two functions for easier user interface. cluster01 and clusterK are not expected to be called directly by the user, except for ease in debugging user-defined clustering functions. 
#' 
#' @aliases cluster01
#' @aliases clusterK
#'
#' @param D the nxn matrix of 0-1 values
#' @param clusterFunction clusterFunction a function that clusters a nxn matrix of dissimilarities/distances. Can also be given character values to indicate use of internal wrapper functions for default methods. See Details for the format of what the function must take as arguments and what format the function must return.
#' @param typeAlg character value of either '01' or 'K' determining whether the function given in clusterFunction should be called by clusterK or cluster01. Only used if clusterFunction is a user-defined function. Otherwise, for methods provided by the package (i.e. by user setting clusterFunction to a character value) clusterD will determine the appropriate input for 'typeAlg' and will ignore user input.
#' @param min.size the minimum number of samples in a cluster. Clusters found below this size will be discarded and samples in the cluster will be given a cluster assignment of "-1" to indicate that they were not clustered.
#' @param orderBy how to order the cluster (either by size or by maximum alpha value)
#' @param format whether to return a list of indices in a cluster or a vector of clustering assignments. List is mainly for compatibility with sequential part. 
#' @param clusterArgs arguments to be passed directly to the clusterFunction, beyond the required input.
#' @param alpha a cutoff value of how much similarity needed for drawing blocks (lower values more strict). 
#' @param findBestK logical, whether should find best K based on average silhouette width
#' @param k single value to be used to determine how many clusters to find, if findBestK=FALSE
#' @param kRange vector of integers. If findBestK=TRUE, this gives the range of k's to look over. Default is k-2 to k+20, subject to those values being greater than 2. Note that default values depend on the input k, so running for different choices of k and findBestK=TRUE can give different answers unless kRange is set to be the same.
#' @param silCutoff Requirement on minimum silhouette width to be included in cluster (only if removeSil=TRUE)
#' @param removeSil logical as to whether remove when silhouette < silCutoff
#' @param ... arguments given to clusterD to be passed to cluster01 or clusterK (depending on the value of typeAlg). Examples include 'k' for clusterK or 'alpha' for cluster01. These should not be the arguments needed by clusterFunction (which should be passed via the argument 'clusterArgs') but the actual arguments of cluster01 or clusterK.
#' 
#' @details cluster01 is for clustering functions that expect as an input D that takes on 0-1 values (e.g. from subclustering). clusterK is for clustering functions that require an input k, the number of clusters, but arbitrary distance/dissimilarity matrix. cluster01 and clusterK are given as separate functions in order to allow the user to provide different clustering functions that expect different types of input and for us to provide different shared processing of the results that is different for these different types of clustering methods (for example, removing low silhouette values is appropriate for clusterK clustering functions rather than cluster01 functions). It is also generally expected that cluster01 algorithms use the 0-1 nature of the input to set criteria as to where to find clusters and therefore do not need a pre-determined 'k'. On the other hand, clusterK functions are assumed to need a predetermined 'k' and are also assumed to cluster all samples to a cluster, and therefore clusterK gives options to exclude poorly clustered samples via silhouette distances.
#'
#' cluster01 required format for input and output for clusterFunction: clusterFunction should be a function that takes (as a minimum) an argument "D" and "alpha". 0-1 clustering algorithms are expected to use the fact that the D input is 0-1 range to find the clusters, rather than a user defined number of clusters; "alpha" is the parameter that tunes the finding of such clusters. For example, a candidate block of samples might be considered a cluster if all values of D are greater than or equal to 1-alpha. The output is a list with each element corresponding to a cluster and the elements of the list corresponding to the indices of the samples that are in the cluster. The list is expected to be in order of 'best clusters' (as defined by the clusterFunction), with first being the best and last being worst.
#'
#' cluster01 methods: "tight" method refers to the method of finding clusters from a subsampling matrix given in the tight algorithm of Tsang and Wong. Arguments for the tight method are 'min.size.core' (default=2), which sets the minimimum number of samples that form a core cluster. "hierarchical" refers to running the hclust algorithm on D and transversing down the tree until getting a block of samples with whose summary of the values  is greater than or equal to 1-alpha. Arguments that can be passed to 'hierarchical' are 'evalClusterMethod' which determines how to summarize the samples' values of D[samples,samples] for comparison to 1-alpha: "minimum" (default) takes the minimum of D[samples,samples] and requires it to be greater than or equal to 1-alpha; "average" requires that each row mean of D[samples,samples] be greater than or equal to 1-alpha. Arguments of hclust can also be passed via clusterArgs to control the hierarchical clustering of D.
#'
#' clusterK required format for input and output for clusterFunction: clusterFunction should be a function that takes as a minimum an argument 'D' and 'k'. The output must be a list similar to that of 'partition.object' of cluster package. Specifically, an element 'clustering' which gives the vector of clusters; and an argument 'silinfo' like that of the partition.object that is a list with silhouette values. Whether these are actually silhouette values is up to the clusterFunction, but they will be used in the following way: silinfo$avg.width will be used to pick the best k (if findBestK=TRUE), silinfo$widths[,"sil_width"] will be used to exclude poorly clustered samples (if removeSil=TRUE), and clusters will be ordered by the average of the values silinfo$widths[,"sil_width"] in each cluster (after removing poorly clustered samples, if removeSil=TRUE). 

#' clusterK methods: "pam" performs pam clustering on the input Dmatrix using \code{\link{pam}} in the cluster package. Arguments to \code{\link{pam}} can be passed via 'clusterArgs', except for the arguments 'x' and 'k' which are given by D and k directly.
#'
#' @return clusterD returns a vector of cluster assignments (if format="vector") or a list of indices for each cluster (if format="list"). Clusters less than min.size are removed. If orderBy="size" the clusters are reordered by the size of the cluster, instead of by the internal ordering of the clusterFunction. 
#'
#' cluster01 and clusterK return a list of indices of the clusters found, which each element of the list corresponding to a cluster and the elements of that list a vector of indices giving the indices of the samples assigned to that cluster. Indices not included in any list are assumed to have not been clustered. The list is assumed to be ordered in terms of the `best' cluster (as defined by the clusterFunction for cluster01 or by average silhoute for clusterK), for example in terms of most internal similarity of the elements, or average silhouette width. 
#'
#' @examples
#' data(simData)
#' subD<-subsampleClustering(simData,k=3,clusterFunction="kmeans",
#' clusterArgs=list(nstart=10),resamp.n=100,samp.p=0.7)
#' 
#' #run hierarchical method for finding blocks, with method of evaluating 
#' #coherence of block set to evalClusterMethod="average", and the hierarchical 
#' #clustering using single linkage:
#' clustSubHier<-clusterD(subD,clusterFunction="hierarchical",alpha=0.1,
#' min.size=5,clusterArgs=list(evalClusterMethod="average",method="single"))
#' 
#' #note passing the wrong arguments for a '01' clusterFunction is caught 
#' #internally and ignored, but without warning:
#' clustSubTight<-clusterD(subD,clusterFunction="tight",alpha=0.1,min.size=5,
#' removeSil=TRUE)
#' 
#' #two twists to pam
#' clustSubPamK<-clusterD(subD,clusterFunction="pam",silCutoff=0,min.size=5,
#' removeSil=TRUE,k=3)
#' clustSubPamBestK<-clusterD(subD,clusterFunction="pam",silCutoff=0,min.size=5,
#' removeSil=TRUE,findBestK=TRUE,kRange=2:10)
#' 
#' #visualize the results of different clusterings
#' library(NMF)
#' clusterDF<-data.frame("hier"=factor(clustSubHier),
#' "tight"=factor(clustSubTight),"PamK"=factor(clustSubPamK),
#' "PamBestK"=factor(clustSubPamBestK))
#' maxNumb<-max(sapply(clusterDF,function(x){max(as.numeric(levels(x)))}))
#' cols<-bigPalette[1:(maxNumb+2)]
#' names(cols)<-as.character(seq(-1,maxNumb,by=1))
#' annColors<-list("hier"=cols,"tight"=cols,"PamK"=cols,"PamBestK"=cols)
#' aheatmap(subD,annCol=clusterDF,annColors=annColors,annLegend=FALSE)


clusterD<-function(D,clusterFunction=c("hierarchical","tight","pam"),typeAlg=c("01","K"),min.size=1, orderBy=c("size","best"),format=c("vector","list"),clusterArgs=NULL,...){
	passedArgs<-list(...)
	orderBy<-match.arg(orderBy)
	format<-match.arg(format)
	clusterFunction<-match.arg(clusterFunction)
	if(!is.function(clusterFunction)) typeAlg<-.checkAlgType(clusterFunction)	
	if(length(passedArgs)>0){ 
		#get rid of wrong args passed because of user confusion between the two
		whRightArgs<-which(names(passedArgs) %in% switch(typeAlg,"01"=.args01,"K"=.argsK))
		if(length(whRightArgs)>0) passedArgs<-passedArgs[whRightArgs]
		else passedArgs<-NULL
	}
	if(typeAlg=="01") {
		res<-do.call("cluster01",c(list(D=D,clusterFunction=clusterFunction,clusterArgs=clusterArgs),passedArgs))
	}
	if(typeAlg=="K") {
		res<-do.call("clusterK",c(list(D=D,clusterFunction=clusterFunction,clusterArgs=clusterArgs),passedArgs))
	}
	N<-nrow(D)

	
	#Now format into desired output
	clusterSize<-sapply(res, length)
    res <- res[clusterSize>=min.size]
	if(length(res)==0){ #No clusters pass
		if(format=="list") return(res)
		else return(rep(-1,nrow(D)))
	} 
	else{
		#now reorders final groups by size
		if(orderBy=="size"){
		  clusterSize<-sapply(res, length) #redo because dropped!
		  res <- res[order(clusterSize,decreasing=TRUE)]			
		}
		names(res)<-as.character(1:length(res))
	
		if(format=="list") return(res)	
		if(format=="vector"){
		
			valMat<-do.call("rbind",mapply(res,names(res),FUN=function(ind,val){cbind(ind,rep(as.numeric(val),length=length(ind)))},SIMPLIFY=FALSE))
			clusterVec<-rep("-1",length=N)
			clusterVec[valMat[,1]]<-valMat[,2]
			clusterVec<-as.numeric(clusterVec)
			names(clusterVec)<-rownames(D)
			return(clusterVec)
		}	
	}
}

.args01<-c("alpha")
#' @rdname clusterD
cluster01<-function(D, clusterFunction=c("hierarchical","tight"), alpha=0.1, clusterArgs=NULL)
{
	if(!is.function(clusterFunction)){
		method<-match.arg(clusterFunction)	
		##These return lists of indices of clusters satisifying alpha criteria
		if(method=="tight") clusterFunction<-.tightClusterDMat
		if(method=="hierarchical") clusterFunction<-.hierClusterDMat
	}
	res<-do.call(clusterFunction,c(list(D=D,alpha=alpha),clusterArgs))
	return(res)
}
.argsK<-c("findBestK","k","kRange","removeSil","silCutoff")
#' @rdname clusterD
clusterK<-function(D,  clusterFunction=c("pam"),findBestK=FALSE, k, kRange,removeSil=FALSE,silCutoff=0,clusterArgs=NULL)
{
	if(!findBestK && missing(k)) stop("If findBestK=FALSE, must provide k")
	if(findBestK){
		if(missing(kRange)){
			if(!missing(k)) kRange<-(k-2):(k+20)
			else kRange<-2:20
		}
		if(any(kRange<2)){ 
			kRange<-kRange[kRange>=2]
			if(length(kRange)==0) stop("Undefined values for kRange; must be greater than or equal to 2")
		}
	}
	##These return lists of indices of clusters satisifying alpha criteria
	if(!is.function(clusterFunction)){
		method<-match.arg(clusterFunction)
		if(method =="pam") clusterFunction<-function(D,k,...){cluster::pam(D,k,diss=TRUE,...)}
	} 
	if(findBestK) ks<-kRange else ks<-k
	kmeansClusters<-lapply(ks,FUN=function(currk){do.call(clusterFunction,c(list(D=D,k=currk),clusterArgs))})	
	if(length(ks)>1){
		whichBest<-which.max(sapply(kmeansClusters, function(z) z$silinfo$avg.width))
		finalCluster<-kmeansClusters[[whichBest]]		
	}
	else finalCluster<-kmeansClusters[[1]]
	sil<-finalCluster$silinfo$widths[,"sil_width"]
	if(removeSil){
		cl<-as.numeric(sil>silCutoff)
		cl[cl==0]<- -1
		cl[cl>0]<-finalCluster$clustering[cl>0]
		sil[cl == -1] <- -Inf #make the -1 cluster the last one in order
	}
	else{ 
		cl<-finalCluster$clustering
	}
	
	#make list of indices and put in order of silhouette width (of positive)
	clList<-tapply(1:length(cl),cl,function(x){x},simplify=FALSE)
	clAveWidth<-tapply(sil,cl,mean,na.rm=TRUE)
	clList[order(clAveWidth,decreasing=TRUE)]
	
	#remove -1 group
	if(removeSil){
		whNotAssign<-which(sapply(clList,function(x){all(cl[x]== -1)}))
		if(length(whNotAssign)>1) stop("Coding error in removing unclustered samples")
		if(length(whNotAssign)>0) clList<-clList[-whNotAssign]
	}
	return(clList)

}


.hierClusterDMat<-function(D,alpha,evalClusterMethod=c("minimum","average"),...)
{
	evalClusterMethod<-match.arg(evalClusterMethod)
	if(is.null(rownames(D))) rownames(D)<-colnames(D)<-as.character(1:nrow(D))
	hDmat<-hclust(dist(D),...)
	dendro<-as.dendrogram(hDmat)
	method<-evalClusterMethod
	############
	#convert into phylo4 (phylobase) object so can traverse tree easily.
	############
	#first into phylo from ape package
	tempPhylo<-try(ape::as.phylo(hDmat),FALSE)
	if(inherits(tempPhylo, "try-error")) stop("the hclust of D cannot be converted to a phylo class (ape package).")
	# require(phylobase)
	phylo4Obj<-try(as(tempPhylo,"phylo4"),FALSE) 
	if(inherits(phylo4Obj, "try-error")) stop("the created phylo object from hclust cannot be converted to a phylo4 class.")
	phylobase::nodeLabels(phylo4Obj)<-paste("Node",1:phylobase::nNodes(phylo4Obj),sep="")
	allTips<-phylobase::getNode(phylo4Obj,  type=c("tip"))
	#each internal node (including root) calculate whether passes value of alpha or not
	nodesToCheck<-phylobase::rootNode(phylo4Obj)
	clusterList<-list()
	
	#was slower with this code:
	# allInternal<-phylobase::getNode(phylo4Obj,  type=c("internal"))
	# allTipsByInternal<-lapply(allInternal,function(currNode){names(phylobase::descendants(phylo4Obj,currNode,"tip"))})
	# allChecks<-sapply(allTipsByInternal,function(currTips){
	# 	if(method=="minimum") check<-all(D[currTips,currTips,drop=FALSE]>=(1-alpha))
	# 	if(method=="average") check<-all(rowMeans(D[currTips,currTips,drop=FALSE])>=(1-alpha))
	# 		return(check)
	# })
#	names(allChecks)<-names(allInternal)
	while(length(nodesToCheck)>0){
		currNode<-nodesToCheck[1]
		nodesToCheck<-nodesToCheck[-1]
		if(currNode%in%allTips){ #block of size 1!
			currTips<-names(currNode)
			check<-TRUE
		}
		else{
			# wh<-match(currNode,allInternal)
			# currTips<-allTipsByInternal[[wh]]
			# check<-allChecks[[wh]]
			currTips<-names(phylobase::descendants(phylo4Obj,currNode,"tip"))
			if(method=="minimum") check<-all(D[currTips,currTips,drop=FALSE]>=(1-alpha))
			if(method=="average") check<-all(rowMeans(D[currTips,currTips,drop=FALSE])>=(1-alpha))
			
		}
		if(check){ #found a block that satisfies
			clusterList<-c(clusterList,list(currTips))	
		} 
		else{ #not satisfy
			childNodes<-phylobase::descendants(phylo4Obj,currNode,"children")
			nodesToCheck<-c(nodesToCheck,childNodes)
		}
	}
#	browser()
	clusterListIndices<-lapply(clusterList,function(tipNames){
		match(tipNames,rownames(D))
	})	
	clusterListIndices<-.orderByAlpha(clusterListIndices,D)
	return(clusterListIndices)
}
.orderByAlpha<-function(res,D)
{
	alphaMax<-unlist(lapply(res, function(x){
		vals<-lower.tri(D[x,x]) #don't grab diag
		1-min(vals) #max(alpha)=1-min(D)
	}))
    res <- res[order(alphaMax, decreasing=TRUE)]					
	
}
.tightClusterDMat <- function(D, alpha, min.size.core=2) 
{
    find.candidates.one <- function(x) {
        tmp <- apply(x >= 1, 1, sum) #how many in each row ==1
		#what if none of them are ==1? will this never happen because have sample of size 1? Depends what diagonal is. 
		if(all(tmp<min.size.core)){ #i.e. only core size groups less than min.size.core (default is 1)
			return(NULL)
		}
		whMax<-which(tmp == max(tmp))
		return(which(x[, whMax[1]] >= 1)) # assumes x is symmetric. Takes largest in size, but arbitrarily picks between them.
    }
    extend.candidate <- function(D, can, alpha ) { 
        can.ex <- which(apply(as.matrix(D[, can] >= 1 - alpha), 1, all)) #find those that are close to those core with 1
        D.temp <- D[can.ex, can.ex]
        if (!is.matrix(D.temp)) {
            D.temp <- as.matrix(D.temp)
            colnames(D.temp) <- names(can.ex)
        }
        D.bad <- apply(as.matrix(D.temp < 1 - alpha), 1,sum)
        while (sum(D.bad) > 0) {
            index <- which(D.bad == max(D.bad))[1]
            D.temp <- D.temp[-index, -index]
            D.bad <- apply(as.matrix(D.temp < 1 - alpha), 
              1, sum)
        }
        return(can.ex[colnames(D.temp)])
    }
	if(is.null(dim(D)) || dim(D)[1]!=dim(D)[2] || any(t(D)!=D)) stop("D must be a symmetric matrix")
	N<-nrow(D)
	colnames(D) <- 1:N
	rownames(D) <- 1:N
    i = 1
    D.temp <- D	
    res <- list()
    while (!is.null(dim(D.temp)) && !is.null(dim(D.temp)) && nrow(D.temp) > 0 & any(D.temp[lower.tri(D.temp)]>1-alpha) & any(D.temp[lower.tri(D.temp)]==1)) {
		#first find those that are always together (resampling =1); pick the largest such group (and if more than 1 of same size will arbitrarily pick one)
        candidate.one <- find.candidates.one(D.temp) 
		if(is.null(candidate.one)){#no more candidates with core always together
			#for now just stop if no core group
        	break 
		}
		#add more on the group if always resamples with the core members >alpha proportion of the time
		candidate <- extend.candidate(D.temp, candidate.one, alpha = alpha) 
        D.temp <- D.temp[-candidate, -candidate]
        res[[i]] <- names(candidate)
        mode(res[[i]]) <- "numeric"
        i = i + 1
    }
	res<-.orderByAlpha(res,D)	
	return(res)

}

.pamClusterDMat<-function(D,alpha,ks, removeSil=TRUE,method=c("pam"))
{
	method<-match.arg(method)
	kmeansClusters<-lapply(ks,FUN=function(k){cluster::pam(D,k)})		
	if(length(ks)>1){
		whichBest<-which.max(sapply(kmeansClusters, function(z) z$silinfo$avg.width))
		finalCluster<-kmeansClusters[[whichBest]]		
	}
	else finalCluster<-kmeansClusters[[1]]
	sil<-finalCluster$silinfo$widths[,"sil_width"]
	if(removeSil){
		cl<-as.numeric(sil>0)
		cl[cl==0]<- -1
		cl[cl>0]<-finalCluster$clustering[cl>0]
		sil[cl== -1] <- -Inf
	}
	else{ 
		cl<-finalCluster$clustering
	}
	#make list of indices and put in order of silhouette width (of positive)
	clList<-tapply(1:length(cl),cl,function(x){x},simplify=FALSE)
	clAveWidth<-tapply(sil,cl,mean)
	clList[order(clAveWidth,decreasing=TRUE)]
	
	#remove -1 group
	if(removeSil){
		whNotAssign<-which(sapply(clList,function(x){all(cl[x]== -1)}))
		if(length(whNotAssign)>1) stop("Coding error in removing unclustered samples")
		if(length(whNotAssign)>0) clList<-clList[-whNotAssign]
	}
	return(clList)
	
}





