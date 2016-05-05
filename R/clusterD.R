#' @title Cluster distance matrix from subsampling
#'
#' @description Given a \code{n x n} matrix of distances, these functions will
#'   try to find the clusters based on the given clustering function. cluster01
#'   and clusterK are internal functions and clusterD is a wrapper around these
#'   two functions for easier user interface. cluster01 and clusterK are not
#'   expected to be called directly by the user, except for ease in debugging
#'   user-defined clustering functions.
#'
#' @aliases cluster01
#' @aliases clusterK
#'
#' @param D either a \code{n x n} matrix of 0-1 values or a \code{p x n} matrix of data.
#' @param clusterFunction clusterFunction a function that clusters a nxn matrix
#'   of dissimilarities/distances. Can also be given character values to
#'   indicate use of internal wrapper functions for default methods. See Details
#'   for the format of what the function must take as arguments and what format
#'   the function must return.
#' @param typeAlg character value of either '01' or 'K' determining whether the
#'   function given in clusterFunction should be called by clusterK or
#'   cluster01. Only used if clusterFunction is a user-defined function.
#'   Otherwise, for methods provided by the package (i.e. by user setting
#'   clusterFunction to a character value) clusterD will determine the
#'   appropriate input for 'typeAlg' and will ignore user input.
#' @param distFunction a distance function to be applied to \code{D}. Only relevant if
#'   input \code{D} is a matrix of data, rather than a distance. See details.
#' @param minSize the minimum number of samples in a cluster. Clusters found
#'   below this size will be discarded and samples in the cluster will be given
#'   a cluster assignment of "-1" to indicate that they were not clustered.
#' @param orderBy how to order the cluster (either by size or by maximum alpha
#'   value).
#' @param format whether to return a list of indices in a cluster or a vector of
#'   clustering assignments. List is mainly for compatibility with sequential
#'   part.
#' @param clusterArgs arguments to be passed directly to the clusterFunction,
#'   beyond the required input.
#' @param alpha a cutoff value of how much similarity needed for drawing blocks
#'   (lower values more strict).
#' @param findBestK logical, whether should find best K based on average
#'   silhouette width.
#' @param k single value to be used to determine how many clusters to find, if
#'   findBestK=FALSE.
#' @param kRange vector of integers. If findBestK=TRUE, this gives the range of
#'   k's to look over. Default is k-2 to k+20, subject to those values being
#'   greater than 2. Note that default values depend on the input k, so running
#'   for different choices of k and findBestK=TRUE can give different answers
#'   unless kRange is set to be the same.
#' @param silCutoff Requirement on minimum silhouette width to be included in
#'   cluster (only if removeSil=TRUE).
#' @param removeSil logical as to whether remove when silhouette < silCutoff
#' @param checkArgs logical as to whether should give warning if arguments given
#'   that don't match clustering choices given. Otherwise, inapplicable
#'   arguments will be ignored without warning.
#' @param ... arguments given to clusterD to be passed to cluster01 or clusterK
#'   (depending on the value of typeAlg). Examples include 'k' for clusterK or
#'   'alpha' for cluster01. These should not be the arguments needed by
#'   clusterFunction (which should be passed via the argument 'clusterArgs') but
#'   the actual arguments of cluster01 or clusterK.
#'
#' @details cluster01 is for clustering functions that expect as an input D that
#'   takes on 0-1 values (e.g. from subclustering). clusterK is for clustering
#'   functions that require an input k, the number of clusters, but arbitrary
#'   distance/dissimilarity matrix. cluster01 and clusterK are given as separate
#'   functions in order to allow the user to provide different clustering
#'   functions that expect different types of input and for us to provide
#'   different shared processing of the results that is different for these
#'   different types of clustering methods (for example, removing low silhouette
#'   values is appropriate for clusterK clustering functions rather than
#'   cluster01 functions). It is also generally expected that cluster01
#'   algorithms use the 0-1 nature of the input to set criteria as to where to
#'   find clusters and therefore do not need a pre-determined 'k'. On the other
#'   hand, clusterK functions are assumed to need a predetermined 'k' and are
#'   also assumed to cluster all samples to a cluster, and therefore clusterK
#'   gives options to exclude poorly clustered samples via silhouette distances.
#'
#'   @details cluster01 required format for input and output for clusterFunction:
#'   clusterFunction should be a function that takes (as a minimum) an argument
#'   "D" and "alpha". 0-1 clustering algorithms are expected to use the fact
#'   that the D input is 0-1 range to find the clusters, rather than a user
#'   defined number of clusters; "alpha" is the parameter that tunes the finding
#'   of such clusters. For example, a candidate block of samples might be
#'   considered a cluster if all values of D are greater than or equal to
#'   1-alpha. The output is a list with each element corresponding to a cluster
#'   and the elements of the list corresponding to the indices of the samples
#'   that are in the cluster. The list is expected to be in order of 'best
#'   clusters' (as defined by the clusterFunction), with first being the best
#'   and last being worst.
#'
#'   @details cluster01 methods: "tight" method refers to the method of finding 
#'     clusters from a subsampling matrix given internally in the tight 
#'     algorithm code of Tsang and Wong. Arguments for the tight method are
#'     'minSize.core' (default=2), which sets the minimimum number of samples
#'     that form a core cluster. "hierarchical01" refers to running the hclust
#'     algorithm on D and transversing down the tree until getting a block of
#'     samples with whose summary of the values  is greater than or equal to
#'     1-alpha. Arguments that can be passed to 'hierarchical' are
#'     'evalClusterMethod' which determines how to summarize the samples' values
#'     of D[samples,samples] for comparison to 1-alpha: "minimum" (default)
#'     takes the minimum of D[samples,samples] and requires it to be greater
#'     than or equal to 1-alpha; "average" requires that each row mean of
#'     D[samples,samples] be greater than or equal to 1-alpha. Arguments of
#'     hclust can also be passed via clusterArgs to control the hierarchical 
#'     clustering of D.
#'
#'   @details clusterK required format for input and output for clusterFunction:
#'   clusterFunction should be a function that takes as a minimum an argument
#'   'D' and 'k'. The output must be a clustering, specified by integer values. 
#'   The function \code{\link{silhouette}} will be used on the clustering to
#'   calculate silhouette scores for each observation.
#'
#' @details clusterK methods: "pam" performs pam clustering on the input 
#'   \code{D} matrix using \code{\link{pam}} in the cluster package. Arguments 
#'   to \code{\link{pam}} can be passed via 'clusterArgs', except for the 
#'   arguments 'x' and 'k' which are given by D and k directly. "hierarchicalK"
#'   performs hierarchical clustering on the input via the \code{\link{hclust}}
#'   and then applies \code{\link{cutree}} with the specified k to obtain
#'   clusters. Arguments to \code{\link{hclust}} can be passed via
#'   \code{clusterArgs}.
#'   @details To provide a distance matrix via the argument \code{distFunction},
#'     the function must be defined to take the distance of the rows of a matrix
#'     (internally, the function will call \code{distFunction(t(x))}. This is to
#'     be compatible with the input for the \code{dist} function.
#'     \code{as.matrix} will be performed on the output of \code{distFunction},
#'     so if the object returned has a \code{as.matrix} method that will convert
#'     the output into a symmetric matrix of distances, this is fine (for
#'     example the class \code{dist} for objects returned by \code{dist} have
#'     such a method). If \code{distFunction=NA}, then a default distance will 
#'     be calculated based on the type of clustering algorithm of 
#'     \code{clusterFunction}. For type "K" the default is to take \code{dist}
#'     as the distance function. For type "01", the default is to take the
#'     absolute value of the correlation between the samples.
#'
#' @return clusterD returns a vector of cluster assignments (if format="vector")
#'   or a list of indices for each cluster (if format="list"). Clusters less
#'   than minSize are removed. If orderBy="size" the clusters are reordered by
#'   the size of the cluster, instead of by the internal ordering of the
#'   clusterFunction.
#'
#'   @return cluster01 and clusterK return a list of indices of the clusters found,
#'   which each element of the list corresponding to a cluster and the elements
#'   of that list a vector of indices giving the indices of the samples assigned
#'   to that cluster. Indices not included in any list are assumed to have not
#'   been clustered. The list is assumed to be ordered in terms of the `best'
#'   cluster (as defined by the clusterFunction for cluster01 or by average
#'   silhoute for clusterK), for example in terms of most internal similarity of
#'   the elements, or average silhouette width.
#'
#' @examples
#' data(simData)
#' cl1<-clusterD(simData,clusterFunction="pam",k=3)
#' cl2<-clusterD(simData,clusterFunction="hierarchical01")
#' cl3<-clusterD(simData,clusterFunction="tight")
#' #change distance to manhattan distance
#' cl4<-clusterD(simData,clusterFunction="pam",k=3,
#'      distFunction=function(x){dist(x,method="manhattan")})
#' 
#' #run hierarchical method for finding blocks, with method of evaluating
#' #coherence of block set to evalClusterMethod="average", and the hierarchical
#' #clustering using single linkage:
#' clustSubHier <- clusterD(simData, clusterFunction="hierarchical01", alpha=0.1,
#' minSize=5, clusterArgs=list(evalClusterMethod="average", method="single"))
#'
#' #do tight
#' clustSubTight <- clusterD(simData, clusterFunction="tight", alpha=0.1,
#' minSize=5)
#'
#' #two twists to pam
#' clustSubPamK <- clusterD(simData, clusterFunction="pam", silCutoff=0, minSize=5,
#' removeSil=TRUE, k=3)
#' clustSubPamBestK <- clusterD(simData, clusterFunction="pam", silCutoff=0,
#' minSize=5, removeSil=TRUE, findBestK=TRUE, kRange=2:10)
#'
#' # note that passing the wrong arguments for an algorithm results in warnings
#' # (which can be turned off with checkArgs=FALSE)
#' clustSubTight_test <- clusterD(simData, clusterFunction="tight", alpha=0.1,
#' minSize=5, removeSil=TRUE)
#' clustSubTight_test2 <- clusterD(simData, clusterFunction="tight", alpha=0.1,
#' clusterArgs=list(evalClusterMethod="average"))
#' @export
#' @importFrom cluster daisy silhouette pam
clusterD<-function(D,clusterFunction=c("hierarchical01","tight","pam","hierarchicalK"),
                   typeAlg=c("01","K"),distFunction=NA,minSize=1, orderBy=c("size","best"),
                   format=c("vector","list"),clusterArgs=NULL,checkArgs=TRUE,...){
	passedArgs<-list(...)
	orderBy<-match.arg(orderBy)
	format<-match.arg(format)
	clusterFunction<-match.arg(clusterFunction)
	if(!is.function(clusterFunction)) typeAlg<-.checkAlgType(clusterFunction)
	if(length(passedArgs)>0){
		#get rid of wrong args passed because of user confusion between the two
		whRightArgs<-which(names(passedArgs) %in% switch(typeAlg,"01"=.args01,"K"=.argsK))
		if(length(whRightArgs)!=length(passedArgs) & checkArgs) warning("Some arguments passed via '...' do not match the choice of typeAlg")
		if(length(whRightArgs)>0) passedArgs<-passedArgs[whRightArgs]
		else passedArgs<-NULL
	}
	### Create distance if needed, and check it.
	#browser()
	if(!all(dim(D)==dim(t(D))) || !all(na.omit(D==t(D)))){
	  x<-D 
	  if(!is.function(distFunction) && is.na(distFunction)){
	    distFunction<-switch(typeAlg,"01"=function(x){abs(cor(t(x)))},"K"=function(x){dist(x)})
	  }
	  D<-try(as.matrix(distFunction(t(x))))	#distances assumed to be of observations on rows
	  if(inherits(D,"try-error")) stop("input distance gives error when applied to x")
	  if(!all(dim(D) == c(ncol(x),ncol(x)))) stop("distance function must result in a ",ncol(x),"by",ncol(x),"matrix of distances")
	  if(!all(D==t(D))) stop("distance function must result in a symmetric matrix")
	  
	}
	if(any(is.na(as.vector(D)))) stop("NA values found in Dbar (could be from too small of subsampling if classifyMethod!='All', see documentation of subsampleClustering)")
	if(any(is.na(D) | is.nan(D) | is.infinite(D))) stop("D matrix contains either NAs, NANs or Infinite values.")
	if(any(D<0)) stop("distance function must give strictly positive values")
	
	####Run clustering:
	if(typeAlg=="01") {
	  if(any(D>1)) stop("distance function must give values between 0 and 1 for clusterFunction", clusterFunction)
		res<-do.call("cluster01",c(list(D=D,clusterFunction=clusterFunction,clusterArgs=clusterArgs,checkArgs=checkArgs),passedArgs))
	}
	if(typeAlg=="K") {
		res<-do.call("clusterK",c(list(D=D,clusterFunction=clusterFunction,clusterArgs=clusterArgs,checkArgs=checkArgs),passedArgs))
	}
	N<-nrow(D)


	#Now format into desired output
	clusterSize<-sapply(res, length)
    if(length(res)>0) res <- res[clusterSize>=minSize]
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
cluster01<-function(D, clusterFunction=c("hierarchical01","tight"), alpha=0.1, clusterArgs=NULL,checkArgs)
{
	if(!is.function(clusterFunction)){
		method<-match.arg(clusterFunction)
		##These return lists of indices of clusters satisifying alpha criteria
		if(method=="tight") clusterFunction<-.tightClusterDMat
		if(method=="hierarchical01") clusterFunction<-.hier01ClusterDMat
	}
	res<-do.call(clusterFunction,c(list(D=D,alpha=alpha,checkArgs=checkArgs),clusterArgs))
	return(res)
}
.hier01ClusterDMat<-function(D,alpha,evalClusterMethod=c("minimum","average"),checkArgs,...)
{
	evalClusterMethod<-match.arg(evalClusterMethod)
	if(is.null(rownames(D))) rownames(D)<-colnames(D)<-as.character(1:nrow(D))
	passedArgs<-list(...)
	hclustArgs<-names(as.list(args(stats::hclust)))
	if(any(!names(passedArgs) %in% hclustArgs)){
	  wh<-which(!names(passedArgs) %in% hclustArgs)
		passedArgs<-passedArgs[-wh]
		if(checkArgs) warning("arguments passed via clusterArgs to hierarchical clustering method not all applicable (should only be arguments to hclust). Will be ignored")
	}
	hDmat<-do.call(stats::hclust,c(list(d=dist(D)),passedArgs))
	method<-evalClusterMethod
	phylo4Obj<-.makePhylobaseTree(hDmat,"hclust")
	# ############
	# #convert into phylo4 (phylobase) object so can traverse tree easily.
	# ############
	# #first into phylo from ape package
	# tempPhylo<-try(ape::as.phylo(hDmat),FALSE)
	# if(inherits(tempPhylo, "try-error")) stop("the hclust of D cannot be converted to a phylo class (ape package).")
	# # require(phylobase)
	# phylo4Obj<-try(as(tempPhylo,"phylo4"),FALSE)
	# if(inherits(phylo4Obj, "try-error")) stop("the created phylo object from hclust cannot be converted to a phylo4 class.")
	# phylobase::nodeLabels(phylo4Obj)<-paste("Node",1:phylobase::nNodes(phylo4Obj),sep="")



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
	if(length(res)>0){
		alphaMax<-unlist(lapply(res, function(x){
			vals<-lower.tri(D[x,x]) #don't grab diag
			1-min(vals) #max(alpha)=1-min(D)
		}))
	    res <- res[order(alphaMax, decreasing=TRUE)]

	}
	else return(res)
}
.tightClusterDMat <- function(D, alpha, minSize.core=2,checkArgs,...)
{
    if(length(list(...))>0 & checkArgs) 	warning("some arguments passed via clusterArgs to tight clustering method are not applicable")
	find.candidates.one <- function(x) {
        tmp <- apply(x >= 1, 1, sum) #how many in each row ==1
		#what if none of them are ==1? will this never happen because have sample of size 1? Depends what diagonal is.
		if(all(tmp<minSize.core)){ #i.e. only core size groups less than minSize.core (default is 1)
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




.argsK<-c("findBestK","k","kRange","removeSil","silCutoff")
#' @rdname clusterD
clusterK<-function(D,  clusterFunction=c("pam","hierarchicalK"),findBestK=FALSE, k, kRange,removeSil=FALSE,silCutoff=0,clusterArgs=NULL,checkArgs)
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
    if(method =="pam") clusterFunction<-function(D,k,checkArgs,...){
      passedArgs<-list(...)
      pamArgs<-names(as.list(args(cluster::pam)))
      if(any(wh<-!names(passedArgs) %in% pamArgs)){
        passedArgs<-passedArgs[-which(wh)]
        if(checkArgs) warning("arguments passed via clusterArgs to pam not all applicable (should only be arguments to pam). Will be ignored")
      }
      do.call(cluster::pam,c(list(x=D,k=k,diss=TRUE,cluster.only=TRUE),passedArgs))
      
    }
    if(method =="hierarchicalK") clusterFunction<-function(D,k,checkArgs,...){
      passedArgs<-list(...)
      hierArgs<-names(as.list(args(stats::hclust)))
      if(any(wh<-!names(passedArgs) %in% hierArgs)){
        passedArgs<-passedArgs[-which(wh)]
        if(checkArgs) warning("arguments passed via clusterArgs to pam not all applicable (should only be arguments to pam). Will be ignored")
      }
#      browser()
      hclustOut<-do.call(stats::hclust,c(list(d=as.dist(D)),passedArgs))
      cutree(hclustOut,k)
    }
  }


  if(findBestK) ks<-kRange else ks<-k
  if(any(ks>= nrow(D))) ks<-ks[ks<nrow(D)]
  #browser()
  clusters<-lapply(ks,FUN=function(currk){do.call(clusterFunction,c(list(D=D,k=currk,checkArgs=checkArgs),clusterArgs))})
  silClusters<-lapply(clusters,function(cl){
    silhouette(cl,dmatrix=D)
  })
  if(length(ks)>1){
    whichBest<-which.max(sapply(silClusters, mean))
    finalCluster<-clusters[[whichBest]]
    sil<-silClusters[[whichBest]][,"sil_width"]
  }
  else{
    finalCluster<-clusters[[1]]
    sil<-silClusters[[1]][,"sil_width"]
  }
  if(removeSil){
    cl<-as.numeric(sil>silCutoff)
    cl[cl==0]<- -1
    cl[cl>0]<-finalCluster[cl>0]
    sil[cl == -1] <- -Inf #make the -1 cluster the last one in order
  }
  else{
    cl<-finalCluster
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


