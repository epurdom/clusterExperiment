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
#' @param x \code{p x n} data matrix on which to run the clustering (samples in
#'   columns).
#' @param diss \code{n x n} data matrix of dissimilarities between the samples
#'   on which to run the clustering
##' @param clusterFunction clusterFunction a function that clusters a nxn matrix
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
#'   silhouette width (only used if clusterFunction of type "K").
#' @param k single value to be used to determine how many clusters to find, if
#'   findBestK=FALSE (only used if clusterFunction of type "K").
#' @param kRange vector of integers. If findBestK=TRUE, this gives the range of
#'   k's to look over. Default is k-2 to k+20, subject to those values being
#'   greater than 2. Note that default values depend on the input k, so running
#'   for different choices of k and findBestK=TRUE can give different answers
#'   unless kRange is set to be the same.
#' @param silCutoff Requirement on minimum silhouette width to be included in
#'   cluster (only if removeSil=TRUE).
#' @param removeSil logical as to whether remove when silhouette < silCutoff
#'   (only used if clusterFunction of type "K")
#' @param checkArgs logical as to whether should give warning if arguments given
#'   that don't match clustering choices given. Otherwise, inapplicable
#'   arguments will be ignored without warning.
#' @param returnD logical as to whether to return the D matrix in output.
#' @param ... arguments given to clusterD to be passed to cluster01 or clusterK
#'   (depending on the value of typeAlg). Examples include 'k' for clusterK or
#'   'alpha' for cluster01. These should not be the arguments needed by
#'   clusterFunction (which should be passed via the argument 'clusterArgs') but
#'   the actual arguments of cluster01 or clusterK.
#' @details To provide a distance matrix via the argument \code{distFunction},
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
#'     (1-cor(x))/2.
#'
#' @details cluster01 and
#'   clusterK are given as separate functions in order to allow the user to
#'   provide different clustering functions that expect different types of input
#'   and for us to provide different shared processing of the results that is
#'   different for these different types of clustering methods (for example,
#'   removing low silhouette values is appropriate for clusterK clustering
#'   functions rather than cluster01 functions). 
#' @return clusterD returns a vector of cluster assignments (if format="vector")
#'   or a list of indices for each cluster (if format="list"). Clusters less
#'   than minSize are removed. If orderBy="size" the clusters are reordered by
#'   the size of the cluster, instead of by the internal ordering of the
#'   clusterFunction.
#'
#' @return cluster01 and clusterK return a list of indices of the clusters found,
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
clusterD<-function(x=NULL, diss=NULL,clusterFunction=c("hierarchical01","tight","pam","hierarchicalK"),
                   typeAlg=c("01","K"),distFunction=NA,minSize=1, orderBy=c("size","best"),
                   format=c("vector","list"),clusterArgs=NULL,checkArgs=TRUE,returnD=FALSE,...){
	input<-.checkXDissInput(x,diss)
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
	#######################
	### Create distance if needed, and check it.
	#######################
	#browser()
	#browser()
	if(input=="X"){
	  if(!is.function(distFunction) && is.na(distFunction)){
	    distFunction<-switch(typeAlg,"01"=function(x){(1-cor(t(x)))/2},"K"=function(x){dist(x)})
	  }
	  D<-try(as.matrix(distFunction(t(x))))	#distances assumed to be of observations on rows
	  if(inherits(D,"try-error")) stop("input distance gives error when applied to x")
	  if(!all(dim(D) == c(ncol(x),ncol(x)))) stop("distance function must result in a ",ncol(x),"by",ncol(x),"matrix of distances")
	  if(!all(D==t(D))) stop("distance function must result in a symmetric matrix")
	  
	}
	else D<-diss
	.checkDistFunction(D)	
	#######################
	####Run clustering:
	#######################
	if(typeAlg=="01") {
	  if(any(D>1)) stop("distance function must give values between 0 and 1 for clusterFunction", clusterFunction)
		res<-do.call("cluster01",c(list(diss=D,clusterFunction=clusterFunction,clusterArgs=clusterArgs,checkArgs=checkArgs),passedArgs))
	}
	if(typeAlg=="K") {
		res<-do.call("clusterK",c(list(diss=D,clusterFunction=clusterFunction,clusterArgs=clusterArgs,checkArgs=checkArgs),passedArgs))
	}

	#######################
	#Now format into desired output
	#######################
	N<-nrow(D)
	clusterSize<-sapply(res, length)
    if(length(res)>0) res <- res[clusterSize>=minSize]
	if(length(res)==0){ #No clusters pass
# 		if(format=="list") return(res)
# 		else return(rep(-1,nrow(D)))
	    if(format=="vector") res<-rep(-1,nrow(D))
	}
	else{
		#now reorders final groups by size
		if(orderBy=="size"){
		  clusterSize<-sapply(res, length) #redo because dropped!
		  res <- res[order(clusterSize,decreasing=TRUE)]
		}
		names(res)<-as.character(1:length(res))

		#if(format=="list") return(res)
		if(format=="vector"){

			valMat<-do.call("rbind",mapply(res,names(res),FUN=function(ind,val){cbind(ind,rep(as.numeric(val),length=length(ind)))},SIMPLIFY=FALSE))
			clusterVec<-rep("-1",length=N)
			clusterVec[valMat[,1]]<-valMat[,2]
			clusterVec<-as.numeric(clusterVec)
			names(clusterVec)<-rownames(D)
		    res<-clusterVec #return(clusterVec)
		}
	}
	if(!returnD) return(res)
	else return(list(result=res,D=D))
}

.args01<-c("alpha")
#' @rdname clusterD
cluster01<-function(diss, clusterFunction=c("hierarchical01","tight"), alpha=0.1, clusterArgs=NULL,checkArgs)
{
  D<-diss
	if(!is.function(clusterFunction)){
		method<-match.arg(clusterFunction)
		##These return lists of indices of clusters satisifying alpha criteria
		if(method=="tight") clusterFunction<-.tightClusterDMat
		if(method=="hierarchical01") clusterFunction<-.hier01ClusterDMat
	}
	res<-do.call(clusterFunction,c(list(D=D,alpha=alpha,checkArgs=checkArgs),clusterArgs))
	return(res)
}
##Need to update this code so converts vector result into lists of indices ...
.hier01ClusterDMat<-function(D,alpha,evalClusterMethod=c("maximum","average"),whichHierDist=c("dist","D"),checkArgs,...)
{
    whichHierDist<-match.arg(whichHierDist)
	evalClusterMethod<-match.arg(evalClusterMethod)
	if(is.null(rownames(D))) rownames(D)<-colnames(D)<-as.character(1:nrow(D))
	passedArgs<-list(...)
	hclustArgs<-names(as.list(args(stats::hclust)))
	if(any(!names(passedArgs) %in% hclustArgs)){
	  wh<-which(!names(passedArgs) %in% hclustArgs)
		passedArgs<-passedArgs[-wh]
		if(checkArgs) warning("arguments passed via clusterArgs to hierarchical clustering method not all applicable (should only be arguments to hclust). Will be ignored")
	}
	S<-round(1-D,10)
	d<-switch(whichHierDist,"dist"=dist(S),"D"=as.dist(D))
	hDmat<-do.call(stats::hclust,c(list(d=d),passedArgs))
	
	method<-evalClusterMethod
	phylo4Obj<-.makePhylobaseTree(hDmat,"hclust")
	allTips<-phylobase::getNode(phylo4Obj,  type=c("tip"))
	#each internal node (including root) calculate whether passes value of alpha or not
	nodesToCheck<-phylobase::rootNode(phylo4Obj)
	clusterList<-list()

	while(length(nodesToCheck)>0){
		currNode<-nodesToCheck[1]
		nodesToCheck<-nodesToCheck[-1]
		if(currNode%in%allTips){ #block of size 1!
			currTips<-names(currNode)
			check<-TRUE
		}
		else{
			currTips<-names(phylobase::descendants(phylo4Obj,currNode,"tip"))
			if(method=="maximum") check<-all(S[currTips,currTips,drop=FALSE]>=(1-alpha))
			if(method=="average") check<-all(rowMeans(S[currTips,currTips,drop=FALSE])>=(1-alpha))

		}
		if(check){ #found a block that satisfies
			clusterList<-c(clusterList,list(currTips))
		}
		else{ #not satisfy
			childNodes<-phylobase::descendants(phylo4Obj,currNode,"children")
			nodesToCheck<-c(nodesToCheck,childNodes)
		}
	}
	clusterListIndices<-lapply(clusterList,function(tipNames){
		match(tipNames,rownames(D))
	})
	clusterListIndices<-.orderByAlpha(clusterListIndices,S)
	return(clusterListIndices)
}
.orderByAlpha<-function(res,S)
{
	if(length(res)>0){
		alphaMax<-unlist(lapply(res, function(x){
			vals<-lower.tri(S[x,x]) #don't grab diag
			1-min(vals) #max(alpha)=1-min(S)
		}))
	    res <- res[order(alphaMax, decreasing=TRUE)]

	}
	else return(res)
}
.tightClusterDMat <- function(D, alpha, minSize.core=2,checkArgs,...)
{
    #previously, D was similarity matrix. To make it match in clusterD, I need it to be D=1-similarity
    #so now convert it back
    S<-1-D #now is similarity matrix...
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
    extend.candidate <- function(S, can, alpha ) {
        can.ex <- which(apply(as.matrix(S[, can] >= 1 - alpha), 1, all)) #find those that are close to those core with 1
        S.temp <- S[can.ex, can.ex]
        if (!is.matrix(S.temp)) {
            S.temp <- as.matrix(S.temp)
            colnames(S.temp) <- names(can.ex)
        }
        S.bad <- apply(as.matrix(S.temp < 1 - alpha), 1,sum)
        while (sum(S.bad) > 0) {
            index <- which(S.bad == max(S.bad))[1]
            S.temp <- S.temp[-index, -index]
            S.bad <- apply(as.matrix(S.temp < 1 - alpha),
              1, sum)
        }
        return(can.ex[colnames(S.temp)])
    }
	if(is.null(dim(S)) || dim(S)[1]!=dim(S)[2] || any(t(S)!=S)) stop("S must be a symmetric matrix")
	N<-nrow(S)
	colnames(S) <- 1:N
	rownames(S) <- 1:N
    i = 1
    S.temp <- S
    res <- list()
    while (!is.null(dim(S.temp)) && !is.null(dim(S.temp)) && nrow(S.temp) > 0 & any(S.temp[lower.tri(S.temp)]>1-alpha) & any(S.temp[lower.tri(S.temp)]==1)) {
		#first find those that are always together (resampling =1); pick the largest such group (and if more than 1 of same size will arbitrarily pick one)
        candidate.one <- find.candidates.one(S.temp)
		if(is.null(candidate.one)){#no more candidates with core always together
			#for now just stop if no core group
        	break
		}
		#add more on the group if always resamples with the core members >alpha proportion of the time
		candidate <- extend.candidate(S.temp, candidate.one, alpha = alpha)
        S.temp <- S.temp[-candidate, -candidate]
        res[[i]] <- names(candidate)
        mode(res[[i]]) <- "numeric"
        i = i + 1
    }
	res<-.orderByAlpha(res,S)
	return(res)

}




.argsK<-c("findBestK","k","kRange","removeSil","silCutoff")
#' @rdname clusterD
clusterK<-function(diss,  clusterFunction=c("pam","hierarchicalK"),findBestK=FALSE, k, kRange,removeSil=FALSE,silCutoff=0,clusterArgs=NULL,checkArgs)
{
  D<-diss
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


