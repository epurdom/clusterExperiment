#' @include internalFunctions.R internalClusterFunctions.R

### input to clustering:
# pam : x or dis
# hier : dis
# kmeans : x
# spectral (SamSPECTRAL for flow cytometry; kernlab for standard; kknn for similarity based on knn rather than kmeans): kernlab is either x or a kernel function
###what rdname should this be? Not S4 methods...

################
##Internal wrapper functions for kmeans and pam
################
.genericClassify<-function(x,centers){
    innerProd<-tcrossprod(t(x),centers) #a n x k matrix of inner-products between them
    distMat<-as.matrix(dist(rbind(t(x),centers)))
    distMat<-distMat[1:ncol(x),(ncol(x)+1):ncol(distMat)]
    apply(distMat,1,which.min)	
}
.getPassedArgs<-function(FUN,passedArgs,checkArgs){
    funArgs<-names(as.list(args(FUN)))
	funName<-tail(as.character(substitute(FUN)),1)
    if(any(wh<-!names(passedArgs) %in% funArgs)){
      passedArgs<-passedArgs[-which(wh)]
      if(checkArgs) warning(.wrongArgsWarning(funName))
    }
	return(passedArgs)
}
.wrongArgsWarning<-function(funName){paste("arguments passed via clusterArgs to the clustering function",funName,"are not all applicable (clusterArgs should only be arguments to,", funName,"). Extra arguments will be ignored")}
##---------
##Kmeans
##---------
#' @importFrom stats kmeans
.kmeansCluster <- function(x,k, checkArgs,cluster.only,...) { 
  passedArgs<-.getPassedArgs(FUN=stats::kmeans,passedArgs=list(...) ,checkArgs=checkArgs)
  out<-do.call(stats::kmeans,c(list(x=t(x),centers=k),passedArgs))
  if(cluster.only) return(out$cluster)
  else return(.kmeansPartitionObject(x,out)) 
} 
.kmeansClassify <- function(x, clusterResult) { 
  centers <- clusterResult$mediods
  suppressWarnings(stats::kmeans(t(x), centers, iter.max = 1, algorithm = "Lloyd")$cluster) #probably uses this so always classifies points to centers
} 
#make partition object same form as pam output
#' @importFrom cluster daisy silhouette
.kmeansPartitionObject<-function(x,kmeansObj){ 
  dissE<-(cluster::daisy(t(x)))^2
  silObj<-try(cluster::silhouette(x=kmeansObj$cluster,dist=dissE))
  silinfo<-list(widths=silObj, clus.avg.widths=summary(silObj)$clus.avg.widths, ave.width=summary(silObj)$avg.width)
  return(list(mediods=kmeansObj$centers, clustering=kmeansObj$cluster, call=NA,silinfo=silinfo, objective=NA, diss=dissE, data=x))
}
.kmeansCF<-clusterFunction(clusterFUN=.kmeansCluster, classifyFUN=.kmeansClassify, inputType="X", inputClassifyType="X", algorithmType="K",outputType="vector")
#internalFunctionCheck(.kmeansCluster,inputType="X",algType="K",outputType="vector")

##---------
##PAM
##---------

#' @importFrom cluster pam
.pamCluster<-function(x,diss,k,checkArgs,cluster.only,...){
      passedArgs<-.getPassedArgs(FUN=cluster::pam,passedArgs=list(...) ,checkArgs=checkArgs)
	  input<-.checkXDissInput(x,diss,checkDiss=FALSE,algType="K")
	  if(input=="X") return(do.call(cluster::pam, c(list(x=t(x),k=k, cluster.only=cluster.only), passedArgs)))
      if(input=="diss" | input=="both") return(do.call(cluster::pam, c(list(x=diss,k=k, diss=TRUE, cluster.only=cluster.only), passedArgs)))
    }
.pamClassify <- function(x, clusterResult) { #x p x n matrix
  .genericClassify(x,clusterResult$medoids)
} 
.pamCF<-clusterFunction(clusterFUN=.pamCluster, classifyFUN=.pamClassify, inputType="either", inputClassifyType="X", algorithmType="K",outputType="vector")

#internalFunctionCheck(.pamCluster,inputType="either",algType="K",outputType="vector")

##---------
##Hiearchical01
##---------

#' @importFrom stats hclust 
#' @importFrom phylobase rootNode getNode descendants
.hier01Cluster<-function(diss,alpha,evalClusterMethod=c("maximum","average"),whichHierDist=c("as.dist","dist"),checkArgs,cluster.only,...)
{
    whichHierDist<-match.arg(whichHierDist)
	evalClusterMethod<-match.arg(evalClusterMethod)
	if(is.null(rownames(diss))) rownames(diss)<-colnames(diss)<-as.character(1:nrow(diss))
	passedArgs<-.getPassedArgs(FUN=stats::hclust,passedArgs=list(...) ,checkArgs=checkArgs)
	S<-round(1-diss,10)
	d<-switch(whichHierDist,"dist"=dist(S),"as.dist"=as.dist(diss))
	hDmat<-do.call(stats::hclust,c(list(d=d),passedArgs))
	
	##Could this be just done by cut with hierarchical cophenic value? Should make it an option to do that. Probably a lot faster...
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
		match(tipNames,rownames(diss))
	})
	clusterListIndices<-.orderByAlpha(clusterListIndices,S)
	##Need to update this code so converts vector result into lists of indices ...
	return(clusterListIndices)
}
.hier01CF<-clusterFunction(clusterFUN=.hier01Cluster, inputType="diss", algorithmType="01",outputType="list")

##---------
##hiearchicalK
##---------
#' @importFrom stats hclust cutree
.hierKCluster<-function(diss,k,checkArgs,cluster.only,...){
	passedArgs<-.getPassedArgs(FUN=stats::hclust,passedArgs=list(...) ,checkArgs=checkArgs)
    hclustOut<-do.call(stats::hclust,c(list(d=as.dist(diss)),passedArgs))
    stats::cutree(hclustOut,k)
}
.hierKCF<-clusterFunction(clusterFUN=.hierKCluster, inputType="diss", algorithmType="K",outputType="vector")

#internalFunctionCheck(.hierKCluster,inputType="diss",algType="K",outputType="vector")


##---------
##Tight
##---------
.tightCluster <- function(diss, alpha, minSize.core=2,checkArgs,cluster.only,...)
{
    #previously, diss was similarity matrix. To make it match all of the code, I need it to be diss=1-similarity so now convert it back
    S<-1-diss #now is similarity matrix...
    if(length(list(...))>0 & checkArgs) 	warning(.wrongArgsWarning("tight"))
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
	##Need to update this code so converts vector result into lists of indices ...
	return(res)

}
.tightCF<-clusterFunction(clusterFUN=.tightCluster, inputType="diss", algorithmType="01",outputType="list")


#########
## Put them together so user/code can access easily
#########
.builtInClusterObjects<-list("pam"=.pamCF,"kmeans"=.kmeansCF,"hierarchical01"=.hier01CF,"hierarchicalK"=.hierKCF,"tight"=.tightCF)

#' @title Built in ClusterFunction options 
#' @description Documents the built-in clustering options that are available in the clusterExperiment package. 
#' @rdname builtInClusteringFunctions
#' @details Built-in clustering methods: The built-in clustering methods, the names of which can be accessed by \code{builtInClusterFunctions} are the following:
#' \itemize{
#'  \item{"pam"}{Based on \code{\link{pam}} in \code{cluster} package. Arguments to that function can be passed via \code{clusterArgs}. }
#'  \item{"kmeans"}{Based on \code{\link{kmeans}} in \code{stats} package. Arguments to that function can be passed via \code{clusterArgs} except for \code{centers} which is reencoded here to be the argument 'k'}
#'  \item{"hierarchical01"}{\code{\link{hclust}} in \code{stats} package is used to build hiearchical clustering. Arguments to that function can be passed via \code{clusterArgs}. The \code{hierarchical01} cuts the hiearchical tree based on the parameter \code{alpha}. It does not use the \code{cutree} function, but instead ... [documentation still needed]   }
#'  \item{"hierarchicalK"}{\code{\link{hclust}} in \code{stats} package is used to build hiearchical clustering and \code{\link{cutree}} is used to cut the tree into \code{k} clusters.}
#'  \item{"tight"}{Based on the algorithm in Tsang and Wong, specifically their method of picking clusters from a co-occurance matrix after subsampling. The clustering encoded here is not the entire tight clustering algorithm, only that single piece that identifies clusters from the co-occurance matrix.  }
#' }
#' @examples
#' builtInClusterFunctions
#' getBuiltInClusterFunction("kmeans")
#' @export
builtInClusterFunctions<-names(.builtInClusterObjects)

#' @rdname builtInClusteringFunctions
#' @aliases getBuiltInClusterFunction
#' @export
setMethod(
  f = "getBuiltInClusterFunction",
  signature = c("character"),
  definition = function(object) {
  	if(!all(object%in%builtInClusterFunctions)) stop("if give character value for 'object' must be one of",paste(builtInClusterFunctions,collapse=","))
  	    m<-match(object,names(.builtInClusterObjects))
  		if(length(m)>1) .builtInClusterObjects[m]
			else .builtInClusterObjects[[m]]
		
    
	    }
)
#' @rdname builtInClusteringFunctions
#' @aliases getBuiltInAlgorithmType
#' @export
setMethod(
  f = "getBuiltInAlgorithmType",
  signature = c("character"),
  definition = function(object) {
	  clObjects<-getBuiltInClusterFunction(object)
	  if(length(clObjects)>1) return(sapply(clObjects,algorithmType))
		  else return(algorithmType(clObjects))
  }
)
#' @rdname builtInClusteringFunctions
#' @export
setMethod(
  f = "getBuiltInAlgorithmType",
  signature = c("factor"),
  definition = function(object) {
	  getBuiltInAlgorithmType(as.character(object))
  }
)
#' @rdname builtInClusteringFunctions
#' @aliases getBuiltInTypeK
#' @export
setMethod(
  f = "getBuiltInTypeK",
  signature = c("ANY"),
  definition = function(object) {
	  allBuiltInTypes<-getBuiltInAlgorithmType(builtInClusterFunctions)
	  return(names(allBuiltInTypes)[allBuiltInTypes=="K"])
  }
)
#' @rdname builtInClusteringFunctions
#' @aliases getBuiltInType01
#' @export
setMethod(
  f = "getBuiltInType01",
  signature = c("ANY"),
  definition = function(object) {
	  allBuiltInTypes<-getBuiltInAlgorithmType(builtInClusterFunctions)
	  return(names(allBuiltInTypes)[allBuiltInTypes=="01"])
  }
	  
)