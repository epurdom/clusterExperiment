##Note to self:
## 5/9/2019: currently, if classifyMethod=="InSample", cluster.only=TRUE, otherwise FALSE. Our default is "All", so generally doing the cluster.only=FALSE. 


#' @include internalFunctions.R internalClusterFunctions.R internalDendroFunctions.R

################
##Internal wrapper functions for kmeans and pam
################
.genericClassify<-function(x,centers){
  if(inherits(x,"DelayedArray") || inherits(centers,"DelayedArray")){
	pracma::distmat(DelayedArray::DelayedArray(t(x)),DelayedArray::DelayedArray(centers))
  }
  else{
    distMat<-pracma::distmat(t(x),centers)
  }
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
##Spectral
##---------

# spectral options (SamSPECTRAL for flow cytometry (function SamSPECTRAL); kernlab for standard ('specc'); kknn for similarity based on knn rather than kmeans): kernlab is either x or a kernel function
#' @importFrom kernlab specc
.speccCluster<-function(x,k,checkArgs,cluster.only,...){
  passedArgs<-.getPassedArgs(FUN=kernlab::specc,passedArgs=list(...) ,checkArgs=checkArgs)
  #add data.matrix here for hdf5Matrix. Not optimized
  out<-try(do.call(kernlab::specc,c(list(x=data.matrix(t(x)), centers=k),passedArgs)))
  if(inherits(out,"try-error"))stop("Spectral clustering failed, probably because k (",k,") was too large relative to the number of samples (",ncol(x),"). k must be less than the number of samples, but how much less is not straightforward.")
  if(cluster.only) return(out@.Data)
  else return(out) 
}
.speccCF<-ClusterFunction(clusterFUN=.speccCluster, classifyFUN=NULL, inputType="X", algorithmType="K",outputType="vector")


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
	#This is hugely computationally expensive. 
  # dissE<-(cluster::daisy(t(x)))^2
  # silObj<-try(cluster::silhouette(x=kmeansObj$cluster,dist=dissE),silent=TRUE)
  # if(!inherits(silObj,"try-error"))
  #   silinfo<-list(widths=silObj, clus.avg.widths=summary(silObj)$clus.avg.widths, ave.width=summary(silObj)$avg.width)
  # else silinfo<-NULL
  return(list(mediods=kmeansObj$centers, clustering=kmeansObj$cluster, call=NA,silinfo=silinfo, objective=NA, diss=dissE, data=x))
}
.kmeansCF<-ClusterFunction(clusterFUN=.kmeansCluster, classifyFUN=.kmeansClassify, inputType="X", inputClassifyType="X", algorithmType="K",outputType="vector")
#internalFunctionCheck(.kmeansCluster,inputType="X",algType="K",outputType="vector")

##---------
##PAM
##---------

#' @importFrom cluster pam
.pamCluster<-function(x,diss,k,checkArgs,cluster.only,...){
  
  passedArgs<-.getPassedArgs(FUN=cluster::pam,passedArgs=list(...) ,checkArgs=checkArgs)
  input<-.checkXDissInput(x,diss,checkDiss=FALSE,algType="K")
  if(input=="X"){
    return(do.call(cluster::pam, c(list(x=t(x),k=k, cluster.only=cluster.only), passedArgs)))  
  } 
  if(input=="diss" | input=="both") return(do.call(cluster::pam, c(list(x=diss,k=k, diss=TRUE, cluster.only=cluster.only), passedArgs)))
}
.pamClassify <- function(x, clusterResult) { #x p x n matrix
  .genericClassify(x,clusterResult$medoids)
} 
.pamCF<-ClusterFunction(clusterFUN=.pamCluster, classifyFUN=.pamClassify, inputType="either", inputClassifyType="X", algorithmType="K",outputType="vector")

#internalFunctionCheck(.pamCluster,inputType="either",algType="K",outputType="vector")

##---------
##clara
##---------

#' @importFrom cluster pam
.claraCluster<-function(x,k,checkArgs,cluster.only,samples=50,keep.data=FALSE,rngR=TRUE,pamLike=TRUE,correct.d=TRUE,medoids.x=FALSE,...){
  passedArgs<-.getPassedArgs(FUN=cluster::clara,passedArgs=list(...) ,checkArgs=checkArgs)
  passedArgs<-c(passedArgs, list(samples=samples, keep.data=keep.data, rngR=rngR, pamLike=pamLike, correct.d=correct.d))
  out<-(do.call(cluster::clara, c(list(x=t(x),k=k), passedArgs)))
  if(cluster.only) return(out$clustering) else return(out)
  
}

.claraCF<-ClusterFunction(clusterFUN=.claraCluster, classifyFUN=.pamClassify, inputType="X", inputClassifyType="X", algorithmType="K",outputType="vector")


##---------
##hiearchicalK
##---------
#' @importFrom stats hclust cutree
.hierKCluster<-function(diss,k,checkArgs,cluster.only,...){
  passedArgs<-.getPassedArgs(FUN=stats::hclust,passedArgs=list(...) ,checkArgs=checkArgs)
  hclustOut<-do.call(stats::hclust,c(list(d=as.dist(diss)),passedArgs))
  stats::cutree(hclustOut,k)
}
.hierKCF<-ClusterFunction(clusterFUN=.hierKCluster, inputType="diss", algorithmType="K",outputType="vector")

#internalFunctionCheck(.hierKCluster,inputType="diss",algType="K",outputType="vector")


##---------
##Hiearchical01
##---------

#' @importFrom stats hclust 
.hier01Cluster<-function(diss,alpha,evalClusterMethod=c("maximum","average"),whichHierDist=c("as.dist","dist"),checkArgs,cluster.only,...)
{
  whichHierDist<-match.arg(whichHierDist)
  evalClusterMethod<-match.arg(evalClusterMethod)
  if(is.null(rownames(diss))) rownames(diss)<-colnames(diss)<-as.character(seq_len(nrow(diss)))
  passedArgs<-.getPassedArgs(FUN=stats::hclust,passedArgs=list(...) ,checkArgs=checkArgs)
  S<-round(1-diss,10)
  d<-switch(whichHierDist,"dist"=dist(S),"as.dist"=as.dist(diss))
  hDmat<-do.call(stats::hclust,c(list(d=d),passedArgs))
  
  ##Could this be just done by cut with hierarchical cophenic value? Should make it an option to do that. Probably a lot faster...
  method<-evalClusterMethod
  phylo4Obj<-.convertToPhyClasses(hDmat,returnClass=c("phylo4"))
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
.hier01CF<-ClusterFunction(clusterFUN=.hier01Cluster, inputType="diss", algorithmType="01",outputType="list")


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
  colnames(S) <- seq_len(N)
  rownames(S) <- seq_len(N)
  i <- 1
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
.tightCF<-ClusterFunction(clusterFUN=.tightCluster, inputType="diss", algorithmType="01",outputType="list")


#########
## Put them together so user/code can access easily
#########
.builtInClusterObjects<-list("pam"=.pamCF,"clara"=.claraCF,"kmeans"=.kmeansCF,"hierarchical01"=.hier01CF,"hierarchicalK"=.hierKCF,"tight"=.tightCF,"spectral"=.speccCF)
.builtInClusterNames<-names(.builtInClusterObjects)

#' @title Built in ClusterFunction options
#' @param object name of built in function.
#' @description Documents the built-in clustering options that are available in
#'   the clusterExperiment package.
#' @rdname builtInClusteringFunctions
#' @details \code{listBuiltInFunctions} will return the character names of the
#'   built-in clustering functions available.
#' @details \code{listBuiltInTypeK} returns the names of the built-in functions
#'   that have type 'K'
#' @details \code{listBuiltInType01} returns the names of the built-in functions
#'   that have type '01'
#' @details \code{getBuiltInFunction} will return the \code{ClusterFunction}
#'   object of a character value that corresponds to a built-in function.
#' @details \code{\link{algorithmType}} and \code{\link{inputType}} will return
#'   the \code{algorithmType} and \code{inputType} of the built-in
#'   clusterFunction corresponding to the character value.
#' @details \strong{Built-in clustering methods:} The built-in clustering
#'   methods, the names of which can be accessed by
#'   \code{listBuiltInFunctions()} are the following:
#' \itemize{ 
#' \item{"pam"}{Based on \code{\link[cluster]{pam}} in
#'   \code{cluster} package. Arguments to that function can be passed via
#'   \code{clusterArgs}. 
#' Input is \code{"either"} (\code{x} or \code{diss}); algorithm type is "K"} 
#' \item{"clara"}{Based on \code{\link[cluster]{clara}} in
#'   \code{cluster} package. Arguments to that function can be passed via
#'   \code{clusterArgs}. Note that we have changed the default arguments of 
#'   that function to match the recommendations in the documentation of
#'  \code{\link[cluster]{clara}} (numerous functions are set to less than optimal 
#'  settings for back-compatiability). Specifically, the following defaults 
#'  are implemented \code{samples=50}, \code{keep.data=FALSE}, 
#' \code{mediods.x=FALSE},\code{rngR=TRUE},
#'  \code{pamLike=TRUE}, \code{correct.d=TRUE}. 
#' Input is \code{"X"}; algorithm type is "K".} 
#' \item{"kmeans"}{Based on \code{\link[stats]{kmeans}} in \code{stats} package.
#' Arguments to that function can be passed via \code{clusterArgs} except for
#' \code{centers} which is reencoded here to be the argument 'k' Input is
#' \code{"X"}; algorithm type is "K"}
#' \item{"hierarchical01"}{\code{\link[stats]{hclust}} in \code{stats} package
#' is used to build hiearchical clustering. Arguments to that function can be
#' passed via \code{clusterArgs}. The \code{hierarchical01} cuts the hiearchical
#' tree based on the parameter \code{alpha}. It does not use the \code{cutree}
#' function, but instead transversing down the tree until getting a block of
#' samples with whose summary of the values  is greater than or equal to
#' 1-alpha. Arguments that can be passed to 'hierarchical01' are
#' 'evalClusterMethod' which determines how to summarize the samples' values of
#' D[samples,samples] for comparison to 1-alpha: "maximum" (default) takes the
#' minimum of D[samples,samples] and requires it to be less than or equal to
#' 1-alpha; "average" requires that each row mean of D[samples,samples] be less
#' than or equal to 1-alpha. Additional arguments of hclust can also be passed
#' via clusterArgs to control the hierarchical clustering of D. Input is
#' \code{"diss"}; algorithm type is "01"}
#' \item{"hierarchicalK"}{\code{\link[stats]{hclust}} in \code{stats} package is used
#'   to build hiearchical clustering and \code{\link{cutree}} is used to cut the
#'   tree into \code{k} clusters.
#' Input is \code{"diss"}; algorithm type is "K"}   
#' \item{"tight"}{Based on the algorithm in
#'   Tsang and Wong, specifically their method of picking clusters from a
#'   co-occurance matrix after subsampling. The clustering encoded here is not
#'   the entire tight clustering algorithm, only that single piece that
#'   identifies clusters from the co-occurance matrix.  
#' Arguments for the tight method are 
#'   'minSize.core' (default=2), which sets the minimimum number of samples that
#'   form a core cluster.
#' Input is \code{"diss"}; algorithm type is "01"} 
#' \item{"spectral"}{\code{\link[kernlab]{specc}} in \code{kernlab} package 
#' is used to perform spectral clustering. Note that spectral clustering can 
#' produce errors if the number of clusters (K) is not sufficiently smaller than 
#' the number of samples (N). K < N is not always sufficient. 
#' Input is \code{"X"}; algorithm type is "K".}
#' }
#' @seealso \code{\link{ClusterFunction}}, \code{\link{algorithmType}},
#'   \code{\link{inputType}}
#' @examples
#' listBuiltInFunctions()
#' algorithmType(c("kmeans","pam","hierarchical01"))
#' inputType(c("kmeans","pam","hierarchical01"))
#' listBuiltInTypeK()
#' listBuiltInType01()
#' @rdname builtInClusteringFunctions
#' @aliases listBuiltInFunctions
#' @return \code{listBuiltInFunctions} returns a character vector of all 
#' the built-in cluster functions' names.
#' @export
listBuiltInFunctions<-function() {
  .builtInClusterNames
  
}
#' @rdname builtInClusteringFunctions
#' @aliases getBuiltInFunction
#' @return \code{getBuiltInFunction} returns the \code{ClusterFunction} 
#' object that corresponds to the character name of a function
#' @export
setMethod(
  f = "getBuiltInFunction",
  signature = c("character"),
  definition = function(object) {
    if(!all(object%in%.builtInClusterNames)) stop("if give character value for a clusterFunction object must be one of",paste(.builtInClusterNames,collapse=","))
    m<-match(object,names(.builtInClusterObjects))
    if(length(m)>1) .builtInClusterObjects[m]
    else .builtInClusterObjects[[m]]
    
    
  }
)

#' @rdname builtInClusteringFunctions
#' @aliases listBuiltInTypeK
#' @return \code{listBuiltInTypeK} returns a character vector of the 
#' names of built-in cluster functions that are of type "K"
#' @export
listBuiltInTypeK<-function() {
  allBuiltInTypes<-algorithmType(.builtInClusterNames)
  return(names(allBuiltInTypes)[allBuiltInTypes=="K"])
}

#' @rdname builtInClusteringFunctions
#' @aliases listBuiltInType01
#' @return \code{listBuiltInType01} returns a character vector of the 
#' names of built-in cluster functions that are of type "01"
#' @export
listBuiltInType01<-function() {
  allBuiltInTypes<-algorithmType(.builtInClusterNames)
  return(names(allBuiltInTypes)[allBuiltInTypes=="01"])
}

