#' @title Cluster distance matrix from subsampling
#'   
#' @description Given input data, this function will try to find the clusters
#'   based on the given ClusterFunction object.
#' @name mainClustering
#'   
#' @param x \code{p x n} data matrix on which to run the clustering (samples in 
#'   columns).
#' @param diss \code{n x n} data matrix of dissimilarities between the samples 
#'   on which to run the clustering
#' @param minSize the minimum number of samples in a cluster. Clusters found 
#'   below this size will be discarded and samples in the cluster will be given 
#'   a cluster assignment of "-1" to indicate that they were not clustered.
#' @param orderBy how to order the cluster (either by size or by maximum alpha 
#'   value). If orderBy="size" the numbering of the clusters are reordered by 
#'   the size of the cluster, instead of by the internal ordering of the 
#'   \code{clusterFUN} defined in the \code{ClusterFunction} object (an internal
#'   ordering is only possible if slot \code{outputType} of the
#'   \code{ClusterFunction} is \code{"list"}).
#' @param format whether to return a list of indices in a cluster or a vector of
#'   clustering assignments. List is mainly for compatibility with sequential 
#'   part.
#' @param clusterArgs arguments to be passed directly to the \code{clusterFUN}
#'   slot of the \code{ClusterFunction} object
#' @param checkArgs logical as to whether should give warning if arguments given
#'   that don't match clustering choices given. Otherwise, inapplicable 
#'   arguments will be ignored without warning.
#' @param returnData logical as to whether to return the \code{diss} or \code{x}
#'   matrix in the output. If \code{FALSE} only the clustering vector is
#'   returned.
#' @param ... arguments passed to the post-processing steps of the clustering.
#'   The available post-processing arguments for a \code{ClusterFunction} object
#'   depend on it's algorithm type and can be found by calling
#'   \code{getPostProcessingArgs}. See details below for documentation.
#' @inheritParams subsampleClustering
#' @inheritParams clusterSingle
#' @details \code{mainClustering} is not meant to be called by the user. It is
#'   only an exported function so as to be able to clearly document the
#'   arguments for \code{mainClustering} which can be passed via the argument
#'   \code{mainClusterArgs} in functions like \code{\link{clusterSingle}} and
#'   \code{\link{clusterMany}}.
#'   
#' @return mainClustering returns a vector of cluster assignments (if
#'   format="vector") or a list of indices for each cluster (if format="list").
#'   Clusters less than minSize are removed.
#'
#' @examples
#' data(simData)
#' cl1<-mainClustering(x=simData,clusterFunction="pam",clusterArgs=list(k=3))
#' cl2<-mainClustering(simData,clusterFunction="hierarchical01",clusterArgs=list(alpha=.1))
#' cl3<-mainClustering(simData,clusterFunction="tight",clusterArgs=list(alpha=.1))
#' #supply a dissimilarity
#' diss<-dist(t(simData),method="manhattan")
#' cl4<-mainClustering(diss=diss,clusterFunction="pam",clusterArgs=list(k=3))
#' 
#' #run hierarchical method for finding blocks, with method of evaluating
#' #coherence of block set to evalClusterMethod="average", and the hierarchical
#' #clustering using single linkage:
#' clustSubHier <- mainClustering(simData, clusterFunction="hierarchical01",
#' minSize=5, clusterArgs=list(alpha=0.1,evalClusterMethod="average", method="single"))
#'
#' #do tight
#' clustSubTight <- mainClustering(simData, clusterFunction="tight", clusterArgs=list(alpha=0.1),
#' minSize=5)
#'
#' #two twists to pam
#' clustSubPamK <- mainClustering(simData, clusterFunction="pam", silCutoff=0, minSize=5,
#' removeSil=TRUE, clusterArgs=list(k=3))
#' clustSubPamBestK <- mainClustering(simData, clusterFunction="pam", silCutoff=0,
#' minSize=5, removeSil=TRUE, findBestK=TRUE, kRange=2:10)
#'
#' # note that passing the wrong arguments for an algorithm results in warnings
#' # (which can be turned off with checkArgs=FALSE)
#' clustSubTight_test <- mainClustering(simData, clusterFunction="tight",
#' clusterArgs=list(alpha=0.1), minSize=5, removeSil=TRUE)
#' clustSubTight_test2 <- mainClustering(simData, clusterFunction="tight",
#' clusterArgs=list(alpha=0.1,evalClusterMethod="average"))
#' @rdname mainClustering
#' @aliases mainClustering,character-method
#' @export
setMethod(
    f = "mainClustering",
    signature = signature(clusterFunction = "character"),
    definition = function(clusterFunction,...){
        mainClustering(getBuiltInFunction(clusterFunction),...)
        
    }
)
#' @rdname mainClustering
#' @export
setMethod(
    f = "mainClustering",
    signature = signature(clusterFunction = "ClusterFunction"),
    definition=function(clusterFunction,inputMatrix, inputType,
                        clusterArgs=NULL,minSize=1, orderBy=c("size","best"),
                        format=c("vector","list"),
                        checkArgs=TRUE,checkDiss=FALSE,returnData=FALSE,...){
        if(missing(inputType)) stop("Internal error: inputType was not passed to mainClustering step")
        orderBy<-match.arg(orderBy)
        format<-match.arg(format)
        #######################
        ### Check arguments.
        #######################
        postProcessArgs<-list(...)
        #remove those based added by clusterSingle/seqCluster
        if("doKPostProcess" %in% names(postProcessArgs)) postProcessArgs<-postProcessArgs[-grep("doKPostProcess",names(postProcessArgs))]
        mainArgs<-c(list(clusterFunction=clusterFunction,
                         clusterArgs=clusterArgs,
                         minSize=minSize, orderBy=orderBy,
						 extraArguments=names(postProcessArgs),
                         format=format,checkArgs=checkArgs),postProcessArgs)
        checkOut<-.checkArgs(inputType=inputType, main=TRUE, subsample=FALSE, sequential=FALSE, mainClusterArgs=mainArgs, subsampleArgs=NULL,warn=TRUE)		
        if(is.character(checkOut)) stop(checkOut)
		mainClusterArgs<-checkOut$mainClusterArgs
		inputType<-mainClusterArgs[["inputType"]]
        doKPostProcess<-mainClusterArgs[["doKPostProcess"]]
        clusterFunction=mainClusterArgs[["clusterFunction"]]
        clusterArgs=mainClusterArgs[["clusterArgs"]]
        minSize=mainClusterArgs[["minSize"]] 
        orderBy=mainClusterArgs[["orderBy"]]
        format=mainClusterArgs[["format"]]
        checkArgs=mainClusterArgs[["checkArgs"]]
        postProcessArgs=mainClusterArgs[mainClusterArgs[["extraArguments"]]]
        
        
        
        #######################
        ####Run clustering:
        #######################
        N <- dim(inputMatrix)[2]
        argsClusterList<-c(clusterArgs, list("checkArgs"=checkArgs, "cluster.only"=TRUE,inputMatrix=inputMatrix,inputType=inputType))
        if(doKPostProcess) {
            if(inputType=="diss") giveDiss<-TRUE
            else giveDiss<-FALSE
            res<-do.call(".postProcessClusterK", c(list(clusterFunction=clusterFunction, clusterArgs=argsClusterList, N=N, orderBy=orderBy, giveDiss=giveDiss), postProcessArgs))
            ###Note to self: .postProcessClusterK returns clusters in list form.
        }
        else{
            res<-do.call(clusterFunction@clusterFUN,argsClusterList)
        }
        
        #######################
        #Now format into desired output, order
        #######################
        ## FIXME: need to get rid of this, unless needed/requested (though it is requested by sequential)
        #this is perhaps not efficient. For now will do this, then consider going back and only converting when, where needed.
        if(clusterFunction@outputType=="vector" & !doKPostProcess){
            res<-.clusterVectorToList(res)
        }
        clusterSize<-sapply(res, length)
        if(length(res)>0) res <- res[clusterSize>=minSize]
        if(length(res)!=0 & orderBy=="size"){ #i.e. there exist clusters found that passed minSize
            clusterSize<-sapply(res, length) #redo because dropped small clusters earlier
            res <- res[order(clusterSize,decreasing=TRUE)]
        }
        if(format=="vector"){
            res<-.clusterListToVector(res,N)
			#FIXME can't I take colnames of a diss too?
            names(res)<-if(inputType!="diss") colnames(inputMatrix) else rownames(inputMatrix)
        }
        if(!returnData) return(res)
        else return(list(results=res,inputMatrix=inputMatrix))
    }
)



#' @rdname mainClustering
#' @aliases getPostProcessingArgs
#' @details Post-processing Arguments: For post-processing the clustering,
#'   currently only type 'K' algorithms have a defined post-processing.
#'   Specifically
#' \itemize{
#'  \item{"findBestK"}{logical, whether should find best K based on average
#'   silhouette width (only used if clusterFunction of type "K").}
#'  \item{"kRange"}{vector of integers to try for k values if findBestK=TRUE. If
#'  \code{k} is given in \code{clusterArgs}, then default is k-2 to k+20,
#'  subject to those values being greater than 2; if not the default is
#'  \code{2:20}. Note that default values depend on the input k, so running for
#'  different choices of k and findBestK=TRUE can give different answers unless
#'  kRange is set to be the same.}
#'  \item{"removeSil"}{logical as to whether remove the assignment of a sample
#'  to a cluster when the sample's silhouette value is less than
#'  \code{silCutoff}}
#'  \item{"silCutoff"}{Cutoff on the minimum silhouette width to be included in
#'   cluster (only used if removeSil=TRUE).}
#' }
#' @export
setMethod(
    f = "getPostProcessingArgs",
    signature = c("ClusterFunction"),
    definition = function(clusterFunction) {
        switch(algorithmType(clusterFunction),"01"=.argsPostCluster01,"K"=.argsPostClusterK)
    }
)

.argsPostCluster01<-c("")
.argsPostClusterK<-c("findBestK","kRange","removeSil","silCutoff")

#' @importFrom cluster silhouette
.postProcessClusterK<-function(clusterFunction,findBestK=FALSE,  kRange,removeSil=FALSE,silCutoff=0,clusterArgs,N,orderBy,giveDiss=FALSE)
{
    # FIXME: this new strategy (where not have diss and x given) makes it not possible to do these post-processing steps unless inputMatrix is a dissimilarity matrix (and so clustering function works on diss)
    doPostProcess<-(findBestK | removeSil ) & giveDiss #whether will calculate silhouette or not; if not, speeds up the function... 
    k<-clusterArgs[["k"]]
    if(!findBestK && is.null(k)) stop("If findBestK=FALSE, must provide k")
    if(!is.null(k)) clusterArgs<-clusterArgs[-which(names(clusterArgs)=="k")]
    if(findBestK){
        if(missing(kRange)){
            if(!is.null(k)) kRange<-(k-2):(k+20)
            else kRange<-2:20
        }
        if(any(kRange<2)){
            kRange<-kRange[kRange>=2]
            if(length(kRange)==0) stop("Undefined values for kRange; must be greater than or equal to 2")
        }
        ks<-kRange 
    }
    else ks<-k
    if(any(ks>= N)) ks<-ks[ks<N]
    clusters<-lapply(ks,FUN=function(currk){
        cl<-do.call(clusterFunction@clusterFUN,c(list(k=currk),clusterArgs))
        if(clusterFunction@outputType=="list") cl<-.clusterListToVector(cl,N=N)
        return(cl)
    })
    if(doPostProcess){
        ###FIXME: Here is problem with post-processing -- has to use the input matrix...
        silClusters<-lapply(clusters,function(cl){
            cluster::silhouette(cl,dmatrix=clusterArgs[["inputMatrix"]])
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
    }
    else{
        cl<-clusters[[1]]
    }
    
    
    #make list of indices and put in order of silhouette width (of positive)
    clList<-tapply(seq_along(cl),cl,function(x){x},simplify=FALSE)
    if(doPostProcess){
        if(orderBy=="best"){
            clAveWidth<-tapply(sil,cl,mean,na.rm=TRUE)
            clList[order(clAveWidth,decreasing=TRUE)]
        }
        #remove -1 group
        if(removeSil){
            whNotAssign<-which(sapply(clList,function(x){all(cl[x]== -1)}))
            if(length(whNotAssign)>1) stop("Coding error in removing unclustered samples")
            if(length(whNotAssign)>0) clList<-clList[-whNotAssign]
        }
    }
    
    return(clList)
    
}


