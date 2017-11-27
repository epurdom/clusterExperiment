#' @title Resampling-based Sequential Ensemble Clustering
#'
#' @description Implementation of the RSEC algorithm (Resampling-based Sequential Ensemble
#' Clustering) for single cell sequencing data. This is a wrapper function
#' around the existing ClusterExperiment workflow that results in the output of
#' RSEC.
#' @param k0s the k0 parameter for sequential clustering (see \code{\link{seqCluster}})
#' @param combineProportion passed to \code{proportion} in \code{\link{combineMany}}
#' @param combineMinSize passed to \code{minSize} in \code{\link{combineMany}}
#' @param dendroReduce passed to \code{dimReduce} in \code{\link{makeDendrogram}}
#' @param dendroNDims passed to \code{nDims} in \code{\link{makeDendrogram}}
#' @param mergeMethod passed to \code{mergeMethod} in \code{\link{mergeClusters}}
#' @param mergeCutoff passed to \code{cutoff} in \code{\link{mergeClusters}}
#' @param rerunClusterMany logical. If the object is a ClusterExperiment object,
#'   determines whether to rerun the clusterMany step. Useful if want to try
#'   different parameters for combining clusters after the clusterMany step,
#'   without the computational costs of the clusterMany step.
#' @return A \code{ClusterExperiment} object is returned containing all of 
#' the clusterings from the steps of RSEC
#' @inheritParams clusterMany,matrix-method
#' @name RSEC
#' @aliases RSEC RSEC-methods RSEC,ClusterExperiment-method RSEC,matrix-method
#' @inheritParams mergeClusters,matrix-method

#' @export
setMethod(
    f = "RSEC",
    signature = signature(x = "matrix"),
    definition = function(x, isCount=FALSE,transFun=NULL,
        dimReduce="PCA",nFilter=NA,
        nPCADims=c(50), k0s=4:15,
        clusterFunction="hierarchical01", #listBuiltInType01(),
        alphas=c(0.1,0.2,0.3),betas=0.9, minSizes=1,
        combineProportion=0.7, combineMinSize=5,
        dendroReduce="mad",dendroNDims=1000,
        mergeMethod="adjP",mergeCutoff=0.05,verbose=FALSE,
        mainClusterArgs=NULL,
        subsampleArgs=NULL,
        seqArgs=NULL,
        ncores=1, random.seed=NULL, run=TRUE
    )
{
    if(dimReduce=="none"){
        nPCADims<-NA
        nFilter<-NA
    }
    if(is.null(seqArgs))seqArgs<-list(verbose=FALSE)  else seqArgs[["verbose"]]<-FALSE #turn off sequential messages
ce<-clusterMany(x,ks=k0s,clusterFunction=clusterFunction,alphas=alphas,betas=betas,minSizes=minSizes,
                    sequential=TRUE,removeSil=FALSE,subsample=TRUE,silCutoff=0,distFunction=NA,
                    isCount=isCount,transFun=transFun,
                    dimReduce=dimReduce,nFilter=nFilter,nPCADims=nPCADims,
                    mainClusterArgs=mainClusterArgs,subsampleArgs=subsampleArgs,
                    seqArgs=seqArgs,ncores=ncores,random.seed=random.seed,run=run)
					
    if(run){
      ce<-.postClusterMany(ce,combineProportion=combineProportion,combineMinSize=combineMinSize,dendroReduce=dendroReduce,dendroNDims=dendroNDims,mergeMethod=mergeMethod,mergeCutoff=mergeCutoff,isCount=isCount)
    }
    return(ce)
})
.methodFormals <- function(f, signature = character()) {
	#to find defaults of RSEC
	#from this conversation:
	#http://r.789695.n4.nabble.com/Q-Get-formal-arguments-of-my-implemented-S4-method-td4702420.html
    fdef <- getGeneric(f)
    method <- selectMethod(fdef, signature)
    genFormals <- base::formals(fdef)
    b <- body(method)
    if(is(b, "{") && is(b[[2]], "<-") && identical(b[[2]][[2]], as.name(".local"))) {
        local <- eval(b[[2]][[3]])
        if(is.function(local))
            return(formals(local))
        warning("Expected a .local assignment to be a function. Corrupted method?")
    }
    genFormals
}
.postClusterMany<-function(ce,...){
    defaultArgs<-.methodFormals("RSEC",signature="matrix")
	passedArgs<-list(...)
	whNotShared<-which(!names(defaultArgs)%in%names(passedArgs) )
	if(length(whNotShared)>0) passedArgs<-c(passedArgs,defaultArgs[whNotShared])
	###CombineMany
	args1<-list()
	if("combineProportion" %in% names(passedArgs)) args1<-c(args1,"proportion"=passedArgs$combineProportion)
	if("combineMinSize" %in% names(passedArgs)) args1<-c(args1,"minSize"=passedArgs$combineMinSize)
		 whClusters<-if("whichClusters" %in% names(passedArgs)) passedArgs$whichClusters else "clusterMany"
  ce<-do.call("combineMany",c(list(x=ce,whichClusters=whClusters),args1))

	##makeDendrogram
  	args1<-list()
  	if("dendroReduce" %in% names(passedArgs)){
		args1<-c(args1,"dimReduce"=passedArgs$dendroReduce)
		if(passedArgs$dendroReduce=="none") passedArgs$dendroNDims<-NA
	}
  	if("dendroNDims" %in% names(passedArgs)) args1<-c(args1,"nDims"=passedArgs$dendroNDims)
		  dendroTry<- try(do.call( "makeDendrogram", c(list(x=ce,ignoreUnassignedVar=TRUE), args1)), silent=TRUE)

		#mergeClusters
  if(!inherits(dendroTry,"try-error")){
    ce<-dendroTry
  	args1<-list()
	if("mergeCutoff" %in% names(passedArgs)) args1<-c(args1,"cutoff"=passedArgs$mergeCutoff)
	if("mergeMethod" %in% names(passedArgs) && passedArgs$mergeMethod!="none"){
		args1<-c(args1,"mergeMethod"=passedArgs$mergeMethod)
      	ce <- do.call( mergeClusters,c(list(x=ce,plot=FALSE,plotInfo="none"), args1, passedArgs[c("isCount")]))
		
	}
	else .mynote("clusters will not be merged because argument 'mergeMethod' was not given (or was equal to 'none')")
  }
  else .mynote("makeDendrogram encountered following error and therefore clusters were not merged:\n", dendroTry)
  return(ce) 
}
#' @export
#' @rdname RSEC
setMethod(
  f = "RSEC",
  signature = signature(x = "SingleCellExperiment"),
  definition = function(x, ...){
    outval <- RSEC(assay(x),  ...)
    retval <- .addBackSEInfo(newObj=outval,oldObj=x)
    return(retval)

  })

#' @export
#' @rdname RSEC
setMethod(
f = "RSEC",
signature = signature(x = "SummarizedExperiment"),
definition = function(x, ...){
	RSEC(as(x,"SingleCellExperiment"),...)

})

#' @export
#' @rdname RSEC
setMethod(
f = "RSEC",
signature = signature(x = "data.frame"),
definition = function(x,...){RSEC(data.matrix(x),...)}
)

#' @export
#' @rdname RSEC
setMethod(
  f = "RSEC",
  signature = signature(x = "ClusterExperiment"),
  definition = function(x, eraseOld=FALSE, rerunClusterMany=FALSE,...){
    if(rerunClusterMany | !"clusterMany" %in% clusterTypes(x)){
      newObj <- RSEC(assay(x),  ...)
      ##Check if pipeline already ran previously and if so increase
      x<-.updateCurrentWorkflow(x,eraseOld,.workflowValues[-1]) #even if didn't make mergeClusters, still update it all
      if(!is.null(x)) retval<-.addNewResult(newObj=newObj,oldObj=x) #make decisions about what to keep.
      else retval<-.addBackSEInfo(newObj=newObj,oldObj=x)
    }
    else{
      retval<-.postClusterMany(x,...)
    }

    return(retval)
  })
