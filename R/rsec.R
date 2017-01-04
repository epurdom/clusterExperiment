#' Resampling-based Sequential Ensemble Clustering
#'
#' Implementation of the RSEC algorithm (Resampling-based Sequential Ensemble
#' Clustering) for single cell sequencing data. This is a wrapper function
#' around the existing clusterExperiment workflow that results in the output of
#' RSEC.
#' @param k0s the k0 parameter for sequential clustering (see \code{\link{seqCluster}})
#' @param combineProportion passed to \code{proportion} in \code{\link{combineMany}}
#' @param combineMinSize passed to \code{minSize} in \code{\link{combineMany}}
#' @param dendroReduce passed to \code{dimReduce} in \code{\link{makeDendrogram}}
#' @param dendroNDims passed to \code{ndims} in \code{\link{makeDendrogram}}
#' @param mergeMethod passed to \code{mergeMethod} in \code{\link{mergeClusters}}
#' @param mergeCutoff passed to \code{cutoff} in \code{\link{mergeClusters}}
#' @inheritParams clusterMany,matrix-method
#' @name RSEC
#' @aliases RSEC RSEC-methods RSEC,ClusterExperiment-method RSEC,matrix-method
#' @inheritParams mergeClusters,matrix-method
#' @export
setMethod(
    f = "RSEC",
    signature = signature(x = "matrix"),
    definition = function(x, isCount=FALSE,transFun=NULL,
        dimReduce="PCA",nVarDims=NA,
        nPCADims=c(50), k0s=4:15,
        clusterFunction=c("tight","hierarchical01"),
        alphas=c(0.1,0.2,0.3),betas=0.9, minSizes=1,
        combineProportion=0.7, combineMinSize=5,
        dendroReduce="mad",dendroNDims=1000,
        mergeMethod="adjP",mergeCutoff=0.05,verbose=FALSE,
        clusterDArgs=NULL,
        subsampleArgs=NULL,
        seqArgs=NULL,
        ncores=1, random.seed=NULL, run=TRUE
    )
{
    if(dimReduce=="none"){
        nPCADims<-NA
        nVarDims<-NA
    }
    if(is.null(seqArgs))seqArgs<-list(verbose=FALSE)  else seqArgs[["verbose"]]<-FALSE #turn off sequential messages
    ce<-clusterMany(x,ks=k0s,clusterFunction=clusterFunction,alphas=alphas,betas=betas,minSizes=minSizes,
                    sequential=TRUE,removeSil=FALSE,subsample=TRUE,silCutoff=0,distFunction=NA,
                    isCount=isCount,transFun=transFun,
                    dimReduce=dimReduce,nVarDims=nVarDims,nPCADims=nPCADims,
                    clusterDArgs=clusterDArgs,subsampleArgs=subsampleArgs,
                    seqArgs=seqArgs,ncores=ncores,random.seed=random.seed,run=run)
    if(run){
      ce<-.postClusterMany(ce,combineProportion,combineMinSize,dendroReduce,dendroNDims,mergeMethod,mergeCutoff,isCount)

#             ce<-combineMany(ce,whichClusters="clusterMany",proportion=combineProportion,minSize=combineMinSize)
#       if(dendroReduce=="none") dendroNDims<-NA
#       dendroTry<-try(makeDendrogram(ce,dimReduce=dendroReduce,ndims=dendroNDims,ignoreUnassignedVar=TRUE),silent=TRUE)
#       if(!inherits(dendroTry,"try-error")){
#         ce<-dendroTry
#         ce<-mergeClusters(ce,mergeMethod=mergeMethod,cutoff=mergeCutoff,plotType="none",isCount=isCount)
#       }
#       else note("makeDendrogram encountered following error and therefore clusters were not merged:\n", dendroTry)
#       return(ce) 
    }
    return(ce)
})
.postClusterMany<-function(ce,combineProportion,combineMinSize,dendroReduce,dendroNDims,mergeMethod,mergeCutoff,isCount){
  ce<-combineMany(ce,whichClusters="clusterMany",proportion=combineProportion,minSize=combineMinSize)
  if(dendroReduce=="none") dendroNDims<-NA
  dendroTry<-try(makeDendrogram(ce,dimReduce=dendroReduce,ndims=dendroNDims,ignoreUnassignedVar=TRUE),silent=TRUE)
  if(!inherits(dendroTry,"try-error")){
    ce<-dendroTry
    ce<-mergeClusters(ce,mergeMethod=mergeMethod,cutoff=mergeCutoff,plotType="none",isCount=isCount)
  }
  else note("makeDendrogram encountered following error and therefore clusters were not merged:\n", dendroTry)
  return(ce) 
}
#' @export
#' @rdname RSEC
setMethod(
  f = "RSEC",
  signature = signature(x = "SummarizedExperiment"),
  definition = function(x, ...){
    outval <- RSEC(assay(x),  ...)
    retval <- .addBackSEInfo(newObj=outval,oldObj=x)
    return(retval)

  })

#' @export
#' @rdname RSEC
setMethod(
  f = "RSEC",
  signature = signature(x = "ClusterExperiment"),
  definition = function(x, eraseOld=FALSE, rerunClusterMany=FALSE,...){
    if(rerunClusterMany){
      newObj <- RSEC(assay(x),  ...)
      ##Check if pipeline already ran previously and if so increase
      x<-.updateCurrentWorkflow(x,eraseOld,.workflowValues[-1]) #even if didn't make mergeClusters, still update it all
      if(!is.null(x)) retval<-.addNewResult(newObj=newObj,oldObj=x) #make decisions about what to keep.
      else retval<-.addBackSEInfo(newObj=newObj,oldObj=x)
    }
    else{
      retval<-.postClusterMany(x)
    }
    validObject(retval)

    return(retval)
  })
