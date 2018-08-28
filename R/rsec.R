#' @title Resampling-based Sequential Ensemble Clustering
#'
#' @description Implementation of the RSEC algorithm (Resampling-based
#'   Sequential Ensemble Clustering) for single cell sequencing data. This is a
#'   wrapper function around the existing ClusterExperiment workflow that
#'   results in the output of RSEC.
#' @param k0s the k0 parameter for sequential clustering (see
#'   \code{\link{seqCluster}})
#' @param consensusProportion passed to \code{proportion} in
#'   \code{\link{makeConsensus}}
#' @param consensusMinSize passed to \code{minSize} in
#'   \code{\link{makeConsensus}}
#' @param dendroReduce passed to \code{reduceMethod} in
#'   \code{\link{makeDendrogram}}
#' @param dendroNDims passed to \code{nDims} in \code{\link{makeDendrogram}}
#' @param mergeMethod passed to \code{mergeMethod} in
#'   \code{\link{mergeClusters}}
#' @param mergeCutoff passed to \code{cutoff} in \code{\link{mergeClusters}}
#' @param mergeLogFCcutoff passed to \code{logFCcutoff} in
#'   \code{\link{mergeClusters}}
#' @param mergeDEMethod passed to \code{DEMethod} argument in
#'   \code{\link{mergeClusters}}. By default, unless otherwise chosen by the
#'   user, if \code{isCount=TRUE}, then \code{mergeDEMethod="limma-voom"}, 
#' otherwise \code{mergeDEMethod="limma"}. These choices are for speed
#'  considerations and the user may want to try \code{mergeDEMethod="edgeR"} 
#' on smaller datasets of counts.
#' @param rerunClusterMany logical. If the object is a ClusterExperiment object,
#'   determines whether to rerun the clusterMany step. Useful if want to try
#'   different parameters for combining clusters after the clusterMany step,
#'   without the computational costs of the clusterMany step.
#' @param stopOnErrors logical. If \code{FALSE}, if RSEC hits an error
#'   \emph{after} the \code{clusterMany} step, it will return the results up to
#'   that point, rather than generating a stop error. The text of error will be
#'   printed as a NOTE. This allows the user to get the results to that point,
#'   so as to not have to rerun the computationally heavy earlier steps. The
#'   \code{TRUE} option is only provided for debugging purposes.
#' @return A \code{\link{ClusterExperiment}} object is returned containing all
#'   of the clusterings from the steps of RSEC
#' @inheritParams clusterMany
#' @name RSEC
#' @aliases RSEC RSEC-methods RSEC,ClusterExperiment-method RSEC,matrix-method RSEC,SingleCellExperiment-method RSEC,SummarizedExperiment-method
#' @inheritParams mergeClusters
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
		passedArgs<-list(...)
		if("isCount" %in% names(passedArgs) & !"mergeDEMethod" %in% names(passedArgs)){
			if(passedArgs$isCount) passedArgs$mergeDEMethod<-"limma-voom"
			else passedArgs$mergeDEMethod<-"limma"
		}
    if(rerunClusterMany | !"clusterMany" %in% clusterTypes(x)){
	  	if(any(c("transFun","isCount") %in% names(list(...))))
	  		warning("The internally saved transformation function of a ClusterExperiment object must be used when given as input and setting 'transFun' or 'isCount' for a 'ClusterExperiment' has no effect other than to set a default for 'mergeDEMethod' (if not set by user).")
      newObj <- do.call("RSEC",c(list(x=as(x,"SingleCellExperiment"),  transFun=transformation(x)),passedArgs))
      ##Check if pipeline already ran previously and if so increase
			x<-.updateCurrentWorkflow(x,eraseOld,newTypeToAdd=.workflowValues[-1],newLabelToAdd=NULL)
			if(!is.null(x)) retval<-.addNewResult(newObj=newObj,oldObj=x) #make decisions about what to keep.
      else{
		  retval<-.addBackSEInfo(newObj=newObj,oldObj=x)

		  }
		  filterStats(retval)<-filterStats(newObj)
		  reducedDims(retval)<-reducedDims(newObj)
    }
    else{
      retval<-do.call(".postClusterMany",c(list(ce=x),passedArgs))
    }

    return(retval)
  })

#' @export
#' @rdname RSEC
setMethod(
  f = "RSEC",
  signature = signature(x = "matrixOrHDF5"),
  definition = function(x, ...){
    return(RSEC(SingleCellExperiment(x),...))

  })

#' @rdname RSEC
#' @param whichAssay numeric or character specifying which assay to use. See
#'   \code{\link[SummarizedExperiment]{assay}} for details.
#' @export
setMethod(
  f = "RSEC",
  signature = signature(x = "SingleCellExperiment"),
  definition = function(x,
		isCount=FALSE,
		transFun=NULL,
    reduceMethod="PCA",
    nFilterDims=defaultNDims(x,reduceMethod,type="filterStats"),
    nReducedDims=defaultNDims(x,reduceMethod,type="reducedDims"),
		k0s=4:15,
		subsample=TRUE,
		sequential=TRUE,
    clusterFunction="hierarchical01", #listBuiltInType01(),
    alphas=c(0.1,0.2,0.3),betas=0.9, minSizes=1,
    consensusProportion=0.7,
    consensusMinSize,
    dendroReduce,
    dendroNDims,
    mergeMethod="adjP",
    mergeCutoff,
    mergeLogFCcutoff,
		mergeDEMethod=if(isCount) "limma-voom" else "limma",
    verbose=FALSE,
    mainClusterArgs=NULL,
    subsampleArgs=NULL,
    seqArgs=NULL, whichAssay = 1,
    ncores=1, random.seed=NULL,
		stopOnErrors=FALSE, run=TRUE
  )
  {
    if(reduceMethod=="none"){
      nReducedDims<-NA
      nFilterDims<-NA
    }
    if(is.null(seqArgs))
      seqArgs<-list(verbose=FALSE)
    else seqArgs[["verbose"]]<-FALSE #turn off sequential messages
    ce<-clusterMany(x,ks=k0s,clusterFunction=clusterFunction,
                    alphas=alphas,betas=betas,minSizes=minSizes,
                    sequential=sequential,removeSil=FALSE,subsample=subsample,
                    silCutoff=0,distFunction=NA,
                    isCount=isCount,transFun=transFun,
										reduceMethod=reduceMethod,nFilterDims=eval(nFilterDims),
                    nReducedDims=eval(nReducedDims),
                    mainClusterArgs=mainClusterArgs,subsampleArgs=subsampleArgs,
                    seqArgs=seqArgs,ncores=ncores,random.seed=random.seed,run=run,
                    whichAssay=whichAssay)

    if(run){
      #first add ones that have default value
      passedArgs<-list(ce=ce,consensusProportion=consensusProportion,
                       mergeMethod=mergeMethod,whichAssay=whichAssay,
											 stopOnErrors=stopOnErrors)
      #add those who will use default value from the function -- is there easier way
      if(!missing(consensusProportion))
        passedArgs<-c(passedArgs,consensusProportion=consensusProportion)
      if(!missing(consensusMinSize))
        passedArgs<-c(passedArgs,consensusMinSize=consensusMinSize)
      if(!missing(dendroReduce))
        passedArgs<-c(passedArgs,dendroReduce=dendroReduce)
      if(!missing(dendroNDims))
        passedArgs<-c(passedArgs,dendroNDims=dendroNDims)
      if(!missing(mergeCutoff))
        passedArgs<-c(passedArgs,mergeCutoff=mergeCutoff)
      if(!missing(mergeLogFCcutoff))
        passedArgs<-c(passedArgs,mergeLogFCcutoff=mergeLogFCcutoff)
			mergeDEMethod<-eval(mergeDEMethod)
			if(!missing(mergeDEMethod))
				passedArgs<-c(passedArgs,mergeDEMethod=eval(mergeDEMethod))
      ce<-do.call(".postClusterMany",passedArgs)
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
.postClusterMany<-function(ce,stopOnErrors=FALSE,...){
  defaultArgs<-.methodFormals("RSEC",signature="SingleCellExperiment")
  #remove those without anything defined
  defaultArgs<-defaultArgs[which(sapply(.methodFormals("RSEC",signature="SingleCellExperiment"),function(x){!isTRUE(x=="")}))]
  passedArgs<-list(...)
  whNotShared<-which(!names(defaultArgs)%in%names(passedArgs) )
  if(length(whNotShared)>0) passedArgs<-c(passedArgs,defaultArgs[whNotShared])
  #------------
  ###CombineMany
  #------------
  args1<-list()
  if("consensusProportion" %in% names(passedArgs)) args1<-c(args1,"proportion"=passedArgs$consensusProportion)
  if("consensusMinSize" %in% names(passedArgs)) args1<-c(args1,"minSize"=passedArgs$consensusMinSize)
  whClusters<-if("whichClusters" %in% names(passedArgs)) passedArgs$whichClusters  	else "clusterMany"
  combineTry<-try(do.call("makeConsensus",c(list(x=ce,whichClusters=whClusters),args1)), silent=TRUE)
  if(!inherits(combineTry,"try-error")){
    ce<-combineTry
    #------------
    ##makeDendrogram
    #------------
    args1<-list()
    if("dendroReduce" %in% names(passedArgs)){
      args1<-c(args1,"reduceMethod"=passedArgs$dendroReduce, "whichAssay"=passedArgs$whichAssay)
      if(passedArgs$dendroReduce=="none") passedArgs$dendroNDims<-NA
    }
    if("dendroNDims" %in% names(passedArgs)) args1<-c(args1,"nDims"=passedArgs$dendroNDims)
    dendroTry<- try(do.call( "makeDendrogram", c(list(x=ce,filterIgnoresUnassigned=TRUE), args1)), silent=TRUE)

    #------------
    #mergeClusters
    #------------
    if(!inherits(dendroTry,"try-error")){
      ce<-dendroTry

      if("mergeMethod" %in% names(passedArgs) && passedArgs$mergeMethod!="none"){
        args1<-list()
        args1 <- c(args1, "whichAssay"=passedArgs$whichAssay)
        args1<-c(args1,"mergeMethod"=passedArgs$mergeMethod)
        if("mergeCutoff" %in% names(passedArgs)) args1<-c(args1,"cutoff"=passedArgs$mergeCutoff)
        if("mergeLogFCCutoff" %in% names(passedArgs)){
          args1<-c(args1,"logFCcutoff="=passedArgs$mergeLogFCCutoff)
        }
				if("mergeDEMethod" %in% names(passedArgs)){
					args1<-c(args1,"DEMethod"=passedArgs$mergeDEMethod)
					mergeTry <- try(do.call( mergeClusters,c(list(x=ce,plot=FALSE,plotInfo="none"), args1 )), silent=TRUE)

				}
				else{
					mergeTry<-"mergeDEMethod argument is missing with no default"
					class(mergeTry)<-"try-error"
				}
        if(!inherits(mergeTry,"try-error")){
          ce<-mergeTry
        }
        else{
        	if(!stopOnErrors).mynote(paste("mergeClusters encountered following error and therefore clusters were not merged:\n", mergeTry))
					else stop(mergeTry)
        }
      }
      else .mynote("clusters will not be merged because argument 'mergeMethod' was not given (or was equal to 'none')")
    }
    else{
    	if(!stopOnErrors).mynote(paste("makeDendrogram encountered following error and therefore clusters were not merged:\n", dendroTry))
			else stop(dendroTry)

    }
  }
  else{
  	if(!stopOnErrors).mynote(paste("makeConsensus encountered following error and therefore clusters from clusterMany were not combined:\n", combineTry))
		else stop(combineTry)

  }
  return(ce)
}
