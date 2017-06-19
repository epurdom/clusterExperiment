#' General wrapper method to cluster the data
#'
#' Given a data matrix, \code{\link{SummarizedExperiment}}, or
#' \code{\link{ClusterExperiment}} object, this function will find clusters,
#' based on a single specification of parameters.
#'
#' @param x the data on which to run the clustering (features in rows).
#' @param diss \code{n x n} data matrix of dissimilarities between the samples
#'   on which to run the clustering (only if \code{subsample=FALSE})
#' @param subsample logical as to whether to subsample via 
#'   \code{\link{subsampleClustering}} to get the distance matrix at each 
#'   iteration; otherwise the distance function will be determined by argument
#'   \code{distFunction} passed in \code{clusterDArgs} (if input a data matrix).
#' @param sequential logical whether to use the sequential strategy (see
#'   details of \code{\link{seqCluster}}).
#' @param clusterDArgs list of additional arguments to be passed to
#'   \code{\link{clusterD}}.
#' @param subsampleArgs list of arguments to be passed to
#'   \code{\link{subsampleClustering}}.
#' @param seqArgs list of additional arguments to be passed to
#'   \code{\link{seqCluster}}.
#' @param isCount logical. Whether the data are in counts, in which case the
#'   default \code{transFun} argument is set as log2(x+1). This is simply a
#'   convenience to the user, and can be overridden by giving an explicit
#'   function to \code{transFun}.
#' @param transFun function A function to use to transform the input data matrix
#'   before clustering.
#' @param dimReduce character A character identifying what type of 
#'   dimensionality reduction to perform before clustering. Options are 
#'   "none","PCA", "var","cv", and "mad". See \code{\link{transform}} for more
#'   details.
#' @param ndims integer An integer identifying how many dimensions to reduce to
#'   in the reduction specified by \code{dimReduce}
#' @param clusterLabel a string used to describe the clustering. By
#'   default it is equal to "clusterSingle", to indicate that this clustering is
#'   the result of a call to \code{clusterSingle}.

#' @param ... arguments to be passed on to the method for signature
#'   \code{matrix}.
#'
#' @details If sequential=TRUE, the sequential clustering controls the 'k'
#'   argument of the underlying clustering so setting 'k=' in the list given to
#'   clusterDArgs or subsampleArgs will not do anything and will produce a
#'   warning to that effect.
#'
#' @return A \code{\link{ClusterExperiment}} object.
#'
#' @seealso \code{\link{clusterMany}} to compare multiple choices of parameters.
#'
#' @name clusterSingle
#'
#' @examples
#' data(simData)
#'
#' \dontrun{
#' #following code takes some time.
#' #use clusterSingle to do sequential clustering
#' #(same as example in seqCluster only using clusterSingle ...)
# ' set.seed(44261)
# ' clustSeqHier_v2 <- clusterSingle(simData, 
# ' sequential=TRUE, subsample=TRUE, subsampleArgs=list(resamp.n=100, samp.p=0.7,
# ' clusterFunction="kmeans", clusterArgs=list(nstart=10)),
# ' seqArgs=list(beta=0.8, k0=5), clusterDArgs=list(minSize=5,clusterFunction="hierarchical01",clusterArgs=list(alpha=0.1)))
#' }
#'
#' #use clusterSingle to do just clustering k=3 with no subsampling
#' clustNothing <- clusterSingle(simData, 
#' subsample=FALSE, sequential=FALSE, clusterDArgs=list(clusterFunction="pam",clusterArgs=list(k=3)))
#' #compare to standard pam
#' cluster::pam(t(simData),k=3,cluster.only=TRUE)
#' @aliases clusterSingle clusterSingle-methods clusterSingle,matrix-method
#'   clusterSingle,ClusterExperiment-method clusterSingle,matrix,missing-method
#'   clusterSingle,matrixOrMissing,matrixOrMissing-method
#' @rdname clusterSingle
#' @export
setMethod(
  f = "clusterSingle",
  signature = signature(x = "missing",diss="matrixOrNULL"),
  definition = function(x, diss,...) {
    	clusterSingle(x=NULL,diss=diss,...)
 
})

#' @rdname clusterSingle
#' @export
setMethod(
  f = "clusterSingle",
  signature = signature(x = "matrixOrNULL",diss="missing"),
  definition = function(x, diss,...) {
    	clusterSingle(x=x,diss=NULL,...)
 
})

#' @rdname clusterSingle
#' @export
setMethod(
  f = "clusterSingle",
  signature = signature(x = "SummarizedExperiment", diss="missing"),
  definition = function(x, ...) {
    outval <- clusterSingle(assay(x),  ...)
    retval <- .addBackSEInfo(newObj=outval,oldObj=x)
    return(retval)
  }
)


#' @rdname clusterSingle
#' @param replaceCoClustering logical. If TRUE, the co-clustering resulting from subsampling is returned in the coClustering object and replaces any existing coClustering object in the slot \code{coClustering}. 
#' @export
setMethod(
  f = "clusterSingle",
  signature = signature(x = "ClusterExperiment", diss="missing"),
  definition = function(x, replaceCoClustering=FALSE,...) {

    outval <- clusterSingle(assay(x),transFun=transformation(x),...)
    retval<-addClusters(x,outval)
	#make most recent clustering the primary cluster
	primaryCluster(retval)<-nClusters(retval)
	if(replaceCoClustering & !is.null(outval@coClustering)) retval@coClustering<-outval@coClustering
	return(retval)
  }
)
#' @rdname clusterSingle
#' @export
setMethod(
  f = "clusterSingle",
  signature = signature(x = "matrixOrNULL",diss="matrixOrNULL"),
  definition = function(x, diss, subsample=TRUE, sequential=FALSE,
      clusterDArgs=NULL, subsampleArgs=NULL, seqArgs=NULL, 
      isCount=FALSE,transFun=NULL, dimReduce=c("none","PCA","var","cv","mad"),
      ndims=NA,clusterLabel="clusterSingle",checkDiss=TRUE) {
	# if(missing(x)) x<-NULL
	# if(missing(diss)) diss<-NULL

    ##########
    ##Check arguments and set defaults as needed
	##Note, some checks are duplicative of internal, but better here, because don't want to find error after already done extensive calculation...
    ##########
    	checkOut<-.checkSubsampleClusterDArgs(x=x,diss=diss,subsample=subsample,sequential=sequential,clusterDArgs=clusterDArgs,subsampleArgs=subsampleArgs,checkDiss=checkDiss)
	if(is.character(checkOut)) stop(checkOut)
	else {
		clusterDArgs<-checkOut$clusterDArgs
		subsampleArgs<-checkOut$subsampleArgs
		input<-checkOut$inputClusterD
	}

    if(sequential){
      if(is.null(seqArgs)) {
		  ##To DO: Question: if missing seqArgs, should we grab k0 from subsampleArgs?
        stop("if sequential=TRUE, must give seqArgs so as to identify k0")
      }
      if(!"k0"%in%names(seqArgs)) {
        stop("seqArgs must contain element 'k0'")
      }
    }
	##########
	## Handle dimensionality reduction:
	##########
	###Don't do this until do the checks, because takes some time.
    if(input %in% c("X")){
      N <- dim(x)[2]
      origX <- x #ngenes x nsamples
      ##########
      ##transformation to data x that will be input to clustering
      ##########
      dimReduce <- match.arg(dimReduce) #should be only 1
      if(length(ndims)>1) {
        stop("clusterSingle only handles one choice of dimensions. If you want to compare multiple choices, try clusterMany")
      }
      if(!is.na(ndims) & dimReduce=="none") {
        warning("specifying ndims has no effect if dimReduce==`none`")
      }
      nPCADims <- ifelse(dimReduce=="PCA", ndims, NA)
      nVarDims <- ifelse(dimReduce %in% c("var","cv","mad"), ndims, NA)
      transObj <- .transData(x, nPCADims=nPCADims, nVarDims=nVarDims,
                             dimReduce=dimReduce, transFun=transFun,
                             isCount=isCount)
      x <- transObj$x
      if(is.null(dim(x)) || NCOL(x)!=NCOL(origX)) {
        stop("Error in the internal transformation of x")
      }
      transFun <- transObj$transFun #need it later to create clusterExperimentObject
      
    }
    else{
      if(dimReduce!="none") stop("dimReduce only applies when diss not given or clusterFunction object doesn't accept the given diss as input")
	  N<-nrow(diss)
    }
    if(input %in% c("both","diss") && !is.null(clusterDArgs) && "distFunction" %in% names(clusterDArgs)){
        if(!is.na(clusterDArgs[["distFunction"]])) stop("if give diss as input to clusterSingle, cannot specify 'distFunction' in clusterDArgs")
    }
	
	
	##########
	## Start running clustering
	##########
    if(sequential){
      outlist <- do.call("seqCluster",
                        c(list(x=x, diss=diss,subsample=subsample,
                               subsampleArgs=subsampleArgs,
                               clusterDArgs=clusterDArgs), seqArgs))
    }
    else{
      ##########
      ##.clusterWrapper just deciphers choices and makes clustering.
      ##########
      finalClusterList <- .clusterWrapper(x=x, diss=diss, 
                                          subsample=subsample,
                                          subsampleArgs=subsampleArgs,
                                          clusterDArgs=clusterDArgs)
      outlist <- list("clustering"=.convertClusterListToVector(finalClusterList$result, N))

    }
    clInfo<-list(list(clusterInfo = outlist$clusterInfo,
                      whyStop = outlist$whyStop,
                      subsample = subsample,
                      sequential = sequential,
                      clusterFunction = clusterFunction,
                      clusterDArgs = clusterDArgs,
                      subsampleArgs = subsampleArgs,
                      seqArgs = seqArgs,
                      dimReduce=dimReduce,
                      ndims=ndims
    ))
    ##########
    ## Convert to clusterExperiment Object
    ##########
    if(input %in% c("X")){
      retval <- clusterExperiment(origX, outlist$clustering,
                                  transformation=transFun,
                                  clusterInfo=clInfo,
                                  clusterTypes="clusterSingle")
      clusterLabels(retval)<-clusterLabel
      if(!sequential & subsample) {
        retval@coClustering<-1-finalClusterList$diss
      }
      validObject(retval)
      return(retval)
    }
    else{
      out<-list(clustering=outlist$clustering,clusterInfo=clInfo,diss=outlist$diss)
    }

  }
)
#wrapper that calls the clusterSampling and clusterD routines in reasonable order.
#called by both seqCluster and clusterSingle
#clusterFunction assumed to be in clusterDArgs and subsampleArgs 
.clusterWrapper <- function(x, diss, subsample, clusterDArgs=NULL, subsampleArgs=NULL) 
{
    if(subsample){
		# ##To DO: Need to revisit why this was here and if add to .checkSubsampleClusterDArgs
		# #not clear why need this check...
		#         if(is.null(subsampleArgs) || is.null(subsampleArgs[["clusterArgs"]]) || !"k" %in% names(subsampleArgs[["clusterArgs"]])) stop("must provide k in 'subsampleArgs' via the 'clusterArgs' argument (or if sequential should have been set by sequential strategy)")
        Dbar<-do.call("subsampleClustering",c(list(x=x),subsampleArgs))
        diss<-1-Dbar #make it a distance.
        x<-NULL
		##To Do: revisit why this is here and if add to .checkSubsampleClusterDArgs
        # if(typeAlg=="K"){
        #     if(is.null(clusterDArgs)) clusterDArgs<-list(k=subsampleArgs[["k"]])
        #     else if(!"k" %in% names(clusterDArgs)) clusterDArgs[["k"]]<-subsampleArgs[["k"]] #either sequential sets this value, or get error in subsampleClustering, so always defined.
        # }
    }
	####To DO:Need to revisit why this was here and if add to .checkSubsampleClusterDArgs
		#     if(typeAlg=="K"){
		# ###Why is this here??? why does findBestK have to be set? Why can't it just be missing???
		# ###Is this simply so I can check that 'k' is defined in clusterDArgs??
		#         findBestK<-FALSE
		#         if(!is.null(clusterDArgs) && "findBestK" %in% names(clusterDArgs)){
		#             findBestK<-clusterDArgs[["findBestK"]]
		#         }
		# ###Ditto. Why is this here??? Aren't there already existing checks??
		#         if(is.null(clusterDArgs) || (!"k" %in% names(clusterDArgs) && !findBestK)) stop("if not type 'K' algorithm, must give k in 'clusterDArgs' (or if sequential should have been set by sequential strategy)")
		#     }
    resList<-do.call("clusterD",c(list(x=x,diss=diss,format="list", returnData=TRUE),clusterDArgs)) 
    return(resList) 
}






