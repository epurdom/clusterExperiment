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
#' @param clusterFunction passed to \code{\link{clusterD}} option
#'   'clusterFunction' to indicate method of clustering, see
#'   \code{\link{clusterD}}.
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
#' set.seed(44261)
#' clustSeqHier_v2 <- clusterSingle(simData, clusterFunction="hierarchical01",
#' sequential=TRUE, subsample=TRUE, subsampleArgs=list(resamp.n=100, samp.p=0.7,
#' clusterFunction="kmeans", clusterArgs=list(nstart=10)),
#' seqArgs=list(beta=0.8, k0=5), clusterDArgs=list(minSize=5))
#' }
#'
#' #use clusterSingle to do just clustering k=3 with no subsampling
#' clustNothing <- clusterSingle(simData, clusterFunction="pam",
#' subsample=FALSE, sequential=FALSE, clusterDArgs=list(k=3))
#' @aliases clusterSingle clusterSingle-methods clusterSingle,matrix-method
#'   clusterSingle,ClusterExperiment-method clusterSingle,matrix,missing-method
#'   clusterSingle,matrixOrMissing,matrixOrMissing-method
#' @rdname clusterSingle
#' @export
setMethod(
  f = "clusterSingle",
  signature = signature(x = "matrixOrMissing",diss="matrixOrMissing",clusterFunction="character"),
  definition = function(x, diss,clusterFunction,...) {
    	clusterSingle(x=x,diss=diss,clusterFunction=getBuiltInClusterFunction(clusterFunction),...)
 
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
#' @export
setMethod(
  f = "clusterSingle",
  signature = signature(x = "ClusterExperiment", diss="missing"),
  definition = function(x, ...) {

    outval <- clusterSingle(assay(x),...)
    
    ## eap: I think we should add it, so I changed it here. You might try a couple of versions.
    retval<-addClusters(outval, x) #should keep primary cluster as most recent, so outval first
    return(retval)
  }
)
#' @rdname clusterSingle
#' @export
setMethod(
  f = "clusterSingle",
  signature = signature(x = "matrixOrMissing",diss="matrixOrMissing",clusterFunction="ClusterFuntion"),
  definition = function(x, diss,clusterFunction,subsample=TRUE, sequential=FALSE,
      clusterDArgs=NULL, subsampleArgs=NULL, seqArgs=NULL, 
      isCount=FALSE,transFun=NULL, dimReduce=c("none","PCA","var","cv","mad"),
      ndims=NA,clusterLabel="clusterSingle",checkDiss=TRUE) {
    if(missing(x)) x<-NULL
    if(missing(diss)) diss<-NULL
	#Following input commands will return only X or Diss because gave the inputType argument...
	input<-.checkXDissInput(x,diss,inputType=inputType(clusterFunction),algType=algorithmType(clusterFunction),checkDiss=checkDiss)
    algType<-algorithmType(clusterFunction)
	
	#check if clusterFunction given by subsampleClustering args
	if(subsample){
		if("clusterFunction" %in% names(subsampleArgs)){
			subsampleCF<-subsampleArgs[["clusterFunction"]
	    	subsampleAlgType<-algorithmType(subsampleCF)
		
			inputSubsample<-.checkXDissInput(x,diss, inputType=inputType(subsampleCF),  algType=algorithmType(subsampleCF), checkDiss=checkDiss) #if algorithm on one is 01 and other isn't, need to check diss again.
			diffSubsampleCF<-TRUE
		}
		else subsampleCF<-NULL
		# else{
		# 	subsampleCF<-clusterFunction
		# 	inputSubsample<-input
		# 	diffSubsampleCF<-FALSE
		# }
	
	}
    if(input %in% c("X")){
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
      N <- dim(x)[2]
      
    }
    else{
      if(dimReduce!="none") stop("dimReduce only applies when diss not given or clusterFunction object doesn't accept the given diss as input")
	  N<-nrow(diss)
    }
    if(input %in% c("both","diss") && !is.null(clusterDArgs) && "distFunction" %in% names(clusterDArgs)){
        if(!is.na(clusterDArgs[["distFunction"]])) stop("if give diss as input to clusterSingle, cannot specify 'distFunction' in clusterDArgs")
    }
	
    ##########
    ##Checks arguments and set defaults as needed:
    ##########
    if(algType == "K"){
      if("findBestK" %in% names(clusterDArgs) & !subsample & sequential){
        if(clusterDArgs[["findBestK"]])
          stop("Cannot do sequential clustering where subsample=FALSE and 'findBestK=TRUE' is passed via clusterDArgs. See help documentation.")
      }
    }
    if(sequential){
	  #if sequential, will run seqCluster
      if(is.null(seqArgs)) {
		  ##To DO:
		  #Question: if missing seqArgs, should we grab k0 from subsampleArgs?
        stop("if sequential=TRUE, must give seqArgs so as to identify k0")
      }
      if(!"k0"%in%names(seqArgs)) {
        stop("seqArgs must contain element 'k0'")
    }
    if(subsample){ 
        if(!is.null(clusterDArgs) && "distFunction" %in% names(clusterDArgs) && !is.na(clusterDArgs[["distFunction"]])){
            warning("if 'subsample=TRUE', 'distFunction' argument in clusterDArgs is ignored.")
            clusterDArgs[["distFunction"]]<-NA
        }
		##------
		##Check generic required args to determine whether should 'borrow' those args from clusterDArgs.
		## only if not sequential because sequential sets k for the subsampling via k0
		##------
		if(!sequential & "clusterArgs" %in% names(clusterDArgs)){
			clusterDReqArgs<-requiredArgs(clusterFunction)
			clusterDReqArgs<-clusterDReqArgs[clusterDReqArgs%in%clusterDArgs[["clusterArgs"]]]
			#find required args, either from subsampleCF if exists or clusterFunction
			if(!is.null(subsampleCF) && algorithmType(subsampleCF)==algorithmType(clusterFunction)){
				subPassedReqArgs<-requiredArgs(subsampleCF) 
			}
			else{ 
				subPassedReqArgs<-requiredArgs(clusterFunction)
			}
			if(!is.null(subsampleArgs) && "clusterArgs" %in% names(subsampleArgs)){
				#check if existing clusterArgs has required names already
				#if not, give them those of clusterD if exist.
				if(!all(subPassedReqArgs %in% names(subsampleArgs[["clusterArgs"]]))) {
					missingArgs<-subPassedReqArgs[!subPassedReqArgs%in%names(subsampleArgs[["clusterArgs"]])]
					missingArgs<-missingArgs[missingArgs%in%clusterDReqArgs]
    				subsampleArgs[["clusterArgs"]][missingArgs]<-clusterDArgs[["clusterArgs"]][clusterDReqArgs]
			    }
			} 	
			else{
				subsampleArgs[["clusterArgs"]]<-clusterDArgs[["clusterArgs"]][clusterDReqArgs]
			}
		}

		##To DO:
		# #Need to consider whether to include these type of warnings or let them be caught later:
		# Messages will be confusing if not done here
		#                             warning("did not give 'k' in 'subsampleArgs'.
		#                     Set to 'k' argument in 'clusterDArgs'")
		#             stop("if not sequential and do subsampling,
		#                  must pass 'k' in subsampleArgs")
		#
		
    }
	##To DO:
	# Similarly, consider whether should do these checks previously had, or just let internal functions.
	# Messages will be confusing if not done here. Maybe best is to make those be good messages!
    # else if(algType=="K" && !is.null(clusterDArgs) && !"k" %in% names(clusterDArgs)){
    #   #if don't specify k, then must have findBestK=TRUE in clusterDArgs;
    #   #is by default, so only need to check that if specified it,
    #   #set it to TRUE
    #   if("findBestK" %in% names(clusterDArgs) && !clusterDArgs[["findBestK"]])
    #     stop("if not sequential and clusterFunction is of type 'K' (e.g. pam)
    #          and findBestK=FALSE in clusterDArgs, must pass 'k' via
    #          clusterDArgs list")
    # }
	
	##########
	## Start running clustering
	##########
    if(sequential){
      outlist <- do.call("seqCluster",
                        c(list(x=x, diss=diss,subsample=subsample,
                               subsampleArgs=subsampleArgs,
                               clusterDArgs=clusterDArgs,
                               clusterFunction=clusterFunction), seqArgs))
    }
    else{
      ##########
      ##.clusterWrapper just deciphers choices and makes clustering.
      ##########
      finalClusterList <- .clusterWrapper(x=x, diss=diss, clusterFunction=clusterFunction,
                                          subsample=subsample,
                                          subsampleArgs=subsampleArgs,
                                          clusterDArgs=clusterDArgs)
      outlist <- list("clustering"=.convertClusterListToVector(finalClusterList$results, N))

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
    if(input %in% c("X","both")){
      retval <- clusterExperiment(origX, outlist$clustering,
                                  transformation=transFun,
                                  clusterInfo=clInfo,
                                  clusterTypes="clusterSingle")
      clusterLabels(retval)<-clusterLabel
      if(!sequential) {
        retval@coClustering<-finalClusterList$D
      }
      validObject(retval)
      return(retval)
    }
    else{
      out<-list(clustering=outlist$clustering,clusterInfo=clInfo)
    }

  }
)
#wrapper that calls the clusterSampling and clusterD routines in reasonable order.
#called by both seqCluster and clusterSingle
.clusterWrapper <- function(x, diss, subsample, clusterFunction,clusterDArgs=NULL,
                            subsampleArgs=NULL) 
{
    if(subsample){
		# ##To DO: Need to revisit why this was here.
		# #not clear why need this check...
		#         if(is.null(subsampleArgs) || is.null(subsampleArgs[["clusterArgs"]]) || !"k" %in% names(subsampleArgs[["clusterArgs"]])) stop("must provide k in 'subsampleArgs' via the 'clusterArgs' argument (or if sequential should have been set by sequential strategy)")
        Dbar<-do.call("subsampleClustering",c(list(x=x),subsampleArgs))
        diss<-1-Dbar #make it a distance.
        x<-NULL
        if(typeAlg=="K"){
            if(is.null(clusterDArgs)) clusterDArgs<-list(k=subsampleArgs[["k"]])
            else if(!"k" %in% names(clusterDArgs)) clusterDArgs[["k"]]<-subsampleArgs[["k"]] #either sequential sets this value, or get error in subsampleClustering, so always defined.
        }
    }
	####To DO:Need to revisit why this was here. 
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
    resList<-do.call("clusterD",c(list(x=x,diss=diss,format="list", clusterFunction=clusterFunction,returnData=TRUE),clusterDArgs)) 
    return(resList) 
}






