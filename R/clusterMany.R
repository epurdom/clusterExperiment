#' Create a matrix of clustering across values of parameters
#'
#' Given a range of parameters, this function will return a matrix with the
#' clustering of the samples across the range, which can be passed to
#' \code{plotClusters} for visualization.
#'
#' @aliases clusterMany
#'
#' @param x the data matrix on which to run the clustering. Can be: matrix (with
#'   genes in rows), a list of datasets overwhich the clusterings should be run,
#'   a \code{SummarizedExperiment} object, or a \code{ClusterExperiment} object.
#' @param ks the range of k values (see details for the meaning of \code{k} for
#'   different choices of other parameters).
#' @param alphas values of alpha to be tried. Only used for clusterFunctions of
#'   type '01'. Determines tightness required in creating clusters from the
#'   dissimilarity matrix. Takes on values in [0,1]. See documentation of
#'   \code{\link{ClusterFunction}}.
#' @param betas values of \code{beta} to be tried in sequential steps. Only used
#'   for \code{sequential=TRUE}. Determines the similarity between two clusters
#'   required in order to deem the cluster stable. Takes on values in [0,1]. See
#'   documentation of \code{\link{seqCluster}}.
#' @param clusterFunction function used for the clustering. Note that unlike in
#'   \code{\link{clusterSingle}}, this must be a character vector of pre-defined
#'   clustering techniques, and can not be a user-defined function. Current
#'   functions can be found by typing \code{listBuiltInFunctions()} into the
#'   command-line.
#' @param minSizes the minimimum size required for a cluster (in the
#'   \code{mainClustering} step). Clusters smaller than this are not kept and samples
#'   are left unassigned.
#' @param distFunction a vector of character strings that are the names of
#'   distance functions found in the global environment. See the help pages of
#'   \code{\link{clusterSingle}} for details about the required format of
#'   distance functions. Currently, this distance function must be applicable
#'   for all clusterFunction types tried. Therefore, it is not possible in
#'   \code{clusterMany} to intermix type "K" and type "01" algorithms if you
#'   also give distances to evaluate via \code{distFunction} unless all
#'   distances give 0-1 values for the distance (and hence are possible for both
#'   type "01" and "K" algorithms).
#' @param nVarDims vector of the number of the most variable features to keep
#'   (when "var", "cv", or "mad" is identified in \code{dimReduce}). If NA is
#'   included, then the full dataset will also be included.
#' @param nPCADims vector of the number of PCs to use (when 'PCA' is identified
#'   in \code{dimReduce}). If NA is included, then the full dataset will also be
#'   included.
#' @param eraseOld logical. Only relevant if input \code{x} is of class
#'   \code{ClusterExperiment}. If TRUE, will erase existing workflow results
#'   (clusterMany as well as mergeClusters and combineMany). If FALSE, existing
#'   workflow results will have "\code{_i}" added to the clusterTypes value,
#'   where \code{i} is one more than the largest such existing workflow
#'   clusterTypes.
#' @param findBestK logical, whether should find best K based on average
#'   silhouette width (only used when clusterFunction of type "K").
#' @param silCutoff Requirement on minimum silhouette width to be included in
#'   cluster (only for combinations where removeSil=TRUE).
#' @param removeSil logical as to whether remove when silhouette < silCutoff
#'   (only used if clusterFunction of type "K")
#' @inheritParams clusterSingle
#' @inheritParams mainClustering
#' @param ncores the number of threads
#' @param random.seed a value to set seed before each run of clusterSingle (so
#'   that all of the runs are run on the same subsample of the data). Note, if
#'   'random.seed' is set, argument 'ncores' should NOT be passed via
#'   subsampleArgs; instead set the argument 'ncores' of clusterMany directly
#'   (which is preferred for improving speed anyway).
#' @param run logical. If FALSE, doesn't run clustering, but just returns matrix
#'   of parameters that will be run, for the purpose of inspection by user (with
#'   rownames equal to the names of the resulting column names of clMat object
#'   that would be returned if \code{run=TRUE}). Even if \code{run=FALSE},
#'   however, the function will create the dimensionality reductions of the data
#'   indicated by the user input.
#' @param ... For signature \code{list}, arguments to be passed on to mclapply
#'   (if ncores>1). For all the other signatures, arguments to be passed to the
#'   method for signature \code{list}.
#' @param verbose logical. If TRUE it will print informative messages.
#' @details Some combinations of these parameters are not feasible. See the
#'   documentation of \code{\link{clusterSingle}} for important information on
#'   how these parameter choices interact.
#' @details While the function allows for multiple values of clusterFunction,
#'   the code does not reuse the same subsampling matrix and try different
#'   clusterFunctions on it. This is because if sequential=TRUE, different
#'   subsample clusterFunctions will create different sets of data to subsample
#'   so it is not possible; if sequential=FALSE, we have not implemented
#'   functionality for this reuse. Setting the \code{random.seed} value,
#'   however, should mean that the subsampled matrix is the same for each, but
#'   there is no gain in computational complexity (i.e. each subsampled
#'   co-occurence matrix is recalculated for each set of parameters).
#'
#' @details The argument \code{ks} is interpreted differently for different
#'   choices of the other parameters. When/if sequential=TRUE, \code{ks} defines
#'   the argument \code{k0} of \code{\link{seqCluster}}. Otherwise, \code{ks}
#'   values are the \code{k} values for \strong{both} the mainClustering and
#'   subsampling step (i.e. assigned to the \code{subsampleArgs} and
#'   \code{mainClusterArgs} that are passed to \code{\link{mainClustering}} and
#'   \code{\link{subsampleClustering}} unless \code{k} is set appropriately in
#'   \code{subsampleArgs}. The passing of these arguments via
#'   \code{subsampleArgs} will only have an effect if `subsample=TRUE`.
#'   Similarly, the passing of \code{mainClusterArgs[["k"]]} will only have an
#'   effect when the clusterFunction argument includes a clustering algorithm of
#'   type "K". When/if "findBestK=TRUE", \code{ks} also defines the
#'   \code{kRange} argument of \code{\link{mainClustering}} unless \code{kRange} is
#'   specified by the user via the \code{mainClusterArgs}; note this means that the
#'   default option of setting \code{kRange} that depends on the input \code{k}
#'   (see \code{\link{mainClustering}}) is not available in \code{clusterMany}, only
#'   in \code{\link{clusterSingle}}.
#' @details If the input is a \code{ClusterExperiment} object, current
#'   implementation is that existing \code{orderSamples},\code{coClustering} or
#'   the many dendrogram slots will be retained.
#' @return If \code{run=TRUE} and the input is either a matrix, a
#'   \code{SummarizedExperiment} object, or a \code{ClusterExperiment} object,
#'   will return a \code{ClusterExperiment} object, where the results are stored
#'   as clusterings with clusterTypes \code{clusterMany}. Depending on
#'   \code{eraseOld} argument above, this will either delete existing such
#'   objects, or change the clusterTypes of existing objects. See argument
#'   \code{eraseOld} above. Arbitrarily the first clustering is set as the
#'   primaryClusteringIndex.
#'
#' @return If \code{run=TRUE} and the input is a list of data sets, a list with
#'   the following objects: \itemize{ \item{\code{clMat}}{ a matrix with each
#'   column corresponding to a clustering and each row to a sample.}
#'   \item{\code{clusterInfo}}{ a list with information regarding clustering
#'   results (only relevant entries for those clusterings with sequential=TRUE)}
#'   \item{\code{paramMatrix}}{ a matrix giving the parameters of each
#'   clustering, where each column is a possible parameter set by the user and
#'   passed to \code{\link{clusterSingle}} and each row of paramMatrix
#'   corresponds to a clustering in \code{clMat}} \item{\code{mainClusterArgs}}{ a
#'   list of (possibly modified) arguments to mainClusterArgs}
#'   \item{\code{seqArgs=seqArgs}}{a list of (possibly modified) arguments to
#'   seqArgs} \item{\code{subsampleArgs}}{a list of (possibly modified)
#'   arguments to subsampleArgs} }
#' @return If \code{run=FALSE} a list similar to that described above, but
#'   without the clustering results.
#'
#' @examples
#' data(simData)
#'
#' #Example: clustering using pam with different dimensions of pca and different
#' #k and whether remove negative silhouette values
#' #check how many and what runs user choices will imply:
#' checkParams <- clusterMany(simData,nPCADims=c(5,10,50),  dimReduce="PCA",
#' clusterFunction="pam",
#' ks=2:4,findBestK=c(TRUE,FALSE),removeSil=c(TRUE,FALSE),run=FALSE)
#' print(head(checkParams$paramMatrix))
#'
#' #Now actually run it
#' cl <- clusterMany(simData,nPCADims=c(5,10,50),  dimReduce="PCA",
#' clusterFunction="pam",ks=2:4,findBestK=c(TRUE,FALSE),removeSil=c(TRUE,FALSE))
#' print(cl)
#' head(colnames(clusterMatrix(cl)))
#'
#' #make names shorter for plotting
#' clMat <- clusterMatrix(cl)
#' colnames(clMat) <- gsub("TRUE", "T", colnames(clMat))
#' colnames(clMat) <- gsub("FALSE", "F", colnames(clMat))
#' colnames(clMat) <- gsub("k=NA,", "", colnames(clMat))
#'
#' par(mar=c(2, 10, 1, 1))
#' plotClusters(clMat, axisLine=-2)
#'
#'
#' \dontrun{
#'	#following code takes around 1+ minutes to run because of the subsampling
#'	#that is redone each time:
#'	system.time(clusterTrack <- clusterMany(simData, ks=2:15,
#'	alphas=c(0.1,0.2,0.3), findBestK=c(TRUE,FALSE), sequential=c(FALSE),
#'	subsample=c(FALSE), removeSil=c(TRUE), clusterFunction="pam",
#'	mainClusterArgs=list(minSize=5, kRange=2:15), ncores=1, random.seed=48120))
#' }
#'
#' @rdname clusterMany
#' @importFrom parallel mclapply
#' @export
# not currently functional:
# @param paramMatrix matrix or data.frame. If given, the algorithm will bypass
#   creating the matrix of possible parameters, and will use the given matrix.
#   There are basically no checks as to whether this matrix is in the right
#   format, and is only intended to be used to feed the results of setting
#   \code{run=FALSE} back into the algorithm (see example).
##Was in @examples:
# #get rid of some of the choices manually
# #note that the supplement arguments could have been changed too, so
# #we give those to clusterMany as well.
# checkParamsMat <- checkParams$paramMatrix[-c(1,2),]
#
# clSmaller <- clusterMany(simData, nPCADims=c(5,10,50),  dimReduce="PCA",
# paramMatrix=checkParamsMat, subsampleArgs=checkParams$subsampleArgs,
# seqArgs=checkParams$seqArgs, mainClusterArgs=checkParams$mainClusterArgs)
#' @export
setMethod(
  f = "clusterMany",
  signature = signature(x = "matrix"),
  definition = function(x,
                        dimReduce="none",nVarDims=NA,nPCADims=NA,
                        transFun=NULL,isCount=FALSE,
                        ...
  ){

	  if(any(dim(x)==0)) stop("x must have non zero dimensions")
    origX <- x
    transObj <- .transData(x, nPCADims=nPCADims, nVarDims=nVarDims,
                           dimReduce=dimReduce, transFun=transFun,
                           isCount=isCount)
    x <- transObj$x
    if(!is.null(dim(x)) && NCOL(x)!=NCOL(origX)) {
      stop("Error in the internal transformation of x")
    }
    transFun <- transObj$transFun #need it later to create clusterExperimentObject

    if(!is.null(dim(x))) {
      x <- list(dataset1=x)
    } #if npcs=NA, then .transData returns a matrix.
    outval <- clusterMany(x, ...)
    ##########
    ## Convert to clusterExperiment Object
    ##########
    if("clMat" %in% names(outval)) {
      retval <- clusterExperiment(origX, outval$clMat,
                                  transformation=transFun,
                                  clusterInfo=outval$clusterInfo,
                                  clusterTypes="clusterMany")
      validObject(retval)
      return(retval)
    }
    else {
      return(outval)
    }
  }

)

#' @rdname clusterMany
#' @export
setMethod(
  f = "clusterMany",
  signature = signature(x = "list"),
  definition = function(x, ks=NA, clusterFunction, alphas=0.1, findBestK=FALSE,
                        sequential=FALSE, removeSil=FALSE, subsample=FALSE,
                        silCutoff=0, distFunction=NA,
                        betas=0.9, minSizes=1,
                        verbose=FALSE,
                        mainClusterArgs=NULL,
                        subsampleArgs=NULL,
                        seqArgs=NULL,
                        ncores=1, random.seed=NULL, run=TRUE,
                        ...
  )
  {
    paramMatrix<-NULL
    data <- x
    if(!is.null(random.seed)){
        if(!is.null(subsampleArgs) && "ncores" %in% names(subsampleArgs)){
            if(subsampleArgs[["ncores"]]>1) stop("setting random.seed will not be reproducible if ncores given to subsampleArgs")
        }
    }
    if(!all(sapply(data, function(y){is.matrix(y) || is.data.frame(y)}))) {
      stop("if data is a list, it must be a list with each element of the list a data.frame or matrix")
    }
    #check all same number of observations: Why do we have this check?? Why does it matter??
    if(!length(unique(sapply(data,NCOL)))==1) {
      stop("All data sets must have the same number of observations")
    }
    if(is.null(names(data))) {
      names(data) <- paste("dataset",1:length(data),sep="")
    }
    dataList<-data
    dataName <- names(dataList)
    if(is.null(paramMatrix)){
      param <- expand.grid(dataset=dataName,
                         k=ks, alpha=alphas, findBestK=findBestK, beta=betas, minSize=minSizes,
                         sequential=sequential, distFunction=distFunction,
                         removeSil=removeSil, subsample=subsample,
                         clusterFunction=clusterFunction, silCutoff=silCutoff)
      ###########
      #Check param matrix:
      #don't vary them across ones that don't matter (i.e. 0-1 versus K);
      #code sets to single value and then will do unique
      #also deals with just in case the user gave duplicated values of something by mistake.
      ###########
	  paramAlgTypes<-algorithmType(param[,"clusterFunction"])
	  if(length(paramAlgTypes)!=nrow(param)) stop("Internal coding error in clusterMany: not getting right number of type of algorithms from param")
      typeK <- which( paramAlgTypes=="K")
      if(length(typeK)>0){
        param[typeK,"alpha"] <- NA #just a nothing value, because doesn't mean anything here
        #--------
		#if findBestK make sure other arguments make sense:
        #--------
        whFindBestK <- which(param[,"findBestK"])
        if(length(whFindBestK)>0){
          #by default make kRange in mainClustering equal to the ks. Note this will be true of ALL
          if(!"kRange" %in% names(mainClusterArgs)) {
            mainClusterArgs[["kRange"]]<-ks
          }
          #if findBestK=TRUE, and sequential=FALSE, then need to set 'k'=NA
          whNoSeq <- which(!param[,"sequential"])
          if(length(intersect(whFindBestK,whNoSeq))>0){
            param[intersect(whFindBestK,whNoSeq),"k"] <- NA
          }
		  #and if subsample=TRUE, then user needs to set k via subsampleArgs
		  ##Might could handle this better by call to .checkSubsampleClusterDArgs
          whNoSeqSub <- which(!param[,"sequential"] & param[,"subsample"])
          if(length(intersect(whFindBestK,whNoSeqSub))>0 &
             is.null(subsampleArgs[["clusterArgs"]]) && is.null(subsampleArgs[["clusterArgs"]][["k"]])){
            stop("must provide k in 'clusterArgs' element of 'subsampleArgs' because there are combinations of findBestK=TRUE, sequential=FALSE and subsample=TRUE. (Note this will set 'k' for all combinations that subsample, not just this parameter combinations)")
          }
        }
      }
      type01 <- which( paramAlgTypes=="01")
      if(length(type01)>0){
        param[type01,"findBestK"] <- FALSE
        param[type01,"removeSil"] <- FALSE
        param[type01,"silCutoff"] <- 0
      }
	  ##Turn off distFunction for those that subsample, because will use that of co-occurance
      whSubsample<-which(param[,"subsample"])
      if(length(whSubsample)>0){
        param[whSubsample,"distFunction"]<-NA
      }
	  ###Check value alpha, beta values
      alpha01 <- which(param[,"alpha"]<0 | param[,"alpha"]>1)
	  if(length(alpha01)>0){
		  warning("alpha value must be in (0,1). Input alphas outside that range are ignored")
		  param[alpha01,"alpha"]<-NA
	  }
      beta01 <- which(param[,"beta"]<0 | param[,"beta"]>1)
	  if(length(beta01)>0){
		  warning("beta value must be in (0,1). Input betas outside that range are ignored")
		  param[beta01,"beta"]<-NA
	  }

      param <- unique(param)

      #####
      #deal with those that are invalid combinations:
	  # Might could handle this better by call to .checkSubsampleClusterDArgs for each parameter combination
	  # Also, if ever reinstate param option, then should apply these checks to that param
      ######
      whInvalid <- which(!param[,"subsample"] & param[,"sequential"]
                         & param[,"findBestK"])
      if(length(whInvalid)>0) {
        param <- param[-whInvalid,]
      }

      whInvalid <- which(!param[,"subsample"] & param[,"sequential"]
                       & param[,"findBestK"])
      if(length(whInvalid)>0) {
        param<-param[-whInvalid,]
      }

      whInvalid <- which(param[,"sequential"] & is.na(param[,"beta"]))
      if(length(whInvalid)>0) {
        param<-param[-whInvalid,]
      }

	  #if type K and not findBestK, need to give the k value.
      whInvalid <- which(is.na(param[,"k"]) & !param[,"findBestK"] & algorithmType(param[,"clusterFunction"])=="K" )
      if(length(whInvalid)>0){
			param<-param[-whInvalid,]
		}

      #####
      #require at least 2 combinations:
      #####
      if(nrow(param)<=1) {
        stop("set of parameters imply only 1 combination. If you wish to run a single clustering, use 'clusterSingle'")
      }

      #give names to the parameter combinations.
      whVary <- which(apply(param,2,function(x){length(unique(x))>1}))
      if(length(whVary)>0) {
        cnames<-apply(param[,whVary,drop=FALSE],1,function(x){
        paste(colnames(param)[whVary],x,sep="=",collapse=",")})
      } else {
        stop("set of parameters imply only 1 combination. If you wish to run a single clustering, use 'clusterSingle'")
      }
      cnames <- gsub("dataset=","",cnames)
      cnames <- gsub("= ","=",cnames)
      cnames[param[,"sequential"]] <- gsub("k=", "k0=", cnames[param[,"sequential"]])
      rownames(param) <- cnames
  } else{ #if paramMatrix!=NULL, have killed off this code for now, because doesn't work.
      if(!run) {
        stop("If paramMatrix is given, run should be TRUE. Otherwise there is no effect.")
      }
      if(is.null(paramMatrix)) {
        stop("invalid input for paramMatrix; must be data.frame or matrix")
      }
      param <- paramMatrix
      if(is.null(rownames(paramMatrix))) {
        stop("input paramMatrix must have row names")
      }
      cnames<-rownames(paramMatrix)
  }

    if(verbose) {
      cat(nrow(param),"parameter combinations,",sum(param[,"sequential"]),"use sequential method.\n")
    }
	if(is.null(mainClusterArgs)) mainClusterArgs<-list(clusterArgs=list())
	if(is.null(subsampleArgs)) subsampleArgs<-list(clusterArgs=list())
    paramFun <- function(i){
      par <- param[i,]
      #make them logical values... otherwise adds a space before the TRUE and doesn't recognize.
      #well, sometimes. Maybe don't need this?
      removeSil <- as.logical(gsub(" ","",par["removeSil"]))
      sequential <- as.logical(gsub(" ","",par["sequential"]))
      subsample <- as.logical(gsub(" ","",par["subsample"]))
      findBestK <- as.logical(gsub(" ","",par["findBestK"]))
      clusterFunction <- as.character(par[["clusterFunction"]])
      distFunction<-if(!is.na(par[["distFunction"]])) as.character(par[["distFunction"]]) else NULL
      if(!is.na(par[["k"]])){
        if(sequential) {
          seqArgs[["k0"]] <- par[["k"]]
        } else{
          #to be safe, set both in case user set one.
          subsampleArgs[["clusterArgs"]][["k"]] <- par[["k"]]
          mainClusterArgs[["clusterArgs"]][["k"]] <- par[["k"]]
        }
      }
      mainClusterArgs[["clusterArgs"]][["alpha"]] <- par[["alpha"]]
      seqArgs[["beta"]] <- par[["beta"]]
      mainClusterArgs[["minSize"]] <- par[["minSize"]]
      mainClusterArgs[["findBestK"]] <- findBestK
      mainClusterArgs[["removeSil"]] <- removeSil
      mainClusterArgs[["silCutoff"]] <- par[["silCutoff"]]
      mainClusterArgs[["checkArgs"]] <- FALSE #turn off printing of warnings that arguments off
	  mainClusterArgs[["clusterFunction"]]<-clusterFunction
      seqArgs[["verbose"]]<-FALSE
      if(!is.null(random.seed)) {
        set.seed(random.seed)
      }
	  ##Note that currently, checkDiss=FALSE, also turns off warnings about arguments
      if(!is.null(distFunction)){
        diss<- allDist[[paste(as.character(par[["dataset"]]),distFunction,sep="--")]]
        clusterSingle(x=dataList[[as.character(par[["dataset"]])]], diss=diss,subsample=subsample, dimReduce="none",
                      mainClusterArgs=mainClusterArgs,
                      subsampleArgs=subsampleArgs, seqArgs=seqArgs,
                      sequential=sequential, transFun=function(x){x},checkDiss=FALSE)       }
      else clusterSingle(x=dataList[[as.character(par[["dataset"]])]], subsample=subsample,
                 mainClusterArgs=mainClusterArgs, dimReduce="none",
                 subsampleArgs=subsampleArgs, seqArgs=seqArgs,
                 sequential=sequential, transFun=function(x){x},checkDiss=FALSE)
	    }
    if(run){
      ##Calculate distances necessary only once
      if(any(!is.na(param[,"distFunction"]))){
        distParam<-unique(param[,c("dataset","distFunction")])
        distParam<-distParam[!is.na(distParam[,"distFunction"]),]
        #browser()
          allDist<-lapply(1:nrow(distParam),function(ii){
            distFun<-as.character(distParam[ii,"distFunction"])
            dataName<-as.character(distParam[ii,"dataset"])
			algCheckType<-if(any(paramAlgTypes=="01")) "01" else "K" #be conservative and check for the 01 type if any of clusterFunctions are 01.
            distMat<-.makeDiss(dataList[[dataName]],distFunction=distFun,checkDiss=TRUE,algType=algCheckType)
			# fun<-get(distFun,envir=globalenv())
			#             distMat<-as.matrix(fun(t(dataList[[dataName]])))
			#             .checkDistFunction(distMat) #check it here!
            return(distMat)
          })
        names(allDist)<-paste(distParam[,"dataset"],distParam[,"distFunction"],sep="--")

      }

      if(verbose) {
        cat("Running Clustering on Parameter Combinations...")
      }
	  #browser()
      if(ncores>1) {
        out <- mclapply(1:nrow(param), FUN=paramFun, mc.cores=ncores, ...)
        nErrors <- which(sapply(out, function(x){inherits(x, "try-error")}))
        if(length(nErrors)>0) {
          stop(length(nErrors)," parameter values (of ",length(out),") hit an error. The first was:\n",out[nErrors[1]])
        }
      } else {
        out <- lapply(1:nrow(param),FUN=paramFun)
      }
      if(verbose) {
        cat("done.\n")
      }
      clMat <- sapply(out, function(x){primaryCluster(x)})

      colnames(clMat) <- unname(cnames)
      pList <- lapply(1:nrow(param), function(i){
        x <- param[i,]
        names(x) <- colnames(param)
        return(x)})
      clInfo <- mapply(pList, out, FUN=function(x, y){
        c(list(choicesParam=x), clusterInfo(y))
      }, SIMPLIFY=FALSE)

      return(list(clMat=clMat, clusterInfo=clInfo, paramMatrix=param,
                  mainClusterArgs=mainClusterArgs, seqArgs=seqArgs,
                  subsampleArgs=subsampleArgs))
    } else{
      if(verbose) {
        cat("Returning Parameter Combinations without running them (to run them choose run=TRUE)\n")
      }
      return(list(paramMatrix=param, mainClusterArgs=mainClusterArgs, seqArgs=seqArgs,
                  subsampleArgs=subsampleArgs))
    }
  }
)

#' @rdname clusterMany
#' @export
setMethod(
  f = "clusterMany",
  signature = signature(x = "ClusterExperiment"),
  definition = function(x, dimReduce="none", nVarDims=NA, nPCADims=NA,
                        eraseOld=FALSE, ...)
  {
    outval<-clusterMany(assay(x), dimReduce=dimReduce, nVarDims=nVarDims,
                        nPCADims=nPCADims, transFun=transformation(x), ...)
    if(is(outval, "ClusterExperiment")) {
      #outval<-.addBackSEInfo(newObj=outval,oldObj=x) #added to '.addNewResult'
      ##Check if clusterMany already ran previously
      x<-.updateCurrentWorkflow(x,eraseOld,"clusterMany")

      if(!is.null(x)) {
        retval <- .addNewResult(newObj=outval,oldObj=x) #make decisions about what to keep.
      } else {
        retval<-.addBackSEInfo(newObj=outval,oldObj=x)
      }

      validObject(retval)
      return(retval)

    } else {
      return(outval)
    }
  }
)


#' @rdname clusterMany
#' @export
setMethod(
  f = "clusterMany",
  signature = signature(x = "SummarizedExperiment"),
  definition = function(x, dimReduce="none", nVarDims=NA, nPCADims=NA,
                        transFun=NULL, isCount=FALSE, ...)
  {
    outval <- clusterMany(assay(x), dimReduce=dimReduce, nVarDims=nVarDims,
                          nPCADims=nPCADims, transFun=transFun, isCount=isCount,
                          ...)
    if(class(outval)=="ClusterExperiment") {
        retval<-.addBackSEInfo(newObj=outval,oldObj=x)

        return(retval)
    }  #need to redo it to make sure get any other part of summarized experiment
    else {
      return(outval)
    }
  }
)
#' @export
#' @rdname clusterMany
setMethod(
f = "clusterMany",
signature = signature(x = "data.frame"),
definition = function(x,...){clusterMany(data.matrix(x),...)}
)



