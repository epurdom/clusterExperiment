#' Create a matrix of clustering across values of parameters
#' 
#' Given a range of parameters, this function will return a matrix with the 
#' clustering of the samples across the range, which can be passed to 
#' \code{plotClusters} for visualization.
#' 
#' @aliases clusterMany
#'   
#' @param x the data matrix on which to run the clustering. Can be object of the
#'   following classes: matrix (with genes in rows), 
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}},
#'   \code{\link[SingleCellExperiment]{SingleCellExperiment}} 
#'   \code{\link{SingleCellFilter}}, or \code{ClusterExperiment}.
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
#' @param nFilter vector of the number of the most variable features to keep 
#'   (when "var", "abscv", or "mad" is identified in \code{dimReduce}). If NA is 
#'   included, then the full dataset will also be included.
#' @param nDimReduce vector of the number of PCs to use (when 'PCA' is identified 
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
#'   \code{kRange} argument of \code{\link{mainClustering}} unless \code{kRange}
#'   is specified by the user via the \code{mainClusterArgs}; note this means
#'   that the default option of setting \code{kRange} that depends on the input
#'   \code{k} (see \code{\link{mainClustering}}) is not available in
#'   \code{clusterMany}, only in \code{\link{clusterSingle}}.
#' @details If the input is a \code{ClusterExperiment} object, current 
#'   implementation is that existing \code{orderSamples},\code{coClustering} or 
#'   the many dendrogram slots will be retained.
#' @return If \code{run=TRUE} and the input is not a list of data matrices, will
#'   return a \code{ClusterExperiment} object, where the results are stored as
#'   clusterings with clusterTypes \code{clusterMany}. Depending on 
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
#'   corresponds to a clustering in \code{clMat}} \item{\code{mainClusterArgs}}{
#'   a list of (possibly modified) arguments to mainClusterArgs} 
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
#' checkParams <- clusterMany(simData,dimReduce="PCA",  
#' nDimReduce=c(5,10,50), clusterFunction="pam", isCount=FALSE,
#' ks=2:4,findBestK=c(TRUE,FALSE),removeSil=c(TRUE,FALSE),run=FALSE)
#' print(head(checkParams$paramMatrix))
#'
#' #Now actually run it
#' cl <- clusterMany(simData,dimReduce="PCA", nDimReduce=c(5,10,50),  isCount=FALSE,
#' clusterFunction="pam",ks=2:4,findBestK=c(TRUE,FALSE),removeSil=c(TRUE,FALSE))
#' print(cl)
#' head(colnames(clusterMatrix(cl)))
#'
#' #make names shorter for plotting
#' clNames <- clusterLabels(cl)
#' clNames <- gsub("TRUE", "T", clNames)
#' clNames <- gsub("FALSE", "F", clNames)
#' clNames <- gsub("k=NA,", "", clNames)
#'
#' par(mar=c(2, 10, 1, 1))
#' plotClusters(cl, axisLine=-2,clusterLabels=clNames)
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
# clSmaller <- clusterMany(simData, nDimReduce=c(5,10,50),  dimReduce="PCA",
# paramMatrix=checkParamsMat, subsampleArgs=checkParams$subsampleArgs,
# seqArgs=checkParams$seqArgs, mainClusterArgs=checkParams$mainClusterArgs)
#' @export 
setMethod(
  f = "clusterMany",
  signature = signature(x = "matrix"),
  definition = function(x,
      dimReduce="none",nDimReduce=NA, transFun=NULL,isCount=FALSE, ...
  ){
	####Basically, matrix version calls makeDimReduce and makeFilterStats and then sends it to the SingleCellFilter version.
	if(missing(dimReduce)) dimReduce<-"none"
	# if(anyNA(nDimReduce)){
# 		if(!"none" in dimReduce) dimReduce<-c(dimReduce,"none")
# 		nDimReduce<-na.omit(nDimReduce)
# 	}
	if(any(dim(x)==0)) stop("x must have non zero dimensions")
	dimReduce<-unique(dimReduce)
	doNone<-any(dimReduce=="none")
	#check can given dimReduce values match built in options.
	dimNam<-dimReduce[dimReduce %in% listBuiltInDimReduce()]
	filtNam<-dimReduce[dimReduce %in% listBuiltInFilterStats()]
	nValid<-length(c(dimNam,filtNam))
	if(doNone) nValid<-nValid+1
	if(!doNone & length(dimNam)==0 & length(filtNam)==0)
		  stop("dimReduce values given are not in built-in dimensionality reduction or built-in filters (and there is no such stored objects if a SingleCellFilter object). Option 'none' also not given, so nothing to do.")
	else if(length(dimReduce)!=nValid)
		warning("Some dimReduce values given are not in built in dimensionality reduction or built in filters (and there is no such stored objects if a SingleCellFilter object). Ignoring options.")
	if(length(dimNam)>0 | length(filtNam)>0){
		if(length(dimNam)>0){
			nDimReduce<-na.omit(nDimReduce)
			if(length(nDimReduce)==0) 
				stop("Must give nDimReduce values if choose a dimReduce option not equal to 'none'")
			maxDims<-max(nDimReduce)
			x<-makeDimReduce(x,dimReduce=dimNam,
				maxDims=maxDims,transFun=transFun,isCount=isCount)			
		}
		if(length(filtNam)>0){
			#Need to think how can pass options to filterData...
		  	x<-makeFilterStats(x,filterStat=filtNam, transFun=transFun,isCount=isCount)  	
		}
	}
	else{
		x<-SingleCellFilter(x)
	}
	return(clusterMany(x,dimReduce=dimReduce,nDimReduce=nDimReduce,transFun=transFun,isCount=isCount,...))
}
)

#' @rdname clusterMany
#' @export
setMethod(
  f = "clusterMany",
  signature = signature(x = "SingleCellFilter"),
  definition = function(x, ks=NA, clusterFunction, 
	  dimReduce="none",nFilter=NA,nDimReduce=NA,
	  alphas=0.1, findBestK=FALSE,
      sequential=FALSE, removeSil=FALSE, subsample=FALSE,
      silCutoff=0, distFunction=NA,
      betas=0.9, minSizes=1,
      transFun=NULL,isCount=FALSE,
      verbose=FALSE,
      mainClusterArgs=NULL,
      subsampleArgs=NULL,
      seqArgs=NULL,
      ncores=1, random.seed=NULL, run=TRUE,
      ...
  )
  {
	inputArgs<-as.list(environment()) #need so can pass all the args, not just the ...
	transFun<-.makeTransFun(transFun=transFun,isCount=isCount)  
    paramMatrix<-NULL
    if(!is.null(random.seed)){
        if(!is.null(subsampleArgs) && "ncores" %in% names(subsampleArgs)){
            if(subsampleArgs[["ncores"]]>1) stop("setting random.seed will not be reproducible if ncores given to subsampleArgs")
        }
    }
	anyFilterSaved<-!is.null(filterStats(x)) && any(dimReduce %in% filterNames(x))
	anyDimSaved<-length(reducedDims(x))>0 && any(dimReduce %in% reducedDimNames(x))
	anyFilter<-!is.null(filterStats(x))
	anyDim<-length(reducedDims(x))>0 
	anyDimBuiltIn<-any(dimReduce %in% listBuiltInDimReduce())
	anyFilterBuiltIn<-any(dimReduce %in% listBuiltInFilterStats())
	if(!all(dimReduce=="none") & !anyFilter & !anyFilterBuiltIn & !anyDim & !anyDimBuiltIn) 
		stop("'dimReduce' does not match any stored or builtin filtering statistics or dimensionality reduction")
	if(!all(dimReduce=="none") & ((!anyFilter & !anyDimSaved & anyFilterBuiltIn) || (!anyDim & !anyFilterSaved & anyDimBuiltIn)) ){
		###This will make it calculate the requested dimReduce values and then send it back to here as a SingleCellFilter object without the args of dimReduce, etc. ...
		## Note that if no Filter saved, and asked for filter
		outval<-do.call(clusterMany,c(list(x=assay(x)),inputArgs[!names(inputArgs)%in%"x"]))
		if(class(outval)=="ClusterExperiment") {
			#lost anything about the meta data, old filtering/dimReduce
			retval<-.addBackSEInfo(newObj=outval,oldObj=x)
			#but now have lost the newly calculated reducedDim etc.!
			if(anyFilterBuiltIn) filterStats(retval)<-filterStats(outval)
			if(anyDimBuiltIn) reducedDims(retval)<-reducedDims(outval)
			return(retval)
		}
		else return(outval)
	}
    else{
		###############
	    #Check inputs of dimReduce slots
		##NOTE: For now, IF there is a reducedDim slot, then will not try 
		##to patch in ones that are missing.
		##This means can list some that want to be calculated. 
		##Either do all of them ahead of time or let all of them be done 
		##during call to clusterMany...
		###############  
	  	doNone<-"none" %in% dimReduce
		if(doNone) dimReduce<-dimReduce[-grep("none",dimReduce)]
		if(length(dimReduce)>0){
			if(any(!dimReduce %in% c(reducedDimNames(x),filterNames(x)))){
				dimReduce<-dimReduce[dimReduce %in%c(reducedDimNames(x),filterNames(x))]
				if(length(dimReduce)>0) warning("Not all of dimReduce value match a reducedDimNames or filterNames of the 'SingleCellFilter' object. Will ignore them:",paste(dimReduce[!dimReduce %in%c(reducedDimNames(x),filterNames(x))],collapse=","))
				else stop("No dimReduce value was given that matches stored reducedDimNames or filterNames of the object.")
			}
			#check if nPCA values
			if(any(dimReduce %in% reducedDimNames(x))){
				maxDimValues<-sapply(reducedDims(x)[dimReduce[dimReduce %in%reducedDimNames(x)]],ncol)
				if(length(na.omit(nDimReduce))>0 && all(na.omit(nDimReduce) > max(maxDimValues))) 
					stop("The values of nDimReduce given are all higher than the maximum components stored in the reducedDims slot of the input object. Run 'makeDimReduce' to get larger number of components.")
				
			}
			
			#check if give nFilter if filterNames and no NA values
			if(any(dimReduce %in% filterNames(x))){
				if(!missing(nFilter) && any(is.na(nFilter))){
					warning("NA values have no meaning for the argument nFilter and will be ignored")		
					nFilter<-na.omit(nFilter)
				}
				if(missing(nFilter) || length(nFilter)==0){
					stop("no valid nFilter values given, but dimReduce values given indicate a filterStat to be used.")
				}
			}
		}
		else{
			nDimReduce<-NA
			nFilter<-NA
			maxDimValues<-NA #indicates that only "none" will be done
		}

		###############
	    #Start creating the combinations
		###############
	    if(is.null(paramMatrix)){
		  if(doNone) dimReduce<-c(dimReduce,"none")
	      param <- expand.grid(dimReduce=dimReduce,
		    nDimReduce=nDimReduce, nFilter=nFilter,k=ks, alpha=alphas, findBestK=findBestK, 
			beta=betas, minSize=minSizes,
	        sequential=sequential, distFunction=distFunction,
	        removeSil=removeSil, subsample=subsample,
	        clusterFunction=clusterFunction,silCutoff=silCutoff)
			
		  ###########
	      #Check param matrix:
	      #don't vary them across ones that don't matter (i.e. 0-1 versus K);
	      #code sets to single value and then will do unique
	      #also deals with just in case the user gave duplicated values of something by mistake.
	      ###########
		  paramAlgTypes<-algorithmType(param[,"clusterFunction"])
		  if(length(paramAlgTypes)!=nrow(param)) stop("Internal coding error in clusterMany: not getting right number of type of algorithms from param")
		  #---
		  #type K fixes
		  #---
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
		  #---
		  # type01 combinations
		  #---
	      type01 <- which( paramAlgTypes=="01")
	      if(length(type01)>0){
	        param[type01,"findBestK"] <- FALSE
	        param[type01,"removeSil"] <- FALSE
	        param[type01,"silCutoff"] <- 0
	      }
		  #---
		  #Turn off distFunction for some combinations
		  #---
		  #those that subsample, because will distance that of co-occurance
	      whSubsample<-which(param[,"subsample"])
	      if(length(whSubsample)>0){
	        param[whSubsample,"distFunction"]<-NA
	      }
		  #those that subsample, because will distance that of co-occurance
	      whDimReduce<-which(param[,"dimReduce"]!="none")
	      if(length(whDimReduce)>0){
	        param[whDimReduce,"distFunction"]<-NA
	      }

		  #---
		  #Check value alpha, beta values
		  #---
	      alpha01 <- which(param[,"alpha"]<0 | param[,"alpha"]>1)
		  if(length(alpha01)>0){
			  stop("alpha value must be in (0,1)")
			  param[alpha01,"alpha"]<-NA
		  }
	      beta01 <- which(param[,"beta"]<0 | param[,"beta"]>1)
		  if(length(beta01)>0){
			  stop("beta value must be in (0,1)")
			  param[beta01,"beta"]<-NA
		  }

  		  #---
		  # deal with nDimReduce NA or larger than the size of the dataset
		  # set it to the maximum value possible.
		  #---
		  whDimReduce<-which(param[,"dimReduce"] %in% reducedDimNames(x))
		  if(length(whDimReduce)>0 && length(na.omit(maxDimValues[whDimReduce]))>0){
			  #if NA, means do the largest possible dimension saved for that method
			  whNADim<-intersect(which(is.na(param[,"nDimReduce"])),whDimReduce)
			  maxDimValues<-maxDimValues[param[whNADim,"dimReduce"]]
			  if(length(whNADim)>0){
				  param[whNADim,"nDimReduce"]<-maxDimValues
			  }
			  if(anyNA(param[whDimReduce,"nDimReduce"])) stop("Internal coding error: didn't get rid of NA dimReduce in checks")
			  whAbove<-intersect(which(param[,"nDimReduce"]>maxDimValues),whDimReduce)
			  if(length(whAbove)>0){
				  param[whAbove,"nDimReduce"]<-maxDimValues[whAbove]
			  }
		  }
		  #now turn to NA is when dimReduce a dim reduce
		  whOther<-which(!param[,"dimReduce"]%in% reducedDimNames(x))
		  if(length(whOther)>0){
			  param[whOther,"nDimReduce"]<-NA
		  }
		  
  		  #---
		  # deal with nFilter NA or larger than the size of the dataset
		  # set it to the maximum value possible.
		  #---
		  whFilter<-which(param[,"dimReduce"] %in% filterNames(x))
		  whTooLarge<-intersect(which(param[,"nFilter"]>NROW(x)),whFilter)
		  if(length(whTooLarge)>0){
			  param[whTooLarge,"dimReduce"]<-"none"
		  }
		  
		  #now turn to NA is when dimReduce a dim reduce
		  whOther<-which(!param[,"dimReduce"]%in% filterNames(x))
		  if(length(whOther)>0){
			  param[whOther,"nFilter"]<-NA
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
			  clFun<-param[,"clusterFunction"]
			  if(any(algorithmType(clFun)=="K")) stop("One of clusterFunctions chosen requires choice of k")
	          else param<-param[-whInvalid,]
			
			}
		  
		  
	      if(any(!is.na(param[,"nFilter"]) & !is.na(param[,"nDimReduce"]))) 
			  stop("Internal error: failed to properly remove inconsistent nFilter, nDimReduce combination.")
	      if(any(is.na(param[,"nFilter"]) & is.na(param[,"nDimReduce"] & !param[,"dimReduce"] %in% "none"))) stop("Internal error: NA in both nFilter, nDimReduce combination without equal to 'none'")
	      #####
	      #require at least 2 combinations:
	      #####
	      if(nrow(param)<=1) {
	        stop("set of parameters imply only 1 combination. If you wish to run a single clustering, use 'clusterSingle'")
	      }

	      #####
	      #give names to the parameter combinations.
	      #####
	      whVary <- which(apply(param,2,function(x){length(unique(x))>1}))
	      if(length(whVary)>0) {
			  #for some reason, this code started changing TRUE/FALSE in to 1/0
	        # cnames<-apply(param[,whVary,drop=FALSE],1,function(x){
	        # paste(colnames(param)[whVary],as.character(x),sep="=",collapse=",")})	
			cnames<-sapply(1:nrow(param),function(ii){
				paste(colnames(param)[whVary],as.character(param[ii,whVary]),sep="=",collapse=",")
			})
	      } else {
	        stop("set of parameters imply only 1 combination. If you wish to run a single clustering, use 'clusterSingle'")
	      }
	      cnames <- gsub("dataset=","",cnames)
	      cnames <- gsub("= ","=",cnames)
	      cnames[param[,"sequential"]] <- gsub("k=", "k0=", cnames[param[,"sequential"]])
		  #should I combine together nDimReduce and nFilter like they were before for the labels?
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
	      cat(nrow(param),"parameter combinations,",sum(param[,"sequential"]),"use sequential method,",sum(param[,"subsample"]),"use subsampling method\n")
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
	      dimReduce<-as.character(par[["dimReduce"]])
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
		  if(dimReduce=="none") 
			  dat<-transformData(x,transFun=transFun) 
		  else if(dimReduce %in% reducedDimNames(x)) 
			  dat<-t(reducedDim(x,dimReduce)[,1:par[["nDimReduce"]]] ) 
		  else if(dimReduce %in% filterNames(x)) 
			  dat<-transformData( filterData(x, type=dimReduce, percentile=par[["nFilter"]]),
		  				transFun=transFun)
		  else stop("Internal error: dimReduce value that not in filterNames or reducedDimNames")
		  #(Note, computational inefficiency: means reordering each time, even if same filter. But not recalculating filter.)
		  if(!is.null(distFunction)){
  			#need to update here when have filter (see below)
	        diss<- allDist[[distFunction]]
	        clusterSingle(x=dat, diss=diss,subsample=subsample, dimReduce="none",
	                      mainClusterArgs=mainClusterArgs,
	                      subsampleArgs=subsampleArgs, seqArgs=seqArgs,
	                      sequential=sequential, transFun=function(x){x},checkDiss=FALSE)       }
	      else
		  clusterSingle(x=dat, subsample=subsample,
	                 mainClusterArgs=mainClusterArgs, dimReduce="none",
	                 subsampleArgs=subsampleArgs, seqArgs=seqArgs,
	                 sequential=sequential, transFun=function(x){x},checkDiss=FALSE) 
		    }
	    if(run){
		##Calculate distances necessary only once
	      if(any(!is.na(param[,"distFunction"]))){
	        distParam<-unique(param[,c("dataset","distFunction")])
	        distParam<-distParam[!is.na(distParam[,"distFunction"]),]
		    ##Assume only take distances on original data (or filtered version of it)
        	#need to update here when have filter
			dat<-assay(x)
			allDist<-lapply(1:nrow(distParam),function(ii){
			distFun<-as.character(distParam[ii,"distFunction"])
				#be conservative and check for the 01 type if any of clusterFunctions are 01.
				algCheckType<-if(any(paramAlgTypes=="01")) "01" else "K" 
				distMat<-.makeDiss(dat, distFunction=distFun, checkDiss=TRUE, algType=algCheckType)
				return(distMat)
			})
			#need to update here when have filter
			##paste(distParam[,"dataset"],distParam[,"distFunction"],sep="--")
	        names(allDist)<-distParam[,"distFunction"]			

	      }

	      if(verbose) {
	        cat("Running Clustering on Parameter Combinations...")
	      }
	  
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
		  outval <- ClusterExperiment(x, clusters=clMat,
			                                    transformation=transFun,
			                                    clusterInfo=clInfo,			                                    clusterTypes="clusterMany",checkTransformAndAssay=FALSE)
	    } else{
	      if(verbose) {
	        cat("Returning Parameter Combinations without running them (to run them choose run=TRUE)\n")
	      }
	      outval<-list(paramMatrix=param, mainClusterArgs=mainClusterArgs, seqArgs=seqArgs,subsampleArgs=subsampleArgs)
	    }
		return(outval)
	}

  }
)

#' @rdname clusterMany
#' @export
setMethod(
  f = "clusterMany",
  signature = signature(x = "ClusterExperiment"),
  definition = function(x, dimReduce="none", nFilter=NA, nDimReduce=NA,
                        eraseOld=FALSE, ...)
  {
  	if(any(c("transFun","isCount") %in% names(list(...)))) 
  		stop("The internally saved transformation function of a ClusterExperiment object must be used when given as input and setting 'transFun' or 'isCount' for a 'ClusterExperiment' is not allowed.")  
    outval<-clusterMany(as(x,"SingleCellFilter"), dimReduce=dimReduce, nFilter=nFilter,
                        nDimReduce=nDimReduce, transFun=transformation(x), ...)
    if(class(outval)=="ClusterExperiment") {
		
	  #outval<-.addBackSEInfo(newObj=outval,oldObj=x) #added to '.addNewResult'
      ##Check if clusterMany already ran previously
      x<-.updateCurrentWorkflow(x,eraseOld,"clusterMany")

      if(!is.null(x)){
		  retval<-.addNewResult(newObj=outval,oldObj=x) #make decisions about what to keep.
		  
	  }
      else{
		  retval<-.addBackSEInfo(newObj=outval,oldObj=x)
	  }
	  filterStats(retval)<-filterStats(outval)
	  reducedDims(retval)<-reducedDims(outval)	  
		  #both above check validity.
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
  definition = function(x, ...){
	  clusterMany(as(x,"SingleCellExperiment"),...)
  }
  )

#' @rdname clusterMany
#' @export
setMethod(
f = "clusterMany",
signature = signature(x = "SingleCellExperiment"),
definition = function(x, ...){
  clusterMany(as(x,"SingleCellFilter"),...)
}
)
  

#' @export
#' @rdname clusterMany
setMethod(
f = "clusterMany",
signature = signature(x = "data.frame"),
definition = function(x,...){clusterMany(data.matrix(x),...)}
)



