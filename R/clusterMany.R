#' Create a matrix of clustering across values of parameters
#'
#' Given a range of parameters, this function will return a matrix with the
#' clustering of the samples across the range, which can be passed to
#' \code{plotClusters} for visualization.
#'
#' @name clusterMany
#' @aliases clusterMany,matrixOrHDF5-method
#' @param x the data matrix on which to run the clustering. Can be object of the
#'   following classes: matrix (with genes in rows),
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}},
#'   \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   or \code{ClusterExperiment}.
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
#' @param clusterFunction function used for the clustering. This must be either 
#'   1) a character vector of built-in clustering techniques, or 2) a
#'   \emph{named} list of \code{\link{ClusterFunction}} objects. Current
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
#' @param nFilterDims vector of the number of the most variable features to keep
#'   (when "var", "abscv", or "mad" is identified in \code{reduceMethod}).
#' @param nReducedDims vector of the number of dimensions to use (when
#'	\code{reduceMethod} gives a dimensionality reduction method).
#' @param eraseOld logical. Only relevant if input \code{x} is of class
#'   \code{ClusterExperiment}. If TRUE, will erase existing workflow results
#'   (clusterMany as well as mergeClusters and makeConsensus). If FALSE, existing
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
#' @param ... For signature \code{matrix}, arguments to be passed on to mclapply
#'   (if ncores>1). For all the other signatures, arguments to be passed to the
#'   method for signature \code{matrix}.
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
#' @return If \code{run=TRUE} will
#'   return a \code{ClusterExperiment} object, where the results are stored as
#'   clusterings with clusterTypes \code{clusterMany}. Depending on
#'   \code{eraseOld} argument above, this will either delete existing such
#'   objects, or change the clusterTypes of existing objects. See argument
#'   \code{eraseOld} above. Arbitrarily the first clustering is set as the
#'   primaryClusteringIndex.
#'
#' @return If \code{run=FALSE} a list with elements:
#' \itemize{
#'   \item{\code{paramMatrix}}{ a matrix giving the parameters of each
#'   clustering, where each column is a possible parameter set by the user and
#'   passed to \code{\link{clusterSingle}} and each row of paramMatrix
#'   corresponds to a clustering in \code{clMat}} 
#'   \item{\code{mainClusterArgs}}{
#'   a list of (possibly modified) arguments to mainClusterArgs}
#'   \item{\code{seqArgs=seqArgs}}{a list of (possibly modified) arguments to
#'   seqArgs} 
#'   \item{\code{subsampleArgs}}{a list of (possibly modified)
#'   arguments to subsampleArgs} }
#'
#' @examples
#' data(simData)
#'
#' #Example: clustering using pam with different dimensions of pca and different
#' #k and whether remove negative silhouette values
#' #check how many and what runs user choices will imply:
#' checkParams <- clusterMany(simData,reduceMethod="PCA", makeMissingDiss=TRUE,
#'    nReducedDims=c(5,10,50), clusterFunction="pam", isCount=FALSE,
#'    ks=2:4,findBestK=c(TRUE,FALSE),removeSil=c(TRUE,FALSE),run=FALSE)
#' print(head(checkParams$paramMatrix))
#'
#' #Now actually run it
#' cl <- clusterMany(simData,reduceMethod="PCA", nReducedDims=c(5,10,50),  isCount=FALSE,
#'    clusterFunction="pam",ks=2:4,findBestK=c(TRUE,FALSE),makeMissingDiss=TRUE, 
#'    removeSil=c(TRUE,FALSE))
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
#'	    alphas=c(0.1,0.2,0.3), findBestK=c(TRUE,FALSE), sequential=c(FALSE),
#'	    subsample=c(FALSE), removeSil=c(TRUE), clusterFunction="pam", 
#'	    makeMissingDiss=TRUE,
#'      mainClusterArgs=list(minSize=5, kRange=2:15), ncores=1, random.seed=48120))
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
# clSmaller <- clusterMany(simData, nReducedDims=c(5,10,50),  reduceMethod="PCA",
# paramMatrix=checkParamsMat, subsampleArgs=checkParams$subsampleArgs,
# seqArgs=checkParams$seqArgs, mainClusterArgs=checkParams$mainClusterArgs)
#' @export
setMethod(
    f = "clusterMany",
    signature = signature(x = "matrixOrHDF5"),
    definition = function(x,
                          reduceMethod="none",nReducedDims=NA, transFun=NULL,isCount=FALSE, ...
    ){
        ####Basically, matrix version calls makeReducedDims and makeFilterStats and then sends it to the SingleCellExperiment version.
        if(missing(reduceMethod)) reduceMethod<-"none"
        # if(anyNA(nReducedDims)){
        # 		if(!"none" in reduceMethod) reduceMethod<-c(reduceMethod,"none")
        # 		nReducedDims<-na.omit(nReducedDims)
        # 	}
        if(any(dim(x)==0)) stop("x must have non zero dimensions")
        reduceMethod<-unique(reduceMethod)
        doNone<-any(reduceMethod=="none")
        #check can given reduceMethod values match built in options.
        dimNam<-reduceMethod[isBuiltInReducedDims(reduceMethod)]
        filtNam<-reduceMethod[isBuiltInFilterStats(reduceMethod)]
        nValid<-length(c(dimNam,filtNam))
        if(doNone) nValid<-nValid+1
        if(!doNone & length(dimNam)==0 & length(filtNam)==0)
            stop("reduceMethod values given are not in built-in dimensionality reduction or built-in filters (and there is no such stored objects if a SingleCellExperiment object). Option 'none' also not given, so nothing to do.")
        else if(length(reduceMethod)!=nValid)
            warning("Some reduceMethod values given are not in built in dimensionality reduction or built in filters (and there is no such stored objects if a SingleCellExperiment object). Ignoring options.")
        if(length(dimNam)>0 | length(filtNam)>0){
            if(length(dimNam)>0){
                nReducedDims<-na.omit(nReducedDims)
                if(length(nReducedDims)==0)
                    stop("Must give nReducedDims values if choose a reduceMethod option not equal to 'none' and not in stored reducedDims slot.")
                maxDims<-max(nReducedDims)
                x<-makeReducedDims(x,reducedDims=dimNam,
                                   maxDims=maxDims,transFun=transFun,isCount=isCount)
            }
            if(length(filtNam)>0){
                #Need to think how can pass options to filterData...
                x<-makeFilterStats(x,filterStat=filtNam, transFun=transFun,isCount=isCount)
            }
        }
        else{
            x<-SingleCellExperiment(x)
        }
        return(clusterMany(x,reduceMethod=reduceMethod,nReducedDims=nReducedDims,transFun=transFun,isCount=isCount,...))
    }
)

#' @rdname clusterMany
#' @param whichAssay numeric or character specifying which assay to use. See
#'   \code{\link[SummarizedExperiment]{assay}} for details.
#' @details The given \code{reduceMethod} values must either be \emph{all}
#'   precalculated filtering/dimensionality reduction stored in the appropriate
#'   location, or must \emph{all} be character values giving a built-in
#'   filtering/dimensionality reduction methods to be calculated. If some of the
#'   filtering/dimensionality methods are already calculated and stored, but not
#'   all, then they will \emph{all} be recalculated (and if they are not all
#'   built-in methods, this will give an error). So to save computational time
#'   with pre-calculated dimensionality reduction, the user must make sure they
#'   are \emph{all} precalculated. Also, user-defined values (i.e. not built-in
#'   functions) cannot be mixed with built-in functions unless they have already
#'   been precalculated (see \code{\link{makeFilterStats}} or
#'   \code{\link{makeReducedDims}}).
#' @export
setMethod(
    f = "clusterMany",
    signature = signature(x = "SingleCellExperiment"),
    definition = function(x, ks=NA, clusterFunction,
                          reduceMethod="none",
                          nFilterDims= defaultNDims(x, 
                              reduceMethod, type="filterStats"),
                          nReducedDims= defaultNDims(x, 
                              reduceMethod, type="reducedDims"),
                          alphas=0.1, findBestK=FALSE,
                          sequential=FALSE, removeSil=FALSE, subsample=FALSE,
                          silCutoff=0, distFunction=NA,
                          betas=0.9, minSizes=1,
                          transFun=NULL,isCount=FALSE,
                          verbose=TRUE,
                          parameterWarnings=FALSE,
                          mainClusterArgs=NULL,
                          subsampleArgs=NULL,
                          seqArgs=NULL,
                          whichAssay=1,
                          makeMissingDiss=FALSE,
                          ncores=1, random.seed=NULL, run=TRUE,
                          ...
    )
    {
        
        #need so can pass all the args, not just the ...
        inputArgs<-as.list(environment()) 
        transFun<-.makeTransFun(transFun=transFun,isCount=isCount)
        paramMatrix<-NULL
        if(!is.null(random.seed)){
            if(!is.null(subsampleArgs) && "ncores" %in% names(subsampleArgs)){
                if(subsampleArgs[["ncores"]]>1) stop("setting random.seed will not be reproducible if ncores given to subsampleArgs")
            }
        }
        
        #######################
        ### Deal with reduceMethods
        #######################
        #issue: have to send reduceMethod, but don't know which are which type
        isExisting<-isReducedDims(x,reduceMethod) | isFilterStats(x,reduceMethod)
        isBuiltIn<- isBuiltInReducedDims(reduceMethod) | isBuiltInFilterStats(reduceMethod)
        isNone<-reduceMethod=="none"
        if(any(isNone)){
            isExisting[isNone]<-TRUE
            isBuiltIn[isNone]<-TRUE
        }
        isBuiltInNotExisting<-isBuiltIn & !isExisting
        
        # anyFilter<-anyValidFilterStats(x)
        # anyDim<-anyValidReducedDims(x)
        # anyFilterSaved<-anyFilter && any(isFilterStats(x,reduceMethod))
        # anyDimSaved<-anyDim && any(isReducedDims(x,reduceMethod))
        # anyDimBuiltIn<-any(isBuiltInReducedDims(reduceMethod))
        # anyFilterBuiltIn<-any(isBuiltInFilterStats(reduceMethod))
        if(!all(isNone | isBuiltIn | isExisting))
            stop("Some values of 'reduceMethod' do not match any stored or built-in filtering statistics or dimensionality reduction")
        if(!all(isNone | isBuiltIn) & !all(isNone | isExisting))
            stop("All values of 'reduceMethod' need to either match an existing (i.e. stored) filtering/dimensionality reduction or they need to all match a built-in function to be calculated")
        
        
        if(all(isBuiltIn) & any(isBuiltInNotExisting) ){
            
            .mynote(paste0("Not all of the methods requested in 'reduceMethod' have been calculated. Will calculate all the methods requested (any pre-existing values -- filtering statistics or dimensionality reductions -- with these names will be recalculated and overwritten): ",paste(reduceMethod,collapse=","),"."))
            
            reduceMethod<-unique(reduceMethod)
            doNone<-any(reduceMethod=="none")
            
            #check can given reduceMethod values match built in options.
            dimNam<-reduceMethod[isBuiltInReducedDims(reduceMethod)]
            filtNam<-reduceMethod[isBuiltInFilterStats(reduceMethod)]
            nValid<-length(c(dimNam,filtNam))
            
            if(doNone) nValid<-nValid+1
            
            if(!doNone & length(dimNam)==0 & length(filtNam)==0)
                stop("reduceMethod values given are not in built-in dimensionality reduction or built-in filters (and there is no such stored objects if a SingleCellExperiment object). Option 'none' also not given, so nothing to do.")
            else if(length(reduceMethod)!=nValid)
                warning("Some reduceMethod values given are not in built in dimensionality reduction or built in filters (and there is no such stored objects if a SingleCellExperiment object). Ignoring options.")
            
            if(length(dimNam)>0 | length(filtNam)>0){
                if(length(dimNam)>0){
                    nReducedDims<-na.omit(nReducedDims)
                    if(length(nReducedDims)==0)
                        stop("Must give nReducedDims values if choose a reduceMethod option not equal to 'none' and not in stored reducedDims slot.")
                    maxDims<-max(nReducedDims)
                    x<-makeReducedDims(x,reducedDims=dimNam, 
                        whichAssay = whichAssay,
                        maxDims=maxDims,transFun=transFun,
                        isCount=isCount)
                }
                if(length(filtNam)>0){
                    #Need to think how can pass options to filterData...
                    x<-makeFilterStats(x,filterStat=filtNam, 
                        whichAssay = whichAssay,
                        transFun=transFun,isCount=isCount)
                }
            }
        }
        
        #----------------
        # Should come here ONLY if all methods needed have been calculated (i.e. built-in function)
        #Check inputs of reduceMethod slots
        ##NOTE: For now, IF there is a reducedDim slot, then will not try
        ##to patch in ones that are missing.
        ##This means can list some that want to be calculated.
        ##Either do all of them ahead of time or let all of them be done
        ##during call to clusterMany...
        #----------------
        doNone<-any(isNone)
        if(doNone) reduceMethod<-reduceMethod[-grep("none",reduceMethod)]
        if(length(reduceMethod)>0){
            if(any(!isReducedDims(x,reduceMethod) & !isFilterStats(x,reduceMethod))){
                stop("Internal Coding Error -- shouldn't have gotten to this point without catching from earlier error check that all necessary statistics are calculated.")
            }
            #check if nReducedDims values
            if(any(isReducedDims(x,reduceMethod))){
                maxDimValues<- 
                    ncolReducedDims(x)[isReducedDims(x,reduceMethod)]
                if(length(na.omit(nReducedDims))>0 && all(na.omit(nReducedDims) > max(maxDimValues)))
                    stop("The values of nReducedDims given are all higher than the maximum components stored in the reducedDims slot of the input object. Run 'makeReducedDims' to get larger number of components.")
                
            }
            
            #check if give nFilterDims if isFilterStats and no NA values
            if(any(isFilterStats(x,reduceMethod))){
                if(!missing(nFilterDims) && any(is.na(nFilterDims))){
                    warning("NA values have no meaning for the argument nFilterDims and will be ignored")
                    nFilterDims<-na.omit(nFilterDims)
                }
                if(length(nFilterDims)==0){
                    stop("no valid nFilterDims values given, but reduceMethod values given indicate a filterStat to be used.")
                }
            }
        }
        else{
            nReducedDims<-NA
            nFilterDims<-NA
            maxDimValues<-NA #indicates that only "none" will be done
        }
        #----------FINISH reduce methods section------
        
        ###############
        # Start creating the parameter combinations
        ###############
        if(is.null(paramMatrix)){
            if(is.list(clusterFunction)){
                .checkFunctionList(clusterFunction)
                clusterFunctionList<-clusterFunction
                clusterFunctionNames<-names(clusterFunction)
            }
            else{
                clusterFunctionList<-getBuiltInFunction(clusterFunction)
                if(length(clusterFunction)==1 & is(clusterFunctionList,"ClusterFunction")){
                    clusterFunctionList<-list(clusterFunctionList)
                    names(clusterFunctionList)<-clusterFunction
                }
                clusterFunctionNames<-clusterFunction
            }
            if(doNone) reduceMethod<-c(reduceMethod,"none")
            param <- expand.grid(clusterFunction=clusterFunctionNames,
                k=ks, 
                alpha=alphas, 
                beta=betas, 
                sequential=sequential, 
                subsample=subsample,
                reduceMethod=reduceMethod,
                nReducedDims=nReducedDims, 
                nFilterDims=nFilterDims,
                minSize=minSizes,
                findBestK=findBestK,
                removeSil=removeSil, 
                silCutoff=silCutoff,
                distFunction=distFunction
                )
            
            if(nrow(param)<=1) {
                stop("set of parameters imply only 1 combination. If you wish to run a single clustering, use 'clusterSingle'")
            }
            #-------------
            # Check param matrix:
            # don't vary them across ones that don't matter (i.e. 0-1 versus K);
            #   -> code sets to single value and then will do unique to remove 
            # also deals with just in case the user gave duplicated 
            #   values of something by mistake.
            #-------------
            
            # small function to pull the list of all functions 
            # from current param matrix            
            cf<-function(param){
                clusterFunctionList[param[,"clusterFunction"]]
            }
            paramAlgTypes<-algorithmType(cf(param))
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
                    ##Might could handle this better by call to .checkArgs
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
            #those that use reducedDims will not use dist
            #but those that filter could use different distances...
						# FIXME: should there be some kind of warning if wind up completely ignoring 'distFunction' that user gave? Ditto on other parameters.
            whDimReduce<-which(param[,"reduceMethod"]!="none" & isReducedDims(x,param[,"reduceMethod"]) )
            if(length(whDimReduce)>0){
                param[whDimReduce,"distFunction"]<-NA
            }
            
            #---
            # deal with nReducedDims NA or larger than the size of the dataset
            # set it to the maximum value possible.
            #---
            whDimReduce<-which(isReducedDims(x,param[,"reduceMethod"]))
            if(length(whDimReduce)>0 && length(na.omit(maxDimValues[whDimReduce]))>0){
                #if NA, means do the largest possible dimension saved for that method
                whNADim<-intersect(which(is.na(param[,"nReducedDims"])),whDimReduce)
                maxDimValuesNA<-maxDimValues[param[whNADim,"reduceMethod"]]
                if(length(whNADim)>0){
                    param[whNADim,"nReducedDims"]<-maxDimValuesNA
                }
                if(anyNA(param[whDimReduce,"nReducedDims"])) stop("Internal coding error: didn't get rid of NA reduceMethod in checks")
                whAbove<-intersect(which(param[,"nReducedDims"] > maxDimValues[param[,"reduceMethod"]]), whDimReduce)
                maxDimValuesAbove<-maxDimValues[param[whAbove,"reduceMethod"]]
                if(length(whAbove)>0){
                    param[whAbove,"nReducedDims"]<-maxDimValuesAbove
                }
            }
            #now turn to NA is when reduceMethod a dim reduce
            whOther<-which(!isReducedDims(x,param[,"reduceMethod"]))
            if(length(whOther)>0){
                param[whOther,"nReducedDims"]<-NA
            }
            
            #---
            # deal with nFilterDims NA or larger than the size of the dataset
            # set it to the maximum value possible.
            #---
            whFilter<-which(isFilterStats(x,param[,"reduceMethod"]))
            whTooLarge<-intersect(which(param[,"nFilterDims"]>NROW(x)),whFilter)
            if(length(whTooLarge)>0){
                param[whTooLarge,"reduceMethod"]<-"none"
            }
            
            #now turn to NA is when reduceMethod a dim reduce
            whOther<-which(!isFilterStats(x,param[,"reduceMethod"]))
            if(length(whOther)>0){
                param[whOther,"nFilterDims"]<-NA
            }
           
            # get rid of duplicates
            param <- unique(param)
            
            #####################
            # Deal with those that are invalid combinations:
            # Also, if ever reinstate param option, then should apply these checks to that param
            #####################
            if(is.null(mainClusterArgs)) 
                mainClusterArgs<-list(clusterArgs=list())
            if(is.null(subsampleArgs)) 
                subsampleArgs<-list(clusterArgs=list())
            paramCheck<-function(paramRow, returnValue, warn){
                totalArgs<- .makeArgsFromParam(paramRow,
                    mainClusterArgs=mainClusterArgs,
                    seqArgs=seqArgs,
                    subsampleArgs=subsampleArgs,
                    clusterFunctionList=clusterFunctionList)
                checkOut<-.checkArgs(inputType="X",
                                 main=TRUE,
                                 subsample=totalArgs$subsample,
                                 sequential=totalArgs$sequential,
                                 mainClusterArgs=totalArgs$mainClusterArgs,
                                 subsampleArgs=totalArgs$subsampleArgs,
                                 seqArgs=totalArgs$seqArgs,
                                 allowMakeDiss=TRUE,
                                 warn = warn) #most of these are because have extra parameters 
                if(returnValue=="logical"){
                    if(is.character(checkOut)) return(FALSE)
                    else return(TRUE)                    
                }
                else{
                    return(checkOut)
                }
            
            }
            #-----
            # Do first run to get rid of invalid selections (i.e. silently)
            # Note! Cannot run apply on param if don't want all character values.
            # e.g. checks<-apply(param,1,paramCheck,returnValue="logical")
            #-----
            checks<-sapply(1:nrow(param), 
                function(i){paramCheck(param[i,],returnValue="logical",
                    warn=FALSE)})
            whInvalid<-which(!checks)
            if(length(whInvalid)>0) {
                if(length(whInvalid)==nrow(param)){
                    checksFull<-lapply(1:nrow(param), 
                        function(i){paramCheck(param[i,],returnValue="full",warn=parameterWarnings)})
                    stop(sprintf("Set of parameters imply %s combinations, none of which are valid (note that if clustering matrix requires calculation of distance matrix, user must now set `makeMissingDiss=TRUE`).  Errors given are:\n\n %s",nrow(param),paste(checksFull[whInvalid],collapse="\n")))
                    
                }
                param <- param[-whInvalid, ,drop=FALSE]
            }
            #-----
            # Do second run to get full warnings that remain
            # & deal with creating distances
            #-----
            checks<-lapply(1:nrow(param), 
                function(i){paramCheck(param[i,],returnValue="full",
                    warn=parameterWarnings)})
            # doDiss means will send diss as inputMatrix
            # doDissPost means will pass diss to mainClusterArgs        
            doDiss<-sapply(checks,function(x){x$doDiss}) | !is.na(distFunction)
            doDissPost<-sapply(checks,function(x){x$doDissPost})
            #in case user asked for distFunction, don't do again:
            doDissPost[doDiss]<-FALSE 
            param$passDistToInput <- doDiss
            param$passDistToMain <- doDissPost
            # missDiss is all places where will need calculate distance
            missDiss<-doDiss | doDissPost
            if(any(missDiss) & !makeMissingDiss){
                stop("Parameter combinations requested require calculation of at least one distance matrix. User must now set `makeMissingDiss=TRUE` to have clusterMany calculated the necessary distances.")
            }
            if(any(missDiss) & makeMissingDiss){
                # give "default" to those not assigned
                param[is.na(distFunction) & missDiss,"distFunction"]<-"default"
                
                ## expand param matrix to include additional information for calculating/passing distances
                subCF<-sapply(checks,function(x){
                    cf<-x$subsampleArgs[["clusterFunction"]]
                    if(!is.null(cf)) algorithmType(cf) else NA})
                mainCF<-sapply(checks,function(x){
                    algorithmType(x$mainClusterArgs[["clusterFunction"]])})                
                # Need clustering algorithm for calculating distances
                param$distAlgType<-NA
                param$distAlgType[missDiss]<-
                    ifelse(param[missDiss,"subsample"],subCF,mainCF)
                param$distAlgType[!doDiss & doDissPost]<-
                    mainCF[!doDiss & doDissPost]
            }
            
            if(any(!is.na(param[,"nFilterDims"]) &
             !is.na(param[,"nReducedDims"]))){
                stop("Internal error: failed to properly remove inconsistent nFilterDims, nReducedDims combination.")
                
            }
            if(any(is.na(param[,"nFilterDims"]) &
                is.na(param[,"nReducedDims"] & 
                !param[,"reduceMethod"] %in% "none"))) stop("Internal error: NA in both nFilterDims, nReducedDims combination without equal to 'none'")
            #####
            #require at least 2 combinations:
            #####
            if(nrow(param)<=1) {
                stop("set of parameters imply only 1 combination, after removing duplicates and invalid combinations. If you wish to run a single clustering, use 'clusterSingle'")
            }
            
            #####
            #give names to the parameter combinations.
            #####
            charParam<-as.matrix(param)
            ## Make default distance more interpretable
            whDefault<-which(charParam[,"distFunction"]=="default")
            if(length(whDefault)>0){
                charParam[whDefault,"distFunction"]<-
                    paste(charParam[whDefault,"distFunction"],
                        charParam[whDefault,"distAlgType"],sep="")           
            }
            ## Remove those that I added
            myAdditions<-c("distAlgType","passDistToMain","passDistToInput")
            wh<-which(names(param) %in% myAdditions)
            if(length(wh)>0) charParam<-charParam[,-wh]
            whVary <- which(apply(charParam,2,function(x){length(unique(x))>1}))
            if(length(whVary)>0) {
                makeLabel<- function(ii){
                    paste(colnames(charParam)[whVary],charParam[ii,whVary],sep="=",collapse=",")
                }
                cnames<-sapply(seq_len(nrow(charParam)),makeLabel)
            } else {
                stop("set of parameters imply only 1 combination. If you wish to run a single clustering, use 'clusterSingle'")
            }
            cnames <- gsub("dataset=","",cnames)
            cnames <- gsub("= ","=",cnames)
            cnames[param[,"sequential"]] <- gsub("k=", "k0=", cnames[param[,"sequential"]])
            #should I combine together nReducedDims and nFilterDims like they were before for the labels?
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
           
        
        ## Uniform way to create and access names of the list of the list of the distances
        createDistNames<-function(paramMatrix,returnNames=TRUE){
            paramNames<-c("reduceMethod", "nFilterDims",
                             "distFunction")
            if(!returnNames) paramNames<-c(paramNames,"distAlgType")
            if(any(! paramNames %in% colnames(paramMatrix)))
                 stop("Internal error: must have names", paste(paramNames,collapse=","))
            distParam<-unique(paramMatrix[, paramNames,drop=FALSE])
            distParam<-distParam[!is.na(distParam[,"distFunction"]), ,drop=FALSE]
            if(returnNames){
                uniqueId<-apply(distParam,1,paste,collapse=",")
                return(uniqueId)
            }
            else return(distParam)
            
        }
        ################
        ## Function that will call clusterSingle for each row of param matrix
        ################
        paramFun <- function(i){
            par <- param[i,]
            totalArgs<-.makeArgsFromParam(par,
                mainClusterArgs=mainClusterArgs,
                seqArgs=seqArgs,
                subsampleArgs=subsampleArgs,
                clusterFunctionList=clusterFunctionList)
            if(!is.null(random.seed)) {
                set.seed(random.seed)
            }
            reduceMethod<-as.character(par[["reduceMethod"]])
            if(reduceMethod=="none")
                dat<-transformData(x,transFun=transFun, whichAssay=whichAssay)
            else if(isReducedDims(x,reduceMethod))
                dat<-t(reducedDim(x,reduceMethod)[,seq_len(par[["nReducedDims"]])] )
            else if(isFilterStats(x,reduceMethod))
                dat<-transformData( filterData(x, 
                    filterStats=reduceMethod, 
                    percentile=par[["nFilterDims"]]),
                    transFun=transFun, 
                    whichAssay=whichAssay)
            else stop("Internal error: reduceMethod value that not in filtering statistics or reducedDimNames")
            #(Note, computational inefficiency: means reordering each time, even if same filter. But not recalculating filter.)
            distFunction<-totalArgs$distFunction
            
            if(!is.null(distFunction)){
                ###FIXME: Error here! Needs to call name that is specific to combination of reduce/dist/algType in allDist[[distFunction]]
                passInput<- as.logical(gsub(" ","",par["passDistToInput"]))
                passMain<- as.logical(gsub(" ","",par["passDistToMain"]))
                distParamMatrix<-cbind(distFunction=as.character(distFunction),
                    nFilterDims=as.character(par[["nFilterDims"]]),
                    reduceMEthod=as.character(par[["reduceMethod"]]))
                dnm<-createDistNames(distParamMatrix,returnNames=TRUE)
                
                if(passInput){
                    out<-clusterSingle(inputMatrix=allDist[[dnm]],
                        inputType="diss",
                        sequential=totalArgs$sequential, 
                        subsample=totalArgs$subsample, 
                        reduceMethod="none",
                        mainClusterArgs=totalArgs$mainClusterArgs,
                        subsampleArgs=totalArgs$subsampleArgs, 
                        seqArgs=totalArgs$seqArgs,
                        transFun=function(x){x},
                        checkDiss=FALSE,
                        makeMissingDiss=FALSE,
                        warnings=parameterWarnings)
                    #because a distance, have to convert it back to a CE object.
                    #Using same internal function as clusterSingle
                    return(.convertOutListToCE(xMatrix=dat,
                            clustering=out$clustering,
                            clusterInfo=out$clusterInfo,
                            coClusterMatrix=NULL, 
                            clusterLabel="clusterSingle",
                            sequential=totalArgs$sequential, 
                            subsample=totalArgs$subsample,
                            transFun=function(x){x},
                            saveSubsamplingMatrix=FALSE, existingCE=NULL))
            
                
                }
                else{
                    if(passMain){
                        return(clusterSingle(inputMatrix=dat,
                            inputType="X",
                            sequential=totalArgs$sequential, 
                            subsample=totalArgs$subsample, 
                            reduceMethod="none",
                            mainClusterArgs=c(totalArgs$mainClusterArgs,
                                list(diss=allDist[[dnm]])),
                            subsampleArgs=totalArgs$subsampleArgs, 
                            seqArgs=totalArgs$seqArgs,
                            transFun=function(x){x},
                            checkDiss=FALSE,
                            makeMissingDiss=FALSE,
                            warnings=parameterWarnings))
                    }
                    else{
                        stop("Internal Error: didn't parse distance passing correctly.")
                    }
                    
                }

            }
            else
                return(clusterSingle(inputMatrix=dat, 
                    inputType="X", 
                    sequential=totalArgs$sequential, 
                    subsample=totalArgs$subsample,
                    reduceMethod="none",
                    mainClusterArgs=totalArgs$mainClusterArgs, 
                    subsampleArgs=totalArgs$subsampleArgs, 
                    seqArgs=totalArgs$seqArgs,
                    transFun=function(x){x},
                    checkDiss=FALSE,
                    makeMissingDiss=FALSE,
                    warnings=parameterWarnings))
        }
        
        ## Run this even if run=FALSE so see message.
        if(any(!is.na(param[,"distFunction"]))){
            ##Get the parameters that imply different datasets.
            uniqueId<-createDistNames(param)
            ## Use to assume only take distances on original data (or filtered version of it). But in fact, it meant previously if needed dissimilarity, would silently calculate it to each call (e.g. if diss->X). So now can get it if makeMissingDiss=TRUE
            if(verbose)
                cat(sprintf("Calculating the %s Requested Distance Matrices needed to run clustering comparisions (if 'distFunction=NA', then needed because of choices of clustering algorithms and 'makeMissingDiss=TRUE').\n",nrow(distParam)))
        }

        if(run){

               

            ##------------
            ##Calculate distances necessary only once
            ##------------
            if(any(!is.na(param[,"distFunction"]))){
                ##Get the parameters that imply different datasets.
                distParam<-createDistNames(param,returnNames=FALSE)
                allDist<-lapply(seq_len(nrow(distParam)),function(ii){
                    distFun<-as.character(distParam[ii,"distFunction"])
                    algCheckType<-distParam[ii,"distAlgType"]
                    redM<-as.character(distParam[ii,"reduceMethod"])
                    if(redM=="none"){
                        dat<-transformData(x,
                            transFun=transFun, whichAssay=whichAssay)
                        
                    }
                    else if(isFilterStats(x,redM)){
                        dat<-transformData( filterData(x, 
                            filterStats=redM,
                            percentile=distParam[ii,"nFilterDims"]),
                            transFun=transFun,whichAssay=whichAssay)
                    }
                    else if(makeMissingDiss & any(missDiss)){
                        ### FIXME: shouldn't be any(missDiss) -- should be only the one for this row!
                        dat<-t(reducedDim(x,type=redM))
                    }
                    else stop("Internal error: should have removed any combinations that included a dimensionality reduction and a distance unless due to makeMissingDiss=TRUE")
                    distMat<-.makeDiss(dat, distFunction=distFun, 
                        checkDiss=TRUE, algType=algCheckType)
                    return(distMat)
                })
            
                nms<-createDistNames(distParam,returnNames=TRUE)
                names(allDist)<-nms
                
                if(verbose) cat("\tdone.\n\n")
            }
            
            ##------------
            ## Run clustering for each row of param
            ##------------
            if(verbose) {
                cat("Running Clustering on Parameter Combinations...\n")
            }
            
            if(ncores>1) {
                out <- mclapply(seq_len(nrow(param)), FUN=paramFun, mc.cores=ncores, ...)
                nErrors <- which(sapply(out, function(x){inherits(x, "try-error")}))
                if(length(nErrors)>0) {
                    stop(length(nErrors)," parameter values (of ",length(out),") hit an error. The first was:\n",out[nErrors[1]])
                }
            } else {
                out <- lapply(seq_len(nrow(param)),FUN=paramFun)
            }
            if(verbose) {
                cat("done.\n")
            }
            clMat <- sapply(out, function(x){primaryCluster(x)})
            
            colnames(clMat) <- unname(cnames)
            pList <- lapply(seq_len(nrow(param)), function(i){
                x <- param[i,]
                names(x) <- colnames(param)
                return(x)})
            clInfo <- mapply(pList, out, FUN=function(x, y){
                c(list(choicesParam=x), clusteringInfo(y))
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
)

#' @rdname clusterMany
#' @export
setMethod(
    f = "clusterMany",
    signature = signature(x = "ClusterExperiment"),
    definition = function(x, reduceMethod="none", nFilterDims=defaultNDims(x,reduceMethod,type="filterStats"), nReducedDims=defaultNDims(x,reduceMethod,type="reducedDims"),
                          eraseOld=FALSE, ...)
    {
        if(any(c("transFun","isCount") %in% names(list(...))))
            stop("The internally saved transformation function of a ClusterExperiment object must be used when given as input and setting 'transFun' or 'isCount' for a 'ClusterExperiment' is not allowed.")
        outval<-clusterMany(as(x,"SingleCellExperiment"), reduceMethod=reduceMethod, nFilterDims=nFilterDims,
                            nReducedDims=nReducedDims, transFun=transformation(x), ...)
        
        if(is(outval,"ClusterExperiment")) {
            
            #outval<-.addBackSEInfo(newObj=outval,oldObj=x) #added to '.addNewResult'
            ##Check if clusterMany already ran previously
            x<-.updateCurrentWorkflow(x,eraseOld,newTypeToAdd="clusterMany",newLabelToAdd=NULL)
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


#' @export
#' @rdname clusterMany
setMethod(
    f = "clusterMany",
    signature = signature(x = "data.frame"),
    definition = function(x,...){clusterMany(data.matrix(x),...)}
)


.makeArgsFromParam<-function(paramRow,
    mainClusterArgs=NULL,
    seqArgs=NULL,
    subsampleArgs=NULL,
    clusterFunctionList=NULL){

    #The following fixes logical values that got messed up
    #otherwise adds a space before the TRUE and doesn't recognize.
    #well, sometimes. 
    removeSil <- as.logical(gsub(" ","",paramRow["removeSil"]))
    sequential <- as.logical(gsub(" ","",paramRow["sequential"]))
    subsample <- as.logical(gsub(" ","",paramRow["subsample"]))
    findBestK <- as.logical(gsub(" ","",paramRow["findBestK"]))
    mainClusterArgs[["findBestK"]] <- findBestK 
    mainClusterArgs[["removeSil"]] <- removeSil
    
    clusterFunctionName <- as.character(paramRow[["clusterFunction"]])
    distFunction <- if(!is.na(paramRow[["distFunction"]])) as.character(paramRow[["distFunction"]]) else NULL
    if(!is.null(clusterFunctionList)){
        clusterFunction<-clusterFunctionList[[clusterFunctionName]]        
    }
    else clusterFunction <- getBuiltInFunction(clusterFunctionName)
    reduceMethod<-as.character(paramRow[["reduceMethod"]])
    if(!is.na(paramRow[["k"]])){
        if(sequential) {
            seqArgs[["k0"]] <- paramRow[["k"]]
        } else{
            #to be safe, set both in case user set one.
            subsampleArgs[["clusterArgs"]][["k"]] <- paramRow[["k"]]
            mainClusterArgs[["clusterArgs"]][["k"]] <- paramRow[["k"]]
        }
    }
    mainClusterArgs[["clusterArgs"]][["alpha"]] <- paramRow[["alpha"]]
    seqArgs[["beta"]] <- paramRow[["beta"]]
    mainClusterArgs[["minSize"]] <- paramRow[["minSize"]]
    mainClusterArgs[["silCutoff"]] <- paramRow[["silCutoff"]]
    mainClusterArgs[["clusterFunction"]]<-clusterFunction
    seqArgs[["verbose"]]<-FALSE
    
    return(list(
        distFunction=distFunction,
        sequential=sequential,
        subsample=subsample,
        mainClusterArgs=mainClusterArgs,
        subsampleArgs=subsampleArgs,
        seqArgs=seqArgs
        ))
}
