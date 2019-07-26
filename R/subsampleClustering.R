#' @useDynLib clusterExperiment, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' Cluster subsamples of the data
#'
#' Given input data, this function will subsample the samples, cluster the
#' subsamples, and return a \code{n x n} matrix with the probability of
#' co-occurance.
#' @name subsampleClustering
#' @param clusterFunction a \code{\link{ClusterFunction}} object that defines
#'   the clustering routine. See \code{\link{ClusterFunction}} for required
#'   format of user-defined clustering routines. User can also give a character
#'   value to the argument \code{clusterFunction} to indicate the use of
#'   clustering routines provided in package. Type
#'   \code{\link{listBuiltInFunctions}} at command prompt to see the built-in
#'   clustering routines. If \code{clusterFunction} is missing, the default is
#'   set to "pam".
#' @param clusterArgs a list of parameter arguments to be passed to the function
#'   defined in the \code{clusterFunction} slot of the \code{ClusterFunction}
#'   object. For any given \code{\link{ClusterFunction}} object, use function
#'   \code{\link{requiredArgs}} to get a list of required arguments for the
#'   object.
#' @param resamp.num the number of subsamples to draw.
#' @param samp.p the proportion of samples to sample for each subsample.
#' @param classifyMethod method for determining which samples should be used in
#'   calculating the co-occurance matrix. "All"= all samples, "OutOfSample"=
#'   those not subsampled, and "InSample"=those in the subsample.  See details
#'   for explanation.
#' @param ncores integer giving the number of cores. If ncores>1, mclapply will
#'   be called.
#' @param ... arguments passed to mclapply (if ncores>1).
#' @inheritParams mainClustering
#' @inheritParams clusterSingle
#'
#' @details \code{subsampleClustering} is not usually called directly by the
#'   user. It is only an exported function so as to be able to clearly document
#'   the arguments for \code{subsampleClustering}  which can be passed via the
#'   argument \code{subsampleArgs} in functions like \code{\link{clusterSingle}}
#'   and \code{\link{clusterMany}}.
#' @details \code{requiredArgs:} The choice of "All" or "OutOfSample" for
#'   \code{requiredArgs} require the classification of arbitrary samples not
#'   originally in the clustering to clusters; this is done via the classifyFUN
#'   provided in the \code{\link{ClusterFunction}} object. If the
#'   \code{\link{ClusterFunction}} object does not have such a function to
#'   define how to classify into a cluster samples not in the subsample that
#'   created the clustering then \code{classifyMethod} must be
#'   \code{"InSample"}. Note that if "All" is chosen, all samples will be
#'   classified into clusters via the classifyFUN, not just those that are
#'   out-of-sample; this could result in different assignments to clusters for
#'   the in-sample samples than their original assignment by the clustering
#'   depending on the classification function. If you do not choose 'All',it is
#'   possible to get NAs in resulting S matrix (particularly if when not enough
#'   subsamples are taken) which can cause errors if you then pass the resulting
#'   D=1-S matrix to \code{\link{mainClustering}}. For this reason the default is
#'   "All".
#' @return A \code{n x n} matrix of co-occurances, i.e. a symmetric matrix with
#'   [i,j] entries equal to the percentage of subsamples where the ith and jth
#'   sample were clustered into the same cluster. The percentage is only out of
#'   those subsamples where the ith and jth samples were both assigned to a
#'   clustering. If \code{classifyMethod=="All"}, this is all subsamples for all
#'   i,j pairs. But if \code{classifyMethod=="InSample"} or
#'   \code{classifyMethod=="OutOfSample"}, then the percentage is only taken on
#'   those subsamples where the ith and jth sample were both in or out of
#'   sample, respectively, relative to the subsample.
#'
#' @examples
#'\dontrun{
#' #takes a bit of time, not run on checks:
#' data(simData)
#' coOccur <- subsampleClustering( inputMatrix=simData, inputType="X",
#' clusterFunction="kmeans",
#' clusterArgs=list(k=3,nstart=10), resamp.n=100, samp.p=0.7)
#'
#' #visualize the resulting co-occurance matrix
#' plotHeatmap(coOccur)
#'}
#' @aliases subsampleClustering,character-method
#' @export
setMethod(
    f = "subsampleClustering",
    signature = signature(clusterFunction = "character"),
    definition = function(clusterFunction,...){
        subsampleClustering(getBuiltInFunction(clusterFunction),...)
        
    }
)

# #' @rdname subsampleClustering
# #' @export
# setMethod(
# f = "subsampleClustering",
# signature = signature(clusterFunction = "missing"),
# definition = function(clusterFunction,...){
# 	subsampleClustering(clusterFunction="pam",...)
# }
# )

#' @rdname subsampleClustering
#' @export
setMethod(
    f = "subsampleClustering",
    signature = signature(clusterFunction = "ClusterFunction"),
    definition=function(clusterFunction, inputMatrix,inputType, clusterArgs=NULL,
                        classifyMethod=c("All","InSample","OutOfSample"),
                        resamp.num = 100, samp.p = 0.7,
						ncores=1,warnings=TRUE,... )
    {
        
        classifyMethod<-match.arg(classifyMethod)
        ###########################
        ######## CHECKS
        ###########################
        if(missing(inputType)) stop("Internal error: inputType was not passed to subsampling step")
        moreArgs<-list(...)
        subsampleArgs<-c(list(clusterFunction=clusterFunction, clusterArgs=clusterArgs, classifyMethod=classifyMethod),moreArgs)
        checkOut<-.checkArgs(inputType=inputType, 
                        main=FALSE, subsample=TRUE, sequential=FALSE,
                        mainClusterArgs=NULL,
                        subsampleArgs=subsampleArgs, 
                        warn=warnings)		
        if(is.character(checkOut)) stop(checkOut)
        else{
            subsampleArgs<-checkOut$subsampleArgs
        }
        classifyMethod<-subsampleArgs[["classifyMethod"]]
        clusterFunction<-subsampleArgs[["clusterFunction"]]
        clusterArgs<-subsampleArgs[["clusterArgs"]]
        inputType<-subsampleArgs[["inputType"]]
        
        #-----
        # Basic parameters, subsamples
        #-----
        N <- dim(inputMatrix)[2]
        subSize <- round(samp.p * N)
        idx<-replicate(resamp.num,sample(seq_len(N),size=subSize))
        #each column a set of indices for the subsample.
        
        ###########################
        # Function that calls the clustering for each subsample
        # Called over a loop (lapply or mclapply)
        ###########################
        perSample<-function(ids){
            ##----
            ##Cluster subsample
            ##----
            #if doing InSample, do cluster.only because will be more efficient, e.g. pam and kmeans.
            argsClusterList <- list(inputType=inputType,
                "cluster.only"=(classifyMethod=="InSample"))
            if(inputType=="diss") 
                argsClusterList<-c(argsClusterList, 
                    list(inputMatrix=inputMatrix[ids,ids,drop=FALSE]))
            else
                argsClusterList<-c(argsClusterList, 
                    list(inputMatrix=inputMatrix[,ids,drop=FALSE]))
            result <- do.call(clusterFunction@clusterFUN,
                              c(argsClusterList,clusterArgs))
            
            ##----
            ##Classify subsample
            ##----
            if(classifyMethod=="All"){
                ##FIXME make this actually only recluster the out-of-sample, but otherwise use the clustering results from the method!
				classX <- do.call(clusterFunction@classifyFUN, list(inputMatrix=inputMatrix,inputType=inputType,clusterResult=result))
            }
            if(classifyMethod=="OutOfSample"){
                classElse <- do.call(clusterFunction@classifyFUN,
                                     list(clusterResult=result,
										 inputMatrix=if(inputType!="diss") inputMatrix[,-ids,drop=FALSE] else inputMatrix[-ids,-ids,drop=FALSE]))
                
                classX <- rep(NA,N)
                classX[-ids] <- classElse
            }
            if(classifyMethod=="InSample"){
                classX <- rep(NA,N)
                if(is.list(result)){
                    if(clusterFunction@outputType=="list"){
                        resultVec <- .clusterListToVector(result,N=length(ids))
                        classX[ids] <- resultVec
                    } else {
                        stop("The clusterFunction given to subsampleClustering returns a list when cluster.only=FALSE but does not have a named element 'clustering' nor outputType='list'")
                    }
                    #			  }
                } else{
                    classX[ids]<-result
                }
            }
            
            #classX is length N
            #classX has NA if method does not classify all of the data.
            return(classX)
        }
        
        if(ncores==1){
            DList<-apply(idx,2,perSample)
        }
        else{
            DList<-parallel::mclapply(seq_len(ncol(idx)), function(nc){ perSample(idx[,nc]) }, mc.cores=ncores,...)
            DList <- simplify2array(DList)
        }
        idnames<-colnames(inputMatrix)
        # DList is a NxB matrix
        return(DList)
    })
