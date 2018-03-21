#' @useDynLib clusterExperiment
#' @importFrom Rcpp sourceCpp
NULL


#' Cluster subsamples of the data
#'
#' Given input data, this function will subsample the samples, cluster the
#' subsamples, and return a \code{n x n} matrix with the probability of
#' co-occurance.
#' @name subsampleClustering
#' @param x the data on which to run the clustering (samples in columns).
#' @param diss a dissimilarity matrix on which to run the clustering.
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
#' @param largeDataset logical indicating whether a more memory-efficient version
#'   should be used because the dataset is large. This is a beta option, and is
#'   in the process of being tested before it becomes the default.
#' @param doGC logical indicating whether frequent calls to gc should be
#'  implemented in children processes (i.e. when ncores>1) to free up memory
#'  for the other processes.
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
#' coOccur <- subsampleClustering(clusterFunction="kmeans", x=simData,
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
  definition=function(clusterFunction, x=NULL,diss=NULL,distFunction=NA,clusterArgs=NULL,
                      classifyMethod=c("All","InSample","OutOfSample"),
                      resamp.num = 100, samp.p = 0.7,ncores=1,checkArgs=TRUE,checkDiss=TRUE,largeDataset=FALSE,doGC=FALSE,which_implementation=c("R", "Csimple", "Cmemory"),... )
  {

    ## Slows down enormously to do rm and gc, so only if largeDataset=TRUE and ncores>1
    #	doGC<-largeDataset & ncores>1
    #######################
    ### Check both types of inputs and create diss if needed, and check it.
    #######################
    input<-.checkXDissInput(x,diss,inputType=clusterFunction@inputType,checkDiss=checkDiss)
    classifyMethod<-match.arg(classifyMethod)
    if(classifyMethod %in% c("All","OutOfSample") && is.null(clusterFunction@classifyFUN)){
      classifyMethod<-"InSample" #silently change it...
    }
    else{
      inputClassify<-.checkXDissInput(x, diss, inputType=clusterFunction@inputClassifyType, checkDiss=FALSE) #don't need to check it twice!
    }
    if((input=="X" & clusterFunction@inputType=="diss") || (classifyMethod!="InSample" && inputClassify=="X" && clusterFunction@inputClassifyType=="diss")){
      diss<-.makeDiss(x,distFunction=distFunction,checkDiss=checkDiss,algType=clusterFunction@algorithmType)
      if(input=="X") input<-"diss"
      if(inputClassify=="X") inputClassify<-"diss"
    }
    #-----
    # Other Checks
    #-----
    reqArgs<-requiredArgs(clusterFunction)
    if(!all(reqArgs %in% names(clusterArgs))) stop(paste("For this clusterFunction algorithm type ('",algorithmType(clusterFunction),"') must supply arguments",reqArgs,"as elements of the list of 'clusterArgs'"))

    #-----
    # Basic parameters, subsamples
    #-----
    if(input %in% c("X","both")) N <- dim(x)[2] else N<-dim(diss)[2]
    subSize <- round(samp.p * N)
    idx<-replicate(resamp.num,sample(1:N,size=subSize)) #each column a set of indices for the subsample.

    #-----
    # Function that calls the clustering for each subsample
    # Called over a loop (lapply or mclapply)
    #-----
    perSample<-function(ids){
      ## Calls rm and gc frequently to free up memory
      ## (mclapply child processes don't know when other processes are using large memory. )

      ##----
      ##Cluster part of subsample
      ##----
      argsClusterList <- .makeDataArgs(dataInput=input,
                                       funInput=clusterFunction@inputType,
                                       xData=x[,ids,drop=FALSE],
                                       dissData=diss[ids,ids,drop=FALSE])

      #if doing InSample, do cluster.only because will be more efficient, e.g. pam and kmeans.
      argsClusterList <- c(argsClusterList,
                           list("checkArgs"=checkArgs,
                                "cluster.only"=(classifyMethod=="InSample")))
      result <- do.call(clusterFunction@clusterFUN,
                        c(argsClusterList,clusterArgs))

      if(doGC){
        rm(argsClusterList)
        gc()
      }

      ##----
      ##Classify part of subsample
      ##----
      if(classifyMethod=="All"){
        argsClassifyList <- .makeDataArgs(dataInput=inputClassify,
                                          funInput=clusterFunction@inputClassifyType,
                                          xData=x, dissData=diss)
        classX <- do.call(clusterFunction@classifyFUN,
                          c(argsClassifyList,list(clusterResult=result)))

        if(doGC){
          rm(argsClassifyList)
          gc()
        }
      }

      if(classifyMethod=="OutOfSample"){
        argsClassifyList <- .makeDataArgs(dataInput=inputClassify,
                                          funInput=clusterFunction@inputClassifyType,
                                          xData=x[,-ids,drop=FALSE],
                                          dissData=diss[-ids,-ids,drop=FALSE])
        classElse <- do.call(clusterFunction@classifyFUN,
                             c(argsClassifyList, list(clusterResult=result)))

        if(doGC){
          rm(argsClassifyList)
          gc()
        }

        classX <- rep(NA,N)
        classX[-ids] <- classElse
      }

      if(classifyMethod=="InSample"){
        classX <- rep(NA,N)

        #methods that do not have
        if(is.list(result)){
          #the next shouldn't happen any more because should be cluster.only=TRUE
          # if("clustering" %in% names(result)){
          # 	classX[ids]<-result$clustering
          # }
          # 		  	  else{
          if(clusterFunction@outputType=="list"){
            resultVec <- .convertClusterListToVector(result,N=length(ids))
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
      if(!largeDataset){

        D <- lapply(seq_len(N), function(i) classX == classX[i])
        return(unlist(D))

      }
      else{

        if(which_implementation == "Csimple") {
          ## we just need to return the vector of cluster labels
          return(classX)
        } else {
          #instead return
          # 1) one vector of length na.omit(classX) of the original indices, where ids in clusters are adjacent in the vector and
          # 2) another vector of length K indicating length of each cluster (allows to decode where the cluster stopes in the above vector),
          # What does this do with NAs? Removes them -- not included.
          clusterIds<-unlist(tapply(1:N,classX,function(x){x},simplify=FALSE),use.names=FALSE)
          clusterLengths<-tapply(1:N,classX,length)
          return(list(clusterIds=clusterIds,clusterLengths=clusterLengths))
        }
      }
    }

    if(ncores==1){

      DList<-apply(idx,2,perSample)

    }
    else{
      DList<-parallel::mclapply(1:ncol(idx), function(nc){ perSample(idx[,nc]) }, mc.cores=ncores,...)
      DList <- simplify2array(DList)
    }
    #N large: get rid of these big matrices from memory
    if(!is.null(diss)){
      idnames<-colnames(diss)
      rm(diss)
    }
    if(!is.null(x)){
      idnames<-colnames(x)
      rm(x)
    }

    if(!largeDataset){

      Dvec <- rowMeans(DList, na.rm = TRUE)
      Dbar <- matrix(Dvec, ncol = N, nrow = N)

    }
    else{
      #############
      #Need to calculate number of times pairs together without building large NxN matrix
      #############

      if(which_implementation == "Csimple") {
        Dbar <- search_pairs(t(DList))
        Dbar<-Dbar+t(Dbar)
        Dbar[is.na(Dbar)]<-0
        diag(Dbar)<-1
      } else {
        if(which_implementation == "R"){
          ###
          # the result of pairList is the upper-triangle of eventual NxN matrix
          ###
          if(ncores==1){
            pairList<-lapply(2:N,function(jj){searchForPairs(jj,clusterList=DList,N=N)})
          }
          else{
            pairList<-parallel::mclapply(2:N,function(jj){searchForPairs(jj,clusterList=DList,N=N)},mc.cores=ncores,...)
          }
          pairList<-unlist(pairList)
        } else {
          matResults<-subsampleLoop(DList, N)
          ord<-order(matResults[,2],matResults[,1])
          pairList<- matResults[ord,3]/matResults[ord,4]
        }

        #Create NxN matrix
        Dbar<-matrix(0,N,N)
        Dbar[upper.tri(Dbar, diag = FALSE)]<-pairList
        Dbar<-Dbar+t(Dbar)
        Dbar[is.na(Dbar)]<-0
        diag(Dbar)<-1

      }

    }

    rownames(Dbar)<-colnames(Dbar)<-idnames
    return(Dbar)
  })
