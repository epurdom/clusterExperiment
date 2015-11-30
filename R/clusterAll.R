
#' General wrap-around for all clustering methods we're trying.
#' 
#' Given a nxp data matrix, this function will find clusters
#' 
#' @param x the data on which to run the clustering (samples in rows)
#' @param subsample logical as to whether to subsample via \code{\link{subsampleClustering}} to get the distance matrix at each iteration; otherwise the distance matrix is dist(x)
#' @param sequential logical whether to use the sequential strategy.
#' @param clusterFunction passed to \code{\link{clusterD}} option 'clusterFunction' to indicate method of clustering, see \code{\link{clusterD}}
#' @param clusterDArgs list of additional arguments to be passed to \code{\link{clusterD}} 
#' @param subsampleArgs list of arguments to be passed to \code{\link{subsampleClustering}}
#' @param seqArgs list of additional arguments to be passed to \code{\link{seqCluster}}
#' 
#' @details If sequential=TRUE, the sequential clustering controls the 'k' argument of the underlying clustering so setting 'k=' in the list given to clusterDArgs or subsampleArgs will not do anything and will produce a warning to that effect. 
#'
#' @return A list with values 
#' \itemize{

#' \item{\code{clustering}}{a vector of length equal to nrows(x) giving the integer-valued cluster ids for each sample. The integer values are assigned in the order that the clusters were found, if sequential=TRUE. "-1" indicates the sample was not clustered.}

#' \item{\code{clusterInfo}}{if sequential=TRUE and clusters were successfully found, a matrix of information regarding the algorithm behavior for each cluster (the starting and stopping K for each cluster, and the number of iterations for each cluster).}

#' \item{\code{whyStop}}{if sequential=TRUE and clusters were successfully found, a character string explaining what triggered the algorithm to stop.}
#' }
#'
#' @examples
#' data(simData)
#' \dontrun{
#' #following code takes some time. 
#' #use clusterAll to do sequential clustering 
#' #(same as example in seqCluster only using clusterAll ...)
#' set.seed(44261)
#' clustSeqHier_v2<-clusterAll(simData,clusterFunction="hierarchical",
#' sequential=TRUE,subsample=TRUE,
#'	subsampleArgs=list(resamp.n=100,samp.p=0.7,clusterFunction="kmeans",
#'	clusterArgs=list(nstart=10)), seqArgs=list(beta=0.8,k0=5),
#'	clusterDArgs=list(minSize=5))
#' }
#' #use clusterAll to do just clustering k=3 with no subsampling
#' clustNothing<-clusterAll(simData,clusterFunction="pam",subsample=FALSE,
#' sequential=FALSE, clusterDArgs=list(k=3))

clusterAll<-function(x,  subsample=TRUE, sequential=FALSE, clusterFunction=c("tight","hierarchical","pam","kmeans"),  clusterDArgs=NULL,subsampleArgs=NULL,seqArgs=NULL) 
{
    if(!is.function(clusterFunction)){
		clusterFunction<-match.arg(clusterFunction)
		if(!subsample & clusterFunction !="pam") stop("If not subsampling, clusterFunction must be 'pam'")
		typeAlg<-.checkAlgType(clusterFunction)
	}
	else{
		if(! "typeAlg" %in% clusterDArgs) stop("if you provide your own clustering algorithm to be passed to clusterD, then you must specify 'typeAlg' in clusterDArgs")
		else typeAlg<-clusterDArgs[["typeAlg"]]
	}
	if(typeAlg == "K"){
		if("findBestK" %in% names(clusterDArgs) & !subsample & sequential){
			if(clusterDArgs[["findBestK"]]) stop("Cannot do sequential clustering where subsample=FALSE and 'findBestK=TRUE' is passed via clusterDArgs. See help documentation.")
		}
		
	}	
    N <- dim(x)[1]
	if(sequential){
		if(is.null(seqArgs)) stop("must give seqArgs so as to identify k0")
		if(!"k0"%in%names(seqArgs)) stop("seqArgs must contain element 'k0'")
		seqOut<-do.call("seqCluster",c(list(x=x,subsample=subsample,subsampleArgs=subsampleArgs,clusterDArgs=clusterDArgs,clusterFunction=clusterFunction),seqArgs))
		#	browser()
		return(seqOut)
	}
	else{
		if(subsample){
			if(is.null(subsampleArgs) || !("k" %in% names(subsampleArgs))){
          if(!is.null(clusterDArgs) & ("k" %in% names(clusterDArgs))){
            #give by default the clusterDArgs to subsampling. 
              warning("did not give 'k' in 'subsampleArgs'. Set to 'k' argument in 'clusterDArgs'")
              if(is.null(subsampleArgs)) subsampleArgs<-list("k"=clusterDArgs[["k"]])
              else subsampleArgs[["k"]]<-clusterDArgs[["k"]]
          }
			  	else	stop("if not sequential and do subsampling, must pass 'k' in subsampleArgs")
			} 
		}
		else if(typeAlg=="K" && !is.null(clusterDArgs) &&  !"k" %in% names(clusterDArgs)){
			#if don't specify k, then must have findBestK=TRUE in clusterDArgs; is by default, so only need to check that if specified it, set it to TRUE
			if("findBestK" %in% names(clusterDArgs) && !clusterDArgs[["findBestK"]]) stop("if not sequential and clusterFunction is of type 'K' (e.g. pam) and findBestK=FALSE in clusterDArgs, must pass 'k' via clusterDArgs list")
				}
		finalClusterList<-.clusterWrapper(x,clusterFunction=clusterFunction, subsample=subsample,  subsampleArgs=subsampleArgs,clusterDArgs=clusterDArgs,typeAlg=typeAlg)
		return(list("clustering"=.convertClusterListToVector(finalClusterList,N)))
	}
}







