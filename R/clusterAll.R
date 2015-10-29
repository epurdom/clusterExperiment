
#' General wrap-around for all clustering methods we're trying.
#' 
#' Given a nxp data matrix, this function will find clusters
#' 
#' @param x the data on which to run the clustering (samples in rows)
#' @param subsample logical as to whether to subsample to get the distance matrix; otherwise the distance matrix is dist(x)
#' @param sequential logical whether to use the sequential strategy.
#' @param clusterMethod method for clustering distance matrix (either from subsampling or from dist(x)); if not subsampling, must be “pam”. Passed to clusterDMat option 'method'.
#' @param clusterArgs list of additional arguments to be passed to clusterDMat
#' @param subsampleArgs list of arguments to be passed to subsamplingClustering
#' 
#' @details If sequential=TRUE, the sequential clustering controls the 'k' argument of the underlying pam/kmeans clustering so setting 'k=' in the list given to clusterArgs will not do anything and will produce a warning to that effect. 
#'
#' @details If subsampleClusterMethod="none" and 'findBestK=FALSE' is passed to subsampleClusterArgs, then each iteration will do pam on dist(x) iterating over k. However, if subsampleClusterMethod="none", you should not set 'findBestK=TRUE' (otherwise clustering dist(x) will be the same for iterating over different k and there is no method implemented to change the choice of how to remove a cluster other than similarity as you change k); an error message will be given. However, if subsampleClusterMethod="pam" (i.e. apply pam to clustering of the distance matrix after subsampling) passing either 'findBestK=TRUE' or 'findBestK=FALSE' will function as expected. Note that the default range of k values for findBestK=TRUE is dependent on the input k; if you want the same range of k values, you should explicitly set kRange via the argument subsampleClusterArgs.
#'
#' @return A list with values 
#' \itemize{

#' \item{\code{clustering}}{a vector of length equal to nrows(x) giving the integer-valued cluster ids for each sample. The integer values are assigned in the order that the clusters were found, if sequential=TRUE. "-1" indicates the sample was not clustered.}

#' \item{\code{clusterInfo}}{if sequential=TRUE and clusters were successfully found, a matrix of information regarding the algorithm behavior for each cluster (the starting and stopping K for each cluster, and the number of iterations for each cluster).}

#' \item{\code{whyStop}}{if sequential=TRUE and clusters were successfully found, a character string explaining what triggered the algorithm to stop.}
#' }
#'
#' @examples
#' 
#use clusterAll to do the same thing...
#'set.seed(44261)
#'clustSeqHier_v2<-clusterAll(z,clusterMethod="hierarchical",sequential=TRUE,subsample=TRUE,
#'	subsampleArgs=list(resamp.n=100,samp.p=0.7),
#'	subsampleClusterArgs=list(nstart=10),
#'	seqArgs=list(beta=0.8,k0=10),
#'	clusterArgs=list(min.size=5))
#' #use clusterAll to do just clustering k=3 with no subsampling
#' clustNothing<-clusterAll(z,clusterMethod="pam",subsample=FALSE,sequential=FALSE
#'	clusterArgs=list(k=3))

clusterAll<-function(x,  subsample=TRUE, sequential=FALSE, clusterMethod=c("tight","hierarchical","pam","kmeans"),  clusterArgs=NULL,subsampleArgs=NULL,seqArgs=NULL) 
{
	#for now, if use pam for subsampleClusterMethod, just use given k.
	
    clusterMethod<-match.arg(clusterMethod)
	if(!subsample & clusterMethod !="pam") stop("If not subsampling, clusterMethod must be 'pam'")
    original.data <- x
    N <- dim(x)[1]
	if(sequential){
		if(is.null(seqArgs)) stop("must give seqArgs so as to identify k0")
		if(!"k0"%in%names(seqArgs)) stop("seqArgs must contain element 'k0'")
		seqOut<-do.call("seqCluster",c(list(x=x,subsample=subsample,subsampleArgs=subsampleArgs,clusterArgs=clusterArgs,clusterMethod=clusterMethod),seqArgs))
		return(seqOut)
	#	browser()
	}
	else{
		if(subsample){
			if(is.null(subsampleArgs) || !("k" %in% names(subsampleArgs))) stop("if not sequential, must pass 'k' in subsampleArgs")
		}
		else if(clusterMethod=="pam" && !is.null(clusterArgs) &&  !"k" %in% names(clusterArgs)){
			if("findBestK" %in% names(clusterArgs) && !clusterArgs[["findBestK"]]) stop("if not sequential and clusterMethod='pam' and findBestK=FALSE in clusterArgs, must pass 'k' via clusterArgs list")
				}
		finalClusterList<-.clusterWrapper(x,clusterMethod=clusterMethod, subsample=subsample,  subsampleArgs=subsampleArgs,clusterArgs=clusterArgs)
		return(list("clustering"=.convertClusterListToVector(finalClusterList,N)))
		#browser()
	}
}







