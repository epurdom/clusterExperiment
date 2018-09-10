#' @title Make hierarchy of set of clusters
#'
#' @aliases makeDendrogram,ClusterExperiment-method
#' @description Makes a dendrogram of a set of clusters based on hclust on the
#'   medoids of the cluster.
#' @param x data to define the medoids from. Matrix and
#'   \code{\link{ClusterExperiment}} supported.
#' @param cluster A numeric vector with cluster assignments. If x is a
#'   ClusterExperiment object, cluster is automatically the primaryCluster(x).
#'   ``-1'' indicates the sample was not assigned to a cluster.
#' @param unassignedSamples how to handle unassigned samples("-1") ; only
#'   relevant for sample clustering. See details.
#' @param reduceMethod character A character identifying what type of
#'   dimensionality reduction to perform before clustering. Can be either a
#'   value stored in either of reducedDims or filterStats slot or a built-in
#'   diminsionality reduction/filtering. The option "coCluster" will use the
#'   co-Clustering matrix stored in the 'coClustering' slot of the
#'   \code{ClusterExperiment} object
#' @param nDims The number of dimensions to keep from \code{reduceMethod}. If
#'   missing calls \code{\link{defaultNDims}}.
#' @param whichCluster an integer index or character string that identifies
#'   which cluster should be used to make the dendrogram. Default is
#'   primaryCluster.
#' @param ... for makeDendrogram, if signature \code{matrix}, arguments passed
#'   to hclust; if signature \code{ClusterExperiment} passed to the method for
#'   signature \code{matrix}. For plotDendrogram, passed to
#'   \code{\link{plot.dendrogram}}.
#' @inheritParams clusterSingle
#' @inheritParams reduceFunctions
#' @details The function returns two dendrograms (as a list if x is a matrix or
#'   in the appropriate slots if x is ClusterExperiment). The cluster dendrogram
#'   is created by applying \code{\link{hclust}} to the medoids of each cluster.
#'   In the sample dendrogram the clusters are again clustered, but now the
#'   samples are also part of the resulting dendrogram. This is done by giving
#'   each sample the value of the medoid of its cluster.
#' @details The argument \code{unassignedSamples} governs what is done with
#'   unassigned samples (defined by a -1 cluster value). If
#'   unassigned=="cluster", then the dendrogram is created by hclust of the
#'   expanded medoid data plus the original unclustered observations. If
#'   \code{unassignedSamples} is "outgroup", then all unassigned samples are put
#'   as an outgroup. If the \code{x} object is a matrix, then
#'   \code{unassignedSamples} can also be "remove", to indicate that samples
#'   with "-1" should be discarded. This is not a permitted option, however,
#'   when \code{x} is a \code{ClusterExperiment} object, because it would return
#'   a dendrogram with less samples than \code{NCOL(x)}, which is not permitted
#'   for the \code{@dendro_samples} slot.
#' @details If any merge information is stored in the input object, it will be
#'   erased by a call to makeDendrogram.
#' @return If x is a matrix, a list with two dendrograms, one in which the
#'   leaves are clusters and one in which the leaves are samples. If x is a
#'   ClusterExperiment object, the dendrograms are saved in the appropriate
#'   slots.
#'
#' @export
#' @seealso makeFilterStats, makeReducedDims
#' @examples
#' data(simData)
#'
#' #create a clustering, for 8 clusters (truth was 3)
#' cl <- clusterSingle(simData, subsample=FALSE,
#' sequential=FALSE, mainClusterArgs=list(clusterFunction="pam", clusterArgs=list(k=8)))
#'
#' #create dendrogram of clusters:
#' hcl <- makeDendrogram(cl)
#' plotDendrogram(hcl)
#' plotDendrogram(hcl, leafType="samples",plotType="colorblock")
#'
#' @name makeDendrogram
#' @rdname makeDendrogram
setMethod(
  f = "makeDendrogram",
  signature = "ClusterExperiment",
  definition = function(x, whichCluster="primaryCluster",reduceMethod="mad",
                        nDims=defaultNDims(x,reduceMethod),filterIgnoresUnassigned=TRUE,
                        unassignedSamples=c("outgroup", "cluster"),
                        whichAssay=1,...)
  {
		passedArgs<-list(...)
					
		checkIgnore<-.depricateArgument(passedArgs=passedArgs,"filterIgnoresUnassigned","ignoreUnassignedVar")
		if(!is.null(checkIgnore)){
			passedArgs<-checkIgnore$passedArgs
			filterIgnoresUnassigned<-checkIgnore$val
		}
    unassignedSamples<-match.arg(unassignedSamples)
    whCl<-.convertSingleWhichCluster(x,whichCluster,passedArgs)
    cl<-clusterMatrix(x)[,whCl]
    ##erase merge information
    if(!is.na(mergeClusterIndex(x)) || !is.na(x@merge_dendrocluster_index)) x<-.eraseMerge(x)

    ########
    ##Transform the data
    ########
    if(length(reduceMethod)>1) stop('makeDendrogram only takes one choice of "reduceMethod" as argument')
    if(reduceMethod!="coCluster"){
      #need to change name of reduceMethod to make it match the
      #clustering information if that option chosen.
      datList<-getReducedData(object=x,whichCluster=whCl,reduceMethod=reduceMethod,
                              nDims=nDims,filterIgnoresUnassigned=TRUE,  whichAssay=whichAssay,returnValue="list")
      x<-datList$objectUpdate
      dat<-datList$dat
      
      outlist <- do.call("makeDendrogram",c(list(
				x=dat, 
				cluster=cl,calculateSample=TRUE,
				unassignedSamples=unassignedSamples),
				passedArgs))
    }
    else{
      if(is.null(x@coClustering)) stop("Cannot choose 'coCluster' if 'coClustering' slot is empty. Run makeConsensus before running 'makeDendrogram' or choose another option for 'reduceMethod'")
      if(is.null(dimnames(x@coClustering))) stop("This ClusterExperiment object was made with an old version of clusterExperiment and did not give dimnames to the coClustering slot.")
     outlist<-do.call("makeDendrogram",c(list(
			  x=as.dist(1-x@coClustering),
				cluster=cl,calculateSample=TRUE,
				unassignedSamples=unassignedSamples),
				passedArgs)) 
    }
		
    x@dendro_clusters <- outlist$clusters
    x@dendro_index<-whCl

		
	    x@dendro_samples <- outlist$samples
	    x@dendro_outbranch<- any(cl<0) & unassignedSamples=="outgroup"
			
		
		ch<-.checkDendrogram(x)
    if(!is.logical(ch)) stop(ch)
			return(x)
  })



#' @rdname makeDendrogram
#' @export
setMethod(
  f = "makeDendrogram",
  signature = "dist",
  definition = function(x, cluster,
                        unassignedSamples=c("outgroup", "cluster", "remove"),
                        calculateSample=TRUE,...) {
    unassigned <- match.arg(unassignedSamples)
    cl <- cluster
    nSamples<-attributes(x)$Size
    if(nSamples != length(cl)) {
      stop("cluster must be the same length as the number of samples")
    }
    if(is.null(attributes(x)$Labels)) {
      attributes(x)$Labels <- as.character(seq_len(nSamples))
    }
    
    clNum<-.convertToNum(cl)
    
    #############
    # Cluster dendrogram
    #############
    whKeep <- which(clNum >= 0) #remove -1, -2
    if(length(whKeep) == 0) stop("all samples have clusterIds<0")
    if(length(unique(cl[whKeep]))==1) stop("Only 1 cluster given. Can not make a dendrogram.")
    clFactor <- factor(cl[whKeep])
    
    #each pair of clusters, need to get median of the distance values
    #do a double by, just returning the values as a vector, and then take the median
    medoids<-do.call("rbind", by(as.matrix(x)[whKeep,whKeep], clFactor, function(z){
      out<-as.vector(by(t(z),clFactor,function(y){median(as.vector(unlist(y)))}))
      names(out)<-levels(clFactor)
      return(out)
    }))
    diag(medoids)<-0 #incase the median of within is not zero...
    rownames(medoids) <- levels(clFactor)
    colnames(medoids) <- levels(clFactor)
    nPerCluster <- table(clFactor)
    clusterD<-as.dendrogram(stats::hclust(as.dist(medoids),members=nPerCluster,...))
    #############
    # Samples dendrogram
    #############
    fullD<-.makeSampleDendro(x,clusterDendro=clusterD, cl=clNum,type=c("dist"), unassignedSamples=unassigned,sampleEdgeLength=0,  outbranchLength=1,calculateSample=calculateSample)
    return(list(samples=fullD,clusters=clusterD))
  })



#' @rdname makeDendrogram
#' @importFrom DelayedArray DelayedArray
#' @export
setMethod(
  f = "makeDendrogram",
  signature = "matrixOrHDF5",
  definition = function(x, cluster,
                        unassignedSamples=c("outgroup", "cluster", "remove"),
                        calculateSample=TRUE,...) {
    unassigned <- match.arg(unassignedSamples)
    cl <- cluster
    if(ncol(x) != length(cl)) {
      stop("cluster must be the same length as the number of samples")
    }
    if(is.null(colnames(x))) {
      colnames(x) <- as.character(seq_len(ncol(x)))
    }
    
    clNum<-.convertToNum(cl)
    
    #############
    # Cluster dendrogram
    #############
    whKeep <- which(clNum >= 0) #remove -1, -2
    if(length(whKeep) == 0) stop("all samples have clusterIds<0")
    if(length(unique(cl[whKeep]))==1) stop("Only 1 cluster given. Can not make a dendrogram.")
    clFactor <- factor(cl[whKeep])
    
    medoids <- do.call("rbind", by(t(x[,whKeep]), clFactor, function(z){apply(z, 2, median)}))
    rownames(medoids) <- levels(clFactor)
    nPerCluster <- table(clFactor)
    clusterD<-as.dendrogram(stats::hclust(dist(medoids)^2,members=nPerCluster,...))
    
    nSamples<-length(clNum)
		fullD<- .makeSampleDendro(x,clusterDendro=clusterD, cl=clNum,type=c("mat"), unassignedSamples=unassigned,sampleEdgeLength=0,  outbranchLength=1,calculateSample=calculateSample)
		
    return(list(samples=fullD,clusters=clusterD))
  })


# 	#wrapper so don't repeat code for class matrix and class dist
# .getUnassignedTree<-function(x,clNum,clusterD,type=c("dist","mat"),whKeep,nSamples,unassigned,calculateSample){
# 	type<-match.arg(type)
#   if(calculateSample){
# 		outlierHclust<-NULL
# 		if(length(whKeep) != nSamples && unassigned != "remove"){
# 		  outlierDat <- if(type=="mat") x[,-whKeep,drop=FALSE] else as.matrix(x)[-whKeep,-whKeep,drop=FALSE]
#
#       if(unassigned=="outgroup"| length(whKeep)==0){
#         #run hclust on unassigned
# 				if(ncol(outlierDat)>5 | length(whKeep)==0){
# 					outlierHclust <- if(type=="mat") stats::hclust(dist(t(outlierDat))^2) else  stats::hclust(as.dist(outlierDat))
# 				}
#       }
#       if(unassigned=="cluster"){
# 				stop("option unassigned='cluster' not yet operational")
# 				#cluster to existing clusters based on this dist matrix -- may need to update predict function for distances???
# 				#but then once get new cluster assignments, then should just update clNum that give to .makeSampleDendro
# 				unassigned<-"remove"
#       }
#     }
# 		if(length(whKeep)==0) fullD<-outlierHclust
# 		else fullD<-.makeSampleDendro(clusterD, clNum, unassignedSamples=unassigned,sampleEdgeLength=0, sampleNames=colnames(x),outbranchHclust=outlierHclust)
#   }
# 	else fullD<-NULL
# 	return(fullD)
# }

