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
#'   erased by a call to mergeDendrogram.
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
				cluster=cl,
				unassignedSamples=unassignedSamples),
				passedArgs))
    }
    else{
      if(is.null(x@coClustering)) stop("Cannot choose 'coCluster' if 'coClustering' slot is empty. Run makeConsensus before running 'makeDendrogram' or choose another option for 'reduceMethod'")
      if(is.null(dimnames(x@coClustering))) stop("This ClusterExperiment object was made with an old version of clusterExperiment and did not give dimnames to the coClustering slot.")
     outlist<-do.call("makeDendrogram",c(list(
			  x=as.dist(1-x@coClustering),
				cluster=cl,
				unassignedSamples=unassignedSamples),
				passedArgs)) 
    }
    x@dendro_samples <- outlist$samples
    x@dendro_clusters <- outlist$clusters
    x@dendro_index<-whCl

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
                        ...) {
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
    
    #make fake dist with just medoids as value per sample (only for whKeep samples):
    fakeData <- do.call("rbind",lapply(levels(clFactor), function(z){
      ind <- which(clFactor == z) #indices of those in this cluster
      med <- medoids[z,] #vector of mediod with other clusters
      medExp<-rep(med,times=nPerCluster)
      mat <- matrix(medExp, nrow=length(ind), ncol=sum(nPerCluster) ,byrow=TRUE)
      rownames(mat) <- rownames(as.matrix(x)[whKeep,])[ind]
      colnames(mat) <- colnames(as.matrix(x)[,whKeep])
      return(mat)
    }))
    if(length(whKeep) != nSamples && unassigned != "remove"){
      #need to add the unassigned samples into fakeData
      outlierDat <- as.matrix(x)[-whKeep,-whKeep,drop=FALSE]
      if(unassigned=="outgroup"){
        #hard to use merge and then get the indices to go back to the same ones
        #cheat and add large amount to the unassigned so that they don't cluster to
        maxAss <- max(c(max(fakeData),max(outlierDat)))
        offDiag<-matrix(maxAss+10^6,nrow=nrow(fakeData),ncol=ncol(outlierDat))
        fakeData <- rbind(cbind(fakeData, offDiag),cbind(t(offDiag),outlierDat))
      }
      
      if(unassigned=="cluster"){
        #add remaining distances to fake data and let them cluster
        offDiag <- as.matrix(x)[whKeep,-whKeep,drop=FALSE]
        #put in order of fake data, then outliers
        offDiag<-offDiag[row.names(fakeData),,drop=FALSE]
        
        fakeData <- rbind(cbind(fakeData, offDiag),cbind(t(offDiag),outlierDat))
      }
    }
    #make sure fakeData in same order as original data so order.dendrogram will work
    sampleNames<-attributes(x)$Labels
    m<-na.omit(match(sampleNames,rownames(fakeData)))
    fakeData<-fakeData[m,m]
    fullD <- as.dendrogram(stats::hclust(as.dist(fakeData)), ...)
    if(length(whKeep) != nSamples && unassigned == "outgroup"){
      #need to get rid of super long outgroup arm
      armLength <- max(attributes(fullD[[1]])$height,
                       attributes(fullD[[2]])$height)
      attributes(fullD)$height <- armLength + .1 * armLength
    }
    #   orderFullD<-dendextend::rotate(fullD,order=colnames(x)[orderSamples[,"index"]])
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
                        ...) {
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
    
    
    #############
    # Samples dendrogram
    #############
    
    #make fake data with just medoids as value per sample:
    dat <- t(x) #make like was in old code
    fakeData <- do.call("rbind", lapply(levels(clFactor), function(z){
      ind <- which(clFactor == z) #indices of those in this cluster
      med <- medoids[z,]
      mat <- matrix(rep(med, length(ind)), nrow=length(ind), byrow=TRUE)
      rownames(mat) <- rownames(dat[whKeep,])[ind]
      return(mat)
    }))
    fakeData <- fakeData[rownames(dat[whKeep,]),] #why do I need this??
    if(length(whKeep) != nrow(dat) && unassigned != "remove"){
      if(unassigned=="outgroup"){
        #hard to use merge and then get the indices to go back to the same ones
        #cheat and add large amount to the unassigned so that they don't cluster to
        
        outlierDat <- dat[-whKeep,,drop=FALSE]
        maxAss <- max(dat[whKeep,,drop=FALSE])
        outlierDat <- outlierDat + maxAss + 10e6
        ############
        ###This a workaround which will hopefully be dealt with in future hdf5:
        ############
        if(inherits(fakeData,"DelayedMatrix")|| inherits(outlierDat,"DelayedMatrix")) fakeData<-rbind(DelayedArray::DelayedArray(fakeData), DelayedArray(outlierDat))
        else fakeData <- rbind(fakeData, outlierDat)
        fakeData <- fakeData[rownames(dat),,drop=FALSE]
        # fullD<-as.dendrogram(stats::hclust(dist(fakeData)))
        # unassD<-as.dendrogram(stats::hclust(dist(dat[-whKeep,])))
        # return(merge(fullD,unassD))
        
      }
      
      if(unassigned=="cluster"){
        #add remaining to fake data and let them cluster
        fakeData <- rbind(fakeData,dat[-whKeep,,drop=FALSE])
        fakeData <- fakeData[rownames(dat),,drop=FALSE]
        #return(as.dendrogram(stats::hclust(dist(fakeData))))
      }
    }
    #make sure fakeData in same order as original data so order.dendrogram will work
    fakeData<-fakeData[na.omit(match(rownames(dat),rownames(fakeData))),]
    fullD <- as.dendrogram(stats::hclust(dist(fakeData)^2), ...)
    if(length(whKeep) != nrow(dat) && unassigned == "outgroup"){
      #need to get rid of super long outgroup arm
      armLength <- max(attributes(fullD[[1]])$height,
                       attributes(fullD[[2]])$height)
      attributes(fullD)$height <- armLength + .1 * armLength
    }
    #   orderFullD<-dendextend::rotate(fullD,order=colnames(x)[orderSamples[,"index"]])
    return(list(samples=fullD,clusters=clusterD))
  })


    ####Past Attempt for samples dendro: make separate hclust and paste together.
    #     #For each cluster, run hclust and convert to phylobase tree
    #     dendPerCl<-by(t(x[,whRm]),clFactor,function(z){.makePhylobaseTree(as.dendrogram(stats::hclust(dist(z)^2)),type="dendro")})
    #     clusterDPhybase<-.makePhylobaseTree(clusterD,type="dendro")
    #     clusterEdge<-phylobase::edges(clusterDPhybase)
    #     #phylobase::rootNode(dendPerCl[[1]]) #gives root node, which needs to be merged into tip of main clusterDend
    #
    #     currMax<-max(clusterEdge)
    #     edgeList<-list()
    #     nn<-phylobase::getNode(clusterDPhybase)
    #     nodeMat<-data.frame("NodeId"=nn,"NodeName"=names(nn),"newNodeId"=nn,"tree"="ClusterTree")
    #     nodeMatCluster<-nodeMat #so have copy not changed in for loop
    #
    #     ###For loop changes node numbers so have proper edge list:
    #     for(kk in 1:length(dendPerCl)){
    #         temp<-phylobase::edges(dendPerCl[[kk]])
    #         tempNew<-temp+currMax+1
    #         nn<-phylobase::getNode(dendPerCl[[kk]])
    #         tempNode<-data.frame("NodeId"=nn,"NodeName"=names(nn),"newNodeId"=nn+currMax+1,"tree"=paste("TipTree",kk,sep=""))
    #
    #         #root node:
    #         rootNode<-phylobase::rootNode(dendPerCl[[kk]])
    #
    #         #remove existing root:
    #         whRootEdge<-which(temp[,"ancestor"]==0)
    #         whRootNodeId<-tempNew[whRootEdge,"descendant"] #new node id. Needs to be replaced with that of cluster
    #         tempNew<-tempNew[-whRootEdge,]
    #
    #         #match to node in clusterDPhybase, and add edge
    #         wh<-which(nodeMatCluster[,"NodeName"]==levels(clFactor)[kk])
    #         clNodeIdBigD<-nodeMatCluster[wh,"NodeId"] #will replace existing root node with this one
    #         tempNew<-apply(tempNew,2,function(x){x[which(x==whRootNodeId)]<-clNodeIdBigD; return(x)})
    #         tempNode$newNodeId[which(tempNode$newNodeId==whRootNodeId)]<-clNodeIdBigD
    #
    #         #update
    #         currMax<-max(tempNew)
    #         edgeList[[kk]]<-tempNew
    #         nodeMat<-rbind(nodeMat,tempNode)
    #     }
    #     whTips<-setdiff(grep("TipTree",nodeMat$tree),grep("Node",nodeMat$NodeName))
    # #make tree for phylo (can't get phylobase to accept my tree and crashes...)
    #     nodeMat$finalNodeId<-NA
    #     nodeMat$finalNodeId[whTips]<-1:length(whTips)
    #     notTips<-c(1:nrow(nodeMat))[-whTips]
    #     nodeMat$finalNodeId[-whTips]<-seq(from=length(whTips)+2,length=length(notTips),by=1) #internal nodes, not root
    #     combineEdge<-rbind(clusterEdge,do.call("rbind",edgeList))
    #     combineEdgeFinal<-apply(combineEdge,2,function(x){
    #         m<-match(x,nodeMat[,"newNodeId"])
    #         return(nodeMat[m,"finalNodeId"])
    #     })
    #     #combineEdgeFinal[combineEdgeFinal[,1]==0,1]<-0
    #     combineEdgeFinal[combineEdge[,1]==0,1]<-length(whTips)+1
    #     phyloObj<-list(edge=combineEdgeFinal,tip.label=as.character(1:))
    #
    #     #phylo4(combineEdgeFinal)


###Past attempt to just reordering, didn't work...
#     orderCluster<-order.dendrogram(clusterD)
#     #For each cluster, run hclust and get order
#     indByCl<-tapply(whRm,clFactor,function(ind){ind})
#     orderPerCl<-lapply(indByCl,function(ind){
#         d<-as.dendrogram(stats::hclust(dist(t(x[,ind]))^2))
#         o<-order.dendrogram(d)
#         return(cbind(orderWithinCluster=o,index=ind))
#     })
#     names(orderPerCl)<-levels(clFactor)
#     orderPerCl<-orderPerCl[orderCluster]
#     ## -1 group
#     if(length(whRm)!=ncol(x)){
#         ind<-(1:ncol(x))[-whRm]
#         order1<-cbind(orderWithinCluster=order.dendrogram(as.dendrogram(stats::hclust(dist(t(x[,-whRm])^2)))),index=ind)
#         orderPerCl<-c(orderPerCl,list(order1))
#     }
#     currMax<-0
#     for(kk in 1:length(orderPerCl)){
#         orderPerCl[[kk]]<-cbind( orderPerCl[[kk]],"externalOrder"=orderPerCl[[kk]][,"orderWithinCluster"]+currMax,group=kk)
#         currMax<-max(orderPerCl[[kk]])
#     }
#     orderSamples<-do.call("rbind",orderPerCl)
#     orderSamples<-orderSamples[order(orderSamples[,"externalOrder"]),]



