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
#' @param whichCluster an integer index or character string that identifies 
#'   which cluster should be used to make the dendrogram. Default is 
#'   primaryCluster.
#' @param ... for makeDendrogram, if signature \code{matrix}, arguments passed 
#'   to hclust; if signature \code{ClusterExperiment} passed to the method for 
#'   signature \code{matrix}. For plotDendrogram, passed to
#'   \code{\link{plot.dendrogram}}.
#' @inheritParams clusterSingle
#' @inheritParams transform
#' @details The function returns two dendrograms (as a list if x is a matrix or
#' in the appropriate slots if x is ClusterExperiment). The cluster dendrogram
#' is created by applying \code{\link{hclust}} to the medoids of each cluster.
#' In the sample dendrogram the clusters are again clustered, but now the
#' samples are also part of the resulting dendrogram. This is done by giving
#' each sample the value of the medoid of its cluster.
#'
#' @details The argument \code{unassignedSamples} governs what is done with
#' unassigned samples (defined by a -1 cluster value). If unassigned=="cluster",
#' then the dendrogram is created by hclust of the expanded medoid data plus the
#' original unclustered observations. If \code{unassignedSamples} is "outgroup",
#' then all unassigned samples are put as an outgroup. If the \code{x} object is
#' a matrix, then \code{unassignedSamples} can also be "remove", to indicate
#' that samples with "-1" should be discarded. This is not a permitted option,
#' however, when \code{x} is a \code{ClusterExperiment} object, because it would
#' return a dendrogram with less samples than \code{NCOL(x)}, which is not
#' permitted for the \code{@dendro_samples} slot.
#'
#'
#' @return If x is a matrix, a list with two dendrograms, one in which the
#' leaves are clusters and one in which the leaves are samples. If x is a
#' ClusterExperiment object, the dendrograms are saved in the appropriate slots.
#'
#' @export
#'
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
#' plotDendrogram(hcl, leafType="samples",labelType="colorblock")
#'
#' @name makeDendrogram
#' @rdname makeDendrogram
setMethod(
  f = "makeDendrogram",
  signature = "ClusterExperiment",
  definition = function(x, whichCluster="primaryCluster",dimReduce=c("none", "PCA", "var","cv","mad"),
                        ndims=NA,ignoreUnassignedVar=TRUE,unassignedSamples=c("outgroup", "cluster"),...)
  {
    unassignedSamples<-match.arg(unassignedSamples)
    if(is.character(whichCluster)) whCl<-.TypeIntoIndices(x,whClusters=whichCluster) else whCl<-whichCluster
    if(length(whCl)!=1) stop("Invalid value for 'whichCluster'. Current value identifies ",length(whCl)," clusterings, but 'whichCluster' must identify only a single clustering.")
    if(!whCl %in% 1:nClusters(x)) stop("Invalid value for 'whichCluster'. Must be integer between 1 and ", nClusters(x))
#    browser()
    cl<-clusterMatrix(x)[,whCl]
	#cl<-convertClusterLegend(x,output="matrixNames")[,whCl]
    ########
    ##Transform the data
    ########
    dimReduce <- match.arg(dimReduce)
    if(length(ndims) > 1) {
      stop("makeDendrogram only handles one choice of dimensions.")
    }
    if(!is.na(ndims) & dimReduce=="none") {
      warning("specifying ndims has no effect if dimReduce==`none`")
    }
    origX <- assay(x)
    nPCADims <- ifelse(dimReduce=="PCA", ndims, NA)
    nVarDims <- ifelse(dimReduce=="var", ndims, NA)
    dimReduceCl<-if(ignoreUnassignedVar) cl else NULL #ifelse doesn't work with NULL
    transObj <- .transData(origX, nPCADims=nPCADims, nVarDims=nVarDims,
                           dimReduce=dimReduce, transFun=transformation(x),clustering=dimReduceCl)
    dat <- transObj$x
	if(is.null(dim(dat)) || NCOL(dat) != NCOL(origX)) {
      stop("Error in the internal transformation of x")
    }
    outlist <- makeDendrogram(x=dat, cluster=cl,unassignedSamples=unassignedSamples, ...)
    x@dendro_samples <- outlist$samples
    x@dendro_clusters <- outlist$clusters
    x@dendro_index<-whCl
	#browser()
	x@dendro_outbranch<- any(cl<0) & unassignedSamples=="outgroup"
    ch<-.checkDendrogram(x)
	if(!is.logical(ch)) stop(ch)
    return(x)
  })

#' @rdname makeDendrogram
#' @export
setMethod(
  f = "makeDendrogram",
  signature = "matrix",
  definition = function(x, cluster,
                        unassignedSamples=c("outgroup", "cluster", "remove"),
                        ...) {
   # browser()
    unassigned <- match.arg(unassignedSamples)
    cl <- cluster
    if(ncol(x) != length(cl)) {
      stop("cluster must be the same length as the number of samples")
    }
    if(is.null(colnames(x))) {
      colnames(x) <- as.character(1:ncol(x))
    }
	clNum<-.convertToNum(cl)
  
    # if(is.factor(cl)) {
    #   warning("cluster is a factor. Converting to numeric, which may not result in valid conversion")
    #   cl<-as.numeric(as.character(cl))
    # }
    #dat <- t(x) #make like was in old code

    #############
    # Cluster dendrogram
    #############
    whRm <- which(clNum >= 0) #remove -1, -2
    if(length(whRm) == 0) stop("all samples have clusterIds<0")
    if(length(unique(cl[whRm]))==1) stop("Only 1 cluster given. Can not make a dendrogram.")
	clFactor <- factor(cl[whRm])
    medoids <- do.call("rbind", by(t(x[,whRm]), clFactor, function(z){apply(z, 2, median)}))
    rownames(medoids) <- levels(clFactor)
    nPerCluster <- table(clFactor)
	#browser()
    clusterD<-as.dendrogram(stats::hclust(dist(medoids)^2,members=nPerCluster,...))
    
    #############
    # Samples dendrogram
    #############
    ####Attempt to make separate hclust and paste together.
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
    #     browser()
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

    #make fake data with just medoids as value per sample:
    dat <- t(x) #make like was in old code
    fakeData <- do.call("rbind", lapply(levels(clFactor), function(x){
        ind <- which(clFactor == x) #indices of those in this cluster
        med <- medoids[x,]
        mat <- matrix(rep(med, length(ind)), nrow=length(ind), byrow=TRUE)
        rownames(mat) <- rownames(dat[whRm,])[ind]
        return(mat)
    }))
    fakeData <- fakeData[rownames(dat[whRm,]),]
    if(length(whRm) != nrow(dat) && unassigned != "remove"){
        if(unassigned=="outgroup"){
            #hard to use merge and then get the indices to go back to the same ones
            #cheat and add large amount to the unassigned so that they don't cluster to
            outlierDat <- dat[-whRm,,drop=FALSE]
            maxAss <- max(dat[whRm,,drop=FALSE])
            outlierDat <- outlierDat + maxAss + 10e6
            fakeData <- rbind(fakeData, outlierDat)
            fakeData <- fakeData[rownames(dat),,drop=FALSE]
            # fullD<-as.dendrogram(stats::hclust(dist(fakeData)))
            # unassD<-as.dendrogram(stats::hclust(dist(dat[-whRm,])))
            # return(merge(fullD,unassD))
            
        }
        
        if(unassigned=="cluster"){
            #add remaining to fake data and let them cluster
            fakeData <- rbind(fakeData,dat[-whRm,,drop=FALSE])
            fakeData <- fakeData[rownames(dat),,drop=FALSE]
            #return(as.dendrogram(stats::hclust(dist(fakeData))))
        }
    }
#	browser()
	#make sure fakeData in same order as original data so order.dendrogram will work
	fakeData<-fakeData[na.omit(match(rownames(dat),rownames(fakeData))),]
	fullD <- as.dendrogram(stats::hclust(dist(fakeData)^2), ...)
    if(length(whRm) != nrow(dat) && unassigned == "outgroup"){
        #need to get rid of super long outgroup arm
        armLength <- max(attributes(fullD[[1]])$height,
                         attributes(fullD[[2]])$height)
        attributes(fullD)$height <- armLength + .1 * armLength
    }
  #   orderFullD<-dendextend::rotate(fullD,order=colnames(x)[orderSamples[,"index"]])
    return(list(samples=fullD,clusters=clusterD))
  })





