#' @title Make hierarchy of set of clusters
#'
#' @description Makes a dendrogram of a set of clusters based on hclust on the
#' mediods of the cluster.
#'
#'
#' @param x data to define the mediods from. Matrix and
#' \code{\link{ClusterExperiment}} supported.
#' @param cluster A numeric vector with cluster assignments. If x is a
#' ClusterExperiment object, cluster is automatically the primaryCluster(x).
#' ``-1'' indicates the sample was not assigned to a cluster.
#' @param unassigned how to handle unassigned samples("-1") ; only relevant
#' for sample clustering. See details.
#' #' @param isCount logical. Whether the data are in counts, in which case the
#' default \code{transFun} argument is set as log(x+1). This is simply a
#' convenience to the user, and can be overridden by giving an explicit function
#' to \code{transFun}.
#' @param transFun function A function to use to transform the input data matrix
#' before clustering.
#' @param dimReduce character A character identifying what type of
#' dimensionality reduction to perform before clustering
#' @param ndims integer An integer identifying how many dimensions to reduce to
#' in the reduction specified by \code{dimReduce}
#' @param ... arguments passed to hclust
#' @details The function returns two dendrograms (as a list if x is a matrix) or
#' in the appropriate slots if x is ClusterExperiment. The cluster dendrogram
#' is created by applying \code{\link{hclust}} to the mediods of each cluster.
#' In the sample dendrogram the clusters are again clustered, but now the
#' samples are also part of the resulting dendrogram. This is done by giving
#' each sample the value of the mediod of its cluster.
#'
#' @details The argument \code{unassignedSamples} governs what is done with
#' unassigned samples (defined by a -1 cluster value). If unassigned=="cluster",
#' then the dendrogram is created by hclust of the expanded mediod data plus the
#' original unclustered observations. If \code{unassignedSamples} is "outgroup",
#' then all unassigned samples are put as an outgroup.
#'
#'
#' @return If x is a matrix, a list with two dendrograms, one in which the
#' leaves are clusters and one in which the leaves are samples. If x is a
#' ClusterExperiment object, the dendrograms are saved in the appropriate slots.

#' @examples
#' data(simData)
#' #create a clustering, for 8 clusters (truth was 3)
#' cl<-clusterSingle(simData,clusterFunction="pam",subsample=FALSE,
#' sequential=FALSE, clusterDArgs=list(k=8))$cl
#' #create dendrogram of clusters:
#' hcl<-clusterHclust(dat=simData,cl,full=FALSE)
#' plot(hcl)
#'
#' #create dendrogram for plotting with data in heatmap:
#' hclData<-clusterHclust(dat=simData,cl,full=TRUE)
#' plotHeatmap(cl,heatData=simCount,clusterData=hclData,colorScale=seqPal5,
#'	annCol=data.frame(PAM8=cl,Truth=trueCluster))
#' @rdname makeDendrogram
setMethod(
  f = "makeDendrogram",
  signature = "ClusterExperiment",
  definition = function(x, dimReduce=c("none", "PCA", "mostVar"),
                        ndims=NA,...)
  {
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
    nVarDims <- ifelse(dimReduce=="mostVar", ndims, NA)
    transObj <- .transData(origX, nPCADims=nPCADims, nVarDims=nVarDims,
                           dimReduce=dimReduce, transFun=transformation(x))
    dat <- transObj$x
    if(is.null(dim(dat)) || NCOL(dat) != NCOL(origX)) {
      stop("Error in the internal transformation of x")
    }
    outlist <- makeDendrogram(x=dat, cluster=primaryCluster(x), ...)
    x@dendro_samples <- outlist$samples
    x@dendro_clusters <- outlist$clusters
    validObject(x)
    return(x)
  })

#' @rdname makeDendrogram
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
    if(is.factor(cl)) {
      warning("cluster is a factor. Converting to numeric, which may not result in valid conversion")
      cl<-as.numeric(as.character(cl))
    }
    dat <- t(x) #make like was in old code

    #############
    whRm <- which(cl >= 0) #remove -1, -2
    if(length(whRm) == 0) stop("all samples have clusterIds<0")
    if(length(unique(cl))==1) stop("Only 1 cluster given. Can not make a dendrogram.")
    clFactor <- factor(cl[whRm])
    mediods <- do.call("rbind", by(dat[whRm,], clFactor,
                                   function(x){apply(x, 2, median)}))
    rownames(mediods) <- levels(clFactor)
    nPerCluster <- table(clFactor)

    ## leaves = samples
    #make fake data with just mediods as value per sample:
    fakeData <- do.call("rbind", lapply(levels(clFactor), function(x){
      ind <- which(clFactor == x) #indices of those in this cluster
      med <- mediods[x,]
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
          return(as.dendrogram(stats::hclust(dist(fakeData))))
        }
      }
      fullD <- as.dendrogram(stats::hclust(dist(fakeData)^2), ...)
      if(length(whRm) != nrow(dat) && unassigned == "outgroup"){
        #need to get rid of super long outgroup arm
        armLength <- max(attributes(fullD[[1]])$height,
                         attributes(fullD[[2]])$height)
        attributes(fullD)$height <- armLength + .1 * armLength
      }

      return(list(samples=fullD,
                  clusters=as.dendrogram(stats::hclust(dist(mediods)^2,
                                                members=nPerCluster,...))))
})

#' @rdname makeDendrogram
#' @inheritParams makeDendrogram
setMethod(
    f = "plotDendrogram",
    signature = "ClusterExperiment",
    definition = function(x,leaves=c("samples","clusters"),
                          main,...)
    {
        leaves<-match.arg(leaves)
        if(missing(main)) main<-ifelse(leaves=="samples","Dendrogram of samples", "Dendrogram of clusters")
        dend<- switch(leaves,"samples"=x@dendro_samples,"clusters"=x@dendro_clusters)
        plot(dend,main=main,...)
    })
