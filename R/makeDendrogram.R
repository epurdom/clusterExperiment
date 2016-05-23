#' @title Make hierarchy of set of clusters
#'
#' @description Makes a dendrogram of a set of clusters based on hclust on the
#' medoids of the cluster.
#'
#' @aliases makeDendrogram
#'
#' @param x data to define the medoids from. Matrix and
#' \code{\link{ClusterExperiment}} supported.
#' @param cluster A numeric vector with cluster assignments. If x is a
#' ClusterExperiment object, cluster is automatically the primaryCluster(x).
#' ``-1'' indicates the sample was not assigned to a cluster.
#' @param unassignedSamples how to handle unassigned samples("-1") ; only relevant
#' for sample clustering. See details.
#' @param dimReduce character A character identifying what type of
#' dimensionality reduction to perform before clustering.
#' @param ndims integer An integer identifying how many dimensions to reduce to
#' in the reduction specified by \code{dimReduce}.
#' @param whichCluster an integer index or character string that identifies
#'   which cluster should be used to make the dendrogram. Default is
#'   primaryCluster.
#' @param ... for makeDendrogram, if signature \code{matrix}, arguments passed
#'   to hclust; if signature \code{ClusterExperiment} passed to the method for
#'   signature \code{matrix}. For plotDendrogram, passed to \code{plot}.
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
#' cl <- clusterSingle(simData, clusterFunction="pam", subsample=FALSE,
#' sequential=FALSE, clusterDArgs=list(k=8))
#'
#' #create dendrogram of clusters:
#' hcl <- makeDendrogram(cl)
#' plotDendrogram(hcl)
#' plotDendrogram(hcl, leaves="samples")
#'
#' @rdname makeDendrogram
setMethod(
  f = "makeDendrogram",
  signature = "ClusterExperiment",
  definition = function(x, whichCluster="primaryCluster",dimReduce=c("none", "PCA", "var"),
                        ndims=NA,unassignedSamples=c("outgroup", "cluster"),...)
  {
    unassignedSamples<-match.arg(unassignedSamples)
    if(is.character(whichCluster)) whCl<-.TypeIntoIndices(x,whClusters=whichCluster) else whCl<-whichCluster
    if(length(whCl)!=1) stop("Invalid value for 'whichCluster'. Current value identifies ",length(whCl)," clusterings, but 'whichCluster' must identify only a single clustering.")
    if(!whCl %in% 1:nClusters(x)) stop("Invalid value for 'whichCluster'. Must be integer between 1 and ", nClusters(x))
#    browser()
    cl<-clusterMatrix(x)[,whCl]
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
    transObj <- .transData(origX, nPCADims=nPCADims, nVarDims=nVarDims,
                           dimReduce=dimReduce, transFun=transformation(x))
    dat <- transObj$x
    if(is.null(dim(dat)) || NCOL(dat) != NCOL(origX)) {
      stop("Error in the internal transformation of x")
    }
    outlist <- makeDendrogram(x=dat, cluster=cl,unassignedSamples=unassignedSamples, ...)
    x@dendro_samples <- outlist$samples
    x@dendro_clusters <- outlist$clusters
    x@dendro_index<-whCl
    validObject(x)
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
    if(is.factor(cl)) {
      warning("cluster is a factor. Converting to numeric, which may not result in valid conversion")
      cl<-as.numeric(as.character(cl))
    }
    dat <- t(x) #make like was in old code

    #############
    whRm <- which(cl >= 0) #remove -1, -2
    if(length(whRm) == 0) stop("all samples have clusterIds<0")
    if(length(unique(cl[whRm]))==1) stop("Only 1 cluster given. Can not make a dendrogram.")
    clFactor <- factor(cl[whRm])
    medoids <- do.call("rbind", by(dat[whRm,], clFactor,
                                   function(x){apply(x, 2, median)}))
    rownames(medoids) <- levels(clFactor)
    nPerCluster <- table(clFactor)

    ## leaves = samples
    #make fake data with just medoids as value per sample:
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
                  clusters=as.dendrogram(stats::hclust(dist(medoids)^2,
                                                members=nPerCluster,...))))
})

#' @rdname makeDendrogram
#' @export
#' @param leaves if "samples" the dendrogram has one leaf per sample, otherwise
#'   it has one per cluster.
#' @param main passed to the \code{plot} function.
#' @param sub passed to the \code{plot} function.
#' @param clusterNames logical. If \code{leaves="clusters"}, then clusters will
#'   be identified with their 'name' value in legend; otherwise the 'clusterIds'
#'   value will be used.
#' @aliases plotDendrogram
#' @details If \code{leaves="clusters"}, the plotting function will work best if
#'   the clusters in the dendrogram correspond to the primary cluster. This is
#'   because the function colors the cluster labels based on the colors of the
#'   clusterIds of the primaryCluster
setMethod(
  f = "plotDendrogram",
  signature = "ClusterExperiment",
  definition = function(x,leaves=c("clusters","samples" ), clusterNames=TRUE,
                        main,sub,...)
  {
    leaves<-match.arg(leaves)
    if(missing(main)) main<-ifelse(leaves=="samples","Dendrogram of samples", "Dendrogram of clusters")
    if(is.null(x@dendro_samples) || is.null(x@dendro_clusters)) stop("No dendrogram is found for this ClusterExperiment Object. Run makeDendrogram first.")
    if(missing(sub)) sub<-paste("Dendrogram made with '",clusterLabels(x)[x@dendro_index],"', cluster index ",x@dendro_index,sep="")
    dend<- switch(leaves,"samples"=x@dendro_samples,"clusters"=x@dendro_clusters)
    labs<-labels(dend)
    if(leaves=="clusters"){
      leg<-clusterLegend(x)[[x@dendro_index]]
      m<-match(labs,leg[,"clusterIds"])
      if(any(is.na(m))) warning("Dendrogram labels do not all match clusterIds of primaryCluster. Dendrogram was not created with current primary cluster, so cannot retreive cluster name or color")
      else{
        #function to change to labels and colors of a node:
        reLabel <- function(n) {
          if(is.leaf(n)) {
            a <- attributes(n)
            m<-match(a$label,leg[,"clusterIds"])
            if(clusterNames) attr(n, "label") <- leg[m,"name"]           #  change the node label
            attr(n, "nodePar") <- c(a$nodePar, list(lab.col = leg[m,"color"],col=leg[m,"color"],pch=19)) #   change the node color
          }
          return(n)
        }
        dend <- dendrapply(dend, reLabel)
      }
    }
    plot(dend,main=main,sub=sub,...)
  })
