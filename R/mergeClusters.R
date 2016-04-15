#' @title Merge clusters based on dendrogram
#'
#' @description Takes an input of hierarchical clusterings of clusters and
#' returns estimates of number of proportion of non-null and merges those below
#' a certain cutoff.
#'
#'
#' @param x data to perform the test on. It can be a matrix or a
#' \code{\link{ClusterExperiment}}.
#' @param cl A numeric vector with cluster assignments to compare to clRef.
#' ``-1'' indicates the sample was not assigned to a cluster.
#' @param dendro If provided, is dendrogram providing hierarchical clustering of
#' clusters in cl; mainly useful to speed up calculations if
#' \code{\link{makeDendrogram}} has already been called. The default for matrix
#' (NULL) is to recalculate it with the given (x, cl) pair.
#' @param mergeMethod method for calculating proportion of non-null that will be
#' used to merge clusters (if 'none', no merging will be done). See details for
#' description of methods.
#' @param cutoff cutoff for merging clusters based on the proportion of
#' non-nulls in the comparison of the clusters (i.e. value between 0, 1, where
#' lower values will make it harder to merge clusters).
#' @param plotType what type of plotting of dendrogram. If 'all', then all the
#' estimates of proportion non-null will be plotted; if 'mergeMethod', then only
#' the value used in the merging is plotted for each node.
#' @param countData logical as to whether input data is a count matrix
#' (in which case log(count+1) will be used for \code{\link{makeDendrogram}}
#' and voom correction will be used in \code{\link{getBestGenes}}). Ignored if
#' input is \code{\link{ClusterExperiment}}.
#' @param ... arguments passed to the \code{\link{plot.phylo}} function of
#' \code{ade4} that plots the dendrogram.
#'
#' @details "JC" refers to the method of Ji and Cai (2007), and implementation
#' of "JC" method is copied from code available on Jiashin Ji's website,
#' December 16, 2015
#' (http://www.stat.cmu.edu/~jiashun/Research/software/NullandProp/). "locfdr"
#' refers to the method of Efron (2004) and is implemented in the package
#' \code{\link{locfdr}}. "MB" refers to the method of Meinshausen and
#' Buhlmann (2005) and is implemented in the package \code{\link{howmany}}.
#' "adjP" refers to the proportion of genes that are found significant based on
#' a FDR adjusted p-values (method "BH") and a cutoff of 0.05.
#'
#' @details If \code{mergeMethod} is not equal to 'none' then the plotting will
#' indicate the merged clusters (assuming \code{plotType} is not 'none').
#' @return If `x` is a matrix, it returns (invisibly) a list with elements
#' \itemize{
#' \item{\code{clustering}}{a vector of length equal to nrows(dat) giving the
#' integer-valued cluster ids for each sample. "-1" indicates the sample was not
#' clustered.}
#' \item{\code{oldClToNew}}{A table of the old cluster labels to the new cluster
#' labels}
#' \item{\code{propDE}}{A table of the proportions that are DE on each node}
#' \item{\code{originalClusterDendro}}{The dendrogram on which the merging was
#' based (based on the original clustering)}
#' }
#' @return If `x` is a \code{\link{ClusterExperiment}}, it returns a new
#' \code{ClusterExperiment} object with an additional clustering based on the
#' merging. This becomes the new primary clustering.
#' @examples
#' data(simData)
#'  #create a clustering, for 8 clusters (truth was 3)
#'  cl<-clusterSingle(simData,clusterFunction="pam",subsample=FALSE,
#'  sequential=FALSE, clusterDArgs=list(k=8))$cl
#' #merge clusters with plotting. Note argument 'use.edge.length' can improve readability
#' mergeResults<-mergeClusters(simData,cl=cl,plot=TRUE,plotType="all",
#' mergeMethod="adjP",use.edge.length=FALSE,countData=FALSE)
#' #compare merged to original on heatmap
#' hclData<-clusterHclust(dat=simData,cl,full=TRUE)
#' plotHeatmap(cl,heatData=simCount,clusterData=hclData,colorScale=seqPal5,
#'	annCol=data.frame(Original=cl,Merged=mergeResults$cl,Truth=trueCluster))
#' @export
#' @importFrom phylobase labels descendant ancestors getNode
#' @importClassesFrom ape phylo
#' @importClassesFrom phylobase phylo4
#' @importFrom ape plot.phylo
#' @rdname mergeClusters
setMethod(f = "mergeClusters",
          signature = signature(x = "matrix"),
          definition = function(x, cl, dendro=NULL,
                          mergeMethod=c("none", "adjP", "locfdr", "MB", "JC"),
                          cutoff=0.1, plotType=c("none", "all", "mergeMethod"),
                          countData=TRUE, ...) {

  if(is.factor(cl)){
    warning("cl is a factor. Converting to numeric, which may not result in valid conversion")
    cl <- as.numeric(as.character(cl))
  }

  if(!is.null(dendro)){
    #check valid
    ncluster <- length(table(cl[cl>0]))
    if(nobs(dendro) != ncluster) {
      warning("Not a valid input dendrogram (not equal to the number of non -1 clusters in cl). Will recompute dendrogram")
      dendro<-NULL
    }
  }

  if(is.null(dendro)) {
    dendro <- ifelse(countData,
                     makeDendrogram(x=log(x+1), cl, leaves="clusters"),
                     makeDendrogram(x=x, cl, leaves="clusters"))
  }

  mergeMethod <- match.arg(mergeMethod)
  plotType <- match.arg(plotType)
  if(plotType=="mergeValue" & mergeMethod=="none") {
    stop("can only plot merge method values if one method is selected")
  }

  #get test-statistics for the contrasts corresponding to each node (and return all)
  sigTable <- getBestGenes(x, cl, type=c("Dendro"), dendro=dendro,
                           returnType=c("Table"), contrastAdj=c("All"),
                           number=nrow(x), p.value=1, voomCorrection=countData)

  #divide table into each node.
  sigByNode <- by(sigTable, sigTable$ContrastName, function(x) {
    mb <- .myTryFunc(pvalues=x$P.Value, FUN=.m1_MB)
    locfdr <- .myTryFunc(tstats=x$t, FUN=.m1_locfdr)
    jc <- .myTryFunc(tstats=x$t, FUN=.m1_JC)
    return(c("adjP"=.m1_adjP(x$adj), "locfdr"=locfdr, "MB"=mb,"JC"=jc))
  })
  newcl <- cl
  phylo4Obj <- .makePhylobaseTree(dendro, "dendro")

  if(mergeMethod != "none"){
    #go up tree and merge clusters
    valsPerNode <- sapply(sigByNode, function(x) {signif(x[[mergeMethod]], 2)})
    nodesBelowCutoff <- names(valsPerNode)[which(valsPerNode<cutoff)] #names of nodes below cutoff

    #find nodes where *all* descendants are below cutoff
    allTipNames <- phylobase::labels(phylo4Obj)[phylobase::getNode(phylo4Obj, type=c("tip"))]
    whToMerge <- sapply(nodesBelowCutoff,function(node){
      desc <- phylobase::descendants(phylo4Obj, node, type = c("all"))
      return(all(names(desc) %in% nodesBelowCutoff | names(desc) %in% allTipNames))
    })
    if(length(whToMerge) > 0){
      nodesToMerge <- nodesBelowCutoff[whToMerge]

      #now find top ones
      whAnc <- sapply(nodesToMerge, function(node){
        anc <- phylobase::ancestors(phylo4Obj, node, type="all")
        return(!any(names(anc) %in% nodesToMerge))
      })
      nodesAtTop <- nodesToMerge[whAnc]

      #make new clusters
      temp <- lapply(nodesAtTop, function(node){
        tips <- phylobase::descendants(phylo4Obj, node, type="tips")
        if(any(!names(tips) %in% as.character(cl))) {
          stop("coding error-- tips don't match values of cl")
        }
        newcl[cl%in% tips] <<- as.numeric(tips[[1]])
      })
      #make consecutive integers
      newcl <- as.numeric(factor(newcl, levels=unique(cl[cl>0])))
      #deal with -1/-2
      newcl[is.na(newcl)] <- cl[is.na(newcl)]
    }
  }

  if(plotType!="none"){
    # phylobase has bug in plotting! submitted to their github
    # move to ape package...
    phyloObj <- as(phylo4Obj, "phylo")

    #####
    #convert names of internal nodes for plotting
    #####
    allInternal <- phyloObj$node
    #match to order of tree
    m <- match(allInternal, names(sigByNode))
    edgeLty <- rep(1, nrow(phyloObj$edge))
    if(mergeMethod != "none" & length(whToMerge) > 0) {
      whMerge <- which(phyloObj$node.label %in% nodesToMerge) #which of nodes merged
      nodeNumbers <- (length(phyloObj$tip) + 1):max(phyloObj$edge)
      whEdge <- which(phyloObj$edge[,1] %in% nodeNumbers[whMerge])
      edgeLty[whEdge] <- 2
    }
    if(plotType == "mergeMethod"){
      phyloObj$node.label <- as.character(valsPerNode)
    }
    if(plotType == "all") {
      phyloObj$node.label <- sapply(sigByNode[m], function(x){
        paste(paste(names(x), signif(x,2), sep=":"), collapse=",\n")})
    }

    ape::plot.phylo(phyloObj, show.node=TRUE, edge.lty=edgeLty, ...)
  }
  nodePropTable <- do.call("rbind", sigByNode)
  nodePropTable <- data.frame("Node"=names(sigByNode),
                              "Contrast"=sigTable$Contrast[match(names(sigByNode), sigTable$ContrastName)],
                              nodePropTable)

  invisible(list(clustering=newcl, oldClToNew=table(Original=cl, New=newcl),
                 propDE=nodePropTable, originalClusterDendro=dendro))
}
)

#' @rdname mergeClusters
setMethod(f = "mergeClusters",
          signature = signature(x = "ClusterExperiment"),
          definition = function(x, ...) {

  if(is.null(x@dendro_clusters)) {
    stop("`makeDendrogram` needs to be called before `mergeClusters`")
  }
  outlist <- mergeClusters(x=transform(x), cl=primaryCluster(x),
                           dendro=x@dendro_clusters,
                           countData=FALSE, ...)

  #add "m" to name of cluster
  idx <- which(outlist$clustering>0)
  cl <- as.numeric(as.factor(outlist$clustering[idx]))
  cl <- paste("m", cl, sep="")
  cl_labels <- as.character(outlist$clustering)
  cl_labels[idx] <- cl

  newObj <- clusterExperiment(x, cl_labels,
                              transformation=transformation(x),
                              clusterType="mergeClusters")
  clusterLabels(newObj) <- "mergeClusters"

  retval <- addClusters(newObj, x)
  retval@dendro_samples <- x@dendro_samples
  return(retval)
}
)

.myTryFunc<-function(FUN,...){
  x<-try(FUN(...))
  if(!inherits(x, "try-error")) return(x)
  else return(NA)
}

#functions for estimating m1/m, the proportion of non-null
.m1_MB<-function(pvalues){
  nCorrect<-max(howmany::lowerbound(howmany::howmany(pvalues))) #the most you can call correctly
  return(nCorrect/length(pvalues))
}
.m1_adjP<-function(adjP){
  sum(adjP<0.05)/length(adjP)
}
.m1_locfdr<-function(tstats){
  locfdrResults<-locfdr::locfdr(tstats,plot=0)#ignore issue of df of t-statistic -- topTable doesn't give it out, and with large samples won't matter.
  p0<-locfdrResults$fp0["mlest","p0"] #estimate proportion null; ignore estimate of variability for now
  return(1-p0)
}
.m1_JC<-function(tstats){
  #copied code from Jianshin's website
  musigma<-try(.EstNull.func(tstats))
  .epsest.func(tstats,musigma$mu,musigma$s) #gives proportion of non-null
}
