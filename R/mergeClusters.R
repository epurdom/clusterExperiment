#' @title Merge clusters based on dendrogram
#'
#' @description Takes an input of hierarchical clusterings of clusters and
#'   returns estimates of number of proportion of non-null and merges those
#'   below a certain cutoff.
#'
#' @aliases mergeClusters
#'
#' @param x data to perform the test on. It can be a matrix or a
#'   \code{\link{ClusterExperiment}}.
#' @param cl A numeric vector with cluster assignments to compare to. ``-1''
#'   indicates the sample was not assigned to a cluster.
#' @param dendro dendrogram providing hierarchical clustering of clusters in cl;
#'   The default for matrix (NULL) is to recalculate it with the given (x, cl)
#'   pair. If x is a \code{\link{ClusterExperiment}} object, the dendrogram in
#'   the slot \code{dendro_clusters} will be used. This means that
#'   \code{\link{makeDendrogram}} needs to be called before
#'   \code{mergeClusters}.
#' @param mergeMethod method for calculating proportion of non-null that will be
#'   used to merge clusters (if 'none', no merging will be done). See details
#'   for description of methods.
#' @param cutoff minimimum value required for NOT merging a cluster, i.e.
#'   two clusters with the proportion of DE below cutoff will be merged.
#'   Must be a value between 0, 1, where
#'   lower values will make it harder to merge clusters.
#' @param plotType what type of plotting of dendrogram. If 'all', then all the
#'   estimates of proportion non-null will be plotted; if 'mergeMethod', then
#'   only the value used in the merging is plotted for each node.
#' @param isCount logical as to whether input data is a count matrix. See details.
#' @param doPlot logical as to whether to plot the dendrogram (overrides
#'   \code{plotType} value). Mainly used for internal coding purposes.
#' @param ... for signature \code{matrix}, arguments passed to the
#'   \code{\link{plot.phylo}} function of \code{ade4} that plots the dendrogram.
#'   For signature \code{ClusterExperiment} arguments passed to the method for
#'   signature \code{matrix}.
#' @inheritParams clusterMany,matrix-method
#'
#' @details If  \code{isCount=TRUE}, and the input is a matrix,
#'   \code{log2(count + 1)} will be used for \code{\link{makeDendrogram}} and the
#'   original data with voom correction will be used in
#'   \code{\link{getBestFeatures}}). If input is
#'   \code{\link{ClusterExperiment}}, then setting \code{isCount=TRUE} also means
#'   that the log2(1+count) will be used as the transformation, like for
#'   the matrix case as well as the voom calculation, and will NOT use the
#'   transformation stored in the object. If FALSE, then transform(x) will be
#'   given to the input and will be used for both \code{makeDendrogram} and
#'   \code{getBestFeatures}, with no voom correction.
#' @details "JC" refers to the method of Ji and Cai (2007), and implementation
#'   of "JC" method is copied from code available on Jiashin Ji's website,
#'   December 16, 2015
#'   (http://www.stat.cmu.edu/~jiashun/Research/software/NullandProp/). "locfdr"
#'   refers to the method of Efron (2004) and is implemented in the package
#'   \code{\link{locfdr}}. "MB" refers to the method of Meinshausen and Buhlmann
#'   (2005) and is implemented in the package \code{\link{howmany}}. "adjP"
#'   refers to the proportion of genes that are found significant based on a FDR
#'   adjusted p-values (method "BH") and a cutoff of 0.05.
#'
#' @details If \code{mergeMethod} is not equal to 'none' then the plotting will
#'   indicate where the clusters will be merged (assuming \code{plotType} is not 'none').
#' @return If `x` is a matrix, it returns (invisibly) a list with elements
#'   \itemize{ \item{\code{clustering}}{ a vector of length equal to ncol(x)
#'   giving the integer-valued cluster ids for each sample. "-1" indicates the
#'   sample was not clustered.} \item{\code{oldClToNew}}{ A table of the old
#'   cluster labels to the new cluster labels.} \item{\code{propDE}}{ A table of
#'   the proportions that are DE on each node.}
#'   \item{\code{originalClusterDendro}}{ The dendrogram on which the merging
#'   was based (based on the original clustering).} }
#' @return If `x` is a \code{\link{ClusterExperiment}}, it returns a new
#'   \code{ClusterExperiment} object with an additional clustering based on the
#'   merging. This becomes the new primary clustering.
#' @examples
#' data(simData)
#'
#' #create a clustering, for 8 clusters (truth was 3)
#' cl<-clusterSingle(simData, clusterFunction="pam", subsample=FALSE,
#' sequential=FALSE, clusterDArgs=list(k=8))
#'
#' #make dendrogram
#' cl <- makeDendrogram(cl)
#'
#' #merge clusters with plotting. Note argument 'use.edge.length' can improve
#' #readability
#' merged <- mergeClusters(cl, plot=TRUE, plotType="all",
#' mergeMethod="adjP", use.edge.length=FALSE)
#'
#' #compare merged to original
#' table(primaryCluster(cl), primaryCluster(merged))
#' @export
#' @importFrom phylobase labels descendants ancestors getNode
#' @importClassesFrom phylobase phylo4
#' @importFrom graphics plot
#' @importFrom ape plot.phylo
#' @importFrom howmany howmany lowerbound
#' @importFrom locfdr locfdr
#' @rdname mergeClusters
setMethod(f = "mergeClusters",
          signature = signature(x = "matrix"),
          definition = function(x, cl, dendro=NULL,
                          mergeMethod=c("none", "adjP", "locfdr", "MB", "JC"),
                          plotType=c("none", "all", "mergeMethod","adjP", "locfdr", "MB", "JC"),
                          cutoff=0.1, doPlot=TRUE,
                          isCount=TRUE, ...) {
  if(is.factor(cl)){
    warning("cl is a factor. Converting to numeric, which may not result in valid conversion")
    cl <- .convertToNum(cl)
  }
  if(!is.null(dendro)){
    #check valid
    ncluster <- length(table(cl[cl>0]))
    if(nobs(dendro) != ncluster) {
      stop("Not a valid input dendrogram (not equal to the number of non -1 clusters in cl).")
    }
  }
  mergeMethod <- match.arg(mergeMethod)
  plotType <- match.arg(plotType)
  if(mergeMethod=="none" & plotType=="none") stop("mergeMethod and plotType both equal 'none'; nothing to be done.")
  if(plotType=="mergeMethod" & mergeMethod=="none") {
    stop("can only plot merge method values if one method is selected")
  }

  #get test-statistics for the contrasts corresponding to each node (and return all)
  sigTable <- getBestFeatures(x, cl, contrastType=c("Dendro"), dendro=dendro,
                               contrastAdj=c("All"),
                              number=nrow(x), p.value=1, isCount=isCount)

  #divide table into each node.
  whMethodCalculate<-if(!mergeMethod=="none") mergeMethod else c()
  if(plotType=="all") whMethodCalculate<-c("adjP", "locfdr", "MB", "JC")
  if(plotType%in% c("adjP", "locfdr", "MB", "JC")) whMethodCalculate<-unique(c(whMethodCalculate,plotType))
  sigByNode <- by(sigTable, sigTable$ContrastName, function(x) {
      mb <-if("MB" %in% whMethodCalculate)  .myTryFunc(pvalues=x$P.Value, FUN=.m1_MB) else NA
      locfdr <-if("locfdr" %in% whMethodCalculate)  .myTryFunc(tstats=x$t, FUN=.m1_locfdr) else NA
      jc <-if("JC" %in% whMethodCalculate)  .myTryFunc(tstats=x$t, FUN=.m1_JC) else NA
      adjP<-if("adjP" %in% whMethodCalculate)  .m1_adjP(x$adj) else NA
      return(c("adjP"=adjP, "locfdr"=locfdr, "MB"=mb,"JC"=jc))
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
    #browser()
    if(length(whToMerge)>0 && length(which(whToMerge)) > 0){
      nodesToMerge <- nodesBelowCutoff[whToMerge]

      #now find top ones
      whAnc <- sapply(nodesToMerge, function(node){
        anc <- phylobase::ancestors(phylo4Obj, node, type="all")
        return(!any(names(anc) %in% nodesToMerge))
      })
     # browser()
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

  #browser()
  nodePropTable <- do.call("rbind", sigByNode)
  annotTable <- data.frame("Node"=names(sigByNode),
                              "Contrast"=sigTable$Contrast[match(names(sigByNode), sigTable$ContrastName)])
    #add merge information:
   if(mergeMethod != "none" && length(whToMerge)>0 && length(which(whToMerge)) > 0){
     logicalMerge<-annotTable$Node %in%nodesToMerge

 } else logicalMerge<-rep(FALSE,length=nrow(annotTable))
  nodePropTable<-data.frame(annotTable,"Merged"=logicalMerge,nodePropTable)
  
  if(mergeMethod=="none"){
    newcl<-NULL #was just the original and nothing changed, so don't return something that makes it look like theres a new clustering
    oldClToNew<-NULL
      }
  else oldClToNew=table(Original=cl, New=newcl)
  out<-list(clustering=newcl, oldClToNew=oldClToNew,
                 propDE=nodePropTable, originalClusterDendro=dendro,mergeMethod=mergeMethod)
  if(doPlot) .plotMerge(dendro,mergeOutput=out,plotType=plotType,mergeMethod=mergeMethod,...)
  invisible(out)
}
)

.plotMerge<-function(dendro,mergeOutput,plotType,mergeMethod,clusterLegendMat=NULL,...){
    sigInfo<-mergeOutput$propDE
    whToMerge<-which(sigInfo$Merged)
    nodesToMerge<-sigInfo$Node[whToMerge]
    methods<-colnames(sigInfo[,-c(1:3)])
    if(plotType!="none"){
        # phylobase has bug in plotting! submitted to their github
        # move to ape package...
        phylo4Obj <- .makePhylobaseTree(dendro, "dendro")
        phyloObj <- as(phylo4Obj, "phylo")
        
        #####
        #convert names of internal nodes for plotting
        #####
        #match to order of tree
        m <- match(phyloObj$node, sigInfo$Node)
        edgeLty <- rep(1, nrow(phyloObj$edge))
        if(mergeMethod != "none" && length(whToMerge) > 0) {
            whMerge <- which(phyloObj$node.label %in% nodesToMerge) #which of nodes merged
            nodeNumbers <- (length(phyloObj$tip) + 1):max(phyloObj$edge)
            whEdge <- which(phyloObj$edge[,1] %in% nodeNumbers[whMerge])
            edgeLty[whEdge] <- 2
        }
        if(plotType == "mergeMethod"){
            if(!mergeMethod %in% methods) stop("mergeMethod not in methods of output")
            phyloObj$node.label <- as.character(sigInfo[m,mergeMethod])
        }
        if(plotType %in% c("all","adjP", "locfdr", "MB", "JC")) {
            meth<-if(plotType=="all") methods else methods[methods%in%plotType]
            phyloObj$node.label <- apply(sigInfo[,meth,drop=FALSE],1, function(x){
                whKp<-which(!is.na(x))
                paste(paste(meth[whKp], signif(x[whKp],2), sep=":"), collapse=",\n")})
            
        }
        #browser()
        ###Add color and name from the object.
        #browser()
        if(!is.null(clusterLegendMat)){
            m<-match(phyloObj$tip.label,clusterLegendMat[,"clusterIds"])
            if(any(is.na(m))) stop("clusterIds do not match dendrogram labels")
            phyloObj$tip.label<-clusterLegendMat[m,"name"]
            tip.color<-clusterLegendMat[m,"color"]
        }
        else tip.color<-"black"
        ape::plot.phylo(phyloObj, show.node=TRUE, edge.lty=edgeLty, tip.color=tip.color,...)
    }
}

#' @rdname mergeClusters
#' @export
#' @param clusterLabel a string used to describe the type of clustering. By
#'   default it is equal to "mergeClusters", to indicate that this clustering is
#'   the result of a call to mergeClusters.
setMethod(f = "mergeClusters",
          signature = signature(x = "ClusterExperiment"),
          definition = function(x, eraseOld=FALSE,isCount=FALSE,
                                mergeMethod="none",plotType="all",clusterLabel="mergeClusters",...) {

  if(is.null(x@dendro_clusters)) {
    stop("`makeDendrogram` needs to be called before `mergeClusters`")
  }
  else{
    cl<-clusterMatrix(x)[,x@dendro_index]
    note("Merging will be done on '",clusterLabels(x)[x@dendro_index],"', with clustering index",x@dendro_index)
  }
  if(isCount) note("If `isCount=TRUE` the data will be transformed with voom() rather than
with the transformation function in the slot `transformation`.
This makes sense only for counts.")
            #browser()
  outlist <- mergeClusters(x=if(!isCount) transform(x) else assay(x),
                           cl=cl,
                           dendro=x@dendro_clusters, plotType=plotType,doPlot=FALSE,
                           isCount=isCount,mergeMethod=mergeMethod, ...)
  if(plotType!="none"){
      .plotMerge(x@dendro_clusters,mergeOutput=outlist,plotType=plotType,mergeMethod=mergeMethod,clusterLegendMat=clusterLegend(x)[[x@dendro_index]])
  }
  
  if(mergeMethod!="none"){#only add a new cluster if there was a mergeMethod. otherwise, mergeClusters just returns original cluster!
    #----
    #add "m" to name of cluster
    #----
    newObj <- clusterExperiment(x, outlist$clustering,
                                transformation=transformation(x),
                                clusterTypes="mergeClusters")
    #add "m" to name of cluster
    newObj<-.addPrefixToClusterNames(newObj,prefix="m",whCluster=1)
    clusterLabels(newObj) <- clusterLabel
    ##Check if pipeline already ran previously and if so increase
    x<-.updateCurrentWorkflow(x,eraseOld,"mergeClusters")
    if(!is.null(x)) retval<-.addNewResult(newObj=newObj,oldObj=x)
    else retval<-.addBackSEInfo(newObj=newObj,oldObj=x)
    invisible(retval)
  }
  else{ #don't do anything, since there was no merging done.
    invisible(x)
  }
}
)

.myTryFunc<-function(FUN,...){
  x<-try(FUN(...),silent=TRUE)
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
  musigma<-try(.EstNull.func(tstats),silent=TRUE)
  .epsest.func(tstats,musigma$mu,musigma$s) #gives proportion of non-null
}
