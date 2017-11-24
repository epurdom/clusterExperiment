#' @title Function for finding best features associated with clusters
#' @description Calls limma on input data to determine features most associated
#'   with found clusters (based on an F-statistic, pairwise comparisons, or
#'   following a tree that clusters the clusters).
#' @aliases getBestFeatures
#' @param x data for the test. Can be a numeric matrix or a
#'   \code{\link{ClusterExperiment}}.
#' @param cluster A numeric vector with cluster assignments.
#'   ``-1'' indicates the sample was not assigned to a cluster.
#' @param contrastType What type of test to do. `F' gives the omnibus
#'   F-statistic, `Dendro' traverses the given dendrogram and does contrasts of
#'   the samples in each side,  `Pairs' does pair-wise contrasts based on the
#'   pairs given in pairMat (if pairMat=NULL, does all pairwise), and
#'   `OneAgainstAll' compares each cluster to the average of all others. Passed
#'   to \code{\link{clusterContrasts}}
#' @param contrastAdj What type of FDR correction to do for contrasts tests
#'   (i.e. if contrastType='Dendro' or 'Pairs').
#' @param isCount logical as to whether input data is count data, in which
#'   case to perform voom correction to data. See details.
#' @param ... options to pass to \code{\link{topTable}} or
#'   \code{\link[limma]{topTableF}} (see \code{\link[limma]{limma}} package)
#' @param normalize.method character value, passed to \code{\link[limma]{voom}} in
#'   \code{\link[limma]{limma}} package. Only used if \code{countData=TRUE}.
#'   Note that the default value is set to "none", which is not the
#'   default value of \code{\link{voom}}.
#' @inheritParams clusterContrasts,ClusterExperiment-method
#' @details getBestFeatures returns the top ranked features corresponding to a
#'   cluster assignment. It uses limma to fit the models, and limma's functions
#'   \code{\link[limma]{topTable}} or \code{\link[limma]{topTableF}} to find the
#'   best features. See the options of these functions to put better control on
#'   what gets returned (e.g. only if significant, only if log-fc is above a
#'   certain amount, etc.). In particular, set `number=` to define how many
#'   significant features to return (where number is per contrast for the
#'   `Pairs` or `Dendro` option)
#'
#' @details When `contrastType` argument implies that the best features should be found
#'   via contrasts (i.e. 'contrastType' is `Pairs` or `Dendro`), then then `contrastAdj`
#'   determines the type of multiple testing correction to perform.
#'   `PerContrast` does FDR correction for each set of contrasts, and does not
#'   guarantee control across all the different contrasts (so probably not the
#'   preferred method). `All` calculates the corrected p-values based on FDR
#'   correction of all of the contrasts tested. `AfterF` controls the FDR based
#'   on a hierarchical scheme that only tests the contrasts in those genes where
#'   the omnibus F statistic is significant. If the user selects `AfterF`, the
#'   user must also supply an option `p.value` to have any effect, and then only
#'   those significant at that p.value level will be returned. Note that
#'   currently the correction for `AfterF` is not guaranteed to control the FDR;
#'   improvements will be added in the future.
#'
#' @details  Note that the default option for \code{\link[limma]{topTable}} is
#'   to not filter based on adjusted p-values (\code{p.value = 1}) and return
#'   only the top 10 most significant (\code{number = 10}) -- these are options
#'   the user can change (these arguments are passed via the \code{...} in
#'   \code{getBestFeatures}). In particular, it only makes sense to set
#'   \code{requireF = TRUE} if \code{p.value} is meaningful (e.g. 0.1 or 0.05);
#'   the default value of \code{p.value = 1} will not result in any effect on
#'   the adjusted p-value otherwise.
#' @details  \code{isCount} triggers whether the "voom" correction will be
#'   performed in \code{limma}. If the input data is a matrix is counts (or a
#'   `ClusterExperiment` object with counts as the primary data before
#'   transformation) this should be set to TRUE and they will be log-transformed
#'   internally by voom for the differential expression analysis in a way that
#'   accounts for the difference in the mean-variance relationships. Otherwise,
#'   dat should be on the correct (log) scale for differential expression
#'   analysis without a need a variance stabilization (e.g. microarray data).
#'   Currently the default is set to FALSE, simply because the isCount has not
#'   been heavily tested. If the But TRUE with \code{x} being counts really
#'   should be the default for RNA-Seq data. If the input data is a
#'   `ClusterExperiment` object, setting `isCount=TRUE` will cause the program
#'   to ignore the internally stored transformation function and instead use
#'   voom with log2(x+0.5). Alternatively, \code{isCount=FALSE} for a
#'   \code{ClusterExperiment} object will cause the DE to be performed with \code{limma}
#'   after transforming the data with the stored transformation. Although some
#'   writing about "voom" seem to suggest that it would be appropriate for
#'   arbitrary transformations, the authors have cautioned against using it for
#'   anything other than count data on mailing lists. For this reason we are not
#'   implementing it for arbitrary transformations at this time (e.g.
#'   log(FPKM+epsilon) transformations).
#'
#' @return A \code{data.frame} in the same format as
#'   \code{\link[limma]{topTable}}, except for the following additional or
#'   changed columns:
#' \itemize{
#'
#' \item{\code{Feature}}{ This is the column called 'ProbeID' by
#' \code{\link{topTable}}}
#'
#' \item{\code{IndexInOriginal}}{ Gives the index of the feature to the original
#' input dataset, \code{x}}
#'
#' \item{\code{Contrast}}{ The contrast that the results corresponds to (if
#' applicable, depends on \code{contrastType} argument)}
#'
#' \item{\code{ContrastName}}{ The name of the contrast that the results
#' corresponds to. For dendrogram searches, this will be the node of the tree of
#' the dendrogram.}
#' }
#'
#' @references Ritchie, ME, Phipson, B, Wu, D, Hu, Y, Law, CW, Shi, W, and Smyth, GK (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research 43, e47. http://nar.oxfordjournals.org/content/43/7/e47
#' @references Law, CW, Chen, Y, Shi, W, and Smyth, GK (2014). Voom: precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biology 15, R29. http://genomebiology.com/2014/15/2/R29
#' @references Smyth, G. K. (2004). Linear models and empirical Bayes methods for assessing differential expression in microarray experiments. Statistical Applications in Genetics and Molecular Biology, Volume 3, Article 3. http://www.statsci.org/smyth/pubs/ebayes.pdf
#' @examples
#' data(simData)
#'
#' #create a clustering, for 8 clusters (truth was 4)
#' cl <- clusterSingle(simData, subsample=FALSE,
#' sequential=FALSE, mainClusterArgs=list(clusterFunction="pam", clusterArgs=list(k=8)))
#'
#' #basic F test, return all, even if not significant:
#' testF <- getBestFeatures(cl, contrastType="F", number=nrow(simData),
#' isCount=FALSE)
#'
#' #Do all pairwise, only return significant, try different adjustments:
#' pairsPerC <- getBestFeatures(cl, contrastType="Pairs", contrastAdj="PerContrast",
#' p.value=0.05, isCount=FALSE)
#' pairsAfterF <- getBestFeatures(cl, contrastType="Pairs", contrastAdj="AfterF",
#' p.value=0.05, isCount=FALSE)
#' pairsAll <- getBestFeatures(cl, contrastType="Pairs", contrastAdj="All",
#' p.value=0.05, isCount=FALSE)
#'
#' #not useful for this silly example, but could look at overlap with Venn
#' allGenes <- paste("Row", 1:nrow(simData),sep="")
#' if(require(limma)){
#'  vennC <- vennCounts(cbind(PerContrast= allGenes %in% pairsPerC$Feature,
#'  AllJoint=allGenes %in% pairsAll$Feature, FHier=allGenes %in%
#'  pairsAfterF$Feature))
#'	vennDiagram(vennC, main="FDR Overlap")
#' }
#'
#' #Do one cluster against all others
#' oneAll <- getBestFeatures(cl, contrastType="OneAgainstAll", contrastAdj="All",
#' p.value=0.05)
#'
#' #Do dendrogram testing
#' hcl <- makeDendrogram(cl)
#' allDendro <- getBestFeatures(hcl, contrastType="Dendro", contrastAdj=c("All"),
#' number=ncol(simData), p.value=0.05)
#'
#' # do DE on counts using voom
#' # compare results to if used simData instead (not on count scale).
#' # Again, not relevant for this silly example, but basic principle useful
#' testFVoom <- getBestFeatures(simCount, primaryCluster(cl), contrastType="F",
#' number=nrow(simData), isCount=TRUE)
#' plot(testF$P.Value[order(testF$Index)],
#' testFVoom$P.Value[order(testFVoom$Index)],log="xy")
#'
#' @export
#' @import limma
#' @importFrom stringr str_pad
#' @rdname getBestFeatures
setMethod(f = "getBestFeatures",
          signature = signature(x = "matrix"),
          definition = function(x, cluster,
                                contrastType=c("F", "Dendro", "Pairs", "OneAgainstAll"),
                                dendro=NULL, pairMat=NULL,
                                contrastAdj=c("All", "PerContrast", "AfterF"),
                                isCount=FALSE, normalize.method="none",...) {

            #... is always sent to topTable, and nothing else
            cl<-cluster
            if(is.factor(cl)) {
              warning("cluster is a factor. Converting to numeric, which may not result in valid conversion")
              cl <- .convertToNum(cl)
            }
            
            dat <- data.matrix(x)
            contrastType <- match.arg(contrastType)
            contrastAdj <- match.arg(contrastAdj)
            returnType <- "Table" 
            
            if(is.null(rownames(dat))) {
              rownames(dat) <- paste("Row", as.character(1:nrow(dat)), sep="")
            }
            
            tmp <- dat
            
            if(any(cl<0)){ #only use those assigned to a cluster to get good genes.
              whNA <- which(cl<0)
              tmp <- tmp[, -whNA]
              cl <- cl[-whNA]
            }
            ###--------
            ### Fix up the names
            ###--------
            pad<-if(length(unique(cl))<100) 2 else 3
            clPretty<-paste("Cl",stringr::str_pad(cl,width=pad,pad="0"),sep="")
            clLevels<-unique(cl[order(clPretty)])
            clPrettyLevels<-unique(clPretty[order(clPretty)])
            #get them ordered nicely.
            clNumFac<-factor(cl,levels=clLevels)
            
            if(contrastType=="Dendro") {
              clPrettyFac<-factor(cl,levels=clLevels,labels=clLevels)
              if(is.null(dendro)) {
                stop("must provide dendro")
              }
              if(!inherits(dendro,"dendrogram") && !inherits(dendro,"phylo4")){
				stop("dendro must be of class 'dendrogram' or 'phylo4'")
              }
            }
            else{
              clPrettyFac<-factor(cl,levels=clLevels,labels=clPrettyLevels)
            }
            
            if(contrastType %in% c("Pairs", "Dendro", "OneAgainstAll")) {
              ###Create fit for running contrasts
              designContr <- model.matrix(~ 0 + clPrettyFac)
              colnames(designContr) <- make.names(levels(clPrettyFac))
              
              if(isCount) {
                v <- voom(tmp, design=designContr, plot=FALSE,
                          normalize.method = normalize.method)
                fitContr <- lmFit(v, designContr)
              } else {
                fitContr <- lmFit(tmp, designContr)
              }
            }
            
            if(contrastType=="F" || contrastAdj=="AfterF") {
              xdat<-data.frame("Cluster"=clPrettyFac)
              designF<-model.matrix(~Cluster,data=xdat)
              #designF <- model.matrix(~clPrettyFac)
              
              if(isCount) {
                v <- voom(tmp, design=designF, plot=FALSE,
                          normalize.method = normalize.method)
                fitF <- lmFit(v, designF)
              } else {
                fitF <- lmFit(tmp, designF)
              }
            } else {
              fitF <- NULL
            }
            
            if(contrastType!="F") contr.result<-clusterContrasts(clNumFac,contrastType=contrastType,dendro=dendro,pairMat=pairMat,outputType = "limma", removeNegative = TRUE)
            tops <- if(contrastType=="F") .getBestFGenes(fitF,...) else .testContrasts(contr.result$contrastMatrix,contrastNames=contr.result$contrastNames,fit=fitContr,fitF=fitF,contrastAdj=contrastAdj,...)
            tops <- data.frame(IndexInOriginal=match(tops$Feature, rownames(tmp)),tops)
            if(returnType=="Index") {
              whGenes <- tops$IndexInOriginal
              names(whGenes) <- tops$Feature
              return(whGenes)
            }
            if(returnType=="Table") {
              return(tops)
            }
          }
)

#' @rdname getBestFeatures
#' @export
setMethod(f = "getBestFeatures",
          signature = signature(x = "ClusterExperiment"),
          definition = function(x,contrastType=c("F", "Dendro", "Pairs", "OneAgainstAll"),
                                isCount=FALSE, ...){
            contrastType <- match.arg(contrastType)
			cl<-primaryCluster(x)
			if(length(unique(cl[cl>0]))==1) stop("only single cluster in clustering -- cannot run getBestFeatures")
            if(contrastType=="Dendro") {
              if(is.null(x@dendro_clusters)) {
                stop("If `contrastType='Dendro'`, `makeDendrogram` must be run before `getBestFeatures`")
              } else {
                if(primaryClusterIndex(x)!= dendroClusterIndex(x)){
					#check if merge from cluster that made dendro
					if(primaryClusterIndex(x)==mergeClusterIndex(x) && x@merge_dendrocluster_index == dendroClusterIndex(x)){
						dendro<-.makeMergeDendrogram(x)
						if(is.null(dendro)) stop("Could not make merge dendrogram")
					}
					else stop("Primary cluster does not match the cluster on which the dendrogram was made. Either replace existing dendrogram with on using the primary cluster (via 'makeDendrogram'), or reset primaryCluster with 'primaryClusterIndex' to be equal to index of 'dendo_index' slot")
				}
                else dendro <- x@dendro_clusters
              }
            }
            
            if(isCount) {
              note(
                "If `isCount=TRUE` the data will be transformed with voom() rather than
with the transformation function in the slot `transformation`.
This makes sense only for counts.")
              dat <- assay(x)
            } else {
              dat <- transformData(x)
            }
            
            getBestFeatures(dat, primaryCluster(x), contrastType=contrastType, dendro=dendro,
                            isCount=isCount, ...)
            
}
)


.getBestFGenes<-function(fit,...){
	## basic limma design
	fit2 <- eBayes(fit)
	#tops<-topTableF(fit2,number=nGenes,genelist=rownames(fit$coef),...)
	tops <- topTableF(fit2,genelist=rownames(fit$coef),...)
	colnames(tops)[colnames(tops)=="ProbeID"]<-"Feature"

	return(tops)
}
.testContrasts<-function(contr.matrix, contrastNames=NULL, fit,fitF,contrastAdj,...){
  ncontr<-ncol(contr.matrix)
  fit2<-contrasts.fit(fit,contr.matrix)
  fit2<-eBayes(fit2)
  args<-list(...)
  if("p.value" %in% names(args)){
    p.value<-args$p.value
  }
  else p.value<-1 #default of topTable
  if("number" %in% names(args)){
    nGenes<-args$number
  }
  else nGenes<-10 #default of topTable
  
  
  #get raw p-values for all(!)
  getRaw<-function(ii){
    if(contrastAdj%in%c("AfterF","All")) {
      tt<-topTable(fit2,coef=ii, number=length(rownames(fit2$coef)),p.value=1,adjust.method="none",genelist=rownames(fit2$coef))
    }
    else{
      tt<-topTable(fit2,coef=ii, genelist=rownames(fit2$coef),...)
    }
    colnames(tt)[colnames(tt)=="ID"]<-"Feature"
    if(nrow(tt)>0){
      tt<-data.frame("Contrast"=unname(colnames(contr.matrix)[ii]),tt,row.names=NULL)
      if(!is.null(contrastNames)){
        tt<-data.frame("ContrastName"=contrastNames[ii],tt,row.names=NULL)
      }
    }
    return(tt)
  }
  tops<-do.call("rbind",lapply(1:ncontr,getRaw))
  if(contrastAdj=="AfterF" & p.value<1){
    #get p-value for F test for all genes, and only consider genes with significant F.
    fitF2<-eBayes(fitF)
    topsF<-topTable(fitF2,genelist=rownames(fit$coef),number=length(rownames(fit$coef)),adjust.method="BH")
    whGenesSigF<-topsF$ProbeID[which(topsF$adj.P.Val < p.value)]
    tops<-tops[tops$Feature %in% whGenesSigF,]
  }
  #do FDR correction on all raw p-values (that remain)
  if(contrastAdj%in%c("AfterF","All")) {
    tops$adj.P.Val<-p.adjust(tops$P.Value,method="BH")
    if(p.value<1) tops<-tops[tops$adj.P.Val<p.value,,drop=FALSE]
    if(nGenes<length(rownames(fit$coef))){
      #just return relevant number per contrast
      tops<-do.call("rbind",by(tops,tops$Contrast,function(x){x[1:min(c(nGenes,nrow(x))),,drop=FALSE]}))
    }
  }
  row.names(tops)<-NULL
  
  return(tops)
  
}

#' @importFrom phylobase descendants nodeLabels subset
.makeMergeDendrogram<-function(object){
	if(is.na(object@dendro_index)) stop("no dendrogram for this ClusterExperiment Object")
  #should this instead just silently return the existing?
	if(is.na(object@merge_index)) stop("no merging was done for this ClusterExperiment Object")
  if(object@merge_dendrocluster_index != object@dendro_index) stop("dendrogram of this object was made from different cluster than that of merge cluster.")
		#test mergeClusters actually subset of the cluster says merged
  whClusterNode<-which(!is.na(object@merge_nodeMerge[,"mergeClusterId"]))
  clusterNode<-object@merge_nodeMerge[whClusterNode,"Node"]
  clusterId<-object@merge_nodeMerge[whClusterNode,"mergeClusterId"]
  phylo4Obj <- .makePhylobaseTree(object@dendro_clusters, "dendro")
  newPhylo4<-phylo4Obj
  if(names(rootNode(phylo4Obj)) %in% clusterNode){
	  stop("coding error -- trying to make dendrogram from merge cluster when only 1 cluster in the clustering.")
  }
  for(node in clusterNode){
    #first remove tips of children nodes so all children of node in question are tips
    
	desc<-phylobase::descendants(newPhylo4, node, type = c("all")) #names are names
    whDescNodes<-which(names(desc) %in% phylobase::nodeLabels(newPhylo4))
    while(length(whDescNodes)>0){
      tipNodeDesc<-unique(unlist(phylobase::descendants(newPhylo4, desc[whDescNodes], type = c("tips")))) #internal ids, not names, and no names to it
      newPhylo4<-phylobase::subset(newPhylo4,tips.exclude=tipNodeDesc,trim.internal =FALSE)
      #redo to check fixed problem
      desc<-phylobase::descendants(newPhylo4, node, type = c("all"))
      whDescNodes<-which(names(desc) %in% phylobase::nodeLabels(newPhylo4))
    }
		
    #should only have tips now
    tipsRemove<-phylobase::descendants(newPhylo4, node, type = c("tips"))
    newPhylo4<-.safePhyloSubset(newPhylo4,tipsRemove=tipsRemove,nodeName=node) #use instead of subset, because run into problems in phylobase in subset when small tree.
  }

  ##################
  #Now need to change tip name to be that of the merge cluster
  #Currently tips should be either 
  #1) Name of combineMany cluster (i.e. integer) which needs to translate to a merge cluster
  #2) Node name which now should be a merge cluster id
  ##################
  newTips<-currTips<-phylobase::tipLabels(newPhylo4) #has *names* as entries
  #browser()

  #Solve 1) First:
  #Find the correspondence between old and new
  #replace the old (i.e. non-nodes) with the new
  corrsp<-getMergeCorrespond(object,by="original") #should be vector with names corresponding to original clusters, entries to merge clusters
  whOldCl<-which(currTips %in% names(corrsp)) #which are cluster names of the original; these are ones that should be 
  if(length(whOldCl)>0){
    newTips[whOldCl]<-corrsp[currTips[whOldCl]]
  }
  ## Solve 2) Now:
  ## should all be clusterNode
  if(!all(clusterNode %in% currTips)) stop("coding error -- some cluster nodes didn't wind up as tips of new tree")
  mClusterNode<-match(clusterNode, currTips)
  newTips[mClusterNode]<-as.character(clusterId)
  names(newTips)<-names(currTips) #this is the internal numbering of the nodes
  phylobase::tipLabels(newPhylo4)<-newTips
  
  #--------
  #just some checks didn't screw up
  #--------
  whMergeNode<-which(currTips %in% clusterNode)
  if(length(intersect(whMergeNode,whOldCl))) stop("coding error -- should be no overlap bwettween merged node in tree tips and old clusters")
  if(length(union(whMergeNode,whOldCl))!= length(currTips)) stop("coding error -- all tips should be either old clusters of merged nodes")
  mCl<-clusterMatrix(object)[,object@merge_index]
  mCl<-unique(mCl[mCl>0])
  if(length(currTips)!= length(mCl)) stop("coding error -- number of tips of new tree not equal to the number of clusters in merged cluster")
  if(length(currTips)!= length(mCl)) stop("coding error -- number of tips of new tree not equal to the number of clusters in merged cluster")
  if(!identical(sort(unname(as.character(mCl))),sort(unname(tipLabels(newPhylo4))))){
	  stop("coding error -- names of new tips of tree do not match cluster ids")
  }
  return(newPhylo4)
  #newPhylo4<-.force.ultrametric(newPhylo4) 
  #convert back to dendrogram class and return
  ##as.dendrogram.phylo in dendextend first converts to hclust with ape function, then dendrogram with their (non-exported) function as.dendrogram.hclust:
  ## as.dendrogram(ape::as.hclust.phylo(object))
  ##Hit a problem from ape that doesn't return matrix in merge entry if only 2 tips, single node. Reported to ape.
  ##Have to restep through and manually fix it
  #previously had to do, but no longer using:
  #@importFrom ape as.hclust.phylo
  #@import dendextend (Couldn't import from because the as.dendrogram.hclust was not exported)
  # xxhclust<-ape::as.hclust.phylo(as(newPhylo4,"phylo"))
  # if(is.null(dim(xxhclust$merge))) xxhclust$merge<-matrix(xxhclust$merge,ncol=2)
  # return(as.dendrogram(xxhclust))
  #return(as.dendrogram(as(newPhylo4,"phylo"))) #as.dendrogram.phylo from dendextend, not exported...

}
 
#' @importFrom phylobase nodeHeight tipLabels edgeLength edges edgeId 
##From http://blog.phytools.org/2017/03/forceultrametric-method-for-ultrametric.html
.force.ultrametric<-function(tree){
	if(!inherits(tree,"phylo4")) stop("tree must be of class phylo4")
	
	allTips<-phylobase::tipLabels(tree)
	depthToTips<-phylobase::nodeHeight(tree,allTips,from="root")
	maxD<-max(depthToTips)
	addValue<-maxD-depthToTips
	allLen<-phylobase::edgeLength(tree)
  edgeMat<-phylobase::edges(tree)
  tipIds<-as.numeric(names(allTips))
  m<-match(tipIds,edgeMat[,2])
  edgeIds<-paste(edgeMat[m,1],edgeMat[m,2],sep="-")

  #check didn't do something stupid:
  checkTipEdges<-phylobase::edgeId(tree,type="tip")
  if(!identical(sort(unname(checkTipEdges)),sort(unname(edgeIds)))) stop("coding error -- didn't correctly get edge ids for tips")

  #replace with new edges:
	allLen[edgeIds]<-allLen[edgeIds]+addValue
	phylobase::edgeLength(tree)<-allLen
	tree
}


