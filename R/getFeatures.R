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
#'   \code{\link{topTableF}} (see \code{\link{limma}} package)
#' @param normalize.method character value, passed to \code{\link{voom}} in
#'   \code{\link{limma}} package. Only used if \code{countData=TRUE}.
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
#'   voom with log2(x+0.5). Alternatively, `isCount=FALSE` for a
#'   `ClusterExperiment` object will cause the DE to be performed with `limma`
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
            cl <- factor(cl)

            if(contrastType=="Dendro") {
              if(is.null(dendro)) {
                stop("must provide dendro")
              }
              if(!inherits(dendro,"dendrogram")) {
                stop("dendro must be of class 'dendrogram'")
              }
            }

            if(contrastType %in% c("Pairs", "Dendro", "OneAgainstAll")) {
              designContr <- model.matrix(~ 0 + cl)
              colnames(designContr) <- make.names(levels(cl))

              if(isCount) {
                v <- voom(tmp, design=designContr, plot=FALSE,
                                 normalize.method = normalize.method)
                fitContr <- lmFit(v, designContr)
              } else {
                fitContr <- lmFit(tmp, designContr)
              }
            }

            if(contrastType=="F" || contrastAdj=="AfterF") {
              designF <- model.matrix(~cl)

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
            #browser()
            if(contrastType!="F") contr.result<-clusterContrasts(cl,contrastType=contrastType,dendro=dendro,pairMat=pairMat,outputType = "limma", removeNegative = TRUE)
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
            if(contrastType=="Dendro") {
              if(is.null(x@dendro_clusters)) {
                stop("If `contrastType='Dendro'`, `makeDendrogram` must be run before `getBestFeatures`")
              } else {
                if(primaryClusterIndex(x)!= dendroClusterIndex(x)) stop("Primary cluster does not match the cluster on which the dendrogram was made. Either replace existing dendrogram with on using the primary cluster (via 'makeDendrogram'), or reset primaryCluster with 'primaryClusterIndex' to be equal to index of 'dendo_index' slot")
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
              dat <- transform(x)
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




