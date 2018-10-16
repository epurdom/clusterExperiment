.demethods<-c("edgeR","limma","limma-voom")
#' @name getBestFeatures
#' @aliases getBestFeatures getBestFeatures,matrixOrHDF5-method
#' @title Function for finding best features associated with clusters
#' @description Calls limma on input data to determine features most associated 
#'   with found clusters (based on an F-statistic, pairwise comparisons, or 
#'   following a tree that clusters the clusters).
#' @aliases getBestFeatures
#' @param x data for the test. Can be a numeric matrix or a 
#'   \code{\link{ClusterExperiment}}.
#' @param cluster A numeric vector with cluster assignments. ``-1'' indicates 
#'   the sample was not assigned to a cluster.
#' @param whichCluster which clustering to use in performing the tests.
#' @param contrastType What type of test to do. `F' gives the omnibus 
#'   F-statistic, `Dendro' traverses the given dendrogram and does contrasts of 
#'   the samples in each side,  `Pairs' does pair-wise contrasts based on the 
#'   pairs given in pairMat (if pairMat=NULL, does all pairwise), and 
#'   `OneAgainstAll' compares each cluster to the average of all others. Passed 
#'   to \code{\link{clusterContrasts}}
#' @param contrastAdj What type of FDR correction to do for contrasts tests 
#'   (i.e. if contrastType='Dendro' or 'Pairs').
#' @param DEMethod character vector describing how the differential expression 
#'   analysis should be performed (replaces previous argument \code{isCount}. 
#'   See details.
#' @param weights weights to use in by edgeR. If \code{x} is a matrix, then
#'   weights should be a matrix of weights, of the same dimensions as \code{x}.
#'   If \code{x} is a \code{ClusterExperiment} object \code{weights} can be a
#'   either a matrix, as previously described, or a character or numeric index
#'   to an assay in \code{x} that contains the weights. We recommend that 
#'   weights be stored as an assay with name \code{"weights"} so that the
#'   weights will also be used with \code{\link{mergeClusters}}, and this is the
#'   default. Setting \code{weights=NULL} ensures that weights will NOT be used,
#'   and only the standard edgeR.
#' @param dgeArgs a list of arguments to pass to \code{\link[edgeR]{DGEList}}
#'   which is the starting point for both \code{edgeR} and \code{limma-voom}
#'   methods of DE. This includes normalization factors/total count values etc.
#' @param ... If \code{x} is a matrix, these are options to pass to 
#'   \code{\link{topTable}} or \code{\link[limma]{topTableF}} (see 
#'   \code{\link[limma]{limma}} package). If \code{x} is a 
#'   \code{ClusterExperiment} object, these arguments can also be those to pass 
#'   to the matrix version.
#' @inheritParams clusterContrasts
#' @details getBestFeatures returns the top ranked features corresponding to a 
#'   cluster assignment. It uses either limma or edgeR to fit the models, and 
#'   limma/edgeR functions \code{\link[limma]{topTable}} or 
#'   \code{\link[limma]{topTableF}} to find the best features. See the options 
#'   of these functions to put better control on what gets returned (e.g. only 
#'   if significant, only if log-fc is above a certain amount, etc.). In 
#'   particular, set `number=` to define how many significant features to return
#'   (where number is per contrast for the `Pairs` or `Dendro` option)
#' @details  \code{DEMethod} triggers what type of differential expression 
#'   analysis will be performed. Three options are available: limma, edgeR, and 
#'   limma with a voom corrections. The last two options are only appropriate 
#'   for count data. If the input \code{x} is a \code{ClusterExperiment} object,
#'   and \code{DEMethod="limma"}, then the data analyzed for DE will be
#'   \emph{after} taking the transformation of the data (as given in the
#'   transformation slot of the object). For the options "limma-voom" and
#'   "edgeR", the transformation slot will be ignored and only the counts data
#'   (as specified by the \code{whichAssay} slot) will be passed to the
#'   programs. Note that for "limma-voom" this implies that the data will be
#'   transformed by voom with the function log2(x+0.5). If \code{weights} is not
#'   \code{NULL}, and \code{DEMethod="edgeR"}, then the function
#'   \code{glmWeightedF} from the \code{zinbwave} package is run; otherwise
#'   \code{glmLRT} from \code{edgeR}.
#' @details Note that the argument \code{DEMethod} replaces the previous option 
#'   \code{isCount}, to decide on the method of DE.
#' @details When `contrastType` argument implies that the best features should 
#'   be found via contrasts (i.e. 'contrastType' is `Pairs` or `Dendro`), then 
#'   then `contrastAdj` determines the type of multiple testing correction to 
#'   perform. `PerContrast` does FDR correction for each set of contrasts, and 
#'   does not guarantee control across all the different contrasts (so probably 
#'   not the preferred method). `All` calculates the corrected p-values based on
#'   FDR correction of all of the contrasts tested. `AfterF` controls the FDR 
#'   based on a hierarchical scheme that only tests the contrasts in those genes
#'   where the omnibus F statistic is significant. If the user selects `AfterF`,
#'   the user must also supply an option `p.value` to have any effect, and then 
#'   only those significant at that p.value level will be returned. Note that 
#'   currently the correction for `AfterF` is not guaranteed to control the FDR;
#'   improvements will be added in the future.
#' @details  Note that the default option for \code{\link[limma]{topTable}} is 
#'   to not filter based on adjusted p-values (\code{p.value = 1}) and return 
#'   only the top 10 most significant (\code{number = 10}) -- these are options 
#'   the user can change (these arguments are passed via the \code{...} in 
#'   \code{getBestFeatures}). In particular, it only makes sense to set 
#'   \code{requireF = TRUE} if \code{p.value} is meaningful (e.g. 0.1 or 0.05); 
#'   the default value of \code{p.value = 1} will not result in any effect on 
#'   the adjusted p-value otherwise.
#' @return A \code{data.frame} in the same format as 
#'   \code{\link[limma]{topTable}}, or \code{\link[edgeR]{topTags}}. The output
#'   differs between these two programs, mainly in the naming of columns.
#'   Furthermore, if weights are used, an additional column \code{padjFilter} is
#'   included as the result of running \code{\link[zinbwave]{glmWeightedF}} with
#'   default option \code{independentFiltering = TRUE}. The following column
#'   names are the same between all of the DE methods.
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
#' the dendrogram. If \code{x} is a \code{ClusterExperiment} object, this name 
#' will make use of the user defined names of the cluster or node in \code{x}.}
#' 
#' \item{\code{InternalName}}{ Only present if \code{x} is a
#' \code{ClusterExperiment} object. In this case this column will give the name
#' of the contrast using the internal ids of the clusters and nodes, not the
#' user-defined names. This provides stability in matching the contrast if the
#' user has changed the names since running \code{getBestFeatures}}
#' 
#' \item{\code{P.Value}}{ The unadjusted p-value (changed from \code{PValue} in
#' \code{topTags})}
#' 
#' \item{\code{adj.P.Val}}{ The adjusted p-value (changed from \code{FDR} or
#' \code{FWER} in \code{topTags})}
#' }
#'
#' @references Ritchie, ME, Phipson, B, Wu, D, Hu, Y, Law, CW, Shi, W, and
#'   Smyth, GK (2015). limma powers differential expression analyses for
#'   RNA-sequencing and microarray studies. Nucleic Acids Research 43, e47.
#'   http://nar.oxfordjournals.org/content/43/7/e47
#' @references Law, CW, Chen, Y, Shi, W, and Smyth, GK (2014). Voom: precision
#'   weights unlock linear model analysis tools for RNA-seq read counts. Genome
#'   Biology 15, R29. http://genomebiology.com/2014/15/2/R29
#' @references Smyth, G. K. (2004). Linear models and empirical Bayes methods
#'   for assessing differential expression in microarray experiments.
#'   Statistical Applications in Genetics and Molecular Biology, Volume 3,
#'   Article 3. http://www.statsci.org/smyth/pubs/ebayes.pdf
#' @seealso \code{\link[edgeR]{glmLRT}} \code{\link[zinbwave]{glmWeightedF}}
#'   \code{\link[limma]{topTable}} \code{\link[edgeR]{topTags}}
#' @examples
#' data(simData)
#'
#' #create a clustering, for 8 clusters (truth was 4)
#' cl <- clusterSingle(simData, subsample=FALSE,
#' sequential=FALSE, mainClusterArgs=list(clusterFunction="pam", clusterArgs=list(k=8)))
#'
#' #basic F test, return all, even if not significant:
#' testF <- getBestFeatures(cl, contrastType="F", number=nrow(simData),
#' DEMethod="limma")
#'
#' #Do all pairwise, only return significant, try different adjustments:
#' pairsPerC <- getBestFeatures(cl, contrastType="Pairs", contrastAdj="PerContrast",
#' p.value=0.05, DEMethod="limma")
#' pairsAfterF <- getBestFeatures(cl, contrastType="Pairs", contrastAdj="AfterF",
#' p.value=0.05, DEMethod="limma")
#' pairsAll <- getBestFeatures(cl, contrastType="Pairs", contrastAdj="All",
#' p.value=0.05, DEMethod="limma")
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
#' p.value=0.05,DEMethod="limma")
#'
#' #Do dendrogram testing
#' hcl <- makeDendrogram(cl)
#' allDendro <- getBestFeatures(hcl, contrastType="Dendro", contrastAdj=c("All"),
#' number=ncol(simData), p.value=0.05,DEMethod="limma")
#'
#' # do DE on counts using voom
#' # compare results to if used simData instead (not on count scale).
#' testFVoom <- getBestFeatures(simCount, primaryCluster(cl), contrastType="F",
#' number=nrow(simData), DEMethod="limma-voom")
#' plot(testF$P.Value[order(testF$Index)],
#' testFVoom$P.Value[order(testFVoom$Index)],log="xy")
#'
#' # do DE on counts using edgeR, compare voom
#' testFEdge <- getBestFeatures(simCount, primaryCluster(cl), contrastType="F",
#' n=nrow(simData), DEMethod="edgeR")
#' plot(testFVoom$P.Value[order(testFVoom$Index)],
#' testFEdge$P.Value[order(testFEdge$Index)],log="xy")
#' @export
#' @import limma
#' @importFrom stringr str_pad
#' @rdname getBestFeatures
setMethod(f = "getBestFeatures",
signature = signature(x = "matrixOrHDF5"),
definition = function(x, cluster,
                      contrastType=c("F", "Dendro", "Pairs", "OneAgainstAll"),
                      dendro=NULL, pairMat=NULL, weights = NULL,
                      contrastAdj=c("All", "PerContrast", "AfterF"),
                      DEMethod=c("edgeR","limma","limma-voom"),
											dgeArgs=NULL,...) {
  
  
  #... is always sent to topTable/topTags, and nothing else
	DEMethod<-match.arg(DEMethod)
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
    rownames(dat) <- paste("Row", as.character(seq_len(nrow(dat))), sep="")
  }
	
  #only use those assigned to a cluster to get good genes.
  tmp <- dat  
  if(any(cl<0)){ 
    whNA <- which(cl<0)
    tmp <- tmp[, -whNA]
    cl <- cl[-whNA]
    if(!is.null(weights)){
      weights<- weights[,-whNA]
    }
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
	
  if(DEMethod %in% c("limma-voom","edgeR")){
  	dge<-do.call(edgeR::DGEList,c(list(counts=tmp),dgeArgs))
    if(!is.null(weights) & DEMethod=="edgeR") dge$weights <- weights
  }
  if(contrastType %in% c("Pairs", "Dendro", "OneAgainstAll")) {
    ###Create fit for running contrasts
    designContr <- model.matrix(~ 0 + clPrettyFac)
    colnames(designContr) <- make.names(levels(clPrettyFac))
    if(DEMethod == "limma-voom"){
      v <- voom(dge, design=designContr, plot=FALSE,
                normalize.method = "none")
      fitContr <- lmFit(v, designContr)
		}
		else if(DEMethod=="edgeR"){
      fitContr <- edgeR::estimateDisp(dge, designContr)
      fitContr <- edgeR::glmFit(fitContr,design = designContr)
		}
		else if(DEMethod=="limma"){
      fitContr <- lmFit(tmp, designContr)
    }
  }  
  if(contrastType=="F" || contrastAdj=="AfterF") {
    xdat<-data.frame("Cluster"=clPrettyFac)
    designF<-model.matrix(~Cluster,data=xdat)
    
    if(DEMethod == "limma-voom"){
      v <- voom(dge, design=designF, plot=FALSE,
                normalize.method = "none")
      fitF <- lmFit(v, designF)
    }
		else if(DEMethod=="edgeR"){
        fitF <- edgeR::estimateDisp(dge, design = designF)
        fitF <- edgeR::glmFit(fitF,design = designF)
    } else if(DEMethod=="limma") {
      fitF <- lmFit(tmp, designF)
    }
  } else {
    fitF <- NULL
  }
  	
  if(contrastType!="F"){
		contr.result<-clusterContrasts(clNumFac,contrastType=contrastType,dendro=dendro,pairMat=pairMat,outputType = "limma", removeUnassigned = TRUE)
  }
  if(contrastType=="F"){
    tops <- .getBestFGenes(fitF, DEMethod=DEMethod,...)
  }
  else{
    tops<-.testContrasts(contr.result$contrastMatrix,
			contrastNames=contr.result$contrastNames,
			fit=fitContr,fitF=fitF,DEMethod=DEMethod,
			contrastAdj=contrastAdj, ...)
  }
  if(contrastType=="Pairs"){
	  tops <- data.frame(ContrastName=tops$Contrast,tops)
  	  
  }
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
#' @param whichAssay numeric or character specifying which assay to use. See
#'   \code{\link[SummarizedExperiment]{assay}} for details.
#' @export
setMethod(
  f = "getBestFeatures",
  signature = signature(x = "ClusterExperiment"),
  definition = 
    function(x, contrastType=c("F", "Dendro", "Pairs", "OneAgainstAll"), whichCluster="primary",
             whichAssay=1,DEMethod, weights=if("weights" %in% assayNames(x)) "weights" else NULL,...)
      {
    whCl<-.convertSingleWhichCluster(x,whichCluster,list(...))
      contrastType <- match.arg(contrastType)
    cl<-clusterMatrix(x)[,whCl]
    if(length(unique(cl[cl>0]))==1) stop("only single cluster in clustering -- cannot run getBestFeatures")
    if(contrastType=="Dendro") {
      if(is.null(x@dendro_clusters)) {
        stop("If `contrastType='Dendro'`, `makeDendrogram` must be run before `getBestFeatures`")
      } else {
        if(whCl!= dendroClusterIndex(x)){
          #check if merge from cluster that made dendro
          if(whCl==mergeClusterIndex(x) && x@merge_dendrocluster_index == dendroClusterIndex(x)){
            dendro<-.makeMergeDendrogram(x)
            if(is.null(dendro)) stop("Could not make merge dendrogram")
          }
          else stop("Cluster given in 'whichCluster' does not match either the cluster on which the dendrogram was made or the merge cluster from this dendrogram. Either replace the existing dendrogram with one based on using your preferred clustering (via 'makeDendrogram'), or set 'whichCluster' to be equal to index of 'dendo_index' or 'merge_index' slot")
        }
        else dendro <- x@dendro_clusters
      }
    }
    passedArgs<-list(...)
    if(DEMethod=="limma") dat<-transformData(x,whichAssay=whichAssay)
    else dat<-assay(x,whichAssay)
    
    if(!is.null(weights) && (is.character(weights) || (is.vector(weights) && is.numeric(weights)))  && length(weights)==1){
    		tops<-getBestFeatures(dat, cl, contrastType=contrastType, dendro=dendro, weights=assay(x, weights),DEMethod=DEMethod,...) 
    }
    else tops<-getBestFeatures(dat, cl, contrastType=contrastType, dendro=dendro, weights=weights,DEMethod=DEMethod,...)


    #### Fix up names
	#### Add column $InternalName -- 
	if(contrastType!="F"){
		wh<-which(colnames(tops)%in%c("IndexInOriginal", "ContrastName"))
		tops<-data.frame(tops[,wh],InternalName=tops$ContrastName,tops[,-wh])
		legMat<-clusterLegend(x)[[whCl]]
		#### Each contrastType needs different parsing
    if(contrastType=="Pairs"){
  	  #ContrastName in form of 'Cl01-Cl02' ...
  	  parseName<-gsub("Cl","",tops$InternalName)
  	  parseName<-strsplit(parseName,"-")
	    #note, incase there is padding in name, but not in clusterIds
  	  firstCl<-as.character(as.numeric(sapply(parseName,.subset2,1)))
  	  secondCl<-as.character(as.numeric(sapply(parseName,.subset2,2)))
  	  m1<-match(firstCl,as.numeric(legMat[,"clusterIds"]))
  	  m2<-match(secondCl,as.numeric(legMat[,"clusterIds"]))
	    if(any(is.na(m1))||any(is.na(m2))) stop("coding error -- cannot match parse cluster id from contrastName (no match to clusterLegend)")
  	  tops$ContrastName<-paste(legMat[m1,"name"],legMat[m2,"name"],sep="-")
    }
    if(contrastType=="OneAgainstAll"){
	    #ContrastName in form of Cl01
	  	parseName<-as.character(as.numeric(gsub("Cl","",tops$InternalName)))
  	  m1<-match(parseName,as.numeric(legMat[,"clusterIds"]))
  	  tops$ContrastName<-legMat[m1,"name"]
    }
    if(contrastType=="Dendro"){
	    #ContrastName in form of InternalNodeId5
	    m <- .matchToDendroData(inputValue=tops$InternalName, dendro, matchColumn="NodeId", returnColumn="NodeIndex")
	    tops$ContrastName<-phylobase::labels(dendro)[m]
    }
		tops$ContrastName<-factor(tops$ContrastName) #make it consistent with the others and the results of matrix version
	}
	return(tops)
  }
)

#' @importFrom zinbwave glmWeightedF
.getBestFGenes<-function(fit, DEMethod, ...){
	if(DEMethod%in% c("limma","limma-voom")){
	  ## basic limma design
	  fit2 <- eBayes(fit)
	  tops <- topTableF(fit2,genelist=rownames(fit$coef),...)
	  colnames(tops)[colnames(tops)=="ProbeID"]<-"Feature"
	}
	else if(DEMethod=="edgeR"){
    if(!is.null(fit$weights)){
      lrt <- zinbwave::glmWeightedF(fit, coef = c(2:ncol(fit$design)) )#Test if all the betas are zero i.e if the groups have the same gene expression
    }
    else{
      lrt <- edgeR::glmLRT(fit, coef = c(2:ncol(fit$design)))
    }
    tops<-edgeR::topTags(lrt,...)
		#Note, fit$coefficients are natural log, while output of topTags is log2 (see topTags documentation). Hence to add intercept, need to multiply by log2(exp(1))
    tops<- data.frame(Feature = rownames(tops), logFC.Intercept=log2(exp(1))*fit$coefficients[rownames(tops),"(Intercept)"],tops) 
	  colnames(tops)[colnames(tops)=="PValue"]<-"P.Value"
	  if("FDR" %in% colnames(tops)) colnames(tops)[colnames(tops)=="FDR"]<-"adj.P.Val"
		if("FWER" %in% colnames(tops)) colnames(tops)[colnames(tops)=="FWER"]<-"adj.P.Val"

  }
  
  return(tops)
}


.testContrasts<-function(contr.matrix, contrastNames=NULL, fit,fitF,contrastAdj, DEMethod,...){
  ncontr<-ncol(contr.matrix)
  if(DEMethod%in%c("limma","limma-voom")){
    fit2<-contrasts.fit(fit,contr.matrix)
    fit2<-eBayes(fit2)
  }
  else if(DEMethod=="edgeR"){
    fit2<-fit
  }
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
    if(DEMethod%in%c("limma","limma-voom")){
      if(contrastAdj%in%c("AfterF","All")) {
        tt<-topTable(fit2,coef=ii, number=length(rownames(fit2$coef)),p.value=1,adjust.method="none",genelist=rownames(fit2$coef))
      }
      else{
        tt<-topTable(fit2,coef=ii, genelist=rownames(fit2$coef),...)
      }
	    colnames(tt)[colnames(tt)=="ID"]<-"Feature"
    }
    if(DEMethod=="edgeR") {
      if(!is.null(fit2$weights)) fit2 <- zinbwave::glmWeightedF(fit2, contrast = contr.matrix[,ii])
      else fit2 <- edgeR::glmLRT(fit2, contrast = contr.matrix[,ii])
      if(contrastAdj%in%c("AfterF","All")) {
        tt<-edgeR::topTags(fit2, n=length(rownames(fit2$coef)),p.value=1,adjust.method="none")
      }
      else{
        tt<-edgeR::topTags(fit2,...)
      }
			if(is.null(rownames(tt))){rownames(tt)<- (seq_len(nGenes))}
      tt<- data.frame("Feature"=rownames(tt), tt)
			colnames(tt)[colnames(tt)=="PValue"]<-"P.Value"
			if("FWER" %in% colnames(tt)) colnames(tt)[colnames(tt)=="FWER"]<-"adj.P.Val"
			if("FDR" %in% colnames(tt)) colnames(tt)[colnames(tt)=="FDR"]<-"adj.P.Val"
    }
    if(nrow(tt)>0){
      tt<-data.frame("Contrast"=unname(colnames(contr.matrix)[ii]),tt,row.names=NULL)
      if(!is.null(contrastNames)){
        tt<-data.frame("ContrastName"=contrastNames[ii],tt,row.names=NULL)
      }
    }
    return(tt)
  }
  tops<-do.call("rbind",lapply(seq_len(ncontr),getRaw))
  if(contrastAdj=="AfterF" & p.value<1){
    #get p-value for F test for all genes, and only consider genes with significant F/LRT.
    if(DEMethod%in%c("limma","limma-voom")){
      fitF2<-eBayes(fitF)
      topsF<-topTable(fitF2,genelist=rownames(fit$coef),number=length(rownames(fit$coef)),adjust.method="BH",p.value=1)
      whGenesSigF<-topsF$ProbeID[which(topsF$adj.P.Val < p.value)]
    }
    else if(DEMethod %in% c("edgeR")) {
	    if(!is.null(fitF$weights)){
	      lrt <- zinbwave::glmWeightedF(fitF, coef = c(2:ncol(fitF$design)) )#Test if all the betas are zero i.e if the groups have the same gene expression
	    }
	    else{
	      lrt <- edgeR::glmLRT(fitF, coef = c(2:ncol(fitF$design)))
	    }
			topsF<-edgeR::topTags(lrt,genelist=rownames(fitF$coef),n=length(rownames(fitF$coef)),adjust.method="BH",p.value=1)
      whGenesSigF<-rownames(topsF)[which(topsF$FDR < p.value)]
    }
    tops<-tops[tops$Feature %in% whGenesSigF,]
  }
  #do FDR correction on all raw p-values (that remain)
  if(contrastAdj%in%c("AfterF","All")) {
    tops$adj.P.Val<-p.adjust(tops$P.Value,method="BH")
    if(p.value<1) tops<-tops[tops$adj.P.Val<p.value,,drop=FALSE]
    if(nGenes<length(rownames(fit$coef))){
      #just return relevant number per contrast
      tops<-do.call("rbind",by(tops,tops$Contrast,function(x){x[seq_len(min(c(nGenes,nrow(x)))),,drop=FALSE]}))
    }
  }
  row.names(tops)<-NULL  
  return(tops)

}

.makeMergeDendrogram<-function(object){
	#this function returns a dendrogram where the tips are the merged clusters, rather than the clusters of the original dendrogram.
  if(is.na(object@dendro_index)) stop("no dendrogram for this ClusterExperiment Object")
  #should this instead just silently return the existing?
  if(is.na(object@merge_index)) stop("no merging was done for this ClusterExperiment Object")
  if(object@merge_dendrocluster_index != object@dendro_index) stop("dendrogram of this object was made from different cluster than that of merge cluster.")
  
  #-----
  #test mergeClusters actually subset of the cluster says merged -- 
  # -- should go back to this with new phylo4d object
  #-----
  whClusterNode<-which(!is.na(object@merge_nodeMerge[,"mergeClusterId"]))
  clusterNode<-object@merge_nodeMerge[whClusterNode,"NodeId"]
  nodePruneIndex <- .matchToDendroData(inputValue=clusterNode,dendro=object@dendro_clusters,matchColumn="NodeId",returnColumn="NodeIndex")
  mergePhylo<-.pruneToNodes(phylo4=object@dendro_clusters,nodesPruned=nodePruneIndex)
  ##################
  #Now need to change tip labels to be that of the merge cluster ids
  ##################  
  tipIndex<-phylobase::getNode(mergePhylo,type="tip")
  clusterIds <- .matchToDendroData(inputValue=tipIndex, dendro=mergePhylo, matchColumn="NodeIndex", returnColumn="ClusterIdMerge")
  phylobase::labels(mergePhylo,type="tip")<-gsub("ClusterId","",clusterIds)
  return(mergePhylo)

}
  # old .makeMergeDendrogram code:
  #
  # phylo4Obj <- .convertToPhyClasses(object@dendro_clusters,"phylo4",convertNodes=TRUE,convertTips=TRUE)
  # newPhylo4<-phylo4Obj
  # if(names(rootNode(phylo4Obj)) %in% clusterNode){
  #   stop("coding error -- trying to make dendrogram from merge cluster when only 1 cluster in the clustering.")
  # }
  #
  #
  # #-----
  # # remove tips
  # # Note that subset function only removes specified tips (i.e. give tips to exclude)
  # # Doesn't allow you to give nodes you want to keep...
  # #-----
  # #After have updated merge, should come back to this...
  # for(node in clusterNode){
  #   #first remove tips of children nodes so all children of node in question are tips
  #
  #   desc<-phylobase::descendants(newPhylo4, node, type = c("all")) #names are names
  #   whDescNodes<-which(names(desc) %in% phylobase::nodeLabels(newPhylo4))
  #   while(length(whDescNodes)>0){
  #     tipNodeDesc<-unique(unlist(phylobase::descendants(newPhylo4, desc[whDescNodes], type = c("tips")))) #internal ids, not names, and no names to it
  #     newPhylo4<-phylobase::subset(newPhylo4,tips.exclude=tipNodeDesc,trim.internal =FALSE)
  #     #redo to check fixed problem
  #     desc<-phylobase::descendants(newPhylo4, node, type = c("all"))
  #     whDescNodes<-which(names(desc) %in% phylobase::nodeLabels(newPhylo4))
  #   }
  #   #should only have tips now
  #   tipsRemove<-phylobase::descendants(newPhylo4, node, type = c("tips"))
  #   newPhylo4<-.safePhyloSubset(newPhylo4,tipsRemove=tipsRemove,nodeName=node) #use instead of subset, because run into problems in phylobase in subset when small tree.
  # }
  #
  # ##################
  # #Now need to change tip name to be that of the merge cluster
  # #Currently tips should be either
  # #1) Name of makeConsensus cluster (i.e. integer) which needs to translate to a merge cluster
  # #2) Node name which now should be a merge cluster id
  # ##################
  # newTips<-currTips<-phylobase::tipLabels(newPhylo4) #has *names* as entries
  # #
  #
  # #Solve 1) First:
  # #Find the correspondence between old and new
  # #replace the old (i.e. non-nodes) with the new
  # corrsp<-getMergeCorrespond(object,by="original") #should be vector with names corresponding to original clusters, entries to merge clusters
  # whOldCl<-which(currTips %in% names(corrsp)) #which are cluster names of the original; these are ones that should be
  # if(length(whOldCl)>0){
  #   newTips[whOldCl]<-corrsp[currTips[whOldCl]]
  # }
  # ## Solve 2) Now:
  # ## should all be clusterNode
  # if(!all(clusterNode %in% currTips)) stop("coding error -- some cluster nodes didn't wind up as tips of new tree")
  # mClusterNode<-match(clusterNode, currTips)
  # newTips[mClusterNode]<-as.character(clusterId)
  # names(newTips)<-names(currTips) #this is the internal numbering of the nodes
  # phylobase::tipLabels(newPhylo4)<-newTips
  #
  # #--------
  # #just some checks didn't screw up
  # #--------
  # whMergeNode<-which(currTips %in% clusterNode)
  # if(length(intersect(whMergeNode,whOldCl))) stop("coding error -- should be no overlap bwettween merged node in tree tips and old clusters")
  # if(length(union(whMergeNode,whOldCl))!= length(currTips)) stop("coding error -- all tips should be either old clusters of merged nodes")
  # mCl<-clusterMatrix(object)[,object@merge_index]
  # mCl<-unique(mCl[mCl>0])
  # if(length(currTips)!= length(mCl)) stop("coding error -- number of tips of new tree not equal to the number of clusters in merged cluster")
  # if(!identical(sort(unname(as.character(mCl))),sort(unname(phylobase::tipLabels(newPhylo4))))){
  #   stop("coding error -- names of new tips of tree do not match cluster ids")
  # }
  # return(newPhylo4)


