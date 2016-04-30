#' @title Function for finding best features associated with clusters
#' @description Calls limma on input data to determine features most associated
#'   with found clusters (based on an F-statistic, pairwise comparisons, or
#'   following a tree that clusters the clusters).
#' @aliases getBestFeatures
#' @param x data for the test. Can be a numeric matrix or a
#'   \code{\link{ClusterExperiment}}.
#' @param cl A numeric vector with cluster assignments.
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
#' @inheritParams clusterContrasts
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
#' cl <- clusterSingle(simData, clusterFunction="pam", subsample=FALSE,
#' sequential=FALSE, clusterDArgs=list(k=8))
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
                                returnType=c("Table", "Index"),
                                contrastAdj=c("All", "PerContrast", "AfterF"),
                                isCount=FALSE, normalize.method="none",...) {

            #... is always sent to topTable, and nothing else
            cl<-cluster
             if(is.factor(cl)) {
              warning("cluster is a factor. Converting to numeric, which may not result in valid conversion")
              cl <- as.numeric(as.character(cl))
            }

            dat <- data.matrix(x)
            contrastType <- match.arg(contrastType)
            contrastAdj <- match.arg(contrastAdj)
            returnType <- match.arg(returnType)

            if(is.null(rownames(dat))) {
              rownames(dat) <- paste("Row", as.character(1:nrow(dat)), sep="")
            }

            tmp <- dat

            if(any(cl== -1)){ #only use those assigned to a cluster to get good genes.
              whNA <- which(cl== -1)
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


            tops <- switch(contrastType,
                           "F"=.getBestFGenes(fitF,...),
                           "Dendro"=.getBestDendroGenes(cl=cl, dendro=dendro,
                                                        contrastAdj=contrastAdj,
                                                        fit=fitContr, fitF=fitF,
                                                        ...),
                           "Pairs"=.getBestPairsGenes(cl=cl, pairMat=pairMat,
                                                      contrastAdj=contrastAdj,
                                                      fit=fitContr, fitF=fitF,
                                                      ...),
                           "OneAgainstAll"=.getBestOne(cl,
                                                       contrastAdj=contrastAdj,
                                                       fit=fitContr, fitF=fitF,
                                                       ...)
            )

            tops <- data.frame(IndexInOriginal=match(tops$Feature, rownames(tmp)),
                               tops)

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
          definition = function(x,
                                contrastType=c("F", "Dendro", "Pairs", "OneAgainstAll"),
                                pairMat=NULL,
                                returnType=c("Table", "Index"),
                                contrastAdj=c("All", "PerContrast", "AfterF"),
                                isCount=FALSE, ...) {

            contrastType <- match.arg(contrastType)

            if(contrastType=="Dendro") {
              if(is.null(x@dendro_clusters)) {
                stop("If `contrastType='Dendro'`, `makeDendrogram` must be run before `getBestFeatures`")
              } else {
              dendro <- x@dendro_clusters
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
                         pairMat=pairMat, returnType=returnType,
                         contrastAdj=contrastAdj, isCount=isCount, ...)

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
.getBestDendroGenes<-function(cl,dendro,...){ #... is given above, not used by this function
	####
	#Convert to object used by phylobase so can navigate easily -- might should make generic function...
	# tempPhylo<-try(dendextend::as.phylo.dendrogram(dendro),FALSE)
	# if(inherits(tempPhylo, "try-error")) stop("the dendrogram object cannot be converted to a phylo class. Check that you gave simple hierarchy of clusters, and not one with fake data per sample")
	# # require(phylobase)
	# phylo4Obj<-try(as(tempPhylo,"phylo4"),FALSE)
	# if(inherits(phylo4Obj, "try-error")) stop("the created phylo object cannot be converted to a phylo4 class. Check that you gave simple hierarchy of clusters, and not one with fake data per sample")
	# phylobase::nodeLabels(phylo4Obj)<-paste("Node",1:phylobase::nNodes(phylo4Obj),sep="")
	phylo4Obj<-.makePhylobaseTree(dendro,type="dendro")

	clChar<-as.character(cl)
	allTipNames<-phylobase::labels(phylo4Obj)[phylobase::getNode(phylo4Obj,  type=c("tip"))]
	if(any(sort(allTipNames)!=sort(unique(clChar)))) stop("tip names of dendro don't match cluster vector values")

	#each internal node (including root) construct contrast between them.
	#(before just tested differences between them)
	allInternal<-phylobase::getNode(phylo4Obj,  type=c("internal"))
	.makeNodeContrast<-function(nodeId){
		children<-phylobase::descendants(phylo4Obj,nodeId,"children") #get immediate descendants
		if(length(children)!=2) stop("More than 2 children for internal node; does not make sense with code")
			#find tips of each child:
		contrAvePerChild<-sapply(children,function(x){
			tips<-phylobase::descendants(phylo4Obj,x,"tip")
			tipNames<-phylobase::labels(phylo4Obj)[tips]
			if(length(tipNames)>1) return(paste("(",paste(paste("X",tipNames,sep=""),collapse="+"),")/",length(tips),sep=""))
				else return(paste("X",tipNames,sep=""))
		})
		return(paste(contrAvePerChild,collapse="-"))
	}
	x<-sapply(allInternal,.makeNodeContrast)
	return(.testContrasts(contrastNames=x, ...))

}

# .testNode<-function(nodeId){
# 	children<-phylobase::descendants(phylo4Obj,nodeId,"children") #get immediate descendants
# 	if(length(children)!=2) stop("More than 2 children for internal node; does not make sense with code")
# 		#find tips of each child:
# 	contrAvePerChild<-sapply(children,function(x){
# 		tips<-phylobase::descendants(phylo4Obj,x,"tip")
# 		tipNames<-phylobase::labels(phylo4Obj)[tips]
# 		paste("(",paste(tipNames,sep="+",")/",length(tips),sep="")
# 	})
#
# 	groupId<-rep(NA,length=length(cl))
# 	sampleIndInEach<-lapply(children,function(x){
# 		tips<-phylobase::descendants(phylo4Obj,x,"tip")
# 		tipNames<-phylobase::labels(phylo4Obj)[tips]
# 		ind<-which(clChar %in% tipNames)
# 		nam<-phylobase::labels(phylo4Obj)[x]
# 		if(x %in% phylobase::getNode(phylo4Obj,type="internal")){
# 			nam<-paste(nam,":",paste(tipNames,collapse="_"),sep="")
# 		}
# 		if(x %in% phylobase::getNode(phylo4Obj,type="tip")){
# 			nam<-paste("Cluster",nam,sep="")
# 		}
# 		groupId[ind]<<-nam ##Assign them labels that correspond to the nodeId
# 		return(ind)
# 		})
# 	vals<-na.omit(unique(groupId))
# 	notIn<-(1:length(cl))[-unlist(sampleIndInEach)]
# 	if(length(notIn)>0){
# 		groupId[notIn]<-"-1"
# 		vals<-c(vals,"-1")
# 	}
# 	nameVals<-make.names(vals)
# 	groupId<-factor(groupId,labels=nameVals,levels=vals)
# 	if(any(is.na(groupId))) stop("Coding error, not all individuals were assigned")
#
# 	## basic limma design
# 	design<-model.matrix(~0+groupId)
# 	colnames(design)<-levels(groupId)
# 	fit<-lmFit(dat,design)
# 	cont<-paste(nameVals[1],nameVals[2],sep="-")
# 	cont.matrix<-makeContrasts(contrasts=cont,levels=fit$design)
# 	fit2<-contrasts.fit(fit,cont.matrix)
# 	fit2<-eBayes(fit2)
# 	tt<-topTable(fit=fit2,coef=1,number=nGenes,genelist=rownames(fit$coef),...)
# 	colnames(tt)[colnames(tt)=="ID"]<-"Feature"
# 	#pretty back up
# 	cont<-gsub("_",",",cont)
# 	cont<-gsub("[\\.]",":",cont)
# 	if(nrow(tt)>0) tt<-data.frame("Contrast"=cont,tt)
# 	return(tt)
#
# }
# #browser()
# tops<-do.call("rbind",lapply(allInternal,.testNode))#,tableArgs=list(...))
#.testNode(allInternal[2])

.getBestPairsGenes<-function(cl,pairMat,...){
	## basic limma design
  #browser()
	if(is.null(pairMat)){ #make pair Mat of all pairwise
		levs<-levels(cl)
		pairMat<-t(apply(expand.grid(levs,levs),1,sort))
		pairMat<-unique(pairMat)
		pairMat<-pairMat[-which(pairMat[,1]==pairMat[,2]),,drop=FALSE]
	}
	if(is.null(dim(pairMat)) || ncol(pairMat)!=2) stop("pairMat must be matrix of 2 columns")
	if(!all(as.character(unique(c(pairMat[,1],pairMat[,2]))) %in% as.character(cl))) stop("Some elements of pairMat do not match cl")
	x <- apply(pairMat,1,function(y){y<-make.names(as.character(y));paste(y[1],y[2],sep="-")})
	return(.testContrasts(contrastNames=x, ...))
}
.getBestOne<-function(cl,...){
	levs<-levels(cl)
	contrVect<-sapply(levs,function(x){
		one<-make.names(x)
		all<-make.names(levs[-which(levs==x)])
		all<-paste("(",paste(all,collapse="+"),")/",length(all),sep="")
		contr<-paste(all,"-",one,sep="")
		return(contr)
	})
	names(contrVect)<-levs
	return(.testContrasts(contrastNames=contrVect, ...))
}

#' @title Function for creating contrasts for a cluster
#' @description Uses cluster to create different types of contrasts to be tested that can then be fed into DE testing programs.
#' @rdname getBestFeatures
#' @aliases clusterContrasts
#' @param cluster Either a vector giving contrasts assignments or a ClusterExperiment object
#' @param contrastType What type of contrast to create.
#'   `Dendro' traverses the given dendrogram and does contrasts of the samples
#'   in each side,  `Pairs' does pair-wise contrasts based on the pairs given in
#'   pairMat (if pairMat=NULL, does all pairwise), and `OneAgainstAll' compares
#'   each cluster to the average of all others.
#' @param dendro The dendrogram to traverse if contrastType="Dendro". Note that this
#'   should be the dendrogram of the clusters, not of the individual samples.
#' @param pairMat matrix giving the pairs of clusters for which to do pair-wise
#'   contrasts (must match to elements of cl). If NULL, will do all pairwise of
#'   the clusters in \code{cluster} (excluding "-1" categories). Each row is a pair
#'   to be compared and must match the names of the clusters in the vector
#'   \code{cluster}.
#' @param removeNegative logical, whether to remove negative valued clusters 
#'   from the design matrix. Appropriate to pick TRUE (default) if design will
#'   be input into linear model on samples that excludes -1.
#' @details The input vector must be numeric clusters, but the external commands
#'   that make the contrast matrix (e.g. \code{\link{makeContrasts}}) require
#'   syntatically valid R names. For this reason, the names of the levels will
#'   be "X1" instead of "1". And negative values (if removeNegative=FALSE) will
#'   be "X.1","X.2", etc.
#' @return If \code{outputType=="limma"}, returns the results of running
#'   \code{\link{makeContrasts}}. This is a matrix with number of columns equal
#'   to the number of contrasts, and rows equal to the number of levels of the
#'   factor that will be fit in a linear model.
#' @examples 
#' data(simData)
#'
#' cl <- clusterMany(simData,nPCADims=c(5,10,50),  dimReduce="PCA",
#' clusterFunction="pam", ks=2:4, findBestK=c(FALSE), removeSil=TRUE,
#' subsample=FALSE)
#' #Pairs:
#' clusterContrasts(cl,contrastType="Pairs")
#' #Dendrogram
#' cl<-makeDendrogram(cl)
#' clusterContrasts(cl,contrastType="Pairs")

setMethod(f = "clusterContrasts",
          signature = "ClusterExperiment",
          definition = function(cluster,contrastType,...){
      if(contrastType=="Dendro"){
        if(is.null(cluster@dendro_clusters)) stop("Must run makeDendrogram before calling clusterContrasts if want contrastType='Dendro'")
        else dendro<-cluster@dendro_clusters
      }
      else dendro<-NULL
    clusterContrasts(primaryCluster(cluster),contrastType=contrastType,dendro=dendro,...)
})
setMethod(f = "clusterContrasts",
          signature = "numeric",
          definition = function(cluster,contrastType=c("Dendro", "Pairs", "OneAgainstAll"),
    dendro=NULL, pairMat=NULL,outputType="limma",removeNegative=TRUE){
              
    if(removeNegative) cl<-cluster[cluster>0] else cl<-cluster
    cl<-factor(cl)
    contrastType<-match.arg(contrastType)
    if(contrastType=="Dendro"){
        if(is.null(dendro)) stop("must provide dendrogram if contrastType='Dendro'")
        ####
        #Convert to object used by phylobase so can navigate easily -- might should make generic function...
        phylo4Obj<-.makePhylobaseTree(dendro,type="dendro")
        clChar<-as.character(cl)
        allTipNames<-phylobase::labels(phylo4Obj)[phylobase::getNode(phylo4Obj,  type=c("tip"))]
        if(any(sort(allTipNames)!=sort(unique(clChar)))) stop("tip names of dendro don't match cluster vector values")
        
        #each internal node (including root) construct contrast between them.
        #(before just tested differences between them)
        allInternal<-phylobase::getNode(phylo4Obj,  type=c("internal"))
        .makeNodeContrast<-function(nodeId){
            children<-phylobase::descendants(phylo4Obj,nodeId,"children") #get immediate descendants
            if(length(children)!=2) stop("More than 2 children for internal node; does not make sense with code")
            #find tips of each child:
            contrAvePerChild<-sapply(children,function(x){
                tips<-phylobase::descendants(phylo4Obj,x,"tip")
                tipNames<-phylobase::labels(phylo4Obj)[tips]
                #should make this code use make.names instead of pasting X...
                if(length(tipNames)>1) return(paste("(",paste(make.names(tipNames),collapse="+"),")/",length(tips),sep=""))
                else return(make.names(tipNames))
            })
            return(paste(contrAvePerChild,collapse="-"))
        }
        contrastNames<-sapply(allInternal,.makeNodeContrast)
        
    }
    if(contrastType=="OneAgainstAll"){
        levs<-levels(cl)
        contrastNames<-sapply(levs,function(x){
            one<-make.names(x)
            all<-make.names(levs[-which(levs==x)])
            all<-paste("(",paste(all,collapse="+"),")/",length(all),sep="")
            contr<-paste(all,"-",one,sep="")
            return(contr)
        })
        names(contrastNames)<-levs
    }
    if(contrastType=="Pairs"){
        if(is.null(pairMat)){ #make pair Mat of all pairwise
            levs<-levels(cl)
            pairMat<-t(apply(expand.grid(levs,levs),1,sort))
            pairMat<-unique(pairMat)
            pairMat<-pairMat[-which(pairMat[,1]==pairMat[,2]),,drop=FALSE]
        }
        if(is.null(dim(pairMat)) || ncol(pairMat)!=2) stop("pairMat must be matrix of 2 columns")
        if(!all(as.character(unique(c(pairMat[,1],pairMat[,2]))) %in% as.character(cl))) stop("Some elements of pairMat do not match cl")
        contrastNames <- apply(pairMat,1,function(y){y<-make.names(as.character(y));paste(y[1],y[2],sep="-")})
    }
#     if(!removeNegative){
#         levnames<-levels(cl)
#         whNeg<-which(cluster<0)
#         if(length(whNeg)>0){
#             levnames[whNeg]<-paste("Neg",cluster[whNeg],sep="")
#         }
#     }
#     if(removeNegative){
#         levnames<-levels(factor(cluster[cluster>0]))
#     }
    levnames<-make.names(levels(cl))
    if(outputType=="limma") return(makeContrasts(contrasts=contrastNames,levels=levnames))
    
})

.testContrasts<-function(contrastNames,fit,fitF,contrastAdj,...){
	ncontr<-length(contrastNames)
	cont.matrix<-makeContrasts(contrasts=contrastNames,levels=fit$design)
	fit2<-contrasts.fit(fit,cont.matrix)
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
	tops<-do.call("rbind",lapply(1:ncontr,function(ii){
		if(contrastAdj%in%c("AfterF","All")) {
			tt<-topTable(fit2,coef=ii, number=length(rownames(fit2$coef)),p.value=1,adjust.method="none",genelist=rownames(fit2$coef))
		}
		else{
			tt<-topTable(fit2,coef=ii, genelist=rownames(fit2$coef),...)
		}
		colnames(tt)[colnames(tt)=="ID"]<-"Feature"
		if(nrow(tt)>0){
			tt<-data.frame("Contrast"=unname(contrastNames[ii]),tt,row.names=NULL)
			if(!is.null(names(contrastNames))){
				 tt<-data.frame("ContrastName"=names(contrastNames)[ii],tt,row.names=NULL)
			}
		}
		return(tt)
		}))
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
