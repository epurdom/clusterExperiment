# getBestFGenes<-function(cl,dat,nGenes=100,correctedY=correctedY){
#         tmp<-dat
#         if(any(cl== -1)){ #only use those assigned to a cluster to get good genes.
#                 whNA<-which(cl== -1)
#                 tmp<-tmp[, -whNA]
#                 cl<-cl[-whNA]
#         }
#         ## F-statistic
#         f <- apply(tmp, 1, function(y) {
#             fit <- lm(y ~ as.factor(cl))
#             summary(fit)$fstatistic[1]
#         })
#         f <- sort(f, decreasing=TRUE)
#         bestGenes<-names(f)[1:nGenes]
#         whGenes<-match(bestGenes,rownames(correctedY))
#         names(whGenes)<-bestGenes
#         return(whGenes)
# }


#' @title Function for finding best genes associated with clusters
#' @description Calls limma on input data to determine genes most associated with found clusters (either pairwise comparisons, or following a tree that clusters the clusters).
#' @param cl A vector with cluster assignments to compare to clRef.``-1'' indicates the sample was not assigned to a cluster.
#' @param dat data for the test, samples in rows
#' @param type What type of test to do. `F' gives the omnibus F-statistic, `Dendro' traverses the given dendrogram and does contrasts of the samples in each side, and `Pairs' does pair-wise contrasts based on the pairs given in pairMat.
#' @param dendro The dendrogram to traverse if type="Dendro". Note that this should be the dendrogram of the clusters, not of the individual samples.
#' @param pairMat matrix giving the pairs of clusters for which to do pair-wise contrasts (must match to elements of cl). If NULL, will do all pairwise of the clusters in \code{cl} (excluding "-1" categories). Each row is a pair to be compared and must match the names of the clusters in the vector \code{cl}. 
#' @param returnType Whether to return the index of genes, or the full table given by topTable or topTableF.
#' @param contrastAdj What type of FDR correction to do for contrasts tests (i.e. type='Dendro' or 'Pairs'). 
#' @param ... options to pass to \code{\link{topTable}} or \code{\link{topTableF}} (see \code{\link{limma}} package)
#' @details getBestGenes returns best genes corresponding to a cluster assignment. It uses limma to fit the models, and limma's functions topTable or topTableF to find the best genes.  See the options of these functions to put better control on what gets returned (e.g. only if significant, only if log-fc is above a certain amount, etc.). In particular, set 'number=' to define how many significant genes to return (where number is per contrast for the 'Pairs' or 'Dendro' option)
#'
#' @details When 'type' argument implies that the best genes should be found via contrasts (i.e. 'type' is 'Dendro' or 'Pairs'), then then 'contrastAdj' determines the type of multiple testing correction to perform. 
#' 'PerContrast' does FDR correction for each set of contrasts, and does not guarantee control across all the different contrasts (so probably not the preferred method). 
#' 'All' calculates the corrected p-values based on FDR correction of all of the contrasts tested. 'AfterF' controls the FDR based on a hierarchical scheme that only tests the contrasts in those genes where the omnibus F statistic is significant.  If the user selects 'AfterF', the user must also supply an option 'p.value' to have any effect, and then only those significant at that p.value level will be returned. 
#'
#' @details  Note that the default option for \code{\link{topTable}} is to not filter based on adjusted p-values (\code{p.value=1}) and return only the top 10 most significant (\code{number=10}) -- these are options the user can change (these arguments are passed via the \code{...} in getBestGenes). In particular, it only makes sense to set \code{requireF=TRUE} if \code{p.value} is meaningful (e.g. 0.1 or 0.05); the default value of \code{p.value=1} will not result in any effect on the adjusted p-value otherwise. 
#'
#' @return A data.frame in the same format as \code{\link{topTable}}, except that the column name \code{ProbeID} is changed to \code{Gene} and a column \code{IndexInOriginal} is provided to give the row index of the gene in the original dataset.

getBestGenes<-function(cl,dat,type=c("F","Dendro","Pairs"),dendro=NULL,pairMat=NULL,returnType=c("Index","Table"),contrastAdj=c("PerContrast","AfterF","All"),...){
	#... is always sent to topTable, and nothing else
#	require(limma)
	dat<-t(data.matrix(dat))
	type<-match.arg(type)
	contrastAdj<-match.arg(contrastAdj)
	returnType<-match.arg(returnType)
	if(is.null(rownames(dat))) rownames(dat)<-paste("Row",as.character(1:nrow(dat)),sep="")
	tmp<-dat
	if(any(cl== -1)){ #only use those assigned to a cluster to get good genes. 
		whNA<-which(cl== -1)
		tmp<-tmp[, -whNA]
		cl<-cl[-whNA]
	}
	cl<-factor(cl)
	if(type=="Dendro"){
		if(is.null(dendro)) stop("must provide dendro")
		if(!inherits(dendro,"dendrogram")) stop("dendro must be of class 'dendrogram'")
	}
	if(type%in% c("Pairs","Dendro")){
	  designContr<-model.matrix(~0+cl)
	  colnames(design)<-make.names(levels(cl))
	  fitContr<-limma::lmFit(tmp,designContr)
	  
	}
	if(type=="F" || contrastAdj=="AfterF"){
	  designF<-model.matrix(~cl)
	  fitF<-limma::lmFit(tmp,designF)
	}
	else fitF<-NULL
	tops<-switch(type,"F"=.getBestFGenes(fitF,...),
		"Dendro"=.getBestDendroGenes(cl=cl,dendro=dendro,contrastAdj=contrastAdj,fit=fitContr,fitF=fitF,...),
		"Pairs"=.getBestPairsGenes(cl=cl,pairMat=pairMat,contrastAdj=contrastAdj,fit=fitContr,fitF=fitF,...)
		)
	tops<-data.frame(IndexInOriginal=match(tops$Gene,rownames(tmp)),tops)
	if(returnType=="Index"){
		whGenes<-tops$IndexInOriginal
		names(whGenes)<-tops$Gene
		return(whGenes)
	}
	if(returnType=="Table"){
		return(tops)
	}
}
.getBestFGenes<-function(fit,...){
	## basic limma design
	fit2<-limma::eBayes(fit)
	#tops<-topTableF(fit2,number=nGenes,genelist=rownames(fit$coef),...)
	tops<-limma::topTable(fit2,genelist=rownames(fit$coef),...)
	colnames(tops)[colnames(tops)=="ProbeID"]<-"Gene"

	return(tops)
}
.getBestDendroGenes<-function(cl,dendro,...){ #... is given above, not used by this function
	####
	#Convert to object used by phylobase so can navigate easily -- might should make generic function...
	tempPhylo<-try(dendextend::as.phylo.dendrogram(dendro),FALSE)
	if(inherits(tempPhylo, "try-error")) stop("the dendrogram object cannot be converted to a phylo class. Check that you gave simple hierarchy of clusters, and not one with fake data per sample")
	# require(phylobase)
	phylo4Obj<-try(as(tempPhylo,"phylo4"),FALSE) 
	if(inherits(phylo4Obj, "try-error")) stop("the created phylo object cannot be converted to a phylo4 class. Check that you gave simple hierarchy of clusters, and not one with fake data per sample")
	phylobase::nodeLabels(phylo4Obj)<-paste("Node",1:phylobase::nNodes(phylo4Obj),sep="")
	
	
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
# 	fit<-limma::lmFit(dat,design)
# 	cont<-paste(nameVals[1],nameVals[2],sep="-")
# 	cont.matrix<-limma::makeContrasts(contrasts=cont,levels=fit$design)
# 	fit2<-limma::contrasts.fit(fit,cont.matrix)
# 	fit2<-limma::eBayes(fit2)
# 	tt<-limma::topTable(fit=fit2,coef=1,number=nGenes,genelist=rownames(fit$coef),...)
# 	colnames(tt)[colnames(tt)=="ID"]<-"Gene"
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
	if(is.null(pairMat)){ #make pair Mat of all pairwise
		levs<-levels(cl)
		pairMat<-t(apply(expand.grid(levs,levs),1,sort))
		pairMat<-unique(pairMat)
		pairMat<-pairMat[-which(pairMat[,1]==pairMat[,2]),]
	}
	if(is.null(dim(pairMat)) || ncol(pairMat)!=2) stop("pairMat must be matrix of 2 columns")
	if(!all(as.character(unique(c(pairMat[,1],pairMat[,2]))) %in% as.character(cl))) stop("Some elements of pairMat do not match cl")
	x <- apply(pairMat,1,function(y){y<-make.names(as.character(y));paste(y[1],y[2],sep="-")})
	return(.testContrasts(contrastNames=x, ...))
}

.testContrasts<-function(contrastNames,fit,fitF,contrastAdj,...){
	ncontr<-length(contrastNames)
	cont.matrix<-limma::makeContrasts(contrasts=contrastNames,levels=fit$design)
	fit2<-limma::contrasts.fit(fit,cont.matrix)
	fit2<-limma::eBayes(fit2) 
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
			tt<-limma::topTable(fit2,coef=ii, number=length(rownames(fit2$coef)),p.value=1,adjust.method="none",genelist=rownames(fit2$coef)) 			
		}
		else{
			tt<-limma::topTable(fit2,coef=ii, genelist=rownames(fit2$coef),...) 			
		}
		colnames(tt)[colnames(tt)=="ID"]<-"Gene"
		if(nrow(tt)>0){
			tt<-data.frame("Contrast"=unname(contrastNames[ii]),tt,row.names=NULL)
			if(!is.null(names(contrastNames))) tt<-data.frame("ContrastName"=names(contrastNames)[ii],tt,row.names=NULL)
		}
		return(tt)
		}))
	if(contrastAdj=="AfterF" & p.value<1){
		#get p-value for F test for all genes, and only consider genes with significant F. 
	  fitF2<-limma::eBayes(fitF)
		topsF<-limma::topTable(fitF2,genelist=rownames(fit$coef),number=length(rownames(fit$coef)),adjust.method="BH")
		whGenesSigF<-topsF$ProbeID[which(topsF$adj.P.Val < p.value)]
		tops<-tops[tops$Gene %in% whGenesSigF,]
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

	return(tops)
	
}