
##Make palette of colors


##Function to do heatmap based on fq, but hierarchical clustering based on correctedY
dualHeatmap<-function(clustVec,fq,correctedY,eps=1,dual=TRUE,clusterCells=TRUE, 
	clusterGenes=TRUE,whGenes=1:500,geneNames=FALSE,cellNames=FALSE,colorScale=seqPal5,
	annCol=NULL,annColors=NULL,alignColors=FALSE,breaks=NA,...){
	#dual=TRUE means use heatmap with fq, hierarch on correctedY; if FALSE both on fq
	#clustVec a vector giving clusters to color the samples with; if missing, then just do all white (NA)
	if(missing(clustVec)) clustVec<-rep(NA,length=ncol(fq)) 
	else{
		if(length(clustVec)!=ncol(fq)) stop("clustVec not of same length as ncol(fq)")
		clustVec[as.character(clustVec)== "-1" ]<- NA #those that don't want to be colored 
	}
	
	#fix up annCol:
	#not sure why this doesn't give back data.frame with factors
	#annCol<-apply(annCol,2,function(x){factor(x)})
	if(is.null(annCol)) annCol<-data.frame(Cluster=clustVec)
	tmpDf<-do.call("data.frame",lapply(1:ncol(annCol),function(ii){factor(annCol[,ii])}))
	names(tmpDf)<-names(annCol)
	annCol<-tmpDf
	if(is.null(annColors)){
		if(alignColors){
			#align the clusters and give them colors from .thisPal
			clMat<-t(data.matrix(annCol)) #converts them all to numbers, required for plotTracking
			#for each column of clDf, get colors from plotTracking
			clMat[clMat== -1]<-max(clMat)+1
			alignObj<-plotTracking(clMat+2,plot=FALSE) #in case any "-1"
			annColors<-mapply(alignObj$groupToColorLegend,annCol,FUN=function(x,fac){
				cols<-x[,"Color"]
				xnam<-levels(fac)[as.numeric(x[,"Original"])-2]
				# print(levels(fac))
	# 			print(as.numeric(x[,"Original"])-2)
				names(cols)<-xnam
				cols[names(cols)=="-1"]<-"white" #unassigned get white
				cols<-cols[order(names(cols))]
				return(cols)
				},SIMPLIFY=FALSE)
		}
		else{#make them have separate colors
			maxPerAnn<-sapply(annCol,function(x){max(as.numeric(x))})
			annColors<-mapply(annCol,c(0,head(cumsum(maxPerAnn),-1)),FUN=function(fac,add){
				cols<-.thisPal[1:nlevels(fac)+add]
				names(cols)<-levels(fac)
				cols[names(cols)=="-1"]<-"white" #unassigned get white
				cols<-cols[order(names(cols))]
				return(cols)
			},SIMPLIFY=FALSE)
		}
		if(ncol(annCol)==1){
			#annCol<-annCol[,1]
			#annColors<-annColors[[1]]
		}
	}
	
	if(dual & clusterCells){
		if(inherits(correctedY, "dendrogram")){
			if("nobs.dendrogram" %in% methods(nobs)){ #not all versions have nobs method for dendrograms; R3.2.0 has it; 3.1.1 doesn't
				if(nobs(correctedY)!=ncol(fq)) stop("correctedY dendrogram is not on same number of observations as fq")
			}
			else{warning("This R version doesn't doesn't allow for checking that dendrogram supplied has correct number observations. If not, surprising downstream errors can occur.")}
			dendroCells<-correctedY	
		} 
		else{
			if(ncol(correctedY)!=ncol(fq)) stop("correctedY matrix not have on same number of observations as fq")
			dendroCells<-as.dendrogram(hclust(dist(t(correctedY[whGenes,]))))
		}
	}
	tmp<-log(data.matrix(fq[whGenes,])+eps) 
	if(!is.logical(cellNames)){
		colnames(tmp)<-cellNames
		cellNames<-TRUE
	}
	if(!is.logical(geneNames)){
		rownames(tmp)<-geneNames
		geneNames<-TRUE
	}
	if(!clusterCells){ #then use clustVec to order them
		ord<-order(clustVec)
		tmp<-tmp[,ord,drop=FALSE]
		clustVec<-clustVec[ord]
		annCol<-annCol[ord,,drop=FALSE]
		clusterCells<-NA
#		if(clusterGenes) dendo<-"row" else dendo<-"none"
	}
	else{
#		if(clusterGenes) dendo<-"both" else dendo<-"col"
	}
	if(!clusterGenes) clusterGenes<-NA
	if(!is.factor(clustVec)) clustVec<-factor(clustVec)
	if(length(breaks)>0 || !is.na(breaks)){ #get arround bug in aheatmap
		#if colors are given, then get back 51, unless give RColorBrewer, in which case get 101! Assume user has to give palette.
		if(length(breaks)==1){
			if(breaks<=1){
				ncols<-51
				if(breaks<1) breaks<-c(seq(min(tmp),quantile(tmp[tmp>0],breaks,na.rm=TRUE),length=ncols),max(tmp))				
				else breaks<-seq(min(tmp),max(tmp),length=ncols+1)
			}
			else{
				warning("Because of bug in aheatmap, breaks should be of length 52 -- otherwise the entire spectrum will not be used. We don't recommend that you set the breaks to a integer number, but let aheatmap determine the breaks")
			}
		}
		else{
			if(length(breaks)!=52) warning("Because of bug in aheatmap, breaks should be of length 52 -- otherwise the entire spectrum will not be used")
		}
	}

	# gplots::heatmap.2(tmp, Colv=if(dual & clusterCells) dendroCells else clusterCells, Rowv=clusterGenes,dendrogram=dendo, scale="none",
	# 	trace="none", key=TRUE, density.info="none", col=colorScale, ColSideColor=col[clustVec],
	# 	cexRow=ifelse(nrow(tmp)<50,1,.6),
	# 	margin=c(ifelse(cellNames, 10,0.1), ifelse(geneNames, 10,0.1)),labRow=rownames(tmp),...)
	out<-NMF::aheatmap(tmp, color = colorScale, scale = "none", Rowv =clusterGenes, Colv = if(dual && !is.na(clusterCells) && clusterCells) dendroCells else clusterCells, 
		 annCol = annCol,annColors=annColors,breaks=breaks,...)
	invisible(list(heatOut=out,annCol=annCol,annColors=annColors,breaks=breaks))
}