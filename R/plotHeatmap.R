
#' Heatmap with two different sources for hierarchical clustering and color scale.
#' 
#' Make heatmap with color scale from one matrix and hiearchical clustering
#' from another. Also color palettes to go with heatmap
#' 
#' 
#' Note that plotHeatmap calles aheatmap under the hood. This allows you to
#' plot multiple heatmaps via par(mfrow=c(2,2)), etc. However, the dendrograms
#' do not resize if you change the size of your plot window in an interactive
#' session of R (this might be a problem for RStudio if you want to pop it out
#' into a large window...).
#' 
#' @name plotHeatmap
#' @aliases plotHeatmap
#' @docType methods
#' @param clusterVector A numeric vector with cluster assignments to show at top of heatmap
#' with cells. ``-1'' indicates the sample was not assigned to a cluster and
#' gets color `white'.
#' @param heatData matrix to define the color scale (i.e. the counts). Will be
#' converted to log-scale internally. Assumes samples on rows and variables on
#' columns.
#' @param clusterData Either a matrix define the hiearchical clustering of
#' samples (e.g. normalized data) or a dendrogram for clustering the samples
#' (only used if dual=TRUE).
#' @param eps Amount to add to heatData so can take log
#' @param dual Logical as to whether should use clusterData for dendrogram;
#' otherwise heatData is used.
#' @param clusterSamples Logical as to whether to do hierarchical clustering of
#' cells.
#' @param clusterVars Logical as to whether to do hiearchical clustering of
#' genes.
#' @param whVars Which genes of heatData matrix to be used. Default assumes top 500,
#' which may be nonsensical if matrix not ordered.
#' @param varNames Logical as to whether show gene names
#' @param sampleNames Logical as to whether show cell names
#' @param colorScale palette of colors for the color scale of heatmap
#' @param annCol data.frame of clusters to show at the top of heatmap (columns
#' are clusters). If NULL, will show clusterVector. Can also be continuous valued data, but see details.
#' @param clusterColors Assignment of colors to the clusters. If NULL, clusters
#' will be assigned colors. If `annCol' should be list of length equal to
#' ncol(annCol) with names equal to the colnames of annCol; each element of the
#' list should be a vector of colors with names corresponding to the levels of
#' the column of annCol. If annCol=NULL (i.e. use clusterVec), then the name of
#' the length-1 list should be `Cluster'
#' @param whMetaDataCont which columns of annCol are continuous data; only used if annCol=NULL (i.e. function \code{plotHeatmap} assigns the clusters' colors)
#' @param alignClusterColors Logical as to whether should align the clusters when
#' colors are assigned (only used if annCol=NULL)
#' @param breaks Either a vector of breaks (should be equal to length 52), or a
#' number between 0 and 1, indicating that the breaks should be equally spaced
#' (on the log scale+eps) upto the `breaks' quantile, see details
#' @param unassignedColor color assigned to cluster values of '-1' ("unassigned")
#' @param missingColor color assigned to cluster values of '-2' ("missing")
#' @param ... passed to aheatmap
#' @details The plotHeatmap function calles \code{\link{aheatmap}} to draw the heatmap. The main point of \code{plotHeatmap} is to 1) allow for two different matrix inputs, one to visualize and one to cluster. 
#' 2) to assign colors to the clusters like in \code{\link{plotClusters}} that lines them up based on their similarity. 
#' The intended purpose is to allow the user to visualize the original count scale of the data (on the log-scale), but create the hierarchical clustering on another, more appropriate dataset for clustering, such as normalized data. Similarly, some of the palettes were developed assuming that the visualization might be on unscaled/uncentered data, rather than the residual from the mean of the gene, and thus palettes need to take on a greater range of relevant values so as to show meaningful comparisons with genes on very different scales.  
#'
#' @details If \code{annCol} contains a column of continuous data, whMetaDataCont should give the index of the column(s); otherwise the annotation data for those columns will be forced into a non-sensical factor (with nlevels equal the the number of samples). 
#'
#' @details If breaks is a numeric value between 0 and 1, then \code{breaks} is assumed to indicate the upper quantile (on the log scale) at which the heatmap color scale should stop. For example, if breaks=0.9, then the breaks will evenly spaced up until the 0.9 upper quantile of the log of the \code{heatData}, and then all values after the 0.9 quantile will be absorbed by the upper-most color bin. This can help to reduce the visual impact of a few highly expressed genes (variables). 
#'
#' @return Returns (invisibly) a list with elements that are passed to aheatmap.
#' \itemize{
#' \item{\code{breaks}}{The breaks used for aheatmap, after adjusting for quantile}
#' \item{\code{annCol}}{the annotation data.frame given to aheatmap}
#' \item{\code{clusterColors}}{the annotation colors given to aheatmap}
#' }
#' @author Elizabeth Purdom
#' @examples
#' 
#' data(simData)
#' data(simData)
#' cl<-rep(1:3,each=100)
#' cl2<-cl
#' changeAssign<-sample(1:length(cl),80)
#' cl2[changeAssign]<-sample(cl[changeAssign])
#' 
#' #simple, minimal, example. Show counts, but cluster on underlying means
#' plotHeatmap(cl,heatData=simCount,clusterData=simData)
#' 
#' #assign cluster colors
#' colors<-bigPalette[20:23]
#' names(colors)<-1:3
#' plotHeatmap(cl,heatData=simCount,clusterData=simData,clusterColors=list(colors))
#' 
#' #show two different clusters
#' anno<-data.frame(cluster1=cl,cluster2=cl2)
#' out<-plotHeatmap(cl,heatData=simCount,clusterData=simData,annCol=anno)
#' #return the values to see format for giving colors to the annotations
#' out$clusterColors
#' 
#' #assign colors to the clusters based on plotClusters algorithm
#' plotHeatmap(cl,heatData=simCount,clusterData=simData,annCol=anno,
#' alignClusterColors=TRUE)
#' 
#' #assign colors manually
#' annoColors<-list(cluster1=c("black","red","green"),
#' cluster2=c("blue","purple","yellow"))
#' plotHeatmap(cl,heatData=simCount,clusterData=simData,annCol=anno,
#' clusterColors=annoColors)
#'
#' #give a continuous valued -- need to indicate columns
#' anno2<-cbind(anno,Cont=c(rnorm(100,0),rnorm(100,2),rnorm(100,3)))
#' plotHeatmap(cl,heatData=simCount,clusterData=simData,annCol=anno2,
#' whMetaDataCont=3)

#' #compare changing breaks quantile on visual effect
#' \dontrun{
#' par(mfrow=c(2,2))
#' plotHeatmap(cl,heatData=simCount,clusterData=simData,colorScale=seqPal1,
#' breaks=1,main="Full length")
#' plotHeatmap(cl,heatData=simCount,clusterData=simData,colorScale=seqPal1,
#' breaks=.99,main="0.99 Quantile Upper Limit")
#' plotHeatmap(cl,heatData=simCount,clusterData=simData,colorScale=seqPal1,
#' breaks=.95,main="0.95 Quantile Upper Limit")
#' plotHeatmap(cl,heatData=simCount,clusterData=simData,colorScale=seqPal1,
#' breaks=.90,main="0.90 Quantile Upper Limit")
#' }
#' 
#' 
#' 
setMethod(
  f = "plotHeatmap",
  signature = signature(x = "ClusterExperiment"),
  definition = function(clusterings,
                        visualize=c("original","transformed","centeredAndScaled"),
                        clusterSamplesWith=c("dendrogram","noClustering","mostVar","PCA"),
                        transformFeatures=c("all","mostVar"),
                        nFeaturesVisualize=500,
                        nFeaturesCluster=500,
                        isCount=isCount, transFun,
  ){	
    clusterOn<-match.arg(clusterOn)
    visualizeOn<-match.arg(visualizeOn)
    if("transformed" %in% c(visualizeOn)) heatData<-transform(x,transFun=transFun,isCount=isCount,dimReduce="none")
    else heatData<-assay(x)
    
    if(clusterOn=="mostVar"){
      clusterData<-transform(x,transFun=transFun,isCount=isCount,dimReduce="PCA",nPCADims=nFeatures)
    }
    else if(clusterOn=="PCA"){
      clusterData<-transform(x,transFun=transFun,isCount=isCount,dimReduce="PCA",nPCADims=nFeatures)
    }    
  })
#assume that the clusterData,heatData has already been reduced in diminsionality in some way.

setMethod(
  f = "plotHeatmap",
  signature = signature(x = "matrix"),
  definition = function(clusterings, heatData,clusterData=heatData, 
      metaData=NULL,whMetaDataCont=NULL,
      clusterSamples=TRUE,showSampleNames=FALSE, 
      clusterFeatures=TRUE,showFeatureNames=FALSE,
      colorScale=if(centerAndScaleFeatures) seqPal3 else seqPal5,
      clusterColors=NULL,alignClusterColors=FALSE,unassignedColor="white",missingColor="grey", breaks=NA,...
  ){	
  

	##########
	##Deal with clusterings...
	##########
	#check clusterings input:
  if(!is.matrix(clusterings) | !is.numeric(clusterings)){ stop("clusterings must be a numeric matrix of cluster values") 
  }
  if(NCOL(origX) != NROW(clusterings)) stop("clusterings must be matrix with same number of rows as columns of x")	
	
	#check metaData input
	
	#add metaData to clusterings
	if(!is.null(clusterings) && !is.null(metaData)){
	  if(!is.null(whMetaDataCont)) whAnnCont<-whMetaDataCont+NCOL(clusterings) else whAnnCont<-NULL
	  clusterings<-cbind(clusterings,metaData)
	} 
	if(is.null(clusterings) && !is.null(metaData)) clusterings<-metaData

	#not sure why this doesn't give back data.frame with factors: annCol<-apply(annCol,2,function(x){factor(x)})
	if(!is.null(clusterings)){
		###Make them explicitly factors
		tmpDf<-do.call("data.frame",lapply(1:ncol(clusterings),function(ii){factor(clusterings[,ii])}))
		names(tmpDf)<-names(clusterings)
		if(!is.null(whMetaDataCont)) tmpDf[,whMetaDataCont]<-metaData[,whMetaDataCont]
		annCol<-tmpDf

		if(is.null(clusterColors)){ #give default colors
			if(is.null(whAnnCont) || length(whAnnCont)<ncol(annCol)){
				if(!is.null(whMetaDataCont)) tmpDf<- annCol[,-whAnnCont,drop=FALSE] else tmpDf<-annCol
				if(alignClusterColors){
					#align the clusters and give them colors from .thisPal
					clMat<-data.matrix(tmpDf) #converts them all to numbers, required for plotClusters
					#for each column of clDf, get colors from plotClusters
					#clMat[clMat== -1]<-max(clMat)+1
					alignObj<-plotClusters(clMat,plot=FALSE) #in case any "-1"; probably don't need now
					
					#make a list of the color alignments
					clusterColors<-mapply(alignObj$groupToColorLegend,tmpDf,FUN=function(x,fac){
						cols<-x[,"Color"]
						xnam<-levels(fac)[as.numeric(x[,"Original"])-2]
						# print(levels(fac))
			# 			print(as.numeric(x[,"Original"])-2)
						names(cols)<-xnam
						cols[names(cols)=="-1"]<-unassignedColor #unassigned get white
						cols[names(cols)=="-2"]<-missingColor #unassigned get white
						cols<-cols[order(names(cols))]
						return(cols)
						},SIMPLIFY=FALSE)
				}
				else{#make them have distinct colors
					maxPerAnn<-sapply(tmpDf,function(x){max(as.numeric(x))})
					clusterColors<-mapply(tmpDf,c(0,head(cumsum(maxPerAnn),-1)),FUN=function(fac,add){
						cols<-.thisPal[1:nlevels(fac)+add]
						names(cols)<-levels(fac)
						cols[names(cols)=="-1"]<-unassignedColor #unassigned 
						cols[names(cols)=="-2"]<-missingColor #
						cols<-cols[order(names(cols))]
						return(cols)
					},SIMPLIFY=FALSE)
				}
			}
		}
	}
	else annCol<-NULL

	##########
	##Deal with numeric matrix for heatmap ...
	##########
	heatData<-data.matrix(heatData) 
	#	tmp<-data.matrix(heatData)
	# 	if(!is.logical(sampleNames)){
	# 		colnames(tmp)<-sampleNames
	# 		sampleNames<-TRUE
	# 	}
	# 	if(!is.logical(varNames)){
	# 		rownames(tmp)<-varNames
	# 		varNames<-TRUE
	# 	}
	###Create the clustering dendrogram:
	if(clusterSamples){
		if(inherits(clusterData, "dendrogram")){
			if("nobs.dendrogram" %in% methods(nobs)){ #not all versions have nobs method for dendrograms; R3.2.0 has it; 3.1.1 doesn't
				if(nobs(clusterData)!=ncol(heatData)) stop("clusterData dendrogram is not on same number of observations as heatData")
			}
			else{warning("This R version doesn't doesn't allow for checking that dendrogram supplied has correct number observations. If not, surprising downstream errors can occur.")}
			dendroCells<-clusterData	
		} 
		else{
		  clusterData<-data.matrix(clusterData)
		  #check valid
			if(ncol(clusterData)!=ncol(heatData)) stop("clusterData matrix not have on same number of observations as heatData")
			dendroCells<-as.dendrogram(hclust(dist(t(clusterData))))
		}
	}
	else{
	  if(!missing(clusterVector)){ #then use clusterVector to order them
	    ord<-order(clusterVector)
	    tmp<-tmp[,ord,drop=FALSE]
	    clusterVector<-clusterVector[ord]
	    annCol<-annCol[ord,,drop=FALSE]
	    clusterSamples<-NA
	  }
	}

	####More things to give to aheatmap
	if(!clusterFeatures) clusterFeatures<-NA
	if(length(breaks)>0 && !is.na(breaks)){ #get arround bug in aheatmap
		#if colors are given, then get back 51, unless give RColorBrewer, in which case get 101! Assume user has to give palette.
		#might not need any more with updated aheatmap.
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

	out<-NMF::aheatmap(tmp, color = colorScale, scale = "none", Rowv =clusterVars, Colv = if(dual && !is.na(clusterSamples) && clusterSamples) dendroCells else clusterSamples, 
		 annCol = annCol,clusterColors=clusterColors,breaks=breaks,...)
		 
	#############
	# add labels to clusters at top of heatmap
	#############
#	browser()
	if(!is.null(annCol)){
		newName<-NMF:::vplayout(NULL) #will be 1 greater (hopefully!) this is fragile. Don't know if it will always work.
		newNameList<-strsplit(newName,"\\.")[[1]]
		oldIndex<-as.numeric(newNameList[[3]])-1
		newNameList[[3]]<-oldIndex
		oldName<-paste(newNameList,collapse=".")
		grid::seekViewport(sprintf("aheatmap-%s",oldName))
		NMF:::vplayout(3,4:5)
		#grid::grid.rect()
		y <- seq(0,1,length=ncol(annCol))
		n<-ncol(annCol)
		y = cumsum(rep(8, n)) - 4 + cumsum(rep(2, n))
#		grid::grid.points(x = grid::unit(rep(0,length(y)),"npc"),y = grid::unit(y[n:1], "bigpts"))
		grid::grid.text(colnames(annCol), x = grid::unit(rep(0.05,length(y)),"npc"),y = grid::unit(y[n:1], "bigpts"), vjust = 0.5, hjust = 0,gp= grid::gpar(fontsize=10))
		grid::upViewport() #close it
		grid::upViewport() #close it
	}

	invisible(list(heatOut=out,annCol=annCol,clusterColors=clusterColors,breaks=breaks))
}

