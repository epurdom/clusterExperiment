#' Functions for Visualizing assigments of samples to clusters across multiple clusterings
#' 
#' Align multiple clusters of the same set of samples and provide a color-coded plot of their shared cluster assignments
#' 
#' @name plotTracking
#' @aliases plotTracking 
#' @docType methods

#' @param clusters A matrix of with each column corresponding to a clustering and
#' each row a sample. The plotting will plot them in order of this matrix,
#' and their order influences the plot greatly.
#' @param orderClusters A predefined order for the clusterings in the plot
#' @param orderSamples A predefined order in which the samples will be plotted. Otherwise the order will be found internally by aligning the clusters (assuming \code{input="clusters"})
#' @param reuseColors Logical. Whether each row should consist of the same set
#' of colors. By default (FALSE) each cluster that the algorithm doesn't
#' identify to the previous rows clusters gets a new color
#' @param matchToTop Logical as to whether all clusters should be aligned to
#' the first row. By default (FALSE) each cluster is aligned to the ordered
#' clusters of the row above it.
#' @param plot Logical as to whether a plot should be produced.
#' @param unassignedColor If ``-1'' in clusters, will be given this color (meant for samples not assigned to cluster).
#' @param missingColor If ``-2'' in clusters, will be given this color (meant for samples that were missing from the clustering, mainly when comparing clusterings run on different sets of samples)
#' @param minRequireColor In aligning colors between rows of clusters, require
#' this percent overlap
#' @param startNewColors logical, indicating whether in aligning colors between rows of clusters, 
#' should the colors restart at beginning of colPalette as long as colors are not in immediately proceeding row
#' (some of the colors at the end of bigPalette are a bit wonky, and so if you
#' have a large clusters matrix, this can be useful)
#' @param colPalette a vector of colors used for the different clusters. Must be as long as the maximum number of clusters found in any single clustering/column given in \code{clusters} or will otherwise return an error.
#' @param input indicate whether the input matrix is matrix of integer assigned clusters, or contains the colors. If \code{input="colors"}, then the object \code{clusters} is a matrix of colors and there is no alignment (this option allows the user to manually adjust the colors and replot, for example)
#' @param clNames names to go with the columns (clusterings) in matrix colorMat
#' @param add whether to add to existing plot
#' @param xCoord values on x-axis at which to plot the rows (samples) 
#' @param ylim vector of limits of y-axis
#' @param tick logical, whether to draw ticks on x-axis for each sample.
#' @param ylab character string for the label of y-axis 
#' @param xlab character string for the label of x-axis 
#' @param axisLine the number of lines in the axis labels on y-axis should be (passed to line=... in the axis call)
#' @param box logical, whether to draw box around the plot
#' @param ... for \code{plotTracking} arguments passed to \code{\link{plot}} if \code{add=FALSE}.


#' @return plotTracking returns (invisibly) the orders and other things that
#' go into making the matrix. Specifically, a list with the following elements
#' \itemize{
#'
#' \item{\code{index}}{ a vector of length equal to ncols(clusters) giving the order of the columns to use to get the original clusters matrix into the order made by plotTracking}
#'
#' \item{\code{colors}}{ matrix of color assignments for each element of original clusters matrix. Matrix is in the same order as original clusters matrix. The matrix \code{colors[index,]} is the matrix that can be given back to plotTracking to recreate the plot (see examples)}
#'
#' \item{\code{aligned}}{ matrix of integer valued cluster assignments that match the colors. This is useful if you want to have cluster identification numbers that are better aligned than that given in the original clusters. Again, the matrix is in same order as original input matrix}
#'
#' \item{\code{groupToColorLegend}}{ list of length equal to the number of columns of input matrix. The elements of the list are matrices, each with three columns named "Original","Aligned", and "Color" giving, respectively, the correspondence between the original cluster ids in \code{clusters}, the aligned cluster ids in \code{aligned}, and the color. }
#' }
#'
#' @author Elizabeth Purdom and Marla Johnson
#' 
#' @seealso consensusClustering package
#' @examples
#' #clustering using pam: try using different dimensions of pca and different k
#' ps<-c(5,10,50)
#' names(ps)<-paste("npc=",ps,sep="")
#' pcaData<-stats::prcomp(simData, center=TRUE, scale=TRUE)
#' cl <- compareChoices(lapply(ps,function(p){pcaData$x[,1:p]}),
#' clusterMethod="pam",ks=2:4,findBestK=c(TRUE,FALSE),subsampleArgs=list("k"=3))
#' colnames(cl$clMat) 
#' #make names shorter for plotting
#' colnames(cl$clMat)<-gsub("TRUE","T",colnames(cl$clMat))
#' colnames(cl$clMat)<-gsub("FALSE","F",colnames(cl$clMat))
#' colnames(cl$clMat)<-gsub("k=NA,","",colnames(cl$clMat))
#' par(mar=c(2,10,1,1))
#' out<-plotTracking(cl$clMat,axisLine=-2)
#' out$groupToColorLegend[1:2]
#' head(out$color[out$index,1:2])
#' 
#' #We can also change the order of the rows. Notice how this dramatically changes the plot!
#' plotTracking(cl$clMat[,c(3:6,1:2,7:ncol(cl$clMat))],axisLine=-2)
#'
#' #notice that the blue and orange are really arbitrarily different 
#' #colors because not next to each other in row.
#' #We can manually change the blue to orange :
#' #first find out color names by looking at the out$colors (but in right order using out$index)
#' out$colors[out$index[1:10],c("npc=5,findBestK=T","npc=5,k=3,findBestK=F")]
#' #change "#1F78B4" to "#FF7F00"
#' newColorMat<-out$colors
#' newColorMat[newColorMat=="#1F78B4"]<-"#FF7F00"
#' #replot by setting 'input="colors"'
#' plotTracking(newColorMat[out$index,],input="colors")
#'
#' #set some arbitrarily to "-1", meaning not clustered.
#' clNew<-apply(cl$clMat,2,function(x){
#'	wh<-sample(1:nrow(cl$clMat),size=10)
#'	x[wh]<- -1
#'	return(x)
#'	})
#' plotTracking(clNew)

#' @rdname ="plotTracking"
setMethod(
  f = "plotTracking",
  signature = signature(clusters = "ClusterCells",whClusters="character"),
  definition = function(clusters, whClusters,...)
  {
    #To Do: should add something for adding additional metadata, perhaps from SummarizedExperiment or user-entered
    plotTracking(allClusters(clusters),...)
  })
#' @rdname ="plotTracking"
setMethod(
  f = "plotTracking",
  signature = signature(clusters = "ClusterCells",whClusters="numeric"),
  definition = function(clusters, whClusters,...)
  {
    #To Do: should add something for adding additional metadata, perhaps from SummarizedExperiment or user-entered
    plotTracking(allClusters(clusters)[,whClusters],...)
  })
#' @rdname ="plotTracking"
setMethod(
  f = "plotTracking",
  signature = signature(clusters = "matrix"),
  definition = function(clusters, orderClusters=NULL,
              orderSamples=NULL,reuseColors=FALSE,matchToTop=FALSE,
              plot=TRUE,unassignedColor="white",missingColor="grey",
              minRequireColor=0.3,startNewColors=FALSE,colPalette=bigPalette,
              input=c("clusters","colors"),clNames=colnames(clusters),
              add=FALSE,xCoord=NULL,ylim=NULL,tick=FALSE,ylab="",xlab="",
              axisLine=0,box=FALSE,...)
{
  #reorder so
  if(is.null(orderSamples)) orderSamples<-1:nrow(clusters)
  else{ if(!all(orderSamples %in% 1:nrow(clusters))) stop("invalid values for orderSamples")}
  if(is.null(orderClusters)) orderClusters<-1:ncol(clusters)
  else{if(!all(orderClusters %in% 1:ncol(clusters))) stop("invalid values for orderClusters")}

  index<-orderSamples #match to old arguments
  clusters<-clusters[,orderClusters] #reorder clusters.
                                                   
  dnames<-dimnames(clusters)
	input<-match.arg(input)

	#arguments to be passed at various calls
	clusterPlotArgs<-list(clNames=clNames,add=add,xCoord=xCoord,ylim=ylim,tick=tick,ylab=ylab,xlab=xlab,axisLine=axisLine,box=box)
	plotTrackArgs<-list(index=index,reuseColors=reuseColors,matchToTop=matchToTop,minRequireColor=minRequireColor,startNewColors=startNewColors,colPalette=colPalette)

	if(input=="clusters"){
		clusters<-t(clusters)
		if(any(as.character(clusters)%in%c("-1","-2"))){
			if(any(apply(clusters,1,function(x){any(is.na(x))}))) stop("clusters should not have 'NA' values; non-clustered samples should get a '-1' or '-2' value depending on why they are not clustered.")
			out<-do.call(".plotTrackingInternal",c(list(clusters=clusters,plot=FALSE),plotTrackArgs,clusterPlotArgs ,list(...)))
			#take out -1
			newColorLeg<-lapply(1:nrow(clusters),function(i){
				leg<-out$groupToColorLegend[[i]]
				if(any(wh<-leg[,"Original"]== -1))
				leg[wh,"Color"]<-unassignedColor
				if(any(wh<-leg[,"Original"]== -2))
				leg[wh,"Color"]<-missingColor
				return(leg)
			})
			names(newColorLeg)<-names(out$groupToColorLegend)
			xColors<-do.call("cbind",lapply(1:nrow(clusters),function(i){
				out$colors[clusters[i,]=="-1",i]<-unassignedColor
				out$colors[clusters[i,]=="-2",i]<-missingColor
				return(out$colors[,i])
			}))
			if(plot) do.call(".clusterTrackingPlot",c(list(colorMat=xColors[out$index,]),clusterPlotArgs,list(...)))
			out$colors<-xColors
			dimnames(out$colors)<-dnames
			out$groupToColorLegend<-newColorLeg
		}
		else{
			out<-do.call(".plotTrackingInternal",c(list(clusters=clusters,plot=plot),plotTrackArgs,clusterPlotArgs ,list(...)))
		}	
		invisible(out)
	}
	else{
		do.call(".clusterTrackingPlot",c(list(colorMat=clusters),clusterPlotArgs,list(...)))		
	}
})

.plotTrackingInternal<-function(clusters, index, reuseColors, matchToTop, plot, minRequireColor, startNewColors, colPalette, ...){
	#clusters is a nmethods x nsamples matrix. The order of the rows determines how the tracking will be done.
	#if given, index is the order of the samples (columns) for plotting. Otherwise determined by program
	#matchToTop allows that the comparison for color assignment be done, not row by row, but always back to the first row.
	##Check that no row has more than length(colPalette) number of unique entries
	tabVals<-apply(clusters,1,function(x){length(table(x))})
	if(any(tabVals> length(colPalette))) stop("Must give colPalette longer than number of clusters. Largest number of distinct clusters in the clusterings:",max(tabVals),". Length of colPalette:",length(colPalette))
	
	pastColorVector<-NULL
	colorM = rbind() #matrix of colors. gives indexes for cluster assignments that give the 'right' color throughout the different rows

	for(i in 1:nrow(clusters)){ 
		#.setClusterColors <- function(past_ct,ct,colorU,colorList,reuseColors,minRequireColor=0.5){
			#description: sets common color of clusters between different K
			#takes previous row's colors and defines color assignments for next row. 
			#assignments are done by:
				#start with largest cluster, assign to old cluster with largest overlap
			#past_ct is the past (i.e. previous row's) assigment to indices;
			#ct is the new cluster ids that need to be aligned to the past_ct
			#colorU is set of possible colors to take from
			#colorList, if given, is three element list giving past assignment
			    #element 1 is colors for each id 
			    #element 2 is highest color index used
			    #element 3 is unique color names; not clear if ever used
			#returns a (updated) colorList in the above format

	#eap: rewrote a little, to make it easier to read. 5/27/2015
		currCl<-as.character(unlist(clusters[i,]))#as.numeric(as.character(unlist(clusters[i,]))) #think I've fixed it so don't need to convert to numeric. 
		if(i == 1) pastCl<-NULL
		else{
			if(!matchToTop) pastCl<-as.character(unlist(clusters[i-1,]))#as.numeric(as.character(unlist(clusters[i-1,])))
			else pastCl<-as.character(unlist(clusters[1,]))#as.numeric(as.character(unlist(clusters[1,])))
		}
		#if(i==2) browser()
		newColorVector<- .setClusterColors(past_ct=pastCl,ct=currCl ,colorU=colPalette,past_colors=pastColorVector,reuseColors=reuseColors,minRequireColor=minRequireColor,startNewColors=startNewColors)
		if(i == 1) colorListTop<-newColorVector
		if(any(is.na(newColorVector))) stop("NA values")
		if(any(sort(table(newColorVector))!=sort(table(unlist(clusters[i,]))))) stop("Coding error, colors to do not have same number as original")
		crossTab<-paste(newColorVector,unlist(clusters[i,]),sep=";")
		if(any(sort(table(crossTab))!=sort(table(unlist(clusters[i,]))))) stop("Coding error, colors to do not have same assignments as original")
		colorM = rbind(colorM,newColorVector) 
		if(!matchToTop) pastColorVector<-newColorVector
			else pastColorVector<-colorM[1,]
	}	
	
	if(is.null(index)){
		tmp<-lapply(1:nrow(colorM),function(i){unlist(colorM[i,])})
		index<-do.call("order",tmp)
	}

	if(plot) .clusterTrackingPlot(t(colorM[,index]),...)
	dimnames(colorM)<-dimnames(clusters)
	allColors<-unique(as.vector(colorM))
	alignCl<-apply(colorM,2,function(x){match(x,allColors)})
	dimnames(alignCl)<-dimnames(clusters)
	groupToColorLegend<-lapply(1:nrow(clusters),function(ii){
		mat<-cbind("Original"=unlist(clusters[ii,]),"Aligned"=unlist(alignCl[ii,]),"Color"=unlist(colorM[ii,]))
		rownames(mat)<-NULL
		(unique(mat))
	})
	names(groupToColorLegend)<-rownames(clusters)
	invisible(list(orderClusters=index,colors=t(colorM),aligned=t(alignCl),groupToColorLegend=groupToColorLegend))

}

.setClusterColors <- function(past_ct,ct,colorU,past_colors,reuseColors,minRequireColor,startNewColors){
	#description: sets common color of clusters between different K
	#takes previous row's colors and defines color assignments for next row. 
	#assignments are done by:
		#start with largest cluster, assign to old cluster with largest overlap
	#past_ct is the past (i.e. previous row's) assigment to indices;
	#ct is the new cluster ids that need to be aligned to the past_ct
	#colorU is set of possible colors to take from
	#colorList, if given, is three element list giving past assignment
	    #element 1 is colors for each id 
	    #element 2 is highest color index used
	    #element 3 is unique color names; not clear if ever used
	#returns a (updated) colorList in the above format
	
	
	colorU<-rep(colorU,30) #so never run out of colors
	if(is.null(past_colors) | is.null(past_ct)){
		#map values of ct to 1:(# unique values)
		ncl<-length(unique(ct))
		v<-1:ncl
		names(v)<-sort(unique(ct))
		clNum<-v[match(ct,names(v))]
		newColors = colorU[clNum]
	}
	else{
		#############
		#check that past_ct and past_colors line up
		#############
		pastTab<-table(past_ct,past_colors) 
		if(nrow(pastTab)!=ncol(pastTab)) stop("past_ct and past_colors have different numbers of categories")
		pastTab<-t(apply(pastTab,1,sort,decreasing=TRUE)) 
		#each row should have only 1 non-zero category, so once order by row, then should all be zero except first
		if(any(colSums(pastTab[,-1,drop=FALSE])>0)) stop("past_ct and past_colors do not imply the same clustering of samples.")

		if(!startNewColors){
			whShared<-which(colorU %in% past_colors)
			m<-max(whShared)
			colorU<-c(colorU[m:length(colorU)],colorU[1:m]) #puts the used by previous rows colors at the end
		}
		colorU<-colorU[!colorU %in% past_colors]
		colorU<-unique(colorU) #make sure unique colors, at least per cluster
		if(length(colorU)==0) stop("not enough colors given -- all are in past_colors")

		colori<-0 #index which color will give next.
		newColors = rep(NA,length=length(past_colors)) #NULL #NA?
		
		####Tabulate the two clusters
		mo=table(past_ct,ct)
		#m=mo/apply(mo,1,sum) #for each old cluster, tells the percentage of the new cluster. Shouldn't it be the other way around?(EAP)
		m=mo #use raw number instead of percentage of cluster.
						# Browse[4]> m
						#        ct
						# past_ct  1  2  3  4  9 10
						#       1 18  7  0  0  0  0
						#       3  4  0 67 45  0  0
						#       9  0  0  1  0 21  1
		
		#EAP: change so do the one with the most overlap first, simplifies code.
		orderClusters<-order(apply(m,2,max),decreasing = TRUE) 
		usedClusters<-c() #index of m of used clusters (pci) if !reuseColors; if reuseColors, used slightly differently
		for(tci in orderClusters){ # for each cluster EAP: in order of size
			
			if(!reuseColors){					
 			   #if(tci==1) browser()
				#pci tells which old cluster has the greatest number of the new cluster
				vals<-unique(m[,tci])
				vals<-vals[vals>0] #remove zeros
				pvals<-vals/sum(m[,tci])
				maxvals<-max(pvals)
				vals<-vals[pvals >= min(c(minRequireColor,maxvals))]#require 30% overlap to take on new cluster value
				if(length(vals)>0){
					vals<-sort(vals,decreasing=TRUE)
					#for each value, determine its matching cluster (if any) then take highest
					matchPci<-sapply(vals,function(v){
						whMatch<-which(m[,tci]==v) 
						whAv<-sapply(whMatch,function(x){!x %in% usedClusters})
						if(any(whAv)) return(whMatch[which(whAv)][1]) #pick the first of the available ones
						else return(NA)
					})
					newColor<-all(is.na(matchPci))
				}
				else{ newColor<-TRUE}
				if( newColor){ #new color
					colori=colori+1
					if(colori > length(colorU)) stop("not enough colors in color list (colorU)")
					newValue<-colorU[colori]	
				}
				else{
					pci <- matchPci[!is.na(matchPci)][1] #take that with highest percentage that matches
					newValue<-unique(past_colors[which(as.character(past_ct)==rownames(m)[pci])])
					usedClusters<-c(usedClusters,pci)
				}
			}
			else{
				#pci tells which old cluster has the greatest number of the new cluster
				mCurr<-m
				if(length(usedClusters)>=nrow(m)){
					colori=colori+1
					newValue<-colorU[colori]
				}
				else{
					if(length(usedClusters)>0 & length(usedClusters)<nrow(m)){
						mCurr<-m[-usedClusters, , drop=FALSE]
					}
					maxC = max(mCurr[,tci])						
					pci = row.names(mCurr)[which(mCurr[,tci] == maxC)]
					newValue<-unique(past_colors[which(as.character(past_ct)==pci)])
				}
				usedClusters<-c(usedClusters,grep(pci,row.names(m)))
			}
			newColors[which(as.character(ct)==colnames(m)[tci])] <-  newValue #EAP: fixed so no longer assumes cluster assignments were numeric integers with no gaps...
			tab<-table(newColors,ct)
			nPerColor<-apply(tab,1,function(x){sum(x!=0)})
			if(any(nPerColor!=1)) stop("color assigned twice within group.")
		}
	}
	if(any(is.na(newColors))) stop("Coding Error: some samples not given a color")
	tabOld<-table(ct)
	tabNew<-table(newColors)
	
	if(length(tabNew)!=length(tabOld) || any(sort(tabNew)!=sort(tabOld))) stop("Coding error, didn't get same number of entries of each")
	return(newColors)
}

.clusterTrackingPlot <- function(colorMat, clNames, add, xCoord, ylim, tick, ylab, xlab, axisLine, box, ...){
  	m<-t(colorMat)
  #description: plots cluster tracking plot
  #input: m - matrix where rows are k, columns are samples, and entries are color assignments
  #names: names that correspond with rows of m
	if(is.null(xCoord)){
		xleft<-seq(0,1-1/ncol(m),by=1/ncol(m))
		xright<-seq(1/ncol(m),1,by=1/ncol(m))
	}
	else{
		d<-unique(round(diff(xCoord),6)) #distance between x values; had to round because got differences due to small numerical storage issues
		if(length(d)!=1) stop("invalid x values given -- not evenly spaced")
		xleft<-xCoord-d/2
		xright<-xCoord+d/2
	}
	if(is.null(ylim)){
		ylim<-c(0,1)
	}
	d<-(ylim[2]-ylim[1])/nrow(m)
	if(nrow(m)>1){ #for safety. theoretically works regardless, but if numerical tolerance problems...
		ybottom<-seq(ylim[2]-d,ylim[1],by=-d)
		ytop<-seq(ylim[2],ylim[1]+d,by=-d)		
	}
	else{
		ybottom<-ylim[1]
		ytop<-ylim[2]
	}
	if(!add) plot(NULL,xlim=range(c(xleft,xright)),ylim=ylim,axes=FALSE,ylab=ylab,xlab=xlab,bty="n",...)
	
	for(i in 1:nrow(m)){
    rect(  xleft=xleft, xright=xright,  ybottom=rep(ybottom[i],ncol(m)) , ytop=rep(ytop[i],ncol(m)), col=m[i,],border=NA,xpd=NA)   
  }
  #hatch lines to indicate samples
  if(tick){
	xl = (xleft+xright)/2
  	axis(1,at=xl,labels=FALSE,tick=TRUE)
#segments(  xl, rep(-0.1,ncol(m)) , xl, rep(0,ncol(m)), col="black")    #** alt white and black color?
  }
	if(!is.null(clNames)){
		y<-(ybottom+ytop)/2
	  	axis(2,at=y,labels=clNames,xpd=NA,las=2,tick=FALSE,line=axisLine)
	} 	
	#draw box around
	if(box){
		rect(xleft=min(xleft), ybottom=min(ybottom), xright=max(xright), ytop=max(ytop))
	}
}


#' Make color annotation for input to dualHeatmap from output of plotTracking
#'
#' @param trackingOut the output from \code{\link{plotTracking}}
#' @rdname plotTracking
makeAnnoColor<-function(trackingOut){
	lapply(trackingOut$groupToColorLegend,function(x){
		z<-as.character(x[,"Color"])
		names(z)<-as.character(as.numeric(x[,"Original"]))
		return(z)
	})
}

