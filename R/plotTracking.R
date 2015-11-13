#' Functions for Visualizing assigments of samples to clusters across multiple clusterings
#' 
#' Align clusters of the same set of samples and provide a plot of their shared
#' values
#' 
#' clusterTrackingPlot is generally called within plotTracking, but is provided
#' as separate function for convenience, for example if you want to manually change a
#' color in the output of plotTracking and replot. 
#' 
#' @name plotTracking
#' @aliases plotTracking clusterTrackingPlot 
#' @docType methods

#' @param clusters A matrix of with each column corresponding to a clustering and
#' each row a sample. The plotting will plot them in order of this matrix,
#' and their order influences the plot greatly.
#' @param index A predefined order in which the samples will be plotted.
#' Otherwise the order is found internally
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
#' @param colPalette a vector of colors used for the different clusters.
#' @return plotTracking provides (invisibly) the orders and other things that
#' go into making the matrix. Specifically, a list with the following elements
#' \itemize{
#'
#' \item{\code{index}}{ a vector of length equal to ncols(clusters) giving the order of the columns to use to get the original clusters matrix into the order made by plotTracking}
#'
#' \item{\code{colors}}{ matrix of color assignments for each element of original clusters matrix. Matrix is in the same order as original clusters matrix. The matrix colors[,index] is the matrix sent to clusterTrackingPlot}
#'
#' \item{\code{aligned}}{ matrix of integer valued cluster assignments that match the colors. this is useful if you want to have cluster identification numbers that are better aligned than that given in the original clusters. Again, matrix is in same order as original clusters matrix}
#'
#' \item{\code{groupToColorLegend}}{ matrix with three columns named "Original","Aligned", and "Color" giving the correspondence between the original cluster ids in clusters, the aligned cluster ids in aligned cluster values, and the color (respectively) }
#' }
#'
#' @author Elizabeth Purdom and Marla Johnson
#' 
#' @examples
#' #clustering using pam: try using different dimensions of pca and different k
#' ps<-c(5,10,50)
#' names(ps)<-paste("npc=",ps,sep="")
#' pcaData<-stats::prcomp(simData, center=TRUE, scale=TRUE)
#' cl <- compareChoices(lapply(ps,function(p){pcaData$x[,1:p]}), clusterMethod="pam",ks=2:4,findBestK=c(TRUE,FALSE))
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
#' #notice that the blue and orange are really arbitrarily different colors because not next to each other in row.
#' #We can manually change the blue to orange :
#' #first find out color names by looking at the out$colors (but in right order using out$index)
#' out$colors[out$index[1:10],c("npc=5,k=NA,findBestK=TRUE","npc=5,k=3,findBestK=FALSE")]
#' #change "#1F78B4" to "#FF7F00"
#' newColorMat<-out$colors
#' newColorMat[newColorMat=="#1F78B4"]<-"#FF7F00"
#' #replot by calling clusterTrackingPlot directly
#' clusterTrackingPlot(newColorMat[out$index,])
#'
#' #We can also change the order of the rows. Notice how this dramatically changes the plot!
#' plotTracking(cl$clMat[,c(3:6,1:2,7:ncol(cl$clMat))],axisLine=-2)

plotTracking<-function(clusters, index=NULL,reuseColors=FALSE,matchToTop=FALSE,plot=TRUE,unassignedColor="white",missingColor="grey",minRequireColor=0.3,startNewColors=FALSE,colPalette=bigPalette,...){
	clusters<-t(clusters)
	if(any(as.character(clusters)%in%c("-1","-2"))){
		if(any(apply(clusters,1,function(x){any(is.na(x))}))) stop("clusters should not have 'NA' values; non-clustered samples should get a '-1' or '-2' value depending on why they are not clustered.")
		#don't think I need this anymore because fixed so not need numeric values.
		# rnames<-row.names(clusters)
		# clusters<-apply(clusters,2,as.numeric)
		# rownames(clusters)<-rnames
		# if(any(is.na(clusters))) stop("could not convert clusters to numeric, but some samples had value '-1'")
		out<-.plotTrackingInternal(clusters, index=index,reuseColors=reuseColors,matchToTop=matchToTop,plot=FALSE,minRequireColor=minRequireColor,startNewColors=startNewColors,colPalette=colPalette,...)
#		browser()
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
		if(plot) clusterTrackingPlot(xColors[out$index,], colnames(out$colors),...)
		out$colors<-xColors
		out$groupToColorLegend<-newColorLeg
	}
	else{
		out<-.plotTrackingInternal(clusters, index=index,reuseColors=reuseColors,matchToTop=matchToTop,plot=plot,minRequireColor=minRequireColor,startNewColors=startNewColors,colPalette=colPalette,...)
	}	
	invisible(out)
}
.plotTrackingInternal<-function(clusters, index=NULL,reuseColors=FALSE,matchToTop=FALSE,plot=TRUE,minRequireColor=0.3,startNewColors=FALSE,colPalette=.thisPal,...){
	#clusters is a nmethods x nsamples matrix. The order of the rows determines how the tracking will be done.
	#if given, index is the order of the samples (columns) for plotting. Otherwise determined by program
	#matchToTop allows that the comparison for color assignment be done, not row by row, but always back to the first row.
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

	if(plot) clusterTrackingPlot(t(colorM[,index]), rownames(clusters),...)
	dimnames(colorM)<-dimnames(clusters)
	allColors<-unique(as.vector(colorM))
	alignCl<-apply(colorM,2,function(x){match(x,allColors)})
	groupToColorLegend<-lapply(1:nrow(clusters),function(ii){
		mat<-cbind("Original"=unlist(clusters[ii,]),"Aligned"=unlist(alignCl[ii,]),"Color"=unlist(colorM[ii,]))
		rownames(mat)<-NULL
		(unique(mat))
	})
	names(groupToColorLegend)<-rownames(clusters)
	invisible(list(index=index,colors=t(colorM),aligned=t(alignCl),groupToColorLegend=groupToColorLegend))

}

.setClusterColors <- function(past_ct,ct,colorU,past_colors,reuseColors,minRequireColor=0.3,startNewColors=FALSE){
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

#' @rdname plotTracking
#' 
#' @param colorMat Matrix where columns are different clusterings, rows are samples,
#' and entries of matrix are colors (i.e. characters that define colors)
#' @param names names to go with the columns (clusterings) in matrix colorMat
#' @param add whether to add to existing plot
#' @param x values on x-axis at which to plot the rows of m
#' @param ylim vector of limits of y-axis
#' @param tick logical, whether to draw ticks on x-axis for each sample.
#' @param ylab character string for the label of y-axis 
#' @param xlab character string for the label of x-axis 
#' @param axisLine the number of lines in the axis labels on y-axis should be (passed to line=... in the axis call)
#' @param box logical, whether to draw box around the plot
#' @param ... for \code{plotTracking} arguments passed to \code{clusterTracking}. For \code{clusterTracking}, arguments passed to \code{\link{plot}} if \code{add=FALSE}.

clusterTrackingPlot <- function(colorMat, names=rownames(m),add=FALSE,x=NULL,ylim=NULL,tick=FALSE,ylab="",xlab="",axisLine=2,box=FALSE,...){
  	m<-t(colorMat)
  #description: plots cluster tracking plot
  #input: m - matrix where rows are k, columns are samples, and entries are color assignments
  #names: names that correspond with rows of m
	if(is.null(x)){
		xleft<-seq(0,1-1/ncol(m),by=1/ncol(m))
		xright<-seq(1/ncol(m),1,by=1/ncol(m))
	}
	else{
		d<-unique(round(diff(x),6)) #distance between x values; had to round because got differences due to small numerical storage issues
		if(length(d)!=1) stop("invalid x values given -- not evenly spaced")
		xleft<-x-d/2
		xright<-x+d/2
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
	if(!is.null(names)){
		y<-(ybottom+ytop)/2
	  	axis(2,at=y,labels=names,xpd=NA,las=2,tick=FALSE,line=axisLine)
	} 	
	#draw box around
	if(box){
		rect(xleft=min(xleft), ybottom=min(ybottom), xright=max(xright), ytop=max(ytop))
	}
}

###After this point, don't export, but didn't want to lose the functions.

.makeFactorsNumbers<-function(mat,MARGIN,sharedFactors=FALSE){
	if(sharedFactors){
		if(is.data.frame(mat)) mat<-as.matrix(mat)
		allcolors<-unique(as.vector(mat))
		out<-apply(mat,MARGIN,function(y){as.numeric(factor(y,levels=allcolors))})
	}
	else out<-apply(mat,MARGIN,function(y){as.numeric(factor(y))})
	if(MARGIN==1) return(t(out)) #apply inverts it when do apply to rows
	else return(out)
}
.plotWithClinical<-function(clusters,clinical,clinicalNumbers=NULL,index=NULL,whClinical=NULL,main="",split=TRUE,mar=c(2.1,7.1,1.1,1.1),plot=TRUE,...){
	cl<-clusters
		order<-plotTracking(cl, ylab="",index=index,tick=FALSE,reuseColors=FALSE,plot=FALSE)	
		if(is.null(clinicalNumbers)) clinicalNumbers<-.makeFactorsNumbers(clinical,MARGIN=2)
		rn<-row.names(cl)
		cn<-colnames(cl)
		cl<-.makeFactorsNumbers(order$colors,MARGIN=1,sharedFactors=TRUE)
	#	cl<-apply(cl,2,as.character)
		##Put them in the middle
		if(split){
			if(is.null(whClinical)){
				nabove<-floor(ncol(clinicalNumbers)/2)
				x<-rbind(tail(t(clinicalNumbers),-nabove),cl,head(t(clinicalNumbers),nabove))
				rownames(x)<-c(tail( colnames(clinical),-nabove),rn,head(colnames(clinical),nabove))
			}
			else{
				nabove<-floor(ncol(clinicalNumbers[,whClinical])/2)
				x<-rbind(tail(t(clinicalNumbers[,whClinical]),-nabove),cl,head(t(clinicalNumbers[,whClinical]),nabove))
				rownames(x)<-c(tail( whClinical,-nabove),rn,head(whClinical,nabove))
			}
		}
		else{
			if(is.null(whClinical)) {
				x<-rbind(t(clinicalNumbers),cl)
				rownames(x)<-c(colnames(clinical),rn)				
			}
			else{
				x<-rbind(t(clinicalNumbers[,whClinical]),cl)
				rownames(x)<-c(whClinical,rn)				
				
			}
		}
		colnames(x)<-cn

		par(mar=mar)
		xcol<-apply(x,2,function(y){.thisPal[y]})
		dimnames(xcol)<-dimnames(x)
		if(plot) clusterTrackingPlot(xcol[,order$index], names=row.names(x),main=main,...)
#		image(x=1:ncol(x),y=1:nrow(x),t(x)[order$index,],col=.thisPal,xlab="",ylab="",xaxt="n",yaxt="n",main=main)
#		axis(2,at=1:nrow(x),labels=rownames(x),las=2,cex=.5)
		invisible(list(x=x,xcol=xcol,order=order$index))
}

.alignClusters<-function(clAlign,clRef=NULL,index=NULL){
	if(!is.null(clRef)){
		if(is.null(dim(clRef))) clRef<-matrix(clRef,nrow=1)
		else if(nrow(clRef)>1) stop("clRef can only be a single alignment at this time")
		colnames(clRef)<-colnames(clAlign)
		if(is.null(row.names(clRef))) rownames(clRef)<-as.character(1:nrow(clRef))
		nrowRef<-nrow(clRef)
		alignCluster<-plotTracking(rbind(clRef,clAlign),index=index,plot=FALSE)
	} 
	else alignCluster<-plotTracking(clAlign,index=index,plot=FALSE)

	if(!is.null(clRef)){ #assume clRef already correct
		oldAssign<-alignCluster$aligned[1:nrowRef, ,drop=FALSE]
		newAssign<-clRef
	}
	else{
		oldAssign<-clAlign
		newAssign<-alignCluster$aligned
		dimnames(newAssign)<-dimnames(clAlign)
	}
	##Make color mapping
	allClId<-sort(unique(as.vector(unlist(newAssign))))
	oldToNew<-lapply(1:nrow(newAssign),function(ii){unique(cbind("New"=unlist(newAssign[ii,]),"Old"=unlist(oldAssign[ii,])))})
	oldToNew<-sapply(oldToNew,function(z){z[match(allClId,z[,"New"]),"Old"]  })
	colnames(oldToNew)<-rownames(newAssign) 
	rownames(oldToNew)<-NULL
	colorMapping<-data.frame("NewClusterId"=as.numeric(allClId),"Color"=.thisPal[allClId],oldToNew,check.names=FALSE)
	
	if(is.null(clRef)){
		return(list(colorMapping=colorMapping,alignedClusters=newAssign,alignedColors=alignCluster$colors))
	}
	else{
		rm(oldAssign)
		rm(newAssign)
		oldAssign<-alignCluster$align[-c(1:nrowRef), , drop=FALSE] #what want to align to
		singlePlotTrackingId<-apply(colorMapping[,-c(1:2),drop=FALSE],1,function(x){unique(na.omit(x))})
		colorMapping<-data.frame(colorMapping[,c(1:2)],"plotTrackingId"=singlePlotTrackingId)
		oldInRef<-na.omit(unique(as.vector(unlist(colorMapping[,-c(1:2)]))))
		# tempMap<-unique(cbind(tempId=lamlAlignCluster[sprintf("Iso: K=%s",kIso),],newId=isoIdAssignments[sprintf("Iso: K=%s",kIso),]))
		#match gene/iso to prop colors already decided on
		notInRef<-unique(as.vector(unlist(oldAssign)))
		notInRef<-setdiff(notInRef,oldInRef)
		notUsed<-setdiff(1:max(colorMapping$NewClusterId),colorMapping$NewClusterId)
		remaining<-c(notUsed,(max(colorMapping$NewClusterId)+1):length(.thisPal))
		if(length(remaining)< length(notInRef)) stop("Not enough in .thisPal to cover all")
		#conver plotTrackingId to final 'NewClusterId'
		tempMap<-rbind(colorMapping[,c(1,3)],cbind("NewClusterId"=remaining[1:length(notInRef)],plotTrackingId=notInRef))
		newAssign<-apply(oldAssign,2,function(x){tempMap[,"NewClusterId"][match(x,tempMap[,"plotTrackingId"])]})
		dimnames(newAssign)<-dimnames(clAlign)
		#make color/cluster mapping
		allClId<-sort(unique(as.vector(newAssign)))
		oldToNew<-lapply(1:nrow(newAssign),function(ii){unique(cbind("New"=unlist(newAssign[ii,]),"Old"=unlist(oldAssign[ii,])))})
		oldToNew<-sapply(oldToNew,function(z){z[match(allClId,z[,"New"]),"Old"]  })
		colnames(oldToNew)<-rownames(newAssign) 
		rownames(oldToNew)<-NULL
		colorMapping<-data.frame("NewClusterId"=allClId,"Color"=.thisPal[allClId],oldToNew,check.names=FALSE)
		alignedColors<-apply(newAssign,2,function(x){.thisPal[x]})
		rownames(alignedColors)<-row.names(clAlign)
		rownames(newAssign)<-row.names(clAlign)
		return(list(colorMapping=colorMapping,alignedClusters=newAssign,alignedColors=alignedColors))

	}
	
}


###Functions to manipulate order of blocks of samples.
.printWithIndex<-function(x,whBlockRow=1,blockOrder=NULL,whRow=NULL,whBlocks=NULL,checkAllIndices=TRUE){
	if(is.null(whRow)) whRow<-1:nrow(x$aligned)
	allCols<-c(1:ncol(x$aligned))
	newx<-rbind(Index=allCols,x$aligned[whRow,])
	# if(is.null(whCol)) whCol<-1:ncol(x$aligned)
	# newIndex<-match(x$index, allCols[whCol])
	newIndex<-x$index
	orderedX<-newx[,newIndex]
	zoom<-orderedX[whBlockRow+1,]
	runs<-rle(zoom)
	blocks<-rep(1:length(runs$lengths),times=runs$length)
	orderedX<-rbind(Blocks=blocks,orderedX)
	if(!is.null(blockOrder)){
		if(length(blockOrder)>length(runs$lengths) | max(blockOrder) > length(runs$lengths)) stop("invalid block orders -- too long")
		blockStarts<-c(1,head(cumsum(runs$lengths),-1)+1)
		blockEnds<-cumsum(runs$lengths)
		if(length(blockStarts)!=length(blockEnds) | length(blockEnds)!=length(runs$lengths)) stop("coding error")
		orderedX<-do.call("cbind",mapply(blockStarts[blockOrder],blockEnds[blockOrder],FUN=function(s,e){orderedX[,seq(s,e)]},SIMPLIFY=FALSE))
	}
	if(checkAllIndices & ! length(unique(orderedX["Index",]))==length(allCols)) stop("whBlocks doesn't capture all indices")
	return(orderedX)
}
.clusterIndInBlock<-function(x,whInd,whBlockRow=1,whBlocks,rightLeftPerBlock=NULL){
	whBlocks<-sort(whBlocks) #don't want to rearrange order of blocks...
	#get indices that are in block:
	orderedX<-.printWithIndex(x,whBlockRow=whBlockRow,blockOrder=NULL)
	blockX<-.printWithIndex(x,whBlockRow=whBlockRow,blockOrder=whBlocks,checkAllIndices=FALSE)[1:2,]
	isFac<-as.numeric(blockX["Index",] %in% whInd)
	blockX<-rbind(blockX,isFac) #make third row in matrix the indicator that have it.
	index<-orderedX["Index",]
	if(!is.null(rightLeftPerBlock)){
		out<-mapply(whBlocks,rightLeftPerBlock,FUN=function(bl,dir){
			x<-blockX[,blockX["Blocks",]==bl]
			isFacX<-which(x[3,]==1)
			newX<-if(dir=="right") cbind(x[,-isFacX],x[,isFacX]) else cbind(x[,isFacX],x[,-isFacX])		
			newIndex<-newX["Index",]
			m<-match(newIndex,index)
			index[sort(m)]<<-newIndex
		})
	}
	else{
		out<-by(t(blockX),blockX["Blocks",],function(x){
			x<-t(x)
			isFacX<-which(x[3,]==1)
			#now replace old matrix index with this one
			newX<-cbind(x[,isFacX],x[,-isFacX])
			newIndex<-newX["Index",]
			m<-match(newIndex,index)
			index[sort(m)]<<-newIndex
			})		
	}
	orderedX["Index",]<-index
	return(orderedX)
}


