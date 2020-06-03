#' @name plotClustersTable
#' @title Plot heatmap of cross-tabs of 2 clusterings
#' @description Plot heatmap of cross-tabulations of two clusterings
#' @aliases plotClustersTable,ClusterExperiment-method
#' @inheritParams getClusterIndex
#' @inheritParams plotHeatmap
#' @param object ClusterExperiment object (or matrix with table result)
#' @param ignoreUnassigned logical as to whether to ignore unassigned clusters
#'   in the plotting. This means they will also be ignored in the calculations
#'   of the proportions (if \code{margin} not NA).
#' @param margin if NA, the actual counts from \code{tableClusters} will be
#'   plotted. Otherwise, \code{\link[base]{prop.table}} will be called and the
#'   argument \code{margin} will be passed to \code{prop.table} to determine
#'   whether proportions should be calculated. If '1', then the proportions in
#'   the rows sum to 1, if '2' the proportions in the columns sum to 1. If
#'   'NULL' then the proportion across the entire matrix will sum to 1. An
#'   additional option has been added so that if you set  \code{margin=0}, the
#'   entry displayed in each cell will be the proportion equal to the size of
#'   the intersection over the size of the union of the clusters (a Jaccard
#'   similarity between the clusters), in which case each entry is a value
#'   between 0 and 1 but no combination of the entries sum to 1.
#' @param cluster logical, whether to cluster the rows and columns of the table.
#'   Passed to arguments \code{clusterFeatures} AND \code{clusterSamples} of
#'   \code{plotHeatmap}.
#' @param clusterLegend list in \code{clusterLegend} format that gives colors
#'   for the clusters tabulated.
#' @param plotType type of plot. If "heatmap", then a heatmap will be created of
#'   the values of the contingency table of the two clusters (calculated as
#'   determined by the argument "margin") using \code{\link{plotHeatmap}}. If
#'   "bubble", then a plot will be created using \code{bubblePlot}, which will
#'   create circles for each cell of the contingencey table whose size
#'   corresponds to the number of samples shared and the color based on the
#'   value of the proportion (as chosen by the argument \code{margin}).
#' @param main title of plot, passed to \code{plotHeatmap} or to the argument
#'   \code{propLabel} in \code{bubblePlot}
#' @param propLabel the label to go with the legend of the color of the
#'   bubbles/circles
#' @param ... arguments passed on to \code{plotHeatmap} or \code{bubblePlot}
#'   depending on choice of \code{plotType}. Note that these functions take
#'   different arguments so that switching from one to the other may not take
#'   all arguments. In particular \code{bubblePlot} calls \code{plot} while
#'   \code{plotHeatmap} calls \code{\link{NMF}{aheatmap}}.
#' @param propTable table of proportions (\code{bubblePlot}))
#' @param sizeTable table of sizes (only for use in \code{bubblePlot} or
#'   \code{plotType="bubble"}). See details.
#' @param gridColor color for grid lines (\code{bubblePlot}))
#' @param cexFactor factor to multiple by to get values of circles. If missing,
#'   finds value automatically, namely by using the maxCex value default.
#'   Overrides value of maxCex. (\code{bubblePlot}))
#' @param maxCex largest value of cex for any point (others will scale
#'   proportionally smaller) (\code{bubblePlot})).
#' @param ylab label for labeling clustering on the y-axis. If NULL, will
#'   determine names. If set to \code{NA} no label for clustering on the y-axis
#'   will be plotted (to turn off legend of the clusterings in heatmap, set
#'   \code{legend=FALSE}).
#' @param xlab label for labeling clustering on the x-axis. If NULL, will
#'   determine names. If set to \code{NA} no label for clustering on the x-axis
#'   will be plotted (to turn off legend of the clusterings in heatmap, set
#'   \code{legend=FALSE}).
#' @param las the value for the las value in the call to \code{\link{axis}} in
#'   labeling the clusters in the bubble plot. Determines whether parallel or
#'   perpindicular labels to the axis (see \code{\link{par}}).
#' @param legend whether to draw legend along top (bubble plot) or the color
#'   legend (heatmap)
#' @param colorScale the color scale for the values of the proportion table
#' @param useNames for \code{tableClusters}, whether the output should be tabled
#'   with names (\code{useNames=TRUE}) or ids (\code{useNames=FALSE})
#' @param tableMethod the type of table calculation to perform. "intersect"
#'   refers to the standard contingency table (\code{\link[base]{table}}), where
#'   each entry of the resulting table is the number of objects in both
#'   clusters. "union" instead gives for each entry the number of objects that
#'   are in the union of both clusters.
#' @seealso \code{\link{plotHeatmap}}
#' @details For \code{plotClustersTable} applied to the class \code{table},
#'   \code{sizeTable} is passed to \code{bubblePlot} to indicate the size of the
#'   circle. If \code{sizeTable=NULL}, then it is assumed that the \code{object}
#'   argument is the table of counts and both the \code{propTable} and
#'   \code{sizeTable} are set to the same value (hence turning off the coloring
#'   of the circle/bubbles). This is equivalent effect to the \code{margin=NA}
#'   option of \code{plotClustersTable} applied to the \code{ClusterExperiment}
#'   class.
#' @details Note that the cluster labels in \code{plotClustersTable} and
#'   \code{tableClusters} are converted to "proper" R names via
#'   \code{make.names}. This is because \code{tableClusters} calls the R
#'   function \code{table}, which makes this conversion
#' @details For \code{plotClustersTable}, \code{whichClusters} should define 2
#'   clusters, while for \code{tableClusters} it can indicate arbitrary number.
#' @seealso \code{\link[base]{table}}
#' @return \code{tableClusters} returns an object of class \code{table} (see
#'   \code{\link[base]{table}}).
#' @return \code{plotClustersTables} returns invisibly the plotted proportion
#'   table. In particular, this is the result of applying
#'   \code{\link[base]{prop.table}} to the results of \code{tableClusters}
#'   (after removing unclustered samples if \code{ignoreUnassigned=TRUE}).
#' @rdname plotClustersTable
#' @author Kelly Street, Elizabeth Purdom
#' @seealso \code{\link[base]{prop.table}}
#' @examples
#' #clustering using pam: try using different dimensions of pca and different k
#' data(simData)
#'
#' cl <- clusterMany(simData, nReducedDims=c(5, 10, 50), reducedDim="PCA",
#' clusterFunction="pam", ks=2:4, findBestK=c(TRUE,FALSE),
#' removeSil=c(TRUE,FALSE), makeMissingDiss=TRUE)
#' #give arbitrary names to clusters for demonstration
#' cl<-renameClusters(cl,value=letters[1:nClusters(cl)[1]],whichCluster=1)
#' tableClusters(cl,whichClusters=1:2)
#' #show options of margin in heatmap format:
#' par(mfrow=c(2,3))
#' plotClustersTable(cl,whichClusters=1:2, margin=NA, legend=FALSE,
#'   ignoreUnassigned=TRUE)
#' plotClustersTable(cl,whichClusters=1:2, margin=0, legend=FALSE,
#'   ignoreUnassigned=TRUE)
#' plotClustersTable(cl,whichClusters=1:2, margin=1, legend=FALSE,
#'   ignoreUnassigned=TRUE)
#' plotClustersTable(cl,whichClusters=1:2, margin=2, legend=FALSE,
#'   ignoreUnassigned=TRUE)
#' plotClustersTable(cl,whichClusters=1:2, margin=NULL, legend=FALSE,
#'   ignoreUnassigned=TRUE)
#'
#' #show options of margin in bubble format:
#' par(mfrow=c(2,3))
#' plotClustersTable(cl,whichClusters=1:2, margin=NA, 
#'    ignoreUnassigned=TRUE, plotType="bubble")
#' plotClustersTable(cl,whichClusters=1:2, margin=0,
#'    ignoreUnassigned=TRUE, plotType="bubble")
#' plotClustersTable(cl,whichClusters=1:2, margin=1,
#'    ignoreUnassigned=TRUE, plotType="bubble")
#' plotClustersTable(cl,whichClusters=1:2, margin=2,
#'    ignoreUnassigned=TRUE, plotType="bubble")
#' plotClustersTable(cl,whichClusters=1:2, margin=NULL,
#'    ignoreUnassigned=TRUE, plotType="bubble")
#' @export
setMethod( 
  f = "plotClustersTable",
  signature = signature(object = "ClusterExperiment"),
  definition = function(object, whichClusters,ignoreUnassigned=FALSE,margin=NA,...){
    whCl<-getClusterIndex(object,whichClusters=whichClusters,noMatch="throwError")
    if(length(whCl)!=2) stop("invalid choice of 'whichClusters' -- must be exactly 2 clusterings chosen.")
	tableAll<-tableClusters(object,whichClusters=whCl,useNames=TRUE,tableMethod="intersect")
	#sort them by their names
	#if ever gets slow, could do this directly so not call tableClusters twice...
	isMargin0<-!is.null(margin) && !is.na(margin) && margin==0
	if(isMargin0) denomTab<-.makeUnion(tableAll)
			
	cL<-clusterLegend(object)[whCl]
	if(ignoreUnassigned){
		rNms<-rownames(tableAll)
		cNms<-colnames(tableAll)
		mat1<-cL[[1]]
		mat2<-cL[[2]]
		#convert names to clusterIds to check for missing
		rNms<-mat1[match(rNms,mat1[,"name"]),"clusterIds"]
		cNms<-mat2[match(cNms,mat2[,"name"]),"clusterIds"]
		
		if("-1" %in% rNms || "-2" %in% rNms){
			wh<-which(! rNms %in% c("-1","-2"))
			if(length(wh)>0){
				tableAll<-tableAll[wh, ,drop=FALSE]
				if(isMargin0) denomTab<-denomTab[wh, ,drop=FALSE]
				#deal with fact that plotHeatmap doesn't fix problem with rowData!
				mat1<-mat1[mat1[,"clusterIds"] %in% rNms[wh], ,drop=FALSE]
				cL[[1]]<-mat1
			}
			else stop("All of the first clustering are unassigned, cannot use ignoreUnassigned=TRUE")
			
		}
		if("-1" %in% cNms || "-2" %in% cNms){
			wh<-which(! cNms %in% c("-1","-2"))
			if(length(wh)>0){
				tableAll<-tableAll[, wh,drop=FALSE]
				if(isMargin0) denomTab<-denomTab[, wh ,drop=FALSE]
			}
			else stop("All of the second clustering are unassigned, cannot use ignoreUnassigned=TRUE")
		}
		
	}
	
	#must be after the ignoreUnassigned so as to get rid of those before calculate
	if(is.null(margin) || !is.na(margin)){
		if(isMargin0){
			plotTable<-tableAll/denomTab
			propLabel="Jaccard Similarity"
		}
		else{
			plotTable<-prop.table(tableAll,margin=margin)
			if(is.null(margin)) propLabel="% overlap (out of total cells)"
			else if(margin==1) propLabel="% overlap (out of row total)"
			else if(margin==2) propLabel="% overlap (out of col total)"
		}
	}
	else{
		plotTable<-tableAll
		propLabel<-"Count of Overlap"
	}
	plotClustersTable(plotTable,clusterLegend=cL,sizeTable=tableAll, main=propLabel,...)
}
)

		
#' @rdname plotClustersTable
#' @export
setMethod( 
  f = "plotClustersTable",
  signature = signature(object = "table"),
  definition = function(object,plotType=c("heatmap","bubble"), 
  		main="",xlab=NULL,ylab=NULL,legend=TRUE,cluster=FALSE,clusterLegend=NULL, sizeTable=NULL, ...){
		plotType<-match.arg(plotType)
		tableAll<-object
		if(any(dim(tableAll)==1)&& plotType=="heatmap") stop("Cannot create heatmap when there is only 1 column or row in the table") #creates error in aheatmap....
			
		
		varNames<-make.names(names(dimnames(tableAll)))
 	 	if(!cluster){
			#determine for each column, which is the cluster in the row that is largest
			#will reorder the columns based on this result.
			rankValues<-rank(sapply(seq_len(ncol(tableAll)),FUN=function(ii){
				whMax<-which.max(tableAll[,ii])
				#return(whMax)
				 if(length(whMax)==0) return(1) #could happen in NaN in entry
				 else return(whMax)
			}),ties.method="first")
			order2<-order(rankValues)
			
		}
		else order2<-seq_len(ncol(tableAll))
		
		#rows=y-axis
		defaultY<-is.null(ylab) || ( is.na(ylab) & plotType=="heatmap")
		rData<-data.frame(rownames(tableAll))
		names(rData)<-if(defaultY) varNames[1] else ylab
		
		#cols=x-axis
		defaultX<-is.null(xlab) || ( is.na(xlab) & plotType=="heatmap")
		cData<-data.frame(colnames(tableAll)[order2])
		names(cData)<-if(defaultX) varNames[2] else xlab
		labelCols<-is.null(xlab) || !is.na(xlab)
		
		if(!is.null(clusterLegend)) names(clusterLegend)<-make.names(names(clusterLegend))
		if(plotType=="heatmap"){
			passedArgs<-list(...)
			if(!"colorScale" %in% names(passedArgs)){
				passedArgs$colorScale<-colorRampPalette(c("white","black"))(12)
			}
			passedArgs<-c(list(data=tableAll[,order2,drop=FALSE],colData=cData, annRow=rData,
					clusterLegend=clusterLegend, main=main,annLegend=legend,
					annotation_names_col =labelCols,
					clusterFeatures=cluster,clusterSamples=cluster)
					,passedArgs)
			 	 	
		
			do.call("plotHeatmap",passedArgs)
		}
		if(plotType=="bubble"){
			if(!is.null(sizeTable))
				bubblePlot(sizeTable=sizeTable[,order2,drop=FALSE],propTable=tableAll[,order2,drop=FALSE],propLabel=main,xlab=xlab,ylab=ylab,legend=legend,...)
			else
				bubblePlot(sizeTable=tableAll[,order2,drop=FALSE], propTable=tableAll[,order2,drop=FALSE],propLabel=main,xlab=xlab,ylab=ylab,legend=legend,...)
			
		} 		
		invisible(tableAll[,order2,drop=FALSE])
  	
  }
	)

#' @rdname plotClustersTable
#' @aliases tableClusters
#' @export
setMethod( 
  f = "tableClusters",
  signature = signature(object = "ClusterExperiment"),
  definition = function(object, whichClusters="primary", useNames=TRUE, tableMethod=c("intersect","union"),...)
  { 
    whichClusters<-getClusterIndex(object,whichClusters=whichClusters,noMatch="throwError")
    tableMethod<-match.arg(tableMethod)
    if(useNames) numCluster<-clusterMatrixNamed(object,whichClusters=whichClusters)
    else numCluster<-clusterMatrix(object)[,whichClusters]
		tabAll<-table(data.frame(numCluster))
    if(tableMethod=="intersect" | length(whichClusters)==1) return(tabAll)
		else return(.makeUnion(tabAll))
})

.makeUnion<-function(tabAll){
	unionTab<-outer(rowSums(tabAll),colSums(tabAll), FUN = "+") - tabAll
	rownames(unionTab)<-rownames(tabAll)
	colnames(unionTab)<-colnames(tabAll)
	names(attributes(unionTab)$dimnames)<-names(attributes(tabAll)$dimnames)
	class(unionTab)<-"table"
	return(unionTab)
}


#' @rdname plotClustersTable

#' @details \code{bubblePlot} is mainly used internally by
#'   \code{plotClustersTable} but is made public for users who want more control
#'   and to allow documentation of the arguments. \code{bubblePlot} plots a
#'   circle for each intersection of two clusters, where the color of the circle
#'   is based on the value in \code{propTable} and the size of the circle is
#'   based on the value in \code{sizeTable}. If \code{propTable}
#'   is equal to \code{sizeTable}, then the \code{propTable} is ignored and the 
#'   coloring of the circles is not performed, only the adjusting of the
#'   size of the circles based on the total size. The size is determined by 
#'   setting the \code{cex} value of the point as
#'   $sqrt(sizeTable[i,j])/sqrt(max(sizeTable))*cexFactor$. 
#' @importFrom grDevices rgb
setMethod( 
  f = "bubblePlot",
  signature = signature(sizeTable="table",propTable = "table"),
	definition=function(propTable,sizeTable,gridColor=rgb(0,0,0,.05),
	maxCex=8,cexFactor, ylab,xlab,propLabel="Value of %",legend=TRUE,
	las=2, colorScale=RColorBrewer::brewer.pal(11,'Spectral')[-6]){

  if(!all(dim(propTable)==dim(sizeTable))) 
	  stop("propTable and sizeTable must be of the same dimensions")
  if(!all(unlist(dimnames(propTable))==unlist(dimnames(sizeTable)))) 
	  stop("propTable and sizeTable must have the same dimnames")
  if(all(na.omit(propTable)==na.omit(sizeTable))) 
	  doProp<-FALSE
  else doProp<-TRUE
  nc.row <- nrow(sizeTable)
  nc.col <- ncol(sizeTable)
  sizeTable<-sizeTable[nc.row:1, ,drop=FALSE]
  propTable<-propTable[nc.row:1, ,drop=FALSE]
  # set up plotting window
  xlim<-c(1,nc.col)
  xlim<-xlim+.1*diff(xlim)*c(-1,1) #increase size 10% around
  ylim<-c(1,nc.row)
  ylim<-ylim+.1*diff(ylim)*c(-1,1) #increase size 10% around
  graphics::plot(0,0,type="n",xlim = xlim, ylim = ylim, ylab="",xlab="",axes=FALSE,frame.plot=FALSE)
   
  # get x-y coords for a grid
  xx <- rep(seq_len(nc.col), each = nc.row)
  yy <- rep(seq_len(nc.row), times = nc.col)   
    
  # set color based on % overlap
  expect.overlap <- min(c(nc.row,nc.col)) / max(c(nc.row,nc.col))
  legend.vals <- pretty(c(0,1),n=50)
  if(doProp){
	  allColors<-.colorby(c(propTable, 0, expect.overlap,legend.vals), colors=colorScale)
		color <- head(allColors,nc.col*nc.row)
	  legend.col<-tail(allColors,length(legend.vals))
	}
	else color<-"grey"
  # put plotting information into data.frame, so we can sort by size (want
  # smaller points plotted over larger points)
  df <- data.frame(xx,yy, color,
      sizeTable = as.numeric(sizeTable), 
      propTable = as.numeric(propTable))[order(as.numeric(sizeTable), decreasing = TRUE), ,drop=FALSE]
  df <- df[df$sizeTable > 0, ,drop=FALSE]
    
  # grid
  abline(v = seq_len(nc.col),col=gridColor)
  abline(h = seq_len(nc.row),col=gridColor)
    
  # points parameters
  cex.pch<-sqrt(df$sizeTable)/sqrt(max(df$sizeTable))
  if(missing(cexFactor)){
	 cexFactor<-maxCex/max(cex.pch)
  }
  cex.pch<-cex.pch*cexFactor
  graphics::points(df$xx,df$yy, cex=cex.pch, col=as.character(df$color), pch=16)
    
  # labels for plots
  axis(1,at=seq_len(nc.col),colnames(sizeTable), adj=1,tick=FALSE,las=las)
  axis(2,at=seq_len(nc.row),rownames(sizeTable), adj=1,tick=FALSE,las=las)
  if(missing(ylab) || is.null(ylab)) ylab<-names(attributes(sizeTable)$dimnames)[1]
  if(missing(xlab)  || is.null(xlab)) xlab<-names(attributes(sizeTable)$dimnames)[2]
  if(!is.na(xlab)) title(xlab=xlab)
  if(!is.na(ylab)) title(ylab=ylab)
		    
  # % overlap legend
  if(legend){
	nboxes<-length(legend.vals)
	usr<-par("usr")
	y<-rep(usr[4]+diff(usr[3:4]*.05), length=nboxes)
	height<-diff(usr[3:4])*.01
	yspace<-diff(usr[3:4])*.01
	if(doProp){
		#legend for color scale
		x<-seq(usr[1],usr[1]+diff(usr[1:2])*.3,length = nboxes)
		width<-(x[2]-x[1])/2
		xspace<-diff(usr[1:2])*.01

		xlabPos<-seq(1,nboxes,length=6)    
		graphics::rect(xright=x-width,xleft=x+width,
			ybottom=y-height,ytop=y+height,
			border=NA,col=legend.col,xpd=NA)
	    graphics::text(x=mean(x), unique(y)+height+yspace, propLabel,xpd=NA,pos=3)
	    graphics::text(x[xlabPos],unique(y)-height-yspace, labels=legend.vals[xlabPos], xpd=NA,pos=1)
	}

	# bubble size legend
	legSizeVals<-pretty(range(as.numeric(sizeTable)),n=5)[-1] #smallest is never needed
	nvals<-length(legSizeVals) #because can't precisely control 'pretty'
	xc<-seq(mean(usr[1:2]),mean(usr[1:2])+diff(usr[1:2])*.5,length = nvals)
	yc<-rep(unique(y), length=nvals)
	graphics::points(xc, yc,
    	cex = sqrt(legSizeVals)/sqrt(max(df$sizeTable))*cexFactor, 
		col=rgb(0,0,0,.4), pch=16,xpd=NA)
    graphics::text(x=xc, y=yc-height-yspace, labels=legSizeVals, xpd=NA,pos=1)
    graphics::text(x=mean(xc), unique(yc)+height+yspace, "# Cells",xpd=NA,pos=3)
  }

})
		
#' @importFrom scales alpha
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
.colorby <- function(x, alpha = 1, colors){
    mypal <- grDevices::colorRampPalette(colors)
    if(class(x) %in% c('character','logical')){
        x <- as.factor(x)
    }
    if(class(x) %in% c('integer','numeric')){
        x <- cut(x, breaks = 100)
    }
    if(is(x,'factor')){
        return(scales::alpha(mypal(length(levels(x)))[x],alpha=alpha))
    }
}
		
