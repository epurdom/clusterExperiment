#' @name plotClustersTable
#' @title Plot heatmap of cross-tabs of 2 clusterings
#' @description Plot heatmap of cross-tabulations of two clusterings
#' @aliases plotClustersTable,ClusterExperiment-method
#' @param object ClusterExperiment object (or matrix with table result)
#' @param ignoreUnassigned logical as to whether to ignore unassigned clusters
#'   in the plotting. This means they will also be ignored in the calculations
#'   of the proportions (if \code{margin} not NA).
#' @param margin if NA, the counts from \code{tableClusters} will be plotted.
#'   Otherwise, \code{\link[base]{prop.table}} will be called and the argument
#'   \code{margin} will be passed to \code{prop.table} to determine whether
#'   proportions should be calculated. If '1', then the proportions in the rows
#'   sum to 1, if '2' the proportions in the columns sum to 1. If 'NULL' then
#'   the proportion across the entire matrix will sum to 1. For a symmetric
#'   version, you can set \code{margin=0}, and the entry displayed in each cell
#'   will be the proportion equal to the size of the intersection over the size
#'   of the union of the clusters (a Jaccard similarity between the clusters),
#'   in which case each entry is a proportion but no combination of the entries 
#'   sum to 1.
#' @param whichClusters which clusters to tabulate. For \code{plotClustersTable}
#'   should be 2 clusters, for \code{tableClusters} can indicate arbitrary
#'   number.
#' @rdname plotClustersTable
#' @seealso \code{\link[base]{prop.table}}
#' @examples
#' #clustering using pam: try using different dimensions of pca and different k
#' data(simData)
#'
#' cl <- clusterMany(simData, nReducedDims=c(5, 10, 50), reducedDim="PCA",
#' clusterFunction="pam", ks=2:4, findBestK=c(TRUE,FALSE),
#' removeSil=c(TRUE,FALSE))
#' #give arbitrary names to clusters for demonstration
#' cl<-renameClusters(cl,value=letters[1:nClusters(cl)[1]],whichCluster=1)
#' tableClusters(cl,whichClusters=1:2)
#' #heatmap of the counts in each entry of table:
#' plotClustersTable(cl,whichClusters=1:2, ignoreUnassigned=TRUE)
#' @export
setMethod( 
  f = "plotClustersTable",
  signature = signature(object = "ClusterExperiment"),
  definition = function(object, whichClusters,ignoreUnassigned=FALSE,margin=NA,...){
    whCl<-.TypeIntoIndices(object,whClusters=whichClusters)
    if(length(whCl)!=2) stop("invalid choice of 'whichClusters' -- must be exactly 2 clusterings chosen.")
		tableAll<-tableClusters(object,whichClusters=whCl,useNames=TRUE,tableMethod="intersect")
		#if ever gets slow, could do this directly so not call tableClusters twice...
		if(!is.na(margin)& margin==0) denomTab<-.makeUnion(tableAll)
				
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
					if(!is.na(margin) & margin==0) denomTab<-denomTab[wh, ,drop=FALSE]
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
					if(!is.na(margin)& margin==0) denomTab<-denomTab[, wh ,drop=FALSE]
				}
				else stop("All of the second clustering are unassigned, cannot use ignoreUnassigned=TRUE")
			}
			
		}
		
		#must be after the ignoreUnassigned so as to get rid of those before calculate
		if(!is.na(margin)){
			if(margin==0){
				plotTable<-tableAll/denomTab
			}
			else{
				plotTable<-prop.table(tableAll,margin=margin)
			}
		}
		else plotTable<-tableAll
		plotClustersTable(plotTable,clusterLegend=cL,sizeTable=tableAll,...)
}
)

		
#' @rdname plotClustersTable
#' @param cluster logical, whether to cluster the rows and columns of the table. Passed
#'  to arguments \code{clusterFeatures} AND \code{clusterSamples} of \code{plotHeatmap}.
#' @param clusterLegend list in \code{clusterLegend} format that gives colors for the
#'  clusters tabulated.
#' @param ... arguments passed on to \code{plotHeatmap} or \code{bubblePlot} 
#'  depending on choice of \code{plotType}
#' @seealso \code{\link{plotHeatmap}}
#' @details Note that the cluster labels in \code{plotClustersTable} and 
#' \code{tableClusters} are converted to "proper" R names via \code{make.names}. This is 
#' because \code{tableClusters} calls the R function \code{table}, which makes this 
#' conversion
#' @inheritParams plotHeatmap
#' conversion.
#' @seealso \code{\link[base]{table}}
#' @export
setMethod( 
  f = "plotClustersTable",
  signature = signature(object = "table"),
  definition = function(object,clusterLegend=NULL,cluster=FALSE,plotType=c("heatmap","bubble"), sizeTable=object, ...){
		plotType<-match.arg(plotType)
		tableAll<-object
		varNames<-make.names(names(dimnames(tableAll)))
 	 	
		if(!cluster){
			rankValues<-rank(sapply(1:ncol(tableAll),FUN=function(ii){
				whMax<-which.max(tableAll[,ii])
	
			}),ties.method="first")
			order2<-order(rankValues)
			
		}
		else order2<-1:ncol(tableAll)
			
		rData<-data.frame(rownames(tableAll))
		names(rData)<-varNames[1]

 	 	cData<-data.frame(colnames(tableAll)[order2])
		names(cData)<-varNames[2]
		if(!is.null(clusterLegend)) names(clusterLegend)<-make.names(names(clusterLegend))
		if(plotType=="heatmap"){
			passedArgs<-list(...)
			if(!"colorScale" %in% names(passedArgs)){
				passedArgs$colorScale<-colorRampPalette(c("white","black"))(12)
			}
			passedArgs<-c(list(data=tableAll[,order2],colData=cData, annRow=rData,
					clusterLegend=clusterLegend,
					clusterFeatures=cluster,clusterSamples=cluster)
					,passedArgs)
			do.call("plotHeatmap",passedArgs)
		}
			if(plotType=="bubble") bubblePlot(propTable=tableAll[,order2],sizeTable=sizeTable[,order2],...)
		invisible(tableAll[,order2])
  	
  }
	)

#' @aliases tableClusters
#' @rdname plotClustersTable
#' @export
setMethod( 
  f = "tableClusters",
  signature = signature(object= "ClusterExperiment",whichClusters="character"),
  definition = function(object, whichClusters,...)
  {
    wh<-.TypeIntoIndices(object,whClusters=whichClusters)
    if(length(wh)==0) stop("invalid choice of 'whichClusters'")
    return(tableClusters(object,whichClusters=wh,...))
    
  })

#' @rdname plotClustersTable
#' @export
setMethod( 
  f = "tableClusters",
  signature = signature(object = "ClusterExperiment",whichClusters="missing"),
  definition = function(object, whichClusters,...)
  {
    tableClusters(object,whichClusters="primaryCluster")
    
  })

#' @rdname plotClustersTable
#' @param useNames for \code{tableClusters}, whether the output should be tabled
#'   with names (\code{useNames=TRUE}) or ids (\code{useNames=FALSE})
#' @param tableMethod the type of table calculation to perform. "intersect" refers to the standard contingency table (\code{\link[base]{table}}), where each entry of the resulting table is the number of objects in both clusters. "union" instead gives for each entry the number of objects that are in the union of both clusters.
#' @export
setMethod( 
  f = "tableClusters",
  signature = signature(object = "ClusterExperiment",whichClusters="numeric"),
  definition = function(object, whichClusters, useNames=TRUE, tableMethod=c("intersect","union"),...)
  { 
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
#' @param propTable table of proportions
#' @param sizeTable table of sizes
#' @param gridColor color for grid lines
#' @param cexFactor factor to multiple by to get values of circles. If missing, finds value automatically, namely by using the maxCex value default. Overrides value of maxCex.
#' @param maxCex largest value of cex for any point (others will scale proportionally smaller). 
#' @param ylab label for y-axis. If missing, uses the name for rows in sizeTable
#' @param xlab label for x-axis. If missing, uses the name for columns in sizeTable
#' @param legend whether to draw legend along top
#' @param colorScale the color scale for the values of the proportion table
#' @details \code{bubblePlot} is mainly used internally by \code{plotClustersTable} but is made public for users who want more control and to allow documentation of the arguments. \code{bubblePlot} plots a circle for each intersection of two clusters, where the color of the circle is based on the value in \code{propTable} and the size of the circle is based on the value in \code{sizeTable}. The size is determined by setting the \code{cex} value of the point as $sqrt(sizeTable[i,j])/sqrt(max(sizeTable))*cexFactor$. 
setMethod( 
  f = "bubblePlot",
  signature = signature(propTable = "table",sizeTable="table"),
	definition=function(propTable,sizeTable,gridColor=rgb(0,0,0,.05),maxCex=8,cexFactor,
		ylab,xlab,legend=TRUE,las=2, colorScale=RColorBrewer::brewer.pal(11,'Spectral')[-6]){
#browser()
	if(!all(dim(propTable)==dim(sizeTable))) stop("propTable and sizeTable must be of the same dimensions")
		if(!all(unlist(dimnames(propTable))==unlist(dimnames(sizeTable)))) stop("propTable and sizeTable must have the same dimnames")
	 nc.row <- nrow(propTable)
	 nc.col <- ncol(propTable)
	propTable<-propTable[nc.row:1,]
	sizeTable<-sizeTable[nc.row:1,]
  # set up plotting window
	xlim<-c(1,nc.col)
	xlim<-xlim+.1*diff(xlim)*c(-1,1) #increase size 10% around
	ylim<-c(1,nc.row)
	ylim<-ylim+.1*diff(ylim)*c(-1,1) #increase size 10% around
  plot(0,0,type="n",xlim = xlim, ylim = ylim, ylab="",xlab="",axes=FALSE,frame.plot=FALSE)
   
  # get x-y coords for a grid
  xx <- rep(1:nc.col, each = nc.row)
  yy <- rep(1:nc.row, times = nc.col)   
    
  # set color based on % overlap
  expect.overlap <- min(c(nc.row,nc.col)) / max(c(nc.row,nc.col))
  legend.vals <- pretty(c(0,1),n=50)
	allColors<-.colorby(c(propTable,0,expect.overlap,legend.vals),colors=colorScale)
  color <- head(allColors,nc.col*nc.row)
  legend.col<-tail(allColors,length(legend.vals))
  # put plotting information into data.frame, so we can sort by size (want
  # smaller points plotted over larger points)
  df <- data.frame(xx,yy, color,
      sizeTable = as.numeric(sizeTable), 
      propTable = as.numeric(propTable))[order(as.numeric(sizeTable), decreasing = TRUE),]
  df <- df[df$sizeTable > 0,]
    
  # grid
  abline(v = 1:nc.col,col=gridColor)
  abline(h = 1:nc.row,col=gridColor)
  # rect(-1,nc.row+1,nc.col+1,nc.row+20, col='white', lty=0) # top
  # rect(-1,-20,nc.col+1,0, col='white', lty=0) # bottom
  # rect(-20,-1,0,nc.row+1, col='white', lty=0) # left
  # rect(nc.col+1,-1,nc.col+20,nc.row+1, col='white', lty=0) # right
    
  # plot points
	cex.pch<-sqrt(df$sizeTable)/sqrt(max(df$sizeTable))
	if(missing(cexFactor)){
		cexFactor<-maxCex/max(cex.pch)
	}
	cex.pch<-cex.pch*cexFactor
  points(df$xx,df$yy, cex=cex.pch, col=as.character(df$color), pch=16)
    
  # labels for plots
  axis(1,at=1:nc.col,colnames(sizeTable), adj=1,tick=FALSE,las=las)
  axis(2,at=1:nc.row,rownames(sizeTable), adj=1,tick=FALSE,las=las)
  if(missing(ylab)) ylab<-names(attributes(sizeTable)$dimnames)[1]
  if(missing(xlab)) xlab<-names(attributes(sizeTable)$dimnames)[2]
	if(!is.null(xlab)) title(xlab=xlab)
	if(!is.null(ylab)) title(ylab=ylab)
		
    
  # % overlap legend
	if(legend){
		#legend for color scale
		nboxes<-length(legend.vals)
		usr<-par("usr")
		x<-seq(usr[1],usr[1]+diff(usr[1:2])*.3,length = nboxes)
		width<-(x[2]-x[1])/2
		y<-rep(usr[4]+diff(usr[3:4]*.05), length=nboxes)
		height<-diff(usr[3:4])*.01
		xspace<-diff(usr[1:2])*.01
		yspace<-diff(usr[3:4])*.01

		xlabPos<-seq(1,nboxes,length=6)    
		rect(xright=x-width,xleft=x+width,
			ybottom=y-height,ytop=y+height,
			border=NA,col=legend.col,xpd=NA)
    text(x=mean(x), unique(y)+height+yspace, "Value of %",xpd=NA,pos=3)
    text(x[xlabPos],unique(y)-height-yspace, labels=legend.vals[xlabPos], xpd=NA,pos=1)

    # bubble size legend
		legSizeVals<-pretty(range(as.numeric(sizeTable)),n=5)[-1] #smallest is never needed
		nvals<-length(legSizeVals) #because can't precisely control 'pretty'
		xc<-seq(mean(usr[1:2]),mean(usr[1:2])+diff(usr[1:2])*.5,length = nvals)
		yc<-rep(unique(y), length=nvals)
		points(xc, yc,
        cex = sqrt(legSizeVals)/sqrt(max(df$sizeTable))*cexFactor, 
				col=rgb(0,0,0,.4), pch=16,xpd=NA)
    text(x=xc, y=yc-height-yspace, labels=legSizeVals, xpd=NA,pos=1)
    text(x=mean(xc), unique(yc)+height+yspace, "# Cells",xpd=NA,pos=3)
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
    if(class(x)=='factor'){
        return(scales::alpha(mypal(length(levels(x)))[x],alpha=alpha))
    }
}
		
