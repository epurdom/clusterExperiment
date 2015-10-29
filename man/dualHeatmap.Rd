\name{dualHeatmap}
\docType{methods}
\alias{dualHeatmap}
\alias{seqPal1}
\alias{seqPal2}
\alias{seqPal3}
\alias{seqPal4}
\alias{seqPal5}
\alias{showPalettes}
\title{Heatmap with two different sources for cluster and color scale.
}

\description{
Make heatmap with color scale from one matrix and hiearchical clustering from another. Also color palettes to go with heatmap
}

\usage{
dualHeatmap(clustVec,fq,correctedY,eps=1,dual=TRUE,clusterCells=TRUE, 
	clusterGenes=TRUE,whGenes=1:500,geneNames=FALSE,cellNames=FALSE,
	colorScale=seqPal5, annCol=NULL,annColors=NULL,
	alignColors=FALSE,breaks=NA,...)
showPalettes()
seqPal2
seqPal3
seqPal4
seqPal5
}

\arguments{
\item{clustVec}{
A vector with cluster assignments to show at top of heatmap with cells. ``-1'' indicates the sample was not assigned to a cluster and gets color `white'.
}
\item{fq}{matrix to define the color scale (i.e. the counts). Will be converted to log-scale internally. Assumes genes on rows and cells on columns.
}
\item{correctedY}{Either a matrix define the hiearchical clustering of samples (i.e. normalized data) or a dendrogram for clustering the samples (only used if dual=TRUE).}
\item{eps}{Amount to add to fq so can take log}
\item{dual}{Logical as to whether should use correctedY for dendrogram; otherwise fq is used.}
\item{clusterCells}{Logical as to whether to do hierarchical clustering of cells.}
\item{clusterGenes}{Logical as to whether to do hiearchical clustering of genes.}
\item{whGenes}{Which genes of fq matrix to be used. Default assumes top 500, which may be nonsensical if matrix not ordered.}
\item{geneNames}{Logical as to whether show gene names}
\item{cellNames}{Logical as to whether show cell names}
\item{colorScale}{palette of colors for the color scale of heatmap}
\item{annCol}{data.frame of clusters to show at the top of heatmap (columns are clusters). If NULL, will show clustVec}
\item{annColors}{Assignment of colors to the clusters. If NULL, clusters will be assigned colors. If `annCol' should be list of length equal to ncol(annCol) with names equal to the colnames of annCol; each element of the list should be a vector of colors with names corresponding to the levels of the column of annCol. If annCol=NULL (i.e. use clusterVec), then the name of the length-1 list should be `Cluster'}
\item{alignColors}{Logical as to whether should align the clusters when colors are assigned}
\item{breaks}{Either a vector of breaks (should be equal to length 52), or a number between 0 and 1, indicating that the breaks should be equally spaced (on the log scale+eps) upto the `breaks' quantile. }
\item{...}{passed to aheatmap}
}

\details{
seqPal1-seqPal4 are palettes for the heatmap. showPalettes() will show you these palettes. 

Note that dualHeatmap calles aheatmap under the hood. This allows you to plot multiple heatmaps via par(mfrow=c(2,2)), etc. However, the dendrograms do not resize if you change the size of your plot window in an interactive session of R (this might be a problem for RStudio if you want to pop it out into a large window...). 
}
\value{
	Returns invisibly the breaks, annCol, and annColors that are sent to aheatmap by the function. 
}
\author{
Elizabeth Purdom
}

\examples{
library(gbcData)
data(fq)
data(correctedY)
data(genesValid)

#show the palette colors
showPalettes()

#simple, minimal, example
cl<-rep(1:10,length=ncol(fq))
dualHeatmap(cl,fq=fq,correctedY=correctedY,whGenes=genesValid$Gene)

#assign cluster colors
colors<-bigPalette[20:30]
names(colors)<-1:10
dualHeatmap(cl,fq=fq,correctedY=correctedY,whGenes=genesValid$Gene,annColors=list(colors))

cl2<-sample(cl)
anno<-data.frame(cluster1=cl,cluster2=cl2)
#return the values to see format for giving colors to the annotations
out<-dualHeatmap(cl,fq=fq,correctedY=correctedY,whGenes=genesValid$Gene,annCol=anno)
out$annColors

#compare the 4 other palettes provided
par(mfrow=c(2,2))
dualHeatmap(cl,fq=fq,correctedY=correctedY,whGenes=genesValid$Gene,colorScale=seqPal2,main="seqPal2")
dualHeatmap(cl,fq=fq,correctedY=correctedY,whGenes=genesValid$Gene,colorScale=seqPal3,main="seqPal3")
dualHeatmap(cl,fq=fq,correctedY=correctedY,whGenes=genesValid$Gene,colorScale=seqPal4,main="seqPal4")
dualHeatmap(cl,fq=fq,correctedY=correctedY,whGenes=genesValid$Gene,colorScale=seqPal5,main="seqPal5")

#compare changing breaks quantile
par(mfrow=c(2,2))
dualHeatmap(cl,fq=fq,correctedY=correctedY,whGenes=genesValid$Gene,colorScale=seqPal1,breaks=1,main="Full length")
dualHeatmap(cl,fq=fq,correctedY=correctedY,whGenes=genesValid$Gene,colorScale=seqPal1,breaks=.99,main="0.99 Quantile Upper Limit")
dualHeatmap(cl,fq=fq,correctedY=correctedY,whGenes=genesValid$Gene,colorScale=seqPal1,breaks=.95,main="0.95 Quantile Upper Limit")
dualHeatmap(cl,fq=fq,correctedY=correctedY,whGenes=genesValid$Gene,colorScale=seqPal1,breaks=.90,main="0.90 Quantile Upper Limit")

}
