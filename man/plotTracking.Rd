\name{plotTracking}
\docType{methods}
\alias{plotTracking}
\alias{showThisPal}
\alias{clusterTrackingPlot}
\alias{alignClusters}
\alias{plotWithClinical}
\alias{printWithIndex}
\alias{clusterIndInBlock}
\alias{bigPalette}

\title{
Functions for Tracking Samples across Clusterings
}

\description{
Align clusters of the same set of samples and provide a plot of their shared values
}

\usage{
plotTracking<-function(clusters, index=NULL,reuseColors=FALSE,matchToTop=FALSE,plot=TRUE,unassignedColor="white",minRequireColor=0.3,startNewColors=FALSE,...)
showThisPal()
clusterTrackingPlot(m, names=rownames(m),add=FALSE,x=NULL,ylim=NULL,tick=FALSE,ylab="",xlab="",...)
plotWithClinical(clusters,clinical,clinicalNumbers=NULL,index=NULL,whClinical=NULL,main="",split=TRUE,mar=c(2.1,7.1,1.1,1.1),plot=TRUE,...)
alignClusters(clAlign,clRef=NULL,index=NULL)
printWithIndex(x,whBlockRow=1,blockOrder=NULL,whRow=NULL,whBlocks=NULL,checkAllIndices=TRUE)
clusterIndInBlock(x,whInd,whBlockRow=1,whBlocks,rightLeftPerBlock=NULL)
bigPalette
}

\arguments{
\item{clusters}{
A matrix of with each row corresponding to a clustering and each column a sample. The plotting will plot them in order of this matrix, and their order influences the plot greatly. 
}
\item{index}{A predefined order in which the samples will be plotted. Otherwise the order is found internally
}
\item{reuseColors}{Logical. Whether each row should consist of the same set of colors. By default (FALSE) each cluster that the algorithm doesn't identify to the previous rows clusters gets a new color}
\item{matchToTop}{Logical as to whether all clusters should be aligned to the first row. By default (FALSE) each cluster is aligned to the ordered clusters of the row above it. }
\item{plot}{Logical as to whether a plot should be produced.}
\item{unassignedColor}{If ``-1'' in clusters, will be given this color.}
\item{minRequireColor}{In aligning colors between rows of clusters, require this percent overlap}
\item{startNewColors}{In aligning colors between rows of clusters, restart at beginning of .thisPal as long as colors not in immediately proceeding (some of the colors at the end of .thisPal are a bit wonky, and so if you have a large clusters matrix, this can be useful)}
\item{m}{Matrix where rows are different clusterings, columns are samples, and entries are colors (i.e. characters that define colors)}
\item{names}{names to go with the rows (clusterings) in matrix m}
\item{clinical}{matrix of other (factors) to be plotted with the plotTracking. Assumed to be rows are samples, columns are clinical factors.}
\item{clinicalNumbers}{Numbers to assign clinical data to colors in thisPal vector (otherwise found automatically)}
\item{whClinical}{which columns (factors) to use. Use to exclude non-factor data (continuous, etc.)}
\item{split}{Logical, whether to split the clinical around the clusters for ease in matching}
\item{clAlign}{clustering to align to clRef}
\item{clRef}{clustering to align against}
\item{x}{output of plotTracking}
\item{whBlockRow}{which row should be used for defining blocks}
\item{blockOrder}{new order for blocks (NULL doesn't change order)}
\item{whRow}{}
\item{whBlocks}{}
\item{checkAllIndices}{Logical, whether throw error if new ordering given in blockOrder doesn't contains all indices}
\item{whInd}{indices that you want to cluster together within blocks}
\item{whBlocks}{which blocks you want to cluster the indices within}
\item{rightLeftPerBlock}{whether cluster them to right or left of block.}
}

\details{
 clusterTrackingPlot is generally called within plotTracking, but is provided as separate function for convenience, for example if you want to change a color and replot.  printWithIndex and clusterIndInBlock are helper functions to try and move around the order in the plotTracking to look prettier. 
 
 bigPalette is a long palette of colors (length 62) used by plotTracking and accompanying functions.
}

\value{
plotTracking provides (invisibly) the orders and other things that go into making the matrix.
}

\author{
Elizabeth Purdom and Marla Johnson
}

\examples{

}
