#' Convert clusterLegend into useful formats
#'
#' Function for converting the information stored in the clusterLegend slot into other useful formats
#' @param object a ClusterExperiment object
#' @param output character value, indicating desired type of conversion.
#' @details convertClusterLegend pulls out information stored in the clusterLegend slot of the object and returns it in useful format.
#' @return If "output=plotAndLegend", \code{convertClusterLegend} will return a list that provides the necessary information to color samples according to cluster and create a legend for it:
#'  \itemize{
#'  \item{"colorVector"}{A vector the same length as the number of samples, assigning a color to each cluster of the primaryCluster of the object}
#'  \item{"legendNames"}{A vector the length of the number of clusters of primaryCluster of the object giving the name of the cluster.}
#'  \item{"legendColors"}{A vector the length of the number of clusters of primaryCluster of the object giving the color of the cluster.}
#' }
#'  @return If "output=aheatmap" a conversion of the clusterLegend to be in the format requested by aheatmap. The column 'name' is used for the names and the column 'color' for the color of the clusters.
#'  @return If "output=matrixNames" or "matrixColors" a matrix the same dimension of clusterMatrix(object), but with the cluster color or cluster name instead of the clusterIds, respectively.
#'  @importFrom RColorBrewer brewer.pal brewer.pal.info
setMethod(
  f = "convertClusterLegend",
  signature = c("ClusterExperiment"),
  definition = function(object,output=c("plotAndLegend","aheatmapFormat","matrixNames","matrixColors")){
    output<-match.arg(output)
    if(output=="aheatmapFormat"){
      outval<-.convertToAheatmap(clusterLegend(object))
    }
    if(output%in% c("matrixNames","matrixColors")){
      outval<-do.call("cbind",lapply(1:nClusters(object),function(ii){
        cl<-clusterMatrix(object)[,ii]
        colMat<-clusterLegend(object)[[ii]]
        m<-match(cl,colMat[,"clusterIds"])
        colReturn<-if(output=="matrixNames") "name" else "color"
        return(colMat[m,colReturn])
      }))

    }
    if(output=="plotAndLegend"){
      cl<-primaryCluster(object)
      colMat<-clusterLegend(object)[[primaryClusterIndex(object)]]
      clColor<-colMat[match(cl,colMat[,"clusterIds"]),"color"]
      legend<-colMat[,"name"]
      color<-colMat[,"color"]
      outval<-list(colorVector=clColor,legendNames=legend,legendColors=color)

    }
    return(outval)

  }
)

.convertToAheatmap<-function(clusterLegend){
    outval<-lapply(clusterLegend,function(x){
      z<-x[,"color"]
      names(z)<-x[,"clusterIds"]
      z<-z[order(names(z))]
      return(z)
    })
    return(outval)

  }



#' @param breaks either vector of breaks, or number of breaks (integer) or a number between 0 and 1 indicating a quantile, between which evenly spaced breaks should be calculated
#' @rdname plottingFunctions
#' @details \code{setBreaks} gives a set of breaks (of length 52) equally spaced between the boundaries of the data. If breaks is between 0 and 1, then the evenly spaced breaks are between these quantiles of the data.
#' @examples 
#' data(simData)
#' setBreaks(.9,simData)
setBreaks<-function(breaks,data){
    if(length(breaks)>0 && !is.na(breaks)){
        #get arround bug in aheatmap
        #if colors are given, then get back 51, unless give RColorBrewer, in which case get 101! Assume user has to give palette.
        #TO DO: might not need any more with updated aheatmap.
        if(length(breaks)==1){
            if(breaks<=1){
                ncols<-51
                if(breaks<1) breaks<-c(seq(min(data),quantile(data[data>0],breaks,na.rm=TRUE),length=ncols),max(data))
                else breaks<-seq(min(data),max(data),length=ncols+1)
            }
            else{
                warning("Because of bug in aheatmap, breaks should be of length 52 -- otherwise the entire spectrum will not be used. We don't recommend that you set the breaks to a integer number, but let aheatmap determine the breaks")
            }
        }
        else{
            if(length(breaks)!=52) warning("Because of bug in aheatmap, breaks should be of length 52 -- otherwise the entire spectrum will not be used")
        }
    }
    return(breaks)

}


.thisPal = c(
	"#A6CEE3",#light blue
	"#1F78B4",#dark blue
	"#B2DF8A",#light green
	"#33A02C",#dark green
	"#FB9A99",#pink
	"#E31A1C",#red
	"#FDBF6F",#light orange
	"#FF7F00",#orange
	"#CAB2D6",#light purple
	"#6A3D9A",#purple
	"#FFFF99",#light yellow
	"#B15928",#brown
	"#bd18ea", #magenta
	"#2ef4ca", #aqua
	"#f4cced", #pink,
	"#05188a", #navy,
	"#f4cc03", #lightorange
	"#e5a25a", #light brown
	"#06f106", #bright green
	"#85848f", #med gray
	"#000000", #black
	"#076f25", #dark green
	"#93cd7f",#lime green
	"#4d0776", #dark purple
	"maroon3",
	"blue",
	"grey"
		)
.thisPal<-.thisPal[c(2,4,6,8,10,12,14,13,1,3,5,7,9,11,16,19,20,15,17,18,21:27)]
.brewers<-RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info[,"category"]=="qual",]
.brewers<-.brewers[c(1:3,6:8,4:5),]
.thisPal<-c(.thisPal,palette(),unlist(sapply(1:nrow(.brewers),function(ii){RColorBrewer::brewer.pal(.brewers[ii,"maxcolors"], rownames(.brewers)[ii])})))
.thisPal<-unique(.thisPal)
.thisPal<-.thisPal[-c(32,34,36,37,40,45:47,49:53,56,62:71,73,75,76,84,90,92 )] #remove because too similar to others
.thisPal<-.thisPal[-34] #very similar to 2
.thisPal<-.thisPal[-31] #very similar to 7

#' @rdname plottingFunctions
bigPalette<-.thisPal

#' @title Various functions useful for plotting
#' @description Most of these functions are called internally by plotting functions, but are exported in case the user finds them useful.
#' @name plottingFunctions
#' @aliases bigPalette showBigPalette
#' @details \code{bigPalette} is a long palette of colors (length 62) used by
#' \code{\link{plotClusters}} and accompanying functions. \code{showBigPalette}
#' creates plot that gives index of each color in bigPalette.
#' @rdname plottingFunctions
#' @details \code{showBigPalette} will plot the \code{bigPalette} functions with their labels and index
showBigPalette<-function(wh=NULL){
  oldMar<-par("mar")
  if(is.null(wh)){
    col<-.thisPal
    wh<-1:length(col)
  }
  else{ col<-.thisPal[wh]}
  par(mar=c(2.1,2.1,2.1,2.1))
	plot(1:length(col),y=1:length(col),pch=19,col=col,cex=3,xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
	text(as.character(wh),x=1:length(col),y=1:length(col),pos=1,xpd=NA)
	text(as.character(col),x=1:length(col),y=1:length(col),pos=1,offset=1.5,xpd=NA)
	par(mar=oldMar)
}
#' @examples 
#' showPalette()


#' @rdname plottingFunctions
#' @details \code{seqPal1}-\code{seqPal4} are palettes for the heatmap. \code{showHeatmapPalettes} will show you these palettes.
#' @examples
#' #show the palette colors
#' showHeatmapPalettes()
#' #compare the palettes on heatmap
#' data(simData)
#' data(simData)
#' cl<-clusterSingle(simData,clusterFunction="pam",subsample=FALSE,
#' sequential=FALSE, clusterDArgs=list(k=8))$cl
#' par(mfrow=c(2,3))
#' plotHeatmap(cl,heatData=simCount,clusterData=simData,colorScale=seqPal1,main="seqPal1")
#' plotHeatmap(cl,heatData=simCount,clusterData=simData,colorScale=seqPal2,main="seqPal2")
#' plotHeatmap(cl,heatData=simCount,clusterData=simData,colorScale=seqPal3,main="seqPal3")
#' plotHeatmap(cl,heatData=simCount,clusterData=simData,colorScale=seqPal4,main="seqPal4")
#' plotHeatmap(cl,heatData=simCount,clusterData=simData,colorScale=seqPal5,main="seqPal5")
showHeatmapPalettes<-function(){
	palettesAll<-list(seqPal1=seqPal1,seqPal2=seqPal2,seqPal3=seqPal3,seqPal4=seqPal4,seqPal5=seqPal5)
	maxLength<-max(sapply(palettesAll,length))
	palettesAllAdj<-lapply(palettesAll,function(x){
		if(length(x)<maxLength) x<-c(x,rep("white",maxLength-length(x)))
			return(x)})
	ll<-list()
	sapply(1:length(palettesAllAdj),function(ii){ll[[2*ii-1]]<<-palettesAllAdj[[ii]]})
	sapply(1:(length(palettesAllAdj)-1),function(ii){ll[[2*ii]]<<-rep("white",length=maxLength)})
	names(ll)[seq(1,length(ll),by=2)]<-names(palettesAll)
	names(ll)[seq(2,length(ll),by=2)]<-rep("",length(palettesAll)-1)
	mat<-do.call("cbind",ll)
	plotClusters(mat,input="colors")
}

#' @rdname plottingFunctions
seqPal5<- colorRampPalette(c("black","navyblue","mediumblue","dodgerblue3","aquamarine4","green4","yellowgreen","yellow"))(16)
#' @rdname plottingFunctions
seqPal2<- colorRampPalette(c("orange","black","blue"))(16)
seqPal2<-(c("yellow","gold2",seqPal2))
seqPal2<-rev(seqPal2)
#' @rdname plottingFunctions
seqPal3<-rev(brewer.pal(11, "RdBu"))
#' @rdname plottingFunctions
seqPal4<-colorRampPalette(c("black","blue","white","red","orange","yellow"))(16)
#' @rdname plottingFunctions
seqPal1<-rev(brewer.pal(11, "Spectral"))
