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
#' Large palette of colors 
#' @name bigPalette
#' @aliases bigPalette showThisPal
#' @details \code{bigPalette} is a long palette of colors (length 62) used by \code{\link{plotTracking}} and
#' accompanying functions. \code{showThisPal} creates plot that gives index of each color in bigPalette.
bigPalette<-.thisPal

#' @rdname bigPalette
showThisPal<-function(){
	plot(1:length(.thisPal),y=1:length(.thisPal),pch=19,col=.thisPal,cex=3)
	text(1:length(.thisPal),x=1:length(.thisPal),y=1:length(.thisPal),pos=1)

}


#' Set of colors useful for heatmap gradients 
#' @name showHeatmapPalettes
#' @aliases showHeatmapPalettes seqPal1 seqPal2 seqPal3 seqPal4 seqPal5
#' @details seqPal1-seqPal4 are palettes for the heatmap. showHeatmapPalettes() will show you
#' these palettes.
#'
#' @examples
#' #show the palette colors
#' showHeatmapPalettes()
#' #compare the palettes on heatmap
#' data(simData)
#' data(simData)
#' cl<-clusterAll(simData,clusterFunction="pam",subsample=FALSE,
#' sequential=FALSE, clusterDArgs=list(k=8))$cl
#' par(mfrow=c(2,3))
#' dualHeatmap(cl,heatData=simCount,clusterData=simData,colorScale=seqPal1,main="seqPal1")
#' dualHeatmap(cl,heatData=simCount,clusterData=simData,colorScale=seqPal2,main="seqPal2")
#' dualHeatmap(cl,heatData=simCount,clusterData=simData,colorScale=seqPal3,main="seqPal3")
#' dualHeatmap(cl,heatData=simCount,clusterData=simData,colorScale=seqPal4,main="seqPal4")
#' dualHeatmap(cl,heatData=simCount,clusterData=simData,colorScale=seqPal5,main="seqPal5")
#' 


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
	plotTracking(mat,input="colors")
}

#' @rdname showHeatmapPalettes
seqPal5<- colorRampPalette(c("black","navyblue","mediumblue","dodgerblue3","aquamarine4","green4","yellowgreen","yellow"))(16)
#' @rdname showHeatmapPalettes
seqPal2<- colorRampPalette(c("orange","black","blue"))(16)
seqPal2<-(c("yellow","gold2",seqPal2))
seqPal2<-rev(seqPal2)
#' @rdname showHeatmapPalettes
seqPal3<-rev(brewer.pal(11, "RdBu"))
#' @rdname showHeatmapPalettes
seqPal4<-colorRampPalette(c("black","blue","white","red","orange","yellow"))(16)
#' @rdname showHeatmapPalettes
seqPal1<-rev(brewer.pal(11, "Spectral"))
