#' Convert clusterLegend into useful formats
#'
#' Function for converting the information stored in the clusterLegend slot into
#' other useful formats.
#'
#' @param object a \code{ClusterExperiment} object.
#' @param output character value, indicating desired type of conversion.
#' @inheritParams getClusterIndex
#' @details convertClusterLegend pulls out information stored in the
#'   \code{clusterLegend} slot of the object and returns it in useful format.
#'
#' @return If \code{output="plotAndLegend"}, \code{"convertClusterLegend"} will
#'   return a list that provides the necessary information to color samples
#'   according to cluster and create a legend for it:
#'  \itemize{
#'  \item{"colorVector"}{ A vector the same length as the number of samples,
#'  assigning a color to each cluster of the primaryCluster of the object.}
#'  \item{"legendNames"}{ A vector the length of the number of clusters of
#'  primaryCluster of the object giving the name of the cluster.}
#'  \item{"legendColors"}{ A vector the length of the number of clusters of
#'  primaryCluster of the object giving the color of the cluster.}
#' }
#' @return If \code{output="aheatmap"} a conversion of the clusterLegend to be
#'   in the format requested by \code{\link[NMF]{aheatmap}}. The column 'name'
#'   is used for the names and the column 'color' for the color of the clusters.
#' @return If \code{output="matrixNames"} or \code{"matrixColors"} a matrix the
#'   same dimension of \code{clusterMatrix(object)}, but with the cluster color
#'   or cluster name instead of the clusterIds, respectively.
#'
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @export
#' @rdname plottingFunctions
#' @aliases convertClusterLegend convertClusterLegend,ClusterExperiment-method
setMethod(
  f = "convertClusterLegend",
  signature = c("ClusterExperiment"),
  definition = function(object,output=c("plotAndLegend","aheatmapFormat","matrixNames","matrixColors"),whichClusters=ifelse(output=="plotAndLegend","primary","all")){
    output<-match.arg(output)
    whichClusters<-getClusterIndex(object,whichClusters=whichClusters,noMatch="throwError")
    if(output=="aheatmapFormat"){
      outval<-.convertToAheatmap(clusterLegend(object)[whichClusters])
    }
    if(output%in% c("matrixNames","matrixColors")){
      outval<-do.call("cbind",lapply(seq_along(whichClusters),function(ii){
        cl<-clusterMatrix(object)[,whichClusters,drop=FALSE][,ii]
        colMat<-clusterLegend(object)[whichClusters][[ii]]
        m<-match(cl,colMat[,"clusterIds"])
        colReturn<-if(output=="matrixNames") "name" else "color"
        return(colMat[m,colReturn])
      }))
      colnames(outval)<-clusterLabels(object)[whichClusters]
      
    }
    if(output=="plotAndLegend"){
      if(length(whichClusters)>1) stop("given whichClusters indicates more than 1 clustering which is not allowed for option 'plotAndLegend'")
      cl<-clusterMatrix(object)[,whichClusters]
      colMat<-clusterLegend(object)[[whichClusters]]
      clColor<-colMat[match(cl,colMat[,"clusterIds"]),"color"]
      legend<-colMat[,"name"]
      color<-colMat[,"color"]
      outval<-list(colorVector=clColor,legendNames=legend,legendColors=color)
      
    }
    return(outval)
    
  }
)

.convertToClusterLegend<-function(acol){
	if(!is.list(acol)){ #single vector
		acol<-list(acol)
	}
	out<-lapply(acol,function(x){
		if(!is.null(names(x))) nms<-names(x)
		else nms<-as.character(seq_along(x))
		cols<-unname(x)
		return(cbind("color"=cols,"name"=nms))
	})
	names(out)<-names(acol)
	return(out)
}

.convertToAheatmap<-function(clusterLegend, names=FALSE){
  outval<-lapply(clusterLegend,function(x){
    if(!is.null(dim(x))){
      z<-x[,"color"]
      if(names) {
        names(z)<-x[,"name"]
      } else {
        names(z)<-x[,"clusterIds"]
      }
      z<-z[order(names(z))]
      return(z)
    }
    else return(x)
  })
  return(outval)
  
}

#' @title Various functions useful for plotting
#' @description Most of these functions are called internally by plotting
#'   functions, but are exported in case the user finds them useful.
#' @name plottingFunctions
#' @aliases bigPalette showPalette
#' @details \code{bigPalette} is a long palette of colors (length 58) used by 
#'   \code{\link{plotClusters}} and accompanying functions. \code{showPalette}
#'   creates plot that gives index of each color in a vector of colors.
#'   \code{massivePalette} is a combination of \code{bigPalette} and the
#'   non-grey colors of \code{\link{colors}()} (length 487). 
#'   \code{massivePalette} is mainly useful for when doing
#'   \code{\link{plotClusters}} of a very large number of clusterings, each with
#'   many clusters, so that the code doesn't run out of colors. However, many of
#'   the colors will be very similar to each other.
#' @param which numeric. Which colors to plot. Must be a numeric vector with
#'   values between 1 and length of \code{colPalette}. If missing, all colors
#'   plotted.
#' @param cex numeric value giving the cex for the text of the plot.
#' @param colPalette a vector of character colors. By default, the palette 
#'   \code{bigPalette} is used
#' @details \code{showPalette} will plot the \code{colPalette} colors with their
#'   labels and index.
#'
#' @export
#'
#' @examples
#' showPalette()
#' showPalette(massivePalette,cex=0.6)
showPalette<-function(colPalette=bigPalette,which=NULL,cex=1){
  oldPar<-par(no.readonly = TRUE)
  wh<-which
  if(is.null(wh)){
    wh<-seq_along(colPalette)
  }
  else{ colPalette<-colPalette[wh]}
  n<-ceiling(sqrt(length(colPalette)))
  nblank<-n^2-length(colPalette)
  xwide<-n
  yup<-n
  x1<-rep(c(seq_len(xwide))-.5,yup)
  x2<-rep(c(seq_len(xwide))+.5,yup)
  xtext<-rep(c(seq_len(xwide)),yup)
  ycolor1<-rep(seq(1,yup*2,by=2)-.5,each=xwide)
  ycolor2<-rep(seq(1,yup*2,by=2)+.5,each=xwide)
  ytext<-rep(seq(2,yup*2,by=2)+.5,each=xwide)

  par(mar=c(0,0,0,0),omi=c(0,0,0,0))
  plot.new()
  plot.window(xlim=c(.5,xwide+.5),ylim=c(.5,(yup*2)+.5))
  rect(x1,ycolor1,x2,ycolor2,col=c(colPalette,rep("white",nblank)),border=FALSE)
  if(length(colPalette)>100){
	  half<-ceiling(length(colPalette)/2)
	 adj.text<-cbind(rep(.5,half*2),rep(c(0,1),times=half))
	 adj.text<-adj.text[seq_along(colPalette),]
  }
  else adj.text<-matrix(c(0.5,0),nrow=length(colPalette),ncol=2,byrow=TRUE)
  for(i in seq_along(colPalette)){
      text(xtext[i],ytext[i]-1,colPalette[i],cex=cex,adj=adj.text[i,])
      if(length(colPalette)<=100) text(xtext[i],ytext[i]-2,wh[i],cex=cex,adj=c(0.5,1))
  }
	par(oldPar)
}

#' @rdname plottingFunctions
#' @export
bigPalette<-c(
	'#E31A1C',
	'#1F78B4',
	'#33A02C',
	'#FF7F00',
	'#6A3D9A',
	'#B15928',
	'#A6CEE3',
	'#bd18ea',
	'cyan',
	'#B2DF8A',
	'#FB9A99',
	"deeppink4",
	'#00B3FFFF',
	'#CAB2D6',
	'#FFFF99',
	'#05188a',
	'#CCFF00FF',
	'cornflowerblue',
	'#f4cc03',
	'black',
	'blueviolet',
	'#4d0776',
	'maroon3',
	'blue',
#	'grey',
	'#E5D8BD',
	'cadetblue4',
	'#e5a25a',
	"lightblue1",
	'#F781BF',
	'#FC8D62',
	'#8DA0CB',
	'#E78AC3',
	'green3',
	'#E7298A',
	'burlywood3',
	'#A6D854',
	"firebrick",
	'#FFFFCC',
	"mediumpurple",
	'#1B9E77',
	'#FFD92F',
	'deepskyblue4',
	"yellow3",
	'#00FFB2FF',
	'#FDBF6F',
	'#FDCDAC',
	"gold3",
	'#F4CAE4',
	'#E6F5C9',
	'#FF00E6FF',
	'#7570B3',
	"goldenrod",
	'#85848f',
	"lightpink3",
	"olivedrab",
#	"plum",
#	"lightskyblue3",
#	"mediumturquoise",
	'cadetblue3'
)

#' @importFrom grDevices colors
.rcolors<-function(){
	set.seed(23589)
	x<-sample(colors()[-c(152:361)])
	return(x)
}

##If want to keep the same colors could:
# if(!exists(".Random.seed", envir = .GlobalEnv)) {
#     message("calling runif(1)"); runif(1) }
# old.R.s <- .Random.seed
# ## will reset everything on exiting this function:
# on.exit(assign(".Random.seed", old.R.s, envir=.GlobalEnv))
# ## set seed for sample() "back compatibly":
# suppressWarnings(RNGversion("3.5.0"))
# set.seed(seed)
# ## return random permutation of "my colors"
# sample(colors()[-c(152:361)])

#' @rdname plottingFunctions
#' @export
massivePalette<-unique(c(bigPalette,.rcolors()))




#' @param breaks either vector of breaks, or number of breaks (integer) or a
#'   number between 0 and 1 indicating a quantile, between which evenly spaced
#'   breaks should be calculated. If missing or NA, will determine evenly spaced
#'   breaks in the range of the data.
#' @param makeSymmetric whether to make the range of the breaks symmetric around zero (only used if not all of the data is non-positive and not all of the data is non-negative)
#' @param returnBreaks logical as to whether to return the vector of breaks. See details.
#' @details if returnBreaks if FALSE, instead of returning the vector of breaks, the function will just return the second smallest and second largest value of the breaks. This is useful for alternatively just setting values of the data matrix larger than these values to this value if breaks was a percentile. This argument is only used if \code{breaks<1}, indicating truncating the breaks for large values of data.
#' @rdname plottingFunctions
#' @details \code{setBreaks} gives a set of breaks (of length 52) equally spaced
#'   between the boundaries of the data. If breaks is between 0 and 1, then the
#'   evenly spaced breaks are between these quantiles of the data.
#' @export
#' @examples
#' setBreaks(data=simData,breaks=.9)
setBreaks<-function(data,breaks=NA,makeSymmetric=FALSE,returnBreaks=TRUE){
  if(all(is.na(data))) stop("data consists only of NA values")
  if(length(unique(na.omit(as.vector(data))))==1){
    warning("data consists of only a single non NA value")
    val<-unique(na.omit(as.vector(data)))
    return(seq(val-1,val+1,length=52))
  }
  isPositive<-all(na.omit(as.vector(data))>=0)
  isNegative<-all(na.omit(as.vector(data))<=0)
  #
  #get arround bug in aheatmap
  #if colors are given, then get back 51 colors, unless give RColorBrewer, in which case get 101! Assume user has to give palette. So breaks has to be +1 of that length
  #TO DO: might not need any more with updated aheatmap.
  ncols<-51
  minData<-min(data,na.rm=TRUE)
  maxData<-max(data,na.rm=TRUE)
  maxAbsData<-max(abs(data),na.rm=TRUE)
  if(!is.vector(breaks)) stop("breaks argument must be a vector")
  if(missing(breaks) || is.na(breaks)){
    #go from minimum to max
    if(makeSymmetric & !isPositive & !isNegative){
      breaks<-seq(-maxAbsData,maxAbsData,length=ncols+1)
      seconds<-c(-maxAbsData,maxAbsData)
    }
    else{
      breaks<-seq(minData,maxData,length=ncols+1)
      seconds<-c(minData,maxData)
    }
    
  }
  else if(length(breaks)>0 && !is.na(breaks)){
    if(length(breaks)==1){
      if(breaks<1){
        if(breaks<0.5) breaks<-1-breaks
        #
        uppQ<-if(isPositive) quantile(data[which(data>0)],breaks,na.rm=TRUE) else quantile(data,breaks,na.rm=TRUE)
        lowQ<-if(isPositive) min(data,na.rm=TRUE) else quantile(data,1-breaks,na.rm=TRUE)
        
        if(makeSymmetric & !isPositive & !isNegative){
          #make breaks symmetric around 0
          absq<-max(abs(c(lowQ,uppQ)),na.rm=TRUE)
          absm<-max(abs(c(min(data,na.rm=TRUE),max(data,na.rm=TRUE))))
          #is largest quantile also max of abs(data)?
          quantAllMax <- if( isTRUE( all.equal(round(absq,5), round(absm,5)))) TRUE else FALSE
          if(!quantAllMax){
            breaks <- c(-absm, seq(-absq,absq,length=ncols-1), absm)
            seconds<-c(-absq,absq)
          }
          else{
            #equally spaced
            breaks <- seq(-absm,absm,length=ncols+1)
            seconds<-c(-absm,absm)
          }
        }
        else{
          #determine if those quantiles are min/max of data
          quantMin <- if( isTRUE( all.equal(round(lowQ,5), round(minData,5)))) TRUE else FALSE
          quantMax<-if( isTRUE( all.equal(round(uppQ,5),round(maxData,5)))) TRUE else FALSE
          
          if(!quantMin & !quantMax){
            breaks <- c(minData, seq(lowQ,uppQ,length=ncols-1), maxData)
            seconds<-c(lowQ,uppQ)
          }
          if(!quantMin & quantMax){
            breaks <- c(minData, seq(lowQ,maxData,length=ncols))
            seconds<-c(lowQ,maxData)
          }
          if(quantMin & !quantMax){
            breaks <- c(seq(minData,uppQ,length=ncols), maxData)
            seconds<-c(minData,uppQ)
          }
          if(quantMin & quantMax){
            breaks<-seq(minData,maxData,length=ncols+1)
            seconds<-c(minData,maxData)
          }
        }
      }
      else{ #breaks is number of breaks
        if(length(breaks)!=52) warning("Because of bug in aheatmap, breaks should be of length 52 -- otherwise the entire spectrum of colors will not be used")
        if(makeSymmetric& !isPositive & !isNegative){
          breaks<-seq(-maxAbsData,maxAbsData,length=breaks)
          seconds<-c(-maxAbsData,maxAbsData)
        }
        else{
          breaks<-seq(minData,maxData,length=breaks)
          seconds<-c(minData,maxData)
        }
      }
    }
  }
  if(!length(unique(breaks))==length(breaks)){
    breaks<-sort(unique(breaks))
    warning("setBreaks did not create unique breaks, and therefore the resulting breaks will not be of exactly the length requested. If you have replicable code that can recreate this warning, we would appreciate submitting it to the issues tracker on github (existing issue: #186)")
    seconds<-c(min(breaks),max(breaks))
  }
  if(returnBreaks) return(breaks)
  else return(seconds)
  
}



#' @rdname plottingFunctions
#'
#' @details \code{seqPal1}-\code{seqPal4} are palettes for the heatmap.
#'   \code{showHeatmapPalettes} will show you these palettes.
#'
#' @export
#'
#' @examples
#'
#' #show the palette colors
#' showHeatmapPalettes()
#'
#' #compare the palettes on heatmap
#' cl <- clusterSingle(simData, subsample=FALSE,
#' sequential=FALSE, mainClusterArgs=list(clusterFunction="pam", clusterArgs=list(k=8)))
#'
#' \dontrun{
#' par(mfrow=c(2,3))
#' plotHeatmap(cl, colorScale=seqPal1, main="seqPal1")
#' plotHeatmap(cl, colorScale=seqPal2, main="seqPal2")
#' plotHeatmap(cl, colorScale=seqPal3, main="seqPal3")
#' plotHeatmap(cl, colorScale=seqPal4, main="seqPal4")
#' plotHeatmap(cl, colorScale=seqPal5, main="seqPal5")
#' par(mfrow=c(1,1))
#' }
#'
showHeatmapPalettes<-function(){
	palettesAll<-list(seqPal1=seqPal1,seqPal2=seqPal2,seqPal3=seqPal3,seqPal4=seqPal4,seqPal5=seqPal5)
	maxLength<-max(sapply(palettesAll,length))
	palettesAllAdj<-lapply(palettesAll,function(x){
		if(length(x)<maxLength) x<-c(x,rep("white",maxLength-length(x)))
			return(x)})
	ll<-list()
	sapply(seq_along(palettesAllAdj),function(ii){ll[[2*ii-1]]<<-palettesAllAdj[[ii]]})
	sapply(seq_len(length(palettesAllAdj)-1),function(ii){ll[[2*ii]]<<-rep("white",length=maxLength)})
	names(ll)[seq(1,length(ll),by=2)]<-names(palettesAll)
	names(ll)[seq(2,length(ll),by=2)]<-rep("",length(palettesAll)-1)
	mat<-do.call("cbind",ll)
	plotClusters(mat,input="colors")
}

#' @rdname plottingFunctions
#' @export
seqPal5<- colorRampPalette(c("black","navyblue","mediumblue","dodgerblue3","aquamarine4","green4","yellowgreen","yellow"))(16)
#' @rdname plottingFunctions
#' @export
seqPal2<-colorRampPalette(c("yellow","orange","black","blue","dodgerblue"))(16)[-c(2,14)]
seqPal2<-rev(seqPal2)
#seqPal2<- colorRampPalette(c("yellow","orange","black","blue","dodgerblue"))(16)
#seqPal2<-(c("yellow","gold2",seqPal2))
#' @rdname plottingFunctions
#' @export
seqPal3<-rev(RColorBrewer::brewer.pal(11, "RdBu"))
#' @rdname plottingFunctions
#' @export
#make it symmetric around white
seqPal4<-colorRampPalette(c("black","blue","white","red","orange","yellow"))(16)[c(1:9,12,13,14,16)]
#' @rdname plottingFunctions
#' @export
seqPal1<-rev(RColorBrewer::brewer.pal(11, "Spectral"))


#' @export
#' @param title title for the clusterLegend plot
#' @param clusterNames vector of names for the clusters; vector should have names
#'  that correspond to the clusterIds in the ClusterExperiment object. If this
#'  argument is missing, will use the names in the "name" column of the clusterLegend
#'  slot of the object.
#' @param add logical. Whether legend should be added to the existing plot.
#' @param location character passed to \code{x} argument of legend indicating 
#'  where to place legend.
#' @param ... arguments passed to legend
#' @inheritParams getClusterIndex
#' @rdname plottingFunctions
#' @aliases plotClusterLegend
setMethod(
  f = "plotClusterLegend",
  signature = c("ClusterExperiment"),
  definition = function(object,whichCluster="primary",clusterNames,title,add=FALSE,location=if(add)"topright" else "center",...){
    whichCluster<-getSingleClusterIndex(object,whichCluster,list(...))
    legMat<-clusterLegend(object)[[whichCluster]]
    if(!missing(clusterNames)){
      if(is.null(names(clusterNames))) stop("clusterNames must be named vector")
      m<-match(legMat[,"clusterIds"],names(clusterNames))
      if(any(is.na(m))){
        whClusters<-which(legMat[,"clusterIds"]>0)
        if(any(is.na(m[whClusters]))) stop("not all clusters are found in names of clusterNames argument")
        else{ #give default names to -1 and -2
          whMiss1<-which(legMat[,"clusterIds"]== -1 & is.na(m))
          if(length(whMiss1)>0) clusterNames<-c(clusterNames,"-1"="Not assigned")
          whMiss2<-which(legMat[,"clusterIds"]== -2 & is.na(m))
          if(length(whMiss2)>0) clusterNames<-c(clusterNames,"-2"="Not clustered")
        }
      }
      m<-match(legMat[,"clusterIds"],names(clusterNames))
      clusterNames<-clusterNames[m]
    }
    else{
      clusterNames<-legMat[,"name"]
    }
    if(missing(title)){
      title<-paste("Clusters of",clusterLabels(object)[whichCluster])
    }
    if(any(as.numeric(legMat[,"clusterIds"])<0)){
      #put -1/-2 last
      isNeg<-as.numeric(legMat[,"clusterIds"])<0
      ord<-c(which(!isNeg),which(isNeg))
    }
    else ord<-seq_len(nrow(legMat))
    if(!add){
			graphics::plot(0,0,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
			graphics::legend(x=location,legend=clusterNames[ord],fill=legMat[ord,"color"],title=title,...)
		}
		else{
			graphics::legend(x=location,legend=clusterNames[ord],fill=legMat[ord,"color"],title=title,...)
		}
    
    
  })


 ###Old definition of bigPalette (.thisPal):
 # bigPalette = c(
 # 	"#A6CEE3",#1 light blue
 # 	"#1F78B4",#2 dark blue
 # 	"#B2DF8A",#3 light green
 # 	"#33A02C",#4 dark green
 # 	"#FB9A99",#5 pink
 # 	"#E31A1C",#6 red
 # 	"#FDBF6F",#7 light orange
 # 	"#FF7F00",#8 orange
 # 	"#CAB2D6",#9 light purple
 # 	"#6A3D9A",#10 purple
 # 	"#FFFF99",#11 light yellow
 # 	"#B15928",#12 brown
 # 	"#bd18ea", #13 magenta
 # 	"#2ef4ca", #14 aqua
 # 	"#f4cced", #15 pink,
 # 	"#05188a", #16 navy,
 # 	"#f4cc03", #17 lightorange
 # 	"#e5a25a", #18 light brown
 # 	"#06f106", #19 bright green
 # 	"#85848f", #20 med gray
 # 	"#000000", #21 black
 # 	"#076f25", #22 dark green
 # 	"#93cd7f",#23 lime green
 # 	"#4d0776", #24 dark purple
 # 	"maroon3",
 # 	"blue",
 # 	"grey"
 # 		)
 # bigPalette<-bigPalette[c(2,4,6,8,10,12,14,13,1,3,5,7,9,11,16,19,20,15,17,18,21:27)]
 # .brewers<-RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info[,"category"]=="qual",]
 # .brewers<-.brewers[c(1:3,6:8,4:5),]
 # bigPalette<-c(bigPalette,palette()[-1],unlist(sapply(1:nrow(.brewers),function(ii){RColorBrewer::brewer.pal(.brewers[ii,"maxcolors"], rownames(.brewers)[ii])})))
 # bigPalette<-unique(bigPalette)
 # bigPalette<-bigPalette[-c(32,34,36,37,40,45:47,49:53,56,62:71,73,75,76,84,90,92 )] #remove because too similar to others
 # bigPalette<-bigPalette[-34] #very similar to 2
 # bigPalette<-bigPalette[-31] #very similar to 7
