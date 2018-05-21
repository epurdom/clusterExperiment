#' @rdname plotContrastHeatmap
#' @aliases plotContrastHeatmap
#' @title Plot heatmaps showing significant genes per contrast
#' @description Plots a heatmap of the data, with the genes grouped based on the
#'   contrast for which they were significant.
#' @param object ClusterExperiment object on which biomarkers were found
#' @param signifTable A \code{data.frame} in format of the result of
#'   \code{\link{getBestFeatures}}. It must minimally contain columns 'Contrast'
#'   and 'IndexInOriginal' giving the grouping and original index of the
#'   features in the \code{assay(object)}
#' @param whichCluster if not NULL, indicates cluster used in making the 
#'   significance table. Used to match to names in \code{clusterLegend(object)} 
#'   (relevant for one-vs-all contrast so that color aligns).
#' @param contrastColors vector of colors to be given to contrasts. Should match
#'   the name of the contrasts in the 'Contrast' column of \code{signifTable} or
#'   'ContrastName', if given.. If  missing, default colors given by match to
#'   the cluster names of \code{whichCluster} (see above), or otherwise given a
#'   default assignment.
#' @param ... Arguments passed to \code{\link{plotHeatmap}}
#' @details If the column 'ContrastName' is given in \code{signifTable}, these
#'   names will be used to describe the contrast in the legend.
#' @details Within each contrast, the genes are sorted by log fold-change if the
#'   column "logFC" is in the \code{signifTable} data.frame
#' @return A heatmap is created. The output of \code{plotHeatmap} is returned.
#' @seealso \code{\link{plotHeatmap}}, \code{\link{makeBlankData}},
#'   \code{\link{getBestFeatures}}
#' @export
#'
#' @examples
#' data(simData)
#'
#' cl <- clusterSingle(simData, subsample=FALSE,
#' sequential=FALSE, mainClusterArgs=list(clusterFunction="pam", clusterArgs=list(k=8)))
#'
#' #Do all pairwise, only return significant, try different adjustments:
#' pairsPerC <- getBestFeatures(cl, contrastType="Pairs", number=5,
#' p.value=0.05, isCount=FALSE)
#' plotContrastHeatmap(cl,pairsPerC)
setMethod(
  f = "plotContrastHeatmap",
  signature = "ClusterExperiment",
  definition = function(object,signifTable,whichCluster=NULL,contrastColors=NULL,...) {
    if(!all(c("IndexInOriginal","Contrast") %in% colnames(signifTable ))) stop("signifTable must have columns 'IndexInOriginal' and 'Contrast'")
    if(!is.numeric(signifTable$IndexInOriginal)) stop("Column 'IndexInOriginal' Must consist of numeric values")
    if(!all(signifTable$IndexInOriginal %in% seq_len(nrow(object)))) stop("Column 'IndexInOriginal' must consist of indices that match the row indices of 'object'")
    #divide by contrast, sort by FC (if exists) and return index
    geneByContrast<-by(signifTable,signifTable$Contrast,function(x){
      if("logFC" %in% names(x)){
        x<-x[order(x$logFC),]
      }
      x$IndexInOriginal
      
    })
    ##Check names of contrastColors
    if(!is.null(contrastColors)){
      if(!all(sort(names(contrastColors))==sort(unique(signifTable$Contrast)))){
        if("ContrastName" %in% colnames(signifTable)){
          if(!all(sort(names(contrastColors))==sort(unique(signifTable$ContrastName)))){
            warning("names of contrastColors do not match 'Contrast' or 'ContrastName' values; will be ignored.")
            contrastColors<-NULL
            
          }
          else contrastMatch<-FALSE		
        }
        else{
          warning("names of contrastColors do not match 'Contrast' values; will be ignored.")
          contrastColors<-NULL
        }
      }
      else contrastMatch<-TRUE
    }
    if(is.null(contrastColors)) contrastMatch<-FALSE
    
    if("ContrastName" %in% colnames(signifTable)){
      #give names to be contrast names
      gpnames<-unique(signifTable[,c("Contrast","ContrastName")])
      if(nrow(gpnames)==length(geneByContrast)){
        m<-match(names(geneByContrast),gpnames[,"Contrast"])
        names(geneByContrast)<-gpnames[m,"ContrastName"]
        #give colors the new names so will match
        if(!is.null(contrastColors) & contrastMatch){
          m<-match(names(contrastColors),gpnames[,"Contrast"])
          names(contrastColors)<-gpnames[m,"ContrastName"]
        }
      }
      
      #fix so order is order of nodes 
      nodeGrep<-grep("Node",names(geneByContrast))
      if(length(nodeGrep)==length(geneByContrast)){
        nsplit<-strsplit(names(geneByContrast),"Node")
        norder<-order(as.numeric(sapply(nsplit,.subset2,2)))
        geneByContrast<-geneByContrast[norder]
      }
    }
    
    
    if(!is.null(whichCluster)){
	    whichCluster<-.convertSingleWhichCluster(object,whichCluster,passedArgs=list(...))
      cl<-clusterMatrix(object)[,whichCluster]
      clMat<-clusterLegend(object)[[whichCluster]]
      clMat<-clMat[which(clMat[,"clusterIds"]>0),] #remove negatives
      pad<-if(length(unique(cl[cl>0]))<100) 2 else 3  
      mCl<-paste("Cl",stringr::str_pad(clMat[,"clusterIds"],width=pad,pad="0"),sep="") #mimic getBestFeatures
      for(ii in 1:length(mCl)){
        names(geneByContrast)<-gsub(mCl[[ii]],clMat[ii,"name"],names(geneByContrast))
      }
      ###If the contrast is one against all, should have name of cluster so give color
      if(is.null(contrastColors) & identical(sort(names(geneByContrast)),unname(sort(clMat[,"name"])))){
        contrastColors<-clMat[,"color"]
        names(contrastColors)<-names(geneByContrast)
      }
    }
    if(is.null(contrastColors)){
      #do it here so don't mess up above default assignment of colors
      contrastColors<-tail(massivePalette,length(geneByContrast)) #least likely to be important colors by accident
      names(contrastColors)<-names(geneByContrast)
    }
    
    plotHeatmap(object,clusterFeaturesData=geneByContrast,clusterLegend=list("Gene Group"=contrastColors),...)
    
  }
)