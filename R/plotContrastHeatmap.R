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
#'   significance table. Used to match to colors in \code{clusterLegend(object)} 
#'   (relevant for one-vs-all contrast so that color aligns). See description of 
#'   argument in \code{\link{getClusterIndex}} for futher details.
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
#' @details Note that if \code{whichCluster} is NOT given (the default) then there is 
#' no automatic match of colors with contrasts based on the information in 
#' \code{object}.
#' @return A heatmap is created. The output of \code{plotHeatmap} is returned.
#' @seealso \code{\link{plotHeatmap}}, \code{\link{makeBlankData}},
#'   \code{\link{getBestFeatures}}
#' @export
#'
#' @examples
#' data(simData)
#'
#' cl <- clusterSingle(simData, subsample=FALSE,
#' sequential=FALSE, 
#' mainClusterArgs=list(clusterFunction="pam", clusterArgs=list(k=8)))
#'
#' #Do all pairwise, only return significant, try different adjustments:
#' pairsPerC <- getBestFeatures(cl, contrastType="Pairs", number=5,
#' p.value=0.05, DEMethod="limma")
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
      if(!all(sort(names(contrastColors)) == sort(unique(signifTable$Contrast))) ){
        if("ContrastName" %in% colnames(signifTable)){
          if(!all(sort(names(contrastColors)) == sort(unique(signifTable$ContrastName))) ){
            warning("names of contrastColors do not match 'Contrast' or 'ContrastName' values; will be ignored.")
            contrastColors<-NULL
            
          }
          else contrastMatch<-FALSE	#FALSE because don't have to fix the names to match ContrastName, since already do...	
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
	  internalNames<-names(geneByContrast) #incase need the Contrast name later in code
      gpnames<-unique(signifTable[,c("Contrast", "ContrastName", "InternalName")])
      if(nrow(gpnames)==length(geneByContrast)){
        m<-match(names(geneByContrast),gpnames[,"Contrast"])
        names(geneByContrast)<-gpnames[m,"ContrastName"]
		internalNames<-gpnames[m,"InternalName"]
		#give colors the new names so will match
        if(!is.null(contrastColors) & contrastMatch){
          m<-match(names(contrastColors),gpnames[,"Contrast"])
          names(contrastColors)<-gpnames[m,"ContrastName"]
        }
      }
      
      #fix so order is order of contrast names
	  order<-order(names(geneByContrast))
	  geneByContrast<-geneByContrast[order]

    }
    
    
    if(!is.null(whichCluster)){
      ## Assign colors to the contrasts 
      ## (only if contrast is oneAgainstAll)
	  	whichCluster<-getSingleClusterIndex(object,
          whichCluster,passedArgs=list(...))
      cl<-clusterMatrix(object)[,whichCluster]
      clMat<-clusterLegend(object)[[whichCluster]]
      clMat<-clMat[which(clMat[,"clusterIds"]>0),] #remove negatives
      ### If the contrast is one against all, should have name of cluster so give color
      ### Note that need to use InternalId, in case the names are not unique?
      internalNames<-gsub("Cl","",internalNames)
			### This is messy. Basically as.numeric will give NAs 
      # This can give warnings for NAs, and if option(warn=2) create error...so have suppressWarnings
      internalNames<-try(sort(suppressWarnings(as.numeric(internalNames))),silent=TRUE)
      if(inherits(internalNames, "try-error") ||
          any(is.na(internalNames))) isOneVAll<-FALSE ## I'm not clear that anything gets here...
      else isOneVAll<-TRUE
      clIds<-unname(sort(as.numeric(clMat[,"clusterIds"])))
      if(is.null(contrastColors) && isOneVAll &&
          identical(internalNames, clIds)){
              contrastColors<-clMat[,"color"]
              names(contrastColors)<-names(geneByContrast)
      }
    }
    if(is.null(contrastColors)){
      #do it here so don't mess up above default assignment of colors
      contrastColors<-tail(massivePalette,length(geneByContrast)) #least likely to be important colors by accident
      names(contrastColors)<-names(geneByContrast)
    }
    
    plotHeatmap(object,
        clusterFeaturesData=geneByContrast,
        clusterLegend=list("Gene Group"=contrastColors),...)
    
  }
)
