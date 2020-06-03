#' @name plotHeatmap
#' @title Heatmap for showing clustering results and more
#' @description Make heatmap with color scale from one matrix and hiearchical
#'   clustering of samples/features from another. Also built in functionality
#'   for showing the clusterings with the heatmap. Builds on
#'   \code{\link[pheatmap]{pheatmap}} function of \code{pheatmap} package.
#' @docType methods
#' @param colData If input to \code{data} is either a
#'   \code{\link{ClusterExperiment}},or \code{SummarizedExperiment} object or
#'   \code{SingleCellExperiment}, then \code{colData} must index the
#'   colData stored as a \code{DataFrame} in \code{colData} slot of the
#'   object. Whether that data is continuous or not will be determined by the
#'   properties of \code{colData} (no user input is needed). If input to
#'   \code{data} is matrix, \code{colData} is a matrix of additional data on
#'   the samples to show above heatmap. In this case, unless indicated by
#'   \code{whColDataCont}, \code{colData} will be converted into factors,
#'   even if numeric. ``-1'' indicates the sample was not assigned to a cluster
#'   and gets color `unassignedColor' and ``-2`` gets the color 'missingColor'.
#' @param data data to use to determine the heatmap. Can be a matrix,
#'   \code{\link{ClusterExperiment}},
#'   \code{\link[SingleCellExperiment]{SingleCellExperiment}} or
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}} object. The
#'   interpretation of parameters depends on the type of the input to
#'   \code{data}.
#' @param whColDataCont Which of the \code{colData} columns are continuous
#'   and should not be converted to counts. \code{NULL} indicates no additional
#'   \code{colData}. Only used if \code{data} input is matrix.
#' @param visualizeData either a character string, indicating what form of the
#'   data should be used for visualizing the data (i.e. for making the
#'   color-scale), or a data.frame/matrix with same number of samples as
#'   \code{assay(data)}. If a new data.frame/matrix, any character arguments to
#'   clusterFeaturesData will be ignored.
#' @param clusterSamplesData If \code{data} is a matrix,
#'   \code{clusterSamplesData} is either a matrix whose columns will be used by
#'   \code{hclust} to define the hiearchical clustering of samples (e.g.
#'   normalized data) or a pre-existing dendrogram (of class
#'  \code{\link[stats]{dendrogram}}) that clusters the samples. If
#'   \code{data} is a \code{ClusterExperiment} object, \code{clusterSamplesData}
#'   should be either character or integers or logical which indicates how (and
#'   whether) the samples should be clustered (or gives indices of the order for
#'   the samples). See details.
#' @inheritParams getClusterIndex
#' @param clusterFeaturesData  If \code{data} is a matrix, either a matrix whose
#'   rows will be used in \code{hclust} to define the hiearchical clustering of
#'   features (e.g. normalized data) or a pre-existing dendrogram that clusters
#'   the features. If \code{data} is a \code{ClusterExperiment} object, the
#'   input should be either character or integers indicating which features
#'   should be used (see details).
#' @param clusterSamples Logical as to whether to do hierarchical clustering of
#'   cells (if FALSE, any input to clusterSamplesData is ignored).
#' @param clusterFeatures Logical as to whether to do hiearchical clustering of
#'   features (if FALSE, any input to clusterFeaturesData is ignored).
#' @param showSampleNames Logical as to whether show sample names.
#' @param showFeatureNames Logical as to whether show feature names.
#' @param colorScale palette of colors for the color scale of the heatmap.
#' @param clusterLegend Assignment of colors to the clusters. If \code{NULL},
#'   \code{colData} columns will be assigned colors internally. See details
#'   for more.
#' @param alignColData Logical as to whether should align the colors of the
#'   \code{colData} (only if \code{clusterLegend} not given and
#'   \code{colData} is not \code{NULL}).
#' @param breaks Either a vector of breaks (should be equal to length 52), or a
#'   number between 0 and 1, indicating that the breaks should be equally spaced
#'   (based on the range in the data) upto the `breaks' quantile, see
#'   \code{\link{setBreaks}}
#' @param unassignedColor color assigned to cluster values of '-1'
#'   ("unassigned").
#' @param missingColor color assigned to cluster values of '-2' ("missing").
#' @param ... for signature \code{matrix}, arguments passed to \code{pheatmap}.
#'   For the other signatures, passed to the method for signature \code{matrix}.
#'   Not all arguments can be passed to \code{pheatmap} effectively, see
#'   details.
#' @param nFeatures integer indicating how many features should be used (if
#'   \code{clusterFeaturesData} is 'var' or 'PCA').
#' @param isSymmetric logical. if TRUE indicates that the input matrix is
#'   symmetric. Useful when plotting a co-clustering matrix or other sample by
#'   sample matrices (e.g., correlation).
#' @param overRideClusterLimit logical. Whether to override the internal limit
#'   that only allows 10 clusterings/annotations. If overridden, may result in
#'   incomprehensible errors from \code{aheatmap}. Only override this if you
#'   have a very large plotting device and want to see if \code{aheatmap} can
#'   render it.
#' @param plot logical indicating whether to plot the heatmap. Mainly useful for
#'   package mantaince to avoid calls to aheatmap on unit tests that take a long
#'   time.
#' @param symmetricBreaks logical as to whether the breaks created for the color
#'   scale should be symmetrical around 0
#' @param capBreaksLegend logical as to whether the legend for the breaks should
#'   be capped. Only relevant if \code{breaks} is a value < 1, in which case if
#'   \code{capBreaksLegend=TRUE}, only the values between the quantiles
#'   requested will show in the color scale legend.
#' @param whichAssay numeric or character specifying which assay to use. See
#'   \code{\link[SummarizedExperiment]{assay}} for details.
#' @inheritParams clusterSingle
#' @details The plotHeatmap function calls \code{\link[pheatmap]{pheatmap}} to draw
#'   the heatmap. The main points of \code{plotHeatmap} are to 1) allow for
#'   different matrix inputs, separating out the color scale visualization and
#'   the clustering of the samples/features. 2) to visualize the clusters and
#'   meta data with the heatmap. The intended use case is to allow the user to
#'   visualize the original count scale of the data (on the log-scale), but
#'   create the hierarchical clustering on another, more appropriate dataset for
#'   clustering, such as normalized data. Similarly, some of the palettes in the
#'   package were developed assuming that the visualization might be on
#'   unscaled/uncentered data, rather than the residual from the mean of the
#'   gene, and thus palettes need to take on a greater range of relevant values
#'   so as to show meaningful comparisons with genes on very different scales.
#' @details If \code{data} is a \code{ClusterExperiment} object,
#'   \code{visualizeData} indicates what kind of transformation should be done
#'   to \code{assay(data)} for calculating the color scale. The features will be
#'   clustered based on these data as well. A different data.frame or matrix can
#'   be given for the visualization. For example, if the
#'   \code{ClusterExperiment} object contains normalized data, but the user
#'   wishes that the color scale be based on the log-counts for easier
#'   interpretation, \code{visualizeData} could be set to be the
#'   \code{log2(counts + 1)}.
#' @details If \code{data} is a \code{ClusterExperiment} object,
#'   \code{clusterSamplesData} can be used to indicate the type of clustering
#'   for the samples. If equal to `dendrogramValue` the dendrogram stored in
#'   \code{data} will be used; if dendrogram is missing, a new one will be
#'   created based on the \code{primaryCluster} of data using
#'   \code{\link{makeDendrogram}}, assuming no errors are created (if errors are
#'   created, then \code{clusterSamplesData} will be set to "primaryCluster").
#'   If \code{clusterSamplesData} is equal to "hclust", then standard
#'   hierachical clustering of the transformed data will be used. If
#'   \code{clusterSamplesData} is equal to 'orderSamplesValue' no clustering of
#'   the samples will be done, and instead the samples will be ordered as in the
#'   slot \code{orderSamples} of \code{data}. If \code{clusterSamplesData} is
#'   equal to 'primaryCluster', again no clustering will be done, and instead
#'   the samples will be ordered based on grouping the samples to match the
#'   primaryCluster of \code{data}; however, if the primaryCluster of
#'   \code{data} is only one cluster or consists soley of -1/-2 values,
#'   \code{clusterSamplesData} will be set to "hclust". If
#'   \code{clusterSamplesData}  is not a character value,
#'   \code{clusterSamplesData} can be a integer valued vector giving the order
#'   of the samples.
#' @details If \code{data} is a matrix, then \code{colData} is a data.frame
#'   of annotation data to be plotted above the heatmap and
#'   \code{whColDataCont} gives the index of the column(s) of this dataset
#'   that should be consider continuous. Otherwise the annotation data for
#'   \code{colData} will be forced into a factor (which will be nonsensical
#'   for continous data). If \code{data} is a \code{ClusterExperiment} object,
#'   \code{colData} should refer to a index or column name of the
#'   \code{colData} slot of \code{data}. In this case \code{colData} will be
#'   added to any choices of clusterings chosen by the \code{whichClusters}
#'   argument (if any). If both clusterings and sample data are chosen, the
#'   clusterings will be shown closest to data (i.e. on bottom).
#' @details If \code{data} is a \code{ClusterExperiment} object,
#'   \code{clusterFeaturesData} is not a dataset, but instead indicates which
#'   features should be shown in the heatmap. In this case
#'   \code{clusterFeatures} can be one of the following: \itemize{ \item{"all"}{
#'   All rows/genes will be shown} \item{character giving dimensionality
#'   reduction}{Should match one of values saved in \code{reducedDims} slot or a
#'   builtin function in \code{listBuiltInReducedDims()}. \code{nFeatures} then
#'   gives the number of dimensions to show. The heatmap will then be of the
#'   dimension reduction vectors} \item{character giving filtering}{ Should
#'   match one of values saved in \code{filterStats} slot or a builtin function
#'   in \code{listBuiltInFilterStats()}. \code{nFeatures} gives the number of
#'   genes to keep after filtering.} \item{character giving gene/row names}{ }
#'   \item{vector of integers giving row indices}{ } \item{a list of indices or
#'   rownames}{This is used to indicate that the features should be grouped
#'   according to the elements of the list, with blank (white) space between
#'   them (see \code{\link{makeBlankData}} for more details). In this case, no
#'   clustering is done of the features.} }
#' @details If \code{breaks} is a numeric value between 0 and 1, then
#'   \code{breaks} is assumed to indicate the upper quantile (on the log scale)
#'   at which the heatmap color scale should stop. For example, if
#'   \code{breaks=0.9}, then the breaks will evenly spaced up until the 0.9
#'   upper quantile of \code{data}, and then all values after the 0.9 quantile
#'   will be absorbed by the upper-most color bin. This can help to reduce the
#'   visual impact of a few highly expressed genes (features). See
#'   \code{\link{setBreaks}} for more details.
#' @details Note that plotHeatmap calls \code{\link[pheatmap]{pheatmap}} under the
#'   hood. 
#' @details \code{clusterLegend} takes the place of argument \code{annColors}
#'   from \code{aheatmap} for giving colors to the annotation on the heatmap.
#'   \code{clusterLegend} should be list of length equal to
#'   \code{ncol(colData)} with names equal to the colnames of
#'   \code{colData}. Each element of the list should be a either the format
#'   requested by \code{\link[pheatmap]{pheatmap}} (a vector of colors with names
#'   corresponding to the levels of the column of \code{colData}), or should
#'   be format of the \code{clusterLegend} slot in a \code{ClusterExperiment}
#'   object. Color assignments to the rows/genes should also be passed via
#'   \code{clusterLegend} (assuming \code{annRow} is an argument passed to
#'   \code{...}). If \code{clusterFeaturesData} is a \emph{named} list
#'   describing groupings of genes then the colors for those groups can be given
#'   in \code{clusterLegend} under the name "Gene Group". 
#' @details Many arguments can be passed on to \code{pheatmap}, however, some
#'   are set internally by \code{plotHeatmap.} In particular, setting the values
#'   of \code{cluster_rows} or \code{cluster_cols} will cause errors.
#'   \code{color} in \code{pheatmap} is replaced by \code{colorScale} in
#'   \code{plotHeatmap.} The \code{annotation_col} to give annotation to the
#'   samples is replaced by the \code{colData} in \code{plotHeatmap}; moreover,
#'   the \code{annotation_colors} option in \code{pheatmap} will also be set
#'   internally to give more vibrant colors than the default in \code{pheatmap}
#'   (for \code{ClusterExperiment} objects, these values can also be set in the
#'   \code{clusterLegend} slot ). 
#' @details Other options given by the user in \code{...} should be passed on to
#'   \code{pheatmap}, though they have not been all tested. There are some
#'   arguments to \code{pheatmap} that we have internally changed the default
#'   values: \code{scale="none"},\code{na_col="white"},\code{border_col=NA}. If
#'   a user sets these values, they will be whatever the user sets, but if these
#'   arguments are not set, our defaults are used instead. All other arguments
#'   to \code{pheatmap} not discussed here are set to the default of
#'   \code{pheatmap}. Useful options of \code{pheatmap} that the user might want
#'   to adjust include \code{treeheight_col=0} or \col{treeheight_row=0} to
#'   suppress plotting of the dendrograms, \code{annotation_legend=FALSE} to
#'   suppress the legend of factors shown beside columns/rows, and
#'   \code{show_rownames=FALSE} or \code{show_colnames=FALSE} to suppress
#'   plotting of row/column labels. If you find annotation or labels are getting
#'   cropped, you may find adjusting \code{cellwidth} or \code{cellheight} is
#'   helpful for giving more space on the sides.
#'   
#' @return Returns (invisibly) a list with elements
#' \itemize{
#' \item{\code{heatmapOut}}{ The output from the final call of
#' \code{\link[pheatmap]{pheatmap}}. This can be useful if you want heatmaps in
#' a grid, as \code{mfrow}/\code{mfcol} will not work. See examples below for
#' how to do this.}
#' \item{\code{colData}}{ the annotation data.frame given to the argument
#' \code{annotation_col} in \code{pheatmap}.}
#' \item{\code{clusterLegend}}{ the (internally updated) \code{clusterLegend}
#' annotation colors.}
#' \item{\code{breaks}}{ The breaks used for \code{pheatmap}, after adjusting
#' for quantile and to match the number of colors.}
#' }
#' @author Elizabeth Purdom
#' @seealso \code{\link[pheatmap]{pheatmap}}, \code{\link{makeBlankData}},
#'   \code{\link{showHeatmapPalettes}}, \code{\link{makeDendrogram}},
#'   \code{\link[stats]{dendrogram}}
#' @export
#' @examples
#' data(simData)
#'
#' cl <- rep(1:3,each=100)
#' cl2 <- cl
#' changeAssign <- sample(1:length(cl), 80)
#' cl2[changeAssign] <- sample(cl[changeAssign])
#' ce <- ClusterExperiment(simCount, cl2, transformation=function(x){log2(x+1)})
#'
#' #simple, minimal, example. Show counts, but cluster on underlying means
#' plotHeatmap(ce)
#'
#' #assign cluster colors
#' clColors <- bigPalette[20:23]
#' names(clColors) <- 1:3
#' plotHeatmap(data=simCount, clusterSamplesData=simData,
#' colData=data.frame(cl), clusterLegend=list(clColors))
#'
#' #show two different clusters
#' anno <- data.frame(cluster1=cl, cluster2=cl2)
#' out <- plotHeatmap(simData, colData=anno)
#'
#' #return the values to see format for giving colors to the annotations
#' head(out$clusterLegend)
#'
#' #assign colors to the clusters based on plotClusters algorithm
#' plotHeatmap(simData, colData=anno, alignColData=TRUE)
#'
#' #assign colors manually
#' annoColors <- list(cluster1=c("black", "red", "green"),
#' cluster2=c("blue","purple","yellow"))
#'
#' plotHeatmap(simData, colData=anno, clusterLegend=annoColors)
#'
#' #give a continuous valued -- need to indicate columns
#' anno2 <- cbind(anno, Cont=c(rnorm(100, 0), rnorm(100, 2), rnorm(100, 3)))
#' plotHeatmap(simData, colData=anno2, whColDataCont=3)
#'
#' # The following compares changing breaks quantile on visual effect
#' # Also shows how to run multi-grid heatmaps by 
#' # saving output and use grid.arrange
#' \dontrun{
#' # Save each of the 4 plots, without plotting, with plot=FALSE
#' out100<-plotHeatmap(simData, colorScale=seqPal1, 
#'     breaks=1, plot=FALSE, main="Full set of breaks")
#' out99<-plotHeatmap(simData,colorScale=seqPal1, 
#'     breaks=.99, plot=FALSE, main="0.99 Quantile Upper Limit")
#' out95<-plotHeatmap(simData,colorScale=seqPal1, 
#'     breaks=.95, plot=FALSE, main="0.95 Quantile Upper Limit")
#' out90<-plotHeatmap(simData, colorScale=seqPal1,
#'     breaks=.90, plot=FALSE, main="0.90 Quantile Upper Limit")
#' # define layout matrix (where will put each plot based on number 1-4)
#' if (require(gridExtra, character.only = TRUE)){
#' library(gridExtra)
#' layoutMatrix=rbind(c(1,2),c(3,4))
#' # Plot them, using the pheatmap output that is returned (see pheatmap)
#' grid.arrange(list(out100$heatmapOut$gtable,
#'    out99$heatmapOut$gtable,
#'    out95$heatmapOut$gtable,
#'    out90$heatmapOut$gtable), 
#'    layout_matrix=layoutMatrix)
#' }
#' }
#'
#' @rdname plotHeatmap
#' @aliases plotHeatmap,SingleCellExperiment-method
#' @importFrom stats hclust dist
#' @importFrom pheatmap pheatmap
setMethod(
  f = "plotHeatmap",
  signature = signature(data = "SingleCellExperiment"),
  definition = function(data, isCount=FALSE,transFun=NULL,...
  ){
    #get transformation function
    transformation<-.makeTransFun(transFun=transFun,isCount=isCount)
    fakeCL<-sample(c(1,2),size=NCOL(data),replace=TRUE)
    fakeCE<-ClusterExperiment(data, fakeCL,transformation=transformation,
                              checkTransformAndAssay=FALSE)
    if("whichClusters" %in% names(list(...))) stop("cannot provide argument 'whichClusters' for input data not of class 'ClusterExperiment'")
    plotHeatmap(fakeCE,whichClusters="none",...)
  })

#' @rdname plotHeatmap
#' @export
setMethod(
  f = "plotHeatmap",
  signature = signature(data = "SummarizedExperiment"),
  definition = function(data, isCount=FALSE,transFun=NULL,...
  ){
    plotHeatmap(as(data,"SingleCellExperiment"),...)
  })


#' @rdname plotHeatmap
#' @export
setMethod(
  f = "plotHeatmap",
  signature = signature(data = "table"),
  definition = function(data,...
  ){
    plotHeatmap(unclass(data),...)
  })

#' @rdname plotHeatmap
#' @param nBlankLines Only applicable if input is \code{ClusterExperiment} object. Indicates the number of lines to put between groups of features if \code{clusterFeaturesData} gives groups of genes (see details and \code{\link{makeBlankData}}).
setMethod(
  f = "plotHeatmap",
  signature = signature(data = "ClusterExperiment"),
  definition = function(data,
                        clusterSamplesData=c("dendrogramValue", "hclust", "orderSamplesValue", "primaryCluster"),
                        clusterFeaturesData="var", nFeatures=NA,
                        visualizeData=c("transformed","centeredAndScaled","original"),
                        whichClusters= c("primary","workflow","all","none"),
                        colData=NULL,clusterFeatures=TRUE, nBlankLines=2,
                        colorScale, whichAssay=1,
                        ...
  ){

    .convertTry<-function(x,tryResult){if(!inherits(tryResult,"try-error")) return(tryResult) else return(x)}
    userList<-list(...)    
		checkIgnore<-.depricateArgument(passedArgs=userList,"colData","sampleData") #06/2018 added in BioC 3.8
		if(!is.null(checkIgnore)){
			userList<-checkIgnore$passedArgs
			colData<-checkIgnore$val
		}
    #########
    ##Determine visualization data and default colorScale based on that
    #########
    externalData<-FALSE
    visualizeData <- .convertTry(visualizeData,
                                 try(match.arg(visualizeData), silent=TRUE))
    if(is.character(visualizeData)){
      if(!visualizeData %in% c("transformed","centeredAndScaled","original")) stop("visualizeData value, '",visualizeData,"',is invalid option")
    }
    else{
      if(!is.data.frame(visualizeData) && !is.matrix(visualizeData)) stop("if visualizeData is not character, must be either data frame or matrix")
      externalData<-TRUE
      if(!ncol(visualizeData)==ncol(assay(data))) stop("if give separate visualizeData, must have same number of sample (columns) as assay(data)")
    }
    if(missing(colorScale)) {
      colorScale <- seqPal5
      if(is.character(visualizeData)) {
        if (visualizeData == "centeredAndScaled") {
          colorScale <- seqPal4
        }
      }
    }

    ####
    ##Transform data and determine which features to use
    ####
    # clusterFeaturesData <- .convertTry(clusterFeaturesData,
    #                                    try(match.arg(clusterFeaturesData),
    #                                        silent=TRUE))

    if(is.list(clusterFeaturesData)){
      groupFeatures<-clusterFeaturesData
      clusterFeaturesData<-unlist(clusterFeaturesData)
    }
    else groupFeatures<-NULL
    if(!externalData){
      if(!clusterFeatures && visualizeData=="original"){
        heatData<-assay(data, whichAssay)
      }
      else{
        if(length(clusterFeaturesData)==1 && isPossibleReducedDims(data,clusterFeaturesData)){
          ##### Dimensionality reduction ####
          if(!isReducedDims(data,clusterFeaturesData)){
            data<-makeReducedDims(data,reducedDims=clusterFeaturesData,maxDims=nFeatures,whichAssay=whichAssay)
          }
          heatData<-t(reducedDim(data,type=clusterFeaturesData))
        }
        else{
          #remaining options subset the data,
          #either by filtering or by user-specified genes
          #These operations will filter the input data object to only the relevant genes
          if(length(clusterFeaturesData)==1 && isPossibleFilterStats(data, clusterFeaturesData) ){
            ##### Filter ####
            if(!isFilterStats(data,clusterFeaturesData)){
              data<-makeFilterStats(data,filterStats=clusterFeaturesData,whichAssay=whichAssay)
            }
            if(is.na(nFeatures)) nFeatures<-min(NROW(data),500)
            data<-filterData(data,filterStats=clusterFeaturesData,percentile=nFeatures)
          }
          else{
            ### Other character values ####
            if(is.character(clusterFeaturesData)){#gene names
              if(length(clusterFeaturesData)==1 && clusterFeaturesData=="all")
                whRows<-seq_len(NROW(data))
              else{
                ### Give specific genes to use ####
                if(is.null(rownames(data))) stop("Cannot give feature names in clusterFeaturesData unless assay(data) has rownames")
                else{
                  whRows<-match(clusterFeaturesData,rownames(data))
                  if(all(is.na(whRows))) stop("None of the feature names in clusterFeaturesData match rownames(assay(data))")
                  if(any(is.na(whRows))){
                    warning("Not all of the feature names in clusterFeaturesData match rownames(assay(data))")
                    whRows<-na.omit(whRows)
                  }
                }
              }
            }
            else{
              ### Numeric row indices ####
              if(any(!clusterFeaturesData %in% seq_len(NROW(data)))) stop("invalid indices for clusterFeaturesData")
              whRows<-clusterFeaturesData
            }
            data<-data[whRows,]
          }
          heatData<-switch(visualizeData,
                           "original"=assay(data, whichAssay),
                           "transformed"=transformData(data, whichAssay),
                           "centeredAndScaled"=t(scale(t(transformData(data, whichAssay)), center=TRUE, scale=TRUE))
          )
        }
      }
    }
    else{
      heatData<-visualizeData
    }


    ######
    #Make colData based on clusterings and columns of colData
    ######
    #---
    #Get clusterings
    #---
    whCl<-getClusterIndex(data,whichClusters=whichClusters,noMatch="silentlyRemove")
    #
    if(length(whCl)>0){
      clusterData<-clusterMatrixNamed(data,whichClusters=whCl)
    }
    else{
      if(any( whichClusters!="none")) warning("given whichClusters value does not match any clusters, none will be plotted")
      clusterData<-NULL
    }
    #---
    #get colData values and subset to those asked for by user
    #---
    sData<-.pullColData(data,colData)
    #identify which numeric
    if(!is.null(sData)) whCont<-which(sapply(seq_len(ncol(sData)),function(ii){is.numeric(sData[,ii])}))
    whColDataCont<-NULL

    if(!is.null(clusterData) & !is.null(sData)){
      colData<-data.frame(clusterData,sData,stringsAsFactors=FALSE,check.names=FALSE)
      if(length(whCont)>0)  whColDataCont<-whCont+ncol(clusterData)
    }
    else{
      if(!is.null(clusterData)) colData<-clusterData
      if(!is.null(sData)){
        colData<-sData
        if(length(whCont)>0) whColDataCont<-whCont
      }
      if(is.null(sData) & is.null(clusterData)) colData<-NULL
    }

    #------
    #check user didn't give something different for colors
    #------
    clLegend<-clusterLegend(data)[whCl] #note, clusterLegend gives names even though not stored internally with @clusterLegend so will match, which plotHeatmap needs
    if(length(clLegend)==0) clLegend<-NULL

    if(!"symmetricBreaks" %in% names(userList) && !externalData && visualizeData %in% c("centeredAndScaled","centered")){
      userList$symmetricBreaks<-TRUE
    }
    userAlign<-"alignColData" %in% names(userList) & !is.null(userList$alignColData)
    userLegend<-"clusterLegend" %in% names(userList) & !is.null(userList$clusterLegend)
    if(userAlign | userLegend){ #if user asks for alignment, don't assign clusterLegend
      if(userLegend){
        userClLegend<-userList[["clusterLegend"]]
        annotNames<-colnames(colData)
        if("annRow" %in% names(userList)) annotNames<-c(annotNames,colnames(userList$annRow))
        if(!is.null(groupFeatures)) annotNames<-c(annotNames,"Gene Group")
        userClLegend<-userClLegend[names(userClLegend) %in% annotNames]
        if(length(userClLegend)==0){
          warning("names of list given by user in clusterLegend do not match clusters nor colData chosen. Will be ignored.")
        }
        else{
          #keep existing clLegend from ClusterExperiment object if not conflict with user input:
          whNotShared<-which(!names(clLegend)%in% names(userClLegend))
          if(length(whNotShared)>0) clLegend<-c(userClLegend,clLegend[whNotShared]) else clLegend<-userClLegend
          clLegend<-.convertToAheatmap(clLegend, names=TRUE)

        }
        userList<-userList[-grep("clusterLegend",names(userList))]
      }
      else{
        if(userAlign){
          al<-userList[["alignColData"]]
          if(al) clLegend<-NULL
        }
      }
    }
    else clLegend<-.convertToAheatmap(clLegend,names=TRUE)
    ######
    #Create clusterSamplesData
    ######
    clusterSamplesData<-.convertTry(clusterSamplesData, try(match.arg(clusterSamplesData), silent=TRUE))
    if(is.logical(clusterSamplesData)) clusterSamples<-clusterSamplesData
    else{
      clusterSamples<-TRUE
      if(is.numeric(clusterSamplesData)){
        heatData<-heatData[,clusterSamplesData,drop=FALSE]
        if(!is.null(colData)) colData<-colData[clusterSamplesData,,drop=FALSE]
        clusterSamplesData<-heatData
        clusterSamples<-FALSE
      }
      else if(is.character(clusterSamplesData)){
        if(clusterSamplesData=="orderSamplesValue"){
          heatData<-heatData[,orderSamples(data),drop=FALSE]
          if(!is.null(colData)) colData<-colData[orderSamples(data), ,drop=FALSE]
          clusterSamplesData<-heatData
          clusterSamples<-FALSE
        }
        else if(clusterSamplesData=="dendrogramValue"){
					if(is.null(data@dendro_samples)){
			      clusterSamplesData <- try( convertToDendrogram(makeDendrogram(data)) ,silent = TRUE) 
	          if(inherits(clusterSamplesData, "try-error")){
	            warning("cannot make dendrogram from 'data' with default makeDendrogram options. Ordering by primary cluster without dendrogram")
	            clusterSamplesData<-"primaryCluster"
	          }
          }
          else{
						#make sure get the sample ids as labels of the tips:
            clusterSamplesData<-try(convertToDendrogram(data),silent=TRUE)
						if(inherits(clusterSamplesData, "try-error")){
	            warning("cannot make dendrogram class from stored dendrograms. Ordering by primary cluster without dendrogram")
	            clusterSamplesData<-"primaryCluster"
	          }
          }
        }
        if(is.character(clusterSamplesData) && clusterSamplesData=="primaryCluster"){
          whUnassign<-which(primaryCluster(data) %in% c(-1,-2))
          if(length(whUnassign)==nSamples(data) || length(unique(primaryCluster(data)[-whUnassign]))==1){
            #in this case, all -1/-2 or same cluster, just do heatmap with hclust
            warning("Cannot order by primary cluster because all one cluster and/or all clustering values are -1/-2. Using standard hiearchical clustering.")
            clusterSamplesData<-"hclust"
          }
          else{
            heatData<-heatData[,order(primaryCluster(data))]
            if(!is.null(colData)) colData<-colData[order(primaryCluster(data)),,drop=FALSE]
            clusterSamplesData<-heatData
            clusterSamples<-FALSE
          }
        }
        if(is.character(clusterSamplesData) && clusterSamplesData=="hclust"){
          #if hclust, then use the visualizeData data
          #unless visualizeData data is original, in which case use transformed (and possibly filtered)
          if(is.character(visualizeData) && visualizeData=="original")
            clusterSamplesData<-transformData(data,whichAssay=whichAssay)
          else clusterSamplesData<-heatData
        }
      }
      else stop("clusterSamplesData must be either character, or vector of indices of samples")
    }

    #################
    ###Deal with grouping of genes
    #################
    if(!is.null(groupFeatures)){
      #convert groupFeatures to indices on new set of data.
      groupFeatures<-lapply(groupFeatures,function(x){match(x,whRows)})
      blankData<-makeBlankData(data=heatData,groupsOfFeatures=groupFeatures,nBlankFeatures=nBlankLines)
      #replace heatData with one with blanks -- won't cluster them now...
      heatData<-data.matrix(blankData$dataWBlanks)
      labRow<-blankData$rowNamesWBlanks
      clusterFeatures<-FALSE

      if(!is.null(names(groupFeatures))){
        annRow<-list("Gene Group"=factor(blankData$featureGroupNamesWBlanks,levels=names(groupFeatures)))
        #show color-coding of gene groupings:
        if("annRow" %in% names(userList)){
          stop("Cannot provide 'annRow' if grouping features")
          # if(is.list(userList$annRow)) userList$annRow<-c(userList$annRow,annRow)
          # else userList$annRow<-c(list(userList$annRow),annRow)
        }
        else userList$annRow<-annRow
        if(!"Gene Group" %in% names(clLegend)){
          nGroups<-length(groupFeatures)
          groupColors<-bigPalette[seq_len(nGroups)]
          names(groupColors)<-levels(annRow[["Gene Group"]])
          clLegend<-c(clLegend, "Gene Group"=list(groupColors))
        }
      }
    }
    else{
      labRow<-rownames(heatData)
    }
    do.call("plotHeatmap",
        c(list(data=heatData,
            clusterSamplesData=clusterSamplesData,
            clusterFeaturesData=heatData, #set it so user doesn't try to pass it and have something weird happen because dimensions wrong, etc.
            colData=colData,whColDataCont=whColDataCont,
            clusterSamples=clusterSamples,labRow=labRow,
            clusterLegend=clLegend,clusterFeatures=clusterFeatures,
            colorScale=colorScale),userList))


  })


#' @rdname plotHeatmap
setMethod(
  f = "plotHeatmap",
  signature = signature(data = "data.frame"),
  definition = function(data,...){
	  plotHeatmap(data.matrix(data),...)
  }
)
#' @rdname plotHeatmap
setMethod(
  f = "plotHeatmap",
  signature = signature(data = "ExpressionSet"),
  definition = function(data,...){
	  plotHeatmap(exprs(data),...)
  }
)

#' @rdname plotHeatmap
setMethod(
  f = "plotHeatmap",
  signature = signature(data = "matrixOrHDF5"),
  definition = function(data,colData=NULL,
                        clusterSamplesData=NULL,
                        clusterFeaturesData=NULL,
                        whColDataCont=NULL,
                        clusterSamples=TRUE,showSampleNames=FALSE,
                        clusterFeatures=TRUE,showFeatureNames=FALSE,
                        colorScale=seqPal5,
                        clusterLegend=NULL,alignColData=FALSE,
                        unassignedColor="white",missingColor="grey",
                        breaks=NA,symmetricBreaks=FALSE,
                        capBreaksLegend=FALSE,
                        isSymmetric=FALSE, overRideClusterLimit=FALSE,
						plot=TRUE,...
  ){
    
    ##########
    ##Deal with numeric matrix for heatmap ...
    ##########
    heatData<-data.matrix(data)

    ##########
    ##Deal with passed arguments to pheatmap ...
    ##########
    pHeatmapArgs<-list(...)  
    # old version was sampleData -> colData
    checkIgnore<-.depricateArgument(passedArgs=pHeatmapArgs,
            newArgName="colData",oldArgName="sampleData") #06/2018 added in BioC 3.8
	if(!is.null(checkIgnore)){
		pHeatmapArgs<-checkIgnore$passedArgs
		colData<-checkIgnore$val
	}
    
    # got rid of labelTracks argument for plotHeatmap, 
    # because now part of pheatmap
    checkIgnore<-.depricateArgument(passedArgs=pHeatmapArgs,
            oldArgName="labelTracks") #05/2020 in development
	if(!is.null(checkIgnore)){
		pHeatmapArgs<-checkIgnore$passedArgs
        pHeatmapArgs<-c(pHeatmapArgs,list(labelTracks=checkIgnore$val))
	}


    # This is function for assigning a value (with my choice being given in argument `value`) to a pheatmap option, but making sure that I don't override a user values.
    # If value=NULL, then just uses pheatmap default (from pHeatmapDefaultArgs)
    pHeatmapDefaultArgs<-as.list(args(pheatmap::pheatmap))
    getHeatmapValue<- function(string,value=NULL){ #note, doesn't work for pulling function 'reorder' so put in manually
      if(string %in% names(pHeatmapArgs)) val<-pHeatmapArgs[[string]]
      else{
        if(is.null(value)) val<-pHeatmapDefaultArgs[[string]]
        else val<-value
      }
      return(val)
    }
    badValues<-c("cluster_rows","cluster_cols","color",
        "annotation_col","annotation_colors")
    replacedValues<-c("clusterSamples","clusterFeatures",
        "colorScale","colData","clusterLegend")
    if(any(badValues %in% names(pHeatmapArgs))) 
        stop("The following arguments to aheatmap cannot be set by the user in plotHeatmap:",paste(badValues,collapse=","),". They are over-ridden by: ",paste(replacedValues,collapse=","))

    ##########
    ### Create the object passed to pheatmap for the clustering:
    ###
    ### RowV/ColV is what is actually passed to pheatmap for the clustering
    ### passed to argument `cluster_rows`/`cluster_columns` (boolean values determining if rows/columns should be clustered or hclust object)
    ##########
    
    #---- samples (Colv):
    if(clusterSamples){ 
        dendroSamples<-.convertDataToDendro(clusterSamplesData,
            N=ncol(heatData),
            dimDirection="columns")
        if(is.null(dendroSamples)) Colv<-TRUE #then just pass the data in heatData
        else Colv<-dendroSamples
    }
    else Colv<-clusterSamples
    
    #---- features (Rowv):
    if(isSymmetric){
      Rowv<-Colv
    }
    else{
      if(clusterFeatures){
          dendroFeatures<-.convertDataToDendro(clusterFeaturesData,
              N=nrow(heatData),
              dimDirection="rows")
          if(is.null(dendroFeatures)) Rowv<-TRUE #then just pass the data in heatData
          else Rowv<-dendroFeatures
      }
      else Rowv<-clusterFeatures
    }
    
    
    ##########
    ##Deal with annotation of samples (colData) ...
    ##########
    if(!is.null(colData)){
        if(!is.matrix(colData) & !is.data.frame(colData)) 
            stop("colData must be a either a matrix or a data.frame")
        if(NCOL(heatData) != NROW(colData)) 
            stop("colData must have same number of rows as columns of heatData")
        colDataOut<-.fixColData(annoteData=colData,
            clusterLegend=clusterLegend,
            whColumnsCont=whColDataCont,
            unassignedColor=unassignedColor,
            overRideClusterLimit=overRideClusterLimit, 
            missingColor=missingColor, 
            alignClusters=alignColData)
        annCol<-colDataOut$annCol
        annColors<-colDataOut$annColors
        clusterLegend<-colDataOut$clusterLegend
    }
    else{ #no colData provided -- just a heatmap with no annotation
      annCol<-NA
      annColors<-NA
    }
    
    #############
    # put into pheatmap
    #############
    
    # Could do this more smoothly now with pheatmap 
    # which has argument legend_breaks...
    # But for now keeping it. 
    capBreaks<-length(breaks)==1 & capBreaksLegend
    breaks<-setBreaks(data=heatData, breaks=breaks, makeSymmetric=symmetricBreaks,returnBreaks=!capBreaks,numberOfColors=length(colorScale))
    if(capBreaks){ #so the legend is not so weird
      if(length(breaks)!=2)
        stop("coding error in new breaks function")
      heatData[which(heatData<breaks[1])]<-breaks[1]
      heatData[which(heatData>breaks[2])]<-breaks[2]
      breaks<-seq(breaks[1],breaks[2],length=length(colorScale)+1)
    }
    #have to have matching names for the annotation (quirk of pheatmap)!
    colnames(heatData)<-rownames(annCol) 
    out<-do.call(pheatmap::pheatmap,c(list(mat=heatData,
        cluster_rows=Rowv,
        cluster_cols=Colv,
        color=colorScale,
        breaks=breaks,
        scale=getHeatmapValue("scale","none"), 
        na_col=getHeatmapValue("na_col","white"),
        border_color=getHeatmapValue("border_color",NA),
        annotation_col=annCol, silent=!plot,
        annotation_colors=annColors),pHeatmapArgs))   

    
    invisible(list(heatmapOut=out,
        colData=annCol,
        clusterLegend=clusterLegend,
        breaks=breaks))
  }
)

#' @param object either hclust, dendrogram, or matrix-like object that will serve to cluster the data shown in heatmap. If matrix, then will assume the ROWS of object are the values to be clustered
#' @param N the size of the dimension in the data shown in heatmap (i.e. not used for clustering, but for visualizing)
#' @param dimDirection one of "rows" or "columns" indicating whether this data will be guiding the clustering of the rows or columns of the data visualized in heatmap (mainly so gives interpretable error messages)
.convertDataToDendro<-function(object, N, dimDirection)
{
    dimDirection<-match.arg(dimDirection,c("rows","columns"))
    nameOfObject<-switch(dimDirection,
        "rows"="clusterFeaturesData",
        "columns"="clusterSamplesData")
    if(inherits(object,"hclust")){
        if(object$merge != (N-1)) stop(paste(nameOfObject,"hclust object is not on same number of", dimDirection," as heatData"))
        return(object)
    }
	if(inherits(object,"dendrogram")){
        if(nobs(object)!=N) stop(paste(nameOfObject,"dendrogram is not on same number of ", dimDirection," as heatData"))
        return(as.hclust(object))
        
    }
    if(is.null(object)) return(NULL)
    if(!is.data.frame(object) & !is.matrix(object) & !inherits(object,"DelayedArray")) stop(nameOfObject," must either be dendrogram, or data.frame/matrix/DelayedArray class")

    #check valid (do this before make it matrix to save on computation)
    if(dimDirection=="rows"){
        if(nrow(object)!=N) stop(nameOfObject, "matrix does not have same number of ", dimDirection," as heatData")
        object<-data.matrix(object)        
    }
    else{
        if(ncol(object)!=N) stop(nameOfObject, "matrix does not have same number of ", dimDirection," as heatData")
        object<-t(data.matrix(object))
        
    }
    ## FIXME!!! Need to figure out exactly what pheatmap would do here. For now, reverted to hclust....
        # dendroSamples<-NMF:::cluster_mat(t(object),param=TRUE,distfun=getHeatmapValue("distfun"),hclustfun=getHeatmapValue("hclustfun"),reorderfun=getHeatmapValue("reorderfun",value=function(d, w) reorder(d, w)))$dendrogram
    return(stats::hclust(stats::dist(object))) #dist finds distances between rows
}


#' @param object either hclust, dendrogram, or matrix-like object that will serve to cluster the data shown in heatmap. If matrix, then will assume the ROWS of object are the values to be clustered
#' @param N the size of the dimension in the data shown in heatmap (i.e. not used for clustering, but for visualizing)
#' @param dimDirection one of "rows" or "columns" indicating whether this data will be guiding the clustering of the rows or columns of the data visualized in heatmap (mainly so gives interpretable error messages)
#' @param whColumnsCont which columns of annoteData should be considered continuous, either vector of integers, or logical valued vector. Otherwise, will be assumed all are factors...
#' @details The point of this function is to convert the annotation data into factors where continuous. This allows integer-valued cluster assignments to not be treated as continuous values. It does this by running .makeColors. This will:
#-------------------
# 1) Make annoteData explicitly factors (except for whColumnsCont variables which are not given to function)
# 2) Make a default clusterLegend, including incorporating and checking user-given clusterLegend
# 3) Make a numeric summary of factors (for if alignSamples==TRUE)
#-------------------

.fixColData<-function(annoteData,
    clusterLegend,whColumnsCont,
    unassignedColor,overRideClusterLimit, 
    missingColor, alignClusters){

    #----
    #check annoteData input:
    #----
    if(NCOL(annoteData)>10){
      # is this a problem with pheatmap??? Not too important perhaps to fix, but...
      if(overRideClusterLimit) warning("More than 10 annotations/clusterings may result in errors in pheatmap. You have >10 but have chosen to override the internal stop by setting overRideClusterLimit=TRUE.")
      else stop("More than 10 annotations/clusterings may not be reliably shown in plotHeatmap. To override this limitation and try for yourself, set overRideClusterLimit=TRUE.")
    }
    
    if(is.data.frame(annoteData)) annoteData<-droplevels(annoteData)
    if(!is.null(whColumnsCont)){
      if(any(logical(whColumnsCont))) whColumnsCont<-which(whColumnsCont)
    }
    if(length(whColumnsCont)>0){
        annCol<-annoteData #want this here so get continuous values that are getting dropped from annoteData
        annoteData<-annoteData[,-whColumnsCont,drop=FALSE]
    } 
    #----
    # Run .makeColors
    #----
    defaultColorLegend<-.makeColors(annoteData,
            colors=massivePalette,
            unassignedColor=unassignedColor,
            missingColor=missingColor, 
            distinctColors=TRUE, 
            matchClusterLegend = clusterLegend, 
            matchTo="name") 
    tmpDfNum<-defaultColorLegend$numClusters
    colnames(tmpDfNum)<-colnames(annoteData)
    #replace annoteData with results of .makeColors and update/create annCol
    annoteData<-defaultColorLegend$facClusters
    if(length(whColumnsCont)>0){
      annCol[,-whColumnsCont]<-annoteData
    }
    else annCol<-annoteData
        
    #-----
    # if aligning cluster colors (alignClusters=TRUE), 
    # then give values and create clusterLegend based on aligned
    # using plotClusters
    #-----
    if(alignClusters & is.null(clusterLegend) &  
        (is.null(whColumnsCont) || length(whColumnsCont)<ncol(annCol))){
          #align the clusters and give them colors
          alignObj<-plotClusters(tmpDfNum ,plot=FALSE,
              unassignedColor=unassignedColor, 
              missingColor=missingColor)
          defaultColorLegend<-.makeColors(annoteData,
              clNumMat=tmpDfNum, colors=massivePalette,
              unassignedColor=unassignedColor,
              missingColor=missingColor, 
              matchClusterLegend=alignObj$clusterLegend,
              matchTo="numIds")
    }
    #preserve those in given clusterLegend that don't match annoteData (could go with features/rows)
		
    if(is.list(clusterLegend)){ #could be single vector, but in that case, will loose them
      whKeep<-names(clusterLegend)[which(!names(clusterLegend) %in%  
          names(defaultColorLegend$colorList))]
      clusterLegend<-c(defaultColorLegend$colorList, clusterLegend[whKeep])
    }
    else clusterLegend<-defaultColorLegend$colorList

    # Convert to aheatmap format, if needed
    annColors<-.convertToAheatmap(clusterLegend, names=TRUE)
    
    ##########################
    # remove any unused level colors to clean up legend,
    # give names to colors if not given, and
    # make them in same order as in annCol factor
    ##########################
    whInAnnColors<-which(names(annColors) %in% colnames(annCol))
    if(!is.null(whColumnsCont) & length(whColumnsCont)>0){
      whInAnnColors<-setdiff(whInAnnColors,whColumnsCont)
    }
    
    #-----
    # Fix up annotation/clusterLegend for the heatmap
    # Current problems with pheatmap:
    #   If extra colors in colorLegend, will plot them regardless of drop_levels. Have submitted bug report:
    #   https://github.com/raivokolde/pheatmap/issues/69
    # This code gets rid of extra levels. Unlike aheatmap, doesn't seem to be a problem to have a single level in factor
    # Also problem if not have enough colors, so regardless want to leave in that check...
    #-----
    prunedList<-lapply(whInAnnColors,function(ii){
      nam<-names(annColors)[[ii]] #the name of variable
      x<-annColors[[ii]] ##list of colors
      levs<-levels(annCol[,nam])
      if(length(x)<length(levs)) stop("number of colors given for ",nam," is less than the number of levels in the data")
      if(is.null(names(x))){
        #if user didn't give names to colors, assign them in order of levels
        x<-x[1:length(levs)] #shorten if needed
        names(x)<-levs
      }
      else{
        if(any(!levs %in% names(x))) stop("colors given for ",nam," do not cover all levels in the data")
        whlevs<-match(levs,names(x))
        x<-x[whlevs]
      }
      return(x)
    })
    names(prunedList)<-names(annColors)[whInAnnColors]
    annColors[whInAnnColors]<-prunedList
    
    return(list(annCol=annCol,
        annColors=annColors,
        clusterLegend=clusterLegend))
}


#' @rdname plotHeatmap
#' @aliases plotCoClustering
#' @param invert logical determining whether the coClustering matrix should be
#'   inverted to be 1-coClustering for plotting. By default, if the diagonal
#'   elements are all zero, invert=TRUE, and otherwise invert=FALSE. If
#'   coClustering matrix is not a 0-1 matrix (e.g. if equal to a distance matrix
#'   output from \code{\link{clusterSingle}}, then the user should manually set
#'   this parameter to FALSE.)
#' @param saveDistance logical. When the \code{coClustering} slot contains
#'   indices of the clusterings or a NxB set of clusterings, the hamming
#'   distance will be calculated before running the plot. This argument
#'   determines whether the \code{ClusterExperiment} object with that distance
#'   in \code{coClustering} slot should be returned (so as to avoid
#'   re-calculating it in the future) or not.
#' @details \code{plotCoClustering} is a convenience function to plot the
#'   heatmap of the co-clustering distance matrix from the \code{coClustering}
#'   slot of a \code{ClusterExperiment} object (either by calculating the
#'   hamming distance of the clusterings stored in the \code{coClustering} slot,
#'   or the distance stored in the \code{coClustering} slot if it has already
#'   been calculated.
#' @export
setMethod(
  f = "plotCoClustering",
  signature = "ClusterExperiment",
  definition = function(data, invert, saveDistance=FALSE,...){
    typeCoCl<-.typeOfCoClustering(data)
    if(typeCoCl=="null") stop("coClustering slot is empty")
    if(typeCoCl=="indices"){
        #calculate the distance
        coClustering(data)<-.clustersHammingDistance(
            t(clusterMatrix(data,whichClusters=data@coClustering)))        
    }
    else{
        #Need to catch if not a symmetric matrix
        #(case of subsampling matrix)
        if(typeCoCl=="clusterings"){
            #calculate the distance
            coClustering(data)<-.clustersHammingDistance(
                t(coClustering(data)))        
            
        }
    }
    if(missing(invert)) 
        invert<-ifelse(all(diag(data@coClustering)==0),TRUE,FALSE)
    if(invert){
        #Make it a similarity matrix (better anyway for the sparse representation)
        coClustering(data) <- 
            as(as(1-data@coClustering,"sparseMatrix"),"symmetricMatrix")
    }
    	
    # Do all this so don't have to erase merge info from data (so can return calculated distance to user)
    fakeCE<-ClusterExperiment(as(data@coClustering,"matrix"),
                              clusterMatrix(data),
                              transformation=function(x){x},
                              checkTransformAndAssay=FALSE


    )
    for(sName in c('clusterMatrix', 'primaryIndex', 'clusterInfo', 
        'clusterTypes', 'dendro_samples', 'dendro_clusters', 'dendro_index', 
        'clusterLegend', 'orderSamples', 'merge_index', 
        'merge_dendrocluster_index', 'merge_method', 'merge_demethod', 
        'merge_cutoff', 'merge_nodeProp', 'merge_nodeMerge','colData')){
        slot(fakeCE, sName)<-slot(data,sName)
    }   
    plotHeatmap(fakeCE,isSymmetric=TRUE,clusterFeaturesData="all",...)
    if(saveDistance) return(data)
    else invisible()
  })
