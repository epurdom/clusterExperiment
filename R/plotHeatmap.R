#' Heatmap for showing clustering results and more
#'
#' Make heatmap with color scale from one matrix and hiearchical clustering of
#' samples/features from another. Also built in functionality for showing the
#' clusterings with the heatmap. Builds on \code{\link[NMF]{aheatmap}} function
#' of \code{NMF} package.
#'
#' @docType methods
#' @param sampleData If input is either a \code{\link{ClusterExperiment}} or
#'   \code{SummarizedExperiment} object, then \code{sampleData} must index the
#'   sampleData stored as a \code{DataFrame} in \code{colData} slot of the
#'   object. Whether that data is continuous or not will be determined by the
#'   properties of \code{colData} (no user input is needed). If input is matrix,
#'   \code{sampleData} is a matrix of additional data on the samples to show
#'   above heatmap. Unless indicated by \code{whSampleDataCont},
#'   \code{sampleData} will be converted into factors, even if numeric. ``-1''
#'   indicates the sample was not assigned to a cluster and gets color
#'   `unassignedColor' and ``-2`` gets the color 'missingColor'.
#' @param data data to use to determine the heatmap. Can be a matrix,
#'   \code{\link{ClusterExperiment}} or
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}} object. The
#'   interpretation of parameters depends on the type of the input.
#' @param whSampleDataCont Which of the \code{sampleData} columns are continuous
#'   and should not be converted to counts. \code{NULL} indicates no additional
#'   \code{sampleData}.
#' @param visualizeData either a character string, indicating what form of the
#'   data should be used for visualizing the data (i.e. for making the
#'   color-scale), or a data.frame/matrix with same dimensions of
#'   \code{assay(data)}.
#' @param clusterSamplesData If \code{data} is a matrix, either a matrix that 
#'   will be used to in \code{hclust} to define the hiearchical clustering of 
#'   samples (e.g. normalized data) or a pre-existing dendrogram that clusters 
#'   the samples. If \code{data} is a \code{ClusterExperiment} object, the input
#'   should be either character or integers or logical. Indicates how (and
#'   whether) the samples should be clustered (or gives indices of the order for
#'   the samples). See details.
#' @param whichClusters character string, or vector of characters or integers,
#'   indicating what clusters should be visualized with the heatmap.
#' @param clusterFeaturesData  If \code{data} is a matrix, either a matrix that
#'   will be used in \code{hclust} to define the hiearchical clustering of
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
#'   \code{sampleData} columns will be assigned colors internally.
#'   \code{clusterLegend} should be list of length equal to
#'   \code{ncol(sampleData)} with names equal to the colnames of
#'   \code{sampleData}. Each element of the list should be a either the format
#'   requested by \code{\link[NMF]{aheatmap}} (a vector of colors with names
#'   corresponding to the levels of the column of \code{sampleData}), or should
#'   be format of \code{ClusterExperiment}.
#' @param alignSampleData Logical as to whether should align the colors of the
#'   \code{sampleData} (only if \code{clusterLegend} not given and
#'   \code{sampleData} is not \code{NULL}).
#' @param breaks Either a vector of breaks (should be equal to length 52), or a
#'   number between 0 and 1, indicating that the breaks should be equally spaced
#'   (based on the range in the data) upto the `breaks' quantile, see
#'   \code{\link{setBreaks}}
#' @param unassignedColor color assigned to cluster values of '-1'
#'   ("unassigned").
#' @param missingColor color assigned to cluster values of '-2' ("missing").
#' @param ... for signature \code{matrix}, arguments passed to \code{aheatmap}. 
#'   For the other signatures, passed to the method for signature \code{matrix}.
#'   Not all arguments can be passed to aheatmap effectively, see details.
#' @param nFeatures integer indicating how many features should be used (if
#'   \code{clusterFeaturesData} is 'var' or 'PCA').
#' @param isSymmetric logical. if TRUE indicates that the input matrix is
#'   symmetric. Useful when plotting a co-clustering matrix or other sample by
#'   sample matrices (e.g., correlation).
#' @param overRideClusterLimit logical. Whether to override the internal limit 
#'   that only allows 10 clusterings/annotations. If overridden, may result in 
#'   incomprehensible errors from aheatmap. Only override this if you have a
#'   very large plotting device and want to see if aheatmap can render it.
#' @inheritParams clusterSingle
#'
#' @details The plotHeatmap function calls \code{\link[NMF]{aheatmap}} to draw
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
#'   \code{data} will be used; if missing, a new one will be created based on
#'   the \code{primaryCluster} of data. If equal to "hclust", then standard
#'   hierachical clustering of the transformed data will be used. If
#'   'orderSamplesValue' no clustering of the samples will be done, and instead
#'   the samples will be ordered as in the slot \code{orderSamples} of
#'   \code{data}. If equal to 'primaryCluster', again no clustering will be
#'   done, and instead the samples will be ordered based on grouping the samples
#'   to match the primaryCluster of \code{data}. If not one of these values,
#'   \code{clusterSamplesData} can be a character vector matching the
#'   clusterLabels (colnames of clusterMatrix).
#' @details If \code{data} is a matrix, then \code{sampleData} is a data.frame
#'   of annotation data to be plotted above the heatmap and
#'   \code{whSampleDataCont} gives the index of the column(s) of this dataset
#'   that should be consider continuous. Otherwise the annotation data for
#'   \code{sampleData} will be forced into a factor (which will be nonsensical
#'   for continous data). If \code{data} is a \code{ClusterExperiment} object,
#'   \code{sampleData} should refer to a index or column name of the
#'   \code{colData} slot of \code{data}. In this case \code{sampleData} will be
#'   added to any choices of clusterings chosen by the \code{whichClusters}
#'   argument (if any). If both clusterings and sample data are chosen, the
#'   clusterings will be shown closest to data (i.e. on bottom).
#' @details If \code{data} is a \code{ClusterExperiment} object,
#'   \code{clusterFeaturesData} is not a dataset, but instead indicates which
#'   features should be shown in the heatmap. "var" selects the
#'   \code{nFeatures} most variable genes (based on
#'   \code{transformation(assay(data))}); "PCA" results in a heatmap of the top
#'   \code{nFeatures} PCAs of the \code{transformation(assay(data))}.
#'   clusterFeaturesData can also be a vector of characters or integers,
#'   indicating the rownames or indices respectively of \code{assay(data)} that
#'   should be shown. For all of these options, the features are clustered based
#'   on the \code{visualizeData} data. Finally, in the \code{ClusterExperiment}
#'   version of \code{plotHeatmap}, \code{clusterFeaturesData} can be a list of
#'   indices or rownames, indicating that the features should be grouped
#'   according to the elements of the list, with blank (white) space between
#'   them (see \code{\link{makeBlankData}} for more details). In this case, no
#'   clustering is done of the features.
#' @details If \code{breaks} is a numeric value between 0 and 1, then
#'   \code{breaks} is assumed to indicate the upper quantile (on the log scale)
#'   at which the heatmap color scale should stop. For example, if
#'   \code{breaks=0.9}, then the breaks will evenly spaced up until the 0.9
#'   upper quantile of \code{data}, and then all values after the
#'   0.9 quantile will be absorbed by the upper-most color bin. This can help to
#'   reduce the visual impact of a few highly expressed genes (features). 
#' @details Note that plotHeatmap calls \code{\link[NMF]{aheatmap}} under the
#'   hood. This allows you to plot multiple heatmaps via
#'   \code{par(mfrow=c(2,2))}, etc. However, the dendrograms do not resize if
#'   you change the size of your plot window in an interactive session of R
#'   (this might be a problem for RStudio if you want to pop it out into a large
#'   window...).
#' @details Many arguments can be passed on to aheatmap, however, some are set 
#'   internally by \code{plotHeatmap.} In particular, setting the values of 
#'   \code{Rowv} or \code{Colv} will cause errors. \code{color} in 
#'   \code{aheatmap} is replaced by \code{colorScale} in \code{plotHeatmap.} The
#'   \code{annCol} to give annotation to the samples is replaced by the
#'   \code{sampleData}; moreover, the \code{annColors} option in \code{aheatmap}
#'   will also be set internally to give more vibrant colors than the default in
#'   \code{aheatmap} (for \code{ClusterExperiment} objects, these values can
#'   also be set in the \code{clusterLegend} slot ). Other options should be
#'   passed on to \code{aheatmap}, though they have not been all tested.
#' 
#' @return Returns (invisibly) a list with elements
#' \itemize{
#' \item{\code{aheatmapOut}}{ The output from the final call of
#' \code{\link[NMF]{aheatmap}}.}
#' \item{\code{sampleData}}{ the annotation data.frame given to the argument
#' \code{annCol} in \code{aheatmap}.}
#' \item{\code{clusterLegend}}{ the annotation colors given to the argument
#' \code{annColors} \code{aheatmap}.}
#' \item{\code{breaks}}{ The breaks used for \code{aheatmap}, after adjusting
#' for quantile.}
#' }
#' @author Elizabeth Purdom
#'
#' @export
#'
#' @examples
#' data(simData)
#'
#' cl <- rep(1:3,each=100)
#' cl2 <- cl
#' changeAssign <- sample(1:length(cl), 80)
#' cl2[changeAssign] <- sample(cl[changeAssign])
#' ce <- clusterExperiment(simCount, cl2, transformation=function(x){log2(x+1)})
#'
#' #simple, minimal, example. Show counts, but cluster on underlying means
#' plotHeatmap(ce)
#'
#' #assign cluster colors
#' colors <- bigPalette[20:23]
#' names(colors) <- 1:3
#' plotHeatmap(data=simCount, clusterSamplesData=simData,
#' sampleData=data.frame(cl), clusterLegend=list(colors))
#'
#' #show two different clusters
#' anno <- data.frame(cluster1=cl, cluster2=cl2)
#' out <- plotHeatmap(simData, sampleData=anno)
#'
#' #return the values to see format for giving colors to the annotations
#' out$clusterLegend
#'
#' #assign colors to the clusters based on plotClusters algorithm
#' plotHeatmap(simData, sampleData=anno, alignSampleData=TRUE)
#'
#' #assign colors manually
#' annoColors <- list(cluster1=c("black", "red", "green"),
#' cluster2=c("blue","purple","yellow"))
#'
#' plotHeatmap(simData, sampleData=anno, clusterLegend=annoColors)
#'
#' #give a continuous valued -- need to indicate columns
#' anno2 <- cbind(anno, Cont=c(rnorm(100, 0), rnorm(100, 2), rnorm(100, 3)))
#' plotHeatmap(simData, sampleData=anno2, whSampleDataCont=3)
#'
#' #compare changing breaks quantile on visual effect
#' \dontrun{
#' par(mfrow=c(2,2))
#' plotHeatmap(simData, colorScale=seqPal1, breaks=1, main="Full length")
#' plotHeatmap(simData,colorScale=seqPal1, breaks=.99, main="0.99 Quantile Upper
#' Limit")
#' plotHeatmap(simData,colorScale=seqPal1, breaks=.95, main="0.95 Quantile Upper
#' Limit")
#' plotHeatmap(simData, colorScale=seqPal1, breaks=.90, main="0.90 Quantile
#' Upper Limit")
#' }
#'
#' @rdname plotHeatmap
#' @aliases plotHeatmap
#' @importFrom stats hclust dist
#' @importFrom NMF aheatmap
setMethod(
    f = "plotHeatmap",
    signature = signature(data = "SummarizedExperiment"),
    definition = function(data, isCount=FALSE,transFun=NULL,...
    ){
      transformation<-.transData(assay(data),isCount=isCount,transFun=transFun,dimReduce="none")$transFun
      fakeCL<-sample(1:2,size=NCOL(data),replace=TRUE)
      fakeCE<-clusterExperiment(data,fakeCL,transformation=transformation)
      if("whichClusters" %in% names(list(...))) stop("cannot provide argument 'whichClusters' for input data of class Summarized Experiment")
      plotHeatmap(fakeCE,whichClusters="none",...)
})
#' @rdname plotHeatmap
setMethod(
  f = "plotHeatmap",
  signature = signature(data = "ClusterExperiment"),
  definition = function(data,
                        clusterSamplesData=c("hclust","dendrogramValue","orderSamplesValue","primaryCluster"),
                        clusterFeaturesData=c("var","all","PCA"), nFeatures=NULL,
                        visualizeData=c("transformed","centeredAndScaled","original"),
                        whichClusters= c("primary","workflow","all","none"),
                        sampleData=NULL,clusterFeatures=TRUE,
                        colorScale,
                       ...
  ){

    .convertTry<-function(x,tryResult){if(!inherits(tryResult,"try-error")) return(tryResult) else return(x)}

    ####
    ##Transform data and determine which features to use
    ####
    clusterFeaturesData <- .convertTry(clusterFeaturesData,
                                       try(match.arg(clusterFeaturesData),
                                           silent=TRUE))

    if(is.list(clusterFeaturesData)){
      groupFeatures<-clusterFeaturesData
      clusterFeaturesData<-unlist(clusterFeaturesData)
    }
    else groupFeatures<-NULL
    if(all(clusterFeaturesData %in% c("var","all","PCA"))){ #
        dimReduce=switch(clusterFeaturesData,
                         "var"="var",
                        "PCA"="PCA",
                        "all"="none")
        if(is.null(nFeatures)) nFeatures<-min(switch(clusterFeaturesData,"var"=500,"all"=nFeatures(data),"PCA"=50),nFeatures(data))
        wh<-1:NROW(data)
    }
    else{
      if(is.character(clusterFeaturesData)){#gene names
        if(is.null(rownames(data))) stop("Cannot give feature names in clusterFeaturesData unless assay(data) has rownames")
        else{
          wh<-match(clusterFeaturesData,rownames(data))
          if(all(is.na(wh))) stop("None of the feature names in clusterFeaturesData match rownames(assay(data))")
          if(any(is.na(wh))){
            warning("Not all of the feature names in clusterFeaturesData match rownames(assay(data))")
            wh<-na.omit(wh)
          }
        }
      }
      else{
          if(any(!clusterFeaturesData %in% 1:NROW(data))) stop("invalid indices for clusterFeaturesData")
          wh<-clusterFeaturesData
      }
      dimReduce<-"none"
    }
    transObj<-.transData(transFun = transformation(data), x=assay(data[wh,]), nPCADims=nFeatures,nVarDims = nFeatures,dimReduce = dimReduce)
    if(dimReduce%in%"PCA") wh<-1:nFeatures
    if(dimReduce=="var") wh<-transObj$whMostVar #give indices that will pull
    #########
    ##Assign visualization data and clusterFeaturesData
    #########
    #browser()
    visualizeData <- .convertTry(visualizeData,
                                 try(match.arg(visualizeData), silent=TRUE))
    if(is.character(visualizeData)){
      if(!visualizeData %in% c("transformed","centeredAndScaled","original")) stop("visualizeData value, '",visualizeData,"',is invalid option")
    }
    if(missing(colorScale)) {
      colorScale <- seqPal5
      if(is.character(visualizeData)) {
        if (visualizeData == "centeredAndScaled") {
          colorScale <- seqPal3
        }
      }
    }

    if(all(clusterFeaturesData=="PCA")) heatData<-transObj$x
    else{
        if(!is.data.frame(visualizeData) && !is.matrix(visualizeData)){
            heatData<-switch(visualizeData,
                            "original"=assay(data[wh,]),
                            "transformed"=transObj$x,
                            "centeredAndScaled"=t(scale(t(transObj$x),center=TRUE,scale=TRUE))
                            )
        }
        else{
            if(!all(dim(visualizeData)==dim(assay(data)))) stop("if give separate visualizeData, must be of same dimensions as assay(data)")
           heatData<-data.matrix(visualizeData)[wh,] #still use the variables identified above.
        }
    }

    ######
    #Make sampleData based on clusterings and columns of colData
    ######
    #Get clusterings
    if(is.character(whichClusters)) whCl<-.TypeIntoIndices(data,whClusters=whichClusters)
    else whCl<-whichClusters
    if(length(whCl)>0){
      if(!is.numeric(whCl)) stop("invalid whichClusters choices")
      if(!all(whCl %in% 1:nClusters(data))) stop("Indices in whichClusters invalid: not in 1 to nClusters(data)")
      clusterData<-clusterMatrixNamed(data)[,whCl,drop=FALSE]
    }
    else{
      if(whichClusters!="none") warning("given whichClusters value does not match any clusters")
      clusterData<-NULL
    }
    clLegend<-clusterLegend(data)[whCl] #note, this gives names even though not stored internally so will match, which plotHeatmap needs
    if(length(clLegend)==0) clLegend<-NULL
    #browser()
    #check user didn't give something different for colors
    userList<-list(...)
    userAlign<-"alignSampleData" %in% names(userList)
    userLegend<-"clusterLegend" %in% names(userList)
    if(userAlign | userLegend){ #if user asks for alignment, don't assign clusterLegend
      if(userLegend){
        userClLegend<-userList[["clusterLegend"]]
        #keep existing clLegend from clusterExperiment object if not conflict with user input:
        whNotShared<-which(!names(clLegend)%in% names(userClLegend))
        if(length(whNotShared)>0) clLegend<-c(userClLegend,clLegend[whNotShared]) else clLegend<-userClLegend
        clLegend<-.convertToAheatmap(clLegend, names=TRUE)
        userList<-userList[-grep("clusterLegend",names(userList))]
      }
      else{
        if(userAlign){
          al<-userList[["alignSampleData"]]
          if(al) clLegend<-NULL
        }
      }
    }

    #get colData values
    sData<-.pullSampleData(data,sampleData)
    #identify which numeric
    if(!is.null(sData)) whCont<-which(sapply(1:ncol(sData),function(ii){is.numeric(sData[,ii])}))
    whSampleDataCont<-NULL

    if(!is.null(clusterData) & !is.null(sData)){
      sampleData<-data.frame(clusterData,sData,stringsAsFactors=FALSE,check.names=FALSE)
      if(length(whCont)>0)  whSampleDataCont<-whCont+ncol(clusterData)
    }
    else{
      if(!is.null(clusterData)) sampleData<-clusterData
      if(!is.null(sData)){
        sampleData<-sData
        if(length(whCont)>0) whSampleDataCont<-whCont
      }
      if(is.null(sData) & is.null(clusterData)) sampleData<-NULL
    }
    
    
    ######
    #Create clusterSamplesData
    ######
    clusterSamplesData<-.convertTry(clusterSamplesData,try(match.arg(clusterSamplesData),silent=TRUE))
    #browser()
    if(is.logical(clusterSamplesData)) clusterSamples<-clusterSamplesData
    else {
      clusterSamples<-TRUE
      #browser()
      if(is.numeric(clusterSamplesData)){
          heatData<-heatData[,clusterSamplesData,drop=FALSE]
          if(!is.null(sampleData)) sampleData<-sampleData[clusterSamplesData,,drop=FALSE]
          clusterSamplesData<-heatData
          clusterSamples<-FALSE
      }
      else if(is.character(clusterSamplesData)){
          if(clusterSamplesData=="orderSamplesValue"){
            heatData<-heatData[,orderSamples(data),drop=FALSE]
            if(!is.null(sampleData)) sampleData<-sampleData[orderSamples(data), ,drop=FALSE]
            clusterSamplesData<-heatData
              clusterSamples<-FALSE
  
          }
          else if(clusterSamplesData=="primaryCluster"){
              heatData<-heatData[,order(primaryCluster(data))]
              if(!is.null(sampleData)) sampleData<-sampleData[order(primaryCluster(data)),,drop=FALSE]
              clusterSamplesData<-heatData
              clusterSamples<-FALSE
          }
          else if(clusterSamplesData=="dendrogramValue"){
              if(is.null(data@dendro_samples)){
                clusterSamplesData<-makeDendrogram(data)@dendro_samples
              }
              else{
                  clusterSamplesData<-data@dendro_samples
              }
          }
          else if(clusterSamplesData=="hclust"){
              #if hclust, then use the visualizeData data, unless visualizeData data is original, in which case use transformed
              clusterSamplesData <- heatData
  
              if(is.character(visualizeData)) {
                if(visualizeData=="original") {
                  transObj$x
                }
              }
          }
      }
      else stop("clusterSamplesData must be either character, or vector of indices of samples")
    }
    if(!is.null(groupFeatures)){
      #convert groupFeatures to indices on new set of data.
      groupFeatures<-lapply(groupFeatures,function(x){match(x,wh)})
      blankData<-makeBlankData(heatData,groupFeatures)
      heatData<-data.matrix(blankData$dataWBlanks)
      labRow<-blankData$rowNamesWBlanks
      #browser()
      clusterFeatures<-FALSE
    }
    else{
      labRow<-rownames(heatData)
    }
    do.call("plotHeatmap",c(list(data=heatData,
                clusterSamplesData=clusterSamplesData,
                clusterFeaturesData=heatData, #set it so user doesn't try to pass it and have something weird happen because dimensions wrong, etc.
                sampleData=sampleData,whSampleDataCont=whSampleDataCont,
                clusterSamples=clusterSamples,labRow=labRow,
                clusterLegend=clLegend,clusterFeatures=clusterFeatures,
                colorScale=colorScale),userList))


  })

#' @rdname plotHeatmap
setMethod(
    f = "plotHeatmap",
    signature = signature(data = "matrix"),
    definition = function(data,sampleData=NULL,
                          clusterSamplesData=NULL,
                          clusterFeaturesData=NULL,
                          whSampleDataCont=NULL,
                          clusterSamples=TRUE,showSampleNames=FALSE,
                          clusterFeatures=TRUE,showFeatureNames=FALSE,
                          colorScale=seqPal5,
                          clusterLegend=NULL,alignSampleData=FALSE,
                          unassignedColor="white",missingColor="grey", breaks=NA,isSymmetric=FALSE, overRideClusterLimit=FALSE,...
    ){


      ##########
      ##Deal with numeric matrix for heatmap ...
      ##########
      heatData<-data.matrix(data)
    aHeatmapArgs<-list(...)
    aHeatmapDefaultArgs<-as.list(args(NMF::aheatmap))
    getHeatmapValue<- function(string,value=NULL){ #note, doesn't work for pulling function 'reorder' so put in manually
        if(string %in% names(aHeatmapArgs)) val<-aHeatmapArgs[[string]]
        else{
            if(is.null(value)) val<-aHeatmapDefaultArgs[[string]]
            else val<-value
        }
        return(val)
    }  
    #browser()
    badValues<-c("Rowv","Colv","color","annCol","annColors")
    if(any(badValues %in% names(aHeatmapArgs))) stop("The following arguments to aheatmap cannot be set by the user in plotHeatmap:",paste(badValues,collapse=","))

    
    
      ###Create the clustering dendrogram:

    if(clusterSamples){
      if(inherits(clusterSamplesData, "dendrogram")){
        if(nobs(clusterSamplesData)!=ncol(heatData)) stop("clusterSamplesData dendrogram is not on same number of observations as heatData")
        dendroSamples<-clusterSamplesData
      }
      else{
          if(is.null(clusterSamplesData)){
              dendroSamples<-NULL
          }
          else{
              ##Call NMF:::cluster_mat so do the same thing:
              
              if(!is.data.frame(clusterSamplesData) & !is.matrix(clusterSamplesData)) stop("clusterSamplesData must either be dendrogram, or data.frame/matrix")
              clusterSamplesData<-data.matrix(clusterSamplesData)
              #check valid
              if(ncol(clusterSamplesData)!=ncol(heatData)) stop("clusterSamplesData matrix does not have on same number of observations as heatData")
#              browser()
              dendroSamples<-NMF:::cluster_mat(t(clusterSamplesData),param=TRUE,distfun=getHeatmapValue("distfun"),hclustfun=getHeatmapValue("hclustfun"),reorderfun=getHeatmapValue("reorderfun",value=function(d, w) reorder(d, w)))$dendrogram
              
              #dendroSamples<-as.dendrogram(stats::hclust(stats::dist(t(clusterSamplesData)))) #dist finds distances between rows
                 
         }
      }
    }
    else{
      clusterSamples<-NA
    }
    if(!is.na(clusterSamples) && clusterSamples && is.null(dendroSamples)) Colv<-TRUE #then just pass the data
    else Colv<-if(!is.na(clusterSamples) && clusterSamples) dendroSamples else clusterSamples
    
    if(isSymmetric){
        Rowv<-Colv
        Colv<-"Rowv"
     }
    else{
        if(clusterFeatures){
            if(inherits(clusterFeaturesData, "dendrogram")){
                if(nobs(clusterFeaturesData)!=ncol(heatData)) stop("clusterFeaturesData dendrogram is not on same number of observations as heatData")
                dendroFeatures<-clusterFeaturesData
            }
            else{
                if(is.null(clusterFeaturesData)){
                    dendroFeatures<-NULL
                }
                else{
                    ##Call NMF:::cluster_mat so do the same thing:
                    
                    if(!is.data.frame(clusterFeaturesData) & !is.matrix(clusterFeaturesData)) stop("clusterFeaturesData must either be dendrogram, or data.frame/matrix")
                    clusterFeaturesData<-data.matrix(clusterFeaturesData)
                    #check valid
                    if(ncol(clusterFeaturesData)!=ncol(heatData)) stop("clusterFeaturesData matrix not have on same number of observations as heatData")
                    dendroFeatures<-NMF:::cluster_mat(clusterFeaturesData,param=TRUE,distfun=getHeatmapValue("distfun"),hclustfun=getHeatmapValue("hclustfun"),reorderfun=getHeatmapValue("reorderfun",value=function(d, w) reorder(d, w)))$dendrogram
                    #                dendroFeatures<-as.dendrogram(stats::hclust(stats::dist(clusterFeaturesData))) #dist finds distances between rows
                    
                }
                
            }
        }
        else{
            clusterFeatures<-NA
        }
        if(!is.na(clusterFeatures) && clusterFeatures && is.null(dendroFeatures)) Rowv<-TRUE #then just pass the data
        else Rowv<-if(!is.na(clusterFeatures) && clusterFeatures) dendroFeatures else clusterFeatures
        
        
    }
      #browser()



      ##########
      ##Deal with annotation of samples (sampleData) ...
      ##########
      #check sampleData input:

      if(!is.null(sampleData)){
        if(!is.matrix(sampleData) & !is.data.frame(sampleData)) stop("sampleData must be a either a matrix or a data.frame")
        if(NCOL(data) != NROW(sampleData)) stop("sampleData must have same number of rows as columns of heatData")
        if(NCOL(sampleData)>10){
            if(overRideClusterLimit) warning("More than 10 annotations/clusterings can result in incomprehensible errors in aheamap. You have >10 but have chosen to override the internal stop by setting overRideClusterLimit=TRUE.")
            else stop("More than 10 annotations/clusterings cannot be reliably shown in plotHeatmap. To override this limitation and try for yourself, set overRideClusterLimit=TRUE.")
        }
                    ###Make sampleData explicitly factors, except for whSampleDataCont
        ###(not sure why this simpler code doesn't give back data.frame with factors: annCol<-apply(annCol,2,function(x){factor(x)}))
        #browser()
        tmpDf<-do.call("data.frame",lapply(1:ncol(sampleData),function(ii){factor(sampleData[,ii])}))
        names(tmpDf)<-colnames(sampleData)
        if(!is.null(whSampleDataCont)) tmpDf[,whSampleDataCont]<-sampleData[,whSampleDataCont]
        annCol<-tmpDf
        #browser()
        convertNames <- TRUE
        if(is.null(clusterLegend)){ #assign default colors
          convertNames <- TRUE
          if(is.null(whSampleDataCont) || length(whSampleDataCont)<ncol(annCol)){ #if not all continuous
            #Pull out the factors to assign clusters
            if(!is.null(whSampleDataCont)) tmpDf<- annCol[,-whSampleDataCont,drop=FALSE] else tmpDf<-annCol
            tmpDfNum<-do.call("cbind",lapply(1:ncol(tmpDf),function(ii){.convertToNum(tmpDf[,ii])}))
            colnames(tmpDfNum)<-colnames(tmpDf)
            if(alignSampleData){
              #align the clusters and give them colors
              alignObj<-plotClusters(tmpDfNum ,plot=FALSE,unassignedColor=unassignedColor, missingColor=missingColor)
              clusterLegend<-lapply(1:ncol(tmpDfNum),function(ii){
                xNum<-tmpDfNum[,ii]
                xOrig<-tmpDf[,ii]
                colMat<-alignObj$clusterLegend[[ii]]
                m<-match(colMat[,"clusterIds"],as.character(xNum))
                colMat<-cbind(colMat,"name"=as.character(xOrig)[m])
                return(colMat)
              })
              #browser()
            }
            else{#give each distinct colors, compared to row before
              maxPerAnn<-apply(tmpDfNum,2,max) #max cluster value (not including -1,-2)
              maxPreviousColor<-c(0,head(cumsum(maxPerAnn),-1))
              pal<-rep(bigPalette,length=sum(maxPerAnn)) #make sure don't run out of colors

              clusterLegend<-lapply(1:ncol(tmpDfNum),FUN=function(ii){
                facInt<-tmpDfNum[,ii]
                facOrig<-tmpDf[,ii]
                add<-maxPreviousColor[[ii]]
                colors<-pal[1:max(facInt)+add]
                cols<-cbind("clusterIds"=levels(factor(facInt[facInt>0])),"color"=colors,"name"=levels(facOrig[facInt>0]))
                if(any(facInt== -1)) cols<-rbind(cols,c("clusterIds"="-1","color"=unassignedColor,"name"="-1") )
                if(any(facInt== -2)) cols<-rbind(cols,c("clusterIds"="-2","color"=unassignedColor,"name"="-2") )
                cols<-cols[order(cols[,"name"]),]
                return(cols)
              })
            }
            names(clusterLegend)<-colnames(tmpDf)

          }
        }
        #browser()
        if(!any(sapply(clusterLegend,function(x){is.null(dim(x))}))) {
          annColors<-.convertToAheatmap(clusterLegend, names=convertNames)
        } else {
          annColors<-clusterLegend #in case give in format wanted by aheatmap to begin with
        }
      }
      else{
        annCol<-NA
        annColors<-NA
      }

      #############
      # put into aheatmap
      #############
      breaks<-setBreaks(data=heatData,breaks=breaks)
      #browser()
      out<-NMF::aheatmap(heatData,
                         Rowv =Rowv,Colv = Colv,
                         color = colorScale, scale = getHeatmapValue("scale","none"),
                         annCol = annCol,annColors=annColors,breaks=breaks,...)

      #############
      # add labels to clusters at top of heatmap
      #############

      if(!any(is.na(annCol))){
        newName<-NMF:::vplayout(NULL) #will be 1 greater (hopefully!) this is fragile. Don't know if it will always work.
        newNameList<-strsplit(newName,"\\.")[[1]]
        oldIndex<-as.numeric(newNameList[[3]])-1
        newNameList[[3]]<-oldIndex
        oldName<-paste(newNameList,collapse=".")
        grid::seekViewport(sprintf("aheatmap-%s",oldName))
        NMF:::vplayout(3,4:5)
        #grid::grid.rect()
        y <- seq(0,1,length=ncol(annCol))
        n<-ncol(annCol)
        y = cumsum(rep(8, n)) - 4 + cumsum(rep(2, n))
        #		grid::grid.points(x = grid::unit(rep(0,length(y)),"npc"),y = grid::unit(y[n:1], "bigpts"))
        grid::grid.text(colnames(annCol), x = grid::unit(rep(0.05,length(y)),"npc"),y = grid::unit(y[n:1], "bigpts"), vjust = 0.5, hjust = 0,gp= grid::gpar(fontsize=10))
        grid::upViewport() #close it
        grid::upViewport() #close it
      }

      invisible(list(aheatmapOut=out,sampleData=annCol,clusterLegend=clusterLegend,breaks=breaks))
    }
)

#' @rdname plotHeatmap
#' @aliases plotCoClustering
#' @param invert logical determining whether the coClustering matrix should be 
#'   inverted to be 1-coClustering for plotting. By default, if the diagonal 
#'   elements are all zero, invert=TRUE, and otherwise invert=FALSE. If
#'   coClustering matrix is not a 0-1 matrix (e.g. if equal to a distance matrix
#'   output from \code{\link{clusterSingle}}, then the user should manually set
#'   this parameter to FALSE.)
#' @details \code{plotCoClustering} is a convenience function to plot the heatmap
#' of the co-clustering matrix stored in the \code{coClustering} slot of a
#' \code{ClusterExperiment} object.
#' @export
setMethod(
  f = "plotCoClustering",
  signature = "ClusterExperiment",
  definition = function(data, invert= ifelse(!is.null(data@coClustering) && all(diag(data@coClustering)==0), TRUE, FALSE), ...){
    if(is.null(data@coClustering)) stop("coClustering slot is empty")
      if(invert) data@coClustering<-1-data@coClustering
    fakeCE<-clusterExperiment(data@coClustering,
                              clusterMatrix(data),
                              transformation=function(x){x},
                              clusterInfo=clusterInfo(data),
                              clusterTypes=clusterTypes(data),
                              orderSamples=orderSamples(data),
                              dendro_samples=data@dendro_samples,
                              dendro_clusters=data@dendro_clusters,
                              dendro_index=data@dendro_index,
                              primaryIndex=data@primaryIndex
                              
                              
    )
    clusterLegend(fakeCE)<-clusterLegend(data)
    colData(fakeCE)<-colData(data)
    plotHeatmap(fakeCE,isSymmetric=TRUE,clusterFeaturesData="all",...)
  })
