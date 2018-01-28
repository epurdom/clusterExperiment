#' Heatmap for showing clustering results and more
#' 
#' Make heatmap with color scale from one matrix and hiearchical clustering of 
#' samples/features from another. Also built in functionality for showing the 
#' clusterings with the heatmap. Builds on \code{\link[NMF]{aheatmap}} function 
#' of \code{NMF} package.
#' 
#' @docType methods
#' @param sampleData If input to \code{data} is either a
#'   \code{\link{ClusterExperiment}},or \code{SummarizedExperiment} object or
#'   \code{SingleCellExperiment}, then \code{sampleData} must index the 
#'   sampleData stored as a \code{DataFrame} in \code{colData} slot of the 
#'   object. Whether that data is continuous or not will be determined by the 
#'   properties of \code{colData} (no user input is needed). If input to
#'   \code{data} is matrix, \code{sampleData} is a matrix of additional data on
#'   the samples to show above heatmap. In this case, unless indicated by
#'   \code{whSampleDataCont}, \code{sampleData} will be converted into factors,
#'   even if numeric. ``-1'' indicates the sample was not assigned to a cluster
#'   and gets color `unassignedColor' and ``-2`` gets the color 'missingColor'.
#' @param data data to use to determine the heatmap. Can be a matrix, 
#'   \code{\link{ClusterExperiment}},
#'   \code{\link[SingleCellExperiment]{SingleCellExperiment}}or 
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}} object. The 
#'   interpretation of parameters depends on the type of the input to
#'   \code{data}.
#' @param whSampleDataCont Which of the \code{sampleData} columns are continuous
#'   and should not be converted to counts. \code{NULL} indicates no additional 
#'   \code{sampleData}. Only used if \code{data} input is matrix.
#' @param visualizeData either a character string, indicating what form of the 
#'   data should be used for visualizing the data (i.e. for making the 
#'   color-scale), or a data.frame/matrix with same number of samples as 
#'   \code{assay(data)}. If a new data.frame/matrix, any character arguments to 
#'   clusterFeaturesData will be ignored.
#' @param clusterSamplesData If \code{data} is a matrix, 
#'   \code{clusterSamplesData} is either a matrix that will be used by 
#'   \code{hclust} to define the hiearchical clustering of samples (e.g. 
#'   normalized data) or a pre-existing dendrogram that clusters the samples. If
#'   \code{data} is a \code{ClusterExperiment} object, \code{clusterSamplesData}
#'   should be either character or integers or logical which indicates how (and 
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
#'   \code{sampleData} columns will be assigned colors internally. See details
#'   for more.
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
#'   Not all arguments can be passed to \code{aheatmap} effectively, see
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
#'   visual impact of a few highly expressed genes (features).
#' @details Note that plotHeatmap calls \code{\link[NMF]{aheatmap}} under the 
#'   hood. This allows you to plot multiple heatmaps via 
#'   \code{par(mfrow=c(2,2))}, etc. However, the dendrograms do not resize if 
#'   you change the size of your plot window in an interactive session of R 
#'   (this might be a problem for RStudio if you want to pop it out into a large
#'   window...). Also, plotting to a pdf adds a blank page; see help pages of 
#'   \code{\link[NMF]{aheatmap}} for how to turn this off.
#' @details \code{clusterLegend} takes the place of argument \code{annColors}
#'   from \code{aheatmap} for giving colors to the annotation on the heatmap.
#'   \code{clusterLegend} should be list of length equal to 
#'   \code{ncol(sampleData)} with names equal to the colnames of 
#'   \code{sampleData}. Each element of the list should be a either the format 
#'   requested by \code{\link[NMF]{aheatmap}} (a vector of colors with names 
#'   corresponding to the levels of the column of \code{sampleData}), or should 
#'   be format of the \code{clusterLegend} slot in a \code{ClusterExperiment}
#'   object. Color assignments to the rows/genes should also be passed via
#'   \code{clusterLegend} (assuming \code{annRow} is an argument passed to
#'   \code{...}). If \code{clusterFeaturesData} is a \emph{named} list
#'   describing groupings of genes then the colors for those groups can be given
#'   in \code{clusterLegend} under the name "Gene Group".
#' @details If you have a factor with many levels, it is important to note that
#'   \code{\link[NMF]{aheatmap}} does not recycle colors across factors in the
#'   \code{sampleData}, and in fact runs out of colors and the remaining levels
#'   get the color white. Thus if you have many factors or many levels in those
#'   factors, you should set their colors via \code{clusterLegend}.
#' @details Many arguments can be passed on to \code{aheatmap}, however, some are set 
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
#' @seealso \code{\link[NMF]{aheatmap}}, \code{\link{makeBlankData}}
#' @export
#'
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
    signature = signature(data = "SingleCellExperiment"),
    definition = function(data, isCount=FALSE,transFun=NULL,...
    ){
     #get transformation function
	 transformation<-.makeTransFun(transFun=transFun,isCount=isCount) 
	 fakeCL<-sample(1:2,size=NCOL(data),replace=TRUE)
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
#' @param nBlankLines Only applicable if input is \code{ClusterExperiment} object. Indicates the number of lines to put between groups of features if \code{clusterFeaturesData} gives groups of genes (see details and \code{\link{makeBlankData}}).  
setMethod(
  f = "plotHeatmap",
  signature = signature(data = "ClusterExperiment"),
  definition = function(data,
    clusterSamplesData=c("dendrogramValue", "hclust", "orderSamplesValue", "primaryCluster"),
    clusterFeaturesData="var", nFeatures=NA,
    visualizeData=c("transformed","centeredAndScaled","original"),
    whichClusters= c("primary","workflow","all","none"),
    sampleData=NULL,clusterFeatures=TRUE, nBlankLines=2,
    colorScale,
   ...
  ){

    .convertTry<-function(x,tryResult){if(!inherits(tryResult,"try-error")) return(tryResult) else return(x)}

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
          colorScale <- seqPal3
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
	    	heatData<-assay(data)
	    }
		else{
			if(length(clusterFeaturesData)==1 && isPossibleReducedDims(data,clusterFeaturesData)){
				##### Dimensionality reduction ####
				if(!isReducedDims(data,clusterFeaturesData)){
					data<-makeReducedDims(data,reducedDims=clusterFeaturesData,maxDims=nFeatures)
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
						data<-makeFilterStats(data,filterStats=clusterFeaturesData)
					}
					if(is.na(nFeatures)) nFeatures<-min(NROW(data),500)
					data<-filterData(data,filterStats=clusterFeaturesData,percentile=nFeatures)
				}
			    else{
					### Other character values ####
					if(is.character(clusterFeaturesData)){#gene names
				        if(length(clusterFeaturesData)==1 && clusterFeaturesData=="all") 
							whRows<-1:NROW(data)
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
					  if(any(!clusterFeaturesData %in% 1:NROW(data))) stop("invalid indices for clusterFeaturesData")
					  whRows<-clusterFeaturesData
					}
					data<-data[whRows,]
				}
				heatData<-switch(visualizeData,
		                    "original"=assay(data),
		                    "transformed"=transformData(data),
		                    "centeredAndScaled"=t(scale(t(transformData(data)), center=TRUE, scale=TRUE))
		        )
			}
		}
	}
	else{
		heatData<-visualizeData
	}
 

    ######
    #Make sampleData based on clusterings and columns of colData
    ######
	#---
    #Get clusterings
	#---
    whCl<-.TypeIntoIndices(data,whClusters=whichClusters)
	#
    if(length(whCl)>0){
      clusterData<-clusterMatrixNamed(data)[,whCl,drop=FALSE]
    }
    else{
      if(any( whichClusters!="none")) warning("given whichClusters value does not match any clusters, none will be plotted")
      clusterData<-NULL
    }
	#---
    #get colData values and subset to those asked for by user
	#---
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
    
    #------
    #check user didn't give something different for colors
    #------
    clLegend<-clusterLegend(data)[whCl] #note, clusterLegend gives names even though not stored internally with @clusterLegend so will match, which plotHeatmap needs
    if(length(clLegend)==0) clLegend<-NULL


    userList<-list(...)
    userAlign<-"alignSampleData" %in% names(userList) & !is.null(userList$alignSampleData)
    userLegend<-"clusterLegend" %in% names(userList) & !is.null(userList$clusterLegend)
	if(userAlign | userLegend){ #if user asks for alignment, don't assign clusterLegend
      if(userLegend){
        userClLegend<-userList[["clusterLegend"]]
		annotNames<-colnames(sampleData)
		if("annRow" %in% names(userList)) annotNames<-c(annotNames,colnames(userList$annRow))
		if(!is.null(groupFeatures)) annotNames<-c(annotNames,"Gene Group")
		userClLegend<-userClLegend[names(userClLegend) %in% annotNames]
		if(length(userClLegend)==0){
			warning("names of list given by user in clusterLegend do not match clusters nor sampleData chosen. Will be ignored.")
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
          al<-userList[["alignSampleData"]]
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
          else if(clusterSamplesData=="dendrogramValue"){
              if(is.null(data@dendro_samples)){
                clusterSamplesData<-try(makeDendrogram(data)@dendro_samples,silent = TRUE)
				if(inherits(clusterSamplesData, "try-error")){
					warning("cannot make dendrogram from 'data' with default makeDendrogram options. Ordering by primary cluster without dendrogram")
					clusterSamplesData<-"primaryCluster"
				}
              }
              else{
                  clusterSamplesData<-data@dendro_samples
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
	              if(!is.null(sampleData)) sampleData<-sampleData[order(primaryCluster(data)),,drop=FALSE]
	              clusterSamplesData<-heatData
	              clusterSamples<-FALSE			  	
			  }
          }
          if(is.character(clusterSamplesData) && clusterSamplesData=="hclust"){
              #if hclust, then use the visualizeData data
			  #unless visualizeData data is original, in which case use transformed (and possibly filtered)
			  if(is.character(visualizeData) && visualizeData=="original") 
				  clusterSamplesData<-transformData(data)
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
	  blankData<-makeBlankData(heatData,groupFeatures,nBlankLines=nBlankLines)
	  #replace heatData with one with blanks -- won't cluster them now...
      heatData<-data.matrix(blankData$dataWBlanks)
      labRow<-blankData$rowNamesWBlanks
      clusterFeatures<-FALSE
	  
	  if(!is.null(names(groupFeatures))){
		  annRow<-list("Gene Group"=factor(blankData$groupNamesWBlanks,levels=names(groupFeatures)))
		  #show color-coding of gene groupings:
		  if("annRow" %in% names(userList)){
			  stop("Cannot provide 'annRow' if grouping features")
			  # if(is.list(userList$annRow)) userList$annRow<-c(userList$annRow,annRow)
			  # else userList$annRow<-c(list(userList$annRow),annRow)
		  }
		  else userList$annRow<-annRow
		  if(!"Gene Group" %in% names(clLegend)){
			  nGroups<-length(groupFeatures)
			  groupColors<-bigPalette[1:nGroups]	
			  names(groupColors)<-levels(annRow[["Gene Group"]])
		  	  clLegend<-c(clLegend, "Gene Group"=list(groupColors))
		  }
	  }
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
    signature = signature(data = "matrix"),
    definition = function(data,sampleData=NULL,
                          clusterSamplesData=NULL,
                          clusterFeaturesData=NULL,
                          whSampleDataCont=NULL,
                          clusterSamples=TRUE,showSampleNames=FALSE,
                          clusterFeatures=TRUE,showFeatureNames=FALSE,
                          colorScale=seqPal5,
                          clusterLegend=NULL,alignSampleData=FALSE,
                          unassignedColor="white",missingColor="grey", breaks=NA,isSymmetric=FALSE, overRideClusterLimit=FALSE, plot=TRUE,...
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
    badValues<-c("Rowv","Colv","color","annCol","annColors")
	replacedValues<-c("clusterSamplesData","clusterFeaturesData","colorScale","sampleData","clusterLegend")
    if(any(badValues %in% names(aHeatmapArgs))) stop("The following arguments to aheatmap cannot be set by the user in plotHeatmap:",paste(badValues,collapse=","),". They are over-ridden by: ",paste(replacedValues,collapse=","))

    
    
      ##########
      ###Create the clustering dendrogram (samples):
      ##########

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
    
        ##########
        ###Create the clustering dendrogram (features):
        ##########
    if(isSymmetric){
        Rowv<-Colv
        Colv<-"Rowv"
     }
    else{
        if(clusterFeatures){
            if(inherits(clusterFeaturesData, "dendrogram")){
                if(nobs(clusterFeaturesData)!=nrow(heatData)) stop("clusterFeaturesData dendrogram is not on same number of observations as heatData")
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
        #-------------------
		###Make sampleData explicitly factors, except for whSampleDataCont
        #-------------------
		
		#--- check that no ordered factors...
        anyOrdered<-sapply(1:ncol(sampleData),function(ii){is.ordered(sampleData[,ii])})
		if(any(anyOrdered)) stop("The function aheatmap in the NMF package that is called to create the heatmap does not currently accept ordered factors (https://github.com/renozao/NMF/issues/83)")
	        ###(not sure why this simpler code doesn't give back data.frame with factors: annCol<-apply(annCol,2,function(x){factor(x)}))
        tmpDf<-do.call("data.frame", lapply(1:ncol(sampleData), function(ii){ factor(sampleData[,ii]) }))
        names(tmpDf)<-colnames(sampleData)
        if(!is.null(whSampleDataCont)){
        	if(any(logical(whSampleDataCont))) whSampleDataCont<-which(whSampleDataCont)
    			if(length(whSampleDataCont)>0) tmpDf[,whSampleDataCont]<-sampleData[,whSampleDataCont]
        } 
        annCol<-tmpDf
        convertNames <- TRUE
		
        ##########
        ##Deal with colors ...
        ##########
		#-----
		#assign default colors
		#-----
        if(is.null(clusterLegend)){ 
          convertNames <- TRUE
		  #------
		  #if not all continuous, assign default colors from palette
		  #------
          if(is.null(whSampleDataCont) || length(whSampleDataCont)<ncol(annCol)){ 
			#Pull out those that are factors to assign clusters
            if(!is.null(whSampleDataCont)) tmpDf<- annCol[,-whSampleDataCont,drop=FALSE] else tmpDf<-annCol
			#make numeric
            tmpDfNum<-do.call("cbind", lapply(1:ncol(tmpDf), function(ii){ .convertToNum(tmpDf[,ii]) }))
            colnames(tmpDfNum)<-colnames(tmpDf)
            if(alignSampleData){
              #align the clusters and give them colors
              alignObj<-plotClusters(tmpDfNum ,plot=FALSE,unassignedColor=unassignedColor, missingColor=missingColor)
			  mkLegend<-function(ii){
                xNum<-tmpDfNum[,ii]
                xOrig<-tmpDf[,ii]
                colMat<-alignObj$clusterLegend[[ii]]
                m<-match(colMat[,"clusterIds"],as.character(xNum))
                colMat<-cbind(colMat[,c("clusterIds","color")],"name"=as.character(xOrig)[m])
                return(colMat)
              }
              clusterLegend<-lapply(1:ncol(tmpDfNum),mkLegend)
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
		#-----
		# Convert to aheatmap format, if needed
		#-----
        if( any(sapply(clusterLegend,function(x){!is.null(dim(x))}))) {
          annColors<-.convertToAheatmap(clusterLegend, names=convertNames)
        } else {
          annColors<-clusterLegend #in case give in format wanted by aheatmap to begin with; actually .convertToAheatmap would probably handle this fine too. Not clear we need catch.
        }
		#-----
		# remove any unused level colors to clean up legend, 
		# give names to colors if not given, and 
		# make them in same order as in annCol factor
		#-----
		whInAnnColors<-which(names(annColors)%in% colnames(annCol))
		if(!is.null(whSampleDataCont) & length(whSampleDataCont)>0){
			whInAnnColors<-setdiff(whInAnnColors,whSampleDataCont)
		} 
		prunedList<-lapply(whInAnnColors,function(ii){
			nam<-names(annColors)[[ii]] #the name of variable
			x<-annColors[[ii]] ##list of colors
			levs<-levels(annCol[,nam])
			if(length(x)<length(levs)) stop("number of colors given for ",nam," is less than the number of levels in the data")
			#Note to self:
			#It appears that if there is only 1 level of a factor, aheatmap 
			#doesn't plot it unless there are colors longer than 1
			#But appears only problem with extra colors are just that the first ones must match the levels -- not being matched (or at least not well) to names of the colors. 
			#So going to leave "extra" colors in place, but put them at the end.
			if(is.null(names(x))){
				#if user didn't give names to colors, assign them in order of levels
				#x<-x[1:length(levs)] #shorten if needed
				names(x)<-levs
			}
			else{
				if(any(!levs %in% names(x))) stop("colors given for ",nam," do not cover all levels in the data")
				whlevs<-match(levs,names(x))
				x<-c(x[whlevs],x[-whlevs])
			}
			return(x)
		})
		names(prunedList)<-names(annColors)[whInAnnColors]
		annColors[whInAnnColors]<-prunedList
		
		#-----
		# Give default colors to any that are missing (because aheatmap runs out of colors...)
		#-----
		
      }
      else{ #no sampleData provided -- just a heatmap with no annotation
        annCol<-NA
        annColors<-NA
      }

      #############
      # put into aheatmap
      #############
      breaks<-setBreaks(data=heatData,breaks=breaks)
	  if(plot){
	      out<-NMF::aheatmap(heatData,
	                         Rowv =Rowv,Colv = Colv,
	                         color = colorScale, scale = getHeatmapValue("scale","none"),
	                         annCol = annCol,annColors=annColors,breaks=breaks,...)

	      #############
	      # add labels to clusters at top of heatmap
	      #############
		  if(!is.null(dim(annCol))){
			  newName<-NMF:::vplayout(NULL) #will be 1 greater (hopefully!) this is fragile. Don't know if it will always work.
			#problems with development version of NMF
	        newNameList<-try(strsplit(newName,"\\.")[[1]])
			if(!inherits(newNameList,"try-error")){
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

	      }
	  }
	  else out<-NULL

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
    fakeCE<-ClusterExperiment(data@coClustering,
                              clusterMatrix(data),
                              transformation=function(x){x},
                              clusterInfo=clusteringInfo(data),
                              clusterTypes=clusterTypes(data),
                              orderSamples=orderSamples(data),
                              dendro_samples=data@dendro_samples,
                              dendro_clusters=data@dendro_clusters,
                              dendro_index=data@dendro_index,
                              dendro_outbranch=data@dendro_outbranch,
                              primaryIndex=data@primaryIndex,
							  checkTransformAndAssay=FALSE
                              
                              
    )
    clusterLegend(fakeCE)<-clusterLegend(data)
    colData(fakeCE)<-colData(data)
    plotHeatmap(fakeCE,isSymmetric=TRUE,clusterFeaturesData="all",...)
  })
