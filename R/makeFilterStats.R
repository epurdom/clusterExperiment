#' @name makeFilterStats
#' @title Calculate filtering statistics
#' @description Function for calculating, per row (gene), built-in statistical 
#'    functions that might be used for filtering. 
#' @param object object from which user wants to calculate per-row statistics
#' @param filterStats character vector of statistics to calculate. 
#' 	  Must be one of the character values given by \code{listBuildInFilterStats()}.
#' @param ... Values passed on the the 'SummarizedExperiment' method.
#' @inheritParams makeReducedDims
#' @inheritParams transformData
#'
#' @examples
#' data(simData)
#' listBuiltInFilterStats()
#' scf<-makeFilterStats(simData,filterStats=c("var","mad"))
#' scf
#' scfFiltered<-filterData(scf,type="mad",percentile=10)
#' scfFiltered
#' assay(scfFiltered)[1:10,1:10]
#' @rdname makeFilterStats
#' @aliases makeFilterStats,SummarizedExperiment-method
#' @export
setMethod(
  f = "makeFilterStats",
  signature = "SummarizedExperiment",
  definition = function(object,filterStats=listBuiltInFilterStats(),transFun=NULL,isCount=FALSE)
{

  ###################
  ##Check user inputs
  ###################
  #check valid options for reducedDims
  validDim<-listBuiltInFilterStats()
  filterStats<-unique(filterStats)
  if(!all(filterStats %in% validDim)){
	  stop("Not all of the filterStats given are valid")
  }
  
  ###################
  ##Clean up data:
  ###################
  #transform data
  x<-transformData(object,transFun=transFun,isCount=isCount)

  ###################
  ##Do loop over filterStats values:
  ###################
  if('abscv' %in% filterStats){
	  #origfilterStats<-filterStats
  #if abscv, add var and mean since calculating anyway and make abscv last one
	  filterStats<-unique(c(filterStats,"var","mean"))
	  doCV<-TRUE
	  whCV<-grep('abscv',filterStats)
	  filterStats<-c(filterStats[-whCV])
  }
  else doCV<-FALSE
  filterStatData<-sapply(filterStats,function(statName){
	  f<-.matchToStats[[statName]]
	  f(x)
  })
  if(doCV){
	  filterStatData<-cbind(filterStatData, "abscv"=sqrt(filterStatData[,"var"])/abs(filterStatData[,"mean"]))
	  #filterStatData<-filterStatData[,origfilterStats] #put it in order, though user shouldn't depend on it.
  }
  filterStats(object)<-filterStatData #should leave in place existing ones, update conflicting ones, and add new ones!
  return(object)

}
)
#' @rdname makeFilterStats
#' @export
setMethod(
  f = "makeFilterStats",
  signature = "matrix",
  definition = function(object,...)
{
	makeFilterStats(SummarizedExperiment(object),...)
}
)

#' @rdname makeFilterStats
#' @export
#' @param whichClusterIgnoreUnassigned indicates clustering that should be used to filter out unassigned samples from the calculations. If \code{NULL} no filtering of samples will be done. See details for more information.
#' @details \code{whichClusterIgnoreUnassigned} is only an option when applied to a
#' \code{ClusterExperiment} classs and indicates that the filtering statistics should be
#' calculated based on samples that are unassigned by the
#'  designated clustering. The name given to the filter in this case is of the form
#'  \code{<filterStats>_<clusterLabel>}, i.e. the clustering label of the clustering is
#'  appended to the standard name for the filtering statistic.
#'
setMethod(
  f = "makeFilterStats",
  signature = "ClusterExperiment",
  definition = function(object,whichClusterIgnoreUnassigned=NULL,filterStats=listBuiltInFilterStats(),...)
{
	if(any(c("transFun","isCount") %in% names(list(...)))) 
		stop("The internally saved transformation function of a ClusterExperiment object must be used when given as input and setting 'transFun' or 'isCount' for a 'ClusterExperiment' is not allowed.")  
	filterStats<-unique(filterStats)
	if(!is.null(whichClusterIgnoreUnassigned)){
		whCluster<-.TypeIntoIndices(object,whichClusterIgnoreUnassigned)
		if(length(whCluster)>1) warning("'whichClusterIgnoreUnassigned' corresponds to multiple clusterings. Ignoring input")
		else if(length(whCluster)==0) warning("'whichClusterIgnoreUnassigned' does not correspond to a clustering. Ignoring input")
		else{
			#give new names to filters to indicate based on clustering.
			newNames<-paste(filterStats,clusterLabels(object)[whCluster],sep="_")
			whDo<-which(!newNames %in% filterNames(object))
			if(length(whDo)>0){
				whAssigned<-which(clusterMatrix(object)[,whCluster]>0)
				if(length(whAssigned)>0){
					out<-makeFilterStats(object[,whAssigned],filterStats=filterStats[whDo],...)
					whNew<-match(filterStats[whDo],filterNames(out))
					filterNames(out)[whNew]<-newNames[whDo]	
				}
				else 
					stop("All samples are unassigned for clustering", clusterLabels(object)[whCluster])
				
			}
		}
	}
	else{
		out<-makeFilterStats(as(object,"SingleCellExperiment"),filterStats=filterStats,transFun=transformation(object),...)
	}
#	
	filterStats(object)<-filterStats(out)
	return(object)
}
)
#' @rdname makeFilterStats
#' @export
listBuiltInFilterStats<-function(){c('var', 'abscv', 'mad','mean','iqr','median')}
#' @importFrom matrixStats rowVars rowMeans2 rowMads rowMedians rowIQRs
.matchToStats<-SimpleList(
	'var'=matrixStats::rowVars,
	'mad'=matrixStats::rowMads,
	'mean'=matrixStats::rowMeans2,
	'iqr'=matrixStats::rowIQRs,
	'median'=matrixStats::rowMedians)


#' @rdname makeFilterStats
#' @aliases filterData
#' @param type The type of filtering statistic to use to filter. 
#' @param cutoff numeric. A value at which to filter the rows (genes) for the test statistic
#' @param percentile numeric. Either a number between 0,1 indicating what percentage of the rows (genes) to keep or an integer value indicated the number of rows (genes) to keep
#' @param absolute whether to take the absolute value of the filter statistic
#' @param keepLarge logical whether to keep rows (genes) with large values of the test statistic or small values of the test statistic. 
#' @details Note that \code{filterData} returns a SingleCellExperiment object. To get the actual data out use either assay or \code{\link{transformData}} if transformed data is desired.
#' @return A SingleCellExperiment object with the rows (genes) removed based on filters
#' @export
#' @importFrom stats quantile
setMethod( "filterData","SingleCellExperiment",
	function(object,type,cutoff,percentile, absolute=FALSE,keepLarge=TRUE){
	stat<-if(absolute) abs(filterStats(object,type,checkValid=TRUE)) else filterStats(object,type,checkValid=TRUE)
	if(missing(cutoff) & missing(percentile)) stop("must provide one of cutoff or percentile")
	if(!missing(cutoff) & !missing(percentile)) stop("can only provide one of cutoff or percentile")
	if(!missing(cutoff)){
		whKeep<- if(keepLarge) which(stat>cutoff) else which(cutoff > stat)
	}
	if(!missing(percentile)){
		if(0<percentile & percentile <1){
			quantile<- quantile(stat,probs=if(keepLarge) percentile else 1-percentile)
			whKeep<-if(keepLarge) which(stat>quantile) else which(stat<quantile)
		}
		else{
			if(is.na(percentile) || percentile>=1){
				if(is.na(percentile) || percentile>NROW(object)){
					warning("the number of most features requested after filtering is either missing or larger than the number of features. Will not do any filtering")
					whKeep<-1:NROW(object)
				}
				else whKeep<- order(stat,decreasing=ifelse(keepLarge,TRUE,FALSE))[1:percentile]
			}
			else stop("Invalid value for percentile. Must be either between 0,1 or a positive integer number to keep")
		}
	}
	object[whKeep,]

})	