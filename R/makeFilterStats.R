#' @name reduceFunctions
#' @title Filtering statistics and Dimensionality Reduction Functions
#' @description Functions for calculating and manipulating either filtering statistics, stored in rowData, or the dimensionality reduction results, stored in reducedDims.
#' @aliases reduceFunctions makeFilterStats makeFilterStats,SummarizedExperiment-method 
#' @param object object from which user wants to calculate per-row statistics
#' @param filterStats character vector of statistics to calculate. 
#' 	  Must be one of the character values given by \code{listBuildInFilterStats()}.
#' @return \code{makeFilterStats} returns a \code{\link[SummarizedExperiment]{SummarizedExperiment}} object with the requested filtering statistics will be added to the \code{DataFrame} in the \code{rowData} slot and given names corresponding to the \code{filterStats} values. Warning: the function will overwrite existing columns in \code{rowData} with the same name. Columns in the \code{rowData} slot with different names should not be affected.
#' @inheritParams transformData
#' @rdname reduceFunctions
#'
#' @examples
#' data(simData)
#' listBuiltInFilterStats()
#' scf<-makeFilterStats(simData,filterStats=c("var","mad"))
#' scf
#' scfFiltered<-filterData(scf,filterStats="mad",percentile=10)
#' scfFiltered
#' assay(scfFiltered)[1:10,1:10]
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
  filterStats<-unique(filterStats)
  if(!all(filterStats %in% listBuiltInFilterStats())){
	  stop("Not all of the filterStats given are valid. Must be one of listBuiltInFilterStats().")
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
#' @rdname reduceFunctions
#' @export
setMethod(
  f = "makeFilterStats",
  signature = "matrix",
  definition = function(object,...)
{
	makeFilterStats(SummarizedExperiment(object),...)
}
)

#' @rdname reduceFunctions
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
		else if(length(whCluster)==0) warning("'whichClusterIgnoreUnassigned' does not correspond to a clustering. Ignoring argument.")
		else{
			#give new names to filters to indicate based on clustering.
			newNames<-.makeClusterFilterStats(filterStats,clusterLabels(object)[whCluster])
			whDo<-which(!isFilterStats(object,newNames))
			if(length(whDo)>0){
				whAssigned<-which(clusterMatrix(object)[,whCluster]>0)
				if(length(whAssigned)>0){
					out<-makeFilterStats(object[,whAssigned],filterStats=filterStats[whDo],...)
					whNew<-match(filterStats[whDo],colnames(rowData(out)))
					colnames(rowData(out))[whNew]<-newNames[whDo]	
				}
				else 
					stop("All samples are unassigned for clustering", clusterLabels(object)[whCluster])
				
			}
			else return(object) #already calculated all of the requested values
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
.makeClusterFilterStats<-function(filterStats,clusterName){
	make.names(paste(filterStats,clusterName,sep="_"))
}
#' @rdname reduceFunctions
#' @export
listBuiltInFilterStats<-function(){c('var', 'abscv', 'mad','mean','iqr','median')}
#' @importFrom matrixStats rowVars rowMeans2 rowMads rowMedians rowIQRs
.matchToStats<-SimpleList(
	'var'=matrixStats::rowVars,
	'mad'=matrixStats::rowMads,
	'mean'=matrixStats::rowMeans2,
	'iqr'=matrixStats::rowIQRs,
	'median'=matrixStats::rowMedians)


#' @rdname reduceFunctions
#' @aliases filterData
#' @param cutoff numeric. A value at which to filter the rows (genes) for the test statistic
#' @param percentile numeric. Either a number between 0,1 indicating what percentage of the rows (genes) to keep or an integer value indicated the number of rows (genes) to keep
#' @param absolute whether to take the absolute value of the filter statistic
#' @param keepLarge logical whether to keep rows (genes) with large values of the test statistic or small values of the test statistic. 
#' @details Note that \code{filterData} returns a SingleCellExperiment object. To get the actual data out use either assay or \code{\link{transformData}} if transformed data is desired.
#' @return \code{filterData} returns a SingleCellExperiment object with the rows (genes) removed based on filters
#' @export
#' @importFrom stats quantile
setMethod( "filterData","SummarizedExperiment",
	function(object,filterStats,cutoff,percentile, absolute=FALSE,keepLarge=TRUE){
	stat<-filterStats(object,filterStats)
	if(!is.null(dim(stat))){
		if(NCOL(stat)==1) stat<-stat[,1]
		else stop("Attempting to filtering on more than one filter statistic")
	}
	if(absolute) stat<-abs(stat)
	if(missing(cutoff) & missing(percentile)) stop("must provide one of cutoff or percentile")
	if(!missing(cutoff) & !missing(percentile)) stop("can only provide one of cutoff or percentile")
	if(!missing(cutoff)){
		whKeep<- if(keepLarge) which(stat>cutoff) else which(cutoff > stat)
	}
	if(!missing(percentile)){
		if(0<percentile & percentile <1){
			quantile<- quantile(stat,probs=if(keepLarge) percentile else 1-percentile,na.rm=TRUE)
			whKeep<-if(keepLarge) which(stat>quantile) else which(stat<quantile)
		}
		else{
			if(is.na(percentile) || percentile>=1){
				if(is.na(percentile) || percentile>NROW(object)){
					warning("the number of most features requested after filtering is either missing or larger than the number of features. Will not do any filtering")
					whKeep<-1:NROW(object)
				}
				else whKeep<- order(stat,decreasing=ifelse(keepLarge,TRUE,FALSE),na.last=NA)[1:percentile]
			}
			else stop("Invalid value for percentile. Must be either between 0,1 or a positive integer number to keep")
		}
	}
	object[whKeep,]

})	

#' @rdname reduceFunctions
#' @param reduceMethod character. A method or methods for reducing the size of the data, either by filtering the rows (genes) or by a dimensionality reduction method. Must either be 1) must match the name of a built-in method, in which case if it is not already existing in the object will be passed to \code{\link{makeFilterStats}} or \code{link{makeReducedDims}}, or 2) must match a stored filtering statistic or dimensionality reduction in the object
#' @param typeToShow character (optional). If given, should be one of "filterStats" or "reducedDims" to indicate of the values in the reduceMethod vector, only show those corresponding to "filterStats" or "reducedDims" options. 
#' @return \code{defaultNDims} returns a numeric vector giving the default dimensions the methods in \code{clusterExperiment} will use for reducing the size of the data. If \code{typeToShow} is missing, the resulting vector will be equal to the length of \code{reduceMethod}. Otherwise, it will be a vector with all the unique valid default values for the \code{typeToShow} (note that different dimensionality reduction methods can have different maximal dimensions, so the result may not be of length one in this case).
#' @aliases defaultNDims defaultNDims,SingleCellExperiment-method
setMethod( "defaultNDims","SingleCellExperiment",function(object,reduceMethod,typeToShow){
	nDims<-rep(NA,length(reduceMethod))
	isFilter<-isBuiltInFilterStats(reduceMethod) || isFilterStats(object,reduceMethod)
	isRed<-isReducedDims(object,reduceMethod ) || isBuiltInReducedDims(reduceMethod)
  	if(isFilter)
		nDims[isBuiltInFilterStats(reduceMethod) || isFilterStats(object,reduceMethod)]<-min(500,NROW(object))
	else if(isReducedDims(object,reduceMethod ))
		nDims[isReducedDims(object,reduceMethod)] <- ncolReducedDims(object)[reduceMethod[isReducedDims(object,reduceMethod )]]
	else if(isBuiltInReducedDims(reduceMethod)) nDims[isBuiltInReducedDims(reduceMethod)]<-min(c(50,dim(object)))
	if(!missing(typeToShow)){
		if(typeToShow=="filterStats") nDims<-unique(nDims[isFilter])
		if(typeToShow=="reducedDims") nDims<-unique(nDims[isRed])
		if(length(nDims)==0) nDims<-NA
	}
	return(nDims)
	
})
setMethod( "defaultNDims","matrix",function(object,...){
	return(defaultNDims(SingleCellExperiment(object),...))
})

#' @rdname reduceFunctions
#' @return \code{filterNames} returns a vector of the columns of the rowData that are considered valid filtering statistics. Currently any numeric column in rowData is a valid filtering statistic. 
#' @aliases filterNames
#' @export
setMethod( "filterNames","SummarizedExperiment",function(object){
	checkValid<-TRUE
	if(!checkValid || ncol(rowData(object))==0) colnames(rowData(object))
	else{
		whValid<- sapply(rowData(object),is.numeric)
		# #for now, do not allow NA values for valid filters.
		# whNa<-which(apply(rowData,2,anyNA))
		# whValid[whNA]<-FALSE
		return(names(rowData(object))[whValid])
		
	}
})	