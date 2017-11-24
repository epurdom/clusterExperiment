#' Helper methods for the SingleCellFilter class
#'
#' This is a collection of helper methods for the SingleCellFilter class.
#' @name SingleCellFilter-methods
#' @aliases SingleCellFilter
#' @rdname SingleCellFilter-methods
#' @export
setGeneric("filterStats", function(object,...) { standardGeneric("filterStats")})
#' @rdname SingleCellFilter-methods
#' @export
setMethod( "filterStats","SingleCellFilter",function(object){object@filterStats})		

#' @rdname SingleCellFilter-methods
#' @export
setGeneric("filterStat", function(object,type,...) { standardGeneric("filterStat")})
#' @rdname SingleCellFilter-methods
#' @export
setMethod( "filterStat",c("SingleCellFilter","character"),
	function(object,type){
		if(is.null(object@filterStats)) stop("There are no filter statistics saved for this object")
		if(!type %in% filterNames(object)) stop(paste("'",type, "' is not the name of a filter statistic held by the object",sep=""))
		return(object@filterStats[,type] )
	})		


#' @rdname SingleCellFilter-methods
#' @export
setGeneric("filterNames", function(object,...) { standardGeneric("filterNames")})
#' @rdname SingleCellFilter-methods
#' @export
setMethod( "filterNames","SingleCellFilter",function(object){colnames(object@filterStats)})		


#' @rdname SingleCellFilter-methods
#' @export
setGeneric("filterData", function(object,...) { standardGeneric("filterData")})

#' @rdname SingleCellFilter-methods
#' @param object a SingleCellFilter object
#' @param type character indicating which filter statistic to use. Must match a single filter statistic in \code{object}
#' @param cutoff numeric. A value at which to filter the rows (genes) for the test statistic
#' @param percentile numeric. Either a number between 0,1 indicating what percentage of the rows (genes) to keep or an integer value indicated the number of rows (genes) to keep
#' @param absolute whether to take the absolute value of the filter statistic
#' @param keepLarge logical whether to keep rows (genes) with large values of the test statistic or small values of the test statistic. 
#' @export
#' @importFrom stats quantile
setMethod( "filterData","SingleCellFilter",
	function(object,type,cutoff,percentile, absolute=FALSE,keepLarge=TRUE){
	stat<-if(absolute) abs(filterStat(object,type)) else filterStat(object,type)
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
			if(percentile>=1){
				whKeep<- order(stat,decreasing=ifelse(keepLarge,FALSE,TRUE))[1:percentile]
			}
			else stop("Invalid value for percentile. Must be either between 0,1 or a positive integer number to keep")
		}
	}
	object[whKeep,]

})		