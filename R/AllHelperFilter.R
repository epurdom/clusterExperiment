setGeneric("filterData", function(object,...) { standardGeneric("filterData")})
setGeneric("filterNames", function(object,...) { standardGeneric("filterNames")})
setGeneric("filterStats", function(object,type,...) { standardGeneric("filterStats")})
setGeneric("filterStats<-", function(object, ..., value) standardGeneric("filterStats<-"))

#' Helper methods for the SingleCellFilter class
#'
#' This is a collection of helper methods for the SingleCellFilter class.
#' @name SingleCellFilter-methods
#' @aliases SingleCellFilter
#' @rdname SingleCellFilter-methods
#' @export
setMethod( "filterStats",c("SingleCellFilter","character"),
	function(object,type){
		if(is.null(object@filterStats)) stop("There are no filter statistics saved for this object")
		if(!all(type %in% filterNames(object))){
			type<-type[type%in%filterNames(object)]
			if(length(type)==0) stop("None of the values of 'type' argument are valid filter names ")	
			else	warning(paste("Not all values of '",type, "' are names of a filter statistic held by the object",sep=""))
		} 
		return(object@filterStats[,type] )
	})	
#' @rdname SingleCellFilter-methods
#' @export
setMethod( "filterStats",c("SingleCellFilter","missing"),
	function(object,type){
		return(object@filterStats)
	})	

#' @rdname SingleCellFilter-methods
#' @details Note that the replacement functions never actually completely replace the slot \code{filterStats}; they update existing filters of the same name and add filters with new names to the existing filters. 
#' @export
setReplaceMethod("filterStats", "SingleCellFilter", function(object, type, ...,value) {
	if(missing(type)){
		if(!is.matrix(value) || is.null(colnames(value))){
			stop("If not indicating a type in the replacement, must give matrix of values with column names")
		}
		type<-colnames(value)
	}
    if(length(type)>1){
    	if(!is.matrix(value) || ncol(value)!=length(type)) 
			stop("If replacing several filter Stats at the same time, new value must be a matrix")
		if(!is.null(colnames(value)) && !all(type==colnames(value))) 
			stop("If replacement value has names must match type")
		colnames(value)<-type
    }
	else if(length(type)==1){
		value<-matrix(value,ncol=1)
		colnames(value)<-type
	}
	fs <- filterStats(object) #existing filters
	if(!is.null(fs)){
	    whTypeExist<-which(type %in% colnames(fs))
		if(length(whTypeExist)>0){
			fs[,type[whTypeExist]] <- value[,whTypeExist,drop=FALSE]
		}
	    whTypeNew<-which(!type %in% colnames(fs))
		if(length(whTypeNew)>0){
			fs<-cbind(fs,value[,whTypeNew,drop=FALSE])
		}	
	    object@filterStats <- fs		
	}
	else object@filterStats<-value
    validObject(object)
    return(object)
})

#' @rdname SingleCellFilter-methods
#' @export
setMethod( "filterNames","SingleCellFilter",function(object){colnames(object@filterStats)})		



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
	stat<-if(absolute) abs(filterStats(object,type)) else filterStats(object,type)
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
				if(percentile>NROW(object)){
					warning("the number of most features requested after filtering is larger than the number of features. Will not do any filtering")
					whKeep<-1:NROW(object)
				}
				else whKeep<- order(stat,decreasing=ifelse(keepLarge,FALSE,TRUE))[1:percentile]
			}
			else stop("Invalid value for percentile. Must be either between 0,1 or a positive integer number to keep")
		}
	}
	object[whKeep,]

})		