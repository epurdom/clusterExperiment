setGeneric("filterData", function(object,...) { standardGeneric("filterData")})
setGeneric("filterNames", function(object,...) { standardGeneric("filterNames")})
# setGeneric("filterNames<-", function(object,value) { standardGeneric("filterNames<-")})
setGeneric("filterStats", function(object,type,...) { standardGeneric("filterStats")})
setGeneric("filterStats<-", function(object, ..., value) standardGeneric("filterStats<-"))

#' @rdname makeFilterStats
#' @aliases filterNames
#' @export
setMethod( "filterNames","SummarizedExperiment",function(object,checkValid=FALSE){
	if(!checkValid) colnames(rowData(object))
	else{
		whValid<-apply(rowData,2,is.numeric)
		return(colnames(rowData(object))[whValid])
	}
})		

# #' @rdname makeFilterStats
# #' @aliases filterNames<-
# #' @export
# setReplaceMethod("filterNames", "SingleCellFilter", function(object,value) {
# 	    if(length(value)!=NCOL(filterStats(object))) stop("value must be a vector of length equal to NCOL(filterStats(object)):",NCOL(filterStats(object)))
# 	    colnames(object@filterStats) <- value
# 		validObject(object)
# 		return(object)
# 	 }
# )


#' @param object A SummarizedExperiment object
#' @param type a type of filter to retrieve. Should match the filter name.
#' @aliases filterStats,SummarizedExperiment,character-method filterStats
#' @rdname makeFilterStats
#' @export
setMethod( "filterStats",c("SummarizedExperiment","character"),
	function(object,type,checkValid=FALSE){
		if(ncol(rowData(object))==0) stop("There are no filter statistics saved for this object")
		fn<-filterNames(object,checkValid=checkValid)
		if(!all(type %in% fn)){
			type<-type[type%in%fn]
			if(length(type)==0){
				if(checkValid) 
					stop("None of the values of 'type' argument are valid filter names ")
				else stop("None of the values of 'type' argument are names in rowData")
			} 
			else{
				if(checkValid) 
					warning(paste("Not all values of '",type, "' are names of a valid filtering statistic held in rowData",sep=""))
				else warning(paste("Not all values of '",type, "' are names of a column of rowData",sep=""))
			}
		}
		return(rowData[,type] )
	})
#' @rdname makeFilterStats
#' @export
setMethod( "filterStats",c("SummarizedExperiment","missing"),
	function(object,type,checkValid=FALSE){
		return(rowData(object)[,filterNames(object,checkValid=checkValid)])
	})

#' @rdname makeFilterStats
#' @details Note that the replacement functions never actually completely
#'   replace the slot \code{filterStats} unless the replacement value is NULL
#'   They update existing filters of the
#'   same name and add filters with new names to the existing filters.
#' @aliases filterStats<-
#' @export
setReplaceMethod("filterStats", "SummarizedExperiment", function(object, type, ...,value) {
	if(missing(type)){
		# if(is.null(value)){
		# 	object@filterStats<-NULL
		# 	return(object)
		# }
		if(!is.matrix(value) || is.null(colnames(value))){
			stop("If not indicating a type in the replacement, must give matrix of values with column names")
		}
		type<-colnames(value)
	}
    if(length(type)>1){
    	if(!is.matrix(value) || ncol(value)!=length(type))
			stop("If replacing several filtering statistics at the same time, new value must be a matrix")
		if(!is.null(colnames(value)) && !all(type==colnames(value)))
			stop("If replacement value has names must match type")
		colnames(value)<-type
    }
	else if(length(type)==1){
		value<-matrix(value,ncol=1)
		colnames(value)<-type
	}
	fs <- filterStats(object,checkValid=FALSE) #all rowData
	if(ncol(fs)>0){
	    whTypeExist<-which(type %in% colnames(fs))
		if(length(whTypeExist)>0){
			fs[,type[whTypeExist]] <- value[,whTypeExist,drop=FALSE]
		}
	    whTypeNew<-which(!type %in% colnames(fs))
		if(length(whTypeNew)>0){
			fs<-cbind(fs,value[,whTypeNew,drop=FALSE])
		}
	    rowData(object) <- fs
	}
	else rowData(object)<-value
    validObject(object)
    return(object)
})

