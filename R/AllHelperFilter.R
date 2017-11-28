setGeneric("filterData", function(object,...) { standardGeneric("filterData")})
setGeneric("filterNames", function(object,...) { standardGeneric("filterNames")})
setGeneric("filterNames<-", function(object,value) { standardGeneric("filterNames<-")})
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
#' @details Note that the replacement functions never actually completely
#'   replace the slot \code{filterStats} unless the replacement value is NULL
#'   They update existing filters of the
#'   same name and add filters with new names to the existing filters.
#' @export
setReplaceMethod("filterStats", "SingleCellFilter", function(object, type, ...,value) {
	if(missing(type)){
		if(is.null(value)){
			object@filterStats<-NULL
			return(object)
		}
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
#' @export
setReplaceMethod("filterNames", "SingleCellFilter", function(object,value) {
	    if(length(value)!=NCOL(filterStats(object))) stop("value must be a vector of length equal to NCOL(filterStats(object)):",NCOL(filterStats(object)))
	    colnames(object@filterStats) <- value
		validObject(object)
		return(object) 
	 }
)
	
	
#' @rdname SingleCellFilter-methods
#' @export
setMethod("[", c("SingleCellFilter", "ANY", "ANY"), function(x, i, j, ..., drop=TRUE) {
    out<-callNextMethod()
	out@filterStats<-filterStats(out)[i, , drop=FALSE]
    
	return(out)
})

##Code from SingleCellExperiment:::
scat <- function(fmt, vals=character(), exdent=2, ...) {
    vals <- ifelse(nzchar(vals), vals, "''")
    lbls <- paste(S4Vectors:::selectSome(vals), collapse=" ")
    txt <- sprintf(fmt, length(vals), lbls)
    cat(strwrap(txt, exdent=exdent, ...), sep="\n")
}

.scf_show <- function(object) {
    callNextMethod()
    scat("filterNames(%d): %s\n", filterNames(object))
}

#' @rdname SingleCellFilter-methods
#' @export
setMethod("show", "SingleCellFilter", .scf_show)
