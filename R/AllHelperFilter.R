####### Note: These are notexported, so the documentation is commented out.

setGeneric("filterStats", function(object,type,...) { standardGeneric("filterStats")})
setGeneric("filterStats<-", function(object, ..., value) standardGeneric("filterStats<-"))

##Note: have set so that all logicals on the filtering and dimReduce are coded here, so that if change things about filterNames, etc. only effects this code.
setGeneric("isFilterStats", function(object, ...) standardGeneric("isFilterStats"))
setGeneric("anyValidFilterStats", function(object, ...) standardGeneric("anyValidFilterStats"))
setGeneric("isReducedDims", function(object, ...) standardGeneric("isReducedDims"))
setGeneric("anyValidReducedDims", function(object, ...) standardGeneric("anyValidReducedDims"))
setGeneric("isBuiltInFilterStats", function(type, ...) standardGeneric("isBuiltInFilterStats"))
setGeneric("isBuiltInReducedDims", function(type, ...) standardGeneric("isBuiltInReducedDims"))
setGeneric("ncolReducedDims",function(object,...) standardGeneric("ncolReducedDims"))
setGeneric("isPossibleReducedDims",function(object,...) standardGeneric("isPossibleReducedDims"))
setGeneric("isPossibleFilterStats",function(object,...) standardGeneric("isPossibleFilterStats"))

setMethod( "isPossibleFilterStats","SingleCellExperiment",function(object,type){
	type %in% c(listBuiltInFilterStats(),filterNames(object))
})
setMethod( "isPossibleReducedDims","SingleCellExperiment",function(object,type){
	type %in% c(listBuiltInReducedDims(),reducedDimNames(object))
})
setMethod( "ncolReducedDims","SingleCellExperiment",function(object){
	sapply(reducedDims(object),ncol)
})
setMethod( "ncolReducedDims","SingleCellExperiment",function(object){
	sapply(reducedDims(object),ncol)
})
setMethod( "isBuiltInReducedDims","character",function(type){
	type %in% listBuiltInReducedDims()
})
setMethod( "isBuiltInFilterStats","character",function(type){
	type %in% listBuiltInFilterStats()
})		

setMethod( "isReducedDims","SingleCellExperiment",function(object,type){
	type %in% reducedDimNames(object)
})
setMethod( "isFilterStats","SummarizedExperiment",function(object,type){
	type %in% filterNames(object)
})		
setMethod( "anyValidFilterStats","SummarizedExperiment",function(object){
	length(filterNames(object))>0

})		
setMethod( "anyValidReducedDims","SummarizedExperiment",function(object){
	length(reducedDimNames(object))>0

})		
	

# #' @rdname reduceFunctions
# #' @aliases filterNames<-
# #' @export
# setReplaceMethod("filterNames", "SummarizedExperiment", function(object,value) {
# 	fs<-filterStats(object)
# 	if(length(value)!=NCOL(fs)) stop("value must be a vector of length equal to NCOL(filterStats(object)):",NCOL(filterStats(object)))
# 	colnames(fs) <- value
# 	filterStats(object)<-fs
# 	validObject(object)
# 	return(object)
# 	}
# )


# #' @param object A SummarizedExperiment object
# #' @param type a type of filter to retrieve. Should match the filter name.
# #' @aliases filterStats,SummarizedExperiment,character-method filterStats
# #' @rdname reduceFunctions
setMethod( 
  f="filterStats",
  c("SummarizedExperiment","character"),
  function(object,type){
    if(ncol(rowData(object))==0) stop("There are no filter statistics saved for this object")
    fn<-filterNames(object)
    if(!all(type %in% fn)){
      type<-type[type%in%fn]
      if(length(type)==0) stop("None of the values of 'type' argument are valid filter names ")
      else warning(paste("Not all values of '",type, "' are names of a valid filtering statistic held in rowData",sep=""))
    }
    return(rowData(object)[,type,drop=FALSE] )
  })
# #' @rdname reduceFunctions
setMethod( "filterStats",c("SummarizedExperiment","missing"),
	function(object,type){
		return(rowData(object)[,filterNames(object),drop=FALSE])
	})

# #' @rdname reduceFunctions
# #' @details Note that the replacement functions never actually completely
# #'   replace the slot \code{filterStats} unless the replacement value is NULL
# #'   They update existing filters of the
# #'   same name and add filters with new names to the existing filters.
# #' @aliases filterStats<-
#' @importFrom S4Vectors DataFrame
setReplaceMethod("filterStats", "SummarizedExperiment", function(object, type, ...,value) {
  isMatrixLike<-is.matrix(value) || class(value)=="DataFrame"
  if(missing(type)){
    # if(is.null(value)){
    # 	object@filterStats<-NULL
    # 	return(object)
    # }
    if(!isMatrixLike || is.null(colnames(value))){
      stop("If not indicating a type in the replacement, must give matrix of values with column names")
    }
    type<-colnames(value)
  }
  if(length(type)>1){
    if(!isMatrixLike || ncol(value)!=length(type))
      stop("If replacing several filtering statistics at the same time, new value must be a matrix")
    if(!is.null(colnames(value)) && !all(type==colnames(value)))
      stop("If replacement value has names must match type")
    colnames(value)<-type
  }
  else if(length(type)==1){
    if(!isMatrixLike) value<-matrix(value,ncol=1)
    colnames(value)<-type
  }
  fs <- rowData(object) #all rowData
  if(NCOL(fs)>0){
    whTypeExist<-which(type %in% colnames(fs))
    if(length(whTypeExist)>0){
      fs[,type[whTypeExist]] <- S4Vectors::DataFrame(value[,whTypeExist,drop=FALSE])
    }
    whTypeNew<-which(!type %in% colnames(fs))
    if(length(whTypeNew)>0){
      fs<-S4Vectors::DataFrame(fs,value[,whTypeNew,drop=FALSE])
    }
    rowData(object) <- fs
  }
  else rowData(object)<-value
  validObject(object)
  return(object)
})

