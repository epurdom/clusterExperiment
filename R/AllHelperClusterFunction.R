#' Helper methods for the ClusterFunction class
#' 
#' This is a collection of helper methods for the ClusterExperiment class.
#' @name ClusterFunction-methods
#' @aliases ClusterFunction-methods
#' @param object input to the method, either a \code{ClusterFunction}
#'   class or a character describing a built-in \code{ClusterFunction} object. 
#'   Can also be a \code{list} of \code{ClusterFunction} objects, in which case 
#'   the list must have names for each function.
#' @details Note that when subsetting the data, the dendrogram information and 
#'   the co-clustering matrix are lost.
#' @return \code{requiredArgs}  returns a list of the required args of a
#'   function (via a call to \code{\link{requiredArgs}})
#' @export
setMethod(
  f = "requiredArgs",
  signature = c("character"),
  definition = function(object) {
	  requiredArgs(getBuiltInFunction(object))    
  }
)

#' @rdname ClusterFunction-methods
#' @param genericOnly logical If TRUE, return only the generic required
#'   arguments (i.e. those required by the algorithm type) and not the arguments
#'   specific to that clustering found in the slot \code{requiredArgs}. If FALSE
#'   both sets of arguments are returned.
#' @export
setMethod(
  f = "requiredArgs",
  signature = c("ClusterFunction"),
  definition = function(object,genericOnly=FALSE) {
    algType<-algorithmType(object)
    if(!genericOnly){
      if(!is.na(object@requiredArgs)) reqArgs<-object@requiredArgs
      else reqArgs<-NULL
      if(algType=="01") return(unique(sort(c(reqArgs,.required01Args))))
      if(algType=="K") return(unique(sort(c(reqArgs,.requiredKArgs))))
    }
    else{
      if(algType=="01") return(unique(sort(.required01Args)))
      if(algType=="K") return(unique(sort(.requiredKArgs)))	  	
    }
  }
)
#' @rdname ClusterFunction-methods
#' @export
setMethod(
  f = "requiredArgs",
  signature = c("list"),
  definition = function(object) {
		.checkFunctionList(object)
  	return(sapply(object,requiredArgs))
	}
	
)
#' @rdname ClusterFunction-methods
#' @aliases requiredArgs
#' @export
setMethod(
  f = "requiredArgs",
  signature = c("character"),
  definition = function(object) {
	  requiredArgs(getBuiltInFunction(object))
  }
)
#' @rdname ClusterFunction-methods
#' @export
setMethod(
  f = "requiredArgs",
  signature = c("character"),
  definition = function(object) {
    clObjects<-getBuiltInFunction(object)
    if(length(clObjects)>1) return(lapply(clObjects,requiredArgs))
    else return(requiredArgs(clObjects))
  }
)
#' @rdname ClusterFunction-methods
#' @export
setMethod(
  f = "requiredArgs",
  signature = c("factor"),
  definition = function(object) {
	  requiredArgs(as.character(object))
  }
)

#' @rdname ClusterFunction-methods
#' @aliases algorithmType
#' @return \code{algorithmType} returns a character value giving the type of 
#' clustering function ("01" or "K")
#' @export
setMethod(
  f = "algorithmType",
  signature = c("ClusterFunction"),
  definition = function(object) {
    object@algorithmType
  }
)
#' @rdname ClusterFunction-methods
#' @export
setMethod(
  f = "algorithmType",
  signature = c("character"),
  definition = function(object) {
    clObjects<-getBuiltInFunction(object)
    if(length(clObjects)>1) return(sapply(clObjects,algorithmType))
    else return(algorithmType(clObjects))
  }
)
#' @rdname ClusterFunction-methods
#' @export
setMethod(
  f = "algorithmType",
  signature = c("factor"),
  definition = function(object) {
	  algorithmType(as.character(object))
  }
)
#' @rdname ClusterFunction-methods
#' @export
setMethod(
  f = "algorithmType",
  signature = c("list"),
  definition = function(object) {
		.checkFunctionList(object)
  	return(sapply(object,algorithmType))
	}
	
)
.checkFunctionList<-function(object){
	if(is.null(names(object))){
		stop("if you give a list of ClusterFunction objects, the list must have names")
	}
	if(!all(sapply(object,class)=="ClusterFunction")) stop("if not giving the names of built in functions, must be list of objects of class 'ClusterFunction'")
}

#' @rdname ClusterFunction-methods
#' @aliases inputType
#' @return \code{inputType} returns a character value giving the input 
#' type of the object
#' @export
setMethod(
  f = "inputType",
  signature = c("ClusterFunction"),
  definition = function(object) {
    object@inputType
  }
)
#' @rdname ClusterFunction-methods
#' @export
setMethod(
  f = "inputType",
  signature = c("list"),
  definition = function(object) {
		.checkFunctionList(object)
  	return(sapply(object,inputType))
	}
	
)
#' @rdname ClusterFunction-methods
#' @export
setMethod(
  f = "inputType",
  signature = c("character"),
  definition = function(object) {
    clObjects<-getBuiltInFunction(object)
    if(length(clObjects)>1) return(sapply(clObjects,inputType))
    else return(inputType(clObjects))
  }
)
#' @rdname ClusterFunction-methods
#' @export
setMethod(
  f = "inputType",
  signature = c("factor"),
  definition = function(object) {
	  inputType(as.character(object))
  }
)
