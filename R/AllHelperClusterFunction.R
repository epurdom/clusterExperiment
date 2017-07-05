#' Helper methods for the ClusterFunction class
#'
#' This is a collection of helper methods for the ClusterExperiment class.
#' @name ClusterFunction-methods
#' @aliases ClusterFunction-methods
#' @param object input to the method, usually either a \code{ClusterFunction} class or a character describing a built-in \code{ClusterFunction} object.
#' @details Note that when subsetting the data, the dendrogram information and
#' the co-clustering matrix are lost.
#' @export
setMethod(
  f = "requiredArgs",
  signature = c("character"),
  definition = function(object) {
	  requiredArgs(getBuiltInFunction(object))    
  }
)

#' @rdname ClusterFunction-methods
#' @param genericOnly logical If TRUE, return only the generic required arguments (i.e. those required by the algorithm type) and not the arguments specific to that clustering found in the slot \code{requiredArgs}. If FALSE both sets of arguments are returned.
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
#' @aliases inputType
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
	  getBuiltInInputType(as.character(object))
  }
)
