#' Helper methods for the ClusterFunction class
#'
#' This is a collection of helper methods for the ClusterExperiment class.
#' @name ClusterFunction-methods
#' @aliases ClusterFunction-methods
#' @details Note that when subsetting the data, the dendrogram information and
#' the co-clustering matrix are lost.
#' @export
setMethod(
  f = "requiredArgs",
  signature = c("character"),
  definition = function(object) {
	  requiredArgs(getBuiltInClusterFunction(object))    
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
  	  if(!onlyGeneric){
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
  signature = c("character"),
  definition = function(object) {
	  requiredArgs(getBuiltInClusterFunction(object))
  }
)
#' @rdname ClusterFunction-methods
#' @export
setMethod(
  f = "algorithmType",
  signature = c("ClusterFunction"),
  definition = function(object) {
	  object@algorithmType
	    }
)

