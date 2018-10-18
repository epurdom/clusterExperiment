setGeneric("combineMany", function(x,...) { standardGeneric("combineMany")})
setGeneric("removeUnclustered",function(x,...) {standardGeneric("removeUnclustered")})

#' @name clusterExperiment-deprecated
#' @aliases clusterExperiment-deprecated combineMany combineMany,ANY-method
#' @title Deprecated functions in package \sQuote{clusterExperiment}
#' @description These functions are provided for compatibility with older versions
#'  of \sQuote{clusterExperiment} only, and will be defunct at the next release.
#' @details The following functions are deprecated and will be made defunct; use
#'   the replacement indicated below:
#'  \itemize{
#'    \item{\code{combineMany}: \code{\link{makeConsensus}}}
#'    \item{\code{removeUnclustered}: \code{\link{removeUnassigned}}}
#'  }
#' @export
#' @param x any object
#' @param ... additional arguments
#' @rdname clusterExperiment-deprecated
setMethod(
  f = "combineMany",
  signature = signature(x = "ANY"),
  definition = function(x,...) {
		.Deprecated("makeConsensus") #06/2018 added in BioC 3.8
		makeConsensus(x,...)
	}
	)
setMethod(
  f = "removeUnclustered",
  signature = signature(x = "ANY"),
  definition = function(x,...) {
		.Deprecated("removeUnassigned") #06/2018 added in BioC 3.8
		removeUnassigned(x,...)
  }
)

.depricateArgument<-function(passedArgs,newArgName,oldArgName){
	if(oldArgName %in% names(passedArgs)){
		warning(paste0("the argument name '",oldArgName,"' has been changed to '",newArgName,"'. Please change future calls to reflect this, as the old argument will be depricated in future releases."))
		val<-passedArgs[[oldArgName]]
		passedArgs<-passedArgs[-which(names(passedArgs)==oldArgName)]
		return(list(passedArgs=passedArgs,val=val))
	}
	else return(NULL)
}
