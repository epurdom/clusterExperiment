setGeneric("combineMany", function(x,...) { standardGeneric("combineMany")})

#' @name clusterExperiment-deprecated
#' @aliases clusterExperiment-deprecated combineMany combineMany,ANY-method
#' @title Deprecated functions in package \sQuote{clusterExperiment}
#' @description These functions are provided for compatibility with older versions
#'  of \sQuote{clusterExperiment} only, and will be defunct at the next release.
#' @details The following functions are deprecated and will be made defunct; use
#'   the replacement indicated below:
#'  \itemize{
#'    \item{\code{combineMany}: \code{\link{makeConsensus}}}
#'  }
#' @export
#' @param x any object
#' @param ... additional arguments
#' @rdname clusterExperiment-deprecated
setMethod(
  f = "combineMany",
  signature = signature(x = "ANY"),
  definition = function(x,...) {
		.Deprecated("makeConsensus")
		makeConsensus(x,...)
	}
	)
