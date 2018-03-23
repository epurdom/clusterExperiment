#' @name transformData
#' @title Transform the original data in a ClusterExperiment object
#'
#' @description Provides the transformed data
#'
#' @param object a matrix, SummarizedExperiment, SingleCellExperiment or ClusterExperiment object.
#' @param transFun a transformation function to be applied to the data. If the transformation 
#' applied to the data creates an error or NA values, then the function will throw an error.
#' If object is of class \code{ClusterExperiment}, the stored transformation will be used 
#' and giving this parameter will result in an error.
#' @param isCount if \code{transFun=NULL}, then \code{isCount=TRUE} will determine the
# transformation as defined by \code{function(x){log2(x+1)}}, and \code{isCount=FALSE} 
#' will give a transformation function \code{function(x){x}}. Ignored if \code{transFun=NULL}.
#' If object is of class \code{ClusterExperiment}, the stored transformation will be used 
#' and giving this parameter will result in an error.
#' @param ... Values passed on the the 'matrix' method.
#' @return A DataFrame defined by \code{assay(x)} suitably transformed
#' @details The data matrix defined by \code{assay(x)} is transformed based on
#'   the transformation function either defined in x (in the case of a 
#' \code{ClusterExperiment} object) or by user given values for other classes. 
#'
#'
#' @examples
#' mat <- matrix(data=rnorm(200), ncol=10)
#' mat[1,1] <- -1 #force a negative value
#' labels <- gl(5, 2)
#' cc <- ClusterExperiment(mat, as.numeric(labels), transformation =
#' function(x){x^2}) #define transformation as x^2
#' z<-transformData(cc) 
#' @aliases transformData,matrix-method
#' @export
setMethod(
  f = "transformData",
  signature = "matrixOrHDF5",
  definition = function(object,transFun=NULL,isCount=FALSE) {
	  transFun<-.makeTransFun(transFun=transFun,isCount=isCount)
	  x <- try(transFun(object), silent=TRUE)
	  if(inherits(x, "try-error"))
	    stop("User-supplied `transFun` produces the following error on the input data matrix:\n",x)
	  if(anyNA(x))
	    stop("User-supplied `transFun` produces NA values")
	  return(x)
  }
)
#' @export
#' @rdname transformData
setMethod(
  f = "transformData",
  signature = "ClusterExperiment",
  definition = function(object,...) {
  	if(any(c("transFun","isCount") %in% names(list(...)))) 
  		stop("The internally saved transformation function of a ClusterExperiment object must be used when given as input and setting 'transFun' or 'isCount' for a 'ClusterExperiment' is not allowed.")  
	  return(transformData(assay(object),transFun=transformation(object)))
  }
)
#' @export
#' @rdname transformData
setMethod(
  f = "transformData",
  signature = "SingleCellExperiment",
  definition = function(object,...) {
	  return(transformData(assay(object),...))
  }
)
#' @export
#' @rdname transformData
setMethod(
  f = "transformData",
  signature = "SummarizedExperiment",
  definition = function(object,...) {
	  return(transformData(as(object,"SingleCellExperiment"),...))
  }
)

#small function to uniformally return transformation function from combination of transFun and isCount
.makeTransFun<-function(transFun=NULL,isCount=FALSE){
  if(is.null(transFun)){
    transFun <- if(isCount) function(x){log2(x+1)} else function(x){x}
  }
  return(transFun)
}






	
